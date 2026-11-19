//! Topological DAG / partial-order Pair-HMM MDL + BOCPD alignment.

use crate::bocpd::{
    MIN_POST_PROB, NEG_INF, base_probs, bocpd_gc_profile, joint_pair_bits, logsumexp2, logsumexp3,
};
use crate::mdl::dp::HmmParams;
use crate::seeding::{MASK_RY20, MASK_W10, canonical_spaced_seeds_with_opts};
use std::collections::{HashMap, HashSet};

/// A topologically ordered node in a sequence DAG for partial-order MDL alignment.
#[derive(Debug, Clone)]
pub struct GraphNode {
    /// Segment identifier from GFA (optionally suffixed with orientation).
    pub name: String,
    /// Oriented DNA sequence.
    pub seq: Vec<u8>,
    /// Indices of predecessor nodes in the topological slice.
    pub preds: Vec<usize>,
}

/// Result of aligning a linear query sequence against a sequence DAG.
#[derive(Debug, Clone, PartialEq)]
pub struct GraphAlignment {
    /// Node indices traversed along the optimal Viterbi path.
    pub path_nodes: Vec<usize>,
    /// Kept base-pair matches `(node_idx, node_pos, query_pos)`.
    pub kept_pairs: Vec<(usize, usize, usize)>,
    /// Number of abstained Horizon Bubbles along the path.
    pub bubbles: usize,
    /// Total Viterbi MDL score (bits).
    pub score_bits: f64,
}

/// `(node_pos_opt, query_pos_opt, sub_bits, post_prob)` column record for DAG traceback.
type DagCol = (Option<(usize, usize)>, Option<usize>, f64, f64);

/// Align `seq_y` against a topologically ordered sequence DAG `nodes` using
/// normalized Pair-HMM MDL + BOCPD DP with Forward-Backward posterior gating.
#[must_use]
#[allow(
    clippy::too_many_lines,
    clippy::similar_names,
    clippy::many_single_char_names
)]
pub fn align_seq_to_dag_mdl(
    nodes: &[GraphNode],
    seq_y: &[u8],
    params: &HmmParams,
    min_block_bits: f64,
) -> GraphAlignment {
    let m = seq_y.len();
    if nodes.is_empty() || m == 0 {
        return GraphAlignment {
            path_nodes: Vec::new(),
            kept_pairs: Vec::new(),
            bubbles: usize::from(m >= 6),
            score_bits: 0.0,
        };
    }

    let num_nodes = nodes.len();
    let mut succs: Vec<Vec<usize>> = vec![Vec::new(); num_nodes];
    for (u, node) in nodes.iter().enumerate() {
        for &p in &node.preds {
            if p < u {
                succs[p].push(u);
            }
        }
    }

    let gc_y = bocpd_gc_profile(seq_y, 1.0 / 120.0, 80);
    let pi_y: Vec<[f64; 4]> = gc_y.into_iter().map(base_probs).collect();
    let pi_nodes: Vec<Vec<[f64; 4]>> = nodes
        .iter()
        .map(|nd| {
            bocpd_gc_profile(&nd.seq, 1.0 / 120.0, 80)
                .into_iter()
                .map(base_probs)
                .collect()
        })
        .collect();

    let id_prior = params.identity_prior.clamp(0.40, 0.98);
    let go_bits = params.gap_open_bits.clamp(-8.0, -1.5);
    let ge_bits = params.gap_extend_bits.clamp(-4.0, -0.2);

    let ln2 = 2.0_f64.ln();
    let p_open = 2.0_f64.powf(go_bits);
    let p_ext = 2.0_f64.powf(ge_bits);
    let mm_bits = (1.0 - 2.0 * p_open).max(1e-6).log2();
    let im_bits = (1.0 - p_ext).max(1e-6).log2();
    let go = go_bits * ln2;
    let ge = ge_bits * ln2;
    let mm_ln = mm_bits * ln2;
    let im_ln = im_bits * ln2;

    let cols = m + 1;
    let mut m_dp: Vec<Vec<f64>> = Vec::with_capacity(num_nodes);
    let mut ix_dp: Vec<Vec<f64>> = Vec::with_capacity(num_nodes);
    let mut iy_dp: Vec<Vec<f64>> = Vec::with_capacity(num_nodes);
    let mut fm: Vec<Vec<f64>> = Vec::with_capacity(num_nodes);
    let mut fix: Vec<Vec<f64>> = Vec::with_capacity(num_nodes);
    let mut fiy: Vec<Vec<f64>> = Vec::with_capacity(num_nodes);
    let mut pred_bt: Vec<Vec<[Option<usize>; 3]>> = Vec::with_capacity(num_nodes);

    for nd in nodes {
        let rows = nd.seq.len() + 1;
        m_dp.push(vec![NEG_INF; rows * cols]);
        ix_dp.push(vec![NEG_INF; rows * cols]);
        iy_dp.push(vec![NEG_INF; rows * cols]);
        fm.push(vec![NEG_INF; rows * cols]);
        fix.push(vec![NEG_INF; rows * cols]);
        fiy.push(vec![NEG_INF; rows * cols]);
        pred_bt.push(vec![[None, None, None]; cols]);
    }

    for u in 0..num_nodes {
        let nu = nodes[u].seq.len();
        if nodes[u].preds.is_empty() {
            m_dp[u][0] = 0.0;
            fm[u][0] = 0.0;
        } else {
            for j in 0..=m {
                let mut best_m = NEG_INF;
                let mut best_ix = NEG_INF;
                let mut best_iy = NEG_INF;
                let mut pm = None;
                let mut pix = None;
                let mut piy = None;
                let mut sum_fm = NEG_INF;
                let mut sum_fix = NEG_INF;
                let mut sum_fiy = NEG_INF;

                for &p in &nodes[u].preds {
                    if p >= u {
                        continue;
                    }
                    let pend = nodes[p].seq.len() * cols + j;
                    if m_dp[p][pend] > best_m {
                        best_m = m_dp[p][pend];
                        pm = Some(p);
                    }
                    if ix_dp[p][pend] > best_ix {
                        best_ix = ix_dp[p][pend];
                        pix = Some(p);
                    }
                    if iy_dp[p][pend] > best_iy {
                        best_iy = iy_dp[p][pend];
                        piy = Some(p);
                    }
                    sum_fm = logsumexp2(sum_fm, fm[p][pend]);
                    sum_fix = logsumexp2(sum_fix, fix[p][pend]);
                    sum_fiy = logsumexp2(sum_fiy, fiy[p][pend]);
                }
                m_dp[u][j] = best_m;
                ix_dp[u][j] = best_ix;
                iy_dp[u][j] = best_iy;
                fm[u][j] = sum_fm;
                fix[u][j] = sum_fix;
                fiy[u][j] = sum_fiy;
                pred_bt[u][j] = [pm, pix, piy];
            }
        }

        for j in 1..=m {
            let v_iy = (m_dp[u][j - 1] + go_bits).max(iy_dp[u][j - 1] + ge_bits);
            if v_iy > iy_dp[u][j] {
                iy_dp[u][j] = v_iy;
                pred_bt[u][j][2] = None;
            }
            fiy[u][j] = logsumexp2(fiy[u][j], logsumexp2(fm[u][j - 1] + go, fiy[u][j - 1] + ge));
        }

        for i in 1..=nu {
            let row = i * cols;
            let prev_row = (i - 1) * cols;
            for j in 0..=m {
                let idx = row + j;
                let up = prev_row + j;
                let v_ix = (m_dp[u][up] + go_bits).max(ix_dp[u][up] + ge_bits);
                ix_dp[u][idx] = v_ix;
                fix[u][idx] = logsumexp2(fm[u][up] + go, fix[u][up] + ge);

                if j > 0 {
                    let left = row + j - 1;
                    let v_iy = (m_dp[u][left] + go_bits).max(iy_dp[u][left] + ge_bits);
                    iy_dp[u][idx] = v_iy;
                    fiy[u][idx] = logsumexp2(fm[u][left] + go, fiy[u][left] + ge);

                    let diag = prev_row + j - 1;
                    let sb = joint_pair_bits(
                        nodes[u].seq[i - 1],
                        seq_y[j - 1],
                        pi_nodes[u][i - 1],
                        pi_y[j - 1],
                        id_prior,
                    );
                    let prev_m = (m_dp[u][diag] + mm_bits)
                        .max(ix_dp[u][diag] + im_bits)
                        .max(iy_dp[u][diag] + im_bits);
                    m_dp[u][idx] = prev_m + sb;
                    let prev_fm = logsumexp3(
                        fm[u][diag] + mm_ln,
                        fix[u][diag] + im_ln,
                        fiy[u][diag] + im_ln,
                    );
                    fm[u][idx] = prev_fm + sb * ln2;
                }
            }
        }
    }

    let mut best_end_node = 0usize;
    let mut best_end_score = NEG_INF;
    let mut z_fwd = NEG_INF;
    for u in 0..num_nodes {
        if succs[u].is_empty() {
            let end_idx = nodes[u].seq.len() * cols + m;
            let sc = m_dp[u][end_idx].max(ix_dp[u][end_idx]).max(iy_dp[u][end_idx]);
            if sc > best_end_score {
                best_end_score = sc;
                best_end_node = u;
            }
            z_fwd = logsumexp2(
                z_fwd,
                logsumexp3(fm[u][end_idx], fix[u][end_idx], fiy[u][end_idx]),
            );
        }
    }

    if z_fwd <= NEG_INF * 0.5 || best_end_score <= NEG_INF * 0.5 {
        return GraphAlignment {
            path_nodes: Vec::new(),
            kept_pairs: Vec::new(),
            bubbles: 1,
            score_bits: best_end_score,
        };
    }

    let mut bm: Vec<Vec<f64>> = Vec::with_capacity(num_nodes);
    let mut bix: Vec<Vec<f64>> = Vec::with_capacity(num_nodes);
    let mut biy: Vec<Vec<f64>> = Vec::with_capacity(num_nodes);
    for nd in nodes {
        let rows = nd.seq.len() + 1;
        bm.push(vec![NEG_INF; rows * cols]);
        bix.push(vec![NEG_INF; rows * cols]);
        biy.push(vec![NEG_INF; rows * cols]);
    }

    for u in (0..num_nodes).rev() {
        let nu = nodes[u].seq.len();
        if succs[u].is_empty() {
            let end_idx = nu * cols + m;
            bm[u][end_idx] = 0.0;
            bix[u][end_idx] = 0.0;
            biy[u][end_idx] = 0.0;
        } else {
            let row = nu * cols;
            for j in 0..=m {
                let mut s_bm = NEG_INF;
                let mut s_bix = NEG_INF;
                let mut s_biy = NEG_INF;
                for &s in &succs[u] {
                    s_bm = logsumexp2(s_bm, bm[s][j]);
                    s_bix = logsumexp2(s_bix, bix[s][j]);
                    s_biy = logsumexp2(s_biy, biy[s][j]);
                }
                bm[u][row + j] = s_bm;
                bix[u][row + j] = s_bix;
                biy[u][row + j] = s_biy;
            }
        }

        for i in (0..=nu).rev() {
            for j in (0..=m).rev() {
                if i == nu && j == m {
                    continue;
                }
                let idx = i * cols + j;
                let m_next = if i < nu && j < m {
                    let sb = joint_pair_bits(
                        nodes[u].seq[i],
                        seq_y[j],
                        pi_nodes[u][i],
                        pi_y[j],
                        id_prior,
                    );
                    bm[u][(i + 1) * cols + j + 1] + sb * ln2
                } else {
                    NEG_INF
                };
                let ix_open = if i < nu {
                    bix[u][(i + 1) * cols + j] + go
                } else {
                    NEG_INF
                };
                let ix_ext = if i < nu {
                    bix[u][(i + 1) * cols + j] + ge
                } else {
                    NEG_INF
                };
                let iy_open = if j < m {
                    biy[u][i * cols + j + 1] + go
                } else {
                    NEG_INF
                };
                let iy_ext = if j < m {
                    biy[u][i * cols + j + 1] + ge
                } else {
                    NEG_INF
                };
                bm[u][idx] = logsumexp2(bm[u][idx], logsumexp3(m_next + mm_ln, ix_open, iy_open));
                bix[u][idx] = logsumexp2(bix[u][idx], logsumexp2(m_next + im_ln, ix_ext));
                biy[u][idx] = logsumexp2(biy[u][idx], logsumexp2(m_next + im_ln, iy_ext));
            }
        }
    }

    let mut raw_cols: Vec<DagCol> = Vec::new();
    let mut path_nodes = vec![best_end_node];
    let mut u = best_end_node;
    let mut i = nodes[u].seq.len();
    let mut j = m;
    let end_idx = i * cols + j;
    let mut state = if m_dp[u][end_idx] >= ix_dp[u][end_idx]
        && m_dp[u][end_idx] >= iy_dp[u][end_idx]
    {
        0usize
    } else if ix_dp[u][end_idx] >= iy_dp[u][end_idx] {
        1usize
    } else {
        2usize
    };

    loop {
        if i == 0 {
            if let Some(p) = pred_bt[u][j][state] {
                u = p;
                i = nodes[u].seq.len();
                path_nodes.push(u);
                continue;
            }
            if j == 0 {
                break;
            }
            raw_cols.push((None, Some(j - 1), ge_bits, 0.0));
            j -= 1;
            let left_m = m_dp[u][j] + go_bits;
            let left_iy = iy_dp[u][j] + ge_bits;
            state = if left_m >= left_iy { 0 } else { 2 };
            continue;
        }

        let cur_idx = i * cols + j;
        if state == 0 && j > 0 {
            let sb = joint_pair_bits(
                nodes[u].seq[i - 1],
                seq_y[j - 1],
                pi_nodes[u][i - 1],
                pi_y[j - 1],
                id_prior,
            );
            let post_p = (fm[u][cur_idx] + bm[u][cur_idx] - z_fwd).min(0.0).exp();
            raw_cols.push((Some((u, i - 1)), Some(j - 1), sb, post_p));
            let diag = (i - 1) * cols + (j - 1);
            let vm = m_dp[u][diag] + mm_bits;
            let vix = ix_dp[u][diag] + im_bits;
            let viy = iy_dp[u][diag] + im_bits;
            state = if vm >= vix && vm >= viy {
                0
            } else if vix >= viy {
                1
            } else {
                2
            };
            i -= 1;
            j -= 1;
        } else if state == 1 {
            raw_cols.push((Some((u, i - 1)), None, ge_bits, 0.0));
            let up = (i - 1) * cols + j;
            let vm = m_dp[u][up] + go_bits;
            let vix = ix_dp[u][up] + ge_bits;
            state = usize::from(vm < vix);
            i -= 1;
        } else if j > 0 {
            raw_cols.push((None, Some(j - 1), ge_bits, 0.0));
            let left = i * cols + (j - 1);
            let vm = m_dp[u][left] + go_bits;
            let viy = iy_dp[u][left] + ge_bits;
            state = if vm >= viy { 0 } else { 2 };
            j -= 1;
        } else {
            break;
        }
    }

    raw_cols.reverse();
    path_nodes.reverse();

    let k = raw_cols.len();
    let mut abstain = vec![false; k];
    let mut idx = 0;
    while idx < k {
        let (cx, cy, _, pp) = raw_cols[idx];
        if cx.is_some() && cy.is_some() && pp >= MIN_POST_PROB {
            let start = idx;
            let mut block_bits = 0.0;
            while idx < k
                && ((raw_cols[idx].0.is_some()
                    && raw_cols[idx].1.is_some()
                    && raw_cols[idx].3 >= MIN_POST_PROB)
                    || (idx + 1 < k
                        && raw_cols[idx].2 > -2.5
                        && raw_cols[idx + 1].3 >= MIN_POST_PROB))
            {
                block_bits += raw_cols[idx].2;
                idx += 1;
            }
            if block_bits < min_block_bits {
                abstain[start..idx].fill(true);
            }
        } else {
            abstain[idx] = true;
            idx += 1;
        }
    }

    let mut kept_pairs = Vec::with_capacity(k);
    let mut bubbles = 0usize;
    idx = 0;
    while idx < k {
        if abstain[idx] {
            let start = idx;
            while idx < k && abstain[idx] {
                idx += 1;
            }
            if idx - start >= 6 {
                bubbles += 1;
            }
        } else {
            if let (Some((nid, np)), Some(yp)) = (raw_cols[idx].0, raw_cols[idx].1) {
                kept_pairs.push((nid, np, yp));
            }
            idx += 1;
        }
    }

    GraphAlignment {
        path_nodes,
        kept_pairs,
        bubbles,
        score_bits: best_end_score,
    }
}

/// Seeded DAG alignment: uses `MASK_W10` + `MASK_RY20` spaced seeds to isolate the
/// connected/reachable DAG subgraph matching `seq_y`, then runs [`align_seq_to_dag_mdl`].
/// Maps returned node indices back to the original `nodes` slice.
#[must_use]
pub fn align_seq_to_dag_seeded(
    nodes: &[GraphNode],
    seq_y: &[u8],
    params: &HmmParams,
    min_block_bits: f64,
) -> GraphAlignment {
    if nodes.len() <= 8 {
        return align_seq_to_dag_mdl(nodes, seq_y, params, min_block_bits);
    }

    let mut seed_hits: HashSet<usize> = HashSet::new();
    for (mask, ry) in [(MASK_W10, false), (MASK_RY20, true)] {
        if let Ok(q_seeds) = canonical_spaced_seeds_with_opts(seq_y, mask, 8, false, ry)
            && !q_seeds.is_empty()
        {
            for (nid, nd) in nodes.iter().enumerate() {
                if let Ok(n_seeds) = canonical_spaced_seeds_with_opts(&nd.seq, mask, 8, false, ry)
                    && n_seeds.keys().any(|k| q_seeds.contains_key(k))
                {
                    seed_hits.insert(nid);
                }
            }
        }
    }

    if seed_hits.is_empty() {
        return align_seq_to_dag_mdl(nodes, seq_y, params, min_block_bits);
    }

    // Expand backward and forward reachability from seeded nodes to preserve full bubble paths.
    let mut succs: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
    for (u, nd) in nodes.iter().enumerate() {
        for &p in &nd.preds {
            if p < u {
                succs[p].push(u);
            }
        }
    }

    let mut active: Vec<bool> = vec![false; nodes.len()];
    for &s in &seed_hits {
        active[s] = true;
    }
    for u in (0..nodes.len()).rev() {
        if active[u] {
            for &p in &nodes[u].preds {
                if p < u {
                    active[p] = true;
                }
            }
        }
    }
    for u in 0..nodes.len() {
        if active[u] {
            for &s in &succs[u] {
                active[s] = true;
            }
        }
    }

    let mut old_to_new: HashMap<usize, usize> = HashMap::new();
    let mut new_to_old: Vec<usize> = Vec::new();
    let mut sub_nodes: Vec<GraphNode> = Vec::new();
    for (old_id, nd) in nodes.iter().enumerate() {
        if !active[old_id] {
            continue;
        }
        let new_id = sub_nodes.len();
        old_to_new.insert(old_id, new_id);
        new_to_old.push(old_id);
        let mapped_preds = nd
            .preds
            .iter()
            .filter_map(|p| old_to_new.get(p).copied())
            .collect();
        sub_nodes.push(GraphNode {
            name: nd.name.clone(),
            seq: nd.seq.clone(),
            preds: mapped_preds,
        });
    }

    let mut aln = align_seq_to_dag_mdl(&sub_nodes, seq_y, params, min_block_bits);
    for nid in &mut aln.path_nodes {
        *nid = new_to_old[*nid];
    }
    for (nid, _, _) in &mut aln.kept_pairs {
        *nid = new_to_old[*nid];
    }
    aln
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dag_mdl_aligns_through_correct_bubble_branch() {
        let nodes = vec![
            GraphNode {
                name: "n0".into(),
                seq: b"ACGTACGTACGTACGTACGT".to_vec(),
                preds: vec![],
            },
            GraphNode {
                name: "n1".into(),
                seq: b"GGGGGGGGGGGGGGGGGGGG".to_vec(),
                preds: vec![0],
            },
            GraphNode {
                name: "n2".into(),
                seq: b"ATATATATATATATATATAT".to_vec(),
                preds: vec![0],
            },
            GraphNode {
                name: "n3".into(),
                seq: b"TGCATGCATGCATGCATGCA".to_vec(),
                preds: vec![1, 2],
            },
        ];
        let query = b"ACGTACGTACGTACGTACGTATATATATATATATATATATTGCATGCATGCATGCATGCA";
        let aln = align_seq_to_dag_mdl(&nodes, query, &HmmParams::default(), 10.0);
        assert_eq!(aln.path_nodes, vec![0, 2, 3]);
        assert_eq!(aln.kept_pairs.len(), 60);
        assert_eq!(aln.bubbles, 0);
        assert!(aln.score_bits > 60.0);
    }

    #[test]
    fn seeded_dag_aligns_ignoring_disconnected_decoy_component() {
        let mut nodes = Vec::new();
        // 6 unrelated decoy nodes
        for idx in 0..6 {
            nodes.push(GraphNode {
                name: format!("decoy_{idx}"),
                seq: b"CCCCCCCCCCCCCCCCCCCCCCCCC".to_vec(),
                preds: if idx == 0 { vec![] } else { vec![idx - 1] },
            });
        }
        // Real 4-node diamond component (indices 6..10)
        nodes.push(GraphNode {
            name: "r0".into(),
            seq: b"ACTGACTGACTGACTGACTGACTGACTG".to_vec(),
            preds: vec![],
        });
        nodes.push(GraphNode {
            name: "r1".into(),
            seq: b"TTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec(),
            preds: vec![6],
        });
        nodes.push(GraphNode {
            name: "r2".into(),
            seq: b"GATTACAGATTACAGATTACAGATTACA".to_vec(),
            preds: vec![6],
        });
        nodes.push(GraphNode {
            name: "r3".into(),
            seq: b"CATGCATGCATGCATGCATGCATGCATG".to_vec(),
            preds: vec![7, 8],
        });
        let query = b"ACTGACTGACTGACTGACTGACTGACTGGATTACAGATTACAGATTACAGATTACACATGCATGCATGCATGCATGCATGCATG";
        let aln = align_seq_to_dag_seeded(&nodes, query, &HmmParams::default(), 10.0);
        assert_eq!(aln.path_nodes, vec![6, 8, 9]);
        assert_eq!(aln.kept_pairs.len(), query.len());
    }
}
