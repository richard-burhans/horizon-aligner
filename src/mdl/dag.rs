//! Topological DAG / partial-order Pair-HMM MDL + BOCPD alignment.

use crate::bocpd::{
    MIN_POST_PROB, NEG_INF, base_probs, bocpd_gc_profile, joint_pair_bits, logsumexp2, logsumexp3,
};
use crate::mdl::dp::HmmParams;
use crate::seeding::{FlatSeedIndex, MASK_RY20, MASK_W10, sdust_mask};
use std::collections::{HashMap, HashSet};

/// Adaptive banded per-node DP grid (`abPOA`, Gao et al. 2021).
#[derive(Clone)]
struct DagBandGrid {
    bounds: Vec<(usize, usize, usize)>, // (j_lo, j_hi, offset) per row 0..=nu
    data: Vec<f64>,
}

impl DagBandGrid {
    fn new(bounds: Vec<(usize, usize, usize)>, total_cells: usize, init: f64) -> Self {
        Self {
            bounds,
            data: vec![init; total_cells],
        }
    }

    #[inline]
    fn row_range(&self, i: usize) -> (usize, usize) {
        let (lo, hi, _) = self.bounds[i];
        (lo, hi)
    }

    #[inline]
    fn get(&self, i: usize, j: usize) -> f64 {
        let (lo, hi, off) = self.bounds[i];
        if j < lo || j > hi {
            NEG_INF
        } else {
            self.data[off + (j - lo)]
        }
    }

    #[inline]
    fn set(&mut self, i: usize, j: usize, val: f64) {
        let (lo, hi, off) = self.bounds[i];
        if j >= lo && j <= hi {
            self.data[off + (j - lo)] = val;
        }
    }
}

fn build_node_bounds(
    nu: usize,
    s_min: usize,
    s_max: usize,
    total_max: usize,
    m: usize,
) -> (Vec<(usize, usize, usize)>, usize) {
    let mut bounds = Vec::with_capacity(nu + 1);
    let mut off = 0usize;
    let use_band = m > 128 && total_max > 0;
    let band_w = (m / 8).clamp(64, 512);
    for i in 0..=nu {
        let (j_lo, j_hi) = if use_band {
            let lo_est = ((s_min + i) * m) / total_max;
            let hi_est = ((s_max + i) * m) / total_max;
            (lo_est.saturating_sub(band_w), (hi_est + band_w).min(m))
        } else {
            (0, m)
        };
        let len = j_hi.saturating_sub(j_lo) + 1;
        bounds.push((j_lo, j_hi, off));
        off += len;
    }
    (bounds, off)
}

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

    let mut start_min = vec![0usize; num_nodes];
    let mut start_max = vec![0usize; num_nodes];
    for u in 0..num_nodes {
        let mut mn = usize::MAX;
        let mut mx = 0usize;
        let mut has_pred = false;
        for &p in &nodes[u].preds {
            if p < u {
                has_pred = true;
                mn = mn.min(start_min[p] + nodes[p].seq.len());
                mx = mx.max(start_max[p] + nodes[p].seq.len());
            }
        }
        if has_pred {
            start_min[u] = mn;
            start_max[u] = mx;
        }
    }
    let total_max = (0..num_nodes)
        .filter(|&u| succs[u].is_empty())
        .map(|u| start_max[u] + nodes[u].seq.len())
        .max()
        .unwrap_or(0);

    let mut fm: Vec<DagBandGrid> = Vec::with_capacity(num_nodes);
    let mut fix: Vec<DagBandGrid> = Vec::with_capacity(num_nodes);
    let mut fiy: Vec<DagBandGrid> = Vec::with_capacity(num_nodes);

    for (u, nd) in nodes.iter().enumerate() {
        let (bounds, cells) =
            build_node_bounds(nd.seq.len(), start_min[u], start_max[u], total_max, m);
        fm.push(DagBandGrid::new(bounds.clone(), cells, NEG_INF));
        fix.push(DagBandGrid::new(bounds.clone(), cells, NEG_INF));
        fiy.push(DagBandGrid::new(bounds, cells, NEG_INF));
    }

    for u in 0..num_nodes {
        let nu = nodes[u].seq.len();
        let (j0_lo, j0_hi) = fm[u].row_range(0);
        if nodes[u].preds.is_empty() {
            fm[u].set(0, 0, 0.0);
        } else {
            for j in j0_lo..=j0_hi {
                let mut sum_fm = NEG_INF;
                let mut sum_fix = NEG_INF;
                let mut sum_fiy = NEG_INF;
                for &p in &nodes[u].preds {
                    if p >= u {
                        continue;
                    }
                    let pend = nodes[p].seq.len();
                    sum_fm = logsumexp2(sum_fm, fm[p].get(pend, j));
                    sum_fix = logsumexp2(sum_fix, fix[p].get(pend, j));
                    sum_fiy = logsumexp2(sum_fiy, fiy[p].get(pend, j));
                }
                fm[u].set(0, j, sum_fm);
                fix[u].set(0, j, sum_fix);
                fiy[u].set(0, j, sum_fiy);
            }
        }

        if nodes[u].preds.is_empty() {
            for j in (j0_lo.max(1))..=j0_hi {
                let cur = fiy[u].get(0, j);
                let from_left = logsumexp2(fm[u].get(0, j - 1) + go, fiy[u].get(0, j - 1) + ge);
                fiy[u].set(0, j, logsumexp2(cur, from_left));
            }
        }

        for i in 1..=nu {
            let (j_lo, j_hi) = fm[u].row_range(i);
            for j in j_lo..=j_hi {
                let v_fix = logsumexp2(fm[u].get(i - 1, j) + go, fix[u].get(i - 1, j) + ge);
                fix[u].set(i, j, v_fix);

                if j > 0 {
                    let v_fiy = logsumexp2(fm[u].get(i, j - 1) + go, fiy[u].get(i, j - 1) + ge);
                    fiy[u].set(i, j, v_fiy);

                    let sb = joint_pair_bits(
                        nodes[u].seq[i - 1],
                        seq_y[j - 1],
                        pi_nodes[u][i - 1],
                        pi_y[j - 1],
                        id_prior,
                    );
                    let prev_fm = logsumexp3(
                        fm[u].get(i - 1, j - 1) + mm_ln,
                        fix[u].get(i - 1, j - 1) + im_ln,
                        fiy[u].get(i - 1, j - 1) + im_ln,
                    );
                    fm[u].set(i, j, prev_fm + sb * ln2);
                }
            }
        }
    }

    let mut z_fwd = NEG_INF;
    for u in 0..num_nodes {
        if succs[u].is_empty() {
            let nu = nodes[u].seq.len();
            z_fwd = logsumexp2(
                z_fwd,
                logsumexp3(fm[u].get(nu, m), fix[u].get(nu, m), fiy[u].get(nu, m)),
            );
        }
    }
    let score_bits = z_fwd / ln2;

    if z_fwd <= NEG_INF * 0.5 {
        return GraphAlignment {
            path_nodes: Vec::new(),
            kept_pairs: Vec::new(),
            bubbles: 1,
            score_bits,
        };
    }

    let mut bm: Vec<DagBandGrid> = fm
        .iter()
        .map(|g| DagBandGrid::new(g.bounds.clone(), g.data.len(), NEG_INF))
        .collect();
    let mut bix: Vec<DagBandGrid> = fm
        .iter()
        .map(|g| DagBandGrid::new(g.bounds.clone(), g.data.len(), NEG_INF))
        .collect();
    let mut biy: Vec<DagBandGrid> = fm
        .iter()
        .map(|g| DagBandGrid::new(g.bounds.clone(), g.data.len(), NEG_INF))
        .collect();

    for u in (0..num_nodes).rev() {
        let nu = nodes[u].seq.len();
        if succs[u].is_empty() {
            bm[u].set(nu, m, 0.0);
            bix[u].set(nu, m, 0.0);
            biy[u].set(nu, m, 0.0);
        } else {
            let (j_lo, j_hi) = bm[u].row_range(nu);
            for j in j_lo..=j_hi {
                let mut s_bm = NEG_INF;
                let mut s_bix = NEG_INF;
                let mut s_biy = NEG_INF;
                for &s in &succs[u] {
                    s_bm = logsumexp2(s_bm, bm[s].get(0, j));
                    s_bix = logsumexp2(s_bix, bix[s].get(0, j));
                    s_biy = logsumexp2(s_biy, biy[s].get(0, j));
                }
                bm[u].set(nu, j, s_bm);
                bix[u].set(nu, j, s_bix);
                biy[u].set(nu, j, s_biy);
            }
        }

        for i in (0..=nu).rev() {
            let (j_lo, j_hi) = bm[u].row_range(i);
            for j in (j_lo..=j_hi).rev() {
                if i == nu && j == m {
                    continue;
                }
                let m_next = if i < nu && j < m {
                    let sb = joint_pair_bits(
                        nodes[u].seq[i],
                        seq_y[j],
                        pi_nodes[u][i],
                        pi_y[j],
                        id_prior,
                    );
                    bm[u].get(i + 1, j + 1) + sb * ln2
                } else {
                    NEG_INF
                };
                let ix_open = if i < nu {
                    bix[u].get(i + 1, j) + go
                } else {
                    NEG_INF
                };
                let ix_ext = if i < nu {
                    bix[u].get(i + 1, j) + ge
                } else {
                    NEG_INF
                };
                let allow_iy = i > 0 || nodes[u].preds.is_empty();
                let iy_open = if allow_iy && j < m {
                    biy[u].get(i, j + 1) + go
                } else {
                    NEG_INF
                };
                let iy_ext = if allow_iy && j < m {
                    biy[u].get(i, j + 1) + ge
                } else {
                    NEG_INF
                };
                let cur_bm = bm[u].get(i, j);
                let cur_bix = bix[u].get(i, j);
                let cur_biy = biy[u].get(i, j);
                bm[u].set(
                    i,
                    j,
                    logsumexp2(cur_bm, logsumexp3(m_next + mm_ln, ix_open, iy_open)),
                );
                bix[u].set(
                    i,
                    j,
                    logsumexp2(cur_bix, logsumexp2(m_next + im_ln, ix_ext)),
                );
                biy[u].set(
                    i,
                    j,
                    logsumexp2(cur_biy, logsumexp2(m_next + im_ln, iy_ext)),
                );
            }
        }
    }

    // Pure Forward-Backward MEA decoding over the DAG:
    // reuse `fm[u]` for MEA cumulative score and `bm[u]` for posterior match probability.
    for u in 0..num_nodes {
        let nu = nodes[u].seq.len();
        let (j0_lo, j0_hi) = fm[u].row_range(0);
        if nodes[u].preds.is_empty() {
            for j in j0_lo..=j0_hi {
                fm[u].set(0, j, 0.0);
            }
        } else {
            for j in j0_lo..=j0_hi {
                let mut best_p = NEG_INF;
                for &p in &nodes[u].preds {
                    if p < u {
                        best_p = best_p.max(fm[p].get(nodes[p].seq.len(), j));
                    }
                }
                fm[u].set(0, j, best_p);
            }
            for j in (j0_lo.max(1))..=j0_hi {
                let left = fm[u].get(0, j - 1);
                if left > fm[u].get(0, j) {
                    fm[u].set(0, j, left);
                }
            }
        }

        for i in 1..=nu {
            let (j_lo, j_hi) = fm[u].row_range(i);
            for j in j_lo..=j_hi {
                let post_m = if j > 0 {
                    (fm[u].get(i, j) + bm[u].get(i, j) - z_fwd).min(0.0).exp()
                } else {
                    0.0
                };
                if j > 0 {
                    bm[u].set(i, j, post_m);
                }
                let d = if j > 0 {
                    fm[u].get(i - 1, j - 1) + (post_m - 0.25)
                } else {
                    NEG_INF
                };
                let up = fm[u].get(i - 1, j);
                let left = if j > 0 { fm[u].get(i, j - 1) } else { NEG_INF };
                fm[u].set(i, j, d.max(up).max(left));
            }
        }
    }

    let mut best_end_node = 0usize;
    let mut best_mea = NEG_INF;
    for u in 0..num_nodes {
        if succs[u].is_empty() {
            let sc = fm[u].get(nodes[u].seq.len(), m);
            if sc > best_mea {
                best_mea = sc;
                best_end_node = u;
            }
        }
    }

    if best_mea <= NEG_INF * 0.5 {
        return GraphAlignment {
            path_nodes: Vec::new(),
            kept_pairs: Vec::new(),
            bubbles: 1,
            score_bits,
        };
    }

    let mut raw_cols: Vec<DagCol> = Vec::new();
    let mut path_nodes = vec![best_end_node];
    let mut u = best_end_node;
    let mut i = nodes[u].seq.len();
    let mut j = m;

    loop {
        if i == 0 {
            let cur = fm[u].get(0, j);
            let mut switched_pred = None;
            for &p in &nodes[u].preds {
                if p < u && (fm[p].get(nodes[p].seq.len(), j) - cur).abs() < 1e-9 {
                    switched_pred = Some(p);
                    break;
                }
            }
            if let Some(p) = switched_pred {
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
            continue;
        }

        let cur = fm[u].get(i, j);
        if cur <= NEG_INF * 0.5 {
            break;
        }
        let post_p = if j > 0 { bm[u].get(i, j) } else { 0.0 };
        let d = if j > 0 {
            fm[u].get(i - 1, j - 1) + (post_p - 0.25)
        } else {
            NEG_INF
        };
        let up = fm[u].get(i - 1, j);
        if j > 0 && (cur - d).abs() < 1e-9 {
            let sb = joint_pair_bits(
                nodes[u].seq[i - 1],
                seq_y[j - 1],
                pi_nodes[u][i - 1],
                pi_y[j - 1],
                id_prior,
            );
            raw_cols.push((Some((u, i - 1)), Some(j - 1), sb, post_p));
            i -= 1;
            j -= 1;
        } else if (cur - up).abs() < 1e-9 {
            raw_cols.push((Some((u, i - 1)), None, ge_bits, 0.0));
            i -= 1;
        } else if j > 0 {
            raw_cols.push((None, Some(j - 1), ge_bits, 0.0));
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
        score_bits,
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

    let y_sdust = sdust_mask(seq_y);
    let mut seed_hits: HashSet<usize> = HashSet::new();
    for (mask, ry) in [(MASK_W10, false), (MASK_RY20, true)] {
        if let Ok(q_idx) =
            FlatSeedIndex::build_with_sdust(seq_y, mask, 8, false, ry, Some(&y_sdust))
            && !q_idx.keys.is_empty()
        {
            for (nid, nd) in nodes.iter().enumerate() {
                if nd.seq.len() < mask.len() {
                    continue;
                }
                if let Ok(n_idx) =
                    FlatSeedIndex::build_with_sdust(&nd.seq, mask, 8, false, ry, None)
                    && q_idx.shared_unique_keys(&n_idx) > 0
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
        let query =
            b"ACTGACTGACTGACTGACTGACTGACTGGATTACAGATTACAGATTACAGATTACACATGCATGCATGCATGCATGCATGCATG";
        let aln = align_seq_to_dag_seeded(&nodes, query, &HmmParams::default(), 10.0);
        assert_eq!(aln.path_nodes, vec![6, 8, 9]);
        assert_eq!(aln.kept_pairs.len(), query.len());
    }

    #[test]
    fn dag_handles_query_insertion_spanning_node_boundary() {
        let nodes = vec![
            GraphNode {
                name: "n0".into(),
                seq: b"ACTGACTGACTGACTGACTGACTGACTGACTG".to_vec(),
                preds: vec![],
            },
            GraphNode {
                name: "n1".into(),
                seq: b"GATTACAGATTACAGATTACAGATTACAGATT".to_vec(),
                preds: vec![0],
            },
        ];
        let query = b"ACTGACTGACTGACTGACTGACTGACTGACTGCCCCCC\
                      GATTACAGATTACAGATTACAGATTACAGATT";
        let aln = align_seq_to_dag_seeded(&nodes, query, &HmmParams::default(), 10.0);
        assert_eq!(aln.path_nodes, vec![0, 1]);
        assert!(aln.kept_pairs.iter().any(|&(n, _, _)| n == 0));
        assert!(aln.kept_pairs.iter().any(|&(n, _, _)| n == 1));
    }
}
