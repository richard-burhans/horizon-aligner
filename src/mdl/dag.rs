//! Topological DAG / partial-order Pair-HMM MDL + BOCPD alignment.

use crate::bocpd::{
    BocpdScratch, MIN_POST_PROB, NEG_INF, base_probs, bocpd_gc_profile_into, joint_pair_bits,
    logsumexp2, logsumexp3,
};
use crate::mdl::dp::{HmmParams, NormalizedHmmParams};

/// Adaptive banded per-node DP grid (`abPOA`, Gao et al. 2021).
#[derive(Clone)]
struct DagBandGrid {
    bounds: Vec<(usize, usize, usize)>, // (j_lo, j_hi, offset) per row 0..=nu
    data: Vec<f64>,
}

impl DagBandGrid {
    fn new(bounds: Vec<(usize, usize, usize)>, total_cells: usize, init: f64) -> Option<Self> {
        let mut data = Vec::new();
        if data.try_reserve_exact(total_cells).is_err() {
            return None;
        }
        data.resize(total_cells, init);
        Some(Self { bounds, data })
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

type NodeBounds = (Vec<(usize, usize, usize)>, usize);

fn build_node_bounds(
    nu: usize,
    s_min: usize,
    s_max: usize,
    total_max: usize,
    m: usize,
) -> Option<NodeBounds> {
    let mut bounds = Vec::with_capacity(nu + 1);
    let mut off = 0usize;
    let use_band = m > 128 && total_max > 0;
    let band_w = (m / 8).clamp(64, 512);
    for i in 0..=nu {
        let (j_lo, j_hi) = if use_band {
            let lo_est = super::dp::band_center(s_min + i, m, total_max);
            let hi_est = super::dp::band_center(s_max + i, m, total_max);
            (lo_est.saturating_sub(band_w), (hi_est + band_w).min(m))
        } else {
            (0, m)
        };
        let len = j_hi.saturating_sub(j_lo).checked_add(1)?;
        bounds.push((j_lo, j_hi, off));
        off = off.checked_add(len)?;
    }
    Some((bounds, off))
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
    let mut unique_preds: Vec<Vec<usize>> = vec![Vec::new(); num_nodes];
    let mut succs: Vec<Vec<usize>> = vec![Vec::new(); num_nodes];
    for (u, node) in nodes.iter().enumerate() {
        for &p in &node.preds {
            if p >= u {
                return GraphAlignment {
                    path_nodes: Vec::new(),
                    kept_pairs: Vec::new(),
                    bubbles: usize::from(m >= 6),
                    score_bits: 0.0,
                };
            }
            if !unique_preds[u].contains(&p) {
                unique_preds[u].push(p);
                succs[p].push(u);
            }
        }
    }

    let mut bocpd_scratch = BocpdScratch::default();
    let mut gc_buf = Vec::new();
    bocpd_gc_profile_into(seq_y, 1.0 / 120.0, 80, &mut bocpd_scratch, &mut gc_buf);
    let pi_y: Vec<[f64; 4]> = gc_buf.iter().copied().map(base_probs).collect();
    let mut pi_nodes: Vec<Vec<[f64; 4]>> = Vec::with_capacity(num_nodes);
    for nd in nodes {
        bocpd_gc_profile_into(&nd.seq, 1.0 / 120.0, 80, &mut bocpd_scratch, &mut gc_buf);
        pi_nodes.push(gc_buf.iter().copied().map(base_probs).collect());
    }

    let NormalizedHmmParams {
        id_prior,
        go_bits,
        ge_bits,
        mm_bits,
        im_bits,
        ..
    } = params.normalized();

    let ln2 = 2.0_f64.ln();
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
        let Some((bounds, cells)) =
            build_node_bounds(nd.seq.len(), start_min[u], start_max[u], total_max, m)
        else {
            return GraphAlignment {
                path_nodes: Vec::new(),
                kept_pairs: Vec::new(),
                bubbles: 1,
                score_bits: 0.0,
            };
        };
        let (Some(g_fm), Some(g_fix), Some(g_fiy)) = (
            DagBandGrid::new(bounds.clone(), cells, NEG_INF),
            DagBandGrid::new(bounds.clone(), cells, NEG_INF),
            DagBandGrid::new(bounds, cells, NEG_INF),
        ) else {
            return GraphAlignment {
                path_nodes: Vec::new(),
                kept_pairs: Vec::new(),
                bubbles: 1,
                score_bits: 0.0,
            };
        };
        fm.push(g_fm);
        fix.push(g_fix);
        fiy.push(g_fiy);
    }

    for u in 0..num_nodes {
        let nu = nodes[u].seq.len();
        let (j0_lo, j0_hi) = fm[u].row_range(0);
        if unique_preds[u].is_empty() {
            fm[u].set(0, 0, 0.0);
        } else {
            for j in j0_lo..=j0_hi {
                let mut sum_fm = NEG_INF;
                let mut sum_fix = NEG_INF;
                let mut sum_fiy = NEG_INF;
                for &p in &unique_preds[u] {
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

        if unique_preds[u].is_empty() {
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

    drop(fix);
    drop(fiy);

    let alloc_like = |template: &[DagBandGrid]| -> Option<Vec<DagBandGrid>> {
        template
            .iter()
            .map(|g| DagBandGrid::new(g.bounds.clone(), g.data.len(), NEG_INF))
            .collect()
    };
    let (Some(mut bm), Some(mut bix), Some(mut biy)) =
        (alloc_like(&fm), alloc_like(&fm), alloc_like(&fm))
    else {
        return GraphAlignment {
            path_nodes: Vec::new(),
            kept_pairs: Vec::new(),
            bubbles: 1,
            score_bits,
        };
    };

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
    drop(bix);
    drop(biy);

    // Pure Forward-Backward MEA decoding over the DAG:
    // reuse `fm[u]` for MEA cumulative score and `bm[u]` for posterior match probability.
    for u in 0..num_nodes {
        let nu = nodes[u].seq.len();
        let (j0_lo, j0_hi) = fm[u].row_range(0);
        if unique_preds[u].is_empty() {
            for j in j0_lo..=j0_hi {
                fm[u].set(0, j, 0.0);
            }
        } else {
            for j in j0_lo..=j0_hi {
                let mut best_p = NEG_INF;
                for &p in &unique_preds[u] {
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
            for &p in &unique_preds[u] {
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

    let mut abstain = Vec::new();
    let (raw_kept, bubbles) = super::dp::filter_mea_pairs_and_count_bubbles(
        &raw_cols,
        &mut abstain,
        MIN_POST_PROB,
        min_block_bits,
    );
    let kept_pairs = raw_kept
        .into_iter()
        .map(|((nid, np), yp)| (nid, np, yp))
        .collect();

    GraphAlignment {
        path_nodes,
        kept_pairs,
        bubbles,
        score_bits,
    }
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
    fn dag_rejects_empty_and_nontopological_preds() {
        let params = HmmParams::default();
        let empty = align_seq_to_dag_mdl(&[], b"ACGT", &params, 10.0);
        assert_eq!(empty.kept_pairs, [] as [(usize, usize, usize); 0]);

        let nodes = vec![GraphNode {
            name: "n0".into(),
            seq: b"ACGTACGTACGT".to_vec(),
            preds: vec![0],
        }];
        let bad_topo = align_seq_to_dag_mdl(&nodes, b"ACGTACGTACGT", &params, 10.0);
        assert_eq!(bad_topo.kept_pairs, [] as [(usize, usize, usize); 0]);
        assert_eq!(bad_topo.bubbles, 1);
    }
}
