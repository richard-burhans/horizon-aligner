//! Banded MDL+BOCPD dynamic programming and X-drop flank extension.

use crate::bocpd::{
    GAP_EXTEND_BITS, GAP_OPEN_BITS, IDENTITY_PRIOR, MIN_POST_PROB, NEG_INF, base_probs,
    bocpd_gc_profile, joint_pair_bits, logsumexp2, logsumexp3,
};

/// Divergence-adaptive Pair-HMM parameters (in probability / bit units).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HmmParams {
    /// Expected nucleotide identity in homologous match columns (`0.50..=0.98`).
    pub identity_prior: f64,
    /// Log2 probability of opening an insertion or deletion (`< -1.0`).
    pub gap_open_bits: f64,
    /// Log2 probability of extending an insertion or deletion (`< 0.0`).
    pub gap_extend_bits: f64,
}

impl Default for HmmParams {
    fn default() -> Self {
        Self {
            identity_prior: IDENTITY_PRIOR,
            gap_open_bits: GAP_OPEN_BITS,
            gap_extend_bits: GAP_EXTEND_BITS,
        }
    }
}

/// Result of aligning one inter-anchor gap.
#[derive(Debug, Clone, PartialEq)]
pub struct GapResult {
    /// Kept base-pair coordinates `(x, y)` relative to the gap start.
    pub kept_pairs: Vec<(usize, usize)>,
    /// Number of abstained Horizon Bubbles in this gap.
    pub bubbles: usize,
    /// Total Viterbi bit-savings score for the gap.
    pub score_bits: f64,
}

struct BandGrid {
    n: usize,
    m: usize,
    band_w: usize,
    stride: usize,
    data: Vec<f64>,
}

impl BandGrid {
    fn new(n: usize, m: usize, band_w: usize, init: f64) -> Self {
        let stride = 2 * band_w + 1;
        Self {
            n,
            m,
            band_w,
            stride,
            data: vec![init; (n + 1) * stride],
        }
    }

    #[inline]
    fn center(&self, i: usize) -> isize {
        #[allow(clippy::cast_possible_wrap)]
        let ci = ((i * self.m) / self.n.max(1)) as isize;
        ci
    }

    #[inline]
    fn slot(&self, i: usize, j: usize) -> Option<usize> {
        #[allow(clippy::cast_possible_wrap)]
        let diff = (j as isize) - self.center(i);
        #[allow(clippy::cast_possible_wrap)]
        let bw = self.band_w as isize;
        if diff.abs() <= bw {
            #[allow(clippy::cast_sign_loss)]
            let col = (diff + bw) as usize;
            Some(i * self.stride + col)
        } else {
            None
        }
    }

    #[inline]
    fn get(&self, i: usize, j: usize) -> f64 {
        self.slot(i, j).map_or(NEG_INF, |idx| self.data[idx])
    }

    #[inline]
    fn set(&mut self, i: usize, j: usize, val: f64) {
        if let Some(idx) = self.slot(i, j) {
            self.data[idx] = val;
        }
    }
}

/// Run banded forward-backward MDL + BOCPD alignment on one gap with default [`HmmParams`].
#[must_use]
pub fn banded_mdl_gap(seq_x: &[u8], seq_y: &[u8], band_w: usize, min_block_bits: f64) -> GapResult {
    banded_mdl_gap_with_params(seq_x, seq_y, band_w, min_block_bits, &HmmParams::default())
}

/// Run banded forward-backward MDL + BOCPD alignment on one gap with explicit [`HmmParams`].
#[must_use]
#[allow(
    clippy::too_many_lines,
    clippy::similar_names,
    clippy::many_single_char_names
)]
pub fn banded_mdl_gap_with_params(
    seq_x: &[u8],
    seq_y: &[u8],
    band_w: usize,
    min_block_bits: f64,
    params: &HmmParams,
) -> GapResult {
    let n = seq_x.len();
    let m = seq_y.len();
    if n == 0 || m == 0 {
        return GapResult {
            kept_pairs: Vec::new(),
            bubbles: usize::from(n.max(m) >= 6),
            score_bits: 0.0,
        };
    }
    if n.min(m) <= 4 {
        let k = n.min(m);
        return GapResult {
            kept_pairs: (0..k).map(|i| (i, i)).collect(),
            bubbles: 0,
            score_bits: 0.0,
        };
    }

    let gc_x = bocpd_gc_profile(seq_x, 1.0 / 120.0, 80);
    let gc_y = bocpd_gc_profile(seq_y, 1.0 / 120.0, 80);
    let pi_x: Vec<[f64; 4]> = gc_x.into_iter().map(base_probs).collect();
    let pi_y: Vec<[f64; 4]> = gc_y.into_iter().map(base_probs).collect();

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

    let mut m_dp = BandGrid::new(n, m, band_w, NEG_INF);
    let mut ix_dp = BandGrid::new(n, m, band_w, NEG_INF);
    let mut iy_dp = BandGrid::new(n, m, band_w, NEG_INF);
    let mut fm = BandGrid::new(n, m, band_w, NEG_INF);
    let mut fix = BandGrid::new(n, m, band_w, NEG_INF);
    let mut fiy = BandGrid::new(n, m, band_w, NEG_INF);
    let mut sub_grid = BandGrid::new(n, m, band_w, 0.0);

    m_dp.set(0, 0, 0.0);
    fm.set(0, 0, 0.0);

    for i in 0..=n {
        let center = (i * m) / n;
        let j_lo = center.saturating_sub(band_w);
        let j_hi = (center + band_w).min(m);
        for j in j_lo..=j_hi {
            if i == 0 && j == 0 {
                continue;
            }
            if i > 0 {
                let v_ix = (m_dp.get(i - 1, j) + go_bits).max(ix_dp.get(i - 1, j) + ge_bits);
                ix_dp.set(i, j, v_ix);
                fix.set(
                    i,
                    j,
                    logsumexp2(fm.get(i - 1, j) + go, fix.get(i - 1, j) + ge),
                );
            }
            if j > 0 {
                let v_iy = (m_dp.get(i, j - 1) + go_bits).max(iy_dp.get(i, j - 1) + ge_bits);
                iy_dp.set(i, j, v_iy);
                fiy.set(
                    i,
                    j,
                    logsumexp2(fm.get(i, j - 1) + go, fiy.get(i, j - 1) + ge),
                );
            }
            if i > 0 && j > 0 {
                let sb = joint_pair_bits(
                    seq_x[i - 1],
                    seq_y[j - 1],
                    pi_x[i - 1],
                    pi_y[j - 1],
                    id_prior,
                );
                sub_grid.set(i - 1, j - 1, sb);
                let prev_m = (m_dp.get(i - 1, j - 1) + mm_bits)
                    .max(ix_dp.get(i - 1, j - 1) + im_bits)
                    .max(iy_dp.get(i - 1, j - 1) + im_bits);
                m_dp.set(i, j, prev_m + sb);
                let prev_fm = logsumexp3(
                    fm.get(i - 1, j - 1) + mm_ln,
                    fix.get(i - 1, j - 1) + im_ln,
                    fiy.get(i - 1, j - 1) + im_ln,
                );
                fm.set(i, j, prev_fm + sb * ln2);
            }
        }
    }

    let z_fwd = logsumexp3(fm.get(n, m), fix.get(n, m), fiy.get(n, m));
    let best_end = m_dp.get(n, m).max(ix_dp.get(n, m)).max(iy_dp.get(n, m));
    if z_fwd <= NEG_INF * 0.5 {
        return GapResult {
            kept_pairs: Vec::new(),
            bubbles: 1,
            score_bits: best_end,
        };
    }

    let mut bm = BandGrid::new(n, m, band_w, NEG_INF);
    let mut bix = BandGrid::new(n, m, band_w, NEG_INF);
    let mut biy = BandGrid::new(n, m, band_w, NEG_INF);
    bm.set(n, m, 0.0);
    bix.set(n, m, 0.0);
    biy.set(n, m, 0.0);

    for i in (0..=n).rev() {
        let center = (i * m) / n;
        let j_lo = center.saturating_sub(band_w);
        let j_hi = (center + band_w).min(m);
        for j in (j_lo..=j_hi).rev() {
            if i == n && j == m {
                continue;
            }
            let m_next = if i < n && j < m {
                bm.get(i + 1, j + 1) + sub_grid.get(i, j) * ln2
            } else {
                NEG_INF
            };
            let ix_open = if i < n {
                bix.get(i + 1, j) + go
            } else {
                NEG_INF
            };
            let ix_ext = if i < n {
                bix.get(i + 1, j) + ge
            } else {
                NEG_INF
            };
            let iy_open = if j < m {
                biy.get(i, j + 1) + go
            } else {
                NEG_INF
            };
            let iy_ext = if j < m {
                biy.get(i, j + 1) + ge
            } else {
                NEG_INF
            };
            bm.set(i, j, logsumexp3(m_next + mm_ln, ix_open, iy_open));
            bix.set(i, j, logsumexp2(m_next + im_ln, ix_ext));
            biy.set(i, j, logsumexp2(m_next + im_ln, iy_ext));
        }
    }

    let mut raw_cols: Vec<(Option<usize>, Option<usize>, f64, f64)> = Vec::with_capacity(n + m);
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        let vm = m_dp.get(i, j);
        let vix = ix_dp.get(i, j);
        let viy = iy_dp.get(i, j);
        let best = vm.max(vix).max(viy);
        if best <= NEG_INF * 0.5 {
            break;
        }
        if i > 0 && j > 0 && (vm - best).abs() < 1e-7 {
            let sb = sub_grid.get(i - 1, j - 1);
            let post_p = (fm.get(i, j) + bm.get(i, j) - z_fwd).min(0.0).exp();
            raw_cols.push((Some(i - 1), Some(j - 1), sb, post_p));
            i -= 1;
            j -= 1;
        } else if i > 0 && (vix - best).abs() < 1e-7 {
            raw_cols.push((Some(i - 1), None, ge_bits, 0.0));
            i -= 1;
        } else if j > 0 {
            raw_cols.push((None, Some(j - 1), ge_bits, 0.0));
            j -= 1;
        } else {
            break;
        }
    }
    raw_cols.reverse();

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
    let mut bubbles = 0;
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
            if let (Some(cx), Some(cy)) = (raw_cols[idx].0, raw_cols[idx].1) {
                kept_pairs.push((cx, cy));
            }
            idx += 1;
        }
    }

    GapResult {
        kept_pairs,
        bubbles,
        score_bits: best_end,
    }
}

/// Result of aligning a full contig pair.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContigAlignment {
    /// Kept 0-based `+` strand `(x, y)` base pairs.
    pub aligned_pairs: Vec<(usize, usize)>,
    /// Number of forward LCBs.
    pub lcbs_fwd: usize,
    /// Number of reverse-complement LCBs.
    pub lcbs_rc: usize,
    /// Number of macro Horizon Bubbles (`> max_gap_bp`).
    pub macro_bubbles: usize,
    /// Number of internal Horizon Bubbles emitted inside DP gaps.
    pub internal_bubbles: usize,
}

/// Extend an alignment outward from an anchor using banded MDL + BOCPD X-drop.
#[must_use]
#[allow(clippy::too_many_lines, clippy::many_single_char_names)]
pub fn mdl_xdrop_extend(
    seq_x: &[u8],
    seq_y: &[u8],
    off_x: usize,
    off_y: usize,
    reverse: bool,
    max_ext: usize,
    band_w: usize,
) -> Vec<(usize, usize)> {
    let x_slice: Vec<u8> = if reverse {
        seq_x.iter().rev().take(max_ext).copied().collect()
    } else {
        seq_x.iter().take(max_ext).copied().collect()
    };
    let y_slice: Vec<u8> = if reverse {
        seq_y.iter().rev().take(max_ext).copied().collect()
    } else {
        seq_y.iter().take(max_ext).copied().collect()
    };
    let n = x_slice.len();
    let m = y_slice.len();
    if n < 12 || m < 12 {
        return Vec::new();
    }

    let gc_x = bocpd_gc_profile(&x_slice, 1.0 / 120.0, 60);
    let gc_y = bocpd_gc_profile(&y_slice, 1.0 / 120.0, 60);
    let pi_x: Vec<[f64; 4]> = gc_x.into_iter().map(base_probs).collect();
    let pi_y: Vec<[f64; 4]> = gc_y.into_iter().map(base_probs).collect();

    let x_drop_bits = 15.0;
    let min_peak_bits = 10.0;
    let stride = 2 * band_w + 1;
    let mut dp = vec![NEG_INF; (n + 1) * stride];
    let mut bt = vec![0_u8; (n + 1) * stride]; // 1=diag, 2=up, 3=left

    let slot = |i: usize, j: usize| -> Option<usize> {
        #[allow(clippy::cast_possible_wrap)]
        let diff = (j as isize) - (i as isize);
        #[allow(clippy::cast_possible_wrap)]
        let bw = band_w as isize;
        if diff.abs() <= bw {
            #[allow(clippy::cast_sign_loss)]
            Some(i * stride + (diff + bw) as usize)
        } else {
            None
        }
    };

    if let Some(s0) = slot(0, 0) {
        dp[s0] = 0.0;
    }
    let mut best_score = 0.0;
    let mut best_cell = (0usize, 0usize);

    for i in 0..=n {
        let mut row_best = NEG_INF;
        let j_lo = i.saturating_sub(band_w);
        let j_hi = (i + band_w).min(m);
        for j in j_lo..=j_hi {
            if i == 0 && j == 0 {
                row_best = 0.0;
                continue;
            }
            let mut cand = NEG_INF;
            let mut src = 0_u8;
            if i > 0
                && j > 0
                && let Some(ps) = slot(i - 1, j - 1)
                && dp[ps] > NEG_INF * 0.5
            {
                let sb = joint_pair_bits(
                    x_slice[i - 1],
                    y_slice[j - 1],
                    pi_x[i - 1],
                    pi_y[j - 1],
                    IDENTITY_PRIOR,
                );
                let v = dp[ps] + sb;
                if v > cand {
                    cand = v;
                    src = 1;
                }
            }
            if i > 0
                && let Some(ps) = slot(i - 1, j)
                && dp[ps] > NEG_INF * 0.5
            {
                let pen = if bt[ps] == 2 {
                    GAP_EXTEND_BITS
                } else {
                    GAP_OPEN_BITS
                };
                let v = dp[ps] + pen;
                if v > cand {
                    cand = v;
                    src = 2;
                }
            }
            if j > 0
                && let Some(ps) = slot(i, j - 1)
                && dp[ps] > NEG_INF * 0.5
            {
                let pen = if bt[ps] == 3 {
                    GAP_EXTEND_BITS
                } else {
                    GAP_OPEN_BITS
                };
                let v = dp[ps] + pen;
                if v > cand {
                    cand = v;
                    src = 3;
                }
            }
            if cand >= best_score - x_drop_bits {
                if let Some(s) = slot(i, j) {
                    dp[s] = cand;
                    bt[s] = src;
                }
                if cand > row_best {
                    row_best = cand;
                }
                if cand > best_score {
                    best_score = cand;
                    best_cell = (i, j);
                }
            }
        }
        if row_best < best_score - x_drop_bits {
            break;
        }
    }

    if best_score < min_peak_bits {
        return Vec::new();
    }

    let mut pairs = Vec::new();
    let (mut ci, mut cj) = best_cell;
    while ci > 0 || cj > 0 {
        let Some(s) = slot(ci, cj) else { break };
        match bt[s] {
            1 => {
                let px = if reverse {
                    off_x + seq_x.len() - ci
                } else {
                    off_x + ci - 1
                };
                let py = if reverse {
                    off_y + seq_y.len() - cj
                } else {
                    off_y + cj - 1
                };
                pairs.push((px, py));
                ci -= 1;
                cj -= 1;
            }
            2 => ci -= 1,
            3 => cj -= 1,
            _ => break,
        }
    }
    pairs
}

/// A topologically ordered node in a sequence DAG for partial-order MDL alignment.
#[derive(Debug, Clone)]
pub struct GraphNode {
    /// Segment identifier from GFA.
    pub name: String,
    /// Segment DNA sequence.
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
    // Track best predecessor node at row 0 for each j and state (0=M, 1=Ix, 2=Iy)
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

    // Forward & Viterbi pass in topological node order
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

        // Row 0 horizontal Iy transitions
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

    // Choose best sink node (or any node with no successors) at j = m
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

    // Backward pass in reverse topological order
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

    // Viterbi traceback across DAG junctions
    let mut raw_cols: Vec<(Option<(usize, usize)>, Option<usize>, f64, f64)> = Vec::new();
    let mut path_nodes = vec![best_end_node];
    let mut u = best_end_node;
    let mut i = nodes[u].seq.len();
    let mut j = m;
    // State: 0=M, 1=Ix, 2=Iy
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
            // Move left along row 0 in state Iy
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
            state = if vm >= vix { 0 } else { 1 };
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

    // Posterior + MDL block gating (same Horizon Bubble abstention rule as linear DP)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bocpd::{GAP_EXTEND_BITS, GAP_OPEN_BITS};

    #[test]
    fn pair_hmm_transitions_are_normalized() {
        let p_open = 2.0_f64.powf(GAP_OPEN_BITS);
        let p_mm = 1.0 - 2.0 * p_open;
        assert!((p_mm + 2.0 * p_open - 1.0).abs() < 1e-12);

        let p_ext = 2.0_f64.powf(GAP_EXTEND_BITS);
        let p_im = 1.0 - p_ext;
        assert!((p_im + p_ext - 1.0).abs() < 1e-12);
    }

    #[test]
    fn dag_mdl_aligns_through_correct_bubble_branch() {
        // Diamond DAG:
        //   node 0 (anchor): ACGTACGTACGTACGTACGT (20bp)
        //   node 1 (alt A):  GGGGGGGGGGGGGGGGGGGG (20bp) <- wrong branch
        //   node 2 (alt B):  ATATATATATATATATATAT (20bp) <- matches query!
        //   node 3 (sink):   TGCATGCATGCATGCATGCA (20bp)
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
}
