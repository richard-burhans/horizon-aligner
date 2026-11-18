//! Banded MDL + BOCPD + KP-14 aligner with Horizon-Bubble abstention.

pub mod egraph;
pub mod hal;

use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use thiserror::Error;

/// Canonical KP-14 spaced-seed mask (span 15, weight 12).
pub const MASK_KP14: &str = "111011011011011";
/// Medium-weight spaced-seed mask (span 14, weight 10) for sparse cross-clade regions.
pub const MASK_W10: &str = "11011011011011";
/// Secondary local rescue mask (span 15, weight 9).
pub const MASK_W9: &str = "110101101011011";

const NEG_INF: f64 = -1e18;
const GAP_OPEN_BITS: f64 = -4.5;
const GAP_EXTEND_BITS: f64 = -0.8;
const IDENTITY_PRIOR: f64 = 0.65;
const MIN_POST_PROB: f64 = 0.15;

/// Errors returned by sequence validation or alignment routines.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum AlignError {
    /// Raised when an empty mask or invalid mask character is supplied.
    #[error("invalid spaced-seed mask: {0}")]
    InvalidMask(String),
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

#[inline]
fn base_to_idx(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

#[inline]
fn base_probs(gc: f64) -> [f64; 4] {
    let p_gc = gc * 0.5;
    let p_at = (1.0 - gc) * 0.5;
    // Order matches A=0, C=1, G=2, T=3.
    [p_at, p_gc, p_gc, p_at]
}

#[inline]
fn joint_pair_bits(a: u8, b: u8, pi_x: [f64; 4], pi_y: [f64; 4], identity: f64) -> f64 {
    let (Some(ia), Some(ib)) = (base_to_idx(a), base_to_idx(b)) else {
        return -2.0;
    };
    let pa = pi_x[ia];
    let pb = pi_y[ib];
    let p_null = pa * pb;
    let p_avg = f64::midpoint(pa, pb);
    let match_term = if ia == ib { identity * p_avg } else { 0.0 };
    let p_joint = match_term + (1.0 - identity) * p_null;
    (p_joint / p_null).log2()
}

#[inline]
fn logsumexp2(a: f64, b: f64) -> f64 {
    let m = a.max(b);
    if m < -1e15 {
        NEG_INF
    } else {
        m + ((a - m).exp() + (b - m).exp()).ln()
    }
}

#[inline]
fn logsumexp3(a: f64, b: f64, c: f64) -> f64 {
    logsumexp2(logsumexp2(a, b), c)
}

/// Estimate local GC probability per base using Bayesian Online Changepoint Detection.
#[must_use]
pub fn bocpd_gc_profile(seq: &[u8], hazard: f64, max_run: usize) -> Vec<f64> {
    let mut run_probs = Vec::with_capacity(max_run);
    let mut alphas = Vec::with_capacity(max_run);
    let mut betas = Vec::with_capacity(max_run);
    run_probs.push(1.0);
    alphas.push(1.0);
    betas.push(1.0);

    let mut next_probs = Vec::with_capacity(max_run);
    let mut next_alphas = Vec::with_capacity(max_run);
    let mut next_betas = Vec::with_capacity(max_run);
    let mut gc_est = Vec::with_capacity(seq.len());

    for &b in seq {
        let is_gc = matches!(b.to_ascii_uppercase(), b'G' | b'C');
        let mut pred_gc: f64 = 0.0;
        for idx in 0..run_probs.len() {
            pred_gc += run_probs[idx] * (alphas[idx] / (alphas[idx] + betas[idx]));
        }
        gc_est.push(pred_gc.clamp(0.02, 0.98));

        next_probs.clear();
        next_alphas.clear();
        next_betas.clear();
        next_probs.push(0.0);
        next_alphas.push(1.0 + if is_gc { 1.0 } else { 0.0 });
        next_betas.push(1.0 + if is_gc { 0.0 } else { 1.0 });

        let mut cp_mass = 0.0;
        let keep_len = run_probs.len().min(max_run.saturating_sub(1));
        for idx in 0..run_probs.len() {
            let a = alphas[idx];
            let bt = betas[idx];
            let pred_x = if is_gc { a / (a + bt) } else { bt / (a + bt) };
            let joint = run_probs[idx] * pred_x;
            cp_mass += joint * hazard;
            if idx < keep_len {
                next_probs.push(joint * (1.0 - hazard));
                next_alphas.push(a + if is_gc { 1.0 } else { 0.0 });
                next_betas.push(bt + if is_gc { 0.0 } else { 1.0 });
            }
        }
        next_probs[0] = cp_mass;

        let norm: f64 = next_probs.iter().sum();
        if norm > 0.0 {
            let inv = 1.0 / norm;
            for p in &mut next_probs {
                *p *= inv;
            }
        } else {
            next_probs.fill(0.0);
            next_probs[0] = 1.0;
        }

        std::mem::swap(&mut run_probs, &mut next_probs);
        std::mem::swap(&mut alphas, &mut next_alphas);
        std::mem::swap(&mut betas, &mut next_betas);
    }

    gc_est
}

/// Extract canonical strand-symmetric spaced-seed positions with at most `max_occ` hits.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if `mask` is empty, exceeds 31 active bits, or contains invalid characters.
pub fn canonical_spaced_seeds(
    seq: &[u8],
    mask: &str,
    max_occ: usize,
) -> Result<HashMap<u64, Vec<usize>>, AlignError> {
    canonical_spaced_seeds_impl(seq, mask, max_occ, false, false)
}

/// Extract canonical spaced-seed positions with optional softmask skipping and 1-bit `R/Y`
/// transition-tolerant encoding.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if `mask` is empty, exceeds active bit limits, or contains invalid characters.
pub fn canonical_spaced_seeds_with_opts(
    seq: &[u8],
    mask: &str,
    max_occ: usize,
    skip_softmasked: bool,
    ry_mode: bool,
) -> Result<HashMap<u64, Vec<usize>>, AlignError> {
    canonical_spaced_seeds_impl(seq, mask, max_occ, skip_softmasked, ry_mode)
}

fn canonical_spaced_seeds_impl(
    seq: &[u8],
    mask: &str,
    max_occ: usize,
    skip_softmasked: bool,
    ry_mode: bool,
) -> Result<HashMap<u64, Vec<usize>>, AlignError> {
    if mask.is_empty() {
        return Err(AlignError::InvalidMask("mask is empty".to_owned()));
    }
    let bytes = mask.as_bytes();
    let span = bytes.len();
    let mut offsets = Vec::new();
    for (idx, &ch) in bytes.iter().enumerate() {
        match ch {
            b'1' => offsets.push(idx),
            b'0' => {}
            _ => {
                return Err(AlignError::InvalidMask(format!(
                    "unexpected char '{}'",
                    ch as char
                )));
            }
        }
    }
    let max_bits = if ry_mode { 62 } else { 31 };
    if offsets.is_empty() || offsets.len() > max_bits {
        return Err(AlignError::InvalidMask(format!(
            "active mask bits must be in 1..={max_bits}"
        )));
    }
    if seq.len() < span {
        return Ok(HashMap::new());
    }

    let rev_offsets: Vec<usize> = offsets.iter().rev().map(|&o| span - 1 - o).collect();
    let mut counts: HashMap<u64, usize> = HashMap::with_capacity(seq.len() / 2);
    let mut positions: HashMap<u64, Vec<usize>> = HashMap::with_capacity(seq.len() / 2);

    for i in 0..=(seq.len() - span) {
        let sub = &seq[i..i + span];
        let mut fwd = 0_u64;
        let mut rc = 0_u64;
        let mut valid = true;
        for (&fo, &ro) in offsets.iter().zip(&rev_offsets) {
            if skip_softmasked && (sub[fo].is_ascii_lowercase() || sub[ro].is_ascii_lowercase()) {
                valid = false;
                break;
            }
            let (Some(bf), Some(br)) = (base_to_idx(sub[fo]), base_to_idx(sub[ro])) else {
                valid = false;
                break;
            };
            if ry_mode {
                fwd = (fwd << 1) | ((bf & 1) as u64);
                rc = (rc << 1) | ((1 - (br & 1)) as u64);
            } else {
                fwd = (fwd << 2) | (bf as u64);
                rc = (rc << 2) | ((3 - br) as u64);
            }
        }
        if !valid {
            continue;
        }
        let key = fwd.min(rc);
        let entry = counts.entry(key).or_insert(0);
        *entry += 1;
        if *entry <= max_occ {
            positions.entry(key).or_default().push(i);
        }
    }

    positions.retain(|k, _| counts.get(k).copied().unwrap_or(usize::MAX) <= max_occ);
    Ok(positions)
}

/// Compute the Longest Increasing Subsequence chain over `(x, y)` anchors.
#[must_use]
pub fn lis_chain(anchors: &[(usize, usize)]) -> Vec<(usize, usize)> {
    if anchors.is_empty() {
        return Vec::new();
    }
    let mut sorted = anchors.to_vec();
    sorted.sort_unstable();

    let mut tails_y: Vec<usize> = Vec::with_capacity(sorted.len());
    let mut tail_idx: Vec<usize> = Vec::with_capacity(sorted.len());
    let mut prev = vec![usize::MAX; sorted.len()];

    for (idx, &(_, y)) in sorted.iter().enumerate() {
        let pos = tails_y.partition_point(|&ty| ty < y);
        if pos == tails_y.len() {
            tails_y.push(y);
            tail_idx.push(idx);
        } else {
            tails_y[pos] = y;
            tail_idx[pos] = idx;
        }
        if pos > 0 {
            prev[idx] = tail_idx[pos - 1];
        }
    }

    let mut chain = Vec::with_capacity(tail_idx.len());
    let mut cur = tail_idx[tail_idx.len() - 1];
    while cur != usize::MAX {
        chain.push(sorted[cur]);
        cur = prev[cur];
    }
    chain.reverse();
    chain
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

/// Run banded forward-backward MDL + BOCPD alignment on one gap.
#[must_use]
#[allow(
    clippy::too_many_lines,
    clippy::similar_names,
    clippy::many_single_char_names
)]
pub fn banded_mdl_gap(seq_x: &[u8], seq_y: &[u8], band_w: usize, min_block_bits: f64) -> GapResult {
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

    let ln2 = 2.0_f64.ln();
    let go = GAP_OPEN_BITS * ln2;
    let ge = GAP_EXTEND_BITS * ln2;

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
                let v_ix =
                    (m_dp.get(i - 1, j) + GAP_OPEN_BITS).max(ix_dp.get(i - 1, j) + GAP_EXTEND_BITS);
                ix_dp.set(i, j, v_ix);
                fix.set(
                    i,
                    j,
                    logsumexp2(fm.get(i - 1, j) + go, fix.get(i - 1, j) + ge),
                );
            }
            if j > 0 {
                let v_iy =
                    (m_dp.get(i, j - 1) + GAP_OPEN_BITS).max(iy_dp.get(i, j - 1) + GAP_EXTEND_BITS);
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
                    IDENTITY_PRIOR,
                );
                sub_grid.set(i - 1, j - 1, sb);
                let prev_m = m_dp
                    .get(i - 1, j - 1)
                    .max(ix_dp.get(i - 1, j - 1))
                    .max(iy_dp.get(i - 1, j - 1));
                m_dp.set(i, j, prev_m + sb);
                let prev_fm = logsumexp3(
                    fm.get(i - 1, j - 1),
                    fix.get(i - 1, j - 1),
                    fiy.get(i - 1, j - 1),
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
            bm.set(i, j, logsumexp3(m_next, ix_open, iy_open));
            bix.set(i, j, logsumexp2(m_next, ix_ext));
            biy.set(i, j, logsumexp2(m_next, iy_ext));
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
            raw_cols.push((Some(i - 1), None, GAP_EXTEND_BITS, 0.0));
            i -= 1;
        } else if j > 0 {
            raw_cols.push((None, Some(j - 1), GAP_EXTEND_BITS, 0.0));
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

/// Cluster `(x, y)` hits on one strand into colinear LCB chains.
#[must_use]
pub fn cluster_strand_lcbs(
    hits: &[(usize, usize)],
    max_gap: usize,
    max_diag_drift: usize,
    min_anchors: usize,
    min_span_bp: usize,
    max_mean_step: f64,
) -> Vec<Vec<(usize, usize)>> {
    if hits.is_empty() {
        return Vec::new();
    }
    let mut sorted = hits.to_vec();
    sorted.sort_unstable();

    let keep = |ch: &[(usize, usize)]| -> bool {
        let span_x = ch[ch.len() - 1].0 - ch[0].0;
        #[allow(clippy::cast_precision_loss)]
        let mean_step = (span_x as f64) / (ch.len() as f64);
        ch.len() >= min_anchors && span_x >= min_span_bp && mean_step <= max_mean_step
    };

    let mut active: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut finished: Vec<Vec<(usize, usize)>> = Vec::new();

    for &(x, y) in &sorted {
        let mut best_idx: Option<usize> = None;
        let mut best_cost = usize::MAX;
        let mut still_active: Vec<Vec<(usize, usize)>> = Vec::with_capacity(active.len());

        for ch in active {
            let (lx, ly) = ch[ch.len() - 1];
            if x.saturating_sub(lx) > max_gap {
                if keep(&ch) {
                    finished.push(ch);
                }
                continue;
            }
            still_active.push(ch);
            if x > lx && y > ly {
                let dx = x - lx;
                let dy = y - ly;
                let drift = dx.abs_diff(dy);
                if dx <= max_gap && dy <= max_gap && drift <= max_diag_drift {
                    let cost = drift + dx;
                    if cost < best_cost {
                        best_cost = cost;
                        best_idx = Some(still_active.len() - 1);
                    }
                }
            }
        }
        active = still_active;
        if let Some(idx) = best_idx {
            active[idx].push((x, y));
        } else {
            active.push(vec![(x, y)]);
        }
    }

    for ch in active {
        if keep(&ch) {
            finished.push(ch);
        }
    }
    finished
}

/// Find local `w=9` spaced-seed anchors with 2-hit diagonal consistency.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if `mask` is invalid.
pub fn local_w9_subchain(
    sub_x: &[u8],
    sub_y: &[u8],
    off_x: usize,
    off_y: usize,
    mask: &str,
    max_diag_diff: usize,
    max_step_bp: usize,
) -> Result<Vec<(usize, usize)>, AlignError> {
    let ix = canonical_spaced_seeds(sub_x, mask, 2)?;
    let iy = canonical_spaced_seeds(sub_y, mask, 2)?;
    let mut hits = Vec::new();
    for (k, pxs) in &ix {
        if let Some(pys) = iy.get(k) {
            for &px in pxs {
                for &py in pys {
                    hits.push((px, py));
                }
            }
        }
    }
    let chain = lis_chain(&hits);
    let mut kept = Vec::new();

    for (idx, &(x, y)) in chain.iter().enumerate() {
        #[allow(clippy::cast_possible_wrap)]
        let diag = (x as isize) - (y as isize);
        #[allow(clippy::cast_possible_wrap)]
        let max_dd = max_diag_diff as isize;
        let has_prev = idx > 0 && {
            let (px, py) = chain[idx - 1];
            #[allow(clippy::cast_possible_wrap)]
            let pdiag = (px as isize) - (py as isize);
            (diag - pdiag).abs() <= max_dd && (x - px) <= max_step_bp
        };
        let has_next = idx + 1 < chain.len() && {
            let (nx, ny) = chain[idx + 1];
            #[allow(clippy::cast_possible_wrap)]
            let ndiag = (nx as isize) - (ny as isize);
            (diag - ndiag).abs() <= max_dd && (nx - x) <= max_step_bp
        };
        if has_prev || has_next {
            kept.push((off_x + x, off_y + y));
        }
    }
    Ok(kept)
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

fn merge_colinear_lcbs(
    lcb_list: &[Vec<(usize, usize)>],
    max_local_rescue_bp: usize,
) -> Vec<Vec<(usize, usize)>> {
    if lcb_list.is_empty() {
        return Vec::new();
    }
    let mut ordered = lcb_list.to_vec();
    ordered.sort_by_key(|c| c[0].0);
    let mut merged = vec![ordered[0].clone()];
    for nxt in ordered.into_iter().skip(1) {
        let last = &merged[merged.len() - 1];
        let (px, py) = last[last.len() - 1];
        let (nx, ny) = nxt[0];
        if nx > px && ny > py {
            let dx = nx - px;
            let dy = ny - py;
            if dx <= max_local_rescue_bp && dy <= max_local_rescue_bp && dx.abs_diff(dy) <= 4000 {
                let m_idx = merged.len() - 1;
                merged[m_idx].extend(nxt);
                continue;
            }
        }
        merged.push(nxt);
    }
    merged
}

/// Precomputed strand-symmetric KP-14 seed index for one contig.
#[derive(Debug, Clone)]
pub struct ContigIndex {
    seq: Vec<u8>,
    seq_rc: Vec<u8>,
    seeds_fwd: HashMap<u64, Vec<usize>>,
    seeds_rc: HashMap<u64, Vec<usize>>,
}

const MASK_RY20: &str = "1110110110111011011011011";
/// Standard `LASTZ` `12-of-19` spaced-seed mask.
pub const MASK_LASTZ12: &str = "1110100110010101111";
/// Complementary `11-of-18` spaced-seed mask for multi-shape divergent seeding.
pub const MASK_LASTZ11: &str = "111010011001010111";

/// Configurable parameters for [`align_contig_pair_with_config`].
#[derive(Debug, Clone)]
pub struct AlignConfig {
    /// Minimum `w=12` LCB anchor count required before running fallback seeding and DP (`0` = disabled).
    pub min_top_lcb_anchors: usize,
    /// Skip lowercase softmasked bases during `w=10` / `R/Y` seeding.
    pub skip_softmasked: bool,
    /// Enable transition-tolerant 1-bit `R/Y` (`MASK_RY20`) fallback seeding for divergent homologs.
    pub use_ry_seeds: bool,
    /// Enable complementary `11-of-18` + `12-of-19` HSP-gated multi-shape seeding.
    pub use_multi_shape_seeds: bool,
    /// Minimum anchor count for `w=10` and `R/Y` fallback LCB clustering.
    pub fallback_min_anchors: usize,
    /// Minimum genomic span (bp) for `w=10` and `R/Y` fallback LCB clustering.
    pub fallback_min_span: usize,
    /// Maximum seed occurrence count per contig for `w=10` / `R/Y` fallback seeding.
    pub fallback_max_occ: usize,
    /// Half-band width for banded MDL+BOCPD DP.
    pub band_w: usize,
    /// Maximum X-drop extension distance (bp) at LCB flanks.
    pub max_ext: usize,
    /// Minimum MDL log-likelihood gain (bits) required to accept a DP gap/flank.
    pub min_delta_bits: f64,
    /// Minimum posterior homology probability required to emit a base pair.
    pub min_post_prob: f64,
}

impl Default for AlignConfig {
    fn default() -> Self {
        Self {
            min_top_lcb_anchors: 0,
            skip_softmasked: false,
            use_ry_seeds: false,
            use_multi_shape_seeds: true,
            fallback_min_anchors: 4,
            fallback_min_span: 50,
            fallback_max_occ: 2,
            band_w: 180,
            max_ext: 1500,
            min_delta_bits: 10.0,
            min_post_prob: MIN_POST_PROB,
        }
    }
}

/// Precompute forward and reverse-complement KP-14 seed maps for a contig.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if `MASK_KP14` is invalid.
pub fn prepare_contig_index(seq: &[u8]) -> Result<ContigIndex, AlignError> {
    prepare_contig_index_opts(seq, false)
}

/// Precompute forward and reverse-complement KP-14 seed maps, optionally skipping lowercase
/// softmasked bases during seeding.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if `MASK_KP14` is invalid.
pub fn prepare_contig_index_opts(
    seq: &[u8],
    skip_softmasked: bool,
) -> Result<ContigIndex, AlignError> {
    let seq_rc_case: Vec<u8> = seq
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => other,
        })
        .collect();
    let seeds_fwd = canonical_spaced_seeds_with_opts(seq, MASK_KP14, 2, skip_softmasked, false)?;
    let seeds_rc =
        canonical_spaced_seeds_with_opts(&seq_rc_case, MASK_KP14, 2, skip_softmasked, false)?;
    Ok(ContigIndex {
        seq: seq.to_vec(),
        seq_rc: revcomp(seq),
        seeds_fwd,
        seeds_rc,
    })
}

/// Run the full bidirectional LCB + 2-hit rescue + parallel MDL/BOCPD pipeline on two contigs.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if internal masks fail validation.
pub fn align_contig_pair(seq_x: &[u8], seq_y: &[u8]) -> Result<ContigAlignment, AlignError> {
    let idx_x = prepare_contig_index(seq_x)?;
    let idx_y = prepare_contig_index(seq_y)?;
    align_contig_pair_indexed(&idx_x, &idx_y)
}

/// Count shared KP-14 seeds (on either strand) between two pre-indexed contigs.
#[must_use]
pub fn shared_seed_count(idx_x: &ContigIndex, idx_y: &ContigIndex) -> usize {
    let fwd = idx_x
        .seeds_fwd
        .keys()
        .filter(|k| idx_y.seeds_fwd.contains_key(*k))
        .count();
    let rc = idx_x
        .seeds_fwd
        .keys()
        .filter(|k| idx_y.seeds_rc.contains_key(*k))
        .count();
    fwd.max(rc)
}

type LcbList = Vec<Vec<(usize, usize)>>;

type AnchorPairs = Vec<(usize, usize)>;

fn intersect_seed_maps(
    ix: &HashMap<u64, Vec<usize>>,
    iy_fwd: &HashMap<u64, Vec<usize>>,
    iy_rc: &HashMap<u64, Vec<usize>>,
) -> (AnchorPairs, AnchorPairs) {
    let mut hits_fwd = Vec::new();
    let mut hits_rc = Vec::new();
    for (k, pxs) in ix {
        if let Some(pys) = iy_fwd.get(k) {
            for &px in pxs {
                for &py in pys {
                    hits_fwd.push((px, py));
                }
            }
        }
        if let Some(pys) = iy_rc.get(k) {
            for &px in pxs {
                for &py in pys {
                    hits_rc.push((px, py));
                }
            }
        }
    }
    (hits_fwd, hits_rc)
}

#[allow(clippy::too_many_arguments)]
fn find_spaced_transition_hits<S: std::hash::BuildHasher>(
    seq_x: &[u8],
    offsets: &[usize],
    rev_offsets: &[usize],
    span: usize,
    iy_fwd: &HashMap<u64, Vec<usize>, S>,
    iy_rc: &HashMap<u64, Vec<usize>, S>,
    skip_softmasked: bool,
    min_transitions: usize,
    max_transitions: usize,
) -> (AnchorPairs, AnchorPairs) {
    if seq_x.len() < span {
        return (Vec::new(), Vec::new());
    }
    let w = offsets.len();
    let mut hits_fwd = Vec::new();
    let mut hits_rc = Vec::new();
    let mut keys = Vec::with_capacity(w * w / 2 + w + 1);
    for i in 0..=(seq_x.len() - span) {
        let sub = &seq_x[i..i + span];
        let mut fwd = 0_u64;
        let mut rc = 0_u64;
        let mut valid = true;
        for (&fo, &ro) in offsets.iter().zip(rev_offsets) {
            if skip_softmasked && (sub[fo].is_ascii_lowercase() || sub[ro].is_ascii_lowercase()) {
                valid = false;
                break;
            }
            let (Some(bf), Some(br)) = (base_to_idx(sub[fo]), base_to_idx(sub[ro])) else {
                valid = false;
                break;
            };
            fwd = (fwd << 2) | (bf as u64);
            rc = (rc << 2) | ((3 - br) as u64);
        }
        if !valid {
            continue;
        }
        keys.clear();
        if min_transitions == 0 {
            keys.push(fwd.min(rc));
        }
        if min_transitions <= 1 && max_transitions >= 1 {
            for b in 0..w {
                let mf = fwd ^ (2_u64 << (2 * (w - 1 - b)));
                let mr = rc ^ (2_u64 << (2 * b));
                keys.push(mf.min(mr));
            }
        }
        if max_transitions >= 2 {
            for a in 0..w {
                let ma_f = fwd ^ (2_u64 << (2 * (w - 1 - a)));
                let ma_r = rc ^ (2_u64 << (2 * a));
                for b in (a + 1)..w {
                    let mf = ma_f ^ (2_u64 << (2 * (w - 1 - b)));
                    let mr = ma_r ^ (2_u64 << (2 * b));
                    keys.push(mf.min(mr));
                }
            }
        }
        for &key in &keys {
            if let Some(pys) = iy_fwd.get(&key) {
                for &py in pys {
                    hits_fwd.push((i, py));
                }
            }
            if let Some(pys) = iy_rc.get(&key) {
                for &py in pys {
                    hits_rc.push((i, py));
                }
            }
        }
    }
    (hits_fwd, hits_rc)
}

const LASTZ12_OFFSETS: [usize; 12] = [0, 1, 2, 4, 7, 8, 11, 13, 15, 16, 17, 18];
const LASTZ12_REV_OFFSETS: [usize; 12] = [0, 1, 2, 3, 5, 7, 10, 11, 14, 16, 17, 18];

/// Find `12-of-19` (`MASK_LASTZ12`) seed hits including 1-transition (`A<->G`, `C<->T`) neighbors.
#[must_use]
pub fn find_lastz12_1transition_hits<S: std::hash::BuildHasher>(
    seq_x: &[u8],
    iy_fwd: &HashMap<u64, Vec<usize>, S>,
    iy_rc: &HashMap<u64, Vec<usize>, S>,
    skip_softmasked: bool,
) -> (AnchorPairs, AnchorPairs) {
    find_spaced_transition_hits(
        seq_x,
        &LASTZ12_OFFSETS,
        &LASTZ12_REV_OFFSETS,
        19,
        iy_fwd,
        iy_rc,
        skip_softmasked,
        0,
        1,
    )
}

/// Find `12-of-19` (`MASK_LASTZ12`) seed hits with 2 simultaneous transitions (`A<->G`, `C<->T`).
#[must_use]
pub fn find_lastz12_2transition_hits<S: std::hash::BuildHasher>(
    seq_x: &[u8],
    iy_fwd: &HashMap<u64, Vec<usize>, S>,
    iy_rc: &HashMap<u64, Vec<usize>, S>,
    skip_softmasked: bool,
) -> (AnchorPairs, AnchorPairs) {
    find_spaced_transition_hits(
        seq_x,
        &LASTZ12_OFFSETS,
        &LASTZ12_REV_OFFSETS,
        19,
        iy_fwd,
        iy_rc,
        skip_softmasked,
        2,
        2,
    )
}

fn find_lastz11_1transition_hits<S: std::hash::BuildHasher>(
    seq_x: &[u8],
    iy_fwd: &HashMap<u64, Vec<usize>, S>,
    iy_rc: &HashMap<u64, Vec<usize>, S>,
    skip_softmasked: bool,
) -> (AnchorPairs, AnchorPairs) {
    const OFFSETS: [usize; 11] = [0, 1, 2, 4, 7, 8, 11, 13, 15, 16, 17];
    const REV_OFFSETS: [usize; 11] = [0, 1, 2, 4, 6, 9, 10, 13, 15, 16, 17];
    find_spaced_transition_hits(
        seq_x,
        &OFFSETS,
        &REV_OFFSETS,
        18,
        iy_fwd,
        iy_rc,
        skip_softmasked,
        0,
        1,
    )
}

fn find_w10_1transition_hits<S: std::hash::BuildHasher>(
    seq_x: &[u8],
    iy_fwd: &HashMap<u64, Vec<usize>, S>,
    iy_rc: &HashMap<u64, Vec<usize>, S>,
    skip_softmasked: bool,
) -> (AnchorPairs, AnchorPairs) {
    const OFFSETS: [usize; 10] = [0, 1, 3, 4, 6, 7, 9, 10, 12, 13];
    find_spaced_transition_hits(
        seq_x,
        &OFFSETS,
        &OFFSETS,
        14,
        iy_fwd,
        iy_rc,
        skip_softmasked,
        0,
        1,
    )
}

/// Promote 2-hit gapped diagonal pairs (`<= max_gap` bp apart, `<= max_ddiag` bp indel shift)
/// into candidate LCBs for banded MDL DP.
fn append_two_hit_diagonal_lcbs(
    target: &mut LcbList,
    hits: &[(usize, usize)],
    max_gap: usize,
    max_ddiag: usize,
) {
    if hits.len() < 2 {
        return;
    }
    let mut sorted = hits.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    let mut promoted: Vec<Vec<(usize, usize)>> = Vec::new();
    for i in 0..sorted.len() {
        let (x1, y1) = sorted[i];
        let mut j = i + 1;
        while j < sorted.len() && sorted[j].0 - x1 <= max_gap {
            let (x2, y2) = sorted[j];
            if x2 > x1 && y2 > y1 && y2 - y1 <= max_gap {
                let d1 = x1.cast_signed() - y1.cast_signed();
                let d2 = x2.cast_signed() - y2.cast_signed();
                if (d1 - d2).unsigned_abs() <= max_ddiag {
                    promoted.push(vec![(x1, y1), (x2, y2)]);
                    break;
                }
            }
            j += 1;
        }
    }
    for lcb in promoted {
        let (sx, ex) = (lcb[0].0, lcb[lcb.len() - 1].0);
        if !target
            .iter()
            .any(|c| sx <= c[c.len() - 1].0 && ex >= c[0].0)
        {
            target.push(lcb);
        }
    }
}

fn codon_to_aa(b0: u8, b1: u8, b2: u8) -> Option<u8> {
    // Standard genetic code indexed by (A=0, C=1, G=2, T=3): 0..20 = amino acids, 255 = stop.
    const TABLE: [u8; 64] = [
        // AAA..ATT
        8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, // CAA..CTT
        13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, // GAA..GTT
        3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, // TAA..TTT
        255, 19, 255, 19, 15, 15, 15, 15, 255, 1, 18, 1, 9, 4, 9, 4,
    ];
    let (Some(i0), Some(i1), Some(i2)) = (base_to_idx(b0), base_to_idx(b1), base_to_idx(b2)) else {
        return None;
    };
    let aa = TABLE[(i0 << 4) | (i1 << 2) | i2];
    if aa == 255 { None } else { Some(aa) }
}

const fn aa_group(ax: u8) -> u8 {
    const GROUP: [u8; 20] = [0, 1, 2, 2, 3, 0, 4, 5, 4, 5, 5, 6, 7, 6, 4, 8, 8, 5, 3, 3];
    if (ax as usize) < 20 {
        GROUP[ax as usize]
    } else {
        9
    }
}

const fn aa_pair_score(ax: u8, ay: u8) -> i32 {
    if ax == ay {
        return 15;
    }
    let gx = aa_group(ax);
    let gy = aa_group(ay);
    if gx < 9 && gx == gy { 6 } else { -6 }
}

fn index_aa_spaced_4of6_frames(
    seq: &[u8],
    max_occ: usize,
    skip_softmasked: bool,
) -> HashMap<u32, Vec<usize>> {
    // Reduced 9-group biochemical alphabet spaced 6-of-8 mask ("11101101", 24 bp span).
    const AA_OFFSETS: [usize; 6] = [0, 1, 2, 4, 5, 7];
    let mut map: HashMap<u32, Vec<usize>> = HashMap::new();
    let mut counts: HashMap<u32, usize> = HashMap::new();
    for frame in 0..3 {
        if seq.len() < frame + 24 {
            continue;
        }
        let mut pos = frame;
        while pos + 24 <= seq.len() {
            let mut key = 0_u32;
            let mut valid = true;
            for &c in &AA_OFFSETS {
                let s = &seq[pos + c * 3..pos + c * 3 + 3];
                if skip_softmasked && s.iter().any(u8::is_ascii_lowercase) {
                    valid = false;
                    break;
                }
                let Some(aa) = codon_to_aa(s[0], s[1], s[2]) else {
                    valid = false;
                    break;
                };
                key = key * 9 + u32::from(aa_group(aa));
            }
            if valid {
                let cnt = counts.entry(key).or_insert(0);
                *cnt += 1;
                if *cnt <= max_occ {
                    map.entry(key).or_default().push(pos);
                } else {
                    map.remove(&key);
                }
            }
            pos += 3;
        }
    }
    map
}

/// Find spaced `4-of-6` translated amino-acid (18 bp) hits across all reading frames.
fn find_translated_aa5_hits(
    seq_x: &[u8],
    seq_y: &[u8],
    max_occ: usize,
    skip_softmasked: bool,
) -> AnchorPairs {
    let ix = index_aa_spaced_4of6_frames(seq_x, max_occ, skip_softmasked);
    let iy = index_aa_spaced_4of6_frames(seq_y, max_occ, skip_softmasked);
    let mut hits = Vec::new();
    for (k, xs) in &ix {
        if let Some(ys) = iy.get(k) {
            for &px in xs {
                for &py in ys {
                    hits.push((px, py));
                }
            }
        }
    }
    hits
}

/// Extend an ungapped codon-level amino-acid HSP from a 5-codon (`15 bp`) seed.
fn extend_aa_ungapped_hsp(
    seq_x: &[u8],
    seq_y: &[u8],
    px: usize,
    py: usize,
    min_score: i32,
) -> Option<(usize, usize, usize, usize, i32)> {
    const AA_XDROP: i32 = 30;
    if px + 15 > seq_x.len() || py + 15 > seq_y.len() {
        return None;
    }
    let mut core = 0i32;
    for c in 0..5 {
        let ax = codon_to_aa(
            seq_x[px + c * 3],
            seq_x[px + c * 3 + 1],
            seq_x[px + c * 3 + 2],
        )?;
        let ay = codon_to_aa(
            seq_y[py + c * 3],
            seq_y[py + c * 3 + 1],
            seq_y[py + c * 3 + 2],
        )?;
        core += aa_pair_score(ax, ay);
    }

    let mut cur = 0i32;
    let mut best_left = 0i32;
    let (mut lx, mut ly) = (px, py);
    let (mut ix, mut iy) = (px, py);
    while ix >= 3 && iy >= 3 {
        ix -= 3;
        iy -= 3;
        let (Some(ax), Some(ay)) = (
            codon_to_aa(seq_x[ix], seq_x[ix + 1], seq_x[ix + 2]),
            codon_to_aa(seq_y[iy], seq_y[iy + 1], seq_y[iy + 2]),
        ) else {
            break;
        };
        cur += aa_pair_score(ax, ay);
        if cur > best_left {
            best_left = cur;
            lx = ix;
            ly = iy;
        } else if best_left - cur > AA_XDROP {
            break;
        }
    }

    cur = 0;
    let mut best_right = 0i32;
    let (mut rx, mut ry) = (px + 15, py + 15);
    ix = px + 15;
    iy = py + 15;
    while ix + 3 <= seq_x.len() && iy + 3 <= seq_y.len() {
        let (Some(ax), Some(ay)) = (
            codon_to_aa(seq_x[ix], seq_x[ix + 1], seq_x[ix + 2]),
            codon_to_aa(seq_y[iy], seq_y[iy + 1], seq_y[iy + 2]),
        ) else {
            break;
        };
        cur += aa_pair_score(ax, ay);
        ix += 3;
        iy += 3;
        if cur > best_right {
            best_right = cur;
            rx = ix;
            ry = iy;
        } else if best_right - cur > AA_XDROP {
            break;
        }
    }

    let total = core + best_left + best_right;
    (total >= min_score).then_some((lx, ly, rx, ry, total))
}

fn append_hsp_gated_aa_lcbs(
    target: &mut LcbList,
    seq_x: &[u8],
    seq_y: &[u8],
    hits: &[(usize, usize)],
    min_aa_score: i32,
) {
    let mut gated = Vec::new();
    let mut strong = Vec::new();
    for &(px, py) in hits {
        if let Some((lx, ly, rx, ry, sc)) =
            extend_aa_ungapped_hsp(seq_x, seq_y, px, py, min_aa_score)
        {
            gated.push((px, py));
            if sc >= min_aa_score + 30 && rx >= lx + 24 && ry >= ly + 24 {
                strong.push(vec![(lx, ly), (rx - 15, ry - 15)]);
            }
        }
    }
    let mut candidates = cluster_strand_lcbs(&gated, 1800, 200, 2, 24, 250.0);
    candidates.extend(strong);
    for lcb in candidates {
        let (sx, ex) = (lcb[0].0, lcb[lcb.len() - 1].0);
        if !target
            .iter()
            .any(|c| sx <= c[c.len() - 1].0 && ex >= c[0].0)
        {
            target.push(lcb);
        }
    }
}

/// Find 2-hit translated amino-acid `5-mer` (15 bp) anchors across all reading frames.
#[allow(dead_code)]
fn append_translated_aa5_lcbs(
    target: &mut LcbList,
    seq_x: &[u8],
    seq_y: &[u8],
    max_occ: usize,
    skip_softmasked: bool,
) {
    let hits = find_translated_aa5_hits(seq_x, seq_y, max_occ, skip_softmasked);
    append_two_hit_diagonal_lcbs(target, &hits, 120, 12);
}

#[inline]
fn hsp_pair_score(a: u8, b: u8) -> i32 {
    let (Some(ia), Some(ib)) = (base_to_idx(a), base_to_idx(b)) else {
        return -15;
    };
    if ia == ib {
        10
    } else if (ia ^ ib) == 2 {
        -3 // transition (A<->G, C<->T)
    } else {
        -12 // transversion
    }
}

/// Extend a 19-bp seed hit `(px, py)` without gaps using transition-aware X-drop;
/// returns `Some((start_x, start_y, end_x, end_y, score))` when `score >= min_score`.
fn extend_ungapped_hsp(
    seq_x: &[u8],
    seq_y: &[u8],
    px: usize,
    py: usize,
    min_score: i32,
) -> Option<(usize, usize, usize, usize, i32)> {
    const XDROP: i32 = 60;
    if px + 19 > seq_x.len() || py + 19 > seq_y.len() {
        return None;
    }
    let mut core = 0i32;
    for off in 0..19 {
        core += hsp_pair_score(seq_x[px + off], seq_y[py + off]);
    }

    let mut cur = 0i32;
    let mut best_left = 0i32;
    let mut lx = px;
    let mut ly = py;
    let mut ix = px;
    let mut iy = py;
    while ix > 0 && iy > 0 {
        ix -= 1;
        iy -= 1;
        cur += hsp_pair_score(seq_x[ix], seq_y[iy]);
        if cur > best_left {
            best_left = cur;
            lx = ix;
            ly = iy;
        } else if best_left - cur > XDROP {
            break;
        }
    }

    cur = 0;
    let mut best_right = 0i32;
    let mut rx = px + 19;
    let mut ry = py + 19;
    ix = px + 19;
    iy = py + 19;
    while ix < seq_x.len() && iy < seq_y.len() {
        cur += hsp_pair_score(seq_x[ix], seq_y[iy]);
        ix += 1;
        iy += 1;
        if cur > best_right {
            best_right = cur;
            rx = ix;
            ry = iy;
        } else if best_right - cur > XDROP {
            break;
        }
    }

    let total = core + best_left + best_right;
    (total >= min_score).then_some((lx, ly, rx, ry, total))
}

fn append_hsp_gated_lastz12_lcbs(
    target: &mut LcbList,
    seq_x: &[u8],
    seq_y: &[u8],
    hits: &[(usize, usize)],
    min_hsp_score: i32,
) {
    let mut filtered_hits: Vec<(usize, usize)> = Vec::new();
    let mut strong_hsps: Vec<Vec<(usize, usize)>> = Vec::new();
    for &(px, py) in hits {
        if let Some((lx, ly, rx, ry, sc)) = extend_ungapped_hsp(seq_x, seq_y, px, py, min_hsp_score)
        {
            filtered_hits.push((px, py));
            if sc >= min_hsp_score + 15 && rx >= lx + 18 && ry >= ly + 18 {
                strong_hsps.push(vec![(lx, ly), (rx - 12, ry - 12)]);
            }
        }
    }
    let mut candidates = cluster_strand_lcbs(&filtered_hits, 1800, 200, 3, 30, 250.0);
    candidates.extend(cluster_strand_lcbs(&filtered_hits, 150, 25, 2, 16, 150.0));
    candidates.extend(strong_hsps);
    for lcb in candidates {
        let (sx, ex) = (lcb[0].0, lcb[lcb.len() - 1].0);
        if !target
            .iter()
            .any(|c| sx <= c[c.len() - 1].0 && ex >= c[0].0)
        {
            target.push(lcb);
        }
    }
}

fn find_w12_raw_lcbs(idx_x: &ContigIndex, idx_y: &ContigIndex) -> (LcbList, LcbList) {
    let (hits_fwd, hits_rc) =
        intersect_seed_maps(&idx_x.seeds_fwd, &idx_y.seeds_fwd, &idx_y.seeds_rc);
    let raw_fwd = cluster_strand_lcbs(&hits_fwd, 3500, 350, 8, 100, 450.0);
    let raw_rc = cluster_strand_lcbs(&hits_rc, 3500, 350, 8, 100, 450.0);
    (raw_fwd, raw_rc)
}

fn append_nonoverlapping_lcbs(
    target: &mut LcbList,
    existing: &[Vec<(usize, usize)>],
    hits: &[(usize, usize)],
    min_anchors: usize,
    min_span: usize,
) {
    let mut candidates = cluster_strand_lcbs(hits, 1800, 200, min_anchors, min_span, 250.0);
    // Tight 2-hit diagonal pass (LASTZ-style): 2 nearby hits within 120bp and <=25bp drift.
    candidates.extend(cluster_strand_lcbs(hits, 120, 25, 2, 16, 120.0));
    for lcb in candidates {
        let (sx, ex) = (lcb[0].0, lcb[lcb.len() - 1].0);
        if !existing
            .iter()
            .chain(target.iter())
            .any(|c| sx <= c[c.len() - 1].0 && ex >= c[0].0)
        {
            target.push(lcb);
        }
    }
}

/// Fast pre-DP score: return the anchor count of the largest `w=12` LCB between two contigs.
#[must_use]
pub fn score_contig_pair_lcb(idx_x: &ContigIndex, idx_y: &ContigIndex) -> usize {
    let (raw_fwd, raw_rc) = find_w12_raw_lcbs(idx_x, idx_y);
    raw_fwd
        .iter()
        .chain(raw_rc.iter())
        .map(Vec::len)
        .max()
        .unwrap_or(0)
}

/// Align two pre-indexed contigs without recomputing KP-14 seed maps.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if internal masks fail validation.
pub fn align_contig_pair_indexed(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
) -> Result<ContigAlignment, AlignError> {
    align_contig_pair_with_config(idx_x, idx_y, &AlignConfig::default())
}

/// Align two pre-indexed contigs, pruning pairs whose largest `w=12` LCB has fewer than
/// `min_top_lcb_anchors` anchors before running fallback seeding or DP.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if internal masks fail validation.
pub fn align_contig_pair_indexed_with_min_lcb(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
    min_top_lcb_anchors: usize,
) -> Result<ContigAlignment, AlignError> {
    let cfg = AlignConfig {
        min_top_lcb_anchors,
        ..AlignConfig::default()
    };
    align_contig_pair_with_config(idx_x, idx_y, &cfg)
}

/// Align two pre-indexed contigs using an explicit [`AlignConfig`].
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if internal masks fail validation.
#[allow(clippy::too_many_lines)]
pub fn align_contig_pair_with_config(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
    cfg: &AlignConfig,
) -> Result<ContigAlignment, AlignError> {
    let span = MASK_KP14.len();
    let max_gap_bp = 2400;
    let max_local_rescue_bp = 40_000;
    let seq_x = &idx_x.seq;
    let seq_y = &idx_y.seq;
    let seq_y_rc = &idx_y.seq_rc;
    let pv_len = seq_y.len();

    let (raw_fwd, raw_rc) = find_w12_raw_lcbs(idx_x, idx_y);

    // Early exit when caller requests a minimum w=12 LCB anchor threshold.
    let top_lcb_anchors = raw_fwd
        .iter()
        .chain(raw_rc.iter())
        .map(Vec::len)
        .max()
        .unwrap_or(0);
    if cfg.min_top_lcb_anchors > 0 && top_lcb_anchors < cfg.min_top_lcb_anchors {
        return Ok(ContigAlignment {
            aligned_pairs: Vec::new(),
            lcbs_fwd: 0,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 0,
        });
    }

    // Sparse-exon fallback: w=10 exact seeds in intervals not covered by w=12 LCBs.
    let ix10 = canonical_spaced_seeds_with_opts(
        seq_x,
        MASK_W10,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
    )?;
    let iy10_fwd = canonical_spaced_seeds_with_opts(
        seq_y,
        MASK_W10,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
    )?;
    let iy10_rc = canonical_spaced_seeds_with_opts(
        seq_y_rc,
        MASK_W10,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
    )?;
    let (hits10_fwd, hits10_rc) = intersect_seed_maps(&ix10, &iy10_fwd, &iy10_rc);
    let mut combined_fwd = raw_fwd.clone();
    let mut combined_rc = raw_rc.clone();
    append_nonoverlapping_lcbs(
        &mut combined_fwd,
        &raw_fwd,
        &hits10_fwd,
        cfg.fallback_min_anchors,
        cfg.fallback_min_span,
    );
    append_nonoverlapping_lcbs(
        &mut combined_rc,
        &raw_rc,
        &hits10_rc,
        cfg.fallback_min_anchors,
        cfg.fallback_min_span,
    );

    // LASTZ-style 12-of-19 + 1-transition hits gated by ungapped X-drop HSP score (>=180).
    let iy_lz_fwd = canonical_spaced_seeds_with_opts(
        seq_y,
        MASK_LASTZ12,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
    )?;
    let iy_lz_rc = canonical_spaced_seeds_with_opts(
        seq_y_rc,
        MASK_LASTZ12,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
    )?;
    let (hits_lz_fwd, hits_lz_rc) =
        find_lastz12_1transition_hits(seq_x, &iy_lz_fwd, &iy_lz_rc, cfg.skip_softmasked);
    append_hsp_gated_lastz12_lcbs(&mut combined_fwd, seq_x, seq_y, &hits_lz_fwd, 180);
    append_hsp_gated_lastz12_lcbs(&mut combined_rc, seq_x, seq_y_rc, &hits_lz_rc, 180);

    if cfg.use_multi_shape_seeds {
        let (double_trans_fwd, double_trans_rc) =
            find_lastz12_2transition_hits(seq_x, &iy_lz_fwd, &iy_lz_rc, cfg.skip_softmasked);
        append_hsp_gated_lastz12_lcbs(&mut combined_fwd, seq_x, seq_y, &double_trans_fwd, 210);
        append_hsp_gated_lastz12_lcbs(&mut combined_rc, seq_x, seq_y_rc, &double_trans_rc, 210);

        let iy_11_fwd = canonical_spaced_seeds_with_opts(
            seq_y,
            MASK_LASTZ11,
            cfg.fallback_max_occ,
            cfg.skip_softmasked,
            false,
        )?;
        let iy_11_rc = canonical_spaced_seeds_with_opts(
            seq_y_rc,
            MASK_LASTZ11,
            cfg.fallback_max_occ,
            cfg.skip_softmasked,
            false,
        )?;
        let (hits_11_fwd, hits_11_rc) =
            find_lastz11_1transition_hits(seq_x, &iy_11_fwd, &iy_11_rc, cfg.skip_softmasked);
        append_hsp_gated_lastz12_lcbs(&mut combined_fwd, seq_x, seq_y, &hits_11_fwd, 210);
        append_hsp_gated_lastz12_lcbs(&mut combined_rc, seq_x, seq_y_rc, &hits_11_rc, 210);

        // HSP-gated 1-bit purine/pyrimidine (R/Y) seeding for multi-transition divergent homologs.
        let ix_ry =
            canonical_spaced_seeds_with_opts(seq_x, MASK_RY20, 2, cfg.skip_softmasked, true)?;
        let iy_ry_fwd =
            canonical_spaced_seeds_with_opts(seq_y, MASK_RY20, 2, cfg.skip_softmasked, true)?;
        let iy_ry_rc =
            canonical_spaced_seeds_with_opts(seq_y_rc, MASK_RY20, 2, cfg.skip_softmasked, true)?;
        let (hits_ry_fwd, hits_ry_rc) = intersect_seed_maps(&ix_ry, &iy_ry_fwd, &iy_ry_rc);
        append_hsp_gated_lastz12_lcbs(&mut combined_fwd, seq_x, seq_y, &hits_ry_fwd, 220);
        append_hsp_gated_lastz12_lcbs(&mut combined_rc, seq_x, seq_y_rc, &hits_ry_rc, 220);

        // HSP-gated w=10 (10-of-14) 1-transition seeding.
        let (hits_w10_fwd, hits_w10_rc) =
            find_w10_1transition_hits(seq_x, &iy10_fwd, &iy10_rc, cfg.skip_softmasked);
        append_hsp_gated_lastz12_lcbs(&mut combined_fwd, seq_x, seq_y, &hits_w10_fwd, 220);
        append_hsp_gated_lastz12_lcbs(&mut combined_rc, seq_x, seq_y_rc, &hits_w10_rc, 220);

        // Amino-acid HSP-gated 6-frame translated spaced 4-of-6 seeding for divergent coding homologs.
        let aa_fwd = find_translated_aa5_hits(seq_x, seq_y, 2, cfg.skip_softmasked);
        let aa_rc = find_translated_aa5_hits(seq_x, seq_y_rc, 2, cfg.skip_softmasked);
        append_hsp_gated_aa_lcbs(&mut combined_fwd, seq_x, seq_y, &aa_fwd, 90);
        append_hsp_gated_aa_lcbs(&mut combined_rc, seq_x, seq_y_rc, &aa_rc, 90);
    }

    // Transition-tolerant R/Y (1-bit purine/pyrimidine) fallback for divergent homologs.
    if cfg.use_ry_seeds {
        let ix_ry =
            canonical_spaced_seeds_with_opts(seq_x, MASK_RY20, 2, cfg.skip_softmasked, true)?;
        let iy_ry_fwd =
            canonical_spaced_seeds_with_opts(seq_y, MASK_RY20, 2, cfg.skip_softmasked, true)?;
        let iy_ry_rc =
            canonical_spaced_seeds_with_opts(seq_y_rc, MASK_RY20, 2, cfg.skip_softmasked, true)?;
        let (hits_ry_fwd, hits_ry_rc) = intersect_seed_maps(&ix_ry, &iy_ry_fwd, &iy_ry_rc);
        let prev_fwd = combined_fwd.clone();
        let prev_rc = combined_rc.clone();
        append_nonoverlapping_lcbs(
            &mut combined_fwd,
            &prev_fwd,
            &hits_ry_fwd,
            cfg.fallback_min_anchors,
            cfg.fallback_min_span,
        );
        append_nonoverlapping_lcbs(
            &mut combined_rc,
            &prev_rc,
            &hits_ry_rc,
            cfg.fallback_min_anchors,
            cfg.fallback_min_span,
        );
    }

    let lcbs_fwd = merge_colinear_lcbs(&combined_fwd, max_local_rescue_bp);
    let lcbs_rc = merge_colinear_lcbs(&combined_rc, max_local_rescue_bp);

    let mut pair_set: HashSet<(usize, usize)> = HashSet::new();
    let mut dp_jobs: Vec<(bool, usize, usize, usize, usize)> = Vec::new();
    let mut xdrop_jobs: Vec<(bool, usize, usize, usize, usize, bool)> = Vec::new();
    let mut macro_bubbles = 0usize;

    for (is_rc, lcb_list, y_work) in [
        (false, &lcbs_fwd, seq_y.as_slice()),
        (true, &lcbs_rc, seq_y_rc.as_slice()),
    ] {
        // Inter-LCB sparse-region rescue using w=10 spaced seeds (with 2-hit w=9 fallback).
        let mut rescued_inter = Vec::new();
        if lcb_list.is_empty() && seq_x.len() <= 2500 && y_work.len() <= 2500 {
            let sub10 = local_w9_subchain(seq_x, y_work, 0, 0, MASK_W10, 200, 1800)?;
            if sub10.len() >= 3 {
                rescued_inter.push(sub10);
            } else {
                let sub9 = local_w9_subchain(seq_x, y_work, 0, 0, MASK_W9, 120, 1200)?;
                if sub9.len() >= 4 {
                    rescued_inter.push(sub9);
                }
            }
        }
        if let Some(first_lcb) = lcb_list.first()
            && seq_x.len() <= 50_000
        {
            let (fx, fy) = first_lcb[0];
            if (100..=2500).contains(&fx) && (100..=2500).contains(&fy) {
                let sub =
                    local_w9_subchain(&seq_x[..fx], &y_work[..fy], 0, 0, MASK_W10, 200, 1800)?;
                if sub.len() >= 3 {
                    rescued_inter.push(sub);
                }
            }
        }
        if let Some(last_lcb) = lcb_list.last()
            && seq_x.len() <= 50_000
        {
            let lx = last_lcb[last_lcb.len() - 1].0 + span;
            let ly = last_lcb[last_lcb.len() - 1].1 + span;
            if seq_x.len() > lx
                && y_work.len() > ly
                && (100..=2500).contains(&(seq_x.len() - lx))
                && (100..=2500).contains(&(y_work.len() - ly))
            {
                let sub =
                    local_w9_subchain(&seq_x[lx..], &y_work[ly..], lx, ly, MASK_W10, 200, 1800)?;
                if sub.len() >= 3 {
                    rescued_inter.push(sub);
                }
            }
        }
        for w in lcb_list.windows(2) {
            let x0 = w[0][w[0].len() - 1].0 + span;
            let y0 = w[0][w[0].len() - 1].1 + span;
            let x1 = w[1][0].0;
            let y1 = w[1][0].1;
            if x1 > x0
                && y1 > y0
                && (100..=120_000).contains(&(x1 - x0))
                && (100..=120_000).contains(&(y1 - y0))
            {
                let sub10 = local_w9_subchain(
                    &seq_x[x0..x1],
                    &y_work[y0..y1],
                    x0,
                    y0,
                    MASK_W10,
                    200,
                    1800,
                )?;
                if sub10.len() >= 3 {
                    rescued_inter.push(sub10);
                } else {
                    let sub9 = local_w9_subchain(
                        &seq_x[x0..x1],
                        &y_work[y0..y1],
                        x0,
                        y0,
                        MASK_W9,
                        120,
                        1200,
                    )?;
                    if sub9.len() >= 4 {
                        rescued_inter.push(sub9);
                    }
                }
            }
        }

        for lcb in lcb_list.iter().chain(&rescued_inter) {
            let mut all_anchors = Vec::with_capacity(lcb.len() * 2);
            all_anchors.push(lcb[0]);
            for w in lcb.windows(2) {
                let x0 = w[0].0 + span;
                let y0 = w[0].1 + span;
                let x1 = w[1].0;
                let y1 = w[1].1;
                if x1 > x0 && y1 > y0 {
                    let lx = x1 - x0;
                    let ly = y1 - y0;
                    if (lx > max_gap_bp || ly > max_gap_bp)
                        && lx <= max_local_rescue_bp
                        && ly <= max_local_rescue_bp
                    {
                        let sub = local_w9_subchain(
                            &seq_x[x0..x1],
                            &y_work[y0..y1],
                            x0,
                            y0,
                            MASK_W9,
                            80,
                            600,
                        )?;
                        all_anchors.extend(sub);
                    }
                }
                all_anchors.push(w[1]);
            }

            for &(ax, ay) in &all_anchors {
                for off in 0..span {
                    let py = if is_rc {
                        pv_len - 1 - (ay + off)
                    } else {
                        ay + off
                    };
                    pair_set.insert((ax + off, py));
                }
            }

            let (fx0, fy0) = all_anchors[0];
            if fx0 >= 12 && fy0 >= 12 {
                xdrop_jobs.push((
                    is_rc,
                    fx0.saturating_sub(600),
                    fx0,
                    fy0.saturating_sub(600),
                    fy0,
                    true,
                ));
            }
            let ex0 = all_anchors[all_anchors.len() - 1].0 + span;
            let ey0 = all_anchors[all_anchors.len() - 1].1 + span;
            if seq_x.len().saturating_sub(ex0) >= 12 && pv_len.saturating_sub(ey0) >= 12 {
                xdrop_jobs.push((
                    is_rc,
                    ex0,
                    (ex0 + 600).min(seq_x.len()),
                    ey0,
                    (ey0 + 600).min(pv_len),
                    false,
                ));
            }

            for w in all_anchors.windows(2) {
                let x0 = w[0].0 + span;
                let y0 = w[0].1 + span;
                let x1 = w[1].0;
                let y1 = w[1].1;
                if x1 <= x0 || y1 <= y0 {
                    continue;
                }
                let lx = x1 - x0;
                let ly = y1 - y0;
                if lx > max_gap_bp || ly > max_gap_bp {
                    macro_bubbles += 1;
                    xdrop_jobs.push((is_rc, x0, x1, y0, y1, false));
                    xdrop_jobs.push((is_rc, x0, x1, y0, y1, true));
                    continue;
                }
                dp_jobs.push((is_rc, x0, x1, y0, y1));
            }
        }
    }

    // Run banded MDL DP gaps in parallel via rayon.
    let dp_results: Vec<(bool, usize, usize, GapResult)> = dp_jobs
        .par_iter()
        .map(|&(is_rc, x0, x1, y0, y1)| {
            let y_work = if is_rc { seq_y_rc } else { seq_y };
            let res = banded_mdl_gap(&seq_x[x0..x1], &y_work[y0..y1], 180, 10.0);
            (is_rc, x0, y0, res)
        })
        .collect();

    let mut internal_bubbles = 0usize;
    for (is_rc, x0, y0, res) in dp_results {
        internal_bubbles += res.bubbles;
        for (dx, dy) in res.kept_pairs {
            let py = if is_rc {
                pv_len - 1 - (y0 + dy)
            } else {
                y0 + dy
            };
            pair_set.insert((x0 + dx, py));
        }
    }

    let xdrop_results: Vec<(bool, Vec<(usize, usize)>)> = xdrop_jobs
        .par_iter()
        .map(|&(is_rc, x0, x1, y0, y1, rev)| {
            let y_work = if is_rc { seq_y_rc } else { seq_y };
            let pairs = mdl_xdrop_extend(&seq_x[x0..x1], &y_work[y0..y1], x0, y0, rev, 1500, 40);
            (is_rc, pairs)
        })
        .collect();

    for (is_rc, pairs) in xdrop_results {
        for (px, py) in pairs {
            let y_orig = if is_rc { pv_len - 1 - py } else { py };
            pair_set.insert((px, y_orig));
        }
    }

    let mut aligned_pairs: Vec<(usize, usize)> = pair_set.into_iter().collect();
    aligned_pairs.sort_unstable();

    Ok(ContigAlignment {
        aligned_pairs,
        lcbs_fwd: lcbs_fwd.len(),
        lcbs_rc: lcbs_rc.len(),
        macro_bubbles,
        internal_bubbles,
    })
}

/// Format aligned LCBs and Horizon Bubbles as minimal GFA v1.1 records.
#[must_use]
pub fn format_gfa_summary(name_x: &str, name_y: &str, aln: &ContigAlignment) -> String {
    let mut out = String::from("H\tVN:Z:1.1\n");
    let total_bubbles = aln.macro_bubbles + aln.internal_bubbles;
    let _ = writeln!(
        out,
        "S\t{name_x}__{name_y}\t*\tLN:i:{}\tLF:i:{}\tLR:i:{}\tHB:i:{total_bubbles}",
        aln.aligned_pairs.len(),
        aln.lcbs_fwd,
        aln.lcbs_rc
    );
    out
}

/// Emit detailed GFA v1.1 block segments, Horizon-Bubble nodes, and links for a contig pair.
#[must_use]
pub fn format_gfa_blocks(name_x: &str, name_y: &str, aln: &ContigAlignment) -> String {
    let mut out = String::from("H\tVN:Z:1.1\n");
    if aln.aligned_pairs.is_empty() {
        return out;
    }

    let runs: Vec<(usize, usize, usize, char)> = partition_colinear_blocks(&aln.aligned_pairs)
        .into_iter()
        .map(|(blk, is_rev)| {
            let strand = if is_rev { '-' } else { '+' };
            (blk[0].0, blk[0].1, blk.len(), strand)
        })
        .collect();

    for (idx, &(bx, by, blen, bst)) in runs.iter().enumerate() {
        let seg_id = format!("{name_x}__{name_y}__blk_{idx}");
        let _ = writeln!(
            out,
            "S\t{seg_id}\t*\tLN:i:{blen}\tXS:i:{bx}\tYS:i:{by}\tST:Z:{bst}"
        );
        if idx > 0 {
            let prev_seg = format!("{name_x}__{name_y}__blk_{}", idx - 1);
            let gap_x = bx.saturating_sub(runs[idx - 1].0 + runs[idx - 1].2);
            if gap_x >= 6 {
                let bub_id = format!("{name_x}__{name_y}__bub_{idx}");
                let _ = writeln!(out, "S\t{bub_id}\t*\tLN:i:{gap_x}\tHB:i:1");
                let _ = writeln!(out, "L\t{prev_seg}\t+\t{bub_id}\t+\t0M");
                let _ = writeln!(out, "L\t{bub_id}\t+\t{seg_id}\t+\t0M");
            } else {
                let _ = writeln!(out, "L\t{prev_seg}\t+\t{seg_id}\t+\t0M");
            }
        }
    }
    out
}

/// Run e-graph duplication/inversion collapse on a sequence of signed block labels.
#[must_use]
pub fn resolve_lcb_duplications(blocks: &[(&str, i8)]) -> Option<egraph::Expr> {
    let mut eg = egraph::EGraph::new();
    let root = egraph::build_seq_expr(&mut eg, blocks)?;
    egraph::saturate_synteny(&mut eg, 4, 800);
    egraph::extract_min_cost(&mut eg, root).map(|(_, expr)| expr)
}

/// Bucket colinear alignment blocks into `bin_bp` reference windows and return
/// deduplicated block labels collapsed by the e-graph duplication/inversion rules.
#[must_use]
pub fn resolve_alignment_duplications(aln: &ContigAlignment, bin_bp: usize) -> Vec<String> {
    let bins: Vec<(String, i8)> = partition_colinear_blocks(&aln.aligned_pairs)
        .into_iter()
        .map(|(blk, is_rev)| {
            let b_idx = blk[0].0 / bin_bp.max(1);
            (format!("B{b_idx}"), if is_rev { -1 } else { 1 })
        })
        .collect();
    let mut dups: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    for win in bins.chunks(12) {
        if win.len() < 2 {
            continue;
        }
        let refs: Vec<(&str, i8)> = win.iter().map(|(s, st)| (s.as_str(), *st)).collect();
        if let Some(expr) = resolve_lcb_duplications(&refs) {
            for d in egraph::collect_dups(&expr) {
                dups.insert(d);
            }
        }
    }
    dups.into_iter().collect()
}

/// Format e-graph-resolved paralog annotations as GFA v1.1 comment/tag lines.
#[must_use]
pub fn format_gfa_paralog_tags(name_x: &str, name_y: &str, dups: &[String]) -> String {
    if dups.is_empty() {
        return String::new();
    }
    format!("# EG:Z:{name_x}__{name_y}\t{}\n", dups.join(","))
}

struct PairwiseMafBlock {
    ref_start: usize,
    ref_end: usize,
    ref_src_size: usize,
    qry_species: String,
    qry_contig: String,
    qry_start: usize,
    qry_strand: char,
    qry_src_size: usize,
    ref_seq: Vec<u8>,
    qry_seq: Vec<u8>,
}

struct ProjectedInterval {
    ref_bases: Vec<u8>,
    qry_bases: Vec<u8>,
    ins_after: Vec<Vec<u8>>,
    qry_start: usize,
    qry_size: usize,
    qry_src_size: usize,
}

fn project_block_interval(
    blk: &PairwiseMafBlock,
    lo: usize,
    hi: usize,
) -> Option<ProjectedInterval> {
    if blk.ref_start > lo || blk.ref_end < hi {
        return None;
    }
    let span = hi - lo;
    let mut ref_pos = blk.ref_start;
    let mut qry_before = 0usize;
    let mut qry_in = 0usize;
    let mut ref_bases = Vec::with_capacity(span);
    let mut qry_bases = Vec::with_capacity(span);
    let mut ins_after = vec![Vec::new(); span];

    for (&rb, &qb) in blk.ref_seq.iter().zip(blk.qry_seq.iter()) {
        if rb != b'-' {
            if (lo..hi).contains(&ref_pos) {
                ref_bases.push(rb);
                qry_bases.push(qb);
                if qb != b'-' {
                    qry_in += 1;
                }
            } else if ref_pos < lo && qb != b'-' {
                qry_before += 1;
            }
            ref_pos += 1;
        } else if ref_pos <= lo {
            if qb != b'-' {
                qry_before += 1;
            }
        } else if ref_pos < hi {
            let prev_idx = ref_pos - lo - 1;
            if qb != b'-' {
                ins_after[prev_idx].push(qb);
                qry_in += 1;
            }
        }
    }

    if ref_bases.len() != span || qry_in == 0 {
        return None;
    }
    Some(ProjectedInterval {
        ref_bases,
        qry_bases,
        ins_after,
        qry_start: blk.qry_start + qry_before,
        qry_size: qry_in,
        qry_src_size: blk.qry_src_size,
    })
}

type ThreadedQueryRow = (String, String, usize, usize, char, usize, String);

fn thread_interval_columns(
    ref_bases: &[u8],
    sp_proj: &[(String, String, char, ProjectedInterval)],
) -> (String, Vec<ThreadedQueryRow>) {
    let span = ref_bases.len();
    let max_ins: Vec<usize> = (0..span)
        .map(|i| {
            sp_proj
                .iter()
                .map(|(_, _, _, p)| p.ins_after[i].len())
                .max()
                .unwrap_or(0)
        })
        .collect();
    let total_cols = span + max_ins.iter().sum::<usize>();
    let mut ref_out = String::with_capacity(total_cols);
    let mut q_outs: Vec<String> = vec![String::with_capacity(total_cols); sp_proj.len()];

    for (i, &rb) in ref_bases.iter().enumerate() {
        ref_out.push(rb as char);
        for (q_idx, (_, _, _, proj)) in sp_proj.iter().enumerate() {
            q_outs[q_idx].push(proj.qry_bases[i] as char);
        }
        for k in 0..max_ins[i] {
            ref_out.push('-');
            for (q_idx, (_, _, _, proj)) in sp_proj.iter().enumerate() {
                let ch = proj.ins_after[i].get(k).copied().unwrap_or(b'-') as char;
                q_outs[q_idx].push(ch);
            }
        }
    }

    let rows = sp_proj
        .iter()
        .zip(q_outs)
        .map(|((sp, qc, st, proj), qseq)| {
            (
                sp.clone(),
                qc.clone(),
                proj.qry_start,
                proj.qry_size,
                *st,
                proj.qry_src_size,
                qseq,
            )
        })
        .collect();
    (ref_out, rows)
}

fn parse_pairwise_maf_blocks(
    query_sp: &str,
    maf_text: &str,
    step_bp: usize,
    by_contig: &mut std::collections::BTreeMap<String, Vec<PairwiseMafBlock>>,
) -> std::collections::BTreeSet<String> {
    type RawRow = (String, usize, usize, char, usize, String);
    let mut cur_rows: Vec<RawRow> = Vec::new();
    let mut tokens: Vec<(String, i8)> = Vec::new();

    let flush = |rows: &mut Vec<RawRow>,
                 by_contig: &mut std::collections::BTreeMap<String, Vec<PairwiseMafBlock>>,
                 tokens: &mut Vec<(String, i8)>| {
        if rows.len() >= 2 {
            let r = &rows[0];
            let q = &rows[1];
            let b_idx = r.1 / step_bp;
            tokens.push((format!("B{b_idx}"), if q.3 == '-' { -1 } else { 1 }));
            by_contig
                .entry(r.0.clone())
                .or_default()
                .push(PairwiseMafBlock {
                    ref_start: r.1,
                    ref_end: r.1 + r.2,
                    ref_src_size: r.4,
                    qry_species: query_sp.to_owned(),
                    qry_contig: q.0.clone(),
                    qry_start: q.1,
                    qry_strand: q.3,
                    qry_src_size: q.4,
                    ref_seq: r.5.as_bytes().to_vec(),
                    qry_seq: q.5.as_bytes().to_vec(),
                });
        }
        rows.clear();
    };

    for line in maf_text.lines() {
        if line.is_empty() || line.starts_with("##") {
            continue;
        }
        if line.starts_with("a ") {
            flush(&mut cur_rows, by_contig, &mut tokens);
        } else if let Some(rest) = line.strip_prefix("s ") {
            let parts: Vec<&str> = rest.split_whitespace().collect();
            if parts.len() >= 6
                && let (Ok(start), Ok(size), Ok(src_size)) = (
                    parts[1].parse::<usize>(),
                    parts[2].parse::<usize>(),
                    parts[4].parse::<usize>(),
                )
            {
                let strand = parts[3].chars().next().unwrap_or('+');
                cur_rows.push((
                    parts[0].to_owned(),
                    start,
                    size,
                    strand,
                    src_size,
                    parts[5].to_owned(),
                ));
            }
        }
    }
    flush(&mut cur_rows, by_contig, &mut tokens);

    let mut sp_dups = std::collections::BTreeSet::new();
    for win in tokens.chunks(12) {
        if win.len() < 2 {
            continue;
        }
        let refs: Vec<(&str, i8)> = win.iter().map(|(s, st)| (s.as_str(), *st)).collect();
        if let Some(expr) = resolve_lcb_duplications(&refs) {
            for d in egraph::collect_dups(&expr) {
                sp_dups.insert(d);
            }
        }
    }
    sp_dups
}

#[allow(clippy::too_many_arguments)]
fn emit_merged_block(
    ref_contig: &str,
    lo: usize,
    hi: usize,
    ref_src_size: usize,
    ref_sub: &str,
    sp_rows: &[ThreadedQueryRow],
    maf_out: &mut String,
    gfa_out: &mut String,
    prev_seg: &mut Option<String>,
) {
    let _ = writeln!(maf_out, "a score={}.0", sp_rows.len() + 1);
    let _ = writeln!(
        maf_out,
        "s {ref_contig} {lo} {} + {ref_src_size} {ref_sub}",
        hi - lo
    );
    for (sp, qc, qs, qz, st, qsrc, qseq) in sp_rows {
        let _ = writeln!(maf_out, "s {sp}.{qc} {qs} {qz} {st} {qsrc} {qseq}");
    }
    maf_out.push('\n');

    let seg_id = format!("ms_{ref_contig}_{lo}_{hi}");
    let _ = writeln!(
        gfa_out,
        "S\t{seg_id}\t*\tLN:i:{}\tSP:i:{}\tSN:Z:{ref_contig}\tSO:i:{lo}\tSR:i:0",
        hi - lo,
        sp_rows.len() + 1
    );
    if let Some(prev) = prev_seg.as_ref() {
        let _ = writeln!(gfa_out, "L\t{prev}\t+\t{seg_id}\t+\t0M");
    }
    *prev_seg = Some(seg_id);
}

/// Merge multiple pairwise MAF texts anchored on `ref_species` using exact reference-interval
/// splitting and base-level threading into multi-species MAF v1 + GFA v1.1 + e-graph paralog tags.
#[must_use]
pub fn merge_pairwise_mafs(
    ref_species: &str,
    pairwise_mafs: &[(String, String)],
    bin_bp: usize,
) -> (String, String, usize, usize) {
    use std::collections::{BTreeMap, BTreeSet};

    let mut by_contig: BTreeMap<String, Vec<PairwiseMafBlock>> = BTreeMap::new();
    let mut dups_by_species: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let step_bp = bin_bp.max(1);

    for (query_sp, maf_text) in pairwise_mafs {
        let sp_dups = parse_pairwise_maf_blocks(query_sp, maf_text, step_bp, &mut by_contig);
        dups_by_species
            .entry(query_sp.clone())
            .or_default()
            .extend(sp_dups);
    }

    let mut maf_out = String::from("##maf version=1 scoring=horizon_merge\n");
    let mut gfa_out = String::from("H\tVN:Z:1.1\n");
    let mut total_blocks = 0usize;
    let mut multi_species_blocks = 0usize;
    let mut prev_seg: Option<String> = None;

    for (ref_contig, mut blocks) in by_contig {
        blocks.sort_unstable_by_key(|b| (b.ref_start, b.ref_end));
        let mut bounds: Vec<usize> = blocks
            .iter()
            .flat_map(|b| [b.ref_start, b.ref_end])
            .collect();
        bounds.sort_unstable();
        bounds.dedup();

        let mut active: Vec<usize> = Vec::new();
        let mut next_idx = 0usize;

        for w in bounds.windows(2) {
            let (lo, hi) = (w[0], w[1]);
            if lo >= hi {
                continue;
            }
            while next_idx < blocks.len() && blocks[next_idx].ref_start <= lo {
                if blocks[next_idx].ref_end >= hi {
                    active.push(next_idx);
                }
                next_idx += 1;
            }
            active.retain(|&idx| blocks[idx].ref_end >= hi);
            if active.is_empty() {
                continue;
            }

            let mut ref_src_size_opt: Option<usize> = None;
            let mut ref_bases_opt: Option<Vec<u8>> = None;
            // Net overlapping alignments per query species by highest match score in [lo, hi).
            let mut best_by_sp: BTreeMap<String, (isize, String, char, ProjectedInterval)> =
                BTreeMap::new();
            for &idx in &active {
                let blk = &blocks[idx];
                if let Some(proj) = project_block_interval(blk, lo, hi) {
                    ref_src_size_opt.get_or_insert(blk.ref_src_size);
                    if ref_bases_opt.is_none() {
                        ref_bases_opt = Some(proj.ref_bases.clone());
                    }
                    let score = interval_match_score(&proj);
                    match best_by_sp.get(&blk.qry_species) {
                        Some((best_sc, _, _, _)) if *best_sc >= score => {}
                        _ => {
                            best_by_sp.insert(
                                blk.qry_species.clone(),
                                (score, blk.qry_contig.clone(), blk.qry_strand, proj),
                            );
                        }
                    }
                }
            }
            let sp_proj: Vec<(String, String, char, ProjectedInterval)> = best_by_sp
                .into_iter()
                .map(|(sp, (_, ctg, strand, proj))| (sp, ctg, strand, proj))
                .collect();
            let (Some(ref_src_size), Some(ref_bases)) = (ref_src_size_opt, ref_bases_opt) else {
                continue;
            };
            let (ref_sub, sp_rows) = thread_interval_columns(&ref_bases, &sp_proj);
            total_blocks += 1;
            if sp_rows.len() >= 2 {
                multi_species_blocks += 1;
            }
            emit_merged_block(
                &ref_contig,
                lo,
                hi,
                ref_src_size,
                &ref_sub,
                &sp_rows,
                &mut maf_out,
                &mut gfa_out,
                &mut prev_seg,
            );
        }
    }

    for (sp, dups) in dups_by_species {
        if !dups.is_empty() {
            let joined: Vec<&str> = dups.iter().map(String::as_str).collect();
            let _ = writeln!(gfa_out, "# EG:Z:{ref_species}__{sp}\t{}", joined.join(","));
        }
    }

    (maf_out, gfa_out, total_blocks, multi_species_blocks)
}

fn interval_match_score(proj: &ProjectedInterval) -> isize {
    let mut score = 0isize;
    for (&rb, &qb) in proj.ref_bases.iter().zip(&proj.qry_bases) {
        if rb != b'-' && qb != b'-' {
            if rb.eq_ignore_ascii_case(&qb) {
                score += 2;
            } else {
                score -= 1;
            }
        }
    }
    score
}

/// Partition sorted aligned pairs into colinear runs of at least 12 bases.
fn partition_colinear_blocks(pairs: &[(usize, usize)]) -> Vec<(Vec<(usize, usize)>, bool)> {
    let mut blocks: Vec<(Vec<(usize, usize)>, bool)> = Vec::new();
    let mut cur: Vec<(usize, usize)> = Vec::new();
    let mut cur_rev: Option<bool> = None;

    for &(x, y) in pairs {
        if let Some(&(px, py)) = cur.last() {
            let dx_ok = x > px && (x - px) <= 20;
            let fwd = dx_ok && y > py && (y - py) <= 20;
            let rev = dx_ok && py > y && (py - y) <= 20;
            let matches = match cur_rev {
                None => fwd || rev,
                Some(false) => fwd,
                Some(true) => rev,
            };
            if matches {
                if cur_rev.is_none() {
                    cur_rev = Some(rev);
                }
                cur.push((x, y));
                continue;
            }
        }
        if cur.len() >= 12 {
            blocks.push((std::mem::take(&mut cur), cur_rev.unwrap_or(false)));
        } else {
            cur.clear();
        }
        cur_rev = None;
        cur.push((x, y));
    }
    if cur.len() >= 12 {
        blocks.push((cur, cur_rev.unwrap_or(false)));
    }
    blocks
}

/// Group aligned `(x, y)` pairs into colinear MAF blocks and emit standard MAF v1 text.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn format_maf_blocks(
    name_x: &str,
    len_x: usize,
    name_y: &str,
    len_y: usize,
    seq_x: &[u8],
    seq_y: &[u8],
    aln: &ContigAlignment,
) -> String {
    let mut out = String::from("##maf version=1 scoring=horizon_mdl\n");
    if aln.aligned_pairs.is_empty() {
        return out;
    }

    for (blk, is_rev) in partition_colinear_blocks(&aln.aligned_pairs) {
        if blk
            .iter()
            .any(|&(x, y)| x >= seq_x.len() || y >= seq_y.len())
        {
            continue;
        }
        let x0 = blk[0].0;
        let x_end = blk[blk.len() - 1].0 + 1;
        let mut ax = String::with_capacity(blk.len() + 8);
        let mut ay = String::with_capacity(blk.len() + 8);

        let mut prev_x = blk[0].0;
        let mut prev_y = blk[0].1;
        ax.push(seq_x[prev_x].to_ascii_uppercase() as char);
        let first_y = seq_y[prev_y].to_ascii_uppercase();
        ay.push(if is_rev {
            revcomp(&[first_y])[0] as char
        } else {
            first_y as char
        });

        for &(x, y) in blk.iter().skip(1) {
            for &bx in &seq_x[(prev_x + 1)..x] {
                ax.push(bx.to_ascii_uppercase() as char);
                ay.push('-');
            }
            if is_rev {
                for &by in seq_y[(y + 1)..prev_y].iter().rev() {
                    ax.push('-');
                    ay.push(revcomp(&[by])[0] as char);
                }
            } else {
                for &by in &seq_y[(prev_y + 1)..y] {
                    ax.push('-');
                    ay.push(by.to_ascii_uppercase() as char);
                }
            }
            ax.push(seq_x[x].to_ascii_uppercase() as char);
            let by = seq_y[y].to_ascii_uppercase();
            ay.push(if is_rev {
                revcomp(&[by])[0] as char
            } else {
                by as char
            });
            prev_x = x;
            prev_y = y;
        }

        let (y_start, y_span, y_strand) = if is_rev {
            let y_lo = blk[blk.len() - 1].1;
            let y_hi = blk[0].1 + 1;
            (len_y - y_hi, y_hi - y_lo, '-')
        } else {
            let y_lo = blk[0].1;
            let y_hi = blk[blk.len() - 1].1 + 1;
            (y_lo, y_hi - y_lo, '+')
        };

        let _ = writeln!(out, "a score={}.0", blk.len());
        let _ = writeln!(out, "s {name_x} {x0} {} + {len_x} {ax}", x_end - x0);
        let _ = writeln!(
            out,
            "s {name_y} {y_start} {y_span} {y_strand} {len_y} {ay}\n"
        );
    }
    out
}

/// Reverse-complement a DNA byte slice (uppercasing ASCII bases).
#[must_use]
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => other,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bocpd_matches_python_oracle() {
        let seq = b"ATATATATATATATATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATAT";
        let gc = bocpd_gc_profile(seq, 1.0 / 120.0, 200);
        assert_eq!(gc.len(), seq.len());
        let expected_first5 = [0.5, 0.333_333_33, 0.250_694_44, 0.201_478_74, 0.168_858_77];
        for (got, exp) in gc.iter().take(5).zip(expected_first5) {
            assert!((got - exp).abs() < 1e-6, "got {got}, expected {exp}");
        }
        let expected_mid5 = [
            0.733_457_91,
            0.776_805_99,
            0.806_661_98,
            0.829_114_1,
            0.846_848_01,
        ];
        for (got, exp) in gc[25..30].iter().zip(expected_mid5) {
            assert!((got - exp).abs() < 1e-6, "got {got}, expected {exp}");
        }
    }

    #[test]
    fn spaced_seeds_are_strand_symmetric_and_match_oracle() {
        let seq = b"ATATATATATATATATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATAT";
        let seeds = canonical_spaced_seeds(seq, MASK_KP14, 2).expect("valid mask");
        assert_eq!(seeds.len(), 28);

        let rc: Vec<u8> = seq
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                other => other,
            })
            .collect();
        let rc_seeds = canonical_spaced_seeds(&rc, MASK_KP14, 2).expect("valid mask");
        assert_eq!(seeds.len(), rc_seeds.len());
        for key in seeds.keys() {
            assert!(rc_seeds.contains_key(key));
        }
    }

    #[test]
    fn lis_chain_finds_longest_colinear_path() {
        let hits = [(10, 100), (20, 50), (30, 120), (40, 130), (50, 90)];
        let chain = lis_chain(&hits);
        assert_eq!(chain.len(), 3);
        for w in chain.windows(2) {
            assert!(w[0].0 < w[1].0 && w[0].1 < w[1].1);
        }
    }

    #[test]
    fn banded_mdl_abstains_on_random_at_rich_pair() {
        let mut state: u64 = 12345;
        let mut next_base = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            let r = (state >> 33) % 100;
            if r < 41 {
                b'A'
            } else if r < 82 {
                b'T'
            } else if r < 91 {
                b'C'
            } else {
                b'G'
            }
        };
        let sx: Vec<u8> = (0..160).map(|_| next_base()).collect();
        let sy: Vec<u8> = (0..160).map(|_| next_base()).collect();
        let res = banded_mdl_gap(&sx, &sy, 48, 12.0);
        assert!(
            res.kept_pairs.len() < 40,
            "expected heavy abstention, got {}",
            res.kept_pairs.len()
        );
        assert!(res.bubbles >= 1);
    }

    #[test]
    fn banded_mdl_recovers_conserved_pair() {
        let sx = b"ATATATATATATATATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATAT";
        let sy = b"ATATATATATATCTATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATAT";
        let res = banded_mdl_gap(sx, sy, 24, 10.0);
        assert!(res.score_bits > 60.0);
        assert_eq!(res.kept_pairs.len(), 60);
        assert_eq!(res.bubbles, 0);
    }

    #[test]
    fn cluster_lcbs_separates_colinear_blocks() {
        let hits = [
            (0, 0),
            (30, 30),
            (60, 60),
            (90, 90),
            (120, 120),
            (150, 150),
            (180, 180),
            (210, 210),
        ];
        let lcbs = cluster_strand_lcbs(&hits, 500, 50, 8, 100, 100.0);
        assert_eq!(lcbs.len(), 1);
        assert_eq!(lcbs[0].len(), 8);
    }

    #[test]
    fn local_w9_requires_two_hit_diagonal_support() {
        let sx = b"AAAAAACCCCCCGGGGGGTTTTTTACGTACGTACGTACGT";
        let sy = b"AAAAAACCCCCCGGGGGGTTTTTTACGTACGTACGTACGT";
        let anchors = local_w9_subchain(sx, sy, 10, 20, MASK_W9, 80, 600).expect("valid mask");
        assert!(anchors.len() >= 2);
        assert!(anchors[0].0 >= 10 && anchors[0].1 >= 20);
    }

    #[test]
    fn mdl_xdrop_extends_conserved_flank() {
        let sx = b"ATATATATGCGCGCGCATATATATGCGCGCGCATATATAT";
        let sy = b"ATATATATGCGCGCGCATATATATGCGCGCGCATATATAT";
        let pairs = mdl_xdrop_extend(sx, sy, 0, 0, false, 100, 16);
        assert!(pairs.len() >= 30);
    }

    #[test]
    fn gfa_summary_contains_header_and_bubble_tag() {
        let aln = ContigAlignment {
            aligned_pairs: vec![(0, 0), (1, 1)],
            lcbs_fwd: 1,
            lcbs_rc: 0,
            macro_bubbles: 2,
            internal_bubbles: 3,
        };
        let gfa = format_gfa_summary("chr1_x", "chr1_y", &aln);
        assert!(gfa.starts_with("H\tVN:Z:1.1\n"));
        assert!(gfa.contains("HB:i:5"));
    }

    #[test]
    fn indexed_alignment_matches_unindexed() {
        let sx = b"ATATATATATATATATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATATATATATATATATATATGCGCGCGCGCGCGCGCGCGC";
        let sy = b"ATATATATATATCTATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATATATATATATATATATATGCGCGCGCGCGCGCGCGCGC";
        let unindexed = align_contig_pair(sx, sy).expect("valid");
        let ix = prepare_contig_index(sx).expect("valid");
        let iy = prepare_contig_index(sy).expect("valid");
        let indexed = align_contig_pair_indexed(&ix, &iy).expect("valid");
        assert_eq!(unindexed, indexed);
    }

    #[test]
    fn shared_seed_count_separates_homologs_from_unrelated() {
        let sx = b"ATATATATATATATATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATAT";
        let sy_same = b"ATATATATATATCTATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATAT";
        let sy_diff = b"CCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCC";
        let ix = prepare_contig_index(sx).expect("valid");
        let iy_same = prepare_contig_index(sy_same).expect("valid");
        let iy_diff = prepare_contig_index(sy_diff).expect("valid");
        assert!(shared_seed_count(&ix, &iy_same) >= 10);
        assert_eq!(shared_seed_count(&ix, &iy_diff), 0);
    }

    #[test]
    fn maf_formatter_emits_valid_blocks() {
        let sx = b"ACGTACGTACGTACGT";
        let sy = b"ACGTACGTACGTACGT";
        let aln = ContigAlignment {
            aligned_pairs: (0..16).map(|i| (i, i)).collect(),
            lcbs_fwd: 1,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let maf = format_maf_blocks("chr1", 16, "chr2", 16, sx, sy, &aln);
        assert!(maf.starts_with("##maf version=1"));
        assert!(maf.contains("s chr1 0 16 + 16 ACGTACGTACGTACGT"));
        assert!(maf.contains("s chr2 0 16 + 16 ACGTACGTACGTACGT"));
    }

    #[test]
    fn gfa_blocks_emits_segments_bubbles_and_links() {
        let mut pairs: Vec<(usize, usize)> = (0..15).map(|i| (i, i)).collect();
        pairs.extend((100..115).map(|i| (i, i)));
        let aln = ContigAlignment {
            aligned_pairs: pairs,
            lcbs_fwd: 1,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 1,
        };
        let gfa = format_gfa_blocks("c1", "c2", &aln);
        assert!(gfa.contains("S\tc1__c2__blk_0"));
        assert!(gfa.contains("S\tc1__c2__bub_1\t*\tLN:i:85\tHB:i:1"));
        assert!(gfa.contains("L\tc1__c2__blk_0\t+\tc1__c2__bub_1\t+\t0M"));
    }

    #[test]
    fn maf_and_gfa_handle_duplicate_x_and_strand_flips() {
        let sx = vec![b'A'; 64];
        let sy = vec![b'A'; 64];
        // Duplicate x=5 and a strand flip mid-stream should not panic.
        let mut pairs: Vec<(usize, usize)> = (0..14).map(|i| (i, 30 - i)).collect();
        pairs.push((5, 5));
        pairs.push((14, 31));
        pairs.sort_unstable();
        let aln = ContigAlignment {
            aligned_pairs: pairs,
            lcbs_fwd: 0,
            lcbs_rc: 1,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let maf = format_maf_blocks("c1", 64, "c2", 64, &sx, &sy, &aln);
        let gfa = format_gfa_blocks("c1", "c2", &aln);
        assert!(maf.contains("##maf version=1"));
        assert!(gfa.contains("H\tVN:Z:1.1"));
    }

    #[test]
    fn min_top_lcb_prunes_small_lcb_pairs() {
        let sx = b"ACGTACGTACGTACGTGCGCGCGCATATATATACGTACGTACGTACGT";
        let sy = b"ACGTACGTACGTACGTGCGCGCGCATATATATACGTACGTACGTACGT";
        let ix = prepare_contig_index(sx).expect("valid");
        let iy = prepare_contig_index(sy).expect("valid");
        let pruned = align_contig_pair_indexed_with_min_lcb(&ix, &iy, 10_000).expect("valid");
        assert_eq!(pruned.aligned_pairs.len(), 0);
    }

    #[test]
    fn resolve_dups_and_merge_pairwise_mafs_work() {
        let maf1 = "##maf version=1\na score=12.0\ns chr1 100 12 + 5000 ACGTACGTACGT\ns q1 10 12 + 5000 ACGTACGTACGT\n\na score=12.0\ns chr1 200 12 + 5000 ACGTACGTACGT\ns q1 200 12 + 5000 ACGTACGTACGT\n";
        let maf2 = "##maf version=1\na score=12.0\ns chr1 105 12 + 5000 CGTACGTACGTA\ns q2 50 12 - 5000 CGTACGTACGTA\n";
        let (merged_maf, merged_gfa, total, multi) = merge_pairwise_mafs(
            "Pf",
            &[
                ("Pr".to_owned(), maf1.to_owned()),
                ("Pv".to_owned(), maf2.to_owned()),
            ],
            1000,
        );
        assert_eq!(total, 4);
        assert_eq!(multi, 1);
        assert!(merged_maf.contains("s chr1 105 7 + 5000 CGTACGT"));
        assert!(merged_maf.contains("s Pr.q1"));
        assert!(merged_maf.contains("s Pv.q2"));
        assert!(merged_gfa.contains("# EG:Z:Pf__Pr\tB0"));
    }

    #[test]
    fn merge_threads_query_insertion_columns() {
        let maf1 =
            "##maf version=1\na score=4.0\ns chr1 10 4 + 100 AC--GT\ns q1 0 6 + 100 ACTTGT\n";
        let maf2 = "##maf version=1\na score=4.0\ns chr1 10 4 + 100 ACGT\ns q2 0 4 + 100 ACGT\n";
        let (merged_maf, _, total, multi) = merge_pairwise_mafs(
            "Pf",
            &[
                ("Pr".to_owned(), maf1.to_owned()),
                ("Pv".to_owned(), maf2.to_owned()),
            ],
            1000,
        );
        assert_eq!(total, 1);
        assert_eq!(multi, 1);
        assert!(merged_maf.contains("s chr1 10 4 + 100 AC--GT"));
        assert!(merged_maf.contains("s Pr.q1 0 6 + 100 ACTTGT"));
        assert!(merged_maf.contains("s Pv.q2 0 4 + 100 AC--GT"));
    }

    #[test]
    fn ry_seeds_recover_transitions_and_softmask_skips_lowercase() {
        // Transitions (A<->G, C<->T) at multiple positions: R/Y seeds still match identically.
        let sx = b"ACGTTGCACTGATCGTAGCTAGCTAAGGCCAT";
        let sy = b"ATGTTGTACTGATTGTAGCTAGCTAAGGCTAT";
        let ry_x = canonical_spaced_seeds_with_opts(sx, MASK_RY20, 2, false, true).expect("valid");
        let ry_y = canonical_spaced_seeds_with_opts(sy, MASK_RY20, 2, false, true).expect("valid");
        assert!(!ry_x.is_empty());
        assert_eq!(ry_x, ry_y);

        let masked = b"acgttgcactgatcgtagctagctaaggccat";
        let skipped =
            canonical_spaced_seeds_with_opts(masked, MASK_RY20, 2, true, true).expect("valid");
        assert_eq!(skipped.len(), 0);
    }

    #[test]
    fn lastz12_1transition_recovers_transition_hit() {
        let sx = b"ACGTTGCACTGATCGTAGCTAGCTAAGGCCAT";
        // Position 1 has C -> T transition; exact 12-of-19 misses, 1-transition recovers (0, 0).
        let sy = b"ATGTTGCACTGATCGTAGCTAGCTAAGGCCAT";
        let iy_fwd = canonical_spaced_seeds(sy, MASK_LASTZ12, 2).expect("valid");
        let iy_rc = canonical_spaced_seeds(&revcomp(sy), MASK_LASTZ12, 2).expect("valid");
        let (fwd, _) = find_lastz12_1transition_hits(sx, &iy_fwd, &iy_rc, false);
        assert!(fwd.contains(&(0, 0)), "expected (0, 0) in {fwd:?}");
    }

    #[test]
    fn lastz12_2transition_recovers_double_transition_hit() {
        let sx = b"ACGTTGCACTGATCGTAGCTAGCTAAGGCCAT";
        // Positions 1 (C->T) and 11 (A->G) both have transitions; 1-transition misses, 2-transition recovers (0, 0).
        let sy = b"ATGTTGCACTGGTCGTAGCTAGCTAAGGCCAT";
        let iy_fwd = canonical_spaced_seeds(sy, MASK_LASTZ12, 2).expect("valid");
        let iy_rc = canonical_spaced_seeds(&revcomp(sy), MASK_LASTZ12, 2).expect("valid");
        let (fwd1, _) = find_lastz12_1transition_hits(sx, &iy_fwd, &iy_rc, false);
        assert!(!fwd1.contains(&(0, 0)), "1-transition should miss (0, 0)");
        let (fwd2, _) = find_lastz12_2transition_hits(sx, &iy_fwd, &iy_rc, false);
        assert!(fwd2.contains(&(0, 0)), "expected (0, 0) in {fwd2:?}");
    }

    #[test]
    fn merge_nets_best_scoring_overlapping_block() {
        let weak =
            "##maf version=1\na score=1.0\ns chr1 0 6 + 100 AAAAAA\ns q_weak 0 6 + 100 TTTTTT\n";
        let strong =
            "##maf version=1\na score=6.0\ns chr1 0 6 + 100 AAAAAA\ns q_best 0 6 + 100 AAAAAA\n";
        let combined = format!("{weak}\n{strong}");
        let (merged_maf, _, total, _) =
            merge_pairwise_mafs("hg38", &[("mm10".to_owned(), combined)], 1000);
        assert_eq!(total, 1);
        assert!(merged_maf.contains("s mm10.q_best"));
        assert!(!merged_maf.contains("s mm10.q_weak"));
    }

    #[test]
    fn translated_aa5_recovers_synonymous_divergent_coding_pair() {
        // Two coding sequences encoding the same 12-aa peptide with synonymous 3rd-codon swaps.
        let sx = b"GCTACTGGTAAAGTTGACTTCCCAATGCAACGTTGG";
        let sy = b"GCCACAGGCAAAGTCGATTTTCCTATGCAGAGGTGG";
        let mut lcbs = Vec::new();
        append_translated_aa5_lcbs(&mut lcbs, sx, sy, 2, false);
        assert!(
            !lcbs.is_empty(),
            "expected translated 5-mer AA LCB for synonymous coding pair"
        );
        let mut aa_gated = Vec::new();
        let aa_hits = find_translated_aa5_hits(sx, sy, 2, false);
        append_hsp_gated_aa_lcbs(&mut aa_gated, sx, sy, &aa_hits, 90);
        assert!(
            !aa_gated.is_empty(),
            "expected AA HSP gate to accept synonymous 12-codon stretch"
        );
    }
}
