//! Banded MDL + BOCPD + KP-14 aligner with Horizon-Bubble abstention.

use std::collections::HashMap;
use thiserror::Error;

/// Canonical KP-14 spaced-seed mask (span 15, weight 12).
pub const MASK_KP14: &str = "111011011011011";

const NEG_INF: f64 = -1e18;
const GAP_OPEN_BITS: f64 = -4.5;
const GAP_EXTEND_BITS: f64 = -0.8;
const IDENTITY_PRIOR: f64 = 0.65;
const MIN_POST_PROB: f64 = 0.35;

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
    if offsets.is_empty() || offsets.len() > 31 {
        return Err(AlignError::InvalidMask(
            "active mask bits must be in 1..=31".to_owned(),
        ));
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
}
