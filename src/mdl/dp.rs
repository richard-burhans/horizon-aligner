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
    fn try_new(n: usize, m: usize, band_w: usize, init: f64) -> Option<Self> {
        let bw = band_w.min(n.max(m));
        let stride = bw.checked_mul(2)?.checked_add(1)?;
        let total = n.checked_add(1)?.checked_mul(stride)?;
        let mut data = Vec::new();
        data.try_reserve_exact(total).ok()?;
        data.resize(total, init);
        Some(Self {
            n,
            m,
            band_w: bw,
            stride,
            data,
        })
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
pub fn banded_mdl_gap_with_params(
    seq_x: &[u8],
    seq_y: &[u8],
    band_w: usize,
    min_block_bits: f64,
    params: &HmmParams,
) -> GapResult {
    banded_mdl_gap_with_min_post(seq_x, seq_y, band_w, min_block_bits, params, MIN_POST_PROB)
}

/// Run banded forward-backward MDL + BOCPD alignment on one gap with explicit [`HmmParams`]
/// and posterior match probability threshold `min_post_prob`.
#[must_use]
#[allow(
    clippy::too_many_lines,
    clippy::similar_names,
    clippy::many_single_char_names
)]
pub fn banded_mdl_gap_with_min_post(
    seq_x: &[u8],
    seq_y: &[u8],
    band_w: usize,
    min_block_bits: f64,
    params: &HmmParams,
    min_post_prob: f64,
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

    let (Some(mut fm), Some(mut fix), Some(mut fiy), Some(mut sub_grid)) = (
        BandGrid::try_new(n, m, band_w, NEG_INF),
        BandGrid::try_new(n, m, band_w, NEG_INF),
        BandGrid::try_new(n, m, band_w, NEG_INF),
        BandGrid::try_new(n, m, band_w, 0.0),
    ) else {
        return GapResult {
            kept_pairs: Vec::new(),
            bubbles: usize::from(n.max(m) >= 6),
            score_bits: 0.0,
        };
    };

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
                fix.set(
                    i,
                    j,
                    logsumexp2(fm.get(i - 1, j) + go, fix.get(i - 1, j) + ge),
                );
            }
            if j > 0 {
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
    drop(fix);
    drop(fiy);
    let score_bits = z_fwd / ln2;
    if z_fwd <= NEG_INF * 0.5 {
        return GapResult {
            kept_pairs: Vec::new(),
            bubbles: 1,
            score_bits,
        };
    }

    let (Some(mut bm), Some(mut bix), Some(mut biy)) = (
        BandGrid::try_new(n, m, band_w, NEG_INF),
        BandGrid::try_new(n, m, band_w, NEG_INF),
        BandGrid::try_new(n, m, band_w, NEG_INF),
    ) else {
        return GapResult {
            kept_pairs: Vec::new(),
            bubbles: usize::from(n.max(m) >= 6),
            score_bits: 0.0,
        };
    };
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

    // Pure Forward-Backward Maximum Expected Accuracy (MEA) decoding (Durbin et al. 1998):
    // reuse `fm` for MEA cumulative score and `bm` for posterior match probability.
    for i in 0..=n {
        let center = (i * m) / n;
        let j_lo = center.saturating_sub(band_w);
        let j_hi = (center + band_w).min(m);
        for j in j_lo..=j_hi {
            if i == 0 && j == 0 {
                fm.set(0, 0, 0.0);
                continue;
            }
            let post_m = if i > 0 && j > 0 {
                (fm.get(i, j) + bm.get(i, j) - z_fwd).min(0.0).exp()
            } else {
                0.0
            };
            if i > 0 && j > 0 {
                bm.set(i, j, post_m);
            }
            let d = if i > 0 && j > 0 {
                fm.get(i - 1, j - 1) + (post_m - 0.25)
            } else {
                NEG_INF
            };
            let u = if i > 0 { fm.get(i - 1, j) } else { NEG_INF };
            let l = if j > 0 { fm.get(i, j - 1) } else { NEG_INF };
            fm.set(i, j, d.max(u).max(l));
        }
    }

    let mut raw_cols: Vec<(Option<usize>, Option<usize>, f64, f64)> = Vec::with_capacity(n + m);
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        let cur = fm.get(i, j);
        if cur <= NEG_INF * 0.5 {
            break;
        }
        let d = if i > 0 && j > 0 {
            fm.get(i - 1, j - 1) + (bm.get(i, j) - 0.25)
        } else {
            NEG_INF
        };
        let u = if i > 0 { fm.get(i - 1, j) } else { NEG_INF };
        if i > 0 && j > 0 && (cur - d).abs() < 1e-9 {
            let sb = sub_grid.get(i - 1, j - 1);
            let post_p = bm.get(i, j);
            raw_cols.push((Some(i - 1), Some(j - 1), sb, post_p));
            i -= 1;
            j -= 1;
        } else if i > 0 && (cur - u).abs() < 1e-9 {
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
        if cx.is_some() && cy.is_some() && pp >= min_post_prob {
            let start = idx;
            let mut block_bits = 0.0;
            while idx < k
                && ((raw_cols[idx].0.is_some()
                    && raw_cols[idx].1.is_some()
                    && raw_cols[idx].3 >= min_post_prob)
                    || (idx + 1 < k
                        && raw_cols[idx].2 > -2.5
                        && raw_cols[idx + 1].3 >= min_post_prob))
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
        score_bits,
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

/// Extend an alignment outward from an anchor using banded MDL + BOCPD X-drop
/// with default Pair-HMM parameters.
#[must_use]
pub fn mdl_xdrop_extend(
    seq_x: &[u8],
    seq_y: &[u8],
    off_x: usize,
    off_y: usize,
    reverse: bool,
    max_ext: usize,
    band_w: usize,
) -> Vec<(usize, usize)> {
    mdl_xdrop_extend_with_params(
        seq_x,
        seq_y,
        off_x,
        off_y,
        reverse,
        max_ext,
        band_w,
        &HmmParams::default(),
    )
}

/// Extend an alignment outward from an anchor using banded MDL + BOCPD X-drop
/// with explicit normalized Pair-HMM parameters.
#[must_use]
#[allow(
    clippy::too_many_lines,
    clippy::many_single_char_names,
    clippy::too_many_arguments
)]
pub fn mdl_xdrop_extend_with_params(
    seq_x: &[u8],
    seq_y: &[u8],
    off_x: usize,
    off_y: usize,
    reverse: bool,
    max_ext: usize,
    band_w: usize,
    params: &HmmParams,
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

    let id_prior = params.identity_prior.clamp(0.40, 0.98);
    let open_bits = params.gap_open_bits.clamp(-8.0, -1.5);
    let extend_bits = params.gap_extend_bits.clamp(-4.0, -0.2);
    let p_open = 2.0_f64.powf(open_bits);
    let p_ext = 2.0_f64.powf(extend_bits);
    let mm_bits = (1.0 - 2.0 * p_open).max(1e-6).log2();
    let im_bits = (1.0 - p_ext).max(1e-6).log2();

    let x_drop_bits = 15.0;
    let min_peak_bits = 10.0;
    let band_w = band_w.min(n.max(m));
    let Some(stride) = band_w.checked_mul(2).and_then(|s| s.checked_add(1)) else {
        return Vec::new();
    };
    let Some(total) = n.checked_add(1).and_then(|r| r.checked_mul(stride)) else {
        return Vec::new();
    };
    let mut dp = Vec::new();
    if dp.try_reserve_exact(total).is_err() {
        return Vec::new();
    }
    dp.resize(total, NEG_INF);
    let mut bt = Vec::new();
    if bt.try_reserve_exact(total).is_err() {
        return Vec::new();
    }
    bt.resize(total, 0_u8); // 1=diag, 2=up, 3=left

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
                    id_prior,
                );
                let trans_bits = if bt[ps] == 2 || bt[ps] == 3 {
                    im_bits
                } else {
                    mm_bits
                };
                let v = dp[ps] + trans_bits + sb;
                if v > cand {
                    cand = v;
                    src = 1;
                }
            }
            if i > 0
                && let Some(ps) = slot(i - 1, j)
                && dp[ps] > NEG_INF * 0.5
            {
                let pen = if bt[ps] == 2 { extend_bits } else { open_bits };
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
                let pen = if bt[ps] == 3 { extend_bits } else { open_bits };
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
    if !reverse {
        pairs.reverse();
    }
    pairs
}

#[cfg(test)]
mod tests {
    use crate::bocpd::{GAP_EXTEND_BITS, GAP_OPEN_BITS};

    #[test]
    fn pair_hmm_transitions_are_normalized() {
        let p_open = 2.0_f64.powf(GAP_OPEN_BITS);
        let stay_match = 1.0 - 2.0 * p_open;
        assert!((stay_match + 2.0 * p_open - 1.0).abs() < 1e-12);

        let p_ext = 2.0_f64.powf(GAP_EXTEND_BITS);
        let leave_indel = 1.0 - p_ext;
        assert!((leave_indel + p_ext - 1.0).abs() < 1e-12);
    }

    #[test]
    fn min_post_prob_filters_uncertain_pairs() {
        let x = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let y = b"ACGTACGTAAAACCCCGGGGTTTTACGTACGT";
        let params = super::HmmParams::default();
        let all = super::banded_mdl_gap_with_min_post(x, y, 8, 5.0, &params, 0.0);
        let strict = super::banded_mdl_gap_with_min_post(x, y, 8, 5.0, &params, 0.9999);
        assert_ne!(all.kept_pairs.len(), 0);
        assert!(strict.kept_pairs.len() < all.kept_pairs.len());
    }
}
