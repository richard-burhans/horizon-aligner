//! Banded MDL+BOCPD dynamic programming for inter-anchor gap alignment.

use crate::bocpd::{
    BocpdScratch, GAP_EXTEND_BITS, GAP_OPEN_BITS, IDENTITY_PRIOR, MIN_POST_PROB, NEG_INF,
    base_probs, base_to_idx_fast, bocpd_gc_profile_into,
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

#[derive(Debug, Clone, Copy)]
pub(crate) struct NormalizedHmmParams {
    pub id_prior: f64,
    pub go_bits: f64,
    pub ge_bits: f64,
    pub p_open: f64,
    pub p_ext: f64,
    pub p_mm: f64,
    pub p_im: f64,
    pub mm_bits: f64,
    pub im_bits: f64,
}

impl HmmParams {
    #[must_use]
    #[allow(clippy::similar_names)]
    pub(crate) fn normalized(&self) -> NormalizedHmmParams {
        let id_raw = if self.identity_prior.is_finite() {
            self.identity_prior
        } else {
            IDENTITY_PRIOR
        };
        let go_raw = if self.gap_open_bits.is_finite() {
            self.gap_open_bits
        } else {
            GAP_OPEN_BITS
        };
        let ge_raw = if self.gap_extend_bits.is_finite() {
            self.gap_extend_bits
        } else {
            GAP_EXTEND_BITS
        };
        let id_prior = id_raw.clamp(0.40, 0.98);
        let go_bits = go_raw.clamp(-8.0, -1.5);
        let ge_bits = ge_raw.clamp(-4.0, -0.2);
        let p_open = 2.0_f64.powf(go_bits);
        let p_ext = 2.0_f64.powf(ge_bits);
        let p_mm = (1.0 - 2.0 * p_open).max(1e-6);
        let p_im = (1.0 - p_ext).max(1e-6);
        NormalizedHmmParams {
            id_prior,
            go_bits,
            ge_bits,
            p_open,
            p_ext,
            p_mm,
            p_im,
            mm_bits: p_mm.log2(),
            im_bits: p_im.log2(),
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

/// Reusable per-worker scratch buffers for [`banded_mdl_gap_with_scratch`].
#[derive(Debug, Default)]
pub(crate) struct DpScratch {
    bocpd: BocpdScratch,
    gc_x: Vec<f64>,
    gc_y: Vec<f64>,
    pi_x: Vec<[f64; 4]>,
    pi_y: Vec<[f64; 4]>,
    enc_x: Vec<u8>,
    enc_y: Vec<u8>,
    scale_f: Vec<f64>,
    fm: Vec<f64>,
    fix: Vec<f64>,
    fiy: Vec<f64>,
    sub_grid: Vec<f64>,
    bm: Vec<f64>,
    bix: Vec<f64>,
    biy: Vec<f64>,
    raw_cols: Vec<(Option<usize>, Option<usize>, f64, f64)>,
    abstain: Vec<bool>,
    pub(super) vit_m: Vec<f32>,
    pub(super) vit_x: Vec<f32>,
    pub(super) vit_y: Vec<f32>,
    pub(super) vit_trace: Vec<u8>,
}

pub(crate) fn prepare_grid<T: Copy>(buf: &mut Vec<T>, total: usize, init: T) -> bool {
    buf.clear();
    if buf.capacity() < total && buf.try_reserve_exact(total).is_err() {
        return false;
    }
    buf.resize(total, init);
    true
}

fn abstained_gap_result(bubbles: usize, score_bits: f64) -> GapResult {
    GapResult {
        kept_pairs: Vec::new(),
        bubbles,
        score_bits,
    }
}

struct BandLayout {
    n: usize,
    m: usize,
    band_w: usize,
    stride: usize,
    total: usize,
}

impl BandLayout {
    fn new(n: usize, m: usize, band_w: usize) -> Option<Self> {
        let bw = band_w.min(n.max(m));
        let stride = bw.checked_mul(2)?.checked_add(1)?;
        let total = n.checked_add(1)?.checked_mul(stride)?;
        Some(Self {
            n,
            m,
            band_w: bw,
            stride,
            total,
        })
    }

    #[inline]
    fn center(&self, i: usize) -> usize {
        band_center(i, self.m, self.n)
    }

    #[inline]
    fn slot(&self, i: usize, j: usize) -> Option<usize> {
        #[allow(clippy::cast_possible_wrap)]
        let diff = (j as isize) - (self.center(i) as isize);
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
    fn get(&self, data: &[f64], i: usize, j: usize) -> f64 {
        self.slot(i, j).map_or(NEG_INF, |idx| data[idx])
    }
}

#[inline]
fn emit_ratio_enc(
    ia: u8,
    ib: u8,
    pi_x: [f64; 4],
    pi_y: [f64; 4],
    identity: f64,
    mismatch_ratio: f64,
) -> f64 {
    if ia >= 4 || ib >= 4 {
        return 1.0;
    }
    if ia != ib {
        return mismatch_ratio;
    }
    let pa = pi_x[ia as usize];
    let pb = pi_y[ib as usize];
    let p_null = pa * pb;
    let p_avg = f64::midpoint(pa, pb);
    identity * p_avg / p_null + mismatch_ratio
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
pub fn banded_mdl_gap_with_min_post(
    seq_x: &[u8],
    seq_y: &[u8],
    band_w: usize,
    min_block_bits: f64,
    params: &HmmParams,
    min_post_prob: f64,
) -> GapResult {
    let mut scratch = DpScratch::default();
    banded_mdl_gap_with_scratch(
        seq_x,
        seq_y,
        band_w,
        min_block_bits,
        params,
        min_post_prob,
        &mut scratch,
    )
}

#[inline]
pub(super) fn band_center(i: usize, m: usize, n: usize) -> usize {
    usize::try_from((i as u128 * m as u128) / (n.max(1) as u128)).unwrap_or(m)
}

#[must_use]
pub(crate) fn dip_floor_bits_for_identity(identity: f64) -> f64 {
    let p_id = if identity.is_finite() {
        identity
    } else {
        IDENTITY_PRIOR
    }
    .clamp(0.25, 0.999);
    crate::bocpd::joint_pair_bits(b'A', b'C', [0.25; 4], [0.25; 4], p_id) - 0.5
}

pub(super) fn filter_mea_pairs_and_count_bubbles<X: Copy>(
    raw_cols: &[(Option<X>, Option<usize>, f64, f64)],
    abstain: &mut Vec<bool>,
    min_post_prob: f64,
    min_block_bits: f64,
    dip_floor_bits: f64,
) -> (Vec<(X, usize)>, usize) {
    let min_post_prob = if min_post_prob.is_nan() {
        0.0
    } else {
        min_post_prob
    };
    let k = raw_cols.len();
    abstain.clear();
    abstain.resize(k, false);
    if min_block_bits.is_nan() || min_block_bits == f64::INFINITY {
        abstain.fill(true);
        let bubbles = usize::from(k >= 6);
        return (Vec::new(), bubbles);
    }
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
                    // Intentionally bridge 1-bp posterior dips (including mildly negative matches)
                    // when the immediately following column is an aligned pair recovering above `min_post_prob`.
                    || (idx + 1 < k
                        && raw_cols[idx].0.is_some()
                        && raw_cols[idx].1.is_some()
                        && raw_cols[idx].2 > dip_floor_bits
                        && raw_cols[idx + 1].0.is_some()
                        && raw_cols[idx + 1].1.is_some()
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
    (kept_pairs, bubbles)
}

/// Run banded forward-backward MDL + BOCPD alignment on one gap reusing `scratch` buffers.
pub(crate) fn banded_mdl_gap_with_scratch(
    seq_x: &[u8],
    seq_y: &[u8],
    band_w: usize,
    min_block_bits: f64,
    params: &HmmParams,
    min_post_prob: f64,
    scratch: &mut DpScratch,
) -> GapResult {
    let n = seq_x.len();
    let m = seq_y.len();
    if n == 0 || m == 0 || !min_block_bits.is_finite() {
        return GapResult {
            kept_pairs: Vec::new(),
            bubbles: usize::from(n.max(m) >= 6),
            score_bits: 0.0,
        };
    }
    if let Some(fast) =
        super::dp_fast::try_diagonal_fast_path(seq_x, seq_y, min_block_bits, min_post_prob)
    {
        return fast;
    }
    if let Some(fast_vit) = super::dp_fast::try_banded_viterbi_fast_path(
        seq_x,
        seq_y,
        band_w,
        min_block_bits,
        min_post_prob,
        scratch,
    ) {
        return fast_vit;
    }
    if n >= 24 && m >= 24 && min_post_prob <= 0.90 {
        let (pref, suff) = super::dp_fast::matching_prefix_suffix_lens(seq_x, seq_y);
        if pref + suff >= 16 && (pref + suff) < n.min(m) {
            let mid_x = &seq_x[pref..(n - suff)];
            let mid_y = &seq_y[pref..(m - suff)];
            #[allow(clippy::cast_precision_loss)]
            let flank_bits = ((pref + suff) as f64) * 1.5;
            let adj_min_bits = (min_block_bits - flank_bits).max(0.0);
            let mid_res = banded_mdl_gap_core(
                mid_x,
                mid_y,
                band_w,
                adj_min_bits,
                params,
                min_post_prob,
                scratch,
            );
            if mid_res.kept_pairs.is_empty()
                || !mid_res.score_bits.is_finite()
                || mid_res.score_bits <= NEG_INF * 0.5
            {
                let mut kept_pairs = Vec::new();
                let mut score_bits = 0.0;
                #[allow(clippy::cast_precision_loss)]
                let pref_bits = (pref as f64) * 1.5;
                #[allow(clippy::cast_precision_loss)]
                let suff_bits = (suff as f64) * 1.5;
                if pref > 0 && pref_bits >= min_block_bits {
                    kept_pairs.extend((0..pref).map(|i| (i, i)));
                    score_bits += pref_bits;
                }
                if suff > 0 && suff_bits >= min_block_bits {
                    kept_pairs.extend((0..suff).map(|k| (n - suff + k, m - suff + k)));
                    score_bits += suff_bits;
                }
                return GapResult {
                    kept_pairs,
                    bubbles: mid_res.bubbles,
                    score_bits,
                };
            }
            let mut kept_pairs = Vec::with_capacity(pref + mid_res.kept_pairs.len() + suff);
            for i in 0..pref {
                kept_pairs.push((i, i));
            }
            for (mx, my) in mid_res.kept_pairs {
                kept_pairs.push((pref + mx, pref + my));
            }
            for k in 0..suff {
                kept_pairs.push((n - suff + k, m - suff + k));
            }
            return GapResult {
                kept_pairs,
                bubbles: mid_res.bubbles,
                score_bits: flank_bits + mid_res.score_bits,
            };
        }
    }
    banded_mdl_gap_core(
        seq_x,
        seq_y,
        band_w,
        min_block_bits,
        params,
        min_post_prob,
        scratch,
    )
}

#[allow(
    clippy::too_many_lines,
    clippy::similar_names,
    clippy::many_single_char_names,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss
)]
fn banded_mdl_gap_core(
    seq_x: &[u8],
    seq_y: &[u8],
    band_w: usize,
    min_block_bits: f64,
    params: &HmmParams,
    min_post_prob: f64,
    scratch: &mut DpScratch,
) -> GapResult {
    let n = seq_x.len();
    let m = seq_y.len();

    bocpd_gc_profile_into(
        seq_x,
        1.0 / 120.0,
        80,
        &mut scratch.bocpd,
        &mut scratch.gc_x,
    );
    bocpd_gc_profile_into(
        seq_y,
        1.0 / 120.0,
        80,
        &mut scratch.bocpd,
        &mut scratch.gc_y,
    );

    scratch.pi_x.clear();
    scratch
        .pi_x
        .extend(scratch.gc_x.iter().copied().map(base_probs));
    scratch.pi_y.clear();
    scratch
        .pi_y
        .extend(scratch.gc_y.iter().copied().map(base_probs));

    scratch.enc_x.clear();
    scratch
        .enc_x
        .extend(seq_x.iter().copied().map(base_to_idx_fast));
    scratch.enc_y.clear();
    scratch
        .enc_y
        .extend(seq_y.iter().copied().map(base_to_idx_fast));

    let NormalizedHmmParams {
        id_prior,
        ge_bits,
        p_open,
        p_ext,
        p_mm,
        p_im,
        ..
    } = params.normalized();
    let mismatch_ratio = 1.0 - id_prior;

    let Some(layout) = BandLayout::new(n, m, band_w) else {
        return GapResult {
            kept_pairs: Vec::new(),
            bubbles: usize::from(n.max(m) >= 6),
            score_bits: 0.0,
        };
    };
    let bw = layout.band_w;
    let stride = layout.stride;
    let total = layout.total;

    if !prepare_grid(&mut scratch.fm, total, 0.0)
        || !prepare_grid(&mut scratch.fix, total, 0.0)
        || !prepare_grid(&mut scratch.fiy, total, 0.0)
        || !prepare_grid(&mut scratch.sub_grid, total, 1.0)
    {
        return abstained_gap_result(usize::from(n.max(m) >= 6), 0.0);
    }
    scratch.scale_f.clear();
    scratch.scale_f.resize(n + 1, 1.0);

    if let Some(s0) = layout.slot(0, 0) {
        scratch.fm[s0] = 1.0;
    }

    let pi_x = &scratch.pi_x;
    let pi_y = &scratch.pi_y;
    let enc_x = &scratch.enc_x;
    let enc_y = &scratch.enc_y;

    let mut log2_scale_sum = 0.0_f64;

    // Scaled linear-space Forward pass (Durbin et al. 1998, Section 3.2):
    // zero exp()/ln()/log2() per cell; one log2() per row.
    for i in 0..=n {
        let center = band_center(i, m, n);
        let j_lo = center.saturating_sub(bw);
        let j_hi = (center + bw).min(m);
        let col_lo = (j_lo as isize - center as isize + bw as isize) as usize;
        let row_cur = i * stride;
        let mut row_max = 0.0_f64;

        if i == 0 {
            for (col, j) in (col_lo..).zip(j_lo..=j_hi) {
                if col > 0 && j > 0 {
                    let v = scratch.fm[row_cur + col - 1] * p_open
                        + scratch.fiy[row_cur + col - 1] * p_ext;
                    scratch.fiy[row_cur + col] = v;
                    if v > row_max {
                        row_max = v;
                    }
                } else if j == 0 {
                    row_max = row_max.max(scratch.fm[row_cur + col]);
                }
            }
        } else {
            let prev_center = band_center(i - 1, m, n);
            let shift = (center - prev_center) as isize;
            let row_prev = (i - 1) * stride;
            let ia = enc_x[i - 1];
            let pix = pi_x[i - 1];

            for (col, j) in (col_lo..).zip(j_lo..=j_hi) {
                let prev_col = col as isize + shift;
                let ix_v = if prev_col >= 0 && (prev_col as usize) < stride {
                    let pc = prev_col as usize;
                    scratch.fm[row_prev + pc] * p_open + scratch.fix[row_prev + pc] * p_ext
                } else {
                    0.0
                };
                scratch.fix[row_cur + col] = ix_v;
                if ix_v > row_max {
                    row_max = ix_v;
                }

                let iy_v = if j > 0 && col > 0 {
                    scratch.fm[row_cur + col - 1] * p_open + scratch.fiy[row_cur + col - 1] * p_ext
                } else {
                    0.0
                };
                scratch.fiy[row_cur + col] = iy_v;
                if iy_v > row_max {
                    row_max = iy_v;
                }

                if j > 0 {
                    let diag_col = prev_col - 1;
                    if diag_col >= 0 && (diag_col as usize) < stride {
                        let dc = diag_col as usize;
                        let er = emit_ratio_enc(
                            ia,
                            enc_y[j - 1],
                            pix,
                            pi_y[j - 1],
                            id_prior,
                            mismatch_ratio,
                        );
                        scratch.sub_grid[row_prev + dc] = er;
                        let m_v = (scratch.fm[row_prev + dc] * p_mm
                            + scratch.fix[row_prev + dc] * p_im
                            + scratch.fiy[row_prev + dc] * p_im)
                            * er;
                        scratch.fm[row_cur + col] = m_v;
                        if m_v > row_max {
                            row_max = m_v;
                        }
                    }
                }
            }
        }

        if row_max > 0.0 {
            scratch.scale_f[i] = row_max;
            let inv = 1.0 / row_max;
            let col_hi = (j_hi as isize - center as isize + bw as isize) as usize;
            for c in col_lo..=col_hi {
                scratch.fm[row_cur + c] *= inv;
                scratch.fix[row_cur + c] *= inv;
                scratch.fiy[row_cur + c] *= inv;
            }
            log2_scale_sum += row_max.log2();
        } else {
            return abstained_gap_result(1, NEG_INF);
        }
    }

    let z_scaled = layout.slot(n, m).map_or(0.0_f64, |snm| {
        scratch.fm[snm] + scratch.fix[snm] + scratch.fiy[snm]
    });
    if !z_scaled.is_finite() || z_scaled <= 0.0 {
        return abstained_gap_result(1, NEG_INF);
    }
    let score_bits = log2_scale_sum + z_scaled.log2();

    if !prepare_grid(&mut scratch.bm, total, 0.0)
        || !prepare_grid(&mut scratch.bix, total, 0.0)
        || !prepare_grid(&mut scratch.biy, total, 0.0)
    {
        return abstained_gap_result(usize::from(n.max(m) >= 6), 0.0);
    }
    if let Some(snm) = layout.slot(n, m) {
        scratch.bm[snm] = 1.0;
        scratch.bix[snm] = 1.0;
        scratch.biy[snm] = 1.0;
    }

    // Scaled linear-space Backward pass using the exact Forward row scales `scale_f[i + 1]`.
    for i in (0..=n).rev() {
        let center = band_center(i, m, n);
        let j_lo = center.saturating_sub(bw);
        let j_hi = (center + bw).min(m);
        let col_lo = (j_lo as isize - center as isize + bw as isize) as usize;
        let col_hi = (j_hi as isize - center as isize + bw as isize) as usize;
        let row_cur = i * stride;
        let has_next_row = i < n;
        let next_shift = if has_next_row {
            let next_center = band_center(i + 1, m, n);
            (next_center - center) as isize
        } else {
            0
        };
        let row_next = (i + 1) * stride;
        let inv_next_scale = if has_next_row {
            1.0 / scratch.scale_f[i + 1]
        } else {
            1.0
        };

        let mut col = col_hi;
        for j in (j_lo..=j_hi).rev() {
            if !(i == n && j == m) {
                let next_col = col as isize - next_shift;
                let m_next = if has_next_row
                    && j < m
                    && next_col + 1 >= 0
                    && ((next_col + 1) as usize) < stride
                {
                    let nc1 = (next_col + 1) as usize;
                    scratch.bm[row_next + nc1] * scratch.sub_grid[row_cur + col]
                } else {
                    0.0
                };
                let (ix_open, ix_ext) =
                    if has_next_row && next_col >= 0 && (next_col as usize) < stride {
                        let bix_v = scratch.bix[row_next + next_col as usize];
                        (bix_v * p_open, bix_v * p_ext)
                    } else {
                        (0.0, 0.0)
                    };
                let (iy_open, iy_ext) = if j < m && col + 1 < stride {
                    let biy_v = scratch.biy[row_cur + col + 1];
                    (biy_v * p_open, biy_v * p_ext)
                } else {
                    (0.0, 0.0)
                };
                scratch.bm[row_cur + col] = m_next * p_mm + ix_open + iy_open;
                scratch.bix[row_cur + col] = m_next * p_im + ix_ext;
                scratch.biy[row_cur + col] = m_next * p_im + iy_ext;
            }
            if col == 0 {
                break;
            }
            col -= 1;
        }
        if has_next_row {
            for c in col_lo..=col_hi {
                scratch.bm[row_cur + c] *= inv_next_scale;
                scratch.bix[row_cur + c] *= inv_next_scale;
                scratch.biy[row_cur + c] *= inv_next_scale;
            }
        }
    }

    let inv_z_scaled = 1.0 / z_scaled;
    // Pure Forward-Backward Maximum Expected Accuracy (MEA) decoding (Durbin et al. 1998):
    // reuse `fm` for MEA cumulative score and `bm` for posterior match probability.
    for i in 0..=n {
        let center = band_center(i, m, n);
        let j_lo = center.saturating_sub(bw);
        let j_hi = (center + bw).min(m);
        let col_lo = (j_lo as isize - center as isize + bw as isize) as usize;
        let row_cur = i * stride;
        let shift = if i > 0 {
            let prev_center = band_center(i - 1, m, n);
            (center - prev_center) as isize
        } else {
            0
        };
        let row_prev = i.saturating_sub(1) * stride;

        for (col, j) in (col_lo..).zip(j_lo..=j_hi) {
            if i == 0 && j == 0 {
                scratch.fm[row_cur + col] = 0.0;
                continue;
            }
            let post_m = if i > 0 && j > 0 {
                let p =
                    (scratch.fm[row_cur + col] * scratch.bm[row_cur + col] * inv_z_scaled).min(1.0);
                scratch.bm[row_cur + col] = p;
                p
            } else {
                0.0
            };
            let prev_col = col as isize + shift;
            let d = if i > 0 && j > 0 && prev_col > 0 && ((prev_col - 1) as usize) < stride {
                scratch.fm[row_prev + (prev_col - 1) as usize] + (post_m - 0.25)
            } else {
                NEG_INF
            };
            let u = if i > 0 && prev_col >= 0 && (prev_col as usize) < stride {
                scratch.fm[row_prev + prev_col as usize]
            } else {
                NEG_INF
            };
            let l = if j > 0 && col > 0 {
                scratch.fm[row_cur + col - 1]
            } else {
                NEG_INF
            };
            scratch.fm[row_cur + col] = d.max(u).max(l);
        }
    }

    scratch.raw_cols.clear();
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        let cur = layout.get(&scratch.fm, i, j);
        if cur <= NEG_INF * 0.5 {
            break;
        }
        let d = if i > 0 && j > 0 {
            layout.get(&scratch.fm, i - 1, j - 1) + (layout.get(&scratch.bm, i, j) - 0.25)
        } else {
            NEG_INF
        };
        let u = if i > 0 {
            layout.get(&scratch.fm, i - 1, j)
        } else {
            NEG_INF
        };
        if i > 0 && j > 0 && (cur - d).abs() < 1e-9 {
            let er = layout.get(&scratch.sub_grid, i - 1, j - 1);
            let sb = er.log2();
            let post_p = layout.get(&scratch.bm, i, j);
            scratch
                .raw_cols
                .push((Some(i - 1), Some(j - 1), sb, post_p));
            i -= 1;
            j -= 1;
        } else if i > 0 && (cur - u).abs() < 1e-9 {
            scratch.raw_cols.push((Some(i - 1), None, ge_bits, 0.0));
            i -= 1;
        } else if j > 0 {
            scratch.raw_cols.push((None, Some(j - 1), ge_bits, 0.0));
            j -= 1;
        } else {
            break;
        }
    }
    scratch.raw_cols.reverse();

    let dip_floor_bits = dip_floor_bits_for_identity(params.identity_prior);
    let (kept_pairs, bubbles) = filter_mea_pairs_and_count_bubbles(
        &scratch.raw_cols,
        &mut scratch.abstain,
        min_post_prob,
        min_block_bits,
        dip_floor_bits,
    );

    GapResult {
        kept_pairs,
        bubbles,
        score_bits,
    }
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

    #[test]
    fn empty_banded_gap_behaves_cleanly() {
        let empty_res = super::banded_mdl_gap(b"", b"ACGT", 8, 5.0);
        assert_eq!(empty_res.kept_pairs, [] as [(usize, usize); 0]);
        assert_eq!(empty_res.bubbles, 0);
    }

    #[test]
    fn dip_bridging_does_not_bridge_across_gap_when_min_post_prob_zero() {
        let mut abstain = Vec::new();
        let cols2 = vec![
            (Some(0usize), Some(0usize), 2.5, 0.9),
            (None, Some(1), -2.0, 0.0), // gap column (must NOT be bridged into)
            (Some(1), Some(2), 2.5, 0.9),
        ];
        let (kept2, _) =
            super::filter_mea_pairs_and_count_bubbles(&cols2, &mut abstain, 0.0, 2.0, -2.5);
        assert_eq!(kept2, vec![(0, 0), (1, 2)]);
        assert!(abstain[1], "gap column must remain abstained");
    }

    #[test]
    fn dip_floor_bits_scales_monotonically_with_identity() {
        let low = super::dip_floor_bits_for_identity(0.55);
        let mid = super::dip_floor_bits_for_identity(0.65);
        let high = super::dip_floor_bits_for_identity(0.90);
        assert!(low.is_finite() && mid.is_finite() && high.is_finite());
        assert!(low < 0.0 && mid < 0.0 && high < 0.0);
        assert!(
            low > mid && mid > high,
            "higher identity prior should penalize single-base mismatches more strongly: {low} > {mid} > {high}"
        );
    }

    #[test]
    fn diagonal_fast_path_recovers_near_identical_gap() {
        let x = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut y = *x;
        y[10] = b'T';
        let res = super::banded_mdl_gap(x, &y, 16, 10.0);
        assert_eq!(res.bubbles, 0);
        assert_eq!(res.kept_pairs.len(), x.len());
        assert_eq!(res.kept_pairs[0], (0, 0));
        assert_eq!(res.kept_pairs[x.len() - 1], (x.len() - 1, x.len() - 1));
    }

    #[test]
    fn prefix_suffix_trimming_handles_internal_indel() {
        let mut x = b"ACGTACGTACGTACGTACGTACGT".to_vec();
        x.extend_from_slice(b"GGG");
        x.extend_from_slice(b"TGCATGCATGCATGCATGCATGCA");
        let mut y = b"ACGTACGTACGTACGTACGTACGT".to_vec();
        y.extend_from_slice(b"TGCATGCATGCATGCATGCATGCA");
        let res = super::banded_mdl_gap(&x, &y, 16, 10.0);
        assert_eq!(res.kept_pairs.len(), y.len());
        assert_eq!(res.kept_pairs[0], (0, 0));
        assert_eq!(res.kept_pairs[y.len() - 1], (x.len() - 1, y.len() - 1));
    }

    #[test]
    fn banded_viterbi_fast_path_recovers_small_indel_with_flank_snps() {
        let mut x = b"AACTGACTGACTGACTGACTGACTGACTGACTG".to_vec();
        x.extend_from_slice(b"TTT");
        x.extend_from_slice(b"GATTACAGATTACAGATTACAGATTACAA");
        let mut y = b"ATCTGACTGACTGACTGACTGACTGACTGACTG".to_vec();
        y.extend_from_slice(b"GATTACAGATTACAGATTACAGATTACAT");
        let res = super::banded_mdl_gap(&x, &y, 16, 10.0);
        assert_eq!(res.bubbles, 0);
        assert_eq!(res.kept_pairs.len(), y.len());
    }
}
