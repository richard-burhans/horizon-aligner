//! Banded MDL+BOCPD dynamic programming for inter-anchor gap alignment.

pub(crate) use super::dp_params::{
    BandLayout, DpScratch, NormalizedHmmParams, band_center, dip_floor_bits_for_identity,
    prepare_grid,
};
pub use super::dp_params::{GapPosteriorTrace, GapResult, HmmParams, MeaTraceColumn};
use crate::bocpd::{MIN_POST_PROB, NEG_INF, base_probs, base_to_idx_fast, bocpd_gc_profile_into};

fn abstained_gap_result(bubbles: usize, score_bits: f64) -> GapResult {
    GapResult {
        kept_pairs: Vec::new(),
        bubbles,
        score_bits,
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
    banded_mdl_gap_with_min_post(
        seq_x,
        seq_y,
        band_w,
        min_block_bits,
        &HmmParams::default(),
        MIN_POST_PROB,
    )
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

/// Run exact Forward-Backward Pair-HMM DP (bypassing heuristic fast paths) and return
/// the full in-band posterior grid alongside the MEA trace and abstention mask.
#[must_use]
pub fn banded_mdl_gap_posterior_trace(
    seq_x: &[u8],
    seq_y: &[u8],
    band_w: usize,
    min_block_bits: f64,
    params: &HmmParams,
    min_post_prob: f64,
) -> GapPosteriorTrace {
    let n = seq_x.len();
    let m = seq_y.len();
    let mut scratch = DpScratch::default();
    let res = banded_mdl_gap_core(
        seq_x,
        seq_y,
        band_w,
        min_block_bits,
        params,
        min_post_prob,
        &mut scratch,
    );
    let mut posterior_cells = Vec::new();
    let mut eff_bw = 0usize;
    if let Some(layout) = BandLayout::new(n, m, band_w) {
        eff_bw = layout.band_w;
        if scratch.bm.len() == layout.total {
            for i in 1..=n {
                let center = layout.center(i);
                let j_lo = center.saturating_sub(layout.band_w).max(1);
                let j_hi = (center + layout.band_w).min(m);
                for j in j_lo..=j_hi {
                    if let Some(slot) = layout.slot(i, j) {
                        posterior_cells.push((i, j, scratch.bm[slot]));
                    }
                }
            }
        }
    }
    let mea_cols = scratch
        .raw_cols
        .iter()
        .enumerate()
        .map(|(k, &(xi, yj, bits, post))| {
            let abst = scratch.abstain.get(k).copied().unwrap_or(false);
            (xi, yj, bits, post, abst)
        })
        .collect();
    GapPosteriorTrace {
        n,
        m,
        band_w: eff_bw,
        score_bits: res.score_bits,
        log2_zb: scratch.last_log2_zb,
        gc_x: scratch.gc_x,
        gc_y: scratch.gc_y,
        posterior_cells,
        mea_cols,
        kept_pairs: res.kept_pairs,
        bubbles: res.bubbles,
        max_block_bits: scratch.last_max_block_bits,
    }
}

pub(super) fn filter_mea_pairs_and_count_bubbles<X: Copy>(
    raw_cols: &[(Option<X>, Option<usize>, f64, f64)],
    abstain: &mut Vec<bool>,
    min_post_prob: f64,
    min_block_bits: f64,
    dip_floor_bits: f64,
) -> (Vec<(X, usize)>, usize) {
    let (kept, bubbles, _) = filter_mea_pairs_with_max_bits(
        raw_cols,
        abstain,
        min_post_prob,
        min_block_bits,
        dip_floor_bits,
    );
    (kept, bubbles)
}

fn filter_mea_pairs_with_max_bits<X: Copy>(
    raw_cols: &[(Option<X>, Option<usize>, f64, f64)],
    abstain: &mut Vec<bool>,
    min_post_prob: f64,
    min_block_bits: f64,
    dip_floor_bits: f64,
) -> (Vec<(X, usize)>, usize, f64) {
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
        return (Vec::new(), bubbles, 0.0);
    }
    let mut max_kept_block_bits = 0.0_f64;
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
            } else {
                max_kept_block_bits = max_kept_block_bits.max(block_bits);
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
    (kept_pairs, bubbles, max_kept_block_bits)
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
        super::dp_fast::try_diagonal_fast_path(seq_x, seq_y, min_block_bits, params, min_post_prob)
    {
        return fast;
    }
    if let Some(fast_vit) = super::dp_fast::try_banded_viterbi_fast_path(
        seq_x,
        seq_y,
        band_w,
        min_block_bits,
        params,
        min_post_prob,
        scratch,
    ) {
        return fast_vit;
    }
    if let Some(trimmed) = try_prefix_suffix_trimmed_gap(
        seq_x,
        seq_y,
        band_w,
        min_block_bits,
        params,
        min_post_prob,
        scratch,
    ) {
        return trimmed;
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

fn try_prefix_suffix_trimmed_gap(
    seq_x: &[u8],
    seq_y: &[u8],
    band_w: usize,
    min_block_bits: f64,
    params: &HmmParams,
    min_post_prob: f64,
    scratch: &mut DpScratch,
) -> Option<GapResult> {
    let n = seq_x.len();
    let m = seq_y.len();
    if n < 24
        || m < 24
        || min_post_prob > 0.90
        || !min_block_bits.is_finite()
        || super::dp_fast::is_fast_path_eligible(seq_x, seq_y, params).is_none()
    {
        return None;
    }
    let (pref, suff) = super::dp_fast::matching_prefix_suffix_lens(seq_x, seq_y);
    if pref + suff < 16 || (pref + suff) >= n.min(m) {
        return None;
    }
    let norm = params.normalized();
    let stay_bits = norm.mm_bits;

    let mid_x = &seq_x[pref..(n - suff)];
    let mid_y = &seq_y[pref..(m - suff)];
    // Run exact Forward-Backward MEA on the untrimmed interior, then stitch the exact-match
    // prefix/suffix columns scored with exact local BOCPD GC profiles (preventing whole-gap GC
    // leakage on homopolymer flanks) and apply `min_block_bits` once across the complete trace.
    let mid_res = banded_mdl_gap_core(mid_x, mid_y, band_w, 0.0, params, min_post_prob, scratch);
    let mid_valid = mid_res.score_bits.is_finite() && mid_res.score_bits > NEG_INF * 0.5;

    let mut stitched: Vec<(Option<usize>, Option<usize>, f64, f64)> =
        Vec::with_capacity(pref + scratch.raw_cols.len() + suff + 1);
    let mut flank_net_bits = 0.0_f64;
    if pref > 0 {
        bocpd_gc_profile_into(
            &seq_x[..pref],
            1.0 / 120.0,
            80,
            &mut scratch.bocpd,
            &mut scratch.gc_x,
        );
        for (i, (&b, &gc)) in seq_x[..pref].iter().zip(scratch.gc_x.iter()).enumerate() {
            let pi = base_probs(gc);
            let sb = crate::bocpd::joint_pair_bits(b, b, pi, pi, norm.id_prior);
            flank_net_bits += sb + stay_bits;
            stitched.push((Some(i), Some(i), sb, 1.0));
        }
    }
    if mid_valid && !scratch.raw_cols.is_empty() {
        for &(cx, cy, sb, pp) in &scratch.raw_cols {
            stitched.push((cx.map(|x| pref + x), cy.map(|y| pref + y), sb, pp));
        }
    } else {
        for ix in pref..(n - suff) {
            stitched.push((Some(ix), None, norm.ge_bits, 0.0));
        }
        for iy in pref..(m - suff) {
            stitched.push((None, Some(iy), norm.ge_bits, 0.0));
        }
    }
    if suff > 0 {
        bocpd_gc_profile_into(
            &seq_x[(n - suff)..],
            1.0 / 120.0,
            80,
            &mut scratch.bocpd,
            &mut scratch.gc_x,
        );
        for (k, (&b, &gc)) in seq_x[(n - suff)..]
            .iter()
            .zip(scratch.gc_x.iter())
            .enumerate()
        {
            let pi = base_probs(gc);
            let sb = crate::bocpd::joint_pair_bits(b, b, pi, pi, norm.id_prior);
            flank_net_bits += sb + stay_bits;
            stitched.push((Some(n - suff + k), Some(m - suff + k), sb, 1.0));
        }
    }

    let dip_floor_bits = dip_floor_bits_for_identity(params.identity_prior);
    let (kept_pairs, bubbles, max_block_bits) = filter_mea_pairs_with_max_bits(
        &stitched,
        &mut scratch.abstain,
        min_post_prob,
        min_block_bits,
        dip_floor_bits,
    );
    scratch.last_max_block_bits = max_block_bits;
    let total_score = if mid_valid {
        mid_res.score_bits + flank_net_bits
    } else {
        flank_net_bits
    };
    Some(GapResult {
        kept_pairs,
        bubbles,
        score_bits: total_score,
    })
}

/// Run exact Forward-Backward DP (bypassing fast paths) and return `(GapResult, max_block_bits)`.
/// Used by CRC calibration so biological decoys are never censored by heuristic fast paths.
pub(crate) fn banded_mdl_gap_core_with_max_bits(
    seq_x: &[u8],
    seq_y: &[u8],
    band_w: usize,
    min_block_bits: f64,
    params: &HmmParams,
    min_post_prob: f64,
    scratch: &mut DpScratch,
) -> (GapResult, f64) {
    scratch.last_max_block_bits = 0.0;
    if seq_x.is_empty() || seq_y.is_empty() || !min_block_bits.is_finite() {
        return (
            GapResult {
                kept_pairs: Vec::new(),
                bubbles: usize::from(seq_x.len().max(seq_y.len()) >= 6),
                score_bits: 0.0,
            },
            0.0,
        );
    }
    let res = banded_mdl_gap_core(
        seq_x,
        seq_y,
        band_w,
        min_block_bits,
        params,
        min_post_prob,
        scratch,
    );
    (res, scratch.last_max_block_bits)
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
    scratch.raw_cols.clear();
    scratch.last_max_block_bits = 0.0;
    scratch.last_log2_zb = NEG_INF;
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
    scratch.log2_scale_f_prefix.clear();
    scratch.log2_scale_f_prefix.resize(n + 1, 0.0);
    scratch.log2_scale_b_suffix.clear();
    scratch.log2_scale_b_suffix.resize(n + 1, 0.0);

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
            let inv = 1.0 / row_max;
            let col_hi = (j_hi as isize - center as isize + bw as isize) as usize;
            for c in col_lo..=col_hi {
                scratch.fm[row_cur + c] *= inv;
                scratch.fix[row_cur + c] *= inv;
                scratch.fiy[row_cur + c] *= inv;
            }
            log2_scale_sum += row_max.log2();
            scratch.log2_scale_f_prefix[i] = log2_scale_sum;
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

    // Independently row-scaled linear-space Backward pass so `bm`/`bix`/`biy` stay in [0, 1]
    // and can never overflow to `+inf` even when Forward and Backward peak on different diagonals.
    let mut log2_b_sum = 0.0_f64;
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
        let mut row_max_b = 0.0_f64;

        let mut col = col_hi;
        for j in (j_lo..=j_hi).rev() {
            if i == n && j == m {
                row_max_b = row_max_b.max(1.0);
            } else {
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
                let bm_v = m_next * p_mm + ix_open + iy_open;
                let bix_v = m_next * p_im + ix_ext;
                let biy_v = m_next * p_im + iy_ext;
                scratch.bm[row_cur + col] = bm_v;
                scratch.bix[row_cur + col] = bix_v;
                scratch.biy[row_cur + col] = biy_v;
                row_max_b = row_max_b.max(bm_v).max(bix_v).max(biy_v);
            }
            if col == 0 {
                break;
            }
            col -= 1;
        }
        if row_max_b > 0.0 && row_max_b.is_finite() {
            let inv_b = 1.0 / row_max_b;
            for c in col_lo..=col_hi {
                scratch.bm[row_cur + c] *= inv_b;
                scratch.bix[row_cur + c] *= inv_b;
                scratch.biy[row_cur + c] *= inv_b;
            }
            log2_b_sum += row_max_b.log2();
        }
        scratch.log2_scale_b_suffix[i] = log2_b_sum;
    }
    scratch.last_log2_zb = layout.slot(0, 0).map_or(NEG_INF, |s0| {
        let bm0 = scratch.bm[s0];
        if bm0 > 0.0 && bm0.is_finite() {
            scratch.log2_scale_b_suffix[0] + bm0.log2()
        } else {
            NEG_INF
        }
    });

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
        let delta_log2 =
            scratch.log2_scale_f_prefix[i] + scratch.log2_scale_b_suffix[i] - score_bits;
        let use_linear_row_factor = delta_log2.is_finite() && delta_log2.abs() < 500.0;
        let row_factor = if use_linear_row_factor {
            delta_log2.exp2()
        } else {
            0.0
        };

        for (col, j) in (col_lo..).zip(j_lo..=j_hi) {
            if i == 0 && j == 0 {
                scratch.fm[row_cur + col] = 0.0;
                continue;
            }
            let post_m = if i > 0 && j > 0 {
                let fm_v = scratch.fm[row_cur + col];
                let bm_v = scratch.bm[row_cur + col];
                let raw = if fm_v > 0.0 && bm_v > 0.0 {
                    if use_linear_row_factor {
                        fm_v * bm_v * row_factor
                    } else if delta_log2.is_finite() {
                        (fm_v.log2() + bm_v.log2() + delta_log2).exp2()
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                let p = if raw.is_finite() {
                    raw.clamp(0.0, 1.0)
                } else {
                    0.0
                };
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
    let (kept_pairs, bubbles, max_block_bits) = filter_mea_pairs_with_max_bits(
        &scratch.raw_cols,
        &mut scratch.abstain,
        min_post_prob,
        min_block_bits,
        dip_floor_bits,
    );
    scratch.last_max_block_bits = max_block_bits;

    GapResult {
        kept_pairs,
        bubbles,
        score_bits,
    }
}

#[cfg(test)]
#[path = "tests_dp.rs"]
mod tests;
