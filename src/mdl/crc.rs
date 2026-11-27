//! Conformal Risk Control (CRC) null bit-savings calibration (Angelopoulos et al. 2024).

use crate::bocpd::MIN_POST_PROB;
use crate::mdl::dp::{DpScratch, HmmParams, banded_mdl_gap_core_with_max_bits};

/// Default CRC target false-positive risk level `alpha`.
pub const DEFAULT_CRC_ALPHA: f64 = 0.07;

/// Default calibration sample count `n_cal`. At `n_cal = 49` and `alpha = 0.07`,
/// `(n_gt + 1) / (n_cal + 1) <= 0.07` holds at `n_gt = 2` (`3/50 = 0.06`), selecting a
/// genuine 3rd-largest order statistic rather than degenerating to a single-sample maximum.
pub const DEFAULT_CRC_N_CAL: usize = 49;

/// Diagnostic summary returned by [`calibrate_null_crc_report`].
#[derive(Debug, Clone, Copy)]
pub struct CrcCalibrationReport {
    /// Selected bit-savings threshold `λ̂ >= base_floor_bits`.
    pub threshold_bits: f64,
    /// Input search-space / base floor before empirical calibration.
    pub base_floor_bits: f64,
    /// Number of valid reversed biological decoy windows scored.
    pub n_cal_bio: usize,
    /// Median uncensored biological decoy `block_bits` (0.0 if synthetic fallback was used).
    pub median_null_bits: f64,
    /// Maximum uncensored biological decoy `block_bits` (0.0 if synthetic fallback was used).
    pub max_null_bits: f64,
    /// True when sequences were too short/N-heavy for biological decoys and synthetic GC decoys were used.
    pub used_synthetic_fallback: bool,
}

use crate::splitmix64;

/// Calibrate a Conformal Risk Control (Angelopoulos et al. 2024) bit-savings
/// threshold `λ̂ >= base_floor_bits` from `n_cal` null/decoy block bit scores
/// such that the finite-sample false-positive bound `(count(score > λ̂) + 1) / (n_cal + 1) <= alpha`
/// holds (adding `+1e-6` so `dp.rs`'s `>=` emit gate rejects scores at the selected order statistic).
#[must_use]
pub fn calibrate_crc_threshold(null_block_bits: &[f64], alpha: f64, base_floor_bits: f64) -> f64 {
    if null_block_bits.is_empty() || !alpha.is_finite() || alpha <= 0.0 {
        return base_floor_bits;
    }
    let mut scores: Vec<f64> = null_block_bits
        .iter()
        .copied()
        .filter(|s| s.is_finite())
        .collect();
    if scores.is_empty() {
        return base_floor_bits;
    }
    scores.sort_by(f64::total_cmp);
    let n = scores.len();
    #[allow(clippy::cast_precision_loss)]
    let denom = (n + 1) as f64;
    for (idx, &s) in scores.iter().enumerate() {
        let candidate = s.max(base_floor_bits);
        // Count calibration scores strictly above `s` (matched to `>=` emit gate via `+1e-6`).
        let n_gt = n - idx - 1;
        #[allow(clippy::cast_precision_loss)]
        let upper_risk = (n_gt as f64 + 1.0) / denom;
        if upper_risk <= alpha {
            return candidate + if s >= base_floor_bits { 1e-6 } else { 0.0 };
        }
    }
    scores[n - 1].max(base_floor_bits) + 1e-6
}

/// Calibrate `λ̂` via `n_cal` deterministic synthetic decoy windows drawn at composition `gc_frac`
/// under `params`, cycling across `[64, 128, 192]` bp windows while never dropping below
/// `base_floor_bits`.
#[must_use]
pub fn calibrate_null_crc_for_gc(
    gc_frac: f64,
    n_cal: usize,
    alpha: f64,
    base_floor_bits: f64,
    params: &HmmParams,
) -> f64 {
    calibrate_null_crc_for_gc_with_min_post(
        gc_frac,
        n_cal,
        alpha,
        base_floor_bits,
        params,
        MIN_POST_PROB,
    )
}

fn calibrate_null_crc_for_gc_with_min_post(
    gc_frac: f64,
    n_cal: usize,
    alpha: f64,
    base_floor_bits: f64,
    params: &HmmParams,
    min_post_prob: f64,
) -> f64 {
    const SYNTH_LENS: [usize; 3] = [64, 128, 192];
    if n_cal == 0 {
        return base_floor_bits;
    }
    let gc = if gc_frac.is_finite() { gc_frac } else { 0.5 }.clamp(0.05, 0.95);
    let mut scores = Vec::with_capacity(n_cal);
    let mut x = Vec::with_capacity(192);
    let mut y = Vec::with_capacity(192);
    let mut scratch = DpScratch::default();
    for draw in 0..n_cal {
        let win_len = SYNTH_LENS[draw % SYNTH_LENS.len()];
        x.resize(win_len, b'A');
        y.resize(win_len, b'A');
        let mut state =
            0x9E37_79B9_7F4A_7C15_u64 ^ ((draw as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9));
        for buf in [&mut x, &mut y] {
            for b in buf.iter_mut() {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                #[allow(clippy::cast_precision_loss)]
                let u = ((state >> 11) as f64) / ((1_u64 << 53) as f64);
                *b = if u < gc * 0.5 {
                    b'G'
                } else if u < gc {
                    b'C'
                } else if u < gc + (1.0 - gc) * 0.5 {
                    b'A'
                } else {
                    b'T'
                };
            }
        }
        let band_w = (win_len / 2).clamp(16, 64);
        let (_, max_block_bits) = banded_mdl_gap_core_with_max_bits(
            &x,
            &y,
            band_w,
            0.0,
            params,
            min_post_prob,
            &mut scratch,
        );
        scores.push(max_block_bits.max(0.0));
    }
    calibrate_crc_threshold(&scores, alpha, base_floor_bits)
}

/// Estimate pooled GC fraction by sampling up to 4,096 positions uniformly strided across the full
/// length of both `seq_x` and `seq_y` (avoiding bias from a single terminal isochore).
#[must_use]
fn estimate_sampled_gc(seq_x: &[u8], seq_y: &[u8]) -> f64 {
    const MAX_SAMPLES_PER_SEQ: usize = 4096;
    let mut gc = 0usize;
    let mut acgt = 0usize;
    for seq in [seq_x, seq_y] {
        if seq.is_empty() {
            continue;
        }
        let step = (seq.len() / MAX_SAMPLES_PER_SEQ).max(1);
        for &b in seq.iter().step_by(step).take(MAX_SAMPLES_PER_SEQ) {
            match b.to_ascii_uppercase() {
                b'G' | b'C' => {
                    gc += 1;
                    acgt += 1;
                }
                b'A' | b'T' => acgt += 1,
                _ => {}
            }
        }
    }
    #[allow(clippy::cast_precision_loss)]
    if acgt > 0 {
        gc as f64 / acgt as f64
    } else {
        0.5
    }
}

/// Calibrate `λ̂` using biological reversed negative windows whose lengths `(lx, ly)` and band
/// width are matched to the empirical `target_windows` (e.g. live inter-LCB `dp_jobs` spans).
#[allow(clippy::too_many_arguments)]
#[must_use]
pub fn calibrate_null_crc_matched(
    seq_x: &[u8],
    seq_y: &[u8],
    target_windows: &[(usize, usize)],
    band_w: usize,
    n_cal: usize,
    alpha: f64,
    base_floor_bits: f64,
    params: &HmmParams,
    min_post_prob: f64,
) -> f64 {
    calibrate_null_crc_report(
        seq_x,
        seq_y,
        target_windows,
        band_w,
        n_cal,
        alpha,
        base_floor_bits,
        params,
        min_post_prob,
    )
    .threshold_bits
}

/// Run empirical biological negative CRC calibration and return full diagnostic statistics.
#[allow(clippy::too_many_arguments)]
fn collect_bio_decoy_scores(
    seq_x: &[u8],
    seq_y: &[u8],
    target_windows: &[(usize, usize)],
    band_w: usize,
    n_cal: usize,
    params: &HmmParams,
    min_post_prob: f64,
) -> Vec<f64> {
    const MIN_WIN: usize = 64;
    const MAX_WIN: usize = 320;
    const MULTI_SCALE: [usize; 4] = [64, 128, 192, 256];

    let max_avail_x = (seq_x.len() / 2).clamp(MIN_WIN, MAX_WIN);
    let max_avail_y = (seq_y.len() / 2).clamp(MIN_WIN, MAX_WIN);
    let seed_mix = (seq_x.len() as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
        ^ (seq_y.len() as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9);

    let mut scores = Vec::with_capacity(n_cal);
    let mut y_rev = Vec::with_capacity(MAX_WIN);
    let mut scratch = DpScratch::default();
    let max_tries = n_cal.saturating_mul(6).max(32);

    for attempt in 0..max_tries {
        if scores.len() >= n_cal {
            break;
        }
        let draw_idx = scores.len();
        let (lx, ly) = if target_windows.is_empty() {
            let w = MULTI_SCALE[draw_idx % MULTI_SCALE.len()];
            (w.min(max_avail_x), w.min(max_avail_y))
        } else {
            let (tx, ty) = target_windows[draw_idx % target_windows.len()];
            (
                tx.clamp(MIN_WIN, max_avail_x),
                ty.clamp(MIN_WIN, max_avail_y),
            )
        };
        let max_x = seq_x.len() - lx;
        let max_y = seq_y.len() - ly;
        let hx = splitmix64(seed_mix ^ ((attempt as u64).wrapping_mul(2) + 1));
        let hy = splitmix64(seed_mix ^ ((attempt as u64).wrapping_mul(2) + 2));
        #[allow(clippy::cast_possible_truncation)]
        let sx = (hx as usize) % (max_x + 1);
        #[allow(clippy::cast_possible_truncation)]
        let mut sy = (hy as usize) % (max_y + 1);
        // Keep biological decoys disjoint/non-collinear on equal-length self-comparisons.
        if seq_x.len() == seq_y.len() && sx.abs_diff(sy) < lx.max(ly) {
            let sep = lx.max(ly);
            sy = if sx >= sep {
                sx - sep
            } else {
                (sx + sep).min(max_y)
            };
        }

        let wx = &seq_x[sx..sx + lx];
        let wy = &seq_y[sy..sy + ly];
        let acgt_x = wx
            .iter()
            .filter(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .count();
        let acgt_y = wy
            .iter()
            .filter(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .count();
        if acgt_x * 2 < lx || acgt_y * 2 < ly {
            continue;
        }

        y_rev.clear();
        y_rev.extend(wy.iter().rev().copied());
        let eff_band = band_w.min(lx.max(ly)).clamp(16, 96);
        let (_, max_block_bits) = banded_mdl_gap_core_with_max_bits(
            wx,
            &y_rev,
            eff_band,
            0.0,
            params,
            min_post_prob,
            &mut scratch,
        );
        scores.push(max_block_bits.max(0.0));
    }
    scores
}

/// Calibrate `λ̂` from biological negative windows matched to empirical gap spans (`target_windows`)
/// and return full null-score diagnostics alongside the final threshold.
///
/// Falls back to synthetic GC-matched decoys only when sequences are too short or contain fewer
/// than `n_cal` valid ACGT windows.
#[allow(clippy::too_many_arguments)]
#[must_use]
pub fn calibrate_null_crc_report(
    seq_x: &[u8],
    seq_y: &[u8],
    target_windows: &[(usize, usize)],
    band_w: usize,
    n_cal: usize,
    alpha: f64,
    base_floor_bits: f64,
    params: &HmmParams,
    min_post_prob: f64,
) -> CrcCalibrationReport {
    const MIN_WIN: usize = 64;

    if n_cal == 0 {
        return CrcCalibrationReport {
            threshold_bits: base_floor_bits,
            base_floor_bits,
            n_cal_bio: 0,
            median_null_bits: 0.0,
            max_null_bits: 0.0,
            used_synthetic_fallback: false,
        };
    }

    let fallback = |n_bio: usize| {
        let obs_gc = estimate_sampled_gc(seq_x, seq_y);
        let gc_crc = calibrate_null_crc_for_gc_with_min_post(
            obs_gc,
            n_cal,
            alpha,
            base_floor_bits,
            params,
            min_post_prob,
        );
        CrcCalibrationReport {
            threshold_bits: gc_crc,
            base_floor_bits,
            n_cal_bio: n_bio,
            median_null_bits: 0.0,
            max_null_bits: 0.0,
            used_synthetic_fallback: true,
        }
    };

    if seq_x.len() < MIN_WIN * 2 || seq_y.len() < MIN_WIN * 2 {
        return fallback(0);
    }

    let scores = collect_bio_decoy_scores(
        seq_x,
        seq_y,
        target_windows,
        band_w,
        n_cal,
        params,
        min_post_prob,
    );
    if scores.len() < n_cal {
        return fallback(scores.len());
    }

    let mut sorted = scores.clone();
    sorted.sort_by(f64::total_cmp);
    let median_null_bits = sorted[sorted.len() / 2];
    let max_null_bits = sorted[sorted.len() - 1];
    let bio_crc = calibrate_crc_threshold(&scores, alpha, base_floor_bits);

    CrcCalibrationReport {
        threshold_bits: bio_crc,
        base_floor_bits,
        n_cal_bio: scores.len(),
        median_null_bits,
        max_null_bits,
        used_synthetic_fallback: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crc_threshold_respects_n_cal_15_and_alpha_0_07() {
        let scores = vec![
            1.2, 2.5, 3.1, 4.0, 4.5, 5.2, 6.0, 6.8, 7.1, 8.0, 8.5, 9.2, 9.8, 11.4, 12.6,
        ];
        // For n_cal = 15, alpha = 0.07: (0 + 1) / 16 = 0.0625 <= 0.07, so lambda_hat > 12.6
        let lam = calibrate_crc_threshold(&scores, 0.07, 10.0);
        assert!(lam > 12.6 && lam < 12.61, "got {lam}");

        let low_scores = vec![1.0; 15];
        assert!((calibrate_crc_threshold(&low_scores, 0.07, 10.0) - 10.0).abs() < 1e-9);

        let params = HmmParams::default();
        let gc_lam = calibrate_null_crc_for_gc(0.18, 15, 0.07, 10.0, &params);
        assert!(gc_lam >= 10.0 && gc_lam.is_finite());
    }

    #[test]
    fn default_n_cal_49_selects_third_largest_order_statistic() {
        // 49 sorted scores 1.0..=49.0 at alpha = 0.07:
        // idx 46 (val = 47.0): n_gt = 2, (2 + 1) / 50 = 0.06 <= 0.07 -> selects 47.0 + 1e-6
        // rather than degenerating to the single maximum (49.0).
        let scores: Vec<f64> = (1..=u32::try_from(DEFAULT_CRC_N_CAL).expect("u32"))
            .map(f64::from)
            .collect();
        let lam = calibrate_crc_threshold(&scores, DEFAULT_CRC_ALPHA, 10.0);
        assert!(
            (lam - 47.000_001).abs() < 1e-9,
            "expected 3rd-largest order statistic 47.000001, got {lam}"
        );
    }

    #[test]
    fn biological_crc_report_calibrates_multiscale_and_matched_windows() {
        let mut state = 0x1234_5678_9ABC_DEF0_u64;
        let bases = *b"ACGT";
        let seq_x: Vec<u8> = (0..4_000)
            .map(|_| {
                state = splitmix64(state);
                bases[(state & 3) as usize]
            })
            .collect();
        let seq_y: Vec<u8> = (0..4_000)
            .map(|_| {
                state = splitmix64(state);
                bases[(state & 3) as usize]
            })
            .collect();
        let params = HmmParams::default();
        let rep = calibrate_null_crc_report(
            &seq_x,
            &seq_y,
            &[(90, 95), (180, 175), (240, 260)],
            180,
            DEFAULT_CRC_N_CAL,
            DEFAULT_CRC_ALPHA,
            10.0,
            &params,
            MIN_POST_PROB,
        );
        assert!(!rep.used_synthetic_fallback);
        assert_eq!(rep.n_cal_bio, DEFAULT_CRC_N_CAL);
        assert!(rep.threshold_bits >= 10.0 && rep.threshold_bits.is_finite());
        assert!(rep.max_null_bits >= rep.median_null_bits);
    }

    #[test]
    fn bio_crc_does_not_censor_palindromic_repeat_decoys() {
        // C1 regression: balanced palindromic repeat windows (exact self-matches under reversal)
        // must be scored by exact Pair-HMM DP rather than censored to 0.0 by heuristic fast paths.
        let motif = b"ACGTACGTACGTACGTTGCATGCATGCATGCA";
        let rep_seq = motif.repeat(20);
        let params = HmmParams::default();
        let rep = calibrate_null_crc_report(
            &rep_seq,
            &rep_seq,
            &[(64, 64)],
            64,
            15,
            DEFAULT_CRC_ALPHA,
            10.0,
            &params,
            MIN_POST_PROB,
        );
        assert!(!rep.used_synthetic_fallback);
        assert!(
            rep.median_null_bits > 20.0 && rep.threshold_bits > 20.0,
            "palindromic decoys must not be censored to 0.0: {rep:?}"
        );
    }
}
