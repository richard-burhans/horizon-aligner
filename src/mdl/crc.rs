//! Conformal Risk Control (CRC) null bit-savings calibration (Angelopoulos et al. 2024).

use crate::bocpd::MIN_POST_PROB;
use crate::mdl::dp::{DpScratch, HmmParams, banded_mdl_gap_with_scratch};

/// Calibrate a Conformal Risk Control (Angelopoulos et al. 2024) bit-savings
/// threshold `λ̂ >= base_floor_bits` from `n_cal` null/decoy block bit scores
/// such that the finite-sample false-positive bound `(count(score >= λ̂) + 1) / (n_cal + 1) <= alpha` holds.
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
        // Number of calibration scores >= candidate (strictly speaking, > s if we set λ̂ = s + eps,
        // or >= s if we include ties). Count scores strictly above s if we add a tiny eps.
        let n_ge = n - idx - 1;
        #[allow(clippy::cast_precision_loss)]
        let upper_risk = (n_ge as f64 + 1.0) / denom;
        if upper_risk <= alpha {
            return candidate + if s >= base_floor_bits { 1e-6 } else { 0.0 };
        }
    }
    scores[n - 1].max(base_floor_bits) + 1e-6
}

/// Calibrate `λ̂` via `n_cal` deterministic decoy windows drawn at composition `gc_frac`
/// under `params`, ensuring `(n_cal * R_hat + 1) / (n_cal + 1) <= alpha` while never
/// dropping below `base_floor_bits`.
#[must_use]
pub fn calibrate_null_crc_for_gc(
    gc_frac: f64,
    n_cal: usize,
    alpha: f64,
    base_floor_bits: f64,
    params: &HmmParams,
) -> f64 {
    if n_cal == 0 {
        return base_floor_bits;
    }
    let gc = if gc_frac.is_finite() { gc_frac } else { 0.5 }.clamp(0.05, 0.95);
    let win_len = 64usize;
    let mut scores = Vec::with_capacity(n_cal);
    let mut x = vec![b'A'; win_len];
    let mut y = vec![b'A'; win_len];
    let mut scratch = DpScratch::default();
    for draw in 0..n_cal {
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
        let res = banded_mdl_gap_with_scratch(
            &x,
            &y,
            16,
            base_floor_bits,
            params,
            MIN_POST_PROB,
            &mut scratch,
        );
        let null_bits = if res.kept_pairs.is_empty() {
            0.0
        } else {
            res.score_bits.max(0.0)
        };
        scores.push(null_bits);
    }
    calibrate_crc_threshold(&scores, alpha, base_floor_bits)
}

#[must_use]
pub(crate) fn estimate_sampled_gc(seq_x: &[u8], seq_y: &[u8]) -> f64 {
    let mut gc = 0usize;
    let mut acgt = 0usize;
    for &b in seq_x.iter().take(4096).chain(seq_y.iter().take(4096)) {
        match b.to_ascii_uppercase() {
            b'G' | b'C' => {
                gc += 1;
                acgt += 1;
            }
            b'A' | b'T' => acgt += 1,
            _ => {}
        }
    }
    #[allow(clippy::cast_precision_loss)]
    if acgt > 0 {
        gc as f64 / acgt as f64
    } else {
        0.5
    }
}

/// Calibrate `λ̂` from `n_cal` empirical biological negative windows drawn from
/// `seq_x` and reversed non-collinear windows of `seq_y`, falling back to
/// [`calibrate_null_crc_for_gc`] when fewer than `n_cal` valid ACGT windows exist.
#[must_use]
pub fn calibrate_null_crc_from_seqs(
    seq_x: &[u8],
    seq_y: &[u8],
    n_cal: usize,
    alpha: f64,
    base_floor_bits: f64,
    params: &HmmParams,
) -> f64 {
    let obs_gc = estimate_sampled_gc(seq_x, seq_y);
    let gc_crc = calibrate_null_crc_for_gc(obs_gc, n_cal, alpha, base_floor_bits, params);
    calibrate_null_crc_from_seqs_with_gc_floor(
        seq_x,
        seq_y,
        n_cal,
        alpha,
        base_floor_bits,
        params,
        gc_crc,
    )
}

#[must_use]
pub(crate) fn calibrate_null_crc_from_seqs_with_gc_floor(
    seq_x: &[u8],
    seq_y: &[u8],
    n_cal: usize,
    alpha: f64,
    base_floor_bits: f64,
    params: &HmmParams,
    gc_crc: f64,
) -> f64 {
    const WIN: usize = 64;
    if n_cal == 0 || seq_x.len() < WIN * 2 || seq_y.len() < WIN * 2 {
        return gc_crc;
    }

    let max_x = seq_x.len() - WIN;
    let max_y = seq_y.len() - WIN;
    let mut scores = Vec::with_capacity(n_cal);
    let mut y_rev = [b'A'; WIN];
    let mut scratch = DpScratch::default();
    let mut tries = 0usize;
    let max_tries = n_cal * 4;
    let mut draw = 0usize;
    while scores.len() < n_cal && tries < max_tries {
        tries += 1;
        draw += 1;
        let sx = (draw.wrapping_mul(max_x / n_cal.max(1) + 97)) % (max_x + 1);
        let sy = (max_y.saturating_sub(sx % (max_y + 1)) + draw.wrapping_mul(131)) % (max_y + 1);
        let wx = &seq_x[sx..sx + WIN];
        let wy = &seq_y[sy..sy + WIN];
        let acgt_x = wx
            .iter()
            .filter(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .count();
        let acgt_y = wy
            .iter()
            .filter(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .count();
        if acgt_x * 2 < WIN || acgt_y * 2 < WIN {
            continue;
        }
        for (dst, &src) in y_rev.iter_mut().zip(wy.iter().rev()) {
            *dst = src;
        }
        let res = banded_mdl_gap_with_scratch(
            wx,
            &y_rev,
            16,
            base_floor_bits,
            params,
            MIN_POST_PROB,
            &mut scratch,
        );
        let null_bits = if res.kept_pairs.is_empty() {
            0.0
        } else {
            res.score_bits.max(0.0)
        };
        scores.push(null_bits);
    }
    if scores.len() < n_cal {
        return gc_crc;
    }
    calibrate_crc_threshold(&scores, alpha, base_floor_bits).max(gc_crc)
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
}
