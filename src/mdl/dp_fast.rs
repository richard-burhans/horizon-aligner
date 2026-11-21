//! Fast-path gap aligners (diagonal, narrow banded Viterbi, exact flank trimming).

use super::dp::{DpScratch, GapResult, band_center, prepare_grid};
use crate::bocpd::base_to_idx_fast;

pub(super) fn try_diagonal_fast_path(
    seq_x: &[u8],
    seq_y: &[u8],
    min_block_bits: f64,
    min_post_prob: f64,
) -> Option<GapResult> {
    const WIN: usize = 24;
    let n = seq_x.len();
    if n != seq_y.len() || n < 16 || min_post_prob > 0.90 || !min_block_bits.is_finite() {
        return None;
    }
    let mut matches = 0usize;
    let mut mismatches = 0usize;
    let mut win_mm = 0usize;
    let mut is_mm = [false; WIN];

    for (i, (&bx, &by)) in seq_x.iter().zip(seq_y).enumerate() {
        let ix = base_to_idx_fast(bx);
        let iy = base_to_idx_fast(by);
        if ix >= 4 || iy >= 4 {
            return None;
        }
        let mm = ix != iy;
        if i >= WIN && is_mm[i % WIN] {
            win_mm -= 1;
        }
        is_mm[i % WIN] = mm;
        if mm {
            mismatches += 1;
            win_mm += 1;
            if win_mm > 2 {
                return None;
            }
        } else {
            matches += 1;
        }
    }
    if matches * 10 < n * 9 {
        return None;
    }
    #[allow(clippy::cast_precision_loss)]
    let score_est = (matches as f64) * 1.6 - (mismatches as f64) * 2.5;
    if score_est < min_block_bits + 2.0 {
        return None;
    }
    Some(GapResult {
        kept_pairs: (0..n).map(|i| (i, i)).collect(),
        bubbles: 0,
        score_bits: score_est,
    })
}

#[allow(
    clippy::too_many_lines,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]
pub(super) fn try_banded_viterbi_fast_path(
    seq_x: &[u8],
    seq_y: &[u8],
    band_w: usize,
    min_block_bits: f64,
    min_post_prob: f64,
    scratch: &mut DpScratch,
) -> Option<GapResult> {
    const MATCH_SCORE: f32 = 1.8;
    const MISMATCH_SCORE: f32 = -2.2;
    const GAP_OPEN: f32 = -4.0;
    const GAP_EXT: f32 = -1.2;
    const NEG: f32 = -1.0e20;

    let n = seq_x.len();
    let m = seq_y.len();
    if n < 32
        || m < 32
        || n > 2048
        || m > 2048
        || n.abs_diff(m) > 12
        || min_post_prob > 0.85
        || !min_block_bits.is_finite()
    {
        return None;
    }
    let bw = n.abs_diff(m) + 4;
    if band_w < bw {
        return None;
    }
    let stride = 2 * bw + 1;
    let rows = n + 1;
    let total = rows * stride;

    if !prepare_grid(&mut scratch.vit_m, total, NEG)
        || !prepare_grid(&mut scratch.vit_x, total, NEG)
        || !prepare_grid(&mut scratch.vit_y, total, NEG)
        || !prepare_grid(&mut scratch.vit_trace, total, 0)
    {
        return None;
    }
    let dp_m = &mut scratch.vit_m;
    let dp_x = &mut scratch.vit_x;
    let dp_y = &mut scratch.vit_y;
    let trace = &mut scratch.vit_trace;

    let slot = |i: usize, j: usize| -> Option<usize> {
        let center = band_center(i, m, n) as isize;
        let diff = (j as isize) - center;
        let bwi = bw as isize;
        if diff.abs() <= bwi {
            Some(i * stride + (diff + bwi) as usize)
        } else {
            None
        }
    };

    let s0 = slot(0, 0)?;
    dp_m[s0] = 0.0;
    for j in 1..=bw.min(m) {
        if let Some(s) = slot(0, j) {
            dp_y[s] = GAP_OPEN + (j as f32 - 1.0) * GAP_EXT;
            trace[s] = 2;
        }
    }
    for i in 1..=bw.min(n) {
        if let Some(s) = slot(i, 0) {
            dp_x[s] = GAP_OPEN + (i as f32 - 1.0) * GAP_EXT;
            trace[s] = 1;
        }
    }

    for i in 1..=n {
        let bx = base_to_idx_fast(seq_x[i - 1]);
        let center = band_center(i, m, n);
        let j_lo = center.saturating_sub(bw).max(1);
        let j_hi = (center + bw).min(m);
        for j in j_lo..=j_hi {
            let Some(cur) = slot(i, j) else {
                continue;
            };
            let by = base_to_idx_fast(seq_y[j - 1]);
            let emit = if bx < 4 && by < 4 {
                if bx == by {
                    MATCH_SCORE
                } else {
                    MISMATCH_SCORE
                }
            } else {
                0.0
            };

            let mut m_pred = 0u8;
            if let Some(prev) = slot(i - 1, j - 1) {
                let pm = dp_m[prev];
                let px = dp_x[prev];
                let py = dp_y[prev];
                let (best_prev, pred) = if pm >= px && pm >= py {
                    (pm, 0u8)
                } else if px >= py {
                    (px, 1u8)
                } else {
                    (py, 2u8)
                };
                dp_m[cur] = best_prev + emit;
                m_pred = pred;
            }
            let mut x_pred = 0u8;
            if let Some(up) = slot(i - 1, j) {
                let from_m = dp_m[up] + GAP_OPEN;
                let from_x = dp_x[up] + GAP_EXT;
                if from_m >= from_x {
                    dp_x[cur] = from_m;
                    x_pred = 0;
                } else {
                    dp_x[cur] = from_x;
                    x_pred = 1;
                }
            }
            let mut y_pred = 0u8;
            if let Some(left) = slot(i, j - 1) {
                let from_m = dp_m[left] + GAP_OPEN;
                let from_y = dp_y[left] + GAP_EXT;
                if from_m >= from_y {
                    dp_y[cur] = from_m;
                    y_pred = 0;
                } else {
                    dp_y[cur] = from_y;
                    y_pred = 2;
                }
            }

            let m_val = dp_m[cur];
            let x_val = dp_x[cur];
            let y_val = dp_y[cur];
            let win = if m_val >= x_val && m_val >= y_val {
                0u8
            } else if x_val >= y_val {
                1u8
            } else {
                2u8
            };
            // Pack winning state + per-state predecessors: bits 0..2=win, 2..4=m_pred, 4..6=x_pred, 6..8=y_pred.
            trace[cur] = win | (m_pred << 2) | (x_pred << 4) | (y_pred << 6);
        }
    }

    let s_end = slot(n, m)?;
    let end_score = dp_m[s_end].max(dp_x[s_end]).max(dp_y[s_end]);
    if !end_score.is_finite() || f64::from(end_score) < min_block_bits + 4.0 {
        return None;
    }

    let mut pairs = Vec::with_capacity(n.min(m));
    let mut i = n;
    let mut j = m;
    let mut st = trace[s_end] & 0b11;
    let mut matches = 0usize;
    let mut cols = 0usize;
    while i > 0 || j > 0 {
        let cur = slot(i, j)?;
        let packed = trace[cur];
        cols += 1;
        if st == 0 && i > 0 && j > 0 {
            let bx = base_to_idx_fast(seq_x[i - 1]);
            let by = base_to_idx_fast(seq_y[j - 1]);
            if bx < 4 && bx == by {
                matches += 1;
            }
            pairs.push((i - 1, j - 1));
            st = (packed >> 2) & 0b11;
            i -= 1;
            j -= 1;
        } else if st == 1 && i > 0 {
            st = (packed >> 4) & 0b11;
            i -= 1;
        } else if j > 0 {
            st = (packed >> 6) & 0b11;
            j -= 1;
        } else {
            return None;
        }
    }
    if cols == 0 || matches * 100 < cols * 88 {
        return None;
    }
    pairs.reverse();
    Some(GapResult {
        kept_pairs: pairs,
        bubbles: 0,
        score_bits: f64::from(end_score),
    })
}

pub(super) fn matching_prefix_suffix_lens(seq_x: &[u8], seq_y: &[u8]) -> (usize, usize) {
    let min_len = seq_x.len().min(seq_y.len());
    let mut pref = 0usize;
    while pref < min_len {
        let ix = base_to_idx_fast(seq_x[pref]);
        let iy = base_to_idx_fast(seq_y[pref]);
        if ix >= 4 || ix != iy {
            break;
        }
        pref += 1;
    }
    let mut suff = 0usize;
    let max_suff = min_len - pref;
    while suff < max_suff {
        let ix = base_to_idx_fast(seq_x[seq_x.len() - 1 - suff]);
        let iy = base_to_idx_fast(seq_y[seq_y.len() - 1 - suff]);
        if ix >= 4 || ix != iy {
            break;
        }
        suff += 1;
    }
    (pref, suff)
}
