//! Fast-path gap aligners (diagonal, narrow banded Viterbi, exact flank trimming).

use super::dp::{DpScratch, GapResult, HmmParams, band_center, prepare_grid};
use crate::bocpd::{
    GAP_EXTEND_BITS, GAP_OPEN_BITS, IDENTITY_PRIOR, base_probs, base_to_idx_fast, joint_pair_bits,
};

/// Composition-aware Pair-HMM bit scores derived strictly from `params.normalized()` and the
/// empirical GC content of `seq_x`/`seq_y`, ensuring fast paths never exceed exact MDL bits.
struct FastPathBits {
    at_match: f64,
    gc_match: f64,
    mismatch: f64,
    stay_match: f64,
    leave_gap: f64,
    gap_open: f64,
    gap_extend: f64,
}

fn fast_path_scores(gc_x: f64, gc_y: f64, params: &HmmParams) -> FastPathBits {
    let pi_x = base_probs(gc_x);
    let pi_y = base_probs(gc_y);
    let norm = params.normalized();
    FastPathBits {
        at_match: joint_pair_bits(b'A', b'A', pi_x, pi_y, norm.id_prior),
        gc_match: joint_pair_bits(b'C', b'C', pi_x, pi_y, norm.id_prior),
        mismatch: joint_pair_bits(b'A', b'C', pi_x, pi_y, norm.id_prior),
        stay_match: norm.mm_bits,
        leave_gap: norm.im_bits,
        gap_open: norm.go_bits,
        gap_extend: norm.ge_bits,
    }
}

/// Fast gap aligners assume near-default Pair-HMM parameters and balanced nucleotide composition
/// so low-complexity (`poly-A`, `(AT)_n`) gaps always route through exact BOCPD+MDL.
/// Returns `Some((gc_x, gc_y))` when eligible, avoiding a redundant base-count scan.
#[must_use]
pub(super) fn is_fast_path_eligible(
    seq_x: &[u8],
    seq_y: &[u8],
    params: &HmmParams,
) -> Option<(f64, f64)> {
    if !params.identity_prior.is_finite()
        || !params.gap_open_bits.is_finite()
        || !params.gap_extend_bits.is_finite()
        || (params.identity_prior - IDENTITY_PRIOR).abs() > 0.08
        || (params.gap_open_bits - GAP_OPEN_BITS).abs() > 0.5
        || (params.gap_extend_bits - GAP_EXTEND_BITS).abs() > 0.25
    {
        return None;
    }
    let mut gcs = [0.5_f64; 2];
    for (slot, seq) in [seq_x, seq_y].into_iter().enumerate() {
        let mut counts = [0usize; 4];
        let mut acgt = 0usize;
        for &b in seq {
            let idx = base_to_idx_fast(b);
            if idx < 4 {
                counts[idx as usize] += 1;
                acgt += 1;
            }
        }
        if acgt < 16 {
            return None;
        }
        let gc = counts[1] + counts[2];
        // Require GC in [32%, 68%], max single-base <= 40%, and min single-base >= 12%.
        if gc * 100 < acgt * 32 || gc * 100 > acgt * 68 {
            return None;
        }
        let max_c = counts.iter().copied().max().unwrap_or(0);
        let min_c = counts.iter().copied().min().unwrap_or(0);
        if max_c * 100 > acgt * 40 || min_c * 100 < acgt * 12 {
            return None;
        }
        #[allow(clippy::cast_precision_loss)]
        {
            gcs[slot] = gc as f64 / acgt as f64;
        }
    }
    Some((gcs[0], gcs[1]))
}

pub(super) fn try_diagonal_fast_path(
    seq_x: &[u8],
    seq_y: &[u8],
    min_block_bits: f64,
    params: &HmmParams,
    min_post_prob: f64,
) -> Option<GapResult> {
    const WIN: usize = 24;
    let n = seq_x.len();
    if n != seq_y.len() || n < 16 || min_post_prob > 0.90 || !min_block_bits.is_finite() {
        return None;
    }
    let (gc_x, gc_y) = is_fast_path_eligible(seq_x, seq_y, params)?;
    let scores = fast_path_scores(gc_x, gc_y, params);
    let mut matches = 0usize;
    let mut emit_sum = 0.0_f64;
    let mut net_score = 0.0_f64;
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
            win_mm += 1;
            if win_mm > 2 {
                return None;
            }
            emit_sum += scores.mismatch;
            net_score += scores.mismatch + scores.stay_match;
        } else {
            matches += 1;
            let mb = if ix == 1 || ix == 2 {
                scores.gc_match
            } else {
                scores.at_match
            };
            emit_sum += mb;
            net_score += mb + scores.stay_match;
        }
    }
    if matches * 10 < n * 9 {
        return None;
    }
    // Both exact net Pair-HMM score and emission sum must clear `min_block_bits + 2.0`.
    let conservative_score = net_score.min(emit_sum);
    if !conservative_score.is_finite() || conservative_score < min_block_bits + 2.0 {
        return None;
    }
    Some(GapResult {
        kept_pairs: (0..n).map(|i| (i, i)).collect(),
        bubbles: 0,
        score_bits: net_score,
    })
}

#[allow(
    clippy::too_many_lines,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation
)]
pub(super) fn try_banded_viterbi_fast_path(
    seq_x: &[u8],
    seq_y: &[u8],
    band_w: usize,
    min_block_bits: f64,
    params: &HmmParams,
    min_post_prob: f64,
    scratch: &mut DpScratch,
) -> Option<GapResult> {
    const WIN: usize = 24;
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
    let (gc_x, gc_y) = is_fast_path_eligible(seq_x, seq_y, params)?;
    let scores = fast_path_scores(gc_x, gc_y, params);
    let at_match = scores.at_match as f32;
    let gc_match = scores.gc_match as f32;
    let mismatch_score = scores.mismatch as f32;
    let stay_match = scores.stay_match as f32;
    let leave_gap = scores.leave_gap as f32;
    let gap_open = scores.gap_open as f32;
    let gap_ext = scores.gap_extend as f32;
    let layout = super::dp::BandLayout::new(n, m, bw)?;

    if !prepare_grid(&mut scratch.vit_m, layout.total, NEG)
        || !prepare_grid(&mut scratch.vit_x, layout.total, NEG)
        || !prepare_grid(&mut scratch.vit_y, layout.total, NEG)
        || !prepare_grid(&mut scratch.vit_trace, layout.total, 0)
    {
        return None;
    }
    let dp_m = &mut scratch.vit_m;
    let dp_x = &mut scratch.vit_x;
    let dp_y = &mut scratch.vit_y;
    let trace = &mut scratch.vit_trace;
    let slot = |i: usize, j: usize| layout.slot(i, j);

    let s0 = slot(0, 0)?;
    dp_m[s0] = 0.0;
    for j in 1..=bw.min(m) {
        if let Some(s) = slot(0, j) {
            dp_y[s] = gap_open + (j as f32 - 1.0) * gap_ext;
            let y_pred = if j > 1 { 2u8 } else { 0u8 };
            trace[s] = 2 | (y_pred << 6);
        }
    }
    for i in 1..=bw.min(n) {
        if let Some(s) = slot(i, 0) {
            dp_x[s] = gap_open + (i as f32 - 1.0) * gap_ext;
            let x_pred = u8::from(i > 1);
            trace[s] = 1 | (x_pred << 4);
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
                    if bx == 1 || bx == 2 {
                        gc_match
                    } else {
                        at_match
                    }
                } else {
                    mismatch_score
                }
            } else {
                0.0
            };

            let mut m_pred = 0u8;
            if let Some(prev) = slot(i - 1, j - 1) {
                let pm = dp_m[prev] + stay_match;
                let px = dp_x[prev] + leave_gap;
                let py = dp_y[prev] + leave_gap;
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
                let from_m = dp_m[up] + gap_open;
                let from_x = dp_x[up] + gap_ext;
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
                let from_m = dp_m[left] + gap_open;
                let from_y = dp_y[left] + gap_ext;
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
    let mut win_nonmatch = 0usize;
    let mut is_nonmatch = [false; WIN];
    let mut cur_block_bits = 0.0_f64;
    let mut cur_block_len = 0usize;
    while i > 0 || j > 0 {
        let cur = slot(i, j)?;
        let packed = trace[cur];
        let col_idx = cols;
        cols += 1;
        let nonmatch = if st == 0 && i > 0 && j > 0 {
            let bx = base_to_idx_fast(seq_x[i - 1]);
            let by = base_to_idx_fast(seq_y[j - 1]);
            if bx >= 4 || by >= 4 {
                // Any ambiguous base routes to exact Forward-Backward MEA.
                return None;
            }
            let is_match = bx == by;
            if is_match {
                matches += 1;
                cur_block_bits += if bx == 1 || bx == 2 {
                    scores.gc_match
                } else {
                    scores.at_match
                };
            } else {
                cur_block_bits += scores.mismatch;
            }
            cur_block_len += 1;
            pairs.push((i - 1, j - 1));
            st = (packed >> 2) & 0b11;
            i -= 1;
            j -= 1;
            !is_match
        } else if st == 1 && i > 0 {
            if cur_block_len > 0 && cur_block_bits < min_block_bits {
                return None;
            }
            cur_block_bits = 0.0;
            cur_block_len = 0;
            st = (packed >> 4) & 0b11;
            i -= 1;
            true
        } else if st == 2 && j > 0 {
            if cur_block_len > 0 && cur_block_bits < min_block_bits {
                return None;
            }
            cur_block_bits = 0.0;
            cur_block_len = 0;
            st = (packed >> 6) & 0b11;
            j -= 1;
            true
        } else {
            return None;
        };
        if col_idx >= WIN && is_nonmatch[col_idx % WIN] {
            win_nonmatch -= 1;
        }
        is_nonmatch[col_idx % WIN] = nonmatch;
        if nonmatch {
            win_nonmatch += 1;
            // Allow at most 4 non-match columns (e.g. one 3-bp codon indel + 1 SNP) per 24 bp;
            // any localized divergence bubble falls through to exact Forward-Backward MEA.
            if win_nonmatch > 4 {
                return None;
            }
        }
    }
    if cur_block_len > 0 && cur_block_bits < min_block_bits {
        return None;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn banded_viterbi_fast_path_handles_leading_deletion_and_insertion_of_len_3() {
        // R1 regression: leading 3 bp deletion in y (x starts with 3 extra bases) must trace back
        // through (i >= 2, 0) instead of aborting with None.
        let base = b"ACGTACGTACGTTGCATGCATGCATGCAACGTACGTACGT";
        let mut x_del = b"CCG".to_vec();
        x_del.extend_from_slice(base);
        let mut scratch = DpScratch::default();
        let res_del = try_banded_viterbi_fast_path(
            &x_del,
            base,
            16,
            10.0,
            &HmmParams::default(),
            0.15,
            &mut scratch,
        )
        .expect("leading deletion of len 3 must succeed in banded Viterbi fast path");
        assert_eq!(res_del.kept_pairs.len(), base.len());
        assert_eq!(res_del.kept_pairs[0], (3, 0));

        let mut y_ins = b"GGT".to_vec();
        y_ins.extend_from_slice(base);
        let res_ins = try_banded_viterbi_fast_path(
            base,
            &y_ins,
            16,
            10.0,
            &HmmParams::default(),
            0.15,
            &mut scratch,
        )
        .expect("leading insertion of len 3 must succeed in banded Viterbi fast path");
        assert_eq!(res_ins.kept_pairs.len(), base.len());
        assert_eq!(res_ins.kept_pairs[0], (0, 3));
    }
}
