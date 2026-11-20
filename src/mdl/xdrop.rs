//! Banded MDL + BOCPD X-drop flank extension.

use super::dp::{HmmParams, NormalizedHmmParams, prepare_grid_f64};
use crate::bocpd::{BocpdScratch, NEG_INF, base_probs, bocpd_gc_profile_into, joint_pair_bits};

/// Reusable per-worker scratch buffers for [`mdl_xdrop_extend_with_scratch`].
#[derive(Debug, Default)]
pub(crate) struct XDropScratch {
    bocpd: BocpdScratch,
    x_rev: Vec<u8>,
    y_rev: Vec<u8>,
    gc_x: Vec<f64>,
    gc_y: Vec<f64>,
    pi_x: Vec<[f64; 4]>,
    pi_y: Vec<[f64; 4]>,
    dp: Vec<f64>,
    bt: Vec<u8>,
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
#[allow(clippy::too_many_arguments)]
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
    let mut scratch = XDropScratch::default();
    mdl_xdrop_extend_with_scratch(
        seq_x,
        seq_y,
        off_x,
        off_y,
        reverse,
        max_ext,
        band_w,
        params,
        &mut scratch,
    )
}

#[allow(
    clippy::too_many_lines,
    clippy::many_single_char_names,
    clippy::too_many_arguments
)]
pub(crate) fn mdl_xdrop_extend_with_scratch(
    seq_x: &[u8],
    seq_y: &[u8],
    off_x: usize,
    off_y: usize,
    reverse: bool,
    max_ext: usize,
    band_w: usize,
    params: &HmmParams,
    scratch: &mut XDropScratch,
) -> Vec<(usize, usize)> {
    let x_slice: &[u8] = if reverse {
        scratch.x_rev.clear();
        scratch
            .x_rev
            .extend(seq_x.iter().rev().take(max_ext).copied());
        &scratch.x_rev
    } else {
        &seq_x[..seq_x.len().min(max_ext)]
    };
    let y_slice: &[u8] = if reverse {
        scratch.y_rev.clear();
        scratch
            .y_rev
            .extend(seq_y.iter().rev().take(max_ext).copied());
        &scratch.y_rev
    } else {
        &seq_y[..seq_y.len().min(max_ext)]
    };
    let n = x_slice.len();
    let m = y_slice.len();
    if n < 12 || m < 12 {
        return Vec::new();
    }

    bocpd_gc_profile_into(
        x_slice,
        1.0 / 120.0,
        60,
        &mut scratch.bocpd,
        &mut scratch.gc_x,
    );
    bocpd_gc_profile_into(
        y_slice,
        1.0 / 120.0,
        60,
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
    let pi_x = &scratch.pi_x;
    let pi_y = &scratch.pi_y;

    let NormalizedHmmParams {
        id_prior,
        go_bits: open_bits,
        ge_bits: extend_bits,
        mm_bits,
        im_bits,
        ..
    } = params.normalized();

    let x_drop_bits = 15.0;
    let min_peak_bits = 10.0;
    let band_w = band_w.min(n.max(m));
    let Some(stride) = band_w.checked_mul(2).and_then(|s| s.checked_add(1)) else {
        return Vec::new();
    };
    let Some(total) = n.checked_add(1).and_then(|r| r.checked_mul(stride)) else {
        return Vec::new();
    };
    if !prepare_grid_f64(&mut scratch.dp, total, NEG_INF) {
        return Vec::new();
    }
    scratch.bt.clear();
    if scratch.bt.capacity() < total && scratch.bt.try_reserve_exact(total).is_err() {
        return Vec::new();
    }
    scratch.bt.resize(total, 0_u8); // 1=diag, 2=up, 3=left
    let dp = &mut scratch.dp;
    let bt = &mut scratch.bt;

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
    use super::super::dp::HmmParams;
    use super::mdl_xdrop_extend_with_params;

    #[test]
    fn reverse_xdrop_behaves_cleanly() {
        let params = HmmParams::default();
        let x = b"TTTTTTTTACGTACGTACGTACGT";
        let y = b"GGGGGGGGACGTACGTACGTACGT";
        let pairs =
            mdl_xdrop_extend_with_params(&x[..16], &y[..16], 100, 200, true, 40, 20, &params);
        assert_ne!(pairs, [] as [(usize, usize); 0]);
        assert!(pairs.windows(2).all(|w| w[0].0 < w[1].0 && w[0].1 < w[1].1));
    }
}
