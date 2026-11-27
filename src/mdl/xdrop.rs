//! Banded MDL + BOCPD X-drop flank extension.

use super::dp::{HmmParams, NormalizedHmmParams, prepare_grid};
use crate::bocpd::{BocpdScratch, NEG_INF, base_probs, bocpd_gc_profile_into, joint_pair_bits};

/// Minimum net MDL bit-savings required after a distal mismatch or gap to retain an X-drop
/// flank tail; trims lucky 1–3 bp random overhangs past SV or inversion breakpoints.
const XDROP_TIP_CONFIRM_BITS: f64 = 8.0;

struct XDropStep {
    pair: Option<(usize, usize)>,
    is_exact: bool,
    cum_bits: f64,
}

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
    min_peak_bits: f64,
    params: &HmmParams,
    scratch: &mut XDropScratch,
) -> Vec<(usize, usize)> {
    if min_peak_bits.is_nan() || (min_peak_bits.is_infinite() && min_peak_bits.is_sign_positive()) {
        return Vec::new();
    }
    let min_peak_bits = min_peak_bits.max(0.0);
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
    let Some(layout) = super::dp::BandLayout::new_zero_centered(n, m, band_w) else {
        return Vec::new();
    };
    if !prepare_grid(&mut scratch.dp, layout.total, NEG_INF)
        || !prepare_grid(&mut scratch.bt, layout.total, 0_u8)
    {
        return Vec::new();
    }
    let dp = &mut scratch.dp;
    let bt = &mut scratch.bt; // 1=diag, 2=up, 3=left
    let band_w = layout.band_w;
    let slot = |i: usize, j: usize| layout.slot(i, j);

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

    // Backtrace from `best_cell` to `(0, 0)` recording aligned pair, exact-match flag, and DP score.
    let mut steps: Vec<XDropStep> = Vec::new();
    let (mut ci, mut cj) = best_cell;
    while ci > 0 || cj > 0 {
        let Some(s) = slot(ci, cj) else { break };
        let cum_bits = dp[s];
        match bt[s] {
            1 => {
                let bx = x_slice[ci - 1];
                let by = y_slice[cj - 1];
                let is_exact = bx.eq_ignore_ascii_case(&by)
                    && matches!(bx.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T');
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
                steps.push(XDropStep {
                    pair: Some((px, py)),
                    is_exact,
                    cum_bits,
                });
                ci -= 1;
                cj -= 1;
            }
            2 => {
                steps.push(XDropStep {
                    pair: None,
                    is_exact: false,
                    cum_bits,
                });
                ci -= 1;
            }
            3 => {
                steps.push(XDropStep {
                    pair: None,
                    is_exact: false,
                    cum_bits,
                });
                cj -= 1;
            }
            _ => break,
        }
    }
    // Put steps into extension order (anchor -> distal tip) and trim weakly confirmed
    // mismatch/gap overhangs past SV or inversion breakpoints.
    steps.reverse();
    let mut end = steps.len();
    while end > 0 {
        let tip_score = steps[end - 1].cum_bits;
        let Some(last_non_match) = (0..end).rev().find(|&k| !steps[k].is_exact) else {
            break;
        };
        let pre_score = if last_non_match > 0 {
            steps[last_non_match - 1].cum_bits
        } else {
            0.0
        };
        if tip_score - pre_score < XDROP_TIP_CONFIRM_BITS {
            end = last_non_match;
        } else {
            break;
        }
    }
    if end == 0 || steps[end - 1].cum_bits < min_peak_bits {
        return Vec::new();
    }

    let mut pairs: Vec<(usize, usize)> = steps[..end].iter().filter_map(|st| st.pair).collect();
    if reverse {
        pairs.reverse();
    }
    pairs
}

#[cfg(test)]
mod tests {
    use super::super::dp::HmmParams;
    use super::{XDropScratch, mdl_xdrop_extend_with_scratch};

    #[test]
    fn reverse_xdrop_behaves_cleanly() {
        let params = HmmParams::default();
        let x = b"TTTTTTTTACGTACGTACGTACGT";
        let y = b"GGGGGGGGACGTACGTACGTACGT";
        let mut scratch = XDropScratch::default();
        let pairs = mdl_xdrop_extend_with_scratch(
            &x[..16],
            &y[..16],
            100,
            200,
            true,
            40,
            20,
            10.0,
            &params,
            &mut scratch,
        );
        assert_ne!(pairs, [] as [(usize, usize); 0]);
        assert!(pairs.windows(2).all(|w| w[0].0 < w[1].0 && w[0].1 < w[1].1));
        for bad_thresh in [f64::INFINITY, f64::NAN] {
            let abstained = mdl_xdrop_extend_with_scratch(
                &x[..16],
                &y[..16],
                100,
                200,
                true,
                40,
                20,
                bad_thresh,
                &params,
                &mut scratch,
            );
            assert_eq!(abstained, [] as [(usize, usize); 0]);
        }
    }

    #[test]
    fn xdrop_trims_breakpoint_overhang_but_keeps_confirmed_snp() {
        let params = HmmParams::default();
        let mut scratch = XDropScratch::default();

        // 48 bp exact conserved flank, then an SV breakpoint at idx 48:
        // 1 mismatch (`A` vs `C`), 2 lucky random matches (`GT` vs `GT`), then 24 bp unrelated DNA.
        // Without tip-overhang trimming, `+1.56 + 1.56 - 1.51 > 0` pulls `best_cell` past idx 48
        // and emits a false SNP at idx 48.
        let flank = b"ACGTACGTACGTTGCATGCATGCAACGTACGTACGTTGCATGCATGCA";
        let mut x_sv = flank.to_vec();
        x_sv.extend_from_slice(b"AGTCCCAAATTTGGGCCCAAATTTGGG");
        let mut y_sv = flank.to_vec();
        y_sv.extend_from_slice(b"CGTTTTGGGAAACCCGGGTTTAAACCC");

        let pairs_sv = mdl_xdrop_extend_with_scratch(
            &x_sv,
            &y_sv,
            1000,
            2000,
            false,
            120,
            20,
            10.0,
            &params,
            &mut scratch,
        );
        assert!(
            !pairs_sv.is_empty(),
            "conserved 48 bp flank must still extend"
        );
        let max_x = pairs_sv.iter().map(|&(px, _)| px).max().unwrap_or(0);
        assert!(
            max_x < 1000 + flank.len(),
            "xdrop must not overhang past the SV breakpoint at {}; got max_x={max_x}",
            1000 + flank.len()
        );

        // Conversely, a genuine SNP inside a conserved flank followed by 24 identical bases
        // MUST be retained past the SNP.
        let mut x_snp = flank.to_vec();
        x_snp.extend_from_slice(b"AACGTACGTACGTTGCATGCATGC");
        let mut y_snp = flank.to_vec();
        y_snp.extend_from_slice(b"CACGTACGTACGTTGCATGCATGC");

        let pairs_snp = mdl_xdrop_extend_with_scratch(
            &x_snp,
            &y_snp,
            1000,
            2000,
            false,
            120,
            20,
            10.0,
            &params,
            &mut scratch,
        );
        let max_x_snp = pairs_snp.iter().map(|&(px, _)| px).max().unwrap_or(0);
        assert!(
            max_x_snp >= 1000 + flank.len() + 20,
            "genuine SNP followed by 24 identical bases must be kept; got max_x_snp={max_x_snp}"
        );
    }
}
