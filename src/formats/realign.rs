//! Bounded HOXD70 affine-gap micro-realignment and colinear run-group materialization.

use crate::seeding::complement_base;

/// Materialized alignment strings and coordinates for one colinear run group.
pub(crate) struct MaterializedRunGroup {
    pub x0: usize,
    pub x_end: usize,
    pub y_start: usize,
    pub y_span: usize,
    pub y_strand: char,
    pub score_bp: usize,
    pub ax: Vec<u8>,
    pub ay: Vec<u8>,
}

/// Bounded Gotoh affine-gap global micro-realignment for short inter-run gaps inside
/// a colinear MAF block (`gap_x.len(), gap_y.len() <= MAX_MICRO_REALIGN_BP`).
/// Uses HOXD70 substitution and affine gap costs (Chiaromonte et al., 2002) so point
/// substitutions and diverged codons stay ungapped while genuine micro-indels collapse
/// into minimal contiguous affine gaps without splitting into paired frameshifts.
const MAX_MICRO_REALIGN_BP: usize = crate::netting::MAX_COLINEAR_GAP_BP;
const MICRO_NEG_INF: i32 = -100_000_000;
const MICRO_GAP_OPEN: i32 = -400;
const MICRO_GAP_EXT: i32 = -30;

/// Standard HOXD70 DNA substitution matrix (Chiaromonte et al., 2002); ambiguous `N`/non-ACGT
/// bases score neutrally (`0`) so masked tracts stay ungapped instead of splitting into skips.
const fn hoxd70_score(a: u8, b: u8) -> i32 {
    match (a.to_ascii_uppercase(), b.to_ascii_uppercase()) {
        (b'A', b'A') | (b'T', b'T') => 91,
        (b'C', b'C') | (b'G', b'G') => 100,
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C') => -31,
        (b'A', b'C') | (b'C', b'A') | (b'G', b'T') | (b'T', b'G') => -114,
        (b'A', b'T') | (b'T', b'A') => -123,
        (b'C', b'G') | (b'G', b'C') => -125,
        _ => 0,
    }
}

fn best_of_three(v0: i32, v1: i32, v2: i32) -> (i32, u8) {
    if v0 >= v1 && v0 >= v2 {
        (v0, 0)
    } else if v1 >= v2 {
        (v1, 1)
    } else {
        (v2, 2)
    }
}

#[allow(clippy::too_many_lines)]
fn align_micro_gap_affine(ref_gap: &[u8], qry_gap: &[u8], ax: &mut Vec<u8>, ay: &mut Vec<u8>) {
    const STRIDE: usize = MAX_MICRO_REALIGN_BP + 1;
    let m = ref_gap.len();
    let n = qry_gap.len();
    if m == 0 {
        for &by in qry_gap {
            ax.push(b'-');
            ay.push(by);
        }
        return;
    }
    if n == 0 {
        for &bx in ref_gap {
            ax.push(bx);
            ay.push(b'-');
        }
        return;
    }
    if m == 1 && n == 1 {
        ax.push(ref_gap[0]);
        ay.push(qry_gap[0]);
        return;
    }
    if m > MAX_MICRO_REALIGN_BP || n > MAX_MICRO_REALIGN_BP {
        for &bx in ref_gap {
            ax.push(bx);
            ay.push(b'-');
        }
        for &by in qry_gap {
            ax.push(b'-');
            ay.push(by);
        }
        return;
    }

    let cols = n + 1;
    let mut prev_m = [MICRO_NEG_INF; STRIDE];
    let mut prev_x = [MICRO_NEG_INF; STRIDE];
    let mut prev_y = [MICRO_NEG_INF; STRIDE];
    let mut cur_m = [MICRO_NEG_INF; STRIDE];
    let mut cur_x = [MICRO_NEG_INF; STRIDE];
    let mut cur_y = [MICRO_NEG_INF; STRIDE];
    // Bit-packed traceback per cell: bits 0..2 = tb_m, bits 2..4 = tb_x, bits 4..6 = tb_y.
    let mut tb = [0_u8; STRIDE * STRIDE];

    prev_m[0] = 0;
    for j in 1..=n {
        prev_y[j] = MICRO_GAP_OPEN + MICRO_GAP_EXT * i32::try_from(j).unwrap_or(0);
        tb[j] = if j > 1 { 2 << 4 } else { 0 };
    }

    for i in 1..=m {
        let bx = ref_gap[i - 1];
        let row_off = i * cols;
        cur_m[0] = MICRO_NEG_INF;
        cur_x[0] = MICRO_GAP_OPEN + MICRO_GAP_EXT * i32::try_from(i).unwrap_or(0);
        cur_y[0] = MICRO_NEG_INF;
        tb[row_off] = u8::from(i > 1) << 2;

        for j in 1..=n {
            let sub = hoxd70_score(bx, qry_gap[j - 1]);
            let (bm, sm) = best_of_three(prev_m[j - 1], prev_x[j - 1], prev_y[j - 1]);
            cur_m[j] = bm + sub;

            let (ref_gap_v, ref_gap_s) = best_of_three(
                prev_m[j] + MICRO_GAP_OPEN + MICRO_GAP_EXT,
                prev_x[j] + MICRO_GAP_EXT,
                prev_y[j] + MICRO_GAP_OPEN + MICRO_GAP_EXT,
            );
            cur_x[j] = ref_gap_v;

            let (qry_gap_v, qry_gap_s) = best_of_three(
                cur_m[j - 1] + MICRO_GAP_OPEN + MICRO_GAP_EXT,
                cur_x[j - 1] + MICRO_GAP_OPEN + MICRO_GAP_EXT,
                cur_y[j - 1] + MICRO_GAP_EXT,
            );
            cur_y[j] = qry_gap_v;
            tb[row_off + j] = sm | (ref_gap_s << 2) | (qry_gap_s << 4);
        }
        std::mem::swap(&mut prev_m, &mut cur_m);
        std::mem::swap(&mut prev_x, &mut cur_x);
        std::mem::swap(&mut prev_y, &mut cur_y);
    }

    let (_, mut state) = best_of_three(prev_m[n], prev_x[n], prev_y[n]);
    let mut rev_x = [0u8; MAX_MICRO_REALIGN_BP * 2];
    let mut rev_y = [0u8; MAX_MICRO_REALIGN_BP * 2];
    let mut rev_len = 0usize;
    let (mut i, mut j) = (m, n);
    while i > 0 || j > 0 {
        let packed = tb[i * cols + j];
        match state {
            0 => {
                state = packed & 0b11;
                i -= 1;
                j -= 1;
                rev_x[rev_len] = ref_gap[i];
                rev_y[rev_len] = qry_gap[j];
            }
            1 => {
                state = (packed >> 2) & 0b11;
                i -= 1;
                rev_x[rev_len] = ref_gap[i];
                rev_y[rev_len] = b'-';
            }
            _ => {
                state = (packed >> 4) & 0b11;
                j -= 1;
                rev_x[rev_len] = b'-';
                rev_y[rev_len] = qry_gap[j];
            }
        }
        rev_len += 1;
    }
    ax.reserve(rev_len);
    ay.reserve(rev_len);
    for k in (0..rev_len).rev() {
        ax.push(rev_x[k]);
        ay.push(rev_y[k]);
    }
}

fn emit_inter_run_gap(
    seq_x: &[u8],
    seq_y: &[u8],
    x_range: std::ops::Range<usize>,
    y_range: std::ops::Range<usize>,
    is_rev: bool,
    ax: &mut Vec<u8>,
    ay: &mut Vec<u8>,
) {
    let ref_len = x_range.len();
    let qry_len = y_range.len();
    if ref_len == 1 && qry_len == 1 {
        ax.push(seq_x[x_range.start].to_ascii_uppercase());
        let raw_y = seq_y[y_range.start].to_ascii_uppercase();
        ay.push(if is_rev {
            complement_base(raw_y)
        } else {
            raw_y
        });
        return;
    }
    if ref_len > 0
        && qry_len > 0
        && ref_len <= MAX_MICRO_REALIGN_BP
        && qry_len <= MAX_MICRO_REALIGN_BP
    {
        let mut ref_up = [0u8; MAX_MICRO_REALIGN_BP];
        for (dst, &bx) in ref_up[..ref_len].iter_mut().zip(&seq_x[x_range]) {
            *dst = bx.to_ascii_uppercase();
        }
        let mut qry_up = [0u8; MAX_MICRO_REALIGN_BP];
        if is_rev {
            for (dst, &by) in qry_up[..qry_len]
                .iter_mut()
                .zip(seq_y[y_range].iter().rev())
            {
                *dst = complement_base(by.to_ascii_uppercase());
            }
        } else {
            for (dst, &by) in qry_up[..qry_len].iter_mut().zip(&seq_y[y_range]) {
                *dst = by.to_ascii_uppercase();
            }
        }
        align_micro_gap_affine(&ref_up[..ref_len], &qry_up[..qry_len], ax, ay);
    } else {
        for &bx in &seq_x[x_range] {
            ax.push(bx.to_ascii_uppercase());
            ay.push(b'-');
        }
        if is_rev {
            for &by in seq_y[y_range].iter().rev() {
                ax.push(b'-');
                ay.push(complement_base(by.to_ascii_uppercase()));
            }
        } else {
            for &by in &seq_y[y_range] {
                ax.push(b'-');
                ay.push(by.to_ascii_uppercase());
            }
        }
    }
}

/// Validate a colinear run group and materialize its uppercase gapped alignment rows.
#[must_use]
pub(crate) fn materialize_run_group(
    len_x: usize,
    len_y: usize,
    seq_x: &[u8],
    seq_y: &[u8],
    grp: &[crate::mdl::AlignedRun],
    is_rev: bool,
) -> Option<MaterializedRunGroup> {
    if grp.is_empty() || len_x < seq_x.len() || len_y < seq_y.len() {
        return None;
    }
    let valid = grp.iter().all(|r| {
        r.len > 0
            && r.is_rc == is_rev
            && r.x_start
                .checked_add(r.len)
                .is_some_and(|e| e <= seq_x.len())
            && if r.is_rc {
                r.y_start < seq_y.len() && r.y_start + 1 >= r.len
            } else {
                r.y_start
                    .checked_add(r.len)
                    .is_some_and(|e| e <= seq_y.len())
            }
    }) && grp.windows(2).all(|w| {
        let prev = w[0];
        let cur = w[1];
        let x_ok = prev.x_start + prev.len <= cur.x_start;
        let y_ok = if is_rev {
            let prev_y_lo = prev.y_start + 1 - prev.len;
            cur.y_start < prev_y_lo
        } else {
            prev.y_start + prev.len <= cur.y_start
        };
        x_ok && y_ok
    });
    if !valid {
        return None;
    }
    let x0 = grp[0].x_start;
    let last = grp[grp.len() - 1];
    let x_end = last.x_start + last.len;
    let score_bp: usize = grp.iter().map(|r| r.len).sum();
    let mut ax = Vec::with_capacity(x_end - x0 + 16);
    let mut ay = Vec::with_capacity(x_end - x0 + 16);

    let mut prev_x_end = x0;
    let mut prev_y_next = grp[0].y_start;
    for (ri, r) in grp.iter().enumerate() {
        if ri > 0 {
            let y_range = if is_rev {
                (r.y_start + 1)..prev_y_next
            } else {
                prev_y_next..r.y_start
            };
            emit_inter_run_gap(
                seq_x,
                seq_y,
                prev_x_end..r.x_start,
                y_range,
                is_rev,
                &mut ax,
                &mut ay,
            );
        }
        for &bx in &seq_x[r.x_start..(r.x_start + r.len)] {
            ax.push(bx.to_ascii_uppercase());
        }
        if is_rev {
            let y_lo = r.y_start + 1 - r.len;
            for &by in seq_y[y_lo..=r.y_start].iter().rev() {
                ay.push(complement_base(by.to_ascii_uppercase()));
            }
            prev_y_next = y_lo;
        } else {
            let y_hi = r.y_start + r.len;
            for &by in &seq_y[r.y_start..y_hi] {
                ay.push(by.to_ascii_uppercase());
            }
            prev_y_next = y_hi;
        }
        prev_x_end = r.x_start + r.len;
    }

    let (y_start, y_span, y_strand) = if is_rev {
        let y_lo = last.y_start + 1 - last.len;
        let y_hi = grp[0].y_start + 1;
        let iv = crate::coord::Interval::new(y_lo, y_hi, '-', len_y);
        (iv.ucsc_start(), y_hi - y_lo, '-')
    } else {
        let y_lo = grp[0].y_start;
        let y_hi = last.y_start + last.len;
        (y_lo, y_hi - y_lo, '+')
    };

    Some(MaterializedRunGroup {
        x0,
        x_end,
        y_start,
        y_span,
        y_strand,
        score_bp,
        ax,
        ay,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mdl::AlignedRun;

    #[test]
    fn materialize_run_group_realigns_micro_gaps_without_balanced_skip_artifacts() {
        // 1-bp point substitution (gx=1, gy=1): must emit 1 mismatch column, ZERO gaps!
        let sx = b"ACGTACGT";
        let sy = b"ACGTTCGT";
        let grp = [
            AlignedRun {
                x_start: 0,
                y_start: 0,
                len: 4,
                is_rc: false,
            },
            AlignedRun {
                x_start: 5,
                y_start: 5,
                len: 3,
                is_rc: false,
            },
        ];
        let mat = materialize_run_group(sx.len(), sy.len(), sx, sy, &grp, false)
            .expect("materialize 1-bp mismatch");
        assert_eq!(mat.ax, b"ACGTACGT");
        assert_eq!(mat.ay, b"ACGTTCGT");

        // 1-bp mismatch + 3-bp in-frame insertion (gx=1, gy=4): must collapse to a single 3-bp gap!
        let sx2 = b"ACGTACGT";
        let sy2 = b"ACGTGGGTCGT";
        let grp2 = [
            AlignedRun {
                x_start: 0,
                y_start: 0,
                len: 4,
                is_rc: false,
            },
            AlignedRun {
                x_start: 5,
                y_start: 8,
                len: 3,
                is_rc: false,
            },
        ];
        let mat2 = materialize_run_group(sx2.len(), sy2.len(), sx2, sy2, &grp2, false)
            .expect("materialize 1-vs-4 micro-gap");
        assert_eq!(mat2.ax.len(), mat2.ay.len());
        assert_eq!(
            mat2.ax
                .iter()
                .filter(|&&b| b != b'-')
                .copied()
                .collect::<Vec<_>>(),
            sx2
        );
        assert_eq!(
            mat2.ay
                .iter()
                .filter(|&&b| b != b'-')
                .copied()
                .collect::<Vec<_>>(),
            sy2
        );
        assert!(!mat2.ay.contains(&b'-'));
        assert_eq!(
            mat2.ax
                .iter()
                .fold(0usize, |acc, &b| acc + usize::from(b == b'-')),
            3
        );
    }

    #[test]
    fn micro_gap_realign_does_not_shatter_diverged_codon_window_into_paired_frameshifts() {
        let mut ax = Vec::new();
        let mut ay = Vec::new();
        align_micro_gap_affine(b"GCTA", b"CTAG", &mut ax, &mut ay);
        assert_eq!(ax, b"GCTA");
        assert_eq!(ay, b"CTAG");
    }

    #[test]
    fn micro_gap_realign_prefers_single_inframe_deletion_over_split_frameshifts() {
        let mut ax = Vec::new();
        let mut ay = Vec::new();
        align_micro_gap_affine(b"GGTTCGATCC", b"GGCGTCC", &mut ax, &mut ay);
        assert_eq!(ax.len(), ay.len());
        assert!(!ax.contains(&b'-'), "deletion must occur only in query row");
        let mut gap_runs = Vec::new();
        let mut cur = 0usize;
        for &b in &ay {
            if b == b'-' {
                cur += 1;
            } else if cur > 0 {
                gap_runs.push(cur);
                cur = 0;
            }
        }
        if cur > 0 {
            gap_runs.push(cur);
        }
        assert_eq!(
            gap_runs,
            vec![3],
            "expected one contiguous 3-bp in-frame deletion, got runs {gap_runs:?} in ay={}",
            String::from_utf8_lossy(&ay)
        );
    }

    #[test]
    fn micro_gap_realign_handles_n_runs_30bp_40bp_overflow_and_rc_gaps() {
        // R5 regression: 24 bp N-run must stay ungapped.
        let mut ax = Vec::new();
        let mut ay = Vec::new();
        let n24 = [b'N'; 24];
        align_micro_gap_affine(&n24, &n24, &mut ax, &mut ay);
        assert_eq!(ax, n24);
        assert_eq!(ay, n24);

        // Exact MAX_MICRO_REALIGN_BP (40x40) stack boundary.
        let seq40 = b"ACGTACGTACGTTGCATGCATGCATGCAACGTACGTACGT";
        assert_eq!(seq40.len(), MAX_MICRO_REALIGN_BP);
        ax.clear();
        ay.clear();
        align_micro_gap_affine(seq40, seq40, &mut ax, &mut ay);
        assert_eq!(ax, seq40);
        assert_eq!(ay, seq40);

        // Overflow (> MAX_MICRO_REALIGN_BP) falls back to clean sequential skip without panic.
        let seq41 = [b'A'; MAX_MICRO_REALIGN_BP + 1];
        ax.clear();
        ay.clear();
        align_micro_gap_affine(&seq41, &seq41, &mut ax, &mut ay);
        assert_eq!(ax.len(), (MAX_MICRO_REALIGN_BP + 1) * 2);
        assert_eq!(ay.len(), (MAX_MICRO_REALIGN_BP + 1) * 2);

        // R4 regression: 30 bp inter-run gap realigns cleanly.
        let mut sx = b"ACGTACGT".to_vec();
        sx.extend_from_slice(b"GCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");
        sx.extend_from_slice(b"TGCATGCA");
        let mut sy = b"ACGTACGT".to_vec();
        sy.extend_from_slice(b"GCTAGCTAGCTAGCTAACGCTAGCTAGC");
        sy.extend_from_slice(b"TGCATGCA");
        let grp = [
            AlignedRun {
                x_start: 0,
                y_start: 0,
                len: 8,
                is_rc: false,
            },
            AlignedRun {
                x_start: 38,
                y_start: 36,
                len: 8,
                is_rc: false,
            },
        ];
        let mat = materialize_run_group(sx.len(), sy.len(), &sx, &sy, &grp, false)
            .expect("materialize 30 bp gap");
        assert!(mat.ax.len() <= 48);

        // Reverse-strand inter-run micro-gap (gx=4, gy=4 with 1 SNP on RC strand).
        let rx = b"ACGTCCGTTGCA";
        // RC of "ACGTCGGTTGCA" = "TGCAACCGACGT"
        let ry = b"TGCAACCGACGT";
        let grp_rc = [
            AlignedRun {
                x_start: 0,
                y_start: 11,
                len: 4,
                is_rc: true,
            },
            AlignedRun {
                x_start: 8,
                y_start: 3,
                len: 4,
                is_rc: true,
            },
        ];
        let mat_rc = materialize_run_group(rx.len(), ry.len(), rx, ry, &grp_rc, true)
            .expect("materialize rc micro-gap");
        assert_eq!(mat_rc.y_strand, '-');
        assert_eq!(mat_rc.ax, b"ACGTCCGTTGCA");
        assert_eq!(mat_rc.ay, b"ACGTCGGTTGCA");
    }
}
