//! Six-frame biochemical amino-acid (5-of-6) spaced seeding.

#[cfg(test)]
use super::append_two_hit_diagonal_lcbs;
use super::{AnchorPairs, LcbList, cluster_strand_lcbs, codon_to_aa};

const AA5_KEY_SPACE: usize = 59_049; // 9^5

struct Aa5Index {
    max_occ: usize,
    counts: Vec<u16>,
    slots: Vec<u32>,
}

impl Aa5Index {
    fn build(seq: &[u8], max_occ: usize, skip_softmasked: bool) -> Self {
        const AA_OFFSETS: [usize; 5] = [0, 1, 2, 4, 5];
        let max_occ = max_occ.min(u16::MAX as usize);
        let Some(slot_len) = (max_occ > 0)
            .then(|| AA5_KEY_SPACE.checked_mul(max_occ))
            .flatten()
        else {
            return Self {
                max_occ: 0,
                counts: Vec::new(),
                slots: Vec::new(),
            };
        };
        if seq.len() < 18 || seq.len() > u32::MAX as usize {
            return Self {
                max_occ,
                counts: Vec::new(),
                slots: Vec::new(),
            };
        }
        let mut counts = vec![0_u16; AA5_KEY_SPACE];
        let mut slots = vec![0_u32; slot_len];
        for frame in 0..3 {
            if seq.len() < frame + 18 {
                continue;
            }
            let mut pos = frame;
            while pos + 18 <= seq.len() {
                let mut key = 0_usize;
                let mut valid = true;
                for &c in &AA_OFFSETS {
                    let s = &seq[pos + c * 3..pos + c * 3 + 3];
                    if skip_softmasked && s.iter().any(u8::is_ascii_lowercase) {
                        valid = false;
                        break;
                    }
                    let Some(aa) = codon_to_aa(s[0], s[1], s[2]) else {
                        valid = false;
                        break;
                    };
                    key = key * 9 + usize::from(aa_group(aa));
                }
                if valid {
                    let c = counts[key] as usize;
                    if c < max_occ {
                        slots[key * max_occ + c] = u32::try_from(pos).unwrap_or(u32::MAX);
                    }
                    counts[key] = counts[key].saturating_add(1);
                }
                pos += 3;
            }
        }
        Self {
            max_occ,
            counts,
            slots,
        }
    }

    fn intersect_into(
        &self,
        other: &Self,
        rev_lens: Option<(usize, usize)>,
        out: &mut AnchorPairs,
    ) {
        if self.counts.is_empty() || other.counts.is_empty() {
            return;
        }
        let mo = self.max_occ;
        for key in 0..AA5_KEY_SPACE {
            let cx = self.counts[key] as usize;
            let cy = other.counts[key] as usize;
            if cx > 0 && cx <= mo && cy > 0 && cy <= other.max_occ {
                let xs = &self.slots[key * mo..key * mo + cx];
                let ys = &other.slots[key * other.max_occ..key * other.max_occ + cy];
                if let Some((lx, ly)) = rev_lens {
                    for &rx in xs {
                        let Some(px) = lx.checked_sub(18 + (rx as usize)) else {
                            continue;
                        };
                        for &ry in ys {
                            let Some(py) = ly.checked_sub(18 + (ry as usize)) else {
                                continue;
                            };
                            out.push((px, py));
                        }
                    }
                } else {
                    for &px in xs {
                        for &py in ys {
                            out.push((px as usize, py as usize));
                        }
                    }
                }
            }
        }
    }
}

const fn aa_group(ax: u8) -> u8 {
    const GROUP: [u8; 20] = [0, 1, 2, 2, 3, 0, 4, 5, 4, 5, 5, 6, 7, 6, 4, 8, 8, 5, 3, 3];
    if (ax as usize) < 20 {
        GROUP[ax as usize]
    } else {
        9
    }
}

const fn aa_pair_score(ax: u8, ay: u8) -> i32 {
    if ax == ay {
        return 15;
    }
    let gx = aa_group(ax);
    let gy = aa_group(ay);
    if gx < 9 && gx == gy { 6 } else { -6 }
}

/// Find spaced `5-of-6` translated amino-acid (18 bp) hits across all 6 reading frames (`+`/`+` and `-`/`-` in forward coordinates).
#[cfg(test)]
pub(crate) fn find_translated_aa5_hits(
    seq_x: &[u8],
    seq_y: &[u8],
    max_occ: usize,
    skip_softmasked: bool,
) -> AnchorPairs {
    let ((x_fwd, y_fwd), (x_rev, y_rev)) = rayon::join(
        || {
            rayon::join(
                || Aa5Index::build(seq_x, max_occ, skip_softmasked),
                || Aa5Index::build(seq_y, max_occ, skip_softmasked),
            )
        },
        || {
            if seq_x.len() >= 18 && seq_y.len() >= 18 {
                let (rc_x, rc_y) = rayon::join(
                    || super::revcomp_preserve_case(seq_x),
                    || super::revcomp_preserve_case(seq_y),
                );
                let (xr, yr) = rayon::join(
                    || Aa5Index::build(&rc_x, max_occ, skip_softmasked),
                    || Aa5Index::build(&rc_y, max_occ, skip_softmasked),
                );
                (Some(xr), Some(yr))
            } else {
                (None, None)
            }
        },
    );
    let mut hits = Vec::new();
    x_fwd.intersect_into(&y_fwd, None, &mut hits);
    if let (Some(xr), Some(yr)) = (x_rev, y_rev) {
        xr.intersect_into(&yr, Some((seq_x.len(), seq_y.len())), &mut hits);
    }
    hits
}

/// Find spaced `5-of-6` translated amino-acid hits for both forward and reverse-complement
/// target orientations while indexing `seq_x` and `rc_x` only once.
pub(crate) fn find_translated_aa5_hits_both(
    seq_x: &[u8],
    seq_y: &[u8],
    seq_y_rc: &[u8],
    max_occ: usize,
    skip_softmasked: bool,
) -> (AnchorPairs, AnchorPairs) {
    let ((x_fwd, x_rev), (y_fwd, yrc_fwd)) = rayon::join(
        || {
            rayon::join(
                || Aa5Index::build(seq_x, max_occ, skip_softmasked),
                || {
                    if seq_x.len() >= 18 {
                        let rc_x = super::revcomp_preserve_case(seq_x);
                        Some(Aa5Index::build(&rc_x, max_occ, skip_softmasked))
                    } else {
                        None
                    }
                },
            )
        },
        || {
            rayon::join(
                || Aa5Index::build(seq_y, max_occ, skip_softmasked),
                || Aa5Index::build(seq_y_rc, max_occ, skip_softmasked),
            )
        },
    );

    let mut hits_fwd = Vec::new();
    x_fwd.intersect_into(&y_fwd, None, &mut hits_fwd);
    if let Some(xr) = &x_rev {
        xr.intersect_into(&yrc_fwd, Some((seq_x.len(), seq_y.len())), &mut hits_fwd);
    }

    let mut hits_rc = Vec::new();
    x_fwd.intersect_into(&yrc_fwd, None, &mut hits_rc);
    if let Some(xr) = &x_rev {
        xr.intersect_into(&y_fwd, Some((seq_x.len(), seq_y_rc.len())), &mut hits_rc);
    }

    (hits_fwd, hits_rc)
}

fn rc_codon_to_aa(seq: &[u8], end_pos: usize) -> Option<u8> {
    if end_pos < 3 || end_pos > seq.len() {
        return None;
    }
    codon_to_aa(
        crate::complement_base(seq[end_pos - 1]),
        crate::complement_base(seq[end_pos - 2]),
        crate::complement_base(seq[end_pos - 3]),
    )
}

#[allow(clippy::too_many_lines)]
fn extend_aa_ungapped_hsp_oriented(
    seq_x: &[u8],
    seq_y: &[u8],
    px: usize,
    py: usize,
    min_score: i32,
    rev_coding: bool,
) -> Option<(usize, usize, usize, usize, i32)> {
    const AA_XDROP: i32 = 30;
    if px + 18 > seq_x.len() || py + 18 > seq_y.len() {
        return None;
    }
    let (req_offsets, wildcard_off): ([usize; 5], usize) = if rev_coding {
        ([0, 1, 3, 4, 5], 6)
    } else {
        ([0, 1, 2, 4, 5], 9)
    };
    let codon_at = |seq: &[u8], start: usize| -> Option<u8> {
        if rev_coding {
            rc_codon_to_aa(seq, start + 3)
        } else {
            codon_to_aa(seq[start], seq[start + 1], seq[start + 2])
        }
    };
    let mut core = 0i32;
    for c in req_offsets {
        let ax = codon_at(seq_x, px + c * 3)?;
        let ay = codon_at(seq_y, py + c * 3)?;
        core += aa_pair_score(ax, ay);
    }
    // Include the wildcard codon if both codons are valid (non-stop); otherwise skip it cleanly.
    if let (Some(ax), Some(ay)) = (
        codon_at(seq_x, px + wildcard_off),
        codon_at(seq_y, py + wildcard_off),
    ) {
        core += aa_pair_score(ax, ay);
    }

    let mut cur = 0i32;
    let mut best_left = 0i32;
    let (mut lx, mut ly) = (px, py);
    let (mut ix, mut iy) = (px, py);
    while ix >= 3 && iy >= 3 {
        ix -= 3;
        iy -= 3;
        let (Some(ax), Some(ay)) = (codon_at(seq_x, ix), codon_at(seq_y, iy)) else {
            break;
        };
        cur += aa_pair_score(ax, ay);
        if cur > best_left {
            best_left = cur;
            lx = ix;
            ly = iy;
        } else if best_left - cur > AA_XDROP {
            break;
        }
    }

    cur = 0;
    let mut best_right = 0i32;
    let (mut rx, mut ry) = (px + 18, py + 18);
    ix = px + 18;
    iy = py + 18;
    while ix + 3 <= seq_x.len() && iy + 3 <= seq_y.len() {
        let (Some(ax), Some(ay)) = (codon_at(seq_x, ix), codon_at(seq_y, iy)) else {
            break;
        };
        cur += aa_pair_score(ax, ay);
        ix += 3;
        iy += 3;
        if cur > best_right {
            best_right = cur;
            rx = ix;
            ry = iy;
        } else if best_right - cur > AA_XDROP {
            break;
        }
    }

    let total = core + best_left + best_right;
    (total >= min_score).then_some((lx, ly, rx, ry, total))
}

/// Extend an ungapped codon-level amino-acid HSP from a 5-of-6 spaced (`18 bp`) seed,
/// evaluating both forward (`+`) and reverse-complement (`-`) coding-frame orientations.
fn extend_aa_ungapped_hsp(
    seq_x: &[u8],
    seq_y: &[u8],
    px: usize,
    py: usize,
    min_score: i32,
) -> Option<(usize, usize, usize, usize, i32)> {
    match (
        extend_aa_ungapped_hsp_oriented(seq_x, seq_y, px, py, min_score, false),
        extend_aa_ungapped_hsp_oriented(seq_x, seq_y, px, py, min_score, true),
    ) {
        (Some(a), Some(b)) => Some(if a.4 >= b.4 { a } else { b }),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    }
}

pub(crate) fn append_hsp_gated_aa_lcbs(
    target: &mut LcbList,
    seq_x: &[u8],
    seq_y: &[u8],
    hits: &[(usize, usize)],
    min_aa_score: i32,
) {
    let mut gated = Vec::new();
    let mut strong = Vec::new();
    let mut sorted_hits = hits.to_vec();
    sorted_hits.sort_unstable_by_key(|&(px, py)| ((px as i128) - (py as i128), px));
    let mut max_rx_by_diag_frame: std::collections::HashMap<(i128, usize), usize> =
        std::collections::HashMap::new();
    for (px, py) in sorted_hits {
        let key = ((px as i128) - (py as i128), px % 3);
        if max_rx_by_diag_frame.get(&key).is_some_and(|&rx| px < rx) {
            gated.push((px, py));
            continue;
        }
        if let Some((lx, ly, rx, ry, sc)) =
            extend_aa_ungapped_hsp(seq_x, seq_y, px, py, min_aa_score)
        {
            max_rx_by_diag_frame
                .entry(key)
                .and_modify(|m| *m = (*m).max(rx))
                .or_insert(rx);
            gated.push((px, py));
            if sc >= min_aa_score + 30 && rx >= lx + 24 && ry >= ly + 24 {
                strong.push(vec![(lx, ly), (rx - 18, ry - 18)]);
            }
        }
    }
    let mut candidates = cluster_strand_lcbs(&gated, 1800, 200, 2, 24, 250.0);
    candidates.extend(strong);
    super::chain::push_disjoint_lcbs(target, &[], candidates);
}

/// Find 2-hit translated amino-acid 5-of-6 (`18 bp`) anchors across all reading frames.
#[cfg(test)]
pub(crate) fn append_translated_aa5_lcbs(
    target: &mut LcbList,
    seq_x: &[u8],
    seq_y: &[u8],
    max_occ: usize,
    skip_softmasked: bool,
) {
    let hits = find_translated_aa5_hits(seq_x, seq_y, max_occ, skip_softmasked);
    append_two_hit_diagonal_lcbs(target, &hits, 120, 12);
}

#[cfg(test)]
mod tests {
    use super::append_hsp_gated_aa_lcbs;
    use crate::seeding::LcbList;

    #[test]
    fn aa_5of6_tolerates_stop_codon_at_wildcard_codon3() {
        let x = b"ATGGCTAGCGATCGATCGATGGCTAGCGATCGATCGATGGCTAGCGATCGATCGATGGCTAGCGATCGATCGATGGCTAGCGATCGATCG".to_vec();
        let mut y = x.clone();
        y[9..12].copy_from_slice(b"TAA");
        let mut lcbs: LcbList = Vec::new();
        append_hsp_gated_aa_lcbs(&mut lcbs, &x, &y, &[(0, 0)], 30);
        assert_ne!(lcbs.len(), 0);
    }

    #[test]
    fn aa_5of6_extends_minus_strand_coding_with_plus_strand_stops() {
        // Antisense of Leu (TTA -> TAA on +) and Ser (TCA -> TGA on +) would hit + stop codons,
        // so only reverse-complement codon translation succeeds.
        let coding_x = b"ATGGCTCTGTTATCATCAGCTCTGTTATCATCAGCTCTGTTATCATCAATGGCT".to_vec();
        let coding_y = b"ATGGCCTTGCTCTCGTCGGCCCTGCTCTCGTCGGCCCTGCTCTCGTCGATGGCT".to_vec();
        let rc_x = crate::revcomp(&coding_x);
        let rc_y = crate::revcomp(&coding_y);
        let (hits_fwd, _) = super::find_translated_aa5_hits_both(&rc_x, &rc_y, &coding_y, 2, false);
        assert_ne!(hits_fwd.len(), 0);
        let mut lcbs: LcbList = Vec::new();
        append_hsp_gated_aa_lcbs(&mut lcbs, &rc_x, &rc_y, &hits_fwd, 30);
        assert!(
            !lcbs.is_empty(),
            "expected minus-strand coding HSP to pass gating"
        );
    }

    #[test]
    fn aa_5of6_minus_strand_tolerates_wildcard_at_fwd_offset_6() {
        let coding = b"ATGGCTCTGTTATCATCAGCTCTGTTATCATCA".to_vec();
        let rc_x = crate::revcomp(&coding);
        let mut rc_y = rc_x.clone();
        // Forward offset 6..9 is the reverse-complement wildcard codon (rc codon index 3).
        rc_y[6..9].copy_from_slice(b"NNN");
        let mut lcbs: LcbList = Vec::new();
        append_hsp_gated_aa_lcbs(&mut lcbs, &rc_x, &rc_y, &[(0, 0)], 20);
        assert!(
            !lcbs.is_empty(),
            "expected minus-strand coding HSP to tolerate wildcard at px+6..9"
        );
    }

    #[test]
    fn aa_5of6_handles_short_and_empty_sequences_without_panic() {
        let rc_empty = crate::revcomp(b"");
        let (fwd, rev) = super::find_translated_aa5_hits_both(b"ACGT", b"", &rc_empty, 4, false);
        assert_eq!(fwd, []);
        assert_eq!(rev, []);

        let aln_short = crate::align_contig_pair(b"ACGTACGT", b"ACGT").expect("short align");
        assert_eq!(aln_short.aligned_pairs, []);

        let aln_empty = crate::align_contig_pair(b"", b"ACGTACGT").expect("empty align");
        assert_eq!(aln_empty.aligned_pairs, []);
    }
}
