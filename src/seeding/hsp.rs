//! DNA ungapped HSP extension and `w=12` LCB helpers.

use super::{AnchorPairs, ContigIndex, LcbList, cluster_strand_lcbs};
use crate::bocpd::base_to_idx;
use rayon::prelude::*;

#[inline]
fn hsp_pair_score(a: u8, b: u8) -> i32 {
    let (Some(ia), Some(ib)) = (base_to_idx(a), base_to_idx(b)) else {
        return -15;
    };
    if ia == ib {
        10
    } else if (ia ^ ib) == 2 {
        -3 // transition (A<->G, C<->T)
    } else {
        -12 // transversion
    }
}

/// Extend a 19-bp seed hit `(px, py)` without gaps using transition-aware X-drop;
/// returns `Some((start_x, start_y, end_x, end_y, score))` when `score >= min_score`.
fn extend_ungapped_hsp(
    seq_x: &[u8],
    seq_y: &[u8],
    px: usize,
    py: usize,
    min_score: i32,
) -> Option<(usize, usize, usize, usize, i32)> {
    const XDROP: i32 = 60;
    if px + 19 > seq_x.len() || py + 19 > seq_y.len() {
        return None;
    }
    let mut core = 0i32;
    for off in 0..19 {
        core += hsp_pair_score(seq_x[px + off], seq_y[py + off]);
    }

    let mut cur = 0i32;
    let mut best_left = 0i32;
    let mut lx = px;
    let mut ly = py;
    let mut ix = px;
    let mut iy = py;
    while ix > 0 && iy > 0 {
        ix -= 1;
        iy -= 1;
        cur += hsp_pair_score(seq_x[ix], seq_y[iy]);
        if cur > best_left {
            best_left = cur;
            lx = ix;
            ly = iy;
        } else if best_left - cur > XDROP {
            break;
        }
    }

    cur = 0;
    let mut best_right = 0i32;
    let mut rx = px + 19;
    let mut ry = py + 19;
    ix = px + 19;
    iy = py + 19;
    while ix < seq_x.len() && iy < seq_y.len() {
        cur += hsp_pair_score(seq_x[ix], seq_y[iy]);
        ix += 1;
        iy += 1;
        if cur > best_right {
            best_right = cur;
            rx = ix;
            ry = iy;
        } else if best_right - cur > XDROP {
            break;
        }
    }

    let total = core + best_left + best_right;
    (total >= min_score).then_some((lx, ly, rx, ry, total))
}

pub(crate) fn append_hsp_gated_lastz12_lcbs(
    target: &mut LcbList,
    seq_x: &[u8],
    seq_y: &[u8],
    hits: &[(usize, usize)],
    min_hsp_score: i32,
) {
    #[allow(clippy::cast_precision_loss)]
    let space_bits = ((seq_x.len().max(1) as f64) * (seq_y.len().max(1) as f64))
        .log2()
        .max(10.0);
    #[allow(clippy::cast_possible_truncation)]
    let strong_margin = (1.8 * space_bits).round() as i32;
    let mut filtered_hits: AnchorPairs = Vec::new();
    let mut strong_hsps: LcbList = Vec::new();
    let eval_hit = |&(px, py): &(usize, usize)| {
        let (lx, ly, rx, ry, sc) = extend_ungapped_hsp(seq_x, seq_y, px, py, min_hsp_score)?;
        let strong = (sc >= min_hsp_score + strong_margin && rx > lx + 19 && ry > ly + 19)
            .then_some(((lx, ly), (rx - 19, ry - 19)));
        Some((px, py, strong))
    };
    let evaluated: Vec<_> = if hits.len() >= 2048 {
        hits.par_iter().filter_map(eval_hit).collect()
    } else {
        hits.iter().filter_map(eval_hit).collect()
    };
    filtered_hits.reserve(evaluated.len());
    for (px, py, strong) in evaluated {
        filtered_hits.push((px, py));
        if let Some((p0, p1)) = strong {
            strong_hsps.push(vec![p0, p1]);
        }
    }
    let mut candidates = cluster_strand_lcbs(&filtered_hits, 1800, 200, 3, 30, 250.0);
    candidates.extend(cluster_strand_lcbs(&filtered_hits, 150, 25, 2, 16, 150.0));
    candidates.extend(strong_hsps);
    let mut intervals = super::IntervalSet::from_lcbs(target.iter());
    for lcb in candidates {
        let (sx, ex) = (lcb[0].0, lcb[lcb.len() - 1].0);
        if intervals.insert_if_disjoint(sx, ex) {
            target.push(lcb);
        }
    }
}

pub(crate) fn find_w12_raw_lcbs(idx_x: &ContigIndex, idx_y: &ContigIndex) -> (LcbList, LcbList) {
    let (hits_fwd, hits_rc) = idx_x.seeds.intersect(&idx_y.seeds);
    let raw_fwd = cluster_strand_lcbs(&hits_fwd, 2500, 250, 7, 70, 220.0);
    let raw_rc = cluster_strand_lcbs(&hits_rc, 2500, 250, 7, 70, 220.0);
    (raw_fwd, raw_rc)
}

pub(crate) fn append_nonoverlapping_lcbs(
    target: &mut LcbList,
    existing: &[Vec<(usize, usize)>],
    hits: &[(usize, usize)],
    min_anchors: usize,
    min_span: usize,
) {
    let mut candidates = cluster_strand_lcbs(hits, 1800, 200, min_anchors, min_span, 250.0);
    // Tight 2-hit diagonal pass (LASTZ-style): 2 nearby hits within 120bp and <=25bp drift.
    candidates.extend(cluster_strand_lcbs(hits, 120, 25, 2, 16, 120.0));
    let mut intervals = super::IntervalSet::from_lcbs(existing.iter().chain(target.iter()));
    for lcb in candidates {
        let (sx, ex) = (lcb[0].0, lcb[lcb.len() - 1].0);
        if intervals.insert_if_disjoint(sx, ex) {
            target.push(lcb);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::append_hsp_gated_lastz12_lcbs;
    use crate::seeding::LcbList;

    #[test]
    fn hsp_gated_lastz12_parallel_branch_handles_2048_plus_hits() {
        let unit = b"ACGTTGCATGCAACTGGTCAGTCGATCGTAGCTAGCTGAT";
        let mut seq = Vec::with_capacity(2400);
        while seq.len() < 2400 {
            seq.extend_from_slice(unit);
        }
        let hits: Vec<(usize, usize)> = (0..2200).map(|i| (i, i)).collect();
        let mut lcbs: LcbList = Vec::new();
        append_hsp_gated_lastz12_lcbs(&mut lcbs, &seq, &seq, &hits, 120);
        assert_ne!(lcbs.len(), 0);
    }
}
