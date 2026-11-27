//! DNA ungapped HSP extension and `w=12` LCB helpers.

use super::{AnchorPairs, ContigIndex, LcbList, cluster_strand_lcbs};
use crate::bocpd::IDENTITY_PRIOR;
use crate::dna::base_to_idx;
use rayon::prelude::*;

/// Integer scaling factor mapping bit scores (`log2` likelihood ratios) into K80 ungapped HSP units.
pub(crate) const HSP_SCALE: f64 = 7.25;
const HSP_KAPPA: f64 = 2.0;
const HSP_XDROP_BITS: f64 = 8.0;

#[inline]
fn k80_hsp_params() -> (i32, i32, i32, i32) {
    let p_id = IDENTITY_PRIOR.clamp(0.25, 0.999);
    let p_mismatch = 1.0 - p_id;
    let p_ts = p_mismatch * (HSP_KAPPA / (HSP_KAPPA + 2.0));
    let p_tv_each = p_mismatch * (1.0 / (HSP_KAPPA + 2.0));
    #[allow(clippy::cast_possible_truncation)]
    let match_sc = ((p_id / 0.25).log2() * HSP_SCALE).round() as i32;
    #[allow(clippy::cast_possible_truncation)]
    let transition_sc = ((p_ts / 0.25).log2() * HSP_SCALE).round() as i32;
    #[allow(clippy::cast_possible_truncation)]
    let tv_sc = ((p_tv_each / 0.25).log2() * HSP_SCALE).round() as i32;
    #[allow(clippy::cast_possible_truncation)]
    let xdrop = (HSP_XDROP_BITS * HSP_SCALE).round() as i32;
    (match_sc, transition_sc, tv_sc, xdrop)
}

#[inline]
fn hsp_pair_score(a: u8, b: u8, match_sc: i32, transition_sc: i32, tv_sc: i32) -> i32 {
    let (Some(ia), Some(ib)) = (base_to_idx(a), base_to_idx(b)) else {
        return tv_sc - 3;
    };
    if ia == ib {
        match_sc
    } else if (ia ^ ib) == 2 {
        transition_sc // transition (A<->G, C<->T)
    } else {
        tv_sc // transversion
    }
}

/// Extend a `seed_span`-bp seed hit `(px, py)` without gaps using transition-aware X-drop;
/// returns `Some((start_x, start_y, end_x, end_y, score))` when `score >= min_score`.
fn extend_ungapped_hsp(
    seq_x: &[u8],
    seq_y: &[u8],
    px: usize,
    py: usize,
    seed_span: usize,
    min_score: i32,
    hsp_params: (i32, i32, i32, i32),
) -> Option<(usize, usize, usize, usize, i32)> {
    let (match_sc, transition_sc, tv_sc, xdrop) = hsp_params;
    let seed_span = seed_span.max(1);
    if px.checked_add(seed_span).is_none_or(|e| e > seq_x.len())
        || py.checked_add(seed_span).is_none_or(|e| e > seq_y.len())
    {
        return None;
    }
    let mut core = 0i32;
    for off in 0..seed_span {
        core += hsp_pair_score(
            seq_x[px + off],
            seq_y[py + off],
            match_sc,
            transition_sc,
            tv_sc,
        );
    }

    let mut cur = 0i64;
    let mut best_left = 0i64;
    let mut lx = px;
    let mut ly = py;
    let mut ix = px;
    let mut iy = py;
    while ix > 0 && iy > 0 {
        ix -= 1;
        iy -= 1;
        cur += i64::from(hsp_pair_score(
            seq_x[ix],
            seq_y[iy],
            match_sc,
            transition_sc,
            tv_sc,
        ));
        if cur > best_left {
            best_left = cur;
            lx = ix;
            ly = iy;
        } else if best_left - cur > i64::from(xdrop) {
            break;
        }
    }

    cur = 0;
    let mut best_right = 0i64;
    let mut rx = px + seed_span;
    let mut ry = py + seed_span;
    ix = px + seed_span;
    iy = py + seed_span;
    while ix < seq_x.len() && iy < seq_y.len() {
        cur += i64::from(hsp_pair_score(
            seq_x[ix],
            seq_y[iy],
            match_sc,
            transition_sc,
            tv_sc,
        ));
        ix += 1;
        iy += 1;
        if cur > best_right {
            best_right = cur;
            rx = ix;
            ry = iy;
        } else if best_right - cur > i64::from(xdrop) {
            break;
        }
    }

    let total_i64 = i64::from(core)
        .saturating_add(best_left)
        .saturating_add(best_right);
    let total = i32::try_from(total_i64.clamp(i64::from(i32::MIN), i64::from(i32::MAX)))
        .unwrap_or(i32::MAX);
    (total_i64 >= i64::from(min_score)).then_some((lx, ly, rx, ry, total))
}

pub(crate) fn append_hsp_gated_lcbs(
    target: &mut LcbList,
    seq_x: &[u8],
    seq_y: &[u8],
    hits: &[(usize, usize)],
    seed_span: usize,
    min_hsp_score: i32,
) {
    let seed_span = seed_span.max(1);
    let hsp_params = k80_hsp_params();
    #[allow(clippy::cast_precision_loss)]
    let space_bits = ((seq_x.len().max(1) as f64) * (seq_y.len().max(1) as f64))
        .log2()
        .max(10.0);
    #[allow(clippy::cast_possible_truncation)]
    let strong_margin = (1.8 * space_bits).round() as i32;
    let mut filtered_hits: AnchorPairs = Vec::new();
    let mut strong_hsps: LcbList = Vec::new();
    let mut sorted_hits = hits.to_vec();
    sorted_hits.sort_unstable_by_key(|&(px, py)| ((px as i128) - (py as i128), px));
    let eval_diag_group = |group: &[(usize, usize)]| {
        let mut local_hits = Vec::new();
        let mut local_strong = Vec::new();
        let mut max_rx = 0usize;
        for &(px, py) in group {
            if px < max_rx {
                local_hits.push((px, py));
                continue;
            }
            if let Some((lx, ly, rx, ry, sc)) =
                extend_ungapped_hsp(seq_x, seq_y, px, py, seed_span, min_hsp_score, hsp_params)
            {
                max_rx = max_rx.max(rx);
                local_hits.push((px, py));
                if sc >= min_hsp_score.saturating_add(strong_margin)
                    && rx > lx + seed_span
                    && ry > ly + seed_span
                {
                    local_strong.push(vec![(lx, ly), (rx - seed_span, ry - seed_span)]);
                }
            }
        }
        (local_hits, local_strong)
    };
    let diag_chunks: Vec<&[(usize, usize)]> = sorted_hits
        .chunk_by(|a, b| ((a.0 as i128) - (a.1 as i128)) == ((b.0 as i128) - (b.1 as i128)))
        .collect();
    let evaluated: Vec<_> = if sorted_hits.len() >= 2048 {
        diag_chunks.par_iter().map(|g| eval_diag_group(g)).collect()
    } else {
        diag_chunks.iter().map(|g| eval_diag_group(g)).collect()
    };
    for (lh, ls) in evaluated {
        filtered_hits.extend(lh);
        strong_hsps.extend(ls);
    }
    let mut candidates = cluster_strand_lcbs(&filtered_hits, 1800, 200, 3, 30, 250.0);
    candidates.extend(cluster_strand_lcbs(&filtered_hits, 150, 25, 2, 16, 150.0));
    candidates.extend(strong_hsps);
    super::chain::push_disjoint_lcbs(target, candidates);
}

pub(crate) fn find_w12_raw_lcbs(idx_x: &ContigIndex, idx_y: &ContigIndex) -> (LcbList, LcbList) {
    let (hits_fwd, hits_rc) = idx_x.seeds.intersect(&idx_y.seeds);
    let raw_fwd = cluster_strand_lcbs(&hits_fwd, 2500, 250, 7, 70, 220.0);
    let raw_rc = cluster_strand_lcbs(&hits_rc, 2500, 250, 7, 70, 220.0);
    (raw_fwd, raw_rc)
}

pub(crate) fn append_nonoverlapping_lcbs(
    target: &mut LcbList,
    hits: &[(usize, usize)],
    min_anchors: usize,
    min_span: usize,
) {
    let mut candidates = cluster_strand_lcbs(hits, 1800, 200, min_anchors, min_span, 250.0);
    // Tight 2-hit diagonal pass (LASTZ-style): 2 nearby hits within 120bp and <=25bp drift.
    candidates.extend(cluster_strand_lcbs(hits, 120, 25, 2, 16, 120.0));
    super::chain::push_disjoint_lcbs(target, candidates);
}

#[cfg(test)]
mod tests {
    use super::append_hsp_gated_lcbs;
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
        append_hsp_gated_lcbs(&mut lcbs, &seq, &seq, &hits, 19, 120);
        assert_ne!(lcbs.len(), 0);
    }

    #[test]
    fn k80_hsp_params_match_expected_ordering_and_signs() {
        let (match_sc, transition_sc, tv_sc, xdrop) = super::k80_hsp_params();
        assert!(match_sc > 0);
        assert!(transition_sc < 0 && tv_sc < transition_sc);
        assert!(xdrop > match_sc);
    }

    #[test]
    fn hsp_threshold_i32_max_does_not_overflow_or_wrap() {
        let seq = b"ACGTACGTACGTTGCATGCATGCATGCAACGTACGTACGT";
        let hits = vec![(0usize, 0usize), (15, 15)];
        let mut out = Vec::new();
        append_hsp_gated_lcbs(&mut out, seq, seq, &hits, 12, i32::MAX);
        assert_eq!(out.len(), 0);
    }
}
