//! Multi-shape fallback seeding and inter-LCB sparse-gap rescue for contig-pair alignment.

use super::AlignConfig;
use crate::AlignError;
use crate::seeding::{
    AnchorPairs, ContigIndex, FlatSeedIndex, LcbList, MASK_LASTZ11, MASK_LASTZ12, MASK_PH11_18,
    MASK_RY20, MASK_W9, MASK_W10, append_hsp_gated_aa_lcbs, append_hsp_gated_lastz12_lcbs,
    append_nonoverlapping_lcbs, append_two_hit_diagonal_lcbs, find_lastz11_1transition_hits_flat,
    find_lastz12_1transition_hits_flat, find_lastz12_2transition_hits_flat,
    find_translated_aa5_hits_both, find_w10_1transition_hits_flat, local_w9_subchain,
    merge_colinear_lcbs,
};
use std::collections::HashMap;
use std::time::{Duration, Instant};

type HitPair = (AnchorPairs, AnchorPairs);

pub(super) struct FallbackOutput {
    pub(super) lcbs_fwd: LcbList,
    pub(super) lcbs_rc: LcbList,
    pub(super) spans_fwd: HashMap<(usize, usize), usize>,
    pub(super) spans_rc: HashMap<(usize, usize), usize>,
    pub(super) t_w10: Duration,
    pub(super) t_lz1: Duration,
    pub(super) t_seed_end: Duration,
}

fn record_spans(map: &mut HashMap<(usize, usize), usize>, lcbs: &[Vec<(usize, usize)>], s: usize) {
    for lcb in lcbs {
        for &pt in lcb {
            map.entry(pt).or_insert(s);
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn collect_fallback_lcbs(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
    raw_fwd: &[Vec<(usize, usize)>],
    raw_rc: &[Vec<(usize, usize)>],
    covered_x: &[bool],
    covered_bp: usize,
    primary_cov: f64,
    min_fallback_cov: f64,
    base_hsp_thresh: i32,
    hsp_bonus: i32,
    aa_bonus: i32,
    cfg: &AlignConfig,
    t_start: Instant,
) -> Result<FallbackOutput, AlignError> {
    let max_local_rescue_bp = super::MAX_LOCAL_RESCUE_BP;

    let filter_uncovered = |hits: &mut Vec<(usize, usize)>| {
        if covered_bp > 0 {
            hits.retain(|&(px, _)| !covered_x.get(px).copied().unwrap_or(false));
        }
    };

    let mut combined_fwd = raw_fwd.to_vec();
    let mut combined_rc = raw_rc.to_vec();
    let mut spans_fwd: HashMap<(usize, usize), usize> = HashMap::new();
    let mut spans_rc: HashMap<(usize, usize), usize> = HashMap::new();

    // Sparse-exon fallback: w=10 exact seeds in intervals not covered by w=12 LCBs
    // (skip when w=12 already covers >= min_fallback_cov of the reference contig).
    let y10_idx_opt = if primary_cov < min_fallback_cov {
        Some(collect_w10_lcbs(
            idx_x,
            idx_y,
            &filter_uncovered,
            &mut combined_fwd,
            &mut combined_rc,
            &mut spans_fwd,
            &mut spans_rc,
            cfg,
        )?)
    } else {
        None
    };
    let t_w10 = t_start.elapsed();

    // LASTZ-style 12-of-19 + 1-transition hits gated by ungapped X-drop HSP score
    // (skip when w=12 already covers >= min_fallback_cov of the reference contig).
    let y_lz12 = if primary_cov < min_fallback_cov {
        Some(collect_lz12_lcbs(
            idx_x,
            idx_y,
            &filter_uncovered,
            &mut combined_fwd,
            &mut combined_rc,
            &mut spans_fwd,
            &mut spans_rc,
            base_hsp_thresh + hsp_bonus,
            cfg,
        )?)
    } else {
        None
    };
    let t_lz1 = t_start.elapsed();

    let mut cached_ry_hits = None;
    if cfg.use_multi_shape_seeds && primary_cov < min_fallback_cov {
        cached_ry_hits = collect_multi_shape_lcbs(
            idx_x,
            idx_y,
            y_lz12.as_ref(),
            y10_idx_opt.as_ref(),
            &filter_uncovered,
            &mut combined_fwd,
            &mut combined_rc,
            &mut spans_fwd,
            &mut spans_rc,
            base_hsp_thresh,
            hsp_bonus,
            aa_bonus,
            cfg,
        )?;
    }

    // Transition-tolerant R/Y (1-bit purine/pyrimidine) fallback for divergent homologs.
    if cfg.use_ry_seeds {
        collect_ry_lcbs(
            idx_x,
            idx_y,
            cached_ry_hits,
            &mut combined_fwd,
            &mut combined_rc,
            &mut spans_fwd,
            &mut spans_rc,
            cfg,
        )?;
    }

    let lcbs_fwd = merge_colinear_lcbs(&combined_fwd, max_local_rescue_bp);
    let lcbs_rc = merge_colinear_lcbs(&combined_rc, max_local_rescue_bp);
    let t_seed_end = t_start.elapsed();

    Ok(FallbackOutput {
        lcbs_fwd,
        lcbs_rc,
        spans_fwd,
        spans_rc,
        t_w10,
        t_lz1,
        t_seed_end,
    })
}

#[allow(clippy::too_many_arguments)]
fn collect_w10_lcbs<F: Fn(&mut AnchorPairs)>(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
    filter_uncovered: &F,
    combined_fwd: &mut LcbList,
    combined_rc: &mut LcbList,
    spans_fwd: &mut HashMap<(usize, usize), usize>,
    spans_rc: &mut HashMap<(usize, usize), usize>,
    cfg: &AlignConfig,
) -> Result<FlatSeedIndex, AlignError> {
    let x10_idx = FlatSeedIndex::build_with_sdust(
        &idx_x.seq,
        MASK_W10,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
        Some(&idx_x.sdust),
    )?;
    let y10_idx = FlatSeedIndex::build_with_sdust(
        &idx_y.seq,
        MASK_W10,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
        Some(&idx_y.sdust),
    )?;
    let (mut hits10_fwd, mut hits10_rc) = x10_idx.intersect(&y10_idx);
    filter_uncovered(&mut hits10_fwd);
    filter_uncovered(&mut hits10_rc);
    let prev_fwd = combined_fwd.len();
    append_nonoverlapping_lcbs(
        combined_fwd,
        &[],
        &hits10_fwd,
        cfg.fallback_min_anchors,
        cfg.fallback_min_span,
    );
    record_spans(spans_fwd, &combined_fwd[prev_fwd..], MASK_W10.len());
    let prev_rc = combined_rc.len();
    append_nonoverlapping_lcbs(
        combined_rc,
        &[],
        &hits10_rc,
        cfg.fallback_min_anchors,
        cfg.fallback_min_span,
    );
    record_spans(spans_rc, &combined_rc[prev_rc..], MASK_W10.len());
    Ok(y10_idx)
}

#[allow(clippy::too_many_arguments)]
fn collect_lz12_lcbs<F: Fn(&mut AnchorPairs)>(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
    filter_uncovered: &F,
    combined_fwd: &mut LcbList,
    combined_rc: &mut LcbList,
    spans_fwd: &mut HashMap<(usize, usize), usize>,
    spans_rc: &mut HashMap<(usize, usize), usize>,
    thresh: i32,
    cfg: &AlignConfig,
) -> Result<FlatSeedIndex, AlignError> {
    let idx = FlatSeedIndex::build_with_sdust(
        &idx_y.seq,
        MASK_LASTZ12,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
        Some(&idx_y.sdust),
    )?;
    let (mut hits_lz_fwd, mut hits_lz_rc) = find_lastz12_1transition_hits_flat(
        &idx_x.seq,
        Some(&idx_x.sdust),
        &idx,
        cfg.skip_softmasked,
    );
    filter_uncovered(&mut hits_lz_fwd);
    filter_uncovered(&mut hits_lz_rc);
    append_hsp_pair_lcbs(
        &idx_x.seq,
        &idx_y.seq,
        &idx_y.seq_rc,
        &hits_lz_fwd,
        &hits_lz_rc,
        thresh,
        MASK_LASTZ12.len(),
        combined_fwd,
        combined_rc,
        spans_fwd,
        spans_rc,
    );
    Ok(idx)
}

#[allow(clippy::too_many_arguments)]
fn append_hsp_pair_lcbs(
    seq_x: &[u8],
    seq_y: &[u8],
    seq_y_rc: &[u8],
    hits_fwd: &[(usize, usize)],
    hits_rc: &[(usize, usize)],
    thresh: i32,
    span_len: usize,
    combined_fwd: &mut LcbList,
    combined_rc: &mut LcbList,
    spans_fwd: &mut HashMap<(usize, usize), usize>,
    spans_rc: &mut HashMap<(usize, usize), usize>,
) {
    let prev_fwd = combined_fwd.len();
    append_hsp_gated_lastz12_lcbs(combined_fwd, seq_x, seq_y, hits_fwd, thresh);
    record_spans(spans_fwd, &combined_fwd[prev_fwd..], span_len);
    let prev_rc = combined_rc.len();
    append_hsp_gated_lastz12_lcbs(combined_rc, seq_x, seq_y_rc, hits_rc, thresh);
    record_spans(spans_rc, &combined_rc[prev_rc..], span_len);
}

#[allow(clippy::too_many_arguments)]
fn collect_multi_shape_lcbs<F: Fn(&mut AnchorPairs)>(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
    y_lz12: Option<&FlatSeedIndex>,
    y10_idx_opt: Option<&FlatSeedIndex>,
    filter_uncovered: &F,
    combined_fwd: &mut LcbList,
    combined_rc: &mut LcbList,
    spans_fwd: &mut HashMap<(usize, usize), usize>,
    spans_rc: &mut HashMap<(usize, usize), usize>,
    base_hsp_thresh: i32,
    hsp_bonus: i32,
    aa_bonus: i32,
    cfg: &AlignConfig,
) -> Result<Option<HitPair>, AlignError> {
    let seq_x = &idx_x.seq;
    let seq_y = &idx_y.seq;
    let seq_y_rc = &idx_y.seq_rc;
    let search_cells = (seq_x.len() as u128) * (seq_y.len() as u128);
    if search_cells < 1_000_000_000
        && let Some(lz_idx) = y_lz12
    {
        let (mut f, mut r) = find_lastz12_2transition_hits_flat(
            seq_x,
            Some(&idx_x.sdust),
            lz_idx,
            cfg.skip_softmasked,
        );
        filter_uncovered(&mut f);
        filter_uncovered(&mut r);
        append_hsp_pair_lcbs(
            seq_x,
            seq_y,
            seq_y_rc,
            &f,
            &r,
            base_hsp_thresh + 30 + hsp_bonus,
            MASK_LASTZ12.len(),
            combined_fwd,
            combined_rc,
            spans_fwd,
            spans_rc,
        );
    }

    let y_lz11 = FlatSeedIndex::build_with_sdust(
        seq_y,
        MASK_LASTZ11,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
        Some(&idx_y.sdust),
    )?;
    let (mut f11, mut r11) =
        find_lastz11_1transition_hits_flat(seq_x, Some(&idx_x.sdust), &y_lz11, cfg.skip_softmasked);
    filter_uncovered(&mut f11);
    filter_uncovered(&mut r11);
    append_hsp_pair_lcbs(
        seq_x,
        seq_y,
        seq_y_rc,
        &f11,
        &r11,
        base_hsp_thresh + 30 + hsp_bonus,
        MASK_LASTZ11.len(),
        combined_fwd,
        combined_rc,
        spans_fwd,
        spans_rc,
    );

    collect_ph11_lcbs(
        idx_x,
        idx_y,
        filter_uncovered,
        combined_fwd,
        combined_rc,
        spans_fwd,
        spans_rc,
        cfg,
    )?;

    let cached_ry_hits = if search_cells < 1_000_000_000 {
        collect_ry_and_w10_1trans_lcbs(
            idx_x,
            idx_y,
            y10_idx_opt,
            filter_uncovered,
            combined_fwd,
            combined_rc,
            spans_fwd,
            spans_rc,
            base_hsp_thresh + 40 + hsp_bonus,
            cfg,
        )?
    } else {
        None
    };

    collect_aa5_lcbs(
        idx_x,
        idx_y,
        filter_uncovered,
        combined_fwd,
        combined_rc,
        spans_fwd,
        spans_rc,
        (base_hsp_thresh / 2) + aa_bonus,
        cfg,
    );

    Ok(cached_ry_hits)
}

#[allow(clippy::too_many_arguments)]
fn collect_ph11_lcbs<F: Fn(&mut AnchorPairs)>(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
    filter_uncovered: &F,
    combined_fwd: &mut LcbList,
    combined_rc: &mut LcbList,
    spans_fwd: &mut HashMap<(usize, usize), usize>,
    spans_rc: &mut HashMap<(usize, usize), usize>,
    cfg: &AlignConfig,
) -> Result<(), AlignError> {
    let x_ph = FlatSeedIndex::build_with_sdust(
        &idx_x.seq,
        MASK_PH11_18,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
        Some(&idx_x.sdust),
    )?;
    let y_ph = FlatSeedIndex::build_with_sdust(
        &idx_y.seq,
        MASK_PH11_18,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
        Some(&idx_y.sdust),
    )?;
    let (mut ph_fwd, mut ph_rc) = x_ph.intersect(&y_ph);
    filter_uncovered(&mut ph_fwd);
    filter_uncovered(&mut ph_rc);
    let pf = combined_fwd.len();
    append_two_hit_diagonal_lcbs(combined_fwd, &ph_fwd, 64, 4);
    record_spans(spans_fwd, &combined_fwd[pf..], MASK_PH11_18.len());
    let pr = combined_rc.len();
    append_two_hit_diagonal_lcbs(combined_rc, &ph_rc, 64, 4);
    record_spans(spans_rc, &combined_rc[pr..], MASK_PH11_18.len());
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn collect_ry_and_w10_1trans_lcbs<F: Fn(&mut AnchorPairs)>(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
    y10_idx_opt: Option<&FlatSeedIndex>,
    filter_uncovered: &F,
    combined_fwd: &mut LcbList,
    combined_rc: &mut LcbList,
    spans_fwd: &mut HashMap<(usize, usize), usize>,
    spans_rc: &mut HashMap<(usize, usize), usize>,
    thresh: i32,
    cfg: &AlignConfig,
) -> Result<Option<HitPair>, AlignError> {
    let x_ry = FlatSeedIndex::build_with_sdust(
        &idx_x.seq,
        MASK_RY20,
        2,
        cfg.skip_softmasked,
        true,
        Some(&idx_x.sdust),
    )?;
    let y_ry = FlatSeedIndex::build_with_sdust(
        &idx_y.seq,
        MASK_RY20,
        2,
        cfg.skip_softmasked,
        true,
        Some(&idx_y.sdust),
    )?;
    let (mut fry, mut rry) = x_ry.intersect(&y_ry);
    filter_uncovered(&mut fry);
    filter_uncovered(&mut rry);
    let cached = cfg.use_ry_seeds.then(|| (fry.clone(), rry.clone()));
    append_hsp_pair_lcbs(
        &idx_x.seq,
        &idx_y.seq,
        &idx_y.seq_rc,
        &fry,
        &rry,
        thresh,
        MASK_RY20.len(),
        combined_fwd,
        combined_rc,
        spans_fwd,
        spans_rc,
    );
    if let Some(y10_idx) = y10_idx_opt {
        let (mut fw10, mut rw10) = find_w10_1transition_hits_flat(
            &idx_x.seq,
            Some(&idx_x.sdust),
            y10_idx,
            cfg.skip_softmasked,
        );
        filter_uncovered(&mut fw10);
        filter_uncovered(&mut rw10);
        append_hsp_pair_lcbs(
            &idx_x.seq,
            &idx_y.seq,
            &idx_y.seq_rc,
            &fw10,
            &rw10,
            thresh,
            MASK_W10.len(),
            combined_fwd,
            combined_rc,
            spans_fwd,
            spans_rc,
        );
    }
    Ok(cached)
}

const AA5_ANCHOR_SPAN_BP: usize = 18;

#[allow(clippy::too_many_arguments)]
fn collect_aa5_lcbs<F: Fn(&mut AnchorPairs)>(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
    filter_uncovered: &F,
    combined_fwd: &mut LcbList,
    combined_rc: &mut LcbList,
    spans_fwd: &mut HashMap<(usize, usize), usize>,
    spans_rc: &mut HashMap<(usize, usize), usize>,
    thresh: i32,
    cfg: &AlignConfig,
) {
    let (mut aa_fwd, mut aa_rc) = find_translated_aa5_hits_both(
        &idx_x.seq,
        &idx_y.seq,
        &idx_y.seq_rc,
        4,
        cfg.skip_softmasked,
    );
    filter_uncovered(&mut aa_fwd);
    filter_uncovered(&mut aa_rc);
    let pf = combined_fwd.len();
    append_hsp_gated_aa_lcbs(combined_fwd, &idx_x.seq, &idx_y.seq, &aa_fwd, thresh);
    record_spans(spans_fwd, &combined_fwd[pf..], AA5_ANCHOR_SPAN_BP);
    let pr = combined_rc.len();
    append_hsp_gated_aa_lcbs(combined_rc, &idx_x.seq, &idx_y.seq_rc, &aa_rc, thresh);
    record_spans(spans_rc, &combined_rc[pr..], AA5_ANCHOR_SPAN_BP);
}

#[allow(clippy::too_many_arguments)]
fn collect_ry_lcbs(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
    cached_ry_hits: Option<HitPair>,
    combined_fwd: &mut LcbList,
    combined_rc: &mut LcbList,
    spans_fwd: &mut HashMap<(usize, usize), usize>,
    spans_rc: &mut HashMap<(usize, usize), usize>,
    cfg: &AlignConfig,
) -> Result<(), AlignError> {
    let (hits_ry_fwd, hits_ry_rc) = if let Some(hits) = cached_ry_hits {
        hits
    } else {
        let x_ry = FlatSeedIndex::build_with_sdust(
            &idx_x.seq,
            MASK_RY20,
            2,
            cfg.skip_softmasked,
            true,
            Some(&idx_x.sdust),
        )?;
        let y_ry = FlatSeedIndex::build_with_sdust(
            &idx_y.seq,
            MASK_RY20,
            2,
            cfg.skip_softmasked,
            true,
            Some(&idx_y.sdust),
        )?;
        x_ry.intersect(&y_ry)
    };
    let prev_fwd = combined_fwd.len();
    append_nonoverlapping_lcbs(
        combined_fwd,
        &[],
        &hits_ry_fwd,
        cfg.fallback_min_anchors,
        cfg.fallback_min_span,
    );
    record_spans(spans_fwd, &combined_fwd[prev_fwd..], MASK_RY20.len());
    let prev_rc = combined_rc.len();
    append_nonoverlapping_lcbs(
        combined_rc,
        &[],
        &hits_ry_rc,
        cfg.fallback_min_anchors,
        cfg.fallback_min_span,
    );
    record_spans(spans_rc, &combined_rc[prev_rc..], MASK_RY20.len());
    Ok(())
}

fn push_rescued(
    rescued_inter: &mut LcbList,
    anchor_spans: &mut HashMap<(usize, usize), usize>,
    sub: Vec<(usize, usize)>,
    mask_len: usize,
) {
    for &pt in &sub {
        anchor_spans.entry(pt).or_insert(mask_len);
    }
    rescued_inter.push(sub);
}

pub(super) fn rescue_inter_lcb_gaps(
    seq_x: &[u8],
    y_work: &[u8],
    lcb_list: &[Vec<(usize, usize)>],
    anchor_spans: &mut HashMap<(usize, usize), usize>,
    default_span: usize,
) -> Result<LcbList, AlignError> {
    let mut rescued_inter = Vec::new();

    if lcb_list.is_empty() && seq_x.len() <= 2500 && y_work.len() <= 2500 {
        let sub10 = local_w9_subchain(seq_x, y_work, 0, 0, MASK_W10, 200, 1800)?;
        if sub10.len() >= 3 {
            push_rescued(&mut rescued_inter, anchor_spans, sub10, MASK_W10.len());
        } else {
            let sub9 = local_w9_subchain(seq_x, y_work, 0, 0, MASK_W9, 120, 1200)?;
            if sub9.len() >= 4 {
                push_rescued(&mut rescued_inter, anchor_spans, sub9, MASK_W9.len());
            }
        }
    }
    if let Some(first_lcb) = lcb_list.first()
        && let Some(&(fx, fy)) = first_lcb.first()
        && seq_x.len() <= 50_000
        && (100..=2500).contains(&fx)
        && (100..=2500).contains(&fy)
    {
        let sub = local_w9_subchain(&seq_x[..fx], &y_work[..fy], 0, 0, MASK_W10, 200, 1800)?;
        if sub.len() >= 3 {
            push_rescued(&mut rescued_inter, anchor_spans, sub, MASK_W10.len());
        }
    }
    if let Some(last_lcb) = lcb_list.last()
        && let Some(&last_pt) = last_lcb.last()
        && seq_x.len() <= 50_000
    {
        let sp = anchor_spans.get(&last_pt).copied().unwrap_or(default_span);
        let lx = last_pt.0 + sp;
        let ly = last_pt.1 + sp;
        if seq_x.len() > lx
            && y_work.len() > ly
            && (100..=2500).contains(&(seq_x.len() - lx))
            && (100..=2500).contains(&(y_work.len() - ly))
        {
            let sub = local_w9_subchain(&seq_x[lx..], &y_work[ly..], lx, ly, MASK_W10, 200, 1800)?;
            if sub.len() >= 3 {
                push_rescued(&mut rescued_inter, anchor_spans, sub, MASK_W10.len());
            }
        }
    }
    for w in lcb_list.windows(2) {
        let (Some(&prev_pt), Some(&next_pt)) = (w[0].last(), w[1].first()) else {
            continue;
        };
        let sp = anchor_spans.get(&prev_pt).copied().unwrap_or(default_span);
        let x0 = prev_pt.0 + sp;
        let y0 = prev_pt.1 + sp;
        let x1 = next_pt.0;
        let y1 = next_pt.1;
        if x1 > x0
            && y1 > y0
            && (100..=25_000).contains(&(x1 - x0))
            && (100..=25_000).contains(&(y1 - y0))
        {
            let sub10 =
                local_w9_subchain(&seq_x[x0..x1], &y_work[y0..y1], x0, y0, MASK_W10, 200, 1800)?;
            let min_sub10 = if x1 - x0 <= 2500 && y1 - y0 <= 2500 {
                2
            } else {
                3
            };
            if sub10.len() >= min_sub10 {
                push_rescued(&mut rescued_inter, anchor_spans, sub10, MASK_W10.len());
            } else {
                let sub9 =
                    local_w9_subchain(&seq_x[x0..x1], &y_work[y0..y1], x0, y0, MASK_W9, 120, 1200)?;
                if sub9.len() >= 4 {
                    push_rescued(&mut rescued_inter, anchor_spans, sub9, MASK_W9.len());
                }
            }
        }
    }
    Ok(rescued_inter)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rescue_inter_lcb_gaps_finds_w10_subchain_on_empty_lcbs() {
        let mut seq = Vec::with_capacity(120);
        let mut state = 0x1234_5678_9abc_def0_u64;
        for _ in 0..120 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            seq.push(b"ACGT"[((state >> 62) & 3) as usize]);
        }
        let mut spans = HashMap::new();
        let rescued = rescue_inter_lcb_gaps(&seq, &seq, &[], &mut spans, 12).expect("rescue");
        assert_ne!(rescued, [] as [Vec<(usize, usize)>; 0]);
        assert!(!spans.is_empty());
    }
}
