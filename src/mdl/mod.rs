//! Banded MDL+BOCPD dynamic programming, X-drop extension, and contig-pair alignment.
//! ## Horizon Bubbles
//! In this layer, areas where the alignment lacks sufficient confidence to align (either because a gap is too large to enter DP, or because the MDL DP abstains due to negative bit-savings) are termed **Horizon Bubbles**.
//! - Unaligned inter-LCB gaps larger than `max_gap_bp` yield *macro bubbles*.
//! - DP dropouts within a gap yield *internal bubbles*, returned via `GapResult`.

use crate::AlignError;
use crate::bocpd::MIN_POST_PROB;
use crate::seeding::{
    ContigIndex, FlatSeedIndex, MASK_LASTZ11, MASK_LASTZ12, MASK_PH11_18, MASK_PRIMARY_11_15,
    MASK_RY20, MASK_W9, MASK_W10, append_hsp_gated_aa_lcbs, append_hsp_gated_lastz12_lcbs,
    append_nonoverlapping_lcbs, append_two_hit_diagonal_lcbs, find_lastz11_1transition_hits_flat,
    find_lastz12_1transition_hits_flat, find_lastz12_2transition_hits_flat,
    find_translated_aa5_hits, find_w10_1transition_hits_flat, find_w12_raw_lcbs, local_w9_subchain,
    merge_colinear_lcbs, prepare_contig_index,
};
use rayon::prelude::*;

mod dag;
mod dp;
pub use dag::*;
pub use dp::*;

/// Karlin-Altschul / MDL random-search threshold coefficient (`~0.5 * log2(N*M)` + overhead).
const KARLIN_MDL_COEFF: f64 = 0.55;

/// Configurable parameters for [`align_contig_pair_with_config`].
#[derive(Debug, Clone)]
pub struct AlignConfig {
    /// Minimum `w=12` LCB anchor count required before running fallback seeding and DP (`0` = disabled).
    pub min_top_lcb_anchors: usize,
    /// Skip lowercase softmasked bases during `w=10` / `R/Y` seeding.
    pub skip_softmasked: bool,
    /// Enable transition-tolerant 1-bit `R/Y` (`MASK_RY20`) fallback seeding for divergent homologs.
    pub use_ry_seeds: bool,
    /// Enable complementary `11-of-18` + `12-of-19` HSP-gated multi-shape seeding.
    pub use_multi_shape_seeds: bool,
    /// Minimum anchor count for `w=10` and `R/Y` fallback LCB clustering.
    pub fallback_min_anchors: usize,
    /// Minimum genomic span (bp) for `w=10` and `R/Y` fallback LCB clustering.
    pub fallback_min_span: usize,
    /// Maximum seed occurrence count per contig for `w=10` / `R/Y` fallback seeding.
    pub fallback_max_occ: usize,
    /// Half-band width for banded MDL+BOCPD DP.
    pub band_w: usize,
    /// Maximum X-drop extension distance (bp) at LCB flanks.
    pub max_ext: usize,
    /// Minimum MDL log-likelihood gain (bits) required to accept a DP gap/flank.
    pub min_delta_bits: f64,
    pub min_fallback_cov: Option<f64>,
    pub hsp_threshold: Option<i32>,
    /// Minimum posterior homology probability required to emit a base pair.
    pub min_post_prob: f64,
}

impl Default for AlignConfig {
    fn default() -> Self {
        Self {
            min_top_lcb_anchors: 0,
            skip_softmasked: false,
            use_ry_seeds: false,
            use_multi_shape_seeds: true,
            fallback_min_anchors: 4,
            fallback_min_span: 50,
            fallback_max_occ: 2,
            band_w: 180,
            max_ext: 1500,
            min_delta_bits: 10.0,
            min_fallback_cov: None,
            hsp_threshold: None,
            min_post_prob: MIN_POST_PROB,
        }
    }
}

///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if internal masks fail validation.
pub fn align_contig_pair(seq_x: &[u8], seq_y: &[u8]) -> Result<ContigAlignment, AlignError> {
    let idx_x = prepare_contig_index(seq_x)?;
    let idx_y = prepare_contig_index(seq_y)?;
    align_contig_pair_indexed(&idx_x, &idx_y)
}

/// Count shared KP-14 seeds (on either strand) between two pre-indexed contigs.
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if internal masks fail validation.
pub fn align_contig_pair_indexed(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
) -> Result<ContigAlignment, AlignError> {
    align_contig_pair_with_config(idx_x, idx_y, &AlignConfig::default())
}

/// Align two pre-indexed contigs, pruning pairs whose largest `w=12` LCB has fewer than
/// `min_top_lcb_anchors` anchors before running fallback seeding or DP.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if internal masks fail validation.
pub fn align_contig_pair_indexed_with_min_lcb(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
    min_top_lcb_anchors: usize,
) -> Result<ContigAlignment, AlignError> {
    let cfg = AlignConfig {
        min_top_lcb_anchors,
        ..AlignConfig::default()
    };
    align_contig_pair_with_config(idx_x, idx_y, &cfg)
}

/// Align two pre-indexed contigs using an explicit [`AlignConfig`].
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if internal masks fail validation.
#[allow(clippy::too_many_lines)]
pub fn align_contig_pair_with_config(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
    cfg: &AlignConfig,
) -> Result<ContigAlignment, AlignError> {
    let span = MASK_PRIMARY_11_15.len();
    let max_gap_bp = 2400;
    let max_local_rescue_bp = 40_000;
    let seq_x = &idx_x.seq;
    let seq_y = &idx_y.seq;
    let seq_y_rc = &idx_y.seq_rc;
    let pv_len = seq_y.len();

    let profile = std::env::var_os("HORIZON_PROFILE").is_some();
    let t_start = std::time::Instant::now();
    let (raw_fwd, raw_rc) = find_w12_raw_lcbs(idx_x, idx_y);
    let t_w12 = t_start.elapsed();

    // Early exit when caller requests a minimum w=12 LCB anchor threshold.
    let top_lcb_anchors = raw_fwd
        .iter()
        .chain(raw_rc.iter())
        .map(Vec::len)
        .max()
        .unwrap_or(0);
    if cfg.min_top_lcb_anchors > 0 && top_lcb_anchors < cfg.min_top_lcb_anchors {
        return Ok(ContigAlignment {
            aligned_pairs: Vec::new(),
            lcbs_fwd: 0,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 0,
        });
    }

    #[allow(clippy::cast_precision_loss)]
    let search_space_bits = ((seq_x.len().max(1) as f64) * (seq_y.len().max(1) as f64))
        .log2()
        .max(10.0);
    // Continuous Karlin-Altschul / MDL scaling with log2(N * M), avoiding fixed size cutoffs.
    #[allow(clippy::cast_possible_truncation)]
    let hsp_bonus = ((search_space_bits - 20.0).max(0.0) * 3.8).round() as i32;
    #[allow(clippy::cast_possible_truncation)]
    let aa_bonus = ((search_space_bits - 20.0).max(0.0) * 1.7).round() as i32;
    let effective_min_bits = cfg.min_delta_bits.max(KARLIN_MDL_COEFF * search_space_bits);

    let mut covered_x = vec![false; seq_x.len()];
    let mut covered_bp = 0usize;
    for lcb in raw_fwd.iter().chain(raw_rc.iter()) {
        if let (Some(first), Some(last)) = (lcb.first(), lcb.last()) {
            let s = first.0.min(seq_x.len());
            let e = (last.0 + span).min(seq_x.len());
            for flag in &mut covered_x[s..e] {
                if !*flag {
                    *flag = true;
                    covered_bp += 1;
                }
            }
        }
    }
    #[allow(clippy::cast_precision_loss)]
    // MDL-derived default for HSP thresholds
    #[allow(clippy::cast_possible_truncation)]
    let base_hsp_thresh = cfg
        .hsp_threshold
        .unwrap_or_else(|| (100.0 + (effective_min_bits * 1.5)).round() as i32);

    // MDL-derived default fallback coverage threshold
    // If the expected bit savings of the primary coverage is less than the effective min bits
    // across the sequence, trigger fallback.
    #[allow(clippy::cast_precision_loss)]
    let min_fallback_cov = cfg.min_fallback_cov.unwrap_or_else(|| {
        let expected_bits_per_base = 2.0; // 2 bits per nucleotide match
        let bit_ratio = effective_min_bits / (seq_x.len().max(1) as f64 * expected_bits_per_base);
        // Map this to a reasonable coverage threshold (e.g. 0.05 to 0.95)
        (bit_ratio * 10.0).clamp(0.05, 0.95)
    });

    #[allow(clippy::cast_precision_loss)]
    let primary_cov = (covered_bp as f64) / (seq_x.len().max(1) as f64);

    let filter_uncovered = |hits: &mut Vec<(usize, usize)>| {
        if covered_bp > 0 {
            hits.retain(|&(px, _)| !covered_x.get(px).copied().unwrap_or(false));
        }
    };

    // Sparse-exon fallback: w=10 exact seeds in intervals not covered by w=12 LCBs.
    let x10_idx = FlatSeedIndex::build_with_sdust(
        seq_x,
        MASK_W10,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
        Some(&idx_x.sdust),
    )?;
    let y10_idx = FlatSeedIndex::build_with_sdust(
        seq_y,
        MASK_W10,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
        Some(&idx_y.sdust),
    )?;
    let (mut hits10_fwd, mut hits10_rc) = x10_idx.intersect(&y10_idx);
    let t_w10 = t_start.elapsed();
    filter_uncovered(&mut hits10_fwd);
    filter_uncovered(&mut hits10_rc);
    let mut combined_fwd = raw_fwd.clone();
    let mut combined_rc = raw_rc.clone();
    append_nonoverlapping_lcbs(
        &mut combined_fwd,
        &raw_fwd,
        &hits10_fwd,
        cfg.fallback_min_anchors,
        cfg.fallback_min_span,
    );
    append_nonoverlapping_lcbs(
        &mut combined_rc,
        &raw_rc,
        &hits10_rc,
        cfg.fallback_min_anchors,
        cfg.fallback_min_span,
    );

    // LASTZ-style 12-of-19 + 1-transition hits gated by ungapped X-drop HSP score
    // (skip when w=12 already covers >= 85% of the target contig).
    let y_lz12 = if primary_cov < min_fallback_cov {
        let idx = FlatSeedIndex::build_with_sdust(
            seq_y,
            MASK_LASTZ12,
            cfg.fallback_max_occ,
            cfg.skip_softmasked,
            false,
            Some(&idx_y.sdust),
        )?;
        let (mut hits_lz_fwd, mut hits_lz_rc) = find_lastz12_1transition_hits_flat(
            seq_x,
            Some(&idx_x.sdust),
            &idx,
            cfg.skip_softmasked,
        );
        filter_uncovered(&mut hits_lz_fwd);
        filter_uncovered(&mut hits_lz_rc);
        append_hsp_gated_lastz12_lcbs(
            &mut combined_fwd,
            seq_x,
            seq_y,
            &hits_lz_fwd,
            base_hsp_thresh + hsp_bonus,
        );
        append_hsp_gated_lastz12_lcbs(
            &mut combined_rc,
            seq_x,
            seq_y_rc,
            &hits_lz_rc,
            base_hsp_thresh + hsp_bonus,
        );
        Some(idx)
    } else {
        None
    };
    let t_lz1 = t_start.elapsed();

    if cfg.use_multi_shape_seeds && primary_cov < min_fallback_cov {
        // Skip O(N * 264) 2-transition scan once search space exceeds ~10^9 cell lookups
        // where single-transition + AA5 + BOCPD/MDL already dominate runtime.
        if (seq_x.len() as u128) * (seq_y.len() as u128) < 1_000_000_000
            && let Some(lz_idx) = &y_lz12
        {
            let (mut double_trans_fwd, mut double_trans_rc) = find_lastz12_2transition_hits_flat(
                seq_x,
                Some(&idx_x.sdust),
                lz_idx,
                cfg.skip_softmasked,
            );
            filter_uncovered(&mut double_trans_fwd);
            filter_uncovered(&mut double_trans_rc);
            append_hsp_gated_lastz12_lcbs(
                &mut combined_fwd,
                seq_x,
                seq_y,
                &double_trans_fwd,
                base_hsp_thresh + 30 + hsp_bonus,
            );
            append_hsp_gated_lastz12_lcbs(
                &mut combined_rc,
                seq_x,
                seq_y_rc,
                &double_trans_rc,
                base_hsp_thresh + 30 + hsp_bonus,
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
        let (mut hits_11_fwd, mut hits_11_rc) = find_lastz11_1transition_hits_flat(
            seq_x,
            Some(&idx_x.sdust),
            &y_lz11,
            cfg.skip_softmasked,
        );
        filter_uncovered(&mut hits_11_fwd);
        filter_uncovered(&mut hits_11_rc);
        append_hsp_gated_lastz12_lcbs(
            &mut combined_fwd,
            seq_x,
            seq_y,
            &hits_11_fwd,
            base_hsp_thresh + 30 + hsp_bonus,
        );
        append_hsp_gated_lastz12_lcbs(
            &mut combined_rc,
            seq_x,
            seq_y_rc,
            &hits_11_rc,
            base_hsp_thresh + 30 + hsp_bonus,
        );

        // PatternHunter II weight-11 spaced seed with 2-hit diagonal gating (Ma et al. 2002).
        let x_ph = FlatSeedIndex::build_with_sdust(
            seq_x,
            MASK_PH11_18,
            cfg.fallback_max_occ,
            cfg.skip_softmasked,
            false,
            Some(&idx_x.sdust),
        )?;
        let y_ph = FlatSeedIndex::build_with_sdust(
            seq_y,
            MASK_PH11_18,
            cfg.fallback_max_occ,
            cfg.skip_softmasked,
            false,
            Some(&idx_y.sdust),
        )?;
        let (mut ph_fwd, mut ph_rc) = x_ph.intersect(&y_ph);
        filter_uncovered(&mut ph_fwd);
        filter_uncovered(&mut ph_rc);
        append_two_hit_diagonal_lcbs(&mut combined_fwd, &ph_fwd, 64, 4);
        append_two_hit_diagonal_lcbs(&mut combined_rc, &ph_rc, 64, 4);

        // Apply R/Y (1-bit) and w=10 1-transition only while random-collision density
        // stays tractable (< 10^9 search space cells).
        let mut cached_ry: Option<(FlatSeedIndex, FlatSeedIndex)> = None;
        if (seq_x.len() as u128) * (seq_y.len() as u128) < 1_000_000_000 {
            // HSP-gated 1-bit purine/pyrimidine (R/Y) seeding for multi-transition divergent homologs.
            let x_ry = FlatSeedIndex::build_with_sdust(
                seq_x,
                MASK_RY20,
                2,
                cfg.skip_softmasked,
                true,
                Some(&idx_x.sdust),
            )?;
            let y_ry = FlatSeedIndex::build_with_sdust(
                seq_y,
                MASK_RY20,
                2,
                cfg.skip_softmasked,
                true,
                Some(&idx_y.sdust),
            )?;
            let (mut hits_ry_fwd, mut hits_ry_rc) = x_ry.intersect(&y_ry);
            if cfg.use_ry_seeds {
                cached_ry = Some((x_ry, y_ry));
            }
            filter_uncovered(&mut hits_ry_fwd);
            filter_uncovered(&mut hits_ry_rc);
            append_hsp_gated_lastz12_lcbs(
                &mut combined_fwd,
                seq_x,
                seq_y,
                &hits_ry_fwd,
                base_hsp_thresh + 40 + hsp_bonus,
            );
            append_hsp_gated_lastz12_lcbs(
                &mut combined_rc,
                seq_x,
                seq_y_rc,
                &hits_ry_rc,
                base_hsp_thresh + 40 + hsp_bonus,
            );

            // HSP-gated w=10 (10-of-14) 1-transition seeding.
            let (mut hits_w10_fwd, mut hits_w10_rc) =
                find_w10_1transition_hits_flat(seq_x, &y10_idx, cfg.skip_softmasked);
            filter_uncovered(&mut hits_w10_fwd);
            filter_uncovered(&mut hits_w10_rc);
            append_hsp_gated_lastz12_lcbs(
                &mut combined_fwd,
                seq_x,
                seq_y,
                &hits_w10_fwd,
                base_hsp_thresh + 40 + hsp_bonus,
            );
            append_hsp_gated_lastz12_lcbs(
                &mut combined_rc,
                seq_x,
                seq_y_rc,
                &hits_w10_rc,
                base_hsp_thresh + 40 + hsp_bonus,
            );
        }

        // Amino-acid HSP-gated 6-frame translated spaced 5-of-6 seeding for divergent coding homologs.
        let mut aa_fwd = find_translated_aa5_hits(seq_x, seq_y, 4, cfg.skip_softmasked);
        let mut aa_rc = find_translated_aa5_hits(seq_x, seq_y_rc, 4, cfg.skip_softmasked);
        filter_uncovered(&mut aa_fwd);
        filter_uncovered(&mut aa_rc);
        append_hsp_gated_aa_lcbs(
            &mut combined_fwd,
            seq_x,
            seq_y,
            &aa_fwd,
            (base_hsp_thresh / 2) + aa_bonus,
        );
        append_hsp_gated_aa_lcbs(
            &mut combined_rc,
            seq_x,
            seq_y_rc,
            &aa_rc,
            (base_hsp_thresh / 2) + aa_bonus,
        );

        // Transition-tolerant R/Y (1-bit purine/pyrimidine) fallback for divergent homologs.
        if cfg.use_ry_seeds {
            let (x_ry, y_ry) = match cached_ry {
                Some(pair) => pair,
                None => (
                    FlatSeedIndex::build_with_sdust(
                        seq_x,
                        MASK_RY20,
                        2,
                        cfg.skip_softmasked,
                        true,
                        Some(&idx_x.sdust),
                    )?,
                    FlatSeedIndex::build_with_sdust(
                        seq_y,
                        MASK_RY20,
                        2,
                        cfg.skip_softmasked,
                        true,
                        Some(&idx_y.sdust),
                    )?,
                ),
            };
            let (hits_ry_fwd, hits_ry_rc) = x_ry.intersect(&y_ry);
            let prev_fwd = combined_fwd.clone();
            let prev_rc = combined_rc.clone();
            append_nonoverlapping_lcbs(
                &mut combined_fwd,
                &prev_fwd,
                &hits_ry_fwd,
                cfg.fallback_min_anchors,
                cfg.fallback_min_span,
            );
            append_nonoverlapping_lcbs(
                &mut combined_rc,
                &prev_rc,
                &hits_ry_rc,
                cfg.fallback_min_anchors,
                cfg.fallback_min_span,
            );
        }
    } else if cfg.use_ry_seeds {
        let x_ry = FlatSeedIndex::build_with_sdust(
            seq_x,
            MASK_RY20,
            2,
            cfg.skip_softmasked,
            true,
            Some(&idx_x.sdust),
        )?;
        let y_ry = FlatSeedIndex::build_with_sdust(
            seq_y,
            MASK_RY20,
            2,
            cfg.skip_softmasked,
            true,
            Some(&idx_y.sdust),
        )?;
        let (hits_ry_fwd, hits_ry_rc) = x_ry.intersect(&y_ry);
        let prev_fwd = combined_fwd.clone();
        let prev_rc = combined_rc.clone();
        append_nonoverlapping_lcbs(
            &mut combined_fwd,
            &prev_fwd,
            &hits_ry_fwd,
            cfg.fallback_min_anchors,
            cfg.fallback_min_span,
        );
        append_nonoverlapping_lcbs(
            &mut combined_rc,
            &prev_rc,
            &hits_ry_rc,
            cfg.fallback_min_anchors,
            cfg.fallback_min_span,
        );
    }

    let lcbs_fwd = merge_colinear_lcbs(&combined_fwd, max_local_rescue_bp);
    let lcbs_rc = merge_colinear_lcbs(&combined_rc, max_local_rescue_bp);
    let t_seed_end = t_start.elapsed();

    let mut pairs_vec: Vec<(usize, usize)> = Vec::new();
    let mut dp_jobs: Vec<(bool, usize, usize, usize, usize)> = Vec::new();
    let mut xdrop_jobs: Vec<(bool, usize, usize, usize, usize, bool)> = Vec::new();
    let mut macro_bubbles = 0usize;

    for (is_rc, lcb_list, y_work) in [
        (false, &lcbs_fwd, seq_y.as_slice()),
        (true, &lcbs_rc, seq_y_rc.as_slice()),
    ] {
        // Inter-LCB sparse-region rescue using w=10 spaced seeds (with 2-hit w=9 fallback).
        let mut rescued_inter = Vec::new();
        if lcb_list.is_empty() && seq_x.len() <= 2500 && y_work.len() <= 2500 {
            let sub10 = local_w9_subchain(seq_x, y_work, 0, 0, MASK_W10, 200, 1800)?;
            if sub10.len() >= 3 {
                rescued_inter.push(sub10);
            } else {
                let sub9 = local_w9_subchain(seq_x, y_work, 0, 0, MASK_W9, 120, 1200)?;
                if sub9.len() >= 4 {
                    rescued_inter.push(sub9);
                }
            }
        }
        if let Some(first_lcb) = lcb_list.first()
            && seq_x.len() <= 50_000
        {
            let (fx, fy) = first_lcb[0];
            if (100..=2500).contains(&fx) && (100..=2500).contains(&fy) {
                let sub =
                    local_w9_subchain(&seq_x[..fx], &y_work[..fy], 0, 0, MASK_W10, 200, 1800)?;
                if sub.len() >= 3 {
                    rescued_inter.push(sub);
                }
            }
        }
        if let Some(last_lcb) = lcb_list.last()
            && seq_x.len() <= 50_000
        {
            let lx = last_lcb[last_lcb.len() - 1].0 + span;
            let ly = last_lcb[last_lcb.len() - 1].1 + span;
            if seq_x.len() > lx
                && y_work.len() > ly
                && (100..=2500).contains(&(seq_x.len() - lx))
                && (100..=2500).contains(&(y_work.len() - ly))
            {
                let sub =
                    local_w9_subchain(&seq_x[lx..], &y_work[ly..], lx, ly, MASK_W10, 200, 1800)?;
                if sub.len() >= 3 {
                    rescued_inter.push(sub);
                }
            }
        }
        for w in lcb_list.windows(2) {
            let x0 = w[0][w[0].len() - 1].0 + span;
            let y0 = w[0][w[0].len() - 1].1 + span;
            let x1 = w[1][0].0;
            let y1 = w[1][0].1;
            if x1 > x0
                && y1 > y0
                && (100..=25_000).contains(&(x1 - x0))
                && (100..=25_000).contains(&(y1 - y0))
            {
                let sub10 = local_w9_subchain(
                    &seq_x[x0..x1],
                    &y_work[y0..y1],
                    x0,
                    y0,
                    MASK_W10,
                    200,
                    1800,
                )?;
                let min_sub10 = if x1 - x0 <= 2500 && y1 - y0 <= 2500 {
                    2
                } else {
                    3
                };
                if sub10.len() >= min_sub10 {
                    rescued_inter.push(sub10);
                } else {
                    let sub9 = local_w9_subchain(
                        &seq_x[x0..x1],
                        &y_work[y0..y1],
                        x0,
                        y0,
                        MASK_W9,
                        120,
                        1200,
                    )?;
                    if sub9.len() >= 4 {
                        rescued_inter.push(sub9);
                    }
                }
            }
        }

        for lcb in lcb_list.iter().chain(&rescued_inter) {
            // General MDL gate: any LCB must earn at least `effective_min_bits` from its
            // anchor base matches minus mismatch penalty (with no size-specific magic branch).
            if lcb.len() < 8 {
                let mut matches = 0usize;
                let mut mismatches = 0usize;
                let mut last_end_x = 0usize;
                for &(ax, ay) in lcb {
                    for off in 0..span {
                        let px = ax + off;
                        let py = ay + off;
                        if px < last_end_x || px >= seq_x.len() || py >= pv_len {
                            continue;
                        }
                        if seq_x[px].eq_ignore_ascii_case(&y_work[py]) {
                            matches += 1;
                        } else {
                            mismatches += 1;
                        }
                    }
                    last_end_x = (ax + span).max(last_end_x);
                }
                #[allow(clippy::cast_precision_loss)]
                let net_bits = 1.6 * (matches as f64) - 2.4 * (mismatches as f64);
                if net_bits < effective_min_bits {
                    continue;
                }
            }
            let mut all_anchors = Vec::with_capacity(lcb.len() * 2);
            all_anchors.push(lcb[0]);
            for w in lcb.windows(2) {
                let x0 = w[0].0 + span;
                let y0 = w[0].1 + span;
                let x1 = w[1].0;
                let y1 = w[1].1;
                if x1 > x0 && y1 > y0 {
                    let lx = x1 - x0;
                    let ly = y1 - y0;
                    if (lx > max_gap_bp || ly > max_gap_bp)
                        && lx <= max_local_rescue_bp
                        && ly <= max_local_rescue_bp
                    {
                        let sub = local_w9_subchain(
                            &seq_x[x0..x1],
                            &y_work[y0..y1],
                            x0,
                            y0,
                            MASK_W9,
                            80,
                            600,
                        )?;
                        all_anchors.extend(sub);
                    }
                }
                all_anchors.push(w[1]);
            }

            for &(ax, ay) in &all_anchors {
                for off in 0..span {
                    if ax + off >= seq_x.len() || ay + off >= pv_len {
                        break;
                    }
                    if !seq_x[ax + off].eq_ignore_ascii_case(&y_work[ay + off]) {
                        continue;
                    }
                    let py = if is_rc {
                        pv_len - 1 - (ay + off)
                    } else {
                        ay + off
                    };
                    pairs_vec.push((ax + off, py));
                }
            }

            let (fx0, fy0) = all_anchors[0];
            if fx0 >= 12 && fy0 >= 12 {
                xdrop_jobs.push((
                    is_rc,
                    fx0.saturating_sub(600),
                    fx0,
                    fy0.saturating_sub(600),
                    fy0,
                    true,
                ));
            }
            let ex0 = all_anchors[all_anchors.len() - 1].0 + span;
            let ey0 = all_anchors[all_anchors.len() - 1].1 + span;
            if seq_x.len().saturating_sub(ex0) >= 12 && pv_len.saturating_sub(ey0) >= 12 {
                xdrop_jobs.push((
                    is_rc,
                    ex0,
                    (ex0 + 600).min(seq_x.len()),
                    ey0,
                    (ey0 + 600).min(pv_len),
                    false,
                ));
            }

            for w in all_anchors.windows(2) {
                let x0 = w[0].0 + span;
                let y0 = w[0].1 + span;
                let x1 = w[1].0;
                let y1 = w[1].1;
                if x1 <= x0 || y1 <= y0 {
                    continue;
                }
                let lx = x1 - x0;
                let ly = y1 - y0;
                if lx > max_gap_bp || ly > max_gap_bp {
                    macro_bubbles += 1;
                    xdrop_jobs.push((is_rc, x0, x1, y0, y1, false));
                    xdrop_jobs.push((is_rc, x0, x1, y0, y1, true));
                    continue;
                }
                dp_jobs.push((is_rc, x0, x1, y0, y1));
            }
        }
    }

    // Estimate divergence-adaptive Pair-HMM parameters empirically from LCB anchors (with Beta prior shrinkage).
    let hmm_params = estimate_hmm_params(seq_x, seq_y, seq_y_rc, &lcbs_fwd, &lcbs_rc, span);

    // Run banded MDL DP gaps in parallel via rayon.
    let t_dp_start = std::time::Instant::now();
    let dp_results: Vec<(bool, usize, usize, GapResult)> = dp_jobs
        .par_iter()
        .map(|&(is_rc, x0, x1, y0, y1)| {
            let y_work = if is_rc { seq_y_rc } else { seq_y };
            let res = banded_mdl_gap_with_min_post(
                &seq_x[x0..x1],
                &y_work[y0..y1],
                cfg.band_w,
                effective_min_bits,
                &hmm_params,
                cfg.min_post_prob,
            );
            (is_rc, x0, y0, res)
        })
        .collect();

    let mut internal_bubbles = 0usize;
    for (is_rc, x0, y0, res) in dp_results {
        internal_bubbles += res.bubbles;
        for (dx, dy) in res.kept_pairs {
            let py = if is_rc {
                pv_len - 1 - (y0 + dy)
            } else {
                y0 + dy
            };
            pairs_vec.push((x0 + dx, py));
        }
    }

    let t_dp = t_dp_start.elapsed();
    let t_xd_start = std::time::Instant::now();
    let xdrop_results: Vec<(bool, Vec<(usize, usize)>)> = xdrop_jobs
        .par_iter()
        .map(|&(is_rc, x0, x1, y0, y1, rev)| {
            let y_work = if is_rc { seq_y_rc } else { seq_y };
            let pairs = mdl_xdrop_extend_with_params(
                &seq_x[x0..x1],
                &y_work[y0..y1],
                x0,
                y0,
                rev,
                cfg.max_ext,
                40,
                &hmm_params,
            );
            (is_rc, pairs)
        })
        .collect();

    for (is_rc, pairs) in xdrop_results {
        for (px, py) in pairs {
            let y_orig = if is_rc { pv_len - 1 - py } else { py };
            pairs_vec.push((px, y_orig));
        }
    }

    let t_xd = t_xd_start.elapsed();
    if profile {
        eprintln!(
            "[profile] w12={:.2}s (fwd={},rc={},cov={:.1}%) w10={:.2}s lz1={:.2}s seed_end={:.2}s (lcbs={}+{}) dp={:.2}s (jobs={}) xdrop={:.2}s (jobs={})",
            t_w12.as_secs_f64(),
            raw_fwd.len(),
            raw_rc.len(),
            primary_cov * 100.0,
            t_w10.saturating_sub(t_w12).as_secs_f64(),
            t_lz1.saturating_sub(t_w10).as_secs_f64(),
            t_seed_end.saturating_sub(t_lz1).as_secs_f64(),
            lcbs_fwd.len(),
            lcbs_rc.len(),
            t_dp.as_secs_f64(),
            dp_jobs.len(),
            t_xd.as_secs_f64(),
            xdrop_jobs.len(),
        );
    }

    pairs_vec.sort_unstable();
    pairs_vec.dedup();
    let aligned_pairs = pairs_vec;

    Ok(ContigAlignment {
        aligned_pairs,
        lcbs_fwd: lcbs_fwd.len(),
        lcbs_rc: lcbs_rc.len(),
        macro_bubbles,
        internal_bubbles,
    })
}

/// Estimate divergence-adaptive Pair-HMM parameters from LCB anchors with Bayesian shrinkage.
fn estimate_hmm_params(
    seq_x: &[u8],
    seq_y: &[u8],
    seq_y_rc: &[u8],
    lcbs_fwd: &[Vec<(usize, usize)>],
    lcbs_rc: &[Vec<(usize, usize)>],
    span: usize,
) -> HmmParams {
    let mut matches = 0usize;
    let mut total_bases = 0usize;
    let mut indel_events = 0usize;
    let mut total_indel_len = 0usize;
    let mut total_span = 0usize;

    for (y_work, lcbs) in [(seq_y, lcbs_fwd), (seq_y_rc, lcbs_rc)] {
        for lcb in lcbs {
            for &(ax, ay) in lcb.iter().take(64) {
                if ax + span <= seq_x.len() && ay + span <= y_work.len() {
                    for off in 0..span {
                        let bx = seq_x[ax + off].to_ascii_uppercase();
                        let by = y_work[ay + off].to_ascii_uppercase();
                        if matches!(bx, b'A' | b'C' | b'G' | b'T')
                            && matches!(by, b'A' | b'C' | b'G' | b'T')
                        {
                            total_bases += 1;
                            if bx == by {
                                matches += 1;
                            }
                        }
                    }
                }
                if total_bases >= 4000 {
                    break;
                }
            }
            for w in lcb.windows(2).take(64) {
                let dx = w[1].0.saturating_sub(w[0].0);
                let dy = w[1].1.saturating_sub(w[0].1);
                if dx > 0 && dy > 0 && dx <= 500 && dy <= 500 {
                    total_span += dx.max(dy);
                    if dx != dy {
                        indel_events += 1;
                        total_indel_len += dx.abs_diff(dy);
                    }
                }
            }
            if total_bases >= 4000 {
                break;
            }
        }
    }

    if total_bases < 60 {
        return HmmParams::default();
    }

    // LCB anchors are enriched for conserved exons; cap upper identity_prior at 0.72
    // so inter-anchor gaps in divergent mammals still receive adequate homology sensitivity.
    #[allow(clippy::cast_precision_loss)]
    let identity_prior =
        ((matches as f64 + 195.0) / (total_bases as f64 + 300.0)).clamp(0.60, 0.72);
    #[allow(clippy::cast_precision_loss)]
    let gap_open_bits = if total_span >= 200 {
        ((indel_events as f64 + 1.0) / (total_span as f64 + 22.0))
            .log2()
            .clamp(-5.5, -4.0)
    } else {
        HmmParams::default().gap_open_bits
    };
    #[allow(clippy::cast_precision_loss)]
    let gap_extend_bits = if indel_events >= 3 && total_indel_len > 0 {
        let mean_len = (total_indel_len as f64 + 2.5) / (indel_events as f64 + 1.0);
        ((mean_len - 1.0) / mean_len)
            .max(0.01)
            .log2()
            .clamp(-1.8, -0.5)
    } else {
        HmmParams::default().gap_extend_bits
    };

    HmmParams {
        identity_prior,
        gap_open_bits,
        gap_extend_bits,
    }
}
