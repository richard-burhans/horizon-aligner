//! Banded MDL+BOCPD dynamic programming, X-drop extension, and contig-pair alignment.
//! ## Horizon Bubbles
//! In this layer, areas where the alignment lacks sufficient confidence to align (either because a gap is too large to enter DP, or because the MDL DP abstains due to negative bit-savings) are termed **Horizon Bubbles**.
//! - Unaligned inter-LCB gaps larger than `max_gap_bp` yield *macro bubbles*.
//! - DP dropouts within a gap yield *internal bubbles*, returned via `GapResult`.

use crate::AlignError;
use crate::bocpd::MIN_POST_PROB;
use crate::seeding::{
    ContigIndex, MASK_PRIMARY_11_15, MASK_W9, find_w12_raw_lcbs, local_w9_subchain,
    prepare_contig_index,
};
use rayon::prelude::*;

mod dag;
mod dag_seed;
mod dp;
mod fallback;
mod xdrop;
pub use dag::*;
pub use dag_seed::*;
pub use dp::*;
pub use xdrop::*;

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
    /// Optional minimum primary-seed coverage fraction below which fallback seeding triggers.
    pub min_fallback_cov: Option<f64>,
    /// Optional minimum ungapped HSP score threshold for LASTZ-style seed gating.
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

/// Align two DNA contigs using default configuration.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if internal masks fail validation.
pub fn align_contig_pair(seq_x: &[u8], seq_y: &[u8]) -> Result<ContigAlignment, AlignError> {
    let idx_x = prepare_contig_index(seq_x)?;
    let idx_y = prepare_contig_index(seq_y)?;
    align_contig_pair_indexed(&idx_x, &idx_y)
}

/// Align two pre-indexed contigs using default configuration.
///
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

struct PrimaryMetrics {
    covered_x: Vec<bool>,
    covered_bp: usize,
    primary_cov: f64,
    effective_min_bits: f64,
    base_hsp_thresh: i32,
    min_fallback_cov: f64,
    hsp_bonus: i32,
    aa_bonus: i32,
}

fn compute_primary_metrics(
    len_x: usize,
    len_y: usize,
    raw_fwd: &[Vec<(usize, usize)>],
    raw_rc: &[Vec<(usize, usize)>],
    span: usize,
    cfg: &AlignConfig,
) -> PrimaryMetrics {
    #[allow(clippy::cast_precision_loss)]
    let search_space_bits = ((len_x.max(1) as f64) * (len_y.max(1) as f64))
        .log2()
        .max(10.0);
    // Continuous Karlin-Altschul / MDL scaling with log2(N * M), avoiding fixed size cutoffs.
    #[allow(clippy::cast_possible_truncation)]
    let hsp_bonus = ((search_space_bits - 20.0).max(0.0) * 3.8).round() as i32;
    #[allow(clippy::cast_possible_truncation)]
    let aa_bonus = ((search_space_bits - 20.0).max(0.0) * 1.7).round() as i32;
    let effective_min_bits = cfg.min_delta_bits.max(KARLIN_MDL_COEFF * search_space_bits);

    let mut covered_x = vec![false; len_x];
    let mut covered_bp = 0usize;
    for lcb in raw_fwd.iter().chain(raw_rc.iter()) {
        if let (Some(first), Some(last)) = (lcb.first(), lcb.last()) {
            let s = first.0.min(len_x);
            let e = (last.0 + span).min(len_x);
            for flag in &mut covered_x[s..e] {
                if !*flag {
                    *flag = true;
                    covered_bp += 1;
                }
            }
        }
    }
    #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
    let base_hsp_thresh = cfg
        .hsp_threshold
        .unwrap_or_else(|| (100.0 + (effective_min_bits * 1.5)).round() as i32);

    // MDL-derived default fallback coverage threshold:
    // If the expected bit savings of the primary coverage is less than the effective min bits
    // across the sequence, trigger fallback.
    #[allow(clippy::cast_precision_loss)]
    let min_fallback_cov = cfg.min_fallback_cov.unwrap_or_else(|| {
        let expected_bits_per_base = 2.0;
        let bit_ratio = effective_min_bits / (len_x.max(1) as f64 * expected_bits_per_base);
        (bit_ratio * 10.0).clamp(0.05, 0.95)
    });

    #[allow(clippy::cast_precision_loss)]
    let primary_cov = (covered_bp as f64) / (len_x.max(1) as f64);

    PrimaryMetrics {
        covered_x,
        covered_bp,
        primary_cov,
        effective_min_bits,
        base_hsp_thresh,
        min_fallback_cov,
        hsp_bonus,
        aa_bonus,
    }
}

/// Align two pre-indexed contigs using an explicit [`AlignConfig`].
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if internal masks fail validation.
pub fn align_contig_pair_with_config(
    idx_x: &ContigIndex,
    idx_y: &ContigIndex,
    cfg: &AlignConfig,
) -> Result<ContigAlignment, AlignError> {
    let span = MASK_PRIMARY_11_15.len();
    let seq_x = &idx_x.seq;
    let seq_y = &idx_y.seq;
    let seq_y_rc = &idx_y.seq_rc;

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

    let PrimaryMetrics {
        covered_x,
        covered_bp,
        primary_cov,
        effective_min_bits,
        base_hsp_thresh,
        min_fallback_cov,
        hsp_bonus,
        aa_bonus,
    } = compute_primary_metrics(seq_x.len(), seq_y.len(), &raw_fwd, &raw_rc, span, cfg);

    let fallback::FallbackOutput {
        lcbs_fwd,
        lcbs_rc,
        mut spans_fwd,
        mut spans_rc,
        t_w10,
        t_lz1,
        t_seed_end,
    } = fallback::collect_fallback_lcbs(
        idx_x,
        idx_y,
        &raw_fwd,
        &raw_rc,
        &covered_x,
        covered_bp,
        primary_cov,
        min_fallback_cov,
        base_hsp_thresh,
        hsp_bonus,
        aa_bonus,
        cfg,
        t_start,
    )?;

    let mut ext = schedule_and_run_extension(
        seq_x,
        seq_y,
        seq_y_rc,
        &lcbs_fwd,
        &lcbs_rc,
        &mut spans_fwd,
        &mut spans_rc,
        span,
        MAX_GAP_BP,
        MAX_LOCAL_RESCUE_BP,
        effective_min_bits,
        cfg,
    )?;

    if profile {
        eprintln!(
            "[profile] w12={:.2}s (fwd={},rc={},cov={:.1}%) w10={:.2}s lz1={:.2}s seed_end={:.2}s (lcbs={}+{}) dp+xd={:.2}s (dp={},xd={})",
            t_w12.as_secs_f64(),
            raw_fwd.len(),
            raw_rc.len(),
            primary_cov * 100.0,
            t_w10.saturating_sub(t_w12).as_secs_f64(),
            t_lz1.saturating_sub(t_w10).as_secs_f64(),
            t_seed_end.saturating_sub(t_lz1).as_secs_f64(),
            lcbs_fwd.len(),
            lcbs_rc.len(),
            ext.t_dp.as_secs_f64(),
            ext.n_dp_jobs,
            ext.n_xdrop_jobs,
        );
    }

    ext.pairs_vec.sort_unstable();
    ext.pairs_vec.dedup();

    Ok(ContigAlignment {
        aligned_pairs: ext.pairs_vec,
        lcbs_fwd: lcbs_fwd.len(),
        lcbs_rc: lcbs_rc.len(),
        macro_bubbles: ext.macro_bubbles,
        internal_bubbles: ext.internal_bubbles,
    })
}

const MAX_GAP_BP: usize = 2400;
pub(super) const MAX_LOCAL_RESCUE_BP: usize = 40_000;
const XDROP_BAND_W: usize = 40;

struct ExtensionOutput {
    pairs_vec: Vec<(usize, usize)>,
    macro_bubbles: usize,
    internal_bubbles: usize,
    t_dp: std::time::Duration,
    n_dp_jobs: usize,
    n_xdrop_jobs: usize,
}

#[allow(clippy::too_many_arguments)]
fn schedule_and_run_extension(
    seq_x: &[u8],
    seq_y: &[u8],
    seq_y_rc: &[u8],
    lcbs_fwd: &[Vec<(usize, usize)>],
    lcbs_rc: &[Vec<(usize, usize)>],
    spans_fwd: &mut std::collections::HashMap<(usize, usize), usize>,
    spans_rc: &mut std::collections::HashMap<(usize, usize), usize>,
    span: usize,
    max_gap_bp: usize,
    max_local_rescue_bp: usize,
    effective_min_bits: f64,
    cfg: &AlignConfig,
) -> Result<ExtensionOutput, AlignError> {
    let mut pairs_vec: Vec<(usize, usize)> = Vec::new();
    let mut dp_jobs: Vec<(bool, usize, usize, usize, usize)> = Vec::new();
    let mut xdrop_jobs: Vec<(bool, usize, usize, usize, usize, bool)> = Vec::new();
    let mut macro_bubbles = 0usize;

    for (is_rc, lcb_list, y_work, anchor_spans) in [
        (false, lcbs_fwd, seq_y, spans_fwd),
        (true, lcbs_rc, seq_y_rc, spans_rc),
    ] {
        let rescued_inter =
            fallback::rescue_inter_lcb_gaps(seq_x, y_work, lcb_list, anchor_spans, span)?;
        macro_bubbles += schedule_strand_lcb_jobs(
            is_rc,
            seq_x,
            y_work,
            lcb_list,
            &rescued_inter,
            anchor_spans,
            span,
            max_gap_bp,
            max_local_rescue_bp,
            cfg.max_ext,
            effective_min_bits,
            &mut pairs_vec,
            &mut dp_jobs,
            &mut xdrop_jobs,
        )?;
    }

    // Estimate divergence-adaptive Pair-HMM parameters empirically from LCB anchors (with Beta prior shrinkage).
    let hmm_params = estimate_hmm_params(seq_x, seq_y, seq_y_rc, lcbs_fwd, lcbs_rc, span);
    // Calibrate finite-sample Conformal Risk Control floor (n_cal = 15, alpha = 0.07) at observed GC.
    let mut gc_cnt = 0usize;
    let mut acgt_cnt = 0usize;
    for &b in seq_x.iter().take(4096).chain(seq_y.iter().take(4096)) {
        match b.to_ascii_uppercase() {
            b'G' | b'C' => {
                gc_cnt += 1;
                acgt_cnt += 1;
            }
            b'A' | b'T' => acgt_cnt += 1,
            _ => {}
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let observed_gc = if acgt_cnt > 0 {
        gc_cnt as f64 / acgt_cnt as f64
    } else {
        0.5
    };
    let crc_min_bits = dp::calibrate_null_crc_for_gc(observed_gc, 15, 0.07, effective_min_bits, &hmm_params);

    let t_dp_start = std::time::Instant::now();
    let internal_bubbles = execute_dp_and_xdrop(
        seq_x,
        seq_y,
        seq_y_rc,
        crc_min_bits,
        &hmm_params,
        cfg,
        &dp_jobs,
        &xdrop_jobs,
        &mut pairs_vec,
    );
    let t_dp = t_dp_start.elapsed();
    Ok(ExtensionOutput {
        pairs_vec,
        macro_bubbles,
        internal_bubbles,
        t_dp,
        n_dp_jobs: dp_jobs.len(),
        n_xdrop_jobs: xdrop_jobs.len(),
    })
}

#[allow(clippy::too_many_arguments)]
fn execute_dp_and_xdrop(
    seq_x: &[u8],
    seq_y: &[u8],
    seq_y_rc: &[u8],
    effective_min_bits: f64,
    hmm_params: &HmmParams,
    cfg: &AlignConfig,
    dp_jobs: &[(bool, usize, usize, usize, usize)],
    xdrop_jobs: &[(bool, usize, usize, usize, usize, bool)],
    pairs_vec: &mut Vec<(usize, usize)>,
) -> usize {
    let pv_len = seq_y.len();
    let dp_results: Vec<(bool, usize, usize, GapResult)> = dp_jobs
        .par_iter()
        .map_init(
            dp::DpScratch::default,
            |scratch, &(is_rc, x0, x1, y0, y1)| {
                let y_work = if is_rc { seq_y_rc } else { seq_y };
                let res = dp::banded_mdl_gap_with_scratch(
                    &seq_x[x0..x1],
                    &y_work[y0..y1],
                    cfg.band_w,
                    effective_min_bits,
                    hmm_params,
                    cfg.min_post_prob,
                    scratch,
                );
                (is_rc, x0, y0, res)
            },
        )
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

    let xdrop_results: Vec<(bool, Vec<(usize, usize)>)> = xdrop_jobs
        .par_iter()
        .map_init(
            xdrop::XDropScratch::default,
            |scratch, &(is_rc, x0, x1, y0, y1, rev)| {
                let y_work = if is_rc { seq_y_rc } else { seq_y };
                let pairs = xdrop::mdl_xdrop_extend_with_scratch(
                    &seq_x[x0..x1],
                    &y_work[y0..y1],
                    x0,
                    y0,
                    rev,
                    cfg.max_ext,
                    XDROP_BAND_W,
                    hmm_params,
                    scratch,
                );
                (is_rc, pairs)
            },
        )
        .collect();

    for (is_rc, pairs) in xdrop_results {
        for (px, py) in pairs {
            let y_orig = if is_rc { pv_len - 1 - py } else { py };
            pairs_vec.push((px, y_orig));
        }
    }
    internal_bubbles
}

#[allow(clippy::too_many_arguments)]
fn schedule_strand_lcb_jobs(
    is_rc: bool,
    seq_x: &[u8],
    y_work: &[u8],
    lcb_list: &[Vec<(usize, usize)>],
    rescued_inter: &[Vec<(usize, usize)>],
    anchor_spans: &std::collections::HashMap<(usize, usize), usize>,
    default_span: usize,
    max_gap_bp: usize,
    max_local_rescue_bp: usize,
    max_ext: usize,
    effective_min_bits: f64,
    pairs_vec: &mut Vec<(usize, usize)>,
    dp_jobs: &mut Vec<(bool, usize, usize, usize, usize)>,
    xdrop_jobs: &mut Vec<(bool, usize, usize, usize, usize, bool)>,
) -> Result<usize, AlignError> {
    let pv_len = y_work.len();
    let span_for = |pt: (usize, usize)| anchor_spans.get(&pt).copied().unwrap_or(default_span);
    let mut macro_bubbles = 0usize;

    for lcb in lcb_list.iter().chain(rescued_inter) {
        if !lcb_passes_mdl_gate(lcb, seq_x, y_work, pv_len, &span_for, effective_min_bits) {
            continue;
        }
        let Some(&first_pt) = lcb.first() else {
            continue;
        };
        let mut all_anchors = Vec::with_capacity(lcb.len() * 2);
        all_anchors.push(first_pt);
        for w in lcb.windows(2) {
            let s0 = span_for(w[0]);
            let x0 = w[0].0 + s0;
            let y0 = w[0].1 + s0;
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
            for off in 0..span_for((ax, ay)) {
                if ax + off >= seq_x.len() || ay + off >= pv_len {
                    break;
                }
                if seq_x[ax + off].eq_ignore_ascii_case(&y_work[ay + off]) {
                    let py = if is_rc {
                        pv_len - 1 - (ay + off)
                    } else {
                        ay + off
                    };
                    pairs_vec.push((ax + off, py));
                }
            }
        }

        let (fx0, fy0) = all_anchors[0];
        if fx0 >= 12 && fy0 >= 12 {
            xdrop_jobs.push((
                is_rc,
                fx0.saturating_sub(max_ext),
                fx0,
                fy0.saturating_sub(max_ext),
                fy0,
                true,
            ));
        }
        let last_a = all_anchors[all_anchors.len() - 1];
        let last_span = span_for(last_a);
        let ex0 = last_a.0 + last_span;
        let ey0 = last_a.1 + last_span;
        if seq_x.len().saturating_sub(ex0) >= 12 && pv_len.saturating_sub(ey0) >= 12 {
            xdrop_jobs.push((
                is_rc,
                ex0,
                ex0.saturating_add(max_ext).min(seq_x.len()),
                ey0,
                ey0.saturating_add(max_ext).min(pv_len),
                false,
            ));
        }

        for w in all_anchors.windows(2) {
            let s0 = span_for(w[0]);
            let x0 = w[0].0 + s0;
            let y0 = w[0].1 + s0;
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
    Ok(macro_bubbles)
}

fn lcb_passes_mdl_gate<F: Fn((usize, usize)) -> usize>(
    lcb: &[(usize, usize)],
    seq_x: &[u8],
    y_work: &[u8],
    pv_len: usize,
    span_for: &F,
    effective_min_bits: f64,
) -> bool {
    let mut gc_cnt_x = 0usize;
    let mut gc_cnt_y = 0usize;
    let mut tot_bases = 0usize;
    let mut last_end_x = 0usize;
    for &(ax, ay) in lcb {
        let a_span = span_for((ax, ay));
        for off in 0..a_span {
            let px = ax + off;
            let py = ay + off;
            if px < last_end_x || px >= seq_x.len() || py >= pv_len {
                continue;
            }
            tot_bases += 1;
            if matches!(seq_x[px], b'G' | b'C' | b'g' | b'c') {
                gc_cnt_x += 1;
            }
            if matches!(y_work[py], b'G' | b'C' | b'g' | b'c') {
                gc_cnt_y += 1;
            }
        }
        last_end_x = (ax + a_span).max(last_end_x);
    }
    #[allow(clippy::cast_precision_loss)]
    let x_gc = if tot_bases > 0 {
        (gc_cnt_x as f64 / tot_bases as f64).clamp(0.05, 0.95)
    } else {
        0.5
    };
    #[allow(clippy::cast_precision_loss)]
    let y_gc = if tot_bases > 0 {
        (gc_cnt_y as f64 / tot_bases as f64).clamp(0.05, 0.95)
    } else {
        0.5
    };
    let pi_x = crate::bocpd::base_probs(x_gc);
    let pi_y = crate::bocpd::base_probs(y_gc);
    let mut net_bits = 0.0_f64;
    let mut n_gaps = 0usize;
    let mut prev_end_x = 0usize;
    let mut prev_end_y = 0usize;
    last_end_x = 0;
    for &(ax, ay) in lcb {
        let a_span = span_for((ax, ay));
        if last_end_x > 0 && (ax > prev_end_x || ay > prev_end_y) {
            n_gaps += 1;
        }
        for off in 0..a_span {
            let px = ax + off;
            let py = ay + off;
            if px < last_end_x || px >= seq_x.len() || py >= pv_len {
                continue;
            }
            net_bits += crate::bocpd::joint_pair_bits(
                seq_x[px],
                y_work[py],
                pi_x,
                pi_y,
                crate::bocpd::IDENTITY_PRIOR,
            );
        }
        last_end_x = (ax + a_span).max(last_end_x);
        prev_end_x = ax + a_span;
        prev_end_y = ay + a_span;
    }
    #[allow(clippy::cast_precision_loss)]
    let gap_penalty = (n_gaps as f64) * crate::bocpd::GAP_OPEN_BITS.abs();
    net_bits - gap_penalty >= effective_min_bits
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

    'outer: for (y_work, lcbs) in [(seq_y, lcbs_fwd), (seq_y_rc, lcbs_rc)] {
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
                break 'outer;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mdl_gate_and_hmm_estimation_smoke() {
        let x = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let y = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let lcb: Vec<(usize, usize)> = (0..10).map(|i| (i * 4, i * 4)).collect();
        assert!(lcb_passes_mdl_gate(&lcb, x, y, y.len(), &|_| 4, 8.0));
        assert!(!lcb_passes_mdl_gate(
            &lcb,
            x,
            b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            40,
            &|_| 4,
            8.0
        ));

        let params = estimate_hmm_params(x, y, y, &[lcb], &[], 4);
        assert!((0.60..=0.72).contains(&params.identity_prior));
    }
}
