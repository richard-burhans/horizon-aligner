//! Banded MDL+BOCPD dynamic programming, X-drop extension, and contig-pair alignment.

use crate::AlignError;
use crate::bocpd::{
    GAP_EXTEND_BITS, GAP_OPEN_BITS, IDENTITY_PRIOR, MIN_POST_PROB, NEG_INF, base_probs,
    bocpd_gc_profile, joint_pair_bits, logsumexp2, logsumexp3,
};
use crate::seeding::{
    ContigIndex, MASK_KP14, MASK_LASTZ11, MASK_LASTZ12, MASK_RY20, MASK_W9, MASK_W10,
    append_hsp_gated_aa_lcbs, append_hsp_gated_lastz12_lcbs, append_nonoverlapping_lcbs,
    canonical_spaced_seeds_with_opts, find_lastz11_1transition_hits, find_lastz12_1transition_hits,
    find_lastz12_2transition_hits, find_translated_aa5_hits, find_w10_1transition_hits,
    find_w12_raw_lcbs, intersect_seed_maps, local_w9_subchain, merge_colinear_lcbs,
    prepare_contig_index,
};
use rayon::prelude::*;
use std::collections::HashSet;

/// Result of aligning one inter-anchor gap.
#[derive(Debug, Clone, PartialEq)]
pub struct GapResult {
    /// Kept base-pair coordinates `(x, y)` relative to the gap start.
    pub kept_pairs: Vec<(usize, usize)>,
    /// Number of abstained Horizon Bubbles in this gap.
    pub bubbles: usize,
    /// Total Viterbi bit-savings score for the gap.
    pub score_bits: f64,
}

struct BandGrid {
    n: usize,
    m: usize,
    band_w: usize,
    stride: usize,
    data: Vec<f64>,
}

impl BandGrid {
    fn new(n: usize, m: usize, band_w: usize, init: f64) -> Self {
        let stride = 2 * band_w + 1;
        Self {
            n,
            m,
            band_w,
            stride,
            data: vec![init; (n + 1) * stride],
        }
    }

    #[inline]
    fn center(&self, i: usize) -> isize {
        #[allow(clippy::cast_possible_wrap)]
        let ci = ((i * self.m) / self.n.max(1)) as isize;
        ci
    }

    #[inline]
    fn slot(&self, i: usize, j: usize) -> Option<usize> {
        #[allow(clippy::cast_possible_wrap)]
        let diff = (j as isize) - self.center(i);
        #[allow(clippy::cast_possible_wrap)]
        let bw = self.band_w as isize;
        if diff.abs() <= bw {
            #[allow(clippy::cast_sign_loss)]
            let col = (diff + bw) as usize;
            Some(i * self.stride + col)
        } else {
            None
        }
    }

    #[inline]
    fn get(&self, i: usize, j: usize) -> f64 {
        self.slot(i, j).map_or(NEG_INF, |idx| self.data[idx])
    }

    #[inline]
    fn set(&mut self, i: usize, j: usize, val: f64) {
        if let Some(idx) = self.slot(i, j) {
            self.data[idx] = val;
        }
    }
}

/// Run banded forward-backward MDL + BOCPD alignment on one gap.
#[must_use]
#[allow(
    clippy::too_many_lines,
    clippy::similar_names,
    clippy::many_single_char_names
)]
pub fn banded_mdl_gap(seq_x: &[u8], seq_y: &[u8], band_w: usize, min_block_bits: f64) -> GapResult {
    let n = seq_x.len();
    let m = seq_y.len();
    if n == 0 || m == 0 {
        return GapResult {
            kept_pairs: Vec::new(),
            bubbles: usize::from(n.max(m) >= 6),
            score_bits: 0.0,
        };
    }
    if n.min(m) <= 4 {
        let k = n.min(m);
        return GapResult {
            kept_pairs: (0..k).map(|i| (i, i)).collect(),
            bubbles: 0,
            score_bits: 0.0,
        };
    }

    let gc_x = bocpd_gc_profile(seq_x, 1.0 / 120.0, 80);
    let gc_y = bocpd_gc_profile(seq_y, 1.0 / 120.0, 80);
    let pi_x: Vec<[f64; 4]> = gc_x.into_iter().map(base_probs).collect();
    let pi_y: Vec<[f64; 4]> = gc_y.into_iter().map(base_probs).collect();

    let ln2 = 2.0_f64.ln();
    let go = GAP_OPEN_BITS * ln2;
    let ge = GAP_EXTEND_BITS * ln2;

    let mut m_dp = BandGrid::new(n, m, band_w, NEG_INF);
    let mut ix_dp = BandGrid::new(n, m, band_w, NEG_INF);
    let mut iy_dp = BandGrid::new(n, m, band_w, NEG_INF);
    let mut fm = BandGrid::new(n, m, band_w, NEG_INF);
    let mut fix = BandGrid::new(n, m, band_w, NEG_INF);
    let mut fiy = BandGrid::new(n, m, band_w, NEG_INF);
    let mut sub_grid = BandGrid::new(n, m, band_w, 0.0);

    m_dp.set(0, 0, 0.0);
    fm.set(0, 0, 0.0);

    for i in 0..=n {
        let center = (i * m) / n;
        let j_lo = center.saturating_sub(band_w);
        let j_hi = (center + band_w).min(m);
        for j in j_lo..=j_hi {
            if i == 0 && j == 0 {
                continue;
            }
            if i > 0 {
                let v_ix =
                    (m_dp.get(i - 1, j) + GAP_OPEN_BITS).max(ix_dp.get(i - 1, j) + GAP_EXTEND_BITS);
                ix_dp.set(i, j, v_ix);
                fix.set(
                    i,
                    j,
                    logsumexp2(fm.get(i - 1, j) + go, fix.get(i - 1, j) + ge),
                );
            }
            if j > 0 {
                let v_iy =
                    (m_dp.get(i, j - 1) + GAP_OPEN_BITS).max(iy_dp.get(i, j - 1) + GAP_EXTEND_BITS);
                iy_dp.set(i, j, v_iy);
                fiy.set(
                    i,
                    j,
                    logsumexp2(fm.get(i, j - 1) + go, fiy.get(i, j - 1) + ge),
                );
            }
            if i > 0 && j > 0 {
                let sb = joint_pair_bits(
                    seq_x[i - 1],
                    seq_y[j - 1],
                    pi_x[i - 1],
                    pi_y[j - 1],
                    IDENTITY_PRIOR,
                );
                sub_grid.set(i - 1, j - 1, sb);
                let prev_m = m_dp
                    .get(i - 1, j - 1)
                    .max(ix_dp.get(i - 1, j - 1))
                    .max(iy_dp.get(i - 1, j - 1));
                m_dp.set(i, j, prev_m + sb);
                let prev_fm = logsumexp3(
                    fm.get(i - 1, j - 1),
                    fix.get(i - 1, j - 1),
                    fiy.get(i - 1, j - 1),
                );
                fm.set(i, j, prev_fm + sb * ln2);
            }
        }
    }

    let z_fwd = logsumexp3(fm.get(n, m), fix.get(n, m), fiy.get(n, m));
    let best_end = m_dp.get(n, m).max(ix_dp.get(n, m)).max(iy_dp.get(n, m));
    if z_fwd <= NEG_INF * 0.5 {
        return GapResult {
            kept_pairs: Vec::new(),
            bubbles: 1,
            score_bits: best_end,
        };
    }

    let mut bm = BandGrid::new(n, m, band_w, NEG_INF);
    let mut bix = BandGrid::new(n, m, band_w, NEG_INF);
    let mut biy = BandGrid::new(n, m, band_w, NEG_INF);
    bm.set(n, m, 0.0);
    bix.set(n, m, 0.0);
    biy.set(n, m, 0.0);

    for i in (0..=n).rev() {
        let center = (i * m) / n;
        let j_lo = center.saturating_sub(band_w);
        let j_hi = (center + band_w).min(m);
        for j in (j_lo..=j_hi).rev() {
            if i == n && j == m {
                continue;
            }
            let m_next = if i < n && j < m {
                bm.get(i + 1, j + 1) + sub_grid.get(i, j) * ln2
            } else {
                NEG_INF
            };
            let ix_open = if i < n {
                bix.get(i + 1, j) + go
            } else {
                NEG_INF
            };
            let ix_ext = if i < n {
                bix.get(i + 1, j) + ge
            } else {
                NEG_INF
            };
            let iy_open = if j < m {
                biy.get(i, j + 1) + go
            } else {
                NEG_INF
            };
            let iy_ext = if j < m {
                biy.get(i, j + 1) + ge
            } else {
                NEG_INF
            };
            bm.set(i, j, logsumexp3(m_next, ix_open, iy_open));
            bix.set(i, j, logsumexp2(m_next, ix_ext));
            biy.set(i, j, logsumexp2(m_next, iy_ext));
        }
    }

    let mut raw_cols: Vec<(Option<usize>, Option<usize>, f64, f64)> = Vec::with_capacity(n + m);
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        let vm = m_dp.get(i, j);
        let vix = ix_dp.get(i, j);
        let viy = iy_dp.get(i, j);
        let best = vm.max(vix).max(viy);
        if best <= NEG_INF * 0.5 {
            break;
        }
        if i > 0 && j > 0 && (vm - best).abs() < 1e-7 {
            let sb = sub_grid.get(i - 1, j - 1);
            let post_p = (fm.get(i, j) + bm.get(i, j) - z_fwd).min(0.0).exp();
            raw_cols.push((Some(i - 1), Some(j - 1), sb, post_p));
            i -= 1;
            j -= 1;
        } else if i > 0 && (vix - best).abs() < 1e-7 {
            raw_cols.push((Some(i - 1), None, GAP_EXTEND_BITS, 0.0));
            i -= 1;
        } else if j > 0 {
            raw_cols.push((None, Some(j - 1), GAP_EXTEND_BITS, 0.0));
            j -= 1;
        } else {
            break;
        }
    }
    raw_cols.reverse();

    let k = raw_cols.len();
    let mut abstain = vec![false; k];
    let mut idx = 0;
    while idx < k {
        let (cx, cy, _, pp) = raw_cols[idx];
        if cx.is_some() && cy.is_some() && pp >= MIN_POST_PROB {
            let start = idx;
            let mut block_bits = 0.0;
            while idx < k
                && ((raw_cols[idx].0.is_some()
                    && raw_cols[idx].1.is_some()
                    && raw_cols[idx].3 >= MIN_POST_PROB)
                    || (idx + 1 < k
                        && raw_cols[idx].2 > -2.5
                        && raw_cols[idx + 1].3 >= MIN_POST_PROB))
            {
                block_bits += raw_cols[idx].2;
                idx += 1;
            }
            if block_bits < min_block_bits {
                abstain[start..idx].fill(true);
            }
        } else {
            abstain[idx] = true;
            idx += 1;
        }
    }

    let mut kept_pairs = Vec::with_capacity(k);
    let mut bubbles = 0;
    idx = 0;
    while idx < k {
        if abstain[idx] {
            let start = idx;
            while idx < k && abstain[idx] {
                idx += 1;
            }
            if idx - start >= 6 {
                bubbles += 1;
            }
        } else {
            if let (Some(cx), Some(cy)) = (raw_cols[idx].0, raw_cols[idx].1) {
                kept_pairs.push((cx, cy));
            }
            idx += 1;
        }
    }

    GapResult {
        kept_pairs,
        bubbles,
        score_bits: best_end,
    }
}

/// Result of aligning a full contig pair.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContigAlignment {
    /// Kept 0-based `+` strand `(x, y)` base pairs.
    pub aligned_pairs: Vec<(usize, usize)>,
    /// Number of forward LCBs.
    pub lcbs_fwd: usize,
    /// Number of reverse-complement LCBs.
    pub lcbs_rc: usize,
    /// Number of macro Horizon Bubbles (`> max_gap_bp`).
    pub macro_bubbles: usize,
    /// Number of internal Horizon Bubbles emitted inside DP gaps.
    pub internal_bubbles: usize,
}

/// Extend an alignment outward from an anchor using banded MDL + BOCPD X-drop.
#[must_use]
#[allow(clippy::too_many_lines, clippy::many_single_char_names)]
pub fn mdl_xdrop_extend(
    seq_x: &[u8],
    seq_y: &[u8],
    off_x: usize,
    off_y: usize,
    reverse: bool,
    max_ext: usize,
    band_w: usize,
) -> Vec<(usize, usize)> {
    let x_slice: Vec<u8> = if reverse {
        seq_x.iter().rev().take(max_ext).copied().collect()
    } else {
        seq_x.iter().take(max_ext).copied().collect()
    };
    let y_slice: Vec<u8> = if reverse {
        seq_y.iter().rev().take(max_ext).copied().collect()
    } else {
        seq_y.iter().take(max_ext).copied().collect()
    };
    let n = x_slice.len();
    let m = y_slice.len();
    if n < 12 || m < 12 {
        return Vec::new();
    }

    let gc_x = bocpd_gc_profile(&x_slice, 1.0 / 120.0, 60);
    let gc_y = bocpd_gc_profile(&y_slice, 1.0 / 120.0, 60);
    let pi_x: Vec<[f64; 4]> = gc_x.into_iter().map(base_probs).collect();
    let pi_y: Vec<[f64; 4]> = gc_y.into_iter().map(base_probs).collect();

    let x_drop_bits = 15.0;
    let min_peak_bits = 10.0;
    let stride = 2 * band_w + 1;
    let mut dp = vec![NEG_INF; (n + 1) * stride];
    let mut bt = vec![0_u8; (n + 1) * stride]; // 1=diag, 2=up, 3=left

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
                    IDENTITY_PRIOR,
                );
                let v = dp[ps] + sb;
                if v > cand {
                    cand = v;
                    src = 1;
                }
            }
            if i > 0
                && let Some(ps) = slot(i - 1, j)
                && dp[ps] > NEG_INF * 0.5
            {
                let pen = if bt[ps] == 2 {
                    GAP_EXTEND_BITS
                } else {
                    GAP_OPEN_BITS
                };
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
                let pen = if bt[ps] == 3 {
                    GAP_EXTEND_BITS
                } else {
                    GAP_OPEN_BITS
                };
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
    pairs
}
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
    let span = MASK_KP14.len();
    let max_gap_bp = 2400;
    let max_local_rescue_bp = 40_000;
    let seq_x = &idx_x.seq;
    let seq_y = &idx_y.seq;
    let seq_y_rc = &idx_y.seq_rc;
    let pv_len = seq_y.len();

    let (raw_fwd, raw_rc) = find_w12_raw_lcbs(idx_x, idx_y);

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

    // Sparse-exon fallback: w=10 exact seeds in intervals not covered by w=12 LCBs.
    let ix10 = canonical_spaced_seeds_with_opts(
        seq_x,
        MASK_W10,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
    )?;
    let iy10_fwd = canonical_spaced_seeds_with_opts(
        seq_y,
        MASK_W10,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
    )?;
    let iy10_rc = canonical_spaced_seeds_with_opts(
        seq_y_rc,
        MASK_W10,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
    )?;
    let (hits10_fwd, hits10_rc) = intersect_seed_maps(&ix10, &iy10_fwd, &iy10_rc);
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

    // LASTZ-style 12-of-19 + 1-transition hits gated by ungapped X-drop HSP score (>=180).
    let iy_lz_fwd = canonical_spaced_seeds_with_opts(
        seq_y,
        MASK_LASTZ12,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
    )?;
    let iy_lz_rc = canonical_spaced_seeds_with_opts(
        seq_y_rc,
        MASK_LASTZ12,
        cfg.fallback_max_occ,
        cfg.skip_softmasked,
        false,
    )?;
    let (hits_lz_fwd, hits_lz_rc) =
        find_lastz12_1transition_hits(seq_x, &iy_lz_fwd, &iy_lz_rc, cfg.skip_softmasked);
    append_hsp_gated_lastz12_lcbs(&mut combined_fwd, seq_x, seq_y, &hits_lz_fwd, 180);
    append_hsp_gated_lastz12_lcbs(&mut combined_rc, seq_x, seq_y_rc, &hits_lz_rc, 180);

    if cfg.use_multi_shape_seeds {
        let (double_trans_fwd, double_trans_rc) =
            find_lastz12_2transition_hits(seq_x, &iy_lz_fwd, &iy_lz_rc, cfg.skip_softmasked);
        append_hsp_gated_lastz12_lcbs(&mut combined_fwd, seq_x, seq_y, &double_trans_fwd, 210);
        append_hsp_gated_lastz12_lcbs(&mut combined_rc, seq_x, seq_y_rc, &double_trans_rc, 210);

        let iy_11_fwd = canonical_spaced_seeds_with_opts(
            seq_y,
            MASK_LASTZ11,
            cfg.fallback_max_occ,
            cfg.skip_softmasked,
            false,
        )?;
        let iy_11_rc = canonical_spaced_seeds_with_opts(
            seq_y_rc,
            MASK_LASTZ11,
            cfg.fallback_max_occ,
            cfg.skip_softmasked,
            false,
        )?;
        let (hits_11_fwd, hits_11_rc) =
            find_lastz11_1transition_hits(seq_x, &iy_11_fwd, &iy_11_rc, cfg.skip_softmasked);
        append_hsp_gated_lastz12_lcbs(&mut combined_fwd, seq_x, seq_y, &hits_11_fwd, 210);
        append_hsp_gated_lastz12_lcbs(&mut combined_rc, seq_x, seq_y_rc, &hits_11_rc, 210);

        // HSP-gated 1-bit purine/pyrimidine (R/Y) seeding for multi-transition divergent homologs.
        let ix_ry =
            canonical_spaced_seeds_with_opts(seq_x, MASK_RY20, 2, cfg.skip_softmasked, true)?;
        let iy_ry_fwd =
            canonical_spaced_seeds_with_opts(seq_y, MASK_RY20, 2, cfg.skip_softmasked, true)?;
        let iy_ry_rc =
            canonical_spaced_seeds_with_opts(seq_y_rc, MASK_RY20, 2, cfg.skip_softmasked, true)?;
        let (hits_ry_fwd, hits_ry_rc) = intersect_seed_maps(&ix_ry, &iy_ry_fwd, &iy_ry_rc);
        append_hsp_gated_lastz12_lcbs(&mut combined_fwd, seq_x, seq_y, &hits_ry_fwd, 220);
        append_hsp_gated_lastz12_lcbs(&mut combined_rc, seq_x, seq_y_rc, &hits_ry_rc, 220);

        // HSP-gated w=10 (10-of-14) 1-transition seeding.
        let (hits_w10_fwd, hits_w10_rc) =
            find_w10_1transition_hits(seq_x, &iy10_fwd, &iy10_rc, cfg.skip_softmasked);
        append_hsp_gated_lastz12_lcbs(&mut combined_fwd, seq_x, seq_y, &hits_w10_fwd, 220);
        append_hsp_gated_lastz12_lcbs(&mut combined_rc, seq_x, seq_y_rc, &hits_w10_rc, 220);

        // Amino-acid HSP-gated 6-frame translated spaced 5-of-6 seeding for divergent coding homologs.
        let aa_fwd = find_translated_aa5_hits(seq_x, seq_y, 4, cfg.skip_softmasked);
        let aa_rc = find_translated_aa5_hits(seq_x, seq_y_rc, 4, cfg.skip_softmasked);
        append_hsp_gated_aa_lcbs(&mut combined_fwd, seq_x, seq_y, &aa_fwd, 90);
        append_hsp_gated_aa_lcbs(&mut combined_rc, seq_x, seq_y_rc, &aa_rc, 90);
    }

    // Transition-tolerant R/Y (1-bit purine/pyrimidine) fallback for divergent homologs.
    if cfg.use_ry_seeds {
        let ix_ry =
            canonical_spaced_seeds_with_opts(seq_x, MASK_RY20, 2, cfg.skip_softmasked, true)?;
        let iy_ry_fwd =
            canonical_spaced_seeds_with_opts(seq_y, MASK_RY20, 2, cfg.skip_softmasked, true)?;
        let iy_ry_rc =
            canonical_spaced_seeds_with_opts(seq_y_rc, MASK_RY20, 2, cfg.skip_softmasked, true)?;
        let (hits_ry_fwd, hits_ry_rc) = intersect_seed_maps(&ix_ry, &iy_ry_fwd, &iy_ry_rc);
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

    let mut pair_set: HashSet<(usize, usize)> = HashSet::new();
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
                && (100..=120_000).contains(&(x1 - x0))
                && (100..=120_000).contains(&(y1 - y0))
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
                    let py = if is_rc {
                        pv_len - 1 - (ay + off)
                    } else {
                        ay + off
                    };
                    pair_set.insert((ax + off, py));
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

    // Run banded MDL DP gaps in parallel via rayon.
    let dp_results: Vec<(bool, usize, usize, GapResult)> = dp_jobs
        .par_iter()
        .map(|&(is_rc, x0, x1, y0, y1)| {
            let y_work = if is_rc { seq_y_rc } else { seq_y };
            let res = banded_mdl_gap(&seq_x[x0..x1], &y_work[y0..y1], 180, 10.0);
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
            pair_set.insert((x0 + dx, py));
        }
    }

    let xdrop_results: Vec<(bool, Vec<(usize, usize)>)> = xdrop_jobs
        .par_iter()
        .map(|&(is_rc, x0, x1, y0, y1, rev)| {
            let y_work = if is_rc { seq_y_rc } else { seq_y };
            let pairs = mdl_xdrop_extend(&seq_x[x0..x1], &y_work[y0..y1], x0, y0, rev, 1500, 40);
            (is_rc, pairs)
        })
        .collect();

    for (is_rc, pairs) in xdrop_results {
        for (px, py) in pairs {
            let y_orig = if is_rc { pv_len - 1 - py } else { py };
            pair_set.insert((px, y_orig));
        }
    }

    let mut aligned_pairs: Vec<(usize, usize)> = pair_set.into_iter().collect();
    aligned_pairs.sort_unstable();

    Ok(ContigAlignment {
        aligned_pairs,
        lcbs_fwd: lcbs_fwd.len(),
        lcbs_rc: lcbs_rc.len(),
        macro_bubbles,
        internal_bubbles,
    })
}
