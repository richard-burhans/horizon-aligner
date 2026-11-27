//! Parameter structs, band geometry, and scratch buffers shared by banded MDL+BOCPD DP,
//! Viterbi fast paths, and X-drop flank extensions.

use crate::bocpd::{BocpdScratch, GAP_EXTEND_BITS, GAP_OPEN_BITS, IDENTITY_PRIOR, NEG_INF};

/// Divergence-adaptive Pair-HMM parameters (in probability / bit units).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HmmParams {
    /// Expected nucleotide identity in homologous match columns (`0.50..=0.98`).
    pub identity_prior: f64,
    /// Log2 probability of opening an insertion or deletion (`< -1.0`).
    pub gap_open_bits: f64,
    /// Log2 probability of extending an insertion or deletion (`< 0.0`).
    pub gap_extend_bits: f64,
}

impl Default for HmmParams {
    fn default() -> Self {
        Self {
            identity_prior: IDENTITY_PRIOR,
            gap_open_bits: GAP_OPEN_BITS,
            gap_extend_bits: GAP_EXTEND_BITS,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct NormalizedHmmParams {
    pub id_prior: f64,
    pub go_bits: f64,
    pub ge_bits: f64,
    pub p_open: f64,
    pub p_ext: f64,
    pub p_mm: f64,
    pub p_im: f64,
    pub mm_bits: f64,
    pub im_bits: f64,
}

impl HmmParams {
    #[must_use]
    #[allow(clippy::similar_names)]
    pub(crate) fn normalized(&self) -> NormalizedHmmParams {
        let id_raw = if self.identity_prior.is_finite() {
            self.identity_prior
        } else {
            IDENTITY_PRIOR
        };
        let go_raw = if self.gap_open_bits.is_finite() {
            self.gap_open_bits
        } else {
            GAP_OPEN_BITS
        };
        let ge_raw = if self.gap_extend_bits.is_finite() {
            self.gap_extend_bits
        } else {
            GAP_EXTEND_BITS
        };
        let id_prior = id_raw.clamp(0.40, 0.98);
        let go_bits = go_raw.clamp(-8.0, -1.5);
        let ge_bits = ge_raw.clamp(-4.0, -0.2);
        let p_open = 2.0_f64.powf(go_bits);
        let p_ext = 2.0_f64.powf(ge_bits);
        let p_mm = (1.0 - 2.0 * p_open).max(1e-6);
        let p_im = (1.0 - p_ext).max(1e-6);
        NormalizedHmmParams {
            id_prior,
            go_bits,
            ge_bits,
            p_open,
            p_ext,
            p_mm,
            p_im,
            mm_bits: p_mm.log2(),
            im_bits: p_im.log2(),
        }
    }
}

/// Result of aligning one inter-anchor gap.
#[derive(Debug, Clone, PartialEq)]
pub struct GapResult {
    /// Kept base-pair coordinates `(x, y)` relative to the gap start.
    pub kept_pairs: Vec<(usize, usize)>,
    /// Number of abstained Horizon Bubbles in this gap.
    pub bubbles: usize,
    /// Total Forward partition-function bit-savings score for the gap.
    pub score_bits: f64,
}

/// One MEA path column: `(x_idx, y_idx, emit_or_gap_bits, post_prob, abstained)`.
pub type MeaTraceColumn = (Option<usize>, Option<usize>, f64, f64, bool);

/// Full Forward-Backward posterior grid and MEA trace for diagnostics and plotting.
#[derive(Debug, Clone, PartialEq)]
pub struct GapPosteriorTrace {
    /// Target sequence length.
    pub n: usize,
    /// Query sequence length.
    pub m: usize,
    /// Effective clamped band half-width.
    pub band_w: usize,
    /// Total Forward-Backward MDL bit-savings (`log2(Z_F)`).
    pub score_bits: f64,
    /// Backward partition-function log2 score (`log2(Z_B)`); matches `score_bits` up to float rounding.
    pub log2_zb: f64,
    /// BOCPD local GC profile along target `X`.
    pub gc_x: Vec<f64>,
    /// BOCPD local GC profile along query `Y`.
    pub gc_y: Vec<f64>,
    /// In-band match posteriors `(i, j, P(M_{i,j} | X, Y))` in 1-based DP coordinates.
    pub posterior_cells: Vec<(usize, usize, f64)>,
    /// MEA path columns `(x_idx, y_idx, emit_or_gap_bits, post_prob, abstained)`.
    pub mea_cols: Vec<MeaTraceColumn>,
    /// Kept base-pair coordinates after posterior and block-bits filtering.
    pub kept_pairs: Vec<(usize, usize)>,
    /// Number of abstained Horizon Bubbles.
    pub bubbles: usize,
    /// Maximum contiguous block bit-savings along the MEA path.
    pub max_block_bits: f64,
}

/// Reusable per-worker scratch buffers for [`super::dp::banded_mdl_gap_with_scratch`].
#[derive(Debug, Default)]
pub(crate) struct DpScratch {
    pub(super) bocpd: BocpdScratch,
    pub(super) gc_x: Vec<f64>,
    pub(super) gc_y: Vec<f64>,
    pub(super) pi_x: Vec<[f64; 4]>,
    pub(super) pi_y: Vec<[f64; 4]>,
    pub(super) enc_x: Vec<u8>,
    pub(super) enc_y: Vec<u8>,
    pub(super) log2_scale_f_prefix: Vec<f64>,
    pub(super) log2_scale_b_suffix: Vec<f64>,
    pub(super) fm: Vec<f64>,
    pub(super) fix: Vec<f64>,
    pub(super) fiy: Vec<f64>,
    pub(super) sub_grid: Vec<f64>,
    pub(super) bm: Vec<f64>,
    pub(super) bix: Vec<f64>,
    pub(super) biy: Vec<f64>,
    pub(super) raw_cols: Vec<(Option<usize>, Option<usize>, f64, f64)>,
    pub(super) abstain: Vec<bool>,
    pub(super) vit_m: Vec<f32>,
    pub(super) vit_x: Vec<f32>,
    pub(super) vit_y: Vec<f32>,
    pub(super) vit_trace: Vec<u8>,
    pub(super) last_max_block_bits: f64,
    pub(super) last_log2_zb: f64,
}

pub(crate) fn prepare_grid<T: Copy>(buf: &mut Vec<T>, total: usize, init: T) -> bool {
    buf.clear();
    if buf.capacity() < total && buf.try_reserve_exact(total).is_err() {
        return false;
    }
    buf.resize(total, init);
    true
}

#[inline]
pub(crate) fn band_center(i: usize, m: usize, n: usize) -> usize {
    usize::try_from((i as u128 * m as u128) / (n.max(1) as u128)).unwrap_or(m)
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct BandLayout {
    pub(super) n: usize,
    pub(super) m: usize,
    pub(super) band_w: usize,
    pub(super) stride: usize,
    pub(super) total: usize,
    zero_centered: bool,
}

impl BandLayout {
    pub(super) fn new(n: usize, m: usize, band_w: usize) -> Option<Self> {
        let bw = band_w.min(n.max(m));
        let stride = bw.checked_mul(2)?.checked_add(1)?;
        let total = n.checked_add(1)?.checked_mul(stride)?;
        Some(Self {
            n,
            m,
            band_w: bw,
            stride,
            total,
            zero_centered: false,
        })
    }

    pub(super) fn new_zero_centered(n: usize, m: usize, band_w: usize) -> Option<Self> {
        let mut layout = Self::new(n, m, band_w)?;
        layout.zero_centered = true;
        Some(layout)
    }

    #[inline]
    pub(super) fn center(&self, i: usize) -> usize {
        if self.zero_centered {
            i
        } else {
            band_center(i, self.m, self.n)
        }
    }

    #[inline]
    pub(super) fn slot(&self, i: usize, j: usize) -> Option<usize> {
        #[allow(clippy::cast_possible_wrap)]
        let diff = (j as isize) - (self.center(i) as isize);
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
    pub(super) fn get(&self, data: &[f64], i: usize, j: usize) -> f64 {
        self.slot(i, j).map_or(NEG_INF, |idx| data[idx])
    }
}

#[must_use]
pub(crate) fn dip_floor_bits_for_identity(identity: f64) -> f64 {
    let p_id = if identity.is_finite() {
        identity
    } else {
        IDENTITY_PRIOR
    }
    .clamp(0.25, 0.999);
    crate::bocpd::joint_pair_bits(b'A', b'C', [0.25; 4], [0.25; 4], p_id) - 0.5
}
