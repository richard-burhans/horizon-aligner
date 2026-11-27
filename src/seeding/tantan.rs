//! `tantan`-style Forward-Backward HMM probabilistic tandem-repeat masker (Frith 2011).

use crate::dna::base_to_idx;

/// Chunk size (bp) for bounded-RAM overlapping-window Forward-Backward decoding (~12 MiB peak).
const TANTAN_CHUNK_BP: usize = 65_536;

/// Configuration for [`tantan_mask`] / [`tantan_posterior_probs`].
#[derive(Debug, Clone, Copy)]
pub struct TantanConfig {
    /// Maximum tandem-repeat period `P` modeled by the HMM (default: `24`).
    pub max_period: usize,
    /// Transition probability from background into a tandem repeat (`r`, default: `0.005`).
    pub repeat_open_prob: f64,
    /// Transition probability from a tandem repeat back to background (`e`, default: `0.05`).
    pub repeat_end_prob: f64,
    /// Per-base identity probability between adjacent repeat units (default: `0.85`).
    pub repeat_match_prob: f64,
    /// Posterior probability threshold above which a base is masked (default: `0.5`).
    pub mask_threshold: f64,
}

impl Default for TantanConfig {
    fn default() -> Self {
        Self {
            max_period: 24,
            repeat_open_prob: 0.005,
            repeat_end_prob: 0.05,
            repeat_match_prob: 0.85,
            mask_threshold: 0.5,
        }
    }
}

struct TantanParams {
    max_p: usize,
    open_per_p: f64,
    stay_bg: f64,
    end_prob: f64,
    stay_rep: f64,
    lr_match: f64,
    lr_mismatch: f64,
}

impl TantanParams {
    fn from_cfg(cfg: &TantanConfig) -> Self {
        let def = TantanConfig::default();
        let max_p = cfg.max_period.max(1);
        let open_prob = if cfg.repeat_open_prob.is_finite() {
            cfg.repeat_open_prob
        } else {
            def.repeat_open_prob
        }
        .clamp(1e-6, 0.5);
        let end_prob = if cfg.repeat_end_prob.is_finite() {
            cfg.repeat_end_prob
        } else {
            def.repeat_end_prob
        }
        .clamp(1e-6, 0.99);
        let p_match = if cfg.repeat_match_prob.is_finite() {
            cfg.repeat_match_prob
        } else {
            def.repeat_match_prob
        }
        .clamp(0.30, 0.999);
        #[allow(clippy::cast_precision_loss)]
        let open_per_p = open_prob / (max_p as f64);
        Self {
            max_p,
            open_per_p,
            stay_bg: 1.0 - open_prob,
            end_prob,
            stay_rep: 1.0 - end_prob,
            lr_match: p_match / 0.25,
            lr_mismatch: ((1.0 - p_match) / 3.0) / 0.25,
        }
    }
}

/// Run exact scaled Forward-Backward HMM over one contiguous slice `seq`.
fn tantan_window_probs(seq: &[u8], p: &TantanParams) -> Vec<f64> {
    let seq_len = seq.len();
    if seq_len == 0 {
        return Vec::new();
    }
    let max_p = p.max_p;
    // Compact 1 B/bp encoding: 0..=3 for ACGT, 4 for N/other (16x smaller than Option<usize>).
    let codes: Vec<u8> = seq
        .iter()
        .map(|&b| base_to_idx(b).map_or(4, |idx| u8::try_from(idx).unwrap_or(4)))
        .collect();
    let emit_lr = |pos: usize, period: usize| -> f64 {
        if pos < period {
            return 0.0;
        }
        let lhs = codes[pos];
        let rhs = codes[pos - period];
        if lhs < 4 && rhs < 4 {
            if lhs == rhs {
                p.lr_match
            } else {
                p.lr_mismatch
            }
        } else {
            1.0
        }
    };

    let n_states = max_p + 1;
    let mut fwd = vec![0.0_f64; seq_len * n_states];
    let mut scales = vec![1.0_f64; seq_len];

    fwd[0] = 1.0;
    for (i, scale_slot) in scales.iter_mut().enumerate().skip(1) {
        let prev = (i - 1) * n_states;
        let cur = i * n_states;

        let mut from_rep_sum = 0.0_f64;
        for per in 1..=max_p {
            from_rep_sum += fwd[prev + per];
        }
        let bg_val = fwd[prev] * p.stay_bg + from_rep_sum * p.end_prob;
        fwd[cur] = bg_val;
        let mut total = bg_val;

        let from_bg = fwd[prev] * p.open_per_p;
        for per in 1..=max_p {
            let val = (from_bg + fwd[prev + per] * p.stay_rep) * emit_lr(i, per);
            fwd[cur + per] = val;
            total += val;
        }
        let sc = if total > 0.0 { total } else { 1.0 };
        *scale_slot = sc;
        let inv = 1.0 / sc;
        for s in 0..n_states {
            fwd[cur + s] *= inv;
        }
    }

    let mut bwd_next = vec![1.0_f64; n_states];
    let mut bwd_cur = vec![0.0_f64; n_states];
    let mut probs = vec![0.0_f64; seq_len];
    probs[seq_len - 1] = (1.0 - fwd[(seq_len - 1) * n_states]).clamp(0.0, 1.0);

    for i in (0..(seq_len - 1)).rev() {
        let b0_next = bwd_next[0];
        let mut to_rep_sum = 0.0_f64;
        // Single pass over period states computes both to_rep_sum and bwd_cur[per].
        for per in 1..=max_p {
            let b_em = bwd_next[per] * emit_lr(i + 1, per);
            to_rep_sum += b_em;
            bwd_cur[per] = p.end_prob * b0_next + p.stay_rep * b_em;
        }
        bwd_cur[0] = p.stay_bg * b0_next + p.open_per_p * to_rep_sum;

        let inv = 1.0 / scales[i + 1];
        for v in &mut bwd_cur {
            *v *= inv;
        }
        let cur = i * n_states;
        let mut norm_z = 0.0_f64;
        for s in 0..n_states {
            norm_z += fwd[cur + s] * bwd_cur[s];
        }
        let p_bg = if norm_z > 0.0 {
            (fwd[cur] * bwd_cur[0]) / norm_z
        } else {
            1.0
        };
        probs[i] = (1.0 - p_bg).clamp(0.0, 1.0);
        std::mem::swap(&mut bwd_next, &mut bwd_cur);
    }
    probs
}

/// Compute Forward-Backward posterior tandem-repeat probability `P(repeat_i | seq)`
/// at every position `i` (Frith 2011).
///
/// Sequences longer than 64 kb are decoded in overlapping windows (`margin >= 8 * max_period`)
/// so peak memory stays bounded at `O(chunk * max_period)` (~12 MiB) even on whole chromosomes.
#[must_use]
pub fn tantan_posterior_probs(seq: &[u8], cfg: &TantanConfig) -> Vec<f64> {
    let seq_len = seq.len();
    if seq_len == 0 {
        return Vec::new();
    }
    let params = TantanParams::from_cfg(cfg);
    if seq_len <= TANTAN_CHUNK_BP {
        return tantan_window_probs(seq, &params);
    }

    let margin = (8 * params.max_p).max(256);
    let mut out = vec![0.0_f64; seq_len];
    let mut start = 0usize;
    while start < seq_len {
        let end = (start + TANTAN_CHUNK_BP).min(seq_len);
        let win_lo = start.saturating_sub(margin);
        let win_hi = (end + margin).min(seq_len);
        let win_probs = tantan_window_probs(&seq[win_lo..win_hi], &params);
        let copy_lo = start - win_lo;
        let copy_hi = end - win_lo;
        out[start..end].copy_from_slice(&win_probs[copy_lo..copy_hi]);
        start = end;
    }
    out
}

/// Mask positions with `P(repeat_i | seq) >= cfg.mask_threshold`.
#[must_use]
pub fn tantan_mask(seq: &[u8], cfg: &TantanConfig) -> Vec<bool> {
    let thresh = if cfg.mask_threshold.is_finite() {
        cfg.mask_threshold
    } else {
        0.5
    };
    tantan_posterior_probs(seq, cfg)
        .into_iter()
        .map(|p| p >= thresh)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::splitmix64;

    #[test]
    fn tantan_masks_imperfect_period_7_tandem_repeat_and_spares_unique_flanks() {
        let mut st = 0x7A7A_1234_5678_9ABC_u64;
        let unique = |n: usize, s: &mut u64| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    *s = splitmix64(*s);
                    b"ACGT"[(*s & 3) as usize]
                })
                .collect()
        };
        let flank_l = unique(80, &mut st);
        let mut repeat = Vec::new();
        for copy in 0..14usize {
            let mut unit = *b"ACGTTGC";
            if copy % 4 == 2 {
                unit[3] = b'A';
            }
            repeat.extend_from_slice(&unit);
        }
        let flank_r = unique(80, &mut st);

        let mut seq = Vec::new();
        seq.extend_from_slice(&flank_l);
        let rep_start = seq.len();
        seq.extend_from_slice(&repeat);
        let rep_end = seq.len();
        seq.extend_from_slice(&flank_r);

        let cfg = TantanConfig::default();
        let mask = tantan_mask(&seq, &cfg);
        assert_eq!(mask.len(), seq.len());

        let masked_flank_l = mask[..60].iter().filter(|&&m| m).count();
        let masked_flank_r = mask[rep_end + 20..].iter().filter(|&&m| m).count();
        let masked_rep = mask[rep_start + 7..rep_end - 7]
            .iter()
            .filter(|&&m| m)
            .count();
        let rep_core = (rep_end - 7) - (rep_start + 7);
        assert!(
            masked_flank_l <= 5 && masked_flank_r <= 5,
            "unique flanks should remain unmasked (L={masked_flank_l}, R={masked_flank_r})"
        );
        assert!(
            masked_rep * 10 >= rep_core * 8,
            "expected >= 80% of period-7 tandem repeat core masked by tantan, got {masked_rep}/{rep_core}"
        );
    }

    #[test]
    fn tantan_handles_edge_inputs_nan_config_and_chunked_boundary_agreement() {
        let nan_cfg = TantanConfig {
            max_period: 0,
            repeat_open_prob: f64::NAN,
            repeat_end_prob: f64::NAN,
            repeat_match_prob: f64::NAN,
            mask_threshold: f64::NAN,
        };
        assert_eq!(tantan_posterior_probs(b"", &nan_cfg), Vec::<f64>::new());
        let single = tantan_posterior_probs(b"A", &nan_cfg);
        assert_eq!(single.len(), 1);
        assert!(single[0].is_finite());

        // Verify chunked decoding (> 64 kb) matches full-sequence HMM across the chunk boundary.
        let mut st = 0xDEAD_BEEF_1111_2222_u64;
        let mut seq: Vec<u8> = (0..(TANTAN_CHUNK_BP + 1024))
            .map(|_| {
                st = splitmix64(st);
                b"ACGT"[(st & 3) as usize]
            })
            .collect();
        // Embed a period-5 repeat straddling the 65,536-bp chunk boundary!
        for (pos, base) in seq
            .iter_mut()
            .enumerate()
            .take(TANTAN_CHUNK_BP + 40)
            .skip(TANTAN_CHUNK_BP - 40)
        {
            *base = b"ACGTT"[pos % 5];
        }
        let cfg = TantanConfig::default();
        let chunked = tantan_posterior_probs(&seq, &cfg);
        let exact = tantan_window_probs(&seq, &TantanParams::from_cfg(&cfg));
        let max_diff = chunked
            .iter()
            .zip(&exact)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_diff < 1e-4,
            "chunked tantan must agree with full-sequence HMM across boundary (max_diff={max_diff:.6})"
        );
    }
}
