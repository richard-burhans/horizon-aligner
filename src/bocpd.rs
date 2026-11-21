//! BOCPD compositional segmentation and MDL log-likelihood scoring primitives.

pub(crate) const NEG_INF: f64 = -1e18;
/// Default affine gap-open penalty (bits).
pub const GAP_OPEN_BITS: f64 = -4.5;
/// Default affine gap-extend penalty (bits).
pub const GAP_EXTEND_BITS: f64 = -0.8;
/// Default expected homologous nucleotide identity prior.
pub const IDENTITY_PRIOR: f64 = 0.65;
pub(crate) const MIN_POST_PROB: f64 = 0.15;

const BASE_IDX_LUT: [u8; 256] = {
    let mut lut = [4u8; 256];
    lut[b'A' as usize] = 0;
    lut[b'a' as usize] = 0;
    lut[b'C' as usize] = 1;
    lut[b'c' as usize] = 1;
    lut[b'G' as usize] = 2;
    lut[b'g' as usize] = 2;
    lut[b'T' as usize] = 3;
    lut[b't' as usize] = 3;
    lut
};

#[inline]
pub(crate) fn base_to_idx_fast(b: u8) -> u8 {
    BASE_IDX_LUT[b as usize]
}

pub(crate) fn base_to_idx(b: u8) -> Option<usize> {
    let idx = BASE_IDX_LUT[b as usize];
    if idx < 4 { Some(idx as usize) } else { None }
}

/// Background nucleotide probabilities `[A, C, G, T]` for a given GC fraction.
#[inline]
#[must_use]
pub fn base_probs(gc: f64) -> [f64; 4] {
    if gc.is_nan() {
        return [0.25; 4];
    }
    let gc = gc.clamp(0.02, 0.98);
    let p_gc = gc * 0.5;
    let p_at = (1.0 - gc) * 0.5;
    // Order matches A=0, C=1, G=2, T=3.
    [p_at, p_gc, p_gc, p_at]
}

/// Pair-HMM log2 likelihood ratio (bits) for aligned bases `(a, b)` under background
/// compositions `pi_x`, `pi_y` and homologous identity prior `identity`.
#[inline]
#[must_use]
pub fn joint_pair_bits(a: u8, b: u8, pi_x: [f64; 4], pi_y: [f64; 4], identity: f64) -> f64 {
    let (Some(ia), Some(ib)) = (base_to_idx(a), base_to_idx(b)) else {
        return 0.0;
    };
    let pa = pi_x[ia];
    let pb = pi_y[ib];
    let p_null = pa * pb;
    if p_null <= 0.0 {
        return 0.0;
    }
    let id = if identity.is_finite() {
        identity.clamp(0.0, 1.0)
    } else {
        IDENTITY_PRIOR
    };
    let p_avg = f64::midpoint(pa, pb);
    let match_term = if ia == ib { id * p_avg } else { 0.0 };
    let p_joint = match_term + (1.0 - id) * p_null;
    (p_joint / p_null).log2()
}

#[inline]
pub(crate) fn logsumexp2(a: f64, b: f64) -> f64 {
    let m = a.max(b);
    if m < -1e15 {
        NEG_INF
    } else {
        m + ((a - m).exp() + (b - m).exp()).ln()
    }
}

#[inline]
pub(crate) fn logsumexp3(a: f64, b: f64, c: f64) -> f64 {
    let m = a.max(b).max(c);
    if m < -1e15 {
        NEG_INF
    } else {
        m + ((a - m).exp() + (b - m).exp() + (c - m).exp()).ln()
    }
}

/// Reusable scratch buffers for [`bocpd_gc_profile_into`].
#[derive(Debug, Default)]
pub struct BocpdScratch {
    run_probs: Vec<f64>,
    alphas: Vec<f64>,
    betas: Vec<f64>,
    next_probs: Vec<f64>,
    next_alphas: Vec<f64>,
    next_betas: Vec<f64>,
}

/// Estimate local GC probability per base using Bayesian Online Changepoint Detection,
/// writing into `gc_est` and reusing `scratch` to avoid per-call allocations.
pub fn bocpd_gc_profile_into(
    seq: &[u8],
    hazard: f64,
    max_run: usize,
    scratch: &mut BocpdScratch,
    gc_est: &mut Vec<f64>,
) {
    scratch.run_probs.clear();
    scratch.alphas.clear();
    scratch.betas.clear();
    scratch.run_probs.push(1.0);
    scratch.alphas.push(1.0);
    scratch.betas.push(1.0);

    gc_est.clear();
    gc_est.reserve(seq.len());

    let hazard = if hazard.is_nan() {
        0.0
    } else {
        hazard.clamp(0.0, 1.0)
    };
    let one_minus_h = 1.0 - hazard;
    let keep_cap = max_run.max(2) - 1;

    for &b in seq {
        let idx_base = BASE_IDX_LUT[b as usize];
        let is_acgt = idx_base < 4;
        let is_gc = idx_base == 1 || idx_base == 2;
        let len = scratch.run_probs.len();
        let mut pred_gc: f64 = 0.0;
        if !is_acgt {
            for idx in 0..len {
                let a = scratch.alphas[idx];
                let bt = scratch.betas[idx];
                pred_gc += scratch.run_probs[idx] * (a / (a + bt));
            }
            gc_est.push(pred_gc.clamp(0.02, 0.98));
            continue;
        }

        scratch.next_probs.clear();
        scratch.next_alphas.clear();
        scratch.next_betas.clear();
        let gc_inc = if is_gc { 1.0 } else { 0.0 };
        let at_inc = 1.0 - gc_inc;
        scratch.next_probs.push(0.0);
        scratch.next_alphas.push(1.0 + gc_inc);
        scratch.next_betas.push(1.0 + at_inc);

        let mut cp_mass = 0.0;
        let keep_len = len.min(keep_cap);
        for idx in 0..len {
            let a = scratch.alphas[idx];
            let bt = scratch.betas[idx];
            let p_gc = a / (a + bt);
            let rp = scratch.run_probs[idx];
            pred_gc += rp * p_gc;
            let pred_x = if is_gc { p_gc } else { 1.0 - p_gc };
            let joint = rp * pred_x;
            cp_mass += joint * hazard;
            if idx < keep_len {
                scratch.next_probs.push(joint * one_minus_h);
                scratch.next_alphas.push(a + gc_inc);
                scratch.next_betas.push(bt + at_inc);
            } else if let Some(last_p) = scratch.next_probs.last_mut() {
                // Fold tail run-length probability into the oldest retained bin instead of dropping mass.
                *last_p += joint * one_minus_h;
            }
        }
        gc_est.push(pred_gc.clamp(0.02, 0.98));
        scratch.next_probs[0] = cp_mass;

        let norm: f64 = scratch.next_probs.iter().sum();
        if norm > 0.0 {
            let inv = 1.0 / norm;
            for p in &mut scratch.next_probs {
                *p *= inv;
            }
        } else {
            scratch.next_probs.fill(0.0);
            scratch.next_probs[0] = 1.0;
        }

        std::mem::swap(&mut scratch.run_probs, &mut scratch.next_probs);
        std::mem::swap(&mut scratch.alphas, &mut scratch.next_alphas);
        std::mem::swap(&mut scratch.betas, &mut scratch.next_betas);
    }
}

/// Estimate local GC probability per base using Bayesian Online Changepoint Detection.
#[must_use]
pub fn bocpd_gc_profile(seq: &[u8], hazard: f64, max_run: usize) -> Vec<f64> {
    let mut scratch = BocpdScratch::default();
    let mut gc_est = Vec::with_capacity(seq.len());
    bocpd_gc_profile_into(seq, hazard, max_run, &mut scratch, &mut gc_est);
    gc_est
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bocpd_ignores_n_bases() {
        let with_n = b"GCGCGCGCGCGCGCGCGCGCNNNNNNNNNNNNNNNNNNNN";
        let gc = bocpd_gc_profile(with_n, 1.0 / 120.0, 200);
        assert_eq!(gc.len(), with_n.len());
        assert!(gc[20] > 0.70);
        assert!((gc[39] - gc[20]).abs() < 1e-12);
    }
}
