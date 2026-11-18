//! BOCPD compositional segmentation and MDL log-likelihood scoring primitives.

pub(crate) const NEG_INF: f64 = -1e18;
pub(crate) const GAP_OPEN_BITS: f64 = -4.5;
pub(crate) const GAP_EXTEND_BITS: f64 = -0.8;
pub(crate) const IDENTITY_PRIOR: f64 = 0.65;
pub(crate) const MIN_POST_PROB: f64 = 0.15;

pub(crate) fn base_to_idx(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

#[inline]
pub(crate) fn base_probs(gc: f64) -> [f64; 4] {
    let p_gc = gc * 0.5;
    let p_at = (1.0 - gc) * 0.5;
    // Order matches A=0, C=1, G=2, T=3.
    [p_at, p_gc, p_gc, p_at]
}

#[inline]
pub(crate) fn joint_pair_bits(a: u8, b: u8, pi_x: [f64; 4], pi_y: [f64; 4], identity: f64) -> f64 {
    let (Some(ia), Some(ib)) = (base_to_idx(a), base_to_idx(b)) else {
        return -2.0;
    };
    let pa = pi_x[ia];
    let pb = pi_y[ib];
    let p_null = pa * pb;
    let p_avg = f64::midpoint(pa, pb);
    let match_term = if ia == ib { identity * p_avg } else { 0.0 };
    let p_joint = match_term + (1.0 - identity) * p_null;
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
    logsumexp2(logsumexp2(a, b), c)
}

/// Estimate local GC probability per base using Bayesian Online Changepoint Detection.
#[must_use]
pub fn bocpd_gc_profile(seq: &[u8], hazard: f64, max_run: usize) -> Vec<f64> {
    let mut run_probs = Vec::with_capacity(max_run);
    let mut alphas = Vec::with_capacity(max_run);
    let mut betas = Vec::with_capacity(max_run);
    run_probs.push(1.0);
    alphas.push(1.0);
    betas.push(1.0);

    let mut next_probs = Vec::with_capacity(max_run);
    let mut next_alphas = Vec::with_capacity(max_run);
    let mut next_betas = Vec::with_capacity(max_run);
    let mut gc_est = Vec::with_capacity(seq.len());

    for &b in seq {
        let is_gc = matches!(b.to_ascii_uppercase(), b'G' | b'C');
        let mut pred_gc: f64 = 0.0;
        for idx in 0..run_probs.len() {
            pred_gc += run_probs[idx] * (alphas[idx] / (alphas[idx] + betas[idx]));
        }
        gc_est.push(pred_gc.clamp(0.02, 0.98));

        next_probs.clear();
        next_alphas.clear();
        next_betas.clear();
        next_probs.push(0.0);
        next_alphas.push(1.0 + if is_gc { 1.0 } else { 0.0 });
        next_betas.push(1.0 + if is_gc { 0.0 } else { 1.0 });

        let mut cp_mass = 0.0;
        let keep_len = run_probs.len().min(max_run.saturating_sub(1));
        for idx in 0..run_probs.len() {
            let a = alphas[idx];
            let bt = betas[idx];
            let pred_x = if is_gc { a / (a + bt) } else { bt / (a + bt) };
            let joint = run_probs[idx] * pred_x;
            cp_mass += joint * hazard;
            if idx < keep_len {
                next_probs.push(joint * (1.0 - hazard));
                next_alphas.push(a + if is_gc { 1.0 } else { 0.0 });
                next_betas.push(bt + if is_gc { 0.0 } else { 1.0 });
            }
        }
        next_probs[0] = cp_mass;

        let norm: f64 = next_probs.iter().sum();
        if norm > 0.0 {
            let inv = 1.0 / norm;
            for p in &mut next_probs {
                *p *= inv;
            }
        } else {
            next_probs.fill(0.0);
            next_probs[0] = 1.0;
        }

        std::mem::swap(&mut run_probs, &mut next_probs);
        std::mem::swap(&mut alphas, &mut next_alphas);
        std::mem::swap(&mut betas, &mut next_betas);
    }

    gc_est
}
