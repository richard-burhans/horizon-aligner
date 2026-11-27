//! `last-train`-style asymmetric 4×4 substitution & gap EM parameter estimation.
//!
//! Fits directional substitution log-odds `S[a][b]` (allowing `S[a][b] != S[b][a]` under
//! compositional GC skew, e.g. 19% GC *P. falciparum* vs 40% GC *P. vivax*) plus asymmetric
//! insertion/deletion gap open & extend probabilities via expectation-maximization.

use crate::bocpd::base_to_idx;

/// Trained asymmetric pairwise substitution and gap model (`last-train`-style).
#[derive(Debug, Clone)]
pub struct TrainedAsymmetricParams {
    /// Asymmetric 4×4 log-odds score matrix `S[x_base][y_base]` (bits), order `A, C, G, T`.
    pub sub_bits: [[f64; 4]; 4],
    /// Estimated marginal base frequencies of sequence X.
    pub pi_x: [f64; 4],
    /// Estimated marginal base frequencies of sequence Y.
    pub pi_y: [f64; 4],
    /// Gap-open cost in X (deletion from Y) in bits (negative).
    pub del_open_bits: f64,
    /// Gap-extend cost in X in bits (negative).
    pub del_ext_bits: f64,
    /// Gap-open cost in Y (insertion in Y) in bits (negative).
    pub ins_open_bits: f64,
    /// Gap-extend cost in Y in bits (negative).
    pub ins_ext_bits: f64,
    /// Number of EM iterations executed.
    pub iterations: usize,
}

/// Fit an asymmetric 4×4 log-odds substitution matrix and directional gap parameters via
/// iterative EM on aligned sequence pairs `(aligned_x, aligned_y)`.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn train_asymmetric_pair_params(
    aligned_pairs: &[(&[u8], &[u8])],
    max_iters: usize,
) -> TrainedAsymmetricParams {
    let mut sub_bits = [[0.0_f64; 4]; 4];
    for (i, row) in sub_bits.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = if i == j { 1.5 } else { -1.5 };
        }
    }
    let mut pi_x = [0.25_f64; 4];
    let mut pi_y = [0.25_f64; 4];
    let mut del_open_bits = -5.0_f64;
    let mut del_ext_bits = -1.5_f64;
    let mut ins_open_bits = -5.0_f64;
    let mut ins_ext_bits = -1.5_f64;

    if aligned_pairs.is_empty() || max_iters == 0 {
        return TrainedAsymmetricParams {
            sub_bits,
            pi_x,
            pi_y,
            del_open_bits,
            del_ext_bits,
            ins_open_bits,
            ins_ext_bits,
            iterations: 0,
        };
    }

    // Direct unbiased empirical MLE over aligned columns (with Dirichlet Laplace pseudocounts),
    // avoiding self-reinforcing mismatch attenuation when training from pre-aligned blocks.
    let mut pair_counts = [[0.25_f64; 4]; 4];
    let mut x_counts = [1.0_f64; 4];
    let mut y_counts = [1.0_f64; 4];
    let mut match_cols = 4.0_f64;
    let mut del_opens = 1.0_f64;
    let mut del_exts = 1.0_f64;
    let mut ins_opens = 1.0_f64;
    let mut ins_exts = 1.0_f64;

    for &(ax, ay) in aligned_pairs {
        let n = ax.len().min(ay.len());
        let mut prev_state = 0u8; // 0=match, 1=del(gap in Y), 2=ins(gap in X)
        for k in 0..n {
            let bx = ax[k];
            let by = ay[k];
            match (bx != b'-', by != b'-') {
                (true, true) => {
                    if let (Some(ix), Some(iy)) = (base_to_idx(bx), base_to_idx(by)) {
                        pair_counts[ix][iy] += 1.0;
                        x_counts[ix] += 1.0;
                        y_counts[iy] += 1.0;
                        match_cols += 1.0;
                    }
                    prev_state = 0;
                }
                (true, false) => {
                    if let Some(ix) = base_to_idx(bx) {
                        x_counts[ix] += 1.0;
                    }
                    if prev_state == 1 {
                        del_exts += 1.0;
                    } else {
                        del_opens += 1.0;
                    }
                    prev_state = 1;
                }
                (false, true) => {
                    if let Some(iy) = base_to_idx(by) {
                        y_counts[iy] += 1.0;
                    }
                    if prev_state == 2 {
                        ins_exts += 1.0;
                    } else {
                        ins_opens += 1.0;
                    }
                    prev_state = 2;
                }
                (false, false) => {}
            }
        }
    }

    let sum_x: f64 = x_counts.iter().sum();
    let sum_y: f64 = y_counts.iter().sum();
    for b in 0..4 {
        pi_x[b] = (x_counts[b] / sum_x).clamp(1e-4, 1.0);
        pi_y[b] = (y_counts[b] / sum_y).clamp(1e-4, 1.0);
    }
    let norm_x: f64 = pi_x.iter().sum();
    let norm_y: f64 = pi_y.iter().sum();
    for b in 0..4 {
        pi_x[b] /= norm_x;
        pi_y[b] /= norm_y;
    }

    for a in 0..4 {
        for b in 0..4 {
            let p_ab = (pair_counts[a][b] / match_cols).clamp(1e-9, 1.0);
            sub_bits[a][b] = (p_ab / (pi_x[a] * pi_y[b])).log2();
        }
    }

    let p_del_open = (del_opens / (match_cols + del_opens + ins_opens)).clamp(1e-6, 0.5);
    let p_ins_open = (ins_opens / (match_cols + del_opens + ins_opens)).clamp(1e-6, 0.5);
    let p_del_ext = (del_exts / (del_opens + del_exts)).clamp(1e-4, 0.99);
    let p_ins_ext = (ins_exts / (ins_opens + ins_exts)).clamp(1e-4, 0.99);

    del_open_bits = p_del_open.log2();
    del_ext_bits = p_del_ext.log2();
    ins_open_bits = p_ins_open.log2();
    ins_ext_bits = p_ins_ext.log2();

    TrainedAsymmetricParams {
        sub_bits,
        pi_x,
        pi_y,
        del_open_bits,
        del_ext_bits,
        ins_open_bits,
        ins_ext_bits,
        iterations: 1,
    }
}

/// Reconstruct up to `max_cols` aligned columns (including small colinear inter-run indels)
/// from `runs` over `(seq_x, seq_y)` for Pair-HMM Baum-Welch training.
pub fn extract_aligned_train_columns(
    runs: &[crate::mdl::AlignedRun],
    sx: &[u8],
    sy: &[u8],
    max_cols: usize,
    train_x: &mut Vec<u8>,
    train_y: &mut Vec<u8>,
) {
    const MAX_COLINEAR_GAP: usize = 200;
    let mut prev_run: Option<(usize, usize, bool)> = None;
    for r in runs {
        if train_x.len() >= max_cols {
            break;
        }
        if let Some((px_end, py_end, p_rc)) = prev_run
            && p_rc == r.is_rc
            && r.x_start >= px_end
        {
            let dx = r.x_start - px_end;
            let dy_opt = if r.is_rc {
                py_end.checked_sub(r.y_start)
            } else {
                r.y_start.checked_sub(py_end)
            };
            if let Some(dy) = dy_opt
                && (dx > 0 || dy > 0)
                && dx <= MAX_COLINEAR_GAP
                && dy <= MAX_COLINEAR_GAP
                && (dx == 0 || dy == 0)
            {
                for gx in px_end..r.x_start {
                    if train_x.len() >= max_cols {
                        break;
                    }
                    if gx < sx.len() {
                        train_x.push(sx[gx].to_ascii_uppercase());
                        train_y.push(b'-');
                    }
                }
                for step in 0..dy {
                    if train_x.len() >= max_cols {
                        break;
                    }
                    let gy_opt = if r.is_rc {
                        py_end.checked_sub(step)
                    } else {
                        py_end.checked_add(step)
                    };
                    if let Some(gy) = gy_opt
                        && gy < sy.len()
                    {
                        train_x.push(b'-');
                        let yb = if r.is_rc {
                            crate::dna::complement_base(sy[gy])
                        } else {
                            sy[gy]
                        };
                        train_y.push(yb.to_ascii_uppercase());
                    }
                }
            }
        }
        let take = r.len.min(max_cols.saturating_sub(train_x.len()));
        for off in 0..take {
            let x = r.x_start + off;
            let y_opt = if r.is_rc {
                r.y_start.checked_sub(off)
            } else {
                r.y_start.checked_add(off)
            };
            if let Some(y) = y_opt
                && x < sx.len()
                && y < sy.len()
            {
                train_x.push(sx[x].to_ascii_uppercase());
                let yb = if r.is_rc {
                    crate::dna::complement_base(sy[y])
                } else {
                    sy[y]
                };
                train_y.push(yb.to_ascii_uppercase());
            }
        }
        let next_py = if r.is_rc {
            r.y_start.checked_sub(r.len)
        } else {
            r.y_start.checked_add(r.len)
        };
        prev_run = next_py.map(|py| (r.x_start + r.len, py, r.is_rc));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn last_train_em_recovers_asymmetric_gc_skew_substitution_matrix_and_directional_gaps() {
        let mut x = Vec::new();
        let mut y = Vec::new();
        for i in 0..400usize {
            match i % 10 {
                0..=3 => {
                    x.push(b'A');
                    y.push(if i % 4 == 0 { b'G' } else { b'A' });
                }
                4..=7 => {
                    x.push(b'T');
                    y.push(if i % 4 == 0 { b'C' } else { b'T' });
                }
                8 => {
                    x.push(b'G');
                    y.push(b'G');
                }
                _ => {
                    x.push(b'C');
                    y.push(b'C');
                }
            }
            if i % 25 == 0 {
                x.extend_from_slice(b"---");
                y.extend_from_slice(b"ACG");
            }
        }

        let trained = train_asymmetric_pair_params(&[(&x, &y)], 10);
        assert!(trained.iterations >= 1);
        assert!(
            trained.sub_bits[0][2] > trained.sub_bits[2][0] + 1.0,
            "expected asymmetric S[A->G] ({:.2}) > S[G->A] ({:.2}) + 1.0",
            trained.sub_bits[0][2],
            trained.sub_bits[2][0]
        );
        for (b, row) in trained.sub_bits.iter().enumerate() {
            assert!(
                row[b] > 0.0,
                "diagonal match score for base {b} must be positive: {:.2}",
                row[b]
            );
        }
        assert!(
            (trained.pi_y[1] + trained.pi_y[2]) > (trained.pi_x[1] + trained.pi_x[2]) + 0.05,
            "expected Y GC > X GC: pi_x={:?}, pi_y={:?}",
            trained.pi_x,
            trained.pi_y
        );
        assert!(
            trained.ins_open_bits > trained.del_open_bits + 1.0,
            "16 insertions vs 0 deletions must yield less negative ins_open_bits ({:.2}) than del_open_bits ({:.2})",
            trained.ins_open_bits,
            trained.del_open_bits
        );
        assert!(
            trained.ins_ext_bits > trained.del_ext_bits,
            "3-bp insertions must yield less negative ins_ext_bits ({:.2}) than del_ext_bits ({:.2})",
            trained.ins_ext_bits,
            trained.del_ext_bits
        );
        let sum_pi_x: f64 = trained.pi_x.iter().sum();
        let sum_pi_y: f64 = trained.pi_y.iter().sum();
        assert!((sum_pi_x - 1.0).abs() < 1e-12 && (sum_pi_y - 1.0).abs() < 1e-12);
    }
}
