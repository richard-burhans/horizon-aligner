//! `SQUINT`-inspired first-principles locus alignment models:
//! 1. Continuous-time HKY85/F81 Markov substitution kernels `P(τ) = exp(τ Q)`.
//! 2. Causal Fixed-Share repeat-null baseline `M_0` codelength estimator.
//! 3. Codon-aware `+1` / `-1` nt frameshift detection and repair.

use crate::bocpd::base_to_idx;
use std::collections::HashMap;

/// Continuous-time HKY85 / F81 DNA substitution model with equilibrium base frequencies
/// `pi = [π_A, π_C, π_G, π_T]` and transition/transversion ratio `kappa` (`kappa = 1.0` -> F81).
#[derive(Debug, Clone, Copy)]
pub struct MarkovKernel {
    pub pi: [f64; 4],
    pub kappa: f64,
}

impl MarkovKernel {
    /// Construct an HKY85 Markov kernel, normalizing `pi` to sum to 1.0.
    #[must_use]
    pub fn new(pi: [f64; 4], kappa: f64) -> Self {
        let mut clamped = [0.0; 4];
        let mut sum = 0.0;
        for (dst, &v) in clamped.iter_mut().zip(&pi) {
            let safe = if v.is_finite() && v > 0.0 { v } else { 0.25 };
            *dst = safe;
            sum += safe;
        }
        for v in &mut clamped {
            *v /= sum;
        }
        let kappa = if kappa.is_finite() && kappa > 0.0 {
            kappa
        } else {
            1.0
        };
        Self { pi: clamped, kappa }
    }

    /// Exact closed-form continuous-time transition matrix `P(τ) = exp(τ Q)` (`4 × 4`,
    /// row-stochastic, satisfying `P(τ_1) P(τ_2) = P(τ_1 + τ_2)` and detailed balance).
    /// Base indexing: `0 = A`, `1 = C`, `2 = G`, `3 = T`.
    #[must_use]
    pub fn transition_matrix(&self, tau: f64) -> [[f64; 4]; 4] {
        let t = tau.max(0.0);
        let [pi_a, pi_c, pi_g, pi_t] = self.pi;
        let pi_r = pi_a + pi_g;
        let pi_y = pi_c + pi_t;
        let rate = 2.0 * (pi_r * pi_y + self.kappa * (pi_a * pi_g + pi_c * pi_t));
        let beta = if rate > 1e-15 { 1.0 / rate } else { 1.0 };

        let e2 = (-beta * t).exp();
        let er = (-beta * t * (1.0 + pi_r * (self.kappa - 1.0))).exp();
        let ey = (-beta * t * (1.0 + pi_y * (self.kappa - 1.0))).exp();

        let is_purine = |b: usize| b == 0 || b == 2;
        let mut mat = [[0.0_f64; 4]; 4];
        for (i, row) in mat.iter_mut().enumerate() {
            let (pi_c_class, ec) = if is_purine(i) { (pi_r, er) } else { (pi_y, ey) };
            for (j, cell) in row.iter_mut().enumerate() {
                if is_purine(i) != is_purine(j) {
                    *cell = self.pi[j] * (1.0 - e2);
                } else if i == j {
                    *cell = self.pi[j]
                        + self.pi[j] * ((1.0 - pi_c_class) / pi_c_class) * e2
                        + ((pi_c_class - self.pi[j]) / pi_c_class) * ec;
                } else {
                    *cell = self.pi[j]
                        * (1.0 + ((1.0 - pi_c_class) / pi_c_class) * e2 - (1.0 / pi_c_class) * ec);
                }
            }
        }
        mat
    }
}

/// Compute per-base causal Fixed-Share repeat-null (`M_0`) codelengths (bits) for `seq`.
///
/// Mixes a Dirichlet-smoothed background expert with a causal local repeat-copy expert
/// (`context_order`-mer suffix predictor) using Herbster & Warmuth Fixed-Share switching rate
/// `alpha ∈ (0, 1)`, so tandem/low-complexity repeats compress cleanly under `M_0`.
#[must_use]
pub fn fixed_share_null_bits(seq: &[u8], alpha: f64, context_order: usize) -> Vec<f64> {
    let n = seq.len();
    if n == 0 {
        return Vec::new();
    }
    let alpha = if alpha.is_finite() { alpha } else { 0.05 }.clamp(1e-4, 0.5);
    let order = context_order.clamp(1, 12);
    let mask: u64 = (1u64 << (2 * order)) - 1;

    let mut w_bg = 0.5_f64;
    let mut w_rep = 0.5_f64;
    let mut bg_counts = [0.5_f64; 4];
    let mut bg_tot = 2.0_f64;
    let mut ctx_counts: HashMap<u64, [f64; 4]> = HashMap::new();
    let mut ctx_key = 0u64;
    let mut valid_ctx = 0usize;
    let mut bits = Vec::with_capacity(n);

    for &raw_b in seq {
        let Some(b) = base_to_idx(raw_b) else {
            bits.push(2.0);
            valid_ctx = 0;
            continue;
        };
        let p_bg = bg_counts[b] / bg_tot;
        let p_rep = if valid_ctx >= order {
            ctx_counts.get(&ctx_key).map_or(0.25, |c| {
                let tot: f64 = c.iter().sum::<f64>() + 2.0;
                (c[b] + 0.5) / tot
            })
        } else {
            0.25
        };

        let p_mix = (w_bg * p_bg + w_rep * p_rep).clamp(1e-12, 1.0);
        bits.push(-p_mix.log2());

        // Causal posterior + Fixed-Share weight sharing step.
        let post_bg = (w_bg * p_bg) / p_mix;
        let post_rep = (w_rep * p_rep) / p_mix;
        w_bg = (1.0 - alpha) * post_bg + alpha * post_rep;
        w_rep = (1.0 - alpha) * post_rep + alpha * post_bg;

        // Strictly causal update AFTER predicting current base.
        bg_counts[b] += 1.0;
        bg_tot += 1.0;
        if valid_ctx >= order {
            let entry = ctx_counts.entry(ctx_key).or_insert([0.0; 4]);
            entry[b] += 1.0;
        }
        ctx_key = ((ctx_key << 2) | (b as u64)) & mask;
        valid_ctx += 1;
    }
    bits
}

/// A detected coding frameshift event (`+1` insertion or `-1` deletion in query relative to reference frame).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FrameshiftEvent {
    pub ref_nt_pos: usize,
    pub query_nt_pos: usize,
    /// `+1` (inserted nucleotide in query) or `-1` (deleted nucleotide in query).
    pub shift_nt: i8,
}

/// Output of [`repair_coding_frameshifts`].
#[derive(Debug, Clone, PartialEq)]
pub struct FrameshiftRepairResult {
    pub aligned_ref: Vec<u8>,
    pub aligned_query: Vec<u8>,
    pub repaired_query_cds: Vec<u8>,
    pub events: Vec<FrameshiftEvent>,
    pub score_bits: f64,
}

fn translate_codon(codon: &[u8]) -> u8 {
    // Index order: A=0, C=1, G=2, T=3 (64 entries).
    const TABLE: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
    if codon.len() < 3 {
        return b'X';
    }
    let (Some(a), Some(b), Some(c)) = (
        base_to_idx(codon[0]),
        base_to_idx(codon[1]),
        base_to_idx(codon[2]),
    ) else {
        return b'X';
    };
    TABLE[(a << 4) | (b << 2) | c]
}

fn score_codon_pair(ref_c: &[u8], qry_c: &[u8]) -> f64 {
    let aa_r = translate_codon(ref_c);
    let aa_q = translate_codon(qry_c);
    if aa_r == b'X' || aa_q == b'X' {
        return 0.0;
    }
    let nt_matches = ref_c
        .iter()
        .zip(qry_c)
        .filter(|&(&a, &b)| base_to_idx(a).is_some() && base_to_idx(a) == base_to_idx(b))
        .count();
    if nt_matches == 3 {
        7.0
    } else if aa_r == aa_q {
        5.5
    } else if aa_q == b'*' && aa_r != b'*' {
        -6.0
    } else {
        match nt_matches {
            2 => 1.5,
            1 => -1.5,
            _ => -4.0,
        }
    }
}

/// Best 2-of-3 match between a 3-nt reference codon and a 2-nt deleted query codon.
/// Returns `(score, repaired_3nt_codon, deleted_offset_in_codon)`.
fn score_del1_codon(ref_c: &[u8], qry_2: &[u8]) -> (f64, [u8; 3], usize) {
    let acgt_eq = |a: u8, b: u8| -> usize {
        usize::from(base_to_idx(a).is_some() && base_to_idx(a) == base_to_idx(b))
    };
    let candidates: [([u8; 3], usize, usize); 3] = [
        // Missing base 0, 1, or 2 in query
        (
            [ref_c[0], qry_2[0], qry_2[1]],
            acgt_eq(ref_c[1], qry_2[0]) + acgt_eq(ref_c[2], qry_2[1]),
            0,
        ),
        (
            [qry_2[0], ref_c[1], qry_2[1]],
            acgt_eq(ref_c[0], qry_2[0]) + acgt_eq(ref_c[2], qry_2[1]),
            1,
        ),
        (
            [qry_2[0], qry_2[1], ref_c[2]],
            acgt_eq(ref_c[0], qry_2[0]) + acgt_eq(ref_c[1], qry_2[1]),
            2,
        ),
    ];
    let mut best = candidates[0];
    for cand in candidates.into_iter().skip(1) {
        if cand.1 > best.1 {
            best = cand;
        }
    }
    let sc = match best.1 {
        2 => 3.5,
        1 => 0.0,
        _ => -3.0,
    };
    (sc, best.0, best.2)
}

/// Align `query_cds` against an in-frame coding `ref_cds` using codon-aware amino-acid scoring
/// with explicit `+1` / `-1` nt frameshift transitions, returning detected frameshifts and an
/// in-frame repaired query CDS.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn repair_coding_frameshifts(
    ref_cds: &[u8],
    query_cds: &[u8],
    fs_penalty_bits: f64,
) -> FrameshiftRepairResult {
    const MAX_FS_DP_CELLS: usize = 16_000_000;
    let n_codons = ref_cds.len() / 3;
    let m = query_cds.len();
    let Some(cells) = (n_codons + 1).checked_mul(m + 1) else {
        return FrameshiftRepairResult {
            aligned_ref: Vec::new(),
            aligned_query: Vec::new(),
            repaired_query_cds: Vec::new(),
            events: Vec::new(),
            score_bits: 0.0,
        };
    };
    if n_codons == 0 || m == 0 || cells > MAX_FS_DP_CELLS {
        return FrameshiftRepairResult {
            aligned_ref: Vec::new(),
            aligned_query: Vec::new(),
            repaired_query_cds: Vec::new(),
            events: Vec::new(),
            score_bits: 0.0,
        };
    }

    let fs_pen = if fs_penalty_bits.is_finite() {
        fs_penalty_bits.abs()
    } else {
        8.0
    }
    .max(1.0);
    let codon_gap_pen = 5.5_f64;
    let neg_inf = -1.0e15_f64;
    let stride = m + 1;
    let mut dp = vec![neg_inf; cells];
    // Ops: 1=codon(3,3), 2=del_codon(3,0), 3=ins_codon(0,3), 4=fs_ins1(0,1), 5=fs_del1(3,2)
    let mut bt = vec![0u8; cells];
    dp[0] = 0.0;

    for c in 0..=n_codons {
        let ref_c = if c > 0 {
            &ref_cds[(c - 1) * 3..c * 3]
        } else {
            &[][..]
        };
        for j in 0..=m {
            if c == 0 && j == 0 {
                continue;
            }
            let mut best = neg_inf;
            let mut op = 0u8;

            if c > 0 && j >= 3 {
                let sc =
                    dp[(c - 1) * stride + (j - 3)] + score_codon_pair(ref_c, &query_cds[j - 3..j]);
                if sc > best {
                    best = sc;
                    op = 1;
                }
            }
            if c > 0 {
                let sc = dp[(c - 1) * stride + j] - codon_gap_pen;
                if sc > best {
                    best = sc;
                    op = 2;
                }
            }
            if j >= 3 {
                let sc = dp[c * stride + (j - 3)] - codon_gap_pen;
                if sc > best {
                    best = sc;
                    op = 3;
                }
            }
            if j >= 1 {
                let sc = dp[c * stride + (j - 1)] - fs_pen;
                if sc > best {
                    best = sc;
                    op = 4;
                }
            }
            if c > 0 && j >= 2 {
                let (d_sc, _, _) = score_del1_codon(ref_c, &query_cds[j - 2..j]);
                let sc = dp[(c - 1) * stride + (j - 2)] + d_sc - fs_pen;
                if sc > best {
                    best = sc;
                    op = 5;
                }
            }
            dp[c * stride + j] = best;
            bt[c * stride + j] = op;
        }
    }

    let score_bits = dp[n_codons * stride + m];
    let mut ops = Vec::new();
    let mut c = n_codons;
    let mut j = m;
    while c > 0 || j > 0 {
        let op = bt[c * stride + j];
        ops.push((op, c, j));
        match op {
            1 => {
                c -= 1;
                j -= 3;
            }
            2 => c -= 1,
            3 => j -= 3,
            4 => j -= 1,
            5 => {
                c -= 1;
                j -= 2;
            }
            _ => break,
        }
    }
    ops.reverse();

    let mut aligned_ref = Vec::new();
    let mut aligned_query = Vec::new();
    let mut repaired_query_cds = Vec::new();
    let mut events = Vec::new();

    for (op, c_end, j_end) in ops {
        match op {
            1 => {
                let rc = &ref_cds[(c_end - 1) * 3..c_end * 3];
                let qc = &query_cds[j_end - 3..j_end];
                aligned_ref.extend_from_slice(rc);
                aligned_query.extend_from_slice(qc);
                repaired_query_cds.extend_from_slice(qc);
            }
            2 => {
                let rc = &ref_cds[(c_end - 1) * 3..c_end * 3];
                aligned_ref.extend_from_slice(rc);
                aligned_query.extend_from_slice(b"---");
            }
            3 => {
                let qc = &query_cds[j_end - 3..j_end];
                aligned_ref.extend_from_slice(b"---");
                aligned_query.extend_from_slice(qc);
                repaired_query_cds.extend_from_slice(qc);
            }
            4 => {
                // +1 nt frameshift insertion in query: omit from repaired in-frame CDS.
                aligned_ref.push(b'-');
                aligned_query.push(query_cds[j_end - 1]);
                events.push(FrameshiftEvent {
                    ref_nt_pos: c_end * 3,
                    query_nt_pos: j_end - 1,
                    shift_nt: 1,
                });
            }
            5 => {
                // -1 nt frameshift deletion in query: emit '-' in aligned_query and restore 3-nt
                // codon frame in repaired_query_cds.
                let rc = &ref_cds[(c_end - 1) * 3..c_end * 3];
                let q2 = &query_cds[j_end - 2..j_end];
                let (_, repaired, del_off) = score_del1_codon(rc, q2);
                let mut q_aln = repaired;
                q_aln[del_off] = b'-';
                aligned_ref.extend_from_slice(rc);
                aligned_query.extend_from_slice(&q_aln);
                repaired_query_cds.extend_from_slice(&repaired);
                events.push(FrameshiftEvent {
                    ref_nt_pos: (c_end - 1) * 3,
                    query_nt_pos: j_end - 2,
                    shift_nt: -1,
                });
            }
            _ => {}
        }
    }

    FrameshiftRepairResult {
        aligned_ref,
        aligned_query,
        repaired_query_cds,
        events,
        score_bits,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn markov_kernel_satisfies_chapman_kolmogorov_stochasticity_and_detailed_balance() {
        let model = MarkovKernel::new([0.35, 0.15, 0.15, 0.35], 3.0);
        let p1 = model.transition_matrix(0.15);
        let p2 = model.transition_matrix(0.25);
        let p12 = model.transition_matrix(0.40);

        for i in 0..4 {
            let row_sum: f64 = p1[i].iter().sum();
            assert!((row_sum - 1.0).abs() < 1e-12, "row {i} sum != 1: {row_sum}");
            for j in 0..4 {
                let lhs = model.pi[i] * p1[i][j];
                let rhs = model.pi[j] * p1[j][i];
                assert!(
                    (lhs - rhs).abs() < 1e-12,
                    "detailed balance violated at ({i},{j})"
                );
                let prod: f64 = (0..4).map(|k| p1[i][k] * p2[k][j]).sum();
                assert!(
                    (prod - p12[i][j]).abs() < 1e-10,
                    "Chapman-Kolmogorov failed at ({i},{j}): {prod} vs {}",
                    p12[i][j]
                );
            }
        }
    }

    #[test]
    fn fixed_share_null_compresses_tandem_repeats_far_below_random_dna() {
        let repeat = b"ACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTTACGTT";
        let mut st = 0x5555_AAAA_1111_2222_u64;
        let unique: Vec<u8> = (0..repeat.len())
            .map(|_| {
                st = crate::splitmix64(st);
                b"ACGT"[(st & 3) as usize]
            })
            .collect();
        let rep_bits: f64 = fixed_share_null_bits(repeat, 0.05, 4).iter().sum();
        let uniq_bits: f64 = fixed_share_null_bits(&unique, 0.05, 4).iter().sum();
        assert!(
            rep_bits < uniq_bits * 0.65,
            "Fixed-Share M0 must compress tandem repeat ({rep_bits:.1} bits) far below random ({uniq_bits:.1} bits)"
        );
    }

    #[test]
    fn repair_coding_frameshifts_detects_and_repairs_plus1_and_minus1_frameshifts() {
        let ref_cds = b"ATGGCTGAACTGAAGCGTACCTTCGACGGTATCCTGGCTAAGCACTGA";
        assert_eq!(ref_cds.len() % 3, 0);

        let mut q_ins = Vec::new();
        q_ins.extend_from_slice(&ref_cds[..18]);
        q_ins.push(b'T');
        q_ins.extend_from_slice(&ref_cds[18..]);

        let res = repair_coding_frameshifts(ref_cds, &q_ins, 8.0);
        assert_eq!(
            res.events.len(),
            1,
            "expected exactly 1 detected +1 frameshift event, got {:?}",
            res.events
        );
        assert_eq!(res.events[0].shift_nt, 1);
        assert_eq!(res.repaired_query_cds.len() % 3, 0);
        assert_eq!(res.repaired_query_cds, ref_cds);

        let mut q_del = Vec::new();
        q_del.extend_from_slice(&ref_cds[..24]);
        q_del.extend_from_slice(&ref_cds[25..]);
        let res_del = repair_coding_frameshifts(ref_cds, &q_del, 8.0);
        assert_eq!(res_del.events.len(), 1);
        assert_eq!(res_del.events[0].shift_nt, -1);
        assert_eq!(res_del.repaired_query_cds.len() % 3, 0);
        let ug_ref: Vec<u8> = res_del
            .aligned_ref
            .iter()
            .copied()
            .filter(|&b| b != b'-')
            .collect();
        let ug_qry: Vec<u8> = res_del
            .aligned_query
            .iter()
            .copied()
            .filter(|&b| b != b'-')
            .collect();
        assert_eq!(ug_ref, ref_cds);
        assert_eq!(ug_qry, q_del, "aligned_query must ungap back to query_cds");
        assert_eq!(score_codon_pair(b"NNN", b"NNN"), 0.0);
    }
}
