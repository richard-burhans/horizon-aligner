//! Provably strand-symmetric Eisenstein polynomial k-mer/spaced-seed hashing over
//! `F_q` (`q = 2^61 - 1`), linear power-sum set sketches with exact symmetric-difference
//! reconciliation, and O(N*H) auto-hinge clade selection.

use crate::AlignError;
use crate::splitmix64;
mod reconcile;
mod sketch;

pub use reconcile::PowerSumSketch;
pub use sketch::{AutoHingePlan, EisensteinSamplerSketch, select_auto_hinges_eisenstein};

/// Mersenne prime `q = 2^61 - 1`. Satisfies `q ≡ 1 (mod 3)`, so a primitive cube root of
/// unity `ω` lies directly in `F_q`.
pub const MERSENNE_61: u64 = (1u64 << 61) - 1;

#[inline]
const fn fq_reduce_u128(z: u128) -> u64 {
    let q = MERSENNE_61 as u128;
    #[allow(clippy::cast_possible_truncation)]
    let lo = (z & q) as u64;
    #[allow(clippy::cast_possible_truncation)]
    let hi = (z >> 61) as u64;
    let s = lo + hi;
    let s = (s & MERSENNE_61) + (s >> 61);
    if s >= MERSENNE_61 { s - MERSENNE_61 } else { s }
}

#[inline]
pub(super) const fn fq_add(a: u64, b: u64) -> u64 {
    let s = a + b;
    if s >= MERSENNE_61 { s - MERSENNE_61 } else { s }
}

#[inline]
pub(super) const fn fq_sub(a: u64, b: u64) -> u64 {
    if a >= b { a - b } else { MERSENNE_61 - (b - a) }
}

#[inline]
pub(super) const fn fq_neg(a: u64) -> u64 {
    if a == 0 { 0 } else { MERSENNE_61 - a }
}

#[inline]
pub(super) const fn fq_mul(a: u64, b: u64) -> u64 {
    fq_reduce_u128((a as u128) * (b as u128))
}

pub(super) const fn fq_pow(mut base: u64, mut exp: u64) -> u64 {
    let mut acc = 1u64;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = fq_mul(acc, base);
        }
        base = fq_mul(base, base);
        exp >>= 1;
    }
    acc
}

#[inline]
pub(super) const fn fq_inv(a: u64) -> u64 {
    fq_pow(a, MERSENNE_61 - 2)
}

const fn find_primitive_cube_root() -> u64 {
    let exp = (MERSENNE_61 - 1) / 3;
    let mut g = 2u64;
    while g < 1_000 {
        let cand = fq_pow(g, exp);
        if cand > 1 && cand < MERSENNE_61 && fq_pow(cand, 3) == 1 {
            return cand;
        }
        g += 1;
    }
    panic!("no primitive cube root found in F_q")
}

/// Primitive cube root of unity `ω ∈ F_q` (`ω^3 = 1`, `ω^2 + ω + 1 = 0`).
pub const EISENSTEIN_OMEGA: u64 = find_primitive_cube_root();

pub(super) const NEG_ONE: u64 = MERSENNE_61 - 1;
const NEG_OMEGA: u64 = MERSENNE_61 - EISENSTEIN_OMEGA;

#[inline]
const fn eisenstein_base(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(1),
        b'T' | b't' => Some(NEG_ONE),
        b'C' | b'c' => Some(EISENSTEIN_OMEGA),
        b'G' | b'g' => Some(NEG_OMEGA),
        _ => None,
    }
}

/// Strand-symmetric polynomial k-mer and spaced-seed hasher over `F_{2^61 - 1}` using the
/// sixth-root Eisenstein map (`A = 1, T = -1, C = ω, G = -ω`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EisensteinHasher {
    pub x: u64,
    pub lambda: u64,
}

impl Default for EisensteinHasher {
    fn default() -> Self {
        Self::from_seed(0x4549_5345_4E53_5445)
    }
}

impl EisensteinHasher {
    /// Construct a hasher with field evaluation point `x` and shift `λ` derived from `seed`.
    #[must_use]
    pub fn from_seed(seed: u64) -> Self {
        let mut s = splitmix64(seed);
        let mut x = s & MERSENNE_61;
        while x < 2 || x == MERSENNE_61 - 1 || x == MERSENNE_61 {
            s = splitmix64(s);
            x = s & MERSENNE_61;
        }
        s = splitmix64(s ^ 0x9E37_79B9_7F4A_7C15);
        let mut lambda = s & MERSENNE_61;
        while lambda == 0 || lambda == MERSENNE_61 {
            s = splitmix64(s);
            lambda = s & MERSENNE_61;
        }
        Self { x, lambda }
    }

    /// Exact Schwartz-Zippel pairwise collision upper bound `2(k - 1) / q` for distinct
    /// double-stranded k-mers.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn collision_bound(k: usize) -> f64 {
        2.0 * (k.saturating_sub(1) as f64) / (MERSENNE_61 as f64)
    }

    #[inline]
    const fn combine_strands(&self, fwd: u64, rev: u64) -> u64 {
        fq_mul(fq_add(fwd, self.lambda), fq_add(rev, self.lambda))
    }

    /// Direct O(k) strand-symmetric hash `S(s) = (F(s) + λ)(R(s) + λ) mod q` for arbitrary `k >= 1`.
    #[must_use]
    pub fn hash_kmer(&self, kmer: &[u8]) -> Option<u64> {
        if kmer.is_empty() {
            return None;
        }
        let mut fwd = 0u64;
        let mut x_pow = 1u64;
        let mut rev = 0u64;
        for &b in kmer {
            let v = eisenstein_base(b)?;
            fwd = fq_add(fwd, fq_mul(v, x_pow));
            x_pow = fq_mul(x_pow, self.x);
            rev = fq_sub(fq_mul(rev, self.x), v);
        }
        Some(self.combine_strands(fwd, rev))
    }

    /// Zero-allocation streaming O(N) strand-symmetric k-mer hash iterator invoking `on_kmer(pos, hash)`
    /// at every valid ACGT window of length `k >= 1`, bit-identical to [`Self::hash_kmer`].
    pub fn for_each_kmer_hash(&self, seq: &[u8], k: usize, mut on_kmer: impl FnMut(usize, u64)) {
        if k == 0 || seq.len() < k {
            return;
        }
        let x_inv = fq_inv(self.x);
        let x_km1 = fq_pow(self.x, (k - 1) as u64);

        let mut fwd = 0u64;
        let mut rev = 0u64;
        let mut x_pow = 1u64;
        let mut run_len = 0usize;

        for (i, &b) in seq.iter().enumerate() {
            let Some(v) = eisenstein_base(b) else {
                fwd = 0;
                rev = 0;
                x_pow = 1;
                run_len = 0;
                continue;
            };
            if run_len < k {
                fwd = fq_add(fwd, fq_mul(v, x_pow));
                x_pow = fq_mul(x_pow, self.x);
                rev = fq_sub(fq_mul(rev, self.x), v);
                run_len += 1;
                if run_len == k {
                    on_kmer(i + 1 - k, self.combine_strands(fwd, rev));
                }
            } else {
                let Some(v_out) = eisenstein_base(seq[i - k]) else {
                    fwd = 0;
                    rev = 0;
                    x_pow = 1;
                    run_len = 0;
                    continue;
                };
                fwd = fq_add(fq_mul(fq_sub(fwd, v_out), x_inv), fq_mul(v, x_km1));
                rev = fq_sub(fq_mul(fq_add(rev, fq_mul(v_out, x_km1)), self.x), v);
                on_kmer(i + 1 - k, self.combine_strands(fwd, rev));
            }
        }
    }

    /// Rolling O(N) strand-symmetric k-mer hash stream `(pos, hash)` for arbitrary `k >= 1`,
    /// bit-identical to [`Self::hash_kmer`] at every valid ACGT window.
    #[must_use]
    pub fn rolling_kmer_hashes(&self, seq: &[u8], k: usize) -> Vec<(usize, u64)> {
        if k == 0 || seq.len() < k {
            return Vec::new();
        }
        let mut out = Vec::with_capacity(seq.len() - k + 1);
        self.for_each_kmer_hash(seq, k, |pos, h| out.push((pos, h)));
        out
    }

    /// Direct strand-symmetric spaced-seed hash over an arbitrary-length mask (`1` = care, `0` = don't-care),
    /// invariant under reverse complement and matching the union of `mask` and its mirror image.
    #[must_use]
    pub fn hash_spaced_seed(&self, window: &[u8], care_offsets: &[usize]) -> Option<u64> {
        let span = window.len();
        if span == 0 || care_offsets.is_empty() {
            return None;
        }
        let mut fwd = 0u64;
        let mut rev = 0u64;
        for &t in care_offsets {
            if t >= span {
                return None;
            }
            let x_t = fq_pow(self.x, t as u64);
            let vf = eisenstein_base(window[t])?;
            let vr = eisenstein_base(window[span - 1 - t])?;
            fwd = fq_add(fwd, fq_mul(vf, x_t));
            rev = fq_sub(rev, fq_mul(vr, x_t));
        }
        Some(self.combine_strands(fwd, rev))
    }

    /// Extract all strand-symmetric spaced-seed hashes `(pos, hash)` across `seq` for arbitrary mask span.
    ///
    /// # Errors
    /// Returns [`AlignError::InvalidMask`] if `mask` is empty, contains characters other than `'0'`/`'1'`,
    /// or has zero care positions.
    pub fn spaced_seed_hashes(
        &self,
        seq: &[u8],
        mask: &str,
    ) -> Result<Vec<(usize, u64)>, AlignError> {
        let span = mask.len();
        if span == 0 {
            return Err(AlignError::InvalidMask(
                "spaced-seed mask cannot be empty".into(),
            ));
        }
        let mut care_offsets = Vec::with_capacity(span);
        for (i, b) in mask.bytes().enumerate() {
            match b {
                b'1' => care_offsets.push(i),
                b'0' => {}
                _ => {
                    return Err(AlignError::InvalidMask(format!(
                        "invalid character '{}' in spaced-seed mask",
                        b as char
                    )));
                }
            }
        }
        if care_offsets.is_empty() {
            return Err(AlignError::InvalidMask(
                "spaced-seed mask must contain at least one '1'".into(),
            ));
        }
        if seq.len() < span {
            return Ok(Vec::new());
        }
        // Precompute x^t for every care offset so evaluation per window is O(weight) multiplications.
        let x_pows: Vec<u64> = care_offsets
            .iter()
            .map(|&t| fq_pow(self.x, t as u64))
            .collect();
        let mut out = Vec::with_capacity(seq.len() - span + 1);
        for pos in 0..=(seq.len() - span) {
            let win = &seq[pos..pos + span];
            let mut fwd = 0u64;
            let mut rev = 0u64;
            let mut ok = true;
            for (&t, &xp) in care_offsets.iter().zip(&x_pows) {
                let (Some(vf), Some(vr)) =
                    (eisenstein_base(win[t]), eisenstein_base(win[span - 1 - t]))
                else {
                    ok = false;
                    break;
                };
                fwd = fq_add(fwd, fq_mul(vf, xp));
                rev = fq_sub(rev, fq_mul(vr, xp));
            }
            if ok {
                out.push((pos, self.combine_strands(fwd, rev)));
            }
        }
        Ok(out)
    }
}

// ---------- Polynomial arithmetic over F_{2^61 - 1}[z] for power-sum set reconciliation ----------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seeding::revcomp;
    use crate::splitmix64;
    use std::collections::HashSet;

    fn random_dna(seed: u64, len: usize) -> Vec<u8> {
        let bases = *b"ACGT";
        let mut st = seed;
        (0..len)
            .map(|_| {
                st = splitmix64(st);
                bases[(st & 3) as usize]
            })
            .collect()
    }

    #[test]
    fn eisenstein_omega_satisfies_cube_root_field_identities() {
        const {
            assert!(EISENSTEIN_OMEGA > 1 && EISENSTEIN_OMEGA < MERSENNE_61);
        }
        let w = u128::from(EISENSTEIN_OMEGA);
        let q = u128::from(MERSENNE_61);
        let w2 = (w * w) % q;
        let w3 = (w2 * w) % q;
        assert_eq!(w3, 1, "ω^3 must equal 1 mod 2^61-1");
        assert_eq!((w2 + w + 1) % q, 0, "ω^2 + ω + 1 must equal 0 mod 2^61-1");
    }

    #[test]
    fn eisenstein_rolling_matches_direct_is_strand_symmetric_and_handles_long_kmers() {
        let hasher = EisensteinHasher::default();
        for &k in &[1usize, 7, 31, 61, 101] {
            let mut seq = random_dna(0xBEEF + k as u64, 300);
            seq[120] = b'N';
            seq[121] = b'N';
            let rc = revcomp(&seq);

            let rolled = hasher.rolling_kmer_hashes(&seq, k);
            assert!(!rolled.is_empty(), "expected rolled hashes for k={k}");
            for &(pos, h) in &rolled {
                let direct = hasher.hash_kmer(&seq[pos..pos + k]).expect("direct");
                assert_eq!(h, direct, "rolling != direct at pos={pos}, k={k}");
                let kmer_rc = revcomp(&seq[pos..pos + k]);
                let h_rc = hasher.hash_kmer(&kmer_rc).expect("rc direct");
                assert_eq!(h, h_rc, "not strand-symmetric at pos={pos}, k={k}");
            }

            let rolled_rc = hasher.rolling_kmer_hashes(&rc, k);
            assert_eq!(rolled.len(), rolled_rc.len());
            let set_fwd: HashSet<u64> = rolled.iter().map(|&(_, h)| h).collect();
            let set_rc: HashSet<u64> = rolled_rc.iter().map(|&(_, h)| h).collect();
            assert_eq!(
                set_fwd, set_rc,
                "whole-sequence revcomp hash sets must match for k={k}"
            );
        }
    }

    #[test]
    fn eisenstein_asymmetric_long_spaced_seed_is_strand_symmetric() {
        let hasher = EisensteinHasher::default();
        // Deliberately asymmetric 45-bp spaced seed (cannot fit in 32-bp packed word!).
        let mask = "111011001110101100101110110100111011010110111";
        let seq = random_dna(0xCAFE_0045, 250);
        let rc = revcomp(&seq);

        let fwd_hashes = hasher.spaced_seed_hashes(&seq, mask).expect("fwd spaced");
        let rc_hashes = hasher.spaced_seed_hashes(&rc, mask).expect("rc spaced");
        assert_eq!(fwd_hashes.len(), seq.len() - mask.len() + 1);
        assert_eq!(fwd_hashes.len(), rc_hashes.len());
        for (&(_, hf), &(_, hr)) in fwd_hashes.iter().zip(rc_hashes.iter().rev()) {
            assert_eq!(
                hf, hr,
                "spaced-seed hash must match exact mirrored revcomp position"
            );
        }
    }

    #[test]
    fn power_sum_sketch_linearly_reconciles_exact_two_sided_kmer_differences() {
        let hasher = EisensteinHasher::default();
        let seq_a = random_dna(0x1111_2222, 400);
        let mut seq_b = seq_a.clone();
        seq_b[200] = if seq_b[200] == b'A' { b'G' } else { b'A' };

        let k = 7usize; // 1 SNP alters <= 7 k-mers per side => <= 14 total diffs
        let degree = 32usize; // 2 * 14 = 28 <= 32
        let sk_a = PowerSumSketch::from_sequence(&hasher, &seq_a, k, degree);
        let sk_b = PowerSumSketch::from_sequence(&hasher, &seq_b, k, degree);

        let kmers_a: HashSet<u64> = hasher
            .rolling_kmer_hashes(&seq_a, k)
            .into_iter()
            .map(|(_, h)| h)
            .collect();
        let kmers_b: HashSet<u64> = hasher
            .rolling_kmer_hashes(&seq_b, k)
            .into_iter()
            .map(|(_, h)| h)
            .collect();
        let mut true_only_a: Vec<u64> = kmers_a.difference(&kmers_b).copied().collect();
        let mut true_only_b: Vec<u64> = kmers_b.difference(&kmers_a).copied().collect();
        true_only_a.sort_unstable();
        true_only_b.sort_unstable();
        assert!(!true_only_a.is_empty() && !true_only_b.is_empty());

        let (mut dec_a, mut dec_b) = sk_a
            .decode_symmetric_diff(&sk_b)
            .expect("exact two-sided power-sum decoding must succeed");
        dec_a.sort_unstable();
        dec_b.sort_unstable();
        assert_eq!(dec_a, true_only_a, "decoded A \\ B mismatch");
        assert_eq!(dec_b, true_only_b, "decoded B \\ A mismatch");

        let rc_a = revcomp(&seq_a);
        let sk_rc = PowerSumSketch::from_sequence(&hasher, &rc_a, k, degree);
        assert_eq!(
            sk_a.decode_symmetric_diff(&sk_rc),
            Some((Vec::new(), Vec::new())),
            "revcomp genome must have zero symmetric difference without canonicalization"
        );
    }

    #[test]
    fn select_auto_hinges_eisenstein_separates_two_clades_and_assigns_tips_to_own_clade() {
        let hasher = EisensteinHasher::default();
        let clade_a = random_dna(0xAAAA_0001, 4_000);
        let clade_b = random_dna(0xBBBB_0002, 4_000);
        let mutate = |src: &[u8], seed: u64| -> Vec<u8> {
            let bases = *b"ACGT";
            let mut st = seed;
            src.iter()
                .map(|&b| {
                    st = splitmix64(st);
                    if st % 100 < 2 {
                        bases[((st >> 8) & 3) as usize]
                    } else {
                        b
                    }
                })
                .collect()
        };
        let a1 = mutate(&clade_a, 1);
        let a2 = mutate(&clade_a, 2);
        let a3 = revcomp(&mutate(&clade_a, 3));
        let b1 = mutate(&clade_b, 4);
        let b2 = mutate(&clade_b, 5);

        let sketches: Vec<(String, EisensteinSamplerSketch)> = [
            ("a1", a1.as_slice()),
            ("a2", a2.as_slice()),
            ("a3", a3.as_slice()),
            ("b1", b1.as_slice()),
            ("b2", b2.as_slice()),
        ]
        .into_iter()
        .map(|(name, seq)| {
            (
                name.to_owned(),
                EisensteinSamplerSketch::from_contigs(&hasher, &[seq], 31, 512),
            )
        })
        .collect();

        let (hinges, tip_map) = select_auto_hinges_eisenstein(&sketches, 2).expect("auto hinges");
        assert_eq!(hinges.len(), 2);
        assert_eq!(tip_map.len(), 3);
        let picked_clade_a = hinges.iter().any(|h| h.starts_with('a'));
        let picked_clade_b = hinges.iter().any(|h| h.starts_with('b'));
        assert!(
            picked_clade_a && picked_clade_b,
            "expected one hinge from clade a and one from clade b, got {hinges:?}"
        );
        for (tip, hinge) in &tip_map {
            assert_eq!(
                tip.as_bytes()[0],
                hinge.as_bytes()[0],
                "tip {tip} must be assigned to its own clade's hinge (got {hinge})"
            );
        }
    }

    #[test]
    fn power_sum_sketch_fails_closed_for_degree_zero_and_oversized_diff() {
        let a = PowerSumSketch::from_unique_hashes(21, 0, &[11, 22]);
        let b = PowerSumSketch::from_unique_hashes(21, 0, &[33, 44]);
        assert_eq!(
            a.decode_symmetric_diff(&b),
            None,
            "degree=0 power-sum sketch must fail closed (never claim disjoint sets are equal)"
        );

        // 4 diffs > degree/2 (degree=4 handles at most 2 diffs) -> fail closed.
        let s1 = PowerSumSketch::from_unique_hashes(21, 4, &[1, 2, 3, 4, 100, 200]);
        let s2 = PowerSumSketch::from_unique_hashes(21, 4, &[1, 2, 3, 4, 300, 400, 500]);
        assert_eq!(s1.decode_symmetric_diff(&s2), None);
    }
}
