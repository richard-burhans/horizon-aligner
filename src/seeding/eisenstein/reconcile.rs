//! Exact two-sided symmetric-difference reconciliation via linear power-sum sketches over `F_q`.

use super::{EisensteinHasher, MERSENNE_61, NEG_ONE, fq_add, fq_inv, fq_mul, fq_neg, fq_sub};
use crate::splitmix64;
use std::collections::HashSet;

fn poly_trim(p: &mut Vec<u64>) {
    while p.last() == Some(&0) {
        p.pop();
    }
}

fn poly_mul_mod(a: &[u64], b: &[u64], modulus: &[u64]) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut prod = vec![0u64; a.len() + b.len() - 1];
    for (i, &ai) in a.iter().enumerate() {
        if ai == 0 {
            continue;
        }
        for (j, &bj) in b.iter().enumerate() {
            if bj != 0 {
                prod[i + j] = fq_add(prod[i + j], fq_mul(ai, bj));
            }
        }
    }
    poly_rem(&mut prod, modulus);
    prod
}

fn poly_rem(a: &mut Vec<u64>, modulus: &[u64]) {
    poly_trim(a);
    let mdeg = modulus.len().saturating_sub(1);
    if modulus.is_empty() || a.len() <= mdeg {
        return;
    }
    let inv_lead = fq_inv(modulus[mdeg]);
    while a.len() > mdeg {
        let shift = a.len() - 1 - mdeg;
        let factor = fq_mul(*a.last().unwrap_or(&0), inv_lead);
        for (k, &mk) in modulus.iter().enumerate() {
            a[shift + k] = fq_sub(a[shift + k], fq_mul(factor, mk));
        }
        poly_trim(a);
    }
}

fn poly_div_exact(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut rem = a.to_vec();
    poly_trim(&mut rem);
    let bdeg = b.len().saturating_sub(1);
    if b.is_empty() || rem.len() <= bdeg {
        return Vec::new();
    }
    let inv_lead = fq_inv(b[bdeg]);
    let mut q = vec![0u64; rem.len() - bdeg];
    while rem.len() > bdeg {
        let shift = rem.len() - 1 - bdeg;
        let factor = fq_mul(*rem.last().unwrap_or(&0), inv_lead);
        q[shift] = factor;
        for (k, &bk) in b.iter().enumerate() {
            rem[shift + k] = fq_sub(rem[shift + k], fq_mul(factor, bk));
        }
        poly_trim(&mut rem);
    }
    poly_trim(&mut q);
    q
}

fn poly_monic(mut p: Vec<u64>) -> Vec<u64> {
    poly_trim(&mut p);
    if let Some(&lead) = p.last()
        && lead != 1
    {
        let inv = fq_inv(lead);
        for c in &mut p {
            *c = fq_mul(*c, inv);
        }
    }
    p
}

fn poly_gcd(mut a: Vec<u64>, mut b: Vec<u64>) -> Vec<u64> {
    poly_trim(&mut a);
    poly_trim(&mut b);
    while !b.is_empty() {
        poly_rem(&mut a, &b);
        std::mem::swap(&mut a, &mut b);
    }
    poly_monic(a)
}

fn poly_pow_mod(mut base: Vec<u64>, mut exp: u64, modulus: &[u64]) -> Vec<u64> {
    let mut res = vec![1u64];
    poly_rem(&mut base, modulus);
    while exp > 0 {
        if exp & 1 == 1 {
            res = poly_mul_mod(&res, &base, modulus);
        }
        base = poly_mul_mod(&base, &base, modulus);
        exp >>= 1;
    }
    res
}

/// Berlekamp-Massey over `F_q`: returns minimal LFSR connection polynomial `C(z) = 1 + C_1 z + ... + C_L z^L`.
fn berlekamp_massey(seq: &[u64]) -> Vec<u64> {
    let mut conn = vec![1u64];
    let mut prev_conn = vec![1u64];
    let mut lfsr_len = 0usize;
    let mut step_shift = 1usize;
    let mut prev_disc = 1u64;

    for n in 0..seq.len() {
        let mut disc = seq[n];
        for i in 1..=lfsr_len {
            if i < conn.len() {
                disc = fq_add(disc, fq_mul(conn[i], seq[n - i]));
            }
        }
        if disc == 0 {
            step_shift += 1;
        } else {
            let coef = fq_mul(disc, fq_inv(prev_disc));
            let old_conn = conn.clone();
            if conn.len() < prev_conn.len() + step_shift {
                conn.resize(prev_conn.len() + step_shift, 0);
            }
            for (i, &bi) in prev_conn.iter().enumerate() {
                conn[i + step_shift] = fq_sub(conn[i + step_shift], fq_mul(coef, bi));
            }
            if 2 * lfsr_len <= n {
                lfsr_len = n + 1 - lfsr_len;
                prev_conn = old_conn;
                prev_disc = disc;
                step_shift = 1;
            } else {
                step_shift += 1;
            }
        }
    }
    conn.truncate(lfsr_len + 1);
    conn
}

/// Find all distinct roots in `F_q` of a monic polynomial `phi(z)` using Cantor-Zassenhaus.
/// Returns `None` if `phi` does not split completely into `deg(phi)` distinct linear factors in `F_q`.
fn cantor_zassenhaus_roots(phi: &[u64]) -> Option<Vec<u64>> {
    let deg = phi.len().saturating_sub(1);
    if deg == 0 {
        return Some(Vec::new());
    }
    // Check that all roots lie in F_q and are distinct: gcd(phi(z), z^q - z) must have degree == deg.
    let mut x_q = poly_pow_mod(vec![0, 1], MERSENNE_61, phi);
    if x_q.len() < 2 {
        x_q.resize(2, 0);
    }
    x_q[1] = fq_sub(x_q[1], 1);
    poly_trim(&mut x_q);
    let split = poly_gcd(phi.to_vec(), x_q);
    if split.len().saturating_sub(1) != deg {
        return None;
    }

    let mut stack = vec![poly_monic(phi.to_vec())];
    let mut roots = Vec::with_capacity(deg);
    let mut rng_state = 0x517C_C1B7_2722_0A95_u64;
    let half_qm1 = (MERSENNE_61 - 1) / 2;

    while let Some(f) = stack.pop() {
        let fdeg = f.len().saturating_sub(1);
        if fdeg == 0 {
            continue;
        }
        if fdeg == 1 {
            roots.push(fq_neg(f[0]));
            continue;
        }
        // Check if 0 is a root directly.
        if f[0] == 0 {
            roots.push(0);
            stack.push(f[1..].to_vec());
            continue;
        }
        let mut split_found = false;
        for _ in 0..64 {
            rng_state = splitmix64(rng_state);
            let a = rng_state & MERSENNE_61;
            if a == MERSENNE_61 {
                continue;
            }
            let mut h = poly_pow_mod(vec![a, 1], half_qm1, &f);
            if h.is_empty() {
                h.push(NEG_ONE);
            } else {
                h[0] = fq_sub(h[0], 1);
            }
            poly_trim(&mut h);
            let g = poly_gcd(f.clone(), h);
            let gdeg = g.len().saturating_sub(1);
            if gdeg > 0 && gdeg < fdeg {
                let other = poly_monic(poly_div_exact(&f, &g));
                stack.push(g);
                stack.push(other);
                split_found = true;
                break;
            }
        }
        if !split_found {
            return None;
        }
    }
    Some(roots)
}

fn poly_eval(p: &[u64], z: u64) -> u64 {
    let mut acc = 0u64;
    for &c in p.iter().rev() {
        acc = fq_add(fq_mul(acc, z), c);
    }
    acc
}

fn poly_deriv_eval(p: &[u64], z: u64) -> u64 {
    let mut acc = 0u64;
    for (idx, &c) in p.iter().enumerate().skip(1).rev() {
        let mult = (idx as u64) % MERSENNE_61;
        acc = fq_add(fq_mul(acc, z), fq_mul(c, mult));
    }
    acc
}

/// Linear degree-`D` power-sum sketch `P_j(A) = Σ_{h ∈ A} h^j (mod q)` over `F_{2^61 - 1}`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PowerSumSketch {
    pub k: usize,
    pub count: u64,
    pub sums: Vec<u64>,
}

impl PowerSumSketch {
    /// Build a degree-`degree` power-sum sketch from a slice of distinct Eisenstein hashes.
    #[must_use]
    pub fn from_unique_hashes(k: usize, degree: usize, hashes: &[u64]) -> Self {
        let mut sums = vec![0u64; degree];
        for &h in hashes {
            let h_red = if h >= MERSENNE_61 { h % MERSENNE_61 } else { h };
            let mut pow = h_red;
            for s in &mut sums {
                *s = fq_add(*s, pow);
                pow = fq_mul(pow, h_red);
            }
        }
        Self {
            k,
            count: hashes.len() as u64,
            sums,
        }
    }

    /// Build a degree-`degree` power-sum sketch directly from a DNA sequence (deduplicating k-mer hashes).
    #[must_use]
    pub fn from_sequence(hasher: &EisensteinHasher, seq: &[u8], k: usize, degree: usize) -> Self {
        let mut unique = HashSet::<u64>::new();
        hasher.for_each_kmer_hash(seq, k, |_, h| {
            unique.insert(h);
        });
        let vals: Vec<u64> = unique.into_iter().collect();
        Self::from_unique_hashes(k, degree, &vals)
    }

    /// Linearly add `other` to `self` modulo `q` (for disjoint k-mer sets).
    #[must_use]
    pub fn add(&self, other: &Self) -> Option<Self> {
        if self.k != other.k || self.sums.len() != other.sums.len() {
            return None;
        }
        let sums = self
            .sums
            .iter()
            .zip(&other.sums)
            .map(|(&a, &b)| fq_add(a, b))
            .collect();
        Some(Self {
            k: self.k,
            count: fq_add(self.count % MERSENNE_61, other.count % MERSENNE_61),
            sums,
        })
    }

    /// Linearly subtract `other` from `self` modulo `q` (encoding `A \ B` vs `B \ A`).
    #[must_use]
    pub fn sub(&self, other: &Self) -> Option<Self> {
        if self.k != other.k || self.sums.len() != other.sums.len() {
            return None;
        }
        let sums = self
            .sums
            .iter()
            .zip(&other.sums)
            .map(|(&a, &b)| fq_sub(a, b))
            .collect();
        Some(Self {
            k: self.k,
            count: fq_sub(self.count % MERSENNE_61, other.count % MERSENNE_61),
            sums,
        })
    }

    /// Decode exact two-sided symmetric difference `(only_in_self, only_in_other)` when
    /// `2 * (|A \ B| + |B \ A|) <= degree`. Fail-closed: verifies every decoded root has
    /// multiplicity `+1` or `-1` and reproduces all `degree` power sums exactly.
    #[must_use]
    pub fn decode_symmetric_diff(&self, other: &Self) -> Option<(Vec<u64>, Vec<u64>)> {
        let diff = self.sub(other)?;
        if diff.sums.is_empty() {
            return None;
        }
        if diff.sums.iter().all(|&s| s == 0) {
            if self.count == other.count {
                return Some((Vec::new(), Vec::new()));
            }
            if self.count.checked_sub(other.count) == Some(1) {
                return Some((vec![0], Vec::new()));
            }
            if other.count.checked_sub(self.count) == Some(1) {
                return Some((Vec::new(), vec![0]));
            }
            return None;
        }

        // Berlekamp-Massey gives C(z) = 1 + C_1 z + ... + C_L z^L; reciprocal Φ(z) = Σ C_{L-i} z^i.
        let conn = berlekamp_massey(&diff.sums);
        let l = conn.len().saturating_sub(1);
        if l == 0 || 2 * l > diff.sums.len() {
            return None;
        }
        let phi: Vec<u64> = conn.iter().rev().copied().collect();
        let roots = cantor_zassenhaus_roots(&phi)?;
        if roots.len() != l {
            return None;
        }

        let mut only_self = Vec::new();
        let mut only_other = Vec::new();
        for &r in &roots {
            if r == 0 || poly_eval(&phi, r) != 0 {
                return None;
            }
            let q_i = poly_div_exact(&phi, &[fq_neg(r), 1]);
            let mut numer = 0u64;
            for (m, &qm) in q_i.iter().enumerate() {
                numer = fq_add(numer, fq_mul(qm, diff.sums[m]));
            }
            let denom = fq_mul(r, poly_deriv_eval(&phi, r));
            if denom == 0 {
                return None;
            }
            let mult = fq_mul(numer, fq_inv(denom));
            if mult == 1 {
                only_self.push(r);
            } else if mult == NEG_ONE {
                only_other.push(r);
            } else {
                return None;
            }
        }

        // Account for potential zero hash if count delta differs by ±1.
        let expected_delta = fq_sub(
            only_self.len() as u64 % MERSENNE_61,
            only_other.len() as u64 % MERSENNE_61,
        );
        if expected_delta != diff.count {
            if fq_add(expected_delta, 1) == diff.count {
                only_self.push(0);
            } else if fq_sub(expected_delta, 1) == diff.count {
                only_other.push(0);
            } else {
                return None;
            }
        }

        // Fail-closed verification: recompute all D power sums and verify exact equality.
        let check_a = Self::from_unique_hashes(self.k, diff.sums.len(), &only_self);
        let check_b = Self::from_unique_hashes(self.k, diff.sums.len(), &only_other);
        let check_diff = check_a.sub(&check_b)?;
        if check_diff.sums != diff.sums {
            return None;
        }
        only_self.sort_unstable();
        only_other.sort_unstable();
        Some((only_self, only_other))
    }
}
