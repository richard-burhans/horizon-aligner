//! Dependency-free DNA alphabet encoding and Watson–Crick complement primitives.

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

/// Fast branchless lookup mapping `A/a -> 0`, `C/c -> 1`, `G/g -> 2`, `T/t -> 3`, and any other byte -> `4`.
#[inline]
#[must_use]
pub const fn base_to_idx_fast(b: u8) -> u8 {
    BASE_IDX_LUT[b as usize]
}

/// Map canonical ASCII DNA bases (`A/C/G/T`, case-insensitive) to `Some(0..=3)` and non-ACGT bytes to `None`.
#[inline]
#[must_use]
pub const fn base_to_idx(b: u8) -> Option<usize> {
    let idx = BASE_IDX_LUT[b as usize];
    if idx < 4 { Some(idx as usize) } else { None }
}

/// Complement a single DNA base (preserving case for `A`/`C`/`G`/`T`).
#[inline]
#[must_use]
pub const fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        other => other,
    }
}

/// Reverse-complement a DNA byte slice while preserving ASCII case (needed for softmasked seeding).
#[must_use]
pub fn revcomp_preserve_case(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_base(b)).collect()
}

/// Reverse-complement a DNA byte slice (uppercasing ASCII bases).
#[must_use]
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| complement_base(b).to_ascii_uppercase())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_alphabet_and_revcomp_work() {
        assert_eq!(base_to_idx(b'A'), Some(0));
        assert_eq!(base_to_idx(b'c'), Some(1));
        assert_eq!(base_to_idx(b'G'), Some(2));
        assert_eq!(base_to_idx(b't'), Some(3));
        assert_eq!(base_to_idx(b'N'), None);
        assert_eq!(base_to_idx_fast(b'N'), 4);
        assert_eq!(revcomp_preserve_case(b"AcGtN"), b"NaCgT");
        assert_eq!(revcomp(b"AcGtN"), b"NACGT");
    }
}
