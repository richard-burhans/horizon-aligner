//! Spaced-seed indexing, transition-tolerant and translated seeding, and LCB chaining.

use crate::AlignError;
use crate::dna::base_to_idx;
use rayon::prelude::*;
use std::collections::HashMap;

pub use translated::aa5_keys;

/// KP-14 spaced seed mask (weight 11, span 15).
pub const MASK_PRIMARY_11_15: &str = "111011011011011";
/// Palindromic weight-15, span-21 spaced seed (`4^15 / 2 = 536.9M` canonical keys) for
/// mega-chromosome (`>= 10 Mbp`) alignments where weight-11 (`2.1M` keys) saturates.
pub const MASK_LARGE_15_21: &str = "111011011010110110111";
/// Optimal non-periodic weight-11 `PatternHunter` spaced seed (Ma et al., 2002).
pub const MASK_PH11_18: &str = "111010010100110111";
/// Secondary `w=10` spaced seed mask for divergent intra-LCB chaining.
pub const MASK_W10: &str = "11011011011011";
/// Local `w=9` spaced seed mask for intra-gap anchor seeding.
pub const MASK_W9: &str = "110101101011011";
/// Transition-tolerant 1-bit `R/Y` spaced-seed mask (weight 18, span 25).
pub(crate) const MASK_RY18: &str = "1110110110111011011011011";
/// Standard `LASTZ` `12-of-19` spaced-seed mask.
pub const MASK_LASTZ12: &str = "1110100110010101111";
/// Complementary `11-of-18` spaced-seed mask for multi-shape divergent seeding.
pub const MASK_LASTZ11: &str = "111010011001010111";

#[inline]
fn is_low_complexity_window(sub: &[u8]) -> bool {
    let mut counts = [0usize; 4];
    let mut total = 0usize;
    for &b in sub {
        if let Some(idx) = base_to_idx(b) {
            counts[idx] += 1;
            total += 1;
        }
    }
    if total == 0 {
        return true;
    }
    #[allow(clippy::cast_precision_loss)]
    let inv_n = 1.0 / (total as f64);
    let mut entropy_bits = 0.0_f64;
    for &c in &counts {
        if c > 0 {
            #[allow(clippy::cast_precision_loss)]
            let p = (c as f64) * inv_n;
            entropy_bits -= p * p.log2();
        }
    }
    entropy_bits <= 1.05
}

/// Linear-time SDUST-style triplet low-complexity masker (`W in {16, 24, 32}`, `T = 2.0`,
/// Morgulis et al., 2006). Masks homopolymers, dinucleotide, trinucleotide, and tetranucleotide
/// simple repeats in `O(N)` time so seeders avoid `O(span)` per-position checks.
#[must_use]
pub fn sdust_mask(seq: &[u8]) -> Vec<bool> {
    let n = seq.len();
    let mut masked = vec![false; n];
    if n < 16 {
        return masked;
    }

    let mut triplets = vec![255_u8; n - 2];
    for (i, slot) in triplets.iter_mut().enumerate() {
        if let (Some(a), Some(b), Some(c)) = (
            base_to_idx(seq[i]),
            base_to_idx(seq[i + 1]),
            base_to_idx(seq[i + 2]),
        ) {
            *slot = u8::try_from((a << 4) | (b << 2) | c).unwrap_or(255);
        }
    }

    for &w in &[16_usize, 24_usize, 32_usize] {
        if n < w {
            continue;
        }
        let nw = w - 2;
        let min_c = 3_u16;
        let thresh = 2_i32 * i32::try_from(w - 3).unwrap_or(i32::MAX);
        let mut counts = [0_u16; 64];
        let mut sum_pairs = 0_i32;

        for &t in &triplets[..nw] {
            if t < 64 {
                let ti = t as usize;
                sum_pairs += i32::from(counts[ti]);
                counts[ti] += 1;
            }
        }
        if sum_pairs > thresh {
            for (p, &t) in triplets[..nw].iter().enumerate() {
                if t < 64 && counts[t as usize] >= min_c {
                    masked[p] = true;
                    masked[p + 1] = true;
                    masked[p + 2] = true;
                }
            }
        }

        for start in 1..=(triplets.len() - nw) {
            let out_t = triplets[start - 1];
            if out_t < 64 {
                let oi = out_t as usize;
                counts[oi] -= 1;
                sum_pairs -= i32::from(counts[oi]);
            }
            let in_t = triplets[start + nw - 1];
            if in_t < 64 {
                let ii = in_t as usize;
                sum_pairs += i32::from(counts[ii]);
                counts[ii] += 1;
            }
            if sum_pairs > thresh {
                for (offset, &t) in triplets[start..start + nw].iter().enumerate() {
                    if t < 64 && counts[t as usize] >= min_c {
                        let p = start + offset;
                        if !masked[p] || !masked[p + 1] || !masked[p + 2] {
                            masked[p] = true;
                            masked[p + 1] = true;
                            masked[p + 2] = true;
                        }
                    }
                }
            }
        }
    }

    masked
}

/// Extract canonical strand-symmetric spaced-seed positions with at most `max_occ` hits.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if `mask` is empty, exceeds 31 active bits, or contains invalid characters.
pub fn canonical_spaced_seeds(
    seq: &[u8],
    mask: &str,
    max_occ: usize,
) -> Result<HashMap<u64, Vec<usize>>, AlignError> {
    canonical_spaced_seeds_impl(seq, mask, max_occ, false, false)
}

#[cfg(test)]
pub(crate) fn canonical_spaced_seeds_with_opts(
    seq: &[u8],
    mask: &str,
    max_occ: usize,
    skip_softmasked: bool,
    ry_mode: bool,
) -> Result<HashMap<u64, Vec<usize>>, AlignError> {
    canonical_spaced_seeds_impl(seq, mask, max_occ, skip_softmasked, ry_mode)
}

pub(crate) fn parse_mask_offsets(
    mask: &str,
    ry_mode: bool,
) -> Result<(usize, Vec<usize>), AlignError> {
    if mask.is_empty() {
        return Err(AlignError::InvalidMask("mask is empty".to_owned()));
    }
    let bytes = mask.as_bytes();
    let span = bytes.len();
    let mut offsets = Vec::new();
    for (idx, &ch) in bytes.iter().enumerate() {
        match ch {
            b'1' => offsets.push(idx),
            b'0' => {}
            _ => {
                return Err(AlignError::InvalidMask(format!(
                    "unexpected char '{}'",
                    ch as char
                )));
            }
        }
    }
    let max_bits = if ry_mode { 62 } else { 31 };
    if offsets.is_empty() || offsets.len() > max_bits {
        return Err(AlignError::InvalidMask(format!(
            "active mask bits must be in 1..={max_bits}"
        )));
    }
    Ok((span, offsets))
}

fn canonical_spaced_seeds_impl(
    seq: &[u8],
    mask: &str,
    max_occ: usize,
    skip_softmasked: bool,
    ry_mode: bool,
) -> Result<HashMap<u64, Vec<usize>>, AlignError> {
    let (span, offsets) = parse_mask_offsets(mask, ry_mode)?;
    if seq.len() < span {
        return Ok(HashMap::new());
    }

    let rev_offsets: Vec<usize> = offsets.iter().map(|&o| span - 1 - o).collect();
    let cap = (seq.len() / 2).min(1 << 20);
    let mut counts: HashMap<u64, usize> = HashMap::with_capacity(cap);
    let mut positions: HashMap<u64, Vec<usize>> = HashMap::with_capacity(cap);

    for i in 0..=(seq.len() - span) {
        let sub = &seq[i..i + span];
        if span >= 12 && is_low_complexity_window(sub) {
            continue;
        }
        let Some((fwd, rc)) =
            encode_spaced_pair(sub, &offsets, &rev_offsets, skip_softmasked, ry_mode)
        else {
            continue;
        };
        let key = fwd.min(rc);
        let entry = counts.entry(key).or_insert(0);
        *entry += 1;
        if *entry <= max_occ {
            positions.entry(key).or_default().push(i);
        } else if *entry == max_occ + 1 {
            positions.remove(&key);
        }
    }

    Ok(positions)
}

mod chain;
pub use chain::{cluster_strand_lcbs, lis_chain, local_w9_subchain};
pub(crate) use chain::{lcb_interval_set, merge_colinear_lcbs};

pub(crate) use flat_index::FlatSeedIndex;

mod sketch;
pub use sketch::{
    ContigIndex, IndexedGenome, adaptive_top_k_for_fragmentation, effective_contig_count,
    merge_contig_sketches, prepare_contig_index, prepare_contig_index_owned,
    select_top_k_for_target_refs, select_top_k_pairs, shared_seed_count, sketch_jaccard,
};

pub(crate) type LcbList = Vec<Vec<(usize, usize)>>;

pub(crate) type AnchorPairs = Vec<(usize, usize)>;

pub(crate) fn encode_spaced_pair(
    sub: &[u8],
    offsets: &[usize],
    rev_offsets: &[usize],
    skip_softmasked: bool,
    ry_mode: bool,
) -> Option<(u64, u64)> {
    let mut fwd = 0_u64;
    let mut rc = 0_u64;
    for (&fo, &ro) in offsets.iter().zip(rev_offsets) {
        if skip_softmasked && (sub[fo].is_ascii_lowercase() || sub[ro].is_ascii_lowercase()) {
            return None;
        }
        let (bf, br) = (base_to_idx(sub[fo])?, base_to_idx(sub[ro])?);
        if ry_mode {
            fwd = (fwd << 1) | ((bf & 1) as u64);
            rc = (rc << 1) | ((1 - (br & 1)) as u64);
        } else {
            fwd = (fwd << 2) | (bf as u64);
            rc = (rc << 2) | ((3 - br) as u64);
        }
    }
    Some((fwd, rc))
}

fn push_transition_keys(
    fwd: u64,
    rc: u64,
    w: usize,
    min_transitions: usize,
    max_transitions: usize,
    keys: &mut Vec<u64>,
) {
    keys.clear();
    if min_transitions == 0 {
        keys.push(fwd);
        keys.push(rc);
    }
    if min_transitions <= 1 && max_transitions >= 1 {
        for b in 0..w {
            let mask = 2_u64 << (2 * (w - 1 - b));
            keys.push(fwd ^ mask);
            keys.push(rc ^ mask);
        }
    }
    if max_transitions >= 2 {
        for a in 0..w {
            let mask_a = 2_u64 << (2 * (w - 1 - a));
            let ma_f = fwd ^ mask_a;
            let ma_r = rc ^ mask_a;
            for b in (a + 1)..w {
                let mask_b = 2_u64 << (2 * (w - 1 - b));
                keys.push(ma_f ^ mask_b);
                keys.push(ma_r ^ mask_b);
            }
        }
    }
    keys.sort_unstable();
    keys.dedup();
}
const LASTZ12_OFFSETS: [usize; 12] = [0, 1, 2, 4, 7, 8, 11, 13, 15, 16, 17, 18];
const LASTZ12_REV_OFFSETS: [usize; 12] = [18, 17, 16, 14, 11, 10, 7, 5, 3, 2, 1, 0];

#[inline]
pub(crate) fn is_masked_window(seq: &[u8], sdust: Option<&[bool]>, i: usize, span: usize) -> bool {
    if let Some(m) = sdust {
        let mid = m.get(i + span / 2).copied().unwrap_or(false);
        let end = span
            .checked_sub(1)
            .and_then(|d| m.get(i + d).copied())
            .unwrap_or(false);
        if m.get(i).copied().unwrap_or(false) || mid || end {
            return true;
        }
    }
    sdust.is_none() && span >= 12 && is_low_complexity_window(&seq[i..i + span])
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn find_spaced_transition_hits_flat(
    seq_x: &[u8],
    sdust_x: Option<&[bool]>,
    offsets: &[usize],
    rev_offsets: &[usize],
    span: usize,
    iy: &FlatSeedIndex,
    skip_softmasked: bool,
    min_transitions: usize,
    max_transitions: usize,
) -> (AnchorPairs, AnchorPairs) {
    const CHUNK_LEN: usize = 32_768;
    debug_assert!(sdust_x.is_none_or(|m| m.len() >= seq_x.len()));
    if seq_x.len() < span {
        return (Vec::new(), Vec::new());
    }
    let n_pos = seq_x.len() - span + 1;
    let scan_range = |start_i: usize, end_i: usize| -> (AnchorPairs, AnchorPairs) {
        let w = offsets.len();
        let mut hits_fwd = Vec::new();
        let mut hits_rc = Vec::new();
        let mut keys = Vec::with_capacity(w * w + 2 * w + 1);
        for i in start_i..end_i {
            if is_masked_window(seq_x, sdust_x, i, span) {
                continue;
            }
            let sub = &seq_x[i..i + span];
            let Some((fwd, rc)) =
                encode_spaced_pair(sub, offsets, rev_offsets, skip_softmasked, false)
            else {
                continue;
            };
            push_transition_keys(fwd, rc, w, min_transitions, max_transitions, &mut keys);
            for &key in &keys {
                for &py in iy.get(key) {
                    let py_u = py as usize;
                    hits_fwd.push((i, py_u));
                    hits_rc.push((i, iy.seq_len - iy.span - py_u));
                }
            }
        }
        (hits_fwd, hits_rc)
    };

    if n_pos < 65_536 {
        return scan_range(0, n_pos);
    }
    let chunks: Vec<(usize, usize)> = (0..n_pos)
        .step_by(CHUNK_LEN)
        .map(|s| (s, (s + CHUNK_LEN).min(n_pos)))
        .collect();
    let per_chunk: Vec<(AnchorPairs, AnchorPairs)> = chunks
        .into_par_iter()
        .map(|(s, e)| scan_range(s, e))
        .collect();
    let total_fwd: usize = per_chunk.iter().map(|(f, _)| f.len()).sum();
    let total_rc: usize = per_chunk.iter().map(|(_, r)| r.len()).sum();
    let mut hits_fwd = Vec::with_capacity(total_fwd);
    let mut hits_rc = Vec::with_capacity(total_rc);
    for (mut f, mut r) in per_chunk {
        hits_fwd.append(&mut f);
        hits_rc.append(&mut r);
    }
    (hits_fwd, hits_rc)
}

pub(crate) fn find_lastz12_1transition_hits_flat(
    seq_x: &[u8],
    sdust_x: Option<&[bool]>,
    iy: &FlatSeedIndex,
    skip_softmasked: bool,
) -> (AnchorPairs, AnchorPairs) {
    find_spaced_transition_hits_flat(
        seq_x,
        sdust_x,
        &LASTZ12_OFFSETS,
        &LASTZ12_REV_OFFSETS,
        19,
        iy,
        skip_softmasked,
        0,
        1,
    )
}

pub(crate) fn find_lastz12_2transition_hits_flat(
    seq_x: &[u8],
    sdust_x: Option<&[bool]>,
    iy: &FlatSeedIndex,
    skip_softmasked: bool,
) -> (AnchorPairs, AnchorPairs) {
    find_spaced_transition_hits_flat(
        seq_x,
        sdust_x,
        &LASTZ12_OFFSETS,
        &LASTZ12_REV_OFFSETS,
        19,
        iy,
        skip_softmasked,
        2,
        2,
    )
}

pub(crate) fn find_lastz11_1transition_hits_flat(
    seq_x: &[u8],
    sdust_x: Option<&[bool]>,
    iy: &FlatSeedIndex,
    skip_softmasked: bool,
) -> (AnchorPairs, AnchorPairs) {
    const OFFSETS: [usize; 11] = [0, 1, 2, 4, 7, 8, 11, 13, 15, 16, 17];
    const REV_OFFSETS: [usize; 11] = [17, 16, 15, 13, 10, 9, 6, 4, 2, 1, 0];
    find_spaced_transition_hits_flat(
        seq_x,
        sdust_x,
        &OFFSETS,
        &REV_OFFSETS,
        18,
        iy,
        skip_softmasked,
        0,
        1,
    )
}

pub(crate) fn find_w10_1transition_hits_flat(
    seq_x: &[u8],
    sdust_x: Option<&[bool]>,
    iy: &FlatSeedIndex,
    skip_softmasked: bool,
) -> (AnchorPairs, AnchorPairs) {
    const OFFSETS: [usize; 10] = [0, 1, 3, 4, 6, 7, 9, 10, 12, 13];
    const REV_OFFSETS: [usize; 10] = [13, 12, 10, 9, 7, 6, 4, 3, 1, 0];
    find_spaced_transition_hits_flat(
        seq_x,
        sdust_x,
        &OFFSETS,
        &REV_OFFSETS,
        14,
        iy,
        skip_softmasked,
        0,
        1,
    )
}

/// Promote 2-hit gapped diagonal pairs (`<= max_gap` bp apart, `<= max_ddiag` bp indel shift)
/// into candidate LCBs for banded MDL DP.
pub(crate) fn append_two_hit_diagonal_lcbs(
    target: &mut LcbList,
    hits: &[(usize, usize)],
    max_gap: usize,
    max_ddiag: usize,
) {
    if hits.len() < 2 {
        return;
    }
    let mut sorted = hits.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    let mut promoted: Vec<((usize, usize), (usize, usize))> = Vec::new();
    for i in 0..sorted.len() {
        let (x1, y1) = sorted[i];
        let mut j = i + 1;
        while j < sorted.len() && sorted[j].0 - x1 <= max_gap {
            let (x2, y2) = sorted[j];
            if x2 > x1 && y2 > y1 && y2 - y1 <= max_gap && (x2 - x1).abs_diff(y2 - y1) <= max_ddiag
            {
                promoted.push(((x1, y1), (x2, y2)));
                break;
            }
            j += 1;
        }
    }
    let mut intervals = lcb_interval_set(target.iter());
    for (p1, p2) in promoted {
        let lo = p1.0.min(p2.0);
        let hi = p1.0.max(p2.0).saturating_add(1);
        if intervals.insert_if_disjoint(lo, hi) {
            target.push(vec![p1, p2]);
        }
    }
}

pub mod cdbg;
pub mod eisenstein;
mod flat_index;
mod hsp;
pub mod tantan;
mod translated;
pub use cdbg::{
    CdbgCollinearBlock, CdbgUnitig, UnitigOccurrence, build_colored_cdbg,
    extract_cdbg_collinear_blocks,
};
pub use eisenstein::{
    EISENSTEIN_OMEGA, EisensteinHasher, EisensteinSamplerSketch, MERSENNE_61, PowerSumSketch,
    select_auto_hinges_eisenstein,
};
pub(crate) use hsp::*;
pub use tantan::{TantanConfig, tantan_mask, tantan_posterior_probs};
pub(crate) use translated::*;

pub use crate::dna::{complement_base, revcomp, revcomp_preserve_case};

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::splitmix64;

    #[test]
    fn non_palindromic_revcomp_seeds_match_exactly() {
        let seq = b"ACGTTAGCAGTCAGTCAGTCAGGTAACCGGTTAACCGATGCTAGCTAG";
        let rc = revcomp(seq);
        let fwd_map = canonical_spaced_seeds(seq, MASK_PRIMARY_11_15, 4).unwrap();
        let rc_map = canonical_spaced_seeds(&rc, MASK_PRIMARY_11_15, 4).unwrap();
        assert_eq!(fwd_map.len(), rc_map.len());
        for k in fwd_map.keys() {
            assert!(rc_map.contains_key(k), "missing rc seed key {k}");
        }

        let idx_f = prepare_contig_index(seq).unwrap();
        let idx_r = prepare_contig_index(&rc).unwrap();
        assert_eq!(shared_seed_count(&idx_f, &idx_r), fwd_map.len());
        let (_hits_fwd, hits_rc) = idx_f.seeds.intersect(&idx_r.seeds);
        assert!(hits_rc.len() >= fwd_map.len());
    }

    #[test]
    fn sdust_masks_simple_repeats_and_preserves_unique_dna() {
        let mut seq = b"ACGTTGCATGCAACTGGTCAGTCGATCGTAGCTAGCTGAT".to_vec(); // 40bp unique
        let rep_start = seq.len();
        seq.extend_from_slice(b"CACACACACACACACACACACACACACACACACACACACA"); // 40bp (CA)20
        let rep_end = seq.len();
        seq.extend_from_slice(b"TGACCTAGGCTAGATTGCCAATGACAGTTCCAGGCTATGC"); // 40bp unique

        let mask = sdust_mask(&seq);
        assert!(mask[..16].iter().all(|&m| !m), "unique prefix was masked");
        assert!(
            mask[rep_start + 4..rep_end - 4].iter().all(|&m| m),
            "dinucleotide repeat core was not masked"
        );
        assert!(
            mask[rep_end + 16..].iter().all(|&m| !m),
            "unique suffix was masked"
        );
    }

    #[test]
    fn flat_seed_index_buckets_and_parallel_transition_hits_match_expected() {
        // Deterministic pseudo-random 70 kb sequence (>= 32_768 keys and >= 65_536 positions).
        let bases = *b"ACGT";
        let mut state: u64 = 0x1234_5678_9abc_def0;
        let mut seq = Vec::with_capacity(70_000);
        for _ in 0..70_000 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            seq.push(bases[((state >> 32) & 3) as usize]);
        }

        let idx = FlatSeedIndex::build_with_sdust(&seq, MASK_LASTZ12, 8, false, false, None)
            .expect("build FlatSeedIndex");
        assert!(idx.keys.len() >= 32_768);
        assert_eq!(idx.buckets.len(), 65_537);
        for (k_i, &k) in idx.keys.iter().enumerate().step_by(257) {
            let got = idx.get(k);
            assert!(got.contains(&idx.pos[k_i]));
        }
        assert_eq!(idx.get(u64::MAX), []);

        let (hits_fwd, hits_rc) = find_lastz12_1transition_hits_flat(&seq, None, &idx, false);
        assert!(hits_fwd.len() >= 60_000);
        assert_eq!(hits_fwd.len(), hits_rc.len());
    }

    #[test]
    fn asymmetric_lastz11_finds_rc_1transition_at_unmirrored_offset() {
        // Construct y where canonical key is rc < fwd, then mutate x at asymmetric MASK_LASTZ11
        // offset 7 (which has no mirror in OFFSETS, corresponding to rc bit 6).
        let y = b"TTTTTTTTACGTACGTAC".to_vec(); // 18 bp; revcomp starts with GTACGTACGTAAAAAAAA (smaller)
        let mut x = y.clone();
        x[7] = if x[7] == b'A' { b'G' } else { b'A' }; // 1 transition at asymmetric offset 7
        let iy = FlatSeedIndex::build_with_sdust(&y, MASK_LASTZ11, 4, false, false, None)
            .expect("build FlatSeedIndex");
        let (hits_fwd, _) = find_lastz11_1transition_hits_flat(&x, None, &iy, false);
        assert!(
            hits_fwd.contains(&(0, 0)),
            "expected asymmetric MASK_LASTZ11 rc 1-transition hit at (0, 0), got {hits_fwd:?}"
        );
    }

    #[test]
    fn hardcoded_mask_offsets_match_parse_mask_offsets() {
        let (span12, off12) = parse_mask_offsets(MASK_LASTZ12, false).expect("valid MASK_LASTZ12");
        assert_eq!(off12, LASTZ12_OFFSETS);
        let rev12: Vec<usize> = off12.iter().map(|&o| span12 - 1 - o).collect();
        assert_eq!(rev12, LASTZ12_REV_OFFSETS);

        let (span11, off11) = parse_mask_offsets(MASK_LASTZ11, false).expect("valid MASK_LASTZ11");
        assert_eq!(off11, [0, 1, 2, 4, 7, 8, 11, 13, 15, 16, 17]);
        let rev11: Vec<usize> = off11.iter().map(|&o| span11 - 1 - o).collect();
        assert_eq!(rev11, [17, 16, 15, 13, 10, 9, 6, 4, 2, 1, 0]);

        let (span10, off10) = parse_mask_offsets(MASK_W10, false).expect("valid MASK_W10");
        assert_eq!(off10, [0, 1, 3, 4, 6, 7, 9, 10, 12, 13]);
        let rev10: Vec<usize> = off10.iter().map(|&o| span10 - 1 - o).collect();
        assert_eq!(rev10, [13, 12, 10, 9, 7, 6, 4, 3, 1, 0]);
    }

    #[test]
    fn shannon_entropy_filters_low_complexity_windows() {
        assert!(is_low_complexity_window(b"AAAAAAAAAAAAAAAA"));
        assert!(is_low_complexity_window(b"ATATATATATATATAT")); // H = 1.0 <= 1.05
        assert!(is_low_complexity_window(b"NNNNNNNNNNNNNNNN"));
        assert!(!is_low_complexity_window(b"ACGTACGTACGTACGT")); // H = 2.0 > 1.05
        assert!(!is_low_complexity_window(b"ACGTTCGAACGTACTG"));
    }

    #[test]
    fn minhash_sketch_estimates_shared_seeds() {
        let seq_a = b"ACGTACGTACGTACGTTTGCATGCATGCATGCATGCAACCGGTTAA".repeat(4);
        let mut seq_b = seq_a.clone();
        seq_b[40..50].copy_from_slice(b"GGGGTTTTAA");
        let seq_c = b"CCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT".repeat(4);

        let ia = prepare_contig_index(&seq_a).expect("ia");
        let ib = prepare_contig_index(&seq_b).expect("ib");
        let ic = prepare_contig_index(&seq_c).expect("ic");

        let sa = ia.minhash_sketch(64);
        let sb = ib.minhash_sketch(64);
        let sc = ic.minhash_sketch(64);
        assert_ne!(sa, [] as [u64; 0]);
        assert!(sa.windows(2).all(|w| w[0] < w[1]));
        assert!(sketch_jaccard(&sa, &sb) > sketch_jaccard(&sa, &sc));
    }

    #[test]
    fn merge_contig_sketches_keeps_k_smallest_unique_hashes() {
        let seq1 = b"ACGTACGTTTGCATGCATGCATGCAACCGGTTAACGTACGTGGGGTTTT".repeat(3);
        let seq2 = b"TGCATGCAAACCGGTTACGTACGTCCCGGGTTTAAACGTACGTACGTAA".repeat(3);
        let i1 = prepare_contig_index(&seq1).expect("i1");
        let i2 = prepare_contig_index(&seq2).expect("i2");

        let merged = merge_contig_sketches(&[("c1".to_owned(), i1), ("c2".to_owned(), i2)], 16);
        assert!(merged.len() <= 16);
        assert_ne!(merged, [] as [u64; 0]);
        assert!(
            merged.windows(2).all(|w| w[0] < w[1]),
            "merged sketch must be strictly sorted and deduplicated"
        );
    }

    #[test]
    fn sketch_jaccard_is_unbiased_on_unequal_sized_sets() {
        assert!((sketch_jaccard(&[], &[1, 2]) - 0.0).abs() < 1e-12);
        assert!((sketch_jaccard(&[1, 3, 5], &[1, 3, 5]) - 1.0).abs() < 1e-12);
        assert!((sketch_jaccard(&[1, 3, 5], &[2, 4, 6]) - 0.0).abs() < 1e-12);

        // Deterministic 50 kb pseudo-random sequence and its 5 kb prefix (true Jaccard ~ 0.10).
        // Under the biased c/(|A|+|B|-c) formula, 10:1 unequal sizes yield ~0.05; Broder's
        // MIN_k(A ∪ B) estimator recovers ~0.10 without downward size bias.
        let mut state = 0x5EED_CAFE_1234_5678_u64;
        let bases = *b"ACGT";
        let long_seq: Vec<u8> = (0..50_000)
            .map(|_| {
                state = splitmix64(state);
                bases[(state & 3) as usize]
            })
            .collect();
        let short_seq = &long_seq[..5_000];
        let idx_long = prepare_contig_index(&long_seq).expect("long");
        let idx_short = prepare_contig_index(short_seq).expect("short");
        let s_long = idx_long.minhash_sketch(512);
        let s_short = idx_short.minhash_sketch(512);
        let j_hat = sketch_jaccard(&s_short, &s_long);
        assert!(
            (0.07..=0.14).contains(&j_hat),
            "expected Broder MIN_k(A ∪ B) Jaccard near 0.10 for 5kb/50kb prefix, got {j_hat:.4}"
        );
    }
}
