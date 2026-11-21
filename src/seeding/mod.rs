//! Spaced-seed indexing, transition-tolerant and translated seeding, and LCB chaining.

use crate::AlignError;
use crate::bocpd::base_to_idx;
use rayon::prelude::*;
use std::collections::HashMap;

pub use translated::aa5_keys;

/// KP-14 spaced seed mask (weight 11, span 15).
pub const MASK_PRIMARY_11_15: &str = "111011011011011";
/// Optimal non-periodic weight-11 `PatternHunter` spaced seed (Ma et al., 2002).
pub const MASK_PH11_18: &str = "111010010100110111";
/// Secondary `w=10` spaced seed mask for divergent intra-LCB chaining.
pub const MASK_W10: &str = "11011011011011";
/// Local `w=9` spaced seed mask for intra-gap anchor seeding.
pub const MASK_W9: &str = "110101101011011";
pub(crate) const MASK_RY20: &str = "1110110110111011011011011";
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

/// Extract canonical spaced-seed positions with optional softmask skipping and 1-bit `R/Y`
/// transition-tolerant encoding.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if `mask` is empty, exceeds active bit limits, or contains invalid characters.
pub fn canonical_spaced_seeds_with_opts(
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
pub(crate) use chain::{IntervalSet, merge_colinear_lcbs};
pub use chain::{cluster_strand_lcbs, lis_chain, local_w9_subchain};

pub(crate) use flat_index::FlatSeedIndex;

/// Precomputed strand-symmetric KP-14 seed index for one contig.
#[derive(Debug, Clone)]
pub struct ContigIndex {
    pub(crate) seq: Vec<u8>,
    pub(crate) sdust: Vec<bool>,
    pub(crate) seeds: FlatSeedIndex,
}

/// Precompute forward and reverse-complement KP-14 seed maps for a contig.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if `MASK_PRIMARY_11_15` is invalid.
pub fn prepare_contig_index(seq: &[u8]) -> Result<ContigIndex, AlignError> {
    prepare_contig_index_opts(seq, false)
}

/// Precompute forward and reverse-complement KP-14 seed maps, optionally skipping lowercase
/// softmasked bases during seeding.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if `MASK_PRIMARY_11_15` is invalid.
pub fn prepare_contig_index_opts(
    seq: &[u8],
    skip_softmasked: bool,
) -> Result<ContigIndex, AlignError> {
    // Weight-11 canonical spaced-seed space is 4^11 / 2 = 2,097,152 keys; once contig length
    // exceeds ~1/21 of the key space (~100 kb), tighten max_occ from 4 to 2 to bound collisions.
    const PRIMARY_SEED_WEIGHT: u32 = 11;
    let canonical_key_space = 4_usize.pow(PRIMARY_SEED_WEIGHT) / 2;
    let max_occ = if seq.len().saturating_mul(21) >= canonical_key_space {
        2
    } else {
        4
    };
    let sdust = sdust_mask(seq);
    let seeds = FlatSeedIndex::build_with_sdust(
        seq,
        MASK_PRIMARY_11_15,
        max_occ,
        skip_softmasked,
        false,
        Some(&sdust),
    )?;
    Ok(ContigIndex {
        seq: seq.to_vec(),
        sdust,
        seeds,
    })
}

/// Count shared KP-14 seeds (on either strand) between two pre-indexed contigs.
#[must_use]
pub fn shared_seed_count(idx_x: &ContigIndex, idx_y: &ContigIndex) -> usize {
    idx_x.seeds.shared_unique_keys(&idx_y.seeds)
}

#[inline]
const fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

impl ContigIndex {
    /// Extract a sorted bottom-`k` `MinHash` sketch of unique KP-14 seed keys in this contig.
    #[must_use]
    pub fn minhash_sketch(&self, k: usize) -> Vec<u64> {
        use std::collections::BinaryHeap;
        if k == 0 || self.seeds.keys.is_empty() {
            return Vec::new();
        }
        let mut heap: BinaryHeap<u64> = BinaryHeap::with_capacity(k);
        let mut prev: Option<u64> = None;
        for &raw in &self.seeds.keys {
            if prev == Some(raw) {
                continue;
            }
            prev = Some(raw);
            let h = splitmix64(raw);
            if heap.len() < k {
                heap.push(h);
            } else if let Some(&top) = heap.peek()
                && h < top
            {
                heap.pop();
                heap.push(h);
            }
        }
        heap.into_sorted_vec()
    }
}

/// Count shared hashes between two sorted bottom-`k` `MinHash` sketches.
#[must_use]
pub fn shared_sketch_count(a: &[u64], b: &[u64]) -> usize {
    let mut i = 0usize;
    let mut j = 0usize;
    let mut shared = 0usize;
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                shared += 1;
                i += 1;
                j += 1;
            }
        }
    }
    shared
}

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
        keys.push(fwd.min(rc));
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

#[allow(clippy::too_many_arguments)]
pub(crate) fn find_spaced_transition_hits<S: std::hash::BuildHasher>(
    seq_x: &[u8],
    offsets: &[usize],
    rev_offsets: &[usize],
    span: usize,
    iy_fwd: &HashMap<u64, Vec<usize>, S>,
    iy_rc: &HashMap<u64, Vec<usize>, S>,
    skip_softmasked: bool,
    min_transitions: usize,
    max_transitions: usize,
) -> (AnchorPairs, AnchorPairs) {
    if seq_x.len() < span {
        return (Vec::new(), Vec::new());
    }
    let w = offsets.len();
    let mut hits_fwd = Vec::new();
    let mut hits_rc = Vec::new();
    let mut keys = Vec::with_capacity(w * w + 2 * w + 1);
    for i in 0..=(seq_x.len() - span) {
        let Some((fwd, rc)) = encode_spaced_pair(
            &seq_x[i..i + span],
            offsets,
            rev_offsets,
            skip_softmasked,
            false,
        ) else {
            continue;
        };
        push_transition_keys(fwd, rc, w, min_transitions, max_transitions, &mut keys);
        for &key in &keys {
            if let Some(pys) = iy_fwd.get(&key) {
                for &py in pys {
                    hits_fwd.push((i, py));
                }
            }
            if let Some(pys) = iy_rc.get(&key) {
                for &py in pys {
                    hits_rc.push((i, py));
                }
            }
        }
    }
    (hits_fwd, hits_rc)
}

const LASTZ12_OFFSETS: [usize; 12] = [0, 1, 2, 4, 7, 8, 11, 13, 15, 16, 17, 18];
const LASTZ12_REV_OFFSETS: [usize; 12] = [18, 17, 16, 14, 11, 10, 7, 5, 3, 2, 1, 0];

/// Find `12-of-19` (`MASK_LASTZ12`) seed hits including 1-transition (`A<->G`, `C<->T`) neighbors.
#[must_use]
pub fn find_lastz12_1transition_hits<S: std::hash::BuildHasher>(
    seq_x: &[u8],
    iy_fwd: &HashMap<u64, Vec<usize>, S>,
    iy_rc: &HashMap<u64, Vec<usize>, S>,
    skip_softmasked: bool,
) -> (AnchorPairs, AnchorPairs) {
    find_spaced_transition_hits(
        seq_x,
        &LASTZ12_OFFSETS,
        &LASTZ12_REV_OFFSETS,
        19,
        iy_fwd,
        iy_rc,
        skip_softmasked,
        0,
        1,
    )
}

/// Find `12-of-19` (`MASK_LASTZ12`) seed hits with 2 simultaneous transitions (`A<->G`, `C<->T`).
#[must_use]
pub fn find_lastz12_2transition_hits<S: std::hash::BuildHasher>(
    seq_x: &[u8],
    iy_fwd: &HashMap<u64, Vec<usize>, S>,
    iy_rc: &HashMap<u64, Vec<usize>, S>,
    skip_softmasked: bool,
) -> (AnchorPairs, AnchorPairs) {
    find_spaced_transition_hits(
        seq_x,
        &LASTZ12_OFFSETS,
        &LASTZ12_REV_OFFSETS,
        19,
        iy_fwd,
        iy_rc,
        skip_softmasked,
        2,
        2,
    )
}

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
    let mut intervals = IntervalSet::from_lcbs(target.iter());
    for (p1, p2) in promoted {
        if intervals.insert_if_disjoint(p1.0, p2.0) {
            target.push(vec![p1, p2]);
        }
    }
}

mod flat_index;
mod hsp;
mod translated;
pub(crate) use hsp::*;
pub(crate) use translated::*;

/// Fast pre-DP score: return the anchor count of the largest `w=12` LCB between two contigs.
#[must_use]
pub fn score_contig_pair_lcb(idx_x: &ContigIndex, idx_y: &ContigIndex) -> usize {
    let (raw_fwd, raw_rc) = find_w12_raw_lcbs(idx_x, idx_y);
    raw_fwd
        .iter()
        .chain(raw_rc.iter())
        .map(Vec::len)
        .max()
        .unwrap_or(0)
}

/// Select the union of top-`K` target-per-query and query-per-target contig pairs between two assemblies.
#[must_use]
pub fn select_top_k_pairs(
    target_idx: &[(String, ContigIndex)],
    query_idx: &[(String, ContigIndex)],
    top_k: usize,
    min_shared_seeds: usize,
) -> Vec<Vec<bool>> {
    let seed_counts: Vec<Vec<usize>> = (0..target_idx.len())
        .into_par_iter()
        .map(|x_i| {
            (0..query_idx.len())
                .map(|y_i| shared_seed_count(&target_idx[x_i].1, &query_idx[y_i].1))
                .collect()
        })
        .collect();

    let shortlist_k = top_k.saturating_mul(4).max(8);
    let mut need_lcb = vec![vec![false; query_idx.len()]; target_idx.len()];
    let mut row_cands: Vec<(usize, usize)> = Vec::new();
    for (x_i, row) in seed_counts.iter().enumerate() {
        row_cands.clear();
        for (y_i, &seeds) in row.iter().enumerate() {
            if seeds >= min_shared_seeds {
                row_cands.push((seeds, y_i));
            }
        }
        row_cands.sort_unstable_by_key(|&(s, y_i)| (std::cmp::Reverse(s), y_i));
        for &(_, y_i) in row_cands.iter().take(shortlist_k) {
            need_lcb[x_i][y_i] = true;
        }
    }
    let mut col_cands: Vec<(usize, usize)> = Vec::new();
    for y_i in 0..query_idx.len() {
        col_cands.clear();
        for (x_i, row) in seed_counts.iter().enumerate() {
            let seeds = row[y_i];
            if seeds >= min_shared_seeds {
                col_cands.push((seeds, x_i));
            }
        }
        col_cands.sort_unstable_by_key(|&(s, x_i)| (std::cmp::Reverse(s), x_i));
        for &(_, x_i) in col_cands.iter().take(shortlist_k) {
            need_lcb[x_i][y_i] = true;
        }
    }

    let mut lcb_jobs: Vec<(usize, usize)> = Vec::new();
    for (x_i, row) in need_lcb.iter().enumerate() {
        for (y_i, &needed) in row.iter().enumerate() {
            if needed {
                lcb_jobs.push((x_i, y_i));
            }
        }
    }
    let lcb_scores: Vec<((usize, usize), usize)> = lcb_jobs
        .into_par_iter()
        .map(|(x_i, y_i)| {
            let lcb = score_contig_pair_lcb(&target_idx[x_i].1, &query_idx[y_i].1);
            ((x_i, y_i), lcb)
        })
        .collect();

    let mut by_x: Vec<Vec<((usize, usize), usize)>> = vec![Vec::new(); target_idx.len()];
    let mut by_y: Vec<Vec<((usize, usize), usize)>> = vec![Vec::new(); query_idx.len()];
    for ((x_i, y_i), lcb) in lcb_scores {
        let seeds = seed_counts[x_i][y_i];
        by_x[x_i].push(((lcb, seeds), y_i));
        by_y[y_i].push(((lcb, seeds), x_i));
    }
    let mut mask = vec![vec![false; query_idx.len()]; target_idx.len()];
    for (x_i, list) in by_x.iter_mut().enumerate() {
        list.sort_unstable_by_key(|a| (std::cmp::Reverse(a.0), a.1));
        for &(_, y_i) in list.iter().take(top_k) {
            mask[x_i][y_i] = true;
        }
    }
    for (y_i, list) in by_y.iter_mut().enumerate() {
        list.sort_unstable_by_key(|a| (std::cmp::Reverse(a.0), a.1));
        for &(_, x_i) in list.iter().take(top_k) {
            mask[x_i][y_i] = true;
        }
    }
    mask
}

/// Complement a single DNA base (preserving case for `A`/`C`/`G`/`T`).
#[must_use]
pub fn complement_base(b: u8) -> u8 {
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
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

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
        assert!(shared_sketch_count(&sa, &sb) > shared_sketch_count(&sa, &sc));
    }
}
