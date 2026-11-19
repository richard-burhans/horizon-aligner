//! Spaced-seed indexing, transition-tolerant and translated seeding, and LCB chaining.

use crate::AlignError;
use crate::bocpd::base_to_idx;
use rayon::prelude::*;
use std::collections::HashMap;

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
    for &b in sub {
        if let Some(idx) = base_to_idx(b) {
            counts[idx] += 1;
        }
    }
    let max_c = counts[0].max(counts[1]).max(counts[2]).max(counts[3]);
    let distinct = counts.iter().filter(|&&c| c > 0).count();
    max_c * 4 >= sub.len() * 3 || (sub.len() >= 14 && distinct <= 2)
}

/// Linear-time SDUST-style triplet low-complexity masker (`W in {16, 32, 64}`, `T = 2.0`,
/// Morgulis et al., 2006). Masks homopolymers, dinucleotide, trinucleotide, and tetranucleotide
/// simple repeats in `O(N)` time so seeders avoid `O(span)` per-position checks.
#[must_use]
pub(crate) fn sdust_mask(seq: &[u8]) -> Vec<bool> {
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
                        masked[p] = true;
                        masked[p + 1] = true;
                        masked[p + 2] = true;
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

fn canonical_spaced_seeds_impl(
    seq: &[u8],
    mask: &str,
    max_occ: usize,
    skip_softmasked: bool,
    ry_mode: bool,
) -> Result<HashMap<u64, Vec<usize>>, AlignError> {
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
    if seq.len() < span {
        return Ok(HashMap::new());
    }

    let rev_offsets: Vec<usize> = offsets.iter().map(|&o| span - 1 - o).collect();
    let mut counts: HashMap<u64, usize> = HashMap::with_capacity(seq.len() / 2);
    let mut positions: HashMap<u64, Vec<usize>> = HashMap::with_capacity(seq.len() / 2);

    for i in 0..=(seq.len() - span) {
        let sub = &seq[i..i + span];
        if span >= 12 && is_low_complexity_window(sub) {
            continue;
        }
        let mut fwd = 0_u64;
        let mut rc = 0_u64;
        let mut valid = true;
        for (&fo, &ro) in offsets.iter().zip(&rev_offsets) {
            if skip_softmasked && (sub[fo].is_ascii_lowercase() || sub[ro].is_ascii_lowercase()) {
                valid = false;
                break;
            }
            let (Some(bf), Some(br)) = (base_to_idx(sub[fo]), base_to_idx(sub[ro])) else {
                valid = false;
                break;
            };
            if ry_mode {
                fwd = (fwd << 1) | ((bf & 1) as u64);
                rc = (rc << 1) | ((1 - (br & 1)) as u64);
            } else {
                fwd = (fwd << 2) | (bf as u64);
                rc = (rc << 2) | ((3 - br) as u64);
            }
        }
        if !valid {
            continue;
        }
        let key = fwd.min(rc);
        let entry = counts.entry(key).or_insert(0);
        *entry += 1;
        if *entry <= max_occ {
            positions.entry(key).or_default().push(i);
        }
    }

    positions.retain(|k, _| counts.get(k).copied().unwrap_or(usize::MAX) <= max_occ);
    Ok(positions)
}

/// Compute the Longest Increasing Subsequence chain over `(x, y)` anchors.
#[must_use]
pub fn lis_chain(anchors: &[(usize, usize)]) -> Vec<(usize, usize)> {
    if anchors.is_empty() {
        return Vec::new();
    }
    let mut sorted = anchors.to_vec();
    sorted.sort_unstable();

    let mut tails_y: Vec<usize> = Vec::with_capacity(sorted.len());
    let mut tail_idx: Vec<usize> = Vec::with_capacity(sorted.len());
    let mut prev = vec![usize::MAX; sorted.len()];

    for (idx, &(_, y)) in sorted.iter().enumerate() {
        let pos = tails_y.partition_point(|&ty| ty < y);
        if pos == tails_y.len() {
            tails_y.push(y);
            tail_idx.push(idx);
        } else {
            tails_y[pos] = y;
            tail_idx[pos] = idx;
        }
        if pos > 0 {
            prev[idx] = tail_idx[pos - 1];
        }
    }

    let mut chain = Vec::with_capacity(tail_idx.len());
    let mut cur = tail_idx[tail_idx.len() - 1];
    while cur != usize::MAX {
        chain.push(sorted[cur]);
        cur = prev[cur];
    }
    chain.reverse();
    chain
}

/// Cluster `(x, y)` hits on one strand into colinear LCB chains.
///
/// Uses `BTreeSet` indices on `(lx, id)` for `max_gap` expiration and `(ly, id)`
/// for $O(\log |\text{active}|)$ local `y`-window range queries instead of
/// scanning all active chains across the chromosome.
#[must_use]
pub fn cluster_strand_lcbs(
    hits: &[(usize, usize)],
    max_gap: usize,
    max_diag_drift: usize,
    min_anchors: usize,
    min_span_bp: usize,
    max_mean_step: f64,
) -> Vec<Vec<(usize, usize)>> {
    use std::collections::BTreeSet;

    if hits.is_empty() {
        return Vec::new();
    }
    let mut sorted = hits.to_vec();
    sorted.sort_unstable();

    let keep = |ch: &[(usize, usize)]| -> bool {
        let span_x = ch[ch.len() - 1].0 - ch[0].0;
        #[allow(clippy::cast_precision_loss)]
        let mean_step = (span_x as f64) / (ch.len() as f64);
        ch.len() >= min_anchors && span_x >= min_span_bp && mean_step <= max_mean_step
    };

    let mut chains: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut tails: Vec<(usize, usize)> = Vec::new();
    let mut by_x: BTreeSet<(usize, usize)> = BTreeSet::new();
    let mut by_y: BTreeSet<(usize, usize)> = BTreeSet::new();
    let mut finished: Vec<Vec<(usize, usize)>> = Vec::new();

    for &(x, y) in &sorted {
        let x_floor = x.saturating_sub(max_gap);
        while let Some(&(lx, id)) = by_x.first() {
            if lx >= x_floor {
                break;
            }
            by_x.remove(&(lx, id));
            let (_, ly) = tails[id];
            by_y.remove(&(ly, id));
            let ch = std::mem::take(&mut chains[id]);
            if keep(&ch) {
                finished.push(ch);
            }
        }

        let y_floor = y.saturating_sub(max_gap);
        let mut best_id: Option<usize> = None;
        let mut best_cost = usize::MAX;
        for &(_, id) in by_y.range((y_floor, 0)..(y, 0)) {
            let (lx, ly) = tails[id];
            if x > lx && y > ly {
                let dx = x - lx;
                let dy = y - ly;
                let drift = dx.abs_diff(dy);
                if dx <= max_gap && dy <= max_gap && drift <= max_diag_drift {
                    let cost = drift + dx;
                    if cost < best_cost {
                        best_cost = cost;
                        best_id = Some(id);
                    }
                }
            }
        }

        if let Some(id) = best_id {
            let (prev_x, prev_y) = tails[id];
            by_x.remove(&(prev_x, id));
            by_y.remove(&(prev_y, id));
            chains[id].push((x, y));
            tails[id] = (x, y);
            by_x.insert((x, id));
            by_y.insert((y, id));
        } else {
            let id = chains.len();
            chains.push(vec![(x, y)]);
            tails.push((x, y));
            by_x.insert((x, id));
            by_y.insert((y, id));
        }
    }

    for (_, id) in by_x {
        let ch = std::mem::take(&mut chains[id]);
        if keep(&ch) {
            finished.push(ch);
        }
    }
    finished
}

/// Find local `w=9` spaced-seed anchors with 2-hit diagonal consistency.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if `mask` is invalid.
pub fn local_w9_subchain(
    sub_x: &[u8],
    sub_y: &[u8],
    off_x: usize,
    off_y: usize,
    mask: &str,
    max_diag_diff: usize,
    max_step_bp: usize,
) -> Result<Vec<(usize, usize)>, AlignError> {
    let ix = canonical_spaced_seeds(sub_x, mask, 2)?;
    let iy = canonical_spaced_seeds(sub_y, mask, 2)?;
    let mut hits = Vec::new();
    for (k, pxs) in &ix {
        if let Some(pys) = iy.get(k) {
            for &px in pxs {
                for &py in pys {
                    hits.push((px, py));
                }
            }
        }
    }
    let chain = lis_chain(&hits);
    let mut kept = Vec::new();

    for (idx, &(x, y)) in chain.iter().enumerate() {
        #[allow(clippy::cast_possible_wrap)]
        let diag = (x as isize) - (y as isize);
        #[allow(clippy::cast_possible_wrap)]
        let max_dd = max_diag_diff as isize;
        let has_prev = idx > 0 && {
            let (px, py) = chain[idx - 1];
            #[allow(clippy::cast_possible_wrap)]
            let pdiag = (px as isize) - (py as isize);
            (diag - pdiag).abs() <= max_dd && (x - px) <= max_step_bp
        };
        let has_next = idx + 1 < chain.len() && {
            let (nx, ny) = chain[idx + 1];
            #[allow(clippy::cast_possible_wrap)]
            let ndiag = (nx as isize) - (ny as isize);
            (diag - ndiag).abs() <= max_dd && (nx - x) <= max_step_bp
        };
        if has_prev || has_next {
            kept.push((off_x + x, off_y + y));
        }
    }
    Ok(kept)
}

pub(crate) fn merge_colinear_lcbs(
    lcb_list: &[Vec<(usize, usize)>],
    max_local_rescue_bp: usize,
) -> Vec<Vec<(usize, usize)>> {
    if lcb_list.is_empty() {
        return Vec::new();
    }
    let mut ordered = lcb_list.to_vec();
    ordered.sort_by_key(|c| c[0].0);
    let mut merged = vec![ordered[0].clone()];
    for nxt in ordered.into_iter().skip(1) {
        let last = &merged[merged.len() - 1];
        let (px, py) = last[last.len() - 1];
        let (nx, ny) = nxt[0];
        if nx > px && ny > py {
            let dx = nx - px;
            let dy = ny - py;
            if dx <= max_local_rescue_bp && dy <= max_local_rescue_bp && dx.abs_diff(dy) <= 4000 {
                let m_idx = merged.len() - 1;
                merged[m_idx].extend(nxt);
                continue;
            }
        }
        merged.push(nxt);
    }
    merged
}

pub(crate) use flat_index::FlatSeedIndex;

/// Precomputed strand-symmetric KP-14 seed index for one contig.
#[derive(Debug, Clone)]
pub struct ContigIndex {
    pub(crate) seq: Vec<u8>,
    pub(crate) seq_rc: Vec<u8>,
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
    let max_occ = if seq.len() >= 100_000 { 2 } else { 4 };
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
        seq_rc: revcomp(seq),
        sdust,
        seeds,
    })
}

/// Count shared KP-14 seeds (on either strand) between two pre-indexed contigs.
#[must_use]
pub fn shared_seed_count(idx_x: &ContigIndex, idx_y: &ContigIndex) -> usize {
    idx_x.seeds.shared_unique_keys(&idx_y.seeds)
}

type LcbList = Vec<Vec<(usize, usize)>>;

type AnchorPairs = Vec<(usize, usize)>;

/// Disjoint closed-interval set `[sx, ex]` backed by `BTreeMap` for `O(log K)` overlap checks.
pub(crate) struct IntervalSet {
    map: std::collections::BTreeMap<usize, usize>,
}

impl IntervalSet {
    pub(crate) fn new() -> Self {
        Self {
            map: std::collections::BTreeMap::new(),
        }
    }

    pub(crate) fn from_lcbs<'a, I>(lcbs: I) -> Self
    where
        I: IntoIterator<Item = &'a Vec<(usize, usize)>>,
    {
        let mut set = Self::new();
        for c in lcbs {
            if let (Some(first), Some(last)) = (c.first(), c.last()) {
                set.insert_merge(first.0, last.0);
            }
        }
        set
    }

    fn insert_merge(&mut self, mut sx: usize, mut ex: usize) {
        while let Some((&s, &e)) = self.map.range(..=ex).next_back() {
            if e < sx {
                break;
            }
            self.map.remove(&s);
            sx = sx.min(s);
            ex = ex.max(e);
        }
        self.map.insert(sx, ex);
    }

    pub(crate) fn insert_if_disjoint(&mut self, sx: usize, ex: usize) -> bool {
        if let Some((_, &e)) = self.map.range(..=ex).next_back()
            && e >= sx
        {
            return false;
        }
        self.map.insert(sx, ex);
        true
    }
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
    let mut keys = Vec::with_capacity(w * w / 2 + w + 1);
    for i in 0..=(seq_x.len() - span) {
        let sub = &seq_x[i..i + span];
        let mut fwd = 0_u64;
        let mut rc = 0_u64;
        let mut valid = true;
        for (&fo, &ro) in offsets.iter().zip(rev_offsets) {
            if skip_softmasked && (sub[fo].is_ascii_lowercase() || sub[ro].is_ascii_lowercase()) {
                valid = false;
                break;
            }
            let (Some(bf), Some(br)) = (base_to_idx(sub[fo]), base_to_idx(sub[ro])) else {
                valid = false;
                break;
            };
            fwd = (fwd << 2) | (bf as u64);
            rc = (rc << 2) | ((3 - br) as u64);
        }
        if !valid {
            continue;
        }
        keys.clear();
        if min_transitions == 0 {
            keys.push(fwd.min(rc));
        }
        if min_transitions <= 1 && max_transitions >= 1 {
            for b in 0..w {
                let mf = fwd ^ (2_u64 << (2 * (w - 1 - b)));
                let mr = rc ^ (2_u64 << (2 * b));
                keys.push(mf.min(mr));
            }
        }
        if max_transitions >= 2 {
            for a in 0..w {
                let ma_f = fwd ^ (2_u64 << (2 * (w - 1 - a)));
                let ma_r = rc ^ (2_u64 << (2 * a));
                for b in (a + 1)..w {
                    let mf = ma_f ^ (2_u64 << (2 * (w - 1 - b)));
                    let mr = ma_r ^ (2_u64 << (2 * b));
                    keys.push(mf.min(mr));
                }
            }
        }
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
    if seq_x.len() < span {
        return (Vec::new(), Vec::new());
    }
    let n_pos = seq_x.len() - span + 1;
    let scan_range = |start_i: usize, end_i: usize| -> (AnchorPairs, AnchorPairs) {
        let w = offsets.len();
        let mut hits_fwd = Vec::new();
        let mut hits_rc = Vec::new();
        let mut keys = Vec::with_capacity(w * w / 2 + w + 1);
        for i in start_i..end_i {
            if let Some(m) = sdust_x
                && (m[i] || m[i + span / 2] || m[i + span - 1])
            {
                continue;
            }
            let sub = &seq_x[i..i + span];
            if sdust_x.is_none() && span >= 12 && is_low_complexity_window(sub) {
                continue;
            }
            let mut fwd = 0_u64;
            let mut rc = 0_u64;
            let mut valid = true;
            for (&fo, &ro) in offsets.iter().zip(rev_offsets) {
                if skip_softmasked && (sub[fo].is_ascii_lowercase() || sub[ro].is_ascii_lowercase())
                {
                    valid = false;
                    break;
                }
                let (Some(bf), Some(br)) = (base_to_idx(sub[fo]), base_to_idx(sub[ro])) else {
                    valid = false;
                    break;
                };
                fwd = (fwd << 2) | (bf as u64);
                rc = (rc << 2) | ((3 - br) as u64);
            }
            if !valid {
                continue;
            }
            keys.clear();
            if min_transitions == 0 {
                keys.push(fwd.min(rc));
            }
            if min_transitions <= 1 && max_transitions >= 1 {
                for b in 0..w {
                    let mf = fwd ^ (2_u64 << (2 * (w - 1 - b)));
                    let mr = rc ^ (2_u64 << (2 * b));
                    keys.push(mf.min(mr));
                }
            }
            if max_transitions >= 2 {
                for a in 0..w {
                    let ma_f = fwd ^ (2_u64 << (2 * (w - 1 - a)));
                    let ma_r = rc ^ (2_u64 << (2 * a));
                    for b in (a + 1)..w {
                        let mf = ma_f ^ (2_u64 << (2 * (w - 1 - b)));
                        let mr = ma_r ^ (2_u64 << (2 * b));
                        keys.push(mf.min(mr));
                    }
                }
            }
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
    let mut promoted: Vec<Vec<(usize, usize)>> = Vec::new();
    for i in 0..sorted.len() {
        let (x1, y1) = sorted[i];
        let mut j = i + 1;
        while j < sorted.len() && sorted[j].0 - x1 <= max_gap {
            let (x2, y2) = sorted[j];
            if x2 > x1 && y2 > y1 && y2 - y1 <= max_gap {
                let d1 = x1.cast_signed() - y1.cast_signed();
                let d2 = x2.cast_signed() - y2.cast_signed();
                if (d1 - d2).unsigned_abs() <= max_ddiag {
                    promoted.push(vec![(x1, y1), (x2, y2)]);
                    break;
                }
            }
            j += 1;
        }
    }
    let mut intervals = IntervalSet::from_lcbs(target.iter());
    for lcb in promoted {
        let (sx, ex) = (lcb[0].0, lcb[lcb.len() - 1].0);
        if intervals.insert_if_disjoint(sx, ex) {
            target.push(lcb);
        }
    }
}

fn codon_to_aa(b0: u8, b1: u8, b2: u8) -> Option<u8> {
    // Standard genetic code indexed by (A=0, C=1, G=2, T=3): 0..20 = amino acids, 255 = stop.
    const TABLE: [u8; 64] = [
        // AAA..ATT
        8, 11, 8, 11, 16, 16, 16, 16, 14, 15, 14, 15, 7, 7, 10, 7, // CAA..CTT
        13, 6, 13, 6, 12, 12, 12, 12, 14, 14, 14, 14, 9, 9, 9, 9, // GAA..GTT
        3, 2, 3, 2, 0, 0, 0, 0, 5, 5, 5, 5, 17, 17, 17, 17, // TAA..TTT
        255, 19, 255, 19, 15, 15, 15, 15, 255, 1, 18, 1, 9, 4, 9, 4,
    ];
    let (Some(i0), Some(i1), Some(i2)) = (base_to_idx(b0), base_to_idx(b1), base_to_idx(b2)) else {
        return None;
    };
    let aa = TABLE[(i0 << 4) | (i1 << 2) | i2];
    if aa == 255 { None } else { Some(aa) }
}

mod flat_index;
mod translated;
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

/// Select top-`K` (or reciprocal-best when `top_k == 1`) contig pairs between two assemblies.
#[must_use]
pub fn select_top_k_pairs(
    target_idx: &[(String, ContigIndex)],
    query_idx: &[(String, ContigIndex)],
    top_k: usize,
    min_shared_seeds: usize,
) -> Vec<Vec<bool>> {
    use rayon::prelude::*;
    let scores: Vec<((usize, usize), (usize, usize))> = (0..target_idx.len())
        .into_par_iter()
        .flat_map(|x_i| {
            (0..query_idx.len()).into_par_iter().map(move |y_i| {
                let seeds = shared_seed_count(&target_idx[x_i].1, &query_idx[y_i].1);
                let lcb = if seeds >= min_shared_seeds {
                    score_contig_pair_lcb(&target_idx[x_i].1, &query_idx[y_i].1)
                } else {
                    0
                };
                ((x_i, y_i), (lcb, seeds))
            })
        })
        .collect();
    let mut by_x: Vec<Vec<((usize, usize), usize)>> = vec![Vec::new(); target_idx.len()];
    let mut by_y: Vec<Vec<((usize, usize), usize)>> = vec![Vec::new(); query_idx.len()];
    for &((x_i, y_i), sc) in &scores {
        if sc.1 >= min_shared_seeds {
            by_x[x_i].push((sc, y_i));
            by_y[y_i].push((sc, x_i));
        }
    }
    let mut mask = vec![vec![false; query_idx.len()]; target_idx.len()];
    for (x_i, list) in by_x.iter_mut().enumerate() {
        list.sort_unstable_by_key(|a| std::cmp::Reverse(a.0));
        for &(_, y_i) in list.iter().take(top_k) {
            mask[x_i][y_i] = true;
        }
    }
    for (y_i, list) in by_y.iter_mut().enumerate() {
        list.sort_unstable_by_key(|a| std::cmp::Reverse(a.0));
        for &(_, x_i) in list.iter().take(top_k) {
            mask[x_i][y_i] = true;
        }
    }
    mask
}

/// Align two pre-indexed contigs without recomputing KP-14 seed maps.
///
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

/// Reverse-complement a DNA byte slice (uppercasing ASCII bases).
#[must_use]
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => other,
        })
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
}
