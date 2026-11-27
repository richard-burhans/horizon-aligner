//! Contig spaced-seed indexing, `MinHash` sketching, and top-k pair selection.

use super::{
    FlatSeedIndex, MASK_PRIMARY_11_15, cluster_strand_lcbs, find_w12_raw_lcbs, sdust_mask,
};
use crate::{AlignError, splitmix64};
use rayon::prelude::*;

/// Precomputed strand-symmetric KP-14 seed index for one contig.
#[derive(Debug, Clone)]
pub struct ContigIndex {
    pub(crate) seq: Vec<u8>,
    pub(crate) sdust: Vec<bool>,
    pub(crate) seeds: FlatSeedIndex,
}

/// Precompute forward and reverse-complement spaced-seed maps for a contig.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if `MASK_PRIMARY_11_15` is invalid.
pub fn prepare_contig_index(seq: &[u8]) -> Result<ContigIndex, AlignError> {
    prepare_contig_index_owned(seq.to_vec(), false)
}

/// Precompute forward and reverse-complement KP-14 seed maps by taking ownership of `seq`
/// without duplicating the raw DNA buffer.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if `MASK_PRIMARY_11_15` is invalid.
pub fn prepare_contig_index_owned(
    seq: Vec<u8>,
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
    let sdust = sdust_mask(&seq);
    let seeds = FlatSeedIndex::build_with_sdust(
        &seq,
        MASK_PRIMARY_11_15,
        max_occ,
        skip_softmasked,
        false,
        Some(&sdust),
    )?;
    Ok(ContigIndex { seq, sdust, seeds })
}

/// Count shared KP-14 seeds (on either strand) between two pre-indexed contigs.
#[must_use]
pub fn shared_seed_count(idx_x: &ContigIndex, idx_y: &ContigIndex) -> usize {
    idx_x.seeds.shared_unique_keys(&idx_y.seeds)
}

/// Build an assembly-wide contig document-frequency table (`df[key]` = number of contigs
/// containing `key`) for fragmented draft assemblies (`>= 32` contigs) so dispersed repeats
/// appearing across many scaffolds are excluded during `--top-k` pre-screening.
fn build_contig_key_df(contigs: &[(String, ContigIndex)]) -> Option<(Vec<u8>, u8)> {
    const MIN_DRAFT_CONTIGS: usize = 32;
    const MAX_DF_TABLE: usize = 1 << 22; // 4^11 = 4,194,304 (4 MiB for weight-11 KP-14 keys)
    if contigs.len() < MIN_DRAFT_CONTIGS {
        return None;
    }
    let mut max_key = 0usize;
    for (_, ix) in contigs {
        if let Some(&last) = ix.seeds.keys.last()
            && let Ok(v) = usize::try_from(last)
            && v < MAX_DF_TABLE
        {
            max_key = max_key.max(v);
        }
    }
    if max_key == 0 {
        return None;
    }
    let mut df = vec![0u8; max_key + 1];
    for (_, ix) in contigs {
        let mut prev = u64::MAX;
        for &k in &ix.seeds.keys {
            if k == prev {
                continue;
            }
            prev = k;
            if let Ok(idx) = usize::try_from(k)
                && idx < df.len()
            {
                df[idx] = df[idx].saturating_add(1);
            }
        }
    }
    #[allow(clippy::cast_possible_truncation)]
    let max_df = (contigs.len() / 20).clamp(4, 32) as u8;
    Some((df, max_df))
}

impl ContigIndex {
    /// Borrow the owned raw DNA sequence stored inside this contig index.
    #[must_use]
    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    /// Consume this index, freeing the `sdust` mask and seed tables while returning the owned DNA.
    #[must_use]
    pub fn into_seq(self) -> Vec<u8> {
        self.seq
    }

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

/// Estimate Jaccard similarity between two sorted bottom-`k` `MinHash` sketches using Broder's
/// unbiased `MIN_k(A ∪ B)` union estimator (Broder, 1997) in `O(k)` time.
///
/// Assuming both sketches were constructed with the same cap `k`, merging up to
/// `k_cap = max(|a|, |b|)` distinct smallest hashes from `a ∪ b` recovers `MIN_{k_cap}(U_A ∪ U_B)`
/// identically for equal-size saturated sketches and asymmetric small-vs-large sets alike.
#[must_use]
pub fn sketch_jaccard(a: &[u64], b: &[u64]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    let k_cap = a.len().max(b.len());
    let mut i = 0usize;
    let mut j = 0usize;
    let mut union_seen = 0usize;
    let mut shared = 0usize;
    while union_seen < k_cap && (i < a.len() || j < b.len()) {
        match (a.get(i), b.get(j)) {
            (Some(&va), Some(&vb)) => match va.cmp(&vb) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    shared += 1;
                    i += 1;
                    j += 1;
                }
            },
            (Some(_), None) => i += 1,
            (None, Some(_)) => j += 1,
            (None, None) => break,
        }
        union_seen += 1;
    }
    if union_seen == 0 {
        0.0
    } else {
        #[allow(clippy::cast_precision_loss)]
        {
            (shared as f64) / (union_seen as f64)
        }
    }
}

/// Fast pre-DP score: return `(total_lcb_anchors, max_lcb_anchors)` across collinear `w=12` LCBs
/// between two contigs. Summing across all collinear LCBs prevents rearranged multi-LCB homologs
/// from losing to single-LCB paralogs, and a tight 3-anchor / 35-bp fallback on short draft
/// contigs (`min(len_x, len_y) <= 4_000`) rescues compact fragmented scaffolds without admitting
/// isolated 20-bp repeat hits.
type ScoredPartnerList = Vec<((usize, usize, usize), usize)>;

fn score_contig_pair_lcb(idx_x: &ContigIndex, idx_y: &ContigIndex) -> (usize, usize) {
    let (raw_fwd, raw_rc) = find_w12_raw_lcbs(idx_x, idx_y);
    let mut total = 0usize;
    let mut best = 0usize;
    for lcb in raw_fwd.iter().chain(raw_rc.iter()) {
        let len = lcb.len();
        total = total.saturating_add(len);
        if len > best {
            best = len;
        }
    }
    if total == 0 && idx_x.seq().len().min(idx_y.seq().len()) <= 4_000 {
        let (hits_fwd, hits_rc) = idx_x.seeds.intersect(&idx_y.seeds);
        let short_fwd = cluster_strand_lcbs(&hits_fwd, 500, 40, 3, 35, 100.0);
        let short_rc = cluster_strand_lcbs(&hits_rc, 500, 40, 3, 35, 100.0);
        for lcb in short_fwd.iter().chain(short_rc.iter()) {
            let len = lcb.len();
            total = total.saturating_add(len);
            if len > best {
                best = len;
            }
        }
    }
    (total, best)
}

/// Select up to `top_k` query contigs per target contig (plus reciprocal best targets for
/// multi-contig assemblies when collinear LCB support exists).
///
/// Crucially, when one side has `<= top_k` contigs (e.g. a single reference chromosome or 2 Mb
/// locus aligned against a 20,000-scaffold draft genome), reciprocal selection on the small side
/// cannot discriminate targets and is skipped so `--top-k K` never degrades to `--all-pairs`.
/// For each shortlisted target contig, candidates are drawn from the union of top raw seed counts
/// and top seed-density counts (`seeds / len_y`) so both large homologous chromosomes and compact
/// homologous draft scaffolds reach LCB scoring.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn select_top_k_pairs(
    target_idx: &[(String, ContigIndex)],
    query_idx: &[(String, ContigIndex)],
    top_k: usize,
    min_shared_seeds: usize,
) -> Vec<Vec<bool>> {
    if target_idx.is_empty() || query_idx.is_empty() || top_k == 0 {
        return vec![vec![false; query_idx.len()]; target_idx.len()];
    }

    let query_df = build_contig_key_df(query_idx);
    let seed_counts: Vec<Vec<usize>> = (0..target_idx.len())
        .into_par_iter()
        .map(|x_i| {
            (0..query_idx.len())
                .map(|y_i| match &query_df {
                    Some((df, max_df)) => target_idx[x_i].1.seeds.shared_non_repetitive_keys(
                        &query_idx[y_i].1.seeds,
                        df,
                        *max_df,
                    ),
                    None => shared_seed_count(&target_idx[x_i].1, &query_idx[y_i].1),
                })
                .collect()
        })
        .collect();

    let multi_both = target_idx.len() > 1 && query_idx.len() > 1;
    let use_by_x = query_idx.len() > top_k || target_idx.len() <= top_k || multi_both;
    let use_by_y = target_idx.len() > top_k || (multi_both && query_idx.len() > top_k);
    let require_lcb_x = multi_both && target_idx.len() > top_k && query_idx.len() <= top_k;
    let require_lcb_y = multi_both && query_idx.len() > top_k;

    let shortlist_k = top_k.saturating_mul(4).max(8);
    let mut need_lcb = vec![vec![false; query_idx.len()]; target_idx.len()];

    if use_by_x {
        let mut row_raw: Vec<(usize, usize)> = Vec::new();
        let mut row_dens: Vec<(u64, usize)> = Vec::new();
        for (x_i, row) in seed_counts.iter().enumerate() {
            row_raw.clear();
            row_dens.clear();
            for (y_i, &seeds) in row.iter().enumerate() {
                if seeds >= min_shared_seeds {
                    row_raw.push((seeds, y_i));
                    let denom = query_idx[y_i].1.seq().len().max(1_000) as u64;
                    let dens = (seeds as u64).saturating_mul(1_000_000) / denom;
                    row_dens.push((dens, y_i));
                }
            }
            row_raw.sort_unstable_by_key(|&(s, y_i)| (std::cmp::Reverse(s), y_i));
            for &(_, y_i) in row_raw.iter().take(shortlist_k) {
                need_lcb[x_i][y_i] = true;
            }
            row_dens.sort_unstable_by_key(|&(d, y_i)| (std::cmp::Reverse(d), y_i));
            for &(_, y_i) in row_dens.iter().take(shortlist_k) {
                need_lcb[x_i][y_i] = true;
            }
        }
    }

    if use_by_y {
        let mut col_raw: Vec<(usize, usize)> = Vec::new();
        let mut col_dens: Vec<(u64, usize)> = Vec::new();
        for y_i in 0..query_idx.len() {
            col_raw.clear();
            col_dens.clear();
            for (x_i, row) in seed_counts.iter().enumerate() {
                let seeds = row[y_i];
                if seeds >= min_shared_seeds {
                    col_raw.push((seeds, x_i));
                    let denom = target_idx[x_i].1.seq().len().max(1_000) as u64;
                    let dens = (seeds as u64).saturating_mul(1_000_000) / denom;
                    col_dens.push((dens, x_i));
                }
            }
            col_raw.sort_unstable_by_key(|&(s, x_i)| (std::cmp::Reverse(s), x_i));
            for &(_, x_i) in col_raw.iter().take(shortlist_k) {
                need_lcb[x_i][y_i] = true;
            }
            col_dens.sort_unstable_by_key(|&(d, x_i)| (std::cmp::Reverse(d), x_i));
            for &(_, x_i) in col_dens.iter().take(shortlist_k) {
                need_lcb[x_i][y_i] = true;
            }
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
    let lcb_scores: Vec<((usize, usize), (usize, usize))> = lcb_jobs
        .into_par_iter()
        .map(|(x_i, y_i)| {
            let lcb = score_contig_pair_lcb(&target_idx[x_i].1, &query_idx[y_i].1);
            ((x_i, y_i), lcb)
        })
        .collect();

    let mut by_x: Vec<ScoredPartnerList> = vec![Vec::new(); target_idx.len()];
    let mut by_y: Vec<ScoredPartnerList> = vec![Vec::new(); query_idx.len()];
    for ((x_i, y_i), (lcb_total, lcb_max)) in lcb_scores {
        let seeds = seed_counts[x_i][y_i];
        if use_by_x && (!require_lcb_x || lcb_total > 0) {
            by_x[x_i].push(((lcb_total, lcb_max, seeds), y_i));
        }
        // In symmetric multi-contig assemblies, only let reciprocal `by_y` attach query scaffolds
        // that have actual collinear LCB support (`lcb_total > 0`) so thousands of unrelated draft
        // scaffolds with random 11-mer collisions do not explode the alignment queue.
        if use_by_y && (!require_lcb_y || lcb_total > 0) {
            by_y[y_i].push(((lcb_total, lcb_max, seeds), x_i));
        }
    }
    let mut mask = vec![vec![false; query_idx.len()]; target_idx.len()];
    if use_by_x {
        for (x_i, list) in by_x.iter_mut().enumerate() {
            list.sort_unstable_by_key(|a| (std::cmp::Reverse(a.0), a.1));
            for &(_, y_i) in list.iter().take(top_k) {
                mask[x_i][y_i] = true;
            }
        }
    }
    if use_by_y {
        for (y_i, list) in by_y.iter_mut().enumerate() {
            list.sort_unstable_by_key(|a| (std::cmp::Reverse(a.0), a.1));
            for &(_, x_i) in list.iter().take(top_k) {
                mask[x_i][y_i] = true;
            }
        }
    }
    mask
}

/// Select up to `top_k` query contigs for each borrowed target `ContigIndex` (strictly
/// one-directional, bounded by `targets.len() * top_k`, zero index cloning).
#[must_use]
pub fn select_top_k_for_target_refs(
    targets: &[&ContigIndex],
    query_idx: &[(String, ContigIndex)],
    top_k: usize,
    min_shared_seeds: usize,
) -> Vec<(usize, usize)> {
    if targets.is_empty() || query_idx.is_empty() {
        return Vec::new();
    }
    let shortlist_k = top_k.saturating_mul(4).max(8);
    let query_df = build_contig_key_df(query_idx);
    let seed_rows: Vec<Vec<usize>> = targets
        .par_iter()
        .map(|&target| {
            (0..query_idx.len())
                .map(|y_i| match &query_df {
                    Some((df, max_df)) => target.seeds.shared_non_repetitive_keys(
                        &query_idx[y_i].1.seeds,
                        df,
                        *max_df,
                    ),
                    None => shared_seed_count(target, &query_idx[y_i].1),
                })
                .collect()
        })
        .collect();

    if top_k == 0 {
        let mut pairs = Vec::new();
        for (row_idx, counts) in seed_rows.into_iter().enumerate() {
            for (y_i, seeds) in counts.into_iter().enumerate() {
                if seeds >= min_shared_seeds {
                    pairs.push((row_idx, y_i));
                }
            }
        }
        return pairs;
    }

    let mut lcb_jobs: Vec<(usize, usize, usize)> = Vec::new();
    let mut row_raw: Vec<(usize, usize)> = Vec::new();
    let mut row_dens: Vec<(u64, usize)> = Vec::new();
    let mut seen_y: Vec<bool> = vec![false; query_idx.len()];
    let mut picked_y: Vec<usize> = Vec::new();
    for (row_idx, counts) in seed_rows.iter().enumerate() {
        row_raw.clear();
        row_dens.clear();
        for (y_i, &seeds) in counts.iter().enumerate() {
            if seeds >= min_shared_seeds {
                row_raw.push((seeds, y_i));
                let denom = query_idx[y_i].1.seq().len().max(1_000) as u64;
                let dens = (seeds as u64).saturating_mul(1_000_000) / denom;
                row_dens.push((dens, y_i));
            }
        }
        row_raw.sort_unstable_by_key(|&(s, y_i)| (std::cmp::Reverse(s), y_i));
        row_dens.sort_unstable_by_key(|&(d, y_i)| (std::cmp::Reverse(d), y_i));
        picked_y.clear();
        for (_, y_i) in row_raw
            .iter()
            .take(shortlist_k)
            .map(|&(s, y)| (s as u64, y))
            .chain(row_dens.iter().take(shortlist_k).copied())
        {
            if !seen_y[y_i] {
                seen_y[y_i] = true;
                picked_y.push(y_i);
                lcb_jobs.push((row_idx, y_i, counts[y_i]));
            }
        }
        for &y_i in &picked_y {
            seen_y[y_i] = false;
        }
    }

    let scored: Vec<(usize, usize, (usize, usize), usize)> = lcb_jobs
        .into_par_iter()
        .map(|(row_idx, y_i, seeds)| {
            let lcb = score_contig_pair_lcb(targets[row_idx], &query_idx[y_i].1);
            (row_idx, y_i, lcb, seeds)
        })
        .collect();

    let mut by_row: Vec<ScoredPartnerList> = vec![Vec::new(); targets.len()];
    for (row_idx, y_i, (lcb_total, lcb_max), seeds) in scored {
        by_row[row_idx].push(((lcb_total, lcb_max, seeds), y_i));
    }

    let mut out = Vec::new();
    for (row_idx, list) in by_row.iter_mut().enumerate() {
        list.sort_unstable_by_key(|a| (std::cmp::Reverse(a.0), a.1));
        for &(_, y_i) in list.iter().take(top_k) {
            out.push((row_idx, y_i));
        }
    }
    out
}

/// Indexed multi-contig genome pairing a species name with its owned [`ContigIndex`]es.
#[derive(Debug, Clone)]
pub struct IndexedGenome {
    /// Species or assembly identifier.
    pub name: String,
    /// Per-contig `(contig_name, ContigIndex)` entries.
    pub contigs: Vec<(String, ContigIndex)>,
}

impl IndexedGenome {
    /// Build a genome-wide bottom-`k` `MinHash` sketch across all contigs.
    #[must_use]
    pub fn minhash_sketch(&self, k: usize) -> Vec<u64> {
        merge_contig_sketches(&self.contigs, k)
    }

    /// Borrow `(contig_name, seq)` slices without duplicating DNA buffers.
    #[must_use]
    pub fn contig_refs(&self) -> Vec<(&str, &[u8])> {
        self.contigs
            .iter()
            .map(|(n, ix)| (n.as_str(), ix.seq()))
            .collect()
    }

    /// Effective contig count `n_eff = (sum L_i)^2 / sum(L_i^2) = total_bp / auN` (Li 2020):
    /// identically `1.0` for single-contig genomes, resilient to minor unplaced scaffolds,
    /// and equal to `K` for `K` equal-length draft contigs.
    #[must_use]
    pub fn effective_contig_count(&self) -> f64 {
        effective_contig_count(self.contigs.iter().map(|(_, ix)| ix.seq().len()))
    }
}

/// Effective contig count `n_eff = (sum L_i)^2 / sum(L_i^2) = total_bp / auN` (Li 2020).
#[must_use]
pub fn effective_contig_count(lengths: impl Iterator<Item = usize>) -> f64 {
    let mut sum_l = 0.0_f64;
    let mut sum_l2 = 0.0_f64;
    for len_bp in lengths {
        #[allow(clippy::cast_precision_loss)]
        let len = len_bp as f64;
        sum_l += len;
        sum_l2 += len * len;
    }
    if sum_l2 <= 0.0 {
        1.0
    } else {
        ((sum_l * sum_l) / sum_l2).max(1.0)
    }
}

/// Scale a base `--top-k` contig-pair budget adaptively when aligning against an ultra-fragmented
/// draft assembly (`n_eff >= 50`), preventing chromosome-scale references from starving thousands
/// of draft scaffolds while leaving chromosome-scale pairs (`n_eff < 50`) completely unchanged.
#[must_use]
pub fn adaptive_top_k_for_fragmentation(
    base_top_k: usize,
    target_idx: &[(String, ContigIndex)],
    query_idx: &[(String, ContigIndex)],
) -> usize {
    if base_top_k == 0 || base_top_k >= 64 || target_idx.is_empty() || query_idx.is_empty() {
        return base_top_k;
    }
    let n_eff_x = effective_contig_count(target_idx.iter().map(|(_, ix)| ix.seq().len()));
    let n_eff_y = effective_contig_count(query_idx.iter().map(|(_, ix)| ix.seq().len()));
    let n_max = n_eff_x.max(n_eff_y);
    let n_min = n_eff_x.min(n_eff_y).max(1.0);
    let ratio = n_max / n_min;
    if n_max < 50.0 || (n_max < 200.0 && ratio < 4.0) {
        return base_top_k;
    }
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    let scaled = ((base_top_k as f64) * ratio.sqrt().max(1.0)).ceil() as usize;
    let max_cap = base_top_k.saturating_mul(16).clamp(32, 64).max(base_top_k);
    scaled.clamp(base_top_k, max_cap)
}

/// Merge per-contig bottom-`k` `MinHash` sketches for a species into a single sorted sketch of up to
/// `k` smallest seed hashes.
#[must_use]
pub fn merge_contig_sketches(indices: &[(String, ContigIndex)], k: usize) -> Vec<u64> {
    use std::collections::BTreeSet;
    let mut set: BTreeSet<u64> = BTreeSet::new();
    for (_, ix) in indices {
        for h in ix.minhash_sketch(k) {
            if set.len() < k {
                set.insert(h);
            } else if let Some(&top) = set.last()
                && h < top
                && set.insert(h)
            {
                let _ = set.pop_last();
            }
        }
    }
    set.into_iter().collect()
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn prng_dna(len: usize, mut seed: u64) -> Vec<u8> {
        const BASES: [u8; 4] = *b"ACGT";
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            out.push(BASES[((seed >> 32) & 3) as usize]);
        }
        out
    }

    #[test]
    fn select_top_k_pairs_bounds_asymmetric_single_target_queries() {
        let target_seq = prng_dna(1_200, 11);
        let target = vec![(
            "chr1".to_string(),
            prepare_contig_index_owned(target_seq.clone(), false).unwrap(),
        )];
        // Build 10 query contigs that all share >=8 seeds with `target`, with contigs 2 and 7
        // carrying long collinear matches.
        let mut queries = Vec::new();
        for i in 0..10usize {
            let mut q = prng_dna(400, 100 + i as u64);
            if i == 2 || i == 7 {
                q.splice(50..350, target_seq[100..400].iter().copied());
            } else {
                // Plant isolated 20 bp chunks so raw shared_seed_count >= 8 without a big LCB.
                for chunk in 0..3 {
                    let src_off = 500 + chunk * 100;
                    let dst_off = chunk * 80;
                    q[dst_off..dst_off + 20].copy_from_slice(&target_seq[src_off..src_off + 20]);
                }
            }
            queries.push((
                format!("q{i}"),
                prepare_contig_index_owned(q, false).unwrap(),
            ));
        }
        let mask = select_top_k_pairs(&target, &queries, 2, 8);
        let selected: Vec<usize> = mask[0]
            .iter()
            .enumerate()
            .filter_map(|(i, &keep)| keep.then_some(i))
            .collect();
        assert_eq!(
            selected,
            vec![2, 7],
            "single-contig target with top_k=2 must select at most 2 homologous queries, not all 10"
        );
    }

    #[test]
    fn select_top_k_for_target_refs_only_visits_requested_targets() {
        let s0 = prng_dna(800, 31);
        let s1 = prng_dna(800, 32);
        let s2 = prng_dna(800, 33);
        let targets = [
            (
                "t0".to_string(),
                prepare_contig_index_owned(s0.clone(), false).unwrap(),
            ),
            (
                "t1".to_string(),
                prepare_contig_index_owned(s1.clone(), false).unwrap(),
            ),
            (
                "t2".to_string(),
                prepare_contig_index_owned(s2.clone(), false).unwrap(),
            ),
        ];
        let queries = vec![
            (
                "q0".to_string(),
                prepare_contig_index_owned(s0, false).unwrap(),
            ),
            (
                "q1".to_string(),
                prepare_contig_index_owned(s1, false).unwrap(),
            ),
            (
                "q2".to_string(),
                prepare_contig_index_owned(s2, false).unwrap(),
            ),
        ];
        let pairs = select_top_k_for_target_refs(&[&targets[1].1], &queries, 1, 8);
        assert_eq!(pairs, vec![(0, 1)]);
    }

    #[test]
    fn select_top_k_pairs_rescues_compact_draft_contig_in_symmetric_multi_contig() {
        // 3 reference chromosomes vs 3 query contigs where q2 is a compact 60 bp homolog of chr2
        // (too short for the 70 bp w12 LCB threshold, rescued by the short-contig fallback).
        let c0 = prng_dna(1_000, 101);
        let c1 = prng_dna(1_000, 102);
        let c2 = prng_dna(1_000, 103);
        let targets = vec![
            (
                "c0".to_owned(),
                prepare_contig_index_owned(c0.clone(), false).unwrap(),
            ),
            (
                "c1".to_owned(),
                prepare_contig_index_owned(c1.clone(), false).unwrap(),
            ),
            (
                "c2".to_owned(),
                prepare_contig_index_owned(c2.clone(), false).unwrap(),
            ),
        ];
        let mut short_draft = prng_dna(120, 203);
        short_draft[30..90].copy_from_slice(&c2[400..460]);
        let queries = vec![
            (
                "q0".to_owned(),
                prepare_contig_index_owned(c0, false).unwrap(),
            ),
            (
                "q1".to_owned(),
                prepare_contig_index_owned(c1, false).unwrap(),
            ),
            (
                "q2_draft".to_owned(),
                prepare_contig_index_owned(short_draft, false).unwrap(),
            ),
        ];
        let mask = select_top_k_pairs(&targets, &queries, 1, 8);
        assert!(
            mask[2][2],
            "compact 60 bp homologous draft contig must attach to its true reference chromosome"
        );
    }

    #[test]
    fn adaptive_top_k_scales_only_for_fragmented_draft_assemblies() {
        let chr_target: Vec<(String, ContigIndex)> = (0..4u64)
            .map(|i| {
                (
                    format!("chr{i}"),
                    prepare_contig_index_owned(prng_dna(500, 10 + i), false).unwrap(),
                )
            })
            .collect();
        let chr_query: Vec<(String, ContigIndex)> = (0..4u64)
            .map(|i| {
                (
                    format!("qchr{i}"),
                    prepare_contig_index_owned(prng_dna(500, 20 + i), false).unwrap(),
                )
            })
            .collect();
        // Chromosome-scale pair: base_top_k unchanged.
        assert_eq!(
            adaptive_top_k_for_fragmentation(1, &chr_target, &chr_query),
            1
        );
        assert_eq!(
            adaptive_top_k_for_fragmentation(0, &chr_target, &chr_query),
            0
        );

        // Fragmented 100-scaffold draft against 4 chromosomes: ratio=25 -> sqrt(25)=5.
        let draft_query: Vec<(String, ContigIndex)> = (0..100u64)
            .map(|i| {
                (
                    format!("scaf{i}"),
                    prepare_contig_index_owned(prng_dna(200, 100 + i), false).unwrap(),
                )
            })
            .collect();
        let scaled = adaptive_top_k_for_fragmentation(1, &chr_target, &draft_query);
        assert!(
            (4..=32).contains(&scaled),
            "expected adaptive top_k in [4, 32] for 100-contig draft vs 4 chrs, got {scaled}"
        );
        // Must never panic when base_top_k > 64.
        assert_eq!(
            adaptive_top_k_for_fragmentation(100, &chr_target, &draft_query),
            100
        );
    }

    #[test]
    fn select_top_k_filters_dispersed_repeats_on_fragmented_draft() {
        // Construct a target chromosome containing a 300-bp dispersed TE repeat plus two 350-bp
        // unique single-copy exons, aligned against a 64-scaffold draft where 62 short scaffolds
        // carry the dispersed TE repeat (high document frequency) and only 2 scaffolds carry the
        // true single-copy homologs. Non-repetitive DF filtering must select the 2 true homologs.
        let mut s = 0x1234_5678_9ABC_DEF0_u64;
        let mut next_dna = |len: usize| -> Vec<u8> {
            (0..len)
                .map(|_| {
                    s = splitmix64(s.wrapping_add(1));
                    b"ACGT"[(s & 3) as usize]
                })
                .collect()
        };
        let te_repeat = next_dna(300);
        let hom1 = next_dna(350);
        let hom2 = next_dna(350);
        let mut target_seq = next_dna(1_000);
        target_seq.extend_from_slice(&te_repeat);
        target_seq.extend_from_slice(&next_dna(500));
        target_seq.extend_from_slice(&hom1);
        target_seq.extend_from_slice(&next_dna(500));
        target_seq.extend_from_slice(&hom2);
        target_seq.extend_from_slice(&next_dna(1_000));

        let target = vec![(
            "chr1".to_owned(),
            prepare_contig_index(&target_seq).expect("target"),
        )];

        let mut queries = Vec::new();
        for i in 0..62 {
            // Short repeat-saturated decoy scaffold (very high raw seed density if unfiltered!)
            let mut dec = next_dna(40);
            dec.extend_from_slice(&te_repeat);
            dec.extend_from_slice(&next_dna(40));
            queries.push((
                format!("te_decoy_{i}"),
                prepare_contig_index(&dec).expect("dec"),
            ));
        }
        queries.push((
            "true_hom1".to_owned(),
            prepare_contig_index(&hom1).expect("h1"),
        ));
        queries.push((
            "true_hom2".to_owned(),
            prepare_contig_index(&hom2).expect("h2"),
        ));

        let mask = select_top_k_pairs(&target, &queries, 2, 8);
        assert!(mask[0][62], "must select true_hom1 over 62 TE decoys");
        assert!(mask[0][63], "must select true_hom2 over 62 TE decoys");
        let decoy_hits = (0..62).filter(|&i| mask[0][i]).count();
        assert_eq!(
            decoy_hits, 0,
            "dispersed-repeat decoys must not displace single-copy homologs"
        );

        let ref_pairs = select_top_k_for_target_refs(&[&target[0].1], &queries, 2, 8);
        let picked: Vec<usize> = ref_pairs.into_iter().map(|(_, y)| y).collect();
        assert!(
            picked.contains(&62) && picked.contains(&63),
            "picked={picked:?}"
        );
    }
}
