//! Reference-free bottom-up TBA clade indexing and blockset merging.

use super::model::{PairwiseTbaBlock, TbaBlock, TbaRow, canonical_pair_key};
use super::poa_refine::refine_shared_tip_insertions_poa;
use super::tree::GuideTree;
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};

type CladeSpan = (usize, usize, usize, usize);
type CladeHit = Option<(usize, usize)>;
type ColTag = Option<(CladeHit, CladeHit)>;

/// Sorted non-overlapping interval index mapping `(species, contig)` forward positions to
/// `(block_idx, row_idx)` in a clade's blockset, plus O(1) ungapped-base-to-column tables.
struct CladeIndex {
    /// Nested `species -> contig -> sorted spans (fwd_start, fwd_end, block_idx, row_idx)`
    /// enabling zero-allocation borrowed lookups.
    by_sp_ctg: HashMap<String, HashMap<String, Vec<CladeSpan>>>,
    /// `row_cols[bi][ri][k]` = alignment column index of the `k`-th ungapped base of `blocks[bi].rows[ri]`.
    row_cols: Vec<Vec<Vec<u32>>>,
}

impl CladeIndex {
    fn build(blocks: &[TbaBlock]) -> Self {
        Self::build_for_species(blocks, None)
    }

    /// Build a [`CladeIndex`] optionally restricted to a single species `only_species` (e.g. a
    /// target hinge during tip projection, avoiding `O(N_species)` `row_cols` allocations when
    /// attaching hundreds of tips onto a multi-species backbone).
    fn build_for_species(blocks: &[TbaBlock], only_species: Option<&str>) -> Self {
        let mut by_sp_ctg: HashMap<String, HashMap<String, Vec<CladeSpan>>> = HashMap::new();
        let mut row_cols: Vec<Vec<Vec<u32>>> = Vec::with_capacity(blocks.len());
        for (bi, blk) in blocks.iter().enumerate() {
            let mut blk_cols = Vec::with_capacity(blk.rows.len());
            for (ri, r) in blk.rows.iter().enumerate() {
                if only_species.is_some_and(|target| r.species != target) {
                    blk_cols.push(Vec::new());
                    continue;
                }
                let (s, e) = r.fwd_interval();
                if e > s {
                    by_sp_ctg
                        .entry(r.species.clone())
                        .or_default()
                        .entry(r.contig.clone())
                        .or_default()
                        .push((s, e, bi, ri));
                }
                let mut cols = Vec::with_capacity(r.size);
                for (c, &b) in r.seq.iter().enumerate() {
                    if b != b'-' {
                        cols.push(u32::try_from(c).unwrap_or(u32::MAX));
                    }
                }
                blk_cols.push(cols);
            }
            row_cols.push(blk_cols);
        }
        for ctg_map in by_sp_ctg.values_mut() {
            for ivs in ctg_map.values_mut() {
                ivs.sort_unstable_by_key(|&(s, _, _, _)| s);
            }
        }
        Self {
            by_sp_ctg,
            row_cols,
        }
    }

    fn lookup(&self, species: &str, contig: &str, fwd_pos: usize) -> Option<(usize, usize)> {
        let ivs = self.by_sp_ctg.get(species)?.get(contig)?;
        let idx = ivs.partition_point(|&(s, _, _, _)| s <= fwd_pos);
        if idx > 0 {
            let (_, e, bi, ri) = ivs[idx - 1];
            if fwd_pos < e {
                return Some((bi, ri));
            }
        }
        None
    }

    /// O(1) lookup of the alignment column index for a single forward position `fwd_pos` on
    /// `blocks[bi].rows[ri]`.
    fn fwd_pos_to_col(
        &self,
        blk: &TbaBlock,
        bi: usize,
        ri: usize,
        fwd_pos: usize,
    ) -> Option<usize> {
        let r = blk.rows.get(ri)?;
        let (r_s, r_e) = r.fwd_interval();
        if fwd_pos < r_s || fwd_pos >= r_e {
            return None;
        }
        let k = if r.strand == '+' {
            fwd_pos - r_s
        } else {
            r_e - 1 - fwd_pos
        };
        self.row_cols
            .get(bi)
            .and_then(|b| b.get(ri))
            .and_then(|c| c.get(k))
            .map(|&col| col as usize)
    }

    /// O(1) lookup of the tightest column range `[col_lo, col_hi)` covering `[fwd_s, fwd_e)` on
    /// `blocks[bi].rows[ri]`.
    fn fwd_interval_to_cols(
        &self,
        blk: &TbaBlock,
        bi: usize,
        ri: usize,
        fwd_s: usize,
        fwd_e: usize,
    ) -> Option<(usize, usize)> {
        if fwd_s >= fwd_e {
            return None;
        }
        let r = blk.rows.get(ri)?;
        let (r_s, r_e) = r.fwd_interval();
        if fwd_s < r_s || fwd_e > r_e {
            return None;
        }
        let (k_lo, k_hi) = if r.strand == '+' {
            (fwd_s - r_s, fwd_e - r_s)
        } else {
            (r_e - fwd_e, r_e - fwd_s)
        };
        let cols = self.row_cols.get(bi)?.get(ri)?;
        let col_lo = *cols.get(k_lo)? as usize;
        let col_hi = (*cols.get(k_hi.checked_sub(1)?)? as usize) + 1;
        Some((col_lo, col_hi))
    }
}

fn slice_pairwise_subspan(
    p: &PairwiseTbaBlock,
    col_lo: usize,
    col_hi: usize,
) -> Option<(TbaRow, TbaRow)> {
    let a_skip = p.row_a.seq[..col_lo].iter().filter(|&&b| b != b'-').count();
    let b_skip = p.row_b.seq[..col_lo].iter().filter(|&&b| b != b'-').count();
    let mut a_seq = Vec::with_capacity(col_hi - col_lo);
    let mut b_seq = Vec::with_capacity(col_hi - col_lo);
    let mut a_bases = 0usize;
    let mut b_bases = 0usize;
    for col in col_lo..col_hi {
        let ab = p.row_a.seq[col];
        let bb = p.row_b.seq[col];
        if ab == b'-' && bb == b'-' {
            continue;
        }
        if ab != b'-' {
            a_bases += 1;
        }
        if bb != b'-' {
            b_bases += 1;
        }
        a_seq.push(ab);
        b_seq.push(bb);
    }
    if a_bases == 0 || b_bases == 0 {
        return None;
    }
    Some((
        TbaRow {
            species: p.row_a.species.clone(),
            contig: p.row_a.contig.clone(),
            start: p.row_a.start.saturating_add(a_skip),
            size: a_bases,
            strand: p.row_a.strand,
            src_size: p.row_a.src_size,
            seq: a_seq,
        },
        TbaRow {
            species: p.row_b.species.clone(),
            contig: p.row_b.contig.clone(),
            start: p.row_b.start.saturating_add(b_skip),
            size: b_bases,
            strand: p.row_b.strand,
            src_size: p.row_b.src_size,
            seq: b_seq,
        },
    ))
}

/// Thread `left_blk` (where `left_blk.rows[a_idx]` has the exact ungapped bases of `sub_a`)
/// with `right_blk` (where `right_blk.rows[b_idx]` has the exact ungapped bases of `sub_b`)
/// using the column alignment `(sub_a, sub_b)`.
#[allow(clippy::too_many_lines)]
fn thread_blocks_via_pair(
    left_blk: &TbaBlock,
    a_idx: usize,
    right_blk: &TbaBlock,
    b_idx: usize,
    sub_a: &TbaRow,
    sub_b: &TbaRow,
    score: f64,
) -> Option<TbaBlock> {
    let n_l = left_blk.n_cols();
    let n_r = right_blk.n_cols();
    let n_p = sub_a.seq.len();
    let est_cols = n_l.max(n_r).max(n_p);

    let mut l_seqs: Vec<Vec<u8>> = (0..left_blk.rows.len())
        .map(|_| Vec::with_capacity(est_cols))
        .collect();
    let mut r_seqs: Vec<Vec<u8>> = (0..right_blk.rows.len())
        .map(|_| Vec::with_capacity(est_cols))
        .collect();

    let (mut cl, mut cr, mut cp) = (0usize, 0usize, 0usize);
    while cl < n_l || cr < n_r || cp < n_p {
        if cl < n_l && left_blk.rows[a_idx].seq[cl] == b'-' {
            for (i, r) in left_blk.rows.iter().enumerate() {
                l_seqs[i].push(r.seq[cl]);
            }
            for rs in &mut r_seqs {
                rs.push(b'-');
            }
            cl += 1;
        } else if cr < n_r && right_blk.rows[b_idx].seq[cr] == b'-' {
            for ls in &mut l_seqs {
                ls.push(b'-');
            }
            for (i, r) in right_blk.rows.iter().enumerate() {
                r_seqs[i].push(r.seq[cr]);
            }
            cr += 1;
        } else if cp < n_p {
            let pa = sub_a.seq[cp];
            let pb = sub_b.seq[cp];
            cp += 1;
            match (pa != b'-', pb != b'-') {
                (true, true) => {
                    if cl >= n_l || cr >= n_r {
                        return None;
                    }
                    for (i, r) in left_blk.rows.iter().enumerate() {
                        l_seqs[i].push(r.seq[cl]);
                    }
                    for (i, r) in right_blk.rows.iter().enumerate() {
                        r_seqs[i].push(r.seq[cr]);
                    }
                    cl += 1;
                    cr += 1;
                }
                (true, false) => {
                    if cl >= n_l {
                        return None;
                    }
                    for (i, r) in left_blk.rows.iter().enumerate() {
                        l_seqs[i].push(r.seq[cl]);
                    }
                    for rs in &mut r_seqs {
                        rs.push(b'-');
                    }
                    cl += 1;
                }
                (false, true) => {
                    if cr >= n_r {
                        return None;
                    }
                    for ls in &mut l_seqs {
                        ls.push(b'-');
                    }
                    for (i, r) in right_blk.rows.iter().enumerate() {
                        r_seqs[i].push(r.seq[cr]);
                    }
                    cr += 1;
                }
                (false, false) => {}
            }
        } else {
            break;
        }
    }

    if cl != n_l || cr != n_r {
        return None;
    }

    let mut rows = Vec::with_capacity(left_blk.rows.len() + right_blk.rows.len());
    for (r, seq) in left_blk.rows.iter().zip(l_seqs) {
        rows.push(TbaRow {
            species: r.species.clone(),
            contig: r.contig.clone(),
            start: r.start,
            size: r.size,
            strand: r.strand,
            src_size: r.src_size,
            seq,
        });
    }
    for (r, seq) in right_blk.rows.iter().zip(r_seqs) {
        rows.push(TbaRow {
            species: r.species.clone(),
            contig: r.contig.clone(),
            start: r.start,
            size: r.size,
            strand: r.strand,
            src_size: r.src_size,
            seq,
        });
    }
    Some(TbaBlock { score, rows })
}

fn make_singleton_block(row: &TbaRow) -> TbaBlock {
    let ungapped: Vec<u8> = row.seq.iter().copied().filter(|&b| b != b'-').collect();
    TbaBlock {
        score: 0.0,
        rows: vec![TbaRow {
            species: row.species.clone(),
            contig: row.contig.clone(),
            start: row.start,
            size: ungapped.len(),
            strand: row.strand,
            src_size: row.src_size,
            seq: ungapped,
        }],
    }
}

/// Merge two clade blocksets (`left` over `left_species`, `right` over `right_species`) using
/// all cross-clade pairwise alignment blocks.
///
/// Preserves all unlinked clade-specific blocks intact while threading cross-clade links into
/// multi-species blocks with single-coverage netting by descending score.
#[must_use]
pub fn tba_merge_two_blocksets<S: std::hash::BuildHasher>(
    left: &[TbaBlock],
    right: &[TbaBlock],
    left_species: &HashSet<&str, S>,
    cross_pairs: &[&PairwiseTbaBlock],
) -> Vec<TbaBlock> {
    tba_merge_two_blocksets_filtered(left, right, left_species, cross_pairs, None)
}

#[allow(clippy::too_many_lines)]
fn tba_merge_two_blocksets_filtered<S: std::hash::BuildHasher>(
    left: &[TbaBlock],
    right: &[TbaBlock],
    left_species: &HashSet<&str, S>,
    cross_pairs: &[&PairwiseTbaBlock],
    left_only_species: Option<&str>,
) -> Vec<TbaBlock> {
    let left_idx = CladeIndex::build_for_species(left, left_only_species);
    let right_idx = CladeIndex::build(right);

    let mut left_avail: Vec<crate::netting::IntervalSet> = left
        .iter()
        .map(|b| crate::netting::IntervalSet::from_single(0, b.n_cols()))
        .collect();
    let mut right_avail: Vec<crate::netting::IntervalSet> = right
        .iter()
        .map(|b| crate::netting::IntervalSet::from_single(0, b.n_cols()))
        .collect();
    let mut left_singleton_used: HashMap<String, HashMap<String, crate::netting::IntervalSet>> =
        HashMap::new();
    let mut right_singleton_used: HashMap<String, HashMap<String, crate::netting::IntervalSet>> =
        HashMap::new();
    let mut iv_scratch: Vec<usize> = Vec::new();

    let mut oriented: Vec<Cow<'_, PairwiseTbaBlock>> = cross_pairs
        .iter()
        .filter_map(|p| p.oriented_for_left_right(left_species))
        .collect();
    oriented.sort_by(|a, b| b.score.total_cmp(&a.score));

    let mut out_blocks: Vec<TbaBlock> = Vec::new();

    for p_cow in &oriented {
        let p: &PairwiseTbaBlock = p_cow.as_ref();
        let n_cols = p.row_a.seq.len();
        let mut ka = 0usize;
        let mut kb = 0usize;
        // Per-column tag: None if already consumed, else Some((left_hit, right_hit)).
        let mut col_tags: Vec<ColTag> = Vec::with_capacity(n_cols);
        for col in 0..n_cols {
            let ab = p.row_a.seq[col];
            let bb = p.row_b.seq[col];
            let l_ok_and_hit = if ab == b'-' {
                Some(None)
            } else {
                let pa = p.row_a.fwd_pos(ka);
                ka += 1;
                if let Some((bi, ri)) = left_idx.lookup(&p.row_a.species, &p.row_a.contig, pa) {
                    if let Some(c_lo) = left_idx.fwd_pos_to_col(&left[bi], bi, ri, pa)
                        && left_avail[bi].contains_pos(c_lo)
                    {
                        Some(Some((bi, ri)))
                    } else {
                        None
                    }
                } else {
                    let used = left_singleton_used
                        .get(&p.row_a.species)
                        .and_then(|m| m.get(&p.row_a.contig))
                        .is_some_and(|u| u.contains_pos(pa));
                    if used { None } else { Some(None) }
                }
            };

            let r_ok_and_hit = if bb == b'-' {
                Some(None)
            } else {
                let pb = p.row_b.fwd_pos(kb);
                kb += 1;
                if let Some((bi, ri)) = right_idx.lookup(&p.row_b.species, &p.row_b.contig, pb) {
                    if let Some(c_lo) = right_idx.fwd_pos_to_col(&right[bi], bi, ri, pb)
                        && right_avail[bi].contains_pos(c_lo)
                    {
                        Some(Some((bi, ri)))
                    } else {
                        None
                    }
                } else {
                    let used = right_singleton_used
                        .get(&p.row_b.species)
                        .and_then(|m| m.get(&p.row_b.contig))
                        .is_some_and(|u| u.contains_pos(pb));
                    if used { None } else { Some(None) }
                }
            };

            match (ab != b'-', bb != b'-', l_ok_and_hit, r_ok_and_hit) {
                (true, true, Some(lh), Some(rh)) => col_tags.push(Some((lh, rh))),
                (true, false, Some(lh), _) => {
                    let prev_r = col_tags
                        .last()
                        .and_then(|t| t.as_ref())
                        .and_then(|&(_, r)| r);
                    col_tags.push(Some((lh, prev_r)));
                }
                (false, true, _, Some(rh)) => {
                    let prev_l = col_tags
                        .last()
                        .and_then(|t| t.as_ref())
                        .and_then(|&(l, _)| l);
                    col_tags.push(Some((prev_l, rh)));
                }
                _ => col_tags.push(None),
            }
        }

        // Partition `col_tags` into maximal runs with identical `Some((l_hit, r_hit))`.
        let mut c = 0usize;
        while c < n_cols {
            let Some((l_hit, r_hit)) = col_tags[c] else {
                c += 1;
                continue;
            };
            let start_c = c;
            while c < n_cols && col_tags[c] == Some((l_hit, r_hit)) {
                c += 1;
            }
            let Some((sub_a, sub_b)) = slice_pairwise_subspan(p, start_c, c) else {
                continue;
            };
            let (a_start, a_end) = sub_a.fwd_interval();
            let (b_start, b_end) = sub_b.fwd_interval();

            let (left_slice, a_idx_in_l, l_col_range) = if let Some((bi, ri)) = l_hit {
                let Some((col_lo, col_hi)) =
                    left_idx.fwd_interval_to_cols(&left[bi], bi, ri, a_start, a_end)
                else {
                    continue;
                };
                // Enforce strict single-coverage: reject if any internal gap column in
                // `[col_lo, col_hi)` was already consumed by a higher-scoring cross-clade link.
                if !left_avail[bi].contains_range(col_lo, col_hi) {
                    continue;
                }
                let Some(sliced) = left[bi].slice_cols(col_lo, col_hi) else {
                    continue;
                };
                let Some(new_ri) = sliced.rows.iter().position(|r| r.species == sub_a.species)
                else {
                    continue;
                };
                let oriented_l = if sliced.rows[new_ri].strand == sub_a.strand {
                    sliced
                } else {
                    sliced.revcomp()
                };
                (oriented_l, new_ri, Some((bi, col_lo, col_hi)))
            } else {
                (make_singleton_block(&sub_a), 0, None)
            };

            let (right_slice, b_idx_in_r, r_col_range) = if let Some((bi, ri)) = r_hit {
                let Some((col_lo, col_hi)) =
                    right_idx.fwd_interval_to_cols(&right[bi], bi, ri, b_start, b_end)
                else {
                    continue;
                };
                if !right_avail[bi].contains_range(col_lo, col_hi) {
                    continue;
                }
                let Some(sliced) = right[bi].slice_cols(col_lo, col_hi) else {
                    continue;
                };
                let Some(new_ri) = sliced.rows.iter().position(|r| r.species == sub_b.species)
                else {
                    continue;
                };
                let oriented_r = if sliced.rows[new_ri].strand == sub_b.strand {
                    sliced
                } else {
                    sliced.revcomp()
                };
                (oriented_r, new_ri, Some((bi, col_lo, col_hi)))
            } else {
                (make_singleton_block(&sub_b), 0, None)
            };

            #[allow(clippy::cast_precision_loss)]
            let pair_score = p.score * ((c - start_c) as f64 / (n_cols.max(1) as f64));
            let merged_score = left_slice.score + right_slice.score + pair_score;
            if let Some(merged) = thread_blocks_via_pair(
                &left_slice,
                a_idx_in_l,
                &right_slice,
                b_idx_in_r,
                &sub_a,
                &sub_b,
                merged_score,
            ) {
                if let Some((bi, lo, hi)) = l_col_range {
                    left_avail[bi].subtract(lo, hi, &mut iv_scratch);
                } else {
                    left_singleton_used
                        .entry(sub_a.species.clone())
                        .or_default()
                        .entry(sub_a.contig.clone())
                        .or_default()
                        .insert(a_start, a_end, &mut iv_scratch);
                }
                if let Some((bi, lo, hi)) = r_col_range {
                    right_avail[bi].subtract(lo, hi, &mut iv_scratch);
                } else {
                    right_singleton_used
                        .entry(sub_b.species.clone())
                        .or_default()
                        .entry(sub_b.contig.clone())
                        .or_default()
                        .insert(b_start, b_end, &mut iv_scratch);
                }
                out_blocks.push(merged);
            }
        }
    }

    // Preserve all unconsumed clade-specific block slices intact (reference-free!).
    for (bi, ranges) in left_avail.into_iter().enumerate() {
        for (lo, hi) in ranges.iter() {
            if let Some(slice) = left[bi].slice_cols(lo, hi)
                && slice.rows.len() >= 2
            {
                out_blocks.push(slice);
            }
        }
    }
    for (bi, ranges) in right_avail.into_iter().enumerate() {
        for (lo, hi) in ranges.iter() {
            if let Some(slice) = right[bi].slice_cols(lo, hi)
                && slice.rows.len() >= 2
            {
                out_blocks.push(slice);
            }
        }
    }

    out_blocks.sort_by(|a, b| {
        let ra = &a.rows[0];
        let rb = &b.rows[0];
        ra.species
            .cmp(&rb.species)
            .then_with(|| ra.contig.cmp(&rb.contig))
            .then_with(|| ra.fwd_interval().0.cmp(&rb.fwd_interval().0))
    });
    out_blocks
}

/// Run full bottom-up reference-free TBA progressive merging along `tree`, consuming pairwise
/// blocks keyed by unordered canonical species pairs `(min(a, b), max(a, b))`.
#[must_use]
pub fn run_tba_tree_merge<S: std::hash::BuildHasher + Sync>(
    tree: &GuideTree,
    pair_blocks: &HashMap<(String, String), Vec<PairwiseTbaBlock>, S>,
) -> Vec<TbaBlock> {
    match tree {
        GuideTree::Leaf(_) => Vec::new(),
        GuideTree::Internal { left, right } => {
            let (left_blocks, right_blocks) = rayon::join(
                || run_tba_tree_merge(left, pair_blocks),
                || run_tba_tree_merge(right, pair_blocks),
            );
            let left_leaves = left.leaves();
            let right_leaves = right.leaves();
            let left_set: HashSet<&str> = left_leaves.iter().copied().collect();

            let mut cross: Vec<&PairwiseTbaBlock> = Vec::new();
            for &ls in &left_leaves {
                for &rs in &right_leaves {
                    let key = canonical_pair_key(ls, rs);
                    if let Some(blks) = pair_blocks.get(&key) {
                        cross.extend(blks.iter());
                    }
                }
            }
            tba_merge_two_blocksets(&left_blocks, &right_blocks, &left_set, &cross)
        }
    }
}

/// Bottom-up reference-free TBA progressive merge along `tree` that loads cross-clade pairwise
/// blocks on demand via `load_pair(sp_a, sp_b)` at each pair's unique LCA internal node and drops
/// them immediately after that node merges. Uses sequential post-order DFS across subtrees with
/// Rayon parallelism inside each node's pair loading, bounding Phase 2 RAM to a single merge
/// frontier.
///
/// # Errors
/// Propagates any error returned by `load_pair`.
pub fn run_tba_tree_merge_lazy<F>(tree: &GuideTree, load_pair: &F) -> anyhow::Result<Vec<TbaBlock>>
where
    F: Fn(&str, &str) -> anyhow::Result<Vec<PairwiseTbaBlock>> + Sync,
{
    match tree {
        GuideTree::Leaf(_) => Ok(Vec::new()),
        GuideTree::Internal { left, right } => {
            let left_blocks = run_tba_tree_merge_lazy(left, load_pair)?;
            let right_blocks = run_tba_tree_merge_lazy(right, load_pair)?;
            let left_leaves = left.leaves();
            let right_leaves = right.leaves();
            let left_set: HashSet<&str> = left_leaves.iter().copied().collect();

            let mut pair_list: Vec<(&str, &str)> =
                Vec::with_capacity(left_leaves.len() * right_leaves.len());
            for &ls in &left_leaves {
                for &rs in &right_leaves {
                    pair_list.push((ls, rs));
                }
            }
            let owned_cross: Vec<Vec<PairwiseTbaBlock>> = pair_list
                .par_iter()
                .map(|&(ls, rs)| load_pair(ls, rs))
                .collect::<anyhow::Result<Vec<_>>>()?;
            let cross_refs: Vec<&PairwiseTbaBlock> =
                owned_cross.iter().flat_map(|v| v.iter()).collect();
            let merged =
                tba_merge_two_blocksets(&left_blocks, &right_blocks, &left_set, &cross_refs);
            drop(cross_refs);
            drop(owned_cross);
            drop(left_blocks);
            drop(right_blocks);
            Ok(merged)
        }
    }
}

/// Bottom-up reference-free TBA progressive merge along `tree` that removes each cross-clade
/// entry from `pair_blocks` as soon as its LCA internal node merges, freeing memory monotonically.
#[must_use]
pub fn run_tba_tree_merge_consuming<S: std::hash::BuildHasher>(
    tree: &GuideTree,
    pair_blocks: &mut HashMap<(String, String), Vec<PairwiseTbaBlock>, S>,
) -> Vec<TbaBlock> {
    match tree {
        GuideTree::Leaf(_) => Vec::new(),
        GuideTree::Internal { left, right } => {
            let left_blocks = run_tba_tree_merge_consuming(left, pair_blocks);
            let right_blocks = run_tba_tree_merge_consuming(right, pair_blocks);
            let left_leaves = left.leaves();
            let right_leaves = right.leaves();
            let left_set: HashSet<&str> = left_leaves.iter().copied().collect();

            let mut owned_cross: Vec<Vec<PairwiseTbaBlock>> = Vec::new();
            for &ls in &left_leaves {
                for &rs in &right_leaves {
                    let key = canonical_pair_key(ls, rs);
                    if let Some(blks) = pair_blocks.remove(&key) {
                        owned_cross.push(blks);
                    }
                }
            }
            let cross_refs: Vec<&PairwiseTbaBlock> =
                owned_cross.iter().flat_map(|v| v.iter()).collect();
            let merged =
                tba_merge_two_blocksets(&left_blocks, &right_blocks, &left_set, &cross_refs);
            drop(cross_refs);
            drop(owned_cross);
            drop(left_blocks);
            drop(right_blocks);
            merged
        }
    }
}

/// Project Tier-2 tip assemblies `(tip, assigned_hinge)` onto an existing reference-free TBA
/// `backbone` blockset using a streaming `load_tip_pair(tip, hinge)` provider so pairwise blocks
/// for only one tip exist in RAM at a time (`O(1)` working set in tip count `N`).
///
/// After threading all tips, shared conspecific tip insertions inside hinge-gap windows are
/// refined in linear `O(N_tips)` time via [`refine_shared_tip_insertions_poa`].
///
/// # Errors
/// Returns an error if any tip species already exists in `backbone`, appears more than once in
/// `tip_assignments`, maps to itself, or if `load_tip_pair` fails.
pub fn tba_project_tips_streaming<F>(
    backbone: Vec<TbaBlock>,
    tip_assignments: &[(String, String)],
    mut load_tip_pair: F,
) -> anyhow::Result<Vec<TbaBlock>>
where
    F: FnMut(&str, &str) -> anyhow::Result<Vec<PairwiseTbaBlock>>,
{
    let mut present: HashSet<String> = backbone
        .iter()
        .flat_map(|b| b.rows.iter().map(|r| r.species.clone()))
        .collect();
    let mut seen_tips: HashSet<&str> = HashSet::with_capacity(tip_assignments.len());
    for (tip, hinge) in tip_assignments {
        anyhow::ensure!(
            tip != hinge,
            "tip `{tip}` cannot map to itself as its own hinge"
        );
        anyhow::ensure!(
            seen_tips.insert(tip.as_str()),
            "duplicate tip `{tip}` in tip assignments"
        );
        anyhow::ensure!(
            !present.contains(tip.as_str()),
            "tip `{tip}` already exists in backbone blockset (cannot duplicate species)"
        );
    }
    for (_, hinge) in tip_assignments {
        anyhow::ensure!(
            !seen_tips.contains(hinge.as_str()),
            "chained tip `{hinge}` cannot serve as a hinge in tip assignments"
        );
        present.insert(hinge.clone());
    }
    let hinge_species = present.clone();

    let mut current = backbone;
    for (tip, hinge) in tip_assignments {
        let blks = load_tip_pair(tip, hinge)?;
        if !blks.is_empty() {
            let left_refs: HashSet<&str> = present.iter().map(String::as_str).collect();
            let cross: Vec<&PairwiseTbaBlock> = blks.iter().collect();
            current = tba_merge_two_blocksets_filtered(
                &current,
                &[],
                &left_refs,
                &cross,
                Some(hinge.as_str()),
            );
        }
        present.insert(tip.clone());
    }
    if tip_assignments.len() >= 2 {
        let _ = refine_shared_tip_insertions_poa(&mut current, &hinge_species);
    }
    Ok(current)
}

/// Project Tier-2 tip assemblies `(tip, assigned_hinge)` onto an existing reference-free TBA
/// `backbone` blockset using only each tip's pairwise alignment against its assigned hinge
/// (`O(N - H)` pairs; two-tier hinge + tip projection / post-hoc incremental extension).
///
/// # Errors
/// Returns an error if any tip species already exists in `backbone`, appears more than once in
/// `tip_assignments`, or maps to itself.
pub fn tba_project_tips_onto_blockset<S: std::hash::BuildHasher>(
    backbone: Vec<TbaBlock>,
    tip_assignments: &[(String, String)],
    pair_blocks: &HashMap<(String, String), Vec<PairwiseTbaBlock>, S>,
) -> anyhow::Result<Vec<TbaBlock>> {
    tba_project_tips_streaming(backbone, tip_assignments, |tip, hinge| {
        let key = canonical_pair_key(hinge, tip);
        Ok(pair_blocks.get(&key).cloned().unwrap_or_default())
    })
}
