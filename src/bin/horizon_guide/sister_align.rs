//! Bounded, slice-aware sister-pair alignment for progressive transitive rescue in `horizon_guide`.
//!
//! When the hinge is a single chromosome or locus (e.g. `ce11 chrI 2 Mb` or `hg38 chr21`) and
//! query species are whole draft assemblies with 10–250 Mb chromosomes, aligning entire sister
//! chromosomes wastes >95% of runtime outside the hinge-mapped intervals that
//! `rescue_transitive_maf` can actually project back to the hinge. Clustering each contig's
//! hinge-mapped intervals within `50 kb` (so a single ectopic hit at a far telomere cannot stretch
//! a bounding box across a whole chromosome), slicing each cluster to `[lo, hi)`, and capping
//! per-slice fan-out to `sister_top_k = top_k.min(2)` eliminates that blowup losslessly on both
//! `+` and `-` strands.

use crate::IndexedSpecies;
use anyhow::Result;
use horizon_aligner::{
    AlignConfig, ContigIndex, align_contig_pair_runs, format_maf_blocks_runs,
    parse_two_row_maf_blocks, prepare_contig_index_owned, select_top_k_for_target_refs,
};
use rayon::prelude::*;
use std::collections::HashMap;

/// Hinge-mapped bounding cluster (0-based half-open `+`-strand coords) on one query contig.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct MappedContigSpan {
    pub(crate) contig_idx: usize,
    pub(crate) lo: usize,
    pub(crate) hi: usize,
    pub(crate) aligned_bp: usize,
}

/// Extract the top hinge-mapped contig interval clusters from `hinge_to_species_maf`.
///
/// Adjacent aligned blocks on the same contig within `50 kb` are merged into one slice cluster;
/// distant ectopic hits stay separate so a stray telomeric hit cannot stretch a slice across an
/// entire 100 Mb chromosome. Stray noise clusters (`< 100 bp`) are dropped when substantial
/// alignment exists.
pub(crate) fn hinge_mapped_contig_spans(
    species_indices: &[(String, ContigIndex)],
    hinge_to_species_maf: &str,
    max_slices: usize,
) -> Vec<MappedContigSpan> {
    const MERGE_GAP_BP: usize = 50_000;
    const PAD_BP: usize = 500;
    const MIN_FRAGMENT_BP: usize = 100;
    const NOISE_FILTER_TOTAL_BP: usize = 500;

    if hinge_to_species_maf.is_empty() || species_indices.is_empty() {
        return Vec::new();
    }
    let name_to_idx: HashMap<&str, usize> = species_indices
        .iter()
        .enumerate()
        .map(|(idx, (name, _))| (name.as_str(), idx))
        .collect();
    let mut intervals_by_ctg: HashMap<usize, Vec<(usize, usize, usize)>> = HashMap::new();
    let mut total_aligned = 0usize;
    for (_, q_row) in parse_two_row_maf_blocks(hinge_to_species_maf) {
        let Some(&ci) = name_to_idx.get(q_row.src.as_str()) else {
            continue;
        };
        let non_gap = q_row.seq.bytes().filter(|&b| b != b'-').count();
        if non_gap == 0 {
            continue;
        }
        let a = q_row.fwd_pos(0);
        let b = q_row.fwd_pos(non_gap - 1);
        let lo = a.min(b);
        let hi = a.max(b).saturating_add(1);
        total_aligned = total_aligned.saturating_add(non_gap);
        intervals_by_ctg
            .entry(ci)
            .or_default()
            .push((lo, hi, non_gap));
    }

    let min_bp = if total_aligned >= NOISE_FILTER_TOTAL_BP {
        MIN_FRAGMENT_BP
    } else {
        1
    };

    let mut spans = Vec::new();
    for (ci, mut ivs) in intervals_by_ctg {
        let ctg_len = species_indices[ci].1.seq().len();
        ivs.sort_unstable_by_key(|&(lo, hi, _)| (lo, hi));
        let mut cur_lo = ivs[0].0;
        let mut cur_hi = ivs[0].1;
        let mut cur_bp = ivs[0].2;
        for &(lo, hi, bp) in &ivs[1..] {
            if lo <= cur_hi.saturating_add(MERGE_GAP_BP) {
                cur_hi = cur_hi.max(hi);
                cur_bp = cur_bp.saturating_add(bp);
            } else {
                if cur_bp >= min_bp {
                    let plo = cur_lo.saturating_sub(PAD_BP);
                    let phi = cur_hi.saturating_add(PAD_BP).min(ctg_len);
                    if plo < phi {
                        spans.push(MappedContigSpan {
                            contig_idx: ci,
                            lo: plo,
                            hi: phi,
                            aligned_bp: cur_bp,
                        });
                    }
                }
                cur_lo = lo;
                cur_hi = hi;
                cur_bp = bp;
            }
        }
        if cur_bp >= min_bp {
            let plo = cur_lo.saturating_sub(PAD_BP);
            let phi = cur_hi.saturating_add(PAD_BP).min(ctg_len);
            if plo < phi {
                spans.push(MappedContigSpan {
                    contig_idx: ci,
                    lo: plo,
                    hi: phi,
                    aligned_bp: cur_bp,
                });
            }
        }
    }

    spans.sort_unstable_by(|a, b| {
        b.aligned_bp
            .cmp(&a.aligned_bp)
            .then_with(|| a.contig_idx.cmp(&b.contig_idx))
            .then_with(|| a.lo.cmp(&b.lo))
    });
    if max_slices > 0 && spans.len() > max_slices {
        spans.truncate(max_slices);
    }
    spans
}

enum SliceIndex<'a> {
    Borrowed(&'a ContigIndex),
    Owned(ContigIndex),
}

impl SliceIndex<'_> {
    #[inline]
    fn as_ref(&self) -> &ContigIndex {
        match self {
            Self::Borrowed(ix) => ix,
            Self::Owned(ix) => ix,
        }
    }
}

struct PreparedSlice<'a> {
    contig_idx: usize,
    offset: usize,
    index: SliceIndex<'a>,
}

fn prepare_mapped_slices<'a>(
    species_indices: &'a [(String, ContigIndex)],
    spans: &[MappedContigSpan],
    skip_softmasked: bool,
) -> Result<Vec<PreparedSlice<'a>>> {
    // Only slice & re-index when the hinge-mapped cluster is < 50% of a >= 4 kb contig;
    // otherwise borrow the existing full contig index directly (zero re-indexing and zero clones).
    const MIN_SLICE_CONTIG_LEN: usize = 4_000;
    spans
        .par_iter()
        .map(|span| {
            let (_, full_idx) = &species_indices[span.contig_idx];
            let full_seq = full_idx.seq();
            let span_len = span.hi - span.lo;
            if full_seq.len() >= MIN_SLICE_CONTIG_LEN && span_len.saturating_mul(2) < full_seq.len()
            {
                let sub_idx = prepare_contig_index_owned(
                    full_seq[span.lo..span.hi].to_vec(),
                    skip_softmasked,
                )?;
                Ok(PreparedSlice {
                    contig_idx: span.contig_idx,
                    offset: span.lo,
                    index: SliceIndex::Owned(sub_idx),
                })
            } else {
                Ok(PreparedSlice {
                    contig_idx: span.contig_idx,
                    offset: 0,
                    index: SliceIndex::Borrowed(full_idx),
                })
            }
        })
        .collect()
}

fn pick_slice_targets(
    slices: &[PreparedSlice<'_>],
    other_indices: &[(String, ContigIndex)],
    sister_top_k: usize,
    min_shared_seeds: usize,
) -> Vec<(usize, usize)> {
    if slices.is_empty() || other_indices.is_empty() {
        return Vec::new();
    }
    let slice_refs: Vec<&ContigIndex> = slices.iter().map(|s| s.index.as_ref()).collect();
    select_top_k_for_target_refs(&slice_refs, other_indices, sister_top_k, min_shared_seeds)
}

/// Align only the hinge-mapped bounding slices of `ref_sp` and `qry_sp`, capping per-slice
/// fan-out to `sister_top_k = top_k.min(2)` and remapping all run coordinates back onto the
/// full contig sequences before MAF formatting.
pub(crate) fn align_species_pair_sliced(
    ref_sp: &IndexedSpecies,
    ref_spans: &[MappedContigSpan],
    qry_sp: &IndexedSpecies,
    qry_spans: &[MappedContigSpan],
    cfg: &AlignConfig,
    top_k: usize,
) -> Result<(String, usize)> {
    const MIN_SHARED_SEEDS: usize = 8;
    if ref_spans.is_empty() && qry_spans.is_empty() {
        return Ok((String::new(), 0));
    }

    let ref_indices = &ref_sp.1;
    let qry_indices = &qry_sp.1;
    let sister_top_k = if top_k == 0 { 0 } else { top_k.min(2) };

    let (ref_slices, qry_slices) = rayon::join(
        || prepare_mapped_slices(ref_indices, ref_spans, cfg.skip_softmasked),
        || prepare_mapped_slices(qry_indices, qry_spans, cfg.skip_softmasked),
    );
    let ref_slices = ref_slices?;
    let qry_slices = qry_slices?;

    let mut qry_slices_by_ci: HashMap<usize, Vec<usize>> = HashMap::new();
    for (sj, sl) in qry_slices.iter().enumerate() {
        qry_slices_by_ci.entry(sl.contig_idx).or_default().push(sj);
    }
    let mut ref_slices_by_ci: HashMap<usize, Vec<usize>> = HashMap::new();
    for (si, sl) in ref_slices.iter().enumerate() {
        ref_slices_by_ci.entry(sl.contig_idx).or_default().push(si);
    }

    // Deduplicated jobs: (ref_contig_idx, Option<ref_slice_idx>, qry_contig_idx, Option<qry_slice_idx>)
    let mut jobs: Vec<(usize, Option<usize>, usize, Option<usize>)> = Vec::new();
    for (si, yi) in pick_slice_targets(&ref_slices, qry_indices, sister_top_k, MIN_SHARED_SEEDS) {
        let xi = ref_slices[si].contig_idx;
        if let Some(q_list) = qry_slices_by_ci.get(&yi) {
            for &sj in q_list {
                jobs.push((xi, Some(si), yi, Some(sj)));
            }
        } else {
            jobs.push((xi, Some(si), yi, None));
        }
    }
    for (sj, xi) in pick_slice_targets(&qry_slices, ref_indices, sister_top_k, MIN_SHARED_SEEDS) {
        let yi = qry_slices[sj].contig_idx;
        if let Some(r_list) = ref_slices_by_ci.get(&xi) {
            for &si in r_list {
                jobs.push((xi, Some(si), yi, Some(sj)));
            }
        } else {
            jobs.push((xi, None, yi, Some(sj)));
        }
    }

    jobs.sort_unstable();
    jobs.dedup();
    if jobs.is_empty() {
        return Ok((String::new(), 0));
    }

    let n_jobs = jobs.len();
    let results: Vec<String> = jobs
        .par_iter()
        .map(|&(xi, opt_si, yi, opt_sj)| {
            let (ix, x_off) = match opt_si {
                Some(si) => (ref_slices[si].index.as_ref(), ref_slices[si].offset),
                None => (&ref_indices[xi].1, 0),
            };
            let (iy, y_off) = match opt_sj {
                Some(sj) => (qry_slices[sj].index.as_ref(), qry_slices[sj].offset),
                None => (&qry_indices[yi].1, 0),
            };
            let mut aln = align_contig_pair_runs(ix, iy, cfg)?;
            if aln.aligned_runs.is_empty() {
                return Ok(String::new());
            }
            if x_off > 0 || y_off > 0 {
                for run in &mut aln.aligned_runs {
                    run.x_start += x_off;
                    run.y_start += y_off;
                }
            }
            let (ref_name, ref_full) = &ref_indices[xi];
            let (qry_name, qry_full) = &qry_indices[yi];
            let sx = ref_full.seq();
            let sy = qry_full.seq();
            Ok(format_maf_blocks_runs(
                ref_name,
                sx.len(),
                qry_name,
                sy.len(),
                sx,
                sy,
                &aln,
            ))
        })
        .collect::<Result<Vec<_>>>()?;

    let mut maf_buf = String::from("##maf version=1 scoring=horizon_mdl\n");
    for maf in results {
        if !maf.is_empty() {
            let body = maf
                .strip_prefix("##maf version=1 scoring=horizon_mdl\n")
                .unwrap_or(&maf);
            maf_buf.push_str(body);
        }
    }
    Ok((maf_buf, n_jobs))
}

#[cfg(test)]
mod tests {
    use super::*;
    use horizon_aligner::complement_base;

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
    fn hinge_mapped_contig_spans_clusters_distant_hits_and_drops_stray_noise() {
        let mk = |name: &str, len: usize| {
            (
                name.to_string(),
                prepare_contig_index_owned(prng_dna(len, 17), false)
                    .expect("valid synthetic contig"),
            )
        };
        let sp = vec![mk("chr1", 200_000), mk("stray", 10_000)];
        let block600 = "A".repeat(600);
        let block150 = "G".repeat(150);
        let noise20 = "C".repeat(20);
        // chr1 has a 600 bp '-' strand hit at fwd [188_400..189_000) and a distant 150 bp hit at
        // fwd [1_000..1_150) (>50 kb apart -> must become TWO small clusters, NOT one 188 kb box!).
        let maf = format!(
            "a score=600\ns hinge 0 600 + 2000 {block600}\ns chr1 11000 600 - 200000 {block600}\n\n\
             a score=150\ns hinge 0 150 + 2000 {block150}\ns chr1 1000 150 + 200000 {block150}\n\n\
             a score=20\ns hinge 0 20 + 2000 {noise20}\ns stray 500 20 + 10000 {noise20}\n\n"
        );
        let spans = hinge_mapped_contig_spans(&sp, &maf, 10);
        assert_eq!(spans.len(), 2);
        assert_eq!(spans[0].contig_idx, 0);
        assert_eq!(spans[0].lo, 188_400 - 500);
        assert_eq!(spans[0].hi, 189_000 + 500);
        assert_eq!(spans[0].aligned_bp, 600);
        assert_eq!(spans[1].contig_idx, 0);
        assert_eq!(spans[1].lo, 1_000 - 500);
        assert_eq!(spans[1].hi, 1_150 + 500);
        assert_eq!(spans[1].aligned_bp, 150);
    }

    #[test]
    fn sliced_sister_alignment_remaps_rc_and_plus_coords_to_full_contigs() {
        let core = prng_dna(400, 99);
        let core_rc: Vec<u8> = core.iter().rev().map(|&b| complement_base(b)).collect();
        let mut full_x = prng_dna(6_000, 1);
        full_x[3_200..3_600].copy_from_slice(&core);
        let mut full_y = prng_dna(6_000, 2);
        full_y[1_500..1_900].copy_from_slice(&core_rc);

        let ref_sp = (
            "spA".to_string(),
            vec![(
                "cA".to_string(),
                prepare_contig_index_owned(full_x, false).expect("idx x"),
            )],
        );
        let qry_sp = (
            "spB".to_string(),
            vec![(
                "cB".to_string(),
                prepare_contig_index_owned(full_y, false).expect("idx y"),
            )],
        );
        let ref_spans = vec![MappedContigSpan {
            contig_idx: 0,
            lo: 3_000,
            hi: 3_800,
            aligned_bp: 400,
        }];
        let qry_spans = vec![MappedContigSpan {
            contig_idx: 0,
            lo: 1_300,
            hi: 2_100,
            aligned_bp: 400,
        }];
        let (maf, n_pairs) = align_species_pair_sliced(
            &ref_sp,
            &ref_spans,
            &qry_sp,
            &qry_spans,
            &AlignConfig::default(),
            2,
        )
        .expect("sliced align");
        assert_eq!(n_pairs, 1);
        let blocks = parse_two_row_maf_blocks(&maf);
        assert!(!blocks.is_empty(), "expected aligned RC block in {maf}");
        let (rx, ry) = &blocks[0];
        assert_eq!(rx.src_size, 6_000);
        assert_eq!(ry.src_size, 6_000);
        assert_eq!(ry.strand, '-');
        assert!(
            rx.start >= 3_180 && rx.start <= 3_220,
            "rx.start={}",
            rx.start
        );
        assert!(
            ry.start >= 4_080 && ry.start <= 4_120,
            "ry.start={}",
            ry.start
        );
    }
}
