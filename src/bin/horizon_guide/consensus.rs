//! Iterative ancestral consensus sequence construction and MAF remapping for `horizon_guide`.

use crate::rescue::maf_ref_coverage;
use crate::store::SpeciesStore;
use crate::{ContigList, ContigRef, contig_refs};
use anyhow::Result;
use horizon_aligner::coord::IntervalSet;
use horizon_aligner::dna::base_to_idx;
use horizon_aligner::{
    AlignConfig, ContigIndex, ProgressLogger, align_contig_pair_runs, format_maf_blocks_runs,
    parse_two_row_maf_blocks, prepare_contig_index_owned, select_top_k_pairs, shared_seed_count,
};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

trait VoteCount: Copy + Ord {
    const ZERO: Self;
    const ONE: Self;
    fn sat_inc(self) -> Self;
}

impl VoteCount for u8 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    #[inline]
    fn sat_inc(self) -> Self {
        self.saturating_add(1)
    }
}

impl VoteCount for u16 {
    const ZERO: Self = 0;
    const ONE: Self = 1;
    #[inline]
    fn sat_inc(self) -> Self {
        self.saturating_add(1)
    }
}

/// Bucket two-row MAF stanzas across batches in a single linear byte scan without allocating
/// `String` sequences for unmatched contigs or `-`-strand blocks.
fn bucket_maf_stanzas_by_batch(
    high_cov_mafs: &[&str],
    ctg_to_batch: &HashMap<&str, (usize, usize)>,
    n_batches: usize,
) -> Vec<Vec<(usize, std::ops::Range<usize>)>> {
    let mut by_batch: Vec<Vec<(usize, std::ops::Range<usize>)>> =
        (0..n_batches).map(|_| Vec::new()).collect();
    for (maf_idx, &maf_text) in high_cov_mafs.iter().enumerate() {
        let bytes = maf_text.as_bytes();
        let mut pos = 0usize;
        let mut stanza_start: Option<usize> = None;
        let mut stanza_batch: Option<usize> = None;
        let mut seen_first_s = false;
        while pos <= bytes.len() {
            let line_start = pos;
            let line_end = bytes[pos..]
                .iter()
                .position(|&b| b == b'\n')
                .map_or(bytes.len(), |rel| pos + rel);
            let line = maf_text[line_start..line_end].trim_end_matches('\r');
            let is_eof = pos == bytes.len();
            if is_eof || line.is_empty() || line.starts_with("a ") || line == "a" {
                if let (Some(st), Some(b_id)) = (stanza_start.take(), stanza_batch.take()) {
                    by_batch[b_id].push((maf_idx, st..line_start));
                }
                seen_first_s = false;
                if !is_eof && (line.starts_with("a ") || line == "a") {
                    stanza_start = Some(line_start);
                }
            } else if stanza_start.is_some() && !seen_first_s && line.starts_with("s ") {
                seen_first_s = true;
                let mut fields = line.split_ascii_whitespace();
                let _ = fields.next(); // "s"
                if let (Some(src), Some(_), Some(_), Some("+")) =
                    (fields.next(), fields.next(), fields.next(), fields.next())
                    && let Some(&(b_id, _)) = ctg_to_batch.get(src)
                {
                    stanza_batch = Some(b_id);
                }
            }
            if is_eof {
                break;
            }
            pos = (line_end + 1).min(bytes.len());
            if line_end == bytes.len() {
                if let (Some(st), Some(b_id)) = (stanza_start.take(), stanza_batch.take()) {
                    by_batch[b_id].push((maf_idx, st..bytes.len()));
                }
                break;
            }
        }
    }
    by_batch
}

fn vote_consensus_batch<C: VoteCount>(
    hinge_contigs: &[ContigRef<'_>],
    batch: &[usize],
    ctg_to_batch: &HashMap<&str, (usize, usize)>,
    high_cov_mafs: &[&str],
    stanzas: &[(usize, std::ops::Range<usize>)],
    out: &mut ContigList,
) {
    const BASES: [u8; 4] = *b"ACGT";
    let mut counts: Vec<Vec<[C; 4]>> = batch
        .iter()
        .map(|&ci| {
            hinge_contigs[ci]
                .1
                .iter()
                .map(|&b| {
                    let mut c = [C::ZERO; 4];
                    if let Some(bi) = base_to_idx(b) {
                        c[bi] = C::ONE;
                    }
                    c
                })
                .collect()
        })
        .collect();
    for (maf_idx, range) in stanzas {
        for (r0, r1) in parse_two_row_maf_blocks(&high_cov_mafs[*maf_idx][range.clone()]) {
            if r0.strand != '+' {
                continue;
            }
            let Some(&(_, local_i)) = ctg_to_batch.get(r0.src.as_str()) else {
                continue;
            };
            let mut rp = r0.start;
            for (rc, qc) in r0.seq.bytes().zip(r1.seq.bytes()) {
                if rc != b'-' {
                    if qc != b'-'
                        && rp < counts[local_i].len()
                        && let Some(bi) = base_to_idx(qc)
                    {
                        counts[local_i][rp][bi] = counts[local_i][rp][bi].sat_inc();
                    }
                    rp += 1;
                }
            }
        }
    }
    for (local_i, &ci) in batch.iter().enumerate() {
        let ref_seq = hinge_contigs[ci].1;
        let mut cons = Vec::with_capacity(ref_seq.len());
        for (pos, &rb) in ref_seq.iter().enumerate() {
            let c = counts[local_i][pos];
            let rb_idx = base_to_idx(rb);
            let mut best_b = rb;
            let mut best_cnt = rb_idx.map_or(C::ZERO, |ri| c[ri]);
            for (bi, &cnt) in c.iter().enumerate() {
                if cnt > best_cnt {
                    best_cnt = cnt;
                    best_b = if rb.is_ascii_lowercase() {
                        BASES[bi].to_ascii_lowercase()
                    } else {
                        BASES[bi]
                    };
                }
            }
            cons.push(best_b);
        }
        out[ci].1 = cons;
    }
}

pub(crate) fn build_ancestral_consensus_contigs(
    hinge_contigs: &[ContigRef<'_>],
    high_cov_mafs: &[&str],
    need_contig: &[bool],
) -> ContigList {
    // Cap per-batch vote buffer length at 16 Mb (~64 MiB for [u8; 4] or ~128 MiB for [u16; 4])
    // while scanning each MAF only once via precomputed per-batch stanza byte ranges.
    const VOTE_BATCH_BP: usize = 16 * 1024 * 1024;

    let mut out: ContigList = hinge_contigs
        .iter()
        .map(|&(name, _)| (name.to_owned(), Vec::new()))
        .collect();

    let needed_indices: Vec<usize> = (0..hinge_contigs.len())
        .filter(|&i| need_contig.get(i).copied().unwrap_or(true))
        .collect();
    if needed_indices.is_empty() {
        return out;
    }

    let mut batches: Vec<std::ops::Range<usize>> = Vec::new();
    let mut ctg_to_batch: HashMap<&str, (usize, usize)> =
        HashMap::with_capacity(needed_indices.len());
    let mut batch_start = 0usize;
    while batch_start < needed_indices.len() {
        let mut batch_end = batch_start;
        let mut batch_bp = 0usize;
        while batch_end < needed_indices.len() {
            let len = hinge_contigs[needed_indices[batch_end]].1.len();
            if batch_end > batch_start && batch_bp.saturating_add(len) > VOTE_BATCH_BP {
                break;
            }
            batch_bp = batch_bp.saturating_add(len);
            batch_end += 1;
        }
        let b_id = batches.len();
        for (local_i, &ci) in needed_indices[batch_start..batch_end].iter().enumerate() {
            ctg_to_batch.insert(hinge_contigs[ci].0, (b_id, local_i));
        }
        batches.push(batch_start..batch_end);
        batch_start = batch_end;
    }

    let stanzas_by_batch = bucket_maf_stanzas_by_batch(high_cov_mafs, &ctg_to_batch, batches.len());
    let use_u8 = high_cov_mafs.len() < usize::from(u8::MAX);
    for (b_id, range) in batches.into_iter().enumerate() {
        let batch = &needed_indices[range];
        if use_u8 {
            vote_consensus_batch::<u8>(
                hinge_contigs,
                batch,
                &ctg_to_batch,
                high_cov_mafs,
                &stanzas_by_batch[b_id],
                &mut out,
            );
        } else {
            vote_consensus_batch::<u16>(
                hinge_contigs,
                batch,
                &ctg_to_batch,
                high_cov_mafs,
                &stanzas_by_batch[b_id],
                &mut out,
            );
        }
    }

    out
}

pub(crate) fn remap_consensus_maf_to_hinge(
    hinge_contigs: &[ContigRef<'_>],
    cons_maf: &str,
    existing_maf: &str,
) -> String {
    let ctg_map: HashMap<&str, &[u8]> = hinge_contigs.iter().copied().collect();
    let mut covered_by_contig: HashMap<String, IntervalSet> = HashMap::new();
    let mut iv_scratch: Vec<usize> = Vec::new();
    for (r0, _) in parse_two_row_maf_blocks(existing_maf) {
        if r0.strand != '+' {
            continue;
        }
        let r0_size = r0.seq.bytes().filter(|&b| b != b'-').count();
        if r0_size > 0 {
            covered_by_contig.entry(r0.src).or_default().insert(
                r0.start,
                r0.start + r0_size,
                &mut iv_scratch,
            );
        }
    }
    let mut out = String::new();
    for (r0, r1) in parse_two_row_maf_blocks(cons_maf) {
        let Some(&ref_seq) = ctg_map.get(r0.src.as_str()) else {
            continue;
        };
        let r0_size = r0.seq.bytes().filter(|&b| b != b'-').count();
        if r0_size == 0 || r0.strand != '+' {
            continue;
        }
        let s0 = r0.start;
        let e0 = s0 + r0_size;
        if let Some(spans) = covered_by_contig.get(&r0.src)
            && spans.overlap_bp(s0, e0) * 10 >= r0_size * 9
        {
            continue;
        }
        let mut remapped_r0 = String::with_capacity(r0.seq.len());
        let mut rp = r0.start;
        for ch in r0.seq.bytes() {
            if ch == b'-' {
                remapped_r0.push('-');
            } else {
                let rb = ref_seq.get(rp).copied().unwrap_or(ch);
                remapped_r0.push(char::from(rb));
                rp += 1;
            }
        }
        let r1_size = r1.seq.bytes().filter(|&b| b != b'-').count();
        covered_by_contig
            .entry(r0.src.clone())
            .or_default()
            .insert(s0, e0, &mut iv_scratch);
        let _ = writeln!(out, "a score=1.0");
        let _ = writeln!(
            out,
            "s {} {} {} {} {} {}",
            r0.src, r0.start, r0_size, r0.strand, r0.src_size, remapped_r0
        );
        let _ = writeln!(
            out,
            "s {} {} {} {} {} {}",
            r1.src, r1.start, r1_size, r1.strand, r1.src_size, r1.seq
        );
        out.push('\n');
    }
    out
}

/// Identify hinge contigs whose aligned coverage in any low-coverage species is `< 80%`
/// so ancestral consensus construction only materializes needed contigs.
pub(crate) fn mark_uncovered_hinge_contigs(
    hinge_contigs: &[ContigRef<'_>],
    low_cov_species: &[usize],
    sp_to_pw_idx: &HashMap<usize, usize>,
    pairwise_mafs: &[(String, String)],
) -> Vec<bool> {
    let mut need_contig = vec![false; hinge_contigs.len()];
    let ctg_map: HashMap<&str, usize> = hinge_contigs
        .iter()
        .enumerate()
        .map(|(i, &(n, _))| (n, i))
        .collect();
    for &k in low_cov_species {
        let pw_idx = sp_to_pw_idx[&k];
        let mut intervals: Vec<Vec<(usize, usize)>> = vec![Vec::new(); hinge_contigs.len()];
        for (r0, _) in parse_two_row_maf_blocks(&pairwise_mafs[pw_idx].1) {
            if r0.strand == '+'
                && let Some(&ci) = ctg_map.get(r0.src.as_str())
            {
                let bp = r0.seq.bytes().filter(|&b| b != b'-').count();
                if bp > 0 {
                    intervals[ci].push((r0.start, r0.start.saturating_add(bp)));
                }
            }
        }
        for (ci, (_, seq)) in hinge_contigs.iter().enumerate() {
            let mut ivs = std::mem::take(&mut intervals[ci]);
            ivs.sort_unstable();
            let mut uniq_bp = 0usize;
            let mut cur_end = 0usize;
            for (s, e) in ivs {
                if e > cur_end {
                    uniq_bp += e - s.max(cur_end);
                    cur_end = e;
                }
            }
            if !seq.is_empty() && uniq_bp * 5 < seq.len() * 4 {
                need_contig[ci] = true;
            }
        }
    }
    need_contig
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn run_ancestral_consensus_rescue(
    store: &SpeciesStore,
    hinge_idx: usize,
    sp_to_pw_idx: &HashMap<usize, usize>,
    pairwise_mafs: &mut [(String, String)],
    cfg: &AlignConfig,
    top_k: usize,
    respect_softmask: bool,
    logger: ProgressLogger,
) -> Result<()> {
    const MIN_SHARED_SEEDS: usize = 8;
    let stage_t0 = Instant::now();
    let hinge_contigs = contig_refs(&store.hinge(hinge_idx).1);
    // Cap per-contig consensus fan-out at 4 (keeping >= 3 for --rbest / small --top-k) so
    // --top-k 20 on draft genomes does not re-run 16 unrelated chromosomes through 7-tier fallback.
    let rescue_top_k = if top_k == 0 { 0 } else { top_k.clamp(3, 4) };
    let mut total_rescued = 0usize;
    logger.log("ancestral rescue: starting (max 2 rounds)");
    for round in 0..2 {
        let round_t0 = Instant::now();
        let covs: Vec<(usize, f64)> = (0..store.len())
            .into_par_iter()
            .filter(|&i| i != hinge_idx)
            .map(|i| {
                let pw = sp_to_pw_idx[&i];
                (i, maf_ref_coverage(&hinge_contigs, &pairwise_mafs[pw].1))
            })
            .collect();
        let mut seed_mafs: Vec<&str> = covs
            .iter()
            .filter(|&&(_, c)| round > 0 || c >= 0.50)
            .map(|&(i, _)| pairwise_mafs[sp_to_pw_idx[&i]].1.as_str())
            .collect();
        if seed_mafs.is_empty() {
            seed_mafs = covs
                .iter()
                .filter(|&&(_, c)| c > 0.0)
                .map(|&(i, _)| pairwise_mafs[sp_to_pw_idx[&i]].1.as_str())
                .collect();
        }
        let low_cov_species: Vec<usize> = covs
            .iter()
            .filter(|&&(_, c)| c < 0.50)
            .map(|&(i, _)| i)
            .collect();
        if seed_mafs.is_empty() || low_cov_species.is_empty() {
            break;
        }
        let need_contig = mark_uncovered_hinge_contigs(
            &hinge_contigs,
            &low_cov_species,
            sp_to_pw_idx,
            pairwise_mafs,
        );
        let cons_contigs: ContigList =
            build_ancestral_consensus_contigs(&hinge_contigs, &seed_mafs, &need_contig)
                .into_iter()
                .enumerate()
                .filter(|&(ci, _)| need_contig[ci])
                .map(|(_, c)| c)
                .collect();
        if cons_contigs.is_empty() {
            break;
        }
        let cons_indices: Vec<(String, ContigIndex)> = cons_contigs
            .into_par_iter()
            .map(|(n, s)| Ok((n, prepare_contig_index_owned(s, respect_softmask)?)))
            .collect::<Result<Vec<_>>>()?;

        let mut cons_maf_by_sp: HashMap<usize, String> = HashMap::new();
        if let Some(indexed) = store.in_memory() {
            // Flatten all (species, cons_contig, qry_contig) pairs into one Rayon pool with live
            // progress logging so uneven query genomes never leave threads idle or go silent.
            let jobs: Vec<(usize, usize, usize)> = low_cov_species
                .par_iter()
                .flat_map(|&k| {
                    let qry_indices = &indexed[k].1;
                    let top_mask = (rescue_top_k > 0).then(|| {
                        select_top_k_pairs(
                            &cons_indices,
                            qry_indices,
                            rescue_top_k,
                            MIN_SHARED_SEEDS,
                        )
                    });
                    let mut sp_jobs = Vec::new();
                    for xi in 0..cons_indices.len() {
                        for yi in 0..qry_indices.len() {
                            let keep = if let Some(mask) = &top_mask {
                                mask[xi][yi]
                            } else {
                                shared_seed_count(&cons_indices[xi].1, &qry_indices[yi].1)
                                    >= MIN_SHARED_SEEDS
                            };
                            if keep {
                                sp_jobs.push((k, xi, yi));
                            }
                        }
                    }
                    sp_jobs
                })
                .collect();

            let total_jobs = jobs.len();
            logger.log(format_args!(
                "ancestral rescue round {}: scheduled {total_jobs} contig pairs across {} low-cov species (top_k={rescue_top_k})",
                round + 1,
                low_cov_species.len(),
            ));
            let report_step = (total_jobs / 10).max(1);
            let completed_jobs = AtomicUsize::new(0);
            let job_results: Vec<(usize, String)> = jobs
                .par_iter()
                .map(|&(k, xi, yi)| {
                    let (nx, ix) = &cons_indices[xi];
                    let (ny, iy) = &indexed[k].1[yi];
                    let aln = align_contig_pair_runs(ix, iy, cfg)?;
                    let done = completed_jobs.fetch_add(1, Ordering::Relaxed) + 1;
                    if done == total_jobs || done.is_multiple_of(report_step) {
                        let pct = (done * 100) / total_jobs.max(1);
                        logger.log(format_args!(
                            "ancestral rescue round {}: {done}/{total_jobs} contig pairs ({pct}%)",
                            round + 1,
                        ));
                    }
                    if aln.aligned_runs.is_empty() {
                        return Ok((k, String::new()));
                    }
                    let sx = ix.seq();
                    let sy = iy.seq();
                    Ok((
                        k,
                        format_maf_blocks_runs(nx, sx.len(), ny, sy.len(), sx, sy, &aln),
                    ))
                })
                .collect::<Result<Vec<_>>>()?;

            for (k, maf) in job_results {
                if !maf.is_empty() {
                    let body = maf
                        .strip_prefix("##maf version=1 scoring=horizon_mdl\n")
                        .unwrap_or(&maf);
                    cons_maf_by_sp
                        .entry(k)
                        .or_insert_with(|| String::from("##maf version=1 scoring=horizon_mdl\n"))
                        .push_str(body);
                }
            }
        } else {
            // Streaming mode: index & align one low-coverage species at a time across all threads.
            logger.log(format_args!(
                "ancestral rescue round {}: streaming {} low-cov species (top_k={rescue_top_k})",
                round + 1,
                low_cov_species.len(),
            ));
            for &k in &low_cov_species {
                let sp_k = store.load_streaming_species(k)?;
                let qry_indices = &sp_k.1;
                let top_mask = (rescue_top_k > 0).then(|| {
                    select_top_k_pairs(&cons_indices, qry_indices, rescue_top_k, MIN_SHARED_SEEDS)
                });
                let mut sp_jobs = Vec::new();
                for xi in 0..cons_indices.len() {
                    for yi in 0..qry_indices.len() {
                        let keep = if let Some(mask) = &top_mask {
                            mask[xi][yi]
                        } else {
                            shared_seed_count(&cons_indices[xi].1, &qry_indices[yi].1)
                                >= MIN_SHARED_SEEDS
                        };
                        if keep {
                            sp_jobs.push((xi, yi));
                        }
                    }
                }
                let sp_results: Vec<String> = sp_jobs
                    .par_iter()
                    .map(|&(xi, yi)| {
                        let (nx, ix) = &cons_indices[xi];
                        let (ny, iy) = &qry_indices[yi];
                        let aln = align_contig_pair_runs(ix, iy, cfg)?;
                        if aln.aligned_runs.is_empty() {
                            return Ok(String::new());
                        }
                        let sx = ix.seq();
                        let sy = iy.seq();
                        Ok(format_maf_blocks_runs(
                            nx,
                            sx.len(),
                            ny,
                            sy.len(),
                            sx,
                            sy,
                            &aln,
                        ))
                    })
                    .collect::<Result<Vec<_>>>()?;
                for maf in sp_results {
                    if !maf.is_empty() {
                        let body = maf
                            .strip_prefix("##maf version=1 scoring=horizon_mdl\n")
                            .unwrap_or(&maf);
                        cons_maf_by_sp
                            .entry(k)
                            .or_insert_with(|| {
                                String::from("##maf version=1 scoring=horizon_mdl\n")
                            })
                            .push_str(body);
                    }
                }
            }
        }

        let remapped_list: Vec<(usize, String)> = low_cov_species
            .par_iter()
            .filter_map(|&k| {
                let cons_maf = cons_maf_by_sp.get(&k)?;
                let pw_idx = sp_to_pw_idx[&k];
                let remapped = remap_consensus_maf_to_hinge(
                    &hinge_contigs,
                    cons_maf,
                    &pairwise_mafs[pw_idx].1,
                );
                Some((pw_idx, remapped))
            })
            .collect();

        let mut round_rescued = 0usize;
        for (pw_idx, remapped) in remapped_list {
            if !remapped.is_empty() {
                round_rescued += remapped.matches("\na score=").count()
                    + usize::from(remapped.starts_with("a score="));
                pairwise_mafs[pw_idx].1.push_str(&remapped);
            }
        }
        total_rescued += round_rescued;
        logger.log(format_args!(
            "ancestral rescue round {}: {} low-cov species, {} consensus contigs -> {round_rescued} blocks rescued ({:.1}s)",
            round + 1,
            low_cov_species.len(),
            cons_indices.len(),
            round_t0.elapsed().as_secs_f64(),
        ));
        if round_rescued == 0 {
            break;
        }
    }
    logger.log(format_args!(
        "ancestral rescue: done ({total_rescued} blocks rescued in {:.1}s)",
        stage_t0.elapsed().as_secs_f64()
    ));
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ancestral_consensus_u16_branch_handles_over_254_species() {
        let hinge = vec![("c1", b"AAAA".as_slice())];
        let maf_c = "a score=4\ns c1 0 4 + 4 AAAA\ns q 0 4 + 4 CCCC\n";
        let maf_g = "a score=4\ns c1 0 4 + 4 AAAA\ns q 0 4 + 4 GGGG\n";
        // 260 total MAFs (> u8::MAX): 140 vote G, 120 vote C -> G must win cleanly via u16 counts.
        let mut mafs = Vec::with_capacity(260);
        mafs.extend(std::iter::repeat_n(maf_c, 120));
        mafs.extend(std::iter::repeat_n(maf_g, 140));
        let cons = build_ancestral_consensus_contigs(&hinge, &mafs, &[true]);
        assert_eq!(cons[0].1, b"GGGG");
    }
}
