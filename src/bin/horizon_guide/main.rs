//! Nystrom/sketch-selected hinge + progressive multi-species WGA/HAL pipeline (`horizon_guide`).

use anyhow::{Context, Result, bail};
use consensus::run_ancestral_consensus_rescue;
use hal_export::{HalBlock, HalMultiBlock, export_guide_hal, runs_to_hal_blocks};
use horizon_aligner::{
    AlignConfig, ContigIndex, PairSchedule, ProgressLogger, adaptive_top_k_for_fragmentation,
    align_contig_pair_runs, format_maf_blocks_runs, init_rayon_threads, merge_contig_sketches,
    parse_species_fasta_spec, parse_threads_arg, select_top_k_pairs, shared_seed_count,
    sketch_jaccard, validate_species_name,
};
use rayon::prelude::*;
use rescue::run_progressive_sister_rescue;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;
use store::{HINGE_SKETCH_K, SpeciesStore};

mod consensus;
mod hal_export;
mod rescue;
mod sister_align;
mod store;

pub(crate) type ContigList = Vec<(String, Vec<u8>)>;
pub(crate) type ContigRef<'a> = (&'a str, &'a [u8]);
pub(crate) type IndexedSpecies = (String, Vec<(String, ContigIndex)>);
pub(crate) type SpeciesContigs = (String, ContigList);

pub(crate) fn contig_refs(indices: &[(String, ContigIndex)]) -> Vec<ContigRef<'_>> {
    indices
        .iter()
        .map(|(n, ix)| (n.as_str(), ix.seq()))
        .collect()
}

pub(crate) fn select_central_hinge(sketches: &[Vec<u64>], effective_contigs: &[f64]) -> usize {
    (0..sketches.len())
        .into_par_iter()
        .map(|i| {
            let total: f64 = (0..sketches.len())
                .filter(|&j| j != i)
                .map(|j| sketch_jaccard(&sketches[i], &sketches[j]))
                .sum();
            // First-principles length-weighted contiguity factor `1 + ln(total_bp / auN)`:
            // identically `1.0` for single-contig genomes, resilient to tiny unplaced scaffolds,
            // and scales smoothly with effective fragmentation without arbitrary tuning weights.
            let n_eff = effective_contigs.get(i).copied().unwrap_or(1.0).max(1.0);
            let weighted = total / (1.0 + n_eff.ln());
            (weighted, i)
        })
        .reduce(
            || (-1.0_f64, usize::MAX),
            |best, cur| match cur.0.total_cmp(&best.0) {
                std::cmp::Ordering::Greater => cur,
                std::cmp::Ordering::Less => best,
                std::cmp::Ordering::Equal => {
                    if cur.1 < best.1 {
                        cur
                    } else {
                        best
                    }
                }
            },
        )
        .1
        .min(sketches.len().saturating_sub(1))
}

const GUIDE_USAGE: &str = "usage: horizon_guide [--progressive] [--respect-softmask] [--hal] [--hinge <species>] [--threads N | -t N] [--quiet | -q] [--rbest | --top-k K | --all-pairs] <out_prefix> <max_contigs> <sp1=fa1> <sp2=fa2>...";

#[allow(clippy::struct_excessive_bools)]
pub(crate) struct GuideCliArgs {
    pub(crate) progressive: bool,
    pub(crate) respect_softmask: bool,
    write_hal: bool,
    quiet: bool,
    pub(crate) hinge: Option<String>,
    threads: Option<usize>,
    top_k: usize,
    pub(crate) out_prefix: String,
    pub(crate) max_contigs: usize,
    pub(crate) specs: Vec<(String, PathBuf)>,
}

#[allow(clippy::too_many_lines)]
fn parse_args() -> Result<Option<GuideCliArgs>> {
    let raw_args: Vec<String> = std::env::args().skip(1).collect();
    let mut progressive = false;
    let mut respect_softmask = false;
    let mut write_hal = false;
    let mut quiet = false;
    let mut hinge: Option<String> = None;
    let mut threads: Option<usize> = None;
    let mut schedule_opt: Option<PairSchedule> = None;
    let mut pos = Vec::new();
    let mut arg_i = 0usize;
    while arg_i < raw_args.len() {
        match raw_args[arg_i].as_str() {
            "-h" | "--help" => {
                println!("{GUIDE_USAGE}");
                return Ok(None);
            }
            "-V" | "--version" => {
                println!("horizon_guide {}", env!("CARGO_PKG_VERSION"));
                return Ok(None);
            }
            "-t" | "--threads" => {
                let raw = raw_args
                    .get(arg_i + 1)
                    .context("--threads requires an integer argument")?;
                threads = Some(parse_threads_arg(raw)?);
                arg_i += 2;
            }
            "--hinge" => {
                let raw = raw_args
                    .get(arg_i + 1)
                    .context("--hinge requires a species name")?;
                validate_species_name(raw)?;
                hinge = Some(raw.clone());
                arg_i += 2;
            }
            "--progressive" => {
                progressive = true;
                arg_i += 1;
            }
            "--respect-softmask" => {
                respect_softmask = true;
                arg_i += 1;
            }
            "--hal" => {
                write_hal = true;
                arg_i += 1;
            }
            "-q" | "--quiet" => {
                quiet = true;
                arg_i += 1;
            }
            "--rbest" => {
                schedule_opt = Some(PairSchedule::TopK(1));
                arg_i += 1;
            }
            "--all-pairs" => {
                schedule_opt = Some(PairSchedule::AllPairs);
                arg_i += 1;
            }
            "--top-k" => {
                let Some(val) = raw_args.get(arg_i + 1) else {
                    bail!("--top-k requires an integer argument");
                };
                let parsed = val
                    .parse::<usize>()
                    .with_context(|| format!("invalid --top-k value `{val}`"))?;
                schedule_opt = Some(if parsed == 0 {
                    PairSchedule::AllPairs
                } else {
                    PairSchedule::TopK(parsed)
                });
                arg_i += 2;
            }
            other if other.starts_with('-') => bail!("unknown flag `{other}`"),
            other => {
                pos.push(other.to_owned());
                arg_i += 1;
            }
        }
    }
    if pos.len() < 4 {
        bail!("{GUIDE_USAGE}");
    }
    let out_prefix = pos[0].clone();
    let max_contigs: usize = pos[1].parse().context("parsing max_contigs")?;
    let top_k = PairSchedule::resolve_top_k(schedule_opt, max_contigs);
    let mut specs = Vec::new();
    let mut seen_species = std::collections::HashSet::new();
    for spec in &pos[2..] {
        let (sp, path) = parse_species_fasta_spec(spec)?;
        anyhow::ensure!(
            seen_species.insert(sp.clone()),
            "duplicate species name: {sp}"
        );
        specs.push((sp, path));
    }
    if let Some(ref hinge_name) = hinge
        && !seen_species.contains(hinge_name)
    {
        let avail: Vec<&str> = specs.iter().map(|(sp, _)| sp.as_str()).collect();
        bail!(
            "--hinge species `{hinge_name}` not found in input specs (available: {})",
            avail.join(", ")
        );
    }
    Ok(Some(GuideCliArgs {
        progressive,
        respect_softmask,
        write_hal,
        quiet,
        hinge,
        threads,
        top_k,
        out_prefix,
        max_contigs,
        specs,
    }))
}

pub(crate) struct GuideAlignOutput {
    pub(crate) pairwise_mafs: Vec<(String, String)>,
    pub(crate) hal_query_indices: Vec<usize>,
    pub(crate) hal_all_blocks: Vec<HalMultiBlock>,
    pub(crate) direct_runs: Vec<(usize, usize, usize, horizon_aligner::RunAlignment)>,
}

fn schedule_hinge_jobs(
    indexed: &[IndexedSpecies],
    hinge_idx: usize,
    top_k: usize,
) -> (Vec<usize>, Vec<(usize, usize, usize)>) {
    const MIN_SHARED_SEEDS: usize = 8;
    let hinge_indices = &indexed[hinge_idx].1;
    let query_species_indices: Vec<usize> =
        (0..indexed.len()).filter(|&i| i != hinge_idx).collect();
    let all_jobs: Vec<(usize, usize, usize)> = query_species_indices
        .par_iter()
        .flat_map(|&sp_i| {
            let qry_indices = &indexed[sp_i].1;
            let eff_top_k = adaptive_top_k_for_fragmentation(top_k, hinge_indices, qry_indices);
            let top_mask = (eff_top_k > 0).then(|| {
                select_top_k_pairs(hinge_indices, qry_indices, eff_top_k, MIN_SHARED_SEEDS)
            });
            let mut sp_jobs = Vec::new();
            for xi in 0..hinge_indices.len() {
                for yi in 0..qry_indices.len() {
                    let keep = if let Some(mask) = &top_mask {
                        mask[xi][yi]
                    } else {
                        shared_seed_count(&hinge_indices[xi].1, &qry_indices[yi].1)
                            >= MIN_SHARED_SEEDS
                    };
                    if keep {
                        sp_jobs.push((sp_i, xi, yi));
                    }
                }
            }
            sp_jobs
        })
        .collect();
    (query_species_indices, all_jobs)
}

#[allow(clippy::too_many_arguments)]
fn run_progressive_alignments(
    store: &SpeciesStore,
    species_sketches: &[Vec<u64>],
    hinge_idx: usize,
    query_species_indices: &[usize],
    all_jobs: &[(usize, usize, usize)],
    cli: &GuideCliArgs,
    cfg: &AlignConfig,
    logger: ProgressLogger,
) -> Result<GuideAlignOutput> {
    let indexed = store
        .in_memory()
        .expect("run_progressive_alignments requires InMemory store");
    let (_, hinge_indices) = &indexed[hinge_idx];
    let total_jobs = all_jobs.len();
    let report_step = (total_jobs / 10).max(1);
    let completed_jobs = AtomicUsize::new(0);

    let job_results: Vec<(usize, String)> = all_jobs
        .par_iter()
        .map(|&(sp_i, xi, yi)| {
            let (_, qry_indices) = &indexed[sp_i];
            let (nx, ix) = &hinge_indices[xi];
            let (ny, iy) = &qry_indices[yi];
            let aln = align_contig_pair_runs(ix, iy, cfg)?;
            let done = completed_jobs.fetch_add(1, Ordering::Relaxed) + 1;
            if done == total_jobs || done.is_multiple_of(report_step) {
                let pct = (done * 100) / total_jobs.max(1);
                logger.log(format_args!(
                    "pairwise: {done}/{total_jobs} contig pairs ({pct}%)"
                ));
            }
            if aln.aligned_runs.is_empty() {
                return Ok((sp_i, String::new()));
            }
            let sx = ix.seq();
            let sy = iy.seq();
            let maf = format_maf_blocks_runs(nx, sx.len(), ny, sy.len(), sx, sy, &aln);
            Ok((sp_i, maf))
        })
        .collect::<Result<Vec<_>>>()?;

    let mut maf_by_sp: Vec<String> =
        vec![String::from("##maf version=1 scoring=horizon_mdl\n"); indexed.len()];
    for (sp_i, maf) in job_results {
        if !maf.is_empty() {
            let body = maf
                .strip_prefix("##maf version=1 scoring=horizon_mdl\n")
                .unwrap_or(&maf);
            maf_by_sp[sp_i].push_str(body);
        }
    }

    let mut pairwise_mafs = Vec::with_capacity(query_species_indices.len());
    let mut hal_query_indices = Vec::new();
    let mut sp_to_pw_idx = HashMap::with_capacity(query_species_indices.len());
    for &i in query_species_indices {
        let qry_sp = &indexed[i].0;
        if cli.write_hal {
            hal_query_indices.push(i);
        }
        sp_to_pw_idx.insert(i, pairwise_mafs.len());
        pairwise_mafs.push((qry_sp.clone(), std::mem::take(&mut maf_by_sp[i])));
    }

    run_progressive_sister_rescue(
        store,
        species_sketches,
        hinge_idx,
        &sp_to_pw_idx,
        &mut pairwise_mafs,
        cfg,
        cli.top_k,
        logger,
    )?;
    run_ancestral_consensus_rescue(
        store,
        hinge_idx,
        &sp_to_pw_idx,
        &mut pairwise_mafs,
        cfg,
        cli.top_k,
        cli.respect_softmask,
        logger,
    )?;

    Ok(GuideAlignOutput {
        pairwise_mafs,
        hal_query_indices,
        hal_all_blocks: Vec::new(),
        direct_runs: Vec::new(),
    })
}

fn run_star_alignments(
    indexed: &[IndexedSpecies],
    hinge_idx: usize,
    query_species_indices: &[usize],
    all_jobs: &[(usize, usize, usize)],
    write_hal: bool,
    cfg: &AlignConfig,
    logger: ProgressLogger,
) -> Result<GuideAlignOutput> {
    let hinge_indices = &indexed[hinge_idx].1;
    let total_jobs = all_jobs.len();
    let report_step = (total_jobs / 10).max(1);
    let completed_jobs = AtomicUsize::new(0);

    let job_results: Vec<(usize, usize, usize, horizon_aligner::RunAlignment)> = all_jobs
        .par_iter()
        .map(|&(sp_i, xi, yi)| {
            let qry_indices = &indexed[sp_i].1;
            let aln = align_contig_pair_runs(&hinge_indices[xi].1, &qry_indices[yi].1, cfg)?;
            let done = completed_jobs.fetch_add(1, Ordering::Relaxed) + 1;
            if done == total_jobs || done.is_multiple_of(report_step) {
                let pct = (done * 100) / total_jobs.max(1);
                logger.log(format_args!(
                    "pairwise: {done}/{total_jobs} contig pairs ({pct}%)"
                ));
            }
            Ok((sp_i, xi, yi, aln))
        })
        .collect::<Result<Vec<_>>>()?;

    let mut direct_runs = Vec::new();
    let mut hblks_by_sp: Vec<Vec<HalBlock>> = vec![Vec::new(); indexed.len()];
    for (sp_i, xi, yi, aln) in job_results {
        if aln.aligned_runs.is_empty() {
            continue;
        }
        if write_hal {
            hblks_by_sp[sp_i].extend(runs_to_hal_blocks(xi, yi, &aln.aligned_runs));
        }
        direct_runs.push((sp_i, xi, yi, aln));
    }

    let mut hal_query_indices = Vec::new();
    let mut hal_all_blocks = Vec::new();
    if write_hal {
        for &i in query_species_indices {
            let q_idx = hal_query_indices.len();
            hal_query_indices.push(i);
            for (xi, yi, rs, qs, len, rev) in std::mem::take(&mut hblks_by_sp[i]) {
                hal_all_blocks.push((q_idx, xi, yi, rs, qs, len, rev));
            }
        }
    }

    Ok(GuideAlignOutput {
        pairwise_mafs: Vec::new(),
        hal_query_indices,
        hal_all_blocks,
        direct_runs,
    })
}

#[allow(clippy::too_many_lines)]
fn run_streaming_alignments(
    store: &SpeciesStore,
    species_sketches: &mut [Vec<u64>],
    hinge_idx: usize,
    use_progressive: bool,
    cli: &GuideCliArgs,
    cfg: &AlignConfig,
    logger: ProgressLogger,
) -> Result<GuideAlignOutput> {
    const MIN_SHARED_SEEDS: usize = 8;
    let hinge_indices = &store.hinge(hinge_idx).1;
    let query_species_indices: Vec<usize> = (0..store.len()).filter(|&i| i != hinge_idx).collect();
    let mut pairwise_mafs = Vec::with_capacity(query_species_indices.len());
    let mut hal_query_indices = Vec::new();
    let mut hal_all_blocks = Vec::new();
    let mut sp_to_pw_idx = HashMap::with_capacity(query_species_indices.len());
    let collect_star_hal = cli.write_hal && !use_progressive;

    for (done_idx, &sp_i) in query_species_indices.iter().enumerate() {
        let sp_t0 = Instant::now();
        let qry_loaded = store.load_streaming_species(sp_i)?;
        let qry_indices = &qry_loaded.1;
        if use_progressive && species_sketches[sp_i].is_empty() {
            species_sketches[sp_i] = merge_contig_sketches(qry_indices, HINGE_SKETCH_K);
        }
        let eff_top_k = adaptive_top_k_for_fragmentation(cli.top_k, hinge_indices, qry_indices);
        let top_mask = (eff_top_k > 0)
            .then(|| select_top_k_pairs(hinge_indices, qry_indices, eff_top_k, MIN_SHARED_SEEDS));
        let mut sp_jobs = Vec::new();
        for xi in 0..hinge_indices.len() {
            for yi in 0..qry_indices.len() {
                let keep = if let Some(mask) = &top_mask {
                    mask[xi][yi]
                } else {
                    shared_seed_count(&hinge_indices[xi].1, &qry_indices[yi].1) >= MIN_SHARED_SEEDS
                };
                if keep {
                    sp_jobs.push((xi, yi));
                }
            }
        }
        let n_jobs = sp_jobs.len();
        let job_results: Vec<(String, Vec<HalBlock>)> = sp_jobs
            .par_iter()
            .map(|&(xi, yi)| {
                let (nx, ix) = &hinge_indices[xi];
                let (ny, iy) = &qry_indices[yi];
                let aln = align_contig_pair_runs(ix, iy, cfg)?;
                if aln.aligned_runs.is_empty() {
                    return Ok((String::new(), Vec::new()));
                }
                let sx = ix.seq();
                let sy = iy.seq();
                let maf = format_maf_blocks_runs(nx, sx.len(), ny, sy.len(), sx, sy, &aln);
                let hblks = if collect_star_hal {
                    runs_to_hal_blocks(xi, yi, &aln.aligned_runs)
                } else {
                    Vec::new()
                };
                Ok((maf, hblks))
            })
            .collect::<Result<Vec<_>>>()?;

        let q_idx = hal_query_indices.len();
        if cli.write_hal {
            hal_query_indices.push(sp_i);
        }
        let mut sp_maf = String::from("##maf version=1 scoring=horizon_mdl\n");
        let mut n_blocks = 0usize;
        for (maf, hblks) in job_results {
            if !maf.is_empty() {
                let body = maf
                    .strip_prefix("##maf version=1 scoring=horizon_mdl\n")
                    .unwrap_or(&maf);
                n_blocks += body.matches("\na ").count() + usize::from(body.starts_with("a "));
                sp_maf.push_str(body);
            }
            for (xi, yi, rs, qs, len, rev) in hblks {
                hal_all_blocks.push((q_idx, xi, yi, rs, qs, len, rev));
            }
        }
        logger.log(format_args!(
            "streaming pairwise: {}/{} `{}` ({} contigs, {n_jobs} pairs, {n_blocks} blocks, {:.1}s)",
            done_idx + 1,
            query_species_indices.len(),
            store.species_name(sp_i),
            qry_indices.len(),
            sp_t0.elapsed().as_secs_f64(),
        ));
        sp_to_pw_idx.insert(sp_i, pairwise_mafs.len());
        pairwise_mafs.push((store.species_name(sp_i).to_owned(), sp_maf));
        drop(qry_loaded);
    }

    if use_progressive {
        run_progressive_sister_rescue(
            store,
            species_sketches,
            hinge_idx,
            &sp_to_pw_idx,
            &mut pairwise_mafs,
            cfg,
            cli.top_k,
            logger,
        )?;
        run_ancestral_consensus_rescue(
            store,
            hinge_idx,
            &sp_to_pw_idx,
            &mut pairwise_mafs,
            cfg,
            cli.top_k,
            cli.respect_softmask,
            logger,
        )?;
    }

    Ok(GuideAlignOutput {
        pairwise_mafs,
        hal_query_indices,
        hal_all_blocks,
        direct_runs: Vec::new(),
    })
}

fn write_guide_maf_gfa(
    out_prefix: &str,
    hinge_sp: &str,
    hinge_idx: usize,
    species_contigs: &[SpeciesContigs],
    direct_runs: &[(usize, usize, usize, horizon_aligner::RunAlignment)],
    pairwise_mafs: &[(String, String)],
    logger: ProgressLogger,
) -> Result<(usize, usize)> {
    logger.log(format_args!("merging MAF/GFA -> {out_prefix}.{{maf,gfa}}"));
    let mut maf_writer =
        horizon_aligner::io_util::open_text_writer(Some(&format!("{out_prefix}.maf")))
            .with_context(|| format!("creating {out_prefix}.maf"))?;
    let mut gfa_writer =
        horizon_aligner::io_util::open_text_writer(Some(&format!("{out_prefix}.gfa")))
            .with_context(|| format!("creating {out_prefix}.gfa"))?;
    let (total_blocks, multi_blocks) = if !direct_runs.is_empty() || pairwise_mafs.is_empty() {
        let (_, hinge_contigs) = &species_contigs[hinge_idx];
        let direct_refs: Vec<horizon_aligner::DirectContigAlignment<'_>> = direct_runs
            .iter()
            .map(|(sp_i, xi, yi, aln)| {
                let (ref_contig, ref_seq) = &hinge_contigs[*xi];
                let (qry_sp, qry_contigs) = &species_contigs[*sp_i];
                let (qry_contig, qry_seq) = &qry_contigs[*yi];
                horizon_aligner::DirectContigAlignment {
                    ref_contig,
                    ref_seq,
                    qry_species: qry_sp,
                    qry_contig,
                    qry_seq,
                    aln,
                }
            })
            .collect();
        horizon_aligner::merge_direct_alignments_prefixed_to_writers(
            hinge_sp,
            &direct_refs,
            1000,
            &mut maf_writer,
            &mut gfa_writer,
        )?
    } else {
        horizon_aligner::merge_pairwise_mafs_prefixed_to_writers(
            hinge_sp,
            pairwise_mafs,
            1000,
            &mut maf_writer,
            &mut gfa_writer,
        )?
    };
    maf_writer.finish()?;
    gfa_writer.finish()?;
    logger.log(format_args!(
        "wrote MAF/GFA ({total_blocks} total blocks, {multi_blocks} multi-species)"
    ));
    Ok((total_blocks, multi_blocks))
}

#[allow(clippy::too_many_lines)]
fn main() -> Result<()> {
    let Some(cli) = parse_args()? else {
        return Ok(());
    };
    init_rayon_threads(cli.threads)?;
    let t0 = Instant::now();
    let logger = ProgressLogger::new("horizon_guide", t0, cli.quiet);

    let (store, hinge_idx, mut species_sketches) = SpeciesStore::open(&cli, logger)?;
    let use_progressive = cli.progressive && store.len() >= 3;
    let cfg = AlignConfig {
        skip_softmasked: cli.respect_softmask,
        ..AlignConfig::default()
    };

    let mode_label = if use_progressive {
        "progressive"
    } else {
        "star"
    };
    let mut align_out = if let Some(indexed) = store.in_memory() {
        // Flatten all (species, ref_contig, qry_contig) alignment jobs into a single Rayon pool
        // so slow or uneven species do not leave worker threads idle at per-species barriers.
        let (query_species_indices, all_jobs) = schedule_hinge_jobs(indexed, hinge_idx, cli.top_k);
        logger.log(format_args!(
            "hinge={}, scheduled {} pairwise contig jobs (mode={mode_label})",
            indexed[hinge_idx].0,
            all_jobs.len()
        ));
        if use_progressive {
            run_progressive_alignments(
                &store,
                &species_sketches,
                hinge_idx,
                &query_species_indices,
                &all_jobs,
                &cli,
                &cfg,
                logger,
            )?
        } else {
            run_star_alignments(
                indexed,
                hinge_idx,
                &query_species_indices,
                &all_jobs,
                cli.write_hal,
                &cfg,
                logger,
            )?
        }
    } else {
        logger.log(format_args!(
            "hinge={}, streaming {} query species (mode={mode_label})",
            store.species_name(hinge_idx),
            store.len().saturating_sub(1),
        ));
        run_streaming_alignments(
            &store,
            &mut species_sketches,
            hinge_idx,
            use_progressive,
            &cli,
            &cfg,
            logger,
        )?
    };

    // Consume heavy KP-14 seed indices and SDUST masks; in streaming mode, raw query FASTAs are
    // only reloaded when `--hal` export actually requires DNA sequences.
    let mut species_contigs = store.into_species_contigs(cli.write_hal)?;

    let hinge_sp = species_contigs[hinge_idx].0.clone();
    let (total_blocks, multi_blocks) = write_guide_maf_gfa(
        &cli.out_prefix,
        &hinge_sp,
        hinge_idx,
        &species_contigs,
        &align_out.direct_runs,
        &align_out.pairwise_mafs,
        logger,
    )?;
    align_out.direct_runs.clear();
    align_out.direct_runs.shrink_to_fit();

    let (hal_tree_desc, hal_bottom_segs) = if cli.write_hal {
        export_guide_hal(
            &cli,
            &hinge_sp,
            hinge_idx,
            &mut species_contigs,
            align_out,
            logger,
        )?
    } else {
        drop(align_out);
        ("skipped".to_owned(), 0)
    };
    println!(
        "hinge={} species={} hal_tree={} hal_bottom_segs={} total_blocks={} multi_species_blocks={} elapsed_s={:.4}",
        hinge_sp,
        species_contigs.len(),
        hal_tree_desc,
        hal_bottom_segs,
        total_blocks,
        multi_blocks,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use crate::consensus::{build_ancestral_consensus_contigs, remap_consensus_maf_to_hinge};
    use crate::rescue::rescue_transitive_maf;
    use horizon_aligner::complement_base;

    #[test]
    fn complement_and_fwd_pos_handle_both_strands() {
        assert_eq!(complement_base(b'A'), b'T');
        assert_eq!(complement_base(b'g'), b'c');
        let fwd_row = horizon_aligner::MafRow {
            src: "hg38.chr1".to_string(),
            start: 10,
            strand: '+',
            src_size: 100,
            seq: "ACGTA".to_string(),
        };
        assert_eq!(fwd_row.fwd_pos(2), 12);

        let rev_row = horizon_aligner::MafRow {
            src: "hg38.chr1".to_string(),
            start: 10,
            strand: '-',
            src_size: 100,
            seq: "ACGTA".to_string(),
        };
        assert_eq!(rev_row.fwd_pos(0), 89);
        assert_eq!(rev_row.fwd_pos(4), 85);
    }

    #[test]
    fn rescue_transitive_maf_projects_uncovered_hinge_bases() {
        let hinge_contigs = vec![("hg38.chr1", b"ACGTACGTACGTACGTACGT".as_slice())];
        let mid_contigs = vec![("chimp.chr1", b"ACGTACGTACGTACGTACGT".as_slice())];
        let hinge_to_mid = "\
a score=20
s hg38.chr1 0 20 + 20 ACGTACGTACGTACGTACGT
s chimp.chr1 0 20 + 20 ACGTACGTACGTACGTACGT
";
        let mid_to_dst = "\
a score=20
s chimp.chr1 0 20 + 20 ACGTACGTACGTACGTACGT
s dog.chr1 0 20 + 20 ACGTACGTACGTACGTACGT
";
        let rescued = rescue_transitive_maf(
            &hinge_contigs,
            &mid_contigs,
            hinge_to_mid,
            mid_to_dst,
            "",
            true,
            10,
            80,
        );
        assert!(rescued.contains("s hg38.chr1 0 20 + 20 ACGTACGTACGTACGTACGT"));
        assert!(rescued.contains("s dog.chr1 0 20 + 20 ACGTACGTACGTACGTACGT"));

        let already_covered = rescue_transitive_maf(
            &hinge_contigs,
            &mid_contigs,
            hinge_to_mid,
            mid_to_dst,
            hinge_to_mid,
            true,
            10,
            80,
        );
        assert_eq!(already_covered, "");
    }

    #[test]
    fn ancestral_consensus_votes_and_remaps_to_hinge() {
        let hinge_contigs = vec![("hg38.chr1", b"ACGTACGTACGTACGTACGT".as_slice())];
        let maf_s1 = "\
a score=20
s hg38.chr1 0 20 + 20 ACGTACGTACGTACGTACGT
s sp1.chr1 0 20 + 20 ACGTACGTACGTACGTACGT
";
        let maf_s2 = "\
a score=20
s hg38.chr1 0 20 + 20 ACGTACGTACGTACGTACGT
s sp2.chr1 0 20 + 20 ACGTACGTACGTACGTACGT
";
        let consensus =
            build_ancestral_consensus_contigs(&hinge_contigs, &[maf_s1, maf_s2], &[true]);
        assert_eq!(consensus.len(), 1);
        assert_eq!(consensus[0].0, "hg38.chr1");
        assert_eq!(consensus[0].1, b"ACGTACGTACGTACGTACGT");

        let cons_to_target = "\
a score=20
s hg38.chr1 0 20 + 20 ACGTACGTACGTACGTACGT
s sp3.chr1 0 20 + 20 ACGTACGTACGTACGTACGT
";
        let remapped = remap_consensus_maf_to_hinge(&hinge_contigs, cons_to_target, "");
        assert!(remapped.contains("s hg38.chr1 0 20 + 20 ACGTACGTACGTACGTACGT"));
        assert!(remapped.contains("s sp3.chr1 0 20 + 20 ACGTACGTACGTACGTACGT"));

        // 100% overlap with existing direct MAF is deduplicated away.
        let deduped = remap_consensus_maf_to_hinge(&hinge_contigs, cons_to_target, maf_s1);
        assert_eq!(deduped, "");

        // Mixed need_contig = [false, true] must not panic on skipped contigs.
        let mixed_contigs = vec![
            ("hg38.chr1", b"ACGTACGTACGTACGTACGT".as_slice()),
            ("hg38.chr2", b"TTTTTTTTTTTTTTTTTTTT".as_slice()),
        ];
        let mixed_cons =
            build_ancestral_consensus_contigs(&mixed_contigs, &[maf_s1, maf_s2], &[false, true]);
        assert_eq!(mixed_cons.len(), 2);
        assert_eq!(mixed_cons[0].1, [] as [u8; 0]);
        assert_eq!(mixed_cons[1].1, b"TTTTTTTTTTTTTTTTTTTT");

        // Softmasked lowercase hinge bases preserve lowercase when overridden by consensus vote.
        let soft_hinge = vec![("hg38.chr1", b"acgt".as_slice())];
        let maf_g = "a score=4\ns hg38.chr1 0 4 + 4 acgt\ns sp1.chr1 0 4 + 4 GGGG\n";
        let soft_cons = build_ancestral_consensus_contigs(&soft_hinge, &[maf_g, maf_g], &[true]);
        assert_eq!(soft_cons[0].1, b"gggg");
    }

    #[test]
    fn select_central_hinge_prefers_contiguous_reference_over_fragmented_outgroup() {
        // Species 0 (1 contig) has slightly lower raw sketch sum than Species 1 (45 contigs),
        // mimicking sacCer3 chrIII (1 contig) vs fragmented sacCas (45 contigs); auN-based
        // effective contig count must pick Species 0 cleanly, while preserving pure Jaccard
        // when all assemblies have n_eff = 1.0 AND staying resilient to minor unplaced scaffolds.
        let s0 = vec![1, 2, 3, 4, 5, 6];
        let s1 = vec![1, 2, 3, 4, 5, 6, 7];
        let s2 = vec![1, 2, 3, 7, 8, 9];
        let sketches = vec![s0, s1, s2];
        assert_eq!(super::select_central_hinge(&sketches, &[1.0, 1.0, 1.0]), 1);
        assert_eq!(
            super::select_central_hinge(&sketches, &[1.0, 45.0, 60.0]),
            0
        );

        let mut ref_lens = vec![100_000usize];
        ref_lens.extend(std::iter::repeat_n(100usize, 100));
        let draft_lens = vec![5_500usize; 20];
        let n_eff_ref = horizon_aligner::effective_contig_count(ref_lens.into_iter());
        let n_eff_draft = horizon_aligner::effective_contig_count(draft_lens.into_iter());
        assert!(n_eff_ref < 1.3 && (n_eff_draft - 20.0).abs() < 1e-6);
        assert_eq!(
            super::select_central_hinge(&sketches, &[n_eff_ref, n_eff_draft, 25.0]),
            0
        );
    }
}
