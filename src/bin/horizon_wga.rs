//! Multi-contig whole-genome aligner CLI using KP-14 + parallel MDL/BOCPD.

use anyhow::{Context, Result, bail};
use horizon_aligner::{
    AlignConfig, ContigIndex, align_contig_pair_with_config, format_gfa_blocks,
    format_gfa_paralog_tags, format_gfa_summary, format_maf_blocks, prepare_contig_index_opts,
    read_fasta_contigs, resolve_alignment_duplications, select_top_k_pairs, shared_seed_count,
};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

struct CliArgs {
    fa_x: PathBuf,
    fa_y: PathBuf,
    max_contigs: usize,
    maf_out: Option<PathBuf>,
    gfa_out: Option<PathBuf>,
    resolve_dups: bool,
    min_shared_seeds: usize,
    min_top_lcb: usize,
    top_k: usize,
    respect_softmask: bool,
    use_ry_seeds: bool,
    hal_info: Option<PathBuf>,
    min_fallback_cov: Option<f64>,
    hsp_threshold: Option<i32>,
}

fn parse_args() -> Result<CliArgs> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut pos = Vec::new();
    let mut maf_out = None;
    let mut gfa_out = None;
    let mut resolve_dups = false;
    let mut min_shared_seeds = 8usize;
    let mut min_top_lcb = 0usize;
    let mut top_k = 0usize;
    let mut respect_softmask = false;
    let mut use_ry_seeds = false;
    let mut hal_info = None;
    let mut min_fallback_cov = None;
    let mut hsp_threshold = None;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--maf" if i + 1 < args.len() => {
                maf_out = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--gfa" if i + 1 < args.len() => {
                gfa_out = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            "--min-shared-seeds" if i + 1 < args.len() => {
                if let Ok(v) = args[i + 1].parse() {
                    min_shared_seeds = v;
                }
                i += 2;
            }
            "--min-top-lcb" if i + 1 < args.len() => {
                if let Ok(v) = args[i + 1].parse() {
                    min_top_lcb = v;
                }
                i += 2;
            }
            "--top-k" if i + 1 < args.len() => {
                if let Ok(v) = args[i + 1].parse() {
                    top_k = v;
                }
                i += 2;
            }
            "--resolve-dups" => {
                resolve_dups = true;
                i += 1;
            }
            "--rbest" => {
                top_k = 1;
                i += 1;
            }
            "--respect-softmask" => {
                respect_softmask = true;
                i += 1;
            }
            "--ry-seeds" => {
                use_ry_seeds = true;
                i += 1;
            }
            "--no-ry-seeds" => {
                use_ry_seeds = false;
                i += 1;
            }
            "--min-fallback-cov" if i + 1 < args.len() => {
                if let Ok(v) = args[i + 1].parse() {
                    min_fallback_cov = Some(v);
                }
                i += 2;
            }
            "--hsp-threshold" if i + 1 < args.len() => {
                if let Ok(v) = args[i + 1].parse() {
                    hsp_threshold = Some(v);
                }
                i += 2;
            }

            "--hal-info" if i + 1 < args.len() => {
                hal_info = Some(PathBuf::from(&args[i + 1]));
                i += 2;
            }
            _ => {
                pos.push(args[i].clone());
                i += 1;
            }
        }
    }

    if pos.len() < 2 && hal_info.is_none() {
        bail!(
            "usage: horizon_wga <target.fa> <query.fa> [max_contigs] [--maf out.maf] [--gfa out.gfa] [--hal-info file.hal] [--top-k K] [--rbest] [--resolve-dups] [--respect-softmask] [--ry-seeds|--no-ry-seeds]"
        );
    }
    let fa_x = pos.first().map_or_default(PathBuf::from);
    let fa_y = pos.get(1).map_or_default(PathBuf::from);
    let max_contigs = pos.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);

    Ok(CliArgs {
        fa_x,
        fa_y,
        max_contigs,
        maf_out,
        gfa_out,
        resolve_dups,
        min_shared_seeds,
        min_top_lcb,
        top_k,
        respect_softmask,
        use_ry_seeds,
        hal_info,
        min_fallback_cov,
        hsp_threshold,
    })
}

#[allow(clippy::too_many_lines)]
fn main() -> Result<()> {
    let cli = parse_args()?;
    if let Some(hal_path) = &cli.hal_info {
        let info = horizon_aligner::hal::inspect_hal_genomes(hal_path)
            .with_context(|| format!("inspecting {}", hal_path.display()))?
            .with_context(|| format!("{} is not a valid HDF5 HAL file", hal_path.display()))?;
        println!(
            "hal_file={} hdf5_v={} offset_bytes={} len_bytes={} root_members={}",
            hal_path.display(),
            info.superblock.version,
            info.superblock.offset_size,
            info.superblock.length_size,
            info.root_names.len()
        );
        if !info.root_names.is_empty() {
            println!(
                "sample_members={}",
                info.root_names[..info.root_names.len().min(10)].join(",")
            );
        }
        return Ok(());
    }
    let t0 = Instant::now();
    let contigs_x = read_fasta_contigs(&cli.fa_x, cli.max_contigs)?;
    let contigs_y = read_fasta_contigs(&cli.fa_y, cli.max_contigs)?;

    let idx_x: Vec<(String, ContigIndex)> = contigs_x
        .par_iter()
        .map(|(n, s)| {
            Ok((
                n.clone(),
                prepare_contig_index_opts(s, cli.respect_softmask)?,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    let idx_y: Vec<(String, ContigIndex)> = contigs_y
        .par_iter()
        .map(|(n, s)| {
            Ok((
                n.clone(),
                prepare_contig_index_opts(s, cli.respect_softmask)?,
            ))
        })
        .collect::<Result<Vec<_>>>()?;

    let cfg = AlignConfig {
        min_top_lcb_anchors: cli.min_top_lcb,
        skip_softmasked: cli.respect_softmask,
        use_ry_seeds: cli.use_ry_seeds,
        min_fallback_cov: cli.min_fallback_cov,
        hsp_threshold: cli.hsp_threshold,
        ..AlignConfig::default()
    };

    let mut total_pairs = 0usize;
    let mut total_bubbles = 0usize;
    let mut pruned_pairs = 0usize;
    let mut summary_out = String::new();
    let mut maf_writer = cli
        .maf_out
        .as_ref()
        .map(|p| File::create(p).map(BufWriter::new))
        .transpose()
        .context("creating MAF output")?;
    let mut gfa_writer = cli
        .gfa_out
        .as_ref()
        .map(|p| File::create(p).map(BufWriter::new))
        .transpose()
        .context("creating GFA output")?;

    let allowed = (cli.top_k > 0)
        .then(|| select_top_k_pairs(&idx_x, &idx_y, cli.top_k, cli.min_shared_seeds));

    for (x_i, (nx, ix)) in idx_x.iter().enumerate() {
        for (y_i, (ny, iy)) in idx_y.iter().enumerate() {
            if let Some(mask) = allowed.as_ref()
                && !mask[x_i][y_i]
            {
                pruned_pairs += 1;
                continue;
            }
            if shared_seed_count(ix, iy) < cli.min_shared_seeds {
                pruned_pairs += 1;
                continue;
            }
            let aln = align_contig_pair_with_config(ix, iy, &cfg)?;
            if aln.aligned_pairs.is_empty() {
                pruned_pairs += 1;
                continue;
            }
            total_pairs += aln.aligned_pairs.len();
            total_bubbles += aln.macro_bubbles + aln.internal_bubbles;
            summary_out.push_str(&format_gfa_summary(nx, ny, &aln));

            if let Some(w) = maf_writer.as_mut() {
                let sx = &contigs_x[x_i].1;
                let sy = &contigs_y[y_i].1;
                let text = format_maf_blocks(nx, sx.len(), ny, sy.len(), sx, sy, &aln);
                w.write_all(text.as_bytes())?;
            }
            if let Some(w) = gfa_writer.as_mut() {
                let text = format_gfa_blocks(nx, ny, &aln);
                w.write_all(text.as_bytes())?;
            }
            if cli.resolve_dups {
                let dups = resolve_alignment_duplications(&aln, 1000);
                let tag_line = format_gfa_paralog_tags(nx, ny, &dups);
                if !tag_line.is_empty() {
                    summary_out.push_str(&tag_line);
                    if let Some(w) = gfa_writer.as_mut() {
                        w.write_all(tag_line.as_bytes())?;
                    }
                }
            }
        }
    }

    if let Some(w) = maf_writer.as_mut() {
        w.flush()?;
    }
    if let Some(w) = gfa_writer.as_mut() {
        w.flush()?;
    }

    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "contigs_x={} contigs_y={} pruned_pairs={} aligned_pairs={} total_bubbles={} elapsed_s={:.4}",
        contigs_x.len(),
        contigs_y.len(),
        pruned_pairs,
        total_pairs,
        total_bubbles,
        elapsed
    );
    if !summary_out.is_empty() {
        print!("{summary_out}");
    }
    Ok(())
}
