//! Multi-contig whole-genome aligner CLI using KP-14 + parallel MDL/BOCPD.

use anyhow::{Context, Result};
use horizon_aligner::{
    AlignConfig, ContigIndex, align_contig_pair_with_config, format_gfa_blocks,
    format_gfa_paralog_tags, format_gfa_summary, format_maf_blocks, prepare_contig_index_opts,
    resolve_alignment_duplications, score_contig_pair_lcb, shared_seed_count,
};
use rayon::prelude::*;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

fn read_fasta_contigs(path: &Path, max_contigs: usize) -> Result<Vec<(String, Vec<u8>)>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut out = Vec::new();
    let mut cur_name: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = cur_name.take() {
                out.push((name, std::mem::take(&mut cur_seq)));
                if out.len() >= max_contigs {
                    return Ok(out);
                }
            }
            let header = rest
                .split_whitespace()
                .next()
                .unwrap_or("contig")
                .to_owned();
            cur_name = Some(header);
        } else if cur_name.is_some() {
            cur_seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some(name) = cur_name {
        out.push((name, cur_seq));
    }
    Ok(out)
}

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
}

fn parse_args() -> CliArgs {
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

    let base = "/home/agent/.gemini/antigravity-cli/brain/c77b5ff5-88ad-43af-973a-c37dbadc186d/scratch/plasmodium_fastas";
    let fa_x = pos.first().map_or_else(
        || Path::new(base).join("GCA_032713075.1__Plasmodium_falciparum.fa"),
        PathBuf::from,
    );
    let fa_y = pos.get(1).map_or_else(
        || Path::new(base).join("GCA_914969965.1__Plasmodium_vivax.fa"),
        PathBuf::from,
    );
    let max_contigs = pos.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);

    CliArgs {
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
    }
}

fn compute_top_k_mask(
    target_idx: &[(String, ContigIndex)],
    query_idx: &[(String, ContigIndex)],
    min_shared_seeds: usize,
    top_k: usize,
) -> Vec<Vec<bool>> {
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

#[allow(clippy::too_many_lines)]
fn main() -> Result<()> {
    let cli = parse_args();
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
        .then(|| compute_top_k_mask(&idx_x, &idx_y, cli.min_shared_seeds, cli.top_k));

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
