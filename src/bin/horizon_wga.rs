//! Multi-contig whole-genome aligner CLI using KP-14 + parallel MDL/BOCPD.

use anyhow::{Context, Result};
use horizon_aligner::{
    ContigIndex, align_contig_pair_indexed_with_min_lcb, format_gfa_blocks,
    format_gfa_paralog_tags, format_gfa_summary, format_maf_blocks, prepare_contig_index,
    resolve_alignment_duplications, shared_seed_count,
};
use rayon::prelude::*;
use std::fs;
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
}

fn parse_args() -> CliArgs {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut pos = Vec::new();
    let mut maf_out = None;
    let mut gfa_out = None;
    let mut resolve_dups = false;
    let mut min_shared_seeds = 8usize;
    let mut min_top_lcb = 0usize;

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
            "--resolve-dups" => {
                resolve_dups = true;
                i += 1;
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
    }
}

fn main() -> Result<()> {
    let cli = parse_args();
    let t0 = Instant::now();
    let contigs_x = read_fasta_contigs(&cli.fa_x, cli.max_contigs)?;
    let contigs_y = read_fasta_contigs(&cli.fa_y, cli.max_contigs)?;

    let idx_x: Vec<(String, ContigIndex)> = contigs_x
        .par_iter()
        .map(|(n, s)| Ok((n.clone(), prepare_contig_index(s)?)))
        .collect::<Result<Vec<_>>>()?;
    let idx_y: Vec<(String, ContigIndex)> = contigs_y
        .par_iter()
        .map(|(n, s)| Ok((n.clone(), prepare_contig_index(s)?)))
        .collect::<Result<Vec<_>>>()?;

    let mut total_pairs = 0usize;
    let mut total_bubbles = 0usize;
    let mut pruned_pairs = 0usize;
    let mut summary_out = String::new();
    let mut maf_buf = String::new();
    let mut gfa_buf = String::new();

    for (x_i, (nx, ix)) in idx_x.iter().enumerate() {
        for (y_i, (ny, iy)) in idx_y.iter().enumerate() {
            if shared_seed_count(ix, iy) < cli.min_shared_seeds {
                pruned_pairs += 1;
                continue;
            }
            let aln = align_contig_pair_indexed_with_min_lcb(ix, iy, cli.min_top_lcb)?;
            if aln.aligned_pairs.is_empty() {
                pruned_pairs += 1;
                continue;
            }
            total_pairs += aln.aligned_pairs.len();
            total_bubbles += aln.macro_bubbles + aln.internal_bubbles;
            summary_out.push_str(&format_gfa_summary(nx, ny, &aln));

            if cli.maf_out.is_some() {
                let sx = &contigs_x[x_i].1;
                let sy = &contigs_y[y_i].1;
                maf_buf.push_str(&format_maf_blocks(nx, sx.len(), ny, sy.len(), sx, sy, &aln));
            }
            if cli.gfa_out.is_some() {
                gfa_buf.push_str(&format_gfa_blocks(nx, ny, &aln));
            }
            if cli.resolve_dups {
                let dups = resolve_alignment_duplications(&aln, 1000);
                let tag_line = format_gfa_paralog_tags(nx, ny, &dups);
                if !tag_line.is_empty() {
                    summary_out.push_str(&tag_line);
                    if cli.gfa_out.is_some() {
                        gfa_buf.push_str(&tag_line);
                    }
                }
            }
        }
    }

    if let Some(maf_path) = &cli.maf_out {
        fs::write(maf_path, &maf_buf).with_context(|| format!("writing {}", maf_path.display()))?;
    }
    if let Some(gfa_path) = &cli.gfa_out {
        fs::write(gfa_path, &gfa_buf).with_context(|| format!("writing {}", gfa_path.display()))?;
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
