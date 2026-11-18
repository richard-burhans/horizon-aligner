//! Multi-contig whole-genome aligner CLI using KP-14 + parallel MDL/BOCPD.

use anyhow::{Context, Result};
use horizon_aligner::{
    ContigIndex, align_contig_pair_indexed, format_gfa_summary, prepare_contig_index,
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

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let (fa_x, fa_y, max_contigs) = if args.len() >= 3 {
        let mc = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1);
        (PathBuf::from(&args[1]), PathBuf::from(&args[2]), mc)
    } else {
        let base = "/home/agent/.gemini/antigravity-cli/brain/c77b5ff5-88ad-43af-973a-c37dbadc186d/scratch/plasmodium_fastas";
        (
            Path::new(base).join("GCA_032713075.1__Plasmodium_falciparum.fa"),
            Path::new(base).join("GCA_914969965.1__Plasmodium_vivax.fa"),
            1,
        )
    };

    let t0 = Instant::now();
    let contigs_x = read_fasta_contigs(&fa_x, max_contigs)?;
    let contigs_y = read_fasta_contigs(&fa_y, max_contigs)?;

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
    let mut gfa_out = String::new();

    for (nx, ix) in &idx_x {
        for (ny, iy) in &idx_y {
            let aln = align_contig_pair_indexed(ix, iy)?;
            if !aln.aligned_pairs.is_empty() {
                total_pairs += aln.aligned_pairs.len();
                total_bubbles += aln.macro_bubbles + aln.internal_bubbles;
                gfa_out.push_str(&format_gfa_summary(nx, ny, &aln));
            }
        }
    }

    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "contigs_x={} contigs_y={} aligned_pairs={} total_bubbles={} elapsed_s={:.4}",
        contigs_x.len(),
        contigs_y.len(),
        total_pairs,
        total_bubbles,
        elapsed
    );
    if !gfa_out.is_empty() {
        print!("{gfa_out}");
    }
    Ok(())
}
