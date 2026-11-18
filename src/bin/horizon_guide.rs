//! Nystrom/sketch-selected hinge + multi-species WGA pipeline (`horizon_wga` + `horizon_merge`).

use anyhow::{Context, Result, bail};
use horizon_aligner::{
    AlignConfig, ContigIndex, align_contig_pair_with_config, format_maf_blocks,
    merge_pairwise_mafs, prepare_contig_index_opts, shared_seed_count,
};
use rayon::prelude::*;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

type ContigList = Vec<(String, Vec<u8>)>;
type IndexedSpecies = (String, ContigList, Vec<(String, ContigIndex)>);

fn read_fasta_contigs(path: &Path, max_contigs: usize) -> Result<ContigList> {
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

fn select_central_hinge(species: &[IndexedSpecies]) -> usize {
    let mut best_idx = 0usize;
    let mut best_score = 0usize;
    for (i, (_, _, idx_i)) in species.iter().enumerate() {
        let mut total = 0usize;
        for (j, (_, _, idx_j)) in species.iter().enumerate() {
            if i == j {
                continue;
            }
            for (_, ix) in idx_i {
                for (_, iy) in idx_j {
                    total += shared_seed_count(ix, iy);
                }
            }
        }
        if total > best_score {
            best_score = total;
            best_idx = i;
        }
    }
    best_idx
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 4 {
        bail!("usage: horizon_guide <out_prefix> <max_contigs> <sp1=fa1> <sp2=fa2>...");
    }
    let t0 = Instant::now();
    let out_prefix = &args[0];
    let max_contigs: usize = args[1].parse().context("parsing max_contigs")?;

    let mut specs: Vec<(String, PathBuf)> = Vec::new();
    for spec in &args[2..] {
        let Some((sp, path)) = spec.split_once('=') else {
            bail!("expected <species>=<fasta>, got {spec}");
        };
        specs.push((sp.to_owned(), PathBuf::from(path)));
    }

    let indexed: Vec<IndexedSpecies> = specs
        .par_iter()
        .map(|(sp, path)| {
            let contigs = read_fasta_contigs(path, max_contigs)?;
            let idx = contigs
                .iter()
                .map(|(n, s)| Ok((n.clone(), prepare_contig_index_opts(s, false)?)))
                .collect::<Result<Vec<_>>>()?;
            Ok((sp.clone(), contigs, idx))
        })
        .collect::<Result<Vec<_>>>()?;

    let hinge_idx = select_central_hinge(&indexed);
    let (hinge_sp, hinge_contigs, hinge_indices) = &indexed[hinge_idx];
    let cfg = AlignConfig::default();

    let mut pairwise_mafs: Vec<(String, String)> = Vec::new();
    for (i, (qry_sp, qry_contigs, qry_indices)) in indexed.iter().enumerate() {
        if i == hinge_idx {
            continue;
        }
        let mut maf_buf = String::from("##maf version=1 scoring=horizon_mdl\n");
        for (xi, (nx, ix)) in hinge_indices.iter().enumerate() {
            for (yi, (ny, iy)) in qry_indices.iter().enumerate() {
                if shared_seed_count(ix, iy) < 8 {
                    continue;
                }
                let aln = align_contig_pair_with_config(ix, iy, &cfg)?;
                if aln.aligned_pairs.is_empty() {
                    continue;
                }
                let sx = &hinge_contigs[xi].1;
                let sy = &qry_contigs[yi].1;
                maf_buf.push_str(&format_maf_blocks(nx, sx.len(), ny, sy.len(), sx, sy, &aln));
            }
        }
        pairwise_mafs.push((qry_sp.clone(), maf_buf));
    }

    let (merged_maf, merged_gfa, total_blocks, multi_blocks) =
        merge_pairwise_mafs(hinge_sp, &pairwise_mafs, 1000);
    fs::write(format!("{out_prefix}.maf"), merged_maf)?;
    fs::write(format!("{out_prefix}.gfa"), merged_gfa)?;
    println!(
        "hinge={} species={} total_blocks={} multi_species_blocks={} elapsed_s={:.4}",
        hinge_sp,
        indexed.len(),
        total_blocks,
        multi_blocks,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}
