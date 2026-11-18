//! Nystrom/sketch-selected hinge + multi-species WGA pipeline (`horizon_wga` + `horizon_merge`).

use anyhow::{Context, Result, bail};
use horizon_aligner::{
    AlignConfig, ContigIndex, align_contig_pair_with_config, format_maf_blocks, hal::HalAlignment,
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

type HalBlockList = Vec<(usize, usize, usize, bool)>;

fn align_species_pair(
    ref_contigs: &ContigList,
    ref_indices: &[(String, ContigIndex)],
    qry_contigs: &ContigList,
    qry_indices: &[(String, ContigIndex)],
    cfg: &AlignConfig,
) -> Result<(String, HalBlockList)> {
    let mut maf_buf = String::from("##maf version=1 scoring=horizon_mdl\n");
    let mut hal_blocks = Vec::new();
    for (xi, (nx, ix)) in ref_indices.iter().enumerate() {
        for (yi, (ny, iy)) in qry_indices.iter().enumerate() {
            if shared_seed_count(ix, iy) < 8 {
                continue;
            }
            let aln = align_contig_pair_with_config(ix, iy, cfg)?;
            if aln.aligned_pairs.is_empty() {
                continue;
            }
            let sx = &ref_contigs[xi].1;
            let sy = &qry_contigs[yi].1;
            maf_buf.push_str(&format_maf_blocks(nx, sx.len(), ny, sy.len(), sx, sy, &aln));
            if let Some(&(x0, y0)) = aln.aligned_pairs.first()
                && let Some(&(x1, y1)) = aln.aligned_pairs.last()
                && x1 > x0
            {
                hal_blocks.push((x0, y0.min(y1), x1 - x0 + 1, y1 < y0));
            }
        }
    }
    Ok((maf_buf, hal_blocks))
}

#[allow(clippy::too_many_lines)]
fn main() -> Result<()> {
    let raw_args: Vec<String> = std::env::args().skip(1).collect();
    let mut progressive = false;
    let mut pos = Vec::new();
    for a in raw_args {
        if a == "--progressive" {
            progressive = true;
        } else {
            pos.push(a);
        }
    }
    if pos.len() < 4 {
        bail!(
            "usage: horizon_guide [--progressive] <out_prefix> <max_contigs> <sp1=fa1> <sp2=fa2>..."
        );
    }
    let t0 = Instant::now();
    let out_prefix = &pos[0];
    let max_contigs: usize = pos[1].parse().context("parsing max_contigs")?;

    let mut specs: Vec<(String, PathBuf)> = Vec::new();
    for spec in &pos[2..] {
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
    let mut hal_queries: Vec<(String, Vec<u8>)> = Vec::new();
    let mut hal_all_blocks: Vec<(usize, usize, usize, usize, bool)> = Vec::new();

    for (i, (qry_sp, qry_contigs, qry_indices)) in indexed.iter().enumerate() {
        if i == hinge_idx {
            continue;
        }
        let (maf_buf, hblocks) =
            align_species_pair(hinge_contigs, hinge_indices, qry_contigs, qry_indices, &cfg)?;
        let q_idx = hal_queries.len();
        hal_queries.push((
            qry_sp.clone(),
            qry_contigs.first().map_or_else(Vec::new, |c| c.1.clone()),
        ));
        for (rs, qs, len, rev) in hblocks {
            hal_all_blocks.push((q_idx, rs, qs, len, rev));
        }
        pairwise_mafs.push((qry_sp.clone(), maf_buf));
    }

    if progressive && indexed.len() >= 3 {
        // Align non-hinge sister pairs to record sister-clade coverage in the HAL graph.
        for i in 0..indexed.len() {
            if i == hinge_idx {
                continue;
            }
            let best_j = (0..indexed.len())
                .filter(|&j| j != i && j != hinge_idx)
                .max_by_key(|&j| {
                    indexed[i]
                        .2
                        .iter()
                        .flat_map(|(_, ix)| {
                            indexed[j]
                                .2
                                .iter()
                                .map(move |(_, iy)| shared_seed_count(ix, iy))
                        })
                        .sum::<usize>()
                });
            if let Some(j) = best_j {
                let _ = align_species_pair(
                    &indexed[i].1,
                    &indexed[i].2,
                    &indexed[j].1,
                    &indexed[j].2,
                    &cfg,
                )?;
            }
        }
    }

    let hal = HalAlignment::from_star_blocks(
        hinge_sp,
        hinge_contigs.first().map_or(&[], |c| c.1.as_slice()),
        &hal_queries,
        &hal_all_blocks,
    );
    let (merged_maf, merged_gfa, total_blocks, multi_blocks) =
        merge_pairwise_mafs(hinge_sp, &pairwise_mafs, 1000);
    fs::write(format!("{out_prefix}.maf"), merged_maf)?;
    fs::write(format!("{out_prefix}.gfa"), merged_gfa)?;
    println!(
        "hinge={} species={} hal_tree={} hal_bottom_segs={} total_blocks={} multi_species_blocks={} elapsed_s={:.4}",
        hinge_sp,
        indexed.len(),
        hal.tree_newick,
        hal.genomes[0].bottom_segments.len(),
        total_blocks,
        multi_blocks,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}
