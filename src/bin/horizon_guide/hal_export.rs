//! HAL block extraction and star/progressive HDF5 `.hal` export for `horizon_guide`.

use crate::{ContigList, GuideAlignOutput, GuideCliArgs, SpeciesContigs};
use anyhow::Result;
use horizon_aligner::{ProgressLogger, hal::HalAlignment, parse_two_row_maf_blocks};
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::PathBuf;

pub(crate) type HalBlock = (usize, usize, usize, usize, usize, bool);
pub(crate) type HalMultiBlock = (usize, usize, usize, usize, usize, usize, bool);

pub(crate) fn runs_to_hal_blocks(
    xi: usize,
    yi: usize,
    runs: &[horizon_aligner::AlignedRun],
) -> Vec<HalBlock> {
    runs.iter()
        .filter(|r| r.len > 0)
        .filter_map(|r| {
            let q_start = if r.is_rc {
                (r.y_start + 1).checked_sub(r.len)?
            } else {
                r.y_start
            };
            Some((xi, yi, r.x_start, q_start, r.len, r.is_rc))
        })
        .collect()
}

pub(crate) fn export_guide_hal(
    cli: &GuideCliArgs,
    hinge_sp: &str,
    hinge_idx: usize,
    species_contigs: &mut [SpeciesContigs],
    align_out: GuideAlignOutput,
    logger: ProgressLogger,
) -> Result<(String, usize)> {
    let GuideAlignOutput {
        pairwise_mafs,
        hal_query_indices,
        hal_all_blocks,
        direct_runs: _,
    } = align_out;
    let hinge_contigs = std::mem::take(&mut species_contigs[hinge_idx].1);
    let hal_queries: Vec<(String, ContigList)> = hal_query_indices
        .into_iter()
        .map(|qi| {
            let sp = species_contigs[qi].0.clone();
            let contigs = std::mem::take(&mut species_contigs[qi].1);
            (sp, contigs)
        })
        .collect();
    let hal_path = PathBuf::from(format!("{}.hal", cli.out_prefix));
    let (desc, nb) = if cli.progressive && hal_queries.len() >= 2 {
        logger.log(format_args!(
            "building/writing HAL -> {} (mode=progressive)",
            hal_path.display()
        ));
        // Consume pairwise_mafs by value so each multi-hundred-MB MAF string is freed
        // immediately after extracting compact HAL block coordinates.
        let prog_blocks =
            extract_hal_blocks_from_pairwise_mafs(&hinge_contigs, &hal_queries, pairwise_mafs);
        HalAlignment::write_multi_contig_progressive_hdf5_hal(
            &hal_path,
            hinge_sp,
            &hinge_contigs,
            &hal_queries,
            &prog_blocks,
        )?
    } else {
        drop(pairwise_mafs);
        logger.log(format_args!(
            "building/writing HAL -> {} (mode=star)",
            hal_path.display()
        ));
        let hal = HalAlignment::from_multi_contig_star_blocks(
            hinge_sp,
            &hinge_contigs,
            &hal_queries,
            &hal_all_blocks,
        );
        drop(hal_all_blocks);
        let desc = hal.tree_newick.clone();
        let nb = hal.genomes[0].bottom_segments.len();
        hal.write_hdf5_hal(&hal_path)?;
        (desc, nb)
    };
    logger.log(format_args!("wrote HAL (bottom_segs={nb}, tree={desc})"));
    Ok((desc, nb))
}

#[allow(clippy::type_complexity)]
fn extract_hal_blocks_from_pairwise_mafs(
    hinge_contigs: &ContigList,
    hal_queries: &[(String, ContigList)],
    pairwise_mafs: Vec<(String, String)>,
) -> Vec<HalMultiBlock> {
    let hinge_idx_by_name: HashMap<&str, usize> = hinge_contigs
        .iter()
        .enumerate()
        .map(|(idx, (name, _))| (name.as_str(), idx))
        .collect();
    let per_query: Vec<Vec<HalMultiBlock>> = hal_queries
        .par_iter()
        .zip(pairwise_mafs.into_par_iter())
        .enumerate()
        .map(|(q_idx, ((_, q_contigs), (_, maf_buf)))| {
            let mut local_out = Vec::new();
            let q_idx_by_name: HashMap<&str, usize> = q_contigs
                .iter()
                .enumerate()
                .map(|(idx, (name, _))| (name.as_str(), idx))
                .collect();
            for (h_row, q_row) in parse_two_row_maf_blocks(&maf_buf) {
                let Some(&xi) = hinge_idx_by_name.get(h_row.src.as_str()) else {
                    continue;
                };
                let Some(&yi) = q_idx_by_name.get(q_row.src.as_str()) else {
                    continue;
                };
                let rev = h_row.strand != q_row.strand;
                let push_run = |out: &mut Vec<HalMultiBlock>,
                                hinge_run_start: usize,
                                qry_run_start: usize,
                                run_len: usize| {
                    let rs = if h_row.strand == '+' {
                        h_row.start + hinge_run_start
                    } else {
                        h_row
                            .src_size
                            .saturating_sub(h_row.start + hinge_run_start + run_len)
                    };
                    let qs = if q_row.strand == '+' {
                        q_row.start + qry_run_start
                    } else {
                        q_row
                            .src_size
                            .saturating_sub(q_row.start + qry_run_start + run_len)
                    };
                    out.push((q_idx, xi, yi, rs, qs, run_len, rev));
                };
                let mut h_off = 0usize;
                let mut q_off = 0usize;
                let mut hinge_run_start = 0usize;
                let mut qry_run_start = 0usize;
                let mut run_len = 0usize;
                for (hb, qb) in h_row.seq.bytes().zip(q_row.seq.bytes()) {
                    if hb != b'-' && qb != b'-' {
                        if run_len == 0 {
                            hinge_run_start = h_off;
                            qry_run_start = q_off;
                        }
                        run_len += 1;
                    } else if run_len > 0 {
                        push_run(&mut local_out, hinge_run_start, qry_run_start, run_len);
                        run_len = 0;
                    }
                    if hb != b'-' {
                        h_off += 1;
                    }
                    if qb != b'-' {
                        q_off += 1;
                    }
                }
                if run_len > 0 {
                    push_run(&mut local_out, hinge_run_start, qry_run_start, run_len);
                }
            }
            local_out
        })
        .collect();
    per_query.into_iter().flatten().collect()
}
