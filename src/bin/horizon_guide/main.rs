//! Nystrom/sketch-selected hinge + progressive multi-species WGA/HAL pipeline (`horizon_guide`).

use anyhow::{Context, Result, bail};
use horizon_aligner::{
    AlignConfig, ContigIndex, align_contig_pair_runs, format_maf_blocks_runs, hal::HalAlignment,
    parse_two_row_maf_blocks, prepare_contig_index_opts, read_fasta_contigs, select_top_k_pairs,
    shared_seed_count, shared_sketch_count,
};
use rayon::prelude::*;
use rescue::{run_ancestral_consensus_rescue, run_progressive_sister_rescue};
use std::collections::HashMap;
use std::path::PathBuf;
use std::time::Instant;

mod rescue;

pub(crate) type ContigList = Vec<(String, Vec<u8>)>;
pub(crate) type IndexedSpecies = (String, ContigList, Vec<(String, ContigIndex)>);

fn compute_species_sketches(species: &[IndexedSpecies]) -> Vec<Vec<u64>> {
    const HINGE_SKETCH_K: usize = 512;
    species
        .par_iter()
        .map(|(_, _, idx_list)| {
            use std::collections::BinaryHeap;
            let mut heap: BinaryHeap<u64> = BinaryHeap::with_capacity(HINGE_SKETCH_K);
            for (_, ix) in idx_list {
                for h in ix.minhash_sketch(HINGE_SKETCH_K) {
                    if heap.len() < HINGE_SKETCH_K {
                        heap.push(h);
                    } else if let Some(&top) = heap.peek()
                        && h < top
                    {
                        heap.pop();
                        heap.push(h);
                    }
                }
            }
            let mut sorted = heap.into_sorted_vec();
            sorted.dedup();
            sorted
        })
        .collect()
}

fn select_central_hinge(sketches: &[Vec<u64>]) -> usize {
    (0..sketches.len())
        .into_par_iter()
        .map(|i| {
            let total: usize = (0..sketches.len())
                .filter(|&j| j != i)
                .map(|j| shared_sketch_count(&sketches[i], &sketches[j]))
                .sum();
            (total, std::cmp::Reverse(i))
        })
        .max()
        .map_or(0, |(_, std::cmp::Reverse(i))| i)
}

type HalBlock = (usize, usize, usize, usize, usize, bool);
type HalBlockList = Vec<HalBlock>;

fn runs_to_hal_blocks(xi: usize, yi: usize, runs: &[horizon_aligner::AlignedRun]) -> Vec<HalBlock> {
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

pub(crate) fn align_species_pair(
    ref_contigs: &ContigList,
    ref_indices: &[(String, ContigIndex)],
    qry_contigs: &ContigList,
    qry_indices: &[(String, ContigIndex)],
    cfg: &AlignConfig,
    top_k: usize,
) -> Result<(String, HalBlockList)> {
    let min_shared_seeds = 8usize;
    let top_mask = if top_k > 0 {
        Some(select_top_k_pairs(
            ref_indices,
            qry_indices,
            top_k,
            min_shared_seeds,
        ))
    } else {
        None
    };

    let pairs: Vec<(usize, usize)> = (0..ref_indices.len())
        .flat_map(|xi| (0..qry_indices.len()).map(move |yi| (xi, yi)))
        .filter(|&(xi, yi)| {
            if let Some(mask) = &top_mask {
                mask[xi][yi]
            } else {
                shared_seed_count(&ref_indices[xi].1, &qry_indices[yi].1) >= min_shared_seeds
            }
        })
        .collect();

    let results: Vec<(String, Vec<HalBlock>)> = pairs
        .par_iter()
        .map(|&(xi, yi)| {
            let (nx, ix) = &ref_indices[xi];
            let (ny, iy) = &qry_indices[yi];
            let aln = align_contig_pair_runs(ix, iy, cfg)?;
            if aln.aligned_runs.is_empty() {
                return Ok((String::new(), Vec::new()));
            }
            let sx = &ref_contigs[xi].1;
            let sy = &qry_contigs[yi].1;
            let maf = format_maf_blocks_runs(nx, sx.len(), ny, sy.len(), sx, sy, &aln);
            let hblks = runs_to_hal_blocks(xi, yi, &aln.aligned_runs);
            Ok((maf, hblks))
        })
        .collect::<Result<Vec<_>>>()?;

    let mut maf_buf = String::from("##maf version=1 scoring=horizon_mdl\n");
    let mut hal_blocks = Vec::new();
    for (maf, hblks) in results {
        if !maf.is_empty() {
            let body = maf
                .strip_prefix("##maf version=1 scoring=horizon_mdl\n")
                .unwrap_or(&maf);
            maf_buf.push_str(body);
        }
        hal_blocks.extend(hblks);
    }
    Ok((maf_buf, hal_blocks))
}

const GUIDE_USAGE: &str = "usage: horizon_guide [--progressive] [--respect-softmask] [--hal] [--rbest | --top-k K | --all-pairs] <out_prefix> <max_contigs> <sp1=fa1> <sp2=fa2>...";

struct GuideCliArgs {
    progressive: bool,
    respect_softmask: bool,
    write_hal: bool,
    top_k: usize,
    out_prefix: String,
    max_contigs: usize,
    specs: Vec<(String, PathBuf)>,
}

fn parse_args() -> Result<Option<GuideCliArgs>> {
    let raw_args: Vec<String> = std::env::args().skip(1).collect();
    let mut progressive = false;
    let mut respect_softmask = false;
    let mut write_hal = false;
    let mut top_k_opt: Option<usize> = None;
    let mut pos = Vec::new();
    let mut arg_i = 0usize;
    while arg_i < raw_args.len() {
        match raw_args[arg_i].as_str() {
            "-h" | "--help" => {
                println!("{GUIDE_USAGE}");
                return Ok(None);
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
            "--rbest" => {
                top_k_opt = Some(1);
                arg_i += 1;
            }
            "--all-pairs" => {
                top_k_opt = Some(0);
                arg_i += 1;
            }
            "--top-k" => {
                let Some(val) = raw_args.get(arg_i + 1) else {
                    bail!("--top-k requires an integer argument");
                };
                let parsed = val
                    .parse::<usize>()
                    .with_context(|| format!("invalid --top-k value `{val}`"))?;
                top_k_opt = Some(parsed);
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
    let top_k = top_k_opt.unwrap_or(usize::from(max_contigs != 1));
    let mut specs = Vec::new();
    let mut seen_species = std::collections::HashSet::new();
    for spec in &pos[2..] {
        let Some((sp, path)) = spec.split_once('=') else {
            bail!("expected <species>=<fasta>, got {spec}");
        };
        anyhow::ensure!(
            !sp.is_empty() && !path.is_empty(),
            "expected non-empty <species>=<fasta>, got {spec}"
        );
        anyhow::ensure!(
            seen_species.insert(sp.to_owned()),
            "duplicate species name: {sp}"
        );
        specs.push((sp.to_owned(), PathBuf::from(path)));
    }
    Ok(Some(GuideCliArgs {
        progressive,
        respect_softmask,
        write_hal,
        top_k,
        out_prefix,
        max_contigs,
        specs,
    }))
}

#[allow(clippy::too_many_lines)]
fn main() -> Result<()> {
    let Some(cli) = parse_args()? else {
        return Ok(());
    };
    let t0 = Instant::now();

    let mut indexed: Vec<IndexedSpecies> = cli
        .specs
        .par_iter()
        .map(|(sp, path)| {
            let contigs = read_fasta_contigs(path, cli.max_contigs)?;
            let idx = contigs
                .par_iter()
                .map(|(n, s)| {
                    Ok((
                        n.clone(),
                        prepare_contig_index_opts(s, cli.respect_softmask)?,
                    ))
                })
                .collect::<Result<Vec<_>>>()?;
            Ok((sp.clone(), contigs, idx))
        })
        .collect::<Result<Vec<_>>>()?;

    let species_sketches = compute_species_sketches(&indexed);
    let hinge_idx = select_central_hinge(&species_sketches);
    let mut pairwise_mafs: Vec<(String, String)> = Vec::new();
    let mut hal_query_indices: Vec<usize> = Vec::new();
    let mut hal_all_blocks: Vec<(usize, usize, usize, usize, usize, usize, bool)> = Vec::new();
    let mut sp_to_pw_idx: HashMap<usize, usize> = HashMap::new();
    let mut direct_runs: Vec<(usize, usize, usize, horizon_aligner::RunAlignment)> = Vec::new();

    {
        let (_, hinge_contigs, hinge_indices) = &indexed[hinge_idx];
        let cfg = AlignConfig {
            skip_softmasked: cli.respect_softmask,
            ..AlignConfig::default()
        };

        // Flatten all (species, ref_contig, qry_contig) alignment jobs into a single Rayon pool
        // so slow or uneven species do not leave worker threads idle at per-species barriers.
        let min_shared_seeds = 8usize;
        let query_species_indices: Vec<usize> =
            (0..indexed.len()).filter(|&i| i != hinge_idx).collect();
        let all_jobs: Vec<(usize, usize, usize)> = query_species_indices
            .par_iter()
            .flat_map(|&sp_i| {
                let qry_indices = &indexed[sp_i].2;
                let top_mask = (cli.top_k > 0).then(|| {
                    select_top_k_pairs(hinge_indices, qry_indices, cli.top_k, min_shared_seeds)
                });
                let mut sp_jobs = Vec::new();
                for xi in 0..hinge_indices.len() {
                    for yi in 0..qry_indices.len() {
                        let keep = if let Some(mask) = &top_mask {
                            mask[xi][yi]
                        } else {
                            shared_seed_count(&hinge_indices[xi].1, &qry_indices[yi].1)
                                >= min_shared_seeds
                        };
                        if keep {
                            sp_jobs.push((sp_i, xi, yi));
                        }
                    }
                }
                sp_jobs
            })
            .collect();

        if cli.progressive && indexed.len() >= 3 {
            let job_results: Vec<(usize, String)> = all_jobs
                .par_iter()
                .map(|&(sp_i, xi, yi)| {
                    let (_, qry_contigs, qry_indices) = &indexed[sp_i];
                    let (nx, ix) = &hinge_indices[xi];
                    let (ny, iy) = &qry_indices[yi];
                    let aln = align_contig_pair_runs(ix, iy, &cfg)?;
                    if aln.aligned_runs.is_empty() {
                        return Ok((sp_i, String::new()));
                    }
                    let sx = &hinge_contigs[xi].1;
                    let sy = &qry_contigs[yi].1;
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

            for i in query_species_indices {
                let qry_sp = &indexed[i].0;
                if cli.write_hal {
                    hal_query_indices.push(i);
                }
                sp_to_pw_idx.insert(i, pairwise_mafs.len());
                pairwise_mafs.push((qry_sp.clone(), std::mem::take(&mut maf_by_sp[i])));
            }

            run_progressive_sister_rescue(
                &indexed,
                &species_sketches,
                hinge_idx,
                &sp_to_pw_idx,
                &mut pairwise_mafs,
                &cfg,
                cli.top_k,
            )?;
            run_ancestral_consensus_rescue(
                &indexed,
                hinge_idx,
                &sp_to_pw_idx,
                &mut pairwise_mafs,
                &cfg,
                cli.top_k,
                cli.respect_softmask,
            )?;
        } else {
            let job_results: Vec<(usize, usize, usize, horizon_aligner::RunAlignment)> = all_jobs
                .par_iter()
                .map(|&(sp_i, xi, yi)| {
                    let qry_indices = &indexed[sp_i].2;
                    let aln =
                        align_contig_pair_runs(&hinge_indices[xi].1, &qry_indices[yi].1, &cfg)?;
                    Ok((sp_i, xi, yi, aln))
                })
                .collect::<Result<Vec<_>>>()?;

            let mut hblks_by_sp: Vec<Vec<HalBlock>> = vec![Vec::new(); indexed.len()];
            for (sp_i, xi, yi, aln) in job_results {
                if aln.aligned_runs.is_empty() {
                    continue;
                }
                if cli.write_hal {
                    hblks_by_sp[sp_i].extend(runs_to_hal_blocks(xi, yi, &aln.aligned_runs));
                }
                direct_runs.push((sp_i, xi, yi, aln));
            }
            if cli.write_hal {
                for i in query_species_indices {
                    let q_idx = hal_query_indices.len();
                    hal_query_indices.push(i);
                    for (xi, yi, rs, qs, len, rev) in std::mem::take(&mut hblks_by_sp[i]) {
                        hal_all_blocks.push((q_idx, xi, yi, rs, qs, len, rev));
                    }
                }
            }
        }
    }

    // Drop heavy KP-14 seed indices, SDUST masks, and RC sequences before merging & HAL export.
    for (_, _, idxs) in &mut indexed {
        idxs.clear();
        idxs.shrink_to_fit();
    }

    let hinge_sp = indexed[hinge_idx].0.clone();
    let out_prefix = &cli.out_prefix;
    let mut maf_writer =
        horizon_aligner::io_util::open_text_writer(Some(&format!("{out_prefix}.maf")))
            .with_context(|| format!("creating {out_prefix}.maf"))?;
    let mut gfa_writer =
        horizon_aligner::io_util::open_text_writer(Some(&format!("{out_prefix}.gfa")))
            .with_context(|| format!("creating {out_prefix}.gfa"))?;
    let (total_blocks, multi_blocks) = if !direct_runs.is_empty() || pairwise_mafs.is_empty() {
        let (_, hinge_contigs, _) = &indexed[hinge_idx];
        let direct_refs: Vec<horizon_aligner::DirectContigAlignment<'_>> = direct_runs
            .iter()
            .map(|(sp_i, xi, yi, aln)| {
                let (ref_contig, ref_seq) = &hinge_contigs[*xi];
                let (qry_sp, qry_contigs, _) = &indexed[*sp_i];
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
            &hinge_sp,
            &direct_refs,
            1000,
            &mut maf_writer,
            &mut gfa_writer,
        )?
    } else {
        horizon_aligner::merge_pairwise_mafs_prefixed_to_writers(
            &hinge_sp,
            &pairwise_mafs,
            1000,
            &mut maf_writer,
            &mut gfa_writer,
        )?
    };
    maf_writer.finish()?;
    gfa_writer.finish()?;
    let (hal_tree_desc, hal_bottom_segs) = if cli.write_hal {
        let hinge_contigs = std::mem::take(&mut indexed[hinge_idx].1);
        let hal_queries: Vec<(String, ContigList)> = hal_query_indices
            .into_iter()
            .map(|qi| {
                let sp = indexed[qi].0.clone();
                let contigs = std::mem::take(&mut indexed[qi].1);
                (sp, contigs)
            })
            .collect();
        let hal = if cli.progressive && hal_queries.len() >= 2 {
            let prog_blocks =
                extract_hal_blocks_from_pairwise_mafs(&hinge_contigs, &hal_queries, &pairwise_mafs);
            HalAlignment::from_multi_contig_progressive_blocks(
                &hinge_sp,
                &hinge_contigs,
                &hal_queries,
                &prog_blocks,
            )
        } else {
            HalAlignment::from_multi_contig_star_blocks(
                &hinge_sp,
                &hinge_contigs,
                &hal_queries,
                &hal_all_blocks,
            )
        };
        let desc = hal.tree_newick.clone();
        let nb = hal.genomes[0].bottom_segments.len();
        hal.write_hdf5_hal(std::path::Path::new(&format!("{out_prefix}.hal")))?;
        (desc, nb)
    } else {
        ("skipped".to_owned(), 0)
    };
    println!(
        "hinge={} species={} hal_tree={} hal_bottom_segs={} total_blocks={} multi_species_blocks={} elapsed_s={:.4}",
        hinge_sp,
        indexed.len(),
        hal_tree_desc,
        hal_bottom_segs,
        total_blocks,
        multi_blocks,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

#[allow(clippy::type_complexity)]
fn extract_hal_blocks_from_pairwise_mafs(
    hinge_contigs: &ContigList,
    hal_queries: &[(String, ContigList)],
    pairwise_mafs: &[(String, String)],
) -> Vec<(usize, usize, usize, usize, usize, usize, bool)> {
    let hinge_idx_by_name: HashMap<&str, usize> = hinge_contigs
        .iter()
        .enumerate()
        .map(|(idx, (name, _))| (name.as_str(), idx))
        .collect();
    let per_query: Vec<Vec<(usize, usize, usize, usize, usize, usize, bool)>> = hal_queries
        .par_iter()
        .zip(pairwise_mafs.par_iter())
        .enumerate()
        .map(|(q_idx, ((_, q_contigs), (_, maf_buf)))| {
            let mut local_out = Vec::new();
            let q_idx_by_name: HashMap<&str, usize> = q_contigs
                .iter()
                .enumerate()
                .map(|(idx, (name, _))| (name.as_str(), idx))
                .collect();
            for (h_row, q_row) in parse_two_row_maf_blocks(maf_buf) {
                let Some(&xi) = hinge_idx_by_name.get(h_row.src.as_str()) else {
                    continue;
                };
                let Some(&yi) = q_idx_by_name.get(q_row.src.as_str()) else {
                    continue;
                };
                let rev = h_row.strand != q_row.strand;
                let push_run = |out: &mut Vec<(usize, usize, usize, usize, usize, usize, bool)>,
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

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use crate::rescue::{
        build_ancestral_consensus_contigs, remap_consensus_maf_to_hinge, rescue_transitive_maf,
    };
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
        let hinge_contigs = vec![("hg38.chr1".to_string(), b"ACGTACGTACGTACGTACGT".to_vec())];
        let mid_contigs = vec![("chimp.chr1".to_string(), b"ACGTACGTACGTACGTACGT".to_vec())];
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
        let hinge_contigs = vec![("hg38.chr1".to_string(), b"ACGTACGTACGTACGTACGT".to_vec())];
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
            ("hg38.chr1".to_string(), b"ACGTACGTACGTACGTACGT".to_vec()),
            ("hg38.chr2".to_string(), b"TTTTTTTTTTTTTTTTTTTT".to_vec()),
        ];
        let mixed_cons =
            build_ancestral_consensus_contigs(&mixed_contigs, &[maf_s1, maf_s2], &[false, true]);
        assert_eq!(mixed_cons.len(), 2);
        assert_eq!(mixed_cons[0].1, [] as [u8; 0]);
        assert_eq!(mixed_cons[1].1, b"TTTTTTTTTTTTTTTTTTTT");
    }
}
