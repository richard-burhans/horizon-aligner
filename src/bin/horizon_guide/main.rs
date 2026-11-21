//! Nystrom/sketch-selected hinge + multi-species WGA pipeline (`horizon_wga` + `horizon_merge`).

use anyhow::{Context, Result, bail};
use horizon_aligner::{
    AlignConfig, ContigIndex, align_contig_pair_runs, format_maf_blocks_runs, hal::HalAlignment,
    merge_pairwise_mafs, prepare_contig_index_opts, read_fasta_contigs, select_top_k_pairs,
    shared_seed_count,
};
use rayon::prelude::*;
use rescue::{run_ancestral_consensus_rescue, run_progressive_sister_rescue};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

mod rescue;

pub(crate) type ContigList = Vec<(String, Vec<u8>)>;
pub(crate) type IndexedSpecies = (String, ContigList, Vec<(String, ContigIndex)>);

fn select_central_hinge(species: &[IndexedSpecies]) -> usize {
    (0..species.len())
        .into_par_iter()
        .map(|i| {
            let idx_i = &species[i].2;
            let total: usize = (0..species.len())
                .into_par_iter()
                .filter(|&j| j != i)
                .map(|j| {
                    let idx_j = &species[j].2;
                    idx_i
                        .par_iter()
                        .map(|(_, ix)| {
                            idx_j
                                .iter()
                                .map(|(_, iy)| shared_seed_count(ix, iy))
                                .sum::<usize>()
                        })
                        .sum::<usize>()
                })
                .sum();
            (total, std::cmp::Reverse(i))
        })
        .max()
        .map_or(0, |(_, std::cmp::Reverse(i))| i)
}

type HalBlock = (usize, usize, usize, usize, usize, bool);
type HalBlockList = Vec<HalBlock>;

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
            let hblks: Vec<HalBlock> = aln
                .aligned_runs
                .iter()
                .filter(|r| r.len > 0)
                .filter_map(|r| {
                    let q_start = if r.is_rc {
                        (r.y_start + 1).checked_sub(r.len)?
                    } else {
                        r.y_start
                    };
                    Some((xi, yi, r.x_start, q_start, r.len, r.is_rc))
                })
                .collect();
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
    for spec in &pos[2..] {
        let Some((sp, path)) = spec.split_once('=') else {
            bail!("expected <species>=<fasta>, got {spec}");
        };
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

    let indexed: Vec<IndexedSpecies> = cli
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

    let hinge_idx = select_central_hinge(&indexed);
    let (hinge_sp, hinge_contigs, hinge_indices) = &indexed[hinge_idx];
    let cfg = AlignConfig {
        skip_softmasked: cli.respect_softmask,
        ..AlignConfig::default()
    };

    let mut pairwise_mafs: Vec<(String, String)> = Vec::new();
    let mut hal_queries: Vec<(String, ContigList)> = Vec::new();
    let mut hal_all_blocks: Vec<(usize, usize, usize, usize, usize, usize, bool)> = Vec::new();
    let mut sp_to_pw_idx: HashMap<usize, usize> = HashMap::new();

    for (i, (qry_sp, qry_contigs, qry_indices)) in indexed.iter().enumerate() {
        if i == hinge_idx {
            continue;
        }
        let (maf_buf, hblocks) = align_species_pair(
            hinge_contigs,
            hinge_indices,
            qry_contigs,
            qry_indices,
            &cfg,
            cli.top_k,
        )?;
        if cli.write_hal {
            let q_idx = hal_queries.len();
            hal_queries.push((qry_sp.clone(), qry_contigs.clone()));
            for (xi, yi, rs, qs, len, rev) in hblocks {
                hal_all_blocks.push((q_idx, xi, yi, rs, qs, len, rev));
            }
        }
        sp_to_pw_idx.insert(i, pairwise_mafs.len());
        pairwise_mafs.push((qry_sp.clone(), maf_buf));
    }

    if cli.progressive && cli.write_hal {
        eprintln!(
            "[horizon_guide] note: --hal exports direct star-tree blocks; --progressive rescue is reflected in .maf/.gfa"
        );
    }

    if cli.progressive && indexed.len() >= 3 {
        run_progressive_sister_rescue(
            &indexed,
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
    }

    let (merged_maf, merged_gfa, total_blocks, multi_blocks) =
        merge_pairwise_mafs(hinge_sp, &pairwise_mafs, 1000);
    let merged_maf = prefix_hinge_species_in_maf(hinge_sp, &merged_maf);
    let out_prefix = &cli.out_prefix;
    fs::write(format!("{out_prefix}.maf"), merged_maf)
        .with_context(|| format!("writing {out_prefix}.maf"))?;
    fs::write(format!("{out_prefix}.gfa"), merged_gfa)
        .with_context(|| format!("writing {out_prefix}.gfa"))?;
    let (hal_tree_desc, hal_bottom_segs) = if cli.write_hal {
        let hal = HalAlignment::from_multi_contig_star_blocks(
            hinge_sp,
            hinge_contigs,
            &hal_queries,
            &hal_all_blocks,
        );
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

fn prefix_hinge_species_in_maf(hinge_sp: &str, maf: &str) -> String {
    let prefix = format!("{hinge_sp}.");
    let mut out = String::with_capacity(maf.len() + 256);
    let mut first_s_in_block = true;
    for line in maf.lines() {
        if line.is_empty() || line.starts_with("a ") || line == "a" {
            first_s_in_block = true;
            out.push_str(line);
            out.push('\n');
        } else if first_s_in_block && let Some(rest) = line.strip_prefix("s ") {
            first_s_in_block = false;
            if let Some((src, tail)) = rest.split_once(' ')
                && !src.starts_with(&prefix)
            {
                out.push_str("s ");
                out.push_str(&prefix);
                out.push_str(src);
                out.push(' ');
                out.push_str(tail);
                out.push('\n');
            } else {
                out.push_str(line);
                out.push('\n');
            }
        } else {
            if line.starts_with("s ") {
                first_s_in_block = false;
            }
            out.push_str(line);
            out.push('\n');
        }
    }
    out
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::prefix_hinge_species_in_maf;
    use crate::rescue::{
        build_ancestral_consensus_contigs, remap_consensus_maf_to_hinge, rescue_transitive_maf,
    };
    use horizon_aligner::complement_base;

    #[test]
    fn prefix_hinge_species_adds_missing_species_prefix() {
        let maf = "##maf version=1\na score=2.0\ns NBAG03001603.1 0 4 + 10 ACGT\ns hg38.chr22 0 4 + 10 ACGT\n\n";
        let got = prefix_hinge_species_in_maf("Pan_troglodytes", maf);
        assert!(got.contains("s Pan_troglodytes.NBAG03001603.1 0 4 + 10 ACGT\n"));
        assert!(got.contains("s hg38.chr22 0 4 + 10 ACGT\n"));
    }

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
