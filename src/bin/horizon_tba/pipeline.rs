use crate::cli::TbaCliArgs;
use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use horizon_aligner::io_util::{AtomicTempFile, open_text_reader, open_text_writer};
use horizon_aligner::tba::{
    GuideTree, PairwiseTbaBlock, ProgressiveAncestralConfig, TbaBlock, TwoTierPlan,
    aligned_bp_similarity_matrix, build_upgma_tree, canonical_pair_key,
    format_ancestors_summary_tsv, parse_multi_tba_reader, parse_newick, parse_pairwise_tba_reader,
    run_progressive_ancestral_tba, run_tba_tree_merge, run_tba_tree_merge_consuming,
    run_tba_tree_merge_lazy, similarity_submatrix, tba_project_tips_onto_blockset,
};
use horizon_aligner::{
    AlignConfig, ContigIndex, ProgressLogger, TantanConfig, adaptive_top_k_for_fragmentation,
    align_contig_pair_runs, format_maf_blocks_runs, merge_contig_sketches,
    prepare_contig_index_owned, read_fasta_contigs, select_top_k_pairs, shared_seed_count,
    sketch_jaccard, tantan_mask,
};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::{Read as _, Write as _};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

const MAX_NEWICK_FILE_BYTES: u64 = 16 * 1024 * 1024;
pub(crate) type IndexedSpecies = (String, Vec<(String, ContigIndex)>);

pub(crate) fn load_explicit_tree(tree_arg: &str, species_names: &[String]) -> Result<GuideTree> {
    let text = if Path::new(tree_arg).is_file() {
        let file = std::fs::File::open(tree_arg)
            .with_context(|| format!("opening Newick tree file `{tree_arg}`"))?;
        let mut buf = String::new();
        let n = file
            .take(MAX_NEWICK_FILE_BYTES + 1)
            .read_to_string(&mut buf)
            .with_context(|| format!("reading Newick tree file `{tree_arg}`"))?;
        anyhow::ensure!(
            u64::try_from(n).unwrap_or(u64::MAX) <= MAX_NEWICK_FILE_BYTES,
            "Newick tree file `{tree_arg}` exceeds {MAX_NEWICK_FILE_BYTES} bytes"
        );
        buf
    } else {
        tree_arg.to_owned()
    };
    let tree = parse_newick(&text)?;
    let tree_set: HashSet<&str> = tree.leaves().into_iter().collect();
    let cli_set: HashSet<&str> = species_names.iter().map(String::as_str).collect();
    anyhow::ensure!(
        tree_set == cli_set,
        "Newick tree leaves ({tree_set:?}) do not match input species ({cli_set:?})"
    );
    Ok(tree)
}

pub(crate) fn find_pairwise_maf_file(
    dir: &Path,
    sp_a: &str,
    sp_b: &str,
) -> Option<(PathBuf, bool)> {
    const SUFFIXES: &[(&str, &str)] = &[
        ("__", ".maf"),
        ("__", ".maf.gz"),
        ("__", ".maf.zst"),
        (".", ".sing.maf"),
        (".", ".sing.maf.gz"),
        (".", ".maf"),
    ];
    for (first, second, swapped) in [(sp_a, sp_b, false), (sp_b, sp_a, true)] {
        for &(sep, ext) in SUFFIXES {
            let p = dir.join(format!("{first}{sep}{second}{ext}"));
            if p.is_file() && std::fs::metadata(&p).is_ok_and(|m| m.len() > 0) {
                return Some((p, swapped));
            }
        }
    }
    None
}

pub(crate) fn compute_sketches(species: &[IndexedSpecies]) -> Vec<Vec<u64>> {
    const K: usize = 512;
    species
        .par_iter()
        .map(|(_, idx_list)| merge_contig_sketches(idx_list, K))
        .collect()
}

pub(crate) fn sketch_similarity_matrix(sketches: &[Vec<u64>]) -> Vec<Vec<f64>> {
    let n = sketches.len();
    let mut sim = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        sim[i][i] = 1.0;
        for j in (i + 1)..n {
            let val = sketch_jaccard(&sketches[i], &sketches[j]);
            sim[i][j] = val;
            sim[j][i] = val;
        }
    }
    sim
}

pub(crate) fn align_species_pair_direct(
    sp_a: &IndexedSpecies,
    sp_b: &IndexedSpecies,
    cfg: &AlignConfig,
    top_k: usize,
    emit_maf_path: Option<&Path>,
    parallel_contigs: bool,
    return_blocks: bool,
) -> Result<Vec<PairwiseTbaBlock>> {
    let min_shared_seeds = 8usize;
    let eff_top_k = adaptive_top_k_for_fragmentation(top_k, &sp_a.1, &sp_b.1);
    let top_mask =
        (eff_top_k > 0).then(|| select_top_k_pairs(&sp_a.1, &sp_b.1, eff_top_k, min_shared_seeds));
    let mut contig_pairs = Vec::new();
    for xi in 0..sp_a.1.len() {
        for yi in 0..sp_b.1.len() {
            let keep = if let Some(mask) = &top_mask {
                mask[xi][yi]
            } else {
                shared_seed_count(&sp_a.1[xi].1, &sp_b.1[yi].1) >= min_shared_seeds
            };
            if keep {
                contig_pairs.push((xi, yi));
            }
        }
    }
    let want_maf = emit_maf_path.is_some();
    let run_one = |&(xi, yi): &(usize, usize)| -> Result<(Vec<PairwiseTbaBlock>, Option<String>)> {
        let (nx, ix) = &sp_a.1[xi];
        let (ny, iy) = &sp_b.1[yi];
        let aln = align_contig_pair_runs(ix, iy, cfg)?;
        if aln.aligned_runs.is_empty() {
            return Ok((Vec::new(), None));
        }
        let sx = ix.seq();
        let sy = iy.seq();
        let blks = if return_blocks {
            PairwiseTbaBlock::from_run_alignment(&sp_a.0, nx, sx, &sp_b.0, ny, sy, &aln)
        } else {
            Vec::new()
        };
        let maf_body = want_maf.then(|| {
            let chunk = format_maf_blocks_runs(nx, sx.len(), ny, sy.len(), sx, sy, &aln);
            chunk
                .strip_prefix("##maf version=1 scoring=horizon_mdl\n")
                .unwrap_or(&chunk)
                .to_owned()
        });
        Ok((blks, maf_body))
    };

    let chunks: Vec<(Vec<PairwiseTbaBlock>, Option<String>)> = if parallel_contigs {
        contig_pairs
            .par_iter()
            .map(run_one)
            .collect::<Result<Vec<_>>>()?
    } else {
        contig_pairs
            .iter()
            .map(run_one)
            .collect::<Result<Vec<_>>>()?
    };

    let mut all_blocks = Vec::new();
    if let Some(out_p) = emit_maf_path {
        let use_gz = out_p
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("gz"));
        let (tmp_guard, buf_w) = AtomicTempFile::create_exclusive(out_p)?;
        if use_gz {
            let mut gz = GzEncoder::new(buf_w, Compression::fast());
            gz.write_all(b"##maf version=1 scoring=horizon_mdl\n")?;
            for (blks, body_opt) in chunks {
                all_blocks.extend(blks);
                if let Some(body) = body_opt {
                    gz.write_all(body.as_bytes())?;
                }
            }
            let mut buf_w = gz.finish()?;
            buf_w.flush()?;
            drop(buf_w);
        } else {
            let mut w = buf_w;
            w.write_all(b"##maf version=1 scoring=horizon_mdl\n")?;
            for (blks, body_opt) in chunks {
                all_blocks.extend(blks);
                if let Some(body) = body_opt {
                    w.write_all(body.as_bytes())?;
                }
            }
            w.flush()?;
            drop(w);
        }
        tmp_guard.persist(out_p)?;
    } else {
        for (blks, _) in chunks {
            all_blocks.extend(blks);
        }
    }
    Ok(all_blocks)
}

pub(crate) fn load_single_pair_from_dir(
    maf_dir: &Path,
    sp_a: &str,
    sp_b: &str,
) -> Result<Vec<PairwiseTbaBlock>> {
    let Some((path, swapped)) = find_pairwise_maf_file(maf_dir, sp_a, sp_b) else {
        return Ok(Vec::new());
    };
    let reader = open_text_reader(&path)?;
    let (first, second) = if swapped { (sp_b, sp_a) } else { (sp_a, sp_b) };
    parse_pairwise_tba_reader(first, second, reader)
}

type PairBlockMap = HashMap<(String, String), Vec<PairwiseTbaBlock>>;

fn load_pair_subset_from_dir(
    cli: &TbaCliArgs,
    maf_dir: &Path,
    sp_pairs: &[(usize, usize)],
    logger: ProgressLogger,
) -> Result<PairBlockMap> {
    anyhow::ensure!(
        maf_dir.is_dir(),
        "--from-mafs directory does not exist: `{}`",
        maf_dir.display()
    );
    logger.log(format_args!(
        "loading pairwise MAFs from `{}` for {} species ({} pairs)",
        maf_dir.display(),
        cli.species_names.len(),
        sp_pairs.len()
    ));
    let loaded: Vec<((String, String), Vec<PairwiseTbaBlock>)> = sp_pairs
        .par_iter()
        .filter_map(|&(i, j)| {
            let sp_a = &cli.species_names[i];
            let sp_b = &cli.species_names[j];
            let (path, swapped) = find_pairwise_maf_file(maf_dir, sp_a, sp_b)?;
            let res = (|| -> Result<_> {
                let reader = open_text_reader(&path)?;
                let (first, second) = if swapped { (sp_b, sp_a) } else { (sp_a, sp_b) };
                let blks = parse_pairwise_tba_reader(first, second, reader)?;
                Ok((canonical_pair_key(sp_a, sp_b), blks))
            })();
            Some(res)
        })
        .collect::<Result<Vec<_>>>()?;
    let pair_map: PairBlockMap = loaded.into_iter().collect();
    anyhow::ensure!(
        !pair_map.is_empty(),
        "no matching pairwise MAF files found in `{}` for species {:?}",
        maf_dir.display(),
        cli.species_names
    );
    let total_in_blocks: usize = pair_map.values().map(Vec::len).sum();
    logger.log(format_args!(
        "loaded {}/{} pairwise MAF files ({total_in_blocks} blocks)",
        pair_map.len(),
        sp_pairs.len()
    ));
    Ok(pair_map)
}

fn load_pairwise_from_dir(
    cli: &TbaCliArgs,
    explicit_tree: Option<GuideTree>,
    maf_dir: &Path,
    sp_pairs: &[(usize, usize)],
    logger: ProgressLogger,
) -> Result<(GuideTree, PairBlockMap)> {
    let pair_map = load_pair_subset_from_dir(cli, maf_dir, sp_pairs, logger)?;
    let tree = if let Some(t) = explicit_tree {
        t
    } else {
        let sim = aligned_bp_similarity_matrix(&cli.species_names, &pair_map);
        build_upgma_tree(&cli.species_names, &sim)?
    };
    Ok((tree, pair_map))
}

pub(crate) fn apply_optional_tantan(mut seq: Vec<u8>, enable: bool) -> Vec<u8> {
    if enable {
        let mask = tantan_mask(&seq, &TantanConfig::default());
        for (b, m) in seq.iter_mut().zip(mask) {
            if m {
                b.make_ascii_lowercase();
            }
        }
    }
    seq
}

fn index_fasta_specs(cli: &TbaCliArgs, logger: ProgressLogger) -> Result<Vec<IndexedSpecies>> {
    let skip_soft = cli.respect_softmask || cli.use_tantan;
    let indexed: Vec<IndexedSpecies> = cli
        .fasta_specs
        .par_iter()
        .map(|(sp, path)| {
            let contigs = read_fasta_contigs(path, cli.max_contigs)?;
            let idx = contigs
                .into_par_iter()
                .map(|(n, s)| {
                    let s = apply_optional_tantan(s, cli.use_tantan);
                    Ok((n, prepare_contig_index_owned(s, skip_soft)?))
                })
                .collect::<Result<Vec<_>>>()?;
            Ok((sp.clone(), idx))
        })
        .collect::<Result<Vec<_>>>()?;

    let total_contigs: usize = indexed.iter().map(|(_, idx)| idx.len()).sum();
    let total_bp: usize = indexed
        .iter()
        .flat_map(|(_, idx)| idx.iter())
        .map(|(_, ix)| ix.seq().len())
        .sum();
    #[allow(clippy::cast_precision_loss)]
    let total_mb = (total_bp as f64) / 1_000_000.0;
    logger.log(format_args!(
        "indexed {} species ({total_contigs} contigs, {total_mb:.2} Mb)",
        indexed.len()
    ));
    Ok(indexed)
}

pub(crate) fn align_scheduled_species_pairs(
    cli: &TbaCliArgs,
    indexed: &mut Vec<IndexedSpecies>,
    sp_pairs: &[(usize, usize)],
    reload_from_dir: bool,
    logger: ProgressLogger,
) -> Result<PairBlockMap> {
    if let Some(ref pdir) = cli.pairwise_dir {
        std::fs::create_dir_all(pdir)
            .with_context(|| format!("creating pairwise dir `{}`", pdir.display()))?;
    }
    let cfg = AlignConfig {
        skip_softmasked: cli.respect_softmask || cli.use_tantan,
        ..AlignConfig::default()
    };
    let spill_to_dir = cli.pairwise_dir.is_some();
    let ext = if cli.pairwise_gz { "maf.gz" } else { "maf" };
    // Avoid nested Rayon `par_iter` working-set amplification on multi-contig genomes:
    // parallelize inside each species pair when genomes have >=3 contigs on average, or
    // across species pairs for single-chromosome panels.
    let multi_contig =
        indexed.iter().map(|(_, c)| c.len()).sum::<usize>() >= indexed.len().saturating_mul(3);
    let total_pairs = sp_pairs.len();
    let report_step = (total_pairs / 10).max(1);
    let done_pairs = AtomicUsize::new(0);

    let run_pair = |&(i, j): &(usize, usize)| -> Result<((String, String), Vec<PairwiseTbaBlock>)> {
        // Canonicalize target/query orientation alphabetically so permuting species input order
        // yields bit-identical pairwise alignments, MAFs, GFAs, and VCFs (Invariant I3).
        let (lo_i, hi_i) = if indexed[i].0 <= indexed[j].0 {
            (i, j)
        } else {
            (j, i)
        };
        let sp_a = &indexed[lo_i].0;
        let sp_b = &indexed[hi_i].0;
        let blks = if let Some(ref pdir) = cli.pairwise_dir
            && let Some((existing, swapped)) = find_pairwise_maf_file(pdir, sp_a, sp_b)
        {
            if spill_to_dir {
                Vec::new()
            } else {
                let reader = open_text_reader(&existing)?;
                let (first, second) = if swapped { (sp_b, sp_a) } else { (sp_a, sp_b) };
                parse_pairwise_tba_reader(first, second, reader)?
            }
        } else {
            let out_maf = cli
                .pairwise_dir
                .as_ref()
                .map(|pdir| pdir.join(format!("{sp_a}__{sp_b}.{ext}")));
            align_species_pair_direct(
                &indexed[lo_i],
                &indexed[hi_i],
                &cfg,
                cli.top_k,
                out_maf.as_deref(),
                multi_contig,
                !spill_to_dir,
            )?
        };
        let done = done_pairs.fetch_add(1, Ordering::Relaxed) + 1;
        if done == total_pairs || done.is_multiple_of(report_step) {
            let pct = (done * 100) / total_pairs.max(1);
            logger.log(format_args!(
                "pairwise: {done}/{total_pairs} species pairs ({pct}%)"
            ));
        }
        Ok((canonical_pair_key(sp_a, sp_b), blks))
    };

    let computed: Vec<((String, String), Vec<PairwiseTbaBlock>)> = if multi_contig {
        sp_pairs.iter().map(run_pair).collect::<Result<Vec<_>>>()?
    } else {
        sp_pairs
            .par_iter()
            .map(run_pair)
            .collect::<Result<Vec<_>>>()?
    };

    // Free genome indexes before Phase 2 merge so ContigIndex and TBA merge buffers never overlap.
    indexed.clear();
    indexed.shrink_to_fit();

    if let Some(ref pdir) = cli.pairwise_dir {
        if reload_from_dir {
            load_pair_subset_from_dir(cli, pdir, sp_pairs, logger)
        } else {
            Ok(HashMap::new())
        }
    } else {
        Ok(computed.into_iter().collect())
    }
}

fn load_backbone_for_extend(
    cli: &TbaCliArgs,
    logger: ProgressLogger,
) -> Result<Option<(Vec<TbaBlock>, HashSet<String>)>> {
    let Some(ref bb_path) = cli.extend_maf_path else {
        return Ok(None);
    };
    anyhow::ensure!(
        cli.tree_arg.is_none(),
        "--tree is not used with --extend-maf (backbone topology is already merged)"
    );
    let reader = open_text_reader(bb_path)
        .with_context(|| format!("opening --extend-maf `{}`", bb_path.display()))?;
    let blks = parse_multi_tba_reader(reader)?;
    anyhow::ensure!(
        !blks.is_empty(),
        "--extend-maf `{}` contains no multi-species MAF blocks",
        bb_path.display()
    );
    let sps: HashSet<String> = blks
        .iter()
        .flat_map(|b| b.rows.iter().map(|r| r.species.clone()))
        .collect();
    logger.log(format_args!(
        "loaded backbone `{}` ({} blocks, {} species)",
        bb_path.display(),
        blks.len(),
        sps.len()
    ));
    Ok(Some((blks, sps)))
}

struct TwoTierInputs {
    hinge_tree: Option<GuideTree>,
    final_tips: Vec<(String, String)>,
    pair_blocks: PairBlockMap,
}

fn resolve_two_tier_inputs(
    cli: &TbaCliArgs,
    plan: &TwoTierPlan,
    explicit_hinge_tree: Option<GuideTree>,
    is_extend: bool,
    sp_to_idx: &HashMap<&str, usize>,
    logger: ProgressLogger,
) -> Result<TwoTierInputs> {
    let (sim, pair_blocks_opt, mut indexed_opt) = if let Some(ref maf_dir) = cli.from_mafs_dir {
        let mut candidate_pairs = plan.schedule_pairs(&plan.explicit_tips, sp_to_idx, !is_extend);
        for tip in &plan.unassigned_tips {
            let it = sp_to_idx[tip.as_str()];
            for h in &plan.hinges {
                let ih = sp_to_idx[h.as_str()];
                candidate_pairs.push((it.min(ih), it.max(ih)));
            }
        }
        candidate_pairs.sort_unstable();
        candidate_pairs.dedup();
        let pair_map = load_pair_subset_from_dir(cli, maf_dir, &candidate_pairs, logger)?;
        let sim = aligned_bp_similarity_matrix(&cli.species_names, &pair_map);
        (sim, Some(pair_map), None)
    } else {
        let indexed = index_fasta_specs(cli, logger)?;
        let need_sim =
            !plan.unassigned_tips.is_empty() || (!is_extend && explicit_hinge_tree.is_none());
        let sim = if need_sim {
            sketch_similarity_matrix(&compute_sketches(&indexed))
        } else {
            vec![vec![0.0; cli.species_names.len()]; cli.species_names.len()]
        };
        (sim, None, Some(indexed))
    };

    let assigned = plan.finalize_assignments(&sim, sp_to_idx)?;
    let final_tips: Vec<(String, String)> = assigned.into_iter().map(|(t, h, _)| (t, h)).collect();
    let hinge_tree = if is_extend {
        None
    } else if let Some(t) = explicit_hinge_tree {
        Some(t)
    } else {
        let h_sim = similarity_submatrix(&sim, &plan.hinges, sp_to_idx);
        Some(build_upgma_tree(&plan.hinges, &h_sim)?)
    };

    let pair_blocks = if let Some(pm) = pair_blocks_opt {
        pm
    } else {
        let mut indexed = indexed_opt.take().unwrap_or_default();
        let sp_pairs = plan.schedule_pairs(&final_tips, sp_to_idx, !is_extend);
        let all_vs_all = cli.species_names.len() * cli.species_names.len().saturating_sub(1) / 2;
        logger.log(format_args!(
            "two-tier {}: {} hinges + {} tips -> scheduling {} pairs (vs {all_vs_all} all-vs-all)",
            if is_extend { "extend" } else { "hinge+tips" },
            plan.hinges.len(),
            final_tips.len(),
            sp_pairs.len()
        ));
        align_scheduled_species_pairs(cli, &mut indexed, &sp_pairs, true, logger)?
    };
    Ok(TwoTierInputs {
        hinge_tree,
        final_tips,
        pair_blocks,
    })
}

pub(crate) fn run_two_tier_or_extend(
    cli: &TbaCliArgs,
    hinges_opt: Option<&[String]>,
    tip_map_opt: Option<&[(String, String)]>,
    logger: ProgressLogger,
) -> Result<(Vec<TbaBlock>, String, usize, usize, usize)> {
    let backbone_opt = load_backbone_for_extend(cli, logger)?;
    let bb_sp_ref = backbone_opt.as_ref().map(|(_, s)| s);
    let plan = TwoTierPlan::resolve(&cli.species_names, hinges_opt, tip_map_opt, bb_sp_ref)?;
    let is_extend = backbone_opt.is_some();
    let explicit_hinge_tree = cli
        .tree_arg
        .as_deref()
        .map(|arg| load_explicit_tree(arg, &plan.hinges))
        .transpose()?;
    if cli.from_mafs_dir.is_none() {
        return crate::stream_extend::run_two_tier_streaming(
            cli,
            &plan,
            explicit_hinge_tree,
            backbone_opt,
            logger,
        );
    }
    let sp_to_idx: HashMap<&str, usize> = cli
        .species_names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let inputs = resolve_two_tier_inputs(
        cli,
        &plan,
        explicit_hinge_tree,
        is_extend,
        &sp_to_idx,
        logger,
    )?;

    let (backbone, tree_desc, total_species) = if let Some((bb_blocks, bb_sps)) = backbone_opt {
        (
            bb_blocks,
            "extend".to_owned(),
            bb_sps.len() + inputs.final_tips.len(),
        )
    } else {
        let hinge_tree = inputs
            .hinge_tree
            .ok_or_else(|| anyhow::anyhow!("missing hinge guide tree in two-tier mode"))?;
        let desc = hinge_tree.to_newick();
        logger.log(format_args!(
            "merging {}-hinge TBA backbone along {desc}",
            plan.hinges.len()
        ));
        (
            run_tba_tree_merge(&hinge_tree, &inputs.pair_blocks),
            desc,
            cli.species_names.len(),
        )
    };

    logger.log(format_args!(
        "projecting {} tips onto {} backbone blocks",
        inputs.final_tips.len(),
        backbone.len()
    ));
    let merged = tba_project_tips_onto_blockset(backbone, &inputs.final_tips, &inputs.pair_blocks)?;
    Ok((
        merged,
        tree_desc,
        total_species,
        plan.hinges.len(),
        inputs.final_tips.len(),
    ))
}

pub(crate) fn run_all_vs_all_tba(
    cli: &TbaCliArgs,
    logger: ProgressLogger,
) -> Result<(Vec<TbaBlock>, String, usize)> {
    // Validate any explicit --tree upfront before reading or indexing FASTA/MAF inputs.
    let explicit_tree = cli
        .tree_arg
        .as_deref()
        .map(|arg| load_explicit_tree(arg, &cli.species_names))
        .transpose()?;
    let n_sp = cli.species_names.len();
    let mut sp_pairs: Vec<(usize, usize)> = Vec::with_capacity(n_sp * n_sp.saturating_sub(1) / 2);
    for i in 0..n_sp {
        for j in (i + 1)..n_sp {
            sp_pairs.push((i, j));
        }
    }
    let (guide_tree, lazy_dir, mut pair_blocks) = if let Some(maf_dir) = &cli.from_mafs_dir {
        if let Some(tree) = explicit_tree {
            anyhow::ensure!(
                maf_dir.is_dir(),
                "--from-mafs directory does not exist: `{}`",
                maf_dir.display()
            );
            let found_any = sp_pairs.iter().any(|&(i, j)| {
                find_pairwise_maf_file(maf_dir, &cli.species_names[i], &cli.species_names[j])
                    .is_some()
            });
            anyhow::ensure!(
                found_any,
                "no matching pairwise MAF files found in `{}` for species {:?}",
                maf_dir.display(),
                cli.species_names
            );
            (tree, Some(maf_dir.clone()), HashMap::new())
        } else {
            let (tree, map) = load_pairwise_from_dir(cli, None, maf_dir, &sp_pairs, logger)?;
            (tree, None, map)
        }
    } else {
        let mut indexed = index_fasta_specs(cli, logger)?;
        let tree = if let Some(t) = explicit_tree {
            t
        } else {
            let sketches = compute_sketches(&indexed);
            let sim = sketch_similarity_matrix(&sketches);
            build_upgma_tree(&cli.species_names, &sim)?
        };
        logger.log(format_args!(
            "guide tree={}, scheduling {} species pairs",
            tree.to_newick(),
            sp_pairs.len()
        ));
        let map = align_scheduled_species_pairs(cli, &mut indexed, &sp_pairs, false, logger)?;
        (tree, cli.pairwise_dir.clone(), map)
    };
    let desc = guide_tree.to_newick();
    let blks = if let Some(ref pdir) = lazy_dir {
        logger.log(format_args!(
            "merging reference-free TBA blocksets lazily from `{}` along {desc}",
            pdir.display()
        ));
        run_tba_tree_merge_lazy(&guide_tree, &|sp_a, sp_b| {
            load_single_pair_from_dir(pdir, sp_a, sp_b)
        })?
    } else {
        logger.log(format_args!(
            "merging reference-free TBA blocksets along {desc}"
        ));
        run_tba_tree_merge_consuming(&guide_tree, &mut pair_blocks)
    };
    Ok((blks, desc, n_sp))
}

pub(crate) fn run_progressive_ancestral_mode(
    cli: &TbaCliArgs,
    logger: ProgressLogger,
) -> Result<(Vec<TbaBlock>, String, String, usize)> {
    anyhow::ensure!(
        cli.from_mafs_dir.is_none(),
        "--progressive-ancestors requires FASTA inputs (cannot be combined with --from-mafs)"
    );
    let explicit_tree = cli
        .tree_arg
        .as_deref()
        .map(|arg| load_explicit_tree(arg, &cli.species_names))
        .transpose()?;

    let mut leaf_genomes: HashMap<String, Vec<(String, Vec<u8>)>> =
        HashMap::with_capacity(cli.fasta_specs.len());
    for (sp, path) in &cli.fasta_specs {
        let ctgs: Vec<(String, Vec<u8>)> = read_fasta_contigs(path, cli.max_contigs)?
            .into_iter()
            .map(|(n, s)| (n, apply_optional_tantan(s, cli.use_tantan)))
            .collect();
        leaf_genomes.insert(sp.clone(), ctgs);
    }

    let guide_tree = if let Some(t) = explicit_tree {
        t
    } else {
        let indexed = index_fasta_specs(cli, logger)?;
        let sketches = compute_sketches(&indexed);
        let sim = sketch_similarity_matrix(&sketches);
        build_upgma_tree(&cli.species_names, &sim)?
    };
    let desc = guide_tree.to_newick();
    logger.log(format_args!(
        "progressive ancestral alignment (Figure 1D): {} leaves, max_outgroups={}, tree={desc}",
        cli.species_names.len(),
        cli.max_outgroups
    ));

    let cfg = ProgressiveAncestralConfig {
        max_outgroups: cli.max_outgroups,
        top_k: cli.top_k,
        respect_softmask: cli.respect_softmask || cli.use_tantan,
        use_transducer_ml: cli.ancestral_ml,
    };
    let res = run_progressive_ancestral_tba(&guide_tree, &leaf_genomes, &cfg)?;

    if let Some(ref anc_dir) = cli.emit_ancestors_dir {
        std::fs::create_dir_all(anc_dir)
            .with_context(|| format!("creating --emit-ancestors dir `{}`", anc_dir.display()))?;
        for (anc_name, ctgs, _) in &res.ancestors {
            let fa_path = anc_dir.join(format!("{anc_name}.fa"));
            let mut w = open_text_writer(Some(&fa_path.to_string_lossy()))?;
            for (ctg_name, seq) in ctgs {
                writeln!(w, ">{ctg_name}")?;
                for chunk in seq.chunks(80) {
                    w.write_all(chunk)?;
                    w.write_all(b"\n")?;
                }
            }
            w.finish()?;
        }
        let tsv = format_ancestors_summary_tsv(&res.ancestors);
        std::fs::write(anc_dir.join("ancestors_summary.tsv"), tsv)?;
        logger.log(format_args!(
            "wrote {} reconstructed ancestral genomes into `{}`",
            res.ancestors.len(),
            anc_dir.display()
        ));
    }

    if cli.write_hal {
        let hal_path = PathBuf::from(format!("{}.hal", cli.out_prefix));
        res.hal
            .write_hdf5_hal(&hal_path)
            .with_context(|| format!("writing progressive HAL `{}`", hal_path.display()))?;
        logger.log(format_args!(
            "wrote progressive HAL `{}` ({} genomes)",
            hal_path.display(),
            res.hal.genomes.len()
        ));
    }

    let mode = format!(
        "tba_progressive_ancestors ancestors={} sub_pairs={}",
        res.ancestors.len(),
        res.subproblem_pairs_aligned
    );
    Ok((res.merged_blocks, mode, desc, cli.species_names.len()))
}
