//! Constant-RAM streaming two-tier (`--hinges` / `--auto-hinges`) and post-hoc (`--extend-maf`)
//! tip-by-tip alignment and projection for `horizon_tba`.
//!
//! Instead of indexing all `N` FASTA assemblies into RAM simultaneously, only the `H` hinges are
//! indexed during the hinge backbone stage, and tips are streamed **one hinge + one tip at a time**
//! (`O(1)` `ContigIndex` RAM regardless of whether `N = 10` or `N = 1,400`), projecting each tip's
//! pairwise alignment onto the backbone and immediately dropping its index and pairwise blocks.

use crate::cli::TbaCliArgs;
use crate::pipeline::{
    IndexedSpecies, align_scheduled_species_pairs, align_species_pair_direct,
    apply_optional_tantan, compute_sketches, find_pairwise_maf_file, load_single_pair_from_dir,
    sketch_similarity_matrix,
};
use anyhow::{Context, Result};
use horizon_aligner::tba::{
    GuideTree, PairwiseTbaBlock, TbaBlock, TwoTierPlan, build_upgma_tree,
    run_tba_tree_merge_consuming, run_tba_tree_merge_lazy, similarity_submatrix,
    tba_project_tips_streaming,
};
use horizon_aligner::{
    AlignConfig, ProgressLogger, merge_contig_sketches, prepare_contig_index_owned,
    read_fasta_contigs, sketch_jaccard,
};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

const SKETCH_K: usize = 512;

fn index_single_species(
    species: &str,
    path: &std::path::Path,
    max_contigs: usize,
    respect_softmask: bool,
    use_tantan: bool,
) -> Result<IndexedSpecies> {
    let contigs = read_fasta_contigs(path, max_contigs)
        .with_context(|| format!("reading FASTA for `{species}` (`{}`)", path.display()))?;
    let skip_soft = respect_softmask || use_tantan;
    let idx = contigs
        .into_par_iter()
        .map(|(n, s)| {
            let s = apply_optional_tantan(s, use_tantan);
            Ok((n, prepare_contig_index_owned(s, skip_soft)?))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok((species.to_owned(), idx))
}

fn index_species_subset(
    species_list: &[String],
    fasta_map: &HashMap<&str, &PathBuf>,
    cli: &TbaCliArgs,
) -> Result<Vec<IndexedSpecies>> {
    species_list
        .iter()
        .map(|sp| {
            let path = fasta_map
                .get(sp.as_str())
                .with_context(|| format!("missing FASTA path for species `{sp}`"))?;
            index_single_species(
                sp,
                path,
                cli.max_contigs,
                cli.respect_softmask,
                cli.use_tantan,
            )
        })
        .collect()
}

/// Assign any unassigned tips by streaming each tip FASTA one at a time against precomputed hinge
/// `MinHash` sketches, keeping peak memory bounded to one tip genome.
fn assign_tips_streaming(
    cli: &TbaCliArgs,
    plan: &TwoTierPlan,
    fasta_map: &HashMap<&str, &PathBuf>,
    hinge_sketches: &[(String, Vec<u64>)],
    logger: ProgressLogger,
) -> Result<Vec<(String, String)>> {
    let mut final_tips: Vec<(String, String)> = plan.explicit_tips.clone();
    if plan.unassigned_tips.is_empty() {
        return Ok(final_tips);
    }
    anyhow::ensure!(
        !hinge_sketches.is_empty(),
        "cannot assign tips without hinge sketches"
    );
    logger.log(format_args!(
        "streaming MinHash assignment for {} unassigned tips across {} hinges",
        plan.unassigned_tips.len(),
        hinge_sketches.len()
    ));
    for tip in &plan.unassigned_tips {
        let path = fasta_map
            .get(tip.as_str())
            .with_context(|| format!("missing FASTA path for unassigned tip `{tip}`"))?;
        let tip_indexed = index_single_species(
            tip,
            path,
            cli.max_contigs,
            cli.respect_softmask,
            cli.use_tantan,
        )?;
        let tip_sketch = merge_contig_sketches(&tip_indexed.1, SKETCH_K);
        drop(tip_indexed);
        let best_hinge = hinge_sketches
            .iter()
            .max_by(|(ha, sa), (hb, sb)| {
                let ja = sketch_jaccard(&tip_sketch, sa);
                let jb = sketch_jaccard(&tip_sketch, sb);
                ja.total_cmp(&jb).then_with(|| hb.cmp(ha))
            })
            .map(|(h, _)| h.clone())
            .context("selecting best hinge for unassigned tip")?;
        final_tips.push((tip.clone(), best_hinge));
    }
    Ok(final_tips)
}

/// Execute two-tier (`--hinges`/`--auto-hinges`) or post-hoc (`--extend-maf`) alignment with
/// streaming tip-by-tip indexing and projection (`O(1)` tip index RAM).
#[allow(clippy::too_many_lines)]
pub(crate) fn run_two_tier_streaming(
    cli: &TbaCliArgs,
    plan: &TwoTierPlan,
    explicit_hinge_tree: Option<GuideTree>,
    backbone_opt: Option<(Vec<TbaBlock>, HashSet<String>)>,
    logger: ProgressLogger,
) -> Result<(Vec<TbaBlock>, String, usize, usize, usize)> {
    let is_extend = backbone_opt.is_some();
    let fasta_map: HashMap<&str, &PathBuf> = cli
        .fasta_specs
        .iter()
        .map(|(sp, p)| (sp.as_str(), p))
        .collect();

    if let Some(ref pdir) = cli.pairwise_dir {
        std::fs::create_dir_all(pdir)
            .with_context(|| format!("creating pairwise dir `{}`", pdir.display()))?;
    }

    let (backbone, tree_desc, total_species, final_tips) = if let Some((bb_blocks, bb_sps)) =
        backbone_opt
    {
        // Post-hoc --extend-maf: only sketch hinges if unassigned tips need hinge selection;
        // otherwise skip upfront hinge indexing completely and load one hinge at a time below!
        let final_tips = if plan.unassigned_tips.is_empty() {
            plan.explicit_tips.clone()
        } else {
            let mut hinge_sketches = Vec::with_capacity(plan.hinges.len());
            for h in &plan.hinges {
                let path = fasta_map
                    .get(h.as_str())
                    .with_context(|| format!("missing FASTA path for hinge `{h}`"))?;
                let hi = index_single_species(
                    h,
                    path,
                    cli.max_contigs,
                    cli.respect_softmask,
                    cli.use_tantan,
                )?;
                let sk = merge_contig_sketches(&hi.1, SKETCH_K);
                drop(hi);
                hinge_sketches.push((h.clone(), sk));
            }
            assign_tips_streaming(cli, plan, &fasta_map, &hinge_sketches, logger)?
        };
        let n_tot = bb_sps.len() + final_tips.len();
        (bb_blocks, "extend".to_owned(), n_tot, final_tips)
    } else {
        // Two-tier hinge+tips: index ONLY the H hinges, build hinge backbone, then free hinge
        // indexes before streaming tips one hinge at a time.
        let mut hinge_indexed = index_species_subset(&plan.hinges, &fasta_map, cli)?;
        let need_sketches = explicit_hinge_tree.is_none() || !plan.unassigned_tips.is_empty();
        let hinge_sketches_raw = if need_sketches {
            compute_sketches(&hinge_indexed)
        } else {
            Vec::new()
        };
        let hinge_tree = if let Some(t) = explicit_hinge_tree {
            t
        } else {
            let sim = sketch_similarity_matrix(&hinge_sketches_raw);
            let h_idx: HashMap<&str, usize> = plan
                .hinges
                .iter()
                .enumerate()
                .map(|(i, s)| (s.as_str(), i))
                .collect();
            let h_sim = similarity_submatrix(&sim, &plan.hinges, &h_idx);
            build_upgma_tree(&plan.hinges, &h_sim)?
        };
        let hinge_sketches: Vec<(String, Vec<u64>)> = if plan.unassigned_tips.is_empty() {
            Vec::new()
        } else {
            plan.hinges
                .iter()
                .cloned()
                .zip(hinge_sketches_raw)
                .collect()
        };

        let n_h = plan.hinges.len();
        let mut hinge_pairs: Vec<(usize, usize)> =
            Vec::with_capacity(n_h * n_h.saturating_sub(1) / 2);
        for i in 0..n_h {
            for j in (i + 1)..n_h {
                hinge_pairs.push((i, j));
            }
        }
        let desc = hinge_tree.to_newick();
        logger.log(format_args!(
            "streaming two-tier: aligning {n_h} hinges ({} hinge pairs) along {desc}",
            hinge_pairs.len()
        ));
        let mut pair_map =
            align_scheduled_species_pairs(cli, &mut hinge_indexed, &hinge_pairs, false, logger)?;
        drop(hinge_indexed);

        let bb = if let Some(ref pdir) = cli.pairwise_dir {
            run_tba_tree_merge_lazy(&hinge_tree, &|a, b| load_single_pair_from_dir(pdir, a, b))?
        } else {
            run_tba_tree_merge_consuming(&hinge_tree, &mut pair_map)
        };
        drop(pair_map);

        let final_tips = assign_tips_streaming(cli, plan, &fasta_map, &hinge_sketches, logger)?;
        (bb, desc, cli.species_names.len(), final_tips)
    };

    // Group tips by assigned hinge so we keep at most ONE hinge index + ONE tip index in RAM.
    let mut grouped_tips = final_tips;
    grouped_tips.sort_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(&b.0)));

    let total_tips = grouped_tips.len();
    let all_vs_all = cli.species_names.len() * cli.species_names.len().saturating_sub(1) / 2;
    let hinge_pair_count = if is_extend {
        0
    } else {
        plan.hinges.len() * plan.hinges.len().saturating_sub(1) / 2
    };
    logger.log(format_args!(
        "streaming {}: {} hinges + {} tips -> {} scheduled pairs (vs {all_vs_all} all-vs-all); projecting onto {} backbone blocks",
        if is_extend { "extend" } else { "hinge+tips" },
        plan.hinges.len(),
        total_tips,
        hinge_pair_count + total_tips,
        backbone.len()
    ));

    let cfg = AlignConfig {
        skip_softmasked: cli.respect_softmask || cli.use_tantan,
        ..AlignConfig::default()
    };
    let ext = if cli.pairwise_gz { "maf.gz" } else { "maf" };
    let mut cached_hinge: Option<IndexedSpecies> = None;
    let mut done_tips = 0usize;
    let report_step = (total_tips / 10).max(1);

    let merged = tba_project_tips_streaming(
        backbone,
        &grouped_tips,
        |tip, hinge| -> Result<Vec<PairwiseTbaBlock>> {
            let blks = if let Some(ref pdir) = cli.pairwise_dir
                && find_pairwise_maf_file(pdir, hinge, tip).is_some()
            {
                load_single_pair_from_dir(pdir, hinge, tip)?
            } else {
                if cached_hinge.as_ref().is_none_or(|(h, _)| h != hinge) {
                    cached_hinge = None;
                    let h_path = fasta_map
                        .get(hinge)
                        .with_context(|| format!("missing FASTA path for hinge `{hinge}`"))?;
                    cached_hinge = Some(index_single_species(
                        hinge,
                        h_path,
                        cli.max_contigs,
                        cli.respect_softmask,
                        cli.use_tantan,
                    )?);
                }
                let hinge_idx = cached_hinge
                    .as_ref()
                    .context("internal error: missing cached hinge index")?;
                let t_path = fasta_map
                    .get(tip)
                    .with_context(|| format!("missing FASTA path for tip `{tip}`"))?;
                let tip_idx = index_single_species(
                    tip,
                    t_path,
                    cli.max_contigs,
                    cli.respect_softmask,
                    cli.use_tantan,
                )?;
                let out_maf = cli
                    .pairwise_dir
                    .as_ref()
                    .map(|pdir| pdir.join(format!("{hinge}__{tip}.{ext}")));
                let pair_blks = align_species_pair_direct(
                    hinge_idx,
                    &tip_idx,
                    &cfg,
                    cli.top_k,
                    out_maf.as_deref(),
                    true,
                    true,
                )?;
                drop(tip_idx);
                pair_blks
            };

            done_tips += 1;
            if done_tips == total_tips || done_tips.is_multiple_of(report_step) {
                let pct = (done_tips * 100) / total_tips.max(1);
                logger.log(format_args!(
                    "streaming tips: {done_tips}/{total_tips} ({pct}%)"
                ));
            }
            Ok(blks)
        },
    )?;
    drop(cached_hinge);

    Ok((
        merged,
        tree_desc,
        total_species,
        plan.hinges.len(),
        total_tips,
    ))
}
