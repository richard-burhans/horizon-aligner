//! Reference-free all-vs-all Threaded Blockset Aligner (`horizon_tba`).

mod cli;
mod pipeline;
mod stream_extend;

use anyhow::{Context, Result};
use cli::parse_args;
use horizon_aligner::io_util::{open_text_reader, open_text_writer};
use horizon_aligner::tba::{
    SourceGenomeContigs, parse_hinge_list, parse_tip_map_reader, write_tba_maf_and_gfa,
    write_tba_maf_and_variation_gfa_with_seqs,
};
use horizon_aligner::{
    EisensteinHasher, EisensteinSamplerSketch, ProgressLogger, init_rayon_threads,
    read_fasta_contigs, select_auto_hinges_eisenstein,
};
#[cfg(test)]
use pipeline::{find_pairwise_maf_file, load_explicit_tree};
use pipeline::{run_all_vs_all_tba, run_progressive_ancestral_mode, run_two_tier_or_extend};
use rayon::prelude::*;
use std::time::Instant;

#[allow(clippy::too_many_lines)]
fn main() -> Result<()> {
    let Some(cli) = parse_args()? else {
        return Ok(());
    };
    let mut hinges_opt = cli
        .hinges_arg
        .as_deref()
        .map(parse_hinge_list)
        .transpose()?;
    let mut tip_map_opt = cli
        .tip_map_path
        .as_deref()
        .map(|p| {
            let r = open_text_reader(p)
                .with_context(|| format!("opening --tip-map file `{}`", p.display()))?;
            parse_tip_map_reader(r)
        })
        .transpose()?;

    init_rayon_threads(cli.threads)?;
    let t0 = Instant::now();
    let logger = ProgressLogger::new("horizon_tba", t0, cli.quiet);

    if let Some(num_h) = cli.auto_hinges {
        anyhow::ensure!(
            hinges_opt.is_none() && tip_map_opt.is_none() && cli.extend_maf_path.is_none(),
            "--auto-hinges cannot be combined with --hinges, --tip-map, or --extend-maf"
        );
        anyhow::ensure!(
            cli.from_mafs_dir.is_none(),
            "--auto-hinges requires FASTA inputs (cannot be combined with --from-mafs)"
        );
        let hasher = EisensteinHasher::default();
        let sketches: Vec<(String, EisensteinSamplerSketch)> = cli
            .fasta_specs
            .par_iter()
            .map(|(sp, path)| {
                let ctgs = read_fasta_contigs(path, cli.max_contigs)?;
                let slices: Vec<&[u8]> = ctgs.iter().map(|(_, s)| s.as_slice()).collect();
                Ok((
                    sp.clone(),
                    EisensteinSamplerSketch::from_contigs(&hasher, &slices, 31, 1024),
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        let (auto_h, auto_tips) = select_auto_hinges_eisenstein(&sketches, num_h)?;
        logger.log(format_args!(
            "Eisenstein auto-hinges selected {} hinges ({:?}) + {} tips",
            auto_h.len(),
            auto_h,
            auto_tips.len()
        ));
        hinges_opt = Some(auto_h);
        tip_map_opt = Some(auto_tips);
    }

    let is_two_tier =
        hinges_opt.is_some() || tip_map_opt.is_some() || cli.extend_maf_path.is_some();

    let (merged_blocks, mode_label, tree_desc, n_sp) = if cli.progressive_ancestors {
        anyhow::ensure!(
            !is_two_tier,
            "--progressive-ancestors cannot be combined with --hinges, --tip-map, or --extend-maf"
        );
        run_progressive_ancestral_mode(&cli, logger)?
    } else if is_two_tier {
        let (blks, desc, tot_sp, n_h, n_t) =
            run_two_tier_or_extend(&cli, hinges_opt.as_deref(), tip_map_opt.as_deref(), logger)?;
        let mode = if cli.extend_maf_path.is_some() {
            format!("tba_extend hinges={n_h} tips={n_t}")
        } else {
            format!("tba_two_tier hinges={n_h} tips={n_t}")
        };
        (blks, mode, desc, tot_sp)
    } else {
        let (blks, desc, n_sp) = run_all_vs_all_tba(&cli, logger)?;
        (blks, "tba".to_owned(), desc, n_sp)
    };

    let out_prefix = &cli.out_prefix;
    let out_ext = if cli.out_gz { ".gz" } else { "" };
    let maf_path = format!("{out_prefix}.maf{out_ext}");
    let gfa_path = format!("{out_prefix}.gfa{out_ext}");
    let mut maf_writer =
        open_text_writer(Some(&maf_path)).with_context(|| format!("creating {maf_path}"))?;
    let mut gfa_writer =
        open_text_writer(Some(&gfa_path)).with_context(|| format!("creating {gfa_path}"))?;
    let (total_blocks, multi_blocks) = if cli.variation_gfa {
        let source_genomes: Vec<SourceGenomeContigs> = if cli.fasta_specs.is_empty() {
            Vec::new()
        } else {
            cli.fasta_specs
                .par_iter()
                .map(|(sp, path)| Ok((sp.clone(), read_fasta_contigs(path, cli.max_contigs)?)))
                .collect::<Result<Vec<_>>>()?
        };
        let (tot, multi, vstats) = write_tba_maf_and_variation_gfa_with_seqs(
            &merged_blocks,
            &source_genomes,
            &mut maf_writer,
            &mut gfa_writer,
        )?;
        logger.log(format_args!(
            "variation GFA: {} segments ({} bp), {} links, {} walks, {} polymorphic runs",
            vstats.segments_emitted,
            vstats.total_segment_bp,
            vstats.links_emitted,
            vstats.walks_emitted,
            vstats.polymorphic_runs
        ));
        (tot, multi)
    } else {
        write_tba_maf_and_gfa(&merged_blocks, &mut maf_writer, &mut gfa_writer)?
    };
    maf_writer.finish()?;
    gfa_writer.finish()?;

    let core_blocks = merged_blocks
        .iter()
        .filter(|b| b.rows.len() == n_sp)
        .count();
    logger.log(format_args!(
        "wrote {maf_path} & {gfa_path} ({total_blocks} total, {multi_blocks} >=3-sp, {core_blocks} {n_sp}-sp core)"
    ));
    println!(
        "mode={mode_label} species={n_sp} tree={tree_desc} total_blocks={total_blocks} multi_species_blocks={multi_blocks} core_blocks={core_blocks} elapsed_s={:.4}",
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn find_pairwise_maf_file_resolves_direct_and_swapped_and_skips_empty() {
        let dir = std::env::temp_dir().join(format!("horizon_tba_find_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("mkdir");

        // Empty file must be ignored so interrupted runs re-compute.
        std::fs::write(dir.join("spA__spB.maf"), "").expect("write empty");
        assert!(find_pairwise_maf_file(&dir, "spA", "spB").is_none());

        // Swapped file resolves with swapped=true.
        std::fs::write(dir.join("spB__spA.maf"), "##maf version=1\n").expect("write swapped");
        let (path, swapped) = find_pairwise_maf_file(&dir, "spA", "spB").expect("found swapped");
        assert!(swapped);
        assert_eq!(
            path.file_name().and_then(|s| s.to_str()),
            Some("spB__spA.maf")
        );

        // Direct non-empty file takes precedence over swapped.
        std::fs::write(dir.join("spA__spB.maf"), "##maf version=1\n").expect("write direct");
        let (path, swapped) = find_pairwise_maf_file(&dir, "spA", "spB").expect("found direct");
        assert!(!swapped);
        assert_eq!(
            path.file_name().and_then(|s| s.to_str()),
            Some("spA__spB.maf")
        );

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn load_explicit_tree_validates_inline_file_and_leaf_set() {
        let species = vec!["sp1".to_owned(), "sp2".to_owned(), "sp3".to_owned()];
        let t = load_explicit_tree("((sp1,sp2),sp3);", &species).expect("inline tree");
        assert_eq!(t.leaves().len(), 3);

        // Mismatched leaves must fail.
        assert!(load_explicit_tree("((sp1,sp2),sp4);", &species).is_err());

        // File-backed Newick loads cleanly.
        let path =
            std::env::temp_dir().join(format!("horizon_tba_tree_{}.nwk", std::process::id()));
        std::fs::write(&path, "((sp1,sp3),sp2);\n").expect("write nwk");
        let t_file =
            load_explicit_tree(path.to_str().expect("utf8 path"), &species).expect("file tree");
        assert_eq!(t_file.to_newick(), "((sp1,sp3),sp2);");
        let _ = std::fs::remove_file(&path);
    }
}
