//! Horizon Aligner: KP-14 spaced-seed + BOCPD/MDL genome aligner and GFA/MAF/HAL toolkit.

pub mod bocpd;
pub mod coord;
pub mod egraph;
pub mod formats;
pub mod gfa;
pub mod hal;
pub mod mdl;
pub mod merge;
pub mod netting;
pub mod paf;
pub mod seeding;

pub use bocpd::bocpd_gc_profile;
pub use formats::{
    MafRow, format_gfa_blocks, format_gfa_blocks_runs, format_gfa_paralog_tags,
    format_gfa_summary, format_gfa_summary_runs, format_maf_blocks, format_maf_blocks_runs,
    matches_sample_prefix, parse_maf_s_line, parse_two_row_maf_blocks, read_fasta_contigs,
    read_vm_hwm_kb, resolve_alignment_duplications_runs, resolve_lcb_duplications,
    strip_sample_prefix,
};
pub use netting::net_one_to_one_runs;
pub use mdl::{
    AlignConfig, AlignedRun, ContigAlignment, DagSeedIndex, GapResult, GraphAlignment, GraphNode,
    HmmParams, RunAlignment, align_contig_pair, align_contig_pair_indexed,
    align_contig_pair_indexed_with_min_lcb, align_contig_pair_runs,
    align_contig_pair_with_config, align_seq_to_dag_mdl, align_seq_to_dag_seeded,
    align_seq_to_dag_seeded_indexed, banded_mdl_gap, banded_mdl_gap_with_params, mdl_xdrop_extend,
    mdl_xdrop_extend_with_params,
};
pub use merge::merge_pairwise_mafs;
pub use seeding::{
    ContigIndex, MASK_LASTZ11, MASK_LASTZ12, MASK_PH11_18, MASK_PRIMARY_11_15, MASK_W9, MASK_W10,
    canonical_spaced_seeds, canonical_spaced_seeds_with_opts, cluster_strand_lcbs, complement_base,
    find_lastz12_1transition_hits, find_lastz12_2transition_hits, lis_chain, local_w9_subchain,
    prepare_contig_index, prepare_contig_index_opts, revcomp, revcomp_preserve_case,
    score_contig_pair_lcb, select_top_k_pairs, shared_seed_count,
};

use thiserror::Error;

/// Parse and validate a `--min-fallback-cov` CLI argument in `0.0..=1.0`.
///
/// # Errors
/// Returns an error if `raw` does not parse as a finite `f64` within `0.0..=1.0`.
pub fn parse_min_fallback_cov_arg(raw: &str) -> anyhow::Result<f64> {
    let val: f64 = raw
        .parse()
        .map_err(|e| anyhow::anyhow!("parsing --min-fallback-cov: {e}"))?;
    anyhow::ensure!(
        val.is_finite() && (0.0..=1.0).contains(&val),
        "--min-fallback-cov must be a finite number in 0.0..=1.0 (got {val})"
    );
    Ok(val)
}

/// Parse and validate a non-negative `--hsp-threshold` CLI argument.
///
/// # Errors
/// Returns an error if `raw` does not parse as a non-negative `i32`.
pub fn parse_hsp_threshold_arg(raw: &str) -> anyhow::Result<i32> {
    let val: i32 = raw
        .parse()
        .map_err(|e| anyhow::anyhow!("parsing --hsp-threshold: {e}"))?;
    anyhow::ensure!(val >= 0, "--hsp-threshold must be non-negative (got {val})");
    Ok(val)
}

/// Errors returned by sequence validation or alignment routines.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum AlignError {
    /// Raised when an empty mask or invalid mask character is supplied.
    #[error("invalid spaced-seed mask: {0}")]
    InvalidMask(String),
    /// Raised when an input sequence or auxiliary mask violates size or alignment invariants.
    #[error("invalid alignment input: {0}")]
    InvalidInput(String),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seeding::{
        MASK_RY20, append_hsp_gated_aa_lcbs, append_translated_aa5_lcbs, find_translated_aa5_hits,
    };

    #[test]
    fn bocpd_matches_python_oracle() {
        let seq = b"ATATATATATATATATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATAT";
        let gc = bocpd_gc_profile(seq, 1.0 / 120.0, 200);
        assert_eq!(gc.len(), seq.len());
        let expected_first5 = [0.5, 0.333_333_33, 0.250_694_44, 0.201_478_74, 0.168_858_77];
        for (got, exp) in gc.iter().take(5).zip(expected_first5) {
            assert!((got - exp).abs() < 1e-6, "got {got}, expected {exp}");
        }
        let expected_mid5 = [
            0.733_457_91,
            0.776_805_99,
            0.806_661_98,
            0.829_114_1,
            0.846_848_01,
        ];
        for (got, exp) in gc[25..30].iter().zip(expected_mid5) {
            assert!((got - exp).abs() < 1e-6, "got {got}, expected {exp}");
        }
    }

    #[test]
    fn spaced_seeds_are_strand_symmetric_and_match_oracle() {
        let seq = b"ATATATATATATATATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATAT";
        let seeds = canonical_spaced_seeds(seq, MASK_PRIMARY_11_15, 2).expect("valid mask");
        // 60bp palindrome: 14 boundary windows on the left half each pair with their RC on the right half.
        assert_eq!(seeds.len(), 14);

        let rc: Vec<u8> = seq
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                other => other,
            })
            .collect();
        let rc_seeds = canonical_spaced_seeds(&rc, MASK_PRIMARY_11_15, 2).expect("valid mask");
        assert_eq!(seeds.len(), rc_seeds.len());
        for key in seeds.keys() {
            assert!(rc_seeds.contains_key(key));
        }
    }

    #[test]
    fn lis_chain_finds_longest_colinear_path() {
        let hits = [(10, 100), (20, 50), (30, 120), (40, 130), (50, 90)];
        let chain = lis_chain(&hits);
        assert_eq!(chain.len(), 3);
        for w in chain.windows(2) {
            assert!(w[0].0 < w[1].0 && w[0].1 < w[1].1);
        }
    }

    #[test]
    fn banded_mdl_abstains_on_random_at_rich_pair() {
        let mut state: u64 = 12345;
        let mut next_base = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            let r = (state >> 33) % 100;
            if r < 41 {
                b'A'
            } else if r < 82 {
                b'T'
            } else if r < 91 {
                b'C'
            } else {
                b'G'
            }
        };
        let sx: Vec<u8> = (0..160).map(|_| next_base()).collect();
        let sy: Vec<u8> = (0..160).map(|_| next_base()).collect();
        let res = banded_mdl_gap(&sx, &sy, 48, 12.0);
        assert!(
            res.kept_pairs.len() <= 40,
            "expected heavy abstention (>=75%), got {}/160",
            res.kept_pairs.len()
        );
        assert!(res.bubbles >= 1);
    }

    #[test]
    fn banded_mdl_recovers_conserved_pair() {
        let sx = b"ATATATATATATATATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATAT";
        let sy = b"ATATATATATATCTATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATAT";
        let res = banded_mdl_gap(sx, sy, 24, 10.0);
        assert!(res.score_bits > 60.0);
        assert_eq!(res.kept_pairs.len(), 60);
        assert_eq!(res.bubbles, 0);
    }

    #[test]
    fn cluster_lcbs_separates_colinear_blocks() {
        let hits = [
            (0, 0),
            (30, 30),
            (60, 60),
            (90, 90),
            (120, 120),
            (150, 150),
            (180, 180),
            (210, 210),
        ];
        let lcbs = cluster_strand_lcbs(&hits, 500, 50, 8, 100, 100.0);
        assert_eq!(lcbs.len(), 1);
        assert_eq!(lcbs[0].len(), 8);
    }

    #[test]
    fn local_w9_requires_two_hit_diagonal_support() {
        let sx = b"AAAAAACCCCCCGGGGGGTTTTTTACGTACGTACGTACGT";
        let sy = b"AAAAAACCCCCCGGGGGGTTTTTTACGTACGTACGTACGT";
        let anchors = local_w9_subchain(sx, sy, 10, 20, MASK_W9, 80, 600).expect("valid mask");
        assert!(anchors.len() >= 2);
        assert!(anchors[0].0 >= 10 && anchors[0].1 >= 20);
    }

    #[test]
    fn mdl_xdrop_extends_conserved_flank() {
        let sx = b"ATATATATGCGCGCGCATATATATGCGCGCGCATATATAT";
        let sy = b"ATATATATGCGCGCGCATATATATGCGCGCGCATATATAT";
        let pairs = mdl_xdrop_extend(sx, sy, 0, 0, false, 100, 16);
        assert!(pairs.len() >= 30);
    }

    #[test]
    fn gfa_summary_contains_header_and_bubble_tag() {
        let aln = ContigAlignment {
            aligned_pairs: vec![(0, 0), (1, 1)],
            lcbs_fwd: 1,
            lcbs_rc: 0,
            macro_bubbles: 2,
            internal_bubbles: 3,
        };
        let gfa = format_gfa_summary("chr1_x", "chr1_y", &aln);
        assert!(gfa.starts_with("H\tVN:Z:1.1\n"));
        assert!(gfa.contains("HB:i:5"));
    }

    #[test]
    fn indexed_alignment_matches_unindexed() {
        let sx = b"ATATATATATATATATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATATATATATATATATATATGCGCGCGCGCGCGCGCGCGC";
        let sy = b"ATATATATATATCTATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATATATATATATATATATATGCGCGCGCGCGCGCGCGCGC";
        let unindexed = align_contig_pair(sx, sy).expect("valid");
        let ix = prepare_contig_index(sx).expect("valid");
        let iy = prepare_contig_index(sy).expect("valid");
        let indexed = align_contig_pair_indexed(&ix, &iy).expect("valid");
        assert_eq!(unindexed, indexed);
    }

    #[test]
    fn shared_seed_count_separates_homologs_from_unrelated() {
        let sx = b"ACGTTGCATGCAACTGGTCAGTCGATCGTAGCTAGCTGATCGATCGATGCTAGCATCGAT";
        let sy_same = b"ACGTTGCATGCAACTGGTCAGTCGATCGTCGCTAGCTGATCGATCGATGCTAGCATCGAT";
        let sy_diff = b"CCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCC";
        let ix = prepare_contig_index(sx).expect("valid");
        let iy_same = prepare_contig_index(sy_same).expect("valid");
        let iy_diff = prepare_contig_index(sy_diff).expect("valid");
        assert!(shared_seed_count(&ix, &iy_same) >= 10);
        assert_eq!(shared_seed_count(&ix, &iy_diff), 0);
    }

    #[test]
    fn maf_formatter_emits_valid_blocks() {
        let sx = b"ACGTACGTACGTACGT";
        let sy = b"ACGTACGTACGTACGT";
        let aln = ContigAlignment {
            aligned_pairs: (0..16).map(|i| (i, i)).collect(),
            lcbs_fwd: 1,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let maf = format_maf_blocks("chr1", 16, "chr2", 16, sx, sy, &aln);
        assert!(maf.starts_with("##maf version=1"));
        assert!(maf.contains("s chr1 0 16 + 16 ACGTACGTACGTACGT"));
        assert!(maf.contains("s chr2 0 16 + 16 ACGTACGTACGTACGT"));
    }

    #[test]
    fn gfa_blocks_emits_segments_bubbles_and_links() {
        let mut pairs: Vec<(usize, usize)> = (0..15).map(|i| (i, i)).collect();
        pairs.extend((100..115).map(|i| (i, i)));
        let aln = ContigAlignment {
            aligned_pairs: pairs,
            lcbs_fwd: 1,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 1,
        };
        let gfa = format_gfa_blocks("c1", "c2", &aln);
        assert!(gfa.contains("S\tc1__c2__blk_0"));
        assert!(gfa.contains("S\tc1__c2__bub_1\t*\tLN:i:85\tHB:i:1"));
        assert!(gfa.contains("L\tc1__c2__blk_0\t+\tc1__c2__bub_1\t+\t0M"));
    }

    #[test]
    fn maf_and_gfa_handle_duplicate_x_and_strand_flips() {
        let sx = vec![b'A'; 64];
        let sy = vec![b'A'; 64];
        // Duplicate x=5 and a strand flip mid-stream should not panic.
        let mut pairs: Vec<(usize, usize)> = (0..14).map(|i| (i, 30 - i)).collect();
        pairs.push((5, 5));
        pairs.push((14, 31));
        pairs.sort_unstable();
        let aln = ContigAlignment {
            aligned_pairs: pairs,
            lcbs_fwd: 0,
            lcbs_rc: 1,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let maf = format_maf_blocks("c1", 64, "c2", 64, &sx, &sy, &aln);
        let gfa = format_gfa_blocks("c1", "c2", &aln);
        assert!(maf.contains("##maf version=1"));
        assert!(gfa.contains("H\tVN:Z:1.1"));
    }

    #[test]
    fn min_top_lcb_prunes_small_lcb_pairs() {
        let sx = b"ACGTACGTACGTACGTGCGCGCGCATATATATACGTACGTACGTACGT";
        let sy = b"ACGTACGTACGTACGTGCGCGCGCATATATATACGTACGTACGTACGT";
        let ix = prepare_contig_index(sx).expect("valid");
        let iy = prepare_contig_index(sy).expect("valid");
        let pruned = align_contig_pair_indexed_with_min_lcb(&ix, &iy, 10_000).expect("valid");
        assert_eq!(pruned.aligned_pairs.len(), 0);
    }

    #[test]
    fn resolve_dups_and_merge_pairwise_mafs_work() {
        let maf1 = "##maf version=1\na score=12.0\ns chr1 100 12 + 5000 ACGTACGTACGT\ns q1 10 12 + 5000 ACGTACGTACGT\n\na score=12.0\ns chr1 200 12 + 5000 ACGTACGTACGT\ns q1 200 12 + 5000 ACGTACGTACGT\n";
        let maf2 = "##maf version=1\na score=12.0\ns chr1 105 12 + 5000 CGTACGTACGTA\ns q2 50 12 - 5000 CGTACGTACGTA\n";
        let (merged_maf, merged_gfa, total, multi) = merge_pairwise_mafs(
            "Pf",
            &[
                ("Pr".to_owned(), maf1.to_owned()),
                ("Pv".to_owned(), maf2.to_owned()),
            ],
            1000,
        );
        assert_eq!(total, 4);
        assert_eq!(multi, 1);
        assert!(merged_maf.contains("s chr1 105 7 + 5000 CGTACGT"));
        assert!(merged_maf.contains("s Pr.q1"));
        assert!(merged_maf.contains("s Pv.q2"));
        assert!(merged_gfa.contains("# EG:Z:Pf__Pr\tB0"));
    }

    #[test]
    fn merge_threads_query_insertion_columns() {
        let maf1 =
            "##maf version=1\na score=4.0\ns chr1 10 4 + 100 AC--GT\ns q1 0 6 + 100 ACTTGT\n";
        let maf2 = "##maf version=1\na score=4.0\ns chr1 10 4 + 100 ACGT\ns q2 0 4 + 100 ACGT\n";
        let (merged_maf, _, total, multi) = merge_pairwise_mafs(
            "Pf",
            &[
                ("Pr".to_owned(), maf1.to_owned()),
                ("Pv".to_owned(), maf2.to_owned()),
            ],
            1000,
        );
        assert_eq!(total, 1);
        assert_eq!(multi, 1);
        assert!(merged_maf.contains("s chr1 10 4 + 100 AC--GT"));
        assert!(merged_maf.contains("s Pr.q1 0 6 + 100 ACTTGT"));
        assert!(merged_maf.contains("s Pv.q2 0 4 + 100 AC--GT"));
    }

    #[test]
    fn ry_seeds_recover_transitions_and_softmask_skips_lowercase() {
        // Transitions (A<->G, C<->T) at multiple positions: R/Y seeds still match identically.
        let sx = b"ACGTTGCACTGATCGTAGCTAGCTAAGGCCAT";
        let sy = b"ATGTTGTACTGATTGTAGCTAGCTAAGGCTAT";
        let ry_x = canonical_spaced_seeds_with_opts(sx, MASK_RY20, 2, false, true).expect("valid");
        let ry_y = canonical_spaced_seeds_with_opts(sy, MASK_RY20, 2, false, true).expect("valid");
        assert!(!ry_x.is_empty());
        assert_eq!(ry_x, ry_y);

        let masked = b"acgttgcactgatcgtagctagctaaggccat";
        let skipped =
            canonical_spaced_seeds_with_opts(masked, MASK_RY20, 2, true, true).expect("valid");
        assert_eq!(skipped.len(), 0);
    }

    #[test]
    fn lastz12_1transition_recovers_transition_hit() {
        let sx = b"ACGTTGCACTGATCGTAGCTAGCTAAGGCCAT";
        // Position 1 has C -> T transition; exact 12-of-19 misses, 1-transition recovers (0, 0).
        let sy = b"ATGTTGCACTGATCGTAGCTAGCTAAGGCCAT";
        let iy_fwd = canonical_spaced_seeds(sy, MASK_LASTZ12, 2).expect("valid");
        let iy_rc = canonical_spaced_seeds(&revcomp(sy), MASK_LASTZ12, 2).expect("valid");
        let (fwd, _) = find_lastz12_1transition_hits(sx, &iy_fwd, &iy_rc, false);
        assert!(fwd.contains(&(0, 0)), "expected (0, 0) in {fwd:?}");
    }

    #[test]
    fn lastz12_2transition_recovers_double_transition_hit() {
        let sx = b"ACGTTGCACTGATCGTAGCTAGCTAAGGCCAT";
        // Positions 1 (C->T) and 11 (A->G) both have transitions; 1-transition misses, 2-transition recovers (0, 0).
        let sy = b"ATGTTGCACTGGTCGTAGCTAGCTAAGGCCAT";
        let iy_fwd = canonical_spaced_seeds(sy, MASK_LASTZ12, 2).expect("valid");
        let iy_rc = canonical_spaced_seeds(&revcomp(sy), MASK_LASTZ12, 2).expect("valid");
        let (fwd1, _) = find_lastz12_1transition_hits(sx, &iy_fwd, &iy_rc, false);
        assert!(!fwd1.contains(&(0, 0)), "1-transition should miss (0, 0)");
        let (fwd2, _) = find_lastz12_2transition_hits(sx, &iy_fwd, &iy_rc, false);
        assert!(fwd2.contains(&(0, 0)), "expected (0, 0) in {fwd2:?}");
    }

    #[test]
    fn merge_nets_best_scoring_overlapping_block() {
        let weak =
            "##maf version=1\na score=1.0\ns chr1 0 6 + 100 AAAAAA\ns q_weak 0 6 + 100 TTTTTT\n";
        let strong =
            "##maf version=1\na score=6.0\ns chr1 0 6 + 100 AAAAAA\ns q_best 0 6 + 100 AAAAAA\n";
        let combined = format!("{weak}\n{strong}");
        let (merged_maf, _, total, _) =
            merge_pairwise_mafs("hg38", &[("mm10".to_owned(), combined)], 1000);
        assert_eq!(total, 1);
        assert!(merged_maf.contains("s mm10.q_best"));
        assert!(!merged_maf.contains("s mm10.q_weak"));
    }

    #[test]
    fn translated_aa5_recovers_synonymous_divergent_coding_pair() {
        // Two coding sequences encoding the same 12-aa peptide with synonymous 3rd-codon swaps.
        let sx = b"GCTACTGGTAAAGTTGACTTCCCAATGCAACGTTGG";
        let sy = b"GCCACAGGCAAAGTCGATTTTCCTATGCAGAGGTGG";
        let mut lcbs = Vec::new();
        append_translated_aa5_lcbs(&mut lcbs, sx, sy, 2, false);
        assert!(
            !lcbs.is_empty(),
            "expected translated 5-mer AA LCB for synonymous coding pair"
        );
        let mut aa_gated = Vec::new();
        let aa_hits = find_translated_aa5_hits(sx, sy, 2, false);
        append_hsp_gated_aa_lcbs(&mut aa_gated, sx, sy, &aa_hits, 90);
        assert!(
            !aa_gated.is_empty(),
            "expected AA HSP gate to accept synonymous 12-codon stretch"
        );
    }
}
