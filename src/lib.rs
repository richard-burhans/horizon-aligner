//! Horizon Aligner: information-theoretic BOCPD/MDL genome aligner and GFA/rGFA/MAF/HAL toolkit.

pub mod bocpd;
pub mod coord;
pub mod dna;
pub mod egraph;
pub mod formats;
pub mod gfa;
pub mod hal;
pub mod io_util;
pub mod mdl;
pub mod merge;
pub mod netting;
pub mod paf;
pub mod paf_filter;
pub mod rgfa;
pub mod seeding;
pub mod tba;

/// Deterministic `splitmix64` mixing step used for `MinHash` sketching and CRC decoy sampling.
#[inline]
#[must_use]
pub const fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

pub use bocpd::bocpd_gc_profile;
pub use coord::{Interval, IntervalSet, Strand};
pub use dna::{base_to_idx, complement_base, revcomp, revcomp_preserve_case};
pub use formats::{
    HnexArchive, HnexBlock, HnexBubble, HnexBubbleSpan, HnexRow, HnexSeqMeta, MafRow,
    ResolvedDuplications, UcscChain, UcscChainBlock, UcscNetFill, UcscNetGap, UcscTargetNet,
    build_ucsc_chains, build_ucsc_nets, format_gfa_blocks_runs, format_gfa_paralog_tags,
    format_gfa_summary_runs, format_maf_blocks_runs, format_maf_blocks_runs_avoiding_gaps,
    format_ucsc_chains, format_ucsc_nets, matches_sample_prefix, parse_maf_s_line,
    parse_two_row_maf_blocks, parse_ucsc_chains, prune_resolved_dup_clusters, read_fasta_contigs,
    read_vm_hwm_kb, resolve_alignment_duplications_mdl, strip_sample_prefix,
    write_ucsc_chain_and_net_files, write_ucsc_chains, write_ucsc_nets,
};
pub use io_util::{AtomicTempFile, ProgressLogger};
pub use mdl::{
    AlignConfig, AlignedRun, DagSeedIndex, FrameshiftEvent, FrameshiftRepairResult,
    GapPosteriorTrace, GapResult, GraphAlignment, GraphNode, HmmParams, MarkovKernel, RunAlignment,
    SparsePosteriorMatrix, TrainedAsymmetricParams, align_contig_pair_runs,
    align_locus_probcons_msa, align_seq_to_dag_mdl, align_seq_to_dag_seeded_indexed,
    banded_mdl_gap, banded_mdl_gap_posterior_trace, estimate_dag_hmm_params, fixed_share_null_bits,
    repair_coding_frameshifts, train_asymmetric_pair_params, triplet_consistency_transform,
};
pub use merge::{
    DirectContigAlignment, merge_direct_alignments_prefixed_to_writers,
    merge_pairwise_maf_paths_to_writers, merge_pairwise_mafs,
    merge_pairwise_mafs_prefixed_to_writers,
};
pub use netting::{
    count_unique_covered_bp, net_global_one_to_one_runs_scored, net_one_to_one_runs_scored,
};
pub use paf::{format_paf_records_runs, format_paf_records_runs_avoiding_gaps};
pub use rgfa::{RgfaGraph, RgfaLink, RgfaSegment};
pub use seeding::{
    CdbgCollinearBlock, CdbgUnitig, ContigIndex, EISENSTEIN_OMEGA, EisensteinHasher,
    EisensteinSamplerSketch, IndexedGenome, MASK_LASTZ11, MASK_LASTZ12, MASK_PH11_18,
    MASK_PRIMARY_11_15, MASK_W9, MASK_W10, MERSENNE_61, PowerSumSketch, TantanConfig,
    UnitigOccurrence, adaptive_top_k_for_fragmentation, build_colored_cdbg, canonical_spaced_seeds,
    cluster_strand_lcbs, effective_contig_count, extract_cdbg_collinear_blocks, lis_chain,
    local_w9_subchain, merge_contig_sketches, prepare_contig_index, prepare_contig_index_owned,
    select_auto_hinges_eisenstein, select_top_k_for_target_refs, select_top_k_pairs,
    shared_seed_count, sketch_jaccard, tantan_mask, tantan_posterior_probs,
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

/// Parse and validate a positive `--threads` CLI argument (`>= 1`).
///
/// # Errors
/// Returns an error if `raw` does not parse as a `usize >= 1`.
pub fn parse_threads_arg(raw: &str) -> anyhow::Result<usize> {
    let val: usize = raw
        .parse()
        .map_err(|e| anyhow::anyhow!("parsing --threads: {e}"))?;
    anyhow::ensure!(val >= 1, "--threads must be >= 1 (got {val})");
    Ok(val)
}

/// Configure the global Rayon thread pool when `threads` is `Some(n)`.
///
/// # Errors
/// Returns an error if the global Rayon thread pool has already been initialized or cannot be built.
pub fn init_rayon_threads(threads: Option<usize>) -> anyhow::Result<()> {
    if let Some(n) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n.max(1))
            .build_global()
            .map_err(|e| anyhow::anyhow!("configuring --threads {n}: {e}"))?;
    }
    Ok(())
}

/// Validate a CLI species identifier at the system boundary so it cannot cause path traversal
/// (e.g., in `--pairwise-dir`) or inject MAF/GFA/Newick record delimiters.
///
/// # Errors
/// Returns an error if `sp` is empty, `.` or `..`, contains path separators, control characters,
/// whitespace, or Newick delimiters.
pub fn validate_species_name(sp: &str) -> anyhow::Result<()> {
    anyhow::ensure!(!sp.is_empty(), "species name cannot be empty");
    anyhow::ensure!(
        sp != "."
            && !sp.contains('/')
            && !sp.contains('\\')
            && !sp.contains("..")
            && !sp.chars().any(|c| {
                c.is_control()
                    || c.is_whitespace()
                    || matches!(c, '(' | ')' | ',' | ':' | ';' | '\'' | '"')
            }),
        "invalid species name `{sp}`: must not contain path separators, control characters, whitespace, or Newick delimiters"
    );
    Ok(())
}

/// Parse a CLI `<species>=<path>` specification, validating the species identifier and non-empty path.
///
/// # Errors
/// Returns an error if `spec` lacks `=`, has an invalid species name, or has an empty path.
pub fn parse_species_fasta_spec(spec: &str) -> anyhow::Result<(String, std::path::PathBuf)> {
    let Some((sp, path)) = spec.split_once('=') else {
        anyhow::bail!("expected <species>=<fasta>, got `{spec}`");
    };
    validate_species_name(sp)?;
    anyhow::ensure!(
        !path.is_empty(),
        "expected non-empty <species>=<fasta>, got `{spec}`"
    );
    Ok((sp.to_owned(), std::path::PathBuf::from(path)))
}

/// Contig-pair scheduling policy shared across `horizon_wga`, `horizon_guide`, and `horizon_tba`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairSchedule {
    /// Retain up to `k` highest-scoring partner contigs per contig (`k = 1` is `--rbest`).
    TopK(usize),
    /// Align every contig pair passing the shared-seed floor (`--all-pairs`).
    AllPairs,
}

impl PairSchedule {
    /// Resolve an optional CLI scheduling override against `max_contigs`, returning the effective
    /// `top_k` integer (`0` = all pairs, `1` = reciprocal best hit when `max_contigs != 1`).
    #[must_use]
    pub const fn resolve_top_k(opt: Option<Self>, max_contigs: usize) -> usize {
        match opt {
            Some(Self::TopK(k)) => k,
            Some(Self::AllPairs) => 0,
            None => {
                if max_contigs == 1 {
                    0
                } else {
                    1
                }
            }
        }
    }
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
    use crate::seeding::{MASK_RY18, append_hsp_gated_aa_lcbs, find_translated_aa5_hits_both};

    #[test]
    fn bocpd_matches_python_oracle() {
        let seq = b"ATATATATATATATATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATAT";
        let gc = bocpd_gc_profile(seq, 1.0 / 120.0, 200);
        assert_eq!(gc.len(), seq.len());
        let expected_first5 = [0.5, 0.334_722_22, 0.252_600_90, 0.203_511_81, 0.170_862_46];
        for (got, exp) in gc.iter().take(5).zip(expected_first5) {
            assert!((got - exp).abs() < 1e-6, "got {got}, expected {exp}");
        }
        let expected_mid5 = [
            0.821_913_23,
            0.849_528_15,
            0.868_348_45,
            0.882_449_15,
            0.893_556_68,
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
        let mut scratch = mdl::xdrop::XDropScratch::default();
        let pairs = mdl::xdrop::mdl_xdrop_extend_with_scratch(
            sx,
            sy,
            0,
            0,
            false,
            100,
            16,
            10.0,
            &HmmParams::default(),
            &mut scratch,
        );
        assert!(pairs.len() >= 30);
    }

    #[test]
    fn gfa_summary_contains_header_and_bubble_tag() {
        let aln = RunAlignment {
            aligned_runs: vec![AlignedRun {
                x_start: 0,
                y_start: 0,
                len: 2,
                is_rc: false,
            }],
            total_pairs: 2,
            lcbs_fwd: 1,
            lcbs_rc: 0,
            macro_bubbles: 2,
            internal_bubbles: 3,
        };
        let gfa = format_gfa_summary_runs("chr1_x", "chr1_y", &aln);
        assert!(gfa.starts_with("H\tVN:Z:1.1\n"));
        assert!(gfa.contains("HB:i:5"));
    }

    #[test]
    fn run_alignment_total_pairs_matches_expanded_runs() {
        let sx = b"ATATATATATATATATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATATATATATATATATATATGCGCGCGCGCGCGCGCGCGC";
        let sy = b"ATATATATATATCTATATATGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATATATATATATATATATATGCGCGCGCGCGCGCGCGCGC";
        let ix = prepare_contig_index(sx).expect("valid");
        let iy = prepare_contig_index(sy).expect("valid");
        let cfg = AlignConfig::default();
        let runs = align_contig_pair_runs(&ix, &iy, &cfg).expect("valid");
        assert_eq!(
            crate::mdl::runs_to_pairs(&runs.aligned_runs).len(),
            runs.total_pairs
        );
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
        let aln = RunAlignment {
            aligned_runs: vec![AlignedRun {
                x_start: 0,
                y_start: 0,
                len: 16,
                is_rc: false,
            }],
            total_pairs: 16,
            lcbs_fwd: 1,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let maf = format_maf_blocks_runs("chr1", 16, "chr2", 16, sx, sy, &aln);
        assert!(maf.starts_with("##maf version=1"));
        assert!(maf.contains("s chr1 0 16 + 16 ACGTACGTACGTACGT"));
        assert!(maf.contains("s chr2 0 16 + 16 ACGTACGTACGTACGT"));
    }

    #[test]
    fn gfa_blocks_emits_segments_bubbles_and_links() {
        let aln = RunAlignment {
            aligned_runs: vec![
                AlignedRun {
                    x_start: 0,
                    y_start: 0,
                    len: 15,
                    is_rc: false,
                },
                AlignedRun {
                    x_start: 100,
                    y_start: 100,
                    len: 15,
                    is_rc: false,
                },
            ],
            total_pairs: 30,
            lcbs_fwd: 1,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 1,
        };
        let gfa = format_gfa_blocks_runs("c1", "c2", &aln);
        assert!(gfa.contains("S\tc1__c2__blk_0"));
        assert!(gfa.contains("S\tc1__c2__bub_1\t*\tLN:i:85\tHB:i:1"));
        assert!(gfa.contains("L\tc1__c2__blk_0\t+\tc1__c2__bub_1\t+\t0M"));
    }

    #[test]
    fn maf_and_gfa_handle_duplicate_x_and_strand_flips() {
        let sx = vec![b'A'; 64];
        let sy = vec![b'A'; 64];
        // Duplicate/overlapping x and a strand flip mid-stream should not panic.
        let aln = RunAlignment {
            aligned_runs: vec![
                AlignedRun {
                    x_start: 0,
                    y_start: 30,
                    len: 14,
                    is_rc: true,
                },
                AlignedRun {
                    x_start: 5,
                    y_start: 5,
                    len: 1,
                    is_rc: false,
                },
                AlignedRun {
                    x_start: 14,
                    y_start: 31,
                    len: 12,
                    is_rc: false,
                },
            ],
            total_pairs: 27,
            lcbs_fwd: 1,
            lcbs_rc: 1,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let maf = format_maf_blocks_runs("c1", 64, "c2", 64, &sx, &sy, &aln);
        let gfa = format_gfa_blocks_runs("c1", "c2", &aln);
        assert!(maf.contains("##maf version=1"));
        assert!(gfa.contains("H\tVN:Z:1.1"));
    }

    #[test]
    fn min_top_lcb_prunes_small_lcb_pairs() {
        let sx = b"ACGTACGTACGTACGTGCGCGCGCATATATATACGTACGTACGTACGT";
        let sy = b"ACGTACGTACGTACGTGCGCGCGCATATATATACGTACGTACGTACGT";
        let ix = prepare_contig_index(sx).expect("valid");
        let iy = prepare_contig_index(sy).expect("valid");
        let cfg = AlignConfig {
            min_top_lcb_anchors: 10_000,
            ..AlignConfig::default()
        };
        let pruned = align_contig_pair_runs(&ix, &iy, &cfg).expect("valid");
        assert_eq!(pruned.total_pairs, 0);
    }

    #[test]
    fn align_contig_pair_runs_handles_short_and_empty_sequences_without_panic() {
        let ref_short = prepare_contig_index(b"ACGTACGT").expect("idx");
        let qry_short = prepare_contig_index(b"ACGT").expect("idx");
        let aln_short =
            align_contig_pair_runs(&ref_short, &qry_short, &AlignConfig::default()).expect("short");
        assert_eq!(aln_short.total_pairs, 0);

        let ref_empty = prepare_contig_index(b"").expect("idx");
        let qry_empty = prepare_contig_index(b"ACGTACGT").expect("idx");
        let aln_empty =
            align_contig_pair_runs(&ref_empty, &qry_empty, &AlignConfig::default()).expect("empty");
        assert_eq!(aln_empty.total_pairs, 0);
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
        let ry_x = seeding::canonical_spaced_seeds_with_opts(sx, MASK_RY18, 2, false, true)
            .expect("valid");
        let ry_y = seeding::canonical_spaced_seeds_with_opts(sy, MASK_RY18, 2, false, true)
            .expect("valid");
        assert!(!ry_x.is_empty());
        assert_eq!(ry_x, ry_y);

        let masked = b"acgttgcactgatcgtagctagctaaggccat";
        let skipped = seeding::canonical_spaced_seeds_with_opts(masked, MASK_RY18, 2, true, true)
            .expect("valid");
        assert_eq!(skipped.len(), 0);
    }

    #[test]
    fn lastz12_1transition_recovers_transition_hit() {
        let sx = b"ACGTTGCACTGATCGTAGCTAGCTAAGGCCAT";
        // Position 1 has C -> T transition; exact 12-of-19 misses, 1-transition recovers (0, 0).
        let sy = b"ATGTTGCACTGATCGTAGCTAGCTAAGGCCAT";
        let y_idx =
            seeding::FlatSeedIndex::build_with_sdust(sy, MASK_LASTZ12, 2, false, false, None)
                .expect("valid");
        let (fwd, _) = seeding::find_lastz12_1transition_hits_flat(sx, None, &y_idx, false);
        assert!(fwd.contains(&(0, 0)), "expected (0, 0) in {fwd:?}");
    }

    #[test]
    fn lastz12_2transition_recovers_double_transition_hit() {
        let sx = b"ACGTTGCACTGATCGTAGCTAGCTAAGGCCAT";
        // Positions 1 (C->T) and 11 (A->G) both have transitions; 1-transition misses, 2-transition recovers (0, 0).
        let sy = b"ATGTTGCACTGGTCGTAGCTAGCTAAGGCCAT";
        let y_idx =
            seeding::FlatSeedIndex::build_with_sdust(sy, MASK_LASTZ12, 2, false, false, None)
                .expect("valid");
        let (fwd1, _) = seeding::find_lastz12_1transition_hits_flat(sx, None, &y_idx, false);
        assert!(!fwd1.contains(&(0, 0)), "1-transition should miss (0, 0)");
        let (fwd2, _) = seeding::find_lastz12_2transition_hits_flat(sx, None, &y_idx, false);
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
        let sy_rc = revcomp_preserve_case(sy);
        let (aa_fwd, _) = find_translated_aa5_hits_both(sx, sy, &sy_rc, 2, false);
        assert!(
            !aa_fwd.is_empty(),
            "expected translated 5-mer AA hits for synonymous coding pair"
        );
        let mut aa_gated = Vec::new();
        append_hsp_gated_aa_lcbs(&mut aa_gated, sx, sy, &aa_fwd, 90);
        assert!(
            !aa_gated.is_empty(),
            "expected AA HSP gate to accept synonymous 12-codon stretch"
        );
    }

    #[test]
    fn cli_arg_parsers_validate_bounds() {
        assert_eq!(parse_min_fallback_cov_arg("0.0").expect("0.0"), 0.0);
        assert_eq!(parse_min_fallback_cov_arg("1.0").expect("1.0"), 1.0);
        assert!(parse_min_fallback_cov_arg("NaN").is_err());
        assert!(parse_min_fallback_cov_arg("-0.1").is_err());
        assert!(parse_min_fallback_cov_arg("1.1").is_err());

        assert_eq!(parse_hsp_threshold_arg("0").expect("0"), 0);
        assert_eq!(parse_hsp_threshold_arg("42").expect("42"), 42);
        assert!(parse_hsp_threshold_arg("-1").is_err());
        assert!(parse_hsp_threshold_arg("abc").is_err());

        assert_eq!(parse_threads_arg("1").expect("1"), 1);
        assert_eq!(parse_threads_arg("8").expect("8"), 8);
        assert!(parse_threads_arg("0").is_err());
        assert!(parse_threads_arg("-2").is_err());
        assert!(parse_threads_arg("abc").is_err());

        for ok_name in ["P_falciparum", "hg38-chr1", "sp1.v2"] {
            assert!(
                validate_species_name(ok_name).is_ok(),
                "expected ok: {ok_name}"
            );
        }
        for bad_name in [
            "", ".", "..", "../evil", "sp/a", "sp\\a", "sp\0a", "sp\na", "sp a", "sp(1)", "sp:1",
            "sp;1",
        ] {
            assert!(
                validate_species_name(bad_name).is_err(),
                "expected rejection for `{bad_name:?}`"
            );
        }
    }
}
