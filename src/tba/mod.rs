//! Reference-free Threaded Blockset Aligner (`TBA`) progressive multi-species merger.
//!
//! Merges all-vs-all pairwise MAF blocks bottom-up along a binary [`GuideTree`] without
//! anchoring on any single reference genome: clade-specific blocks absent from other clades
//! are preserved intact while cross-clade pairwise links thread multi-species columns together.

mod hnex;
mod merge;
mod model;
mod poa_refine;
pub mod progressive;
pub mod tree;
pub mod two_tier;
pub mod variation_gfa;
pub mod writer;

pub use merge::{
    run_tba_tree_merge, run_tba_tree_merge_consuming, run_tba_tree_merge_lazy,
    tba_merge_two_blocksets, tba_project_tips_onto_blockset, tba_project_tips_streaming,
};
pub use model::{
    PairwiseTbaBlock, TbaBlock, TbaRow, canonical_pair_key, for_each_multi_tba_block_reader,
    parse_multi_tba_reader, parse_pairwise_tba_reader,
};
pub use poa_refine::refine_shared_tip_insertions_poa;
pub use progressive::{
    AncestralBranchLengths, AncestralNodeResult, AncestralNodeStats, ProgressiveAncestralConfig,
    ProgressiveAncestralResult, format_ancestors_summary_tsv, reconstruct_ancestral_contigs,
    reconstruct_ancestral_contigs_transducer, run_progressive_ancestral_tba, select_node_outgroups,
};
pub use tree::{
    GuideTree, aligned_bp_similarity_matrix, build_bionj_tree, build_upgma_tree, parse_newick,
    similarity_submatrix,
};
pub use two_tier::{TwoTierPlan, parse_hinge_list, parse_tip_map_reader};
pub use variation_gfa::{
    SourceGenomeContigs, VariationGfaBuilder, VariationGfaStats, write_tba_maf_and_variation_gfa,
    write_tba_maf_and_variation_gfa_with_seqs, write_variation_gfa,
    write_variation_gfa_from_maf_reader,
};
pub use writer::write_tba_maf_and_gfa;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::formats::MafRow;
    use std::collections::{HashMap, HashSet};

    #[allow(clippy::too_many_arguments)]
    fn mk_pair(
        sp_a: &str,
        start_a: usize,
        strand_a: char,
        seq_a: &str,
        sp_b: &str,
        start_b: usize,
        strand_b: char,
        seq_b: &str,
    ) -> PairwiseTbaBlock {
        let ra = MafRow {
            src: format!("{sp_a}.chr1"),
            start: start_a,
            strand: strand_a,
            src_size: 100,
            seq: seq_a.to_owned(),
        };
        let rb = MafRow {
            src: format!("{sp_b}.chr1"),
            start: start_b,
            strand: strand_b,
            src_size: 100,
            seq: seq_b.to_owned(),
        };
        PairwiseTbaBlock::from_maf_rows(sp_a, ra, sp_b, rb).expect("valid pair")
    }

    #[test]
    fn tba_preserves_clade_specific_blocks_and_fuses_4way_core() {
        // Tree: ((lav1,lav2),(vin1,vin2));
        // [0..10): pan-genus 4-way conserved core
        // [20..30): lav-only block (absent from vin1/vin2)
        // [40..50): vin-only block (absent from lav1/lav2)
        let tree = parse_newick("((lav1,lav2),(vin1,vin2));").expect("tree");
        let mut pairs: HashMap<(String, String), Vec<PairwiseTbaBlock>> = HashMap::new();
        pairs.insert(
            canonical_pair_key("lav1", "lav2"),
            vec![
                mk_pair("lav1", 0, '+', "ACGTACGTAC", "lav2", 0, '+', "ACGTACGTAC"),
                mk_pair("lav1", 20, '+', "GGGGGGGGGG", "lav2", 20, '+', "GGGGGGGGGG"),
            ],
        );
        pairs.insert(
            canonical_pair_key("vin1", "vin2"),
            vec![
                mk_pair("vin1", 0, '+', "ACGTACGTAC", "vin2", 0, '+', "ACGTACGTAC"),
                mk_pair("vin1", 40, '+', "CCCCCCCCCC", "vin2", 40, '+', "CCCCCCCCCC"),
            ],
        );
        pairs.insert(
            canonical_pair_key("lav1", "vin1"),
            vec![mk_pair(
                "lav1",
                0,
                '+',
                "ACGTACGTAC",
                "vin1",
                0,
                '+',
                "ACGTACGTAC",
            )],
        );

        let merged = run_tba_tree_merge(&tree, &pairs);
        assert_eq!(merged.len(), 3, "expected 4-way core + lav-only + vin-only");
        let depths: Vec<usize> = {
            let mut d: Vec<usize> = merged.iter().map(|b| b.rows.len()).collect();
            d.sort_unstable();
            d
        };
        assert_eq!(depths, vec![2, 2, 4]);
        // Confirm vin-only clade block survived even though lav1/lav2 lack it completely!
        assert!(merged.iter().any(|b| {
            b.rows.len() == 2
                && b.rows.iter().any(|r| r.species == "vin1" && r.start == 40)
                && b.rows.iter().any(|r| r.species == "vin2" && r.start == 40)
        }));
    }

    #[test]
    fn tba_threads_insertions_and_minus_strand_accurately() {
        let tree = parse_newick("((sp1,sp2),sp3);").expect("tree");
        let mut pairs: HashMap<(String, String), Vec<PairwiseTbaBlock>> = HashMap::new();
        // sp1 vs sp2 has a 2bp insertion in sp2.
        pairs.insert(
            canonical_pair_key("sp1", "sp2"),
            vec![mk_pair(
                "sp1",
                0,
                '+',
                "ACGT--ACGT",
                "sp2",
                0,
                '+',
                "ACGTTTACGT",
            )],
        );
        // sp1 vs sp3 aligns on minus strand without the insertion.
        pairs.insert(
            canonical_pair_key("sp1", "sp3"),
            vec![mk_pair(
                "sp1", 0, '+', "ACGTACGT", "sp3", 10, '-', "ACGTACGT",
            )],
        );
        let merged = run_tba_tree_merge(&tree, &pairs);
        assert_eq!(merged.len(), 1);
        let blk = &merged[0];
        assert_eq!(blk.rows.len(), 3);
        assert_eq!(blk.n_cols(), 10);
        let r3 = blk.rows.iter().find(|r| r.species == "sp3").expect("sp3");
        assert_eq!(r3.strand, '-');
        assert_eq!(r3.start, 10);
        assert_eq!(r3.size, 8);
        assert_eq!(r3.seq, b"ACGT--ACGT");
    }

    #[test]
    fn from_run_alignment_matches_maf_roundtrip_on_plus_and_minus_strands() {
        use crate::formats::format_maf_blocks_runs;
        use crate::mdl::{AlignedRun, RunAlignment};

        let seq_a = b"ACGTACGTGGGGTTTTACGTACGT";
        let seq_b = b"ACGTACGTTTGGGGTTTTACGTACGT";
        let aln_fwd = RunAlignment {
            aligned_runs: vec![
                AlignedRun {
                    x_start: 0,
                    y_start: 0,
                    len: 8,
                    is_rc: false,
                },
                AlignedRun {
                    x_start: 8,
                    y_start: 10,
                    len: 16,
                    is_rc: false,
                },
            ],
            total_pairs: 24,
            lcbs_fwd: 1,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let direct_fwd = PairwiseTbaBlock::from_run_alignment(
            "spA", "spA.chr1", seq_a, "spB", "chr1", seq_b, &aln_fwd,
        );
        let maf_fwd = format_maf_blocks_runs(
            "chr1",
            seq_a.len(),
            "chr1",
            seq_b.len(),
            seq_a,
            seq_b,
            &aln_fwd,
        );
        let parsed_fwd = parse_pairwise_tba_reader("spA", "spB", maf_fwd.as_bytes()).expect("fwd");
        assert_eq!(direct_fwd.len(), parsed_fwd.len());
        assert_eq!(direct_fwd[0].row_a, parsed_fwd[0].row_a);
        assert_eq!(direct_fwd[0].row_b, parsed_fwd[0].row_b);

        // Minus-strand run alignment (for is_rc=true, y_start is the top descending coordinate).
        let seq_b_rc = crate::seeding::revcomp(seq_a);
        let aln_rev = RunAlignment {
            aligned_runs: vec![AlignedRun {
                x_start: 2,
                y_start: 21,
                len: 18,
                is_rc: true,
            }],
            total_pairs: 18,
            lcbs_fwd: 0,
            lcbs_rc: 1,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let direct_rev = PairwiseTbaBlock::from_run_alignment(
            "spA", "chr1", seq_a, "spB", "chr2", &seq_b_rc, &aln_rev,
        );
        let maf_rev = format_maf_blocks_runs(
            "chr1",
            seq_a.len(),
            "chr2",
            seq_b_rc.len(),
            seq_a,
            &seq_b_rc,
            &aln_rev,
        );
        let parsed_rev = parse_pairwise_tba_reader("spA", "spB", maf_rev.as_bytes()).expect("rev");
        assert_eq!(direct_rev.len(), parsed_rev.len());
        assert_eq!(direct_rev[0].row_a, parsed_rev[0].row_a);
        assert_eq!(direct_rev[0].row_b, parsed_rev[0].row_b);
        assert_eq!(direct_rev[0].row_b.strand, '-');
    }

    #[test]
    fn split_cross_clade_pairwise_block_scales_score_proportionally() {
        // Left clade has two 10-bp blocks covering spA [0..10) and [10..20) with score 10.0 each.
        // Cross-clade pairwise block spans spA [0..20) <-> spC [0..20) with score 200.0.
        // Splitting p across the two 10-bp sub-runs must apportion 100.0 + 100.0 = 200.0, NOT 400.0.
        let make_row = |sp: &str, start: usize, seq: &str| TbaRow {
            species: sp.to_owned(),
            contig: "chr1".to_owned(),
            start,
            size: seq.len(),
            strand: '+',
            src_size: 100,
            seq: seq.as_bytes().to_vec(),
        };
        let left = vec![
            TbaBlock {
                score: 10.0,
                rows: vec![
                    make_row("spA", 0, "ACGTACGTAC"),
                    make_row("spB", 0, "ACGTACGTAC"),
                ],
            },
            TbaBlock {
                score: 10.0,
                rows: vec![
                    make_row("spA", 10, "TGCATGCATG"),
                    make_row("spB", 10, "TGCATGCATG"),
                ],
            },
        ];
        let right: Vec<TbaBlock> = Vec::new();
        let pairs = [PairwiseTbaBlock {
            score: 200.0,
            row_a: make_row("spA", 0, "ACGTACGTACTGCATGCATG"),
            row_b: make_row("spC", 0, "ACGTACGTACTGCATGCATG"),
        }];
        let left_sp = HashSet::from(["spA", "spB"]);
        let pair_refs: Vec<&PairwiseTbaBlock> = pairs.iter().collect();
        let merged = tba_merge_two_blocksets(&left, &right, &left_sp, &pair_refs);
        assert_eq!(merged.len(), 2);
        let total_score: f64 = merged.iter().map(|b| b.score).sum();
        assert!(
            (total_score - 220.0).abs() < 1e-6,
            "expected proportional total score 10+10+200=220.0, got {total_score} ({:?})",
            merged.iter().map(|b| b.score).collect::<Vec<_>>()
        );
    }

    #[test]
    fn project_tips_onto_backbone_threads_core_accessory_inverted_and_rejects_duplicates() {
        // Hinge backbone over (h1, h2) with a shared core [0..10) and an h2-only block absent from h1.
        // Tip t1a (+ strand) and tip t1b (- strand) project onto h1 covering BOTH the core [0..10)
        // AND an h1-specific accessory locus [30..40).
        let tree = parse_newick("(h1,h2);").expect("tree");
        let mut pairs: HashMap<(String, String), Vec<PairwiseTbaBlock>> = HashMap::new();
        pairs.insert(
            canonical_pair_key("h1", "h2"),
            vec![mk_pair(
                "h1",
                0,
                '+',
                "ACGTACGTAC",
                "h2",
                0,
                '+',
                "ACGTACGTAC",
            )],
        );
        let backbone = run_tba_tree_merge(&tree, &pairs);
        assert_eq!(backbone.len(), 1);

        // Add only tip <-> h1 pairs (zero tip <-> tip or tip <-> h2 pairs).
        pairs.insert(
            canonical_pair_key("h1", "t1a"),
            vec![
                mk_pair("h1", 0, '+', "ACGTACGTAC", "t1a", 5, '+', "ACGTACGTAC"),
                mk_pair("h1", 30, '+', "GGGGTTTTAA", "t1a", 35, '+', "GGGGTTTTAA"),
            ],
        );
        pairs.insert(
            canonical_pair_key("h1", "t1b"),
            vec![
                mk_pair("h1", 0, '+', "ACGTACGTAC", "t1b", 12, '-', "ACGTACGTAC"),
                mk_pair("h1", 30, '+', "GGGGTTTTAA", "t1b", 42, '-', "GGGGTTTTAA"),
            ],
        );
        let tips = vec![
            ("t1a".to_owned(), "h1".to_owned()),
            ("t1b".to_owned(), "h1".to_owned()),
        ];
        let extended =
            tba_project_tips_onto_blockset(backbone.clone(), &tips, &pairs).expect("project tips");
        assert_eq!(
            extended.len(),
            2,
            "expected 4-way core + 3-way h1 accessory"
        );

        let core = extended
            .iter()
            .find(|b| b.rows.len() == 4)
            .expect("4-way core block");
        let t1b_core = core.rows.iter().find(|r| r.species == "t1b").expect("t1b");
        assert_eq!(t1b_core.strand, '-');
        assert_eq!(t1b_core.start, 12);

        let acc = extended
            .iter()
            .find(|b| b.rows.len() == 3)
            .expect("3-way h1+t1a+t1b accessory block");
        let acc_sps: HashSet<&str> = acc.rows.iter().map(|r| r.species.as_str()).collect();
        assert_eq!(acc_sps, HashSet::from(["h1", "t1a", "t1b"]));

        // Fail-closed: rejecting a tip whose species already exists in the backbone.
        let dup_err =
            tba_project_tips_onto_blockset(backbone, &[("h2".to_owned(), "h1".to_owned())], &pairs);
        assert!(dup_err.is_err());
    }

    #[test]
    fn parse_multi_tba_reader_roundtrips_post_hoc_extension() {
        let tree = parse_newick("(h1,h2);").expect("tree");
        let mut pairs: HashMap<(String, String), Vec<PairwiseTbaBlock>> = HashMap::new();
        pairs.insert(
            canonical_pair_key("h1", "h2"),
            vec![mk_pair(
                "h1",
                0,
                '+',
                "ACGT--ACGT",
                "h2",
                0,
                '+',
                "ACGTTTACGT",
            )],
        );
        let backbone = run_tba_tree_merge(&tree, &pairs);
        let mut maf_buf = Vec::new();
        let mut gfa_buf = Vec::new();
        write_tba_maf_and_gfa(&backbone, &mut maf_buf, &mut gfa_buf).expect("write backbone");

        let parsed = parse_multi_tba_reader(maf_buf.as_slice()).expect("parse backbone");
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].rows.len(), 2);
        assert_eq!(parsed[0].n_cols(), 10);

        pairs.insert(
            canonical_pair_key("h1", "new_tip"),
            vec![mk_pair(
                "h1", 0, '+', "ACGTACGT", "new_tip", 4, '+', "ACGTACGT",
            )],
        );
        let extended = tba_project_tips_onto_blockset(
            parsed,
            &[("new_tip".to_owned(), "h1".to_owned())],
            &pairs,
        )
        .expect("post-hoc extend");
        assert_eq!(extended.len(), 1);
        assert_eq!(extended[0].rows.len(), 3);
        let tip_row = extended[0]
            .rows
            .iter()
            .find(|r| r.species == "new_tip")
            .expect("new_tip");
        assert_eq!(tip_row.seq, b"ACGT--ACGT");
    }

    #[test]
    fn run_tba_tree_merge_lazy_loads_each_pair_exactly_once_and_matches_in_memory() {
        use std::sync::Mutex;

        let tree = parse_newick("((lav1,lav2),(vin1,vin2));").expect("tree");
        let mut pairs: HashMap<(String, String), Vec<PairwiseTbaBlock>> = HashMap::new();
        pairs.insert(
            canonical_pair_key("lav1", "lav2"),
            vec![
                mk_pair("lav1", 0, '+', "ACGTACGTAC", "lav2", 0, '+', "ACGTACGTAC"),
                mk_pair("lav1", 20, '+', "GGGGGGGGGG", "lav2", 20, '+', "GGGGGGGGGG"),
            ],
        );
        pairs.insert(
            canonical_pair_key("vin1", "vin2"),
            vec![
                mk_pair("vin1", 0, '+', "ACGTACGTAC", "vin2", 0, '+', "ACGTACGTAC"),
                mk_pair("vin1", 40, '+', "CCCCCCCCCC", "vin2", 40, '+', "CCCCCCCCCC"),
            ],
        );
        pairs.insert(
            canonical_pair_key("lav1", "vin1"),
            vec![mk_pair(
                "lav1",
                0,
                '+',
                "ACGTACGTAC",
                "vin1",
                0,
                '+',
                "ACGTACGTAC",
            )],
        );

        let expected = run_tba_tree_merge(&tree, &pairs);
        let load_counts: Mutex<HashMap<(String, String), usize>> = Mutex::new(HashMap::new());
        let lazy = run_tba_tree_merge_lazy(&tree, &|sp_a, sp_b| {
            let key = canonical_pair_key(sp_a, sp_b);
            *load_counts
                .lock()
                .expect("lock")
                .entry(key.clone())
                .or_insert(0) += 1;
            Ok(pairs.get(&key).cloned().unwrap_or_default())
        })
        .expect("lazy merge");

        assert_eq!(lazy.len(), expected.len());
        for (a, b) in lazy.iter().zip(&expected) {
            assert_eq!(a.rows.len(), b.rows.len());
            for (ra, rb) in a.rows.iter().zip(&b.rows) {
                assert_eq!(
                    (&ra.species, &ra.contig, ra.start, ra.strand, &ra.seq),
                    (&rb.species, &rb.contig, rb.start, rb.strand, &rb.seq)
                );
            }
        }
        let counts = load_counts.into_inner().expect("counts");
        assert_eq!(
            counts.len(),
            6,
            "every unordered species pair in a 4-leaf tree must be queried once"
        );
        assert!(
            counts.values().all(|&c| c == 1),
            "every species pair must be loaded at most once at its LCA node: {counts:?}"
        );
    }

    #[test]
    fn run_tba_tree_merge_consuming_drains_all_pairs_and_matches_borrowed_merge() {
        let tree = parse_newick("((lav1,lav2),(vin1,vin2));").expect("tree");
        let mut pairs: HashMap<(String, String), Vec<PairwiseTbaBlock>> = HashMap::new();
        pairs.insert(
            canonical_pair_key("lav1", "lav2"),
            vec![mk_pair(
                "lav1",
                0,
                '+',
                "ACGTACGTAC",
                "lav2",
                0,
                '+',
                "ACGTACGTAC",
            )],
        );
        pairs.insert(
            canonical_pair_key("vin1", "vin2"),
            vec![mk_pair(
                "vin1",
                0,
                '+',
                "ACGTACGTAC",
                "vin2",
                0,
                '+',
                "ACGTACGTAC",
            )],
        );
        pairs.insert(
            canonical_pair_key("lav1", "vin1"),
            vec![mk_pair(
                "lav1",
                0,
                '+',
                "ACGTACGTAC",
                "vin1",
                0,
                '+',
                "ACGTACGTAC",
            )],
        );
        let expected = run_tba_tree_merge(&tree, &pairs);
        let consumed = run_tba_tree_merge_consuming(&tree, &mut pairs);
        assert!(
            pairs.is_empty(),
            "all pairwise entries must be freed immediately after their LCA node merges"
        );
        assert_eq!(consumed.len(), expected.len());
        assert_eq!(consumed[0].rows.len(), expected[0].rows.len());
    }

    #[test]
    fn select_node_outgroups_picks_nearest_non_descendants_by_tree_distance() {
        let tree = parse_newick("(((A,B),C),(D,E));").expect("tree");
        let ab: HashSet<&str> = ["A", "B"].into_iter().collect();
        let og_ab = progressive::select_node_outgroups(&tree, &ab, 2);
        assert_eq!(og_ab.len(), 2);
        assert_eq!(
            og_ab[0], "C",
            "closest phylogenetic outgroup to (A,B) must be C (distance 3), got {og_ab:?}"
        );
        assert!(
            og_ab[1] == "D" || og_ab[1] == "E",
            "second outgroup must be D or E, got {og_ab:?}"
        );

        let de: HashSet<&str> = ["D", "E"].into_iter().collect();
        let og_de = progressive::select_node_outgroups(&tree, &de, 1);
        assert_eq!(
            og_de,
            vec!["C".to_owned()],
            "closest outgroup to (D,E) must be C (distance 4 vs A/B distance 5)"
        );

        let all: HashSet<&str> = ["A", "B", "C", "D", "E"].into_iter().collect();
        assert!(
            progressive::select_node_outgroups(&tree, &all, 2).is_empty(),
            "root clade has no outgroups"
        );
    }

    #[test]
    fn reconstruct_ancestral_contigs_polarizes_substitutions_deletions_and_insertions() {
        let mk_row = |sp: &str, start: usize, src_size: usize, seq: &str| TbaRow {
            species: sp.to_owned(),
            contig: "chr1".to_owned(),
            start,
            size: seq.bytes().filter(|&b| b != b'-').count(),
            strand: '+',
            src_size,
            seq: seq.as_bytes().to_vec(),
        };

        // Block 1: polarized deletion in sp2 ("TA" present in sp1 and og -> kept in Anc0)
        //          + polarized substitution at col 7 (sp1=A, sp2=G, og=G -> Anc0=G).
        let blk_del_sub = TbaBlock {
            score: 80.0,
            rows: vec![
                mk_row("sp1", 0, 100, "ACGTACGA"),
                mk_row("sp2", 0, 100, "ACG--CGG"),
                mk_row("og", 0, 100, "ACGTACGG"),
            ],
        };
        // Block 2: polarized insertion in sp1 ("TTTT" absent from sp2 and og -> excluded from Anc0).
        let blk_ins = TbaBlock {
            score: 60.0,
            rows: vec![
                mk_row("sp1", 10, 100, "ACGTTTTTCGT"),
                mk_row("sp2", 10, 100, "ACGT----CGT"),
                mk_row("og", 10, 100, "ACGT----CGT"),
            ],
        };
        // Block 3: whole-block deletion in sp2 (present in sp1 + og -> retained in Anc0).
        let blk_whole_del = TbaBlock {
            score: 50.0,
            rows: vec![
                mk_row("sp1", 30, 100, "GGCCAATT"),
                mk_row("og", 30, 100, "GGCCAATT"),
            ],
        };
        // Block 4: lineage-specific insertion in sp1 only (absent from sp2 and og -> excluded).
        let blk_child_only = TbaBlock {
            score: 40.0,
            rows: vec![mk_row("sp1", 50, 100, "CCCCCCCC")],
        };

        let out = progressive::reconstruct_ancestral_contigs(
            "Anc0",
            "sp1",
            "sp2",
            &["og".to_owned()],
            &[blk_del_sub, blk_ins, blk_whole_del, blk_child_only],
        );

        let concat_anc: String = out
            .anc_contigs
            .iter()
            .map(|(_, s)| String::from_utf8_lossy(s).into_owned())
            .collect();
        assert!(
            concat_anc.contains("ACGTACGG"),
            "polarized deletion + substitution must yield ACGTACGG, got {concat_anc}"
        );
        assert!(
            concat_anc.contains("ACGTCGT"),
            "polarized insertion TTTT must be stripped to ACGTCGT, got {concat_anc}"
        );
        assert!(
            concat_anc.contains("GGCCAATT"),
            "outgroup-supported deletion block must be retained in Anc0, got {concat_anc}"
        );
        assert!(
            !concat_anc.contains("CCCCCCCC"),
            "sp1-only insertion block must be excluded from Anc0, got {concat_anc}"
        );
        assert_eq!(out.stats.outgroup_polarized_deletion_bp, 10); // 2 bp in blk1 + 8 bp in blk3
        assert_eq!(out.stats.lineage_insertion_excluded_bp, 12); // 4 bp in blk2 + 8 bp in blk4
        assert!(
            !out.anc_to_left.is_empty() && !out.anc_to_right.is_empty(),
            "parent->child pairwise links must be emitted for both children"
        );
    }

    #[test]
    fn progressive_ancestral_tba_runs_end_to_end_and_builds_valid_hal() {
        let mut state: u64 = 0xDEAD_BEEF_1234_5678;
        let bases = *b"ACGT";
        let anc_seq: Vec<u8> = (0..1_800)
            .map(|_| {
                state = crate::splitmix64(state);
                bases[(state & 3) as usize]
            })
            .collect();
        let mutate = |src: &[u8], seed: u64, pct: u64| -> Vec<u8> {
            let mut st = seed;
            src.iter()
                .map(|&b| {
                    st = crate::splitmix64(st);
                    if st % 100 < pct {
                        bases[((st >> 8) & 3) as usize]
                    } else {
                        b
                    }
                })
                .collect()
        };
        let sp1 = mutate(&anc_seq, 11, 2);
        let sp2 = mutate(&anc_seq, 22, 2);
        let sp3 = mutate(&anc_seq, 33, 3);
        let sp4 = mutate(&anc_seq, 44, 4);

        let tree = parse_newick("((sp1,sp2),(sp3,sp4));").expect("tree");
        let leaf_genomes: HashMap<String, Vec<(String, Vec<u8>)>> =
            [("sp1", sp1), ("sp2", sp2), ("sp3", sp3), ("sp4", sp4)]
                .into_iter()
                .map(|(k, v)| (k.to_owned(), vec![("chr1".to_owned(), v)]))
                .collect();

        let res = progressive::run_progressive_ancestral_tba(
            &tree,
            &leaf_genomes,
            &progressive::ProgressiveAncestralConfig {
                max_outgroups: 2,
                top_k: 0,
                respect_softmask: false,
                use_transducer_ml: false,
            },
        )
        .expect("progressive ancestral tba");

        assert_eq!(
            res.ancestors.len(),
            3,
            "4-leaf binary tree has 3 internal ancestors"
        );
        assert!(
            res.merged_blocks.iter().any(|b| b.rows.len() == 4),
            "expected at least one 4-species core block"
        );
        assert_eq!(
            res.hal.genomes.len(),
            7,
            "HAL must contain 4 leaves + 3 internal ancestors"
        );
        let hal_bytes = res.hal.to_hdf5_bytes();
        assert!(crate::hal::is_hdf5_hal(&hal_bytes));
    }

    #[test]
    fn transducer_ml_ancestor_uses_branch_lengths_to_resolve_asymmetric_substitutions_and_indels() {
        let mk_row = |sp: &str, seq: &str| TbaRow {
            species: sp.to_owned(),
            contig: "chr1".to_owned(),
            start: 0,
            size: seq.bytes().filter(|&b| b != b'-').count(),
            strand: '+',
            src_size: 100,
            seq: seq.as_bytes().to_vec(),
        };
        // Short branch to sp2 (t=0.01), very long branch to sp1 (t=0.80), NO outgroup:
        // Fitch parsimony breaks A vs G tie arbitrarily toward 'A' (sp1), whereas branch-length-aware
        // Felsenstein ML MUST strongly prefer 'G' (sp2, near-zero branch length!) AND keep sp2's
        // base when sp1 has a gap on a long branch (high probability of deletion on long branch)!
        let blk = TbaBlock {
            score: 50.0,
            rows: vec![mk_row("sp1", "AAAA----AAAA"), mk_row("sp2", "GGGGTTTTGGGG")],
        };
        let branches = progressive::AncestralBranchLengths {
            left_t: 0.80,
            right_t: 0.01,
            outgroup_t: 0.30,
            outgroup_stem_t: 0.30,
            outgroup_tip_t: 0.0,
            ins_rate: 0.05,
            del_rate: 0.08,
        };
        let res = progressive::reconstruct_ancestral_contigs_transducer(
            "Anc0",
            "sp1",
            "sp2",
            &[],
            &[blk],
            branches,
        );
        let anc_seq: String = res
            .anc_contigs
            .iter()
            .map(|(_, s)| String::from_utf8_lossy(s).into_owned())
            .collect();
        assert_eq!(
            anc_seq, "GGGGTTTTGGGG",
            "branch-length-aware ML transducer must favor near-zero branch sp2 over long-branch sp1, got {anc_seq}"
        );
    }

    #[test]
    fn poa_refines_shared_conspecific_tip_insertions_into_unified_columns_and_preserves_sequences()
    {
        // Hinge h1 lacks a 64 bp insertion that 3 conspecific tips (t1, t2, t3) share with minor
        // SNPs (~3% divergence). Projecting via only h1<->tip pairs initially staggers them into
        // 3 * 64 = 192 unaligned columns; linear O(N_tips) DAG/POA refinement must collapse them
        // into ~64 aligned multi-tip columns while preserving every row's exact ungapped DNA.
        let flank_l = "ACGTACGTACGTACGTACGT";
        let flank_r = "TGCATGCATGCATGCATGCA";
        let ins_base = "GATTACAGATTACAGATTACAGATTACACTGACTGACTGACTGACTGACGTACGTTGCATGCAT";
        assert_eq!(ins_base.len(), 64);
        let mut ins_t2 = ins_base.as_bytes().to_vec();
        ins_t2[10] = b'A';
        ins_t2[45] = b'C';
        let ins_t2_str = String::from_utf8(ins_t2).expect("utf8");
        let mut ins_t3 = ins_base.as_bytes().to_vec();
        ins_t3[22] = b'G';
        let ins_t3_str = String::from_utf8(ins_t3).expect("utf8");

        let h1_gapped = format!("{flank_l}{}{flank_r}", "-".repeat(64));
        let t1_seq = format!("{flank_l}{ins_base}{flank_r}");
        let t2_seq = format!("{flank_l}{ins_t2_str}{flank_r}");
        let t3_seq = format!("{flank_l}{ins_t3_str}{flank_r}");

        let tree = parse_newick("(h1,h2);").expect("tree");
        let mut pairs: HashMap<(String, String), Vec<PairwiseTbaBlock>> = HashMap::new();
        let h_seq = format!("{flank_l}{flank_r}");
        pairs.insert(
            canonical_pair_key("h1", "h2"),
            vec![mk_pair("h1", 0, '+', &h_seq, "h2", 0, '+', &h_seq)],
        );
        let backbone = run_tba_tree_merge(&tree, &pairs);

        pairs.insert(
            canonical_pair_key("h1", "t1"),
            vec![mk_pair("h1", 0, '+', &h1_gapped, "t1", 0, '+', &t1_seq)],
        );
        pairs.insert(
            canonical_pair_key("h1", "t2"),
            vec![mk_pair("h1", 0, '+', &h1_gapped, "t2", 0, '+', &t2_seq)],
        );
        pairs.insert(
            canonical_pair_key("h1", "t3"),
            vec![mk_pair("h1", 0, '+', &h1_gapped, "t3", 0, '+', &t3_seq)],
        );

        let tips = vec![
            ("t1".to_owned(), "h1".to_owned()),
            ("t2".to_owned(), "h1".to_owned()),
            ("t3".to_owned(), "h1".to_owned()),
        ];
        let extended = tba_project_tips_onto_blockset(backbone, &tips, &pairs).expect("project");
        assert_eq!(extended.len(), 1);
        let blk = &extended[0];
        assert_eq!(blk.rows.len(), 5, "h1 + h2 + t1 + t2 + t3");

        // Staggered unaligned width would be 40 + 192 = 232 cols.
        // POA-aligned width must be close to 40 + 64 = 104 cols (<= 112).
        assert!(
            blk.n_cols() <= 112,
            "expected shared 64bp insertion aligned across t1/t2/t3 into ~104 cols, got {}",
            blk.n_cols()
        );

        // Verify ungapped DNA is 100% preserved for every single row.
        let ungapped = |sp: &str| -> String {
            let r = blk.rows.iter().find(|r| r.species == sp).expect(sp);
            assert_eq!(r.seq.len(), blk.n_cols());
            let ug: Vec<u8> = r.seq.iter().copied().filter(|&b| b != b'-').collect();
            assert_eq!(ug.len(), r.size);
            String::from_utf8(ug).expect("utf8")
        };
        assert_eq!(ungapped("h1"), h_seq);
        assert_eq!(ungapped("h2"), h_seq);
        assert_eq!(ungapped("t1"), t1_seq);
        assert_eq!(ungapped("t2"), t2_seq);
        assert_eq!(ungapped("t3"), t3_seq);

        // Verify that t1, t2, t3 actually share aligned bases inside the hinge-gap columns!
        let hinge_row = blk.rows.iter().find(|r| r.species == "h1").expect("h1");
        let tip_a = blk.rows.iter().find(|r| r.species == "t1").expect("t1");
        let tip_b = blk.rows.iter().find(|r| r.species == "t2").expect("t2");
        let tip_c = blk.rows.iter().find(|r| r.species == "t3").expect("t3");
        let shared_ins_cols = (0..blk.n_cols())
            .filter(|&c| {
                hinge_row.seq[c] == b'-'
                    && tip_a.seq[c] != b'-'
                    && tip_b.seq[c] != b'-'
                    && tip_c.seq[c] != b'-'
            })
            .count();
        assert!(
            shared_ins_cols >= 60,
            "expected >= 60 aligned 3-tip columns inside shared insertion, got {shared_ins_cols}"
        );
    }

    #[test]
    fn poa_keeps_non_homologous_tip_insertions_unaligned_under_mdl_gate() {
        let flank_l = "ACGTACGTACGTACGTACGT";
        let flank_r = "TGCATGCATGCATGCATGCA";
        // Two completely unrelated 64bp insertions (pure AT-rich vs pure GC-rich homopolymers/dinucs)
        let ins_at = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let ins_gc = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
        let h_seq = format!("{flank_l}{flank_r}");
        let h1_gapped = format!("{flank_l}{}{flank_r}", "-".repeat(64));
        let t1_seq = format!("{flank_l}{ins_at}{flank_r}");
        let t2_seq = format!("{flank_l}{ins_gc}{flank_r}");

        let tree = parse_newick("(h1,h2);").expect("tree");
        let mut pairs: HashMap<(String, String), Vec<PairwiseTbaBlock>> = HashMap::new();
        pairs.insert(
            canonical_pair_key("h1", "h2"),
            vec![mk_pair("h1", 0, '+', &h_seq, "h2", 0, '+', &h_seq)],
        );
        let backbone = run_tba_tree_merge(&tree, &pairs);
        pairs.insert(
            canonical_pair_key("h1", "t1"),
            vec![mk_pair("h1", 0, '+', &h1_gapped, "t1", 0, '+', &t1_seq)],
        );
        pairs.insert(
            canonical_pair_key("h1", "t2"),
            vec![mk_pair("h1", 0, '+', &h1_gapped, "t2", 0, '+', &t2_seq)],
        );
        let tips = vec![
            ("t1".to_owned(), "h1".to_owned()),
            ("t2".to_owned(), "h1".to_owned()),
        ];
        let extended = tba_project_tips_onto_blockset(backbone, &tips, &pairs).expect("project");
        let blk = &extended[0];
        let hinge_row = blk.rows.iter().find(|r| r.species == "h1").expect("h1");
        let tip_a = blk.rows.iter().find(|r| r.species == "t1").expect("t1");
        let tip_b = blk.rows.iter().find(|r| r.species == "t2").expect("t2");
        let forced_cols = (0..blk.n_cols())
            .filter(|&c| hinge_row.seq[c] == b'-' && tip_a.seq[c] != b'-' && tip_b.seq[c] != b'-')
            .count();
        assert_eq!(
            forced_cols, 0,
            "non-homologous AT vs GC tip insertions must remain unaligned under MDL gate, got {forced_cols} forced cols"
        );
    }
}
