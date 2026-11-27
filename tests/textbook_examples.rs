//! R9-RUST-VERIFY: Live Rust integration tests verifying worked examples and numerical
//! exercises in "Mathematics and Algorithms of the Horizon Suite" directly against
//! `horizon_aligner::*` (complementing `scripts/figures/verify_textbook_examples.py`).

use horizon_aligner::bocpd::{IDENTITY_PRIOR, base_probs, bocpd_gc_profile, joint_pair_bits};
use horizon_aligner::coord::Interval;
use horizon_aligner::egraph::{
    EGraph, ENode, Expr, MdlCostModel, build_seq_expr, collect_dups, extract_min_mdl_cost,
    saturate_synteny,
};
use horizon_aligner::hal::{
    HalAlignment, read_hdf5_superblock, read_hdf5_v0_dataset_bytes_seek,
    read_hdf5_v0_group_entries_seek,
};
use horizon_aligner::mdl::{
    DEFAULT_CRC_ALPHA, DEFAULT_CRC_N_CAL, GraphNode, HmmParams, align_seq_to_dag_mdl,
    banded_mdl_gap, banded_mdl_gap_with_min_post, calibrate_crc_threshold,
};
use horizon_aligner::merge::merge_pairwise_mafs;
use horizon_aligner::seeding::{
    MASK_PRIMARY_11_15, canonical_spaced_seeds, lis_chain, prepare_contig_index, revcomp,
    sdust_mask, sketch_jaccard,
};
use horizon_aligner::tba::{GuideTree, build_upgma_tree};
use std::collections::HashMap;
use std::io::Cursor;

#[test]
fn ch1_coordinates_and_composition_aware_scoring() {
    // Worked Example 1.1: [10, 35) on 100 bp contig
    let fwd = Interval::new(10, 35, '+', 100);
    let rev = Interval::new(10, 35, '-', 100);
    assert_eq!(fwd.span(), 25);
    assert_eq!(fwd.ucsc_start(), 10);
    assert_eq!(rev.ucsc_start(), 65);
    // Exact involution: 100 - (65 + 25) = 10
    assert_eq!(100 - (rev.ucsc_start() + rev.span()), 10);

    // Exercise 1.1: [420, 600) on 2,500 bp contig, base k = 17
    let ex1 = Interval::new(420, 600, '-', 2500);
    assert_eq!(ex1.ucsc_start(), 1900);
    assert_eq!(ex1.span(), 180);
    assert_eq!(2500 - 1 - (ex1.ucsc_start() + 17), 582);

    // Worked Example 1.2: 19% GC expected null bit-savings
    let pi = base_probs(0.19);
    let aa = joint_pair_bits(b'A', b'A', pi, pi, IDENTITY_PRIOR);
    let gg = joint_pair_bits(b'G', b'G', pi, pi, IDENTITY_PRIOR);
    let mm = joint_pair_bits(b'A', b'C', pi, pi, IDENTITY_PRIOR);
    assert!((aa - 0.967).abs() < 0.005, "aa={aa}");
    assert!((gg - 2.846).abs() < 0.005, "gg={gg}");
    assert!((mm - (-1.515)).abs() < 0.005, "mm={mm}");
    let p_at = (1.0 - 0.19) / 2.0;
    let p_gc = 0.19 / 2.0;
    let exp_null = 2.0 * p_at * p_at * aa
        + 2.0 * p_gc * p_gc * gg
        + (1.0 - 2.0 * p_at * p_at - 2.0 * p_gc * p_gc) * mm;
    assert!((exp_null - (-0.622)).abs() < 0.005, "exp_null={exp_null}");

    // Exercise 1.2: Gibbs' inequality strictly negative across non-degenerate GC in [0.02, 0.98]
    for gc_pct in [2, 10, 19, 41, 50, 75, 98] {
        let g = f64::from(gc_pct) / 100.0;
        let p = base_probs(g);
        let bases = *b"ACGT";
        let mut e = 0.0;
        for (i, &a) in bases.iter().enumerate() {
            for (j, &b) in bases.iter().enumerate() {
                e += p[i] * p[j] * joint_pair_bits(a, b, p, p, IDENTITY_PRIOR);
            }
        }
        assert!(
            e < 0.0,
            "expected strictly negative null expectation at gc={g}, got {e}"
        );
    }

    // Exercise 1.3: GC crossover at g = 0.50 -> 1.5607 bits
    let p50 = base_probs(0.50);
    let aa50 = joint_pair_bits(b'A', b'A', p50, p50, IDENTITY_PRIOR);
    let gg50 = joint_pair_bits(b'G', b'G', p50, p50, IDENTITY_PRIOR);
    assert!((aa50 - gg50).abs() < 1e-12);
    assert!((aa50 - 1.5607).abs() < 0.0001, "aa50={aa50}");
}

#[test]
fn ch2_sdust_spaced_seed_symmetry_and_lis_chain() {
    // Exercise 2.1: 16 bp dinucleotide repeat ATATAT... is masked by SDUST
    let mask = sdust_mask(b"ATATATATATATATAT");
    assert_eq!(mask.len(), 16);
    assert!(mask.iter().all(|&m| m), "expected all 16 bp masked");

    // Chapter 2: reverse-complementation preserves canonical spaced-seed keys, plain reversal does not
    let seq = b"ACGTTGCATAGGCTA";
    let rc_seq = revcomp(seq);
    let rev_only: Vec<u8> = seq.iter().rev().copied().collect();
    let k_fwd: Vec<u64> = canonical_spaced_seeds(seq, MASK_PRIMARY_11_15, 4)
        .expect("fwd")
        .into_keys()
        .collect();
    let k_rc: Vec<u64> = canonical_spaced_seeds(&rc_seq, MASK_PRIMARY_11_15, 4)
        .expect("rc")
        .into_keys()
        .collect();
    let k_rev: Vec<u64> = canonical_spaced_seeds(&rev_only, MASK_PRIMARY_11_15, 4)
        .expect("rev")
        .into_keys()
        .collect();
    assert_eq!(
        k_fwd, k_rc,
        "reverse-complement must preserve canonical key"
    );
    assert_ne!(k_fwd, k_rev, "plain reversal must not swap fwd and rc");

    // Worked Example 2.1: patience-sorting LIS trace
    let pts = vec![(10, 12), (10, 20), (18, 15), (25, 22), (30, 18), (38, 28)];
    let chain = lis_chain(&pts);
    assert_eq!(chain, vec![(10, 12), (18, 15), (30, 18), (38, 28)]);

    // Exercise 2.3: descending-y tie breaker prevents duplicate x = 10 in chain
    let ex23 = lis_chain(&[(10, 12), (10, 20), (15, 25)]);
    assert_eq!(ex23.len(), 2);
    assert!(ex23.windows(2).all(|w| w[0].0 < w[1].0 && w[0].1 < w[1].1));
}

#[test]
fn ch3_bocpd_worked_example_and_n_handling() {
    // Worked Example 3.1: h = 0.1 on "GAA" (exact Adams & MacKay clean (1, 1) reset)
    let prof = bocpd_gc_profile(b"GAA", 0.1, 80);
    assert!((prof[0] - 0.50).abs() < 1e-12, "g1={}", prof[0]);
    assert!((prof[1] - 0.65).abs() < 1e-12, "g2={}", prof[1]);
    assert!((prof[2] - (67.0 / 140.0)).abs() < 1e-12, "g3={}", prof[2]);

    // Exercise 3.1: 500 bp N scaffold gap preserves GC estimate without collapsing to 0.02
    let mut seq = b"GCGCGCGCGCGCGCGCGCGC".to_vec();
    seq.extend(std::iter::repeat_n(b'N', 500));
    let n_prof = bocpd_gc_profile(&seq, 1.0 / 120.0, 80);
    assert!(
        n_prof.last().copied().unwrap_or(0.0) > 0.70,
        "N gap must not drag GC estimate down to 0.02 clamp: {:?}",
        n_prof.last()
    );
}

#[test]
fn ch4_pair_hmm_worked_example_4_1() {
    let prefix = b"AATTGCATATCG";
    let pref_res = banded_mdl_gap(prefix, prefix, 180, 10.0);
    assert!(
        (pref_res.score_bits - 17.8527).abs() < 0.01,
        "prefix score_bits={}",
        pref_res.score_bits
    );
    assert_eq!(pref_res.kept_pairs.len(), 12);
    assert_eq!(pref_res.bubbles, 0);

    // Isolated 6 bp divergent gap: score_bits = -9.4230, abstained by 10.0-bit floor (NOT <0.15 posterior)
    let gap_abstain = banded_mdl_gap(b"GGGGCC", b"CCTTTT", 180, 10.0);
    assert!(
        (gap_abstain.score_bits - (-9.4230)).abs() < 0.005,
        "gap score_bits={}",
        gap_abstain.score_bits
    );
    assert_eq!(gap_abstain.kept_pairs, [] as [(usize, usize); 0]);
    assert_eq!(gap_abstain.bubbles, 1);

    // Lowering the bit floor to -100.0 while keeping min_post_prob = 0.15 retains all 6 diagonal pairs
    let gap_post_only = banded_mdl_gap_with_min_post(
        b"GGGGCC",
        b"CCTTTT",
        180,
        -100.0,
        &HmmParams::default(),
        0.15,
    );
    assert_eq!(gap_post_only.kept_pairs.len(), 6);
    assert_eq!(gap_post_only.bubbles, 0);

    // Full unanchored 22 bp DP retains all 22 pairs with positive total score_bits ~ +15.34
    let full = banded_mdl_gap(
        b"AATTGCATATCGGGGGCCATCG",
        b"AATTGCATATCGCCTTTTATCG",
        180,
        10.0,
    );
    assert!(
        (full.score_bits - 15.3365).abs() < 0.01,
        "full={}",
        full.score_bits
    );
    assert_eq!(full.kept_pairs.len(), 22);
    assert_eq!(full.bubbles, 0);

    // Chapter 4 fast-path eligibility: 16 bp poly-A never bypasses BOCPD+MDL abstention
    let poly_a = [b'A'; 16];
    let poly_res = banded_mdl_gap(&poly_a, &poly_a, 16, 20.0);
    assert!(
        poly_res.kept_pairs.is_empty(),
        "low-complexity poly-A must route through BOCPD+MDL abstention"
    );
}

#[test]
fn ch5_crc_worked_example_and_exercises() {
    // Worked Example 5.1: 15 sorted decoy scores, alpha = 0.07, floor = 10.0 -> 12.6 + 1e-6
    let scores = [
        1.2, 2.5, 3.1, 4.0, 4.5, 5.2, 6.0, 6.8, 7.1, 8.0, 8.5, 9.2, 9.8, 11.4, 12.6,
    ];
    let lam = calibrate_crc_threshold(&scores, 0.07, 10.0);
    assert!((lam - 12.600_001).abs() < 1e-9, "lam={lam}");

    // Exercise 5.1: n_cal = 10, alpha = 0.05 -> exits loop and returns max(scores).max(floor) + 1e-6
    let small = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.5];
    let lam_small = calibrate_crc_threshold(&small, 0.05, 10.0);
    assert!(
        (lam_small - 11.500_001).abs() < 1e-9,
        "lam_small={lam_small}"
    );

    // Production default n_cal = 49, alpha = 0.07 selects 3rd-largest order statistic
    let prod_scores: Vec<f64> = (1..=u32::try_from(DEFAULT_CRC_N_CAL).expect("u32"))
        .map(f64::from)
        .collect();
    let lam_prod = calibrate_crc_threshold(&prod_scores, DEFAULT_CRC_ALPHA, 10.0);
    assert!((lam_prod - 47.000_001).abs() < 1e-9, "lam_prod={lam_prod}");
}

#[test]
fn ch6_egraph_saturation_and_mdl_cost_model() {
    let add_blk = |eg: &mut EGraph, name: &str| {
        eg.add(ENode {
            op: format!("blk:{name}"),
            children: vec![],
        })
    };
    let add_inv = |eg: &mut EGraph, x: usize| {
        eg.add(ENode {
            op: "inv".to_owned(),
            children: vec![x],
        })
    };
    let add_seq = |eg: &mut EGraph, x: usize, y: usize| {
        eg.add(ENode {
            op: "seq".to_owned(),
            children: vec![x, y],
        })
    };

    // Figure 6.1: seq(inv(inv(A)), A) saturates to dup(A)
    let mut eg = EGraph::default();
    let a = add_blk(&mut eg, "A");
    let inv_a = add_inv(&mut eg, a);
    let inv_inv_a = add_inv(&mut eg, inv_a);
    let root = add_seq(&mut eg, inv_inv_a, a);
    saturate_synteny(&mut eg, 10, 512);
    let mut model = MdlCostModel {
        num_blocks: 16,
        block_stats: HashMap::from([("A".to_owned(), (500, 0.05))]),
        ..MdlCostModel::default()
    };
    let (cost, expr) = extract_min_mdl_cost(&eg, root, &model).expect("extract");
    assert!(
        matches!(expr, Expr::Dup(_)),
        "expected Dup(A), got {expr:?}"
    );
    // Worked Example 6.1: K=16, L=500, d=0.05 -> 1189.3 bits
    assert!((cost - 1189.32).abs() < 0.05, "cost={cost}");

    // Exercise 6.1: seq(inv(seq(inv(A), B)), inv(B)) produces NO dup
    let mut eg61 = EGraph::default();
    let ba = add_blk(&mut eg61, "A");
    let bb = add_blk(&mut eg61, "B");
    let ia = add_inv(&mut eg61, ba);
    let ib = add_inv(&mut eg61, bb);
    let inner = add_seq(&mut eg61, ia, bb);
    let inv_inner = add_inv(&mut eg61, inner);
    let root61 = add_seq(&mut eg61, inv_inner, ib);
    saturate_synteny(&mut eg61, 10, 512);
    let (_, expr61) =
        extract_min_mdl_cost(&eg61, root61, &MdlCostModel::from_egraph(&eg61)).expect("extract61");
    assert!(
        collect_dups(&expr61).is_empty(),
        "expected no dup, got {expr61:?}"
    );

    // Exercise 6.2: L = 100 crossover for K = 2 near d* = 0.7134, and Dup retains margin at K = 3, 4
    let mut eg_aa = EGraph::default();
    let root_aa = build_seq_expr(&mut eg_aa, &[("A", 1), ("A", 1)]).expect("seq(A,A)");
    saturate_synteny(&mut eg_aa, 10, 512);

    model.num_blocks = 2;
    model.block_stats = HashMap::from([("A".to_owned(), (100, 0.7130))]);
    let (_, e_below) = extract_min_mdl_cost(&eg_aa, root_aa, &model).expect("k2 below");
    assert!(matches!(e_below, Expr::Dup(_)));

    model.block_stats = HashMap::from([("A".to_owned(), (100, 0.7138))]);
    let (_, e_above) = extract_min_mdl_cost(&eg_aa, root_aa, &model).expect("k2 above");
    assert!(matches!(e_above, Expr::Seq(_, _)));

    for k in [3usize, 4usize] {
        model.num_blocks = k;
        model.block_stats = HashMap::from([("A".to_owned(), (100, 0.75))]);
        let (_, e_sat) = extract_min_mdl_cost(&eg_aa, root_aa, &model).expect("k>=3 sat");
        assert!(
            matches!(e_sat, Expr::Dup(_)),
            "expected Dup at K={k}, d=0.75"
        );
    }
}

#[test]
fn ch7_column_threading_minhash_and_upgma() {
    // Worked Example 7.1: Ref=ACGTA, Sp1=ACGTTTA (ins TT after G), Sp2=ACACTA (ins C after A)
    let maf_sp1 =
        "##maf version=1\n\na score=10.0\ns chr1 0 5 + 100 ACG--TA\ns s1   0 7 + 100 ACGTTTA\n\n";
    let maf_sp2 =
        "##maf version=1\n\na score=10.0\ns chr1 0 5 + 100 ACG-TA\ns s2   0 6 + 100 ACACTA\n\n";
    let (merged_maf, _, total_blocks, multi_blocks) = merge_pairwise_mafs(
        "Ref",
        &[
            ("Sp1".to_owned(), maf_sp1.to_owned()),
            ("Sp2".to_owned(), maf_sp2.to_owned()),
        ],
        1000,
    );
    assert_eq!(total_blocks, 1);
    assert_eq!(multi_blocks, 1);
    assert!(
        merged_maf.contains("ACG--TA"),
        "missing Ref row:\n{merged_maf}"
    );
    assert!(
        merged_maf.contains("ACGTTTA"),
        "missing Sp1 row:\n{merged_maf}"
    );
    assert!(
        merged_maf.contains("ACAC-TA"),
        "missing Sp2 row:\n{merged_maf}"
    );

    // Exercise 7.1: Broder's unbiased sketch_jaccard on identical & 10:1 unequal nested sets
    let mut state: u64 = 0x517c_c1b7_2722_0a95;
    let mut rand_dna = |len: usize| -> Vec<u8> {
        const BASES: [u8; 4] = *b"ACGT";
        (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                BASES[(state & 3) as usize]
            })
            .collect()
    };
    let seq_small = rand_dna(5_000);
    let mut seq_large = seq_small.clone();
    seq_large.extend_from_slice(&rand_dna(45_000));
    let sk_small = prepare_contig_index(&seq_small)
        .expect("idx small")
        .minhash_sketch(512);
    let sk_large = prepare_contig_index(&seq_large)
        .expect("idx large")
        .minhash_sketch(512);
    assert!(!sk_small.is_empty() && !sk_large.is_empty());
    assert!((sketch_jaccard(&sk_small, &sk_small) - 1.0).abs() < 1e-12);
    let j_unequal = sketch_jaccard(&sk_small, &sk_large);
    assert!(
        (j_unequal - 0.10).abs() < 0.04,
        "expected unbiased Jaccard ~0.10 on 5 kb inside 50 kb, got {j_unequal:.4}"
    );

    // Exercise 7.2: UPGMA hand merge on A, B, C, D distances -> ((A,B),(C,D))
    let names = vec![
        "A".to_owned(),
        "B".to_owned(),
        "C".to_owned(),
        "D".to_owned(),
    ];
    let dist = [
        [0.00, 0.10, 0.40, 0.42],
        [0.10, 0.00, 0.44, 0.46],
        [0.40, 0.44, 0.00, 0.20],
        [0.42, 0.46, 0.20, 0.00],
    ];
    let sim: Vec<Vec<f64>> = dist
        .iter()
        .map(|row| row.iter().map(|&d| 1.0 - d).collect())
        .collect();
    let tree = build_upgma_tree(&names, &sim).expect("upgma");
    match tree {
        GuideTree::Internal { left, right } => {
            let mut l = left.leaves();
            let mut r = right.leaves();
            l.sort_unstable();
            r.sort_unstable();
            assert!(
                (l == vec!["A", "B"] && r == vec!["C", "D"])
                    || (l == vec!["C", "D"] && r == vec!["A", "B"]),
                "unexpected UPGMA split: {l:?} vs {r:?}"
            );
        }
        GuideTree::Leaf(_) => panic!("expected internal root"),
    }
}

struct CountingSeek<R> {
    inner: R,
    seeks: usize,
}
impl<R: std::io::Read> std::io::Read for CountingSeek<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.inner.read(buf)
    }
}
impl<R: std::io::Seek> std::io::Seek for CountingSeek<R> {
    fn seek(&mut self, pos: std::io::SeekFrom) -> std::io::Result<u64> {
        self.seeks += 1;
        self.inner.seek(pos)
    }
}

#[test]
fn ch8_dag_topological_invariant_and_hal_layout() {
    // Exercise 8.1: cyclic / non-topological predecessor p >= u returns empty alignment immediately
    let cyclic_nodes = vec![
        GraphNode {
            name: "u0".to_owned(),
            seq: b"ACGTACGTACGT".to_vec(),
            preds: vec![1], // p = 1 >= u = 0 violates topological order
        },
        GraphNode {
            name: "u1".to_owned(),
            seq: b"ACGTACGTACGT".to_vec(),
            preds: vec![0],
        },
    ];
    let aln = align_seq_to_dag_mdl(&cyclic_nodes, b"ACGTACGTACGT", &HmmParams::default(), 10.0);
    assert_eq!(aln.kept_pairs, [] as [(usize, usize, usize); 0]);

    // Chapter 8 & Exercise 8.3: HAL 2-level /<sp>/DNA_ARRAY layout, empty /Genomes,
    // 33-byte TOP_ARRAY, (24 + 9*C)-byte BOTTOM_ARRAY, +1 libHal sentinel records, and >=20 seeks
    let hal = HalAlignment::from_star_blocks(
        "Ref",
        b"ACGTACGTACGTACGTACGT",
        &[("Qry".to_owned(), b"ACGTACGTACGTACGTACGT".to_vec())],
        &[(0, 0, 0, 20, false)],
    );
    let bytes = hal.to_hdf5_bytes();
    let sb = read_hdf5_superblock(&bytes).expect("superblock");

    let mut cur = CountingSeek {
        inner: Cursor::new(bytes.as_slice()),
        seeks: 1, // +1 for initial superblock read at offset 0
    };
    let root =
        read_hdf5_v0_group_entries_seek(&mut cur, sb.root_group_addr.expect("root_group_addr"))
            .expect("root io")
            .expect("root entries");
    let find_addr = |entries: &[horizon_aligner::hal::Hdf5GroupEntry], name: &str| {
        entries.iter().find(|e| e.name == name).map_or_else(
            || panic!("missing HDF5 entry {name}"),
            |e| e.object_header_addr,
        )
    };
    let genomes_addr = find_addr(&root, "Genomes");
    let ref_addr = find_addr(&root, "Ref");
    let qry_addr = find_addr(&root, "Qry");
    let _ver_addr = find_addr(&root, "Verison");

    let ref_entries = read_hdf5_v0_group_entries_seek(&mut cur, ref_addr)
        .expect("ref io")
        .expect("ref entries");
    let dna_addr = find_addr(&ref_entries, "DNA_ARRAY");
    let dna = read_hdf5_v0_dataset_bytes_seek(&mut cur, dna_addr)
        .expect("dna io")
        .expect("dna bytes");
    assert_eq!(dna.len(), 20 / 2, "4 bits/base -> 10 bytes for 20 bp");
    assert!(
        cur.seeks >= 20,
        "2-level contiguous /Ref/DNA_ARRAY traversal requires >= 20 seeks, got {}",
        cur.seeks
    );

    let genomes_group = read_hdf5_v0_group_entries_seek(&mut cur, genomes_addr)
        .expect("genomes io")
        .expect("genomes entries");
    assert!(
        genomes_group.is_empty(),
        "/Genomes is an empty placeholder group"
    );

    let qry_entries = read_hdf5_v0_group_entries_seek(&mut cur, qry_addr)
        .expect("qry io")
        .expect("qry entries");
    let top_q = read_hdf5_v0_dataset_bytes_seek(&mut cur, find_addr(&qry_entries, "TOP_ARRAY"))
        .expect("top io")
        .expect("top bytes");
    // 1 active top segment + 1 mandatory libHal sentinel = 2 * 33 = 66 bytes
    assert_eq!(top_q.len(), 2 * 33);

    let bot_r = read_hdf5_v0_dataset_bytes_seek(&mut cur, find_addr(&ref_entries, "BOTTOM_ARRAY"))
        .expect("bot io")
        .expect("bot bytes");
    // C = 1 child -> 24 + 9*1 = 33 bytes/record; 1 active + 1 sentinel = 66 bytes
    assert_eq!(bot_r.len(), 2 * (24 + 9));
}
