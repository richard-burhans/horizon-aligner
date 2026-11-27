use crate::bocpd::{GAP_EXTEND_BITS, GAP_OPEN_BITS};

#[test]
fn pair_hmm_transitions_are_normalized() {
    let p_open = 2.0_f64.powf(GAP_OPEN_BITS);
    let stay_match = 1.0 - 2.0 * p_open;
    assert!((stay_match + 2.0 * p_open - 1.0).abs() < 1e-12);

    let p_ext = 2.0_f64.powf(GAP_EXTEND_BITS);
    let leave_indel = 1.0 - p_ext;
    assert!((leave_indel + p_ext - 1.0).abs() < 1e-12);
}

#[test]
fn min_post_prob_filters_uncertain_pairs() {
    let x = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
    let y = b"ACGTACGTAAAACCCCGGGGTTTTACGTACGT";
    let params = super::HmmParams::default();
    let all = super::banded_mdl_gap_with_min_post(x, y, 8, 5.0, &params, 0.0);
    let strict = super::banded_mdl_gap_with_min_post(x, y, 8, 5.0, &params, 0.9999);
    assert_ne!(all.kept_pairs.len(), 0);
    assert!(strict.kept_pairs.len() < all.kept_pairs.len());
}

#[test]
fn empty_banded_gap_behaves_cleanly() {
    let empty_res = super::banded_mdl_gap(b"", b"ACGT", 8, 5.0);
    assert_eq!(empty_res.kept_pairs, [] as [(usize, usize); 0]);
    assert_eq!(empty_res.bubbles, 0);
}

#[test]
fn dip_bridging_does_not_bridge_across_gap_when_min_post_prob_zero() {
    let mut abstain = Vec::new();
    let cols2 = vec![
        (Some(0usize), Some(0usize), 2.5, 0.9),
        (None, Some(1), -2.0, 0.0), // gap column (must NOT be bridged into)
        (Some(1), Some(2), 2.5, 0.9),
    ];
    let (kept2, _) =
        super::filter_mea_pairs_and_count_bubbles(&cols2, &mut abstain, 0.0, 2.0, -2.5);
    assert_eq!(kept2, vec![(0, 0), (1, 2)]);
    assert!(abstain[1], "gap column must remain abstained");
}

#[test]
fn dip_floor_bits_scales_monotonically_with_identity() {
    let low = super::dip_floor_bits_for_identity(0.55);
    let mid = super::dip_floor_bits_for_identity(0.65);
    let high = super::dip_floor_bits_for_identity(0.90);
    assert!(low.is_finite() && mid.is_finite() && high.is_finite());
    assert!(low < 0.0 && mid < 0.0 && high < 0.0);
    assert!(
        low > mid && mid > high,
        "higher identity prior should penalize single-base mismatches more strongly: {low} > {mid} > {high}"
    );
}

#[test]
fn diagonal_fast_path_recovers_near_identical_gap() {
    let x = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let mut y = *x;
    y[10] = b'T';
    let res = super::banded_mdl_gap(x, &y, 16, 10.0);
    assert_eq!(res.bubbles, 0);
    assert_eq!(res.kept_pairs.len(), x.len());
    assert_eq!(res.kept_pairs[0], (0, 0));
    assert_eq!(res.kept_pairs[x.len() - 1], (x.len() - 1, x.len() - 1));
}

#[test]
fn prefix_suffix_trimming_handles_internal_indel() {
    let mut x = b"ACGTACGTACGTACGTACGTACGT".to_vec();
    x.extend_from_slice(b"GGG");
    x.extend_from_slice(b"TGCATGCATGCATGCATGCATGCA");
    let mut y = b"ACGTACGTACGTACGTACGTACGT".to_vec();
    y.extend_from_slice(b"TGCATGCATGCATGCATGCATGCA");
    let res = super::banded_mdl_gap(&x, &y, 16, 10.0);
    assert_eq!(res.kept_pairs.len(), y.len());
    assert_eq!(res.kept_pairs[0], (0, 0));
    assert_eq!(res.kept_pairs[y.len() - 1], (x.len() - 1, y.len() - 1));
}

#[test]
fn banded_viterbi_fast_path_recovers_small_indel_with_flank_snps() {
    let mut x = b"AACTGACTGACTGACTGACTGACTGACTGACTG".to_vec();
    x.extend_from_slice(b"TTT");
    x.extend_from_slice(b"GATTACAGATTACAGATTACAGATTACAA");
    let mut y = b"ATCTGACTGACTGACTGACTGACTGACTGACTG".to_vec();
    y.extend_from_slice(b"GATTACAGATTACAGATTACAGATTACAT");
    let res = super::banded_mdl_gap(&x, &y, 16, 10.0);
    assert_eq!(res.bubbles, 0);
    assert_eq!(res.kept_pairs.len(), y.len());
}

#[test]
fn posterior_trace_emits_live_grid_and_horizon_bubble() {
    let x = b"GCTAGCATCGAATTGGCCCCCAAAACCCCAAACGATGCATTCGGAATT";
    let y = b"GCTAGCATCGAATTGGTTTTGGGGTTTTGGGGCGATGCATTCGGAATT";
    let params = super::HmmParams::default();
    let trace = super::banded_mdl_gap_posterior_trace(x, y, 8, 10.0, &params, 0.15);
    assert_eq!(trace.n, 48);
    assert_eq!(trace.m, 48);
    assert_eq!(trace.band_w, 8);
    assert_eq!(trace.bubbles, 1);
    assert!(trace.score_bits > 20.0);
    assert_ne!(trace.posterior_cells.len(), 0);
    assert!(
        trace
            .posterior_cells
            .iter()
            .all(|&(_, _, p)| (0.0..=1.0).contains(&p))
    );
    assert!(
        trace.mea_cols.iter().any(|c| c.4),
        "expected abstained Horizon Bubble columns"
    );
    assert!(
        trace.mea_cols.iter().any(|c| !c.4),
        "expected kept flank columns"
    );
}

#[test]
fn low_complexity_poly_a_and_dinuc_gaps_never_bypass_bocpd_abstention() {
    // C2 regression: 16 bp poly-A yields ~0.745 bits/bp (~11.9 bits) under BOCPD, so it must
    // abstain at min_block_bits = 20.0 rather than hitting diagonal/Viterbi fast paths.
    let poly_a = [b'A'; 16];
    let res_a = super::banded_mdl_gap(&poly_a, &poly_a, 16, 20.0);
    assert!(
        res_a.kept_pairs.is_empty(),
        "16 bp poly-A bypassed BOCPD: {res_a:?}"
    );

    let at_dinuc = b"AT".repeat(20);
    let res_at = super::banded_mdl_gap(&at_dinuc, &at_dinuc, 16, 45.0);
    assert!(
        res_at.kept_pairs.is_empty(),
        "40 bp (AT)_n bypassed BOCPD: {res_at:?}"
    );

    let mut poly_a_ins = vec![b'A'; 38];
    poly_a_ins.extend_from_slice(b"AA");
    let res_ins = super::banded_mdl_gap(&poly_a_ins, &[b'A'; 38], 16, 50.0);
    assert!(
        res_ins.kept_pairs.is_empty(),
        "poly-A indel bypassed BOCPD: {res_ins:?}"
    );
}

#[test]
fn backward_scaling_does_not_overflow_or_emit_nan_posteriors() {
    // C3 regression: wide band + steep gap_extend_bits + shifted homologous blocks must keep
    // every posterior finite in [0, 1] (never turning 0.0 * inf = NaN into 1.0).
    let mut x = b"ACGT".repeat(125);
    x.extend_from_slice(&b"CGTA".repeat(250));
    let mut y = b"TGCA".repeat(250);
    y.extend_from_slice(&b"ACGT".repeat(125));
    let steep = super::HmmParams {
        identity_prior: 0.70,
        gap_open_bits: -6.0,
        gap_extend_bits: -4.0,
    };
    let trace = super::banded_mdl_gap_posterior_trace(&x, &y, 250, 10.0, &steep, 0.15);
    assert_ne!(trace.posterior_cells.len(), 0);
    assert!(
        trace
            .posterior_cells
            .iter()
            .all(|&(_, _, p)| p.is_finite() && (0.0..=1.0).contains(&p)),
        "found non-finite or out-of-range posterior"
    );
}

#[test]
fn viterbi_fast_path_abstains_localized_nonhomologous_bubble_and_respects_crc_threshold() {
    // C1 regression part 1: 450 bp identical + 100 bp completely unrelated (0% identity) + 450 bp identical.
    // Global identity is 90% >= 88%, but the 100 bp interior is a true Horizon Bubble and MUST abstain.
    let flank = b"ACGTACGTACGTTGCATGCATGCATGCAACGTACGTACGTTGCATGCATG".repeat(9); // 450 bp balanced
    assert_eq!(flank.len(), 450);
    let mid_x = b"ACGTTGCACTAGTAGC".repeat(6); // 96 bp
    let mid_y: Vec<u8> = mid_x
        .iter()
        .map(|&b| match b {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        })
        .collect(); // 96 bp, 0% identity to mid_x
    assert!(mid_x.iter().zip(&mid_y).all(|(a, b)| a != b));
    let mut x = Vec::with_capacity(996);
    x.extend_from_slice(&flank);
    x.extend_from_slice(&mid_x);
    x.extend_from_slice(&flank);
    let mut y = Vec::with_capacity(996);
    y.extend_from_slice(&flank);
    y.extend_from_slice(&mid_y);
    y.extend_from_slice(&flank);

    let res = super::banded_mdl_gap(&x, &y, 180, 10.0);
    assert!(
        res.bubbles >= 1,
        "localized 96 bp 0%-identity bubble must not bypass abstention: bubbles={}",
        res.bubbles
    );
    let mid_emitted = res
        .kept_pairs
        .iter()
        .filter(|&&(i, _)| (460..536).contains(&i))
        .count();
    assert!(
        mid_emitted < 20,
        "non-homologous 96 bp bubble interior should abstain, but emitted {mid_emitted} pairs"
    );

    // C1 regression part 2: a 65%-GC 32 bp identical gap (eligible for fast paths) must abstain
    // whenever `min_block_bits` exceeds exact `banded_mdl_gap_core`'s `max_block_bits` (even though
    // old hardcoded Viterbi `+1.8` scored 57.6 - 4.0 = 53.6 bits).
    let short_gc = b"CGCGCGCGCGCGCGCGCGCGAATTAAGCGCCG"; // 32 bp, 62.5% GC, fast-path eligible
    let mut scratch = super::DpScratch::default();
    let (_, exact_max_bits) = super::banded_mdl_gap_core_with_max_bits(
        short_gc,
        short_gc,
        16,
        0.0,
        &super::HmmParams::default(),
        0.15,
        &mut scratch,
    );
    let threshold = exact_max_bits + 0.5;
    assert!(
        threshold < 53.0,
        "sanity: exact max_block_bits ({exact_max_bits:.2}) must be below old Viterbi's 53.6"
    );
    let res_short = super::banded_mdl_gap(short_gc, short_gc, 16, threshold);
    assert!(
        res_short.kept_pairs.is_empty(),
        "fast path emitted {} pairs at threshold={threshold:.2} > exact_max_bits={exact_max_bits:.2}",
        res_short.kept_pairs.len()
    );
}

#[test]
fn prefix_suffix_trimming_does_not_emit_disconnected_subthreshold_island() {
    // C2 regression: 16 bp identical prefix + 84 bp non-homologous interior containing a
    // disconnected 6 bp exact match at 50..56 (~8.5 bits) under a 30.0-bit CRC gate.
    // Trimming must NOT subtract flank bits to let the disconnected 8.5-bit island pass.
    let pref = b"ACGTTGCACGTAGCTA"; // 16 bp balanced
    let mut x = pref.to_vec();
    let mut y = pref.to_vec();
    let mut noise_x = b"ACGTTGCACTAGTAGCACGTTGCACTAGTAGCAC".to_vec(); // 34 bp (16..50)
    let mut noise_y = b"TGCAAACTGATCAGTCTGCAAACTGATCAGTCTG".to_vec(); // 34 bp, 0% id
    noise_x.extend_from_slice(b"GCTACG"); // 6 bp disconnected island (50..56)
    noise_y.extend_from_slice(b"GCTACG");
    noise_x.extend_from_slice(b"ACGTTGCACTAGTAGCACGTTGCACTAGTAGCACGTTGCACTAG"); // 44 bp (56..100)
    noise_y.extend_from_slice(b"TGCAAACTGATCAGTCTGCAAACTGATCAGTCTGCAAACTGATC");
    x.extend_from_slice(&noise_x);
    y.extend_from_slice(&noise_y);
    assert_eq!(x.len(), 100);
    assert_eq!(y.len(), 100);

    let res = super::banded_mdl_gap(&x, &y, 32, 30.0);
    assert!(
        res.kept_pairs.is_empty(),
        "disconnected 6 bp island (~8.5 bits) leaked past 30.0-bit CRC gate: {:?}",
        res.kept_pairs
    );
}

#[test]
fn viterbi_fast_path_respects_per_block_threshold_across_indels() {
    // 64 bp balanced sequence split at position 32 by a 1-bp deletion: two ungapped blocks of
    // 32 bp (~49.9 bits) and 31 bp (~48.4 bits). At min_block_bits = 55.0, exact DP abstains
    // both halves; the Viterbi fast path must NOT sum across the indel to bypass min_block_bits.
    let x = b"ACGTACGTACGTTGCATGCATGCATGCAACGTACGTACGTACGTTGCATGCATGCATGCAACGT";
    let mut y = x[..32].to_vec();
    y.extend_from_slice(&x[33..]);
    let mut scratch = super::DpScratch::default();
    let (exact, _) = super::banded_mdl_gap_core_with_max_bits(
        x,
        &y,
        16,
        55.0,
        &super::HmmParams::default(),
        0.15,
        &mut scratch,
    );
    assert_eq!(exact.kept_pairs, [] as [(usize, usize); 0]);

    let fast = super::banded_mdl_gap(x, &y, 16, 55.0);
    assert!(
        fast.kept_pairs.is_empty(),
        "Viterbi fast path bypassed per-block min_block_bits across an indel: {:?}",
        fast.kept_pairs
    );
}

#[test]
fn prefix_suffix_trimming_uses_local_bocpd_on_homopolymer_flank() {
    // Whole-gap composition is balanced (A=30%, T=12%, C=29%, G=29%, GC=58%), so it passes
    // is_fast_path_eligible, while the 30-bp prefix is pure poly-A. Exact local BOCPD scores
    // A^30 at ~25 bits (not ~53.5 bits from whole-gap 58% GC), so a 35-bit CRC gate must abstain.
    let mut x = vec![b'A'; 30];
    let mut y = vec![b'A'; 30];
    // Two non-homologous 70-bp tails with exact counts (12 T, 29 C, 29 G) and no shifted homopolymers.
    let tail_x = b"CGTCGCTCGTCGCTCGTCGCTCGTCGCTCGTCGCTCGTCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG";
    let tail_y = b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGGTCGCTCGTCGCTCGTCGCTCGTCGCTCGTCGCTCGT";
    // Flip tail_y into an unrelated deterministic permutation with identical (12 T, 29 C, 29 G) counts:
    let mut rng = 0x1234_5678_9abc_def0_u64;
    let mut perm_x = tail_x.to_vec();
    let mut perm_y = tail_y.to_vec();
    for i in (1..perm_x.len()).rev() {
        rng = crate::splitmix64(rng);
        let j = usize::try_from(rng % (i as u64 + 1)).unwrap_or(0);
        perm_x.swap(i, j);
    }
    for i in (1..perm_y.len()).rev() {
        rng = crate::splitmix64(rng);
        let j = usize::try_from(rng % (i as u64 + 1)).unwrap_or(0);
        perm_y.swap(i, j);
    }
    x.extend_from_slice(&perm_x);
    y.extend_from_slice(&perm_y);
    assert_eq!(x.len(), 100);
    assert_eq!(y.len(), 100);

    let res = super::banded_mdl_gap(&x, &y, 32, 35.0);
    assert!(
        res.kept_pairs.is_empty(),
        "whole-gap GC leaked into poly-A prefix and bypassed 35-bit gate: {} pairs",
        res.kept_pairs.len()
    );
}
