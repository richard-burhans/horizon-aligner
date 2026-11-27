//! Phase 0 regression test suite (`VALIDATION_PLAN.md` P0.1-P0.6).

use horizon_aligner::coord::IntervalSet;
use horizon_aligner::io_util::{open_text_reader, open_text_writer};
use horizon_aligner::mdl::{HmmParams, banded_mdl_gap_posterior_trace};
use horizon_aligner::tba::{
    TbaBlock, TbaRow, VariationGfaBuilder, format_ancestors_summary_tsv,
    reconstruct_ancestral_contigs, reconstruct_ancestral_contigs_transducer,
};
use horizon_aligner::{
    AlignedRun, RgfaGraph, RunAlignment, format_maf_blocks_runs_avoiding_gaps,
    format_paf_records_runs_avoiding_gaps,
};
use std::collections::{HashMap, HashSet};
use std::io::{Read, Write};
use std::path::PathBuf;
use std::process::Command;

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("horizon_p0_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("mkdir");
    dir
}

fn rand_dna(seed: &mut u64, len: usize) -> Vec<u8> {
    let bases = *b"ACGT";
    (0..len)
        .map(|_| {
            *seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            bases[((*seed >> 62) & 3) as usize]
        })
        .collect()
}

fn mutate_jc69(seed: &mut u64, seq: &[u8], branch_t: f64) -> Vec<u8> {
    let p_change = 0.75 * (1.0 - (-4.0 * branch_t / 3.0).exp());
    let bases = *b"ACGT";
    seq.iter()
        .map(|&b| {
            *seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            #[allow(clippy::cast_precision_loss)]
            let u = ((*seed >> 11) as f64) / ((1u64 << 53) as f64);
            if u < p_change {
                *seed = seed
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                let alt = bases[((*seed >> 62) & 3) as usize];
                if alt == b {
                    bases[(bases.iter().position(|&x| x == b).unwrap_or(0) + 1) % 4]
                } else {
                    alt
                }
            } else {
                b
            }
        })
        .collect()
}

/// P0.1 (`R11-ONE-TO-ONE-OVERLAP`): runs separated by a <=40 bp gap claimed by another contig
/// must NOT be bridged in MAF or PAF output, guaranteeing exact 0 overlapping bp.
#[test]
fn p0_1_one_to_one_gap_bridging_emits_zero_overlapping_bp() {
    let mut seed = 0xCAFE_BABE_1234_5678;
    let seq_x = rand_dna(&mut seed, 200);
    let seq_y = seq_x.clone();

    // Contig pair A claims [0..50) and [70..120); contig pair B claims [50..70) (20 bp <= 40 bp).
    let runs_a = vec![
        AlignedRun {
            x_start: 0,
            y_start: 0,
            len: 50,
            is_rc: false,
        },
        AlignedRun {
            x_start: 70,
            y_start: 70,
            len: 50,
            is_rc: false,
        },
    ];
    let mut forbidden_x = IntervalSet::default();
    let mut forbidden_y = IntervalSet::default();
    let mut scratch = Vec::new();
    forbidden_x.insert(50, 70, &mut scratch);
    forbidden_y.insert(50, 70, &mut scratch);

    let aln = RunAlignment {
        aligned_runs: runs_a,
        total_pairs: 100,
        lcbs_fwd: 1,
        lcbs_rc: 0,
        macro_bubbles: 0,
        internal_bubbles: 0,
    };

    let maf = format_maf_blocks_runs_avoiding_gaps(
        "ref.c1",
        seq_x.len(),
        "qry.q1",
        seq_y.len(),
        &seq_x,
        &seq_y,
        &aln,
        Some(&forbidden_x),
        Some(&forbidden_y),
    );
    let mut maf_blocks = 0usize;
    for line in maf.lines() {
        if let Some(rest) = line.strip_prefix("s ") {
            let f: Vec<&str> = rest.split_whitespace().collect();
            let start: usize = f[1].parse().expect("start");
            let size: usize = f[2].parse().expect("size");
            let end = start + size;
            assert!(
                end <= 50 || start >= 70,
                "MAF bridged into forbidden [50..70): [{start}..{end})"
            );
            if f[0].starts_with("ref.") {
                maf_blocks += 1;
            }
        }
    }
    assert_eq!(
        maf_blocks, 2,
        "must split into 2 non-overlapping MAF blocks"
    );

    let paf = format_paf_records_runs_avoiding_gaps(
        "c1",
        seq_x.len(),
        "q1",
        seq_y.len(),
        &seq_x,
        &seq_y,
        &aln,
        Some(&forbidden_x),
        Some(&forbidden_y),
    );
    let mut paf_rows = 0usize;
    for line in paf.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let qs: usize = f[2].parse().expect("qs");
        let qe: usize = f[3].parse().expect("qe");
        let ts: usize = f[7].parse().expect("ts");
        let te: usize = f[8].parse().expect("te");
        assert!(
            te <= 50 || ts >= 70,
            "PAF bridged target [50..70): [{ts}..{te})"
        );
        assert!(
            qe <= 50 || qs >= 70,
            "PAF bridged query [50..70): [{qs}..{qe})"
        );
        paf_rows += 1;
    }
    assert_eq!(paf_rows, 2, "must split into 2 non-overlapping PAF records");
}

/// P0.2: `horizon_gfa stats` and `horizon_tba --from-mafs` fail closed (`rc != 0`) on garbage.
#[test]
fn p0_2_gfa_and_tba_fail_closed_on_garbage() {
    let dir = temp_dir("fail_closed");
    let bad_gfa = dir.join("garbage.gfa");
    std::fs::write(&bad_gfa, "THIS IS NOT A VALID GFA FILE\n").expect("write bad gfa");

    let gfa_out = Command::new(env!("CARGO_BIN_EXE_horizon_gfa"))
        .args(["stats", bad_gfa.to_str().expect("utf8")])
        .output()
        .expect("run horizon_gfa");
    assert!(
        !gfa_out.status.success(),
        "horizon_gfa stats must fail closed on garbage input"
    );

    let bad_maf_dir = dir.join("bad_mafs");
    std::fs::create_dir_all(&bad_maf_dir).expect("mkdir");
    std::fs::write(bad_maf_dir.join("spA__spB.maf"), "GARBAGE MAF LINE\n").expect("write bad maf");
    let out_prefix = dir.join("out");
    let tba_out = Command::new(env!("CARGO_BIN_EXE_horizon_tba"))
        .args([
            "--from-mafs",
            bad_maf_dir.to_str().expect("utf8"),
            "--tree",
            "(spA,spB);",
            out_prefix.to_str().expect("utf8"),
            "spA",
            "spB",
        ])
        .output()
        .expect("run horizon_tba");
    assert!(
        !tba_out.status.success(),
        "horizon_tba --from-mafs must fail closed on malformed MAF"
    );
    let tba_err = String::from_utf8_lossy(&tba_out.stderr);
    assert!(
        tba_err.contains("malformed MAF"),
        "expected MAF parse error, got: {tba_err}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// P0.3: `VariationGFA` and rGFA emit globally unique `SR:i:` per `SN:Z:`.
#[test]
fn p0_3_unique_sr_rank_per_sn_in_variation_gfa_and_rgfa() {
    // Variation GFA: sp1 present in block 1, absent in block 2 (only sp2 + sp3).
    let mk = |sp: &str, s: usize, seq: &[u8]| TbaRow {
        species: sp.to_owned(),
        contig: "chr1".to_owned(),
        start: s,
        size: seq.iter().filter(|&&b| b != b'-').count(),
        strand: '+',
        src_size: 100,
        seq: seq.to_vec(),
    };
    let blocks = vec![
        TbaBlock {
            score: 10.0,
            rows: vec![mk("sp1", 0, b"ACGT"), mk("sp2", 0, b"ACGG")],
        },
        TbaBlock {
            score: 10.0,
            rows: vec![mk("sp2", 4, b"TTTT"), mk("sp3", 0, b"TTTA")],
        },
    ];
    let mut gfa = Vec::new();
    let mut builder = VariationGfaBuilder::new(&mut gfa).expect("builder");
    for b in &blocks {
        builder.push_block(b).expect("push");
    }
    builder.finish().expect("finish");
    let gfa_str = String::from_utf8(gfa).expect("utf8");
    let mut sn_to_sr: HashMap<String, HashSet<String>> = HashMap::new();
    for line in gfa_str.lines() {
        if line.starts_with("S\t") {
            let mut sn = None;
            let mut sr = None;
            for tok in line.split('\t').skip(3) {
                if let Some(v) = tok.strip_prefix("SN:Z:") {
                    sn = Some(v.to_owned());
                } else if let Some(v) = tok.strip_prefix("SR:i:") {
                    sr = Some(v.to_owned());
                }
            }
            if let (Some(n), Some(r)) = (sn, sr) {
                sn_to_sr.entry(n).or_default().insert(r);
            }
        }
    }
    for (sn, ranks) in &sn_to_sr {
        assert_eq!(
            ranks.len(),
            1,
            "variation GFA SN `{sn}` has conflicting SR ranks: {ranks:?}"
        );
    }

    // rGFA: ref and query share contig name `chr1` with a 100 bp insertion.
    let mut seed = 0x1234_9876;
    let left = rand_dna(&mut seed, 400);
    let ins = rand_dna(&mut seed, 100);
    let right = rand_dna(&mut seed, 400);
    let mut ref_seq = left.clone();
    ref_seq.extend_from_slice(&right);
    let mut qry_seq = left;
    qry_seq.extend_from_slice(&ins);
    qry_seq.extend_from_slice(&right);

    let mut graph = RgfaGraph::from_reference(&[("chr1", &ref_seq)], 50);
    let grafted = graph.graft_assembly("chr1", &qry_seq);
    assert!(grafted >= 1, "expected >=1 grafted SV");
    let mut rgfa_sn_sr: HashMap<String, HashSet<u32>> = HashMap::new();
    for seg in &graph.segments {
        rgfa_sn_sr
            .entry(seg.origin_name.clone())
            .or_default()
            .insert(seg.rank);
    }
    for (sn, ranks) in &rgfa_sn_sr {
        assert_eq!(
            ranks.len(),
            1,
            "rGFA SN `{sn}` has conflicting SR ranks: {ranks:?}"
        );
    }
}

/// P0.4: BGZF `.vcf.gz` output, `PanSN` diploid merging, unaligned singleton insertion splicing,
/// and `>= 200 bp` shared insertion / deletion calling in `horizon_tba` -> `horizon_resolve --vcf`.
#[test]
#[allow(clippy::too_many_lines)]
fn p0_4_bgzf_pansn_diploid_and_large_indels_end_to_end() {
    let dir = temp_dir("vcf_giab");
    let mut seed = 0x9988_7766_5544_3322;
    let flank_a = rand_dna(&mut seed, 500);
    let shared_ins_300 = rand_dna(&mut seed, 300);
    let flank_b = rand_dna(&mut seed, 500);
    let single_ins_300 = rand_dna(&mut seed, 300);
    let flank_c = rand_dna(&mut seed, 500);
    let del_300 = rand_dna(&mut seed, 300);
    let flank_d = rand_dna(&mut seed, 500);

    // Ref:     A +     B +     C + del_300 + D
    // HG002#1: A + ins + B + s_ins + C + del_300 + D  (hom 300bp ins, het 300bp singleton ins)
    // HG002#2: A + ins + B +     C +         D        (hom 300bp ins, het 300bp del)
    let mut ref_seq = Vec::new();
    ref_seq.extend_from_slice(&flank_a);
    ref_seq.extend_from_slice(&flank_b);
    ref_seq.extend_from_slice(&flank_c);
    ref_seq.extend_from_slice(&del_300);
    ref_seq.extend_from_slice(&flank_d);

    let mut h1_seq = Vec::new();
    h1_seq.extend_from_slice(&flank_a);
    h1_seq.extend_from_slice(&shared_ins_300);
    h1_seq.extend_from_slice(&flank_b);
    h1_seq.extend_from_slice(&single_ins_300);
    h1_seq.extend_from_slice(&flank_c);
    h1_seq.extend_from_slice(&del_300);
    h1_seq.extend_from_slice(&flank_d);

    let mut h2_seq = Vec::new();
    h2_seq.extend_from_slice(&flank_a);
    h2_seq.extend_from_slice(&shared_ins_300);
    h2_seq.extend_from_slice(&flank_b);
    h2_seq.extend_from_slice(&flank_c);
    h2_seq.extend_from_slice(&flank_d);

    let write_fa = |name: &str, seq: &[u8]| {
        let p = dir.join(format!("{name}.fa"));
        std::fs::write(&p, format!(">chr1\n{}\n", String::from_utf8_lossy(seq))).expect("write fa");
        p
    };
    let ref_fa = write_fa("Ref", &ref_seq);
    let h1_fa = write_fa("HG002_1", &h1_seq);
    let h2_fa = write_fa("HG002_2", &h2_seq);
    let tba_prefix = dir.join("panel");

    let tba = Command::new(env!("CARGO_BIN_EXE_horizon_tba"))
        .args([
            "--tree",
            "(Ref,(HG002#1,HG002#2));",
            "--variation-gfa",
            tba_prefix.to_str().expect("utf8"),
            "0",
            &format!("Ref={}", ref_fa.display()),
            &format!("HG002#1={}", h1_fa.display()),
            &format!("HG002#2={}", h2_fa.display()),
        ])
        .output()
        .expect("run horizon_tba");
    assert!(
        tba.status.success(),
        "horizon_tba failed: {}",
        String::from_utf8_lossy(&tba.stderr)
    );

    let gfa_path = dir.join("panel.gfa");
    let resolved_gfa = dir.join("resolved.gfa");
    let vcf_gz = dir.join("panel.vcf.gz");
    let res = Command::new(env!("CARGO_BIN_EXE_horizon_resolve"))
        .args([
            gfa_path.to_str().expect("utf8"),
            "--out",
            resolved_gfa.to_str().expect("utf8"),
            "--vcf",
            vcf_gz.to_str().expect("utf8"),
            "--vcf-ref",
            "Ref",
        ])
        .output()
        .expect("run horizon_resolve");
    assert!(
        res.status.success(),
        "horizon_resolve failed: {}",
        String::from_utf8_lossy(&res.stderr)
    );

    // Verify true BGZF magic (`FEXTRA` + `BC` subfield) + 28-byte BGZF EOF block.
    let raw_bytes = std::fs::read(&vcf_gz).expect("read vcf.gz");
    assert!(raw_bytes.len() > 28, "BGZF output too short");
    assert_eq!(&raw_bytes[0..4], &[0x1f, 0x8b, 0x08, 0x04], "BGZF header");
    assert_eq!(&raw_bytes[12..14], b"BC", "BGZF BC subfield");
    assert_eq!(
        &raw_bytes[raw_bytes.len() - 28..],
        &[
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
        ],
        "BGZF EOF block"
    );

    let mut vcf_text = String::new();
    open_text_reader(&vcf_gz)
        .expect("open bgzf")
        .read_to_string(&mut vcf_text)
        .expect("read bgzf");

    // Verify diploid PanSN sample column (`HG002`, NOT `HG002#1`/`HG002#2`).
    let header = vcf_text
        .lines()
        .find(|l| l.starts_with("#CHROM"))
        .expect("#CHROM header");
    let cols: Vec<&str> = header.split('\t').collect();
    let hg_idx = cols
        .iter()
        .position(|&c| c == "HG002")
        .expect("diploid HG002 column");
    assert!(
        !cols.contains(&"HG002#1") && !cols.contains(&"HG002#2"),
        "PanSN haplotypes must merge into diploid HG002 column: {cols:?}"
    );

    let mut saw_hom_300_ins = false;
    let mut saw_het_300_ins = false;
    let mut saw_het_300_del = false;
    for line in vcf_text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        let ref_len = f[3].len();
        let alt_len = f[4].len();
        let gt = f[hg_idx];
        if alt_len >= ref_len + 280 && gt == "1|1" {
            saw_hom_300_ins = true;
        }
        if alt_len >= ref_len + 280 && (gt == "1|0" || gt == "0|1") {
            saw_het_300_ins = true;
        }
        if ref_len >= alt_len + 280 && (gt == "1|0" || gt == "0|1") {
            saw_het_300_del = true;
        }
    }
    assert!(
        saw_hom_300_ins,
        "must call >=200bp shared homozygous 300bp insertion as 1|1:\n{vcf_text}"
    );
    assert!(
        saw_het_300_ins,
        "must call singleton unaligned 300bp insertion as 1|0:\n{vcf_text}"
    );
    assert!(
        saw_het_300_del,
        "must call >=200bp heterozygous 300bp deletion as 0|1:\n{vcf_text}"
    );

    // Also verify BGZF round-trip on multi-block payload (>130 KiB).
    let big_path = dir.join("multi.vcf.gz");
    let big_payload = "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\n".repeat(4000);
    {
        let mut w = open_text_writer(Some(big_path.to_str().expect("utf8"))).expect("open");
        w.write_all(big_payload.as_bytes()).expect("write");
        w.finish().expect("finish");
    }
    let mut got_big = String::new();
    open_text_reader(&big_path)
        .expect("read")
        .read_to_string(&mut got_big)
        .expect("read_to_string");
    assert_eq!(got_big, big_payload);

    let _ = std::fs::remove_dir_all(&dir);
}

/// P0.5 (`R11-ANCESTRAL-ML-OUTGROUP`): exact 2-level outgroup clade peeling beats Fitch on BOTH
/// correlated sister outgroups `(t_stem=0.20, t_tip=0.01)` AND nearly-independent outgroups
/// `(t_stem=0.01, t_tip=0.20)`.
#[test]
fn p0_5_ancestral_ml_outgroup_peeling_beats_fitch_correlated_and_independent() {
    let n = 10_000usize;
    let mut seed = 0x55AA_1122_3344_7788;
    let anc_true = rand_dna(&mut seed, n);

    let mk = |sp: &str, seq: &[u8]| TbaRow {
        species: sp.to_owned(),
        contig: "chr1".to_owned(),
        start: 0,
        size: seq.len(),
        strand: '+',
        src_size: seq.len(),
        seq: seq.to_vec(),
    };
    let err_rate = |pred: &[u8]| -> f64 {
        assert_eq!(pred.len(), anc_true.len());
        let mm = pred.iter().zip(&anc_true).filter(|(a, b)| a != b).count();
        #[allow(clippy::cast_precision_loss)]
        {
            mm as f64 / anc_true.len() as f64
        }
    };

    let ogs = vec!["O1".to_owned(), "O2".to_owned(), "O3".to_owned()];

    // Case A: Correlated sister outgroups (tL=0.05, tR=0.40, t_stem=0.20, t_tip=0.01).
    let l_a = mutate_jc69(&mut seed, &anc_true, 0.05);
    let r_a = mutate_jc69(&mut seed, &anc_true, 0.40);
    let og_mrca_a = mutate_jc69(&mut seed, &anc_true, 0.20);
    let o1_a = mutate_jc69(&mut seed, &og_mrca_a, 0.01);
    let o2_a = mutate_jc69(&mut seed, &og_mrca_a, 0.01);
    let o3_a = mutate_jc69(&mut seed, &og_mrca_a, 0.01);
    let blk_a = TbaBlock {
        score: 1000.0,
        rows: vec![
            mk("L", &l_a),
            mk("R", &r_a),
            mk("O1", &o1_a),
            mk("O2", &o2_a),
            mk("O3", &o3_a),
        ],
    };
    let fitch_a =
        reconstruct_ancestral_contigs("Anc0", "L", "R", &ogs, std::slice::from_ref(&blk_a));
    let branches_a = horizon_aligner::tba::AncestralBranchLengths {
        left_t: 0.05,
        right_t: 0.40,
        outgroup_t: 0.21,
        outgroup_stem_t: 0.20,
        outgroup_tip_t: 0.01,
        ins_rate: 0.05,
        del_rate: 0.05,
    };
    let ml_a =
        reconstruct_ancestral_contigs_transducer("Anc0", "L", "R", &ogs, &[blk_a], branches_a);
    let err_fitch_a = err_rate(&fitch_a.anc_contigs[0].1);
    let err_ml_a = err_rate(&ml_a.anc_contigs[0].1);
    assert!(
        err_ml_a < err_fitch_a,
        "correlated-outgroup ML ({err_ml_a:.4}) must strictly beat Fitch ({err_fitch_a:.4})"
    );
    let tsv = format_ancestors_summary_tsv(&[(
        "Anc0".to_owned(),
        ml_a.anc_contigs.clone(),
        ml_a.stats.clone(),
    )]);
    assert!(tsv.contains("left_t\tright_t\toutgroup_t"));
    assert!(tsv.contains("0.050000\t0.400000\t0.210000"));

    // Case B: Nearly-independent outgroups (tL=0.05, tR=0.40, t_stem=0.01, t_tip=0.20).
    let l_b = mutate_jc69(&mut seed, &anc_true, 0.05);
    let r_b = mutate_jc69(&mut seed, &anc_true, 0.40);
    let og_mrca_b = mutate_jc69(&mut seed, &anc_true, 0.01);
    let o1_b = mutate_jc69(&mut seed, &og_mrca_b, 0.20);
    let o2_b = mutate_jc69(&mut seed, &og_mrca_b, 0.20);
    let o3_b = mutate_jc69(&mut seed, &og_mrca_b, 0.20);
    let blk_b = TbaBlock {
        score: 1000.0,
        rows: vec![
            mk("L", &l_b),
            mk("R", &r_b),
            mk("O1", &o1_b),
            mk("O2", &o2_b),
            mk("O3", &o3_b),
        ],
    };
    let fitch_b =
        reconstruct_ancestral_contigs("Anc0", "L", "R", &ogs, std::slice::from_ref(&blk_b));
    let branches_b = horizon_aligner::tba::AncestralBranchLengths {
        left_t: 0.05,
        right_t: 0.40,
        outgroup_t: 0.21,
        outgroup_stem_t: 0.01,
        outgroup_tip_t: 0.20,
        ins_rate: 0.05,
        del_rate: 0.05,
    };
    let ml_b =
        reconstruct_ancestral_contigs_transducer("Anc0", "L", "R", &ogs, &[blk_b], branches_b);
    let err_fitch_b = err_rate(&fitch_b.anc_contigs[0].1);
    let err_ml_b = err_rate(&ml_b.anc_contigs[0].1);
    assert!(
        err_ml_b < err_fitch_b && err_ml_b < 0.025,
        "independent-outgroup ML ({err_ml_b:.4}) must beat Fitch ({err_fitch_b:.4}) and stay <2.5%"
    );
}

/// P0.6: Forward-Backward partition-function parity `|log2(Z_F) - log2(Z_B)| <= 1e-6 * max(1, |Z_F|)`.
#[test]
fn p0_6_pair_hmm_forward_backward_partition_parity() {
    let mut seed = 0xDEAD_BEEF_CAFE_0001;
    let params = HmmParams::default();
    let seq1 = rand_dna(&mut seed, 140);
    let mut seq2 = mutate_jc69(&mut seed, &seq1, 0.15);
    seq2.splice(60..60, b"ACGTACGT".iter().copied());
    let mut gc_shift = b"ATATATATATATATATATATATATATATATATATATATAT".to_vec();
    gc_shift.extend_from_slice(b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC");
    let gc_shift_mut = mutate_jc69(&mut seed, &gc_shift, 0.20);

    for (x, y) in [
        (&seq1[..], &seq1[..]),
        (&seq1[..], &seq2[..]),
        (&gc_shift[..], &gc_shift_mut[..]),
    ] {
        let trace = banded_mdl_gap_posterior_trace(x, y, 40, 10.0, &params, 0.5);
        assert!(trace.score_bits.is_finite() && trace.log2_zb.is_finite());
        let diff = (trace.score_bits - trace.log2_zb).abs();
        let tol = 1e-6 * trace.score_bits.abs().max(1.0);
        assert!(
            diff <= tol,
            "Forward/Backward partition mismatch: log2(Z_F)={}, log2(Z_B)={}, diff={diff} > tol={tol}",
            trace.score_bits,
            trace.log2_zb
        );
    }
}
