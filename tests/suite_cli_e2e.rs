//! End-to-end integration tests for `horizon_guide --progressive --hal` + `horizon_wga --hal-info`,
//! and the GFA pangenome pipeline (`horizon_clean` -> `horizon_resolve` -> `horizon_lift`).

mod common;

use common::{mutate_dna, synthetic_dna};
use std::process::Command;

#[test]
fn horizon_guide_progressive_hal_roundtrips_through_hal_info() {
    let dir = std::env::temp_dir().join(format!("horizon_guide_e2e_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("mkdir");

    let anc = synthetic_dna(0xCAFE_BABE, 2_000);
    let sp1 = mutate_dna(&anc, 11, 2);
    let sp2 = mutate_dna(&anc, 22, 2);
    let sp3 = mutate_dna(&anc, 33, 3);

    let fa1 = dir.join("sp1.fa");
    let fa2 = dir.join("sp2.fa");
    let fa3 = dir.join("sp3.fa");
    std::fs::write(&fa1, format!(">chr1\n{}\n", String::from_utf8_lossy(&sp1))).expect("fa1");
    std::fs::write(&fa2, format!(">chr1\n{}\n", String::from_utf8_lossy(&sp2))).expect("fa2");
    std::fs::write(&fa3, format!(">chr1\n{}\n", String::from_utf8_lossy(&sp3))).expect("fa3");

    let out_prefix = dir.join("guide_out");
    let guide = Command::new(env!("CARGO_BIN_EXE_horizon_guide"))
        .args([
            "--progressive",
            "--hal",
            "-q",
            "-t",
            "2",
            out_prefix.to_str().expect("prefix"),
            "1",
            &format!("sp1={}", fa1.display()),
            &format!("sp2={}", fa2.display()),
            &format!("sp3={}", fa3.display()),
        ])
        .output()
        .expect("run horizon_guide");
    assert!(
        guide.status.success(),
        "horizon_guide failed: {}",
        String::from_utf8_lossy(&guide.stderr)
    );
    let hal_path = dir.join("guide_out.hal");
    assert!(hal_path.is_file(), "expected progressive .hal file");

    let info = Command::new(env!("CARGO_BIN_EXE_horizon_wga"))
        .args(["--hal-info", hal_path.to_str().expect("hal")])
        .output()
        .expect("run horizon_wga --hal-info");
    assert!(
        info.status.success(),
        "horizon_wga --hal-info failed: {}",
        String::from_utf8_lossy(&info.stderr)
    );
    let info_stdout = String::from_utf8_lossy(&info.stdout);
    assert!(
        info_stdout.contains("Anc0") && info_stdout.contains("Anc1"),
        "expected progressive ancestral genomes in --hal-info: {info_stdout}"
    );
    for sp in ["sp1", "sp2", "sp3"] {
        assert!(
            info_stdout.contains(sp),
            "expected leaf species {sp} in --hal-info: {info_stdout}"
        );
    }

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn gfa_pipeline_clean_resolve_lift_interoperates_end_to_end() {
    let dir = std::env::temp_dir().join(format!("horizon_gfa_pipe_e2e_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("mkdir");

    let in_gfa = dir.join("in.gfa");
    let clean_gfa = dir.join("clean.gfa");
    let resolved_gfa = dir.join("resolved.gfa");
    let bed_path = dir.join("regions.bed");

    let gfa_content = "\
H\tVN:Z:1.1
S\ts1\tACGTACGTACGTACGTACGT
S\ts2\tGGGGTTTTCCCCAAAAGGGG
S\ts3\tTGCATGCATGCATGCATGCA
L\ts1\t+\ts2\t+\t0M
L\ts2\t+\ts2\t+\t0M
L\ts2\t+\ts3\t+\t0M
W\tspA\t0\tchr1\t0\t60\t>s1>s2>s3
W\tspB\t0\tchr1\t0\t80\t>s1>s2>s2>s3
";
    std::fs::write(&in_gfa, gfa_content).expect("write gfa");
    std::fs::write(&bed_path, "chr1\t0\t20\tfeat_s1\t0\t+\n").expect("write bed");

    let clean = Command::new(env!("CARGO_BIN_EXE_horizon_clean"))
        .args(["--gfa", in_gfa.to_str().expect("in_gfa")])
        .output()
        .expect("run horizon_clean");
    assert!(
        clean.status.success(),
        "horizon_clean failed: {}",
        String::from_utf8_lossy(&clean.stderr)
    );
    std::fs::write(&clean_gfa, &clean.stdout).expect("write clean_gfa");

    let resolve = Command::new(env!("CARGO_BIN_EXE_horizon_resolve"))
        .args([
            clean_gfa.to_str().expect("clean_gfa"),
            "--collapse-dups",
            "--tag-walks",
        ])
        .output()
        .expect("run horizon_resolve");
    assert!(
        resolve.status.success(),
        "horizon_resolve failed: {}",
        String::from_utf8_lossy(&resolve.stderr)
    );
    std::fs::write(&resolved_gfa, &resolve.stdout).expect("write resolved_gfa");
    let resolved_text = String::from_utf8_lossy(&resolve.stdout);
    assert!(
        resolved_text.contains("W\tspB\t0\tchr1\t0\t60\t>s1>s2>s3")
            && resolved_text.contains("EG:Z:"),
        "expected horizon_resolve to collapse >s2>s2 tandem dup and emit EG:Z: tag, got:\n{resolved_text}"
    );

    let lift = Command::new(env!("CARGO_BIN_EXE_horizon_lift"))
        .args([
            "--gfa",
            resolved_gfa.to_str().expect("resolved_gfa"),
            "--bed",
            bed_path.to_str().expect("bed"),
            "--from",
            "spA",
            "--to",
            "spB",
        ])
        .output()
        .expect("run horizon_lift");
    assert!(
        lift.status.success(),
        "horizon_lift failed: {}",
        String::from_utf8_lossy(&lift.stderr)
    );
    let lift_stdout = String::from_utf8_lossy(&lift.stdout);
    assert!(
        lift_stdout.contains("feat_s1") && lift_stdout.contains("EXACT"),
        "expected EXACT liftover of feat_s1 across clean->resolve->lift pipeline: {lift_stdout}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn horizon_resolve_vcf_and_horizon_clean_mapping_filters_cli_e2e() {
    let dir = std::env::temp_dir().join(format!(
        "horizon_mc_features_e2e_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("mkdir");

    let in_gfa = dir.join("pangenome.gfa");
    let resolved_gfa = dir.join("resolved.gfa");
    let vcf_out = dir.join("calls.vcf");
    let mapping_gfa = dir.join("mapping_ready.gfa");

    let big_bubble = "CGTA".repeat(20); // 80 bp HB:i:1
    let gfa_content = format!(
        "H\tVN:Z:1.1\n\
         S\ta\tACGTACGTAC\n\
         S\tb_ref\tA\n\
         S\tb_alt\tG\n\
         S\tc\tTTTTTCCCCC\n\
         S\tpriv_ins\tGGGGAAAA\n\
         S\tbig_hb\t{big_bubble}\tHB:i:1\n\
         S\td\tACACACACAC\n\
         L\ta\t+\tb_ref\t+\t0M\n\
         L\ta\t+\tb_alt\t+\t0M\n\
         L\tb_ref\t+\tc\t+\t0M\n\
         L\tb_alt\t+\tc\t+\t0M\n\
         L\tc\t+\td\t+\t0M\n\
         L\tc\t+\tpriv_ins\t+\t0M\n\
         L\tpriv_ins\t+\td\t+\t0M\n\
         L\td\t+\tbig_hb\t+\t0M\n\
         W\tCHM13\t0\tchr1\t0\t31\t>a>b_ref>c>d\n\
         W\tHG002\t1\tchr1\t0\t39\t>a>b_alt>c>priv_ins>d\n\
         W\tHG002\t2\tchr1\t0\t111\t>a>b_ref>c>d>big_hb\n"
    );
    std::fs::write(&in_gfa, gfa_content).expect("write gfa");

    // 1) horizon_resolve --vcf --vcf-ref CHM13
    let resolve = Command::new(env!("CARGO_BIN_EXE_horizon_resolve"))
        .args([
            in_gfa.to_str().expect("in_gfa"),
            "--out",
            resolved_gfa.to_str().expect("resolved_gfa"),
            "--vcf",
            vcf_out.to_str().expect("vcf_out"),
            "--vcf-ref",
            "CHM13",
        ])
        .output()
        .expect("run horizon_resolve --vcf");
    assert!(
        resolve.status.success(),
        "horizon_resolve --vcf failed: {}",
        String::from_utf8_lossy(&resolve.stderr)
    );
    let vcf_text = std::fs::read_to_string(&vcf_out).expect("read vcf");
    assert!(vcf_text.contains("##fileformat=VCFv4.2\n"));
    assert!(
        vcf_text.contains("chr1\t11\tsnarl_1\tA\tG\t.\tPASS\t"),
        "expected SNP at chr1:11, got:\n{vcf_text}"
    );
    assert!(
        vcf_text.contains("chr1\t21\tsnarl_2\tC\tCGGGGAAAA\t.\tPASS\t"),
        "expected left-anchored insertion at chr1:21, got:\n{vcf_text}"
    );
    assert!(
        vcf_text.contains("\t0\t1|0\n"),
        "expected CHM13=0 and HG002=1|0 phased genotype, got:\n{vcf_text}"
    );

    // 2) horizon_clean --min-haplotypes 2 --clip-bubbles-bp 50
    let clean = Command::new(env!("CARGO_BIN_EXE_horizon_clean"))
        .args([
            "--gfa",
            in_gfa.to_str().expect("in_gfa"),
            "--out",
            mapping_gfa.to_str().expect("mapping_gfa"),
            "--min-cp-bits",
            "100.0",
            "--min-haplotypes",
            "2",
            "--clip-bubbles-bp",
            "50",
        ])
        .output()
        .expect("run horizon_clean");
    assert!(
        clean.status.success(),
        "horizon_clean failed: {}",
        String::from_utf8_lossy(&clean.stderr)
    );
    let clean_text = std::fs::read_to_string(&mapping_gfa).expect("read mapping_gfa");
    assert!(clean_text.contains("S\ta\t"));
    assert!(clean_text.contains("S\tc\t"));
    assert!(clean_text.contains("S\td\t"));
    assert!(
        !clean_text.contains("priv_ins"),
        "private 1-haplotype segment must be pruned:\n{clean_text}"
    );
    assert!(
        !clean_text.contains("big_hb"),
        "80bp HB:i:1 bubble must be clipped at 50bp threshold:\n{clean_text}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}
