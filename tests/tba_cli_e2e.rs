//! End-to-end integration tests for `horizon_tba` (E2E FASTAs, `--pairwise-dir`, and `--from-mafs`).

mod common;

use common::{mutate_dna, synthetic_dna};
use std::process::Command;

#[test]
fn horizon_tba_e2e_pairwise_dir_and_from_mafs_agree() {
    let dir = std::env::temp_dir().join(format!("horizon_tba_e2e_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let pw_dir = dir.join("pairwise");
    std::fs::create_dir_all(&dir).expect("mkdir");

    // Create 3 homologous synthetic 2kb genomes so all 3 species align into a 3-way core block.
    let anc = synthetic_dna(0x1234_5678, 2_000);
    let sp1 = mutate_dna(&anc, 1, 2);
    let sp2 = mutate_dna(&anc, 2, 2);
    let sp3 = mutate_dna(&anc, 3, 3);

    let fa1 = dir.join("sp1.fa");
    let fa2 = dir.join("sp2.fa");
    let fa3 = dir.join("sp3.fa");
    std::fs::write(&fa1, format!(">chr1\n{}\n", String::from_utf8_lossy(&sp1))).expect("fa1");
    std::fs::write(&fa2, format!(">chr1\n{}\n", String::from_utf8_lossy(&sp2))).expect("fa2");
    std::fs::write(&fa3, format!(">chr1\n{}\n", String::from_utf8_lossy(&sp3))).expect("fa3");

    let tba_bin = env!("CARGO_BIN_EXE_horizon_tba");
    let out_e2e = dir.join("out_e2e");
    let status = Command::new(tba_bin)
        .args([
            "-q",
            "--tree",
            "((sp1,sp2),sp3);",
            "--pairwise-dir",
            pw_dir.to_str().expect("pw_dir"),
            out_e2e.to_str().expect("out_e2e"),
            "1",
            &format!("sp1={}", fa1.display()),
            &format!("sp2={}", fa2.display()),
            &format!("sp3={}", fa3.display()),
        ])
        .output()
        .expect("run horizon_tba e2e");
    assert!(
        status.status.success(),
        "horizon_tba E2E failed: {}",
        String::from_utf8_lossy(&status.stderr)
    );
    let e2e_stdout = String::from_utf8_lossy(&status.stdout);
    assert!(
        e2e_stdout.contains("mode=tba species=3"),
        "stdout: {e2e_stdout}"
    );
    assert!(
        e2e_stdout.contains("core_blocks=") && !e2e_stdout.contains("core_blocks=0"),
        "expected >= 1 3-species core block: {e2e_stdout}"
    );
    assert!(pw_dir.join("sp1__sp2.maf").is_file());
    assert!(pw_dir.join("sp1__sp3.maf").is_file());
    assert!(pw_dir.join("sp2__sp3.maf").is_file());

    // Now run --from-mafs on the saved pairwise directory and verify identical block output.
    let out_from = dir.join("out_from");
    let status_from = Command::new(tba_bin)
        .args([
            "-q",
            "--tree",
            "((sp1,sp2),sp3);",
            "--from-mafs",
            pw_dir.to_str().expect("pw_dir"),
            out_from.to_str().expect("out_from"),
            "sp1",
            "sp2",
            "sp3",
        ])
        .output()
        .expect("run horizon_tba --from-mafs");
    assert!(
        status_from.status.success(),
        "horizon_tba --from-mafs failed: {}",
        String::from_utf8_lossy(&status_from.stderr)
    );

    let maf_e2e = std::fs::read_to_string(dir.join("out_e2e.maf")).expect("read e2e maf");
    let maf_from = std::fs::read_to_string(dir.join("out_from.maf")).expect("read from maf");
    assert_eq!(
        maf_e2e, maf_from,
        "direct E2E and --from-mafs merge-only mode must produce byte-identical MAF output"
    );

    // Verify the emitted GFA has valid W-line walks for all 3 species and parses via horizon_gfa.
    let gfa_e2e = std::fs::read_to_string(dir.join("out_e2e.gfa")).expect("read e2e gfa");
    for sp in ["sp1", "sp2", "sp3"] {
        assert!(
            gfa_e2e.contains(&format!("W\t{sp}\t0\tchr1\t")),
            "expected GFA W walk for {sp} in:\n{gfa_e2e}"
        );
    }
    let gfa_bin = env!("CARGO_BIN_EXE_horizon_gfa");
    let gfa_stats = Command::new(gfa_bin)
        .args(["stats", dir.join("out_e2e.gfa").to_str().expect("gfa")])
        .output()
        .expect("run horizon_gfa stats");
    assert!(
        gfa_stats.status.success(),
        "horizon_gfa stats failed on TBA GFA: {}",
        String::from_utf8_lossy(&gfa_stats.stderr)
    );
    let gfa_stdout = String::from_utf8_lossy(&gfa_stats.stdout);
    assert!(
        gfa_stdout.contains("walks=3"),
        "expected walks=3 in horizon_gfa stats: {gfa_stdout}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

fn run_bin_stdout(bin: &str, args: &[&str], label: &str) -> String {
    let res = Command::new(bin).args(args).output().expect(label);
    assert!(
        res.status.success(),
        "{label} failed: {}",
        String::from_utf8_lossy(&res.stderr)
    );
    String::from_utf8_lossy(&res.stdout).into_owned()
}

fn build_two_tier_fixture(dir: &std::path::Path) -> [String; 4] {
    let anc = synthetic_dna(0xCAFE_BABE, 2_000);
    let h1 = mutate_dna(&anc, 11, 2);
    let h2 = mutate_dna(&anc, 12, 5);
    let t1 = mutate_dna(&h1, 13, 1);
    let t2 = mutate_dna(&h2, 14, 1);
    let write_fa = |name: &str, seq: &[u8]| {
        let p = dir.join(format!("{name}.fa"));
        std::fs::write(&p, format!(">chr1\n{}\n", String::from_utf8_lossy(seq))).expect("write fa");
        format!("{name}={}", p.display())
    };
    [
        write_fa("h1", &h1),
        write_fa("h2", &h2),
        write_fa("t1", &t1),
        write_fa("t2", &t2),
    ]
}

#[test]
fn horizon_tba_two_tier_hinges_and_post_hoc_extend_maf_work_end_to_end() {
    let dir = std::env::temp_dir().join(format!("horizon_tba_twotier_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let pw_2tier = dir.join("pw_2tier");
    std::fs::create_dir_all(&dir).expect("mkdir");
    let s = build_two_tier_fixture(&dir);

    let tba_bin = env!("CARGO_BIN_EXE_horizon_tba");
    let pw_s = pw_2tier.to_string_lossy();
    let out_2t = dir.join("out_2tier").to_string_lossy().into_owned();
    let stdout_2t = run_bin_stdout(
        tba_bin,
        &[
            "-q",
            "--hinges",
            "h1,h2",
            "--pairwise-dir",
            &pw_s,
            &out_2t,
            "1",
            &s[0],
            &s[1],
            &s[2],
            &s[3],
        ],
        "two-tier",
    );
    assert!(
        stdout_2t.contains("mode=tba_two_tier hinges=2 tips=2 species=4")
            && !stdout_2t.contains("core_blocks=0"),
        "unexpected two-tier stdout: {stdout_2t}"
    );
    let n_pw = std::fs::read_dir(&pw_2tier)
        .expect("read pw dir")
        .filter_map(Result::ok)
        .count();
    assert_eq!(n_pw, 3, "expected exactly 3 two-tier pairwise MAFs");

    // Step 2: Build 2-hinge backbone alone, then add t1 and t2 POST-HOC via --extend-maf.
    let out_bb = dir.join("out_backbone").to_string_lossy().into_owned();
    run_bin_stdout(tba_bin, &["-q", &out_bb, "1", &s[0], &s[1]], "backbone");
    let tip_map = dir.join("tips.tsv");
    std::fs::write(&tip_map, "tip\thinge\nt1\th1\nt2\th2\n").expect("write tip map");
    let bb_maf = dir.join("out_backbone.maf").to_string_lossy().into_owned();
    let tip_map_s = tip_map.to_string_lossy().into_owned();
    let out_ext = dir.join("out_extended").to_string_lossy().into_owned();
    let stdout_ext = run_bin_stdout(
        tba_bin,
        &[
            "-q",
            "--extend-maf",
            &bb_maf,
            "--tip-map",
            &tip_map_s,
            &out_ext,
            "1",
            &s[0],
            &s[1],
            &s[2],
            &s[3],
        ],
        "extend-maf",
    );
    assert!(
        stdout_ext.contains("mode=tba_extend hinges=2 tips=2 species=4")
            && !stdout_ext.contains("core_blocks=0"),
        "unexpected extend stdout: {stdout_ext}"
    );
    let maf_2tier = std::fs::read_to_string(dir.join("out_2tier.maf")).expect("read 2tier");
    let maf_ext = std::fs::read_to_string(dir.join("out_extended.maf")).expect("read extended");
    assert_eq!(maf_2tier, maf_ext, "one-shot and post-hoc MAF parity");

    // Step 3: Post-hoc --extend-maf combined with --from-mafs matches byte-for-byte.
    let out_fm = dir
        .join("out_extended_from_mafs")
        .to_string_lossy()
        .into_owned();
    run_bin_stdout(
        tba_bin,
        &[
            "-q",
            "--extend-maf",
            &bb_maf,
            "--from-mafs",
            &pw_s,
            "--tip-map",
            &tip_map_s,
            &out_fm,
            "h1",
            "h2",
            "t1",
            "t2",
        ],
        "extend-maf --from-mafs",
    );
    let maf_fm =
        std::fs::read_to_string(dir.join("out_extended_from_mafs.maf")).expect("read ext_fm");
    assert_eq!(maf_ext, maf_fm);

    let gfa_p = dir.join("out_extended.gfa").to_string_lossy().into_owned();
    let gfa_stdout = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_gfa"),
        &["stats", &gfa_p],
        "horizon_gfa stats",
    );
    assert!(gfa_stdout.contains("walks=4"));
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn horizon_tba_pairwise_gz_lazy_merge_matches_uncompressed_maf() {
    let dir = std::env::temp_dir().join(format!("horizon_tba_pw_gz_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let pw_plain = dir.join("pw_plain");
    let pw_gz = dir.join("pw_gz");
    std::fs::create_dir_all(&dir).expect("mkdir");

    let anc = synthetic_dna(0x2468_ACE0, 2_000);
    let sp1 = mutate_dna(&anc, 21, 2);
    let sp2 = mutate_dna(&anc, 22, 2);
    let sp3 = mutate_dna(&anc, 23, 3);

    let fa1 = dir.join("sp1.fa");
    let fa2 = dir.join("sp2.fa");
    let fa3 = dir.join("sp3.fa");
    std::fs::write(&fa1, format!(">chr1\n{}\n", String::from_utf8_lossy(&sp1))).expect("fa1");
    std::fs::write(&fa2, format!(">chr1\n{}\n", String::from_utf8_lossy(&sp2))).expect("fa2");
    std::fs::write(&fa3, format!(">chr1\n{}\n", String::from_utf8_lossy(&sp3))).expect("fa3");

    let tba_bin = env!("CARGO_BIN_EXE_horizon_tba");
    let s1 = format!("sp1={}", fa1.display());
    let s2 = format!("sp2={}", fa2.display());
    let s3 = format!("sp3={}", fa3.display());

    let out_plain = dir.join("out_plain");
    run_bin_stdout(
        tba_bin,
        &[
            "-q",
            "--tree",
            "((sp1,sp2),sp3);",
            "--pairwise-dir",
            pw_plain.to_str().expect("pw_plain"),
            out_plain.to_str().expect("out_plain"),
            "1",
            &s1,
            &s2,
            &s3,
        ],
        "tba plain pairwise-dir",
    );

    let out_gz = dir.join("out_gz");
    run_bin_stdout(
        tba_bin,
        &[
            "-q",
            "--tree",
            "((sp1,sp2),sp3);",
            "--pairwise-dir",
            pw_gz.to_str().expect("pw_gz"),
            "--pairwise-gz",
            out_gz.to_str().expect("out_gz"),
            "1",
            &s1,
            &s2,
            &s3,
        ],
        "tba gzip pairwise-dir",
    );

    for pair in ["sp1__sp2.maf.gz", "sp1__sp3.maf.gz", "sp2__sp3.maf.gz"] {
        let p = pw_gz.join(pair);
        assert!(p.is_file(), "expected gzipped pairwise MAF {}", p.display());
        let bytes = std::fs::read(&p).expect("read gz");
        assert!(
            bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b,
            "expected gzip magic header in {}",
            p.display()
        );
    }

    let maf_plain = std::fs::read_to_string(dir.join("out_plain.maf")).expect("plain maf");
    let maf_gz = std::fs::read_to_string(dir.join("out_gz.maf")).expect("gz maf");
    assert_eq!(
        maf_plain, maf_gz,
        "lazy gzip --pairwise-dir must produce byte-identical MAF to plain --pairwise-dir"
    );

    let out_from_gz = dir.join("out_from_gz");
    run_bin_stdout(
        tba_bin,
        &[
            "-q",
            "--tree",
            "((sp1,sp2),sp3);",
            "--from-mafs",
            pw_gz.to_str().expect("pw_gz"),
            out_from_gz.to_str().expect("out_from_gz"),
            "sp1",
            "sp2",
            "sp3",
        ],
        "tba --from-mafs gz",
    );
    let maf_from_gz = std::fs::read_to_string(dir.join("out_from_gz.maf")).expect("from gz maf");
    assert_eq!(maf_plain, maf_from_gz);

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn horizon_tba_progressive_ancestors_emits_maf_gfa_ancestors_and_hal() {
    let dir = std::env::temp_dir().join(format!("horizon_tba_prog_anc_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let anc_dir = dir.join("ancestors");
    std::fs::create_dir_all(&dir).expect("mkdir");

    let anc = synthetic_dna(0x9876_5432, 2_000);
    let sp1 = mutate_dna(&anc, 31, 2);
    let sp2 = mutate_dna(&anc, 32, 2);
    let sp3 = mutate_dna(&anc, 33, 3);
    let sp4 = mutate_dna(&anc, 34, 4);

    let mut specs = Vec::new();
    for (name, seq) in [("sp1", sp1), ("sp2", sp2), ("sp3", sp3), ("sp4", sp4)] {
        let p = dir.join(format!("{name}.fa"));
        std::fs::write(&p, format!(">chr1\n{}\n", String::from_utf8_lossy(&seq)))
            .expect("write fa");
        specs.push(format!("{name}={}", p.display()));
    }

    let out_prefix = dir.join("prog_out");
    let stdout = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_tba"),
        &[
            "-q",
            "--progressive-ancestors",
            "--max-outgroups",
            "2",
            "--emit-ancestors",
            anc_dir.to_str().expect("anc_dir"),
            "--hal",
            "--tree",
            "((sp1,sp2),(sp3,sp4));",
            out_prefix.to_str().expect("out_prefix"),
            "1",
            &specs[0],
            &specs[1],
            &specs[2],
            &specs[3],
        ],
        "horizon_tba --progressive-ancestors",
    );
    assert!(
        stdout.contains("mode=tba_progressive_ancestors")
            && stdout.contains("core_blocks=")
            && !stdout.contains("core_blocks=0"),
        "expected progressive ancestors summary with >=1 4-sp core block: {stdout}"
    );
    assert!(dir.join("prog_out.maf").is_file());
    assert!(dir.join("prog_out.gfa").is_file());
    assert!(dir.join("prog_out.hal").is_file());
    assert!(anc_dir.join("ancestors_summary.tsv").is_file());

    let gfa_stats = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_gfa"),
        &["stats", dir.join("prog_out.gfa").to_str().expect("gfa")],
        "horizon_gfa stats progressive",
    );
    assert!(gfa_stats.contains("walks=4"), "gfa_stats: {gfa_stats}");

    let hal_info = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_wga"),
        &[
            "--hal-info",
            dir.join("prog_out.hal").to_str().expect("hal"),
        ],
        "horizon_wga --hal-info progressive",
    );
    assert!(
        hal_info.contains("hdf5_v=0")
            && hal_info.contains("root_members=11")
            && hal_info.contains("Anc0,Anc1,Anc2"),
        "expected valid HDF5 HAL with reconstructed Anc0..Anc2 genomes, got:\n{hal_info}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn horizon_tba_auto_hinges_eisenstein_runs_end_to_end() {
    let dir = std::env::temp_dir().join(format!(
        "horizon_tba_auto_hinges_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("mkdir");

    let anc = synthetic_dna(0x1357_9BDF, 2_000);
    let clade_a = mutate_dna(&anc, 1, 2);
    let clade_b = mutate_dna(&anc, 2, 12);
    let seqs = [
        ("a1", mutate_dna(&clade_a, 11, 1)),
        ("a2", mutate_dna(&clade_a, 12, 1)),
        ("b1", mutate_dna(&clade_b, 21, 1)),
        ("b2", mutate_dna(&clade_b, 22, 1)),
    ];
    let mut specs = Vec::new();
    for (sp, seq) in &seqs {
        let fa = dir.join(format!("{sp}.fa"));
        std::fs::write(&fa, format!(">chr1\n{}\n", String::from_utf8_lossy(seq))).expect("write");
        specs.push(format!("{sp}={}", fa.display()));
    }

    let out_prefix = dir.join("auto_out");
    let stdout = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_tba"),
        &[
            "-q",
            "--auto-hinges",
            "2",
            out_prefix.to_str().expect("out_prefix"),
            "1",
            &specs[0],
            &specs[1],
            &specs[2],
            &specs[3],
        ],
        "horizon_tba --auto-hinges 2",
    );
    assert!(
        stdout.contains("mode=tba_two_tier hinges=2 tips=2 species=4")
            && !stdout.contains("core_blocks=0"),
        "expected Eisenstein --auto-hinges 2 two-tier alignment with core blocks: {stdout}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
#[allow(clippy::too_many_lines)]
fn horizon_tba_variation_gfa_out_gz_nested_vcf_qc_from_maf_and_hal_e2e() {
    let dir = std::env::temp_dir().join(format!(
        "horizon_vargfa_hal_e2e_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("mkdir");

    let anc = synthetic_dna(0xCAFE_1122, 1_800);
    let mut sp2 = mutate_dna(&anc, 2, 1);
    let mut sp3 = mutate_dna(&anc, 3, 2);
    // Introduce a 30-bp deletion in sp3 spanning an internal SNP in sp2 (and a flanking SNP so
    // the spanning deletion breaks Tarjan bridges across sp2's internal SNP, triggering LV=1).
    sp2[610] = if sp2[610] == b'A' { b'G' } else { b'A' };
    sp3.drain(600..630);
    let sp4 = anc.clone();

    let mut specs = Vec::new();
    for (name, seq) in [("sp1", &anc), ("sp2", &sp2), ("sp3", &sp3), ("sp4", &sp4)] {
        let fa = dir.join(format!("{name}.fa"));
        std::fs::write(&fa, format!(">chr1\n{}\n", String::from_utf8_lossy(seq))).expect("write");
        specs.push(format!("{name}={}", fa.display()));
    }

    let out_prefix = dir.join("var_out");
    let stdout = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_tba"),
        &[
            "-q",
            "--variation-gfa",
            "--out-gz",
            "--tantan",
            "--progressive-ancestors",
            "--hal",
            "--tree",
            "(((sp1,sp4),sp2),sp3);",
            out_prefix.to_str().expect("out_prefix"),
            "1",
            &specs[0],
            &specs[1],
            &specs[2],
            &specs[3],
        ],
        "horizon_tba --variation-gfa --out-gz --tantan --progressive-ancestors --hal",
    );
    assert!(stdout.contains("mode=tba_progressive_ancestors"));
    let maf_gz = dir.join("var_out.maf.gz");
    let gfa_gz = dir.join("var_out.gfa.gz");
    let hal_file = dir.join("var_out.hal");
    assert!(maf_gz.is_file() && gfa_gz.is_file() && hal_file.is_file());

    // 1. horizon_gfa qc on gzipped variation GFA, plus MIXED_PW, P_only, and --max-segments 1 budget
    let qc_out = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_gfa"),
        &["qc", gfa_gz.to_str().expect("gfa.gz")],
        "horizon_gfa qc",
    );
    assert!(
        qc_out.contains("dialect=W_only") && qc_out.contains("placeholder_segments=0"),
        "qc_out:\n{qc_out}"
    );
    let qc_budget = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_gfa"),
        &[
            "qc",
            gfa_gz.to_str().expect("gfa.gz"),
            "--max-segments",
            "1",
        ],
        "horizon_gfa qc --max-segments 1",
    );
    assert!(
        qc_budget.contains("skipped_over_budget"),
        "expected skipped_over_budget under --max-segments 1:\n{qc_budget}"
    );
    let mixed_gfa = dir.join("mixed.gfa");
    std::fs::write(
        &mixed_gfa,
        "H\tVN:Z:1.1\nS\t1\tACGT\nS\t2\tTTTT\nL\t1\t+\t2\t+\t0M\nP\ts1#0#c1\t1+,2+\t*\nW\ts2\t0\tc1\t0\t8\t>1>2\n",
    )
    .expect("write mixed");
    let qc_mixed = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_gfa"),
        &["qc", mixed_gfa.to_str().expect("mixed")],
        "horizon_gfa qc mixed",
    );
    assert!(
        qc_mixed.contains("dialect=MIXED_PW"),
        "expected MIXED_PW:\n{qc_mixed}"
    );
    let ponly_gfa = dir.join("ponly.gfa");
    std::fs::write(
        &ponly_gfa,
        "H\tVN:Z:1.0\nS\t1\tACGT\nS\t2\tTTTT\nL\t1\t+\t2\t+\t0M\nP\ts1#0#c1\t1+,2+\t*\n",
    )
    .expect("write ponly");
    let qc_ponly = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_gfa"),
        &["qc", ponly_gfa.to_str().expect("ponly")],
        "horizon_gfa qc ponly",
    );
    assert!(
        qc_ponly.contains("dialect=P_only"),
        "expected P_only:\n{qc_ponly}"
    );

    // 2. horizon_gfa from-maf reproduces variation GFA directly from var_out.maf.gz
    let from_maf_gfa = dir.join("from_maf.gfa");
    let _ = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_gfa"),
        &[
            "from-maf",
            maf_gz.to_str().expect("maf.gz"),
            "--out",
            from_maf_gfa.to_str().expect("from_maf"),
        ],
        "horizon_gfa from-maf",
    );
    assert!(from_maf_gfa.is_file());

    // 3. horizon_resolve --vcf --nested-snarls calls both LV=0 and LV=1 variants
    let vcf_out = dir.join("nested.vcf");
    let resolved_out = dir.join("resolved.gfa");
    let _ = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_resolve"),
        &[
            gfa_gz.to_str().expect("gfa.gz"),
            "--out",
            resolved_out.to_str().expect("resolved"),
            "--vcf",
            vcf_out.to_str().expect("vcf"),
            "--vcf-ref",
            "sp1",
            "--nested-snarls",
        ],
        "horizon_resolve --nested-snarls",
    );
    let vcf_text = std::fs::read_to_string(&vcf_out).expect("read vcf");
    assert!(
        vcf_text.contains("##fileformat=VCFv4.2")
            && vcf_text.lines().any(|l| !l.starts_with('#'))
            && vcf_text.contains(";LV=1;PS=snarl_"),
        "expected LV=1 nested snarl record with PS= parent link inside 30-bp spanning deletion:\n{vcf_text}"
    );

    // 4. horizon_hal stats + liftover + to-maf + to-gfa on the progressive HAL
    let hal_stats = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_hal"),
        &["stats", hal_file.to_str().expect("hal")],
        "horizon_hal stats",
    );
    assert!(
        hal_stats.contains("sp1") && hal_stats.contains("sp3"),
        "hal_stats:\n{hal_stats}"
    );
    let bed_in = dir.join("q.bed");
    let bed_out = dir.join("lifted.bed");
    std::fs::write(&bed_in, "chr1\t100\t400\treg1\t0\t+\n").expect("write bed");
    let _ = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_hal"),
        &[
            "liftover",
            hal_file.to_str().expect("hal"),
            "--src",
            "sp1",
            "--dst",
            "sp3",
            "--bed",
            bed_in.to_str().expect("bed_in"),
            "--out",
            bed_out.to_str().expect("bed_out"),
        ],
        "horizon_hal liftover",
    );
    assert!(
        !std::fs::read_to_string(&bed_out)
            .expect("read lifted bed")
            .trim()
            .is_empty(),
        "expected non-empty multi-level HAL liftover"
    );
    let hal_gfa = dir.join("from_hal.gfa");
    let _ = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_hal"),
        &[
            "to-gfa",
            hal_file.to_str().expect("hal"),
            "--ref",
            "sp1",
            "--out",
            hal_gfa.to_str().expect("hal_gfa"),
        ],
        "horizon_hal to-gfa",
    );
    let hal_gfa_text = std::fs::read_to_string(&hal_gfa).expect("read hal_gfa");
    assert!(
        hal_gfa_text.contains("S\t") && hal_gfa_text.contains("W\tsp1\t"),
        "expected variation GFA from horizon_hal to-gfa:\n{hal_gfa_text}"
    );

    // 5. horizon_wga --tantan --train-hmm (verify trained gap parameters on stderr)
    let wga_res = Command::new(env!("CARGO_BIN_EXE_horizon_wga"))
        .args([
            dir.join("sp1.fa").to_str().expect("sp1"),
            dir.join("sp3.fa").to_str().expect("sp3"),
            "--tantan",
            "--train-hmm",
        ])
        .output()
        .expect("run horizon_wga");
    assert!(
        wga_res.status.success(),
        "horizon_wga failed: {}",
        String::from_utf8_lossy(&wga_res.stderr)
    );
    let wga_out = String::from_utf8_lossy(&wga_res.stdout);
    let wga_err = String::from_utf8_lossy(&wga_res.stderr);
    assert!(wga_out.contains("aligned_pairs="), "wga_out:\n{wga_out}");
    assert!(
        wga_err.contains("[horizon_wga --train-hmm]")
            && wga_err.contains("del=(")
            && wga_err.contains("ins=("),
        "expected trained HMM gap parameters on stderr:\n{wga_err}"
    );

    // 6. horizon_locus --repair-frameshifts on a synthetic CDS with a +1 nt frameshift
    let cds_fa = dir.join("cds.fa");
    let ref_cds = "ATGGCTAGCGAACGTACCGTTCTGGCTAGCGAACGTACCGTTCTGGCTAGCGAACGTACCGTTCTG";
    let mut fs_cds = ref_cds.to_owned();
    fs_cds.insert(24, 'A'); // +1 nt frameshift
    std::fs::write(&cds_fa, format!(">ref_cds\n{ref_cds}\n>fs_cds\n{fs_cds}\n")).expect("cds");
    let loc_pref = dir.join("loc_out");
    let _ = run_bin_stdout(
        env!("CARGO_BIN_EXE_horizon_locus"),
        &[
            "--fasta",
            cds_fa.to_str().expect("cds"),
            "--out-prefix",
            loc_pref.to_str().expect("loc_pref"),
            "--ref",
            "ref_cds",
            "--repair-frameshifts",
        ],
        "horizon_locus --repair-frameshifts",
    );
    let fs_tsv = std::fs::read_to_string(dir.join("loc_out.frameshifts.tsv")).expect("fs tsv");
    assert!(
        fs_tsv.contains("fs_cds"),
        "expected detected frameshift in frameshifts.tsv:\n{fs_tsv}"
    );
    let repaired_fa = std::fs::read_to_string(dir.join("loc_out.repaired.fa")).expect("repaired");
    assert!(
        repaired_fa.contains(">fs_cds"),
        "expected repaired sequence header in loc_out.repaired.fa:\n{repaired_fa}"
    );
    let fs_repaired_seq: String = repaired_fa
        .split(">fs_cds")
        .nth(1)
        .expect("fs_cds record")
        .lines()
        .skip(1)
        .take_while(|l| !l.starts_with('>'))
        .collect();
    assert_eq!(
        fs_repaired_seq.len() % 3,
        0,
        "repaired CDS must be frame-aligned (multiple of 3 bp): {fs_repaired_seq}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}
