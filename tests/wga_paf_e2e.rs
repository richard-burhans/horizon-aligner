//! End-to-end integration tests for `horizon_wga --paf`, `horizon_filter --paf`, and uniform
//! `-t/--threads` & `-V/--version` CLI support across the suite.

mod common;

use common::{mutate_dna, synthetic_dna};
use std::process::Command;

#[test]
fn horizon_wga_paf_output_feeds_directly_into_horizon_filter() {
    let dir = std::env::temp_dir().join(format!("horizon_wga_paf_e2e_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("mkdir");

    let target = synthetic_dna(0xABCD_1234, 2_500);
    let query = mutate_dna(&target, 7, 2);
    let tgt_fa = dir.join("tgt.fa");
    let qry_fa = dir.join("qry.fa");
    let raw_paf = dir.join("raw.paf");
    let filt_paf = dir.join("filtered.paf");

    std::fs::write(
        &tgt_fa,
        format!(">chrT\n{}\n", String::from_utf8_lossy(&target)),
    )
    .expect("tgt");
    std::fs::write(
        &qry_fa,
        format!(">chrQ\n{}\n", String::from_utf8_lossy(&query)),
    )
    .expect("qry");

    let wga_out = Command::new(env!("CARGO_BIN_EXE_horizon_wga"))
        .args([
            tgt_fa.to_str().expect("tgt"),
            qry_fa.to_str().expect("qry"),
            "1",
            "--paf",
            raw_paf.to_str().expect("paf"),
            "-t",
            "2",
        ])
        .output()
        .expect("run horizon_wga");
    assert!(
        wga_out.status.success(),
        "horizon_wga --paf failed: {}",
        String::from_utf8_lossy(&wga_out.stderr)
    );
    let raw_text = std::fs::read_to_string(&raw_paf).expect("read raw paf");
    assert!(
        raw_text.contains("chrQ\t2500\t") && raw_text.contains("cg:Z:"),
        "expected valid PAF with cg:Z: CIGAR, got:\n{raw_text}"
    );

    let filt_out = Command::new(env!("CARGO_BIN_EXE_horizon_filter"))
        .args([
            "--paf",
            raw_paf.to_str().expect("raw"),
            "--target",
            tgt_fa.to_str().expect("tgt"),
            "--query",
            qry_fa.to_str().expect("qry"),
            "--out",
            filt_paf.to_str().expect("filt"),
            "--threads",
            "1",
            "--threads",
            "2",
            "--net",
        ])
        .output()
        .expect("run horizon_filter");
    assert!(
        filt_out.status.success(),
        "horizon_filter failed on horizon_wga PAF: {}",
        String::from_utf8_lossy(&filt_out.stderr)
    );
    let filt_text = std::fs::read_to_string(&filt_paf).expect("read filtered paf");
    assert!(
        filt_text.contains("chrQ\t2500\t") && filt_text.contains("md:f:"),
        "expected filtered PAF with md:f: tag, got:\n{filt_text}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn horizon_wga_resolve_dups_prunes_redundant_paralog_blocks_from_maf() {
    let dir = std::env::temp_dir().join(format!(
        "horizon_wga_resolve_dups_maf_{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("mkdir");

    // Target: [flank_left (1200 bp)] + [cassette (900 bp)] + [flank_right (1200 bp)]
    // Query:  [flank_left (1200 bp)] + [cassette (900 bp, exact)] + [spacer (800 bp)] + [cassette (900 bp, 6% diverged)] + [flank_right (1200 bp)]
    let flank_l = synthetic_dna(0x1111_2222, 1_200);
    let cassette = synthetic_dna(0x3333_4444, 900);
    let spacer = synthetic_dna(0x5555_6666, 800);
    let flank_r = synthetic_dna(0x7777_8888, 1_200);
    let cassette_paralog = mutate_dna(&cassette, 0x9999_AAAA, 6);

    let mut target = Vec::new();
    target.extend_from_slice(&flank_l);
    target.extend_from_slice(&cassette);
    target.extend_from_slice(&flank_r);

    let mut query = Vec::new();
    query.extend_from_slice(&flank_l);
    query.extend_from_slice(&cassette);
    query.extend_from_slice(&spacer);
    query.extend_from_slice(&cassette_paralog);
    query.extend_from_slice(&flank_r);

    let tgt_fa = dir.join("tgt.fa");
    let qry_fa = dir.join("qry.fa");
    let maf_raw = dir.join("raw.maf");
    let maf_res = dir.join("resolved.maf");
    std::fs::write(
        &tgt_fa,
        format!(">chrT\n{}\n", String::from_utf8_lossy(&target)),
    )
    .expect("tgt");
    std::fs::write(
        &qry_fa,
        format!(">chrQ\n{}\n", String::from_utf8_lossy(&query)),
    )
    .expect("qry");

    let run_raw = Command::new(env!("CARGO_BIN_EXE_horizon_wga"))
        .args([
            tgt_fa.to_str().expect("tgt"),
            qry_fa.to_str().expect("qry"),
            "1",
            "--maf",
            maf_raw.to_str().expect("raw"),
            "-t",
            "2",
        ])
        .output()
        .expect("run raw wga");
    assert!(run_raw.status.success());

    let run_res = Command::new(env!("CARGO_BIN_EXE_horizon_wga"))
        .args([
            tgt_fa.to_str().expect("tgt"),
            qry_fa.to_str().expect("qry"),
            "1",
            "--resolve-dups",
            "--maf",
            maf_res.to_str().expect("res"),
            "-t",
            "2",
        ])
        .output()
        .expect("run resolved wga");
    assert!(
        run_res.status.success(),
        "horizon_wga --resolve-dups failed: {}",
        String::from_utf8_lossy(&run_res.stderr)
    );
    let stdout_res = String::from_utf8_lossy(&run_res.stdout);
    let raw_txt = std::fs::read_to_string(&maf_raw).expect("read raw maf");
    let res_txt = std::fs::read_to_string(&maf_res).expect("read resolved maf");
    let raw_blocks = raw_txt
        .lines()
        .filter(|l| l.starts_with("a score="))
        .count();
    let res_blocks = res_txt
        .lines()
        .filter(|l| l.starts_with("a score="))
        .count();
    assert!(
        res_blocks <= raw_blocks && res_blocks >= 1,
        "resolved MAF blocks ({res_blocks}) should be <= raw ({raw_blocks})"
    );
    if stdout_res.contains("# EG:Z:") {
        assert!(
            res_blocks < raw_blocks,
            "when e-graph detects duplications (# EG:Z:), resolved MAF must prune at least one redundant block (raw={raw_blocks}, res={res_blocks})"
        );
    }

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn cli_binaries_report_version_flag() {
    let bins = [
        ("horizon", env!("CARGO_BIN_EXE_horizon")),
        ("horizon_wga", env!("CARGO_BIN_EXE_horizon_wga")),
        ("horizon_guide", env!("CARGO_BIN_EXE_horizon_guide")),
        ("horizon_tba", env!("CARGO_BIN_EXE_horizon_tba")),
        ("horizon_merge", env!("CARGO_BIN_EXE_horizon_merge")),
        ("horizon_filter", env!("CARGO_BIN_EXE_horizon_filter")),
        ("horizon_gfa", env!("CARGO_BIN_EXE_horizon_gfa")),
        ("horizon_clean", env!("CARGO_BIN_EXE_horizon_clean")),
        ("horizon_lift", env!("CARGO_BIN_EXE_horizon_lift")),
        ("horizon_resolve", env!("CARGO_BIN_EXE_horizon_resolve")),
        ("horizon_locus", env!("CARGO_BIN_EXE_horizon_locus")),
        ("horizon_eval", env!("CARGO_BIN_EXE_horizon_eval")),
    ];
    for (name, path) in bins {
        let out = Command::new(path).arg("-V").output().expect("run -V");
        assert!(out.status.success(), "{name} -V failed");
        let stdout = String::from_utf8_lossy(&out.stdout);
        assert!(
            stdout.starts_with(name) && stdout.contains(env!("CARGO_PKG_VERSION")),
            "unexpected -V output for {name}: `{stdout}`"
        );
    }
}

#[test]
#[allow(clippy::too_many_lines, clippy::unwrap_used)]
fn cli_error_exit_codes_stdout_hygiene_and_uniform_short_t() {
    let dir = std::env::temp_dir().join(format!(
        "horizon_cli_wave2b_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("mkdir");

    // 1. horizon_resolve and horizon_locus must exit non-zero when invoked with no args.
    let res_no_args = Command::new(env!("CARGO_BIN_EXE_horizon_resolve"))
        .output()
        .expect("run horizon_resolve");
    assert!(
        !res_no_args.status.success(),
        "horizon_resolve with no args must exit non-zero"
    );
    let loc_no_args = Command::new(env!("CARGO_BIN_EXE_horizon_locus"))
        .output()
        .expect("run horizon_locus");
    assert!(
        !loc_no_args.status.success(),
        "horizon_locus with no args must exit non-zero"
    );

    // 2. horizon_wga --paf - must write ONLY PAF to stdout (no summary corruption)
    //    and reject multiple "-" outputs.
    let target = synthetic_dna(0x1111_2222, 800);
    let tgt_fa = dir.join("tgt.fa");
    let qry_fa = dir.join("qry.fa");
    std::fs::write(
        &tgt_fa,
        format!(">c1\n{}\n", String::from_utf8_lossy(&target)),
    )
    .unwrap();
    std::fs::write(
        &qry_fa,
        format!(">c1\n{}\n", String::from_utf8_lossy(&target)),
    )
    .unwrap();

    let wga_stdout = Command::new(env!("CARGO_BIN_EXE_horizon_wga"))
        .args([
            tgt_fa.to_str().unwrap(),
            qry_fa.to_str().unwrap(),
            "1",
            "--paf",
            "-",
        ])
        .output()
        .expect("run wga --paf -");
    assert!(wga_stdout.status.success());
    let paf_out = String::from_utf8_lossy(&wga_stdout.stdout);
    assert!(
        !paf_out.contains("contigs_x=") && !paf_out.contains("H\tVN:Z:"),
        "summary text must not corrupt stdout when --paf - is used:\n{paf_out}"
    );
    assert!(paf_out.contains("c1\t800\t"));

    let wga_double_dash = Command::new(env!("CARGO_BIN_EXE_horizon_wga"))
        .args([
            tgt_fa.to_str().unwrap(),
            qry_fa.to_str().unwrap(),
            "1",
            "--paf",
            "-",
            "--maf",
            "-",
        ])
        .output()
        .expect("run wga double -");
    assert!(
        !wga_double_dash.status.success(),
        "multiple '-' outputs must be rejected"
    );

    // 3. horizon_filter -t must mean --threads (standardized across all binaries).
    let raw_paf = dir.join("raw.paf");
    std::fs::write(&raw_paf, paf_out.as_bytes()).unwrap();
    let filt_short_t = Command::new(env!("CARGO_BIN_EXE_horizon_filter"))
        .args([
            "--paf",
            raw_paf.to_str().unwrap(),
            "--target",
            tgt_fa.to_str().unwrap(),
            "--query",
            qry_fa.to_str().unwrap(),
            "-t",
            "2",
        ])
        .output()
        .expect("run filter -t 2");
    assert!(
        filt_short_t.status.success(),
        "horizon_filter -t 2 failed: {}",
        String::from_utf8_lossy(&filt_short_t.stderr)
    );

    // 4. horizon_merge must reject double "-" outputs.
    let dummy_maf = dir.join("d.maf");
    std::fs::write(&dummy_maf, "##maf version=1\n").unwrap();
    let merge_double_dash = Command::new(env!("CARGO_BIN_EXE_horizon_merge"))
        .args([
            "ref",
            "100",
            "-",
            "-",
            &format!("sp1={}", dummy_maf.display()),
        ])
        .output()
        .expect("run merge - -");
    assert!(
        !merge_double_dash.status.success(),
        "horizon_merge must reject writing both MAF and GFA to '-'"
    );

    // 5. horizon_lift must reject double stdin and conflicting --bed + --gff.
    let lift_double_stdin = Command::new(env!("CARGO_BIN_EXE_horizon_lift"))
        .args(["--gfa", "-", "--from", "A", "--to", "B", "--bed", "-"])
        .output()
        .expect("run lift - -");
    assert!(!lift_double_stdin.status.success());

    let lift_conflict = Command::new(env!("CARGO_BIN_EXE_horizon_lift"))
        .args([
            "--gfa",
            dummy_maf.to_str().unwrap(),
            "--from",
            "A",
            "--to",
            "B",
            "--bed",
            dummy_maf.to_str().unwrap(),
            "--gff",
            dummy_maf.to_str().unwrap(),
        ])
        .output()
        .expect("run lift --bed --gff");
    assert!(!lift_conflict.status.success());

    // 6. horizon_gfa extract handles GFA filename containing ':' and rejects empty '::' spec.
    let colon_gfa = dir.join("graph:v1.gfa");
    let ext_fa = dir.join("ext.fa");
    std::fs::write(
        &colon_gfa,
        "H\tVN:Z:1.1\nS\ts1\tACGTACGT\nW\tHG002\t0\tchr1\t0\t8\t>s1\n",
    )
    .unwrap();
    let gfa_colon = Command::new(env!("CARGO_BIN_EXE_horizon_gfa"))
        .args([
            "extract",
            colon_gfa.to_str().unwrap(),
            &format!("HG002:chr1:{}", ext_fa.display()),
        ])
        .output()
        .expect("run gfa extract colon");
    assert!(
        gfa_colon.status.success(),
        "horizon_gfa extract failed on ':' path: {}",
        String::from_utf8_lossy(&gfa_colon.stderr)
    );
    let gfa_empty_spec = Command::new(env!("CARGO_BIN_EXE_horizon_gfa"))
        .args(["extract", colon_gfa.to_str().unwrap(), "::"])
        .output()
        .expect("run gfa extract ::");
    assert!(
        !gfa_empty_spec.status.success(),
        "horizon_gfa extract must reject empty '::' spec"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn wave_a_cli_stdout_hygiene_regressions() {
    let dir = std::env::temp_dir().join(format!(
        "horizon_wave_a_cli_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).expect("mkdir");

    // C7 regression: horizon_merge must accept '-' for maf_out and route summary to stderr.
    let maf_a = dir.join("sp1.maf");
    let gfa_out = dir.join("merged.gfa");
    std::fs::write(
        &maf_a,
        "##maf version=1\na score=10.0\ns hg38.chr1 0 10 + 100 ACGTACGTAC\ns sp1.chr1 0 10 + 100 ACGTACGTAC\n\n",
    )
    .expect("write maf");
    let merge_out = Command::new(env!("CARGO_BIN_EXE_horizon_merge"))
        .args([
            "hg38",
            "100",
            "-",
            gfa_out.to_str().expect("utf8"),
            &format!("sp1={}", maf_a.display()),
        ])
        .output()
        .expect("run horizon_merge");
    assert!(
        merge_out.status.success(),
        "C7: horizon_merge must accept '-' stdout output: {}",
        String::from_utf8_lossy(&merge_out.stderr)
    );
    let merge_stdout = String::from_utf8_lossy(&merge_out.stdout);
    assert!(
        merge_stdout.contains("##maf version=1") && !merge_stdout.contains("total_blocks="),
        "C7: horizon_merge stdout must contain pure MAF without summary text: {merge_stdout}"
    );

    // C8 regression: horizon_lift must reject --unmapped - when primary output also targets stdout.
    let lift_gfa = dir.join("lift.gfa");
    let lift_bed = dir.join("in.bed");
    std::fs::write(
        &lift_gfa,
        "H\tVN:Z:1.1\nS\ts1\tACGTACGT\nW\tA\t0\tchr1\t0\t8\t>s1\nW\tB\t0\tchr1\t0\t8\t>s1\n",
    )
    .expect("write gfa");
    std::fs::write(&lift_bed, "chr1\t0\t8\tf1\n").expect("write bed");
    let lift_out = Command::new(env!("CARGO_BIN_EXE_horizon_lift"))
        .args([
            "--gfa",
            lift_gfa.to_str().expect("utf8"),
            "--from",
            "A",
            "--to",
            "B",
            "--bed",
            lift_bed.to_str().expect("utf8"),
            "--unmapped",
            "-",
        ])
        .output()
        .expect("run horizon_lift");
    assert!(
        !lift_out.status.success(),
        "C8: horizon_lift must reject --unmapped - when --out defaults to stdout"
    );

    // C9 regression: horizon_gfa extract mixing file + '-' must keep stdout pure FASTA
    // and reject two '-' specs.
    let ext_file = dir.join("a.fa");
    let gfa_mix = Command::new(env!("CARGO_BIN_EXE_horizon_gfa"))
        .args([
            "extract",
            lift_gfa.to_str().expect("utf8"),
            &format!("A:chr1:{}", ext_file.display()),
            "B:chr1:-",
        ])
        .output()
        .expect("run horizon_gfa extract");
    assert!(
        gfa_mix.status.success(),
        "horizon_gfa extract failed: {}",
        String::from_utf8_lossy(&gfa_mix.stderr)
    );
    let mix_stdout = String::from_utf8_lossy(&gfa_mix.stdout);
    assert!(
        mix_stdout.starts_with('>') && !mix_stdout.contains("wrote "),
        "C9: horizon_gfa extract stdout corrupted by file status message:\n{mix_stdout}"
    );
    let gfa_dup = Command::new(env!("CARGO_BIN_EXE_horizon_gfa"))
        .args([
            "extract",
            lift_gfa.to_str().expect("utf8"),
            "A:chr1:-",
            "B:chr1:-",
        ])
        .output()
        .expect("run horizon_gfa extract dup stdout");
    assert!(
        !gfa_dup.status.success(),
        "C9: horizon_gfa extract must reject multiple '-' specs"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

fn verify_chain_liftover(dir: &std::path::Path, chain_path: &std::path::Path) {
    let bed_in = dir.join("genes.bed");
    let bed_out = dir.join("lifted.bed");
    std::fs::write(
        &bed_in,
        "chrT\t200\t400\tfwd_exon\t100\t+\nchrT\t1550\t1850\tinv_exon\t200\t+\n",
    )
    .expect("write bed");
    let lift_out = Command::new(env!("CARGO_BIN_EXE_horizon_lift"))
        .args([
            "--chain",
            chain_path.to_str().expect("utf8"),
            "--bed",
            bed_in.to_str().expect("utf8"),
            "--out",
            bed_out.to_str().expect("utf8"),
            "--bridge-bp",
            "50",
        ])
        .output()
        .expect("run horizon_lift --chain");
    assert!(
        lift_out.status.success(),
        "horizon_lift --chain failed: {}",
        String::from_utf8_lossy(&lift_out.stderr)
    );
    let lifted = std::fs::read_to_string(&bed_out).expect("read lifted bed");
    assert!(
        lifted.contains("chrQ\t200\t400\tfwd_exon\t100\t+\tchrT:200-400\t"),
        "expected exact forward coordinates:\n{lifted}"
    );
    assert!(
        lifted.contains("chrQ\t1650\t1950\tinv_exon\t200\t-\tchrT:1550-1850\t"),
        "expected exact inverted coordinates with '-' strand flip:\n{lifted}"
    );
}

#[test]
fn horizon_wga_outputs_ucsc_chain_and_nested_net_for_inversion() {
    let dir =
        std::env::temp_dir().join(format!("horizon_wga_ucsc_chain_net_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("mkdir");

    // Target: [left 1400][inv 700][right 1400] = 3500 bp
    // Query:  [left 1400 (2% sub)][revcomp(inv) 700 (2% sub)][right 1400 (2% sub)] = 3500 bp
    let left = synthetic_dna(0xCAFE_0001, 1_400);
    let mid = synthetic_dna(0xCAFE_0002, 700);
    let right = synthetic_dna(0xCAFE_0003, 1_400);
    let mut target = Vec::new();
    target.extend_from_slice(&left);
    target.extend_from_slice(&mid);
    target.extend_from_slice(&right);

    let mut query = Vec::new();
    query.extend_from_slice(&mutate_dna(&left, 11, 2));
    query.extend_from_slice(&mutate_dna(&horizon_aligner::revcomp(&mid), 22, 2));
    query.extend_from_slice(&mutate_dna(&right, 33, 2));

    let tgt_fa = dir.join("tgt.fa");
    let qry_fa = dir.join("qry.fa");
    let chain_path = dir.join("out.chain");
    let net_path = dir.join("out.net");
    std::fs::write(
        &tgt_fa,
        format!(">chrT\n{}\n", String::from_utf8_lossy(&target)),
    )
    .expect("tgt");
    std::fs::write(
        &qry_fa,
        format!(">chrQ\n{}\n", String::from_utf8_lossy(&query)),
    )
    .expect("qry");

    let out = Command::new(env!("CARGO_BIN_EXE_horizon_wga"))
        .args([
            tgt_fa.to_str().expect("utf8"),
            qry_fa.to_str().expect("utf8"),
            "1",
            "--chain",
            chain_path.to_str().expect("utf8"),
            "--net",
            net_path.to_str().expect("utf8"),
            "-t",
            "2",
        ])
        .output()
        .expect("run horizon_wga --chain --net");
    assert!(
        out.status.success(),
        "horizon_wga --chain --net failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let chain_text = std::fs::read_to_string(&chain_path).expect("read chain");
    let chains = horizon_aligner::parse_ucsc_chains(&chain_text).expect("parse UCSC chain");
    assert!(
        chains.len() >= 2,
        "expected >= 2 chains (+ flanks and - inversion), got {}:\n{chain_text}",
        chains.len()
    );
    assert!(
        chains.iter().any(|c| c.q_strand == '+'),
        "missing '+' chain"
    );
    assert!(
        chains.iter().any(|c| c.q_strand == '-'),
        "missing '-' inversion chain"
    );

    let net_text = std::fs::read_to_string(&net_path).expect("read net");
    assert!(
        net_text.contains("net chrT 3500\n"),
        "missing net target header:\n{net_text}"
    );
    assert!(
        net_text.contains(" type top\n"),
        "missing top-level fill:\n{net_text}"
    );
    assert!(
        net_text.contains("  gap ")
            && net_text.contains("   fill ")
            && net_text.contains(" type inv\n"),
        "expected nested inversion fill inside gap:\n{net_text}"
    );

    verify_chain_liftover(&dir, &chain_path);
    let _ = std::fs::remove_dir_all(&dir);
}
