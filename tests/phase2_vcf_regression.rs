//! Phase 2 VCF snarl-nesting regression tests (`INFO/PS` parent linkage and orphaned `LV=1` promotion).

use horizon_aligner::io_util::open_text_reader;
use std::collections::HashMap;
use std::io::Read;
use std::path::PathBuf;
use std::process::Command;

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("horizon_p2_{tag}_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("mkdir");
    dir
}

const NESTED_AND_ORPHAN_GFA: &str = "\
H\tVN:Z:1.1
S\ts1\tAAACCCGGGTTT
S\ts2\tGGTTAA
S\ts3\tA
S\ts4\tG
S\ts5\tCCAATT
S\ts6\tTTTTAAAACCCC
S\ts7\tCCCCCCCCCC
S\ts8\tGGAATT
S\ts9\tC
S\ts10\tT
S\ts11\tTTGGAA
S\ts12\tACGTACGTACGT
S\ts13\tGGGGGG
L\ts1\t+\ts2\t+\t0M
L\ts2\t+\ts3\t+\t0M
L\ts2\t+\ts4\t+\t0M
L\ts3\t+\ts5\t+\t0M
L\ts4\t+\ts5\t+\t0M
L\ts5\t+\ts6\t+\t0M
L\ts1\t+\ts7\t+\t0M
L\ts7\t+\ts6\t+\t0M
L\ts6\t+\ts8\t+\t0M
L\ts8\t+\ts9\t+\t0M
L\ts8\t+\ts10\t+\t0M
L\ts9\t+\ts11\t+\t0M
L\ts10\t+\ts11\t+\t0M
L\ts11\t+\ts12\t+\t0M
L\ts6\t+\ts13\t+\t0M
L\ts13\t+\ts12\t+\t0M
W\tRef\t0\tchr1\t0\t62\t>s1>s2>s3>s5>s6>s8>s9>s11>s12
W\tHG002\t1\tchr1\t0\t50\t>s1>s2>s4>s5>s6>s8>s10>s11
W\tHG002\t2\tchr1\t0\t34\t>s1>s7>s6
";

/// Verify that `horizon_resolve --vcf --nested-snarls`:
/// 1. Declares `##INFO=<ID=PS,...>` in the VCF header.
/// 2. Attaches `PS=snarl_<parent_id>` to every emitted `LV=1` record enclosed by an `LV=0` snarl.
/// 3. Promotes orphaned `LV=1` bubbles (whose enclosing `LV=0` snarl was unemitted) to `LV=0` so
///    `vcfbub -l 0` never silently drops them, guaranteeing zero unlinked `LV=1` records.
#[test]
fn phase2_nested_snarl_ps_linkage_and_orphan_promotion() {
    let dir = temp_dir("nested_ps");
    let gfa_path = dir.join("nested.gfa");
    let out_gfa = dir.join("resolved.gfa");
    let out_vcf = dir.join("nested.vcf.gz");
    std::fs::write(&gfa_path, NESTED_AND_ORPHAN_GFA).expect("write gfa");

    let res = Command::new(env!("CARGO_BIN_EXE_horizon_resolve"))
        .args([
            gfa_path.to_str().expect("utf8"),
            "--out",
            out_gfa.to_str().expect("utf8"),
            "--vcf",
            out_vcf.to_str().expect("utf8"),
            "--vcf-ref",
            "Ref",
            "--nested-snarls",
        ])
        .output()
        .expect("run horizon_resolve");
    assert!(
        res.status.success(),
        "horizon_resolve failed: {}",
        String::from_utf8_lossy(&res.stderr)
    );

    let mut vcf = String::new();
    open_text_reader(&out_vcf)
        .expect("open vcf.gz")
        .read_to_string(&mut vcf)
        .expect("read vcf.gz");

    assert!(
        vcf.contains("##INFO=<ID=PS,Number=1,Type=String,"),
        "VCF header must declare INFO/PS:\n{vcf}"
    );

    let mut lv0_spans: HashMap<String, (usize, usize)> = HashMap::new();
    let mut saw_linked_lv1_snp = false;
    let mut saw_promoted_orphan_snp = false;

    for line in vcf.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        let pos: usize = f[1].parse().expect("pos");
        let id = f[2];
        let ref_a = f[3];
        let alt_a = f[4];
        let info = f[7];

        let mut lv: Option<&str> = None;
        let mut ps: Option<&str> = None;
        for kv in info.split(';') {
            if let Some(v) = kv.strip_prefix("LV=") {
                lv = Some(v);
            } else if let Some(v) = kv.strip_prefix("PS=") {
                ps = Some(v);
            }
        }
        match lv {
            Some("0") => {
                assert!(ps.is_none(), "LV=0 record must not carry PS=: {line}");
                lv0_spans.insert(id.to_owned(), (pos, pos + ref_a.len()));
                // Region B orphaned SNP (C -> T inside s8..s11) must be promoted to LV=0.
                if ref_a.ends_with('C') && alt_a.ends_with('T') && pos >= 30 {
                    saw_promoted_orphan_snp = true;
                }
            }
            Some("1") => {
                let parent = ps.unwrap_or_else(|| {
                    panic!("every LV=1 record MUST carry PS=<parent_id>: {line}\nFull VCF:\n{vcf}")
                });
                let (p_start, p_end) = lv0_spans.get(parent).unwrap_or_else(|| {
                    panic!("PS={parent} must reference an earlier LV=0 record: {line}\n{vcf}")
                });
                assert!(
                    pos >= *p_start && pos < *p_end,
                    "LV=1 pos {pos} must lie inside parent {parent} [{p_start}, {p_end}): {line}"
                );
                // Region A nested SNP (A -> G inside s2..s5) must be linked to Region A LV=0 snarl.
                if ref_a.ends_with('A') && alt_a.ends_with('G') {
                    saw_linked_lv1_snp = true;
                }
            }
            other => panic!("unexpected LV={other:?} in line: {line}"),
        }
    }

    assert!(
        saw_linked_lv1_snp,
        "expected nested Region A SNP (A->G) as LV=1 with PS= parent link:\n{vcf}"
    );
    assert!(
        saw_promoted_orphan_snp,
        "expected orphaned Region B SNP (C->T) promoted to LV=0:\n{vcf}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}
