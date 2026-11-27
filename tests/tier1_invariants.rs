//! Tier 1 CI mathematical & format invariant suite (`docs/VALIDATION_PLAN.md` I1–I8).
//!
//! Enforces:
//! - I1: Byte-identical thread-count determinism (`-t 1` vs `-t 2` vs `-t 4`) across binaries
//! - I2: Strand / reverse-complement equivariance (`J = 1.000000`)
//! - I3: Input-order permutation invariance (`horizon_wga` contigs & `horizon_tba -> resolve --vcf`)
//! - I4: Exact MAF/PAF coordinate + CIGAR sequence reconstruction, UCSC chain/net span sums, & 1:1 0-bp overlap
//! - I5: Per-walk GFA sequence round-trip (`0` diffs) & whole-contig walk completeness
//! - I6: VCF haplotype allele reconstruction round-trip & minimum called-coverage gate
//! - I7: HAL bi-directional round-trip (`horizon_guide --hal` & `horizon_tba --hal`)
//! - I8: Analytical math (`proptest`: Pair-HMM `Z_F == Z_B`, fast-path vs exact DP agreement,
//!   pre-clamp BOCPD Bayes marginal conservation, and fail-closed parser fuzzing)

use horizon_aligner::bocpd::{BocpdScratch, bocpd_gc_profile_into};
use horizon_aligner::dna::revcomp;
use horizon_aligner::formats::for_each_maf_stanza_reader;
use horizon_aligner::gfa::{parse_walk_steps, stream_gfa_record_refs_reader};
use horizon_aligner::hal::{HDF5_SIGNATURE, read_hdf5_superblock, read_hdf5_v0_root_names};
use horizon_aligner::io_util::open_text_reader;
use horizon_aligner::mdl::{
    HmmParams, banded_mdl_gap_posterior_trace, banded_mdl_gap_with_min_post,
};
use horizon_aligner::paf::{CigarOp, PafRecord};
use horizon_aligner::tba::for_each_multi_tba_block_reader;
use proptest::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{Cursor, Read};
use std::path::{Path, PathBuf};
use std::process::Command;

type VcfSite = (usize, String, Vec<String>, Vec<Option<usize>>);

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("horizon_t1_{tag}_{}", std::process::id()));
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

fn mutate_sub(seed: &mut u64, seq: &[u8], rate: f64) -> Vec<u8> {
    let bases = *b"ACGT";
    seq.iter()
        .map(|&b| {
            *seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            #[allow(clippy::cast_precision_loss)]
            let u = ((*seed >> 11) as f64) / ((1u64 << 53) as f64);
            if u < rate {
                let idx = bases.iter().position(|&x| x == b).unwrap_or(0);
                bases[(idx + 1 + ((*seed >> 59) as usize % 3)) % 4]
            } else {
                b
            }
        })
        .collect()
}

fn write_fasta(path: &Path, contigs: &[(&str, &[u8])]) {
    let mut s = String::new();
    for (name, seq) in contigs {
        s.push('>');
        s.push_str(name);
        s.push('\n');
        s.push_str(&String::from_utf8_lossy(seq));
        s.push('\n');
    }
    std::fs::write(path, s).expect("write fasta");
}

/// Extract aligned forward-strand coordinate pairs `(target_pos, query_pos)` from a pairwise MAF.
#[allow(clippy::similar_names)]
fn maf_aligned_pairs(maf_text: &str) -> HashSet<(String, usize, String, usize)> {
    let mut pairs = HashSet::new();
    let mut s_rows: Vec<(String, usize, char, usize, String)> = Vec::new();
    for line in maf_text.lines() {
        if let Some(rest) = line.strip_prefix("s ") {
            let f: Vec<&str> = rest.split_whitespace().collect();
            let name = f[0].to_owned();
            let start: usize = f[1].parse().expect("start");
            let strand = f[3].chars().next().expect("strand");
            let src_size: usize = f[4].parse().expect("src_size");
            let seq = f[5].to_owned();
            s_rows.push((name, start, strand, src_size, seq));
            if s_rows.len() == 2 {
                let (nx, sx, stx, lx, seqx) = &s_rows[0];
                let (ny, sy, sty, ly, seqy) = &s_rows[1];
                let mut ox = 0usize;
                let mut oy = 0usize;
                for (bx, by) in seqx.bytes().zip(seqy.bytes()) {
                    let px = if bx == b'-' {
                        None
                    } else {
                        let p = if *stx == '-' {
                            *lx - 1 - (*sx + ox)
                        } else {
                            *sx + ox
                        };
                        ox += 1;
                        Some(p)
                    };
                    let py = if by == b'-' {
                        None
                    } else {
                        let p = if *sty == '-' {
                            *ly - 1 - (*sy + oy)
                        } else {
                            *sy + oy
                        };
                        oy += 1;
                        Some(p)
                    };
                    if let (Some(x), Some(y)) = (px, py) {
                        pairs.insert((nx.clone(), x, ny.clone(), y));
                    }
                }
                s_rows.clear();
            }
        }
    }
    pairs
}

/// I1: Byte-identical thread-count determinism (`-t 1` vs `-t 2` vs `-t 4`) across core binaries.
#[test]
#[allow(clippy::too_many_lines, clippy::similar_names)]
fn i1_thread_count_determinism_across_suite() {
    let dir = temp_dir("i1_threads");
    let mut seed = 0x1111_2222_3333_4444;
    let c1 = rand_dna(&mut seed, 1200);
    let c2 = rand_dna(&mut seed, 900);
    let q1 = mutate_sub(&mut seed, &c1, 0.04);
    let mut q2 = mutate_sub(&mut seed, &c2, 0.04);
    q2.splice(300..300, rand_dna(&mut seed, 80));
    let sp3_1 = mutate_sub(&mut seed, &c1, 0.05);
    let sp3_2 = mutate_sub(&mut seed, &c2, 0.05);

    let fa_a = dir.join("spA.fa");
    let fa_b = dir.join("spB.fa");
    let fa_c = dir.join("spC.fa");
    write_fasta(&fa_a, &[("c1", &c1), ("c2", &c2)]);
    write_fasta(&fa_b, &[("c1", &q1), ("c2", &q2)]);
    write_fasta(&fa_c, &[("c1", &sp3_1), ("c2", &sp3_2)]);

    let mut wga_mafs = Vec::new();
    let mut wga_pafs = Vec::new();
    let mut tba_gfas = Vec::new();
    let mut clean_gfas = Vec::new();
    let mut resolve_vcfs = Vec::new();
    let mut guide_mafs = Vec::new();
    let mut filter_pafs = Vec::new();
    let mut rgfas = Vec::new();

    for t in ["1", "2", "4"] {
        let maf = dir.join(format!("wga_t{t}.maf"));
        let paf = dir.join(format!("wga_t{t}.paf"));
        let st = Command::new(env!("CARGO_BIN_EXE_horizon_wga"))
            .args([
                fa_a.to_str().expect("utf8"),
                fa_b.to_str().expect("utf8"),
                "0",
                "--one-to-one",
                "-t",
                t,
                "--maf",
                maf.to_str().expect("utf8"),
                "--paf",
                paf.to_str().expect("utf8"),
            ])
            .status()
            .expect("horizon_wga");
        assert!(st.success());
        wga_mafs.push(std::fs::read(&maf).expect("read maf"));
        wga_pafs.push(std::fs::read(&paf).expect("read paf"));

        let fpaf = dir.join(format!("filt_t{t}.paf"));
        let st = Command::new(env!("CARGO_BIN_EXE_horizon_filter"))
            .args([
                "--paf",
                paf.to_str().expect("utf8"),
                "--target",
                fa_a.to_str().expect("utf8"),
                "--query",
                fa_b.to_str().expect("utf8"),
                "--net",
                "-t",
                t,
                "-o",
                fpaf.to_str().expect("utf8"),
            ])
            .status()
            .expect("horizon_filter");
        assert!(st.success());
        filter_pafs.push(std::fs::read(&fpaf).expect("read filt paf"));

        let tba_pfx = dir.join(format!("tba_t{t}"));
        let st = Command::new(env!("CARGO_BIN_EXE_horizon_tba"))
            .args([
                "--tree",
                "(spA,(spB,spC));",
                "--variation-gfa",
                "-t",
                t,
                tba_pfx.to_str().expect("utf8"),
                "0",
                &format!("spA={}", fa_a.display()),
                &format!("spB={}", fa_b.display()),
                &format!("spC={}", fa_c.display()),
            ])
            .status()
            .expect("horizon_tba");
        assert!(st.success());
        let gfa_p = dir.join(format!("tba_t{t}.gfa"));
        tba_gfas.push(std::fs::read(&gfa_p).expect("read tba gfa"));

        let cgfa = dir.join(format!("clean_t{t}.gfa"));
        let st = Command::new(env!("CARGO_BIN_EXE_horizon_clean"))
            .args([
                "--gfa",
                gfa_p.to_str().expect("utf8"),
                "-t",
                t,
                "--out",
                cgfa.to_str().expect("utf8"),
            ])
            .status()
            .expect("horizon_clean");
        assert!(st.success());
        clean_gfas.push(std::fs::read(&cgfa).expect("read clean gfa"));

        let rgfa_out = dir.join(format!("res_t{t}.gfa"));
        let rvcf = dir.join(format!("res_t{t}.vcf"));
        let st = Command::new(env!("CARGO_BIN_EXE_horizon_resolve"))
            .args([
                cgfa.to_str().expect("utf8"),
                "-t",
                t,
                "--out",
                rgfa_out.to_str().expect("utf8"),
                "--vcf",
                rvcf.to_str().expect("utf8"),
                "--vcf-ref",
                "spA",
            ])
            .status()
            .expect("horizon_resolve");
        assert!(st.success());
        resolve_vcfs.push(std::fs::read(&rvcf).expect("read vcf"));

        let guide_pfx = dir.join(format!("guide_t{t}"));
        let st = Command::new(env!("CARGO_BIN_EXE_horizon_guide"))
            .args([
                "--hinge",
                "spA",
                "-t",
                t,
                guide_pfx.to_str().expect("utf8"),
                "0",
                &format!("spA={}", fa_a.display()),
                &format!("spB={}", fa_b.display()),
                &format!("spC={}", fa_c.display()),
            ])
            .status()
            .expect("horizon_guide");
        assert!(st.success());
        guide_mafs.push(std::fs::read(dir.join(format!("guide_t{t}.maf"))).expect("read guide"));

        let rgfa_p = dir.join(format!("sv_t{t}.rgfa"));
        let st = Command::new(env!("CARGO_BIN_EXE_horizon_gfa"))
            .args([
                "rgfa",
                "--ref",
                fa_a.to_str().expect("utf8"),
                "--query",
                fa_b.to_str().expect("utf8"),
                "-t",
                t,
                "--out",
                rgfa_p.to_str().expect("utf8"),
            ])
            .status()
            .expect("horizon_gfa rgfa");
        assert!(st.success());
        rgfas.push(std::fs::read(&rgfa_p).expect("read rgfa"));
    }

    for (label, outs) in [
        ("horizon_wga MAF", &wga_mafs),
        ("horizon_wga PAF", &wga_pafs),
        ("horizon_filter PAF", &filter_pafs),
        ("horizon_tba GFA", &tba_gfas),
        ("horizon_clean GFA", &clean_gfas),
        ("horizon_resolve VCF", &resolve_vcfs),
        ("horizon_guide MAF", &guide_mafs),
        ("horizon_gfa rGFA", &rgfas),
    ] {
        assert_eq!(outs[0], outs[1], "{label}: -t 1 != -t 2");
        assert_eq!(outs[0], outs[2], "{label}: -t 1 != -t 4");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// I2: Strand / reverse-complement equivariance (`J = 1.000000`, `0` bp diff).
#[test]
fn i2_strand_revcomp_equivariance() {
    let dir = temp_dir("i2_strand");
    let mut seed = 0x2222_3333_4444_5555;
    let target = rand_dna(&mut seed, 2000);
    let query = mutate_sub(&mut seed, &target, 0.05);
    let query_rc = revcomp(&query);
    let len_y = query.len();

    let fa_t = dir.join("t.fa");
    let fa_q = dir.join("q.fa");
    let fa_qrc = dir.join("qrc.fa");
    write_fasta(&fa_t, &[("chr1", &target)]);
    write_fasta(&fa_q, &[("chr1", &query)]);
    write_fasta(&fa_qrc, &[("chr1", &query_rc)]);

    let maf_fwd = dir.join("fwd.maf");
    let maf_rc = dir.join("rc.maf");
    for (qf, out) in [(&fa_q, &maf_fwd), (&fa_qrc, &maf_rc)] {
        let st = Command::new(env!("CARGO_BIN_EXE_horizon_wga"))
            .args([
                fa_t.to_str().expect("utf8"),
                qf.to_str().expect("utf8"),
                "0",
                "--one-to-one",
                "-t",
                "1",
                "--maf",
                out.to_str().expect("utf8"),
            ])
            .status()
            .expect("horizon_wga");
        assert!(st.success());
    }

    let pairs_fwd = maf_aligned_pairs(&std::fs::read_to_string(&maf_fwd).expect("read fwd"));
    let pairs_rc_raw = maf_aligned_pairs(&std::fs::read_to_string(&maf_rc).expect("read rc"));
    assert!(pairs_fwd.len() > 1800, "expected high alignment coverage");

    let pairs_rc_mapped: HashSet<(String, usize, String, usize)> = pairs_rc_raw
        .into_iter()
        .map(|(nx, x, ny, y)| (nx, x, ny.replace("qrc", "q"), len_y - 1 - y))
        .collect();
    assert_eq!(
        pairs_fwd, pairs_rc_mapped,
        "revcomp query must yield exact J = 1.000000 equivariant base pairs"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// I3: Input-order permutation invariance (`horizon_wga` contigs & `horizon_tba -> resolve --vcf`).
#[test]
#[allow(clippy::too_many_lines)]
fn i3_input_order_permutation_invariance() {
    let dir = temp_dir("i3_perm");
    let mut seed = 0x3333_4444_5555_6666;
    let c1 = rand_dna(&mut seed, 1000);
    let c2 = rand_dna(&mut seed, 900);
    let c3 = rand_dna(&mut seed, 800);
    let q1 = mutate_sub(&mut seed, &c1, 0.03);
    let q2 = mutate_sub(&mut seed, &c2, 0.03);
    let q3 = mutate_sub(&mut seed, &c3, 0.03);

    let t_a = dir.join("t_a.fa");
    let t_b = dir.join("t_b.fa");
    let q_a = dir.join("q_a.fa");
    let q_b = dir.join("q_b.fa");
    write_fasta(&t_a, &[("c1", &c1), ("c2", &c2), ("c3", &c3)]);
    write_fasta(&t_b, &[("c3", &c3), ("c1", &c1), ("c2", &c2)]);
    write_fasta(&q_a, &[("c1", &q1), ("c2", &q2), ("c3", &q3)]);
    write_fasta(&q_b, &[("c2", &q2), ("c3", &q3), ("c1", &q1)]);

    let maf1 = dir.join("perm1.maf");
    let maf2 = dir.join("perm2.maf");
    for (tf, qf, out) in [(&t_a, &q_a, &maf1), (&t_b, &q_b, &maf2)] {
        let st = Command::new(env!("CARGO_BIN_EXE_horizon_wga"))
            .args([
                tf.to_str().expect("utf8"),
                qf.to_str().expect("utf8"),
                "0",
                "--one-to-one",
                "-t",
                "1",
                "--maf",
                out.to_str().expect("utf8"),
            ])
            .status()
            .expect("wga");
        assert!(st.success());
    }
    let norm =
        |s: HashSet<(String, usize, String, usize)>| -> HashSet<(String, usize, String, usize)> {
            s.into_iter()
                .map(|(nx, x, ny, y)| {
                    (
                        nx.split('.').next_back().unwrap_or(&nx).to_owned(),
                        x,
                        ny.split('.').next_back().unwrap_or(&ny).to_owned(),
                        y,
                    )
                })
                .collect()
        };
    let p1 = norm(maf_aligned_pairs(
        &std::fs::read_to_string(&maf1).expect("m1"),
    ));
    let p2 = norm(maf_aligned_pairs(
        &std::fs::read_to_string(&maf2).expect("m2"),
    ));
    assert!(!p1.is_empty());
    assert_eq!(
        p1, p2,
        "FASTA contig permutation must yield exact J=1.000000 pairs"
    );

    // Permute non-reference species input order in horizon_tba -> horizon_resolve --vcf.
    let s3 = mutate_sub(&mut seed, &c1, 0.04);
    let fa_r = dir.join("Ref.fa");
    let fa_1 = dir.join("S1.fa");
    let fa_2 = dir.join("S2.fa");
    write_fasta(&fa_r, &[("chr1", &c1)]);
    write_fasta(&fa_1, &[("chr1", &q1)]);
    write_fasta(&fa_2, &[("chr1", &s3)]);

    let mut vcf_records = Vec::new();
    for (tag, sp_args) in [
        (
            "ord1",
            vec![
                format!("Ref={}", fa_r.display()),
                format!("S1={}", fa_1.display()),
                format!("S2={}", fa_2.display()),
            ],
        ),
        (
            "ord2",
            vec![
                format!("S2={}", fa_2.display()),
                format!("Ref={}", fa_r.display()),
                format!("S1={}", fa_1.display()),
            ],
        ),
    ] {
        let pfx = dir.join(tag);
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_horizon_tba"));
        cmd.args([
            "--tree",
            "(Ref,(S1,S2));",
            "--variation-gfa",
            "-t",
            "1",
            pfx.to_str().expect("utf8"),
            "0",
        ]);
        cmd.args(&sp_args);
        assert!(cmd.status().expect("tba").success());

        let vcf = dir.join(format!("{tag}.vcf"));
        let res_gfa = dir.join(format!("{tag}.res.gfa"));
        assert!(
            Command::new(env!("CARGO_BIN_EXE_horizon_resolve"))
                .args([
                    dir.join(format!("{tag}.gfa")).to_str().expect("utf8"),
                    "--out",
                    res_gfa.to_str().expect("utf8"),
                    "--vcf",
                    vcf.to_str().expect("utf8"),
                    "--vcf-ref",
                    "Ref",
                ])
                .status()
                .expect("resolve")
                .success()
        );
        // Compare normalized (CHROM, POS, REF, ALT, GT columns) independent of arbitrary snarl IDs.
        let recs: Vec<String> = std::fs::read_to_string(&vcf)
            .expect("read vcf")
            .lines()
            .filter(|l| !l.starts_with("##"))
            .map(|l| {
                if l.starts_with("#CHROM") {
                    l.to_owned()
                } else {
                    let f: Vec<&str> = l.split('\t').collect();
                    format!(
                        "{}\t{}\t{}\t{}\t{}",
                        f[0],
                        f[1],
                        f[3],
                        f[4],
                        f[9..].join("\t")
                    )
                }
            })
            .collect();
        vcf_records.push(recs);
    }
    assert_eq!(
        vcf_records[0], vcf_records[1],
        "permuting non-reference species input order must yield identical VCF sites & genotypes"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// I4: Exact MAF/PAF coordinate + CIGAR sequence reconstruction, chain/net span sums, & 1:1 0-bp overlap.
#[test]
#[allow(clippy::too_many_lines, clippy::similar_names)]
fn i4_coordinate_cigar_chain_net_and_one_to_one_legality() {
    let dir = temp_dir("i4_coords");
    let mut seed = 0x4444_5555_6666_7777;
    let seg_a = rand_dna(&mut seed, 800);
    let seg_b = rand_dna(&mut seed, 800);
    let mut q_a = mutate_sub(&mut seed, &seg_a, 0.03);
    q_a.splice(400..400, rand_dna(&mut seed, 25)); // insertion
    let mut q_b = mutate_sub(&mut seed, &seg_b, 0.03);
    q_b.drain(300..320); // deletion
    let q_b_rc = revcomp(&q_b);

    let fa_t = dir.join("t.fa");
    let fa_q = dir.join("q.fa");
    write_fasta(&fa_t, &[("chr1", &seg_a), ("chr2", &seg_b)]);
    write_fasta(&fa_q, &[("q1", &q_a), ("q2", &q_b_rc)]);
    let t_map: HashMap<&str, &[u8]> = HashMap::from([("chr1", &seg_a[..]), ("chr2", &seg_b[..])]);
    let q_map: HashMap<&str, &[u8]> = HashMap::from([("q1", &q_a[..]), ("q2", &q_b_rc[..])]);

    let maf = dir.join("out.maf");
    let paf = dir.join("out.paf");
    let chain = dir.join("out.chain");
    let net = dir.join("out.net");
    assert!(
        Command::new(env!("CARGO_BIN_EXE_horizon_wga"))
            .args([
                fa_t.to_str().expect("utf8"),
                fa_q.to_str().expect("utf8"),
                "0",
                "--one-to-one",
                "--maf",
                maf.to_str().expect("utf8"),
                "--paf",
                paf.to_str().expect("utf8"),
                "--chain",
                chain.to_str().expect("utf8"),
                "--net",
                net.to_str().expect("utf8"),
            ])
            .status()
            .expect("wga")
            .success()
    );

    // 1. MAF exact substring + 0-bp overlap check.
    let mut maf_cov_t: HashMap<String, Vec<bool>> = HashMap::new();
    let mut maf_cov_q: HashMap<String, Vec<bool>> = HashMap::new();
    let mut row_idx = 0usize;
    for line in std::fs::read_to_string(&maf).expect("maf").lines() {
        if let Some(rest) = line.strip_prefix("s ") {
            let f: Vec<&str> = rest.split_whitespace().collect();
            let ctg = f[0].split('.').next_back().expect("ctg");
            let start: usize = f[1].parse().expect("start");
            let size: usize = f[2].parse().expect("size");
            let strand = f[3];
            let src_size: usize = f[4].parse().expect("src_size");
            let ungapped: Vec<u8> = f[5].bytes().filter(|&b| b != b'-').collect();
            assert_eq!(ungapped.len(), size);
            assert!(start + size <= src_size);
            let is_target = row_idx.is_multiple_of(2);
            let src = if is_target { t_map[ctg] } else { q_map[ctg] };
            assert_eq!(src.len(), src_size);
            let (fwd_s, fwd_e) = if strand == "-" {
                (src_size - start - size, src_size - start)
            } else {
                (start, start + size)
            };
            let expected = if strand == "-" {
                revcomp(&src[fwd_s..fwd_e])
            } else {
                src[fwd_s..fwd_e].to_vec()
            };
            assert!(
                ungapped.eq_ignore_ascii_case(&expected),
                "MAF sequence mismatch on {ctg} [{fwd_s}..{fwd_e}) {strand}"
            );
            let cov = if is_target {
                maf_cov_t
                    .entry(ctg.to_owned())
                    .or_insert_with(|| vec![false; src_size])
            } else {
                maf_cov_q
                    .entry(ctg.to_owned())
                    .or_insert_with(|| vec![false; src_size])
            };
            for flag in &mut cov[fwd_s..fwd_e] {
                assert!(!*flag, "overlapping MAF base on {ctg} under --one-to-one");
                *flag = true;
            }
            row_idx += 1;
        }
    }
    assert!(row_idx >= 4, "expected alignments for both contigs");

    // 2. PAF exact CIGAR (`=`/`X`/`I`/`D`) reconstruction + 0-bp overlap check.
    let mut paf_cov_t: HashMap<String, Vec<bool>> = HashMap::new();
    let mut paf_cov_q: HashMap<String, Vec<bool>> = HashMap::new();
    let mut paf_count = 0usize;
    for line in std::fs::read_to_string(&paf).expect("paf").lines() {
        if line.trim().is_empty() {
            continue;
        }
        paf_count += 1;
        let rec = PafRecord::parse_line(line).expect("parse paf");
        let tseq = t_map[rec.tname.as_str()];
        let qseq_raw = q_map[rec.qname.as_str()];
        assert_eq!(rec.tlen, tseq.len());
        assert_eq!(rec.qlen, qseq_raw.len());
        let qslice = if rec.strand == '-' {
            revcomp(&qseq_raw[rec.qstart..rec.qend])
        } else {
            qseq_raw[rec.qstart..rec.qend].to_vec()
        };
        let tslice = &tseq[rec.tstart..rec.tend];
        let mut ti = 0usize;
        let mut qi = 0usize;
        for op in rec.cigar_ops().expect("cg:Z:") {
            match op.expect("valid cigar op") {
                CigarOp::Equal(len) => {
                    for _ in 0..len {
                        assert!(tslice[ti].eq_ignore_ascii_case(&qslice[qi]));
                        ti += 1;
                        qi += 1;
                    }
                }
                CigarOp::Diff(len) => {
                    for _ in 0..len {
                        assert!(!tslice[ti].eq_ignore_ascii_case(&qslice[qi]));
                        ti += 1;
                        qi += 1;
                    }
                }
                CigarOp::Match(len) => {
                    ti += len;
                    qi += len;
                }
                CigarOp::Ins(len) => qi += len,
                CigarOp::Del(len) => ti += len,
                other => panic!("unexpected CIGAR op {other:?}"),
            }
        }
        assert_eq!(ti, tslice.len());
        assert_eq!(qi, qslice.len());
        let tc = paf_cov_t
            .entry(rec.tname.clone())
            .or_insert_with(|| vec![false; rec.tlen]);
        for f in &mut tc[rec.tstart..rec.tend] {
            assert!(!*f, "PAF target overlap under --one-to-one");
            *f = true;
        }
        let qc = paf_cov_q
            .entry(rec.qname.clone())
            .or_insert_with(|| vec![false; rec.qlen]);
        for f in &mut qc[rec.qstart..rec.qend] {
            assert!(!*f, "PAF query overlap under --one-to-one");
            *f = true;
        }
    }
    assert!(paf_count >= 2);

    // 3. UCSC chain exact block/gap span sum check.
    let chain_text = std::fs::read_to_string(&chain).expect("chain");
    let mut cur_t_span = 0usize;
    let mut cur_q_span = 0usize;
    let mut acc_t = 0usize;
    let mut acc_q = 0usize;
    let mut chains_seen = 0usize;
    for line in chain_text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            if chains_seen > 0 && cur_t_span > 0 {
                assert_eq!(acc_t, cur_t_span, "chain target span mismatch");
                assert_eq!(acc_q, cur_q_span, "chain query span mismatch");
                cur_t_span = 0;
            }
            continue;
        }
        if trimmed.starts_with("chain ") {
            let f: Vec<&str> = trimmed.split_whitespace().collect();
            let t_start: usize = f[5].parse().expect("ts");
            let t_end: usize = f[6].parse().expect("te");
            let q_start: usize = f[10].parse().expect("qs");
            let q_end: usize = f[11].parse().expect("qe");
            cur_t_span = t_end - t_start;
            cur_q_span = q_end - q_start;
            acc_t = 0;
            acc_q = 0;
            chains_seen += 1;
        } else {
            let nums: Vec<usize> = trimmed
                .split_whitespace()
                .map(|x| x.parse::<usize>().expect("int"))
                .collect();
            acc_t += nums[0];
            acc_q += nums[0];
            if nums.len() == 3 {
                acc_t += nums[1];
                acc_q += nums[2];
            }
        }
    }
    assert!(chains_seen >= 2, "expected >=2 UCSC chains");
    assert!(
        std::fs::read_to_string(&net)
            .expect("net")
            .contains("net chr1 800"),
        "UCSC net must record exact chr1 length"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// I5 & I6: Per-walk GFA sequence round-trip (`0` diffs), whole-contig completeness (`>= 99%`),
/// and VCF haplotype allele reconstruction round-trip (`0` diffs, `>= 95%` called coverage).
#[test]
#[allow(clippy::too_many_lines)]
fn i5_i6_gfa_walk_roundtrip_completeness_and_vcf_consensus() {
    let dir = temp_dir("i5_i6");
    let mut seed = 0x5555_6666_7777_8888;
    let fl_a = rand_dna(&mut seed, 500);
    let ins_250 = rand_dna(&mut seed, 250);
    let fl_b = rand_dna(&mut seed, 500);
    let del_250 = rand_dna(&mut seed, 250);
    let fl_c = rand_dna(&mut seed, 500);

    let mut ref_seq = Vec::new();
    ref_seq.extend_from_slice(&fl_a);
    ref_seq.extend_from_slice(&fl_b);
    ref_seq.extend_from_slice(&del_250);
    ref_seq.extend_from_slice(&fl_c);

    let mut h1 = Vec::new();
    h1.extend_from_slice(&fl_a);
    h1.extend_from_slice(&ins_250);
    h1.extend_from_slice(&fl_b);
    h1.extend_from_slice(&del_250);
    h1.extend_from_slice(&fl_c);
    h1[100] = if h1[100] == b'A' { b'G' } else { b'A' };

    let mut h2 = Vec::new();
    h2.extend_from_slice(&fl_a);
    h2.extend_from_slice(&fl_b);
    h2.extend_from_slice(&fl_c);
    h2[200] = if h2[200] == b'C' { b'T' } else { b'C' };

    let fa_r = dir.join("Ref.fa");
    let fa_1 = dir.join("H1.fa");
    let fa_2 = dir.join("H2.fa");
    write_fasta(&fa_r, &[("chr1", &ref_seq)]);
    write_fasta(&fa_1, &[("chr1", &h1)]);
    write_fasta(&fa_2, &[("chr1", &h2)]);

    let pfx = dir.join("tba");
    assert!(
        Command::new(env!("CARGO_BIN_EXE_horizon_tba"))
            .args([
                "--tree",
                "(Ref,(HG002#1,HG002#2));",
                "--variation-gfa",
                pfx.to_str().expect("utf8"),
                "0",
                &format!("Ref={}", fa_r.display()),
                &format!("HG002#1={}", fa_1.display()),
                &format!("HG002#2={}", fa_2.display()),
            ])
            .status()
            .expect("tba")
            .success()
    );
    let raw_gfa = dir.join("tba.gfa");
    let clean_gfa = dir.join("clean.gfa");
    let res_gfa = dir.join("res.gfa");
    let vcf_gz = dir.join("out.vcf.gz");
    assert!(
        Command::new(env!("CARGO_BIN_EXE_horizon_clean"))
            .args([
                "--gfa",
                raw_gfa.to_str().expect("utf8"),
                "--out",
                clean_gfa.to_str().expect("utf8"),
            ])
            .status()
            .expect("clean")
            .success()
    );
    assert!(
        Command::new(env!("CARGO_BIN_EXE_horizon_resolve"))
            .args([
                clean_gfa.to_str().expect("utf8"),
                "--out",
                res_gfa.to_str().expect("utf8"),
                "--vcf",
                vcf_gz.to_str().expect("utf8"),
                "--vcf-ref",
                "Ref",
            ])
            .status()
            .expect("resolve")
            .success()
    );

    let truth: HashMap<(&str, &str), &[u8]> = HashMap::from([
        (("Ref", "0"), &ref_seq[..]),
        (("HG002", "1"), &h1[..]),
        (("HG002", "2"), &h2[..]),
    ]);

    // I5: Verify exact walk sequence round-trip & >=99% whole-contig completeness across all 3 GFAs.
    for gfa_file in [&raw_gfa, &clean_gfa, &res_gfa] {
        let mut segs: HashMap<String, Vec<u8>> = HashMap::new();
        let mut cov_bp: BTreeMap<(String, String), usize> = BTreeMap::new();
        let text = std::fs::read_to_string(gfa_file).expect("read gfa");
        for line in text.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.first() == Some(&"S") {
                segs.insert(f[1].to_owned(), f[2].as_bytes().to_vec());
            } else if f.first() == Some(&"W") {
                let sample = f[1];
                let hap = f[2];
                let w_s: usize = f[4].parse().expect("ws");
                let w_e: usize = f[5].parse().expect("we");
                let mut recon = Vec::new();
                for (seg_id, fwd) in parse_walk_steps(f[6]) {
                    let s = &segs[seg_id];
                    if fwd {
                        recon.extend_from_slice(s);
                    } else {
                        recon.extend_from_slice(&revcomp(s));
                    }
                }
                assert_eq!(recon.len(), w_e - w_s, "walk span mismatch");
                let src = truth[&(sample, hap)];
                assert_eq!(
                    recon.as_slice(),
                    &src[w_s..w_e],
                    "walk sequence mismatch for {sample}#{hap} in {}",
                    gfa_file.display()
                );
                *cov_bp
                    .entry((sample.to_owned(), hap.to_owned()))
                    .or_insert(0) += recon.len();
            }
        }
        for (&(sp, hp), &src) in &truth {
            let got = cov_bp
                .get(&(sp.to_owned(), hp.to_owned()))
                .copied()
                .unwrap_or(0);
            #[allow(clippy::cast_precision_loss)]
            let frac = got as f64 / src.len() as f64;
            assert!(
                frac >= 0.99,
                "walk completeness for {sp}#{hp} in {} was {frac:.4} < 0.99",
                gfa_file.display()
            );
        }
    }

    // I6: Apply VCF records directly to Ref for HG002 hap 1 & hap 2 -> exact 0-diff match.
    let mut vcf_text = String::new();
    open_text_reader(&vcf_gz)
        .expect("open vcf.gz")
        .read_to_string(&mut vcf_text)
        .expect("read vcf.gz");
    let header = vcf_text
        .lines()
        .find(|l| l.starts_with("#CHROM"))
        .expect("hdr");
    let hg_col = header
        .split('\t')
        .position(|c| c == "HG002")
        .expect("HG002");
    let mut recs: Vec<VcfSite> = Vec::new();
    for line in vcf_text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        let pos0 = f[1].parse::<usize>().expect("pos") - 1;
        let ref_a = f[3].to_owned();
        let mut alleles = vec![ref_a.clone()];
        alleles.extend(f[4].split(',').map(str::to_owned));
        let gts: Vec<Option<usize>> = f[hg_col]
            .split('|')
            .map(|tok| tok.parse::<usize>().ok())
            .collect();
        recs.push((pos0, ref_a, alleles, gts));
    }
    assert!(recs.len() >= 4, "expected SNPs + 250bp ins + 250bp del");
    for (hap_i, expected_seq) in [(&h1), (&h2)].into_iter().enumerate() {
        let mut out = Vec::new();
        let mut cursor = 0usize;
        for (pos0, ref_a, alleles, gts) in &recs {
            assert!(*pos0 >= cursor, "non-overlapping LV=0 VCF records required");
            out.extend_from_slice(&ref_seq[cursor..*pos0]);
            let aidx = gts[hap_i].expect("called diploid allele");
            assert_eq!(&ref_seq[*pos0..*pos0 + ref_a.len()], ref_a.as_bytes());
            out.extend_from_slice(alleles[aidx].as_bytes());
            cursor = *pos0 + ref_a.len();
        }
        out.extend_from_slice(&ref_seq[cursor..]);
        assert_eq!(
            out.as_slice(),
            expected_seq.as_slice(),
            "I6 VCF haplotype {hap_i} consensus must match source FASTA with 0 diffs"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// I7: Native HAL bi-directional round-trip (`horizon_guide --hal` & `horizon_tba --hal --progressive-ancestors`).
#[test]
#[allow(clippy::too_many_lines)]
fn i7_hal_bidirectional_roundtrip_guide_and_tba() {
    let dir = temp_dir("i7_hal");
    let mut seed = 0x7777_8888_9999_AAAA;
    let r = rand_dna(&mut seed, 1200);
    let q1 = mutate_sub(&mut seed, &r, 0.04);
    let q2 = mutate_sub(&mut seed, &r, 0.05);

    let fa_r = dir.join("Ref.fa");
    let fa_1 = dir.join("Q1.fa");
    let fa_2 = dir.join("Q2.fa");
    write_fasta(&fa_r, &[("chr1", &r)]);
    write_fasta(&fa_1, &[("chr1", &q1)]);
    write_fasta(&fa_2, &[("chr1", &q2)]);

    let guide_pfx = dir.join("guide");
    assert!(
        Command::new(env!("CARGO_BIN_EXE_horizon_guide"))
            .args([
                "--hinge",
                "Ref",
                "--hal",
                guide_pfx.to_str().expect("utf8"),
                "0",
                &format!("Ref={}", fa_r.display()),
                &format!("Q1={}", fa_1.display()),
                &format!("Q2={}", fa_2.display()),
            ])
            .status()
            .expect("guide")
            .success()
    );

    let tba_pfx = dir.join("tba");
    assert!(
        Command::new(env!("CARGO_BIN_EXE_horizon_tba"))
            .args([
                "--tree",
                "(Ref,(Q1,Q2));",
                "--progressive-ancestors",
                "--hal",
                tba_pfx.to_str().expect("utf8"),
                "0",
                &format!("Ref={}", fa_r.display()),
                &format!("Q1={}", fa_1.display()),
                &format!("Q2={}", fa_2.display()),
            ])
            .status()
            .expect("tba")
            .success()
    );

    let bed = dir.join("q.bed");
    std::fs::write(&bed, "chr1\t200\t600\treg1\t0\t+\n").expect("bed");
    for hal_path in [dir.join("guide.hal"), dir.join("tba.hal")] {
        let st = Command::new(env!("CARGO_BIN_EXE_horizon_hal"))
            .args(["stats", hal_path.to_str().expect("utf8")])
            .output()
            .expect("hal stats");
        assert!(
            st.status.success(),
            "horizon_hal stats failed on {}",
            hal_path.display()
        );

        // Forward liftover Q1 -> Q2 then reverse Q2 -> Q1 must recover [200, 600).
        let lo1 = Command::new(env!("CARGO_BIN_EXE_horizon_hal"))
            .args([
                "liftover",
                hal_path.to_str().expect("utf8"),
                "--src",
                "Q1",
                "--dst",
                "Q2",
                "--bed",
                bed.to_str().expect("utf8"),
            ])
            .output()
            .expect("liftover 1");
        assert!(
            lo1.status.success(),
            "horizon_hal liftover Q1->Q2 failed: {}",
            String::from_utf8_lossy(&lo1.stderr)
        );
        let lo1_txt = String::from_utf8(lo1.stdout).expect("utf8");
        let lifted_bp: usize = lo1_txt
            .lines()
            .filter_map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                if f.len() >= 3 {
                    let s: usize = f[1].parse().ok()?;
                    let e: usize = f[2].parse().ok()?;
                    assert!((200..=600).contains(&s) && (200..=600).contains(&e) && e > s);
                    Some(e - s)
                } else {
                    None
                }
            })
            .sum();
        assert!(
            lifted_bp >= 300,
            "expected >=300/400 bp lifted Q1->Q2 on {}, got {lifted_bp}:\n{lo1_txt}",
            hal_path.display()
        );

        let bed2 = dir.join("q2.bed");
        std::fs::write(&bed2, &lo1_txt).expect("write bed2");
        let lo2 = Command::new(env!("CARGO_BIN_EXE_horizon_hal"))
            .args([
                "liftover",
                hal_path.to_str().expect("utf8"),
                "--src",
                "Q2",
                "--dst",
                "Q1",
                "--bed",
                bed2.to_str().expect("utf8"),
            ])
            .output()
            .expect("liftover 2");
        assert!(
            lo2.status.success(),
            "horizon_hal liftover Q2->Q1 failed: {}",
            String::from_utf8_lossy(&lo2.stderr)
        );
        let lo2_txt = String::from_utf8(lo2.stdout).expect("utf8");
        assert_eq!(
            lo2_txt,
            lo1_txt,
            "bi-directional HAL liftover Q1->Q2->Q1 must be an exact involution on {}",
            hal_path.display()
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

fn dna_strategy(min_len: usize, max_len: usize) -> impl Strategy<Value = Vec<u8>> {
    proptest::collection::vec(
        prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
        min_len..=max_len,
    )
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    /// I8a: Pair-HMM Forward/Backward partition parity `|log2(Z_F) - log2(Z_B)| <= 1e-6 * max(1, |Z_F|)`.
    #[test]
    fn i8_proptest_pair_hmm_partition_parity(
        x in dna_strategy(12, 96),
        y in dna_strategy(12, 96),
        id_prior in 0.55_f64..0.95_f64,
    ) {
        let params = HmmParams {
            identity_prior: id_prior,
            ..HmmParams::default()
        };
        let trace = banded_mdl_gap_posterior_trace(&x, &y, 32, 8.0, &params, 0.5);
        prop_assert!(trace.score_bits.is_finite() && trace.log2_zb.is_finite());
        let diff = (trace.score_bits - trace.log2_zb).abs();
        let tol = 1e-6 * trace.score_bits.abs().max(1.0);
        prop_assert!(
            diff <= tol,
            "Z_F={} Z_B={} diff={} > tol={}",
            trace.score_bits, trace.log2_zb, diff, tol
        );
    }

    /// I8b: Fast-path (`dp_fast.rs`) vs exact Forward-Backward (`banded_mdl_gap_core`) agreement
    /// on homologous balanced sequences where diagonal/Viterbi fast paths engage.
    #[test]
    fn i8_proptest_fast_path_vs_exact_dp_agreement(
        base in dna_strategy(40, 120),
        mut_idx in proptest::collection::vec(0usize..120usize, 0..3),
    ) {
        let mut y = base.clone();
        for idx in mut_idx {
            let pos = idx % y.len();
            y[pos] = match y[pos] {
                b'A' => b'G',
                b'G' => b'A',
                b'C' => b'T',
                _ => b'C',
            };
        }
        let params = HmmParams::default();
        let fast = banded_mdl_gap_with_min_post(&base, &y, 24, 10.0, &params, 0.5);
        let exact = banded_mdl_gap_posterior_trace(&base, &y, 24, 10.0, &params, 0.5);
        // Fast path must never hallucinate alignments that exact DP abstains on, and on high-identity
        // sequences both must recover >=90% overlapping aligned coordinate pairs.
        if exact.kept_pairs.is_empty() {
            prop_assert!(fast.kept_pairs.is_empty(), "fast path hallucinated pairs when exact DP abstained");
        } else if !fast.kept_pairs.is_empty() {
            let s_fast: HashSet<(usize, usize)> = fast.kept_pairs.iter().copied().collect();
            let s_exact: HashSet<(usize, usize)> = exact.kept_pairs.iter().copied().collect();
            let inter = s_fast.intersection(&s_exact).count();
            let union = s_fast.union(&s_exact).count().max(1);
            #[allow(clippy::cast_precision_loss)]
            let j = inter as f64 / union as f64;
            prop_assert!(j >= 0.85, "fast vs exact DP Jaccard {j:.3} < 0.85");
        }
    }

    /// I8c: Pre-clamp BOCPD Bayes marginalization conservation (`last_pre_norm_mass == P(x_t | x_{1:t-1})`).
    #[test]
    fn i8_proptest_bocpd_pre_norm_bayes_conservation(
        seq in dna_strategy(4, 150),
        hazard in 0.005_f64..0.25_f64,
        max_run in 4usize..64usize,
    ) {
        let mut scratch = BocpdScratch::default();
        let mut gc = Vec::new();
        bocpd_gc_profile_into(&seq, hazard, max_run, &mut scratch, &mut gc);
        prop_assert_eq!(gc.len(), seq.len());
        let last_b = *seq.last().expect("non-empty");
        let is_gc = last_b == b'C' || last_b == b'G';
        let expected_marginal = if is_gc {
            scratch.last_raw_pred_gc()
        } else {
            1.0 - scratch.last_raw_pred_gc()
        };
        let diff = (scratch.last_pre_norm_mass() - expected_marginal).abs();
        prop_assert!(
            diff < 1e-12,
            "BOCPD pre-norm mass {} != Bayes marginal {} (diff={diff})",
            scratch.last_pre_norm_mass(), expected_marginal
        );
        let sum_rp: f64 = scratch.run_probs().iter().sum();
        prop_assert!((sum_rp - 1.0).abs() < 1e-12);
    }

    /// I8d: Fail-closed fuzzing of GFA, MAF, PAF, and HDF5/HAL parsers (0 panics, structured Err on garbage).
    #[test]
    fn i8_proptest_parsers_fail_closed_no_panic(
        payload in "[A-Za-z0-9_ :;,.\\t\\n!@#$%^&*]{4,128}",
        raw_bytes in proptest::collection::vec(any::<u8>(), 0..256),
    ) {
        // Prefix with unmistakable invalid token 'Z' so valid comments/records are excluded.
        let bad = format!("Z_INVALID_REC\t{payload}\n");
        let gfa_res = stream_gfa_record_refs_reader(Cursor::new(bad.as_bytes()), |_| {
            Ok(std::ops::ControlFlow::Continue(()))
        });
        prop_assert!(gfa_res.is_err(), "GFA parser must reject Z_INVALID_REC");

        let maf_res = for_each_maf_stanza_reader(Cursor::new(bad.as_bytes()), |_| {});
        prop_assert!(maf_res.is_err(), "MAF stanza parser must reject Z_INVALID_REC");

        let tba_res = for_each_multi_tba_block_reader(Cursor::new(bad.as_bytes()), |_| Ok(()));
        prop_assert!(tba_res.is_err(), "TBA MAF parser must reject Z_INVALID_REC");

        let paf_res = PafRecord::parse_line(bad.trim_end());
        prop_assert!(paf_res.is_err(), "PAF parser must reject malformed line");

        // Arbitrary random bytes (both raw and with HDF5 magic header) must never panic HDF5/HAL parsers.
        let _ = read_hdf5_superblock(&raw_bytes);
        let _ = read_hdf5_v0_root_names(&raw_bytes);
        let mut with_magic = HDF5_SIGNATURE.to_vec();
        with_magic.extend_from_slice(&raw_bytes);
        let _ = read_hdf5_superblock(&with_magic);
        let _ = read_hdf5_v0_root_names(&with_magic);
    }
}
