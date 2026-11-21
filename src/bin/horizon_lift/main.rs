//! Graph-aware BED/GFF3 coordinate liftover across GFA v1.0 (`P`) and v1.1 (`W`) walks.

mod features;
mod lift;

use anyhow::{Context, Result, bail};
use features::parse_features;
use horizon_aligner::gfa::{SpoolableGfaInput, open_output_writer};
use lift::{LiftConfig, run_lift};
use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

struct LiftCliArgs {
    gfa_path: String,
    from_sample: String,
    to_sample: String,
    input_path: String,
    is_gff: bool,
    out_path: Option<String>,
    unmapped_path: Option<String>,
    min_frac: f64,
    bridge_bp: usize,
    to_chroms: HashSet<String>,
}

const USAGE: &str = "usage: horizon_lift --gfa <file.gfa[.gz]|-> --from <sample> --to <sample> (--bed <in.bed> | --gff <in.gff3>) [--out <out.bed>] [--unmapped <unmapped.bed>] [--min-frac 0.0] [--bridge-bp N] [--to-chrom C1,C2]";

fn parse_args() -> Result<Option<LiftCliArgs>> {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut gfa_path: Option<String> = None;
    let mut from_sample: Option<String> = None;
    let mut to_sample: Option<String> = None;
    let mut input_path: Option<String> = None;
    let mut is_gff = false;
    let mut out_path: Option<String> = None;
    let mut unmapped_path: Option<String> = None;
    let mut min_frac = 0.0_f64;
    let mut bridge_bp = 0usize;
    let mut to_chroms: HashSet<String> = HashSet::new();

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                println!("{USAGE}");
                return Ok(None);
            }
            "--gfa" => {
                gfa_path = Some(args.get(i + 1).context("missing value for --gfa")?.clone());
                i += 2;
            }
            "--from" => {
                from_sample = Some(args.get(i + 1).context("missing value for --from")?.clone());
                i += 2;
            }
            "--to" => {
                to_sample = Some(args.get(i + 1).context("missing value for --to")?.clone());
                i += 2;
            }
            "--to-chrom" => {
                let raw = args.get(i + 1).context("missing value for --to-chrom")?;
                for c in raw.split(',').map(str::trim).filter(|s| !s.is_empty()) {
                    to_chroms.insert(c.to_string());
                }
                i += 2;
            }
            "--bed" => {
                input_path = Some(args.get(i + 1).context("missing value for --bed")?.clone());
                is_gff = false;
                i += 2;
            }
            "--gff" | "--gff3" => {
                input_path = Some(args.get(i + 1).context("missing value for --gff")?.clone());
                is_gff = true;
                i += 2;
            }
            "--out" | "-o" => {
                out_path = Some(args.get(i + 1).context("missing value for --out")?.clone());
                i += 2;
            }
            "--unmapped" => {
                unmapped_path = Some(
                    args.get(i + 1)
                        .context("missing value for --unmapped")?
                        .clone(),
                );
                i += 2;
            }
            "--min-frac" => {
                min_frac = args
                    .get(i + 1)
                    .context("missing value for --min-frac")?
                    .parse()
                    .context("parsing --min-frac")?;
                anyhow::ensure!(
                    min_frac.is_finite() && (0.0..=1.0).contains(&min_frac),
                    "--min-frac must be finite and in 0.0..=1.0, got {min_frac}"
                );
                i += 2;
            }
            "--bridge-bp" => {
                bridge_bp = args
                    .get(i + 1)
                    .context("missing value for --bridge-bp")?
                    .parse()
                    .context("parsing --bridge-bp")?;
                i += 2;
            }
            _ => bail!("unknown argument '{}'", args[i]),
        }
    }

    let (Some(gfa_path), Some(from_sample), Some(to_sample), Some(input_path)) =
        (gfa_path, from_sample, to_sample, input_path)
    else {
        bail!("{USAGE}");
    };

    Ok(Some(LiftCliArgs {
        gfa_path,
        from_sample,
        to_sample,
        input_path,
        is_gff,
        out_path,
        unmapped_path,
        min_frac,
        bridge_bp,
        to_chroms,
    }))
}

fn main() -> Result<()> {
    let Some(cli) = parse_args()? else {
        return Ok(());
    };

    let t0 = Instant::now();
    let spooled = SpoolableGfaInput::new(&cli.gfa_path)?;
    let features = parse_features(Path::new(&cli.input_path), cli.is_gff)?;
    let cfg = LiftConfig {
        gfa_path: spooled.path(),
        from_sample: &cli.from_sample,
        to_sample: &cli.to_sample,
        to_chroms: &cli.to_chroms,
        features: &features,
        min_frac: cli.min_frac,
        bridge_bp: cli.bridge_bp,
    };

    let mut unmapped_writer = if let Some(ref up) = cli.unmapped_path {
        Some(BufWriter::new(File::create(up)?))
    } else {
        None
    };

    let mut w = open_output_writer(cli.out_path.as_deref())?;
    let summary = run_lift(&cfg, &mut w, unmapped_writer.as_mut())?;
    w.flush()?;
    if let Some(mut uw) = unmapped_writer {
        uw.flush()?;
    }

    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    eprintln!(
        "horizon_lift: features={} exact={} partial={} split={} unmapped={} elapsed={:.2}s peak_rss_mb={:.1}",
        features.len(),
        summary.lifted_exact,
        summary.lifted_partial,
        summary.lifted_split,
        summary.unmapped,
        t0.elapsed().as_secs_f64(),
        hwm_mb
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::features::InputFeature;
    use crate::lift::flip_strand;
    use horizon_aligner::strip_sample_prefix;

    fn test_tmp(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("horizon_lift_{}_{name}", std::process::id()))
    }

    #[test]
    fn lifts_forward_inverted_split_and_unmapped_intervals() {
        // Construct a small GFA with GRCh38 (P line) and CHM13 (W line)
        // Segments: s1=10bp, s2=10bp, s3=10bp (ref-only insertion), s4=10bp
        // GRCh38: s1+, s2+, s3+, s4+ -> total 40bp (0..10, 10..20, 20..30, 30..40)
        // CHM13:  >s1 <s2 >s4       -> total 30bp (0..10=s1+, 10..20=s2-, 20..30=s4+)
        let gfa_path = test_tmp("test.gfa");
        let gfa = "\
H\tVN:Z:1.1
S\ts1\tAAAAAAAAAA
S\ts2\tCCCCCCCCCC
S\ts3\tGGGGGGGGGG
S\ts4\tTTTTTTTTTT
P\tGRCh38.chr22\ts1+,s2+,s3+,s4+\t0M,0M,0M
W\tCHM13\t0\tchr22\t100\t130\t>s1<s2>s4
";
        std::fs::write(&gfa_path, gfa).expect("write test GFA");

        let features = vec![
            // 1. Exact in s1: [2, 8) -> CHM13 [102, 108) '+'
            InputFeature {
                chrom: "chr22".into(),
                start: 2,
                end: 8,
                name: "f_fwd".into(),
                score: "100".into(),
                strand: '+',
            },
            // 2. Inverted in s2: [12, 17) (local 2..7 on s2+) -> CHM13 s2- at [110..120), rev [113, 118) '-'
            InputFeature {
                chrom: "chr22".into(),
                start: 12,
                end: 17,
                name: "f_inv".into(),
                score: "100".into(),
                strand: '+',
            },
            // 3. Split across s1..s4 (spanning s2- and s3 gap): [5, 35)
            InputFeature {
                chrom: "chr22".into(),
                start: 5,
                end: 35,
                name: "f_split".into(),
                score: "100".into(),
                strand: '+',
            },
            // 4. Unmapped inside s3: [22, 28)
            InputFeature {
                chrom: "chr22".into(),
                start: 22,
                end: 28,
                name: "f_unmapped".into(),
                score: "100".into(),
                strand: '+',
            },
        ];

        let cfg = LiftConfig {
            gfa_path: gfa_path.to_str().expect("valid utf8 temp path"),
            from_sample: "GRCh38",
            to_sample: "CHM13",
            to_chroms: &HashSet::new(),
            features: &features,
            min_frac: 0.0,
            bridge_bp: 0,
        };
        let mut out = Vec::new();
        let mut unmapped = Vec::new();
        let summary = run_lift(&cfg, &mut out, Some(&mut unmapped)).expect("run_lift succeeds");
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(summary.lifted_exact, 2);
        assert_eq!(summary.lifted_split, 1);
        assert_eq!(summary.unmapped, 1);

        let out_str = String::from_utf8(out).expect("valid utf8 output");
        let lines: Vec<&str> = out_str.lines().collect();
        assert!(lines[0].starts_with("CHM13.chr22\t102\t108\tf_fwd\t100\t+"));
        assert!(lines[1].starts_with("CHM13.chr22\t113\t118\tf_inv\t100\t-"));
        assert!(
            lines
                .iter()
                .any(|l| l.contains("f_split") && l.contains("SPLIT_1/3"))
        );

        let unm_str = String::from_utf8(unmapped).expect("valid utf8 unmapped");
        assert!(unm_str.contains("f_unmapped"));
    }

    #[test]
    fn disambiguates_multicopy_target_with_flanking_context_jaccard() {
        // Segment `dup` (10bp) appears twice in CHM13:
        //   copy 1 at [10..20) flanked by d1, d2 (distantly located from source context)
        //   copy 2 at [60..70) flanked by s1, s2, s4, s5 (matching source context!)
        // Even with a single-segment feature inside `dup` (so last_target is None),
        // Jaccard context window must pick copy 2 [62..68) over copy 1 [12..18).
        let gfa_path = test_tmp("jaccard_test.gfa");
        let gfa = "\
H\tVN:Z:1.1
S\ts1\tAAAAAAAAAA
S\ts2\tCCCCCCCCCC
S\tdup\tGGGGGGGGGG
S\ts4\tTTTTTTTTTT
S\ts5\tACACACACAC
S\td1\tTGTGTGTGTG
S\td2\tCTAGCTAGCT
P\tGRCh38.chr1\ts1+,s2+,dup+,s4+,s5+\t0M,0M,0M,0M
W\tCHM13\t0\tchr1\t0\t80\t>d1>dup>d2>d1>s1>s2>dup>s4>s5
";
        std::fs::write(&gfa_path, gfa).expect("write test GFA");

        let features = vec![InputFeature {
            chrom: "chr1".into(),
            start: 22,
            end: 28,
            name: "f_dup".into(),
            score: "100".into(),
            strand: '+',
        }];
        let cfg = LiftConfig {
            gfa_path: gfa_path.to_str().expect("valid utf8 temp path"),
            from_sample: "GRCh38",
            to_sample: "CHM13",
            to_chroms: &HashSet::new(),
            features: &features,
            min_frac: 0.0,
            bridge_bp: 0,
        };
        let mut out = Vec::new();
        let summary = run_lift(&cfg, &mut out, None::<&mut Vec<u8>>).expect("run_lift succeeds");
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(summary.lifted_exact, 1);
        let out_str = String::from_utf8(out).expect("valid utf8 output");
        assert!(
            out_str.starts_with("CHM13.chr1\t62\t68\tf_dup\t100\t+"),
            "expected copy 2 at [62, 68), got: {out_str}"
        );
        assert!(out_str.contains(";DUP"), "expected ;DUP tag in: {out_str}");
    }

    #[test]
    fn lifts_across_multiple_walks_sharing_seq_id() {
        let gfa_path = test_tmp("multiwalk_test.gfa");
        let gfa = "\
H\tVN:Z:1.1
S\ts1\tACGTACGTAC
S\ts2\tTGCATGCATG
W\tGRCh38\t0\tchr1\t0\t10\t>s1
W\tGRCh38\t0\tchr1\t100\t110\t>s2
W\tCHM13\t0\tchr1\t0\t20\t>s1>s2
";
        std::fs::write(&gfa_path, gfa).expect("write test GFA");

        let features = vec![InputFeature {
            chrom: "chr1".into(),
            start: 102,
            end: 108,
            name: "f_walk2".into(),
            score: "50".into(),
            strand: '+',
        }];
        let cfg = LiftConfig {
            gfa_path: gfa_path.to_str().expect("valid utf8 temp path"),
            from_sample: "GRCh38",
            to_sample: "CHM13",
            to_chroms: &HashSet::new(),
            features: &features,
            min_frac: 0.0,
            bridge_bp: 0,
        };
        let mut out = Vec::new();
        let summary = run_lift(&cfg, &mut out, None::<&mut Vec<u8>>).expect("run_lift succeeds");
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(summary.lifted_exact, 1);
        let out_str = String::from_utf8(out).expect("valid utf8 output");
        assert!(
            out_str.starts_with("CHM13.chr1\t12\t18\tf_walk2\t50\t+"),
            "expected [12, 18) via second walk, got: {out_str}"
        );
    }

    #[test]
    fn p_path_exact_and_colon_prefix_do_not_duplicate_liftover() {
        assert_eq!(strip_sample_prefix("GRCh38.chr1", "GRCh38"), "chr1");
        assert_eq!(strip_sample_prefix("GRCh38#chr1", "GRCh38"), "chr1");
        assert_eq!(strip_sample_prefix("GRCh38:chr1", "GRCh38"), "chr1");
        assert_eq!(strip_sample_prefix("GRCh38", "GRCh38"), "GRCh38");

        let gfa_path = test_tmp("ppath_dedup_test.gfa");
        let gfa = "\
H\tVN:Z:1.0
S\ts1\tACGTACGTAC
P\tGRCh38\ts1+\t*
P\tCHM13:chr1\ts1+\t*
";
        std::fs::write(&gfa_path, gfa).expect("write test GFA");
        let features = vec![InputFeature {
            chrom: "GRCh38".into(),
            start: 2,
            end: 8,
            name: "f_exact".into(),
            score: "100".into(),
            strand: '+',
        }];
        let cfg = LiftConfig {
            gfa_path: gfa_path.to_str().expect("valid utf8 temp path"),
            from_sample: "GRCh38",
            to_sample: "CHM13",
            to_chroms: &HashSet::new(),
            features: &features,
            min_frac: 0.5,
            bridge_bp: 0,
        };
        let mut out = Vec::new();
        let summary = run_lift(&cfg, &mut out, None::<&mut Vec<u8>>).expect("run_lift succeeds");
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(summary.lifted_exact, 1);
        assert_eq!(summary.lifted_split, 0);
        let out_str = String::from_utf8(out).expect("valid utf8 output");
        assert_eq!(
            out_str.lines().count(),
            1,
            "expected 1 line, got: {out_str}"
        );
    }

    #[test]
    fn parses_gff3_and_diploid_from_walks_without_inflating_frac() {
        let gff_path = test_tmp("gff3_test.gff3");
        std::fs::write(
            &gff_path,
            "##gff-version 3\nchr1\t.\tgene\t3\t8\t.\t+\t.\tID=g1\n",
        )
        .expect("write test GFF3");
        let feats = parse_features(&gff_path, true).expect("parse GFF3");
        let _ = std::fs::remove_file(&gff_path);
        assert_eq!(feats.len(), 1);
        assert_eq!(feats[0].start, 2);
        assert_eq!(feats[0].end, 8);
        assert_eq!(feats[0].name, "g1");

        let gfa_path = test_tmp("diploid_test.gfa");
        let gfa = "\
H\tVN:Z:1.1
S\ts1\tACGTACGTAC
W\tHG002\t1\tchr1\t0\t10\t>s1
W\tHG002\t2\tchr1\t0\t10\t>s1
W\tCHM13\t0\tchr1\t0\t10\t>s1
";
        std::fs::write(&gfa_path, gfa).expect("write test diploid GFA");
        let cfg = LiftConfig {
            gfa_path: gfa_path.to_str().expect("valid utf8"),
            from_sample: "HG002",
            to_sample: "CHM13",
            to_chroms: &HashSet::new(),
            features: &feats,
            min_frac: 0.9,
            bridge_bp: 0,
        };
        let mut out = Vec::new();
        let summary = run_lift(&cfg, &mut out, None::<&mut Vec<u8>>).expect("run_lift");
        let _ = std::fs::remove_file(&gfa_path);
        assert_eq!(summary.unmapped, 0);
        assert_eq!(summary.lifted_exact, 2);
        assert_eq!(summary.lifted_split, 0);
        let out_str = String::from_utf8(out).expect("utf8");
        assert!(
            out_str.contains("\t1.0000\tEXACT"),
            "expected diploid walks to emit EXACT per walk, got: {out_str}"
        );
    }

    #[test]
    fn unstranded_feature_does_not_coalesce_across_inversion() {
        let gfa_path = test_tmp("unstranded_inv_test.gfa");
        let gfa = "\
H\tVN:Z:1.1
S\ts1\tACGTACGTAC
S\ts2\tTGCATGCATG
W\tGRCh38\t0\tchr1\t0\t20\t>s1>s2
W\tCHM13\t0\tchr1\t0\t20\t>s1<s2
";
        std::fs::write(&gfa_path, gfa).expect("write test GFA");
        let features = vec![InputFeature {
            chrom: "chr1".into(),
            start: 5,
            end: 15,
            name: "f_dot".into(),
            score: "0".into(),
            strand: '.',
        }];
        let cfg = LiftConfig {
            gfa_path: gfa_path.to_str().expect("valid utf8"),
            from_sample: "GRCh38",
            to_sample: "CHM13",
            to_chroms: &HashSet::new(),
            features: &features,
            min_frac: 0.0,
            bridge_bp: 5,
        };
        let mut out = Vec::new();
        let summary = run_lift(&cfg, &mut out, None::<&mut Vec<u8>>).expect("run_lift");
        let _ = std::fs::remove_file(&gfa_path);
        assert_eq!(summary.lifted_split, 1);
        assert_eq!(summary.lifted_exact, 0);
        let out_str = String::from_utf8(out).expect("utf8");
        assert_eq!(
            out_str.lines().count(),
            2,
            "expected 2 split pieces, got: {out_str}"
        );
    }

    #[test]
    fn flip_strand_preserves_unstranded_dot() {
        assert_eq!(flip_strand('+'), '-');
        assert_eq!(flip_strand('-'), '+');
        assert_eq!(flip_strand('.'), '.');
    }
}
