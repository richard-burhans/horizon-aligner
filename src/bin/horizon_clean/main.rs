//! Two-pass streaming GFA v1.0/v1.1 cleaner using BOCPD + MDL junction scoring and Horizon Bubbles.

mod engine;

use anyhow::{Context, Result, bail};
use engine::{CleanConfig, clean_gfa};
use horizon_aligner::gfa::{SpoolableGfaInput, open_output_writer};
use std::env;
use std::time::Instant;

struct CliOptions {
    gfa: String,
    out_path: Option<String>,
    threads: Option<usize>,
    cfg: CleanConfig,
}

fn parse_cli_args() -> Result<Option<CliOptions>> {
    const USAGE: &str = "usage: horizon_clean --gfa <in.gfa[.gz]|-> [--out <clean.gfa>] [--flank 128] [--min-cp-bits 12.0] [--prune-unsupported] [--min-haplotypes <K>] [--clip-bubbles-bp <N>] [-t|--threads <N>]";
    let args: Vec<String> = env::args().skip(1).collect();
    let mut gfa_path: Option<String> = None;
    let mut out_path: Option<String> = None;
    let mut threads: Option<usize> = None;
    let mut flank_bp = 128usize;
    let mut min_cp_bits = 12.0_f64;
    let mut prune_unsupported = false;
    let mut min_haplotypes = 0usize;
    let mut clip_bubbles_bp = 0usize;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                println!("{USAGE}");
                return Ok(None);
            }
            "-V" | "--version" => {
                println!("horizon_clean {}", env!("CARGO_PKG_VERSION"));
                return Ok(None);
            }
            "--gfa" | "-i" => {
                gfa_path = Some(args.get(i + 1).context("missing value for --gfa")?.clone());
                i += 2;
            }
            "--out" | "-o" => {
                out_path = Some(args.get(i + 1).context("missing value for --out")?.clone());
                i += 2;
            }
            "-t" | "--threads" => {
                threads = Some(horizon_aligner::parse_threads_arg(
                    args.get(i + 1).context("missing value for --threads")?,
                )?);
                i += 2;
            }
            "--flank" => {
                let val = args.get(i + 1).context("missing value for --flank")?;
                flank_bp = val.parse().context("parsing --flank")?;
                anyhow::ensure!(
                    flank_bp >= 12,
                    "--flank must be at least 12 bp (got {flank_bp})"
                );
                i += 2;
            }
            "--min-cp-bits" => {
                let val = args.get(i + 1).context("missing value for --min-cp-bits")?;
                min_cp_bits = val.parse().context("parsing --min-cp-bits")?;
                anyhow::ensure!(
                    min_cp_bits.is_finite() && min_cp_bits >= 0.0,
                    "--min-cp-bits must be a finite non-negative number (got {min_cp_bits})"
                );
                i += 2;
            }
            "--prune-unsupported" => {
                prune_unsupported = true;
                i += 1;
            }
            "--min-haplotypes" => {
                let val = args
                    .get(i + 1)
                    .context("missing value for --min-haplotypes")?;
                min_haplotypes = val.parse().context("parsing --min-haplotypes")?;
                i += 2;
            }
            "--clip-bubbles-bp" => {
                let val = args
                    .get(i + 1)
                    .context("missing value for --clip-bubbles-bp")?;
                clip_bubbles_bp = val.parse().context("parsing --clip-bubbles-bp")?;
                i += 2;
            }
            _ => bail!("unknown argument '{}'", args[i]),
        }
    }

    let Some(gfa) = gfa_path else {
        bail!("{USAGE}");
    };
    Ok(Some(CliOptions {
        gfa,
        out_path,
        threads,
        cfg: CleanConfig {
            flank_bp,
            min_cp_bits,
            prune_unsupported,
            min_haplotypes,
            clip_bubbles_bp,
        },
    }))
}

fn main() -> Result<()> {
    let Some(opts) = parse_cli_args()? else {
        return Ok(());
    };
    horizon_aligner::init_rayon_threads(opts.threads)?;
    let t0 = Instant::now();
    let spooled = SpoolableGfaInput::new(&opts.gfa)?;
    let effective_gfa = spooled.path();

    let mut w = open_output_writer(opts.out_path.as_deref())?;
    let stats = clean_gfa(effective_gfa, &opts.cfg, &mut w)?;
    w.finish()?;

    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    eprintln!(
        "horizon_clean: total_links={} kept={} severed={} dropped={} bubbles={} split_walks={} dropped_low_hap_segs={} clipped_bubbles={} elapsed={:.2}s peak_rss_mb={:.1}",
        stats.total_links,
        stats.kept_links,
        stats.severed_links,
        stats.dropped_links,
        stats.emitted_bubbles,
        stats.split_walks,
        stats.dropped_low_hap_segs,
        stats.clipped_bubbles,
        t0.elapsed().as_secs_f64(),
        hwm_mb
    );
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn severs_discordant_unsupported_link_and_preserves_valid_walk() {
        let gfa_path =
            std::env::temp_dir().join(format!("horizon_clean_test_{}.gfa", std::process::id()));
        // s1 (60bp AT-rich) -> s2 (60bp AT-rich): continuous, walk-supported
        // s1 (60bp AT-rich) -> s3 (60bp GC-rich): unsupported + strong compositional break
        let s1 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let s2 = "TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA";
        let s3 = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\ts2\t{s2}\n\
             S\ts3\t{s3}\n\
             L\ts1\t+\ts2\t+\t0M\n\
             L\ts1\t+\ts3\t+\t0M\n\
             W\tCHM13\t0\tchr22\t0\t120\t>s1>s2\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 60,
            min_cp_bits: 12.0,
            prune_unsupported: false,
            min_haplotypes: 0,
            clip_bubbles_bp: 0,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(stats.total_links, 2);
        assert_eq!(stats.kept_links, 1);
        assert_eq!(stats.severed_links, 1);
        assert_eq!(stats.emitted_bubbles, 1);

        let out_str = String::from_utf8(out).unwrap();
        assert!(out_str.contains("L\ts1\t+\ts2\t+\t0M"));
        assert!(!out_str.contains("L\ts1\t+\ts3\t+\t0M"));
        assert!(out_str.contains("HB:i:1"));
        assert!(out_str.contains("W\tCHM13\t0\tchr22\t0\t120\t>s1>s2"));
    }

    #[test]
    fn severs_bidirected_link_pair_symmetrically() {
        let gfa_path =
            std::env::temp_dir().join(format!("horizon_clean_sym_test_{}.gfa", std::process::id()));
        let s1 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let s2 = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\ts2\t{s2}\n\
             L\ts1\t+\ts2\t+\t0M\n\
             L\ts2\t-\ts1\t-\t0M\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 60,
            min_cp_bits: 12.0,
            prune_unsupported: false,
            min_haplotypes: 0,
            clip_bubbles_bp: 0,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(stats.total_links, 2);
        assert_eq!(stats.kept_links, 0);
        assert_eq!(stats.severed_links, 2);
        let out_str = String::from_utf8(out).unwrap();
        assert!(!out_str.contains("\nL\t"));
    }

    #[test]
    fn splits_path_and_walk_across_severed_link() {
        let gfa_path = std::env::temp_dir().join(format!(
            "horizon_clean_split_test_{}.gfa",
            std::process::id()
        ));
        let s1 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let s2 = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\ts2\t{s1}\n\
             S\ts3\t{s2}\n\
             S\ts4\t{s2}\n\
             L\ts1\t+\ts2\t+\t0M\n\
             L\ts2\t+\ts3\t+\t0M\n\
             L\ts2\t+\ts4\t+\t0M\n\
             P\tp1\ts1+,s2+,s3+\t0M,0M\n\
             W\tCHM13\t0\tchr1\t0\t180\t>s1>s2>s4\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 60,
            min_cp_bits: 12.0,
            prune_unsupported: false,
            min_haplotypes: 0,
            clip_bubbles_bp: 0,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(stats.split_walks, 2);
        let out_str = String::from_utf8(out).unwrap();
        assert!(out_str.contains("P\tp1_part1\ts1+,s2+\t0M\n"));
        assert!(out_str.contains("P\tp1_part2\ts3+\t*\n"));
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t0\t120\t>s1>s2\n"));
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t120\t180\t>s4\n"));
    }

    #[test]
    fn preserves_walk_without_explicit_l_lines() {
        let gfa_path = std::env::temp_dir().join(format!(
            "horizon_clean_no_l_test_{}.gfa",
            std::process::id()
        ));
        let s1 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\ts2\t{s1}\n\
             W\tCHM13\t0\tchr1\t0\t120\t>s1>s2\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 60,
            min_cp_bits: 12.0,
            prune_unsupported: true,
            min_haplotypes: 0,
            clip_bubbles_bp: 0,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(stats.split_walks, 0);
        let out_str = String::from_utf8(out).unwrap();
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t0\t120\t>s1>s2\n"));
    }

    #[test]
    fn prune_unsupported_severs_low_cp_unsupported_link() {
        let gfa_path = std::env::temp_dir().join(format!(
            "horizon_clean_prune_test_{}.gfa",
            std::process::id()
        ));
        let s1 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\ts2\t{s1}\n\
             L\ts1\t+\ts2\t+\t0M\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 60,
            min_cp_bits: 100.0,
            prune_unsupported: true,
            min_haplotypes: 0,
            clip_bubbles_bp: 0,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(stats.total_links, 1);
        assert_eq!(stats.kept_links, 0);
        assert_eq!(stats.severed_links, 1);
        let out_str = String::from_utf8(out).unwrap();
        assert!(!out_str.contains("\nL\t"));
        assert!(out_str.contains("HB:i:1"));
    }

    #[test]
    fn walk_split_accounts_for_link_overlap_cigar() {
        let gfa_path =
            std::env::temp_dir().join(format!("horizon_clean_ov_test_{}.gfa", std::process::id()));
        let s1 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let s2 = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
        let s3 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\ts2\t{s2}\n\
             S\ts3\t{s3}\n\
             L\ts1\t+\ts2\t+\t10M\n\
             L\ts2\t+\ts3\t+\t10M\n\
             W\tCHM13\t0\tchr1\t0\t160\t>s1>s2>s3\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 60,
            min_cp_bits: 5.0,
            prune_unsupported: false,
            min_haplotypes: 0,
            clip_bubbles_bp: 0,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(stats.split_walks, 1);
        let out_str = String::from_utf8(out).unwrap();
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t0\t60\t>s1\n"));
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t50\t110\t>s2\n"));
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t100\t160\t>s3\n"));
    }

    #[test]
    fn min_haplotypes_filters_private_segments_links_and_splits_walks() {
        let gfa_path = std::env::temp_dir().join(format!(
            "horizon_clean_minhap_{}_{:?}.gfa",
            std::process::id(),
            std::thread::current().id()
        ));
        let s1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 60 bp
        let s_priv = "TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAA"; // 40 bp (only HG002#1)
        let s3 = "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"; // 60 bp
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\ts_priv\t{s_priv}\n\
             S\ts3\t{s3}\n\
             L\ts1\t+\ts_priv\t+\t0M\n\
             L\ts_priv\t+\ts3\t+\t0M\n\
             L\ts1\t+\ts3\t+\t0M\n\
             W\tCHM13\t0\tchr1\t0\t120\t>s1>s3\n\
             W\tHG002\t1\tchr1\t0\t160\t>s1>s_priv>s3\n\
             W\tHG002\t2\tchr1\t0\t120\t>s1>s3\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 60,
            min_cp_bits: 100.0,
            prune_unsupported: false,
            min_haplotypes: 2,
            clip_bubbles_bp: 0,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        let out_str = String::from_utf8(out).unwrap();
        assert_eq!(
            stats.dropped_low_hap_segs, 1,
            "expected s_priv dropped:\n{out_str}"
        );
        assert!(out_str.contains("S\ts1\t"));
        assert!(out_str.contains("S\ts3\t"));
        assert!(
            !out_str.contains("S\ts_priv\t"),
            "private segment must be dropped"
        );
        assert!(
            !out_str.contains("s_priv"),
            "no links or walks may reference dropped s_priv"
        );
        assert!(out_str.contains("L\ts1\t+\ts3\t+\t0M\n"));
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t0\t120\t>s1>s3\n"));
        assert!(out_str.contains("W\tHG002\t2\tchr1\t0\t120\t>s1>s3\n"));
        // HG002#1 walk must split cleanly around dropped 40bp s_priv with exact genomic coordinates!
        assert!(
            out_str.contains("W\tHG002\t1\tchr1\t0\t60\t>s1\n")
                && out_str.contains("W\tHG002\t1\tchr1\t100\t160\t>s3\n"),
            "expected HG002#1 walk split around dropped s_priv at [0..60) and [100..160), got:\n{out_str}"
        );
    }

    #[test]
    fn clip_bubbles_bp_drops_large_horizon_bubbles_and_keeps_small_ones() {
        let gfa_path = std::env::temp_dir().join(format!(
            "horizon_clean_clip_{}_{:?}.gfa",
            std::process::id(),
            std::thread::current().id()
        ));
        let s1 = "ACGT".repeat(30); // 120 bp normal segment (must NOT be clipped)
        let b_small = "TGCA".repeat(12); // 48 bp HB:i:1 (< 100 bp, must be kept)
        let b_big = "CGTA".repeat(30); // 120 bp HB:i:1 (>= 100 bp, must be clipped)
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\tb_big\t{b_big}\tHB:i:1\n\
             S\tb_small\t{b_small}\tHB:i:1\n\
             L\ts1\t+\tb_big\t+\t0M\n\
             L\tb_big\t+\tb_small\t+\t0M\n\
             W\tCHM13\t0\tchr1\t0\t288\t>s1>b_big>b_small\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 48,
            min_cp_bits: 100.0,
            prune_unsupported: false,
            min_haplotypes: 0,
            clip_bubbles_bp: 100,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        let out_str = String::from_utf8(out).unwrap();
        assert_eq!(
            stats.clipped_bubbles, 1,
            "expected b_big clipped:\n{out_str}"
        );
        assert!(out_str.contains("S\ts1\t"));
        assert!(out_str.contains("S\tb_small\t"));
        assert!(
            !out_str.contains("b_big"),
            "clipped large HB:i:1 bubble must be completely removed"
        );
        assert!(
            out_str.contains("W\tCHM13\t0\tchr1\t0\t120\t>s1\n")
                && out_str.contains("W\tCHM13\t0\tchr1\t240\t288\t>b_small\n"),
            "expected walk split cleanly around clipped 120bp bubble, got:\n{out_str}"
        );
    }

    #[test]
    fn default_zero_thresholds_preserve_private_segments_and_large_bubbles_losslessly() {
        let gfa_path = std::env::temp_dir().join(format!(
            "horizon_clean_lossless_{}_{:?}.gfa",
            std::process::id(),
            std::thread::current().id()
        ));
        let s1 = "ACGT".repeat(30);
        let b_big = "CGTA".repeat(30);
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\tb_big\t{b_big}\tHB:i:1\n\
             L\ts1\t+\tb_big\t+\t0M\n\
             W\tCHM13\t0\tchr1\t0\t240\t>s1>b_big\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 60,
            min_cp_bits: 100.0,
            prune_unsupported: false,
            min_haplotypes: 0,
            clip_bubbles_bp: 0,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        let out_str = String::from_utf8(out).unwrap();
        assert_eq!(stats.dropped_low_hap_segs, 0);
        assert_eq!(stats.clipped_bubbles, 0);
        assert!(out_str.contains("S\tb_big\t"));
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t0\t240\t>s1>b_big\n"));
    }
}
