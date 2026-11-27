//! Unified `horizon` subcommand dispatcher for the Horizon genome alignment and pangenome suite.

use std::env;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

const SUBCOMMANDS: &[(&str, &str, &str)] = &[
    (
        "wga",
        "horizon_wga",
        "Pairwise WGA, DAG/POA graph alignment (--graph), and HAL info (--hal-info)",
    ),
    (
        "merge",
        "horizon_merge",
        "O(N log N) sweep-line multi-species MAF + rGFA merger",
    ),
    (
        "guide",
        "horizon_guide",
        "Progressive guide-tree WGA with transitive rescue and HAL export",
    ),
    (
        "tba",
        "horizon_tba",
        "Reference-free TBA (all-vs-all, two-tier --hinges/--auto-hinges, --progressive-ancestors, --from-mafs)",
    ),
    (
        "gfa",
        "horizon_gfa",
        "Streaming GFA v1.0/v1.1 stats, list-walks, chromosome extract, and incremental rgfa",
    ),
    (
        "hal",
        "horizon_hal",
        "Pure-Rust streaming multi-level HAL stats, liftover, to-maf, and to-gfa (zero libhdf5)",
    ),
    (
        "lift",
        "horizon_lift",
        "Graph-aware BED/GFF3 coordinate liftover with micro-bubble coalescing",
    ),
    (
        "clean",
        "horizon_clean",
        "BOCPD + MDL GFA junction scorer, bubble cutter, and haplotype/bubble filter",
    ),
    (
        "resolve",
        "horizon_resolve",
        "MDL-calibrated k-copy e-graph resolver for GFA walks and phased snarl VCFv4.2 export",
    ),
    (
        "filter",
        "horizon_filter",
        "Reference-free MDL, biological CRC, & syntenic-netting PAF alignment filter",
    ),
    (
        "locus",
        "horizon_locus",
        "Regional multi-FASTA locus aligner, BioNJ tree, --hnex, & ProbCons --consistency MSA",
    ),
    (
        "eval",
        "horizon_eval",
        "Pairwise alignment recall/precision evaluator against reference MAF",
    ),
];

fn resolve_binary_name(subcommand: &str) -> Option<&'static str> {
    SUBCOMMANDS
        .iter()
        .find(|&&(cmd, _, _)| cmd == subcommand)
        .map(|&(_, bin, _)| bin)
}

fn find_sibling_or_path(bin_name: &str, current_exe: Option<&Path>) -> PathBuf {
    let exe_name = format!("{bin_name}{}", env::consts::EXE_SUFFIX);
    if let Some(exe) = current_exe
        && let Some(dir) = exe.parent()
    {
        let sibling = dir.join(&exe_name);
        if sibling.is_file() {
            return sibling;
        }
    }
    PathBuf::from(exe_name)
}

fn write_usage<W: std::io::Write>(mut w: W) {
    let _ = writeln!(w, "horizon {}", env!("CARGO_PKG_VERSION"));
    let _ = writeln!(
        w,
        "Information-theoretic genome alignment and pangenome graph suite.\n"
    );
    let _ = writeln!(w, "USAGE:");
    let _ = writeln!(w, "    horizon <SUBCOMMAND> [ARGS...]\n");
    let _ = writeln!(w, "SUBCOMMANDS:");
    for &(cmd, bin, desc) in SUBCOMMANDS {
        let _ = writeln!(w, "    {cmd:<11} {desc} ({bin})");
    }
    let _ = writeln!(
        w,
        "\nRun `horizon <SUBCOMMAND> --help` for subcommand options."
    );
}

fn main() -> ExitCode {
    let mut args = env::args_os();
    let _program = args.next();
    let Some(first) = args.next() else {
        write_usage(std::io::stderr());
        return ExitCode::from(2);
    };
    let cmd = first.to_string_lossy();
    match cmd.as_ref() {
        "-h" | "--help" => {
            write_usage(std::io::stdout());
            ExitCode::SUCCESS
        }
        "help" => {
            if let Some(sub_os) = args.next() {
                let sub = sub_os.to_string_lossy();
                let Some(bin_name) = resolve_binary_name(sub.as_ref()) else {
                    eprintln!("Unknown subcommand: {sub}\n");
                    write_usage(std::io::stderr());
                    return ExitCode::from(2);
                };
                let exe = env::current_exe().ok();
                let target = find_sibling_or_path(bin_name, exe.as_deref());
                return match Command::new(&target).arg("--help").status() {
                    Ok(status) => status
                        .code()
                        .and_then(|c| u8::try_from(c).ok())
                        .map_or(ExitCode::FAILURE, ExitCode::from),
                    Err(err) => {
                        eprintln!("Failed to launch `{}`: {err}", target.display());
                        ExitCode::FAILURE
                    }
                };
            }
            write_usage(std::io::stdout());
            ExitCode::SUCCESS
        }
        "-V" | "--version" | "version" => {
            println!("horizon {}", env!("CARGO_PKG_VERSION"));
            ExitCode::SUCCESS
        }
        sub => {
            let Some(bin_name) = resolve_binary_name(sub) else {
                eprintln!("Unknown subcommand: {sub}\n");
                write_usage(std::io::stderr());
                return ExitCode::from(2);
            };
            let exe = env::current_exe().ok();
            let target = find_sibling_or_path(bin_name, exe.as_deref());
            match Command::new(&target).args(args).status() {
                Ok(status) => status
                    .code()
                    .and_then(|c| u8::try_from(c).ok())
                    .map_or(ExitCode::FAILURE, ExitCode::from),
                Err(err) => {
                    eprintln!("Failed to launch `{}`: {err}", target.display());
                    ExitCode::FAILURE
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolves_known_subcommands_and_rejects_unknown() {
        assert_eq!(resolve_binary_name("wga"), Some("horizon_wga"));
        assert_eq!(resolve_binary_name("gfa"), Some("horizon_gfa"));
        assert_eq!(resolve_binary_name("hal"), Some("horizon_hal"));
        assert_eq!(resolve_binary_name("lift"), Some("horizon_lift"));
        assert_eq!(resolve_binary_name("clean"), Some("horizon_clean"));
        assert_eq!(resolve_binary_name("resolve"), Some("horizon_resolve"));
        assert_eq!(resolve_binary_name("filter"), Some("horizon_filter"));
        assert_eq!(resolve_binary_name("locus"), Some("horizon_locus"));
        assert_eq!(resolve_binary_name("merge"), Some("horizon_merge"));
        assert_eq!(resolve_binary_name("guide"), Some("horizon_guide"));
        assert_eq!(resolve_binary_name("tba"), Some("horizon_tba"));
        assert_eq!(resolve_binary_name("eval"), Some("horizon_eval"));
        assert_eq!(resolve_binary_name("nonexistent"), None);
    }

    #[test]
    fn find_sibling_or_path_prefers_existing_sibling_and_falls_back() {
        let dir = std::env::temp_dir().join(format!("horizon_bin_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create temp dir");
        let exe = dir.join("horizon");
        let sibling = dir.join("horizon_wga");
        std::fs::write(&sibling, b"").expect("write dummy sibling");

        assert_eq!(find_sibling_or_path("horizon_wga", Some(&exe)), sibling);
        assert_eq!(
            find_sibling_or_path("horizon_missing", Some(&exe)),
            PathBuf::from("horizon_missing")
        );
        assert_eq!(
            find_sibling_or_path("horizon_wga", None),
            PathBuf::from("horizon_wga")
        );

        let _ = std::fs::remove_file(&sibling);
        let _ = std::fs::remove_dir(&dir);
    }
}
