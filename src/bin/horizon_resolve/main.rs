//! Two-pass streaming GFA v1.0/v1.1 snarl/bubble e-graph resolver using `MdlCostModel`.

mod resolve;
mod snarl;

use anyhow::{Context, Result, bail};
use horizon_aligner::gfa::{SpoolableGfaInput, open_output_writer};
use std::env;
use std::time::Instant;

use resolve::pass2_rewrite;
use snarl::pass1_index;

pub(crate) struct CliArgs {
    pub(crate) gfa_in: String,
    pub(crate) out_path: Option<String>,
    pub(crate) collapse_dups: bool,
    pub(crate) tag_walks: bool,
    pub(crate) max_bubble_steps: usize,
    pub(crate) min_anchor_bp: usize,
    pub(crate) divergence: Option<f64>,
}

fn parse_args() -> Result<Option<CliArgs>> {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        println!(
            "Usage: horizon_resolve <input.gfa[.gz]> [--out <resolved.gfa>] [--collapse-dups] [--tag-walks] [--max-bubble-steps <N>] [--min-anchor-bp <N>] [--divergence <F>]\n\
             \n\
             Decomposes GFA v1.0/v1.1 P/W paths into bounded snarl/bubble windows between shared anchors,\n\
             runs bounded synteny e-graph saturation with MdlCostModel extraction, and rewrites\n\
             simplified walks/paths back into valid GFA."
        );
        return Ok(None);
    }

    let mut gfa_in: Option<String> = None;
    let mut out_path: Option<String> = None;
    let mut collapse_dups = false;
    let mut tag_walks = false;
    let mut max_bubble_steps = 16usize;
    let mut min_anchor_bp = 200usize;
    let mut divergence: Option<f64> = None;

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "--out" | "-o" => {
                i += 1;
                out_path = Some(args.get(i).context("missing value for --out")?.clone());
            }
            "--collapse-dups" => collapse_dups = true,
            "--tag-walks" => tag_walks = true,
            "--max-bubble-steps" => {
                i += 1;
                max_bubble_steps = args
                    .get(i)
                    .context("missing value for --max-bubble-steps")?
                    .parse()
                    .context("invalid --max-bubble-steps")?;
                anyhow::ensure!(
                    max_bubble_steps >= 2,
                    "--max-bubble-steps must be at least 2 (got {max_bubble_steps})"
                );
            }
            "--min-anchor-bp" => {
                i += 1;
                min_anchor_bp = args
                    .get(i)
                    .context("missing value for --min-anchor-bp")?
                    .parse()
                    .context("invalid --min-anchor-bp")?;
            }
            "--divergence" => {
                i += 1;
                let d: f64 = args
                    .get(i)
                    .context("missing value for --divergence")?
                    .parse()
                    .context("invalid --divergence")?;
                anyhow::ensure!(
                    d.is_finite() && d > 0.0,
                    "--divergence must be a finite positive number, got {d}"
                );
                divergence = Some(d.clamp(1e-6, 0.75));
            }
            other if (!other.starts_with('-') || other == "-") && gfa_in.is_none() => {
                gfa_in = Some(other.to_string());
            }
            other => bail!("unknown argument: {other}"),
        }
        i += 1;
    }

    let gfa_in = gfa_in.context("missing <input.gfa[.gz]>")?;
    Ok(Some(CliArgs {
        gfa_in,
        out_path,
        collapse_dups,
        tag_walks,
        max_bubble_steps: max_bubble_steps.max(2),
        min_anchor_bp,
        divergence,
    }))
}

fn main() -> Result<()> {
    let Some(cli) = parse_args()? else {
        return Ok(());
    };
    let t0 = Instant::now();
    let spooled = SpoolableGfaInput::new(&cli.gfa_in)?;
    let effective_gfa = spooled.path();
    let graph = pass1_index(effective_gfa, cli.min_anchor_bp)?;
    let anchor_count = graph.is_anchor.iter().filter(|&&a| a).count();

    let mut writer = open_output_writer(cli.out_path.as_deref())?;
    let stats = pass2_rewrite(effective_gfa, &mut writer, &graph, &cli)?;
    writer.finish()?;

    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    eprintln!(
        "segments={} anchors={} bridges={} components={} paths={} walks={} windows={} simplified={} dups_detected={} dups_collapsed={} invs_cancelled={} mdl_bits_saved={:.1} elapsed_s={:.3} peak_rss_mb={:.1}",
        graph.id_to_seg.len(),
        anchor_count,
        graph.num_bridges,
        graph.num_components,
        stats.paths_seen,
        stats.walks_seen,
        stats.windows_evaluated,
        stats.windows_simplified,
        stats.dups_detected,
        stats.dups_collapsed,
        stats.invs_cancelled,
        stats.total_mdl_bits_saved,
        t0.elapsed().as_secs_f64(),
        hwm_mb,
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolves_tandem_dup_and_cancels_double_inversion_in_gfa_walk() -> Result<()> {
        let dir = std::env::temp_dir().join(format!("horizon_resolve_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir)?;
        let gfa_in = dir.join("in.gfa");
        let gfa_out = dir.join("out.gfa");
        std::fs::write(
            &gfa_in,
            "H\tVN:Z:1.1\n\
             S\ta\tACTGACTGACTGACTGACTG\n\
             S\tb\tGATTACAGATTACAGATTAC\n\
             S\tc\tCCCCCCCCCCCCCCCCCCCC\n\
             L\ta\t+\tb\t+\t0M\n\
             L\tb\t+\tb\t+\t0M\n\
             L\tb\t+\tc\t+\t0M\n\
             W\tHG002\t1\tchr6\t0\t80\t>a>b>b>c\n\
             W\tHG002\t2\tchr7\t0\t60\t>a<b<b>c\n",
        )?;

        let cli = CliArgs {
            gfa_in: gfa_in.to_string_lossy().to_string(),
            out_path: Some(gfa_out.to_string_lossy().to_string()),
            collapse_dups: true,
            tag_walks: true,
            max_bubble_steps: 16,
            min_anchor_bp: 100,
            divergence: None,
        };
        let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
        assert!(graph.num_bridges >= 2);
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&gfa_out)?);
        let stats = pass2_rewrite(&cli.gfa_in, &mut writer, &graph, &cli)?;
        assert_eq!(stats.walks_seen, 2);
        assert_eq!(stats.dups_detected, 2);
        assert_eq!(stats.dups_collapsed, 2);
        assert!(stats.total_mdl_bits_saved > 0.0);

        let out_text = std::fs::read_to_string(&gfa_out)?;
        assert!(out_text.contains("S\ta\tACTGACTGACTGACTGACTG"));
        assert!(out_text.contains("S\tb\tGATTACAGATTACAGATTAC"));
        assert!(out_text.contains("S\tc\tCCCCCCCCCCCCCCCCCCCC"));
        assert!(out_text.contains("W\tHG002\t1\tchr6\t0\t60\t>a>b>c\tEG:Z:"));
        assert!(out_text.contains("W\tHG002\t2\tchr7\t0\t60\t>a<b>c\tEG:Z:"));
        assert!(out_text.contains("L\ta\t+\tb\t+\t0M"));
        assert!(out_text.contains("L\tb\t+\tc\t+\t0M"));
        assert!(!out_text.contains("L\tb\t+\tb\t+\t0M"));
        let _ = std::fs::remove_dir_all(&dir);
        Ok(())
    }
}
