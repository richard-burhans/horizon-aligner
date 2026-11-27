//! Two-pass streaming GFA v1.0/v1.1 snarl/bubble e-graph resolver using `MdlCostModel`.

mod resolve;
mod snarl;
#[cfg(test)]
mod tests_nested;
mod vcf;
mod vcf_nested;

use anyhow::{Context, Result, bail};
use horizon_aligner::gfa::{SpoolableGfaInput, open_output_writer};
use std::env;
use std::time::Instant;

use resolve::pass2_rewrite;
use snarl::pass1_index;
use vcf::export_snarl_vcf;

pub(crate) struct CliArgs {
    pub(crate) gfa_in: String,
    pub(crate) out_path: Option<String>,
    pub(crate) vcf_path: Option<String>,
    pub(crate) vcf_ref: Option<String>,
    pub(crate) nested_snarls: bool,
    pub(crate) collapse_dups: bool,
    pub(crate) tag_walks: bool,
    pub(crate) max_bubble_steps: usize,
    pub(crate) min_anchor_bp: usize,
    pub(crate) divergence: Option<f64>,
    pub(crate) threads: Option<usize>,
}

#[allow(clippy::too_many_lines)]
fn parse_args() -> Result<Option<CliArgs>> {
    let args: Vec<String> = env::args().skip(1).collect();
    let usage = "Usage: horizon_resolve <input.gfa[.gz]> [--out <resolved.gfa>] [--vcf <out.vcf[.gz]>] [--vcf-ref <sample>] [--nested-snarls] [--collapse-dups] [--tag-walks] [--max-bubble-steps <N>] [--min-anchor-bp <N>] [--divergence <F>] [-t|--threads <N>]\n\
                 \n\
                 Decomposes GFA v1.0/v1.1 P/W paths into bounded snarl/bubble windows between shared anchors,\n\
                 runs bounded synteny e-graph saturation with MdlCostModel extraction, rewrites\n\
                 simplified walks/paths back into valid GFA, and optionally projects Tarjan snarls (plus optional\n\
                 LV=1 nested reference-anchored bubbles via --nested-snarls) to VCFv4.2.";
    if args.iter().any(|a| a == "-h" || a == "--help") {
        println!("{usage}");
        return Ok(None);
    }
    if args.is_empty() {
        bail!("{usage}");
    }

    let mut gfa_in: Option<String> = None;
    let mut out_path: Option<String> = None;
    let mut vcf_path: Option<String> = None;
    let mut vcf_ref: Option<String> = None;
    let mut nested_snarls = false;
    let mut collapse_dups = false;
    let mut tag_walks = false;
    let mut max_bubble_steps = 16usize;
    let mut min_anchor_bp = 200usize;
    let mut divergence: Option<f64> = None;
    let mut threads: Option<usize> = None;

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-V" | "--version" => {
                println!("horizon_resolve {}", env!("CARGO_PKG_VERSION"));
                return Ok(None);
            }
            "--out" | "-o" => {
                i += 1;
                out_path = Some(args.get(i).context("missing value for --out")?.clone());
            }
            "--vcf" => {
                i += 1;
                vcf_path = Some(args.get(i).context("missing value for --vcf")?.clone());
            }
            "--vcf-ref" => {
                i += 1;
                vcf_ref = Some(args.get(i).context("missing value for --vcf-ref")?.clone());
            }
            "--nested-snarls" => nested_snarls = true,
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
            "-t" | "--threads" => {
                i += 1;
                threads = Some(horizon_aligner::parse_threads_arg(
                    args.get(i).context("missing value for --threads")?,
                )?);
            }
            other if (!other.starts_with('-') || other == "-") && gfa_in.is_none() => {
                gfa_in = Some(other.to_string());
            }
            other => bail!("unknown argument: {other}"),
        }
        i += 1;
    }

    let gfa_in = gfa_in.context("missing <input.gfa[.gz]>")?;
    let gfa_is_stdout = out_path.as_deref().is_none_or(|p| p == "-");
    anyhow::ensure!(
        !(gfa_is_stdout && vcf_path.as_deref() == Some("-")),
        "cannot write both resolved GFA and --vcf to stdout ('-')"
    );
    anyhow::ensure!(
        vcf_ref.is_none() || vcf_path.is_some(),
        "--vcf-ref requires --vcf <out.vcf[.gz]>"
    );
    anyhow::ensure!(
        !nested_snarls || vcf_path.is_some(),
        "--nested-snarls requires --vcf <out.vcf[.gz]>"
    );
    Ok(Some(CliArgs {
        gfa_in,
        out_path,
        vcf_path,
        vcf_ref,
        nested_snarls,
        collapse_dups,
        tag_walks,
        max_bubble_steps,
        min_anchor_bp,
        divergence,
        threads,
    }))
}

fn main() -> Result<()> {
    let Some(cli) = parse_args()? else {
        return Ok(());
    };
    horizon_aligner::init_rayon_threads(cli.threads)?;
    let t0 = Instant::now();
    let spooled = SpoolableGfaInput::new(&cli.gfa_in)?;
    let effective_gfa = spooled.path();
    let graph = pass1_index(effective_gfa, cli.min_anchor_bp)?;
    let anchor_count = graph.is_anchor.iter().filter(|&&a| a).count();

    let mut writer = open_output_writer(cli.out_path.as_deref())?;
    let mut vcf_capture = cli
        .vcf_path
        .as_ref()
        .map(|_| resolve::VcfCapture::with_segments(graph.num_segs()));
    let stats = pass2_rewrite(
        effective_gfa,
        &mut writer,
        &graph,
        &cli,
        vcf_capture.as_mut(),
    )?;
    writer.finish()?;

    let mut vcf_summary = String::new();
    if let (Some(vcf_path), Some(capture)) = (cli.vcf_path.as_deref(), vcf_capture.as_ref()) {
        let mut vcf_writer = open_output_writer(Some(vcf_path))?;
        let vstats = export_snarl_vcf(
            &graph,
            &capture.seg_seqs,
            &capture.walks,
            cli.vcf_ref.as_deref(),
            cli.nested_snarls,
            cli.max_bubble_steps.max(64),
            &mut vcf_writer,
        )?;
        vcf_writer.finish()?;
        vcf_summary = format!(
            " vcf_sites={} vcf_lv1={} vcf_snps={} vcf_mnvs={} vcf_indels={} vcf_multi={} vcf_samples={}",
            vstats.sites_emitted,
            vstats.nested_sites,
            vstats.snps,
            vstats.mnvs,
            vstats.indels,
            vstats.multiallelic,
            vstats.samples
        );
    }

    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    eprintln!(
        "segments={} anchors={} bridges={} components={} paths={} walks={} windows={} simplified={} dups_detected={} dups_collapsed={} invs_factored={} mdl_bits_saved={:.1}{} elapsed_s={:.3} peak_rss_mb={:.1}",
        graph.num_segs(),
        anchor_count,
        graph.num_bridges,
        graph.num_components,
        stats.paths_seen,
        stats.walks_seen,
        stats.windows_evaluated,
        stats.windows_simplified,
        stats.dups_detected,
        stats.dups_collapsed,
        stats.invs_factored,
        stats.total_mdl_bits_saved,
        vcf_summary,
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
            vcf_path: None,
            vcf_ref: None,
            nested_snarls: false,
            collapse_dups: true,
            tag_walks: true,
            max_bubble_steps: 16,
            min_anchor_bp: 100,
            divergence: None,
            threads: None,
        };
        let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
        assert!(graph.num_bridges >= 2);
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&gfa_out)?);
        let stats = pass2_rewrite(&cli.gfa_in, &mut writer, &graph, &cli, None)?;
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

    #[test]
    fn tag_walks_and_mdl_savings_work_without_collapse_dups() -> Result<()> {
        let dir =
            std::env::temp_dir().join(format!("horizon_resolve_nocollapse_{}", std::process::id()));
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
             W\tHG002\t1\tchr6\t0\t80\t>a>b>b>c\n",
        )?;

        let cli = CliArgs {
            gfa_in: gfa_in.to_string_lossy().to_string(),
            out_path: Some(gfa_out.to_string_lossy().to_string()),
            vcf_path: None,
            vcf_ref: None,
            nested_snarls: false,
            collapse_dups: false,
            tag_walks: true,
            max_bubble_steps: 16,
            min_anchor_bp: 100,
            divergence: None,
            threads: None,
        };
        let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&gfa_out)?);
        let stats = pass2_rewrite(&cli.gfa_in, &mut writer, &graph, &cli, None)?;
        assert_eq!(stats.dups_detected, 1);
        assert_eq!(stats.dups_collapsed, 0);
        assert_eq!(stats.windows_simplified, 1);
        assert!(stats.total_mdl_bits_saved > 0.0);

        let out_text = std::fs::read_to_string(&gfa_out)?;
        assert!(out_text.contains("W\tHG002\t1\tchr6\t0\t80\t>a>b>b>c\tEG:Z:"));
        let _ = std::fs::remove_dir_all(&dir);
        Ok(())
    }

    #[test]
    fn overlap_walk_subtracts_link_overlaps_and_dup_inv_does_not_inflate_invs_cancelled()
    -> Result<()> {
        let dir = std::env::temp_dir().join(format!(
            "horizon_resolve_ov_{}_{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir)?;
        let gfa_in = dir.join("in.gfa");
        let gfa_out = dir.join("out.gfa");
        std::fs::write(
            &gfa_in,
            "H\tVN:Z:1.1\n\
             S\ta\tACTGACTGACTGACTGACTG\n\
             S\tb\tGATTACAGATTACAGATTAC\n\
             S\tc\tCCCCCCCCCCCCCCCCCCCC\n\
             L\ta\t+\tb\t+\t5M\n\
             L\tb\t+\tb\t+\t5M\n\
             L\tb\t+\tc\t+\t5M\n\
             W\tHG002\t1\tchr6\t0\t65\t>a>b>b>c\n\
             W\tHG002\t2\tchr7\t0\t80\t>a<b<b>c\n",
        )?;
        let cli = CliArgs {
            gfa_in: gfa_in.to_string_lossy().to_string(),
            out_path: Some(gfa_out.to_string_lossy().to_string()),
            vcf_path: None,
            vcf_ref: None,
            nested_snarls: false,
            collapse_dups: true,
            tag_walks: true,
            max_bubble_steps: 16,
            min_anchor_bp: 100,
            divergence: None,
            threads: None,
        };
        let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&gfa_out)?);
        let stats = pass2_rewrite(&cli.gfa_in, &mut writer, &graph, &cli, None)?;
        drop(writer);
        let out_text = std::fs::read_to_string(&gfa_out)?;
        let _ = std::fs::remove_dir_all(&dir);

        assert_eq!(
            stats.invs_factored, 0,
            "collapsing <b<b into dup(inv(b)) is a duplication, not an inversion factoring"
        );
        assert!(
            out_text.contains("W\tHG002\t1\tchr6\t0\t50\t>a>b>c\tEG:Z:"),
            "expected overlap-adjusted end_pos=50 (60 - 5 - 5), got:\n{out_text}"
        );
        Ok(())
    }

    #[test]
    fn collapses_three_copy_tandem_array_into_rep_3_in_gfa_walk() -> Result<()> {
        let dir = std::env::temp_dir().join(format!(
            "horizon_resolve_rep3_{}_{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
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
             W\tHG002\t1\tchr6\t0\t100\t>a>b>b>b>c\n",
        )?;
        let cli = CliArgs {
            gfa_in: gfa_in.to_string_lossy().to_string(),
            out_path: Some(gfa_out.to_string_lossy().to_string()),
            vcf_path: None,
            vcf_ref: None,
            nested_snarls: false,
            collapse_dups: true,
            tag_walks: true,
            max_bubble_steps: 16,
            min_anchor_bp: 100,
            divergence: None,
            threads: None,
        };
        let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&gfa_out)?);
        let stats = pass2_rewrite(&cli.gfa_in, &mut writer, &graph, &cli, None)?;
        drop(writer);
        let out_text = std::fs::read_to_string(&gfa_out)?;
        let _ = std::fs::remove_dir_all(&dir);

        assert_eq!(
            stats.dups_detected, 2,
            "3-copy repeat collapses 2 extra copies"
        );
        assert_eq!(stats.dups_collapsed, 2);
        assert!(
            out_text.contains("rep(b,3)") && out_text.contains("W\tHG002\t1\tchr6\t0\t60\t>a>b>c"),
            "expected rep(b,3) tag and 60 bp collapsed walk, got:\n{out_text}"
        );
        Ok(())
    }

    #[test]
    fn vcf_export_emits_snp_indel_multiallelic_phased_diploid_and_revcomp_genotypes() -> Result<()>
    {
        let dir = std::env::temp_dir().join(format!(
            "horizon_resolve_vcf_{}_{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir)?;
        let gfa_in = dir.join("in.gfa");
        // Snarl 1 (a -> c): multi-allelic SNP (ref=A, alt1=G, alt2=T) at 1-based POS=11
        // Snarl 2 (c -> e): indel (ref="" direct c->e, alt="GGGA" via ins) -> anchor base C from c at 1-based POS=21, REF=C, ALT=CGGGA
        // Snarl 3 (e -> f): invariant (every walk goes e -> f directly) -> must NOT emit a VCF row
        std::fs::write(
            &gfa_in,
            "H\tVN:Z:1.1\n\
             S\ta\tACGTACGTAC\n\
             S\tb_ref\tA\n\
             S\tb_g\tG\n\
             S\tb_t\tT\n\
             S\tc\tTTTTTCCCCC\n\
             S\tins\tGGGA\n\
             S\te\tACACACACAC\n\
             S\tf\tGTGTGTGTGT\n\
             L\ta\t+\tb_ref\t+\t0M\n\
             L\ta\t+\tb_g\t+\t0M\n\
             L\ta\t+\tb_t\t+\t0M\n\
             L\tb_ref\t+\tc\t+\t0M\n\
             L\tb_g\t+\tc\t+\t0M\n\
             L\tb_t\t+\tc\t+\t0M\n\
             L\tc\t+\te\t+\t0M\n\
             L\tc\t+\tins\t+\t0M\n\
             L\tins\t+\te\t+\t0M\n\
             L\te\t+\tf\t+\t0M\n\
             W\tCHM13\t0\tchr1\t0\t41\t>a>b_ref>c>e>f\n\
             W\tHG002\t1\tchr1\t0\t45\t>a>b_g>c>ins>e>f\n\
             W\tHG002\t2\tchr1\t0\t41\t>a>b_t>c>e>f\n\
             W\tHG003\t0\tchr1\t0\t45\t<f<e<ins<c<b_g<a\n",
        )?;
        let cli = CliArgs {
            gfa_in: gfa_in.to_string_lossy().to_string(),
            out_path: None,
            vcf_path: Some("out.vcf".to_owned()),
            vcf_ref: Some("CHM13".to_owned()),
            nested_snarls: false,
            collapse_dups: false,
            tag_walks: false,
            max_bubble_steps: 16,
            min_anchor_bp: 100,
            divergence: None,
            threads: None,
        };
        let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
        let mut gfa_buf = Vec::new();
        let mut cap = resolve::VcfCapture::with_segments(graph.num_segs());
        let _ = pass2_rewrite(&cli.gfa_in, &mut gfa_buf, &graph, &cli, Some(&mut cap))?;

        let mut vcf_buf = Vec::new();
        let vstats = export_snarl_vcf(
            &graph,
            &cap.seg_seqs,
            &cap.walks,
            Some("CHM13"),
            false,
            64,
            &mut vcf_buf,
        )?;
        let _ = std::fs::remove_dir_all(&dir);

        let vcf_text = String::from_utf8(vcf_buf)?;
        assert_eq!(
            vstats.sites_emitted, 2,
            "only the 2 polymorphic snarls should emit:\n{vcf_text}"
        );
        assert_eq!(vstats.snps, 1);
        assert_eq!(vstats.indels, 1);
        assert_eq!(vstats.multiallelic, 1);
        assert_eq!(vstats.samples, 3);
        assert!(vcf_text.contains("##fileformat=VCFv4.2\n"));
        assert!(vcf_text.contains(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCHM13\tHG002\tHG003\n"
        ));

        let data_lines: Vec<&str> = vcf_text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(
            data_lines.len(),
            2,
            "expected 2 VCF records, got:\n{vcf_text}"
        );

        // Site 1: multi-allelic SNP at POS 11 (after 10bp anchor `a`), REF=A, ALT=G,T
        let f1: Vec<&str> = data_lines[0].split('\t').collect();
        assert_eq!(f1[0], "chr1");
        assert_eq!(f1[1], "11");
        assert_eq!(f1[3], "A");
        assert_eq!(f1[4], "G,T");
        assert_eq!(f1[8], "GT");
        assert_eq!(f1[9], "0"); // CHM13 (haploid ref)
        assert_eq!(f1[10], "1|2"); // HG002 (phased diploid G|T)
        assert_eq!(f1[11], "1"); // HG003 (reverse-complement traversal oriented to ref -> G)

        // Site 2: insertion at POS 21 (last base 'C' of 10bp anchor `c` ending at 10+1+10=21), REF=C, ALT=CGGGA
        let f2: Vec<&str> = data_lines[1].split('\t').collect();
        assert_eq!(f2[0], "chr1");
        assert_eq!(f2[1], "21");
        assert_eq!(f2[3], "C");
        assert_eq!(f2[4], "CGGGA");
        assert_eq!(f2[9], "0"); // CHM13
        assert_eq!(f2[10], "1|0"); // HG002 (hap1 has ins, hap2 is ref)
        assert_eq!(f2[11], "1"); // HG003 (reverse-complement traversal has ins!)
        Ok(())
    }

    #[test]
    fn vcf_export_errors_cleanly_when_requested_ref_sample_missing() -> Result<()> {
        let dir = std::env::temp_dir().join(format!(
            "horizon_resolve_vcf_missing_{}_{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir)?;
        let gfa_in = dir.join("in.gfa");
        std::fs::write(
            &gfa_in,
            "H\tVN:Z:1.1\n\
             S\ta\tACGTACGTAC\n\
             W\tCHM13\t0\tchr1\t0\t10\t>a\n",
        )?;
        let cli = CliArgs {
            gfa_in: gfa_in.to_string_lossy().to_string(),
            out_path: None,
            vcf_path: Some("out.vcf".to_owned()),
            vcf_ref: Some("NO_SUCH_SAMPLE".to_owned()),
            nested_snarls: false,
            collapse_dups: false,
            tag_walks: false,
            max_bubble_steps: 16,
            min_anchor_bp: 100,
            divergence: None,
            threads: None,
        };
        let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
        let mut gfa_buf = Vec::new();
        let mut cap = resolve::VcfCapture::with_segments(graph.num_segs());
        let _ = pass2_rewrite(&cli.gfa_in, &mut gfa_buf, &graph, &cli, Some(&mut cap))?;
        let mut vcf_buf = Vec::new();
        let err = export_snarl_vcf(
            &graph,
            &cap.seg_seqs,
            &cap.walks,
            Some("NO_SUCH_SAMPLE"),
            false,
            64,
            &mut vcf_buf,
        );
        let _ = std::fs::remove_dir_all(&dir);
        let msg = format!("{:#}", err.expect_err("expected missing ref sample error"));
        assert!(
            msg.contains("NO_SUCH_SAMPLE"),
            "error should name the missing sample, got: {msg}"
        );
        Ok(())
    }

    #[test]
    fn vcf_export_skips_star_segments_trims_closing_overlap_and_marks_multivisit_missing()
    -> Result<()> {
        let dir = std::env::temp_dir().join(format!(
            "horizon_resolve_vcf_edge_{}_{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir)?;
        let gfa_in = dir.join("in.gfa");
        // Snarl 1 (a -> c): non-0M 2M closing overlap into `c` (`GG` shared at end of interior).
        // Snarl 2 (c -> e): alt haplotype visits a `*` placeholder bubble (must become `.`, not a deletion!).
        // HG004 traverses (a -> c) twice with different alleles -> must emit `.` (Conflict).
        std::fs::write(
            &gfa_in,
            "H\tVN:Z:1.1\n\
             S\ta\tACGTACGTAC\n\
             S\tx_ref\tTTGG\n\
             S\tx_alt\tCCGG\n\
             S\tc\tGGTTTCCCCC\n\
             S\tbub\t*\tLN:i:50\tHB:i:1\n\
             S\te\tACACACACAC\n\
             S\tf\tGTGTGTGTGT\n\
             L\ta\t+\tx_ref\t+\t0M\n\
             L\ta\t+\tx_alt\t+\t0M\n\
             L\tx_ref\t+\tc\t+\t2M\n\
             L\tx_alt\t+\tc\t+\t2M\n\
             L\tc\t+\te\t+\t0M\n\
             L\tc\t+\tbub\t+\t0M\n\
             L\tbub\t+\te\t+\t0M\n\
             L\te\t+\tf\t+\t0M\n\
             W\tCHM13\t0\tchr1\t0\t42\t>a>x_ref>c>e>f\n\
             W\tHG002\t0\tchr1\t0\t92\t>a>x_alt>c>bub>e>f\n\
             W\tHG004\t0\tchr1_p1\t0\t42\t>a>x_ref>c>e>f\n\
             W\tHG004\t0\tchr1_p2\t0\t42\t>a>x_alt>c>e>f\n",
        )?;
        let cli = CliArgs {
            gfa_in: gfa_in.to_string_lossy().to_string(),
            out_path: None,
            vcf_path: Some("out.vcf".to_owned()),
            vcf_ref: Some("CHM13".to_owned()),
            nested_snarls: false,
            collapse_dups: false,
            tag_walks: false,
            max_bubble_steps: 16,
            min_anchor_bp: 100,
            divergence: None,
            threads: None,
        };
        let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
        let mut gfa_buf = Vec::new();
        let mut cap = resolve::VcfCapture::with_segments(graph.num_segs());
        let _ = pass2_rewrite(&cli.gfa_in, &mut gfa_buf, &graph, &cli, Some(&mut cap))?;
        let mut vcf_buf = Vec::new();
        let vstats = export_snarl_vcf(
            &graph,
            &cap.seg_seqs,
            &cap.walks,
            Some("CHM13"),
            false,
            64,
            &mut vcf_buf,
        )?;
        let _ = std::fs::remove_dir_all(&dir);
        let vcf_text = String::from_utf8(vcf_buf)?;
        // Only Snarl 1 emits (Snarl 2 has only missing `*` alt -> invariant -> skipped).
        assert_eq!(vstats.sites_emitted, 1, "vcf:\n{vcf_text}");
        let rec: Vec<&str> = vcf_text
            .lines()
            .find(|l| !l.starts_with('#'))
            .expect("record")
            .split('\t')
            .collect();
        // 2M closing overlap `GG` trimmed from `TTGG`/`CCGG` -> REF=`TT`, ALT=`CC`.
        assert_eq!(rec[3], "TT", "closing 2M overlap must be trimmed from REF");
        assert_eq!(rec[4], "CC", "closing 2M overlap must be trimmed from ALT");
        assert_eq!(rec[9], "0"); // CHM13
        assert_eq!(rec[10], "1"); // HG002
        assert_eq!(rec[11], ".", "multi-visit HG004 must emit missing GT '.'");
        Ok(())
    }

    #[test]
    fn vcf_export_preserves_reference_coordinates_under_collapse_dups() -> Result<()> {
        let dir = std::env::temp_dir().join(format!(
            "horizon_resolve_vcf_dup_pos_{}_{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir)?;
        let gfa_in = dir.join("in.gfa");
        // Reference walk has tandem duplication >b>b (10 + 10 bp) BEFORE a downstream SNP (r vs q)
        // at 1-based position 10 (a) + 20 (b+b) + 10 (c) + 1 = 41. Even when `--collapse-dups`
        // rewrites the GFA walk to >a>b>c>r>d, the VCF POS must remain 41 (never shifted to 31).
        std::fs::write(
            &gfa_in,
            "H\tVN:Z:1.1\n\
             S\ta\tAAAAAAAAAA\n\
             S\tb\tCCCCCCCCCC\n\
             S\tc\tGGGGGGGGGG\n\
             S\tr\tA\n\
             S\tq\tT\n\
             S\td\tTTTTTTTTTT\n\
             L\ta\t+\tb\t+\t0M\n\
             L\tb\t+\tb\t+\t0M\n\
             L\tb\t+\tc\t+\t0M\n\
             L\tc\t+\tr\t+\t0M\n\
             L\tc\t+\tq\t+\t0M\n\
             L\tr\t+\td\t+\t0M\n\
             L\tq\t+\td\t+\t0M\n\
             P\tREF#1#chr1\ta+,b+,b+,c+,r+,d+\t*\n\
             P\tALT#1#chr1\ta+,b+,b+,c+,q+,d+\t*\n",
        )?;
        let cli = CliArgs {
            gfa_in: gfa_in.to_string_lossy().to_string(),
            out_path: None,
            vcf_path: Some("out.vcf".to_owned()),
            vcf_ref: Some("REF".to_owned()),
            nested_snarls: false,
            collapse_dups: true,
            tag_walks: false,
            max_bubble_steps: 16,
            min_anchor_bp: 5,
            divergence: Some(0.01),
            threads: None,
        };
        let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
        let mut gfa_buf = Vec::new();
        let mut cap = resolve::VcfCapture::with_segments(graph.num_segs());
        let stats = pass2_rewrite(&cli.gfa_in, &mut gfa_buf, &graph, &cli, Some(&mut cap))?;
        assert!(
            stats.dups_collapsed >= 1,
            "expected tandem duplication to collapse in GFA rewrite"
        );
        let mut vcf_buf = Vec::new();
        let vstats = export_snarl_vcf(
            &graph,
            &cap.seg_seqs,
            &cap.walks,
            Some("REF"),
            false,
            64,
            &mut vcf_buf,
        )?;
        let _ = std::fs::remove_dir_all(&dir);
        let vcf_text = String::from_utf8(vcf_buf)?;
        assert_eq!(vstats.sites_emitted, 1, "vcf:\n{vcf_text}");
        let rec: Vec<&str> = vcf_text
            .lines()
            .find(|l| !l.starts_with('#'))
            .expect("record")
            .split('\t')
            .collect();
        assert_eq!(
            rec[1], "41",
            "downstream SNP POS must stay at uncollapsed reference coord 41, got:\n{vcf_text}"
        );
        assert_eq!(rec[3], "A");
        assert_eq!(rec[4], "T");
        Ok(())
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn vcf_export_nested_snarls_calls_internal_snps_and_indels_inside_spanning_deletion()
    -> Result<()> {
        let dir = std::env::temp_dir().join(format!(
            "horizon_resolve_nested_{}_{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir)?;
        let gfa_in = dir.join("in.gfa");
        // Spanning deletion in DEL_HAP jumps directly a -> z, breaking global Tarjan bridges
        // across m1..m2..m3. Inside that spanning deletion:
        // - Internal biallelic SNP (r_snp=A vs a_snp=G, both shared by >=2 haplotypes!) at POS=16
        // - Internal 2-bp insertion (ins=TT) between m2 and m3 at POS=21 (anchor base G from m2)
        // - REV_HAP traverses the entire locus in reverse-complement orientation with a_snp + ins
        std::fs::write(
            &gfa_in,
            "H\tVN:Z:1.1\n\
             S\ta\tACGTACGTAC\n\
             S\tm1\tTTTTT\n\
             S\tr_snp\tA\n\
             S\ta_snp\tG\n\
             S\tm2\tCCCCG\n\
             S\tins\tTT\n\
             S\tm3\tAAAAA\n\
             S\tz\tGTGTGTGTGT\n\
             L\ta\t+\tm1\t+\t0M\n\
             L\ta\t+\tz\t+\t0M\n\
             L\tm1\t+\tr_snp\t+\t0M\n\
             L\tm1\t+\ta_snp\t+\t0M\n\
             L\tr_snp\t+\tm2\t+\t0M\n\
             L\ta_snp\t+\tm2\t+\t0M\n\
             L\tm2\t+\tm3\t+\t0M\n\
             L\tm2\t+\tins\t+\t0M\n\
             L\tins\t+\tm3\t+\t0M\n\
             L\tm3\t+\tz\t+\t0M\n\
             W\tREF\t0\tchr1\t0\t36\t>a>m1>r_snp>m2>m3>z\n\
             W\tREF2\t0\tchr1\t0\t36\t>a>m1>r_snp>m2>m3>z\n\
             W\tALT\t0\tchr1\t0\t38\t>a>m1>a_snp>m2>ins>m3>z\n\
             W\tREV\t0\tchr1\t0\t38\t<z<m3<ins<m2<a_snp<m1<a\n\
             W\tDEL\t0\tchr1\t0\t20\t>a>z\n",
        )?;
        let cli = CliArgs {
            gfa_in: gfa_in.to_string_lossy().to_string(),
            out_path: None,
            vcf_path: Some("out.vcf".to_owned()),
            vcf_ref: Some("REF".to_owned()),
            nested_snarls: true,
            collapse_dups: false,
            tag_walks: false,
            max_bubble_steps: 16,
            min_anchor_bp: 100,
            divergence: None,
            threads: None,
        };
        let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
        let mut gfa_buf = Vec::new();
        let mut cap = resolve::VcfCapture::with_segments(graph.num_segs());
        let _ = pass2_rewrite(&cli.gfa_in, &mut gfa_buf, &graph, &cli, Some(&mut cap))?;

        // Default (LV=0 only) collapses the whole spanning deletion into 1 top-level snarl.
        let mut vcf_lv0 = Vec::new();
        let st0 = export_snarl_vcf(
            &graph,
            &cap.seg_seqs,
            &cap.walks,
            Some("REF"),
            false,
            64,
            &mut vcf_lv0,
        )?;
        assert_eq!(st0.sites_emitted, 1);
        assert_eq!(st0.nested_sites, 0);

        // With --nested-snarls, both LV=0 deletion AND LV=1 internal SNP + insertion emit!
        let mut vcf_nested = Vec::new();
        let st1 = export_snarl_vcf(
            &graph,
            &cap.seg_seqs,
            &cap.walks,
            Some("REF"),
            true,
            64,
            &mut vcf_nested,
        )?;
        let _ = std::fs::remove_dir_all(&dir);
        let text = String::from_utf8(vcf_nested)?;
        assert_eq!(st1.sites_emitted, 3, "expected 1 LV=0 + 2 LV=1:\n{text}");
        assert_eq!(st1.nested_sites, 2);
        let rows: Vec<Vec<&str>> = text
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').collect())
            .collect();
        // Samples sorted: ALT(9), DEL(10), REF(11), REF2(12), REV(13)
        // Row 0 (POS 10, LV=0): spanning deletion anchored at last base 'C' of `a`
        assert_eq!(rows[0][1], "10");
        assert!(rows[0][7].ends_with(";LV=0"), "row0: {:?}", rows[0]);
        // Row 1 (POS 16, LV=1): internal SNP A -> G (ALT=1, DEL=., REF=0, REF2=0, REV=1)
        assert_eq!(rows[1][1], "16");
        assert_eq!(rows[1][3], "A");
        assert_eq!(rows[1][4], "G");
        assert!(
            rows[1][7].ends_with(";LV=1;PS=snarl_1"),
            "row1: {:?}",
            rows[1]
        );
        assert_eq!(rows[1][9], "1"); // ALT
        assert_eq!(rows[1][10], "."); // DEL (spans deletion -> missing)
        assert_eq!(rows[1][11], "0"); // REF
        assert_eq!(rows[1][12], "0"); // REF2
        assert_eq!(rows[1][13], "1"); // REV (reverse-complement walk)
        // Row 2 (POS 21, LV=1): internal insertion G -> GTT (ALT=1, DEL=., REF=0, REF2=0, REV=1)
        assert_eq!(rows[2][1], "21");
        assert_eq!(rows[2][3], "G");
        assert_eq!(rows[2][4], "GTT");
        assert!(
            rows[2][7].ends_with(";LV=1;PS=snarl_1"),
            "row2: {:?}",
            rows[2]
        );
        assert_eq!(rows[2][9], "1");
        assert_eq!(rows[2][10], ".");
        assert_eq!(rows[2][11], "0");
        assert_eq!(rows[2][12], "0");
        assert_eq!(rows[2][13], "1");
        Ok(())
    }
}
