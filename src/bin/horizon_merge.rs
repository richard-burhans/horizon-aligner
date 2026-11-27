//! Multi-species MAF + GFA merger with e-graph duplication resolution.

use anyhow::{Context, Result, bail};
use horizon_aligner::io_util::open_text_writer;
use horizon_aligner::{
    init_rayon_threads, merge_pairwise_maf_paths_to_writers, parse_threads_arg,
    validate_species_name,
};
use std::path::PathBuf;
use std::time::Instant;

fn main() -> Result<()> {
    let raw_args: Vec<String> = std::env::args().skip(1).collect();
    let mut threads: Option<usize> = None;
    let mut args: Vec<String> = Vec::with_capacity(raw_args.len());
    let mut idx = 0usize;
    while idx < raw_args.len() {
        match raw_args[idx].as_str() {
            "-h" | "--help" => {
                println!(
                    "usage: horizon_merge [-t N | --threads N] <ref_species> <bin_bp> <maf_out> <gfa_out> <sp=path.maf>..."
                );
                return Ok(());
            }
            "-V" | "--version" => {
                println!("horizon_merge {}", env!("CARGO_PKG_VERSION"));
                return Ok(());
            }
            "-t" | "--threads" => {
                let raw = raw_args
                    .get(idx + 1)
                    .context("missing value for --threads")?;
                threads = Some(parse_threads_arg(raw)?);
                idx += 2;
            }
            "--" => {
                args.extend(raw_args[(idx + 1)..].iter().cloned());
                break;
            }
            other if other.starts_with('-') && other != "-" => bail!("unknown flag `{other}`"),
            other => {
                args.push(other.to_owned());
                idx += 1;
            }
        }
    }
    if args.len() < 5 {
        bail!(
            "usage: horizon_merge [-t N | --threads N] <ref_species> <bin_bp> <maf_out> <gfa_out> <sp=path.maf>..."
        );
    }
    init_rayon_threads(threads)?;

    let ref_species = &args[0];
    validate_species_name(ref_species)?;
    let bin_bp: usize = args[1].parse().context("parsing bin_bp")?;
    anyhow::ensure!(bin_bp > 0, "bin_bp must be > 0");
    let maf_out_path = &args[2];
    let gfa_out_path = &args[3];
    anyhow::ensure!(
        !(maf_out_path == "-" && gfa_out_path == "-"),
        "cannot write both <maf_out> and <gfa_out> to stdout ('-')"
    );
    let writes_stdout = maf_out_path == "-" || gfa_out_path == "-";

    let mut inputs: Vec<(String, PathBuf)> = Vec::new();
    let mut seen_species = std::collections::HashSet::new();
    for spec in &args[4..] {
        let Some((sp, path)) = spec.split_once('=') else {
            bail!("expected <species>=<maf_path>, got {spec}");
        };
        validate_species_name(sp)?;
        anyhow::ensure!(
            !path.is_empty(),
            "expected non-empty <species>=<maf_path>, got {spec}"
        );
        anyhow::ensure!(
            sp != ref_species,
            "query species `{sp}` cannot match reference species `{ref_species}`"
        );
        anyhow::ensure!(
            seen_species.insert(sp.to_owned()),
            "duplicate species name: {sp}"
        );
        inputs.push((sp.to_owned(), PathBuf::from(path)));
    }

    let t0 = Instant::now();
    let mut maf_writer =
        open_text_writer(Some(maf_out_path)).with_context(|| format!("creating {maf_out_path}"))?;
    let mut gfa_writer =
        open_text_writer(Some(gfa_out_path)).with_context(|| format!("creating {gfa_out_path}"))?;
    let (total_blocks, multi_blocks) = merge_pairwise_maf_paths_to_writers(
        ref_species,
        &inputs,
        bin_bp,
        &mut maf_writer,
        &mut gfa_writer,
    )?;
    maf_writer.finish()?;
    gfa_writer.finish()?;

    let summary = format!(
        "ref={} species={} total_blocks={} multi_species_blocks={} elapsed_s={:.4}",
        ref_species,
        inputs.len() + 1,
        total_blocks,
        multi_blocks,
        t0.elapsed().as_secs_f64()
    );
    if writes_stdout {
        eprintln!("{summary}");
    } else {
        println!("{summary}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use horizon_aligner::merge_pairwise_mafs;

    #[test]
    fn merge_pairwise_mafs_merges_two_species() {
        let maf_a = "a score=10\ns hg38.chr1 0 10 + 100 ACGTACGTAC\ns panTro6.chr1 0 10 + 100 ACGTACGTAC\n\n";
        let maf_b =
            "a score=10\ns hg38.chr1 0 10 + 100 ACGTACGTAC\ns mm39.chr1 0 10 + 100 ACGTACGTAC\n\n";
        let inputs = vec![
            ("panTro6".to_string(), maf_a.to_string()),
            ("mm39".to_string(), maf_b.to_string()),
        ];
        let (maf_out, gfa_out, total_blocks, multi_blocks) =
            merge_pairwise_mafs("hg38", &inputs, 100);
        assert!(total_blocks >= 1);
        assert!(multi_blocks >= 1);
        assert!(maf_out.contains("hg38.chr1"));
        assert!(gfa_out.contains("H\tVN:Z:"));
    }
}
