//! Multi-species MAF + GFA merger with e-graph duplication resolution.

use anyhow::{Context, Result, bail};
use horizon_aligner::merge_pairwise_mafs;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.iter().any(|a| a == "-h" || a == "--help") {
        println!(
            "usage: horizon_merge <ref_species> <bin_bp> <maf_out> <gfa_out> <sp=path.maf>..."
        );
        return Ok(());
    }
    if args.len() < 5 {
        bail!("usage: horizon_merge <ref_species> <bin_bp> <maf_out> <gfa_out> <sp=path.maf>...");
    }

    let ref_species = &args[0];
    let bin_bp: usize = args[1].parse().context("parsing bin_bp")?;
    let maf_out_path = PathBuf::from(&args[2]);
    let gfa_out_path = PathBuf::from(&args[3]);

    let mut inputs: Vec<(String, String)> = Vec::new();
    for spec in &args[4..] {
        let Some((sp, path)) = spec.split_once('=') else {
            bail!("expected <species>=<maf_path>, got {spec}");
        };
        let text = fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
        inputs.push((sp.to_owned(), text));
    }

    let t0 = Instant::now();
    let (maf_out, gfa_out, total_blocks, multi_blocks) =
        merge_pairwise_mafs(ref_species, &inputs, bin_bp);
    fs::write(&maf_out_path, maf_out)
        .with_context(|| format!("writing {}", maf_out_path.display()))?;
    fs::write(&gfa_out_path, gfa_out)
        .with_context(|| format!("writing {}", gfa_out_path.display()))?;

    println!(
        "ref={} species={} total_blocks={} multi_species_blocks={} elapsed_s={:.4}",
        ref_species,
        inputs.len() + 1,
        total_blocks,
        multi_blocks,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

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
