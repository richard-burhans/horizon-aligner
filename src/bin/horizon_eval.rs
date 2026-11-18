//! Evaluate `horizon-aligner` recall and abstention against Galaxy LASTZ `rbest_maf` on Pf vs Pv Chr1.

use anyhow::{Context, Result};
use horizon_aligner::{align_contig_pair_with_config, AlignConfig, prepare_contig_index_opts};
use std::collections::HashSet;
use std::fs;
use std::path::Path;
use std::time::Instant;

fn read_named_contig(path: &Path, target: &str) -> Result<Vec<u8>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut seq = Vec::new();
    let mut active = false;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if active {
                break;
            }
            let id = rest.split_whitespace().next().unwrap_or("");
            if id == target {
                active = true;
            }
        } else if active {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    Ok(seq)
}

#[allow(clippy::similar_names)]
fn parse_maf_pairs(path: &Path, cx: &str, cy: &str) -> Result<HashSet<(usize, usize)>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut out = HashSet::new();
    let mut cur: Vec<Vec<String>> = Vec::new();
    for line in text.lines() {
        if line.starts_with("s ") {
            cur.push(line.split_whitespace().map(str::to_owned).collect());
        } else if line.trim().is_empty() && !cur.is_empty() {
            if cur.len() == 2 && cur[0][1].contains(cx) && cur[1][1].contains(cy) {
                let sx: usize = cur[0][2].parse()?;
                let stx = &cur[0][4];
                let src_lx: usize = cur[0][5].parse()?;
                let ax = cur[0][6].as_bytes();

                let sy: usize = cur[1][2].parse()?;
                let sty = &cur[1][4];
                let src_ly: usize = cur[1][5].parse()?;
                let ay = cur[1][6].as_bytes();

                let (mut ix, mut iy) = (0usize, 0usize);
                for (&ca, &cb) in ax.iter().zip(ay) {
                    if ca != b'-' && cb != b'-' {
                        let px = if stx == "+" {
                            sx + ix
                        } else {
                            src_lx - 1 - (sx + ix)
                        };
                        let py = if sty == "+" {
                            sy + iy
                        } else {
                            src_ly - 1 - (sy + iy)
                        };
                        out.insert((px, py));
                    }
                    if ca != b'-' {
                        ix += 1;
                    }
                    if cb != b'-' {
                        iy += 1;
                    }
                }
            }
            cur.clear();
        }
    }
    Ok(out)
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 6 {
        anyhow::bail!("usage: horizon_eval <ref_name> <target_name> <ref.fa> <target.fa> <truth.maf>");
    }
    let ref_name = &args[1];
    let target_name = &args[2];
    let pf = read_named_contig(Path::new(&args[3]), ref_name)?;
    let pv = read_named_contig(Path::new(&args[4]), target_name)?;
    let maf_all = parse_maf_pairs(Path::new(&args[5]), ref_name, target_name)?;


    let mut min_fallback_cov = None;
    let mut hsp_threshold = None;
    let mut i = 6;
    while i < args.len() {
        if args[i] == "--min-fallback-cov" && i + 1 < args.len() {
            min_fallback_cov = Some(args[i + 1].parse()?);
            i += 2;
        } else if args[i] == "--hsp-threshold" && i + 1 < args.len() {
            hsp_threshold = Some(args[i + 1].parse()?);
            i += 2;
        } else {
            i += 1;
        }
    }

    let t0 = Instant::now();
    let idx_pf = prepare_contig_index_opts(&pf, false)?;
    let idx_pv = prepare_contig_index_opts(&pv, false)?;
    
    let cfg = AlignConfig {
        min_fallback_cov,
        hsp_threshold,
        ..AlignConfig::default()
    };
    let aln = align_contig_pair_with_config(&idx_pf, &idx_pv, &cfg)?;

    let elapsed = t0.elapsed().as_secs_f64();

    let proto_set: HashSet<(usize, usize)> = aln.aligned_pairs.iter().copied().collect();
    let overlap = proto_set.intersection(&maf_all).count();
    #[allow(clippy::cast_precision_loss)]
    let recall = (overlap as f64) / (maf_all.len().max(1) as f64);
    println!(
        "elapsed={elapsed:.3}s proto_pairs={} maf_pairs={} overlap={overlap} raw_maf_recall={:.1}% bubbles={}",
        proto_set.len(),
        maf_all.len(),
        recall * 100.0,
        aln.macro_bubbles + aln.internal_bubbles
    );
    Ok(())
}
