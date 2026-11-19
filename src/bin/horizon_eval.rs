//! Evaluate `horizon-aligner` recall, precision, and abstention against truth MAF (supports single- and multi-contig FASTAs).

use anyhow::{Context, Result};
use horizon_aligner::{align_contig_pair_with_config, prepare_contig_index_opts, AlignConfig};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;
use std::time::Instant;

fn read_named_contigs(path: &Path, target: &str) -> Result<HashMap<String, Vec<u8>>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut contigs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut active_id: Option<String> = None;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            let id = rest.split_whitespace().next().unwrap_or("");
            if id == target || id.starts_with(&format!("{}.", target)) {
                active_id = Some(id.to_string());
                contigs.entry(id.to_string()).or_default();
            } else {
                active_id = None;
            }
        } else if let Some(ref id) = active_id {
            if let Some(seq) = contigs.get_mut(id) {
                seq.extend_from_slice(line.trim().as_bytes());
            }
        }
    }
    Ok(contigs)
}

#[allow(clippy::similar_names)]
fn parse_maf_pairs_by_contig(
    path: &Path,
    cx: &str,
    cy: &str,
) -> Result<HashMap<(String, String), HashSet<(usize, usize)>>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut out: HashMap<(String, String), HashSet<(usize, usize)>> = HashMap::new();
    let mut cur: Vec<Vec<String>> = Vec::new();
    for line in text.lines() {
        if line.starts_with("s ") {
            cur.push(line.split_whitespace().map(str::to_owned).collect());
        } else if line.trim().is_empty() && !cur.is_empty() {
            if cur.len() == 2 && cur[0][1].contains(cx) && cur[1][1].contains(cy) {
                let ref_id = cur[0][1].clone();
                let tgt_id = cur[1][1].clone();
                let sx: usize = cur[0][2].parse()?;
                let stx = &cur[0][4];
                let src_lx: usize = cur[0][5].parse()?;
                let ax = cur[0][6].as_bytes();

                let sy: usize = cur[1][2].parse()?;
                let sty = &cur[1][4];
                let src_ly: usize = cur[1][5].parse()?;
                let ay = cur[1][6].as_bytes();

                let pair_set = out.entry((ref_id, tgt_id)).or_default();
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
                        pair_set.insert((px, py));
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
        anyhow::bail!(
            "usage: horizon_eval <ref_name> <target_name> <ref.fa> <target.fa> <truth.maf>"
        );
    }
    let ref_name = &args[1];
    let target_name = &args[2];
    let ref_contigs = read_named_contigs(Path::new(&args[3]), ref_name)?;
    let tgt_contigs = read_named_contigs(Path::new(&args[4]), target_name)?;
    let maf_by_pair = parse_maf_pairs_by_contig(Path::new(&args[5]), ref_name, target_name)?;

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

    let cfg = AlignConfig {
        min_fallback_cov,
        hsp_threshold,
        ..AlignConfig::default()
    };

    let t0 = Instant::now();
    let mut ref_idx_cache = HashMap::new();
    let mut tgt_idx_cache = HashMap::new();
    for (id, seq) in &ref_contigs {
        ref_idx_cache.insert(id.clone(), prepare_contig_index_opts(seq, false)?);
    }
    for (id, seq) in &tgt_contigs {
        tgt_idx_cache.insert(id.clone(), prepare_contig_index_opts(seq, false)?);
    }

    let mut total_proto = 0usize;
    let mut total_maf = 0usize;
    let mut total_overlap = 0usize;
    let mut total_bubbles = 0usize;

    let mut pairs: Vec<((String, String), HashSet<(usize, usize)>)> =
        maf_by_pair.into_iter().collect();
    pairs.sort_by(|a, b| a.0.cmp(&b.0));

    for ((ref_id, tgt_id), maf_set) in &pairs {
        total_maf += maf_set.len();
        let (Some(idx_r), Some(idx_t)) = (ref_idx_cache.get(ref_id), tgt_idx_cache.get(tgt_id))
        else {
            continue;
        };
        let aln = align_contig_pair_with_config(idx_r, idx_t, &cfg)?;
        let proto_set: HashSet<(usize, usize)> = aln.aligned_pairs.iter().copied().collect();
        total_proto += proto_set.len();
        total_overlap += proto_set.intersection(maf_set).count();
        total_bubbles += aln.macro_bubbles + aln.internal_bubbles;
    }

    let elapsed = t0.elapsed().as_secs_f64();
    let recall = if total_maf == 0 {
        0.0
    } else {
        100.0 * total_overlap as f64 / total_maf as f64
    };
    let precision = if total_proto == 0 {
        0.0
    } else {
        100.0 * total_overlap as f64 / total_proto as f64
    };
    println!(
        "elapsed={:.3}s proto_pairs={} maf_pairs={} overlap={} raw_maf_recall={:.1}% raw_maf_precision={:.1}% bubbles={}",
        elapsed, total_proto, total_maf, total_overlap, recall, precision, total_bubbles
    );
    Ok(())
}
