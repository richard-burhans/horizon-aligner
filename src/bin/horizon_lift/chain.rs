//! UCSC `.chain[.gz]` coordinate liftover for `horizon_lift --chain`.
//!
//! Maps 0-based half-open target BED/GFF3 features onto forward (`+`) query coordinates via
//! ungapped chain blocks, handling reverse-strand (`-`) coordinate conversion, strand flipping,
//! `--min-frac`, `--bridge-bp` coalescing, `--to-chrom` filtering, and `--unmapped` reporting.

use crate::features::InputFeature;
use crate::lift::{LiftSummary, flip_strand};
use anyhow::{Context, Result};
use horizon_aligner::{UcscChain, parse_ucsc_chains};
use std::collections::{HashMap, HashSet};
use std::io::{Read as _, Write};
use std::path::Path;

pub(crate) struct ChainLiftConfig<'a> {
    pub chain_path: &'a Path,
    pub to_chroms: &'a HashSet<String>,
    pub features: &'a [InputFeature],
    pub min_frac: f64,
    pub bridge_bp: usize,
}

#[derive(Debug, Clone)]
struct MappedPiece {
    chrom: String,
    start: usize,
    end: usize,
    strand: char,
    mapped_bp: usize,
}

fn lift_feature_on_chain(
    feat: &InputFeature,
    ch: &UcscChain,
    bridge_bp: usize,
) -> (Vec<MappedPiece>, usize) {
    let mut raw_spans: Vec<(usize, usize, usize)> = Vec::new();
    let mut total_mapped = 0usize;

    for blk in &ch.blocks {
        let b_end = blk.t_start + blk.size;
        if b_end <= feat.start {
            continue;
        }
        if blk.t_start >= feat.end {
            break;
        }
        let s = blk.t_start.max(feat.start);
        let e = b_end.min(feat.end);
        if e <= s {
            continue;
        }
        let len = e - s;
        let off = s - blk.t_start;
        let q_strand_s = blk.q_start + off;
        let q_strand_e = q_strand_s + len;
        let (q_lo, q_hi) = ch.strand_to_fwd(q_strand_s, q_strand_e);
        raw_spans.push((q_lo, q_hi, len));
        total_mapped += len;
    }

    if raw_spans.is_empty() {
        return (Vec::new(), 0);
    }
    raw_spans.sort_unstable_by_key(|&(s, e, _)| (s, e));
    let out_strand = if ch.q_strand == '-' {
        flip_strand(feat.strand)
    } else {
        feat.strand
    };

    let mut pieces: Vec<MappedPiece> = Vec::with_capacity(raw_spans.len());
    for (s, e, bp) in raw_spans {
        if let Some(last) = pieces.last_mut()
            && s <= last.end.saturating_add(bridge_bp)
        {
            last.end = last.end.max(e);
            last.mapped_bp += bp;
            continue;
        }
        pieces.push(MappedPiece {
            chrom: ch.q_name.clone(),
            start: s,
            end: e,
            strand: out_strand,
            mapped_bp: bp,
        });
    }
    (pieces, total_mapped)
}

fn write_lifted_pieces<W: Write>(
    feat: &InputFeature,
    total_bp: usize,
    mapped_bp: usize,
    pieces: &[MappedPiece],
    summary: &mut LiftSummary,
    out: &mut W,
) -> std::io::Result<()> {
    #[allow(clippy::cast_precision_loss)]
    let frac = mapped_bp as f64 / total_bp as f64;
    if pieces.len() == 1 {
        if mapped_bp == total_bp {
            summary.lifted_exact += 1;
        } else {
            summary.lifted_partial += 1;
        }
    } else {
        summary.lifted_split += 1;
    }
    for (idx, p) in pieces.iter().enumerate() {
        let status = if pieces.len() > 1 {
            format!("SPLIT_{}/{}", idx + 1, pieces.len())
        } else if mapped_bp == total_bp {
            "EXACT".to_string()
        } else {
            "PARTIAL".to_string()
        };
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}:{}-{}\t{:.4}\t{}",
            p.chrom,
            p.start,
            p.end,
            feat.name,
            feat.score,
            p.strand,
            feat.chrom,
            feat.start,
            feat.end,
            frac,
            status,
        )?;
    }
    Ok(())
}

pub(crate) fn run_chain_lift<W: Write, U: Write>(
    cfg: &ChainLiftConfig<'_>,
    out: &mut W,
    mut unmapped_out: Option<&mut U>,
) -> Result<LiftSummary> {
    let mut reader = horizon_aligner::io_util::open_text_reader(cfg.chain_path)
        .with_context(|| format!("opening UCSC chain '{}'", cfg.chain_path.display()))?;
    let mut text = String::new();
    reader.read_to_string(&mut text)?;
    let chains = parse_ucsc_chains(&text)?;

    let mut by_target: HashMap<&str, Vec<&UcscChain>> = HashMap::new();
    for ch in &chains {
        if !cfg.to_chroms.is_empty() && !cfg.to_chroms.contains(&ch.q_name) {
            continue;
        }
        by_target.entry(ch.t_name.as_str()).or_default().push(ch);
    }
    for list in by_target.values_mut() {
        list.sort_unstable_by(|a, b| b.score.cmp(&a.score).then_with(|| a.id.cmp(&b.id)));
    }

    let mut summary = LiftSummary::default();
    for feat in cfg.features {
        let total_bp = feat.end.saturating_sub(feat.start);
        if total_bp == 0 {
            summary.unmapped += 1;
            continue;
        }
        let mut best: Option<(usize, i64, Vec<MappedPiece>)> = None;
        let mut best_failed_frac = 0.0_f64;

        if let Some(cands) = by_target.get(feat.chrom.as_str()) {
            for &ch in cands {
                if ch.t_end <= feat.start || ch.t_start >= feat.end {
                    continue;
                }
                let (pieces, mapped_bp) = lift_feature_on_chain(feat, ch, cfg.bridge_bp);
                if pieces.is_empty() {
                    continue;
                }
                #[allow(clippy::cast_precision_loss)]
                let frac = mapped_bp as f64 / total_bp as f64;
                if frac < cfg.min_frac {
                    best_failed_frac = best_failed_frac.max(frac);
                    continue;
                }
                let replace = match best {
                    None => true,
                    Some((prev_bp, prev_score, _)) => {
                        mapped_bp > prev_bp || (mapped_bp == prev_bp && ch.score > prev_score)
                    }
                };
                if replace {
                    best = Some((mapped_bp, ch.score, pieces));
                }
            }
        }

        if let Some((mapped_bp, _score, pieces)) = best {
            write_lifted_pieces(feat, total_bp, mapped_bp, &pieces, &mut summary, out)?;
        } else {
            summary.unmapped += 1;
            if let Some(u) = unmapped_out.as_deref_mut() {
                writeln!(
                    u,
                    "{}\t{}\t{}\t{}\t{}\t{}\tUNMAPPED({:.3})",
                    feat.chrom,
                    feat.start,
                    feat.end,
                    feat.name,
                    feat.score,
                    feat.strand,
                    best_failed_frac
                )?;
            }
        }
    }

    Ok(summary)
}

#[cfg(test)]
mod tests {
    use super::*;
    use horizon_aligner::{UcscChainBlock, format_ucsc_chains};

    #[test]
    fn chain_lift_handles_forward_reverse_partial_and_unmapped() {
        // Chain 1 (+): chrT[10..60) -> chrQ[100..150) (block 10..30 -> 100..120, gap 10t/5q, block 40..60 -> 125..145)
        // Chain 2 (-): chrT[100..140) -> chrQ fwd [200..240) on qSize=300 => '-' coords [60..100)
        let chains = vec![
            UcscChain {
                id: 1,
                score: 4000,
                t_name: "chrT".into(),
                t_size: 200,
                t_start: 10,
                t_end: 60,
                q_name: "chrQ".into(),
                q_size: 300,
                q_strand: '+',
                q_start: 100,
                q_end: 145,
                blocks: vec![
                    UcscChainBlock {
                        t_start: 10,
                        q_start: 100,
                        size: 20,
                    },
                    UcscChainBlock {
                        t_start: 40,
                        q_start: 125,
                        size: 20,
                    },
                ],
            },
            UcscChain {
                id: 2,
                score: 3500,
                t_name: "chrT".into(),
                t_size: 200,
                t_start: 100,
                t_end: 140,
                q_name: "chrQ".into(),
                q_size: 300,
                q_strand: '-',
                q_start: 60,
                q_end: 100,
                blocks: vec![UcscChainBlock {
                    t_start: 100,
                    q_start: 60,
                    size: 40,
                }],
            },
        ];
        let dir = std::env::temp_dir().join(format!("horizon_lift_chain_{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("mkdir");
        let ch_path = dir.join("test.chain");
        std::fs::write(&ch_path, format_ucsc_chains(&chains)).expect("write chain");

        let mk = |s: usize, e: usize, name: &str| InputFeature {
            chrom: "chrT".into(),
            start: s,
            end: e,
            name: name.into(),
            score: "0".into(),
            strand: '+',
        };
        let features = vec![
            mk(12, 22, "fwd_exact"),
            mk(105, 115, "inv_exact"),
            mk(20, 50, "bridged_partial"),
            mk(160, 180, "unmapped_gap"),
        ];
        let cfg = ChainLiftConfig {
            chain_path: &ch_path,
            to_chroms: &HashSet::new(),
            features: &features,
            min_frac: 0.5,
            bridge_bp: 10,
        };
        let mut out = Vec::new();
        let mut unm = Vec::new();
        let summary = run_chain_lift(&cfg, &mut out, Some(&mut unm)).expect("run chain lift");
        assert_eq!(summary.lifted_exact, 2);
        assert_eq!(summary.lifted_partial, 1);
        assert_eq!(summary.unmapped, 1);

        let out_s = String::from_utf8(out).expect("utf8");
        assert!(out_s.contains("chrQ\t102\t112\tfwd_exact\t0\t+\tchrT:12-22\t1.0000\tEXACT\n"));
        // t[105..115) => q_rc[65..75) => q_fwd[300-75..300-65) = [225..235) '-'
        assert!(out_s.contains("chrQ\t225\t235\tinv_exact\t0\t-\tchrT:105-115\t1.0000\tEXACT\n"));
        // t[20..50) overlaps [20..30)->[110..120) and [40..50)->[125..135), bridged (gap 5 <= 10) => [110..135)
        assert!(
            out_s.contains("chrQ\t110\t135\tbridged_partial\t0\t+\tchrT:20-50\t0.6667\tPARTIAL\n")
        );
        let unm_s = String::from_utf8(unm).expect("utf8");
        assert!(unm_s.contains("chrT\t160\t180\tunmapped_gap\t0\t+\tUNMAPPED(0.000)\n"));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
