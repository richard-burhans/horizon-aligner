//! Graph-aware BED/GFF3 coordinate liftover across GFA v1.0 (`P`) and v1.1 (`W`) walks.

use anyhow::{Context, Result, bail};
use horizon_aligner::coord::Interval;
use horizon_aligner::gfa::{parse_path_steps, parse_walk_steps, stream_gfa};
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

#[derive(Debug, Clone)]
struct InputFeature {
    chrom: String,
    start: usize, // 0-based half-open
    end: usize,
    name: String,
    score: String,
    strand: char,
}

#[derive(Debug, Clone)]
struct PathRecord {
    display_chrom: String,
    match_keys: Vec<String>,
    start_offset: usize,
    steps: Vec<(String, bool)>,
}

#[derive(Debug, Clone)]
struct StepSpan {
    seg: String,
    is_fwd: bool,
    path_start: usize,
    path_end: usize,
}

#[derive(Debug, Clone)]
struct TargetOcc {
    target_chrom: String,
    target_start: usize,
    target_end: usize,
    is_fwd: bool,
    path_idx: usize,
    step_idx: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LiftedPiece {
    chrom: String,
    start: usize,
    end: usize,
    strand: char,
    src_start: usize,
    src_end: usize,
    duplicated: bool,
}

fn read_vm_hwm_kb() -> Option<usize> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            let num = rest.split_whitespace().next()?.parse::<usize>().ok()?;
            return Some(num);
        }
    }
    None
}

fn parse_features(path: &Path, is_gff: bool) -> Result<Vec<InputFeature>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut out = Vec::new();
    for (idx, line_res) in reader.lines().enumerate() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with("track") {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if is_gff {
            if cols.len() < 9 {
                continue;
            }
            let chrom = cols[0].to_string();
            let start_1: usize = cols[3]
                .parse()
                .with_context(|| format!("invalid GFF3 start at line {}", idx + 1))?;
            let end_1: usize = cols[4]
                .parse()
                .with_context(|| format!("invalid GFF3 end at line {}", idx + 1))?;
            let start = start_1.saturating_sub(1);
            let end = end_1;
            let score = cols[5].to_string();
            let strand = cols[6].chars().next().unwrap_or('+');
            let strand = if strand == '-' { '-' } else { '+' };
            let name = cols[8]
                .split(';')
                .find_map(|kv| kv.strip_prefix("ID=").or_else(|| kv.strip_prefix("Name=")))
                .unwrap_or(cols[2])
                .to_string();
            out.push(InputFeature {
                chrom,
                start,
                end,
                name,
                score,
                strand,
            });
        } else {
            if cols.len() < 3 {
                continue;
            }
            let chrom = cols[0].to_string();
            let start: usize = cols[1]
                .parse()
                .with_context(|| format!("invalid BED start at line {}", idx + 1))?;
            let end: usize = cols[2]
                .parse()
                .with_context(|| format!("invalid BED end at line {}", idx + 1))?;
            let name = cols.get(3).copied().unwrap_or(".").to_string();
            let score = cols.get(4).copied().unwrap_or("0").to_string();
            let strand = cols
                .get(5)
                .and_then(|s| s.chars().next())
                .map_or('+', |c| if c == '-' { '-' } else { '+' });
            out.push(InputFeature {
                chrom,
                start,
                end,
                name,
                score,
                strand,
            });
        }
    }
    Ok(out)
}

fn matches_sample(name: &str, prefix: &str) -> bool {
    name == prefix
        || name.starts_with(&format!("{prefix}."))
        || name.starts_with(&format!("{prefix}#"))
}

fn strip_sample_prefix<'a>(name: &'a str, prefix: &str) -> &'a str {
    if let Some(rest) = name.strip_prefix(&format!("{prefix}.")) {
        rest
    } else if let Some(rest) = name.strip_prefix(&format!("{prefix}#")) {
        rest
    } else {
        name
    }
}

fn flip_strand(strand: char) -> char {
    if strand == '-' { '+' } else { '-' }
}

fn context_jaccard(
    src_steps: &[StepSpan],
    src_pos: usize,
    tgt_steps: &[String],
    tgt_pos: usize,
) -> f64 {
    const CTX_RADIUS: usize = 4;
    let s_lo = src_pos.saturating_sub(CTX_RADIUS);
    let s_hi = (src_pos + CTX_RADIUS + 1).min(src_steps.len());
    let t_lo = tgt_pos.saturating_sub(CTX_RADIUS);
    let t_hi = (tgt_pos + CTX_RADIUS + 1).min(tgt_steps.len());
    let mut src_counts: HashMap<&str, usize> = HashMap::with_capacity(s_hi - s_lo);
    for st in &src_steps[s_lo..s_hi] {
        *src_counts.entry(st.seg.as_str()).or_insert(0) += 1;
    }
    let mut tgt_counts: HashMap<&str, usize> = HashMap::with_capacity(t_hi - t_lo);
    for seg in &tgt_steps[t_lo..t_hi] {
        *tgt_counts.entry(seg.as_str()).or_insert(0) += 1;
    }
    let mut inter = 0usize;
    let mut union = 0usize;
    for (k, &sc) in &src_counts {
        let tc = tgt_counts.get(k).copied().unwrap_or(0);
        inter += sc.min(tc);
        union += sc.max(tc);
    }
    for (k, &tc) in &tgt_counts {
        if !src_counts.contains_key(k) {
            union += tc;
        }
    }
    if union == 0 {
        0.0
    } else {
        #[allow(clippy::cast_precision_loss)]
        {
            inter as f64 / union as f64
        }
    }
}

/// Project a single source interval `[src_s, src_e)` on `src_steps` onto `target_index`.
#[allow(clippy::too_many_lines)]
fn lift_feature(
    feat: &InputFeature,
    src_steps: &[StepSpan],
    seg_lens: &HashMap<String, usize>,
    target_index: &HashMap<String, Vec<TargetOcc>>,
    target_step_names: &[Vec<String>],
    bridge_bp: usize,
) -> Vec<LiftedPiece> {
    if feat.start >= feat.end || src_steps.is_empty() {
        return Vec::new();
    }

    // Binary search for first StepSpan with path_end > feat.start
    let first_idx = src_steps.partition_point(|s| s.path_end <= feat.start);
    let mut raw_pieces: Vec<LiftedPiece> = Vec::new();
    let mut last_target: Option<(String, usize, usize, char)> = None;

    for (rel_idx, step) in src_steps[first_idx..].iter().enumerate() {
        let src_pos = first_idx + rel_idx;
        if step.path_start >= feat.end {
            break;
        }
        let ov_s = feat.start.max(step.path_start);
        let ov_e = feat.end.min(step.path_end);
        if ov_s >= ov_e {
            continue;
        }
        let seg_len = match seg_lens.get(&step.seg) {
            Some(&l) if l > 0 => l,
            _ => continue,
        };
        let step_iv = Interval::new(
            ov_s - step.path_start,
            ov_e - step.path_start,
            if step.is_fwd { '+' } else { '-' },
            seg_len,
        );
        // Convert to forward segment-local coordinates [seg_s, seg_e)
        let (seg_s, seg_e) = if step.is_fwd {
            (step_iv.start, step_iv.end)
        } else {
            let rev_s = step_iv.ucsc_start();
            (rev_s, rev_s + step_iv.span())
        };

        let Some(occs) = target_index.get(&step.seg) else {
            continue;
        };
        if occs.is_empty() {
            continue;
        }
        let duplicated = occs.len() > 1;
        // Prefer highest flanking-context Jaccard similarity (Guarracino et al., 2022),
        // breaking ties by linear proximity to last_target.
        let chosen = if occs.len() == 1 {
            &occs[0]
        } else {
            occs.iter()
                .max_by(|a, b| {
                    let ja = context_jaccard(
                        src_steps,
                        src_pos,
                        &target_step_names[a.path_idx],
                        a.step_idx,
                    );
                    let jb = context_jaccard(
                        src_steps,
                        src_pos,
                        &target_step_names[b.path_idx],
                        b.step_idx,
                    );
                    let dist = |occ: &TargetOcc| -> usize {
                        if let Some((ref prev_chr, prev_s, prev_e, _)) = last_target {
                            if &occ.target_chrom != prev_chr {
                                usize::MAX
                            } else if occ.target_start >= prev_e {
                                occ.target_start - prev_e
                            } else {
                                prev_s.saturating_sub(occ.target_end)
                            }
                        } else {
                            0
                        }
                    };
                    ja.total_cmp(&jb).then_with(|| dist(b).cmp(&dist(a)))
                })
                .unwrap_or(&occs[0])
        };

        let (t_s, t_e) = if chosen.is_fwd {
            (chosen.target_start + seg_s, chosen.target_start + seg_e)
        } else {
            (
                chosen.target_start + (seg_len - seg_e),
                chosen.target_start + (seg_len - seg_s),
            )
        };
        let out_strand = if chosen.is_fwd == step.is_fwd {
            feat.strand
        } else {
            flip_strand(feat.strand)
        };
        let piece = LiftedPiece {
            chrom: chosen.target_chrom.clone(),
            start: t_s,
            end: t_e,
            strand: out_strand,
            src_start: ov_s,
            src_end: ov_e,
            duplicated,
        };
        last_target = Some((piece.chrom.clone(), piece.start, piece.end, piece.strand));

        // Coalesce contiguous (or within bridge_bp on the same chrom & strand) pieces.
        if let Some(prev) = raw_pieces.last_mut() {
            let same_ctx = prev.chrom == piece.chrom && prev.strand == piece.strand;
            let same_rel_dir = chosen.is_fwd == step.is_fwd;
            if same_ctx
                && same_rel_dir
                && piece.start >= prev.end
                && piece.start - prev.end <= bridge_bp
                && piece.src_start >= prev.src_end
                && piece.src_start - prev.src_end <= bridge_bp
            {
                prev.end = piece.end;
                prev.src_end = piece.src_end;
                prev.duplicated |= piece.duplicated;
                continue;
            }
            if same_ctx
                && !same_rel_dir
                && prev.start >= piece.end
                && prev.start - piece.end <= bridge_bp
                && piece.src_start >= prev.src_end
                && piece.src_start - prev.src_end <= bridge_bp
            {
                prev.start = piece.start;
                prev.src_end = piece.src_end;
                prev.duplicated |= piece.duplicated;
                continue;
            }
        }
        raw_pieces.push(piece);
    }

    raw_pieces
}

struct LiftConfig<'a> {
    gfa_path: &'a str,
    from_sample: &'a str,
    to_sample: &'a str,
    features: &'a [InputFeature],
    min_frac: f64,
    bridge_bp: usize,
}

struct LiftSummary {
    lifted_exact: usize,
    lifted_partial: usize,
    lifted_split: usize,
    unmapped: usize,
}

#[allow(clippy::too_many_lines)]
fn run_lift<W1: Write, W2: Write>(
    cfg: &LiftConfig<'_>,
    out: &mut W1,
    mut unmapped_out: Option<&mut W2>,
) -> Result<LiftSummary> {
    let wanted_chroms: HashSet<String> = cfg.features.iter().map(|f| f.chrom.clone()).collect();
    let source_paths = std::cell::RefCell::new(Vec::<PathRecord>::new());
    let target_paths = std::cell::RefCell::new(Vec::<PathRecord>::new());

    // Pass 1: collect matching source and target P/W records.
    stream_gfa(
        cfg.gfa_path,
        |_| Ok(()),
        |_| Ok(()),
        |p| {
            if matches_sample(&p.name, cfg.from_sample) {
                let short = strip_sample_prefix(&p.name, cfg.from_sample).to_string();
                if wanted_chroms.is_empty()
                    || wanted_chroms.contains(&p.name)
                    || wanted_chroms.contains(&short)
                {
                    let steps = parse_path_steps(&p.segment_names)
                        .into_iter()
                        .map(|(s, fwd)| (s.to_string(), fwd))
                        .collect();
                    source_paths.borrow_mut().push(PathRecord {
                        display_chrom: p.name.clone(),
                        match_keys: vec![p.name.clone(), short],
                        start_offset: 0,
                        steps,
                    });
                }
            }
            if matches_sample(&p.name, cfg.to_sample) {
                let short = strip_sample_prefix(&p.name, cfg.to_sample).to_string();
                if wanted_chroms.is_empty()
                    || wanted_chroms.contains(&p.name)
                    || wanted_chroms.contains(&short)
                {
                    let steps = parse_path_steps(&p.segment_names)
                        .into_iter()
                        .map(|(s, fwd)| (s.to_string(), fwd))
                        .collect();
                    target_paths.borrow_mut().push(PathRecord {
                        display_chrom: p.name.clone(),
                        match_keys: vec![p.name, short],
                        start_offset: 0,
                        steps,
                    });
                }
            }
            Ok(())
        },
        |w| {
            if matches_sample(&w.sample, cfg.from_sample) {
                let full = format!("{}.{}", w.sample, w.seq_id);
                if wanted_chroms.is_empty()
                    || wanted_chroms.contains(&w.seq_id)
                    || wanted_chroms.contains(&full)
                {
                    let steps = parse_walk_steps(&w.walk)
                        .into_iter()
                        .map(|(s, fwd)| (s.to_string(), fwd))
                        .collect();
                    source_paths.borrow_mut().push(PathRecord {
                        display_chrom: full.clone(),
                        match_keys: vec![w.seq_id.clone(), full],
                        start_offset: w.start,
                        steps,
                    });
                }
            }
            if matches_sample(&w.sample, cfg.to_sample) {
                let full = format!("{}.{}", w.sample, w.seq_id);
                if wanted_chroms.is_empty()
                    || wanted_chroms.contains(&w.seq_id)
                    || wanted_chroms.contains(&full)
                {
                    let steps = parse_walk_steps(&w.walk)
                        .into_iter()
                        .map(|(s, fwd)| (s.to_string(), fwd))
                        .collect();
                    target_paths.borrow_mut().push(PathRecord {
                        display_chrom: full.clone(),
                        match_keys: vec![w.seq_id, full],
                        start_offset: w.start,
                        steps,
                    });
                }
            }
            Ok(())
        },
    )?;

    let source_paths = source_paths.into_inner();
    let target_paths = target_paths.into_inner();

    if source_paths.is_empty() {
        bail!("no source paths/walks matched '{}'", cfg.from_sample);
    }
    if target_paths.is_empty() {
        bail!("no target paths/walks matched '{}'", cfg.to_sample);
    }

    let mut needed_segs: HashSet<String> = HashSet::new();
    for p in source_paths.iter().chain(target_paths.iter()) {
        for (s, _) in &p.steps {
            needed_segs.insert(s.clone());
        }
    }

    // Pass 2: collect segment lengths only for needed segments (keeps RSS small on 12 GB HPRC).
    let mut seg_lens: HashMap<String, usize> = HashMap::with_capacity(needed_segs.len());
    stream_gfa(
        cfg.gfa_path,
        |s| {
            if needed_segs.contains(&s.name) {
                let len = if s.sequence == b"*" {
                    0
                } else {
                    s.sequence.len()
                };
                seg_lens.insert(s.name, len);
            }
            Ok(())
        },
        |_| Ok(()),
        |_| Ok(()),
        |_| Ok(()),
    )?;

    // Build source StepSpan tables keyed by match_keys.
    let mut source_tables: HashMap<String, (Vec<StepSpan>, usize)> = HashMap::new();
    for p in source_paths {
        let mut cur = p.start_offset;
        let mut spans = Vec::with_capacity(p.steps.len());
        for (seg, is_fwd) in p.steps {
            let slen = seg_lens.get(&seg).copied().unwrap_or(0);
            spans.push(StepSpan {
                seg,
                is_fwd,
                path_start: cur,
                path_end: cur + slen,
            });
            cur += slen;
        }
        for key in p.match_keys {
            source_tables
                .entry(key)
                .or_insert_with(|| (spans.clone(), cur));
        }
    }

    // Build target segment -> TargetOcc index.
    let mut target_index: HashMap<String, Vec<TargetOcc>> =
        HashMap::with_capacity(needed_segs.len());
    let mut target_step_names: Vec<Vec<String>> = Vec::with_capacity(target_paths.len());
    for (path_idx, p) in target_paths.into_iter().enumerate() {
        let mut cur = p.start_offset;
        let mut names = Vec::with_capacity(p.steps.len());
        for (step_idx, (seg, is_fwd)) in p.steps.into_iter().enumerate() {
            let slen = seg_lens.get(&seg).copied().unwrap_or(0);
            names.push(seg.clone());
            target_index.entry(seg).or_default().push(TargetOcc {
                target_chrom: p.display_chrom.clone(),
                target_start: cur,
                target_end: cur + slen,
                is_fwd,
                path_idx,
                step_idx,
            });
            cur += slen;
        }
        target_step_names.push(names);
    }

    let mut summary = LiftSummary {
        lifted_exact: 0,
        lifted_partial: 0,
        lifted_split: 0,
        unmapped: 0,
    };

    for feat in cfg.features {
        let Some((src_steps, _src_len)) = source_tables.get(&feat.chrom) else {
            summary.unmapped += 1;
            if let Some(u) = unmapped_out.as_deref_mut() {
                writeln!(
                    u,
                    "{}\t{}\t{}\t{}\t{}\t{}\tNO_SOURCE_PATH",
                    feat.chrom, feat.start, feat.end, feat.name, feat.score, feat.strand
                )?;
            }
            continue;
        };

        let pieces = lift_feature(
            feat,
            src_steps,
            &seg_lens,
            &target_index,
            &target_step_names,
            cfg.bridge_bp,
        );
        let mapped_bp: usize = pieces.iter().map(|p| p.src_end - p.src_start).sum();
        let total_bp = feat.end.saturating_sub(feat.start).max(1);
        #[allow(clippy::cast_precision_loss)]
        let frac = mapped_bp as f64 / total_bp as f64;

        if pieces.is_empty() || frac < cfg.min_frac {
            summary.unmapped += 1;
            if let Some(u) = unmapped_out.as_deref_mut() {
                writeln!(
                    u,
                    "{}\t{}\t{}\t{}\t{}\t{}\tUNMAPPED({:.3})",
                    feat.chrom, feat.start, feat.end, feat.name, feat.score, feat.strand, frac
                )?;
            }
            continue;
        }

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
            let dup_tag = if p.duplicated { ";DUP" } else { "" };
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}:{}-{}\t{:.4}\t{}{}",
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
                dup_tag
            )?;
        }
    }

    Ok(summary)
}

#[allow(clippy::too_many_lines)]
fn main() -> Result<()> {
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

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--gfa" if i + 1 < args.len() => {
                gfa_path = Some(args[i + 1].clone());
                i += 2;
            }
            "--from" if i + 1 < args.len() => {
                from_sample = Some(args[i + 1].clone());
                i += 2;
            }
            "--to" if i + 1 < args.len() => {
                to_sample = Some(args[i + 1].clone());
                i += 2;
            }
            "--bed" if i + 1 < args.len() => {
                input_path = Some(args[i + 1].clone());
                is_gff = false;
                i += 2;
            }
            "--gff" | "--gff3" if i + 1 < args.len() => {
                input_path = Some(args[i + 1].clone());
                is_gff = true;
                i += 2;
            }
            "--out" | "-o" if i + 1 < args.len() => {
                out_path = Some(args[i + 1].clone());
                i += 2;
            }
            "--unmapped" if i + 1 < args.len() => {
                unmapped_path = Some(args[i + 1].clone());
                i += 2;
            }
            "--min-frac" if i + 1 < args.len() => {
                min_frac = args[i + 1].parse().context("parsing --min-frac")?;
                i += 2;
            }
            "--bridge-bp" if i + 1 < args.len() => {
                bridge_bp = args[i + 1].parse().context("parsing --bridge-bp")?;
                i += 2;
            }
            _ => bail!("unknown argument '{}'", args[i]),
        }
    }

    let (Some(gfa), Some(from), Some(to), Some(inp)) =
        (gfa_path, from_sample, to_sample, input_path)
    else {
        bail!(
            "usage: horizon_lift --gfa <file.gfa[.gz]|-> --from <sample> --to <sample> (--bed <in.bed> | --gff <in.gff3>) [--out <out.bed>] [--unmapped <unmapped.bed>] [--min-frac 0.0]"
        );
    };

    let t0 = Instant::now();
    let features = parse_features(Path::new(&inp), is_gff)?;
    let cfg = LiftConfig {
        gfa_path: &gfa,
        from_sample: &from,
        to_sample: &to,
        features: &features,
        min_frac,
        bridge_bp,
    };

    let mut unmapped_writer = if let Some(ref up) = unmapped_path {
        Some(BufWriter::new(File::create(up)?))
    } else {
        None
    };

    let summary = if let Some(ref op) = out_path {
        let mut w = BufWriter::new(File::create(op)?);
        let s = run_lift(&cfg, &mut w, unmapped_writer.as_mut())?;
        w.flush()?;
        s
    } else {
        let stdout = io::stdout();
        let mut w = BufWriter::new(stdout.lock());
        let s = run_lift(&cfg, &mut w, unmapped_writer.as_mut())?;
        w.flush()?;
        s
    };
    if let Some(mut uw) = unmapped_writer {
        uw.flush()?;
    }

    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
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

    #[test]
    fn lifts_forward_inverted_split_and_unmapped_intervals() {
        // Construct a small GFA with GRCh38 (P line) and CHM13 (W line)
        // Segments: s1=10bp, s2=10bp, s3=10bp (ref-only insertion), s4=10bp
        // GRCh38: s1+, s2+, s3+, s4+ -> total 40bp (0..10, 10..20, 20..30, 30..40)
        // CHM13:  >s1 <s2 >s4       -> total 30bp (0..10=s1+, 10..20=s2-, 20..30=s4+)
        let gfa_path = std::env::temp_dir().join("horizon_lift_test.gfa");
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
        let gfa_path = std::env::temp_dir().join("horizon_lift_jaccard_test.gfa");
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
}
