//! Core GFA walk coordinate projection and micro-bubble coalescing for `horizon_lift`.

use crate::features::InputFeature;
use anyhow::{Result, bail};
use horizon_aligner::coord::Interval;
use horizon_aligner::gfa::{GfaRecord, parse_path_steps, parse_walk_steps, stream_gfa_records};
use horizon_aligner::strip_sample_prefix;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::ops::ControlFlow;

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
    rel_fwd: bool,
    src_start: usize,
    src_end: usize,
    duplicated: bool,
}

pub(crate) fn flip_strand(strand: char) -> char {
    match strand {
        '+' => '-',
        '-' => '+',
        other => other,
    }
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
    let mut src_counts: Vec<(&str, usize)> = Vec::with_capacity(s_hi - s_lo);
    for st in &src_steps[s_lo..s_hi] {
        let key = st.seg.as_str();
        if let Some((_, c)) = src_counts.iter_mut().find(|(k, _)| *k == key) {
            *c += 1;
        } else {
            src_counts.push((key, 1));
        }
    }
    let mut tgt_counts: Vec<(&str, usize)> = Vec::with_capacity(t_hi - t_lo);
    for seg in &tgt_steps[t_lo..t_hi] {
        let key = seg.as_str();
        if let Some((_, c)) = tgt_counts.iter_mut().find(|(k, _)| *k == key) {
            *c += 1;
        } else {
            tgt_counts.push((key, 1));
        }
    }
    let mut inter = 0usize;
    let mut union = 0usize;
    for &(k, sc) in &src_counts {
        let tc = tgt_counts
            .iter()
            .find(|(tk, _)| *tk == k)
            .map_or(0, |&(_, c)| c);
        inter += sc.min(tc);
        union += sc.max(tc);
    }
    for &(k, tc) in &tgt_counts {
        if !src_counts.iter().any(|(sk, _)| *sk == k) {
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
            occs.iter()
                .map(|occ| {
                    let j = context_jaccard(
                        src_steps,
                        src_pos,
                        &target_step_names[occ.path_idx],
                        occ.step_idx,
                    );
                    (occ, j, dist(occ))
                })
                .max_by(|(_, ja, da), (_, jb, db)| ja.total_cmp(jb).then_with(|| db.cmp(da)))
                .map_or(&occs[0], |(occ, _, _)| occ)
        };

        let (t_s, t_e) = if chosen.is_fwd {
            (chosen.target_start + seg_s, chosen.target_start + seg_e)
        } else {
            (
                chosen.target_start + (seg_len - seg_e),
                chosen.target_start + (seg_len - seg_s),
            )
        };
        let rel_fwd = chosen.is_fwd == step.is_fwd;
        let out_strand = if rel_fwd {
            feat.strand
        } else {
            flip_strand(feat.strand)
        };
        let piece = LiftedPiece {
            chrom: chosen.target_chrom.clone(),
            start: t_s,
            end: t_e,
            strand: out_strand,
            rel_fwd,
            src_start: ov_s,
            src_end: ov_e,
            duplicated,
        };
        last_target = Some((piece.chrom.clone(), piece.start, piece.end, piece.strand));

        // Coalesce contiguous (or within bridge_bp on the same chrom, strand, and relative orientation) pieces.
        if let Some(prev) = raw_pieces.last_mut() {
            let same_ctx = prev.chrom == piece.chrom
                && prev.strand == piece.strand
                && prev.rel_fwd == piece.rel_fwd;
            if same_ctx
                && piece.rel_fwd
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
                && !piece.rel_fwd
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

pub(crate) struct LiftConfig<'a> {
    pub(crate) gfa_path: &'a str,
    pub(crate) from_sample: &'a str,
    pub(crate) to_sample: &'a str,
    pub(crate) to_chroms: &'a HashSet<String>,
    pub(crate) features: &'a [InputFeature],
    pub(crate) min_frac: f64,
    pub(crate) bridge_bp: usize,
}

pub(crate) struct LiftSummary {
    pub(crate) lifted_exact: usize,
    pub(crate) lifted_partial: usize,
    pub(crate) lifted_split: usize,
    pub(crate) unmapped: usize,
}

struct LiftPaths {
    source_paths: Vec<PathRecord>,
    target_paths: Vec<PathRecord>,
}

fn record_lift_path(
    p: horizon_aligner::gfa::GfaPath,
    cfg: &LiftConfig<'_>,
    wanted_chroms: &HashSet<String>,
    source_paths: &mut Vec<PathRecord>,
    target_paths: &mut Vec<PathRecord>,
) {
    if horizon_aligner::matches_sample_prefix(&p.name, cfg.from_sample) {
        let short = strip_sample_prefix(&p.name, cfg.from_sample).to_string();
        if wanted_chroms.is_empty()
            || wanted_chroms.contains(&p.name)
            || wanted_chroms.contains(&short)
        {
            let steps = parse_path_steps(&p.segment_names)
                .into_iter()
                .map(|(s, fwd)| (s.to_string(), fwd))
                .collect();
            let match_keys = if short == p.name {
                vec![p.name.clone()]
            } else {
                vec![p.name.clone(), short]
            };
            source_paths.push(PathRecord {
                display_chrom: p.name.clone(),
                match_keys,
                start_offset: 0,
                steps,
            });
        }
    }
    if horizon_aligner::matches_sample_prefix(&p.name, cfg.to_sample) {
        let short = strip_sample_prefix(&p.name, cfg.to_sample);
        if cfg.to_chroms.is_empty()
            || cfg.to_chroms.contains(&p.name)
            || cfg.to_chroms.contains(short)
        {
            let steps = parse_path_steps(&p.segment_names)
                .into_iter()
                .map(|(s, fwd)| (s.to_string(), fwd))
                .collect();
            target_paths.push(PathRecord {
                display_chrom: p.name,
                match_keys: Vec::new(),
                start_offset: 0,
                steps,
            });
        }
    }
}

fn record_lift_walk(
    w: &horizon_aligner::gfa::GfaWalk,
    cfg: &LiftConfig<'_>,
    wanted_chroms: &HashSet<String>,
    source_paths: &mut Vec<PathRecord>,
    target_paths: &mut Vec<PathRecord>,
) {
    if horizon_aligner::matches_sample_prefix(&w.sample, cfg.from_sample) {
        let full = format!("{}.{}", w.sample, w.seq_id);
        let hap_full = if w.haplotype == "0" || w.haplotype.is_empty() {
            None
        } else {
            Some(format!("{}#{}.{}", w.sample, w.haplotype, w.seq_id))
        };
        if wanted_chroms.is_empty()
            || wanted_chroms.contains(&w.seq_id)
            || wanted_chroms.contains(&full)
            || hap_full.as_ref().is_some_and(|h| wanted_chroms.contains(h))
        {
            let steps = parse_walk_steps(&w.walk)
                .into_iter()
                .map(|(s, fwd)| (s.to_string(), fwd))
                .collect();
            let mut match_keys = vec![w.seq_id.clone(), full.clone()];
            if let Some(ref hf) = hap_full {
                match_keys.push(hf.clone());
            }
            source_paths.push(PathRecord {
                display_chrom: hap_full.unwrap_or(full),
                match_keys,
                start_offset: w.start,
                steps,
            });
        }
    }
    if horizon_aligner::matches_sample_prefix(&w.sample, cfg.to_sample) {
        let disp = if w.haplotype == "0" || w.haplotype.is_empty() {
            format!("{}.{}", w.sample, w.seq_id)
        } else {
            format!("{}#{}.{}", w.sample, w.haplotype, w.seq_id)
        };
        if cfg.to_chroms.is_empty()
            || cfg.to_chroms.contains(&w.seq_id)
            || cfg.to_chroms.contains(&disp)
        {
            let steps = parse_walk_steps(&w.walk)
                .into_iter()
                .map(|(s, fwd)| (s.to_string(), fwd))
                .collect();
            target_paths.push(PathRecord {
                display_chrom: disp,
                match_keys: Vec::new(),
                start_offset: w.start,
                steps,
            });
        }
    }
}

fn collect_lift_paths(cfg: &LiftConfig<'_>, wanted_chroms: &HashSet<String>) -> Result<LiftPaths> {
    let mut source_paths = Vec::<PathRecord>::new();
    let mut target_paths = Vec::<PathRecord>::new();

    stream_gfa_records(cfg.gfa_path, |rec| {
        match rec {
            GfaRecord::Path(p) => {
                record_lift_path(p, cfg, wanted_chroms, &mut source_paths, &mut target_paths);
            }
            GfaRecord::Walk(w) => {
                record_lift_walk(&w, cfg, wanted_chroms, &mut source_paths, &mut target_paths);
            }
            _ => {}
        }
        Ok(ControlFlow::Continue(()))
    })?;

    if source_paths.is_empty() {
        bail!("no source paths/walks matched '{}'", cfg.from_sample);
    }
    if target_paths.is_empty() {
        bail!("no target paths/walks matched '{}'", cfg.to_sample);
    }
    Ok(LiftPaths {
        source_paths,
        target_paths,
    })
}

struct LiftTables {
    seg_lens: HashMap<String, usize>,
    source_walks: Vec<(Vec<StepSpan>, usize, usize)>,
    source_tables: HashMap<String, Vec<usize>>,
    target_index: HashMap<String, Vec<TargetOcc>>,
    target_step_names: Vec<Vec<String>>,
}

fn build_lift_tables(gfa_path: &str, paths: LiftPaths) -> Result<LiftTables> {
    let mut needed_segs: HashSet<String> = HashSet::new();
    for p in paths.source_paths.iter().chain(paths.target_paths.iter()) {
        for (s, _) in &p.steps {
            needed_segs.insert(s.clone());
        }
    }

    // Pass 2: collect segment lengths only for needed segments (keeps RSS small on 12 GB HPRC).
    let mut seg_lens: HashMap<String, usize> = HashMap::with_capacity(needed_segs.len());
    stream_gfa_records(gfa_path, |rec| {
        if let GfaRecord::Segment(s) = rec
            && needed_segs.contains(&s.name)
        {
            seg_lens.insert(s.name.clone(), s.seq_len());
        }
        Ok(ControlFlow::Continue(()))
    })?;

    let mut source_walks: Vec<(Vec<StepSpan>, usize, usize)> =
        Vec::with_capacity(paths.source_paths.len());
    let mut source_tables: HashMap<String, Vec<usize>> = HashMap::new();
    for p in paths.source_paths {
        let start_off = p.start_offset;
        let mut cur = start_off;
        let mut spans = Vec::with_capacity(p.steps.len());
        for (seg, is_fwd) in p.steps {
            let slen = seg_lens.get(&seg).copied().unwrap_or(0);
            let next_cur = cur.saturating_add(slen);
            spans.push(StepSpan {
                seg,
                is_fwd,
                path_start: cur,
                path_end: next_cur,
            });
            cur = next_cur;
        }
        let walk_idx = source_walks.len();
        source_walks.push((spans, start_off, cur));
        for key in p.match_keys {
            source_tables.entry(key).or_default().push(walk_idx);
        }
    }

    let mut target_index: HashMap<String, Vec<TargetOcc>> =
        HashMap::with_capacity(needed_segs.len());
    let mut target_step_names: Vec<Vec<String>> = Vec::with_capacity(paths.target_paths.len());
    for (path_idx, p) in paths.target_paths.into_iter().enumerate() {
        let mut cur = p.start_offset;
        let mut names = Vec::with_capacity(p.steps.len());
        for (step_idx, (seg, is_fwd)) in p.steps.into_iter().enumerate() {
            let slen = seg_lens.get(&seg).copied().unwrap_or(0);
            let next_cur = cur.saturating_add(slen);
            names.push(seg.clone());
            target_index.entry(seg).or_default().push(TargetOcc {
                target_chrom: p.display_chrom.clone(),
                target_start: cur,
                target_end: next_cur,
                is_fwd,
                path_idx,
                step_idx,
            });
            cur = next_cur;
        }
        target_step_names.push(names);
    }

    Ok(LiftTables {
        seg_lens,
        source_walks,
        source_tables,
        target_index,
        target_step_names,
    })
}

fn merged_mapped_bp(pieces: &[LiftedPiece]) -> usize {
    if pieces.is_empty() {
        return 0;
    }
    let mut covered: Vec<(usize, usize)> = pieces
        .iter()
        .filter(|p| p.src_end > p.src_start)
        .map(|p| (p.src_start, p.src_end))
        .collect();
    covered.sort_unstable();
    let mut total = 0usize;
    let mut cur_s = 0usize;
    let mut cur_e = 0usize;
    for (s, e) in covered {
        if s > cur_e {
            total += cur_e - cur_s;
            cur_s = s;
            cur_e = e;
        } else {
            cur_e = cur_e.max(e);
        }
    }
    total + (cur_e - cur_s)
}

pub(crate) fn run_lift<W1: Write, W2: Write>(
    cfg: &LiftConfig<'_>,
    out: &mut W1,
    mut unmapped_out: Option<&mut W2>,
) -> Result<LiftSummary> {
    let wanted_chroms: HashSet<String> = cfg.features.iter().map(|f| f.chrom.clone()).collect();
    let paths = collect_lift_paths(cfg, &wanted_chroms)?;
    let tables = build_lift_tables(cfg.gfa_path, paths)?;

    let mut summary = LiftSummary {
        lifted_exact: 0,
        lifted_partial: 0,
        lifted_split: 0,
        unmapped: 0,
    };

    for feat in cfg.features {
        let Some(src_walk_idxs) = tables.source_tables.get(&feat.chrom) else {
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
        let total_bp = feat.end.saturating_sub(feat.start).max(1);
        let mut mapped_any = false;
        let mut best_failed_frac = 0.0_f64;

        for &w_idx in src_walk_idxs {
            let (ref src_steps, s0, s1) = tables.source_walks[w_idx];
            if feat.end <= s0 || feat.start >= s1 {
                continue;
            }
            let pieces = lift_feature(
                feat,
                src_steps,
                &tables.seg_lens,
                &tables.target_index,
                &tables.target_step_names,
                cfg.bridge_bp,
            );
            let mapped_bp = merged_mapped_bp(&pieces);
            #[allow(clippy::cast_precision_loss)]
            let frac = mapped_bp as f64 / total_bp as f64;
            if pieces.is_empty() || frac < cfg.min_frac {
                best_failed_frac = best_failed_frac.max(frac);
                continue;
            }
            mapped_any = true;
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

        if !mapped_any {
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
