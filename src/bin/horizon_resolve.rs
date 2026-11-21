//! Two-pass streaming GFA v1.0/v1.1 snarl/bubble e-graph resolver using `MdlCostModel`.

use anyhow::{Context, Result, bail};
use horizon_aligner::egraph::{
    EGraph, Expr, MdlCostModel, build_seq_expr, extract_min_mdl_cost, saturate_synteny,
};
use horizon_aligner::gfa::{
    GfaRecord, SpoolableGfaInput, open_output_writer, parse_path_steps, parse_walk_steps,
    stream_gfa_records,
};
use std::collections::{HashMap, HashSet};
use std::env;
use std::io::Write;
use std::ops::ControlFlow;
use std::time::Instant;

struct CliArgs {
    gfa_in: String,
    out_path: Option<String>,
    collapse_dups: bool,
    tag_walks: bool,
    max_bubble_steps: usize,
    min_anchor_bp: usize,
    divergence: Option<f64>,
}

#[derive(Default)]
struct ResolveStats {
    paths_seen: usize,
    walks_seen: usize,
    windows_evaluated: usize,
    windows_simplified: usize,
    dups_detected: usize,
    dups_collapsed: usize,
    invs_cancelled: usize,
    total_mdl_bits_saved: f64,
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

fn flatten_expr(expr: &Expr, flip: bool, collapse_dups: bool, out: &mut Vec<(String, bool)>) {
    match expr {
        Expr::Block(lbl) => out.push((lbl.clone(), !flip)),
        Expr::Inv(inner) => flatten_expr(inner, !flip, collapse_dups, out),
        Expr::Seq(a, b) => {
            if flip {
                flatten_expr(b, true, collapse_dups, out);
                flatten_expr(a, true, collapse_dups, out);
            } else {
                flatten_expr(a, false, collapse_dups, out);
                flatten_expr(b, false, collapse_dups, out);
            }
        }
        Expr::Dup(inner) => {
            flatten_expr(inner, flip, collapse_dups, out);
            if !collapse_dups {
                flatten_expr(inner, flip, collapse_dups, out);
            }
        }
    }
}

fn format_expr(expr: &Expr) -> String {
    match expr {
        Expr::Block(lbl) => lbl.clone(),
        Expr::Inv(inner) => format!("inv({})", format_expr(inner)),
        Expr::Dup(inner) => format!("dup({})", format_expr(inner)),
        Expr::Seq(a, b) => format!("seq({},{})", format_expr(a), format_expr(b)),
    }
}

fn count_ops(expr: &Expr) -> (usize, usize) {
    match expr {
        Expr::Block(_) => (0, 0),
        Expr::Inv(inner) => {
            let (d, inv) = count_ops(inner);
            (d, inv + 1)
        }
        Expr::Dup(inner) => {
            let (d, inv) = count_ops(inner);
            (d + 1, inv)
        }
        Expr::Seq(a, b) => {
            let (d1, i1) = count_ops(a);
            let (d2, i2) = count_ops(b);
            (d1 + d2, i1 + i2)
        }
    }
}

fn resolve_window(
    steps: &[(&str, bool)],
    seg_lens: &HashMap<String, usize>,
    divergence: f64,
    collapse_dups: bool,
    stats: &mut ResolveStats,
) -> (Vec<(String, bool)>, Option<String>) {
    if steps.len() <= 1 {
        return (
            steps.iter().map(|&(s, f)| (s.to_string(), f)).collect(),
            None,
        );
    }

    stats.windows_evaluated += 1;
    let blocks: Vec<(&str, i8)> = steps
        .iter()
        .map(|&(s, f)| (s, if f { 1 } else { -1 }))
        .collect();

    let mut eg = EGraph::new();
    let Some(root) = build_seq_expr(&mut eg, &blocks) else {
        return (
            steps.iter().map(|&(s, f)| (s.to_string(), f)).collect(),
            None,
        );
    };

    let mut model = MdlCostModel::from_egraph(&eg);
    for &(seg, _) in steps {
        let len_bp = seg_lens.get(seg).copied().unwrap_or(100);
        model
            .block_stats
            .entry(seg.to_string())
            .or_insert((len_bp, divergence));
    }

    let raw_cost = extract_min_mdl_cost(&eg, root, &model).map_or(0.0, |(c, _)| c);
    saturate_synteny(&mut eg, 6, 4096);
    let Some((opt_cost, expr)) = extract_min_mdl_cost(&eg, root, &model) else {
        return (
            steps.iter().map(|&(s, f)| (s.to_string(), f)).collect(),
            None,
        );
    };

    let (dups, _) = count_ops(&expr);
    stats.dups_detected += dups;

    let mut resolved = Vec::with_capacity(steps.len());
    flatten_expr(&expr, false, collapse_dups, &mut resolved);

    let invs_before = steps.iter().filter(|&&(_, f)| !f).count();
    let invs_after_uncollapsed = if collapse_dups {
        let mut uncollapsed = Vec::with_capacity(steps.len());
        flatten_expr(&expr, false, false, &mut uncollapsed);
        uncollapsed.iter().filter(|(_, f)| !*f).count()
    } else {
        resolved.iter().filter(|(_, f)| !*f).count()
    };
    if invs_before > invs_after_uncollapsed {
        stats.invs_cancelled += invs_before - invs_after_uncollapsed;
    }

    let changed = resolved.len() != steps.len()
        || resolved
            .iter()
            .zip(steps.iter())
            .any(|((a, fa), &(b, fb))| a != b || *fa != fb);

    if changed {
        stats.windows_simplified += 1;
        if collapse_dups && resolved.len() < steps.len() {
            stats.dups_collapsed += dups;
        }
        if raw_cost > opt_cost {
            stats.total_mdl_bits_saved += raw_cost - opt_cost;
        }
        (resolved, Some(format_expr(&expr)))
    } else {
        (
            steps.iter().map(|&(s, f)| (s.to_string(), f)).collect(),
            None,
        )
    }
}

fn estimate_window_divergence(steps: &[(&str, bool)], anchors: &HashSet<String>) -> f64 {
    if steps.is_empty() {
        return 0.02;
    }
    let mut non_anchors = 0usize;
    let mut rev_steps = 0usize;
    for &(seg, fwd) in steps {
        if !anchors.contains(seg) {
            non_anchors += 1;
        }
        if !fwd {
            rev_steps += 1;
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let n = steps.len() as f64;
    #[allow(clippy::cast_precision_loss)]
    let non_anchor_frac = non_anchors as f64 / n;
    #[allow(clippy::cast_precision_loss)]
    let inv_frac = rev_steps as f64 / n;
    (0.02 + 0.08 * non_anchor_frac + 0.05 * inv_frac).clamp(0.01, 0.25)
}

fn decompose_and_resolve(
    steps: &[(&str, bool)],
    anchors: &HashSet<String>,
    seg_lens: &HashMap<String, usize>,
    max_bubble_steps: usize,
    divergence_override: Option<f64>,
    collapse_dups: bool,
    stats: &mut ResolveStats,
) -> (Vec<(String, bool)>, Vec<String>) {
    let mut out_steps = Vec::with_capacity(steps.len());
    let mut tags = Vec::new();
    let mut i = 0usize;
    while i < steps.len() {
        let mut j = i + 1;
        let mut saw_non_anchor = !anchors.contains(steps[i].0);
        while j < steps.len() && j - i < max_bubble_steps {
            let is_anchor = anchors.contains(steps[j].0);
            let is_tandem = steps[j].0 == steps[j - 1].0;
            if is_anchor && saw_non_anchor && !is_tandem && j - i >= 2 {
                break;
            }
            if !is_anchor {
                saw_non_anchor = true;
            }
            j += 1;
        }
        let win_div = divergence_override
            .unwrap_or_else(|| estimate_window_divergence(&steps[i..j], anchors));
        let (win_steps, win_tag) =
            resolve_window(&steps[i..j], seg_lens, win_div, collapse_dups, stats);
        out_steps.extend(win_steps);
        if let Some(t) = win_tag {
            tags.push(t);
        }
        i = j;
    }
    (out_steps, tags)
}

fn pass1_index(
    gfa_path: &str,
    min_anchor_bp: usize,
) -> Result<(HashMap<String, usize>, HashSet<String>)> {
    let mut seg_lens: HashMap<String, usize> = HashMap::new();
    let mut seg_visits: HashMap<String, usize> = HashMap::new();

    stream_gfa_records(gfa_path, |rec| {
        match rec {
            GfaRecord::Segment(s) => {
                let len = s.seq_len();
                seg_lens.insert(s.name, len);
            }
            GfaRecord::Link(_) => {}
            GfaRecord::Path(p) => {
                let mut seen_in_path = HashSet::new();
                for (seg, _) in parse_path_steps(&p.segment_names) {
                    if seen_in_path.insert(seg) {
                        *seg_visits.entry(seg.to_string()).or_insert(0) += 1;
                    }
                }
            }
            GfaRecord::Walk(w) => {
                let mut seen_in_walk = HashSet::new();
                for (seg, _) in parse_walk_steps(&w.walk) {
                    if seen_in_walk.insert(seg) {
                        *seg_visits.entry(seg.to_string()).or_insert(0) += 1;
                    }
                }
            }
        }
        Ok(ControlFlow::Continue(()))
    })?;

    let mut anchors = HashSet::new();
    for (seg, visits) in seg_visits {
        let len = seg_lens.get(&seg).copied().unwrap_or(0);
        if visits >= 2 || len >= min_anchor_bp {
            anchors.insert(seg);
        }
    }
    Ok((seg_lens, anchors))
}

struct ResolveCtx<'a> {
    seg_lens: &'a HashMap<String, usize>,
    anchors: &'a HashSet<String>,
    cli: &'a CliArgs,
}

fn intern_resolve_seg(
    seg_to_id: &mut HashMap<String, u32>,
    id_to_seg: &mut Vec<String>,
    name: &str,
) -> u32 {
    if let Some(&id) = seg_to_id.get(name) {
        id
    } else {
        let id = u32::try_from(id_to_seg.len()).expect("segment count fits u32");
        let owned = name.to_owned();
        seg_to_id.insert(owned.clone(), id);
        id_to_seg.push(owned);
        id
    }
}

fn record_interned_edges(
    steps: &[(String, bool)],
    seg_to_id: &mut HashMap<String, u32>,
    id_to_seg: &mut Vec<String>,
    used: &mut HashSet<(u32, bool, u32, bool)>,
) {
    for win in steps.windows(2) {
        let u = intern_resolve_seg(seg_to_id, id_to_seg, &win[0].0);
        let v = intern_resolve_seg(seg_to_id, id_to_seg, &win[1].0);
        let uf = win[0].1;
        let vf = win[1].1;
        used.insert((u, uf, v, vf));
        used.insert((v, !vf, u, !uf));
    }
}

fn emit_resolved_path<W: Write>(
    p: &horizon_aligner::gfa::GfaPath,
    ctx: &ResolveCtx<'_>,
    seg_to_id: &mut HashMap<String, u32>,
    id_to_seg: &mut Vec<String>,
    used_edges: &mut HashSet<(u32, bool, u32, bool)>,
    stats: &mut ResolveStats,
    out: &mut W,
) -> Result<()> {
    stats.paths_seen += 1;
    let steps = parse_path_steps(&p.segment_names);
    let (resolved, _) = decompose_and_resolve(
        &steps,
        ctx.anchors,
        ctx.seg_lens,
        ctx.cli.max_bubble_steps,
        ctx.cli.divergence,
        ctx.cli.collapse_dups,
        stats,
    );
    record_interned_edges(&resolved, seg_to_id, id_to_seg, used_edges);
    let mut seg_str = String::new();
    for (idx, (s, fwd)) in resolved.iter().enumerate() {
        if idx > 0 {
            seg_str.push(',');
        }
        seg_str.push_str(s);
        seg_str.push(if *fwd { '+' } else { '-' });
    }
    let same_steps = resolved.len() == steps.len()
        && resolved
            .iter()
            .zip(steps.iter())
            .all(|((a, fa), &(b, fb))| a == b && *fa == fb);
    let overlaps = if same_steps {
        p.overlaps.as_str()
    } else {
        "*"
    };
    writeln!(out, "P\t{}\t{}\t{}", p.name, seg_str, overlaps)?;
    Ok(())
}

fn emit_resolved_walk<W: Write>(
    w: &horizon_aligner::gfa::GfaWalk,
    ctx: &ResolveCtx<'_>,
    seg_to_id: &mut HashMap<String, u32>,
    id_to_seg: &mut Vec<String>,
    used_edges: &mut HashSet<(u32, bool, u32, bool)>,
    stats: &mut ResolveStats,
    out: &mut W,
) -> Result<()> {
    stats.walks_seen += 1;
    let steps = parse_walk_steps(&w.walk);
    let (resolved, tags) = decompose_and_resolve(
        &steps,
        ctx.anchors,
        ctx.seg_lens,
        ctx.cli.max_bubble_steps,
        ctx.cli.divergence,
        ctx.cli.collapse_dups,
        stats,
    );
    record_interned_edges(&resolved, seg_to_id, id_to_seg, used_edges);
    let mut walk_str = String::new();
    let mut resolved_bp = 0usize;
    let mut all_known = true;
    for (s, fwd) in &resolved {
        walk_str.push(if *fwd { '>' } else { '<' });
        walk_str.push_str(s);
        let slen = ctx.seg_lens.get(s.as_str()).copied().unwrap_or(0);
        if slen == 0 {
            all_known = false;
        } else {
            resolved_bp += slen;
        }
    }
    let same_steps = resolved.len() == steps.len()
        && resolved
            .iter()
            .zip(steps.iter())
            .all(|((a, fa), &(b, fb))| a == b && *fa == fb);
    let end_pos = if !same_steps && all_known && resolved_bp > 0 {
        w.start.saturating_add(resolved_bp)
    } else {
        w.end
    };
    let extra_tags: String = if let Some(ref raw_tags) = w.tags {
        let kept: Vec<&str> = raw_tags
            .split('\t')
            .filter(|t| !(ctx.cli.tag_walks && !tags.is_empty() && t.starts_with("EG:Z:")))
            .collect();
        if kept.is_empty() {
            String::new()
        } else {
            format!("\t{}", kept.join("\t"))
        }
    } else {
        String::new()
    };
    if ctx.cli.tag_walks && !tags.is_empty() {
        writeln!(
            out,
            "W\t{}\t{}\t{}\t{}\t{}\t{}{}\tEG:Z:{}",
            w.sample,
            w.haplotype,
            w.seq_id,
            w.start,
            end_pos,
            walk_str,
            extra_tags,
            tags.join(";")
        )?;
    } else {
        writeln!(
            out,
            "W\t{}\t{}\t{}\t{}\t{}\t{}{}",
            w.sample, w.haplotype, w.seq_id, w.start, end_pos, walk_str, extra_tags
        )?;
    }
    Ok(())
}

/// Pass 2: rewrite `P`/`W` records and emit only `L` links still traversed by resolved walks.
///
/// Note: buffers `L` lines in `raw_links` so unused links can be pruned after all `P`/`W` lines
/// are seen; memory scales with the number of `L` records in the input GFA.
fn pass2_rewrite<W: Write>(
    gfa_path: &str,
    out: &mut W,
    seg_lens: &HashMap<String, usize>,
    anchors: &HashSet<String>,
    cli: &CliArgs,
) -> Result<ResolveStats> {
    let ctx = ResolveCtx {
        seg_lens,
        anchors,
        cli,
    };
    let mut stats = ResolveStats::default();
    let mut raw_links: Vec<horizon_aligner::gfa::GfaLink> = Vec::new();
    let mut seg_to_id: HashMap<String, u32> = HashMap::new();
    let mut id_to_seg: Vec<String> = Vec::new();
    let mut used_edges: HashSet<(u32, bool, u32, bool)> = HashSet::new();

    writeln!(out, "H\tVN:Z:1.1")?;
    stream_gfa_records(gfa_path, |rec| {
        match rec {
            GfaRecord::Segment(s) => {
                let seq_str = String::from_utf8_lossy(&s.sequence);
                if let Some(ref tags) = s.tags {
                    writeln!(out, "S\t{}\t{}\t{}", s.name, seq_str, tags)?;
                } else {
                    writeln!(out, "S\t{}\t{}", s.name, seq_str)?;
                }
            }
            GfaRecord::Link(l) => {
                raw_links.push(l);
            }
            GfaRecord::Path(p) => {
                emit_resolved_path(
                    &p,
                    &ctx,
                    &mut seg_to_id,
                    &mut id_to_seg,
                    &mut used_edges,
                    &mut stats,
                    out,
                )?;
            }
            GfaRecord::Walk(w) => {
                emit_resolved_walk(
                    &w,
                    &ctx,
                    &mut seg_to_id,
                    &mut id_to_seg,
                    &mut used_edges,
                    &mut stats,
                    out,
                )?;
            }
        }
        Ok(ControlFlow::Continue(()))
    })?;

    let has_paths_or_walks = stats.paths_seen + stats.walks_seen > 0;
    let mut seen_edges: HashSet<(u32, bool, u32, bool)> = HashSet::new();
    for l in raw_links {
        let uf = l.from_orient != b'-';
        let vf = l.to_orient != b'-';
        let u = intern_resolve_seg(&mut seg_to_id, &mut id_to_seg, &l.from);
        let v = intern_resolve_seg(&mut seg_to_id, &mut id_to_seg, &l.to);
        let fwd_edge = (u, uf, v, vf);
        if !has_paths_or_walks || used_edges.contains(&fwd_edge) {
            writeln!(
                out,
                "L\t{}\t{}\t{}\t{}\t{}",
                l.from, l.from_orient as char, l.to, l.to_orient as char, l.overlap
            )?;
            seen_edges.insert((v, !vf, u, !uf));
            seen_edges.insert(fwd_edge);
        }
    }
    let mut missing: Vec<(&str, bool, &str, bool, u32, u32)> = used_edges
        .into_iter()
        .filter(|e| !seen_edges.contains(e))
        .map(|(u, uf, v, vf)| (id_to_seg[u as usize].as_str(), uf, id_to_seg[v as usize].as_str(), vf, u, v))
        .collect();
    missing.sort_unstable();
    for (u_name, uf, v_name, vf, u, v) in missing {
        let rc = (v, !vf, u, !uf);
        if !seen_edges.contains(&(u, uf, v, vf)) && !seen_edges.contains(&rc) {
            writeln!(
                out,
                "L\t{}\t{}\t{}\t{}\t0M",
                u_name,
                if uf { '+' } else { '-' },
                v_name,
                if vf { '+' } else { '-' }
            )?;
            seen_edges.insert((u, uf, v, vf));
            seen_edges.insert(rc);
        }
    }

    out.flush()?;
    Ok(stats)
}

fn main() -> Result<()> {
    let Some(cli) = parse_args()? else {
        return Ok(());
    };
    let t0 = Instant::now();
    let spooled = SpoolableGfaInput::new(&cli.gfa_in)?;
    let effective_gfa = spooled.path();
    let (seg_lens, anchors) = pass1_index(effective_gfa, cli.min_anchor_bp)?;

    let mut writer = open_output_writer(cli.out_path.as_deref())?;
    let stats = pass2_rewrite(effective_gfa, &mut writer, &seg_lens, &anchors, &cli)?;
    writer.flush()?;

    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    eprintln!(
        "segments={} anchors={} paths={} walks={} windows={} simplified={} dups_detected={} dups_collapsed={} invs_cancelled={} mdl_bits_saved={:.1} elapsed_s={:.3} peak_rss_mb={:.1}",
        seg_lens.len(),
        anchors.len(),
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
        let (seg_lens, anchors) = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&gfa_out)?);
        let stats = pass2_rewrite(&cli.gfa_in, &mut writer, &seg_lens, &anchors, &cli)?;
        assert_eq!(stats.walks_seen, 2);
        assert_eq!(stats.dups_detected, 2);
        assert_eq!(stats.dups_collapsed, 2);
        assert!(stats.total_mdl_bits_saved > 0.0);

        let out_text = std::fs::read_to_string(&gfa_out)?;
        assert!(out_text.contains("W\tHG002\t1\tchr6\t0\t60\t>a>b>c\tEG:Z:"));
        assert!(out_text.contains("W\tHG002\t2\tchr7\t0\t60\t>a<b>c\tEG:Z:"));
        assert!(out_text.contains("L\ta\t+\tb\t+\t0M"));
        assert!(out_text.contains("L\tb\t+\tc\t+\t0M"));
        assert!(!out_text.contains("L\tb\t+\tb\t+\t0M"));
        let _ = std::fs::remove_dir_all(&dir);
        Ok(())
    }
}
