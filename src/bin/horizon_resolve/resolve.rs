use anyhow::Result;
use horizon_aligner::egraph::{
    EGraph, Expr, MdlCostModel, build_seq_expr, extract_min_mdl_cost, saturate_synteny,
};
use horizon_aligner::gfa::{
    GfaLink, GfaPath, GfaRecord, GfaWalk, parse_path_steps, parse_walk_steps, stream_gfa_records,
};
use std::collections::HashSet;
use std::io::Write;
use std::ops::ControlFlow;

use crate::CliArgs;
use crate::snarl::PackedGraph;

#[derive(Default)]
pub(crate) struct ResolveStats {
    pub(crate) paths_seen: usize,
    pub(crate) walks_seen: usize,
    pub(crate) windows_evaluated: usize,
    pub(crate) windows_simplified: usize,
    pub(crate) dups_detected: usize,
    pub(crate) dups_collapsed: usize,
    pub(crate) invs_cancelled: usize,
    pub(crate) total_mdl_bits_saved: f64,
}

fn flatten_expr<'a>(
    expr: &'a Expr,
    flip: bool,
    collapse_dups: bool,
    out: &mut Vec<(&'a str, bool)>,
) {
    match expr {
        Expr::Block(lbl) => out.push((lbl.as_str(), !flip)),
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
    steps: &[(u32, bool)],
    graph: &PackedGraph,
    divergence: f64,
    collapse_dups: bool,
    stats: &mut ResolveStats,
) -> (Vec<(u32, bool)>, Option<String>) {
    if steps.len() <= 1 {
        return (steps.to_vec(), None);
    }

    stats.windows_evaluated += 1;
    let blocks: Vec<(&str, i8)> = steps
        .iter()
        .map(|&(id, f)| {
            (
                graph.id_to_seg[id as usize].as_str(),
                if f { 1 } else { -1 },
            )
        })
        .collect();

    let mut eg = EGraph::new();
    let Some(root) = build_seq_expr(&mut eg, &blocks) else {
        return (steps.to_vec(), None);
    };

    let mut model = MdlCostModel::from_egraph(&eg);
    for &(id, _) in steps {
        let seg_name = &graph.id_to_seg[id as usize];
        let len_bp = graph.seg_lens[id as usize].max(100);
        model
            .block_stats
            .entry(seg_name.clone())
            .or_insert((len_bp, divergence));
    }

    let raw_cost = extract_min_mdl_cost(&eg, root, &model).map_or(0.0, |(c, _)| c);
    saturate_synteny(&mut eg, 6, 4096);
    let Some((opt_cost, expr)) = extract_min_mdl_cost(&eg, root, &model) else {
        return (steps.to_vec(), None);
    };

    let (dups, _) = count_ops(&expr);
    stats.dups_detected += dups;

    let mut resolved_names = Vec::with_capacity(steps.len());
    flatten_expr(&expr, false, collapse_dups, &mut resolved_names);
    let resolved: Vec<(u32, bool)> = resolved_names
        .iter()
        .filter_map(|(name, f)| graph.get_id(name).map(|id| (id, *f)))
        .collect();
    if resolved.len() != resolved_names.len() {
        return (steps.to_vec(), None);
    }

    let invs_before = steps.iter().filter(|&&(_, f)| !f).count();
    let invs_after_uncollapsed = if collapse_dups {
        let mut uncollapsed = Vec::with_capacity(steps.len());
        flatten_expr(&expr, false, false, &mut uncollapsed);
        uncollapsed.iter().filter(|(_, f)| !*f).count()
    } else {
        resolved.iter().filter(|&&(_, f)| !f).count()
    };
    if invs_before > invs_after_uncollapsed {
        stats.invs_cancelled += invs_before - invs_after_uncollapsed;
    }

    let changed = resolved.len() != steps.len() || resolved != steps;

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
        (steps.to_vec(), None)
    }
}

fn estimate_window_divergence(steps: &[(u32, bool)], is_anchor: &[bool]) -> f64 {
    if steps.is_empty() {
        return 0.02;
    }
    let mut non_anchors = 0usize;
    let mut rev_steps = 0usize;
    for &(seg, fwd) in steps {
        if !is_anchor.get(seg as usize).copied().unwrap_or(false) {
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
    steps: &[(u32, bool)],
    graph: &PackedGraph,
    max_bubble_steps: usize,
    divergence_override: Option<f64>,
    collapse_dups: bool,
    stats: &mut ResolveStats,
) -> (Vec<(u32, bool)>, Vec<String>) {
    let mut out_steps = Vec::with_capacity(steps.len());
    let mut tags = Vec::new();
    let mut i = 0usize;
    while i < steps.len() {
        let mut j = i + 1;
        let mut saw_non_anchor = !graph
            .is_anchor
            .get(steps[i].0 as usize)
            .copied()
            .unwrap_or(false);
        while j < steps.len() && j - i < max_bubble_steps {
            let is_anchor = graph
                .is_anchor
                .get(steps[j].0 as usize)
                .copied()
                .unwrap_or(false);
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
            .unwrap_or_else(|| estimate_window_divergence(&steps[i..j], &graph.is_anchor));
        let (win_steps, win_tag) =
            resolve_window(&steps[i..j], graph, win_div, collapse_dups, stats);
        out_steps.extend(win_steps);
        if let Some(t) = win_tag {
            tags.push(t);
        }
        i = j;
    }
    (out_steps, tags)
}

struct ResolveCtx<'a> {
    graph: &'a PackedGraph,
    cli: &'a CliArgs,
}

fn record_interned_edges(steps: &[(u32, bool)], used: &mut HashSet<(u32, bool, u32, bool)>) {
    for win in steps.windows(2) {
        let (u, uf) = win[0];
        let (v, vf) = win[1];
        used.insert((u, uf, v, vf));
        used.insert((v, !vf, u, !uf));
    }
}

fn emit_resolved_path<W: Write>(
    p: &GfaPath,
    ctx: &ResolveCtx<'_>,
    used_edges: &mut HashSet<(u32, bool, u32, bool)>,
    stats: &mut ResolveStats,
    out: &mut W,
) -> Result<()> {
    stats.paths_seen += 1;
    let raw_steps = parse_path_steps(&p.segment_names);
    let steps: Vec<(u32, bool)> = raw_steps
        .iter()
        .filter_map(|&(s, f)| ctx.graph.get_id(s).map(|id| (id, f)))
        .collect();
    let (resolved, _) = decompose_and_resolve(
        &steps,
        ctx.graph,
        ctx.cli.max_bubble_steps,
        ctx.cli.divergence,
        ctx.cli.collapse_dups,
        stats,
    );
    record_interned_edges(&resolved, used_edges);
    let mut seg_str = String::new();
    for (idx, &(id, fwd)) in resolved.iter().enumerate() {
        if idx > 0 {
            seg_str.push(',');
        }
        seg_str.push_str(&ctx.graph.id_to_seg[id as usize]);
        seg_str.push(if fwd { '+' } else { '-' });
    }
    let overlaps = if resolved == steps {
        p.overlaps.as_str()
    } else {
        "*"
    };
    writeln!(out, "P\t{}\t{}\t{}", p.name, seg_str, overlaps)?;
    Ok(())
}

fn emit_resolved_walk<W: Write>(
    w: &GfaWalk,
    ctx: &ResolveCtx<'_>,
    used_edges: &mut HashSet<(u32, bool, u32, bool)>,
    stats: &mut ResolveStats,
    out: &mut W,
) -> Result<()> {
    stats.walks_seen += 1;
    let raw_steps = parse_walk_steps(&w.walk);
    let steps: Vec<(u32, bool)> = raw_steps
        .iter()
        .filter_map(|&(s, f)| ctx.graph.get_id(s).map(|id| (id, f)))
        .collect();
    let (resolved, tags) = decompose_and_resolve(
        &steps,
        ctx.graph,
        ctx.cli.max_bubble_steps,
        ctx.cli.divergence,
        ctx.cli.collapse_dups,
        stats,
    );
    record_interned_edges(&resolved, used_edges);
    let mut walk_str = String::new();
    let mut resolved_bp = 0usize;
    let mut all_known = true;
    for &(id, fwd) in &resolved {
        walk_str.push(if fwd { '>' } else { '<' });
        walk_str.push_str(&ctx.graph.id_to_seg[id as usize]);
        let slen = ctx.graph.seg_lens[id as usize];
        if slen == 0 {
            all_known = false;
        } else {
            resolved_bp += slen;
        }
    }
    let same_steps = resolved == steps;
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
pub(crate) fn pass2_rewrite<W: Write>(
    gfa_path: &str,
    out: &mut W,
    graph: &PackedGraph,
    cli: &CliArgs,
) -> Result<ResolveStats> {
    let ctx = ResolveCtx { graph, cli };
    let mut stats = ResolveStats::default();
    let mut raw_links: Vec<GfaLink> = Vec::new();
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
                emit_resolved_path(&p, &ctx, &mut used_edges, &mut stats, out)?;
            }
            GfaRecord::Walk(w) => {
                emit_resolved_walk(&w, &ctx, &mut used_edges, &mut stats, out)?;
            }
        }
        Ok(ControlFlow::Continue(()))
    })?;

    let has_paths_or_walks = stats.paths_seen + stats.walks_seen > 0;
    let mut seen_edges: HashSet<(u32, bool, u32, bool)> = HashSet::new();
    for l in raw_links {
        let uf = l.from_orient != b'-';
        let vf = l.to_orient != b'-';
        let (Some(u), Some(v)) = (graph.get_id(&l.from), graph.get_id(&l.to)) else {
            continue;
        };
        let fwd_edge = (u, uf, v, vf);
        if (!has_paths_or_walks || used_edges.contains(&fwd_edge))
            && !seen_edges.contains(&fwd_edge)
        {
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
        .map(|(u, uf, v, vf)| {
            (
                graph.id_to_seg[u as usize].as_str(),
                uf,
                graph.id_to_seg[v as usize].as_str(),
                vf,
                u,
                v,
            )
        })
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
