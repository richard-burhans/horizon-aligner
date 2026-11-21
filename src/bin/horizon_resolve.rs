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

/// Compact `u32`-interned biedged pangenome graph in CSR format with 2-edge-connected
/// snarl/bubble boundary annotations.
struct PackedGraph {
    seg_to_id: HashMap<String, u32>,
    id_to_seg: Vec<String>,
    seg_lens: Vec<usize>,
    is_anchor: Vec<bool>,
    num_bridges: usize,
    num_components: usize,
}

impl PackedGraph {
    fn get_id(&self, name: &str) -> Option<u32> {
        self.seg_to_id.get(name).copied()
    }

    fn intern(&mut self, name: &str, len: usize) -> u32 {
        if let Some(&id) = self.seg_to_id.get(name) {
            if len > 0 && self.seg_lens[id as usize] == 0 {
                self.seg_lens[id as usize] = len;
            }
            return id;
        }
        let id = u32::try_from(self.id_to_seg.len()).expect("segment count fits u32");
        let owned = name.to_owned();
        self.seg_to_id.insert(owned.clone(), id);
        self.id_to_seg.push(owned);
        self.seg_lens.push(len);
        id
    }
}

#[inline]
const fn exit_side(seg: u32, fwd: bool) -> u32 {
    2 * seg + if fwd { 1 } else { 0 }
}

#[inline]
const fn entry_side(seg: u32, fwd: bool) -> u32 {
    2 * seg + if fwd { 0 } else { 1 }
}

fn tarjan_bridges(
    num_sides: usize,
    total_edges: usize,
    head: &[usize],
    csr_to: &[u32],
    csr_eid: &[usize],
) -> (Vec<bool>, usize) {
    let mut tin = vec![0u32; num_sides];
    let mut low = vec![0u32; num_sides];
    let mut is_bridge = vec![false; total_edges];
    let mut timer = 0u32;
    let mut num_bridges = 0usize;
    let mut stack: Vec<(usize, usize, usize)> = Vec::new();

    for root in 0..num_sides {
        if tin[root] != 0 {
            continue;
        }
        timer += 1;
        tin[root] = timer;
        low[root] = timer;
        stack.push((root, usize::MAX, head[root]));

        while let Some((u, peid, idx)) = stack.last_mut() {
            let u_node = *u;
            let parent_eid = *peid;
            let end = head[u_node + 1];
            if *idx < end {
                let cur = *idx;
                *idx += 1;
                let v = csr_to[cur] as usize;
                let eid = csr_eid[cur];
                if eid == parent_eid {
                    continue;
                }
                if tin[v] != 0 {
                    low[u_node] = low[u_node].min(tin[v]);
                } else {
                    timer += 1;
                    tin[v] = timer;
                    low[v] = timer;
                    stack.push((v, eid, head[v]));
                }
            } else {
                let finished_u = u_node;
                let finished_peid = parent_eid;
                stack.pop();
                if let Some(&(p, _, _)) = stack.last() {
                    low[p] = low[p].min(low[finished_u]);
                    if low[finished_u] > tin[p] && !is_bridge[finished_peid] {
                        is_bridge[finished_peid] = true;
                        num_bridges += 1;
                    }
                }
            }
        }
    }
    (is_bridge, num_bridges)
}

/// Build CSR adjacency for the biedged graph (black segment edges + gray link edges)
/// and run iterative Tarjan bridge-finding to identify snarl boundary anchors.
fn compute_biedged_snarls(
    num_segs: usize,
    seg_lens: &[usize],
    seg_visits: &[u32],
    gray_edges: &[(u32, u32)],
    min_anchor_bp: usize,
) -> (Vec<bool>, usize, usize) {
    let num_sides = num_segs * 2;
    if num_sides == 0 {
        return (Vec::new(), 0, 0);
    }

    let total_edges = num_segs + gray_edges.len();
    let mut deg = vec![1u32; num_sides];
    for &(u, v) in gray_edges {
        deg[u as usize] += 1;
        deg[v as usize] += 1;
    }

    let mut head = vec![0usize; num_sides + 1];
    for s in 0..num_sides {
        head[s + 1] = head[s] + deg[s] as usize;
    }
    let mut csr_to = vec![0u32; head[num_sides]];
    let mut csr_eid = vec![0usize; head[num_sides]];
    let mut cursor = head[..num_sides].to_vec();

    for v in 0..num_segs {
        let s0 = 2 * v;
        let s1 = s0 + 1;
        let p0 = cursor[s0];
        cursor[s0] += 1;
        csr_to[p0] = u32::try_from(s1).expect("side fits u32");
        csr_eid[p0] = v;

        let p1 = cursor[s1];
        cursor[s1] += 1;
        csr_to[p1] = u32::try_from(s0).expect("side fits u32");
        csr_eid[p1] = v;
    }

    for (idx, &(u, v)) in gray_edges.iter().enumerate() {
        let eid = num_segs + idx;
        let pu = cursor[u as usize];
        cursor[u as usize] += 1;
        csr_to[pu] = v;
        csr_eid[pu] = eid;

        let pv = cursor[v as usize];
        cursor[v as usize] += 1;
        csr_to[pv] = u;
        csr_eid[pv] = eid;
    }

    let (is_bridge, num_bridges) = tarjan_bridges(num_sides, total_edges, &head, &csr_to, &csr_eid);

    // Count 2-edge-connected components by BFS ignoring bridge edges.
    let mut comp = vec![u32::MAX; num_sides];
    let mut num_components = 0usize;
    let mut queue = Vec::new();
    for start in 0..num_sides {
        if comp[start] != u32::MAX {
            continue;
        }
        let cid = u32::try_from(num_components).expect("component fits u32");
        num_components += 1;
        comp[start] = cid;
        queue.push(start);
        let mut q_head = 0usize;
        while q_head < queue.len() {
            let u = queue[q_head];
            q_head += 1;
            for cur in head[u]..head[u + 1] {
                let eid = csr_eid[cur];
                if is_bridge[eid] {
                    continue;
                }
                let v = csr_to[cur] as usize;
                if comp[v] == u32::MAX {
                    comp[v] = cid;
                    queue.push(v);
                }
            }
        }
        queue.clear();
    }

    // A segment is a snarl boundary anchor if:
    // - its black edge is a bridge AND (shared by >=2 paths, or >= min_anchor_bp, or a graph tip), OR
    // - it is shared by >=2 paths AND >= min_anchor_bp.
    let mut is_anchor = vec![false; num_segs];
    for v in 0..num_segs {
        let black_bridge = is_bridge[v];
        let is_tip = deg[2 * v] == 1 || deg[2 * v + 1] == 1;
        let multi_visit = seg_visits[v] >= 2;
        let long_seg = seg_lens[v] >= min_anchor_bp;
        if (black_bridge && (multi_visit || long_seg || is_tip)) || (multi_visit && long_seg) {
            is_anchor[v] = true;
        }
    }

    (is_anchor, num_bridges, num_components)
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

fn record_walk_steps(
    steps: &[(&str, bool)],
    graph: &mut PackedGraph,
    seg_visits: &mut Vec<u32>,
    visit_stamp: &mut Vec<u32>,
    stamp: u32,
    gray_edges: &mut Vec<(u32, u32)>,
) {
    let mut prev: Option<(u32, bool)> = None;
    for &(seg, fwd) in steps {
        let id = graph.intern(seg, 0);
        let idx = id as usize;
        if idx >= seg_visits.len() {
            seg_visits.resize(idx + 1, 0);
            visit_stamp.resize(idx + 1, 0);
        }
        if visit_stamp[idx] != stamp {
            visit_stamp[idx] = stamp;
            seg_visits[idx] += 1;
        }
        if let Some((pu, pf)) = prev {
            let su = exit_side(pu, pf);
            let sv = entry_side(id, fwd);
            gray_edges.push(if su <= sv { (su, sv) } else { (sv, su) });
        }
        prev = Some((id, fwd));
    }
}

fn pass1_index(gfa_path: &str, min_anchor_bp: usize) -> Result<PackedGraph> {
    let mut graph = PackedGraph {
        seg_to_id: HashMap::new(),
        id_to_seg: Vec::new(),
        seg_lens: Vec::new(),
        is_anchor: Vec::new(),
        num_bridges: 0,
        num_components: 0,
    };
    let mut seg_visits: Vec<u32> = Vec::new();
    let mut visit_stamp: Vec<u32> = Vec::new();
    let mut cur_stamp = 0u32;
    let mut gray_edges: Vec<(u32, u32)> = Vec::new();

    stream_gfa_records(gfa_path, |rec| {
        match rec {
            GfaRecord::Segment(s) => {
                let len = s.seq_len();
                let id = graph.intern(&s.name, len);
                let idx = id as usize;
                if idx >= seg_visits.len() {
                    seg_visits.resize(idx + 1, 0);
                    visit_stamp.resize(idx + 1, 0);
                }
            }
            GfaRecord::Link(l) => {
                let u = graph.intern(&l.from, 0);
                let v = graph.intern(&l.to, 0);
                let max_idx = (u.max(v) as usize) + 1;
                if max_idx > seg_visits.len() {
                    seg_visits.resize(max_idx, 0);
                    visit_stamp.resize(max_idx, 0);
                }
                let su = exit_side(u, l.from_orient != b'-');
                let sv = entry_side(v, l.to_orient != b'-');
                gray_edges.push(if su <= sv { (su, sv) } else { (sv, su) });
            }
            GfaRecord::Path(p) => {
                cur_stamp += 1;
                let steps = parse_path_steps(&p.segment_names);
                record_walk_steps(
                    &steps,
                    &mut graph,
                    &mut seg_visits,
                    &mut visit_stamp,
                    cur_stamp,
                    &mut gray_edges,
                );
            }
            GfaRecord::Walk(w) => {
                cur_stamp += 1;
                let steps = parse_walk_steps(&w.walk);
                record_walk_steps(
                    &steps,
                    &mut graph,
                    &mut seg_visits,
                    &mut visit_stamp,
                    cur_stamp,
                    &mut gray_edges,
                );
            }
        }
        Ok(ControlFlow::Continue(()))
    })?;

    gray_edges.sort_unstable();
    gray_edges.dedup();
    let (is_anchor, num_bridges, num_components) = compute_biedged_snarls(
        graph.id_to_seg.len(),
        &graph.seg_lens,
        &seg_visits,
        &gray_edges,
        min_anchor_bp,
    );
    graph.is_anchor = is_anchor;
    graph.num_bridges = num_bridges;
    graph.num_components = num_components;
    Ok(graph)
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
    p: &horizon_aligner::gfa::GfaPath,
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
    w: &horizon_aligner::gfa::GfaWalk,
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
fn pass2_rewrite<W: Write>(
    gfa_path: &str,
    out: &mut W,
    graph: &PackedGraph,
    cli: &CliArgs,
) -> Result<ResolveStats> {
    let ctx = ResolveCtx { graph, cli };
    let mut stats = ResolveStats::default();
    let mut raw_links: Vec<horizon_aligner::gfa::GfaLink> = Vec::new();
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

fn main() -> Result<()> {
    let Some(cli) = parse_args()? else {
        return Ok(());
    };
    let t0 = Instant::now();
    let spooled = SpoolableGfaInput::new(&cli.gfa_in)?;
    let effective_gfa = spooled.path();
    let graph = pass1_index(effective_gfa, cli.min_anchor_bp)?;
    let anchor_count = graph.is_anchor.iter().filter(|&&a| a).count();

    let mut writer = open_output_writer(cli.out_path.as_deref())?;
    let stats = pass2_rewrite(effective_gfa, &mut writer, &graph, &cli)?;
    writer.flush()?;

    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    eprintln!(
        "segments={} anchors={} bridges={} components={} paths={} walks={} windows={} simplified={} dups_detected={} dups_collapsed={} invs_cancelled={} mdl_bits_saved={:.1} elapsed_s={:.3} peak_rss_mb={:.1}",
        graph.id_to_seg.len(),
        anchor_count,
        graph.num_bridges,
        graph.num_components,
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
        let graph = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;
        assert!(graph.num_bridges >= 2);
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&gfa_out)?);
        let stats = pass2_rewrite(&cli.gfa_in, &mut writer, &graph, &cli)?;
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
