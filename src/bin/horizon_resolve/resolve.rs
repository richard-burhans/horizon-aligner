use anyhow::Result;
use horizon_aligner::egraph::{
    EGraph, Expr, MdlCostModel, build_seq_expr, extract_min_mdl_cost, saturate_synteny,
};
use horizon_aligner::gfa::{
    GfaPathRef, GfaRecordRef, GfaWalkRef, parse_path_steps, parse_walk_steps,
    stream_gfa_record_refs,
};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::ops::ControlFlow;

use crate::CliArgs;
use crate::snarl::{PackedEdge, PackedGraph, PackedStep};
use crate::vcf::CapturedWalk;

/// Contiguous byte-arena storing GFA segment DNA sequences without per-segment heap allocations.
#[derive(Default)]
pub(crate) struct SegSeqArena {
    data: Vec<u8>,
    /// Start byte offset in `data`; `u64::MAX` marks absent / `*` placeholder segments.
    offsets: Vec<u64>,
    lens: Vec<u32>,
}

impl SegSeqArena {
    pub(crate) fn with_segments(num_segs: usize) -> Self {
        Self {
            data: Vec::new(),
            offsets: vec![u64::MAX; num_segs],
            lens: vec![0; num_segs],
        }
    }

    pub(crate) fn insert(&mut self, id: u32, seq: &[u8]) {
        let idx = id as usize;
        if idx >= self.offsets.len() {
            return;
        }
        let off = self.data.len() as u64;
        self.data.extend_from_slice(seq);
        self.offsets[idx] = off;
        self.lens[idx] = u32::try_from(seq.len()).unwrap_or(u32::MAX);
    }

    #[inline]
    pub(crate) fn has_seq(&self, id: u32) -> bool {
        self.offsets
            .get(id as usize)
            .is_some_and(|&o| o != u64::MAX)
    }

    #[inline]
    pub(crate) fn get(&self, id: u32) -> Option<&[u8]> {
        let &off = self.offsets.get(id as usize)?;
        if off == u64::MAX {
            return None;
        }
        let start = usize::try_from(off).ok()?;
        let len = self.lens[id as usize] as usize;
        self.data.get(start..start.saturating_add(len))
    }

    #[inline]
    pub(crate) fn seq_len(&self, id: u32) -> usize {
        self.lens.get(id as usize).copied().unwrap_or(0) as usize
    }
}

/// Optional pass-2 capture buffer populated only when `--vcf` export is enabled.
#[derive(Default)]
pub(crate) struct VcfCapture {
    pub(crate) seg_seqs: SegSeqArena,
    pub(crate) walks: Vec<CapturedWalk>,
}

impl VcfCapture {
    pub(crate) fn with_segments(num_segs: usize) -> Self {
        Self {
            seg_seqs: SegSeqArena::with_segments(num_segs),
            walks: Vec::new(),
        }
    }
}

#[derive(Default)]
pub(crate) struct ResolveStats {
    pub(crate) paths_seen: usize,
    pub(crate) walks_seen: usize,
    pub(crate) windows_evaluated: usize,
    pub(crate) windows_simplified: usize,
    pub(crate) dups_detected: usize,
    pub(crate) dups_collapsed: usize,
    pub(crate) invs_factored: usize,
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
        Expr::Repeat(inner, k) => {
            let copies = if collapse_dups { 1 } else { *k };
            for _ in 0..copies {
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
        Expr::Repeat(inner, k) => format!("rep({},{k})", format_expr(inner)),
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
        Expr::Repeat(inner, k) => {
            let (d, inv) = count_ops(inner);
            (d + k.saturating_sub(1), inv)
        }
        Expr::Seq(a, b) => {
            let (d1, i1) = count_ops(a);
            let (d2, i2) = count_ops(b);
            (d1 + d2, i1 + i2)
        }
    }
}

/// Build the closed-form optimal factored expression for a duplicate-free window containing
/// consecutive `-`-strand runs (each run of length `L >= 2` factors into `inv(seq(rev(...)))`).
fn format_distinct_factored_expr(steps: &[PackedStep], graph: &PackedGraph) -> String {
    let mut parts: Vec<String> = Vec::new();
    let mut i = 0usize;
    while i < steps.len() {
        if steps[i].fwd() {
            parts.push(graph.seg_name_cow(steps[i].id()).into_owned());
            i += 1;
        } else {
            let start = i;
            while i < steps.len() && !steps[i].fwd() {
                i += 1;
            }
            let run = &steps[start..i];
            let mut acc = graph.seg_name_cow(run[run.len() - 1].id()).into_owned();
            for &step in run[..run.len() - 1].iter().rev() {
                acc = format!("seq({acc},{})", graph.seg_name_cow(step.id()));
            }
            parts.push(format!("inv({acc})"));
        }
    }
    let mut expr = parts[0].clone();
    for p in &parts[1..] {
        expr = format!("seq({expr},{p})");
    }
    expr
}

fn resolve_window(
    steps: &[PackedStep],
    graph: &PackedGraph,
    divergence: f64,
    collapse_dups: bool,
    want_tags: bool,
    stats: &mut ResolveStats,
) -> (Vec<PackedStep>, Option<String>) {
    if steps.len() <= 1 {
        return (steps.to_vec(), None);
    }

    stats.windows_evaluated += 1;

    // Fast path: when all segment IDs in the window are distinct, no duplication/repeat rule
    // (`seq(x,x)=dup(x)`, `rep(x,k)`) can ever fire. Any consecutive `-`-strand runs of length
    // `L >= 2` factor in exact closed form into `inv(seq(rev(...)))`, saving `(L - 1) * inv_op_bits`
    // without materializing a 4,096-node e-graph per window on inverted whole-genome walks.
    let has_dup_seg =
        (0..steps.len()).any(|a| ((a + 1)..steps.len()).any(|b| steps[a].id() == steps[b].id()));
    if !has_dup_seg {
        let mut invs_factored = 0usize;
        let mut run_len = 0usize;
        for &step in steps {
            if step.fwd() {
                if run_len >= 2 {
                    invs_factored += run_len - 1;
                }
                run_len = 0;
            } else {
                run_len += 1;
            }
        }
        if run_len >= 2 {
            invs_factored += run_len - 1;
        }
        if invs_factored > 0 {
            stats.windows_simplified += 1;
            stats.invs_factored += invs_factored;
            #[allow(clippy::cast_precision_loss)]
            {
                stats.total_mdl_bits_saved +=
                    (invs_factored as f64) * MdlCostModel::default().inv_op_bits;
            }
            let tag = want_tags.then(|| format_distinct_factored_expr(steps, graph));
            return (steps.to_vec(), tag);
        }
        return (steps.to_vec(), None);
    }

    let names: Vec<std::borrow::Cow<'_, str>> =
        steps.iter().map(|&s| graph.seg_name_cow(s.id())).collect();
    let blocks: Vec<(&str, i8)> = steps
        .iter()
        .zip(&names)
        .map(|(&s, name)| (name.as_ref(), if s.fwd() { 1 } else { -1 }))
        .collect();

    let mut eg = EGraph::new();
    let Some(root) = build_seq_expr(&mut eg, &blocks) else {
        return (steps.to_vec(), None);
    };

    let mut model = MdlCostModel::from_egraph(&eg);
    for (&step, name) in steps.iter().zip(&names) {
        let len_bp = graph.seg_len(step.id()).max(1);
        model
            .block_stats
            .entry(name.as_ref().to_owned())
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
    let expected_len = resolved_names.len();
    let resolved: Vec<PackedStep> = resolved_names
        .into_iter()
        .filter_map(|(name, f)| graph.get_id(name).map(|id| PackedStep::new(id, f)))
        .collect();
    if resolved.len() != expected_len || (resolved.is_empty() && !steps.is_empty()) {
        return (steps.to_vec(), None);
    }

    // Count factored inversion operators (`inv(seq(b, a))` replacing separate `inv(a), inv(b)`),
    // expanding `dup(x)` by 2x so reductions from duplication collapse are not double-counted.
    let invs_before = steps.iter().filter(|&&s| !s.fwd()).count();
    let invs_factored = invs_before.saturating_sub(count_expanded_invs(&expr));
    stats.invs_factored += invs_factored;

    let structurally_simplified = dups > 0 || invs_factored > 0 || raw_cost > opt_cost + 1e-6;
    let steps_changed = resolved.len() != steps.len() || resolved != steps;

    if structurally_simplified || steps_changed {
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

fn count_expanded_invs(expr: &Expr) -> usize {
    match expr {
        Expr::Block(_) => 0,
        Expr::Inv(a) => 1 + count_expanded_invs(a),
        Expr::Dup(a) => 2 * count_expanded_invs(a),
        Expr::Repeat(a, k) => k.saturating_mul(count_expanded_invs(a)),
        Expr::Seq(a, b) => count_expanded_invs(a) + count_expanded_invs(b),
    }
}

fn estimate_window_divergence(steps: &[PackedStep], is_anchor: &[bool]) -> f64 {
    if steps.is_empty() {
        return 0.02;
    }
    let mut non_anchors = 0usize;
    let mut rev_steps = 0usize;
    for &step in steps {
        if !is_anchor.get(step.id() as usize).copied().unwrap_or(false) {
            non_anchors += 1;
        }
        if !step.fwd() {
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
    steps: &[PackedStep],
    graph: &PackedGraph,
    max_bubble_steps: usize,
    divergence_override: Option<f64>,
    collapse_dups: bool,
    want_tags: bool,
    stats: &mut ResolveStats,
) -> (Vec<PackedStep>, Vec<String>) {
    let mut out_steps = Vec::with_capacity(steps.len());
    let mut tags = Vec::new();
    let mut i = 0usize;
    while i < steps.len() {
        let mut j = i + 1;
        let mut saw_non_anchor = !graph
            .is_anchor
            .get(steps[i].id() as usize)
            .copied()
            .unwrap_or(false);
        while j < steps.len() && j - i < max_bubble_steps {
            let is_anchor = graph
                .is_anchor
                .get(steps[j].id() as usize)
                .copied()
                .unwrap_or(false);
            let is_tandem = steps[j].id() == steps[j - 1].id();
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
        let (win_steps, win_tag) = resolve_window(
            &steps[i..j],
            graph,
            win_div,
            collapse_dups,
            want_tags,
            stats,
        );
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

fn record_interned_edges(steps: &[PackedStep], used: &mut HashSet<PackedEdge>) {
    for win in steps.windows(2) {
        used.insert(PackedEdge::directed(win[0], win[1]).canonical());
    }
}

fn emit_resolved_path<W: Write>(
    p: GfaPathRef<'_>,
    ctx: &ResolveCtx<'_>,
    used_edges: &mut HashSet<PackedEdge>,
    stats: &mut ResolveStats,
    vcf_capture: Option<&mut VcfCapture>,
    out: &mut W,
) -> Result<()> {
    stats.paths_seen += 1;
    let raw_steps = parse_path_steps(p.segment_names);
    let steps: Vec<PackedStep> = raw_steps
        .iter()
        .filter_map(|&(s, f)| ctx.graph.get_id(s).map(|id| PackedStep::new(id, f)))
        .collect();
    let (resolved, _) = decompose_and_resolve(
        &steps,
        ctx.graph,
        ctx.cli.max_bubble_steps,
        ctx.cli.divergence,
        ctx.cli.collapse_dups,
        false,
        stats,
    );
    record_interned_edges(&resolved, used_edges);
    let same_steps = resolved == steps;
    if let Some(cap) = vcf_capture {
        // Capture the original uncollapsed path steps so VCF POS/REF coordinates stay exact on
        // the linear reference assembly even when `--collapse-dups` simplifies the rewritten GFA.
        let mut total_bp = 0usize;
        for &step in &steps {
            total_bp = total_bp.saturating_add(ctx.graph.seg_len(step.id()));
        }
        for win in steps.windows(2) {
            if let Some(&ov) = ctx
                .graph
                .link_overlaps
                .get(&PackedEdge::directed(win[0], win[1]))
            {
                total_bp = total_bp.saturating_sub(ov as usize);
            }
        }
        let (sample, hap, seq_id) = parse_path_sample_fields(p.name);
        cap.walks.push(CapturedWalk {
            sample,
            haplotype: hap,
            seq_id,
            start: 0,
            end: total_bp,
            steps,
        });
    }
    let mut seg_str = String::with_capacity(resolved.len().saturating_mul(6));
    for (idx, &step) in resolved.iter().enumerate() {
        if idx > 0 {
            seg_str.push(',');
        }
        ctx.graph.push_name(step.id(), &mut seg_str);
        seg_str.push(if step.fwd() { '+' } else { '-' });
    }
    let overlaps = if same_steps { p.overlaps } else { "*" };
    writeln!(out, "P\t{}\t{}\t{}", p.name, seg_str, overlaps)?;
    Ok(())
}

/// Split PanSN-style `sample#hap#contig` GFA `P`-line names (falling back to `(name, "0", name)`).
fn parse_path_sample_fields(name: &str) -> (String, String, String) {
    let (sample, hap, contig) = horizon_aligner::gfa::parse_path_sample_fields(name);
    (sample.to_owned(), hap.to_owned(), contig.to_owned())
}

fn emit_resolved_walk<W: Write>(
    w: GfaWalkRef<'_>,
    ctx: &ResolveCtx<'_>,
    used_edges: &mut HashSet<PackedEdge>,
    stats: &mut ResolveStats,
    vcf_capture: Option<&mut VcfCapture>,
    out: &mut W,
) -> Result<()> {
    stats.walks_seen += 1;
    let raw_steps = parse_walk_steps(w.walk);
    let steps: Vec<PackedStep> = raw_steps
        .iter()
        .filter_map(|&(s, f)| ctx.graph.get_id(s).map(|id| PackedStep::new(id, f)))
        .collect();
    let (resolved, tags) = decompose_and_resolve(
        &steps,
        ctx.graph,
        ctx.cli.max_bubble_steps,
        ctx.cli.divergence,
        ctx.cli.collapse_dups,
        ctx.cli.tag_walks,
        stats,
    );
    record_interned_edges(&resolved, used_edges);
    let mut walk_str = String::with_capacity(resolved.len().saturating_mul(6));
    let mut resolved_bp = 0usize;
    let mut all_known = true;
    for &step in &resolved {
        walk_str.push(if step.fwd() { '>' } else { '<' });
        ctx.graph.push_name(step.id(), &mut walk_str);
        let slen = ctx.graph.seg_len(step.id());
        if slen == 0 {
            all_known = false;
        } else {
            resolved_bp = resolved_bp.saturating_add(slen);
        }
    }
    for win in resolved.windows(2) {
        if let Some(&ov) = ctx
            .graph
            .link_overlaps
            .get(&PackedEdge::directed(win[0], win[1]))
        {
            resolved_bp = resolved_bp.saturating_sub(ov as usize);
        }
    }
    let same_steps = resolved == steps;
    let end_pos = if !same_steps && all_known && resolved_bp > 0 {
        w.start.saturating_add(resolved_bp)
    } else {
        w.end
    };
    if let Some(cap) = vcf_capture {
        // Normalize PanSN `sample#hap` when encoded inside `w.sample` (e.g. `HG002#1` with hap `0`)
        // so diploid haplotypes merge into a single phased VCF sample column (`HG002`).
        let (norm_sample, norm_hap) = if let Some((s, h)) = w.sample.split_once('#')
            && !s.is_empty()
            && !h.is_empty()
            && h.bytes().all(|b| b.is_ascii_digit())
            && w.haplotype == "0"
        {
            (s.to_owned(), h.to_owned())
        } else {
            (w.sample.to_owned(), w.haplotype.to_owned())
        };
        // Preserve original uncollapsed walk steps/bounds so VCF coordinates never drift when
        // `--collapse-dups` collapses tandem repeats on the reference walk.
        cap.walks.push(CapturedWalk {
            sample: norm_sample,
            haplotype: norm_hap,
            seq_id: w.seq_id.to_owned(),
            start: w.start,
            end: w.end,
            steps,
        });
    }
    let extra_tags: String = if let Some(raw_tags) = w.tags {
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

/// Write a single GFA `L` line without allocating intermediate `String` names.
fn write_gfa_link<W: Write>(
    out: &mut W,
    graph: &PackedGraph,
    edge: PackedEdge,
    overlap: &str,
) -> Result<()> {
    let u = edge.src();
    let v = edge.dst();
    out.write_all(b"L\t")?;
    graph.write_name(u.id(), out)?;
    out.write_all(if u.fwd() { b"\t+\t" } else { b"\t-\t" })?;
    graph.write_name(v.id(), out)?;
    out.write_all(if v.fwd() { b"\t+\t" } else { b"\t-\t" })?;
    out.write_all(overlap.as_bytes())?;
    out.write_all(b"\n")?;
    Ok(())
}

type MissingEdge<'a> = (
    std::borrow::Cow<'a, str>,
    bool,
    std::borrow::Cow<'a, str>,
    bool,
    PackedEdge,
);

/// Pass 2: rewrite `P`/`W` records using zero-copy GFA streaming and emit only `L` links
/// still traversed by resolved walks (stored compactly as 8-byte `PackedEdge` entries).
#[allow(clippy::too_many_lines)]
pub(crate) fn pass2_rewrite<W: Write>(
    gfa_path: impl AsRef<std::path::Path>,
    out: &mut W,
    graph: &PackedGraph,
    cli: &CliArgs,
    mut vcf_capture: Option<&mut VcfCapture>,
) -> Result<ResolveStats> {
    let gfa_path = gfa_path.as_ref();
    let ctx = ResolveCtx { graph, cli };
    let mut stats = ResolveStats::default();
    let mut raw_links: Vec<PackedEdge> = Vec::new();
    let mut non_zero_overlaps: HashMap<usize, String> = HashMap::new();
    let mut used_edges: HashSet<PackedEdge> = HashSet::new();

    writeln!(out, "H\tVN:Z:1.1")?;
    stream_gfa_record_refs(gfa_path, |rec| {
        match rec {
            GfaRecordRef::Segment(s) => {
                if let Some(cap) = vcf_capture.as_deref_mut()
                    && let Some(id) = graph.get_id(s.name)
                    && !s.sequence.is_empty()
                    && s.sequence != b"*"
                {
                    cap.seg_seqs.insert(id, s.sequence);
                }
                out.write_all(b"S\t")?;
                out.write_all(s.name.as_bytes())?;
                out.write_all(b"\t")?;
                out.write_all(s.sequence)?;
                if let Some(tags) = s.tags {
                    out.write_all(b"\t")?;
                    out.write_all(tags.as_bytes())?;
                }
                out.write_all(b"\n")?;
            }
            GfaRecordRef::Link(l) => {
                if let (Some(u), Some(v)) = (graph.get_id(l.from), graph.get_id(l.to)) {
                    let uf = l.from_orient != b'-';
                    let vf = l.to_orient != b'-';
                    let link_idx = raw_links.len();
                    raw_links.push(PackedEdge::directed(
                        PackedStep::new(u, uf),
                        PackedStep::new(v, vf),
                    ));
                    if l.overlap != "0M" {
                        non_zero_overlaps.insert(link_idx, l.overlap.to_owned());
                    }
                }
            }
            GfaRecordRef::Path(p) => {
                emit_resolved_path(
                    p,
                    &ctx,
                    &mut used_edges,
                    &mut stats,
                    vcf_capture.as_deref_mut(),
                    out,
                )?;
            }
            GfaRecordRef::Walk(w) => {
                emit_resolved_walk(
                    w,
                    &ctx,
                    &mut used_edges,
                    &mut stats,
                    vcf_capture.as_deref_mut(),
                    out,
                )?;
            }
        }
        Ok(ControlFlow::Continue(()))
    })?;

    let has_paths_or_walks = stats.paths_seen + stats.walks_seen > 0;
    let mut seen_edges: HashSet<PackedEdge> = HashSet::with_capacity(used_edges.len());
    for (link_idx, edge) in raw_links.into_iter().enumerate() {
        let canon = edge.canonical();
        if (!has_paths_or_walks || used_edges.contains(&canon)) && seen_edges.insert(canon) {
            let ov = non_zero_overlaps
                .get(&link_idx)
                .map_or("0M", String::as_str);
            write_gfa_link(out, graph, edge, ov)?;
        }
    }
    drop(non_zero_overlaps);

    let mut missing: Vec<MissingEdge<'_>> = used_edges
        .into_iter()
        .filter(|e| !seen_edges.contains(e))
        .map(|edge| {
            // Pick the lexicographically smaller directed orientation to stay deterministic
            // and bit-identical with bidirectional missing-edge sorting.
            let rc = edge.rev();
            let u_a = graph.seg_name_cow(edge.src().id());
            let v_a = graph.seg_name_cow(edge.dst().id());
            let key_a = (
                u_a.as_ref(),
                edge.src().fwd(),
                v_a.as_ref(),
                edge.dst().fwd(),
            );
            let key_b = (v_a.as_ref(), rc.src().fwd(), u_a.as_ref(), rc.dst().fwd());
            if key_a <= key_b {
                (u_a, edge.src().fwd(), v_a, edge.dst().fwd(), edge)
            } else {
                (v_a, rc.src().fwd(), u_a, rc.dst().fwd(), rc)
            }
        })
        .collect();
    missing.sort_unstable_by(|a, b| {
        a.0.as_ref()
            .cmp(b.0.as_ref())
            .then_with(|| a.1.cmp(&b.1))
            .then_with(|| a.2.as_ref().cmp(b.2.as_ref()))
            .then_with(|| a.3.cmp(&b.3))
    });
    for (_, _, _, _, edge) in missing {
        write_gfa_link(out, graph, edge, "0M")?;
    }

    out.flush()?;
    Ok(stats)
}
