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
use std::io::{BufRead, Write};
use std::ops::ControlFlow;
use std::time::Instant;

struct CliArgs {
    gfa_in: String,
    out_path: Option<String>,
    collapse_dups: bool,
    tag_walks: bool,
    max_bubble_steps: usize,
    min_anchor_bp: usize,
    divergence: f64,
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
    let mut divergence = 0.02_f64;

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
                divergence = args
                    .get(i)
                    .context("missing value for --divergence")?
                    .parse()
                    .context("invalid --divergence")?;
            }
            other if (!other.starts_with('-') || other == "-") && gfa_in.is_none() => {
                gfa_in = Some(other.to_string());
            }
            other => bail!("unknown argument: {other}"),
        }
        i += 1;
    }

    let gfa_in = gfa_in.context("missing <input.gfa[.gz]>")?;
    anyhow::ensure!(
        divergence.is_finite() && divergence > 0.0,
        "--divergence must be a finite positive number, got {divergence}"
    );
    Ok(Some(CliArgs {
        gfa_in,
        out_path,
        collapse_dups,
        tag_walks,
        max_bubble_steps: max_bubble_steps.max(2),
        min_anchor_bp,
        divergence: divergence.clamp(1e-6, 0.75),
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

    let raw_cost = extract_min_mdl_cost(&mut eg, root, &model).map_or(0.0, |(c, _)| c);
    saturate_synteny(&mut eg, 6, 4096);
    let Some((opt_cost, expr)) = extract_min_mdl_cost(&mut eg, root, &model) else {
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

fn decompose_and_resolve(
    steps: &[(&str, bool)],
    anchors: &HashSet<String>,
    seg_lens: &HashMap<String, usize>,
    max_bubble_steps: usize,
    divergence: f64,
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
        let (win_steps, win_tag) =
            resolve_window(&steps[i..j], seg_lens, divergence, collapse_dups, stats);
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

/// Pass 2: rewrite `P`/`W` records and emit only `L` links still traversed by resolved walks.
///
/// Note: buffers `L` lines in `raw_l_lines` so unused links can be pruned after all `P`/`W` lines
/// are seen; memory scales with the number of `L` records in the input GFA.
#[allow(clippy::too_many_lines)]
fn pass2_rewrite<W: Write>(
    gfa_path: &str,
    out: &mut W,
    seg_lens: &HashMap<String, usize>,
    anchors: &HashSet<String>,
    cli: &CliArgs,
) -> Result<ResolveStats> {
    let mut reader = horizon_aligner::gfa::open_line_reader(gfa_path)?;
    let mut stats = ResolveStats::default();
    let mut line = String::new();
    let mut raw_l_lines: Vec<String> = Vec::new();
    let mut used_edges: HashSet<(String, bool, String, bool)> = HashSet::new();

    let record_edges = |steps: &[(String, bool)],
                        used: &mut HashSet<(String, bool, String, bool)>| {
        for win in steps.windows(2) {
            let (ref u, uf) = win[0];
            let (ref v, vf) = win[1];
            used.insert((u.clone(), uf, v.clone(), vf));
            used.insert((v.clone(), !vf, u.clone(), !uf));
        }
    };

    while reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.starts_with("L\t") {
            raw_l_lines.push(trimmed.to_string());
            line.clear();
            continue;
        } else if trimmed.starts_with("P\t") {
            let cols: Vec<&str> = trimmed.split('\t').collect();
            if cols.len() >= 4 {
                stats.paths_seen += 1;
                let steps = parse_path_steps(cols[2]);
                let (resolved, _) = decompose_and_resolve(
                    &steps,
                    anchors,
                    seg_lens,
                    cli.max_bubble_steps,
                    cli.divergence,
                    cli.collapse_dups,
                    &mut stats,
                );
                record_edges(&resolved, &mut used_edges);
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
                    cols[3].to_string()
                } else {
                    "*".to_string()
                };
                if cols.len() > 4 {
                    writeln!(
                        out,
                        "P\t{}\t{}\t{}\t{}",
                        cols[1],
                        seg_str,
                        overlaps,
                        cols[4..].join("\t")
                    )?;
                } else {
                    writeln!(out, "P\t{}\t{}\t{}", cols[1], seg_str, overlaps)?;
                }
                line.clear();
                continue;
            }
        } else if trimmed.starts_with("W\t") {
            let cols: Vec<&str> = trimmed.split('\t').collect();
            if cols.len() >= 7 {
                stats.walks_seen += 1;
                let steps = parse_walk_steps(cols[6]);
                let (resolved, tags) = decompose_and_resolve(
                    &steps,
                    anchors,
                    seg_lens,
                    cli.max_bubble_steps,
                    cli.divergence,
                    cli.collapse_dups,
                    &mut stats,
                );
                record_edges(&resolved, &mut used_edges);
                let mut walk_str = String::new();
                let mut resolved_bp = 0usize;
                let mut all_known = true;
                for (s, fwd) in &resolved {
                    walk_str.push(if *fwd { '>' } else { '<' });
                    walk_str.push_str(s);
                    let slen = seg_lens.get(s.as_str()).copied().unwrap_or(0);
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
                let end_str = if !same_steps
                    && all_known
                    && resolved_bp > 0
                    && let Ok(start_pos) = cols[4].parse::<usize>()
                {
                    start_pos.saturating_add(resolved_bp).to_string()
                } else {
                    cols[5].to_string()
                };
                let extra_tags: String = if cols.len() > 7 {
                    let kept: Vec<&str> = cols[7..]
                        .iter()
                        .copied()
                        .filter(|t| !(cli.tag_walks && !tags.is_empty() && t.starts_with("EG:Z:")))
                        .collect();
                    if kept.is_empty() {
                        String::new()
                    } else {
                        format!("\t{}", kept.join("\t"))
                    }
                } else {
                    String::new()
                };
                if cli.tag_walks && !tags.is_empty() {
                    writeln!(
                        out,
                        "W\t{}\t{}\t{}\t{}\t{}\t{}{}\tEG:Z:{}",
                        cols[1],
                        cols[2],
                        cols[3],
                        cols[4],
                        end_str,
                        walk_str,
                        extra_tags,
                        tags.join(";")
                    )?;
                } else {
                    writeln!(
                        out,
                        "W\t{}\t{}\t{}\t{}\t{}\t{}{}",
                        cols[1], cols[2], cols[3], cols[4], end_str, walk_str, extra_tags
                    )?;
                }
                line.clear();
                continue;
            }
        }
        writeln!(out, "{trimmed}")?;
        line.clear();
    }

    let has_paths_or_walks = stats.paths_seen + stats.walks_seen > 0;
    let mut seen_edges: HashSet<(String, bool, String, bool)> = HashSet::new();
    for l_line in &raw_l_lines {
        let cols: Vec<&str> = l_line.split('\t').collect();
        if cols.len() >= 6 {
            let uf = cols[2] != "-";
            let vf = cols[4] != "-";
            let fwd_edge = (cols[1].to_string(), uf, cols[3].to_string(), vf);
            if !has_paths_or_walks || used_edges.contains(&fwd_edge) {
                writeln!(out, "{l_line}")?;
                seen_edges.insert((cols[3].to_string(), !vf, cols[1].to_string(), !uf));
                seen_edges.insert(fwd_edge);
            }
        } else {
            writeln!(out, "{l_line}")?;
        }
    }
    let mut missing: Vec<(String, bool, String, bool)> = used_edges
        .into_iter()
        .filter(|e| !seen_edges.contains(e))
        .collect();
    missing.sort();
    for (u, uf, v, vf) in missing {
        let rc = (v.clone(), !vf, u.clone(), !uf);
        if !seen_edges.contains(&(u.clone(), uf, v.clone(), vf)) && !seen_edges.contains(&rc) {
            writeln!(
                out,
                "L\t{}\t{}\t{}\t{}\t0M",
                u,
                if uf { '+' } else { '-' },
                v,
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
            divergence: 0.01,
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
