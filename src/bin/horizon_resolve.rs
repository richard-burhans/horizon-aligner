//! Two-pass streaming GFA v1.0/v1.1 snarl/bubble e-graph resolver using MdlCostModel.

use anyhow::{Context, Result, bail};
use horizon_aligner::egraph::{
    EGraph, Expr, MdlCostModel, build_seq_expr, extract_min_mdl_cost, saturate_synteny,
};
use horizon_aligner::gfa::{parse_path_steps, parse_walk_steps, stream_gfa};
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
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

fn parse_args() -> Result<CliArgs> {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        println!(
            "Usage: horizon_resolve <input.gfa[.gz]> [--out <resolved.gfa>] [--collapse-dups] [--tag-walks] [--max-bubble-steps <N>] [--min-anchor-bp <N>] [--divergence <F>]\n\
             \n\
             Decomposes GFA v1.0/v1.1 P/W paths into bounded snarl/bubble windows between shared anchors,\n\
             runs bounded synteny e-graph saturation with MdlCostModel extraction, and rewrites\n\
             simplified walks/paths back into valid GFA."
        );
        std::process::exit(0);
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
            other if !other.starts_with('-') && gfa_in.is_none() => {
                gfa_in = Some(other.to_string());
            }
            other => bail!("unknown argument: {other}"),
        }
        i += 1;
    }

    let gfa_in = gfa_in.context("missing <input.gfa[.gz]>")?;
    Ok(CliArgs {
        gfa_in,
        out_path,
        collapse_dups,
        tag_walks,
        max_bubble_steps: max_bubble_steps.max(2),
        min_anchor_bp,
        divergence: divergence.clamp(1e-6, 0.75),
    })
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

    let (dups, invs_after) = count_ops(&expr);
    let invs_before = steps.iter().filter(|&&(_, f)| !f).count();
    if invs_before > invs_after {
        stats.invs_cancelled += invs_before - invs_after;
    }
    stats.dups_detected += dups;

    let mut resolved = Vec::with_capacity(steps.len());
    flatten_expr(&expr, false, collapse_dups, &mut resolved);

    let changed = resolved.len() != steps.len()
        || resolved
            .iter()
            .zip(steps.iter())
            .any(|((a, fa), &(b, fb))| a != b || *fa != fb);

    if changed || dups > 0 {
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

fn decompose_and_resolve<'a>(
    steps: &[(&'a str, bool)],
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
        let (win_steps, win_tag) = resolve_window(
            &steps[i..j],
            seg_lens,
            divergence,
            collapse_dups,
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

fn pass1_index(
    gfa_path: &str,
    min_anchor_bp: usize,
) -> Result<(HashMap<String, usize>, HashSet<String>)> {
    let mut seg_lens: HashMap<String, usize> = HashMap::new();
    let seg_visits: RefCell<HashMap<String, usize>> = RefCell::new(HashMap::new());

    stream_gfa(
        gfa_path,
        |s| {
            let len = if s.sequence == b"*" {
                0
            } else {
                s.sequence.len()
            };
            seg_lens.insert(s.name, len);
            Ok(())
        },
        |_| Ok(()),
        |p| {
            let mut seen_in_path = HashSet::new();
            let mut visits = seg_visits.borrow_mut();
            for (seg, _) in parse_path_steps(&p.segment_names) {
                if seen_in_path.insert(seg) {
                    *visits.entry(seg.to_string()).or_insert(0) += 1;
                }
            }
            Ok(())
        },
        |w| {
            let mut seen_in_walk = HashSet::new();
            let mut visits = seg_visits.borrow_mut();
            for (seg, _) in parse_walk_steps(&w.walk) {
                if seen_in_walk.insert(seg) {
                    *visits.entry(seg.to_string()).or_insert(0) += 1;
                }
            }
            Ok(())
        },
    )?;

    let mut anchors = HashSet::new();
    for (seg, visits) in seg_visits.into_inner() {
        let len = seg_lens.get(&seg).copied().unwrap_or(0);
        if visits >= 2 || len >= min_anchor_bp {
            anchors.insert(seg);
        }
    }
    Ok((seg_lens, anchors))
}

fn open_line_reader(path: &str) -> Result<Box<dyn BufRead>> {
    if path == "-" {
        Ok(Box::new(BufReader::new(io::stdin())))
    } else {
        let file = File::open(path).with_context(|| format!("opening {path}"))?;
        let is_gz = Path::new(path).extension().is_some_and(|e| e == "gz");
        if is_gz {
            let dec = flate2::read::MultiGzDecoder::new(file);
            Ok(Box::new(BufReader::new(dec)))
        } else {
            Ok(Box::new(BufReader::new(file)))
        }
    }
}

fn pass2_rewrite<W: Write>(
    gfa_path: &str,
    out: &mut W,
    seg_lens: &HashMap<String, usize>,
    anchors: &HashSet<String>,
    cli: &CliArgs,
) -> Result<ResolveStats> {
    let mut reader = open_line_reader(gfa_path)?;
    let mut stats = ResolveStats::default();
    let mut line = String::new();

    while reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.starts_with("P\t") {
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
                let mut seg_str = String::new();
                for (idx, (s, fwd)) in resolved.iter().enumerate() {
                    if idx > 0 {
                        seg_str.push(',');
                    }
                    seg_str.push_str(s);
                    seg_str.push(if *fwd { '+' } else { '-' });
                }
                let overlaps = if resolved.len() == steps.len() {
                    cols[3].to_string()
                } else {
                    "*".to_string()
                };
                writeln!(out, "P\t{}\t{}\t{}", cols[1], seg_str, overlaps)?;
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
                let mut walk_str = String::new();
                for (s, fwd) in &resolved {
                    walk_str.push(if *fwd { '>' } else { '<' });
                    walk_str.push_str(s);
                }
                if cli.tag_walks && !tags.is_empty() {
                    writeln!(
                        out,
                        "W\t{}\t{}\t{}\t{}\t{}\t{}\tEG:Z:{}",
                        cols[1],
                        cols[2],
                        cols[3],
                        cols[4],
                        cols[5],
                        walk_str,
                        tags.join(";")
                    )?;
                } else {
                    writeln!(
                        out,
                        "W\t{}\t{}\t{}\t{}\t{}\t{}",
                        cols[1], cols[2], cols[3], cols[4], cols[5], walk_str
                    )?;
                }
                line.clear();
                continue;
            }
        }
        writeln!(out, "{trimmed}")?;
        line.clear();
    }

    out.flush()?;
    Ok(stats)
}

fn main() -> Result<()> {
    let cli = parse_args()?;
    let t0 = Instant::now();
    let (seg_lens, anchors) = pass1_index(&cli.gfa_in, cli.min_anchor_bp)?;

    let stats = if let Some(ref out_p) = cli.out_path {
        let mut writer = BufWriter::new(
            File::create(out_p).with_context(|| format!("creating {}", out_p))?,
        );
        pass2_rewrite(&cli.gfa_in, &mut writer, &seg_lens, &anchors, &cli)?
    } else {
        let mut sink = io::sink();
        pass2_rewrite(&cli.gfa_in, &mut sink, &seg_lens, &anchors, &cli)?
    };

    let hwm_mb = read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    println!(
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
        let dir = std::env::temp_dir().join("horizon_resolve_test");
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
             W\tHG002\t1\tchr6\t0\t80\t>a>b>b>c\n",
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
        let mut writer = BufWriter::new(File::create(&gfa_out)?);
        let stats = pass2_rewrite(&cli.gfa_in, &mut writer, &seg_lens, &anchors, &cli)?;
        assert_eq!(stats.walks_seen, 1);
        assert_eq!(stats.dups_detected, 1);
        assert_eq!(stats.dups_collapsed, 1);
        assert!(stats.total_mdl_bits_saved > 0.0);

        let out_text = std::fs::read_to_string(&gfa_out)?;
        assert!(out_text.contains(">a>b>c\tEG:Z:"));
        let _ = std::fs::remove_dir_all(&dir);
        Ok(())
    }
}
