//! Two-pass streaming GFA v1.0/v1.1 cleaner using BOCPD + MDL junction scoring and Horizon Bubbles.

use anyhow::{Context, Result, bail};
use horizon_aligner::bocpd::bocpd_gc_profile;
use horizon_aligner::gfa::{
    GfaLink, SpoolableGfaInput, parse_path_steps, parse_walk_steps, stream_gfa,
};
use horizon_aligner::seeding::revcomp;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::time::Instant;

#[derive(Debug, Clone)]
struct SegFlank {
    len: usize,
    prefix: Vec<u8>,
    suffix: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct DirectedEdge {
    from: String,
    from_fwd: bool,
    to: String,
    to_fwd: bool,
}

impl DirectedEdge {
    fn new(from: &str, from_fwd: bool, to: &str, to_fwd: bool) -> Self {
        Self {
            from: from.to_string(),
            from_fwd,
            to: to.to_string(),
            to_fwd,
        }
    }

    fn rev(&self) -> Self {
        Self {
            from: self.to.clone(),
            from_fwd: !self.to_fwd,
            to: self.from.clone(),
            to_fwd: !self.from_fwd,
        }
    }
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

fn oriented_3p_flank(flank: &SegFlank, is_fwd: bool) -> Vec<u8> {
    if is_fwd {
        flank.suffix.clone()
    } else {
        revcomp(&flank.prefix)
    }
}

fn oriented_5p_flank(flank: &SegFlank, is_fwd: bool) -> Vec<u8> {
    if is_fwd {
        flank.prefix.clone()
    } else {
        revcomp(&flank.suffix)
    }
}

fn gc_counts(seq: &[u8]) -> (f64, f64) {
    let mut gc = 0.0_f64;
    let mut at = 0.0_f64;
    for &b in seq {
        match b.to_ascii_uppercase() {
            b'G' | b'C' => gc += 1.0,
            b'A' | b'T' => at += 1.0,
            _ => {}
        }
    }
    (gc, gc + at)
}

fn bernoulli_ll_bits(k: f64, n: f64) -> f64 {
    if n <= 0.0 || k <= 0.0 || k >= n {
        return 0.0;
    }
    let p = (k / n).clamp(1e-6, 1.0 - 1e-6);
    k * p.log2() + (n - k) * (1.0 - p).log2()
}

/// Score a GFA link junction using MDL two-regime changepoint gain (bits) gated by BOCPD GC shift.
/// Positive bits indicate a statistically significant compositional break at the junction.
fn score_junction_cp_bits(flank_u: &[u8], flank_v: &[u8]) -> f64 {
    if flank_u.len() < 12 || flank_v.len() < 12 {
        return 0.0;
    }
    let (ku, nu) = gc_counts(flank_u);
    let (kv, nv) = gc_counts(flank_v);
    let n_tot = nu + nv;
    if nu < 10.0 || nv < 10.0 || n_tot < 24.0 {
        return 0.0;
    }
    let ll_null = bernoulli_ll_bits(ku + kv, n_tot);
    let ll_alt = bernoulli_ll_bits(ku, nu) + bernoulli_ll_bits(kv, nv);
    let mdl_penalty = 1.5 * n_tot.log2();
    let cp_bits = (ll_alt - ll_null) - mdl_penalty;
    if cp_bits <= 0.0 {
        return cp_bits;
    }

    let mut concat = Vec::with_capacity(flank_u.len() + flank_v.len());
    concat.extend_from_slice(flank_u);
    concat.extend_from_slice(flank_v);
    let gc_prof = bocpd_gc_profile(&concat, 1.0 / 120.0, 80);
    let split = flank_u.len();
    let left_gc = gc_prof[split - 1];
    let right_gc = gc_prof[(split + 15).min(gc_prof.len() - 1)];
    let bocpd_jump = (right_gc - left_gc).abs();
    if bocpd_jump < 0.10 { 0.0 } else { cp_bits }
}

struct CleanConfig {
    flank_bp: usize,
    min_cp_bits: f64,
    prune_unsupported: bool,
}

struct CleanStats {
    total_links: usize,
    kept_links: usize,
    severed_links: usize,
    emitted_bubbles: usize,
    split_walks: usize,
}

#[allow(clippy::too_many_lines)]
#[allow(clippy::type_complexity)]
fn clean_gfa<W: Write>(gfa_path: &str, cfg: &CleanConfig, out: &mut W) -> Result<CleanStats> {
    let seg_flanks = std::cell::RefCell::new(HashMap::<String, SegFlank>::new());
    let links = std::cell::RefCell::new(Vec::<GfaLink>::new());
    let walk_support = std::cell::RefCell::new(HashMap::<DirectedEdge, usize>::new());

    // Pass 1: index segment flanks (bounded by cfg.flank_bp), links, and walk edge support.
    stream_gfa(
        gfa_path,
        |s| {
            let (len, has_seq) = if s.sequence == b"*" {
                let ln = s
                    .tags
                    .iter()
                    .find_map(|t| {
                        t.strip_prefix("LN:i:")
                            .and_then(|v| v.parse::<usize>().ok())
                    })
                    .unwrap_or(0);
                (ln, false)
            } else {
                (s.sequence.len(), true)
            };
            let prefix = if has_seq && len > 0 {
                s.sequence[..len.min(cfg.flank_bp)].to_vec()
            } else {
                Vec::new()
            };
            let suffix = if has_seq && len > 0 {
                s.sequence[len.saturating_sub(cfg.flank_bp)..].to_vec()
            } else {
                Vec::new()
            };
            seg_flanks.borrow_mut().insert(
                s.name,
                SegFlank {
                    len,
                    prefix,
                    suffix,
                },
            );
            Ok(())
        },
        |l| {
            links.borrow_mut().push(l);
            Ok(())
        },
        |p| {
            let steps = parse_path_steps(&p.segment_names);
            let mut sup = walk_support.borrow_mut();
            for win in steps.windows(2) {
                let e = DirectedEdge::new(win[0].0, win[0].1, win[1].0, win[1].1);
                let rev = e.rev();
                *sup.entry(e.clone()).or_insert(0) += 1;
                if rev != e {
                    *sup.entry(rev).or_insert(0) += 1;
                }
            }
            Ok(())
        },
        |w| {
            let steps = parse_walk_steps(&w.walk);
            let mut sup = walk_support.borrow_mut();
            for win in steps.windows(2) {
                let e = DirectedEdge::new(win[0].0, win[0].1, win[1].0, win[1].1);
                let rev = e.rev();
                *sup.entry(e.clone()).or_insert(0) += 1;
                if rev != e {
                    *sup.entry(rev).or_insert(0) += 1;
                }
            }
            Ok(())
        },
    )?;

    let seg_flanks = seg_flanks.into_inner();
    let links = links.into_inner();
    let walk_support = walk_support.into_inner();

    // Score links symmetrically per canonical bidirected edge pair so L u+ v+ and L v- u-
    // never receive asymmetric sever/keep decisions.
    let mut kept_edges: HashSet<DirectedEdge> = HashSet::new();
    let mut severed_info: HashMap<DirectedEdge, f64> = HashMap::new();
    let mut canon_decisions: HashMap<DirectedEdge, (bool, f64)> = HashMap::new();
    let mut kept_links = 0usize;
    let mut severed_links = 0usize;

    for l in &links {
        let edge = DirectedEdge::new(&l.from, l.from_orient != b'-', &l.to, l.to_orient != b'-');
        let rev_edge = edge.rev();
        let canon = if (&edge.from, edge.from_fwd, &edge.to, edge.to_fwd)
            <= (
                &rev_edge.from,
                rev_edge.from_fwd,
                &rev_edge.to,
                rev_edge.to_fwd,
            ) {
            edge.clone()
        } else {
            rev_edge.clone()
        };

        let (should_sever, cp_bits) = *canon_decisions.entry(canon.clone()).or_insert_with(|| {
            let sup = walk_support
                .get(&canon)
                .copied()
                .unwrap_or(0)
                .max(walk_support.get(&canon.rev()).copied().unwrap_or(0));
            let cp = match (seg_flanks.get(&canon.from), seg_flanks.get(&canon.to)) {
                (Some(fu), Some(fv)) => {
                    let u3 = oriented_3p_flank(fu, canon.from_fwd);
                    let v5 = oriented_5p_flank(fv, canon.to_fwd);
                    score_junction_cp_bits(&u3, &v5)
                }
                _ => 0.0,
            };
            let sever = (cfg.prune_unsupported && sup == 0) || (sup <= 1 && cp >= cfg.min_cp_bits);
            (sever, cp)
        });

        if should_sever {
            severed_links += 1;
            severed_info.insert(edge, cp_bits);
            severed_info.insert(rev_edge, cp_bits);
        } else {
            kept_links += 1;
            kept_edges.insert(edge);
            kept_edges.insert(rev_edge);
        }
    }

    // Pass 2: stream GFA and emit cleaned S, L, Horizon Bubble, and repaired P/W records.
    writeln!(out, "H\tVN:Z:1.1")?;
    let mut emitted_bubbles = 0usize;
    let split_walks = std::cell::Cell::new(0usize);
    let out_cell = std::cell::RefCell::new(out);

    stream_gfa(
        gfa_path,
        |s| {
            let seq_str = String::from_utf8_lossy(&s.sequence);
            if let Some(ref tags) = s.tags {
                writeln!(
                    out_cell.borrow_mut(),
                    "S\t{}\t{}\t{}",
                    s.name,
                    seq_str,
                    tags
                )?;
            } else {
                writeln!(out_cell.borrow_mut(), "S\t{}\t{}", s.name, seq_str)?;
            }
            Ok(())
        },
        |l| {
            let edge =
                DirectedEdge::new(&l.from, l.from_orient != b'-', &l.to, l.to_orient != b'-');
            if let Some(&cp_bits) = severed_info.get(&edge) {
                emitted_bubbles += 1;
                writeln!(
                    out_cell.borrow_mut(),
                    "S\t{}__{}__hb_{}\t*\tLN:i:0\tHB:i:1\tCP:f:{:.2}",
                    l.from,
                    l.to,
                    emitted_bubbles,
                    cp_bits
                )?;
            } else {
                writeln!(
                    out_cell.borrow_mut(),
                    "L\t{}\t{}\t{}\t{}\t{}",
                    l.from,
                    l.from_orient as char,
                    l.to,
                    l.to_orient as char,
                    l.overlap
                )?;
            }
            Ok(())
        },
        |p| {
            let steps = parse_path_steps(&p.segment_names);
            if steps.len() <= 1 {
                writeln!(
                    out_cell.borrow_mut(),
                    "P\t{}\t{}\t{}",
                    p.name,
                    p.segment_names,
                    p.overlaps
                )?;
                return Ok(());
            }
            let mut chunks: Vec<Vec<(&str, bool)>> = vec![vec![steps[0]]];
            for win in steps.windows(2) {
                let e = DirectedEdge::new(win[0].0, win[0].1, win[1].0, win[1].1);
                if kept_edges.contains(&e) {
                    if let Some(last) = chunks.last_mut() {
                        last.push(win[1]);
                    }
                } else {
                    chunks.push(vec![win[1]]);
                }
            }
            if chunks.len() == 1 {
                writeln!(
                    out_cell.borrow_mut(),
                    "P\t{}\t{}\t{}",
                    p.name,
                    p.segment_names,
                    p.overlaps
                )?;
            } else {
                split_walks.set(split_walks.get() + 1);
                let orig_ovs: Vec<&str> = if p.overlaps == "*" || p.overlaps.is_empty() {
                    Vec::new()
                } else {
                    p.overlaps.split(',').collect()
                };
                let mut cur_start = 0usize;
                for (ci, chunk) in chunks.iter().enumerate() {
                    let segs: Vec<String> = chunk
                        .iter()
                        .map(|(s, fwd)| format!("{s}{}", if *fwd { '+' } else { '-' }))
                        .collect();
                    let ovs = if chunk.len() > 1 {
                        if orig_ovs.len() >= cur_start + chunk.len() - 1 {
                            orig_ovs[cur_start..cur_start + chunk.len() - 1].join(",")
                        } else {
                            vec!["0M"; chunk.len() - 1].join(",")
                        }
                    } else {
                        "*".to_string()
                    };
                    cur_start += chunk.len();
                    writeln!(
                        out_cell.borrow_mut(),
                        "P\t{}_part{}\t{}\t{}",
                        p.name,
                        ci + 1,
                        segs.join(","),
                        ovs
                    )?;
                }
            }
            Ok(())
        },
        |w| {
            let steps = parse_walk_steps(&w.walk);
            if steps.len() <= 1 {
                writeln!(
                    out_cell.borrow_mut(),
                    "W\t{}\t{}\t{}\t{}\t{}\t{}",
                    w.sample,
                    w.haplotype,
                    w.seq_id,
                    w.start,
                    w.end,
                    w.walk
                )?;
                return Ok(());
            }
            let mut chunks: Vec<(usize, usize, Vec<(&str, bool)>)> = Vec::new();
            let mut cur_start = w.start;
            let first_len = seg_flanks.get(steps[0].0).map_or(0, |f| f.len);
            let mut cur_end = cur_start + first_len;
            let mut cur_steps = vec![steps[0]];

            for win in steps.windows(2) {
                let step_len = seg_flanks.get(win[1].0).map_or(0, |f| f.len);
                let e = DirectedEdge::new(win[0].0, win[0].1, win[1].0, win[1].1);
                if kept_edges.contains(&e) {
                    cur_steps.push(win[1]);
                    cur_end += step_len;
                } else {
                    chunks.push((cur_start, cur_end, std::mem::take(&mut cur_steps)));
                    cur_start = cur_end;
                    cur_end = cur_start + step_len;
                    cur_steps.push(win[1]);
                }
            }
            chunks.push((cur_start, cur_end, cur_steps));

            if chunks.len() == 1 {
                writeln!(
                    out_cell.borrow_mut(),
                    "W\t{}\t{}\t{}\t{}\t{}\t{}",
                    w.sample,
                    w.haplotype,
                    w.seq_id,
                    w.start,
                    w.end,
                    w.walk
                )?;
            } else {
                split_walks.set(split_walks.get() + 1);
                for (cs, ce, chunk) in chunks {
                    let mut walk_str = String::new();
                    for (s, fwd) in chunk {
                        walk_str.push(if fwd { '>' } else { '<' });
                        walk_str.push_str(s);
                    }
                    writeln!(
                        out_cell.borrow_mut(),
                        "W\t{}\t{}\t{}\t{}\t{}\t{}",
                        w.sample,
                        w.haplotype,
                        w.seq_id,
                        cs,
                        ce,
                        walk_str
                    )?;
                }
            }
            Ok(())
        },
    )?;

    Ok(CleanStats {
        total_links: links.len(),
        kept_links,
        severed_links,
        emitted_bubbles,
        split_walks: split_walks.get(),
    })
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut gfa_path: Option<String> = None;
    let mut out_path: Option<String> = None;
    let mut flank_bp = 128usize;
    let mut min_cp_bits = 12.0_f64;
    let mut prune_unsupported = false;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--gfa" | "-i" => {
                let val = args.get(i + 1).context("missing value for --gfa")?;
                gfa_path = Some(val.clone());
                i += 2;
            }
            "--out" | "-o" => {
                let val = args.get(i + 1).context("missing value for --out")?;
                out_path = Some(val.clone());
                i += 2;
            }
            "--flank" => {
                let val = args.get(i + 1).context("missing value for --flank")?;
                flank_bp = val.parse().context("parsing --flank")?;
                i += 2;
            }
            "--min-cp-bits" => {
                let val = args.get(i + 1).context("missing value for --min-cp-bits")?;
                min_cp_bits = val.parse().context("parsing --min-cp-bits")?;
                i += 2;
            }
            "--prune-unsupported" => {
                prune_unsupported = true;
                i += 1;
            }
            _ => bail!("unknown argument '{}'", args[i]),
        }
    }

    let Some(gfa) = gfa_path else {
        bail!(
            "usage: horizon_clean --gfa <in.gfa[.gz]|-> [--out <clean.gfa>] [--flank 128] [--min-cp-bits 12.0] [--prune-unsupported]"
        );
    };

    let t0 = Instant::now();
    let cfg = CleanConfig {
        flank_bp,
        min_cp_bits,
        prune_unsupported,
    };

    let spooled = SpoolableGfaInput::new(&gfa)?;
    let effective_gfa = spooled.path();

    let stats = if let Some(ref op) = out_path {
        let mut w = BufWriter::new(File::create(op)?);
        let s = clean_gfa(effective_gfa, &cfg, &mut w)?;
        w.flush()?;
        s
    } else {
        let stdout = io::stdout();
        let mut w = BufWriter::new(stdout.lock());
        let s = clean_gfa(effective_gfa, &cfg, &mut w)?;
        w.flush()?;
        s
    };

    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    eprintln!(
        "horizon_clean: total_links={} kept={} severed={} bubbles={} split_walks={} elapsed={:.2}s peak_rss_mb={:.1}",
        stats.total_links,
        stats.kept_links,
        stats.severed_links,
        stats.emitted_bubbles,
        stats.split_walks,
        t0.elapsed().as_secs_f64(),
        hwm_mb
    );
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn severs_discordant_unsupported_link_and_preserves_valid_walk() {
        let gfa_path = std::env::temp_dir().join("horizon_clean_test.gfa");
        // s1 (60bp AT-rich) -> s2 (60bp AT-rich): continuous, walk-supported
        // s1 (60bp AT-rich) -> s3 (60bp GC-rich): unsupported + strong compositional break
        let s1 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let s2 = "TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA";
        let s3 = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\ts2\t{s2}\n\
             S\ts3\t{s3}\n\
             L\ts1\t+\ts2\t+\t0M\n\
             L\ts1\t+\ts3\t+\t0M\n\
             W\tCHM13\t0\tchr22\t0\t120\t>s1>s2\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 60,
            min_cp_bits: 12.0,
            prune_unsupported: false,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(stats.total_links, 2);
        assert_eq!(stats.kept_links, 1);
        assert_eq!(stats.severed_links, 1);
        assert_eq!(stats.emitted_bubbles, 1);

        let out_str = String::from_utf8(out).unwrap();
        assert!(out_str.contains("L\ts1\t+\ts2\t+\t0M"));
        assert!(!out_str.contains("L\ts1\t+\ts3\t+\t0M"));
        assert!(out_str.contains("HB:i:1"));
        assert!(out_str.contains("W\tCHM13\t0\tchr22\t0\t120\t>s1>s2"));
    }

    #[test]
    fn severs_bidirected_link_pair_symmetrically() {
        let gfa_path = std::env::temp_dir().join("horizon_clean_sym_test.gfa");
        let s1 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let s2 = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\ts2\t{s2}\n\
             L\ts1\t+\ts2\t+\t0M\n\
             L\ts2\t-\ts1\t-\t0M\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 60,
            min_cp_bits: 12.0,
            prune_unsupported: false,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(stats.total_links, 2);
        assert_eq!(stats.kept_links, 0);
        assert_eq!(stats.severed_links, 2);
        let out_str = String::from_utf8(out).unwrap();
        assert!(!out_str.contains("\nL\t"));
    }
}
