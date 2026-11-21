//! Two-pass streaming GFA v1.0/v1.1 cleaner using BOCPD + MDL junction scoring and Horizon Bubbles.

use anyhow::{Context, Result, bail};
use horizon_aligner::bocpd::{BocpdScratch, bocpd_gc_profile_into};
use horizon_aligner::gfa::{
    GfaRecord, SpoolableGfaInput, open_output_writer, parse_path_steps, parse_walk_steps,
    stream_gfa_records,
};
use horizon_aligner::seeding::revcomp;
use std::collections::HashMap;
use std::env;
use std::io::Write;
use std::ops::ControlFlow;
use std::time::Instant;

#[derive(Debug, Clone)]
struct SegFlank {
    len: usize,
    prefix: Vec<u8>,
    suffix: Vec<u8>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct DirectedEdge {
    from: u32,
    from_fwd: bool,
    to: u32,
    to_fwd: bool,
}

impl DirectedEdge {
    fn new(from: u32, from_fwd: bool, to: u32, to_fwd: bool) -> Self {
        Self {
            from,
            from_fwd,
            to,
            to_fwd,
        }
    }

    fn rev(self) -> Self {
        Self {
            from: self.to,
            from_fwd: !self.to_fwd,
            to: self.from,
            to_fwd: !self.from_fwd,
        }
    }
}

fn intern_seg(map: &mut HashMap<String, u32>, next_id: &mut u32, name: &str) -> u32 {
    if let Some(&id) = map.get(name) {
        return id;
    }
    let id = *next_id;
    *next_id = next_id
        .checked_add(1)
        .expect("segment interner exceeded u32::MAX segments");
    map.insert(name.to_string(), id);
    id
}

fn oriented_3p_flank(flank: &SegFlank, is_fwd: bool) -> std::borrow::Cow<'_, [u8]> {
    if is_fwd {
        std::borrow::Cow::Borrowed(&flank.suffix)
    } else {
        std::borrow::Cow::Owned(revcomp(&flank.prefix))
    }
}

fn oriented_5p_flank(flank: &SegFlank, is_fwd: bool) -> std::borrow::Cow<'_, [u8]> {
    if is_fwd {
        std::borrow::Cow::Borrowed(&flank.prefix)
    } else {
        std::borrow::Cow::Owned(revcomp(&flank.suffix))
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
fn score_junction_cp_bits_with_scratch(
    flank_u: &[u8],
    flank_v: &[u8],
    scratch: &mut BocpdScratch,
    concat: &mut Vec<u8>,
    gc_prof: &mut Vec<f64>,
) -> f64 {
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

    concat.clear();
    concat.extend_from_slice(flank_u);
    concat.extend_from_slice(flank_v);
    bocpd_gc_profile_into(concat, 1.0 / 120.0, 80, scratch, gc_prof);
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

fn parse_cigar_overlap_bp(cigar: &str) -> usize {
    if cigar.is_empty() || cigar == "*" {
        return 0;
    }
    let mut ov = 0usize;
    let mut num = 0usize;
    for b in cigar.bytes() {
        if b.is_ascii_digit() {
            num = num.saturating_mul(10).saturating_add(usize::from(b - b'0'));
        } else {
            if matches!(b, b'M' | b'=' | b'X') {
                ov = ov.saturating_add(num);
            }
            num = 0;
        }
    }
    ov
}

struct CleanPass1 {
    seg_ids: HashMap<String, u32>,
    seg_flanks: HashMap<u32, SegFlank>,
    edges: Vec<DirectedEdge>,
    walk_support: HashMap<DirectedEdge, usize>,
    link_overlaps: HashMap<DirectedEdge, usize>,
}

fn record_walk_support(
    steps: &[(&str, bool)],
    seg_ids: &mut HashMap<String, u32>,
    next_seg_id: &mut u32,
    walk_support: &mut HashMap<DirectedEdge, usize>,
) {
    for win in steps.windows(2) {
        let u = intern_seg(seg_ids, next_seg_id, win[0].0);
        let v = intern_seg(seg_ids, next_seg_id, win[1].0);
        let e = DirectedEdge::new(u, win[0].1, v, win[1].1);
        let rev = e.rev();
        *walk_support.entry(e).or_insert(0) += 1;
        if rev != e {
            *walk_support.entry(rev).or_insert(0) += 1;
        }
    }
}

fn collect_clean_pass1(gfa_path: &str, flank_bp: usize) -> Result<CleanPass1> {
    let mut seg_ids = HashMap::<String, u32>::new();
    let mut next_seg_id = 0_u32;
    let mut seg_flanks = HashMap::<u32, SegFlank>::new();
    let mut edges = Vec::<DirectedEdge>::new();
    let mut walk_support = HashMap::<DirectedEdge, usize>::new();
    let mut link_overlaps = HashMap::<DirectedEdge, usize>::new();

    stream_gfa_records(gfa_path, |rec| {
        match rec {
            GfaRecord::Segment(s) => {
                let id = intern_seg(&mut seg_ids, &mut next_seg_id, &s.name);
                let has_seq = s.sequence != b"*";
                let len = s.seq_len();
                let prefix = if has_seq && len > 0 {
                    s.sequence[..len.min(flank_bp)].to_vec()
                } else {
                    Vec::new()
                };
                let suffix = if has_seq && len > 0 {
                    s.sequence[len.saturating_sub(flank_bp)..].to_vec()
                } else {
                    Vec::new()
                };
                seg_flanks.insert(
                    id,
                    SegFlank {
                        len,
                        prefix,
                        suffix,
                    },
                );
            }
            GfaRecord::Link(l) => {
                let u = intern_seg(&mut seg_ids, &mut next_seg_id, &l.from);
                let v = intern_seg(&mut seg_ids, &mut next_seg_id, &l.to);
                let e = DirectedEdge::new(u, l.from_orient != b'-', v, l.to_orient != b'-');
                let ov = parse_cigar_overlap_bp(&l.overlap);
                if ov > 0 {
                    link_overlaps.entry(e).or_insert(ov);
                    link_overlaps.entry(e.rev()).or_insert(ov);
                }
                edges.push(e);
            }
            GfaRecord::Path(p) => {
                let steps = parse_path_steps(&p.segment_names);
                record_walk_support(&steps, &mut seg_ids, &mut next_seg_id, &mut walk_support);
            }
            GfaRecord::Walk(w) => {
                let steps = parse_walk_steps(&w.walk);
                record_walk_support(&steps, &mut seg_ids, &mut next_seg_id, &mut walk_support);
            }
        }
        Ok(ControlFlow::Continue(()))
    })?;

    Ok(CleanPass1 {
        seg_ids,
        seg_flanks,
        edges,
        walk_support,
        link_overlaps,
    })
}

fn score_links_symmetrically(
    pass1: &CleanPass1,
    cfg: &CleanConfig,
) -> (HashMap<DirectedEdge, f64>, usize, usize) {
    let mut severed_info: HashMap<DirectedEdge, f64> = HashMap::new();
    let mut canon_decisions: HashMap<DirectedEdge, (bool, f64)> = HashMap::new();
    let mut kept_links = 0usize;
    let mut severed_links = 0usize;
    let mut scratch = BocpdScratch::default();
    let mut concat = Vec::new();
    let mut gc_prof = Vec::new();

    for &edge in &pass1.edges {
        let rev_edge = edge.rev();
        let canon = if (edge.from, edge.from_fwd, edge.to, edge.to_fwd)
            <= (
                rev_edge.from,
                rev_edge.from_fwd,
                rev_edge.to,
                rev_edge.to_fwd,
            ) {
            edge
        } else {
            rev_edge
        };

        let (should_sever, cp_bits) = *canon_decisions.entry(canon).or_insert_with(|| {
            let sup = pass1
                .walk_support
                .get(&canon)
                .copied()
                .unwrap_or(0)
                .max(pass1.walk_support.get(&canon.rev()).copied().unwrap_or(0));
            let cp = match (
                pass1.seg_flanks.get(&canon.from),
                pass1.seg_flanks.get(&canon.to),
            ) {
                (Some(fu), Some(fv)) => {
                    let u3 = oriented_3p_flank(fu, canon.from_fwd);
                    let v5 = oriented_5p_flank(fv, canon.to_fwd);
                    score_junction_cp_bits_with_scratch(
                        &u3,
                        &v5,
                        &mut scratch,
                        &mut concat,
                        &mut gc_prof,
                    )
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
        }
    }
    (severed_info, kept_links, severed_links)
}

fn emit_clean_path<W: Write>(
    p: &horizon_aligner::gfa::GfaPath,
    seg_ids: &HashMap<String, u32>,
    severed_info: &HashMap<DirectedEdge, f64>,
    out: &mut W,
) -> Result<bool> {
    let steps = parse_path_steps(&p.segment_names);
    if steps.len() <= 1 {
        writeln!(out, "P\t{}\t{}\t{}", p.name, p.segment_names, p.overlaps)?;
        return Ok(false);
    }
    let mut chunks: Vec<Vec<(&str, bool)>> = vec![vec![steps[0]]];
    for win in steps.windows(2) {
        let is_severed = seg_ids
            .get(win[0].0)
            .zip(seg_ids.get(win[1].0))
            .is_some_and(|(&u, &v)| {
                let e = DirectedEdge::new(u, win[0].1, v, win[1].1);
                severed_info.contains_key(&e)
            });
        if is_severed {
            chunks.push(vec![win[1]]);
        } else if let Some(last) = chunks.last_mut() {
            last.push(win[1]);
        }
    }
    if chunks.len() == 1 {
        writeln!(out, "P\t{}\t{}\t{}", p.name, p.segment_names, p.overlaps)?;
        return Ok(false);
    }
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
            out,
            "P\t{}_part{}\t{}\t{}",
            p.name,
            ci + 1,
            segs.join(","),
            ovs
        )?;
    }
    Ok(true)
}

fn emit_clean_walk<W: Write>(
    w: &horizon_aligner::gfa::GfaWalk,
    pass1: &CleanPass1,
    severed_info: &HashMap<DirectedEdge, f64>,
    out: &mut W,
) -> Result<bool> {
    type WalkChunk<'a> = (usize, usize, Vec<(&'a str, bool)>);
    let seg_ids = &pass1.seg_ids;
    let seg_flanks = &pass1.seg_flanks;
    let emit_walk = |out: &mut dyn Write, cs: usize, ce: usize, walk_str: &str| -> Result<()> {
        if let Some(ref tags) = w.tags {
            writeln!(
                out,
                "W\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                w.sample, w.haplotype, w.seq_id, cs, ce, walk_str, tags
            )?;
        } else {
            writeln!(
                out,
                "W\t{}\t{}\t{}\t{}\t{}\t{}",
                w.sample, w.haplotype, w.seq_id, cs, ce, walk_str
            )?;
        }
        Ok(())
    };
    let steps = parse_walk_steps(&w.walk);
    if steps.len() <= 1 {
        emit_walk(out, w.start, w.end, &w.walk)?;
        return Ok(false);
    }
    let mut chunks: Vec<WalkChunk<'_>> = Vec::new();
    let mut cur_start = w.start;
    let first_len = seg_ids
        .get(steps[0].0)
        .and_then(|id| seg_flanks.get(id))
        .map_or(0, |f| f.len);
    let mut cur_end = cur_start.saturating_add(first_len);
    let mut cur_steps = vec![steps[0]];

    for win in steps.windows(2) {
        let step_len = seg_ids
            .get(win[1].0)
            .and_then(|id| seg_flanks.get(id))
            .map_or(0, |f| f.len);
        let edge_opt = seg_ids
            .get(win[0].0)
            .zip(seg_ids.get(win[1].0))
            .map(|(&u, &v)| DirectedEdge::new(u, win[0].1, v, win[1].1));
        let ov_bp = edge_opt
            .and_then(|e| pass1.link_overlaps.get(&e).copied())
            .unwrap_or(0);
        let is_severed = edge_opt.is_some_and(|e| severed_info.contains_key(&e));
        if is_severed {
            chunks.push((cur_start, cur_end, std::mem::take(&mut cur_steps)));
            cur_start = cur_end.saturating_sub(ov_bp);
            cur_end = cur_start.saturating_add(step_len);
            cur_steps.push(win[1]);
        } else {
            cur_steps.push(win[1]);
            cur_end = cur_end.saturating_add(step_len.saturating_sub(ov_bp));
        }
    }
    chunks.push((cur_start, cur_end, cur_steps));

    if chunks.len() == 1 {
        emit_walk(out, w.start, w.end, &w.walk)?;
        return Ok(false);
    }
    for (cs, ce, chunk) in chunks {
        let mut walk_str = String::new();
        for (s, fwd) in chunk {
            walk_str.push(if fwd { '>' } else { '<' });
            walk_str.push_str(s);
        }
        emit_walk(out, cs, ce, &walk_str)?;
    }
    Ok(true)
}

fn clean_gfa<W: Write>(gfa_path: &str, cfg: &CleanConfig, out: &mut W) -> Result<CleanStats> {
    let pass1 = collect_clean_pass1(gfa_path, cfg.flank_bp)?;
    let (severed_info, kept_links, severed_links) = score_links_symmetrically(&pass1, cfg);

    writeln!(out, "H\tVN:Z:1.1")?;
    let mut emitted_bubbles = 0usize;
    let mut split_walks = 0usize;

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
                let edge_cp = pass1
                    .seg_ids
                    .get(&l.from)
                    .zip(pass1.seg_ids.get(&l.to))
                    .and_then(|(&u, &v)| {
                        let edge =
                            DirectedEdge::new(u, l.from_orient != b'-', v, l.to_orient != b'-');
                        severed_info.get(&edge).copied()
                    });
                if let Some(cp_bits) = edge_cp {
                    emitted_bubbles += 1;
                    writeln!(
                        out,
                        "S\t{}__{}__hb_{}\t*\tLN:i:0\tHB:i:1\tCP:f:{:.2}",
                        l.from, l.to, emitted_bubbles, cp_bits
                    )?;
                } else {
                    writeln!(
                        out,
                        "L\t{}\t{}\t{}\t{}\t{}",
                        l.from, l.from_orient as char, l.to, l.to_orient as char, l.overlap
                    )?;
                }
            }
            GfaRecord::Path(p) => {
                if emit_clean_path(&p, &pass1.seg_ids, &severed_info, out)? {
                    split_walks += 1;
                }
            }
            GfaRecord::Walk(w) => {
                if emit_clean_walk(&w, &pass1, &severed_info, out)? {
                    split_walks += 1;
                }
            }
        }
        Ok(ControlFlow::Continue(()))
    })?;

    Ok(CleanStats {
        total_links: pass1.edges.len(),
        kept_links,
        severed_links,
        emitted_bubbles,
        split_walks,
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
            "-h" | "--help" => {
                println!(
                    "usage: horizon_clean --gfa <in.gfa[.gz]|-> [--out <clean.gfa>] [--flank 128] [--min-cp-bits 12.0] [--prune-unsupported]"
                );
                return Ok(());
            }
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
                anyhow::ensure!(
                    flank_bp >= 12,
                    "--flank must be at least 12 bp (got {flank_bp})"
                );
                i += 2;
            }
            "--min-cp-bits" => {
                let val = args.get(i + 1).context("missing value for --min-cp-bits")?;
                min_cp_bits = val.parse().context("parsing --min-cp-bits")?;
                anyhow::ensure!(
                    min_cp_bits.is_finite() && min_cp_bits >= 0.0,
                    "--min-cp-bits must be a finite non-negative number (got {min_cp_bits})"
                );
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

    let mut w = open_output_writer(out_path.as_deref())?;
    let stats = clean_gfa(effective_gfa, &cfg, &mut w)?;
    w.finish()?;

    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
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
        let gfa_path =
            std::env::temp_dir().join(format!("horizon_clean_test_{}.gfa", std::process::id()));
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
        let gfa_path =
            std::env::temp_dir().join(format!("horizon_clean_sym_test_{}.gfa", std::process::id()));
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

    #[test]
    fn splits_path_and_walk_across_severed_link() {
        let gfa_path = std::env::temp_dir().join(format!(
            "horizon_clean_split_test_{}.gfa",
            std::process::id()
        ));
        let s1 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let s2 = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\ts2\t{s1}\n\
             S\ts3\t{s2}\n\
             S\ts4\t{s2}\n\
             L\ts1\t+\ts2\t+\t0M\n\
             L\ts2\t+\ts3\t+\t0M\n\
             L\ts2\t+\ts4\t+\t0M\n\
             P\tp1\ts1+,s2+,s3+\t0M,0M\n\
             W\tCHM13\t0\tchr1\t0\t180\t>s1>s2>s4\n"
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

        assert_eq!(stats.split_walks, 2);
        let out_str = String::from_utf8(out).unwrap();
        assert!(out_str.contains("P\tp1_part1\ts1+,s2+\t0M\n"));
        assert!(out_str.contains("P\tp1_part2\ts3+\t*\n"));
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t0\t120\t>s1>s2\n"));
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t120\t180\t>s4\n"));
    }

    #[test]
    fn preserves_walk_without_explicit_l_lines() {
        let gfa_path = std::env::temp_dir().join(format!(
            "horizon_clean_no_l_test_{}.gfa",
            std::process::id()
        ));
        let s1 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\ts2\t{s1}\n\
             W\tCHM13\t0\tchr1\t0\t120\t>s1>s2\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 60,
            min_cp_bits: 12.0,
            prune_unsupported: true,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(stats.split_walks, 0);
        let out_str = String::from_utf8(out).unwrap();
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t0\t120\t>s1>s2\n"));
    }

    #[test]
    fn prune_unsupported_severs_low_cp_unsupported_link() {
        let gfa_path = std::env::temp_dir().join(format!(
            "horizon_clean_prune_test_{}.gfa",
            std::process::id()
        ));
        let s1 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\ts2\t{s1}\n\
             L\ts1\t+\ts2\t+\t0M\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 60,
            min_cp_bits: 100.0,
            prune_unsupported: true,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(stats.total_links, 1);
        assert_eq!(stats.kept_links, 0);
        assert_eq!(stats.severed_links, 1);
        let out_str = String::from_utf8(out).unwrap();
        assert!(!out_str.contains("\nL\t"));
        assert!(out_str.contains("HB:i:1"));
    }

    #[test]
    fn walk_split_accounts_for_link_overlap_cigar() {
        let gfa_path =
            std::env::temp_dir().join(format!("horizon_clean_ov_test_{}.gfa", std::process::id()));
        let s1 = "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let s2 = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
        let s3 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let gfa = format!(
            "H\tVN:Z:1.1\n\
             S\ts1\t{s1}\n\
             S\ts2\t{s2}\n\
             S\ts3\t{s3}\n\
             L\ts1\t+\ts2\t+\t10M\n\
             L\ts2\t+\ts3\t+\t10M\n\
             W\tCHM13\t0\tchr1\t0\t160\t>s1>s2>s3\n"
        );
        std::fs::write(&gfa_path, gfa).unwrap();

        let cfg = CleanConfig {
            flank_bp: 60,
            min_cp_bits: 5.0,
            prune_unsupported: false,
        };
        let mut out = Vec::new();
        let stats = clean_gfa(gfa_path.to_str().unwrap(), &cfg, &mut out).unwrap();
        let _ = std::fs::remove_file(&gfa_path);

        assert_eq!(stats.split_walks, 1);
        let out_str = String::from_utf8(out).unwrap();
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t0\t60\t>s1\n"));
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t50\t110\t>s2\n"));
        assert!(out_str.contains("W\tCHM13\t0\tchr1\t100\t160\t>s3\n"));
    }
}
