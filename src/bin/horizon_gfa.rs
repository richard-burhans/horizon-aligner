//! Streaming GFA v1.0/v1.1 utility (`stats`, `list-walks`, `extract`) using `horizon_aligner::gfa`.

use anyhow::{bail, Context, Result};
use horizon_aligner::gfa::stream_gfa;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

fn rev_comp_into(seq: &[u8], out: &mut Vec<u8>) {
    out.reserve(seq.len());
    for &b in seq.iter().rev() {
        out.push(match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => other,
        });
    }
}

/// Parse a GFA v1.0 path `segment_names` field such as `14+,15-` into `(segment_name, is_forward)` steps.
fn parse_path_steps(seg_names: &str) -> Vec<(&str, bool)> {
    seg_names
        .split(',')
        .filter_map(|tok| {
            let tok = tok.trim();
            if tok.len() < 2 {
                return None;
            }
            let (name, orient) = tok.split_at(tok.len() - 1);
            match orient {
                "+" => Some((name, true)),
                "-" => Some((name, false)),
                _ => None,
            }
        })
        .collect()
}

/// Parse a GFA v1.1 walk string such as `>12<34>56` into `(segment_name, is_forward)` steps.
fn parse_walk_steps(walk: &str) -> Vec<(&str, bool)> {
    let bytes = walk.as_bytes();
    let mut steps = Vec::new();
    let mut i = 0usize;
    while i < bytes.len() {
        let dir = bytes[i];
        if dir != b'>' && dir != b'<' {
            i += 1;
            continue;
        }
        let is_fwd = dir == b'>';
        let start = i + 1;
        let mut end = start;
        while end < bytes.len() && bytes[end] != b'>' && bytes[end] != b'<' {
            end += 1;
        }
        if end > start {
            steps.push((&walk[start..end], is_fwd));
        }
        i = end;
    }
    steps
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

fn cmd_stats(gfa_path: &str) -> Result<()> {
    let t0 = Instant::now();
    let mut s_count: u64 = 0;
    let mut s_bases: u64 = 0;
    let mut l_count: u64 = 0;
    let mut p_count: u64 = 0;
    let mut w_count: u64 = 0;

    stream_gfa(
        gfa_path,
        |s| {
            s_count += 1;
            s_bases += s.sequence.len() as u64;
            Ok(())
        },
        |_| {
            l_count += 1;
            Ok(())
        },
        |_| {
            p_count += 1;
            Ok(())
        },
        |_| {
            w_count += 1;
            Ok(())
        },
    )?;

    let elapsed = t0.elapsed().as_secs_f64();
    let hwm_mb = read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    println!(
        "gfa={} elapsed={:.2}s peak_rss_mb={:.1} segments={} segment_bp={} links={} paths={} walks={}",
        gfa_path, elapsed, hwm_mb, s_count, s_bases, l_count, p_count, w_count
    );
    Ok(())
}

fn cmd_list_walks(gfa_path: &str, filter: Option<&str>, limit: usize) -> Result<()> {
    let mut shown_p = 0usize;
    let mut shown_w = 0usize;
    stream_gfa(
        gfa_path,
        |_| Ok(()),
        |_| Ok(()),
        |p| {
            if shown_p < limit && filter.map_or(true, |f| p.name.contains(f)) {
                println!("P\t{}\tsteps_bytes={}", p.name, p.segment_names.len());
                shown_p += 1;
            }
            Ok(())
        },
        |w| {
            if shown_w < limit
                && filter.map_or(true, |f| w.sample.contains(f) || w.seq_id.contains(f))
            {
                println!(
                    "W\t{}\t{}\t{}\t{}\t{}",
                    w.sample, w.haplotype, w.seq_id, w.start, w.end
                );
                shown_w += 1;
            }
            Ok(())
        },
    )?;
    Ok(())
}

struct ExtractSpec {
    sample_prefix: String,
    seq_id: String,
    out_fa: String,
}

fn cmd_extract(gfa_path: &str, specs: &[ExtractSpec], max_bp: Option<usize>) -> Result<()> {
    let t0 = Instant::now();
    // Pass 1: find matching paths (P) or walks (W) and collect only the needed segment IDs.
    let mut matched_p: Vec<(usize, String, Vec<(String, bool)>)> = Vec::new();
    let mut matched_w: Vec<(usize, String, Vec<(String, bool)>)> = Vec::new();

    stream_gfa(
        gfa_path,
        |_| Ok(()),
        |_| Ok(()),
        |p| {
            if let Some(idx) = specs.iter().position(|spec| {
                p.name == format!("{}.{}", spec.sample_prefix, spec.seq_id)
                    || (p.name.starts_with(&spec.sample_prefix) && p.name.ends_with(&spec.seq_id))
            }) {
                let steps: Vec<(String, bool)> = parse_path_steps(&p.segment_names)
                    .into_iter()
                    .map(|(s, fwd)| (s.to_string(), fwd))
                    .collect();
                matched_p.push((idx, p.name, steps));
            }
            Ok(())
        },
        |w| {
            if let Some(idx) = specs
                .iter()
                .position(|spec| w.sample.starts_with(&spec.sample_prefix) && w.seq_id == spec.seq_id)
            {
                let steps: Vec<(String, bool)> = parse_walk_steps(&w.walk)
                    .into_iter()
                    .map(|(s, fwd)| (s.to_string(), fwd))
                    .collect();
                let header = format!("{}.{}:{}-{}", w.sample, w.seq_id, w.start, w.end);
                matched_w.push((idx, header, steps));
            }
            Ok(())
        },
    )?;

    let mut matched_paths = matched_p;
    matched_paths.extend(matched_w);
    let mut needed_segs: HashSet<String> = HashSet::new();
    for (_, _, steps) in &matched_paths {
        for (s, _) in steps {
            needed_segs.insert(s.clone());
        }
    }

    if matched_paths.is_empty() {
        bail!("no paths or walks matched the requested sample/seq_id specs");
    }

    // Pass 2: stream S lines and store only the segments referenced by matched paths/walks.
    let mut seg_map: HashMap<String, Vec<u8>> = HashMap::with_capacity(needed_segs.len());
    stream_gfa(
        gfa_path,
        |s| {
            if needed_segs.contains(&s.name) {
                seg_map.insert(s.name, s.sequence);
            }
            Ok(())
        },
        |_| Ok(()),
        |_| Ok(()),
        |_| Ok(()),
    )?;

    // Write one FASTA per spec (first matching longest/first record).
    for (idx, spec) in specs.iter().enumerate() {
        let mut seq_out: Vec<u8> = Vec::new();
        let mut header_id = String::new();
        if let Some((_, hdr, steps)) = matched_paths
            .iter()
            .filter(|(w_idx, _, _)| *w_idx == idx)
            .max_by_key(|(_, _, steps)| steps.len())
        {
            header_id = hdr.clone();
            for (seg, is_fwd) in steps {
                let seg_seq = seg_map
                    .get(seg)
                    .with_context(|| format!("missing segment '{}' in GFA", seg))?;
                if *is_fwd {
                    seq_out.extend_from_slice(seg_seq);
                } else {
                    rev_comp_into(seg_seq, &mut seq_out);
                }
                if let Some(cap) = max_bp {
                    if seq_out.len() >= cap {
                        seq_out.truncate(cap);
                        break;
                    }
                }
            }
        }

        let mut writer = BufWriter::new(File::create(&spec.out_fa)?);
        writeln!(writer, ">{}", header_id)?;
        for chunk in seq_out.chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
        writer.flush()?;
        println!(
            "wrote {} ({} bp, id={})",
            spec.out_fa,
            seq_out.len(),
            header_id
        );
    }

    let hwm_mb = read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    println!(
        "extract_done elapsed={:.2}s peak_rss_mb={:.1} cached_segments={}",
        t0.elapsed().as_secs_f64(),
        hwm_mb,
        seg_map.len()
    );
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        bail!(
            "usage:\n  horizon_gfa stats <file.gfa[.gz]|->\n  horizon_gfa list-walks <file.gfa[.gz]> [filter] [limit]\n  horizon_gfa extract <file.gfa[.gz]> [--max-bp N] <sample_prefix:seq_id:out.fa>..."
        );
    }
    match args[1].as_str() {
        "stats" => cmd_stats(&args[2]),
        "list-walks" => {
            let filter = args.get(3).map(String::as_str);
            let limit = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(50);
            cmd_list_walks(&args[2], filter, limit)
        }
        "extract" => {
            let gfa_path = &args[2];
            let mut max_bp = None;
            let mut specs = Vec::new();
            let mut i = 3;
            while i < args.len() {
                if args[i] == "--max-bp" && i + 1 < args.len() {
                    max_bp = Some(args[i + 1].parse()?);
                    i += 2;
                    continue;
                }
                let parts: Vec<&str> = args[i].split(':').collect();
                if parts.len() != 3 {
                    bail!("extract spec must be sample_prefix:seq_id:out.fa, got {}", args[i]);
                }
                specs.push(ExtractSpec {
                    sample_prefix: parts[0].to_string(),
                    seq_id: parts[1].to_string(),
                    out_fa: parts[2].to_string(),
                });
                i += 1;
            }
            cmd_extract(gfa_path, &specs, max_bp)
        }
        other => bail!("unknown subcommand: {}", other),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_walk_steps_and_reverse_complements() {
        let steps = parse_walk_steps(">10<22>3");
        assert_eq!(steps, vec![("10", true), ("22", false), ("3", true)]);

        let mut out = Vec::new();
        rev_comp_into(b"ACGt", &mut out);
        assert_eq!(out, b"aCGT");
    }
}
