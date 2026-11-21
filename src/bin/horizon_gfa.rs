//! Streaming GFA v1.0/v1.1 utility (`stats`, `list-walks`, `extract`) using `horizon_aligner::gfa`.

use anyhow::{Context, Result, bail};
use horizon_aligner::gfa::{
    GfaRecord, parse_path_steps, parse_walk_steps, stream_gfa, stream_gfa_records,
};
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

fn rev_comp_into(seq: &[u8], out: &mut Vec<u8>) {
    out.reserve(seq.len());
    for &b in seq.iter().rev() {
        out.push(horizon_aligner::complement_base(b));
    }
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
            s_bases += s.seq_len() as u64;
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
    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    println!(
        "gfa={gfa_path} elapsed={elapsed:.2}s peak_rss_mb={hwm_mb:.1} segments={s_count} segment_bp={s_bases} links={l_count} paths={p_count} walks={w_count}"
    );
    Ok(())
}

fn cmd_list_walks(gfa_path: &str, filter: Option<&str>, limit: usize) -> Result<()> {
    let mut shown_p = 0usize;
    let mut shown_w = 0usize;
    stream_gfa_records(gfa_path, |rec| {
        match rec {
            GfaRecord::Path(p) if shown_p < limit && filter.is_none_or(|f| p.name.contains(f)) => {
                println!("P\t{}\tsteps_bytes={}", p.name, p.segment_names.len());
                shown_p += 1;
            }
            GfaRecord::Walk(w)
                if shown_w < limit
                    && filter.is_none_or(|f| w.sample.contains(f) || w.seq_id.contains(f)) =>
            {
                println!(
                    "W\t{}\t{}\t{}\t{}\t{}",
                    w.sample, w.haplotype, w.seq_id, w.start, w.end
                );
                shown_w += 1;
            }
            _ => {}
        }
        if shown_p >= limit && shown_w >= limit {
            Ok(std::ops::ControlFlow::Break(()))
        } else {
            Ok(std::ops::ControlFlow::Continue(()))
        }
    })?;
    Ok(())
}

struct ExtractSpec {
    sample_prefix: String,
    seq_id: String,
    out_fa: String,
}

type MatchedPath = (usize, String, Vec<(String, bool)>);

#[allow(clippy::too_many_lines)]
fn cmd_extract(gfa_path: &str, specs: &[ExtractSpec], max_bp: Option<usize>) -> Result<()> {
    let t0 = Instant::now();
    let spooled = horizon_aligner::gfa::SpoolableGfaInput::new(gfa_path)?;
    let eff_path = spooled.path();
    // Pass 1: find matching paths (P) or walks (W) and collect only the needed segment IDs.
    let mut matched_p: Vec<MatchedPath> = Vec::new();
    let mut matched_w: Vec<MatchedPath> = Vec::new();

    stream_gfa(
        eff_path,
        |_| Ok(()),
        |_| Ok(()),
        |p| {
            if let Some(idx) = specs.iter().position(|spec| {
                horizon_aligner::matches_sample_prefix(&p.name, &spec.sample_prefix) && {
                    let rest = horizon_aligner::strip_sample_prefix(&p.name, &spec.sample_prefix);
                    rest == spec.seq_id
                        || rest
                            .strip_suffix(spec.seq_id.as_str())
                            .is_some_and(|prefix| {
                                prefix.ends_with('.')
                                    || prefix.ends_with('#')
                                    || prefix.ends_with(':')
                            })
                }
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
            if let Some(idx) = specs.iter().position(|spec| {
                horizon_aligner::matches_sample_prefix(&w.sample, &spec.sample_prefix)
                    && w.seq_id == spec.seq_id
            }) {
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
    if matched_paths.is_empty() {
        bail!("no paths or walks matched the requested sample/seq_id specs");
    }

    let mut selected_paths = Vec::with_capacity(specs.len());
    for (idx, spec) in specs.iter().enumerate() {
        let Some((_, hdr, steps)) = matched_paths
            .iter()
            .filter(|(w_idx, _, _)| *w_idx == idx)
            .max_by_key(|(_, _, steps)| steps.len())
        else {
            bail!(
                "no path or walk matched spec {}:{}",
                spec.sample_prefix,
                spec.seq_id
            );
        };
        selected_paths.push((hdr.clone(), steps.clone()));
    }

    let mut needed_segs: HashSet<String> = HashSet::new();
    for (_, steps) in &selected_paths {
        for (s, _) in steps {
            needed_segs.insert(s.clone());
        }
    }

    // Pass 2: stream S lines and store only the segments referenced by selected paths/walks.
    let mut seg_map: HashMap<String, Vec<u8>> = HashMap::with_capacity(needed_segs.len());
    stream_gfa(
        eff_path,
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

    // Write one FASTA per spec (longest matching record).
    for (spec, (header_id, steps)) in specs.iter().zip(selected_paths.iter()) {
        let mut seq_out: Vec<u8> = Vec::new();
        if max_bp != Some(0) {
            for (seg, is_fwd) in steps {
                let seg_seq = seg_map
                    .get(seg)
                    .with_context(|| format!("missing segment '{seg}' in GFA"))?;
                if seg_seq == b"*" {
                    bail!(
                        "segment '{seg}' has placeholder sequence '*' and cannot be extracted to FASTA"
                    );
                }
                if *is_fwd {
                    seq_out.extend_from_slice(seg_seq);
                } else {
                    rev_comp_into(seg_seq, &mut seq_out);
                }
                if let Some(cap) = max_bp
                    && seq_out.len() >= cap
                {
                    seq_out.truncate(cap);
                    break;
                }
            }
        }

        let mut writer = BufWriter::new(File::create(&spec.out_fa)?);
        writeln!(writer, ">{header_id}")?;
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

    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    println!(
        "extract_done elapsed={:.2}s peak_rss_mb={:.1} cached_segments={}",
        t0.elapsed().as_secs_f64(),
        hwm_mb,
        seg_map.len()
    );
    Ok(())
}

const USAGE: &str = "usage:\n  horizon_gfa stats <file.gfa[.gz]|->\n  horizon_gfa list-walks <file.gfa[.gz]> [filter] [limit]\n  horizon_gfa extract <file.gfa[.gz]> [--max-bp N] <sample_prefix:seq_id:out.fa>...";

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.iter().skip(1).any(|a| a == "-h" || a == "--help") {
        println!("{USAGE}");
        return Ok(());
    }
    if args.len() < 3 {
        bail!("{USAGE}");
    }
    match args[1].as_str() {
        "stats" => {
            if args.len() != 3 {
                bail!("stats accepts only <file.gfa[.gz]|->");
            }
            cmd_stats(&args[2])
        }
        "list-walks" => {
            if args.len() > 5 {
                bail!("list-walks accepts at most <file.gfa[.gz]> [filter] [limit]");
            }
            let filter = args.get(3).map(String::as_str);
            let limit = match args.get(4) {
                Some(s) => s.parse().with_context(|| format!("invalid limit: {s}"))?,
                None => 50,
            };
            cmd_list_walks(&args[2], filter, limit)
        }
        "extract" => {
            let gfa_path = &args[2];
            let mut max_bp = None;
            let mut specs = Vec::new();
            let mut i = 3;
            while i < args.len() {
                if args[i] == "--max-bp" {
                    let val = args
                        .get(i + 1)
                        .ok_or_else(|| anyhow::anyhow!("missing value for --max-bp"))?;
                    max_bp = Some(val.parse().context("parsing --max-bp")?);
                    i += 2;
                    continue;
                }
                let Some((sample_prefix, rest)) = args[i].split_once(':') else {
                    bail!(
                        "extract spec must be sample_prefix:seq_id:out.fa, got {}",
                        args[i]
                    );
                };
                let Some((seq_id, out_fa)) = rest.rsplit_once(':') else {
                    bail!(
                        "extract spec must be sample_prefix:seq_id:out.fa, got {}",
                        args[i]
                    );
                };
                specs.push(ExtractSpec {
                    sample_prefix: sample_prefix.to_string(),
                    seq_id: seq_id.to_string(),
                    out_fa: out_fa.to_string(),
                });
                i += 1;
            }
            if specs.is_empty() {
                bail!("extract requires at least one sample_prefix:seq_id:out.fa spec");
            }
            cmd_extract(gfa_path, &specs, max_bp)
        }
        other => bail!("unknown subcommand: {other}"),
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

    #[test]
    fn stats_list_walks_and_extract_work_with_max_bp_and_star_error() {
        let tmp = std::env::temp_dir();
        let pid = std::process::id();
        let gfa_path = tmp.join(format!("horizon_gfa_unit_test_{pid}.gfa"));
        let out_fa = tmp.join(format!("horizon_gfa_unit_out_{pid}.fa"));
        let gfa = "\
H\tVN:Z:1.1
S\ts1\tAAACCC
S\ts2\tGGGTTT
S\ts3\t*\tLN:i:10
L\ts1\t+\ts2\t+\t0M
W\tHG002\t0\tchr1\t0\t12\t>s1<s2
P\tref#0#chr1\ts1+,s2-\t0M
";
        std::fs::write(&gfa_path, gfa).expect("write gfa");
        let gfa_str = gfa_path.to_str().expect("utf8 path");

        cmd_stats(gfa_str).expect("cmd_stats");
        cmd_list_walks(gfa_str, Some("HG002"), 10).expect("cmd_list_walks");

        let spec = ExtractSpec {
            sample_prefix: "HG002".to_string(),
            seq_id: "chr1".to_string(),
            out_fa: out_fa.to_str().expect("utf8 out path").to_string(),
        };
        cmd_extract(gfa_str, &[spec], Some(9)).expect("cmd_extract");
        let extracted = std::fs::read_to_string(&out_fa).expect("read out_fa");
        assert_eq!(extracted, ">HG002.chr1:0-12\nAAACCCAAA\n");

        let star_gfa_path = tmp.join(format!("horizon_gfa_star_test_{pid}.gfa"));
        std::fs::write(
            &star_gfa_path,
            "H\tVN:Z:1.1\nS\ts3\t*\tLN:i:10\nW\tHG002\t0\tchr1\t0\t10\t>s3\n",
        )
        .expect("write star gfa");
        let star_spec = ExtractSpec {
            sample_prefix: "HG002".to_string(),
            seq_id: "chr1".to_string(),
            out_fa: out_fa.to_str().expect("utf8 out path").to_string(),
        };
        assert!(
            cmd_extract(
                star_gfa_path.to_str().expect("utf8 star path"),
                &[star_spec],
                None
            )
            .is_err()
        );

        let _ = std::fs::remove_file(&gfa_path);
        let _ = std::fs::remove_file(&star_gfa_path);
        let _ = std::fs::remove_file(&out_fa);
    }
}
