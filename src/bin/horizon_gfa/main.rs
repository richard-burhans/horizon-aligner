//! Streaming GFA v1.0/v1.1 utility (`stats`, `qc`, `list-walks`, `extract`, `rgfa`, `from-maf`) using `horizon_aligner::gfa`.

mod qc;

use anyhow::{Context, Result, bail};
use horizon_aligner::gfa::{
    GfaRecordRef, parse_path_steps, parse_walk_steps, stream_gfa_record_refs,
};
use std::collections::{HashMap, HashSet};
use std::env;
use std::io::Write;
use std::time::Instant;

fn cmd_stats(gfa_path: &str) -> Result<()> {
    let t0 = Instant::now();
    let mut s_count: u64 = 0;
    let mut s_bases: u64 = 0;
    let mut l_count: u64 = 0;
    let mut p_count: u64 = 0;
    let mut w_count: u64 = 0;

    stream_gfa_record_refs(gfa_path, |rec| {
        match rec {
            GfaRecordRef::Segment(s) => {
                s_count += 1;
                s_bases += s.seq_len() as u64;
            }
            GfaRecordRef::Link(_) => l_count += 1,
            GfaRecordRef::Path(_) => p_count += 1,
            GfaRecordRef::Walk(_) => w_count += 1,
        }
        Ok(std::ops::ControlFlow::Continue(()))
    })?;

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
    stream_gfa_record_refs(gfa_path, |rec| {
        match rec {
            GfaRecordRef::Path(p)
                if shown_p < limit && filter.is_none_or(|f| p.name.contains(f)) =>
            {
                println!("P\t{}\tsteps_bytes={}", p.name, p.segment_names.len());
                shown_p += 1;
            }
            GfaRecordRef::Walk(w)
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
        if shown_p + shown_w >= limit {
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

type MatchedPath = (usize, String, Vec<(String, bool)>, usize);

#[allow(clippy::too_many_lines)]
fn cmd_extract(gfa_path: &str, specs: &[ExtractSpec], max_bp: Option<usize>) -> Result<()> {
    let t0 = Instant::now();
    let spooled = horizon_aligner::gfa::SpoolableGfaInput::new(gfa_path)?;
    let eff_path = spooled.path();
    // Pass 1: find matching paths (P) or walks (W) without allocating segment sequences.
    let stdout_specs = specs.iter().filter(|s| s.out_fa == "-").count();
    anyhow::ensure!(
        stdout_specs <= 1,
        "extract accepts at most one spec writing to stdout ('-')"
    );
    let any_stdout = stdout_specs > 0;
    let mut matched_p: Vec<MatchedPath> = Vec::new();
    let mut matched_w: Vec<MatchedPath> = Vec::new();
    let mut seg_lens: HashMap<String, usize> = HashMap::new();

    stream_gfa_record_refs(eff_path, |rec| {
        match rec {
            GfaRecordRef::Segment(s) => {
                seg_lens.insert(s.name.to_owned(), s.seq_len());
            }
            GfaRecordRef::Path(p) => {
                let matching_specs: Vec<usize> = specs
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, spec)| {
                        let matched =
                            horizon_aligner::matches_sample_prefix(p.name, &spec.sample_prefix)
                                && {
                                    let rest = horizon_aligner::strip_sample_prefix(
                                        p.name,
                                        &spec.sample_prefix,
                                    );
                                    rest == spec.seq_id
                                        || rest.strip_suffix(spec.seq_id.as_str()).is_some_and(
                                            |prefix| {
                                                prefix.ends_with('.')
                                                    || prefix.ends_with('#')
                                                    || prefix.ends_with(':')
                                            },
                                        )
                                };
                        matched.then_some(idx)
                    })
                    .collect();
                if !matching_specs.is_empty() {
                    let steps: Vec<(String, bool)> = parse_path_steps(p.segment_names)
                        .into_iter()
                        .map(|(s, fwd)| (s.to_string(), fwd))
                        .collect();
                    for idx in matching_specs {
                        matched_p.push((idx, p.name.to_owned(), steps.clone(), 0));
                    }
                }
            }
            GfaRecordRef::Walk(w) => {
                let matching_specs: Vec<usize> = specs
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, spec)| {
                        let matched =
                            horizon_aligner::matches_sample_prefix(w.sample, &spec.sample_prefix)
                                && w.seq_id == spec.seq_id;
                        matched.then_some(idx)
                    })
                    .collect();
                if !matching_specs.is_empty() {
                    let steps: Vec<(String, bool)> = parse_walk_steps(w.walk)
                        .into_iter()
                        .map(|(s, fwd)| (s.to_string(), fwd))
                        .collect();
                    let header = format!("{}.{}:{}-{}", w.sample, w.seq_id, w.start, w.end);
                    let walk_span = w.end.saturating_sub(w.start);
                    for idx in matching_specs {
                        matched_w.push((idx, header.clone(), steps.clone(), walk_span));
                    }
                }
            }
            GfaRecordRef::Link(_) => {}
        }
        Ok(std::ops::ControlFlow::Continue(()))
    })?;

    let mut matched_paths = matched_p;
    matched_paths.extend(matched_w);
    if matched_paths.is_empty() {
        bail!("no paths or walks matched the requested sample/seq_id specs");
    }

    let mut selected_paths = Vec::with_capacity(specs.len());
    for (idx, spec) in specs.iter().enumerate() {
        let Some((_, hdr, steps, _)) = matched_paths
            .iter()
            .filter(|(w_idx, _, _, _)| *w_idx == idx)
            .max_by_key(|(_, _, steps, walk_span)| {
                let bp = if *walk_span > 0 {
                    *walk_span
                } else {
                    steps
                        .iter()
                        .map(|(s, _)| seg_lens.get(s).copied().unwrap_or(0))
                        .sum()
                };
                (bp, steps.len())
            })
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

    // Pass 2: zero-copy stream S lines and allocate ONLY segments referenced by selected walks.
    let mut seg_map: HashMap<String, Vec<u8>> = HashMap::with_capacity(needed_segs.len());
    stream_gfa_record_refs(eff_path, |rec| {
        if let GfaRecordRef::Segment(s) = rec
            && needed_segs.contains(s.name)
        {
            seg_map.insert(s.name.to_owned(), s.sequence.to_vec());
        }
        Ok(std::ops::ControlFlow::Continue(()))
    })?;

    // Stream one FASTA per spec in 80-bp chunks without buffering full walk sequence in RAM.
    for (spec, (header_id, steps)) in specs.iter().zip(selected_paths.iter()) {
        // Validate every referenced segment exists and is non-placeholder before opening output.
        if max_bp != Some(0) {
            let mut rem_check = max_bp.unwrap_or(usize::MAX);
            for (seg, _) in steps {
                if rem_check == 0 {
                    break;
                }
                let seg_seq = seg_map
                    .get(seg)
                    .with_context(|| format!("missing segment '{seg}' in GFA"))?;
                if seg_seq == b"*" {
                    bail!(
                        "segment '{seg}' has placeholder sequence '*' and cannot be extracted to FASTA"
                    );
                }
                rem_check = rem_check.saturating_sub(seg_seq.len());
            }
        }

        let mut writer = horizon_aligner::io_util::open_text_writer(Some(&spec.out_fa))
            .with_context(|| format!("creating extract output '{}'", spec.out_fa))?;
        writeln!(writer, ">{header_id}")?;

        let mut total_written = 0usize;
        let mut line_buf = [0u8; 80];
        let mut line_pos = 0usize;
        if max_bp != Some(0) {
            let cap = max_bp.unwrap_or(usize::MAX);
            'steps: for (seg, is_fwd) in steps {
                let seg_seq = &seg_map[seg];
                if *is_fwd {
                    for &b in seg_seq {
                        if total_written >= cap {
                            break 'steps;
                        }
                        line_buf[line_pos] = b;
                        line_pos += 1;
                        total_written += 1;
                        if line_pos == 80 {
                            writer.write_all(&line_buf)?;
                            writer.write_all(b"\n")?;
                            line_pos = 0;
                        }
                    }
                } else {
                    for &b in seg_seq.iter().rev() {
                        if total_written >= cap {
                            break 'steps;
                        }
                        line_buf[line_pos] = horizon_aligner::seeding::complement_base(b);
                        line_pos += 1;
                        total_written += 1;
                        if line_pos == 80 {
                            writer.write_all(&line_buf)?;
                            writer.write_all(b"\n")?;
                            line_pos = 0;
                        }
                    }
                }
            }
            if line_pos > 0 {
                writer.write_all(&line_buf[..line_pos])?;
                writer.write_all(b"\n")?;
            }
        }
        writer.finish()?;
        if any_stdout {
            eprintln!(
                "wrote {} ({} bp, id={})",
                spec.out_fa, total_written, header_id
            );
        } else {
            println!(
                "wrote {} ({} bp, id={})",
                spec.out_fa, total_written, header_id
            );
        }
    }

    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    if any_stdout {
        eprintln!(
            "extract_done elapsed={:.2}s peak_rss_mb={:.1} cached_segments={}",
            t0.elapsed().as_secs_f64(),
            hwm_mb,
            seg_map.len()
        );
    } else {
        println!(
            "extract_done elapsed={:.2}s peak_rss_mb={:.1} cached_segments={}",
            t0.elapsed().as_secs_f64(),
            hwm_mb,
            seg_map.len()
        );
    }
    Ok(())
}

const USAGE: &str = "usage:\n  horizon_gfa stats <file.gfa[.gz]|->\n  horizon_gfa qc <file.gfa[.gz]|-> [--max-segments N]\n  horizon_gfa list-walks <file.gfa[.gz]> [filter] [limit]\n  horizon_gfa extract <file.gfa[.gz]> [--max-bp N] <sample_prefix:seq_id:out.fa>...\n  horizon_gfa rgfa --ref <ref.fa> --query <q.fa> [--query <q2.fa>...] [--min-sv-bp 50] [--respect-softmask] [-t|--threads N] [--out <out.rgfa>]\n  horizon_gfa from-maf <in.maf[.gz]|-> [--out <out.gfa[.gz]>]";

fn cmd_from_maf(args: &[String]) -> Result<()> {
    let t0 = Instant::now();
    let mut maf_path: Option<String> = None;
    let mut out_path: Option<String> = None;
    let mut i = 2usize;
    while i < args.len() {
        match args[i].as_str() {
            "--out" | "-o" => {
                out_path = Some(args.get(i + 1).context("missing value for --out")?.clone());
                i += 2;
            }
            other if other.starts_with('-') && other != "-" => {
                bail!("unknown flag in from-maf: {other}");
            }
            other => {
                anyhow::ensure!(
                    maf_path.is_none(),
                    "from-maf accepts only one <in.maf[.gz]|-> argument, got extra `{other}`"
                );
                maf_path = Some(other.to_owned());
                i += 1;
            }
        }
    }
    let maf_path = maf_path.context("from-maf requires <in.maf[.gz]|->")?;
    let reader = horizon_aligner::io_util::open_text_reader(std::path::Path::new(&maf_path))
        .with_context(|| format!("opening MAF input `{maf_path}`"))?;
    let mut writer = horizon_aligner::gfa::open_output_writer(out_path.as_deref())?;
    let stats = horizon_aligner::tba::write_variation_gfa_from_maf_reader(reader, &mut writer)?;
    writer.finish()?;
    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    eprintln!(
        "horizon_gfa from-maf: blocks={} segments={} segment_bp={} links={} walks={} poly_runs={} elapsed={:.2}s peak_rss_mb={:.1}",
        stats.blocks_processed,
        stats.segments_emitted,
        stats.total_segment_bp,
        stats.links_emitted,
        stats.walks_emitted,
        stats.polymorphic_runs,
        t0.elapsed().as_secs_f64(),
        hwm_mb
    );
    Ok(())
}

fn cmd_rgfa(args: &[String]) -> Result<()> {
    let mut ref_fa: Option<String> = None;
    let mut query_fas: Vec<String> = Vec::new();
    let mut out_path: Option<String> = None;
    let mut min_sv_bp = 50usize;
    let mut respect_softmask = false;
    let mut threads: Option<usize> = None;
    let mut i = 2usize;
    while i < args.len() {
        match args[i].as_str() {
            "--ref" => {
                ref_fa = Some(args.get(i + 1).context("missing value for --ref")?.clone());
                i += 2;
            }
            "--query" | "-q" => {
                query_fas.push(
                    args.get(i + 1)
                        .context("missing value for --query")?
                        .clone(),
                );
                i += 2;
            }
            "--min-sv-bp" => {
                min_sv_bp = args
                    .get(i + 1)
                    .context("missing value for --min-sv-bp")?
                    .parse()
                    .context("parsing --min-sv-bp")?;
                anyhow::ensure!(min_sv_bp >= 1, "--min-sv-bp must be >= 1 (got {min_sv_bp})");
                i += 2;
            }
            "--respect-softmask" => {
                respect_softmask = true;
                i += 1;
            }
            "-t" | "--threads" => {
                let raw = args.get(i + 1).context("missing value for --threads")?;
                threads = Some(horizon_aligner::parse_threads_arg(raw)?);
                i += 2;
            }
            "--out" | "-o" => {
                out_path = Some(args.get(i + 1).context("missing value for --out")?.clone());
                i += 2;
            }
            other => bail!("unknown flag in rgfa: {other}"),
        }
    }
    let ref_fa = ref_fa.context("rgfa requires --ref <ref.fa>")?;
    horizon_aligner::init_rayon_threads(threads)?;
    let ref_ctgs = horizon_aligner::read_fasta_contigs(std::path::Path::new(&ref_fa), 0)?;
    let ref_refs: Vec<(&str, &[u8])> = ref_ctgs
        .iter()
        .map(|(n, s)| (n.as_str(), s.as_slice()))
        .collect();
    let mut graph = horizon_aligner::RgfaGraph::from_reference(&ref_refs, min_sv_bp)
        .with_respect_softmask(respect_softmask);
    let mut total_svs = 0usize;
    for q_path in &query_fas {
        let q_ctgs = horizon_aligner::read_fasta_contigs(std::path::Path::new(q_path), 0)?;
        let q_refs: Vec<(&str, &[u8])> = q_ctgs
            .iter()
            .map(|(n, s)| (n.as_str(), s.as_slice()))
            .collect();
        total_svs += graph.graft_assembly_contigs(&q_refs);
    }
    let mut writer = horizon_aligner::gfa::open_output_writer(out_path.as_deref())?;
    graph.write_rgfa(&mut writer)?;
    writer.finish()?;
    eprintln!(
        "horizon_gfa rgfa: segments={} links={} svs_grafted={} ranks={}",
        graph.segments.len(),
        graph.links.len(),
        total_svs,
        graph.next_rank
    );
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if let Some(first) = args.get(1).map(String::as_str) {
        if first == "-h" || first == "--help" {
            println!("{USAGE}");
            return Ok(());
        }
        if first == "-V" || first == "--version" {
            println!("horizon_gfa {}", env!("CARGO_PKG_VERSION"));
            return Ok(());
        }
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
        "qc" => qc::cmd_qc(&args),
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
            let mut gfa_path: Option<&str> = None;
            let mut max_bp = None;
            let mut specs = Vec::new();
            let mut i = 2;
            while i < args.len() {
                if args[i] == "--max-bp" {
                    let val = args
                        .get(i + 1)
                        .ok_or_else(|| anyhow::anyhow!("missing value for --max-bp"))?;
                    max_bp = Some(val.parse().context("parsing --max-bp")?);
                    i += 2;
                    continue;
                }
                if args[i].starts_with('-') && args[i] != "-" {
                    bail!("unknown flag in extract: {}", args[i]);
                }
                if gfa_path.is_none() {
                    gfa_path = Some(&args[i]);
                    i += 1;
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
                anyhow::ensure!(
                    !sample_prefix.is_empty() && !seq_id.is_empty() && !out_fa.is_empty(),
                    "extract spec fields must be non-empty (sample_prefix:seq_id:out.fa), got {}",
                    args[i]
                );
                specs.push(ExtractSpec {
                    sample_prefix: sample_prefix.to_string(),
                    seq_id: seq_id.to_string(),
                    out_fa: out_fa.to_string(),
                });
                i += 1;
            }
            let gfa_path = gfa_path
                .ok_or_else(|| anyhow::anyhow!("missing <file.gfa> argument for extract"))?;
            if specs.is_empty() {
                bail!("extract requires at least one sample_prefix:seq_id:out.fa spec");
            }
            cmd_extract(gfa_path, &specs, max_bp)
        }
        "rgfa" => cmd_rgfa(&args),
        "from-maf" => cmd_from_maf(&args),
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

        let out: Vec<u8> = b"ACGt"
            .iter()
            .rev()
            .map(|&b| horizon_aligner::complement_base(b))
            .collect();
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

    #[test]
    fn from_maf_emits_extractable_variation_gfa() {
        let tmp = std::env::temp_dir();
        let pid = std::process::id();
        let maf_path = tmp.join(format!("horizon_gfa_from_maf_{pid}.maf"));
        let gfa_path = tmp.join(format!("horizon_gfa_from_maf_{pid}.gfa"));
        let out_fa = tmp.join(format!("horizon_gfa_from_maf_{pid}.fa"));
        let maf = "\
##maf version=1
a score=20.0
s sp1.chr1 0 9 + 100 ACGTACGTT
s sp2.chr1 0 8 + 100 ACG-GCGTT

";
        std::fs::write(&maf_path, maf).expect("write maf");
        let args = vec![
            "horizon_gfa".to_owned(),
            "from-maf".to_owned(),
            maf_path.to_str().expect("utf8").to_owned(),
            "--out".to_owned(),
            gfa_path.to_str().expect("utf8").to_owned(),
        ];
        cmd_from_maf(&args).expect("cmd_from_maf");

        let spec = ExtractSpec {
            sample_prefix: "sp2".to_owned(),
            seq_id: "chr1".to_owned(),
            out_fa: out_fa.to_str().expect("utf8").to_owned(),
        };
        cmd_extract(gfa_path.to_str().expect("utf8"), &[spec], None).expect("extract sp2");
        let extracted = std::fs::read_to_string(&out_fa).expect("read extracted");
        assert_eq!(extracted, ">sp2.chr1:0-8\nACGGCGTT\n");

        let _ = std::fs::remove_file(&maf_path);
        let _ = std::fs::remove_file(&gfa_path);
        let _ = std::fs::remove_file(&out_fa);
    }
}
