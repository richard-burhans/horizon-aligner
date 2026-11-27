//! Pure-Rust streaming HAL toolkit (`horizon_hal` / `horizon hal`).
//!
//! Zero-`libhdf5` static binary supporting multi-GiB HDF5 v0/v1 `.hal` alignments:
//! - `stats`: guide tree + per-genome sequence/bp/segment counts (matching `halStats`)
//! - `liftover`: multi-level BED coordinate liftover across arbitrary phylogenetic tree paths
//! - `to-maf`: multi-level MAF v1 projection for a full reference genome or genomic locus
//! - `to-gfa`: multi-level base-level GFA v1.1 variation graph projection (`S` + `L` `0M` + `W`)

use anyhow::{Context, Result, bail};
use horizon_aligner::hal::StreamingHalFile;
use horizon_aligner::io_util::{open_text_reader, open_text_writer};
use horizon_aligner::tba::{TbaBlock, TbaRow, VariationGfaBuilder};
use std::collections::HashSet;
use std::env;
use std::io::{BufRead, Write};
use std::path::Path;
use std::time::Instant;

const USAGE: &str = "\
usage:
  horizon_hal stats <file.hal> [--chroms <genome>]
  horizon_hal liftover <file.hal> --src <genome> --dst <genome> --bed <in.bed> [--out <out.bed>] [--cache-chunks <N>]
  horizon_hal to-maf <file.hal> --ref <genome> [--seq <name>] [--start <N>] [--length <N>] [--targets <g1,g2,...>] [--out <out.maf[.gz]>] [--cache-chunks <N>]
  horizon_hal to-gfa <file.hal> --ref <genome> [--seq <name>] [--start <N>] [--length <N>] [--targets <g1,g2,...>] [--out <out.gfa[.gz]>] [--cache-chunks <N>]";

fn next_val<'a>(args: &'a [String], i: usize, flag: &str) -> Result<&'a str> {
    args.get(i + 1)
        .map(String::as_str)
        .with_context(|| format!("missing value for {flag}"))
}

fn cmd_stats(args: &[String]) -> Result<()> {
    let hal_path = args.first().context("stats requires <file.hal>")?;
    let mut chroms_genome: Option<String> = None;
    let mut i = 1usize;
    while i < args.len() {
        match args[i].as_str() {
            "--chroms" => {
                chroms_genome = Some(next_val(args, i, "--chroms")?.to_owned());
                i += 2;
            }
            other => bail!("unknown flag `{other}` for stats"),
        }
    }
    let mut hal = StreamingHalFile::open(Path::new(hal_path))
        .with_context(|| format!("opening HAL `{hal_path}`"))?;
    if let Some(gname) = chroms_genome {
        let seqs = hal
            .genome_sequences(&gname)?
            .with_context(|| format!("genome `{gname}` not found in HAL"))?;
        for s in seqs {
            println!("{}\t{}", s.name, s.length);
        }
        return Ok(());
    }
    println!("{}", hal.tree_newick());
    println!();
    println!(
        "GenomeName\tNumChildren\tLength\tNumSequences\tNumTopSegments\tNumBottomSegments\tParent"
    );
    for st in hal.stats()? {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            st.name,
            st.num_children,
            st.total_bp,
            st.num_sequences,
            st.num_top_segments,
            st.num_bottom_segments,
            st.parent.as_deref().unwrap_or("-")
        );
    }
    Ok(())
}

fn cmd_liftover(args: &[String]) -> Result<()> {
    let hal_path = args.first().context("liftover requires <file.hal>")?;
    let mut src_genome: Option<String> = None;
    let mut dst_genome: Option<String> = None;
    let mut bed_path: Option<String> = None;
    let mut out_path: Option<String> = None;
    let mut cache_chunks: Option<usize> = None;

    let mut i = 1usize;
    while i < args.len() {
        match args[i].as_str() {
            "--src" => {
                src_genome = Some(next_val(args, i, "--src")?.to_owned());
                i += 2;
            }
            "--dst" => {
                dst_genome = Some(next_val(args, i, "--dst")?.to_owned());
                i += 2;
            }
            "--bed" => {
                bed_path = Some(next_val(args, i, "--bed")?.to_owned());
                i += 2;
            }
            "--out" | "-o" => {
                out_path = Some(next_val(args, i, "--out")?.to_owned());
                i += 2;
            }
            "--cache-chunks" => {
                cache_chunks = Some(
                    next_val(args, i, "--cache-chunks")?
                        .parse()
                        .context("parsing --cache-chunks")?,
                );
                i += 2;
            }
            other => bail!("unknown flag `{other}` for liftover"),
        }
    }
    let src = src_genome.context("liftover requires --src <genome>")?;
    let dst = dst_genome.context("liftover requires --dst <genome>")?;
    let bed = bed_path.context("liftover requires --bed <in.bed>")?;

    let t0 = Instant::now();
    let mut hal = StreamingHalFile::open(Path::new(hal_path))
        .with_context(|| format!("opening HAL `{hal_path}`"))?;
    if let Some(cap) = cache_chunks {
        hal.set_chunk_cache_cap(cap);
    }
    let reader =
        open_text_reader(Path::new(&bed)).with_context(|| format!("opening BED `{bed}`"))?;
    let mut writer = open_text_writer(out_path.as_deref())?;

    let mut in_rows = 0usize;
    let mut out_rows = 0usize;
    for line_res in reader.lines() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track ")
            || trimmed.starts_with("browser ")
        {
            continue;
        }
        let cols: Vec<&str> = trimmed.split_ascii_whitespace().collect();
        if cols.len() < 3 {
            bail!("malformed BED row (expected >= 3 columns): `{trimmed}`");
        }
        let seq = cols[0];
        let start: u64 = cols[1]
            .parse()
            .with_context(|| format!("invalid BED start `{}`", cols[1]))?;
        let end: u64 = cols[2]
            .parse()
            .with_context(|| format!("invalid BED end `{}`", cols[2]))?;
        let name = cols.get(3).copied().unwrap_or(".");
        let score = cols.get(4).copied().unwrap_or("0");
        let in_rev = cols.get(5).copied() == Some("-");
        in_rows += 1;

        let hits = hal.liftover_interval(&src, seq, start, end, &dst)?;
        for hit in hits {
            let final_strand = if in_rev {
                if hit.strand == '+' { '-' } else { '+' }
            } else {
                hit.strand
            };
            writeln!(
                writer,
                "{}\t{}\t{}\t{name}\t{score}\t{final_strand}",
                hit.seq_name, hit.start, hit.end
            )?;
            out_rows += 1;
        }
    }
    writer.finish()?;
    eprintln!(
        "[horizon_hal liftover] {src} -> {dst}: in_intervals={in_rows} lifted_intervals={out_rows} elapsed_s={:.3}",
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

struct ProjectArgs {
    hal_path: String,
    ref_genome: String,
    seq: Option<String>,
    start: u64,
    length: Option<u64>,
    targets: Option<Vec<String>>,
    out_path: Option<String>,
    cache_chunks: Option<usize>,
}

fn parse_project_args(args: &[String], subcmd: &str) -> Result<ProjectArgs> {
    let hal_path = args
        .first()
        .with_context(|| format!("{subcmd} requires <file.hal>"))?
        .clone();
    let mut ref_genome: Option<String> = None;
    let mut seq: Option<String> = None;
    let mut start = 0u64;
    let mut length: Option<u64> = None;
    let mut targets: Option<Vec<String>> = None;
    let mut out_path: Option<String> = None;
    let mut cache_chunks: Option<usize> = None;

    let mut i = 1usize;
    while i < args.len() {
        match args[i].as_str() {
            "--ref" => {
                ref_genome = Some(next_val(args, i, "--ref")?.to_owned());
                i += 2;
            }
            "--seq" => {
                seq = Some(next_val(args, i, "--seq")?.to_owned());
                i += 2;
            }
            "--start" => {
                start = next_val(args, i, "--start")?
                    .parse()
                    .context("parsing --start")?;
                i += 2;
            }
            "--length" => {
                length = Some(
                    next_val(args, i, "--length")?
                        .parse()
                        .context("parsing --length")?,
                );
                i += 2;
            }
            "--targets" => {
                let raw = next_val(args, i, "--targets")?;
                targets = Some(
                    raw.split(',')
                        .map(str::trim)
                        .filter(|s| !s.is_empty())
                        .map(ToOwned::to_owned)
                        .collect(),
                );
                i += 2;
            }
            "--out" | "-o" => {
                out_path = Some(next_val(args, i, "--out")?.to_owned());
                i += 2;
            }
            "--cache-chunks" => {
                cache_chunks = Some(
                    next_val(args, i, "--cache-chunks")?
                        .parse()
                        .context("parsing --cache-chunks")?,
                );
                i += 2;
            }
            other => bail!("unknown flag `{other}` for {subcmd}"),
        }
    }
    let ref_genome = ref_genome.with_context(|| format!("{subcmd} requires --ref <genome>"))?;
    anyhow::ensure!(
        seq.is_some() || (start == 0 && length.is_none()),
        "--start and --length require --seq <contig>"
    );
    Ok(ProjectArgs {
        hal_path,
        ref_genome,
        seq,
        start,
        length,
        targets,
        out_path,
        cache_chunks,
    })
}

fn cmd_to_maf(args: &[String]) -> Result<()> {
    let p = parse_project_args(args, "to-maf")?;
    let t0 = Instant::now();
    let mut hal = StreamingHalFile::open(Path::new(&p.hal_path))
        .with_context(|| format!("opening HAL `{}`", p.hal_path))?;
    if let Some(cap) = p.cache_chunks {
        hal.set_chunk_cache_cap(cap);
    }
    let mut writer = open_text_writer(p.out_path.as_deref())?;
    let seq_slice = p.seq.as_deref().map(|s| (s, p.start, p.length));
    let stats = hal.project_to_maf(&p.ref_genome, seq_slice, p.targets.as_deref(), &mut writer)?;
    writer.finish()?;
    eprintln!(
        "[horizon_hal to-maf] ref={} blocks={} ref_bp={} target_rows={} elapsed_s={:.3}",
        p.ref_genome,
        stats.blocks_emitted,
        stats.ref_bp_emitted,
        stats.target_rows_emitted,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

fn cmd_to_gfa(args: &[String]) -> Result<()> {
    let p = parse_project_args(args, "to-gfa")?;
    let t0 = Instant::now();
    let mut hal = StreamingHalFile::open(Path::new(&p.hal_path))
        .with_context(|| format!("opening HAL `{}`", p.hal_path))?;
    if let Some(cap) = p.cache_chunks {
        hal.set_chunk_cache_cap(cap);
    }
    let mut gfa_writer = open_text_writer(p.out_path.as_deref())?;
    let mut builder = VariationGfaBuilder::new(&mut gfa_writer)?;
    let seq_slice = p.seq.as_deref().map(|s| (s, p.start, p.length));
    // Stream projected blocks directly into VariationGfaBuilder with zero intermediate MAF buffer,
    // deduplicating per-species rows within each block to preserve TBA single-occurrence invariants.
    let mstats =
        hal.for_each_projected_block(&p.ref_genome, seq_slice, p.targets.as_deref(), |blk| {
            let mut rows = Vec::with_capacity(1 + blk.target_rows.len());
            let mut seen_sp = HashSet::with_capacity(1 + blk.target_rows.len());
            seen_sp.insert(blk.ref_row.species.clone());
            rows.push(TbaRow {
                species: blk.ref_row.species,
                contig: blk.ref_row.contig,
                start: usize::try_from(blk.ref_row.start).unwrap_or(0),
                size: usize::try_from(blk.ref_row.size).unwrap_or(0),
                strand: blk.ref_row.strand,
                src_size: usize::try_from(blk.ref_row.src_size).unwrap_or(0),
                seq: blk.ref_row.seq,
            });
            for t in blk.target_rows {
                if seen_sp.insert(t.species.clone()) {
                    rows.push(TbaRow {
                        species: t.species,
                        contig: t.contig,
                        start: usize::try_from(t.start).unwrap_or(0),
                        size: usize::try_from(t.size).unwrap_or(0),
                        strand: t.strand,
                        src_size: usize::try_from(t.src_size).unwrap_or(0),
                        seq: t.seq,
                    });
                }
            }
            if rows.len() >= 2 {
                #[allow(clippy::cast_precision_loss)]
                let score = rows.len() as f64;
                builder.push_block(&TbaBlock { score, rows })?;
            }
            Ok(())
        })?;
    let vstats = builder.finish()?;
    gfa_writer.finish()?;
    eprintln!(
        "[horizon_hal to-gfa] ref={} maf_blocks={} segments={} links={} walks={} bp={} elapsed_s={:.3}",
        p.ref_genome,
        mstats.blocks_emitted,
        vstats.segments_emitted,
        vstats.links_emitted,
        vstats.walks_emitted,
        vstats.total_segment_bp,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() || matches!(args[0].as_str(), "-h" | "--help" | "help") {
        println!("{USAGE}");
        return Ok(());
    }
    if matches!(args[0].as_str(), "-V" | "--version") {
        println!("horizon_hal {}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    match args[0].as_str() {
        "stats" => cmd_stats(&args[1..]),
        "liftover" => cmd_liftover(&args[1..]),
        "to-maf" => cmd_to_maf(&args[1..]),
        "to-gfa" => cmd_to_gfa(&args[1..]),
        other => bail!("unknown subcommand `{other}`\n\n{USAGE}"),
    }
}
