//! Multi-contig whole-genome aligner CLI using KP-14 + parallel MDL/BOCPD.

use anyhow::{Context, Result, bail};
use horizon_aligner::gfa::{GfaRecord, stream_gfa_records};
use horizon_aligner::{
    AlignConfig, ContigIndex, DagSeedIndex, GraphNode, HmmParams, align_contig_pair_runs,
    align_seq_to_dag_seeded_indexed, format_gfa_blocks_runs, format_gfa_paralog_tags,
    format_gfa_summary_runs, format_maf_blocks_runs, net_one_to_one_runs,
    prepare_contig_index_opts, read_fasta_contigs, resolve_alignment_duplications_runs, revcomp,
    select_top_k_pairs, shared_seed_count,
};
use rayon::prelude::*;
use std::collections::{HashMap, VecDeque};
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

#[allow(clippy::struct_excessive_bools)]
struct CliArgs {
    fa_x: PathBuf,
    fa_y: PathBuf,
    max_contigs: usize,
    maf_out: Option<PathBuf>,
    gfa_out: Option<PathBuf>,
    resolve_dups: bool,
    one_to_one: bool,
    min_shared_seeds: usize,
    min_top_lcb: usize,
    top_k: usize,
    respect_softmask: bool,
    use_ry_seeds: bool,
    hal_info: Option<PathBuf>,
    graph_path: Option<PathBuf>,
    min_fallback_cov: Option<f64>,
    hsp_threshold: Option<i32>,
}

#[allow(clippy::too_many_lines)]
fn parse_args() -> Result<Option<CliArgs>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut pos = Vec::new();
    let mut maf_out = None;
    let mut gfa_out = None;
    let mut resolve_dups = false;
    let mut one_to_one = false;
    let mut min_shared_seeds = 8usize;
    let mut min_top_lcb = 0usize;
    let mut top_k = 0usize;
    let mut respect_softmask = false;
    let mut use_ry_seeds = false;
    let mut hal_info = None;
    let mut graph_path = None;
    let mut min_fallback_cov = None;
    let mut hsp_threshold = None;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                println!(
                    "usage: horizon_wga <target.fa> <query.fa> [max_contigs] [--maf out.maf] [--gfa out.gfa] [--min-shared-seeds N] [--min-top-lcb N] [--top-k K | --rbest] [--resolve-dups] [--one-to-one] [--respect-softmask] [--ry-seeds | --no-ry-seeds] [--min-fallback-cov F] [--hsp-threshold N]\n       (max_contigs defaults to 1; pass 0 to align all contigs)\n       horizon_wga --hal-info <file.hal>\n       horizon_wga --graph <graph.gfa> <query.fa> [max_contigs]"
                );
                return Ok(None);
            }
            "--maf" => {
                maf_out = Some(PathBuf::from(
                    args.get(i + 1).context("missing value for --maf")?,
                ));
                i += 2;
            }
            "--gfa" => {
                gfa_out = Some(PathBuf::from(
                    args.get(i + 1).context("missing value for --gfa")?,
                ));
                i += 2;
            }
            "--min-shared-seeds" => {
                min_shared_seeds = args
                    .get(i + 1)
                    .context("missing value for --min-shared-seeds")?
                    .parse()
                    .context("parsing --min-shared-seeds")?;
                i += 2;
            }
            "--min-top-lcb" => {
                min_top_lcb = args
                    .get(i + 1)
                    .context("missing value for --min-top-lcb")?
                    .parse()
                    .context("parsing --min-top-lcb")?;
                i += 2;
            }
            "--top-k" => {
                top_k = args
                    .get(i + 1)
                    .context("missing value for --top-k")?
                    .parse()
                    .context("parsing --top-k")?;
                i += 2;
            }
            "--resolve-dups" => {
                resolve_dups = true;
                i += 1;
            }
            "--one-to-one" => {
                one_to_one = true;
                i += 1;
            }
            "--rbest" => {
                top_k = 1;
                i += 1;
            }
            "--respect-softmask" => {
                respect_softmask = true;
                i += 1;
            }
            "--ry-seeds" => {
                use_ry_seeds = true;
                i += 1;
            }
            "--no-ry-seeds" => {
                use_ry_seeds = false;
                i += 1;
            }
            "--min-fallback-cov" => {
                let raw = args
                    .get(i + 1)
                    .context("missing value for --min-fallback-cov")?;
                min_fallback_cov = Some(horizon_aligner::parse_min_fallback_cov_arg(raw)?);
                i += 2;
            }
            "--hsp-threshold" => {
                let raw = args
                    .get(i + 1)
                    .context("missing value for --hsp-threshold")?;
                hsp_threshold = Some(horizon_aligner::parse_hsp_threshold_arg(raw)?);
                i += 2;
            }
            "--hal-info" => {
                hal_info = Some(PathBuf::from(
                    args.get(i + 1).context("missing value for --hal-info")?,
                ));
                i += 2;
            }
            "--graph" => {
                graph_path = Some(PathBuf::from(
                    args.get(i + 1).context("missing value for --graph")?,
                ));
                i += 2;
            }
            other if other.starts_with('-') => {
                bail!("unknown flag '{other}'");
            }
            _ => {
                pos.push(args[i].clone());
                i += 1;
            }
        }
    }

    if hal_info.is_some() && graph_path.is_some() {
        bail!("--hal-info and --graph cannot be combined");
    }
    if pos.len() < 2 && hal_info.is_none() && (graph_path.is_none() || pos.is_empty()) {
        bail!(
            "usage: horizon_wga <target.fa> <query.fa> [max_contigs] [--maf out.maf] [--gfa out.gfa] [--graph target.gfa] [--hal-info file.hal] [--min-shared-seeds N] [--min-top-lcb N] [--top-k K | --rbest] [--resolve-dups] [--one-to-one] [--respect-softmask] [--ry-seeds | --no-ry-seeds] [--min-fallback-cov F] [--hsp-threshold N]"
        );
    }
    let (fa_x, fa_y, max_contigs) = if hal_info.is_some() {
        if !pos.is_empty() {
            bail!("--hal-info does not accept positional arguments");
        }
        (PathBuf::new(), PathBuf::new(), 1)
    } else if graph_path.is_some() {
        if pos.len() > 2 {
            bail!(
                "--graph mode accepts <query.fa> [max_contigs], got {} positional args",
                pos.len()
            );
        }
        if maf_out.is_some() || gfa_out.is_some() {
            bail!("--maf and --gfa outputs are not supported in --graph mode");
        }
        let mc = match pos.get(1) {
            Some(s) => s.parse().context("parsing max_contigs")?,
            None => 1,
        };
        (PathBuf::new(), PathBuf::from(&pos[0]), mc)
    } else {
        if pos.len() > 3 {
            bail!(
                "expected <target.fa> <query.fa> [max_contigs], got {} positional args",
                pos.len()
            );
        }
        let mc = match pos.get(2) {
            Some(s) => s.parse().context("parsing max_contigs")?,
            None => 1,
        };
        (
            pos.first().map_or_default(PathBuf::from),
            pos.get(1).map_or_default(PathBuf::from),
            mc,
        )
    };

    Ok(Some(CliArgs {
        fa_x,
        fa_y,
        max_contigs,
        maf_out,
        gfa_out,
        resolve_dups,
        one_to_one,
        min_shared_seeds,
        min_top_lcb,
        top_k,
        respect_softmask,
        use_ry_seeds,
        hal_info,
        graph_path,
        min_fallback_cov,
        hsp_threshold,
    }))
}

#[allow(clippy::too_many_lines)]
fn main() -> Result<()> {
    let Some(cli) = parse_args()? else {
        return Ok(());
    };
    if let Some(hal_path) = &cli.hal_info {
        let info = horizon_aligner::hal::inspect_hal_genomes(hal_path)
            .with_context(|| format!("inspecting {}", hal_path.display()))?
            .with_context(|| format!("{} is not a valid HDF5 HAL file", hal_path.display()))?;
        println!(
            "hal_file={} hdf5_v={} offset_bytes={} len_bytes={} root_members={}",
            hal_path.display(),
            info.superblock.version,
            info.superblock.offset_size,
            info.superblock.length_size,
            info.root_names.len()
        );
        if !info.root_names.is_empty() {
            println!(
                "sample_members={}",
                info.root_names[..info.root_names.len().min(10)].join(",")
            );
        }
        return Ok(());
    }
    if let Some(graph_p) = &cli.graph_path {
        let t0 = Instant::now();
        let dag_nodes = load_gfa_dag(graph_p)?;
        let dag_index = DagSeedIndex::build(&dag_nodes);
        let queries = read_fasta_contigs(&cli.fa_y, cli.max_contigs)?;
        let mut total_pairs = 0usize;
        let mut total_bubbles = 0usize;
        for (qname, qseq) in &queries {
            let dag_params = estimate_dag_hmm_params(&dag_nodes, qseq);
            let aln =
                align_seq_to_dag_seeded_indexed(&dag_nodes, &dag_index, qseq, &dag_params, 10.0);
            total_pairs += aln.kept_pairs.len();
            total_bubbles += aln.bubbles;
            let path_str = aln.path_nodes.iter().fold(String::new(), |mut acc, &idx| {
                let _ = write!(acc, ">{}", dag_nodes[idx].name);
                acc
            });
            println!(
                "{}\t{}\t{}\tkept_bp={}\tscore_bits={:.1}\tbubbles={}",
                qname,
                qseq.len(),
                path_str,
                aln.kept_pairs.len(),
                aln.score_bits,
                aln.bubbles
            );
        }
        println!(
            "graph_nodes={} query_contigs={} aligned_pairs={} total_bubbles={} elapsed_s={:.4}",
            dag_nodes.len(),
            queries.len(),
            total_pairs,
            total_bubbles,
            t0.elapsed().as_secs_f64()
        );
        return Ok(());
    }
    let t0 = Instant::now();
    let contigs_x = read_fasta_contigs(&cli.fa_x, cli.max_contigs)?;
    let contigs_y = read_fasta_contigs(&cli.fa_y, cli.max_contigs)?;

    let idx_x: Vec<(String, ContigIndex)> = contigs_x
        .par_iter()
        .map(|(n, s)| {
            Ok((
                n.clone(),
                prepare_contig_index_opts(s, cli.respect_softmask)?,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    let idx_y: Vec<(String, ContigIndex)> = contigs_y
        .par_iter()
        .map(|(n, s)| {
            Ok((
                n.clone(),
                prepare_contig_index_opts(s, cli.respect_softmask)?,
            ))
        })
        .collect::<Result<Vec<_>>>()?;

    let cfg = AlignConfig {
        min_top_lcb_anchors: cli.min_top_lcb,
        skip_softmasked: cli.respect_softmask,
        use_ry_seeds: cli.use_ry_seeds,
        min_fallback_cov: cli.min_fallback_cov,
        hsp_threshold: cli.hsp_threshold,
        ..AlignConfig::default()
    };

    let mut total_pairs = 0usize;
    let mut total_bubbles = 0usize;
    let mut pruned_pairs = 0usize;
    let mut summary_out = String::new();
    let mut maf_writer = cli
        .maf_out
        .as_ref()
        .map(|p| File::create(p).map(BufWriter::new))
        .transpose()
        .context("creating MAF output")?;
    if let Some(w) = maf_writer.as_mut() {
        w.write_all(b"##maf version=1 scoring=horizon_mdl\n")?;
    }
    let mut gfa_writer = cli
        .gfa_out
        .as_ref()
        .map(|p| File::create(p).map(BufWriter::new))
        .transpose()
        .context("creating GFA output")?;
    if let Some(w) = gfa_writer.as_mut() {
        w.write_all(b"H\tVN:Z:1.1\n")?;
    }

    let allowed = (cli.top_k > 0)
        .then(|| select_top_k_pairs(&idx_x, &idx_y, cli.top_k, cli.min_shared_seeds));

    let mut candidate_pairs: Vec<(usize, usize)> = Vec::new();
    for (x_i, (_, ix)) in idx_x.iter().enumerate() {
        for (y_i, (_, iy)) in idx_y.iter().enumerate() {
            if let Some(mask) = allowed.as_ref() {
                if !mask[x_i][y_i] {
                    pruned_pairs += 1;
                    continue;
                }
            } else if shared_seed_count(ix, iy) < cli.min_shared_seeds {
                pruned_pairs += 1;
                continue;
            }
            candidate_pairs.push((x_i, y_i));
        }
    }

    let aligned_results: Vec<Result<_>> = candidate_pairs
        .into_par_iter()
        .map(|(x_i, y_i)| {
            let mut aln = align_contig_pair_runs(&idx_x[x_i].1, &idx_y[y_i].1, &cfg)?;
            if cli.one_to_one && !aln.aligned_runs.is_empty() {
                aln.aligned_runs = net_one_to_one_runs(
                    &aln.aligned_runs,
                    contigs_x[x_i].1.len(),
                    contigs_y[y_i].1.len(),
                );
                aln.total_pairs = aln.aligned_runs.iter().map(|r| r.len).sum();
            }
            Ok((x_i, y_i, aln))
        })
        .collect();

    let mut target_cov_by_contig: Vec<Vec<u64>> = contigs_x
        .iter()
        .map(|(_, s)| vec![0_u64; s.len().div_ceil(64)])
        .collect();
    let mut query_cov_by_contig: Vec<Vec<u64>> = contigs_y
        .iter()
        .map(|(_, s)| vec![0_u64; s.len().div_ceil(64)])
        .collect();

    for item in aligned_results {
        let (x_i, y_i, aln) = item?;
        let nx = idx_x[x_i].0.as_str();
        let ny = idx_y[y_i].0.as_str();
        if aln.aligned_runs.is_empty() {
            pruned_pairs += 1;
            continue;
        }
        total_pairs += aln.total_pairs;
        let x_bits = &mut target_cov_by_contig[x_i];
        let y_bits = &mut query_cov_by_contig[y_i];
        for r in &aln.aligned_runs {
            for px in r.x_start..(r.x_start + r.len) {
                if let Some(w) = x_bits.get_mut(px >> 6) {
                    *w |= 1_u64 << (px & 63);
                }
            }
            if r.is_rc {
                let y_lo = r.y_start.saturating_add(1).saturating_sub(r.len);
                for py in y_lo..=r.y_start {
                    if let Some(w) = y_bits.get_mut(py >> 6) {
                        *w |= 1_u64 << (py & 63);
                    }
                }
            } else {
                for py in r.y_start..(r.y_start + r.len) {
                    if let Some(w) = y_bits.get_mut(py >> 6) {
                        *w |= 1_u64 << (py & 63);
                    }
                }
            }
        }
        total_bubbles += aln.macro_bubbles + aln.internal_bubbles;
        let gfa_sum = format_gfa_summary_runs(nx, ny, &aln);
        if summary_out.is_empty() {
            summary_out.push_str(&gfa_sum);
        } else {
            summary_out.push_str(gfa_sum.strip_prefix("H\tVN:Z:1.1\n").unwrap_or(&gfa_sum));
        }

        if let Some(w) = maf_writer.as_mut() {
            let sx = &contigs_x[x_i].1;
            let sy = &contigs_y[y_i].1;
            let text = format_maf_blocks_runs(nx, sx.len(), ny, sy.len(), sx, sy, &aln);
            let body = text
                .strip_prefix("##maf version=1 scoring=horizon_mdl\n")
                .unwrap_or(&text);
            w.write_all(body.as_bytes())?;
        }
        if let Some(w) = gfa_writer.as_mut() {
            let text = format_gfa_blocks_runs(nx, ny, &aln);
            let body = text.strip_prefix("H\tVN:Z:1.1\n").unwrap_or(&text);
            w.write_all(body.as_bytes())?;
        }
        if cli.resolve_dups {
            let dups = resolve_alignment_duplications_runs(&aln, 1000);
            let tag_line = format_gfa_paralog_tags(nx, ny, &dups);
            if !tag_line.is_empty() {
                summary_out.push_str(&tag_line);
                if let Some(w) = gfa_writer.as_mut() {
                    w.write_all(tag_line.as_bytes())?;
                }
            }
        }
    }

    if let Some(w) = maf_writer.as_mut() {
        w.flush()?;
    }
    if let Some(w) = gfa_writer.as_mut() {
        w.flush()?;
    }

    let total_unique_x: usize = target_cov_by_contig
        .iter()
        .map(|bits| bits.iter().map(|w| w.count_ones() as usize).sum::<usize>())
        .sum();
    let total_unique_y: usize = query_cov_by_contig
        .iter()
        .map(|bits| bits.iter().map(|w| w.count_ones() as usize).sum::<usize>())
        .sum();

    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "contigs_x={} contigs_y={} pruned_pairs={} aligned_pairs={} unique_x_bp={} unique_y_bp={} total_bubbles={} elapsed_s={:.4}",
        contigs_x.len(),
        contigs_y.len(),
        pruned_pairs,
        total_pairs,
        total_unique_x,
        total_unique_y,
        total_bubbles,
        elapsed
    );
    if !summary_out.is_empty() {
        print!("{summary_out}");
    }
    Ok(())
}

#[allow(clippy::too_many_lines)]
fn estimate_dag_hmm_params(nodes: &[GraphNode], seq_y: &[u8]) -> HmmParams {
    use horizon_aligner::complement_base;
    use horizon_aligner::seeding::{MASK_W10, canonical_spaced_seeds_with_opts};
    let Ok(q_seeds) = canonical_spaced_seeds_with_opts(seq_y, MASK_W10, 4, false, false) else {
        return HmmParams::default();
    };
    let span = MASK_W10.len();
    let mut matches = 0usize;
    let mut total = 0usize;
    let mut indel_events = 0usize;
    let mut total_indel_len = 0usize;
    let mut total_span = 0usize;
    for nd in nodes {
        let Ok(n_seeds) = canonical_spaced_seeds_with_opts(&nd.seq, MASK_W10, 4, false, false)
        else {
            continue;
        };
        let mut keys: Vec<u64> = n_seeds.keys().copied().collect();
        keys.sort_unstable();
        let mut node_hits: Vec<(usize, usize)> = Vec::new();
        for k in keys {
            if let (Some(n_pos), Some(q_pos)) = (n_seeds.get(&k), q_seeds.get(&k)) {
                for &np in n_pos {
                    for &qp in q_pos {
                        if np + span <= nd.seq.len() && qp + span <= seq_y.len() {
                            node_hits.push((np, qp));
                            let mut fwd_m = 0usize;
                            let mut rc_m = 0usize;
                            let mut valid_f = 0usize;
                            let mut valid_r = 0usize;
                            for off in 0..span {
                                let bx = nd.seq[np + off].to_ascii_uppercase();
                                let by_f = seq_y[qp + off].to_ascii_uppercase();
                                let by_r = complement_base(seq_y[qp + span - 1 - off])
                                    .to_ascii_uppercase();
                                if matches!(bx, b'A' | b'C' | b'G' | b'T') {
                                    if matches!(by_f, b'A' | b'C' | b'G' | b'T') {
                                        valid_f += 1;
                                        if bx == by_f {
                                            fwd_m += 1;
                                        }
                                    }
                                    if matches!(by_r, b'A' | b'C' | b'G' | b'T') {
                                        valid_r += 1;
                                        if bx == by_r {
                                            rc_m += 1;
                                        }
                                    }
                                }
                            }
                            if fwd_m * valid_r >= rc_m * valid_f {
                                total += valid_f;
                                matches += fwd_m;
                            } else {
                                total += valid_r;
                                matches += rc_m;
                            }
                        }
                        if total >= 2000 {
                            break;
                        }
                    }
                    if total >= 2000 {
                        break;
                    }
                }
            }
            if total >= 2000 {
                break;
            }
        }
        node_hits.sort_unstable();
        for w in node_hits.windows(2).take(64) {
            let dx = w[1].0.saturating_sub(w[0].0);
            let dy = w[1].1.abs_diff(w[0].1);
            if dx > 0 && dy > 0 && dx <= 500 && dy <= 500 {
                total_span += dx.max(dy);
                if dx != dy {
                    indel_events += 1;
                    total_indel_len += dx.abs_diff(dy);
                }
            }
        }
        if total >= 2000 {
            break;
        }
    }
    if total < 60 {
        return HmmParams::default();
    }
    let defaults = HmmParams::default();
    #[allow(clippy::cast_precision_loss)]
    let identity_prior = ((matches as f64 + 195.0) / (total as f64 + 300.0)).clamp(0.60, 0.72);
    #[allow(clippy::cast_precision_loss)]
    let gap_open_bits = if total_span >= 200 {
        ((indel_events as f64 + 1.0) / (total_span as f64 + 22.0))
            .log2()
            .clamp(-5.5, -4.0)
    } else {
        defaults.gap_open_bits
    };
    #[allow(clippy::cast_precision_loss)]
    let gap_extend_bits = if indel_events >= 3 && total_indel_len > 0 {
        let mean_len = (total_indel_len as f64 + 2.5) / (indel_events as f64 + 1.0);
        ((mean_len - 1.0) / mean_len)
            .max(0.01)
            .log2()
            .clamp(-1.8, -0.5)
    } else {
        defaults.gap_extend_bits
    };
    HmmParams {
        identity_prior,
        gap_open_bits,
        gap_extend_bits,
    }
}

type OrientedSegs = Vec<(String, Vec<u8>)>;
type DirectedEdges = Vec<(usize, usize)>;

#[allow(clippy::too_many_lines)]
fn load_gfa_dag(path: &std::path::Path) -> Result<Vec<GraphNode>> {
    let mut segs = Vec::<(String, Vec<u8>)>::new();
    let mut name_to_raw = HashMap::<String, usize>::new();
    let mut raw_edges = Vec::<(String, u8, String, u8)>::new();

    stream_gfa_records(path, |rec| {
        match rec {
            GfaRecord::Segment(s) if s.sequence != b"*" => {
                let idx = segs.len();
                name_to_raw.insert(s.name.clone(), idx);
                segs.push((s.name, s.sequence));
            }
            GfaRecord::Link(l) => {
                raw_edges.push((l.from, l.from_orient, l.to, l.to_orient));
            }
            _ => {}
        }
        Ok(std::ops::ControlFlow::Continue(()))
    })?;

    let has_rev = raw_edges
        .iter()
        .any(|(_, fo, _, to)| *fo == b'-' || *to == b'-');

    let (oriented_segs, directed_edges): (OrientedSegs, DirectedEdges) = if has_rev {
        let mut o_segs = Vec::with_capacity(segs.len() * 2);
        for (name, seq) in &segs {
            o_segs.push((format!("{name}+"), seq.clone()));
            o_segs.push((format!("{name}-"), revcomp(seq)));
        }
        let orient_idx = |raw_i: usize, orient: u8| -> usize {
            if orient == b'-' {
                raw_i * 2 + 1
            } else {
                raw_i * 2
            }
        };
        let flip = |o: u8| -> u8 { if o == b'-' { b'+' } else { b'-' } };
        let mut d_edges = Vec::with_capacity(raw_edges.len() * 2);
        for (u_name, uo, v_name, vo) in &raw_edges {
            if let (Some(&u_raw), Some(&v_raw)) = (name_to_raw.get(u_name), name_to_raw.get(v_name))
            {
                let u = orient_idx(u_raw, *uo);
                let v = orient_idx(v_raw, *vo);
                if u != v {
                    d_edges.push((u, v));
                }
                let rv = orient_idx(v_raw, flip(*vo));
                let ru = orient_idx(u_raw, flip(*uo));
                if rv != ru {
                    d_edges.push((rv, ru));
                }
            }
        }
        (o_segs, d_edges)
    } else {
        let mut d_edges = Vec::with_capacity(raw_edges.len());
        for (u_name, _, v_name, _) in &raw_edges {
            if let (Some(&u), Some(&v)) = (name_to_raw.get(u_name), name_to_raw.get(v_name))
                && u != v
            {
                d_edges.push((u, v));
            }
        }
        (segs, d_edges)
    };

    let n = oriented_segs.len();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut indeg: Vec<usize> = vec![0; n];

    for (u, v) in directed_edges {
        adj[u].push(v);
    }
    for nbrs in &mut adj {
        nbrs.sort_unstable();
        nbrs.dedup();
        for &v in nbrs.iter() {
            indeg[v] += 1;
        }
    }

    let mut q: VecDeque<usize> = (0..n).filter(|&i| indeg[i] == 0).collect();
    let mut order = Vec::with_capacity(n);
    while let Some(u) = q.pop_front() {
        order.push(u);
        for &v in &adj[u] {
            indeg[v] = indeg[v].saturating_sub(1);
            if indeg[v] == 0 {
                q.push_back(v);
            }
        }
    }
    // Append any remaining nodes (if cycles existed) in input order
    if order.len() < n {
        eprintln!(
            "warning: cyclic GFA detected ({} of {} oriented nodes in cycles); back-edges will be dropped to form a DAG",
            n - order.len(),
            n
        );
        let seen: std::collections::HashSet<usize> = order.iter().copied().collect();
        for i in 0..n {
            if !seen.contains(&i) {
                order.push(i);
            }
        }
    }

    let mut raw_to_topo = vec![0usize; n];
    for (topo_i, &raw_i) in order.iter().enumerate() {
        raw_to_topo[raw_i] = topo_i;
    }

    let mut nodes: Vec<GraphNode> = order
        .iter()
        .map(|&raw_i| GraphNode {
            name: oriented_segs[raw_i].0.clone(),
            seq: oriented_segs[raw_i].1.clone(),
            preds: Vec::new(),
        })
        .collect();

    for (u_raw, vs) in adj.iter().enumerate() {
        let u_topo = raw_to_topo[u_raw];
        for &v_raw in vs {
            let v_topo = raw_to_topo[v_raw];
            if u_topo < v_topo {
                nodes[v_topo].preds.push(u_topo);
            }
        }
    }
    Ok(nodes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn load_gfa_dag_handles_bidirected_minus_links() -> Result<()> {
        let dir = std::env::temp_dir().join(format!("horizon_wga_dag_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir)?;
        let gfa = dir.join("inv.gfa");
        // Seg a -> Seg b(-) -> Seg c(+)
        std::fs::write(
            &gfa,
            "H\tVN:Z:1.0\n\
             S\ta\tACTGACTGACTGACTGACTG\n\
             S\tb\tCCCCCCCCCCCCCCCCAAAA\n\
             S\tc\tGATTACAGATTACAGATTAC\n\
             L\ta\t+\tb\t-\t0M\n\
             L\tb\t-\tc\t+\t0M\n",
        )?;
        let nodes = load_gfa_dag(&gfa)?;
        assert_eq!(nodes.len(), 6);
        let b_rc = revcomp(b"CCCCCCCCCCCCCCCCAAAA");
        let mut query = b"ACTGACTGACTGACTGACTG".to_vec();
        query.extend_from_slice(&b_rc);
        query.extend_from_slice(b"GATTACAGATTACAGATTAC");
        let aln =
            horizon_aligner::align_seq_to_dag_seeded(&nodes, &query, &HmmParams::default(), 10.0);
        let path_names: Vec<&str> = aln
            .path_nodes
            .iter()
            .map(|&i| nodes[i].name.as_str())
            .collect();
        assert_eq!(path_names, vec!["a+", "b-", "c+"]);
        assert_eq!(aln.kept_pairs.len(), query.len());
        let _ = std::fs::remove_dir_all(&dir);
        Ok(())
    }

    #[test]
    fn estimate_dag_hmm_params_adapts_identity_and_gaps() {
        let short_nodes = vec![GraphNode {
            name: "s+".to_string(),
            seq: b"ACGTACGTACGT".to_vec(),
            preds: vec![],
        }];
        let def = HmmParams::default();
        let fallback = estimate_dag_hmm_params(&short_nodes, b"ACGTACGTACGT");
        assert!((fallback.identity_prior - def.identity_prior).abs() < 1e-9);

        let bases = *b"ACGT";
        let mut ref_seq = Vec::with_capacity(400);
        for i in 0..400 {
            ref_seq.push(bases[(i * 7 + (i / 11) * 3) % 4]);
        }
        // Introduce 4 small deletions in query so dx != dy across seed pairs
        let mut query = Vec::new();
        for (idx, &b) in ref_seq.iter().enumerate() {
            if matches!(idx, 70 | 140 | 210 | 280 | 281 | 282) {
                continue;
            }
            query.push(b);
        }
        let nodes = vec![GraphNode {
            name: "n0+".to_string(),
            seq: ref_seq,
            preds: vec![],
        }];
        let adapted = estimate_dag_hmm_params(&nodes, &query);
        assert!(adapted.identity_prior >= 0.65);
        assert!(adapted.gap_open_bits >= -5.5 && adapted.gap_open_bits <= -4.0);
        assert!(adapted.gap_extend_bits >= -1.8 && adapted.gap_extend_bits <= -0.5);
    }
}
