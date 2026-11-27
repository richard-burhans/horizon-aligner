//! Multi-contig whole-genome aligner CLI using KP-14 + parallel MDL/BOCPD.

mod graph;

use anyhow::{Context, Result, bail};
use graph::load_gfa_dag;
use horizon_aligner::{
    AlignConfig, ContigIndex, DagSeedIndex, TantanConfig, adaptive_top_k_for_fragmentation,
    align_contig_pair_runs, align_seq_to_dag_seeded_indexed, estimate_dag_hmm_params,
    format_gfa_blocks_runs, format_gfa_paralog_tags, format_gfa_summary_runs, init_rayon_threads,
    net_one_to_one_runs_scored, parse_threads_arg, prepare_contig_index_owned,
    prune_resolved_dup_clusters, read_fasta_contigs, resolve_alignment_duplications_mdl,
    select_top_k_pairs, shared_seed_count, tantan_mask, train_asymmetric_pair_params,
};
use rayon::prelude::*;
use std::fmt::Write as _;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

const WGA_USAGE: &str = "usage: horizon_wga <target.fa> <query.fa> [max_contigs] [--maf out.maf] [--gfa out.gfa] [--paf out.paf] [--chain out.chain] [--net out.net] [--threads N | -t N] [--graph target.gfa] [--hal-info file.hal] [--min-shared-seeds N] [--min-top-lcb N] [--top-k K | --rbest] [--resolve-dups] [--one-to-one] [--respect-softmask] [--tantan] [--train-hmm] [--ry-seeds | --no-ry-seeds] [--no-aa5] [--min-fallback-cov F] [--hsp-threshold N]\n       (max_contigs defaults to 1; pass 0 to align all contigs)";

#[allow(clippy::struct_excessive_bools)]
struct CliArgs {
    fa_x: PathBuf,
    fa_y: PathBuf,
    max_contigs: usize,
    maf_out: Option<PathBuf>,
    gfa_out: Option<PathBuf>,
    paf_out: Option<PathBuf>,
    chain_out: Option<PathBuf>,
    net_out: Option<PathBuf>,
    threads: Option<usize>,
    resolve_dups: bool,
    one_to_one: bool,
    min_shared_seeds: usize,
    min_top_lcb: usize,
    top_k: usize,
    respect_softmask: bool,
    use_tantan: bool,
    train_hmm: bool,
    use_ry_seeds: bool,
    use_aa5: bool,
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
    let mut paf_out = None;
    let mut chain_out = None;
    let mut net_out = None;
    let mut threads = None;
    let mut resolve_dups = false;
    let mut one_to_one = false;
    let mut min_shared_seeds = 8usize;
    let mut min_top_lcb = 0usize;
    let mut top_k = 0usize;
    let mut respect_softmask = false;
    let mut use_tantan = false;
    let mut train_hmm = false;
    let mut use_ry_seeds = false;
    let mut use_aa5 = true;
    let mut hal_info = None;
    let mut graph_path = None;
    let mut min_fallback_cov = None;
    let mut hsp_threshold = None;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                println!("{WGA_USAGE}");
                return Ok(None);
            }
            "-V" | "--version" => {
                println!("horizon_wga {}", env!("CARGO_PKG_VERSION"));
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
            "--paf" => {
                paf_out = Some(PathBuf::from(
                    args.get(i + 1).context("missing value for --paf")?,
                ));
                i += 2;
            }
            "--chain" => {
                chain_out = Some(PathBuf::from(
                    args.get(i + 1).context("missing value for --chain")?,
                ));
                i += 2;
            }
            "--net" => {
                net_out = Some(PathBuf::from(
                    args.get(i + 1).context("missing value for --net")?,
                ));
                i += 2;
            }
            "-t" | "--threads" => {
                let raw = args.get(i + 1).context("missing value for --threads")?;
                threads = Some(parse_threads_arg(raw)?);
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
            "--tantan" => {
                use_tantan = true;
                i += 1;
            }
            "--train-hmm" => {
                train_hmm = true;
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
            "--no-aa5" => {
                use_aa5 = false;
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
    let has_aln_out = maf_out.is_some()
        || gfa_out.is_some()
        || paf_out.is_some()
        || chain_out.is_some()
        || net_out.is_some();
    let (fa_x, fa_y, max_contigs) = if hal_info.is_some() {
        if !pos.is_empty() {
            bail!("--hal-info does not accept positional arguments");
        }
        if has_aln_out {
            bail!(
                "--maf, --gfa, --paf, --chain, and --net outputs are not supported in --hal-info mode"
            );
        }
        (PathBuf::new(), PathBuf::new(), 1)
    } else if graph_path.is_some() {
        if pos.is_empty() || pos.len() > 2 {
            bail!(
                "--graph mode accepts <query.fa> [max_contigs], got {} positional args",
                pos.len()
            );
        }
        if has_aln_out {
            bail!(
                "--maf, --gfa, --paf, --chain, and --net outputs are not supported in --graph mode"
            );
        }
        let mc = match pos.get(1) {
            Some(s) => s.parse().context("parsing max_contigs")?,
            None => 1,
        };
        (PathBuf::new(), PathBuf::from(&pos[0]), mc)
    } else {
        if pos.len() < 2 || pos.len() > 3 {
            bail!("{WGA_USAGE}");
        }
        let mc = match pos.get(2) {
            Some(s) => s.parse().context("parsing max_contigs")?,
            None => 1,
        };
        (PathBuf::from(&pos[0]), PathBuf::from(&pos[1]), mc)
    };

    Ok(Some(CliArgs {
        fa_x,
        fa_y,
        max_contigs,
        maf_out,
        gfa_out,
        paf_out,
        chain_out,
        net_out,
        threads,
        resolve_dups,
        one_to_one,
        min_shared_seeds,
        min_top_lcb,
        top_k,
        respect_softmask,
        use_tantan,
        train_hmm,
        use_ry_seeds,
        use_aa5,
        hal_info,
        graph_path,
        min_fallback_cov,
        hsp_threshold,
    }))
}

fn main() -> Result<()> {
    let Some(cli) = parse_args()? else {
        return Ok(());
    };
    init_rayon_threads(cli.threads)?;
    if let Some(hal_path) = &cli.hal_info {
        return run_hal_info(hal_path);
    }
    if let Some(graph_p) = &cli.graph_path {
        return run_graph_mode(&cli, graph_p);
    }
    run_pairwise_wga(&cli)
}

fn run_hal_info(hal_path: &std::path::Path) -> Result<()> {
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
    Ok(())
}

fn run_graph_mode(cli: &CliArgs, graph_p: &std::path::Path) -> Result<()> {
    anyhow::ensure!(
        !cli.use_tantan && !cli.train_hmm,
        "--tantan and --train-hmm are only supported in pairwise alignment mode (not --graph)"
    );
    let t0 = Instant::now();
    let dag_nodes = load_gfa_dag(graph_p)?;
    let dag_index = DagSeedIndex::build(&dag_nodes);
    let queries = read_fasta_contigs(&cli.fa_y, cli.max_contigs)?;
    let results: Vec<(&str, usize, _)> = queries
        .par_iter()
        .map(|(qname, qseq)| {
            let dag_params = estimate_dag_hmm_params(&dag_nodes, &dag_index, qseq);
            let aln =
                align_seq_to_dag_seeded_indexed(&dag_nodes, &dag_index, qseq, &dag_params, 10.0);
            (qname.as_str(), qseq.len(), aln)
        })
        .collect();
    let mut total_pairs = 0usize;
    let mut total_bubbles = 0usize;
    for (qname, qlen, aln) in results {
        total_pairs += aln.kept_pairs.len();
        total_bubbles += aln.bubbles;
        let path_str = aln.path_nodes.iter().fold(String::new(), |mut acc, &idx| {
            let _ = write!(acc, ">{}", dag_nodes[idx].name);
            acc
        });
        println!(
            "{}\t{}\t{}\tkept_bp={}\tscore_bits={:.1}\tbubbles={}",
            qname,
            qlen,
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
    Ok(())
}

fn apply_optional_tantan(mut seq: Vec<u8>, enable: bool) -> Vec<u8> {
    if enable {
        let mask = tantan_mask(&seq, &TantanConfig::default());
        for (b, masked) in seq.iter_mut().zip(mask) {
            if masked {
                b.make_ascii_lowercase();
            }
        }
    }
    seq
}

#[allow(clippy::too_many_lines)]
fn run_pairwise_wga(cli: &CliArgs) -> Result<()> {
    let t0 = Instant::now();
    let skip_soft = cli.respect_softmask || cli.use_tantan;
    let idx_x: Vec<(String, ContigIndex)> = read_fasta_contigs(&cli.fa_x, cli.max_contigs)?
        .into_par_iter()
        .map(|(n, s)| {
            let s = apply_optional_tantan(s, cli.use_tantan);
            Ok((n, prepare_contig_index_owned(s, skip_soft)?))
        })
        .collect::<Result<Vec<_>>>()?;
    let idx_y: Vec<(String, ContigIndex)> = read_fasta_contigs(&cli.fa_y, cli.max_contigs)?
        .into_par_iter()
        .map(|(n, s)| {
            let s = apply_optional_tantan(s, cli.use_tantan);
            Ok((n, prepare_contig_index_owned(s, skip_soft)?))
        })
        .collect::<Result<Vec<_>>>()?;

    let cfg = AlignConfig {
        min_top_lcb_anchors: cli.min_top_lcb,
        skip_softmasked: skip_soft,
        use_ry_seeds: cli.use_ry_seeds,
        use_aa5: cli.use_aa5,
        min_fallback_cov: cli.min_fallback_cov,
        hsp_threshold: cli.hsp_threshold,
        ..AlignConfig::default()
    };

    let stdout_outputs = [
        &cli.maf_out,
        &cli.gfa_out,
        &cli.paf_out,
        &cli.chain_out,
        &cli.net_out,
    ]
    .iter()
    .filter(|opt| opt.as_deref() == Some(std::path::Path::new("-")))
    .count();
    anyhow::ensure!(
        stdout_outputs <= 1,
        "at most one of --maf, --gfa, --paf, --chain, or --net may write to stdout ('-')"
    );
    let writes_stdout = stdout_outputs > 0;

    let mut total_pairs = 0usize;
    let mut total_bubbles = 0usize;
    let mut pruned_pairs = 0usize;
    let mut summary_out = String::new();
    let mut maf_writer = cli
        .maf_out
        .as_deref()
        .map(|p| horizon_aligner::io_util::open_text_writer(Some(p.to_string_lossy().as_ref())))
        .transpose()
        .context("creating MAF output")?;
    if let Some(w) = maf_writer.as_mut() {
        w.write_all(b"##maf version=1 scoring=horizon_mdl\n")?;
    }
    let mut gfa_writer = cli
        .gfa_out
        .as_deref()
        .map(|p| horizon_aligner::io_util::open_text_writer(Some(p.to_string_lossy().as_ref())))
        .transpose()
        .context("creating GFA output")?;
    if let Some(w) = gfa_writer.as_mut() {
        w.write_all(b"H\tVN:Z:1.1\n")?;
    }
    let mut paf_writer = cli
        .paf_out
        .as_deref()
        .map(|p| horizon_aligner::io_util::open_text_writer(Some(p.to_string_lossy().as_ref())))
        .transpose()
        .context("creating PAF output")?;

    let eff_top_k = adaptive_top_k_for_fragmentation(cli.top_k, &idx_x, &idx_y);
    let allowed = (eff_top_k > 0)
        .then(|| select_top_k_pairs(&idx_x, &idx_y, eff_top_k, cli.min_shared_seeds));

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

    let mut aligned_results: Vec<(usize, usize, horizon_aligner::RunAlignment)> = candidate_pairs
        .into_par_iter()
        .map(|(x_i, y_i)| {
            let mut aln = align_contig_pair_runs(&idx_x[x_i].1, &idx_y[y_i].1, &cfg)?;
            if cli.one_to_one && !aln.aligned_runs.is_empty() {
                aln.aligned_runs = net_one_to_one_runs_scored(
                    &aln.aligned_runs,
                    idx_x[x_i].1.seq(),
                    idx_y[y_i].1.seq(),
                );
                aln.total_pairs = aln.aligned_runs.iter().map(|r| r.len).sum();
            }
            Ok((x_i, y_i, aln))
        })
        .collect::<Result<Vec<_>>>()?;

    // When aligning multi-contig assemblies with --one-to-one, enforce species-wide global
    // 1-to-1 orthology across ALL contigs in O(R log R) time/memory via `src/netting.rs`
    // (eliminating cross-contig paralog double-mapping, e.g., yeast HML/MAT/HMR cassettes).
    let (total_unique_x, total_unique_y) = if cli.one_to_one && (idx_x.len() > 1 || idx_y.len() > 1)
    {
        let seqs_x: Vec<&[u8]> = idx_x.iter().map(|(_, ix)| ix.seq()).collect();
        let seqs_y: Vec<&[u8]> = idx_y.iter().map(|(_, iy)| iy.seq()).collect();
        let mut pair_runs: Vec<(usize, usize, Vec<horizon_aligner::AlignedRun>)> = aligned_results
            .iter_mut()
            .map(|(xi, yi, aln)| (*xi, *yi, std::mem::take(&mut aln.aligned_runs)))
            .collect();
        let unique_cov =
            horizon_aligner::net_global_one_to_one_runs_scored(&mut pair_runs, &seqs_x, &seqs_y);
        for ((_, _, aln), (_, _, kept_runs)) in aligned_results.iter_mut().zip(pair_runs) {
            aln.total_pairs = kept_runs.iter().map(|r| r.len).sum();
            aln.aligned_runs = kept_runs;
        }
        unique_cov
    } else {
        let refs: Vec<(usize, usize, &[horizon_aligner::AlignedRun])> = aligned_results
            .iter()
            .map(|(xi, yi, aln)| (*xi, *yi, aln.aligned_runs.as_slice()))
            .collect();
        horizon_aligner::count_unique_covered_bp(&refs, idx_x.len(), idx_y.len())
    };

    if cli.train_hmm {
        const MAX_TRAIN_COLS: usize = 250_000;
        let mut train_x = Vec::new();
        let mut train_y = Vec::new();
        for &(x_i, y_i, ref aln) in &aligned_results {
            if train_x.len() >= MAX_TRAIN_COLS {
                break;
            }
            horizon_aligner::mdl::extract_aligned_train_columns(
                &aln.aligned_runs,
                idx_x[x_i].1.seq(),
                idx_y[y_i].1.seq(),
                MAX_TRAIN_COLS,
                &mut train_x,
                &mut train_y,
            );
        }
        let trained = train_asymmetric_pair_params(&[(&train_x, &train_y)], 5);
        eprintln!(
            "[horizon_wga --train-hmm] cols={} pi_x=[{:.3},{:.3},{:.3},{:.3}] pi_y=[{:.3},{:.3},{:.3},{:.3}] sub_diag=[{:.2},{:.2},{:.2},{:.2}] del=({:.2},{:.2}) ins=({:.2},{:.2})",
            train_x.len(),
            trained.pi_x[0],
            trained.pi_x[1],
            trained.pi_x[2],
            trained.pi_x[3],
            trained.pi_y[0],
            trained.pi_y[1],
            trained.pi_y[2],
            trained.pi_y[3],
            trained.sub_bits[0][0],
            trained.sub_bits[1][1],
            trained.sub_bits[2][2],
            trained.sub_bits[3][3],
            trained.del_open_bits,
            trained.del_ext_bits,
            trained.ins_open_bits,
            trained.ins_ext_bits,
        );
    }

    if cli.chain_out.is_some() || cli.net_out.is_some() {
        let nx: Vec<&str> = idx_x.iter().map(|(n, _)| n.as_str()).collect();
        let sx: Vec<&[u8]> = idx_x.iter().map(|(_, ix)| ix.seq()).collect();
        let ny: Vec<&str> = idx_y.iter().map(|(n, _)| n.as_str()).collect();
        let sy: Vec<&[u8]> = idx_y.iter().map(|(_, iy)| iy.seq()).collect();
        let refs: Vec<(usize, usize, &[horizon_aligner::AlignedRun])> = aligned_results
            .iter()
            .map(|(xi, yi, aln)| (*xi, *yi, aln.aligned_runs.as_slice()))
            .collect();
        horizon_aligner::write_ucsc_chain_and_net_files(
            &refs,
            &nx,
            &sx,
            &ny,
            &sy,
            cli.chain_out.as_deref(),
            cli.net_out.as_deref(),
        )?;
    }

    // When --one-to-one is active, collect the global claimed target/query intervals across all
    // kept runs so MAF/PAF colinear block grouping never bridges a <=40 bp gap claimed by another
    // contig pair (or diagonal), preserving strict 0-bp target/query overlap across all outputs.
    let (claimed_x, claimed_y): (
        Option<Vec<horizon_aligner::IntervalSet>>,
        Option<Vec<horizon_aligner::IntervalSet>>,
    ) = if cli.one_to_one {
        let mut cx: Vec<horizon_aligner::IntervalSet> = (0..idx_x.len())
            .map(|_| horizon_aligner::IntervalSet::new())
            .collect();
        let mut cy: Vec<horizon_aligner::IntervalSet> = (0..idx_y.len())
            .map(|_| horizon_aligner::IntervalSet::new())
            .collect();
        let mut scratch = Vec::new();
        for &(x_i, y_i, ref aln) in &aligned_results {
            for r in &aln.aligned_runs {
                if r.len == 0 {
                    continue;
                }
                cx[x_i].insert(r.x_start, r.x_start.saturating_add(r.len), &mut scratch);
                if r.is_rc {
                    if let Some(y_lo) = r
                        .y_start
                        .checked_add(1)
                        .and_then(|y1| y1.checked_sub(r.len))
                    {
                        cy[y_i].insert(y_lo, r.y_start + 1, &mut scratch);
                    }
                } else {
                    cy[y_i].insert(r.y_start, r.y_start.saturating_add(r.len), &mut scratch);
                }
            }
        }
        (Some(cx), Some(cy))
    } else {
        (None, None)
    };

    for (x_i, y_i, aln) in aligned_results {
        let nx = idx_x[x_i].0.as_str();
        let ny = idx_y[y_i].0.as_str();
        if aln.aligned_runs.is_empty() {
            pruned_pairs += 1;
            continue;
        }
        total_pairs += aln.total_pairs;
        total_bubbles += aln.macro_bubbles + aln.internal_bubbles;
        let gfa_sum = format_gfa_summary_runs(nx, ny, &aln);
        if summary_out.is_empty() {
            summary_out.push_str(&gfa_sum);
        } else {
            summary_out.push_str(gfa_sum.strip_prefix("H\tVN:Z:1.1\n").unwrap_or(&gfa_sum));
        }

        let sx = idx_x[x_i].1.seq();
        let sy = idx_y[y_i].1.seq();
        let resolved_dups = if cli.resolve_dups {
            resolve_alignment_duplications_mdl(&aln, sx, sy)
        } else {
            horizon_aligner::ResolvedDuplications::default()
        };
        let fx = claimed_x.as_ref().map(|v| &v[x_i]);
        let fy = claimed_y.as_ref().map(|v| &v[y_i]);

        if let Some(w) = maf_writer.as_mut() {
            let pruned_maf_aln;
            let maf_aln = if cli.resolve_dups && !resolved_dups.dup_cluster_ids.is_empty() {
                pruned_maf_aln = prune_resolved_dup_clusters(&aln, &resolved_dups, sx, sy);
                &pruned_maf_aln
            } else {
                &aln
            };
            let text = horizon_aligner::format_maf_blocks_runs_avoiding_gaps(
                nx,
                sx.len(),
                ny,
                sy.len(),
                sx,
                sy,
                maf_aln,
                fx,
                fy,
            );
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
        if let Some(w) = paf_writer.as_mut() {
            let text = horizon_aligner::paf::format_paf_records_runs_avoiding_gaps(
                nx,
                sx.len(),
                ny,
                sy.len(),
                sx,
                sy,
                &aln,
                fx,
                fy,
            );
            w.write_all(text.as_bytes())?;
        }
        if cli.resolve_dups {
            let tag_line = format_gfa_paralog_tags(nx, ny, &resolved_dups.dup_labels);
            if !tag_line.is_empty() {
                summary_out.push_str(&tag_line);
                if let Some(w) = gfa_writer.as_mut() {
                    w.write_all(tag_line.as_bytes())?;
                }
            }
        }
    }

    if let Some(w) = maf_writer {
        w.finish()?;
    }
    if let Some(w) = gfa_writer {
        w.finish()?;
    }
    if let Some(w) = paf_writer {
        w.finish()?;
    }

    let elapsed = t0.elapsed().as_secs_f64();
    if writes_stdout {
        eprintln!(
            "contigs_x={} contigs_y={} pruned_pairs={} aligned_pairs={} unique_x_bp={} unique_y_bp={} total_bubbles={} elapsed_s={:.4}",
            idx_x.len(),
            idx_y.len(),
            pruned_pairs,
            total_pairs,
            total_unique_x,
            total_unique_y,
            total_bubbles,
            elapsed
        );
    } else {
        println!(
            "contigs_x={} contigs_y={} pruned_pairs={} aligned_pairs={} unique_x_bp={} unique_y_bp={} total_bubbles={} elapsed_s={:.4}",
            idx_x.len(),
            idx_y.len(),
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
    }
    Ok(())
}
