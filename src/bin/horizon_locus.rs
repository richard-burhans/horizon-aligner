//! Regional `impg`-style multi-FASTA locus aligner, distance-matrix builder, and MAF/GFA emitter.
//!
//! Designed for VGP `O3` workflows (25 kb ILS windows and BUSCO gene intervals): takes a
//! multi-FASTA of homologous sequences, runs pairwise SDUST + BOCPD/MDL + `Aa5Index` alignment
//! across all pairs, and outputs `<out_prefix>.similarity.tsv`, `<out_prefix>.maf`, and
//! `<out_prefix>.gfa` (with Horizon Bubbles and `# EG:Z:` paralog annotations).

use anyhow::{Context, Result, bail};
use horizon_aligner::formats::{format_maf_blocks_runs, read_fasta_contigs};
use horizon_aligner::mdl::{AlignConfig, RunAlignment, align_contig_pair_runs};
use horizon_aligner::merge::merge_pairwise_mafs;
use horizon_aligner::seeding::{ContigIndex, complement_base, prepare_contig_index_opts};
use rayon::prelude::*;
use std::env;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;
use std::time::Instant;

struct LocusCli {
    fasta_path: String,
    out_prefix: String,
    focal_ref: Option<String>,
    threads: Option<usize>,
    respect_softmask: bool,
}

struct PairMetrics {
    i: usize,
    j: usize,
    aligned_bp: usize,
    matches: usize,
    jaccard_cov: f64,
    p_dist: f64,
    jc69_dist: f64,
    macro_bubbles: usize,
    internal_bubbles: usize,
}

fn count_run_matches(seq_x: &[u8], seq_y: &[u8], aln: &RunAlignment) -> usize {
    let mut matches = 0usize;
    for r in &aln.aligned_runs {
        for off in 0..r.len {
            let x = r.x_start + off;
            let Some(y) = (if r.is_rc {
                r.y_start.checked_sub(off)
            } else {
                Some(r.y_start + off)
            }) else {
                continue;
            };
            if x < seq_x.len() && y < seq_y.len() {
                let yb = if r.is_rc {
                    complement_base(seq_y[y])
                } else {
                    seq_y[y]
                };
                if seq_x[x].eq_ignore_ascii_case(&yb) {
                    matches += 1;
                }
            }
        }
    }
    matches
}

fn jc69_distance(p_dist: f64) -> f64 {
    if p_dist <= 0.0 {
        return 0.0;
    }
    let capped = p_dist.clamp(0.0, 0.74);
    let inner = (1.0 - (4.0 / 3.0) * capped).max(1e-6);
    (-0.75 * inner.ln()).max(0.0)
}

fn unique_covered_bp(intervals: &mut [(usize, usize)]) -> usize {
    if intervals.is_empty() {
        return 0;
    }
    intervals.sort_unstable();
    let mut total = 0usize;
    let (mut cur_s, mut cur_e) = intervals[0];
    for &(s, e) in &intervals[1..] {
        if s <= cur_e {
            cur_e = cur_e.max(e);
        } else {
            total += cur_e.saturating_sub(cur_s);
            cur_s = s;
            cur_e = e;
        }
    }
    total + cur_e.saturating_sub(cur_s)
}

fn compute_pair_metrics(
    i: usize,
    j: usize,
    seq_x: &[u8],
    seq_y: &[u8],
    aln: &RunAlignment,
) -> PairMetrics {
    let aligned_bp = aln.total_pairs;
    let matches = count_run_matches(seq_x, seq_y, aln);
    let mut x_ivs: Vec<(usize, usize)> = Vec::with_capacity(aln.aligned_runs.len());
    let mut y_ivs: Vec<(usize, usize)> = Vec::with_capacity(aln.aligned_runs.len());
    for r in &aln.aligned_runs {
        if r.len == 0 {
            continue;
        }
        let xe = (r.x_start + r.len).min(seq_x.len());
        if r.x_start < xe {
            x_ivs.push((r.x_start, xe));
        }
        let (ys, ye) = if r.is_rc {
            (
                (r.y_start + 1).saturating_sub(r.len),
                (r.y_start + 1).min(seq_y.len()),
            )
        } else {
            (r.y_start, (r.y_start + r.len).min(seq_y.len()))
        };
        if ys < ye {
            y_ivs.push((ys, ye));
        }
    }
    let cov_x = unique_covered_bp(&mut x_ivs);
    let cov_y = unique_covered_bp(&mut y_ivs);
    let total_len = seq_x.len() + seq_y.len();
    #[allow(clippy::cast_precision_loss)]
    let jaccard_cov = if total_len > 0 {
        (cov_x + cov_y) as f64 / total_len as f64
    } else {
        0.0
    };
    #[allow(clippy::cast_precision_loss)]
    let p_dist = if aligned_bp > 0 {
        1.0 - (matches as f64 / aligned_bp as f64)
    } else {
        1.0
    };
    let jc69_dist = if aligned_bp > 0 {
        jc69_distance(p_dist)
    } else {
        f64::INFINITY
    };
    PairMetrics {
        i,
        j,
        aligned_bp,
        matches,
        jaccard_cov,
        p_dist,
        jc69_dist,
        macro_bubbles: aln.macro_bubbles,
        internal_bubbles: aln.internal_bubbles,
    }
}

fn run_locus(cli: &LocusCli) -> Result<()> {
    let pool = cli
        .threads
        .map(|threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads.max(1))
                .build()
                .context("configuring rayon thread pool")
        })
        .transpose()?;
    if let Some(pool) = pool {
        pool.install(|| run_locus_inner(cli))
    } else {
        run_locus_inner(cli)
    }
}

fn run_locus_inner(cli: &LocusCli) -> Result<()> {
    let t0 = Instant::now();
    let contigs = read_fasta_contigs(Path::new(&cli.fasta_path), 0)
        .with_context(|| format!("reading locus FASTA '{}'", cli.fasta_path))?;
    if contigs.len() < 2 {
        bail!(
            "need at least 2 sequences in '{}' (found {})",
            cli.fasta_path,
            contigs.len()
        );
    }

    let indexed: Vec<ContigIndex> = contigs
        .par_iter()
        .map(|(_, seq)| prepare_contig_index_opts(seq, cli.respect_softmask))
        .collect::<std::result::Result<Vec<_>, _>>()?;

    let ref_idx = if let Some(ref wanted) = cli.focal_ref {
        contigs
            .iter()
            .position(|(name, _)| name == wanted)
            .or_else(|| {
                contigs
                    .iter()
                    .position(|(name, _)| horizon_aligner::matches_sample_prefix(name, wanted))
            })
            .with_context(|| format!("focal reference '{wanted}' not found in FASTA"))?
    } else {
        contigs
            .iter()
            .enumerate()
            .max_by_key(|(_, (_, seq))| seq.len())
            .map_or(0, |(idx, _)| idx)
    };

    let cfg = AlignConfig {
        skip_softmasked: cli.respect_softmask,
        ..AlignConfig::default()
    };

    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..contigs.len() {
        for j in (i + 1)..contigs.len() {
            pairs.push((i, j));
        }
    }

    let results: Vec<(PairMetrics, Option<(String, String)>)> = pairs
        .par_iter()
        .map(
            |&(i, j)| -> Result<(PairMetrics, Option<(String, String)>)> {
                if j == ref_idx {
                    let aln = align_contig_pair_runs(&indexed[j], &indexed[i], &cfg)?;
                    let metrics = compute_pair_metrics(i, j, &contigs[j].1, &contigs[i].1, &aln);
                    let maf = format_maf_blocks_runs(
                        &contigs[j].0,
                        contigs[j].1.len(),
                        &contigs[i].0,
                        contigs[i].1.len(),
                        &contigs[j].1,
                        &contigs[i].1,
                        &aln,
                    );
                    Ok((metrics, Some((contigs[i].0.clone(), maf))))
                } else {
                    let aln = align_contig_pair_runs(&indexed[i], &indexed[j], &cfg)?;
                    let metrics = compute_pair_metrics(i, j, &contigs[i].1, &contigs[j].1, &aln);
                    let focal_maf = if i == ref_idx {
                        let maf = format_maf_blocks_runs(
                            &contigs[i].0,
                            contigs[i].1.len(),
                            &contigs[j].0,
                            contigs[j].1.len(),
                            &contigs[i].1,
                            &contigs[j].1,
                            &aln,
                        );
                        Some((contigs[j].0.clone(), maf))
                    } else {
                        None
                    };
                    Ok((metrics, focal_maf))
                }
            },
        )
        .collect::<Result<Vec<_>>>()?;

    let (total_blocks, multi_blocks) =
        write_locus_outputs(&cli.out_prefix, &contigs, ref_idx, results)?;
    let ref_name = &contigs[ref_idx].0;
    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    eprintln!(
        "[horizon_locus] seqs={} pairs={} focal={} maf_blocks={} multi_blocks={} elapsed={:.2}s peak_rss_mb={:.1}",
        contigs.len(),
        pairs.len(),
        ref_name,
        total_blocks,
        multi_blocks,
        t0.elapsed().as_secs_f64(),
        hwm_mb
    );
    Ok(())
}

fn write_locus_outputs(
    out_prefix: &str,
    contigs: &[(String, Vec<u8>)],
    ref_idx: usize,
    results: Vec<(PairMetrics, Option<(String, String)>)>,
) -> Result<(usize, usize)> {
    let mut sim_tsv = String::from(
        "seq_a\tseq_b\tlen_a\tlen_b\taligned_bp\tmatches\tjaccard_cov\tp_dist\tjc69_dist\tmacro_bubbles\tinternal_bubbles\n",
    );
    let mut pairwise_mafs: Vec<(String, String)> = Vec::new();
    for (m, focal_opt) in results {
        let (ref name_a, ref seq_a) = contigs[m.i];
        let (ref name_b, ref seq_b) = contigs[m.j];
        let jc_str = if m.jc69_dist.is_finite() {
            format!("{:.6}", m.jc69_dist)
        } else {
            "NA".to_owned()
        };
        let _ = writeln!(
            sim_tsv,
            "{name_a}\t{name_b}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{jc_str}\t{}\t{}",
            seq_a.len(),
            seq_b.len(),
            m.aligned_bp,
            m.matches,
            m.jaccard_cov,
            m.p_dist,
            m.macro_bubbles,
            m.internal_bubbles
        );
        if let Some(pw) = focal_opt {
            pairwise_mafs.push(pw);
        }
    }

    let ref_name = &contigs[ref_idx].0;
    let (merged_maf, merged_gfa, total_blocks, multi_blocks) =
        merge_pairwise_mafs(ref_name, &pairwise_mafs, 500);

    fs::write(format!("{out_prefix}.similarity.tsv"), sim_tsv)
        .with_context(|| format!("writing {out_prefix}.similarity.tsv"))?;
    fs::write(format!("{out_prefix}.maf"), merged_maf)
        .with_context(|| format!("writing {out_prefix}.maf"))?;
    fs::write(format!("{out_prefix}.gfa"), merged_gfa)
        .with_context(|| format!("writing {out_prefix}.gfa"))?;
    Ok((total_blocks, multi_blocks))
}

fn parse_args() -> Result<Option<LocusCli>> {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        println!(
            "usage: horizon_locus --fasta <locus.fa[.gz]> --out-prefix <prefix> [--ref <seq_name>] [--threads <N>] [--respect-softmask]"
        );
        return Ok(None);
    }
    let mut fasta_path: Option<String> = None;
    let mut out_prefix: Option<String> = None;
    let mut focal_ref: Option<String> = None;
    let mut threads: Option<usize> = None;
    let mut respect_softmask = false;

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "--fasta" | "-i" => {
                i += 1;
                fasta_path = Some(args.get(i).context("missing value for --fasta")?.clone());
            }
            "--out-prefix" | "-o" => {
                i += 1;
                out_prefix = Some(
                    args.get(i)
                        .context("missing value for --out-prefix")?
                        .clone(),
                );
            }
            "--ref" => {
                i += 1;
                focal_ref = Some(args.get(i).context("missing value for --ref")?.clone());
            }
            "--threads" | "-t" => {
                i += 1;
                threads = Some(
                    args.get(i)
                        .context("missing value for --threads")?
                        .parse()
                        .context("invalid --threads")?,
                );
            }
            "--respect-softmask" => respect_softmask = true,
            other => bail!("unknown argument: {other}"),
        }
        i += 1;
    }

    Ok(Some(LocusCli {
        fasta_path: fasta_path.context("missing required --fasta")?,
        out_prefix: out_prefix.context("missing required --out-prefix")?,
        focal_ref,
        threads,
        respect_softmask,
    }))
}

fn main() -> Result<()> {
    let Some(cli) = parse_args()? else {
        return Ok(());
    };
    run_locus(&cli)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locus_pipeline_produces_similarity_maf_and_gfa() {
        let pid = std::process::id();
        let fa_path = std::env::temp_dir().join(format!("horizon_locus_test_{pid}.fa"));
        let out_prefix = std::env::temp_dir().join(format!("horizon_locus_out_{pid}"));

        // Deterministic non-repetitive 160 bp sequence.
        let mut base_bytes = Vec::with_capacity(160);
        let mut x: u64 = 0x1234_5678_9abc_def0;
        for _ in 0..160 {
            x = x
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            base_bytes.push(b"ACGT"[((x >> 62) & 3) as usize]);
        }
        let base = String::from_utf8(base_bytes).expect("ascii");
        let rc_bytes: Vec<u8> = base.bytes().rev().map(complement_base).collect();
        let rc_base = String::from_utf8(rc_bytes).expect("ascii");
        let mut mut_b = base.as_bytes().to_vec();
        mut_b[40] = if mut_b[40] == b'A' { b'T' } else { b'A' };
        mut_b[100] = if mut_b[100] == b'C' { b'G' } else { b'C' };
        let fa = format!(
            ">sp1\n{base}\n>sp2\n{}\n>sp3\n{base}\n>sp4_rc\n{rc_base}\n",
            String::from_utf8_lossy(&mut_b)
        );
        fs::write(&fa_path, fa).expect("write fa");

        let cli = LocusCli {
            fasta_path: fa_path.to_string_lossy().into_owned(),
            out_prefix: out_prefix.to_string_lossy().into_owned(),
            focal_ref: Some("sp1".to_owned()),
            threads: Some(2),
            respect_softmask: false,
        };
        run_locus(&cli).expect("run_locus");

        let tsv = fs::read_to_string(format!("{}.similarity.tsv", cli.out_prefix)).expect("tsv");
        let maf = fs::read_to_string(format!("{}.maf", cli.out_prefix)).expect("maf");
        let gfa = fs::read_to_string(format!("{}.gfa", cli.out_prefix)).expect("gfa");

        let cli_nonfirst_ref = LocusCli {
            fasta_path: fa_path.to_string_lossy().into_owned(),
            out_prefix: out_prefix.to_string_lossy().into_owned(),
            focal_ref: Some("sp3".to_owned()),
            threads: Some(2),
            respect_softmask: false,
        };
        run_locus(&cli_nonfirst_ref).expect("second run_locus with threads=Some(2) and j==ref_idx");
        let maf_sp3 =
            fs::read_to_string(format!("{}.maf", cli_nonfirst_ref.out_prefix)).expect("maf sp3");
        assert!(
            maf_sp3.contains("s sp3 ") && maf_sp3.contains("s sp1.sp1 "),
            "maf_sp3:\n{maf_sp3}"
        );

        let _ = fs::remove_file(&fa_path);
        let _ = fs::remove_file(format!("{}.similarity.tsv", cli.out_prefix));
        let _ = fs::remove_file(format!("{}.maf", cli.out_prefix));
        let _ = fs::remove_file(format!("{}.gfa", cli.out_prefix));

        assert!(tsv.contains("sp1\tsp2\t160\t160\t"), "tsv:\n{tsv}");
        assert!(
            tsv.contains("sp1\tsp3\t160\t160\t160\t160\t1.000000\t0.000000\t0.000000"),
            "tsv:\n{tsv}"
        );
        assert!(
            tsv.contains("sp1\tsp4_rc\t160\t160\t160\t160\t1.000000\t0.000000\t0.000000"),
            "tsv:\n{tsv}"
        );
        assert!(maf.contains("##maf version=1"));
        assert!(gfa.contains("H\tVN:Z:1.1"));
    }
}
