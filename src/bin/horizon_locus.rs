//! Regional `impg`-style multi-FASTA locus aligner, `BioNJ` tree builder, and MAF/GFA/HNEX emitter.
//!
//! Designed for VGP `O3` workflows (25 kb ILS windows and BUSCO gene intervals): takes a
//! multi-FASTA of homologous sequences, runs pairwise SDUST + BOCPD/MDL + `Aa5Index` alignment
//! across all pairs, and outputs `<out_prefix>.similarity.tsv`, `<out_prefix>.bionj.nwk`,
//! `<out_prefix>.maf`, `<out_prefix>.gfa`, optional lossless `<out_prefix>.hnex` (`--hnex`),
//! and optional banded `ProbCons` consistency MSA (`--consistency`).

use anyhow::{Context, Result, bail};
use horizon_aligner::formats::{format_maf_blocks_runs, read_fasta_contigs};
use horizon_aligner::mdl::{AlignConfig, RunAlignment, align_contig_pair_runs};
use horizon_aligner::merge::merge_pairwise_mafs;
use horizon_aligner::seeding::{ContigIndex, complement_base, prepare_contig_index_owned};
use rayon::prelude::*;
use std::env;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;
use std::time::Instant;

#[allow(clippy::struct_excessive_bools)]
struct LocusCli {
    fasta_path: String,
    out_prefix: String,
    focal_ref: Option<String>,
    threads: Option<usize>,
    respect_softmask: bool,
    emit_hnex: bool,
    consistency: bool,
    repair_frameshifts: bool,
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

type PairRunResult = (PairMetrics, Option<(String, String)>, Option<(usize, bool)>);

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

fn compute_pair_metrics(
    i: usize,
    j: usize,
    seq_x: &[u8],
    seq_y: &[u8],
    aln: &RunAlignment,
) -> PairMetrics {
    let aligned_bp = aln.total_pairs;
    let matches = count_run_matches(seq_x, seq_y, aln);
    let mut cov_x_set = horizon_aligner::coord::IntervalSet::new();
    let mut cov_y_set = horizon_aligner::coord::IntervalSet::new();
    let mut scratch: Vec<usize> = Vec::new();
    for r in &aln.aligned_runs {
        if r.len == 0 {
            continue;
        }
        let xe = (r.x_start + r.len).min(seq_x.len());
        if r.x_start < xe {
            cov_x_set.insert(r.x_start, xe, &mut scratch);
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
            cov_y_set.insert(ys, ye, &mut scratch);
        }
    }
    let cov_x = cov_x_set.covered_bp();
    let cov_y = cov_y_set.covered_bp();
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

#[allow(clippy::too_many_lines)]
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
    let mut seen_names = std::collections::HashSet::with_capacity(contigs.len());
    for (name, _) in &contigs {
        anyhow::ensure!(
            !name.trim().is_empty(),
            "empty sequence name in '{}'",
            cli.fasta_path
        );
        anyhow::ensure!(
            seen_names.insert(name.as_str()),
            "duplicate sequence name `{name}` in '{}'",
            cli.fasta_path
        );
    }

    let indexed: Vec<(String, ContigIndex)> = contigs
        .into_par_iter()
        .map(|(name, seq)| Ok((name, prepare_contig_index_owned(seq, cli.respect_softmask)?)))
        .collect::<Result<Vec<_>>>()?;

    let ref_idx = if let Some(ref wanted) = cli.focal_ref {
        indexed
            .iter()
            .position(|(name, _)| name == wanted)
            .or_else(|| {
                indexed
                    .iter()
                    .position(|(name, _)| horizon_aligner::matches_sample_prefix(name, wanted))
            })
            .with_context(|| format!("focal reference '{wanted}' not found in FASTA"))?
    } else {
        indexed
            .iter()
            .enumerate()
            .max_by_key(|(_, (_, ix))| ix.seq().len())
            .map_or(0, |(idx, _)| idx)
    };

    let cfg = AlignConfig {
        skip_softmasked: cli.respect_softmask,
        ..AlignConfig::default()
    };

    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..indexed.len() {
        for j in (i + 1)..indexed.len() {
            pairs.push((i, j));
        }
    }

    let majority_rc = |aln: &RunAlignment| -> bool {
        let mut fwd_bp = 0usize;
        let mut rc_bp = 0usize;
        for r in &aln.aligned_runs {
            if r.is_rc {
                rc_bp += r.len;
            } else {
                fwd_bp += r.len;
            }
        }
        rc_bp > fwd_bp
    };

    let results: Vec<PairRunResult> = pairs
        .par_iter()
        .map(|&(i, j)| -> Result<PairRunResult> {
            if j == ref_idx {
                let aln = align_contig_pair_runs(&indexed[j].1, &indexed[i].1, &cfg)?;
                let sx = indexed[j].1.seq();
                let sy = indexed[i].1.seq();
                let metrics = compute_pair_metrics(i, j, sx, sy, &aln);
                let is_rc = majority_rc(&aln);
                let maf = format_maf_blocks_runs(
                    &indexed[j].0,
                    sx.len(),
                    &indexed[i].0,
                    sy.len(),
                    sx,
                    sy,
                    &aln,
                );
                Ok((metrics, Some((indexed[i].0.clone(), maf)), Some((i, is_rc))))
            } else {
                let aln = align_contig_pair_runs(&indexed[i].1, &indexed[j].1, &cfg)?;
                let sx = indexed[i].1.seq();
                let sy = indexed[j].1.seq();
                let metrics = compute_pair_metrics(i, j, sx, sy, &aln);
                let (focal_maf, focal_orient) = if i == ref_idx {
                    let is_rc = majority_rc(&aln);
                    let maf = format_maf_blocks_runs(
                        &indexed[i].0,
                        sx.len(),
                        &indexed[j].0,
                        sy.len(),
                        sx,
                        sy,
                        &aln,
                    );
                    (Some((indexed[j].0.clone(), maf)), Some((j, is_rc)))
                } else {
                    (None, None)
                };
                Ok((metrics, focal_maf, focal_orient))
            }
        })
        .collect::<Result<Vec<_>>>()?;

    let (total_blocks, multi_blocks) = write_locus_outputs(
        &cli.out_prefix,
        &indexed,
        ref_idx,
        results,
        cli.emit_hnex,
        cli.consistency,
        cli.repair_frameshifts,
    )?;
    let ref_name = &indexed[ref_idx].0;
    #[allow(clippy::cast_precision_loss)]
    let hwm_mb = horizon_aligner::read_vm_hwm_kb().map_or(0.0, |kb| kb as f64 / 1024.0);
    eprintln!(
        "[horizon_locus] seqs={} pairs={} focal={} maf_blocks={} multi_blocks={} elapsed={:.2}s peak_rss_mb={:.1}",
        indexed.len(),
        pairs.len(),
        ref_name,
        total_blocks,
        multi_blocks,
        t0.elapsed().as_secs_f64(),
        hwm_mb
    );
    Ok(())
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn write_locus_outputs(
    out_prefix: &str,
    indexed: &[(String, ContigIndex)],
    ref_idx: usize,
    results: Vec<PairRunResult>,
    emit_hnex: bool,
    consistency: bool,
    repair_frameshifts: bool,
) -> Result<(usize, usize)> {
    let n_seqs = indexed.len();
    let mut dist_mat = vec![vec![0.0_f64; n_seqs]; n_seqs];
    let mut seq_is_rc = vec![false; n_seqs];
    let mut sim_tsv = String::from(
        "seq_a\tseq_b\tlen_a\tlen_b\taligned_bp\tmatches\tjaccard_cov\tp_dist\tjc69_dist\tmacro_bubbles\tinternal_bubbles\n",
    );
    let mut pairwise_mafs: Vec<(String, String)> = Vec::new();
    for (m, focal_opt, orient_opt) in results {
        let (ref name_a, ref ix_a) = indexed[m.i];
        let (ref name_b, ref ix_b) = indexed[m.j];
        if let Some((idx, is_rc)) = orient_opt {
            seq_is_rc[idx] = is_rc;
        }
        let d_val = if m.jc69_dist.is_finite() {
            m.jc69_dist
        } else {
            2.5
        };
        dist_mat[m.i][m.j] = d_val;
        dist_mat[m.j][m.i] = d_val;
        let jc_str = if m.jc69_dist.is_finite() {
            format!("{:.6}", m.jc69_dist)
        } else {
            "NA".to_owned()
        };
        let _ = writeln!(
            sim_tsv,
            "{name_a}\t{name_b}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{jc_str}\t{}\t{}",
            ix_a.seq().len(),
            ix_b.seq().len(),
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

    let ref_name = &indexed[ref_idx].0;
    let (merged_maf, merged_gfa, total_blocks, multi_blocks) =
        merge_pairwise_mafs(ref_name, &pairwise_mafs, 500);

    let seq_names: Vec<String> = indexed.iter().map(|(n, _)| n.clone()).collect();
    let bionj_newick = horizon_aligner::tba::build_bionj_tree(&seq_names, &dist_mat)
        .context("building BioNJ guide tree")?
        .to_newick();
    fs::write(
        format!("{out_prefix}.bionj.nwk"),
        format!("{bionj_newick}\n"),
    )
    .with_context(|| format!("writing {out_prefix}.bionj.nwk"))?;

    fs::write(format!("{out_prefix}.similarity.tsv"), sim_tsv)
        .with_context(|| format!("writing {out_prefix}.similarity.tsv"))?;
    fs::write(format!("{out_prefix}.maf"), &merged_maf)
        .with_context(|| format!("writing {out_prefix}.maf"))?;
    fs::write(format!("{out_prefix}.gfa"), merged_gfa)
        .with_context(|| format!("writing {out_prefix}.gfa"))?;
    if emit_hnex {
        let tba_blocks = horizon_aligner::tba::parse_multi_tba_reader(merged_maf.as_bytes())
            .context("parsing merged MAF for --hnex")?;
        let archive = horizon_aligner::formats::HnexArchive::from_tba_blocks(
            &tba_blocks,
            Some(bionj_newick),
            Vec::new(),
        );
        let mut hnex_buf = Vec::new();
        archive.write_to(&mut hnex_buf)?;
        fs::write(format!("{out_prefix}.hnex"), hnex_buf)
            .with_context(|| format!("writing {out_prefix}.hnex"))?;
    }
    if consistency {
        // Orient every locus sequence to the focal reference strand before ProbCons MSA so
        // reverse-complement homologs align accurately rather than on the opposite strand.
        let oriented_owned: Vec<Vec<u8>> = indexed
            .iter()
            .enumerate()
            .map(|(idx, (_, ix))| {
                if seq_is_rc[idx] {
                    horizon_aligner::seeding::revcomp(ix.seq())
                } else {
                    ix.seq().to_vec()
                }
            })
            .collect();
        let seq_slices: Vec<&[u8]> = oriented_owned.iter().map(Vec::as_slice).collect();
        let msa = horizon_aligner::mdl::align_locus_probcons_msa(
            &seq_slices,
            &horizon_aligner::mdl::HmmParams::default(),
            2,
        );
        let mut msa_fa = String::new();
        for (idx, ((name, _), row)) in indexed.iter().zip(&msa).enumerate() {
            let strand_note = if seq_is_rc[idx] { " strand=-" } else { "" };
            let _ = writeln!(
                msa_fa,
                ">{name}{strand_note}\n{}",
                String::from_utf8_lossy(row)
            );
        }
        fs::write(format!("{out_prefix}.probcons.fa"), msa_fa)
            .with_context(|| format!("writing {out_prefix}.probcons.fa"))?;
    }
    if repair_frameshifts {
        let ref_seq = indexed[ref_idx].1.seq();
        let ref_rem = ref_seq.len() % 3;
        let ref_cds = &ref_seq[..ref_seq.len() - ref_rem];
        let mut repaired_fa = String::new();
        let mut fs_tsv =
            String::from("seq_name\tstrand\tref_nt_pos\tquery_nt_pos\tshift_nt\tscore_bits\n");
        let mut total_events = 0usize;
        let mut skipped_budget = 0usize;
        let _ = writeln!(
            repaired_fa,
            ">{ref_name}\n{}",
            String::from_utf8_lossy(ref_seq)
        );
        for (idx, (name, ix)) in indexed.iter().enumerate() {
            if idx == ref_idx {
                continue;
            }
            let is_rc = seq_is_rc[idx];
            let strand_ch = if is_rc { '-' } else { '+' };
            let strand_note = if is_rc { " strand=-" } else { "" };
            let oriented = if is_rc {
                horizon_aligner::seeding::revcomp(ix.seq())
            } else {
                ix.seq().to_vec()
            };
            // Avoid spurious trailing +1/+2 frameshift calls when both reference and query carry
            // the same non-codon trailing overhang.
            let q_trim = if ref_rem > 0 && oriented.len() % 3 == ref_rem && oriented.len() > ref_rem
            {
                oriented.len() - ref_rem
            } else {
                oriented.len()
            };
            let rep =
                horizon_aligner::mdl::repair_coding_frameshifts(ref_cds, &oriented[..q_trim], 8.0);
            if rep.repaired_query_cds.is_empty() && !ref_cds.is_empty() && q_trim > 0 {
                skipped_budget += 1;
                eprintln!(
                    "[horizon_locus --repair-frameshifts] warning: skipped `{name}` (len={q_trim} exceeds codon DP cell cap)"
                );
            }
            total_events += rep.events.len();
            for ev in &rep.events {
                let _ = writeln!(
                    fs_tsv,
                    "{name}\t{strand_ch}\t{}\t{}\t{}\t{:.2}",
                    ev.ref_nt_pos, ev.query_nt_pos, ev.shift_nt, rep.score_bits
                );
            }
            let out_seq = if rep.repaired_query_cds.is_empty() {
                &oriented
            } else {
                &rep.repaired_query_cds
            };
            let _ = writeln!(
                repaired_fa,
                ">{name}{strand_note}\n{}",
                String::from_utf8_lossy(out_seq)
            );
        }
        fs::write(format!("{out_prefix}.repaired.fa"), repaired_fa)
            .with_context(|| format!("writing {out_prefix}.repaired.fa"))?;
        fs::write(format!("{out_prefix}.frameshifts.tsv"), fs_tsv)
            .with_context(|| format!("writing {out_prefix}.frameshifts.tsv"))?;
        eprintln!(
            "[horizon_locus --repair-frameshifts] focal={ref_name} events={total_events} skipped={skipped_budget} out={out_prefix}.{{repaired.fa,frameshifts.tsv}}"
        );
    }
    Ok((total_blocks, multi_blocks))
}

fn parse_args() -> Result<Option<LocusCli>> {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() {
        bail!(
            "usage: horizon_locus --fasta <locus.fa[.gz]> --out-prefix <prefix> [--ref <seq_name>] [-t N | --threads N] [--respect-softmask] [--hnex] [--consistency] [--repair-frameshifts]"
        );
    }
    let mut fasta_path: Option<String> = None;
    let mut out_prefix: Option<String> = None;
    let mut focal_ref: Option<String> = None;
    let mut threads: Option<usize> = None;
    let mut respect_softmask = false;
    let mut emit_hnex = false;
    let mut consistency = false;
    let mut repair_frameshifts = false;

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                println!(
                    "usage: horizon_locus --fasta <locus.fa[.gz]> --out-prefix <prefix> [--ref <seq_name>] [--threads <N>] [--respect-softmask] [--hnex] [--consistency] [--repair-frameshifts]"
                );
                return Ok(None);
            }
            "-V" | "--version" => {
                println!("horizon_locus {}", env!("CARGO_PKG_VERSION"));
                return Ok(None);
            }
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
                let raw = args.get(i).context("missing value for --threads")?;
                threads = Some(horizon_aligner::parse_threads_arg(raw)?);
            }
            "--respect-softmask" => respect_softmask = true,
            "--hnex" => emit_hnex = true,
            "--consistency" => consistency = true,
            "--repair-frameshifts" => repair_frameshifts = true,
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
        emit_hnex,
        consistency,
        repair_frameshifts,
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
            emit_hnex: true,
            consistency: true,
            repair_frameshifts: false,
        };
        run_locus(&cli).expect("run_locus");

        let tsv = fs::read_to_string(format!("{}.similarity.tsv", cli.out_prefix)).expect("tsv");
        let maf = fs::read_to_string(format!("{}.maf", cli.out_prefix)).expect("maf");
        let gfa = fs::read_to_string(format!("{}.gfa", cli.out_prefix)).expect("gfa");
        let hnex = fs::read_to_string(format!("{}.hnex", cli.out_prefix)).expect("hnex");
        let bionj = fs::read_to_string(format!("{}.bionj.nwk", cli.out_prefix)).expect("bionj");
        let probcons =
            fs::read_to_string(format!("{}.probcons.fa", cli.out_prefix)).expect("probcons");
        assert!(hnex.starts_with("HNEX\tVN:1.0\n"), "hnex:\n{hnex}");
        let parsed_hnex =
            horizon_aligner::formats::HnexArchive::read_from(hnex.as_bytes()).expect("parse hnex");
        assert!(
            !parsed_hnex.blocks.is_empty(),
            "hnex archive must contain aligned blocks"
        );
        assert!(bionj.trim().ends_with(';'), "bionj:\n{bionj}");
        assert!(probcons.contains(">sp1\n"), "probcons:\n{probcons}");
        assert!(
            probcons.contains(">sp4_rc strand=-\n"),
            "sp4_rc must be reverse-complemented to focal strand before ProbCons MSA:\n{probcons}"
        );
        // Verify sp1 and oriented sp4_rc align with exact 100% identity (same sequence).
        let mut seqs_by_hdr = std::collections::HashMap::new();
        let mut cur_hdr = String::new();
        for line in probcons.lines() {
            if let Some(h) = line.strip_prefix('>') {
                cur_hdr = h.split_whitespace().next().unwrap_or("").to_owned();
            } else if !line.is_empty() {
                seqs_by_hdr.insert(cur_hdr.clone(), line.to_owned());
            }
        }
        assert_eq!(seqs_by_hdr.len(), 4);
        assert_eq!(seqs_by_hdr["sp1"], seqs_by_hdr["sp4_rc"]);

        let cli_nonfirst_ref = LocusCli {
            fasta_path: fa_path.to_string_lossy().into_owned(),
            out_prefix: out_prefix.to_string_lossy().into_owned(),
            focal_ref: Some("sp3".to_owned()),
            threads: Some(2),
            respect_softmask: false,
            emit_hnex: false,
            consistency: false,
            repair_frameshifts: false,
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
        let _ = fs::remove_file(format!("{}.hnex", cli.out_prefix));
        let _ = fs::remove_file(format!("{}.bionj.nwk", cli.out_prefix));
        let _ = fs::remove_file(format!("{}.probcons.fa", cli.out_prefix));

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
