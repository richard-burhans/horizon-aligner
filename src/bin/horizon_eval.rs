//! Evaluate `horizon-aligner` recall, precision, and abstention against truth MAF (supports single- and multi-contig FASTAs).

use anyhow::{Context, Result};
use horizon_aligner::formats::for_each_maf_stanza_reader;
use horizon_aligner::{
    AlignConfig, align_contig_pair_runs, init_rayon_threads, parse_hsp_threshold_arg,
    parse_min_fallback_cov_arg, parse_threads_arg, prepare_contig_index_owned,
};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::time::Instant;

const USAGE: &str = "usage: horizon_eval <ref_name> <target_name> <ref.fa> <target.fa> <truth.maf> [-t N | --threads N] [--min-fallback-cov F] [--hsp-threshold N]";

fn read_named_contigs(path: &Path, target: &str) -> Result<HashMap<String, Vec<u8>>> {
    let all = horizon_aligner::read_fasta_contigs(path, 0)?;
    let has_match = all
        .iter()
        .any(|(id, _)| horizon_aligner::matches_sample_prefix(id, target));
    if has_match {
        Ok(all
            .into_iter()
            .filter(|(id, _)| horizon_aligner::matches_sample_prefix(id, target))
            .collect())
    } else {
        Ok(all.into_iter().collect())
    }
}

fn resolve_contig_index<'a, T>(
    cache: &'a HashMap<String, T>,
    maf_id: &str,
    sample_name: &str,
) -> Option<&'a T> {
    cache.get(maf_id).or_else(|| {
        let stripped = horizon_aligner::strip_sample_prefix(maf_id, sample_name);
        (stripped != maf_id).then(|| cache.get(stripped)).flatten()
    })
}

/// Canonical merged diagonal run of aligned `(x, y)` coordinate pairs:
/// - `is_rc == false`: `y = x - diag` for `x in x_start..x_end` (`diag = x_start - y_start`)
/// - `is_rc == true`:  `y = diag - x` for `x in x_start..x_end` (`diag = x_start + y_start`)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct DiagRun {
    is_rc: bool,
    diag: i128,
    x_start: usize,
    x_end: usize,
}

fn normalize_diag_runs(mut runs: Vec<DiagRun>) -> Vec<DiagRun> {
    runs.retain(|r| r.x_end > r.x_start);
    runs.sort_unstable_by_key(|r| (r.is_rc, r.diag, r.x_start, r.x_end));
    let mut merged: Vec<DiagRun> = Vec::with_capacity(runs.len());
    for r in runs {
        if let Some(last) = merged.last_mut()
            && last.is_rc == r.is_rc
            && last.diag == r.diag
            && r.x_start <= last.x_end
        {
            last.x_end = last.x_end.max(r.x_end);
        } else {
            merged.push(r);
        }
    }
    merged
}

fn total_diag_pairs(runs: &[DiagRun]) -> usize {
    runs.iter().map(|r| r.x_end - r.x_start).sum()
}

fn count_diag_overlap(a: &[DiagRun], b: &[DiagRun]) -> usize {
    let (mut i, mut j, mut ov) = (0usize, 0usize, 0usize);
    while i < a.len() && j < b.len() {
        let ra = a[i];
        let rb = b[j];
        match (ra.is_rc, ra.diag).cmp(&(rb.is_rc, rb.diag)) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                ov += ra
                    .x_end
                    .min(rb.x_end)
                    .saturating_sub(ra.x_start.max(rb.x_start));
                if ra.x_end <= rb.x_end {
                    i += 1;
                } else {
                    j += 1;
                }
            }
        }
    }
    ov
}

type MafPairMap = HashMap<(String, String), Vec<DiagRun>>;
type MafEntry = ((String, String), Vec<DiagRun>);

#[allow(clippy::similar_names)]
fn parse_maf_pairs_by_contig(path: &Path, cx: &str, cy: &str) -> Result<MafPairMap> {
    let reader = horizon_aligner::io_util::open_text_reader(path)?;
    let mut out: MafPairMap = HashMap::new();
    let mut mismatch_err: Option<anyhow::Error> = None;
    for_each_maf_stanza_reader(reader, |rows| {
        if mismatch_err.is_some() {
            return;
        }
        for (i, (rx, _)) in rows.iter().enumerate() {
            if !horizon_aligner::matches_sample_prefix(&rx.src, cx) {
                continue;
            }
            for (j, (ry, _)) in rows.iter().enumerate() {
                if i == j || !horizon_aligner::matches_sample_prefix(&ry.src, cy) {
                    continue;
                }
                let ax = rx.seq.as_bytes();
                let ay = ry.seq.as_bytes();
                if ax.len() != ay.len() {
                    mismatch_err = Some(anyhow::anyhow!(
                        "MAF block alignment length mismatch ({} vs {})",
                        ax.len(),
                        ay.len()
                    ));
                    return;
                }
                let runs = out.entry((rx.src.clone(), ry.src.clone())).or_default();
                let is_rc = rx.strand != ry.strand;
                let (mut ix, mut iy) = (0usize, 0usize);
                let mut cur_run: Option<(usize, usize, usize, usize)> = None;
                let flush_run =
                    |runs: &mut Vec<DiagRun>, state: Option<(usize, usize, usize, usize)>| {
                        if let Some((px0, py0, px_last, len)) = state {
                            let x_start = px0.min(px_last);
                            let diag = if is_rc {
                                (px0 as i128) + (py0 as i128)
                            } else {
                                (px0 as i128) - (py0 as i128)
                            };
                            runs.push(DiagRun {
                                is_rc,
                                diag,
                                x_start,
                                x_end: x_start + len,
                            });
                        }
                    };
                for (&ca, &cb) in ax.iter().zip(ay) {
                    if ca != b'-' && cb != b'-' {
                        let px = rx.fwd_pos(ix);
                        let py = ry.fwd_pos(iy);
                        if let Some((_, _, ref mut px_last, ref mut len)) = cur_run {
                            *px_last = px;
                            *len += 1;
                        } else {
                            cur_run = Some((px, py, px, 1));
                        }
                    } else if cur_run.is_some() {
                        flush_run(runs, cur_run.take());
                    }
                    if ca != b'-' {
                        ix += 1;
                    }
                    if cb != b'-' {
                        iy += 1;
                    }
                }
                flush_run(runs, cur_run);
            }
        }
    })?;
    if let Some(err) = mismatch_err {
        return Err(err);
    }
    for runs in out.values_mut() {
        *runs = normalize_diag_runs(std::mem::take(runs));
    }
    Ok(out)
}

struct EvalCliArgs {
    ref_name: String,
    target_name: String,
    ref_fa: std::path::PathBuf,
    tgt_fa: std::path::PathBuf,
    truth_maf: std::path::PathBuf,
    threads: Option<usize>,
    min_fallback_cov: Option<f64>,
    hsp_threshold: Option<i32>,
}

fn parse_eval_args() -> Result<Option<EvalCliArgs>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut pos = Vec::new();
    let mut threads = None;
    let mut min_fallback_cov = None;
    let mut hsp_threshold = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                println!("{USAGE}");
                return Ok(None);
            }
            "-V" | "--version" => {
                println!("horizon_eval {}", env!("CARGO_PKG_VERSION"));
                return Ok(None);
            }
            "-t" | "--threads" => {
                let raw = args.get(i + 1).context("missing value for --threads")?;
                threads = Some(parse_threads_arg(raw)?);
                i += 2;
            }
            "--min-fallback-cov" => {
                let raw = args
                    .get(i + 1)
                    .context("missing value for --min-fallback-cov")?;
                min_fallback_cov = Some(parse_min_fallback_cov_arg(raw)?);
                i += 2;
            }
            "--hsp-threshold" => {
                let raw = args
                    .get(i + 1)
                    .context("missing value for --hsp-threshold")?;
                hsp_threshold = Some(parse_hsp_threshold_arg(raw)?);
                i += 2;
            }
            s if s.starts_with('-') && s != "-" => anyhow::bail!("unknown argument '{s}'"),
            _ => {
                pos.push(args[i].as_str());
                i += 1;
            }
        }
    }
    if pos.len() != 5 {
        anyhow::bail!("{USAGE}");
    }
    Ok(Some(EvalCliArgs {
        ref_name: pos[0].to_string(),
        target_name: pos[1].to_string(),
        ref_fa: pos[2].into(),
        tgt_fa: pos[3].into(),
        truth_maf: pos[4].into(),
        threads,
        min_fallback_cov,
        hsp_threshold,
    }))
}

fn collect_needed_contig_names<'a>(
    maf_by_pair: &'a MafPairMap,
    sample_name: &'a str,
    use_ref: bool,
) -> HashSet<&'a str> {
    maf_by_pair
        .keys()
        .flat_map(|(r, t)| {
            let id = if use_ref { r.as_str() } else { t.as_str() };
            [id, horizon_aligner::strip_sample_prefix(id, sample_name)]
        })
        .collect()
}

type ContigPairTruth<'a> = (
    &'a horizon_aligner::ContigIndex,
    &'a horizon_aligner::ContigIndex,
    &'a [DiagRun],
);

fn align_and_score_pairs(
    aligned_pairs: Vec<ContigPairTruth<'_>>,
    cfg: &AlignConfig,
) -> Result<(usize, usize, usize, usize)> {
    let pair_stats: Vec<Result<(usize, usize, usize, usize)>> = aligned_pairs
        .into_par_iter()
        .map(|(idx_r, idx_t, maf_runs)| {
            let aln = align_contig_pair_runs(idx_r, idx_t, cfg)?;
            let proto_runs = normalize_diag_runs(
                aln.aligned_runs
                    .iter()
                    .filter(|r| r.len > 0)
                    .map(|r| DiagRun {
                        is_rc: r.is_rc,
                        diag: if r.is_rc {
                            (r.x_start as i128) + (r.y_start as i128)
                        } else {
                            (r.x_start as i128) - (r.y_start as i128)
                        },
                        x_start: r.x_start,
                        x_end: r.x_start + r.len,
                    })
                    .collect(),
            );
            let overlap = count_diag_overlap(&proto_runs, maf_runs);
            Ok((
                total_diag_pairs(maf_runs),
                total_diag_pairs(&proto_runs),
                overlap,
                aln.macro_bubbles + aln.internal_bubbles,
            ))
        })
        .collect();

    let mut total_maf = 0usize;
    let mut total_proto = 0usize;
    let mut total_overlap = 0usize;
    let mut total_bubbles = 0usize;
    for res in pair_stats {
        let (maf_len, proto_len, overlap, bubbles) = res?;
        total_maf += maf_len;
        total_proto += proto_len;
        total_overlap += overlap;
        total_bubbles += bubbles;
    }
    Ok((total_maf, total_proto, total_overlap, total_bubbles))
}

fn main() -> Result<()> {
    let Some(cli) = parse_eval_args()? else {
        return Ok(());
    };
    init_rayon_threads(cli.threads)?;
    let ref_contigs = read_named_contigs(&cli.ref_fa, &cli.ref_name)?;
    let tgt_contigs = read_named_contigs(&cli.tgt_fa, &cli.target_name)?;
    let maf_by_pair = parse_maf_pairs_by_contig(&cli.truth_maf, &cli.ref_name, &cli.target_name)?;

    let cfg = AlignConfig {
        min_fallback_cov: cli.min_fallback_cov,
        hsp_threshold: cli.hsp_threshold,
        ..AlignConfig::default()
    };

    let t0 = Instant::now();
    let needed_refs = collect_needed_contig_names(&maf_by_pair, &cli.ref_name, true);
    let needed_tgts = collect_needed_contig_names(&maf_by_pair, &cli.target_name, false);
    let ref_idx_cache: HashMap<String, _> = ref_contigs
        .into_par_iter()
        .filter(|(id, _)| needed_refs.contains(id.as_str()))
        .map(|(id, seq)| Ok((id, prepare_contig_index_owned(seq, false)?)))
        .collect::<Result<_>>()?;
    let tgt_idx_cache: HashMap<String, _> = tgt_contigs
        .into_par_iter()
        .filter(|(id, _)| needed_tgts.contains(id.as_str()))
        .map(|(id, seq)| Ok((id, prepare_contig_index_owned(seq, false)?)))
        .collect::<Result<_>>()?;

    let mut pairs: Vec<MafEntry> = maf_by_pair.into_iter().collect();
    pairs.sort_by(|a, b| a.0.cmp(&b.0));

    let mut aligned_pairs = Vec::with_capacity(pairs.len());
    for ((ref_id, tgt_id), maf_runs) in &pairs {
        let (Some(idx_r), Some(idx_t)) = (
            resolve_contig_index(&ref_idx_cache, ref_id, &cli.ref_name),
            resolve_contig_index(&tgt_idx_cache, tgt_id, &cli.target_name),
        ) else {
            eprintln!(
                "horizon_eval: warning: skipping MAF pair ({ref_id}, {tgt_id}) missing from FASTA inputs ({} pairs)",
                total_diag_pairs(maf_runs)
            );
            continue;
        };
        aligned_pairs.push((idx_r, idx_t, maf_runs.as_slice()));
    }
    anyhow::ensure!(
        pairs.is_empty() || !aligned_pairs.is_empty(),
        "no truth MAF contig pairs matched the provided FASTA inputs"
    );

    let (total_maf, total_proto, total_overlap, total_bubbles) =
        align_and_score_pairs(aligned_pairs, &cfg)?;

    let elapsed = t0.elapsed().as_secs_f64();
    #[allow(clippy::cast_precision_loss)]
    let recall = if total_maf == 0 {
        0.0
    } else {
        100.0 * total_overlap as f64 / total_maf as f64
    };
    #[allow(clippy::cast_precision_loss)]
    let precision = if total_proto == 0 {
        0.0
    } else {
        100.0 * total_overlap as f64 / total_proto as f64
    };
    println!(
        "elapsed={elapsed:.3}s proto_pairs={total_proto} maf_pairs={total_maf} overlap={total_overlap} raw_maf_recall={recall:.1}% raw_maf_precision={precision:.1}% bubbles={total_bubbles}"
    );
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use std::fs;

    fn diag_contains(runs: &[DiagRun], x: usize, y: usize) -> bool {
        runs.iter().any(|r| {
            x >= r.x_start
                && x < r.x_end
                && if r.is_rc {
                    (x as i128) + (y as i128) == r.diag
                } else {
                    (x as i128) - (y as i128) == r.diag
                }
        })
    }

    #[test]
    fn parses_named_contigs_and_maf_pairs() {
        let dir = std::env::temp_dir();
        let pid = std::process::id();
        let fa_path = dir.join(format!("horizon_eval_test_{pid}.fa"));
        let maf_path = dir.join(format!("horizon_eval_test_{pid}.maf"));

        fs::write(&fa_path, ">hg38.chr22\nACGT\nTGCA\n>other\nAAAA\n").unwrap();
        let contigs = read_named_contigs(&fa_path, "hg38").unwrap();
        assert_eq!(
            contigs.get("hg38.chr22").map(Vec::as_slice),
            Some(&b"ACGTTGCA"[..])
        );
        assert!(!contigs.contains_key("other"));

        fs::write(
            &maf_path,
            "\
a score=10
s hg38.chr22 10 3 + 100 ACG
s dog.chr1 5 3 - 50 TGC

a score=10
s hg38.chr22 20 3 - 100 ACG
s dog.chr1 10 3 + 50 ACG

a score=10
s hg38.chr22 30 3 - 100 ACG
s dog.chr1 15 3 - 50 ACG
",
        )
        .unwrap();
        let maf = parse_maf_pairs_by_contig(&maf_path, "hg38", "dog").unwrap();
        let pairs = maf
            .get(&("hg38.chr22".to_string(), "dog.chr1".to_string()))
            .unwrap();
        assert_eq!(total_diag_pairs(pairs), 9);
        // (+, -)
        assert!(diag_contains(pairs, 10, 44));
        assert!(diag_contains(pairs, 11, 43));
        assert!(diag_contains(pairs, 12, 42));
        // (-, +): hg38 fwd positions 79, 78, 77 paired with dog 10, 11, 12 (is_rc = true)
        assert!(diag_contains(pairs, 79, 10));
        assert!(diag_contains(pairs, 78, 11));
        assert!(diag_contains(pairs, 77, 12));
        // (-, -): hg38 fwd positions 69, 68, 67 paired with dog 34, 33, 32 (is_rc = false)
        assert!(diag_contains(pairs, 69, 34));
        assert!(diag_contains(pairs, 68, 33));
        assert!(diag_contains(pairs, 67, 32));

        let _ = fs::remove_file(&fa_path);
        let _ = fs::remove_file(&maf_path);
    }

    #[test]
    fn read_named_contigs_falls_back_to_unprefixed_headers() {
        let dir = std::env::temp_dir();
        let pid = std::process::id();
        let fa_path = dir.join(format!("horizon_eval_unprefixed_{pid}.fa"));
        fs::write(&fa_path, ">chr1\nACGT\n").unwrap();
        let contigs = read_named_contigs(&fa_path, "hg38").unwrap();
        let _ = fs::remove_file(&fa_path);
        assert_eq!(contigs.get("chr1").map(Vec::as_slice), Some(&b"ACGT"[..]));
        assert!(resolve_contig_index(&contigs, "hg38.chr1", "hg38").is_some());
    }

    #[test]
    fn parse_maf_pairs_discards_malformed_stanza_and_handles_paralogs_and_self_align() {
        let dir = std::env::temp_dir();
        let pid = std::process::id();
        let maf_path = dir.join(format!("horizon_eval_paralog_{pid}.maf"));
        // Stanza 1: malformed second s-line -> entire stanza must be ignored.
        // Stanza 2: hg38.chr1 aligned to TWO dog paralogs -> both pairs extracted.
        // Stanza 3: self-species (hg38.chr1 vs hg38.chr2) -> must NOT self-pair chr1 with chr1.
        fs::write(
            &maf_path,
            "a score=10\ns hg38.chr1 0 3 + 100 ACG\ns dog.chr1 0 999 + 10 ACG\n\n\
             a score=20\ns hg38.chr1 10 3 + 100 ACG\ns dog.chr1 0 3 + 50 ACG\ns dog.chr2 5 3 + 50 ACG\n\n\
             a score=30\ns hg38.chr1 20 3 + 100 ACG\ns hg38.chr2 30 3 + 100 ACG\n\n",
        )
        .unwrap();

        let maf_hd = parse_maf_pairs_by_contig(&maf_path, "hg38", "dog").unwrap();
        let p1 = maf_hd
            .get(&("hg38.chr1".to_string(), "dog.chr1".to_string()))
            .unwrap();
        assert!(
            !diag_contains(p1, 0, 0),
            "malformed stanza must be rejected"
        );
        assert!(diag_contains(p1, 10, 0));
        let p2 = maf_hd
            .get(&("hg38.chr1".to_string(), "dog.chr2".to_string()))
            .expect("second paralog row must be extracted");
        assert!(diag_contains(p2, 10, 5));

        let maf_self = parse_maf_pairs_by_contig(&maf_path, "hg38", "hg38").unwrap();
        assert!(
            !maf_self.contains_key(&("hg38.chr1".to_string(), "hg38.chr1".to_string())),
            "cx == cy must not self-pair a row with itself"
        );
        assert!(
            maf_self
                .get(&("hg38.chr1".to_string(), "hg38.chr2".to_string()))
                .is_some_and(|s| diag_contains(s, 20, 30))
        );
        let _ = fs::remove_file(&maf_path);
    }

    #[test]
    fn diag_run_overlap_matches_exact_pair_intersection_across_strands() {
        let a = normalize_diag_runs(vec![
            DiagRun {
                is_rc: false,
                diag: 10,
                x_start: 20,
                x_end: 30,
            },
            DiagRun {
                is_rc: false,
                diag: 10,
                x_start: 25,
                x_end: 35,
            },
            DiagRun {
                is_rc: true,
                diag: 100,
                x_start: 10,
                x_end: 20,
            },
        ]);
        let b = normalize_diag_runs(vec![
            DiagRun {
                is_rc: false,
                diag: 10,
                x_start: 28,
                x_end: 40,
            },
            DiagRun {
                is_rc: true,
                diag: 100,
                x_start: 15,
                x_end: 25,
            },
            DiagRun {
                is_rc: true,
                diag: 99,
                x_start: 10,
                x_end: 20,
            },
        ]);
        // Forward diag 10: [20,35) ∩ [28,40) = 7 bp; RC diag 100: [10,20) ∩ [15,25) = 5 bp.
        assert_eq!(total_diag_pairs(&a), 25);
        assert_eq!(total_diag_pairs(&b), 32);
        assert_eq!(count_diag_overlap(&a, &b), 12);
    }
}
