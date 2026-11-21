//! Evaluate `horizon-aligner` recall, precision, and abstention against truth MAF (supports single- and multi-contig FASTAs).

use anyhow::{Context, Result};
use horizon_aligner::{
    AlignConfig, MafRow, align_contig_pair_with_config, parse_hsp_threshold_arg, parse_maf_s_line,
    parse_min_fallback_cov_arg, prepare_contig_index_opts,
};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::BufRead as _;
use std::path::Path;
use std::time::Instant;

const USAGE: &str = "usage: horizon_eval <ref_name> <target_name> <ref.fa> <target.fa> <truth.maf> [--min-fallback-cov F] [--hsp-threshold N]";

fn read_named_contigs(path: &Path, target: &str) -> Result<HashMap<String, Vec<u8>>> {
    let all = horizon_aligner::read_fasta_contigs(path, usize::MAX)?;
    let matched: HashMap<String, Vec<u8>> = all
        .iter()
        .filter(|(id, _)| horizon_aligner::matches_sample_prefix(id, target))
        .map(|(id, seq)| (id.clone(), seq.clone()))
        .collect();
    if matched.is_empty() {
        Ok(all.into_iter().collect())
    } else {
        Ok(matched)
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

type MafPairMap = HashMap<(String, String), HashSet<(usize, usize)>>;
type MafEntry = ((String, String), HashSet<(usize, usize)>);

#[allow(clippy::similar_names)]
fn parse_maf_pairs_by_contig(path: &Path, cx: &str, cy: &str) -> Result<MafPairMap> {
    let mut reader = horizon_aligner::io_util::open_text_reader(path)?;
    let mut out: MafPairMap = HashMap::new();
    let mut cur: Vec<MafRow> = Vec::new();
    let flush = |cur: &mut Vec<MafRow>, out: &mut MafPairMap| -> Result<()> {
        let rx_row = cur
            .iter()
            .find(|r| horizon_aligner::matches_sample_prefix(&r.src, cx));
        let ry_row = cur
            .iter()
            .find(|r| horizon_aligner::matches_sample_prefix(&r.src, cy));
        if let (Some(rx), Some(ry)) = (rx_row, ry_row) {
            let ax = rx.seq.as_bytes();
            let ay = ry.seq.as_bytes();
            anyhow::ensure!(
                ax.len() == ay.len(),
                "MAF block alignment length mismatch ({} vs {})",
                ax.len(),
                ay.len()
            );

            let pair_set = out.entry((rx.src.clone(), ry.src.clone())).or_default();
            let (mut ix, mut iy) = (0usize, 0usize);
            for (&ca, &cb) in ax.iter().zip(ay) {
                if ca != b'-' && cb != b'-' {
                    let px = rx.fwd_pos(ix);
                    let py = ry.fwd_pos(iy);
                    pair_set.insert((px, py));
                }
                if ca != b'-' {
                    ix += 1;
                }
                if cb != b'-' {
                    iy += 1;
                }
            }
        }
        cur.clear();
        Ok(())
    };
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed == "a"
            || trimmed
                .strip_prefix('a')
                .is_some_and(|r| r.starts_with(char::is_whitespace))
        {
            if !cur.is_empty() {
                flush(&mut cur, &mut out)?;
            }
        } else if let Some((row, _)) = parse_maf_s_line(trimmed) {
            cur.push(row);
        }
        line.clear();
    }
    if !cur.is_empty() {
        flush(&mut cur, &mut out)?;
    }
    Ok(out)
}

struct EvalCliArgs {
    ref_name: String,
    target_name: String,
    ref_fa: std::path::PathBuf,
    tgt_fa: std::path::PathBuf,
    truth_maf: std::path::PathBuf,
    min_fallback_cov: Option<f64>,
    hsp_threshold: Option<i32>,
}

fn parse_eval_args() -> Result<Option<EvalCliArgs>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut pos = Vec::new();
    let mut min_fallback_cov = None;
    let mut hsp_threshold = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                println!("{USAGE}");
                return Ok(None);
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
            s if s.starts_with('-') => anyhow::bail!("unknown argument '{s}'"),
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
        min_fallback_cov,
        hsp_threshold,
    }))
}

fn main() -> Result<()> {
    let Some(cli) = parse_eval_args()? else {
        return Ok(());
    };
    let ref_contigs = read_named_contigs(&cli.ref_fa, &cli.ref_name)?;
    let tgt_contigs = read_named_contigs(&cli.tgt_fa, &cli.target_name)?;
    let maf_by_pair = parse_maf_pairs_by_contig(&cli.truth_maf, &cli.ref_name, &cli.target_name)?;

    let cfg = AlignConfig {
        min_fallback_cov: cli.min_fallback_cov,
        hsp_threshold: cli.hsp_threshold,
        ..AlignConfig::default()
    };

    let t0 = Instant::now();
    let needed_refs: HashSet<&str> = maf_by_pair
        .keys()
        .flat_map(|(r, _)| [r.as_str(), horizon_aligner::strip_sample_prefix(r, &cli.ref_name)])
        .collect();
    let needed_tgts: HashSet<&str> = maf_by_pair
        .keys()
        .flat_map(|(_, t)| [t.as_str(), horizon_aligner::strip_sample_prefix(t, &cli.target_name)])
        .collect();
    let ref_idx_cache: HashMap<String, _> = ref_contigs
        .par_iter()
        .filter(|(id, _)| needed_refs.contains(id.as_str()))
        .map(|(id, seq)| Ok((id.clone(), prepare_contig_index_opts(seq, false)?)))
        .collect::<Result<_>>()?;
    let tgt_idx_cache: HashMap<String, _> = tgt_contigs
        .par_iter()
        .filter(|(id, _)| needed_tgts.contains(id.as_str()))
        .map(|(id, seq)| Ok((id.clone(), prepare_contig_index_opts(seq, false)?)))
        .collect::<Result<_>>()?;

    let mut total_proto = 0usize;
    let mut total_maf = 0usize;
    let mut total_overlap = 0usize;
    let mut total_bubbles = 0usize;

    let mut pairs: Vec<MafEntry> = maf_by_pair.into_iter().collect();
    pairs.sort_by(|a, b| a.0.cmp(&b.0));

    let mut aligned_pairs = Vec::with_capacity(pairs.len());
    for ((ref_id, tgt_id), maf_set) in &pairs {
        let (Some(idx_r), Some(idx_t)) = (
            resolve_contig_index(&ref_idx_cache, ref_id, &cli.ref_name),
            resolve_contig_index(&tgt_idx_cache, tgt_id, &cli.target_name),
        ) else {
            eprintln!(
                "horizon_eval: warning: skipping MAF pair ({ref_id}, {tgt_id}) missing from FASTA inputs ({} pairs)",
                maf_set.len()
            );
            continue;
        };
        aligned_pairs.push((idx_r, idx_t, maf_set));
    }
    anyhow::ensure!(
        pairs.is_empty() || !aligned_pairs.is_empty(),
        "no truth MAF contig pairs matched the provided FASTA inputs"
    );

    let pair_stats: Vec<Result<(usize, usize, usize, usize)>> = aligned_pairs
        .into_par_iter()
        .map(|(idx_r, idx_t, maf_set)| {
            let mut aln = align_contig_pair_with_config(idx_r, idx_t, &cfg)?;
            aln.aligned_pairs.sort_unstable();
            aln.aligned_pairs.dedup();
            let overlap = aln
                .aligned_pairs
                .iter()
                .filter(|pair| maf_set.contains(pair))
                .count();
            Ok((
                maf_set.len(),
                aln.aligned_pairs.len(),
                overlap,
                aln.macro_bubbles + aln.internal_bubbles,
            ))
        })
        .collect();

    for res in pair_stats {
        let (maf_len, proto_len, overlap, bubbles) = res?;
        total_maf += maf_len;
        total_proto += proto_len;
        total_overlap += overlap;
        total_bubbles += bubbles;
    }

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
            "a score=10\ns hg38.chr22 10 3 + 100 ACG\ns dog.chr1 5 3 - 50 TGC\n\n",
        )
        .unwrap();
        let maf = parse_maf_pairs_by_contig(&maf_path, "hg38", "dog").unwrap();
        let pairs = maf
            .get(&("hg38.chr22".to_string(), "dog.chr1".to_string()))
            .unwrap();
        assert!(pairs.contains(&(10, 44)));
        assert!(pairs.contains(&(11, 43)));
        assert!(pairs.contains(&(12, 42)));

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
}
