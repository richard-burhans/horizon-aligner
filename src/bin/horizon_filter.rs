//! Streaming MDL + BOCPD + SDUST filter for PAF alignments (e.g. `FastGA` / `wfmash` outputs).

use anyhow::{Context, Result, bail};
use horizon_aligner::bocpd::{
    GAP_EXTEND_BITS, GAP_OPEN_BITS, IDENTITY_PRIOR, base_probs, joint_pair_bits,
};
use horizon_aligner::formats::read_fasta_contigs;
use horizon_aligner::io_util::open_text_writer;
use horizon_aligner::mdl::{HmmParams, calibrate_null_crc_from_seqs};
use horizon_aligner::paf::{CigarOp, PafRecord, stream_paf_records};
use horizon_aligner::seeding::{complement_base, sdust_mask};
use std::collections::HashMap;
use std::env;
use std::ops::ControlFlow;
use std::time::Instant;

const KARLIN_MDL_COEFF: f64 = 0.55;

struct FilterConfig {
    paf_path: String,
    query_fa: Option<String>,
    target_fa: Option<String>,
    out_path: Option<String>,
    min_bits_override: Option<f64>,
    min_nonrepeat_bp: usize,
    min_aln_len: usize,
    keep_all: bool,
    max_occ: Option<usize>,
}

struct IndexedSeq {
    seq: Vec<u8>,
    masked: Vec<bool>,
    gc: f64,
}

impl IndexedSeq {
    fn from_dna(seq: Vec<u8>) -> Self {
        let masked = sdust_mask(&seq);
        let mut gc_cnt = 0usize;
        let mut acgt = 0usize;
        for &b in &seq {
            match b.to_ascii_uppercase() {
                b'G' | b'C' => {
                    gc_cnt += 1;
                    acgt += 1;
                }
                b'A' | b'T' => acgt += 1,
                _ => {}
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let gc = if acgt > 0 {
            (gc_cnt as f64 / acgt as f64).clamp(0.05, 0.95)
        } else {
            0.5
        };
        Self { seq, masked, gc }
    }
}

struct ScoreResult {
    mdl_bits: f64,
    nonrepeat_bp: usize,
    threshold_bits: f64,
}

fn load_fasta_map(path: &str) -> Result<HashMap<String, IndexedSeq>> {
    let contigs = read_fasta_contigs(std::path::Path::new(path), 0)
        .with_context(|| format!("reading FASTA '{path}'"))?;
    let mut map = HashMap::with_capacity(contigs.len());
    for (name, seq) in contigs {
        map.insert(name, IndexedSeq::from_dna(seq));
    }
    Ok(map)
}

fn lookup_contig<'a>(map: &'a HashMap<String, IndexedSeq>, name: &str) -> Option<&'a IndexedSeq> {
    if let Some(entry) = map.get(name) {
        return Some(entry);
    }
    // Also match bare contig name before `:start-end` range suffixes (common in subrange FASTAs).
    if let Some((prefix, _)) = name.split_once(':')
        && let Some(entry) = map.get(prefix)
    {
        return Some(entry);
    }
    None
}

fn score_from_counts_only(rec: &PafRecord, threshold_bits: f64) -> ScoreResult {
    let mut eq_bp = 0usize;
    let mut diff_bp = 0usize;
    let mut gap_opens = 0usize;
    let mut gap_ext_bp = 0usize;

    if let Some(cigar_iter) = rec.cigar_ops() {
        let mut valid = true;
        let mut has_match_op = false;
        for op_res in cigar_iter {
            let Ok(op) = op_res else {
                valid = false;
                break;
            };
            match op {
                CigarOp::Equal(n) => eq_bp += n,
                CigarOp::Match(_) => has_match_op = true,
                CigarOp::Diff(n) => diff_bp += n,
                CigarOp::Ins(n) | CigarOp::Del(n) if n > 0 => {
                    gap_opens += 1;
                    gap_ext_bp += n.saturating_sub(1);
                }
                _ => {}
            }
        }
        if !valid || has_match_op || (eq_bp == 0 && diff_bp == 0) {
            eq_bp = rec.matches;
            diff_bp = rec.aln_len.saturating_sub(rec.matches);
        }
    } else {
        eq_bp = rec.matches;
        diff_bp = rec.aln_len.saturating_sub(rec.matches);
    }

    let pi = base_probs(0.5);
    let match_gain = joint_pair_bits(b'A', b'A', pi, pi, IDENTITY_PRIOR);
    let mismatch_cost = joint_pair_bits(b'A', b'C', pi, pi, IDENTITY_PRIOR);
    #[allow(clippy::cast_precision_loss)]
    let mdl_bits = (eq_bp as f64) * match_gain
        + (diff_bp as f64) * mismatch_cost
        + (gap_opens as f64) * GAP_OPEN_BITS
        + (gap_ext_bp as f64) * GAP_EXTEND_BITS;
    ScoreResult {
        mdl_bits,
        nonrepeat_bp: eq_bp + diff_bp,
        threshold_bits,
    }
}

#[allow(clippy::too_many_arguments)]
fn score_match_run(
    n: usize,
    is_rc: bool,
    qs: &IndexedSeq,
    ts: &IndexedSeq,
    pi_q: [f64; 4],
    pi_t: [f64; 4],
    q_pos: &mut usize,
    t_pos: &mut usize,
    mdl_bits: &mut f64,
    nonrepeat_bp: &mut usize,
) -> bool {
    for _ in 0..n {
        let qi = if is_rc {
            if *q_pos == 0 {
                return false;
            }
            *q_pos -= 1;
            *q_pos
        } else {
            if *q_pos >= qs.seq.len() {
                return false;
            }
            let cur = *q_pos;
            *q_pos += 1;
            cur
        };
        if *t_pos >= ts.seq.len() {
            return false;
        }
        let ti = *t_pos;
        *t_pos += 1;

        // Down-weight SDUST-masked simple/tandem repeats so low-complexity false
        // positives cannot accumulate positive MDL support.
        if qs.masked[qi] && ts.masked[ti] {
            continue;
        }
        *nonrepeat_bp += 1;
        let qb = if is_rc {
            complement_base(qs.seq[qi])
        } else {
            qs.seq[qi]
        };
        *mdl_bits += joint_pair_bits(qb, ts.seq[ti], pi_q, pi_t, IDENTITY_PRIOR);
    }
    true
}

#[allow(clippy::too_many_lines)]
fn score_paf_record(
    rec: &PafRecord,
    q_map: Option<&HashMap<String, IndexedSeq>>,
    t_map: Option<&HashMap<String, IndexedSeq>>,
    min_bits_override: Option<f64>,
    crc_cache: &mut HashMap<String, HashMap<String, f64>>,
) -> ScoreResult {
    #[allow(clippy::cast_precision_loss)]
    let search_space_bits = ((rec.qlen.max(1) as f64) * (rec.tlen.max(1) as f64))
        .log2()
        .max(10.0);
    let base_floor = (KARLIN_MDL_COEFF * search_space_bits).max(10.0);
    let mut threshold_bits = min_bits_override.unwrap_or(base_floor);

    let (Some(qm), Some(tm)) = (q_map, t_map) else {
        return score_from_counts_only(rec, threshold_bits);
    };
    let (Some(qs), Some(ts)) = (lookup_contig(qm, &rec.qname), lookup_contig(tm, &rec.tname))
    else {
        return score_from_counts_only(rec, threshold_bits);
    };
    if min_bits_override.is_none() {
        let crc_floor = if let Some(&cached) = crc_cache
            .get(rec.qname.as_str())
            .and_then(|m| m.get(rec.tname.as_str()))
        {
            cached
        } else {
            let val = calibrate_null_crc_from_seqs(
                &qs.seq,
                &ts.seq,
                15,
                0.07,
                base_floor,
                &HmmParams::default(),
            );
            crc_cache
                .entry(rec.qname.clone())
                .or_default()
                .insert(rec.tname.clone(), val);
            val
        };
        threshold_bits = threshold_bits.max(crc_floor);
    }
    let reject_corrupt_cigar = ScoreResult {
        mdl_bits: 0.0,
        nonrepeat_bp: 0,
        threshold_bits,
    };
    let Some(cigar_iter) = rec.cigar_ops() else {
        return score_from_counts_only(rec, threshold_bits);
    };
    if rec.qend > qs.seq.len()
        || rec.tend > ts.seq.len()
        || rec.qstart > rec.qend
        || rec.tstart > rec.tend
    {
        return reject_corrupt_cigar;
    }

    let pi_q = base_probs(qs.gc);
    let pi_t = base_probs(ts.gc);
    let is_rc = rec.strand == '-';
    let mut q_pos = if is_rc { rec.qend } else { rec.qstart };
    let mut t_pos = rec.tstart;
    let mut mdl_bits = 0.0_f64;
    let mut nonrepeat_bp = 0usize;

    for op_res in cigar_iter {
        let Ok(op) = op_res else {
            return reject_corrupt_cigar;
        };
        match op {
            CigarOp::Match(n) | CigarOp::Equal(n) | CigarOp::Diff(n) => {
                if !score_match_run(
                    n,
                    is_rc,
                    qs,
                    ts,
                    pi_q,
                    pi_t,
                    &mut q_pos,
                    &mut t_pos,
                    &mut mdl_bits,
                    &mut nonrepeat_bp,
                ) {
                    return reject_corrupt_cigar;
                }
            }
            CigarOp::Ins(n) | CigarOp::SoftClip(n) => {
                if is_rc {
                    let Some(next_q) = q_pos.checked_sub(n) else {
                        return reject_corrupt_cigar;
                    };
                    q_pos = next_q;
                } else {
                    let Some(next_q) = q_pos.checked_add(n).filter(|&p| p <= qs.seq.len()) else {
                        return reject_corrupt_cigar;
                    };
                    q_pos = next_q;
                }
                if matches!(op, CigarOp::Ins(_)) && n > 0 {
                    #[allow(clippy::cast_precision_loss)]
                    {
                        mdl_bits += GAP_OPEN_BITS + (n.saturating_sub(1) as f64) * GAP_EXTEND_BITS;
                    }
                }
            }
            CigarOp::Del(n) | CigarOp::RefSkip(n) => {
                let Some(next_t) = t_pos.checked_add(n).filter(|&p| p <= ts.seq.len()) else {
                    return reject_corrupt_cigar;
                };
                t_pos = next_t;
                if matches!(op, CigarOp::Del(_)) && n > 0 {
                    #[allow(clippy::cast_precision_loss)]
                    {
                        mdl_bits += GAP_OPEN_BITS + (n.saturating_sub(1) as f64) * GAP_EXTEND_BITS;
                    }
                }
            }
            CigarOp::HardClip(_) | CigarOp::Pad(_) => {}
        }
    }

    ScoreResult {
        mdl_bits,
        nonrepeat_bp,
        threshold_bits,
    }
}

fn net_top_k_records(
    mut scored: Vec<(PafRecord, ScoreResult)>,
    max_occ: usize,
) -> Vec<(PafRecord, ScoreResult)> {
    scored.sort_unstable_by(|a, b| b.1.mdl_bits.total_cmp(&a.1.mdl_bits));
    let mut occ_q: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    let mut occ_t: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    let mut kept = Vec::new();

    let overlap_count = |ivs: &[(usize, usize)], s: usize, e: usize| -> usize {
        let len = e.saturating_sub(s);
        if len == 0 {
            return usize::MAX;
        }
        // `ivs` is sorted by interval end (`oe`), so any interval with `oe <= s` cannot overlap.
        let start_idx = ivs.partition_point(|&(_, oe)| oe <= s);
        ivs[start_idx..]
            .iter()
            .filter(|&&(os, oe)| {
                let ov = os.max(s);
                let oe_min = oe.min(e);
                oe_min > ov && (oe_min - ov) * 2 >= len
            })
            .count()
    };

    for (rec, res) in scored {
        let q_ivs = occ_q.get(&rec.qname).map_or(&[][..], Vec::as_slice);
        let t_ivs = occ_t.get(&rec.tname).map_or(&[][..], Vec::as_slice);
        if overlap_count(q_ivs, rec.qstart, rec.qend) >= max_occ
            || overlap_count(t_ivs, rec.tstart, rec.tend) >= max_occ
        {
            continue;
        }
        let q_vec = occ_q.entry(rec.qname.clone()).or_default();
        let q_pos = q_vec.partition_point(|&(_, oe)| oe < rec.qend);
        q_vec.insert(q_pos, (rec.qstart, rec.qend));
        let t_vec = occ_t.entry(rec.tname.clone()).or_default();
        let t_pos = t_vec.partition_point(|&(_, oe)| oe < rec.tend);
        t_vec.insert(t_pos, (rec.tstart, rec.tend));
        kept.push((rec, res));
    }
    kept.sort_unstable_by(|a, b| {
        a.0.qname
            .cmp(&b.0.qname)
            .then(a.0.qstart.cmp(&b.0.qstart))
            .then(a.0.tstart.cmp(&b.0.tstart))
    });
    kept
}

fn run_filter(cfg: &FilterConfig) -> Result<()> {
    let t0 = Instant::now();
    let q_map = cfg.query_fa.as_deref().map(load_fasta_map).transpose()?;
    let same_fa = cfg.target_fa.is_some() && cfg.target_fa.as_deref() == cfg.query_fa.as_deref();
    let t_map = match cfg.target_fa.as_deref() {
        Some(_) if same_fa => None,
        Some(path) => Some(load_fasta_map(path)?),
        None => None,
    };
    let t_ref = if same_fa {
        q_map.as_ref()
    } else {
        t_map.as_ref()
    };

    let mut out = open_text_writer(cfg.out_path.as_deref())?;
    let mut total = 0usize;
    let mut kept = 0usize;
    let mut filtered_mdl = 0usize;
    let mut filtered_len = 0usize;
    let mut net_buf: Vec<(PafRecord, ScoreResult)> = Vec::new();
    let mut crc_cache: HashMap<String, HashMap<String, f64>> = HashMap::new();

    stream_paf_records(&cfg.paf_path, |rec| {
        total += 1;
        let res = score_paf_record(
            &rec,
            q_map.as_ref(),
            t_ref,
            cfg.min_bits_override,
            &mut crc_cache,
        );
        let pass_len = rec.aln_len >= cfg.min_aln_len && res.nonrepeat_bp >= cfg.min_nonrepeat_bp;
        let pass_mdl = res.mdl_bits >= res.threshold_bits;
        let pass = pass_len && pass_mdl;

        if !pass {
            if pass_len {
                filtered_mdl += 1;
            } else {
                filtered_len += 1;
            }
        }

        if cfg.max_occ.is_some() {
            if pass {
                net_buf.push((rec, res));
            }
            return Ok(ControlFlow::Continue(()));
        }

        if pass || cfg.keep_all {
            if pass {
                kept += 1;
            }
            let hb = u8::from(!pass);
            rec.write_line_with_mdl(&mut out, res.mdl_bits, hb)?;
        }
        Ok(ControlFlow::Continue(()))
    })?;

    if let Some(k) = cfg.max_occ {
        let netted = net_top_k_records(net_buf, k);
        kept = netted.len();
        for (rec, res) in netted {
            rec.write_line_with_mdl(&mut out, res.mdl_bits, 0)?;
        }
    }

    out.finish()?;
    eprintln!(
        "[horizon_filter] total={} kept={} filtered_mdl={} filtered_len={} elapsed={:.2}s",
        total,
        kept,
        filtered_mdl,
        filtered_len,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

fn next_arg_val<'a>(args: &'a [String], i: usize, flag: &str) -> Result<&'a str> {
    args.get(i + 1)
        .map(String::as_str)
        .with_context(|| format!("missing value for {flag}"))
}

fn parse_args() -> Result<Option<FilterConfig>> {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut paf_path: Option<String> = None;
    let mut query_fa: Option<String> = None;
    let mut target_fa: Option<String> = None;
    let mut out_path: Option<String> = None;
    let mut min_bits_override: Option<f64> = None;
    let mut min_nonrepeat_bp = 0usize;
    let mut min_aln_len = 0usize;
    let mut keep_all = false;
    let mut max_occ: Option<usize> = None;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                println!(
                    "usage: horizon_filter --paf <in.paf[.gz]|-> [--query <q.fa>] [--target <t.fa>] [--seqs <all.fa>] [--out <out.paf>] [--min-bits <f64>] [--min-nonrepeat-bp <usize>] [--min-len <usize>] [--net] [--max-occ <usize>] [--keep-all]"
                );
                return Ok(None);
            }
            "--paf" | "-i" => {
                paf_path = Some(next_arg_val(&args, i, "--paf")?.to_owned());
                i += 2;
            }
            "--query" | "-q" => {
                query_fa = Some(next_arg_val(&args, i, "--query")?.to_owned());
                i += 2;
            }
            "--target" | "-t" => {
                target_fa = Some(next_arg_val(&args, i, "--target")?.to_owned());
                i += 2;
            }
            "--seqs" | "-s" => {
                let s = next_arg_val(&args, i, "--seqs")?.to_owned();
                query_fa = Some(s.clone());
                target_fa = Some(s);
                i += 2;
            }
            "--out" | "-o" => {
                out_path = Some(next_arg_val(&args, i, "--out")?.to_owned());
                i += 2;
            }
            "--min-bits" => {
                let val: f64 = next_arg_val(&args, i, "--min-bits")?
                    .parse()
                    .context("parsing --min-bits")?;
                anyhow::ensure!(val.is_finite(), "--min-bits must be finite, got {val}");
                min_bits_override = Some(val);
                i += 2;
            }
            "--min-nonrepeat-bp" => {
                min_nonrepeat_bp = next_arg_val(&args, i, "--min-nonrepeat-bp")?
                    .parse()
                    .context("parsing --min-nonrepeat-bp")?;
                i += 2;
            }
            "--min-len" => {
                min_aln_len = next_arg_val(&args, i, "--min-len")?
                    .parse()
                    .context("parsing --min-len")?;
                i += 2;
            }
            "--keep-all" => {
                keep_all = true;
                i += 1;
            }
            "--net" => {
                max_occ = Some(1);
                i += 1;
            }
            "--max-occ" => {
                let k: usize = next_arg_val(&args, i, "--max-occ")?
                    .parse()
                    .context("parsing --max-occ")?;
                anyhow::ensure!(k >= 1, "--max-occ must be >= 1");
                max_occ = Some(k);
                i += 2;
            }
            other => bail!("unknown argument '{other}'"),
        }
    }

    anyhow::ensure!(
        !(keep_all && max_occ.is_some()),
        "--keep-all cannot be combined with --net or --max-occ"
    );
    let paf_path = paf_path.context("missing required --paf <in.paf[.gz]|->")?;
    Ok(Some(FilterConfig {
        paf_path,
        query_fa,
        target_fa,
        out_path,
        min_bits_override,
        min_nonrepeat_bp,
        min_aln_len,
        keep_all,
        max_occ,
    }))
}

fn main() -> Result<()> {
    let Some(cfg) = parse_args()? else {
        return Ok(());
    };
    run_filter(&cfg)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mdl_filter_rejects_simple_repeat_and_accepts_complex_match() {
        // 80bp complex random-like homologous region vs 80bp ATATAT... dinucleotide repeat
        let complex =
            b"ACGTAGTGCTAGCATCGATCGTAGCTAGCTGATCGATCGATCGATCGTAGCTAGCTAGCATCGATCGTAGCTAGCTAGCA";
        let repeat =
            b"ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT";
        let mut q_map = HashMap::new();
        q_map.insert("q_c".to_owned(), IndexedSeq::from_dna(complex.to_vec()));
        q_map.insert("q_r".to_owned(), IndexedSeq::from_dna(repeat.to_vec()));

        let mut t_map = HashMap::new();
        t_map.insert("t_c".to_owned(), IndexedSeq::from_dna(complex.to_vec()));
        t_map.insert("t_r".to_owned(), IndexedSeq::from_dna(repeat.to_vec()));

        let paf_c =
            PafRecord::parse_line("q_c\t80\t0\t80\t+\tt_c\t80\t0\t80\t80\t80\t60\tcg:Z:80=")
                .expect("paf");
        let paf_r =
            PafRecord::parse_line("q_r\t80\t0\t80\t+\tt_r\t80\t0\t80\t80\t80\t60\tcg:Z:80=")
                .expect("paf");

        let mut crc_cache = HashMap::new();
        let res_c = score_paf_record(&paf_c, Some(&q_map), Some(&t_map), None, &mut crc_cache);
        let res_r = score_paf_record(&paf_r, Some(&q_map), Some(&t_map), None, &mut crc_cache);
        assert!(res_c.mdl_bits > res_c.threshold_bits);
        assert_eq!(res_c.nonrepeat_bp, 80);
        assert!(res_r.mdl_bits <= 0.0, "repeat should be SDUST-masked");
        assert_eq!(res_r.nonrepeat_bp, 0);

        // Empirical biological CRC calibration over >=128 bp sequences stays >= base_floor.
        let bio_x = complex.repeat(4);
        let bio_y = bio_x.iter().rev().copied().collect::<Vec<u8>>();
        let crc_floor =
            calibrate_null_crc_from_seqs(&bio_x, &bio_y, 15, 0.07, 10.0, &HmmParams::default());
        assert!(crc_floor >= 10.0);
    }

    #[test]
    fn counts_only_does_not_treat_match_op_as_100pct_exact() {
        let paf_m =
            PafRecord::parse_line("q\t100\t0\t100\t+\tt\t100\t0\t100\t50\t100\t60\tcg:Z:100M")
                .expect("paf");
        let paf_eq =
            PafRecord::parse_line("q\t100\t0\t100\t+\tt\t100\t0\t100\t100\t100\t60\tcg:Z:100=")
                .expect("paf");
        let res_m = score_from_counts_only(&paf_m, 10.0);
        let res_eq = score_from_counts_only(&paf_eq, 10.0);
        assert!(
            res_m.mdl_bits < res_eq.mdl_bits,
            "100M with 50/100 matches must score strictly lower than 100="
        );
    }

    #[test]
    fn corrupt_or_oob_cigar_with_fasta_is_rejected_not_counted() {
        let complex =
            b"ACGTAGTGCTAGCATCGATCGTAGCTAGCTGATCGATCGATCGATCGTAGCTAGCTAGCATCGATCGTAGCTAGCTAGCA";
        let mut q_map = HashMap::new();
        q_map.insert("q".to_owned(), IndexedSeq::from_dna(complex.to_vec()));
        let mut t_map = HashMap::new();
        t_map.insert("t".to_owned(), IndexedSeq::from_dna(complex.to_vec()));

        let paf_oob = PafRecord::parse_line("q\t80\t0\t80\t+\tt\t80\t0\t80\t80\t80\t60\tcg:Z:200=")
            .expect("paf");
        let mut crc_cache = HashMap::new();
        let res = score_paf_record(
            &paf_oob,
            Some(&q_map),
            Some(&t_map),
            Some(10.0),
            &mut crc_cache,
        );
        assert_eq!(res.mdl_bits, 0.0);
        assert_eq!(res.nonrepeat_bp, 0);
    }
}
