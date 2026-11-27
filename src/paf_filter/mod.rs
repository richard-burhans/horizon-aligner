//! Streaming SDUST + MDL + operator-matched biological CRC scoring and top-k interval netting
//! for PAF alignments (used by `horizon_filter`).

use crate::bocpd::{GAP_EXTEND_BITS, GAP_OPEN_BITS, IDENTITY_PRIOR, base_probs, joint_pair_bits};
use crate::formats::read_fasta_contigs;
use crate::io_util::open_text_writer;
use crate::mdl::{HmmParams, calibrate_crc_threshold, calibrate_null_crc_for_gc};
use crate::paf::{CigarOp, PafRecord, stream_paf_records};
use crate::seeding::{complement_base, sdust_mask};
use crate::splitmix64;
use anyhow::{Context, Result};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::ops::ControlFlow;

mod net;
use net::net_top_k_partition_with_mode;
pub use net::{net_top_k_records, net_top_k_records_with_mode};

const KARLIN_MDL_COEFF: f64 = 0.55;

/// Ranking objective used when netting overlapping PAF records (`--net` / `--max-occ`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum NetScoreMode {
    /// Rank by cumulative MDL bit-savings (linear in aligned length; default).
    #[default]
    Mdl,
    /// Identity-weighted sublinear-length ranking (`ani^2 * ln(length)`) with symmetric two-axis reciprocal-best netting,
    /// favoring high-identity orthologs over longer low-identity paralogs/palindromes while
    /// retaining MDL/CRC pre-filtering.
    LogLengthAni,
}

/// Configuration for [`run_paf_filter`].
#[derive(Debug, Clone)]
pub struct PafFilterConfig {
    /// Input PAF path (plain, `.gz`, `.zst`, or `"-"` for stdin).
    pub paf_path: String,
    /// Optional query FASTA path for sequence-aware SDUST+MDL scoring.
    pub query_fa: Option<String>,
    /// Optional target FASTA path for sequence-aware SDUST+MDL scoring.
    pub target_fa: Option<String>,
    /// Optional output PAF path (`None` or `Some("-")` writes to stdout).
    pub out_path: Option<String>,
    /// Optional explicit MDL bit-savings floor overriding Karlin/CRC calibration.
    pub min_bits_override: Option<f64>,
    /// Conformal Risk Control alpha level for reversed biological decoy calibration.
    pub crc_alpha: f64,
    /// Number of reversed biological decoy windows (`0` disables CRC).
    pub crc_n_cal: usize,
    /// Minimum unmasked non-repeat aligned base pairs required to keep a record.
    pub min_nonrepeat_bp: usize,
    /// Minimum PAF alignment block length (`aln_len`) required to keep a record.
    pub min_aln_len: usize,
    /// Emit failing records tagged with `hb:i:1` instead of dropping them.
    pub keep_all: bool,
    /// Optional maximum overlapping multiplicity (`Some(1)` = strict 1-to-1 netting).
    pub max_occ: Option<usize>,
    /// Ranking objective applied during `--net` / `--max-occ` interval netting.
    pub net_score: NetScoreMode,
}

/// Summary statistics returned by [`run_paf_filter`].
#[derive(Debug, Clone, Copy, Default)]
pub struct PafFilterStats {
    /// Total PAF records scanned.
    pub total: usize,
    /// Records retained after MDL/length filtering and optional netting.
    pub kept: usize,
    /// Records dropped by the MDL bit-savings floor.
    pub filtered_mdl: usize,
    /// Records dropped by `--min-len` or `--min-nonrepeat-bp`.
    pub filtered_len: usize,
}

/// Sequence plus precomputed SDUST repeat mask and global GC content for PAF scoring.
#[derive(Debug, Clone)]
pub struct IndexedSeq {
    /// Raw ASCII DNA bytes.
    pub seq: Vec<u8>,
    /// Per-base SDUST low-complexity mask.
    pub masked: Vec<bool>,
    /// Global GC fraction clamped to `[0.05, 0.95]`.
    pub gc: f64,
}

impl IndexedSeq {
    /// Index a DNA sequence with SDUST masking and global GC estimation.
    #[must_use]
    pub fn from_dna(seq: Vec<u8>) -> Self {
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

/// Scoring output for a single PAF record.
#[derive(Debug, Clone, Copy)]
pub struct ScoreResult {
    /// Total MDL bit-savings under the PAF scoring operator.
    pub mdl_bits: f64,
    /// Number of unmasked non-repeat aligned columns scored.
    pub nonrepeat_bp: usize,
    /// Effective bit-savings acceptance threshold applied to this record.
    pub threshold_bits: f64,
}

/// Load a FASTA file and index all contigs in parallel with Rayon.
///
/// # Errors
/// Returns an error if the FASTA file cannot be opened or read.
pub fn load_fasta_map(path: &str) -> Result<HashMap<String, IndexedSeq>> {
    let contigs = read_fasta_contigs(std::path::Path::new(path), 0)
        .with_context(|| format!("reading FASTA '{path}'"))?;
    let map: HashMap<String, IndexedSeq> = contigs
        .into_par_iter()
        .map(|(name, seq)| (name, IndexedSeq::from_dna(seq)))
        .collect();
    Ok(map)
}

/// Resolve a contig name (including `contig:start-end` subrange names) against `map`, returning
/// `(indexed_contig, start_offset)`.
#[must_use]
pub fn lookup_contig<'a, S: std::hash::BuildHasher>(
    map: &'a HashMap<String, IndexedSeq, S>,
    name: &str,
) -> Option<(&'a IndexedSeq, usize)> {
    if let Some(entry) = map.get(name) {
        return Some((entry, 0));
    }
    if let Some((prefix, range)) = name.rsplit_once(':')
        && let Some((s_str, e_str)) = range.split_once('-')
        && let (Ok(start), Ok(end)) = (s_str.parse::<usize>(), e_str.parse::<usize>())
        && start <= end
        && let Some(entry) = map.get(prefix)
        && end <= entry.seq.len()
    {
        return Some((entry, start));
    }
    None
}

/// Score a PAF record using only its match/length/CIGAR counts (sequence-free mode).
#[must_use]
pub fn score_from_counts_only(rec: &PafRecord, threshold_bits: f64) -> ScoreResult {
    let mut eq_bp = 0usize;
    let mut diff_bp = 0usize;
    let mut gap_opens = 0usize;
    let mut gap_ext_bp = 0usize;
    let mut gap_bp = 0usize;
    let mut ref_skip_bp = 0usize;

    if let Some(cigar_iter) = rec.cigar_ops() {
        let mut valid = true;
        let mut has_match_op = false;
        for op_res in cigar_iter {
            let Ok(op) = op_res else {
                valid = false;
                break;
            };
            match op {
                CigarOp::Equal(n) => {
                    let Some(next) = eq_bp.checked_add(n) else {
                        valid = false;
                        break;
                    };
                    eq_bp = next;
                }
                CigarOp::Match(_) => has_match_op = true,
                CigarOp::Diff(n) => {
                    let Some(next) = diff_bp.checked_add(n) else {
                        valid = false;
                        break;
                    };
                    diff_bp = next;
                }
                CigarOp::Ins(n) | CigarOp::Del(n) if n > 0 => {
                    let (Some(no), Some(ne), Some(ng)) = (
                        gap_opens.checked_add(1),
                        gap_ext_bp.checked_add(n.saturating_sub(1)),
                        gap_bp.checked_add(n),
                    ) else {
                        valid = false;
                        break;
                    };
                    gap_opens = no;
                    gap_ext_bp = ne;
                    gap_bp = ng;
                }
                CigarOp::RefSkip(n) if n > 0 => {
                    let Some(ns) = ref_skip_bp.checked_add(n) else {
                        valid = false;
                        break;
                    };
                    ref_skip_bp = ns;
                }
                _ => {}
            }
        }
        if !valid {
            gap_opens = 0;
            gap_ext_bp = 0;
            eq_bp = rec.matches;
            diff_bp = rec.aln_len.saturating_sub(rec.matches);
        } else if has_match_op || (eq_bp == 0 && diff_bp == 0) {
            eq_bp = rec.matches;
            diff_bp = rec
                .aln_len
                .saturating_sub(rec.matches)
                .saturating_sub(gap_bp)
                .saturating_sub(ref_skip_bp);
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
        nonrepeat_bp: eq_bp.saturating_add(diff_bp),
        threshold_bits,
    }
}

#[allow(clippy::too_many_arguments)]
fn score_match_run(
    n: usize,
    is_rc: bool,
    qs: &IndexedSeq,
    ts: &IndexedSeq,
    q_lo: usize,
    q_hi: usize,
    t_hi: usize,
    pi_q: [f64; 4],
    pi_t: [f64; 4],
    q_pos: &mut usize,
    t_pos: &mut usize,
    mdl_bits: &mut f64,
    nonrepeat_bp: &mut usize,
) -> bool {
    for _ in 0..n {
        let qi = if is_rc {
            if *q_pos <= q_lo {
                return false;
            }
            *q_pos -= 1;
            *q_pos
        } else {
            if *q_pos >= q_hi {
                return false;
            }
            let cur = *q_pos;
            *q_pos += 1;
            cur
        };
        if *t_pos >= t_hi {
            return false;
        }
        let ti = *t_pos;
        *t_pos += 1;

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

#[allow(clippy::too_many_arguments)]
fn paf_operator_local_sw(
    qx: &[u8],
    q_mask: &[bool],
    ty_rev: &[u8],
    t_mask_rev: &[bool],
    pi_q: [f64; 4],
    pi_t: [f64; 4],
    band_w: usize,
    prev_h: &mut [f64],
    prev_e: &mut [f64],
    cur_h: &mut [f64],
    cur_e: &mut [f64],
) -> f64 {
    let n = qx.len();
    let m = ty_rev.len();
    prev_h[..=m].fill(0.0);
    prev_e[..=m].fill(f64::NEG_INFINITY);
    let mut best = 0.0_f64;

    for i in 1..=n {
        let j_lo = i.saturating_sub(band_w).max(1);
        let j_hi = (i + band_w).min(m);
        cur_h[j_lo - 1] = 0.0;
        cur_e[j_lo - 1] = f64::NEG_INFINITY;
        let mut f_gap = f64::NEG_INFINITY;
        let qb = qx[i - 1];
        let qm = q_mask[i - 1];
        for j in j_lo..=j_hi {
            let sub = if qm && t_mask_rev[j - 1] {
                0.0
            } else {
                joint_pair_bits(qb, ty_rev[j - 1], pi_q, pi_t, IDENTITY_PRIOR)
            };
            let diag = prev_h[j - 1] + sub;
            let e_val = (prev_h[j] + GAP_OPEN_BITS).max(prev_e[j] + GAP_EXTEND_BITS);
            f_gap = (cur_h[j - 1] + GAP_OPEN_BITS).max(f_gap + GAP_EXTEND_BITS);
            let h = diag.max(e_val).max(f_gap).max(0.0);
            cur_h[j] = h;
            cur_e[j] = e_val;
            if h > best {
                best = h;
            }
        }
        prev_h[j_lo..=j_hi].copy_from_slice(&cur_h[j_lo..=j_hi]);
        prev_e[j_lo..=j_hi].copy_from_slice(&cur_e[j_lo..=j_hi]);
    }
    best
}

/// Calibrate the PAF scoring operator's CRC bit-savings floor against multi-scale reversed
/// biological decoy windows (`64, 128, 256, 512` bp).
#[must_use]
pub fn calibrate_paf_operator_crc(
    qs: &IndexedSeq,
    ts: &IndexedSeq,
    n_cal: usize,
    alpha: f64,
    base_floor: f64,
) -> f64 {
    const MIN_WIN: usize = 64;
    const MAX_WIN: usize = 512;
    const SCALES: [usize; 4] = [64, 128, 256, 512];
    const BAND_W: usize = 32;

    if n_cal == 0 || !alpha.is_finite() || alpha <= 0.0 {
        return base_floor;
    }

    if qs.seq.len() < MIN_WIN * 2 || ts.seq.len() < MIN_WIN * 2 {
        let avg_gc = f64::midpoint(qs.gc, ts.gc);
        return calibrate_null_crc_for_gc(avg_gc, n_cal, alpha, base_floor, &HmmParams::default());
    }

    let pi_q = base_probs(qs.gc);
    let pi_t = base_probs(ts.gc);
    let max_avail = (qs.seq.len().min(ts.seq.len()) / 2).clamp(MIN_WIN, MAX_WIN);
    let seed_mix = (qs.seq.len() as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
        ^ (ts.seq.len() as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9);

    let mut scores = Vec::with_capacity(n_cal);
    let mut t_rev = Vec::with_capacity(MAX_WIN);
    let mut t_mask_rev = Vec::with_capacity(MAX_WIN);
    let mut prev_h = vec![0.0_f64; MAX_WIN + 2];
    let mut prev_e = vec![f64::NEG_INFINITY; MAX_WIN + 2];
    let mut cur_h = vec![0.0_f64; MAX_WIN + 2];
    let mut cur_e = vec![f64::NEG_INFINITY; MAX_WIN + 2];

    let max_tries = n_cal.saturating_mul(6).max(32);
    for attempt in 0..max_tries {
        if scores.len() >= n_cal {
            break;
        }
        let win = SCALES[scores.len() % SCALES.len()].min(max_avail);
        let max_q = qs.seq.len() - win;
        let max_t = ts.seq.len() - win;
        let hq = splitmix64(seed_mix ^ ((attempt as u64).wrapping_mul(2) + 1));
        let ht = splitmix64(seed_mix ^ ((attempt as u64).wrapping_mul(2) + 2));
        #[allow(clippy::cast_possible_truncation)]
        let sq = (hq as usize) % (max_q + 1);
        #[allow(clippy::cast_possible_truncation)]
        let mut st = (ht as usize) % (max_t + 1);
        if qs.seq.len() == ts.seq.len() && sq.abs_diff(st) < win {
            st = if sq >= win {
                sq - win
            } else {
                (sq + win).min(max_t)
            };
        }
        let wq = &qs.seq[sq..sq + win];
        let wt = &ts.seq[st..st + win];
        let acgt_q = wq
            .iter()
            .filter(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .count();
        let acgt_t = wt
            .iter()
            .filter(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .count();
        if acgt_q * 2 < win || acgt_t * 2 < win {
            continue;
        }
        t_rev.clear();
        t_rev.extend(wt.iter().rev().copied());
        t_mask_rev.clear();
        t_mask_rev.extend(ts.masked[st..st + win].iter().rev().copied());

        let sw_bits = paf_operator_local_sw(
            wq,
            &qs.masked[sq..sq + win],
            &t_rev,
            &t_mask_rev,
            pi_q,
            pi_t,
            BAND_W,
            &mut prev_h,
            &mut prev_e,
            &mut cur_h,
            &mut cur_e,
        );
        scores.push(sw_bits);
    }

    if scores.len() < n_cal {
        let avg_gc = f64::midpoint(qs.gc, ts.gc);
        return calibrate_null_crc_for_gc(avg_gc, n_cal, alpha, base_floor, &HmmParams::default());
    }
    calibrate_crc_threshold(&scores, alpha, base_floor)
}

#[allow(clippy::too_many_lines)]
/// Score a single PAF record against optional indexed FASTAs under an explicit bit-savings floor.
#[must_use]
pub fn score_paf_record_with_floor<S: std::hash::BuildHasher>(
    rec: &PafRecord,
    q_map: Option<&HashMap<String, IndexedSeq, S>>,
    t_map: Option<&HashMap<String, IndexedSeq, S>>,
    threshold_bits: f64,
) -> ScoreResult {
    let (Some(qm), Some(tm)) = (q_map, t_map) else {
        return score_from_counts_only(rec, threshold_bits);
    };
    let (Some((qs, q_off)), Some((ts, t_off))) =
        (lookup_contig(qm, &rec.qname), lookup_contig(tm, &rec.tname))
    else {
        return score_from_counts_only(rec, threshold_bits);
    };
    let reject_corrupt_cigar = ScoreResult {
        mdl_bits: 0.0,
        nonrepeat_bp: 0,
        threshold_bits,
    };
    let Some(cigar_iter) = rec.cigar_ops() else {
        return score_from_counts_only(rec, threshold_bits);
    };
    let (Some(q_lo), Some(q_hi), Some(t_lo), Some(t_hi)) = (
        q_off.checked_add(rec.qstart),
        q_off.checked_add(rec.qend),
        t_off.checked_add(rec.tstart),
        t_off.checked_add(rec.tend),
    ) else {
        return reject_corrupt_cigar;
    };
    if q_hi > qs.seq.len() || t_hi > ts.seq.len() || q_lo > q_hi || t_lo > t_hi {
        return reject_corrupt_cigar;
    }

    let pi_q = base_probs(qs.gc);
    let pi_t = base_probs(ts.gc);
    let is_rc = rec.strand == '-';
    let mut q_pos = if is_rc { q_hi } else { q_lo };
    let mut t_pos = t_lo;
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
                    q_lo,
                    q_hi,
                    t_hi,
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
                    let Some(next_q) = q_pos.checked_sub(n).filter(|&p| p >= q_lo) else {
                        return reject_corrupt_cigar;
                    };
                    q_pos = next_q;
                } else {
                    let Some(next_q) = q_pos.checked_add(n).filter(|&p| p <= q_hi) else {
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
                let Some(next_t) = t_pos.checked_add(n).filter(|&p| p <= t_hi) else {
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

    let expected_q = if is_rc { q_lo } else { q_hi };
    if q_pos != expected_q || t_pos != t_hi {
        return reject_corrupt_cigar;
    }

    ScoreResult {
        mdl_bits,
        nonrepeat_bp,
        threshold_bits,
    }
}

fn base_floor_for_record(rec: &PafRecord) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let search_space_bits = ((rec.qlen.max(1) as f64) * (rec.tlen.max(1) as f64))
        .log2()
        .max(10.0);
    (KARLIN_MDL_COEFF * search_space_bits).max(10.0)
}

type CrcCache = HashMap<String, HashMap<String, f64>>;

#[allow(clippy::too_many_arguments)]
fn flush_paf_batch(
    batch: &mut Vec<PafRecord>,
    cfg: &PafFilterConfig,
    q_map: Option<&HashMap<String, IndexedSeq>>,
    t_ref: Option<&HashMap<String, IndexedSeq>>,
    crc_cache: &mut CrcCache,
    net_buf: &mut Vec<(PafRecord, ScoreResult)>,
    out: &mut crate::io_util::TextWriter,
    stats: &mut PafFilterStats,
) -> Result<()> {
    if batch.is_empty() {
        return Ok(());
    }
    if cfg.min_bits_override.is_none()
        && cfg.crc_n_cal > 0
        && let (Some(qm), Some(tm)) = (q_map, t_ref)
    {
        let mut missing: Vec<(String, String, f64)> = Vec::new();
        let mut seen: HashSet<(&str, &str)> = HashSet::new();
        for rec in batch.iter() {
            let cached = crc_cache
                .get(rec.qname.as_str())
                .is_some_and(|m| m.contains_key(rec.tname.as_str()));
            if !cached
                && seen.insert((rec.qname.as_str(), rec.tname.as_str()))
                && lookup_contig(qm, &rec.qname).is_some()
                && lookup_contig(tm, &rec.tname).is_some()
            {
                missing.push((
                    rec.qname.clone(),
                    rec.tname.clone(),
                    base_floor_for_record(rec),
                ));
            }
        }
        if !missing.is_empty() {
            let computed: Vec<(String, String, f64)> = missing
                .into_par_iter()
                .filter_map(|(qn, tn, base_floor)| {
                    let (qs, _) = lookup_contig(qm, &qn)?;
                    let (ts, _) = lookup_contig(tm, &tn)?;
                    let val = calibrate_paf_operator_crc(
                        qs,
                        ts,
                        cfg.crc_n_cal,
                        cfg.crc_alpha,
                        base_floor,
                    );
                    Some((qn, tn, val))
                })
                .collect();
            for (qn, tn, val) in computed {
                crc_cache.entry(qn).or_default().insert(tn, val);
            }
        }
    }

    let scores: Vec<ScoreResult> = batch
        .par_iter()
        .map(|rec| {
            let base_floor = base_floor_for_record(rec);
            let mut threshold_bits = cfg.min_bits_override.unwrap_or(base_floor);
            if cfg.min_bits_override.is_none()
                && cfg.crc_n_cal > 0
                && let Some(&crc_floor) = crc_cache
                    .get(rec.qname.as_str())
                    .and_then(|m| m.get(rec.tname.as_str()))
            {
                threshold_bits = threshold_bits.max(crc_floor);
            }
            score_paf_record_with_floor(rec, q_map, t_ref, threshold_bits)
        })
        .collect();

    for (rec, res) in batch.drain(..).zip(scores) {
        let pass_len = rec.aln_len >= cfg.min_aln_len && res.nonrepeat_bp >= cfg.min_nonrepeat_bp;
        let pass_mdl = res.mdl_bits >= res.threshold_bits;
        let pass = pass_len && pass_mdl;

        if !pass {
            if pass_len {
                stats.filtered_mdl += 1;
            } else {
                stats.filtered_len += 1;
            }
        }

        if cfg.max_occ.is_some() {
            if pass {
                net_buf.push((rec, res));
            } else if cfg.keep_all {
                rec.write_line_with_mdl(out, res.mdl_bits, 1)?;
            }
            continue;
        }

        if pass || cfg.keep_all {
            if pass {
                stats.kept += 1;
            }
            let hb = u8::from(!pass);
            rec.write_line_with_mdl(out, res.mdl_bits, hb)?;
        }
    }
    Ok(())
}

/// Stream, score, filter, and optionally net a PAF file according to `cfg`.
///
/// # Errors
/// Returns an error if reading inputs or writing the filtered PAF fails.
pub fn run_paf_filter(cfg: &PafFilterConfig) -> Result<PafFilterStats> {
    const BATCH_SIZE: usize = 1024;
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
    let mut stats = PafFilterStats::default();
    let mut net_buf: Vec<(PafRecord, ScoreResult)> = Vec::new();
    let mut crc_cache: CrcCache = HashMap::new();
    let mut batch: Vec<PafRecord> = Vec::with_capacity(BATCH_SIZE);

    stream_paf_records(&cfg.paf_path, |rec| {
        stats.total += 1;
        batch.push(rec);
        if batch.len() >= BATCH_SIZE {
            flush_paf_batch(
                &mut batch,
                cfg,
                q_map.as_ref(),
                t_ref,
                &mut crc_cache,
                &mut net_buf,
                &mut out,
                &mut stats,
            )?;
        }
        Ok(ControlFlow::Continue(()))
    })?;
    flush_paf_batch(
        &mut batch,
        cfg,
        q_map.as_ref(),
        t_ref,
        &mut crc_cache,
        &mut net_buf,
        &mut out,
        &mut stats,
    )?;

    if let Some(k) = cfg.max_occ {
        let (netted, dropped) = net_top_k_partition_with_mode(net_buf, k, cfg.net_score);
        stats.kept = netted.len();
        for (rec, res) in netted {
            rec.write_line_with_mdl(&mut out, res.mdl_bits, 0)?;
        }
        if cfg.keep_all {
            for (rec, res) in dropped {
                rec.write_line_with_mdl(&mut out, res.mdl_bits, 1)?;
            }
        }
    }

    out.finish()?;
    Ok(stats)
}
