//! Top-k interval netting and identity-weighted reciprocal-best filtering for scored PAF records.

use super::{NetScoreMode, ScoreResult};
use crate::paf::PafRecord;
use std::collections::HashMap;

#[derive(Default)]
struct BinnedIntervals {
    ivs: Vec<(usize, usize)>,
    bins: HashMap<usize, Vec<usize>>,
    large_ivs: Vec<usize>,
    seen_stamp: Vec<u32>,
    cur_stamp: u32,
    scratch_spans: Vec<(usize, usize)>,
}

impl BinnedIntervals {
    const BIN_SHIFT: usize = 16;
    const MAX_BINS_PER_IV: usize = 1024;

    fn next_stamp(&mut self) -> u32 {
        if self.seen_stamp.len() < self.ivs.len() {
            self.seen_stamp.resize(self.ivs.len(), 0);
        }
        self.cur_stamp = self.cur_stamp.wrapping_add(1);
        if self.cur_stamp == 0 {
            self.seen_stamp.fill(0);
            self.cur_stamp = 1;
        }
        self.cur_stamp
    }

    fn overlap_count(&mut self, s: usize, e: usize, max_occ: usize) -> usize {
        let len = e.saturating_sub(s);
        if len == 0 {
            return usize::MAX;
        }
        let half_len = len.div_ceil(2);
        let b0 = s >> Self::BIN_SHIFT;
        let b1 = (e - 1) >> Self::BIN_SHIFT;
        let stamp = self.next_stamp();
        let mut count = 0usize;

        let mut check_id = |id: usize| -> bool {
            if self.seen_stamp[id] == stamp {
                return false;
            }
            self.seen_stamp[id] = stamp;
            let (os, oe) = self.ivs[id];
            let ov = os.max(s);
            let oe_min = oe.min(e);
            if oe_min > ov && (oe_min - ov) >= half_len {
                count += 1;
                return count >= max_occ;
            }
            false
        };

        for &id in &self.large_ivs {
            if check_id(id) {
                return count;
            }
        }
        if b1.saturating_sub(b0) >= Self::MAX_BINS_PER_IV {
            for id in 0..self.ivs.len() {
                if check_id(id) {
                    return count;
                }
            }
            return count;
        }
        for bin in b0..=b1 {
            let Some(ids) = self.bins.get(&bin) else {
                continue;
            };
            for &id in ids {
                if check_id(id) {
                    return count;
                }
            }
        }
        count
    }

    /// Return `true` if the union of already-inserted intervals covers at least `threshold` bases
    /// of `[s, e)`. Short-circuits on any single interval covering `>= threshold` with zero alloc.
    fn covered_bp_at_least(&mut self, s: usize, e: usize, threshold: usize) -> bool {
        if threshold == 0 {
            return true;
        }
        if e <= s || (e - s) < threshold {
            return false;
        }
        let b0 = s >> Self::BIN_SHIFT;
        let b1 = (e - 1) >> Self::BIN_SHIFT;
        let stamp = self.next_stamp();
        self.scratch_spans.clear();

        let mut push_id = |id: usize| -> bool {
            if self.seen_stamp[id] == stamp {
                return false;
            }
            self.seen_stamp[id] = stamp;
            let (os, oe) = self.ivs[id];
            let ov = os.max(s);
            let oe_min = oe.min(e);
            if oe_min > ov {
                if oe_min - ov >= threshold {
                    return true;
                }
                self.scratch_spans.push((ov, oe_min));
            }
            false
        };

        for &id in &self.large_ivs {
            if push_id(id) {
                return true;
            }
        }
        if b1.saturating_sub(b0) >= Self::MAX_BINS_PER_IV {
            for id in 0..self.ivs.len() {
                if push_id(id) {
                    return true;
                }
            }
        } else {
            for bin in b0..=b1 {
                if let Some(ids) = self.bins.get(&bin) {
                    for &id in ids {
                        if push_id(id) {
                            return true;
                        }
                    }
                }
            }
        }
        if self.scratch_spans.is_empty() {
            return false;
        }
        self.scratch_spans.sort_unstable();
        let mut total = 0usize;
        let (mut cur_s, mut cur_e) = self.scratch_spans[0];
        for &(ns, ne) in self.scratch_spans.iter().skip(1) {
            if ns <= cur_e {
                cur_e = cur_e.max(ne);
            } else {
                total += cur_e - cur_s;
                if total >= threshold {
                    return true;
                }
                cur_s = ns;
                cur_e = ne;
            }
        }
        total + (cur_e - cur_s) >= threshold
    }

    fn insert(&mut self, s: usize, e: usize) {
        if e <= s {
            return;
        }
        let id = self.ivs.len();
        self.ivs.push((s, e));
        let b0 = s >> Self::BIN_SHIFT;
        let b1 = (e - 1) >> Self::BIN_SHIFT;
        if b1.saturating_sub(b0) >= Self::MAX_BINS_PER_IV {
            self.large_ivs.push(id);
        } else {
            for bin in b0..=b1 {
                self.bins.entry(bin).or_default().push(id);
            }
        }
    }
}

fn paf_tie_key(
    r: &PafRecord,
) -> (
    &str,
    usize,
    usize,
    &str,
    usize,
    usize,
    char,
    usize,
    usize,
    u8,
) {
    (
        r.qname.as_str(),
        r.qstart,
        r.qend,
        r.tname.as_str(),
        r.tstart,
        r.tend,
        r.strand,
        r.matches,
        r.aln_len,
        r.mapq,
    )
}

pub(super) type ScoredPafVec = Vec<(PafRecord, ScoreResult)>;

#[allow(clippy::cast_precision_loss)]
fn net_rank_score(rec: &PafRecord, res: &ScoreResult, mode: NetScoreMode) -> f64 {
    match mode {
        NetScoreMode::Mdl => res.mdl_bits,
        NetScoreMode::LogLengthAni => {
            let aln_len = rec.aln_len.max(1) as f64;
            let ani = (rec.matches.min(rec.aln_len) as f64 / aln_len).clamp(0.0, 1.0);
            let eff_len = rec.aln_len.max(2) as f64;
            // Identity-weighted sublinear-length ranking matching sweepga's `log-length-ani`,
            // ensuring high-identity orthologs outrank longer low-identity palindromes.
            ani * ani * eff_len.ln()
        }
    }
}

fn single_axis_pass(
    scored: &[(PafRecord, ScoreResult)],
    max_occ: usize,
    on_query: bool,
) -> Vec<bool> {
    let mut occ: HashMap<&str, BinnedIntervals> = HashMap::new();
    let mut keep = vec![false; scored.len()];
    for (idx, (rec, _)) in scored.iter().enumerate() {
        // Reject degenerate zero-length alignments on either axis before mutating axis state so
        // an invalid half-empty record cannot shadow valid alignments on the opposite axis.
        if rec.qend <= rec.qstart || rec.tend <= rec.tstart {
            continue;
        }
        let (name, start, end) = if on_query {
            (rec.qname.as_str(), rec.qstart, rec.qend)
        } else {
            (rec.tname.as_str(), rec.tstart, rec.tend)
        };
        let len = end - start;
        let shadowed = if max_occ == 1 {
            occ.get_mut(name)
                .is_some_and(|b| b.covered_bp_at_least(start, end, len.div_ceil(2)))
        } else {
            occ.get_mut(name)
                .is_some_and(|b| b.overlap_count(start, end, max_occ) >= max_occ)
        };
        if !shadowed {
            occ.entry(name).or_default().insert(start, end);
            keep[idx] = true;
        }
    }
    keep
}

pub(super) fn net_top_k_partition_with_mode(
    mut scored: ScoredPafVec,
    max_occ: usize,
    mode: NetScoreMode,
) -> (ScoredPafVec, ScoredPafVec) {
    if max_occ == 0 {
        scored.sort_unstable_by(|a, b| paf_tie_key(&a.0).cmp(&paf_tie_key(&b.0)));
        return (Vec::new(), scored);
    }
    scored.sort_unstable_by(|a, b| {
        net_rank_score(&b.0, &b.1, mode)
            .total_cmp(&net_rank_score(&a.0, &a.1, mode))
            .then_with(|| b.1.mdl_bits.total_cmp(&a.1.mdl_bits))
            .then_with(|| paf_tie_key(&a.0).cmp(&paf_tie_key(&b.0)))
    });

    let mut kept = Vec::new();
    let mut dropped = Vec::new();

    if mode == NetScoreMode::LogLengthAni {
        // Symmetric two-axis reciprocal-best netting (query-best ∩ target-best): a deleted-region
        // paralog that loses on target cannot fall through to steal an unrelated free target copy.
        let keep_q = single_axis_pass(&scored, max_occ, true);
        let keep_t = single_axis_pass(&scored, max_occ, false);
        for (idx, pair) in scored.into_iter().enumerate() {
            if keep_q[idx] && keep_t[idx] {
                kept.push(pair);
            } else {
                dropped.push(pair);
            }
        }
    } else {
        let mut occ_q: HashMap<String, BinnedIntervals> = HashMap::new();
        let mut occ_t: HashMap<String, BinnedIntervals> = HashMap::new();
        for (rec, res) in scored {
            let q_len = rec.qend.saturating_sub(rec.qstart);
            let t_len = rec.tend.saturating_sub(rec.tstart);
            let shadowed = q_len == 0
                || t_len == 0
                || occ_q
                    .get_mut(&rec.qname)
                    .is_some_and(|idx| idx.overlap_count(rec.qstart, rec.qend, max_occ) >= max_occ)
                || occ_t
                    .get_mut(&rec.tname)
                    .is_some_and(|idx| idx.overlap_count(rec.tstart, rec.tend, max_occ) >= max_occ);
            if shadowed {
                dropped.push((rec, res));
                continue;
            }
            occ_q
                .entry(rec.qname.clone())
                .or_default()
                .insert(rec.qstart, rec.qend);
            occ_t
                .entry(rec.tname.clone())
                .or_default()
                .insert(rec.tstart, rec.tend);
            kept.push((rec, res));
        }
    }

    kept.sort_unstable_by(|a, b| paf_tie_key(&a.0).cmp(&paf_tie_key(&b.0)));
    dropped.sort_unstable_by(|a, b| paf_tie_key(&a.0).cmp(&paf_tie_key(&b.0)));
    (kept, dropped)
}

/// Net scored PAF records deterministically using cumulative MDL ranking, retaining at most
/// `max_occ` records overlapping `>= 50%` on either query or target.
#[must_use]
pub fn net_top_k_records(
    scored: Vec<(PafRecord, ScoreResult)>,
    max_occ: usize,
) -> Vec<(PafRecord, ScoreResult)> {
    net_top_k_partition_with_mode(scored, max_occ, NetScoreMode::Mdl).0
}

/// Net scored PAF records deterministically under the requested [`NetScoreMode`] ranking,
/// retaining at most `max_occ` records overlapping `>= 50%` on either query or target.
#[must_use]
pub fn net_top_k_records_with_mode(
    scored: Vec<(PafRecord, ScoreResult)>,
    max_occ: usize,
    mode: NetScoreMode,
) -> Vec<(PafRecord, ScoreResult)> {
    net_top_k_partition_with_mode(scored, max_occ, mode).0
}
