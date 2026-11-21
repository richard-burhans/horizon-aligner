//! Transitive sister rescue and ancestral consensus rescue for `horizon_guide`.

use crate::{ContigList, IndexedSpecies, align_species_pair};
use anyhow::Result;
use horizon_aligner::{
    AlignConfig, ContigIndex, complement_base, parse_two_row_maf_blocks, prepare_contig_index_opts,
    shared_sketch_count,
};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fmt::Write as _;

#[derive(Clone, Copy)]
struct MidHingeRun {
    m_start: usize,
    m_end: usize,
    h_cidx: usize,
    h_start: usize,
    mid_rev: bool,
}

fn build_mid_to_hinge_runs(
    hinge_idx_by_name: &HashMap<&str, usize>,
    mid_idx_by_name: &HashMap<&str, usize>,
    n_mid_contigs: usize,
    hinge_to_mid_maf: &str,
) -> Vec<Vec<MidHingeRun>> {
    let mut by_mid: Vec<Vec<MidHingeRun>> = vec![Vec::new(); n_mid_contigs];
    for (h_row, m_row) in parse_two_row_maf_blocks(hinge_to_mid_maf) {
        let Some(&h_cidx) = hinge_idx_by_name.get(h_row.src.as_str()) else {
            continue;
        };
        let Some(&m_cidx) = mid_idx_by_name.get(m_row.src.as_str()) else {
            continue;
        };
        let mid_rev = (h_row.strand == '-') ^ (m_row.strand == '-');
        let m_is_rev = m_row.strand == '-';
        let mut hk = 0usize;
        let mut mk = 0usize;
        let mut run: Option<(usize, usize, usize)> = None; // (m_start_pos, h_start_pos, len)
        let flush_run = |run: Option<(usize, usize, usize)>, out: &mut Vec<MidHingeRun>| {
            if let Some((ms, hs, len)) = run {
                if m_is_rev {
                    let m_lo = ms + 1 - len;
                    let h_at_m_lo = if mid_rev { hs + len - 1 } else { hs + 1 - len };
                    out.push(MidHingeRun {
                        m_start: m_lo,
                        m_end: ms + 1,
                        h_cidx,
                        h_start: h_at_m_lo,
                        mid_rev,
                    });
                } else {
                    out.push(MidHingeRun {
                        m_start: ms,
                        m_end: ms + len,
                        h_cidx,
                        h_start: hs,
                        mid_rev,
                    });
                }
            }
        };
        for (hb, mb) in h_row.seq.bytes().zip(m_row.seq.bytes()) {
            if hb != b'-' && mb != b'-' {
                let hp = h_row.fwd_pos(hk);
                let mp = m_row.fwd_pos(mk);
                match run {
                    Some((ms, hs, len))
                        if (if m_is_rev {
                            ms.checked_sub(len) == Some(mp)
                        } else {
                            ms + len == mp
                        }) && (if h_row.strand == '-' {
                            hs.checked_sub(len) == Some(hp)
                        } else {
                            hs + len == hp
                        }) =>
                    {
                        run = Some((ms, hs, len + 1));
                    }
                    other => {
                        flush_run(other, &mut by_mid[m_cidx]);
                        run = Some((mp, hp, 1));
                    }
                }
            } else if let Some(prev) = run.take() {
                flush_run(Some(prev), &mut by_mid[m_cidx]);
            }
            if hb != b'-' {
                hk += 1;
            }
            if mb != b'-' {
                mk += 1;
            }
        }
        flush_run(run, &mut by_mid[m_cidx]);
    }
    for runs in &mut by_mid {
        runs.sort_by_key(|r| r.m_start);
    }
    by_mid
}

fn lookup_mid_to_hinge(
    runs: &[MidHingeRun],
    mp: usize,
    hinge_contigs: &ContigList,
) -> Option<(usize, usize, u8, bool)> {
    let idx = runs.partition_point(|r| r.m_start <= mp);
    for r in runs[..idx].iter().rev() {
        if mp < r.m_end {
            let d = mp - r.m_start;
            let hp = if r.mid_rev {
                r.h_start - d
            } else {
                r.h_start + d
            };
            let hb = hinge_contigs
                .get(r.h_cidx)
                .and_then(|(_, seq)| seq.get(hp))
                .copied()
                .unwrap_or(b'N')
                .to_ascii_uppercase();
            return Some((r.h_cidx, hp, hb, r.mid_rev));
        }
    }
    None
}

fn build_covered_dst(
    hinge_contigs: &ContigList,
    hinge_idx_by_name: &HashMap<&str, usize>,
    dst_direct_maf: &str,
) -> Vec<Vec<u64>> {
    let mut covered_dst: Vec<Vec<u64>> = hinge_contigs
        .iter()
        .map(|(_, seq)| vec![0u64; seq.len().div_ceil(64)])
        .collect();
    for (h_row, _) in parse_two_row_maf_blocks(dst_direct_maf) {
        let Some(&h_cidx) = hinge_idx_by_name.get(h_row.src.as_str()) else {
            continue;
        };
        let ctg_len = hinge_contigs[h_cidx].1.len();
        let mut hk = 0usize;
        for hb in h_row.seq.bytes() {
            if hb != b'-' {
                let hp = h_row.fwd_pos(hk);
                if hp < ctg_len {
                    covered_dst[h_cidx][hp >> 6] |= 1u64 << (hp & 63);
                }
                hk += 1;
            }
        }
    }
    covered_dst
}

#[inline]
fn is_bitset_covered(bits: &[u64], pos: usize) -> bool {
    bits.get(pos >> 6)
        .is_some_and(|&word| (word & (1u64 << (pos & 63))) != 0)
}

/// Project `src -> mid -> dst` via `hinge_to_mid_maf` and `mid_to_dst_maf`, emitting rescued
/// `hinge <-> dst` MAF blocks for hinge positions not already covered in `dst_direct_maf`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn rescue_transitive_maf(
    hinge_contigs: &ContigList,
    mid_contigs: &ContigList,
    hinge_to_mid_maf: &str,
    mid_to_dst_maf: &str,
    dst_direct_maf: &str,
    mid_is_first_in_pair: bool,
    min_run_len: usize,
    min_ident_pct: usize,
) -> String {
    let hinge_idx_by_name: HashMap<&str, usize> = hinge_contigs
        .iter()
        .enumerate()
        .map(|(idx, (name, _))| (name.as_str(), idx))
        .collect();
    let mid_idx_by_name: HashMap<&str, usize> = mid_contigs
        .iter()
        .enumerate()
        .map(|(idx, (name, _))| (name.as_str(), idx))
        .collect();

    let mid_to_hinge = build_mid_to_hinge_runs(
        &hinge_idx_by_name,
        &mid_idx_by_name,
        mid_contigs.len(),
        hinge_to_mid_maf,
    );
    if mid_to_hinge.iter().all(Vec::is_empty) {
        return String::new();
    }

    // Track hinge positions already covered by direct hinge <-> dst alignment.
    let covered_dst = build_covered_dst(hinge_contigs, &hinge_idx_by_name, dst_direct_maf);

    // Collect projected (hinge_pos, dst_pos, h_base, dst_base_oriented) points grouped by
    // (h_cidx, dst_src, dst_src_size, dst_rel_rev).
    #[allow(clippy::type_complexity)]
    let mut projected: HashMap<(usize, String, usize, bool), Vec<(usize, usize, u8, u8)>> =
        HashMap::new();
    for (r1, r2) in parse_two_row_maf_blocks(mid_to_dst_maf) {
        let (m_row, d_row) = if mid_is_first_in_pair {
            (r1, r2)
        } else {
            (r2, r1)
        };
        let Some(&m_cidx) = mid_idx_by_name.get(m_row.src.as_str()) else {
            continue;
        };
        let mut mk = 0usize;
        let mut dk = 0usize;
        for (mb, db) in m_row.seq.bytes().zip(d_row.seq.bytes()) {
            if mb != b'-' && db != b'-' {
                let mp = m_row.fwd_pos(mk);
                let dp = d_row.fwd_pos(dk);
                if let Some((h_cidx, hp, hb, mid_rev)) =
                    lookup_mid_to_hinge(&mid_to_hinge[m_cidx], mp, hinge_contigs)
                    && !is_bitset_covered(&covered_dst[h_cidx], hp)
                {
                    let dst_rel_rev = mid_rev ^ (m_row.strand == '-') ^ (d_row.strand == '-');
                    let db_oriented = if (d_row.strand == '-') == dst_rel_rev {
                        db
                    } else {
                        complement_base(db)
                    };
                    projected
                        .entry((h_cidx, d_row.src.clone(), d_row.src_size, dst_rel_rev))
                        .or_default()
                        .push((hp, dp, hb, db_oriented));
                }
            }
            if mb != b'-' {
                mk += 1;
            }
            if db != b'-' {
                dk += 1;
            }
        }
    }

    let mut out = String::new();
    let mut keys: Vec<_> = projected.keys().cloned().collect();
    keys.sort();
    for key in keys {
        let (h_cidx, d_src, d_src_size, dst_rel_rev) = key;
        let Some(pts) = projected.remove(&(h_cidx, d_src.clone(), d_src_size, dst_rel_rev)) else {
            continue;
        };
        let (h_name, h_seq) = &hinge_contigs[h_cidx];
        emit_rescued_runs(
            &mut out,
            h_name,
            h_seq.len(),
            &d_src,
            d_src_size,
            dst_rel_rev,
            pts,
            min_run_len,
            min_ident_pct,
        );
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn emit_rescued_runs(
    out: &mut String,
    h_name: &str,
    h_src_size: usize,
    d_src: &str,
    d_src_size: usize,
    dst_rel_rev: bool,
    mut pts: Vec<(usize, usize, u8, u8)>,
    min_run_len: usize,
    min_ident_pct: usize,
) {
    pts.sort_unstable_by_key(|p| (p.0, p.1));
    pts.dedup_by_key(|p| p.0);

    let mut run_start = 0usize;
    while run_start < pts.len() {
        let mut run_end = run_start + 1;
        while run_end < pts.len() {
            let prev = pts[run_end - 1];
            let cur = pts[run_end];
            let expected_dp = if dst_rel_rev {
                prev.1.checked_sub(1)
            } else {
                Some(prev.1 + 1)
            };
            if cur.0 == prev.0 + 1 && Some(cur.1) == expected_dp {
                run_end += 1;
            } else {
                break;
            }
        }
        let run = &pts[run_start..run_end];
        let matches = run
            .iter()
            .filter(|p| p.2.eq_ignore_ascii_case(&p.3))
            .count();
        let mismatches = run.len() - matches;
        let net_score = matches.cast_signed() * 10 - mismatches.cast_signed() * 12;
        if run.len() >= min_run_len
            && matches * 100 >= run.len() * min_ident_pct
            && net_score >= 100
        {
            let h_start = run[0].0;
            let size = run.len();
            let d_min = run.iter().map(|p| p.1).min().unwrap_or(0);
            let d_max = run.iter().map(|p| p.1).max().unwrap_or(0);
            let d_start = if dst_rel_rev {
                d_src_size.saturating_sub(1 + d_max)
            } else {
                d_min
            };
            let d_strand = if dst_rel_rev { '-' } else { '+' };
            let h_aln: String = run.iter().map(|p| p.2 as char).collect();
            let d_aln: String = run.iter().map(|p| p.3 as char).collect();
            let _ = writeln!(out, "a score={size}.0");
            let _ = writeln!(out, "s {h_name} {h_start} {size} + {h_src_size} {h_aln}");
            let _ = writeln!(
                out,
                "s {d_src} {d_start} {size} {d_strand} {d_src_size} {d_aln}\n"
            );
        }
        run_start = run_end;
    }
}

fn maf_ref_coverage(hinge_contigs: &ContigList, maf_text: &str) -> f64 {
    let total_bp: usize = hinge_contigs.iter().map(|(_, s)| s.len()).sum();
    if total_bp == 0 {
        return 0.0;
    }
    let ctg_map: HashMap<&str, usize> = hinge_contigs
        .iter()
        .enumerate()
        .map(|(i, (n, _))| (n.as_str(), i))
        .collect();
    let mut intervals: Vec<Vec<(usize, usize)>> = vec![Vec::new(); hinge_contigs.len()];
    for (r0, _) in parse_two_row_maf_blocks(maf_text) {
        if r0.strand != '+' {
            continue;
        }
        if let Some(&ci) = ctg_map.get(r0.src.as_str()) {
            let r0_size = r0.seq.bytes().filter(|&b| b != b'-').count();
            let s = r0.start;
            let e = (r0.start + r0_size).min(hinge_contigs[ci].1.len());
            if s < e {
                intervals[ci].push((s, e));
            }
        }
    }
    let cov_bp: usize = intervals
        .iter_mut()
        .map(|ivs| {
            if ivs.is_empty() {
                return 0;
            }
            ivs.sort_unstable();
            let mut sum = 0usize;
            let (mut cur_s, mut cur_e) = ivs[0];
            for &(s, e) in &ivs[1..] {
                if s <= cur_e {
                    cur_e = cur_e.max(e);
                } else {
                    sum += cur_e - cur_s;
                    cur_s = s;
                    cur_e = e;
                }
            }
            sum + (cur_e - cur_s)
        })
        .sum();
    #[allow(clippy::cast_precision_loss)]
    let frac = (cov_bp as f64) / (total_bp as f64);
    frac
}

fn base_to_idx(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

pub(crate) fn build_ancestral_consensus_contigs(
    hinge_contigs: &ContigList,
    high_cov_mafs: &[&str],
    need_contig: &[bool],
) -> ContigList {
    const BASES: [u8; 4] = *b"ACGT";
    let ctg_map: HashMap<&str, usize> = hinge_contigs
        .iter()
        .enumerate()
        .filter(|&(i, _)| need_contig.get(i).copied().unwrap_or(true))
        .map(|(i, (n, _))| (n.as_str(), i))
        .collect();
    let mut counts: Vec<Vec<[u32; 4]>> = hinge_contigs
        .iter()
        .enumerate()
        .map(|(i, (_, seq))| {
            if !need_contig.get(i).copied().unwrap_or(true) {
                return Vec::new();
            }
            seq.iter()
                .map(|&b| {
                    let mut c = [0u32; 4];
                    if let Some(bi) = base_to_idx(b) {
                        c[bi] = 1;
                    }
                    c
                })
                .collect()
        })
        .collect();

    for maf_text in high_cov_mafs {
        for (r0, r1) in parse_two_row_maf_blocks(maf_text) {
            if r0.strand != '+' {
                continue;
            }
            let Some(&ci) = ctg_map.get(r0.src.as_str()) else {
                continue;
            };
            let mut rp = r0.start;
            for (rc, qc) in r0.seq.bytes().zip(r1.seq.bytes()) {
                if rc != b'-' {
                    if qc != b'-'
                        && rp < counts[ci].len()
                        && let Some(bi) = base_to_idx(qc)
                    {
                        counts[ci][rp][bi] += 1;
                    }
                    rp += 1;
                }
            }
        }
    }

    hinge_contigs
        .iter()
        .enumerate()
        .map(|(ci, (name, ref_seq))| {
            if counts[ci].is_empty() {
                return (name.clone(), Vec::new());
            }
            let mut cons = Vec::with_capacity(ref_seq.len());
            for (pos, &rb) in ref_seq.iter().enumerate() {
                let c = counts[ci][pos];
                let rb_upper = rb.to_ascii_uppercase();
                let rb_idx = BASES.iter().position(|&b| b == rb_upper);
                let mut best_b = rb;
                let mut best_cnt = rb_idx.map_or(0, |ri| c[ri]);
                for (bi, &cnt) in c.iter().enumerate() {
                    if cnt > best_cnt {
                        best_cnt = cnt;
                        best_b = BASES[bi];
                    }
                }
                cons.push(best_b);
            }
            (name.clone(), cons)
        })
        .collect()
}

pub(crate) fn remap_consensus_maf_to_hinge(
    hinge_contigs: &ContigList,
    cons_maf: &str,
    existing_maf: &str,
) -> String {
    let ctg_map: HashMap<&str, &[u8]> = hinge_contigs
        .iter()
        .map(|(n, s)| (n.as_str(), s.as_slice()))
        .collect();
    let mut covered_by_contig: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    for (r0, _) in parse_two_row_maf_blocks(existing_maf) {
        if r0.strand != '+' {
            continue;
        }
        let r0_size = r0.seq.bytes().filter(|&b| b != b'-').count();
        if r0_size > 0 {
            covered_by_contig
                .entry(r0.src)
                .or_default()
                .push((r0.start, r0.start + r0_size));
        }
    }
    let mut out = String::new();
    for (r0, r1) in parse_two_row_maf_blocks(cons_maf) {
        let Some(&ref_seq) = ctg_map.get(r0.src.as_str()) else {
            continue;
        };
        let r0_size = r0.seq.bytes().filter(|&b| b != b'-').count();
        if r0_size == 0 {
            continue;
        }
        if r0.strand != '+' {
            continue;
        }
        let s0 = r0.start;
        let e0 = s0 + r0_size;
        if let Some(spans) = covered_by_contig.get(&r0.src) {
            let mut clipped: Vec<(usize, usize)> = spans
                .iter()
                .filter(|&&(a, b)| a < e0 && b > s0)
                .map(|&(a, b)| (a.max(s0), b.min(e0)))
                .collect();
            clipped.sort_unstable();
            let mut ov = 0usize;
            let mut cur: Option<(usize, usize)> = None;
            for (a, b) in clipped {
                match cur {
                    Some((cs, ce)) if a <= ce => cur = Some((cs, ce.max(b))),
                    Some((cs, ce)) => {
                        ov += ce - cs;
                        cur = Some((a, b));
                    }
                    None => cur = Some((a, b)),
                }
            }
            if let Some((cs, ce)) = cur {
                ov += ce - cs;
            }
            if ov * 10 >= r0_size * 9 {
                continue;
            }
        }
        let mut remapped_r0 = String::with_capacity(r0.seq.len());
        let mut rp = r0.start;
        for ch in r0.seq.bytes() {
            if ch == b'-' {
                remapped_r0.push('-');
            } else {
                let rb = ref_seq.get(rp).copied().unwrap_or(ch);
                remapped_r0.push(char::from(rb));
                rp += 1;
            }
        }
        let r1_size = r1.seq.bytes().filter(|&b| b != b'-').count();
        covered_by_contig
            .entry(r0.src.clone())
            .or_default()
            .push((s0, e0));
        let _ = writeln!(out, "a score=1.0");
        let _ = writeln!(
            out,
            "s {} {} {} {} {} {}",
            r0.src, r0.start, r0_size, r0.strand, r0.src_size, remapped_r0
        );
        let _ = writeln!(
            out,
            "s {} {} {} {} {} {}",
            r1.src, r1.start, r1_size, r1.strand, r1.src_size, r1.seq
        );
        out.push('\n');
    }
    out
}

/// Build a bottom-up UPGMA guide-tree ordering of rescue pairs `(rep_a, rep_b)` over non-hinge
/// species using `MinHash` sketch Jaccard similarity, skipping pairs where both clades already
/// have high hinge coverage.
fn build_upgma_rescue_pairs(
    n_species: usize,
    hinge_idx: usize,
    sketches: &[Vec<u64>],
    high_cov: &[bool],
) -> Vec<(usize, usize)> {
    let mut clades: Vec<(usize, Vec<usize>, bool)> = (0..n_species)
        .filter(|&i| i != hinge_idx)
        .map(|i| (i, vec![i], high_cov[i]))
        .collect();
    let mut pairs = Vec::new();
    while clades.len() >= 2 {
        let mut best: Option<(f64, usize, usize)> = None;
        for a in 0..clades.len() {
            for b in (a + 1)..clades.len() {
                let rep_a = clades[a].0;
                let rep_b = clades[b].0;
                let shared = shared_sketch_count(&sketches[rep_a], &sketches[rep_b]);
                let union_est =
                    (sketches[rep_a].len() + sketches[rep_b].len()).saturating_sub(shared);
                #[allow(clippy::cast_precision_loss)]
                let sim = if union_est == 0 {
                    0.0
                } else {
                    (shared as f64) / (union_est as f64)
                };
                if best.is_none_or(|(best_sim, _, _)| sim > best_sim) {
                    best = Some((sim, a, b));
                }
            }
        }
        let Some((_, a, b)) = best else {
            break;
        };
        let (ca, cb) = (clades.swap_remove(b), clades.swap_remove(a));
        if !(ca.2 && cb.2) {
            pairs.push((ca.0, cb.0));
        }
        let merged_high = ca.2 && cb.2;
        let mut members = ca.1;
        members.extend(cb.1);
        // Prefer a low-coverage representative so it isn't stranded by a high-coverage sibling,
        // breaking ties by larger sketch coverage (or lower index).
        let rep = match (ca.2, cb.2) {
            (false, true) => ca.0,
            (true, false) => cb.0,
            _ => {
                if sketches[ca.0].len() >= sketches[cb.0].len() {
                    ca.0
                } else {
                    cb.0
                }
            }
        };
        clades.push((rep, members, merged_high));
    }
    pairs
}

fn append_nonoverlapping_rescued_maf(
    hinge_idx_by_name: &HashMap<&str, usize>,
    covered: &mut [Vec<(usize, usize)>],
    target_maf: &mut String,
    rescued_maf: &str,
) {
    if rescued_maf.is_empty() {
        return;
    }
    for (h_row, d_row) in parse_two_row_maf_blocks(rescued_maf) {
        let Some(&h_cidx) = hinge_idx_by_name.get(h_row.src.as_str()) else {
            continue;
        };
        let size = h_row.seq.bytes().filter(|&b| b != b'-').count();
        if size == 0 {
            continue;
        }
        let s = h_row.start;
        let e = s + size;
        let ivs = &mut covered[h_cidx];
        // `ivs` is kept sorted by start; check the preceding and next intervals in O(log B).
        let pos = ivs.partition_point(|&(cs, _)| cs < e);
        let overlaps = pos > 0 && ivs[pos - 1].1 > s;
        if overlaps {
            continue;
        }
        let ins = ivs.partition_point(|&(cs, _)| cs < s);
        ivs.insert(ins, (s, e));
        let _ = writeln!(target_maf, "a score={size}.0");
        let _ = writeln!(
            target_maf,
            "s {} {} {size} {} {} {}",
            h_row.src, h_row.start, h_row.strand, h_row.src_size, h_row.seq
        );
        let d_size = d_row.seq.bytes().filter(|&b| b != b'-').count();
        let _ = writeln!(
            target_maf,
            "s {} {} {d_size} {} {} {}\n",
            d_row.src, d_row.start, d_row.strand, d_row.src_size, d_row.seq
        );
    }
}

pub(crate) fn run_progressive_sister_rescue(
    indexed: &[IndexedSpecies],
    sketches: &[Vec<u64>],
    hinge_idx: usize,
    sp_to_pw_idx: &HashMap<usize, usize>,
    pairwise_mafs: &mut [(String, String)],
    cfg: &AlignConfig,
    top_k: usize,
) -> Result<()> {
    let hinge_contigs = &indexed[hinge_idx].1;
    let hinge_idx_by_name: HashMap<&str, usize> = hinge_contigs
        .iter()
        .enumerate()
        .map(|(idx, (name, _))| (name.as_str(), idx))
        .collect();
    let high_cov: Vec<bool> = (0..indexed.len())
        .into_par_iter()
        .map(|i| {
            if i == hinge_idx {
                false
            } else {
                let pw = sp_to_pw_idx[&i];
                maf_ref_coverage(hinge_contigs, &pairwise_mafs[pw].1) >= 0.85
            }
        })
        .collect();
    let sister_pairs = build_upgma_rescue_pairs(indexed.len(), hinge_idx, sketches, &high_cov);

    let rescued_chunks: Vec<(usize, String, usize, String)> = sister_pairs
        .par_iter()
        .map(|&(i, j)| {
            let (sister_maf, _) = align_species_pair(
                &indexed[i].1,
                &indexed[i].2,
                &indexed[j].1,
                &indexed[j].2,
                cfg,
                top_k,
            )?;
            let pw_i = sp_to_pw_idx[&i];
            let pw_j = sp_to_pw_idx[&j];
            let (rescued_j, rescued_i) = rayon::join(
                || {
                    rescue_transitive_maf(
                        hinge_contigs,
                        &indexed[i].1,
                        &pairwise_mafs[pw_i].1,
                        &sister_maf,
                        &pairwise_mafs[pw_j].1,
                        true,
                        12,
                        75,
                    )
                },
                || {
                    rescue_transitive_maf(
                        hinge_contigs,
                        &indexed[j].1,
                        &pairwise_mafs[pw_j].1,
                        &sister_maf,
                        &pairwise_mafs[pw_i].1,
                        false,
                        12,
                        75,
                    )
                },
            );
            Ok((pw_j, rescued_j, pw_i, rescued_i))
        })
        .collect::<Result<Vec<_>>>()?;

    let mut covered_by_pw: Vec<Vec<Vec<(usize, usize)>>> =
        vec![vec![Vec::new(); hinge_contigs.len()]; pairwise_mafs.len()];
    for (pw_j, rescued_j, pw_i, rescued_i) in rescued_chunks {
        append_nonoverlapping_rescued_maf(
            &hinge_idx_by_name,
            &mut covered_by_pw[pw_j],
            &mut pairwise_mafs[pw_j].1,
            &rescued_j,
        );
        append_nonoverlapping_rescued_maf(
            &hinge_idx_by_name,
            &mut covered_by_pw[pw_i],
            &mut pairwise_mafs[pw_i].1,
            &rescued_i,
        );
    }
    Ok(())
}

pub(crate) fn mark_uncovered_hinge_contigs(
    hinge_contigs: &ContigList,
    low_cov_species: &[usize],
    sp_to_pw_idx: &HashMap<usize, usize>,
    pairwise_mafs: &[(String, String)],
) -> Vec<bool> {
    let mut need_contig = vec![false; hinge_contigs.len()];
    let ctg_map: HashMap<&str, usize> = hinge_contigs
        .iter()
        .enumerate()
        .map(|(i, (n, _))| (n.as_str(), i))
        .collect();
    for &k in low_cov_species {
        let pw_idx = sp_to_pw_idx[&k];
        let mut intervals: Vec<Vec<(usize, usize)>> = vec![Vec::new(); hinge_contigs.len()];
        for (r0, _) in parse_two_row_maf_blocks(&pairwise_mafs[pw_idx].1) {
            if r0.strand == '+'
                && let Some(&ci) = ctg_map.get(r0.src.as_str())
            {
                let bp = r0.seq.bytes().filter(|&b| b != b'-').count();
                if bp > 0 {
                    intervals[ci].push((r0.start, r0.start.saturating_add(bp)));
                }
            }
        }
        for (ci, (_, seq)) in hinge_contigs.iter().enumerate() {
            let mut ivs = std::mem::take(&mut intervals[ci]);
            ivs.sort_unstable();
            let mut uniq_bp = 0usize;
            let mut cur_end = 0usize;
            for (s, e) in ivs {
                if e > cur_end {
                    uniq_bp += e - s.max(cur_end);
                    cur_end = e;
                }
            }
            if !seq.is_empty() && uniq_bp * 5 < seq.len() * 4 {
                need_contig[ci] = true;
            }
        }
    }
    need_contig
}

pub(crate) fn run_ancestral_consensus_rescue(
    indexed: &[IndexedSpecies],
    hinge_idx: usize,
    sp_to_pw_idx: &HashMap<usize, usize>,
    pairwise_mafs: &mut [(String, String)],
    cfg: &AlignConfig,
    top_k: usize,
    respect_softmask: bool,
) -> Result<()> {
    let hinge_contigs = &indexed[hinge_idx].1;
    let rescue_top_k = if top_k == 0 { 0 } else { top_k.max(3) };
    for round in 0..2 {
        let covs: Vec<(usize, f64)> = (0..indexed.len())
            .into_par_iter()
            .filter(|&i| i != hinge_idx)
            .map(|i| {
                let pw = sp_to_pw_idx[&i];
                (i, maf_ref_coverage(hinge_contigs, &pairwise_mafs[pw].1))
            })
            .collect();
        let seed_mafs: Vec<&str> = covs
            .iter()
            .filter(|&&(_, c)| round > 0 || c >= 0.50)
            .map(|&(i, _)| pairwise_mafs[sp_to_pw_idx[&i]].1.as_str())
            .collect();
        let low_cov_species: Vec<usize> = covs
            .iter()
            .filter(|&&(_, c)| c < 0.50)
            .map(|&(i, _)| i)
            .collect();
        if seed_mafs.is_empty() || low_cov_species.is_empty() {
            break;
        }
        let need_contig = mark_uncovered_hinge_contigs(
            hinge_contigs,
            &low_cov_species,
            sp_to_pw_idx,
            pairwise_mafs,
        );
        let cons_contigs: ContigList =
            build_ancestral_consensus_contigs(hinge_contigs, &seed_mafs, &need_contig)
                .into_iter()
                .enumerate()
                .filter(|&(ci, _)| need_contig[ci])
                .map(|(_, c)| c)
                .collect();
        if cons_contigs.is_empty() {
            break;
        }
        let cons_indices: Vec<(String, ContigIndex)> = cons_contigs
            .par_iter()
            .map(|(n, s)| Ok((n.clone(), prepare_contig_index_opts(s, respect_softmask)?)))
            .collect::<Result<Vec<_>>>()?;
        let remapped_list: Vec<(usize, String)> = low_cov_species
            .par_iter()
            .map(|&k| {
                let (cons_maf, _) = align_species_pair(
                    &cons_contigs,
                    &cons_indices,
                    &indexed[k].1,
                    &indexed[k].2,
                    cfg,
                    rescue_top_k,
                )?;
                let pw_idx = sp_to_pw_idx[&k];
                let remapped = remap_consensus_maf_to_hinge(
                    hinge_contigs,
                    &cons_maf,
                    &pairwise_mafs[pw_idx].1,
                );
                Ok((pw_idx, remapped))
            })
            .collect::<Result<Vec<_>>>()?;

        let mut any_added = false;
        for (pw_idx, remapped) in remapped_list {
            if !remapped.is_empty() {
                pairwise_mafs[pw_idx].1.push_str(&remapped);
                any_added = true;
            }
        }
        if !any_added {
            break;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn upgma_prefers_low_coverage_representative_and_dedups_rescued_maf() {
        let sketches = vec![
            vec![1, 2, 3],
            vec![10, 11, 12, 13, 14],
            vec![10, 11, 12, 13],
            vec![20, 21, 22, 23, 24, 25],
        ];
        let high_cov = vec![false, true, false, true];
        let pairs = build_upgma_rescue_pairs(4, 0, &sketches, &high_cov);
        assert_eq!(pairs.len(), 2);
        assert!(
            pairs.iter().all(|&(a, b)| a == 2 || b == 2),
            "low-coverage species 2 should remain the representative across UPGMA levels: {pairs:?}"
        );

        let mut hinge_idx = HashMap::new();
        hinge_idx.insert("hg38.chr1", 0usize);
        let mut covered = vec![Vec::new()];
        let mut target = String::new();
        let rescued_a = "a score=20.0\ns hg38.chr1 0 20 + 100 ACGTACGTACGTACGTACGT\ns dog.chr1 0 20 + 100 ACGTACGTACGTACGTACGT\n\n";
        let rescued_b = "a score=20.0\ns hg38.chr1 10 20 + 100 ACGTACGTACGTACGTACGT\ns dog.chr1 10 20 + 100 ACGTACGTACGTACGTACGT\n\n";
        append_nonoverlapping_rescued_maf(&hinge_idx, &mut covered, &mut target, rescued_a);
        append_nonoverlapping_rescued_maf(&hinge_idx, &mut covered, &mut target, rescued_b);
        assert_eq!(
            target.matches("\na score=").count() + usize::from(target.starts_with("a score=")),
            1
        );
    }
}
