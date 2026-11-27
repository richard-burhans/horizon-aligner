//! Transitive sister rescue and ancestral consensus rescue for `horizon_guide`.

use crate::sister_align::{align_species_pair_sliced, hinge_mapped_contig_spans};
use crate::store::SpeciesStore;
use crate::{ContigRef, IndexedSpecies, contig_refs};
use anyhow::Result;
use horizon_aligner::{
    AlignConfig, ProgressLogger, complement_base, parse_two_row_maf_blocks, sketch_jaccard,
};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::time::Instant;

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
) -> Vec<(Vec<MidHingeRun>, Vec<usize>)> {
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
    by_mid
        .into_iter()
        .map(|mut runs| {
            runs.sort_by_key(|r| r.m_start);
            let mut prefix_max_end = Vec::with_capacity(runs.len());
            let mut max_end = 0usize;
            for r in &runs {
                max_end = max_end.max(r.m_end);
                prefix_max_end.push(max_end);
            }
            (runs, prefix_max_end)
        })
        .collect()
}

fn lookup_mid_to_hinge(
    runs: &[MidHingeRun],
    prefix_max_end: &[usize],
    mp: usize,
    hinge_contigs: &[ContigRef<'_>],
) -> Option<(usize, usize, u8, bool)> {
    let idx = runs.partition_point(|r| r.m_start <= mp);
    for i in (0..idx).rev() {
        if prefix_max_end[i] <= mp {
            break;
        }
        let r = &runs[i];
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
    hinge_contigs: &[ContigRef<'_>],
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
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn rescue_transitive_maf(
    hinge_contigs: &[ContigRef<'_>],
    mid_contigs: &[ContigRef<'_>],
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
        .map(|(idx, &(name, _))| (name, idx))
        .collect();
    let mid_idx_by_name: HashMap<&str, usize> = mid_contigs
        .iter()
        .enumerate()
        .map(|(idx, &(name, _))| (name, idx))
        .collect();

    let mid_to_hinge = build_mid_to_hinge_runs(
        &hinge_idx_by_name,
        &mid_idx_by_name,
        mid_contigs.len(),
        hinge_to_mid_maf,
    );
    if mid_to_hinge.iter().all(|(runs, _)| runs.is_empty()) {
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
        let mut cached_key: Option<(usize, bool)> = None;
        let mut cached_pts: Vec<(usize, usize, u8, u8)> = Vec::new();
        for (mb, db) in m_row.seq.bytes().zip(d_row.seq.bytes()) {
            if mb != b'-' && db != b'-' {
                let mp = m_row.fwd_pos(mk);
                let dp = d_row.fwd_pos(dk);
                if let Some((h_cidx, hp, hb, mid_rev)) = lookup_mid_to_hinge(
                    &mid_to_hinge[m_cidx].0,
                    &mid_to_hinge[m_cidx].1,
                    mp,
                    hinge_contigs,
                ) && !is_bitset_covered(&covered_dst[h_cidx], hp)
                {
                    let dst_rel_rev = mid_rev ^ (m_row.strand == '-') ^ (d_row.strand == '-');
                    let db_oriented = if (d_row.strand == '-') == dst_rel_rev {
                        db
                    } else {
                        complement_base(db)
                    };
                    if cached_key != Some((h_cidx, dst_rel_rev))
                        && let Some((prev_h, prev_rev)) = cached_key.replace((h_cidx, dst_rel_rev))
                        && !cached_pts.is_empty()
                    {
                        projected
                            .entry((prev_h, d_row.src.clone(), d_row.src_size, prev_rev))
                            .or_default()
                            .append(&mut cached_pts);
                    }
                    cached_pts.push((hp, dp, hb, db_oriented));
                }
            }
            if mb != b'-' {
                mk += 1;
            }
            if db != b'-' {
                dk += 1;
            }
        }
        if let Some((prev_h, prev_rev)) = cached_key
            && !cached_pts.is_empty()
        {
            projected
                .entry((prev_h, d_row.src.clone(), d_row.src_size, prev_rev))
                .or_default()
                .append(&mut cached_pts);
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

pub(crate) fn maf_ref_coverage(hinge_contigs: &[ContigRef<'_>], maf_text: &str) -> f64 {
    let total_bp: usize = hinge_contigs.iter().map(|(_, s)| s.len()).sum();
    if total_bp == 0 {
        return 0.0;
    }
    let ctg_map: HashMap<&str, usize> = hinge_contigs
        .iter()
        .enumerate()
        .map(|(i, &(n, _))| (n, i))
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

fn build_upgma_rescue_pairs(
    n_species: usize,
    hinge_idx: usize,
    sketches: &[Vec<u64>],
    high_cov: &[bool],
) -> Vec<(usize, usize)> {
    let mut sim_mat = vec![vec![0.0_f64; n_species]; n_species];
    for i in 0..n_species {
        sim_mat[i][i] = 1.0;
        for j in (i + 1)..n_species {
            let s = sketch_jaccard(&sketches[i], &sketches[j]);
            sim_mat[i][j] = s;
            sim_mat[j][i] = s;
        }
    }

    let mut clades: Vec<(usize, Vec<usize>, bool)> = (0..n_species)
        .filter(|&i| i != hinge_idx)
        .map(|i| (i, vec![i], high_cov[i]))
        .collect();
    let mut pairs = Vec::new();
    while clades.len() >= 2 {
        let mut best: Option<(f64, usize, usize)> = None;
        for a in 0..clades.len() {
            for b in (a + 1)..clades.len() {
                let mut sum = 0.0_f64;
                let mut cnt = 0usize;
                for &i in &clades[a].1 {
                    for &j in &clades[b].1 {
                        sum += sim_mat[i][j];
                        cnt += 1;
                    }
                }
                #[allow(clippy::cast_precision_loss)]
                let sim = if cnt == 0 { 0.0 } else { sum / (cnt as f64) };
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

pub(crate) fn append_nonoverlapping_rescued_maf(
    hinge_idx_by_name: &HashMap<&str, usize>,
    covered: &mut [Vec<(usize, usize)>],
    target_maf: &mut String,
    rescued_maf: &str,
) -> usize {
    if rescued_maf.is_empty() {
        return 0;
    }
    let mut added = 0usize;
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
        added += 1;
    }
    added
}

/// Order UPGMA rescue pairs outward by BFS from already high-coverage species (`high_cov`) and
/// partition them into vertex-disjoint dependency waves. Bridge pairs linking a reached species
/// to an unreached sister commit rescued blocks before downstream multi-hop pairs read that
/// sister's MAF — working identically for `(A_hi, (B_lo, C_lo))`, `((A_hi, B_lo), C_lo)`, and
/// `((A_hi, B_lo), (C_lo, D_lo))` topologies.
fn build_upgma_rescue_waves(
    pairs_bottom_up: &[(usize, usize)],
    high_cov: &[bool],
) -> Vec<Vec<(usize, usize)>> {
    if pairs_bottom_up.is_empty() {
        return Vec::new();
    }
    let max_node = pairs_bottom_up
        .iter()
        .flat_map(|&(u, v)| [u, v])
        .max()
        .unwrap_or(0);
    let mut reached = vec![false; (max_node + 1).max(high_cov.len())];
    for (i, &hc) in high_cov.iter().enumerate() {
        reached[i] = hc;
    }

    let mut remaining: Vec<(usize, usize)> = pairs_bottom_up.to_vec();
    let mut ordered: Vec<(usize, usize)> = Vec::with_capacity(pairs_bottom_up.len());
    while !remaining.is_empty() {
        let mut frontier = Vec::new();
        let mut next_rem = Vec::new();
        for &(u, v) in remaining.iter().rev() {
            if reached[u] || reached[v] {
                frontier.push((u, v));
            } else {
                next_rem.push((u, v));
            }
        }
        next_rem.reverse();
        if frontier.is_empty() {
            // Disconnected component with no high-coverage seed: fall back to top-down UPGMA order.
            if let Some(edge) = next_rem.pop() {
                reached[edge.0] = true;
                reached[edge.1] = true;
                ordered.push(edge);
            }
        } else {
            for &(u, v) in &frontier {
                reached[u] = true;
                reached[v] = true;
            }
            ordered.extend(frontier);
        }
        remaining = next_rem;
    }

    let mut pair_wave = vec![0usize; ordered.len()];
    let mut max_wave = 0usize;
    for k in 0..ordered.len() {
        let (u, v) = ordered[k];
        let mut w = 0usize;
        for prev in 0..k {
            let (pu, pv) = ordered[prev];
            if pu == u || pu == v || pv == u || pv == v {
                w = w.max(pair_wave[prev] + 1);
            }
        }
        pair_wave[k] = w;
        max_wave = max_wave.max(w);
    }
    let mut waves = vec![Vec::new(); max_wave + 1];
    for (pair, w) in ordered.into_iter().zip(pair_wave) {
        waves[w].push(pair);
    }
    waves
}

#[allow(clippy::too_many_arguments)]
fn rescue_sister_pair_refs(
    sp_i: &IndexedSpecies,
    sp_j: &IndexedSpecies,
    hinge_contigs: &[ContigRef<'_>],
    pw_i: usize,
    pw_j: usize,
    pairwise_mafs: &[(String, String)],
    cfg: &AlignConfig,
    top_k: usize,
) -> Result<(usize, String, usize, String, usize)> {
    let max_mapped = if top_k == 0 { 0 } else { top_k.clamp(4, 20) };
    let spans_i = hinge_mapped_contig_spans(&sp_i.1, &pairwise_mafs[pw_i].1, max_mapped);
    let spans_j = hinge_mapped_contig_spans(&sp_j.1, &pairwise_mafs[pw_j].1, max_mapped);
    if spans_i.is_empty() && spans_j.is_empty() {
        return Ok((pw_j, String::new(), pw_i, String::new(), 0));
    }
    let (sister_maf, n_pairs) =
        align_species_pair_sliced(sp_i, &spans_i, sp_j, &spans_j, cfg, top_k)?;
    let mid_i = contig_refs(&sp_i.1);
    let mid_j = contig_refs(&sp_j.1);
    let (rescued_j, rescued_i) = rayon::join(
        || {
            if spans_i.is_empty() {
                String::new()
            } else {
                rescue_transitive_maf(
                    hinge_contigs,
                    &mid_i,
                    &pairwise_mafs[pw_i].1,
                    &sister_maf,
                    &pairwise_mafs[pw_j].1,
                    true,
                    12,
                    75,
                )
            }
        },
        || {
            if spans_j.is_empty() {
                String::new()
            } else {
                rescue_transitive_maf(
                    hinge_contigs,
                    &mid_j,
                    &pairwise_mafs[pw_j].1,
                    &sister_maf,
                    &pairwise_mafs[pw_i].1,
                    false,
                    12,
                    75,
                )
            }
        },
    );
    Ok((pw_j, rescued_j, pw_i, rescued_i, n_pairs))
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn run_progressive_sister_rescue(
    store: &SpeciesStore,
    sketches: &[Vec<u64>],
    hinge_idx: usize,
    sp_to_pw_idx: &HashMap<usize, usize>,
    pairwise_mafs: &mut [(String, String)],
    cfg: &AlignConfig,
    top_k: usize,
    logger: ProgressLogger,
) -> Result<()> {
    let stage_t0 = Instant::now();
    let hinge_contigs = contig_refs(&store.hinge(hinge_idx).1);
    let hinge_idx_by_name: HashMap<&str, usize> = hinge_contigs
        .iter()
        .enumerate()
        .map(|(idx, &(name, _))| (name, idx))
        .collect();
    let high_cov: Vec<bool> = (0..store.len())
        .into_par_iter()
        .map(|i| {
            if i == hinge_idx {
                false
            } else {
                let pw = sp_to_pw_idx[&i];
                maf_ref_coverage(&hinge_contigs, &pairwise_mafs[pw].1) >= 0.85
            }
        })
        .collect();
    let sister_pairs = build_upgma_rescue_pairs(store.len(), hinge_idx, sketches, &high_cov);
    if sister_pairs.is_empty() {
        logger.log("sister rescue: skipped (all species already >=85% covered)");
        return Ok(());
    }
    let rescue_waves = build_upgma_rescue_waves(&sister_pairs, &high_cov);
    let n_waves = rescue_waves.len();
    logger.log(format_args!(
        "sister rescue: starting {} UPGMA pairs across {n_waves} waves",
        sister_pairs.len(),
    ));

    let mut covered_by_pw: Vec<Vec<Vec<(usize, usize)>>> =
        vec![vec![Vec::new(); hinge_contigs.len()]; pairwise_mafs.len()];
    let mut rescued_blocks = 0usize;
    for (wave_idx, wave) in rescue_waves.into_iter().enumerate() {
        let wave_t0 = Instant::now();
        let rescued_chunks: Vec<(usize, String, usize, String, usize)> =
            if let Some(indexed) = store.in_memory() {
                wave.par_iter()
                    .map(|&(i, j)| {
                        rescue_sister_pair_refs(
                            &indexed[i],
                            &indexed[j],
                            &hinge_contigs,
                            sp_to_pw_idx[&i],
                            sp_to_pw_idx[&j],
                            pairwise_mafs,
                            cfg,
                            top_k,
                        )
                    })
                    .collect::<Result<Vec<_>>>()?
            } else {
                // Streaming mode: load only the two species in each pair sequentially (each indexed
                // & aligned across all Rayon threads) so peak RAM stays O(2 species), never O(N).
                let mut chunks = Vec::with_capacity(wave.len());
                for &(i, j) in &wave {
                    let pw_i = sp_to_pw_idx[&i];
                    let pw_j = sp_to_pw_idx[&j];
                    let has_i = pairwise_mafs[pw_i].1.contains("\na ")
                        || pairwise_mafs[pw_i].1.starts_with("a ");
                    let has_j = pairwise_mafs[pw_j].1.contains("\na ")
                        || pairwise_mafs[pw_j].1.starts_with("a ");
                    if !has_i && !has_j {
                        chunks.push((pw_j, String::new(), pw_i, String::new(), 0));
                        continue;
                    }
                    let sp_i = store.load_streaming_species(i)?;
                    let sp_j = store.load_streaming_species(j)?;
                    chunks.push(rescue_sister_pair_refs(
                        &sp_i,
                        &sp_j,
                        &hinge_contigs,
                        pw_i,
                        pw_j,
                        pairwise_mafs,
                        cfg,
                        top_k,
                    )?);
                }
                chunks
            };

        let mut wave_blocks = 0usize;
        let mut wave_pairs = 0usize;
        for (pw_j, rescued_j, pw_i, rescued_i, n_pairs) in rescued_chunks {
            wave_pairs += n_pairs;
            wave_blocks += append_nonoverlapping_rescued_maf(
                &hinge_idx_by_name,
                &mut covered_by_pw[pw_j],
                &mut pairwise_mafs[pw_j].1,
                &rescued_j,
            );
            wave_blocks += append_nonoverlapping_rescued_maf(
                &hinge_idx_by_name,
                &mut covered_by_pw[pw_i],
                &mut pairwise_mafs[pw_i].1,
                &rescued_i,
            );
        }
        rescued_blocks += wave_blocks;
        logger.log(format_args!(
            "sister rescue: wave {}/{n_waves} ({} sp-pairs, {wave_pairs} contig pairs) -> +{wave_blocks} blocks ({:.1}s)",
            wave_idx + 1,
            wave.len(),
            wave_t0.elapsed().as_secs_f64()
        ));
    }
    logger.log(format_args!(
        "sister rescue: done ({rescued_blocks} blocks rescued in {:.1}s)",
        stage_t0.elapsed().as_secs_f64()
    ));
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

    #[test]
    fn upgma_rescue_averages_all_pairs_across_clade_members() {
        // Species 0 is hinge. Non-hinge: 1, 2, 3, 4.
        // 1 and 2 merge first (high overlap). Species 1 (rep after merge since low_cov=false
        // for 1 and true for 2) shares 0 hashes with 3 and 1 hash with 4, whereas species 2
        // shares 5 hashes with 3 and 0 hashes with 4.
        // Single-rep (using only rep=1) would merge {1,2} with 4; true all-pairs UPGMA averages
        // sim(1,3)+sim(2,3) > sim(1,4)+sim(2,4) and merges {1,2} with 3 before 4.
        let sketches = vec![
            vec![999],
            vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 50],
            vec![1, 2, 3, 4, 5, 6, 7, 8, 20, 21, 22, 23, 24],
            vec![20, 21, 22, 23, 24, 30, 31, 32, 33, 34],
            vec![50, 60, 61, 62, 63, 64, 65, 66, 67, 68],
        ];
        let high_cov = vec![true, false, true, false, false];
        let pairs = build_upgma_rescue_pairs(5, 0, &sketches, &high_cov);
        assert_eq!(pairs.len(), 3);
        assert!(
            (pairs[0] == (1, 2) || pairs[0] == (2, 1)),
            "first merge must join sister pair 1 & 2: {pairs:?}"
        );
        assert!(
            (pairs[1] == (1, 3) || pairs[1] == (3, 1)),
            "second merge must use all-pairs clade average to join {{1,2}} with 3 before 4: {pairs:?}"
        );
    }

    #[test]
    fn upgma_rescue_waves_schedule_bridge_before_dependent_clade_pair() {
        // Topology 1: (A_hi=1, (B_lo=2, C_lo=3)) -> bottom-up [(2, 3), (2, 1)].
        // Wave 0 must run bridge (2, 1) before Wave 1 runs (2, 3).
        let hc_1 = [false, true, false, false];
        let waves = build_upgma_rescue_waves(&[(2usize, 3usize), (2usize, 1usize)], &hc_1);
        assert_eq!(waves, vec![vec![(2, 1)], vec![(2, 3)]]);

        // Topology 2: ((A_hi=1, B_lo=2), C_lo=3) -> bottom-up [(2, 1), (2, 3)].
        // Must STILL schedule (2, 1) in Wave 0 and (2, 3) in Wave 1 (simple reverse breaks this!).
        let waves2 = build_upgma_rescue_waves(&[(2usize, 1usize), (2usize, 3usize)], &hc_1);
        assert_eq!(waves2, vec![vec![(2, 1)], vec![(2, 3)]]);

        // Topology 3: ((A_hi=1, B_lo=2), (C_lo=3, D_lo=4)) -> bottom-up [(2, 1), (3, 4), (2, 3)].
        // Must propagate A -> B -> C -> D across 3 sequential waves.
        let hc_3 = [false, true, false, false, false];
        let waves3 = build_upgma_rescue_waves(&[(2, 1), (3, 4), (2, 3)], &hc_3);
        assert_eq!(waves3, vec![vec![(2, 1)], vec![(2, 3)], vec![(3, 4)]]);

        // Disjoint bridge pairs (1, 2) and (3, 4) with high_cov {1, 3} parallelize in wave 0.
        let hc_disjoint = [false, true, false, true, false];
        let disjoint = build_upgma_rescue_waves(&[(1, 2), (3, 4)], &hc_disjoint);
        assert_eq!(disjoint.len(), 1);
        assert_eq!(disjoint[0].len(), 2);
    }
}
