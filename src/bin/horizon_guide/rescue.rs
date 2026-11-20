//! Transitive sister rescue and ancestral consensus rescue for `horizon_guide`.

use crate::{ContigList, IndexedSpecies, align_species_pair};
use anyhow::Result;
use horizon_aligner::{
    AlignConfig, ContigIndex, complement_base, parse_two_row_maf_blocks, prepare_contig_index_opts,
    shared_seed_count,
};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;

type MidToHingeMap = HashMap<(usize, usize), (usize, usize, u8, bool)>;

fn build_mid_to_hinge_map(
    hinge_idx_by_name: &HashMap<&str, usize>,
    mid_idx_by_name: &HashMap<&str, usize>,
    hinge_to_mid_maf: &str,
) -> MidToHingeMap {
    let mut mid_to_hinge: MidToHingeMap = HashMap::new();
    for (h_row, m_row) in parse_two_row_maf_blocks(hinge_to_mid_maf) {
        let Some(&h_cidx) = hinge_idx_by_name.get(h_row.src.as_str()) else {
            continue;
        };
        let Some(&m_cidx) = mid_idx_by_name.get(m_row.src.as_str()) else {
            continue;
        };
        let mid_rev = (h_row.strand == '-') ^ (m_row.strand == '-');
        let mut hk = 0usize;
        let mut mk = 0usize;
        for (hb, mb) in h_row.seq.bytes().zip(m_row.seq.bytes()) {
            if hb != b'-' && mb != b'-' {
                let hp = h_row.fwd_pos(hk);
                let mp = m_row.fwd_pos(mk);
                let hb_fwd = if h_row.strand == '-' {
                    complement_base(hb)
                } else {
                    hb
                };
                mid_to_hinge.insert((m_cidx, mp), (h_cidx, hp, hb_fwd, mid_rev));
            }
            if hb != b'-' {
                hk += 1;
            }
            if mb != b'-' {
                mk += 1;
            }
        }
    }
    mid_to_hinge
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

    let mid_to_hinge =
        build_mid_to_hinge_map(&hinge_idx_by_name, &mid_idx_by_name, hinge_to_mid_maf);
    if mid_to_hinge.is_empty() {
        return String::new();
    }

    // Track hinge positions already covered by direct hinge <-> dst alignment.
    let mut covered_dst: HashSet<(usize, usize)> = HashSet::new();
    for (h_row, _) in parse_two_row_maf_blocks(dst_direct_maf) {
        let Some(&h_cidx) = hinge_idx_by_name.get(h_row.src.as_str()) else {
            continue;
        };
        let mut hk = 0usize;
        for hb in h_row.seq.bytes() {
            if hb != b'-' {
                covered_dst.insert((h_cidx, h_row.fwd_pos(hk)));
                hk += 1;
            }
        }
    }

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
                if let Some(&(h_cidx, hp, hb, mid_rev)) = mid_to_hinge.get(&(m_cidx, mp))
                    && !covered_dst.contains(&(h_cidx, hp))
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
    pts.sort_unstable_by_key(|p| p.0);
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
    let mut covered: Vec<Vec<bool>> = hinge_contigs
        .iter()
        .map(|(_, s)| vec![false; s.len()])
        .collect();
    for (r0, _) in parse_two_row_maf_blocks(maf_text) {
        if r0.strand != '+' {
            continue;
        }
        if let Some(&ci) = ctg_map.get(r0.src.as_str()) {
            let r0_size = r0.seq.bytes().filter(|&b| b != b'-').count();
            let end = (r0.start + r0_size).min(covered[ci].len());
            for flag in &mut covered[ci][r0.start.min(end)..end] {
                *flag = true;
            }
        }
    }
    let cov_bp: usize = covered
        .iter()
        .map(|v| v.iter().filter(|&&b| b).count())
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
            let mut ov = 0usize;
            for &(a, b) in spans {
                if a < e0 && b > s0 {
                    ov += b.min(e0) - a.max(s0);
                }
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

pub(crate) fn run_progressive_sister_rescue(
    indexed: &[IndexedSpecies],
    hinge_idx: usize,
    sp_to_pw_idx: &HashMap<usize, usize>,
    pairwise_mafs: &mut [(String, String)],
    cfg: &AlignConfig,
    top_k: usize,
) -> Result<()> {
    let hinge_contigs = &indexed[hinge_idx].1;
    let mut aligned_pairs: HashSet<(usize, usize)> = HashSet::new();
    for i in 0..indexed.len() {
        if i == hinge_idx {
            continue;
        }
        let best_j = (0..indexed.len())
            .into_par_iter()
            .filter(|&j| j != i && j != hinge_idx)
            .map(|j| {
                let sc: usize = indexed[i]
                    .2
                    .par_iter()
                    .map(|(_, ix)| {
                        indexed[j]
                            .2
                            .iter()
                            .map(|(_, iy)| shared_seed_count(ix, iy))
                            .sum::<usize>()
                    })
                    .sum();
                (sc, std::cmp::Reverse(j))
            })
            .max()
            .map(|(_, std::cmp::Reverse(j))| j);
        if let Some(j) = best_j {
            let pair_key = (i.min(j), i.max(j));
            if !aligned_pairs.insert(pair_key) {
                continue;
            }
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
            let rescued_j = rescue_transitive_maf(
                hinge_contigs,
                &indexed[i].1,
                &pairwise_mafs[pw_i].1,
                &sister_maf,
                &pairwise_mafs[pw_j].1,
                true,
                12,
                75,
            );
            if !rescued_j.is_empty() {
                pairwise_mafs[pw_j].1.push_str(&rescued_j);
            }
            let rescued_i = rescue_transitive_maf(
                hinge_contigs,
                &indexed[j].1,
                &pairwise_mafs[pw_j].1,
                &sister_maf,
                &pairwise_mafs[pw_i].1,
                false,
                12,
                75,
            );
            if !rescued_i.is_empty() {
                pairwise_mafs[pw_i].1.push_str(&rescued_i);
            }
        }
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
        let mut any_added = false;
        for k in low_cov_species {
            let (cons_maf, _) = align_species_pair(
                &cons_contigs,
                &cons_indices,
                &indexed[k].1,
                &indexed[k].2,
                cfg,
                rescue_top_k,
            )?;
            let pw_idx = sp_to_pw_idx[&k];
            let remapped =
                remap_consensus_maf_to_hinge(hinge_contigs, &cons_maf, &pairwise_mafs[pw_idx].1);
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
