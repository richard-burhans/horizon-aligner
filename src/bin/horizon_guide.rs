//! Nystrom/sketch-selected hinge + multi-species WGA pipeline (`horizon_wga` + `horizon_merge`).

use anyhow::{Context, Result, bail};
use horizon_aligner::{
    AlignConfig, ContigIndex, MafRow, align_contig_pair_with_config, format_maf_blocks,
    hal::HalAlignment, merge_pairwise_mafs, parse_two_row_maf_blocks, prepare_contig_index_opts,
    read_fasta_contigs, select_top_k_pairs, shared_seed_count,
};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

type ContigList = Vec<(String, Vec<u8>)>;
type IndexedSpecies = (String, ContigList, Vec<(String, ContigIndex)>);

fn select_central_hinge(species: &[IndexedSpecies]) -> usize {
    (0..species.len())
        .into_par_iter()
        .map(|i| {
            let idx_i = &species[i].2;
            let total: usize = (0..species.len())
                .into_par_iter()
                .filter(|&j| j != i)
                .map(|j| {
                    let idx_j = &species[j].2;
                    idx_i
                        .par_iter()
                        .map(|(_, ix)| {
                            idx_j
                                .iter()
                                .map(|(_, iy)| shared_seed_count(ix, iy))
                                .sum::<usize>()
                        })
                        .sum::<usize>()
                })
                .sum();
            (total, std::cmp::Reverse(i))
        })
        .max()
        .map_or(0, |(_, std::cmp::Reverse(i))| i)
}

type HalBlock = (usize, usize, usize, bool);
type HalBlockList = Vec<HalBlock>;

fn align_species_pair(
    ref_contigs: &ContigList,
    ref_indices: &[(String, ContigIndex)],
    qry_contigs: &ContigList,
    qry_indices: &[(String, ContigIndex)],
    cfg: &AlignConfig,
    top_k: usize,
) -> Result<(String, HalBlockList)> {
    let min_shared_seeds = 8usize;
    let top_mask = if top_k > 0 {
        Some(select_top_k_pairs(
            ref_indices,
            qry_indices,
            top_k,
            min_shared_seeds,
        ))
    } else {
        None
    };

    let pairs: Vec<(usize, usize)> = (0..ref_indices.len())
        .flat_map(|xi| (0..qry_indices.len()).map(move |yi| (xi, yi)))
        .filter(|&(xi, yi)| {
            if let Some(mask) = &top_mask {
                mask[xi][yi]
            } else {
                shared_seed_count(&ref_indices[xi].1, &qry_indices[yi].1) >= min_shared_seeds
            }
        })
        .collect();

    let results: Vec<(String, Option<HalBlock>)> = pairs
        .par_iter()
        .map(|&(xi, yi)| {
            let (nx, ix) = &ref_indices[xi];
            let (ny, iy) = &qry_indices[yi];
            let aln = align_contig_pair_with_config(ix, iy, cfg)?;
            if aln.aligned_pairs.is_empty() {
                return Ok((String::new(), None));
            }
            let sx = &ref_contigs[xi].1;
            let sy = &qry_contigs[yi].1;
            let maf = format_maf_blocks(nx, sx.len(), ny, sy.len(), sx, sy, &aln);
            let hblk = if let Some(&(x0, y0)) = aln.aligned_pairs.first()
                && let Some(&(x1, y1)) = aln.aligned_pairs.last()
                && x1 > x0
            {
                Some((x0, y0.min(y1), x1 - x0 + 1, y1 < y0))
            } else {
                None
            };
            Ok((maf, hblk))
        })
        .collect::<Result<Vec<_>>>()?;

    let mut maf_buf = String::from("##maf version=1 scoring=horizon_mdl\n");
    let mut hal_blocks = Vec::new();
    for (maf, hblk) in results {
        if !maf.is_empty() {
            maf_buf.push_str(&maf);
        }
        if let Some(b) = hblk {
            hal_blocks.push(b);
        }
    }
    Ok((maf_buf, hal_blocks))
}

fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        other => other,
    }
}

fn fwd_pos(row: &MafRow, nongap_k: usize) -> usize {
    if row.strand == '-' {
        row.src_size - 1 - (row.start + nongap_k)
    } else {
        row.start + nongap_k
    }
}

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
                let hp = fwd_pos(&h_row, hk);
                let mp = fwd_pos(&m_row, mk);
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
fn rescue_transitive_maf(
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
                covered_dst.insert((h_cidx, fwd_pos(&h_row, hk)));
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
                let mp = fwd_pos(&m_row, mk);
                let dp = fwd_pos(&d_row, dk);
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
                d_src_size - 1 - d_max
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

fn build_ancestral_consensus_contigs(
    hinge_contigs: &ContigList,
    high_cov_mafs: &[&str],
) -> ContigList {
    const BASES: [u8; 4] = *b"ACGT";
    let ctg_map: HashMap<&str, usize> = hinge_contigs
        .iter()
        .enumerate()
        .map(|(i, (n, _))| (n.as_str(), i))
        .collect();
    let mut counts: Vec<Vec<[u32; 4]>> = hinge_contigs
        .iter()
        .map(|(_, seq)| {
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
                let mut best_b = rb;
                let mut best_cnt = 0u32;
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

fn remap_consensus_maf_to_hinge(hinge_contigs: &ContigList, cons_maf: &str) -> String {
    let ctg_map: HashMap<&str, &[u8]> = hinge_contigs
        .iter()
        .map(|(n, s)| (n.as_str(), s.as_slice()))
        .collect();
    let mut out = String::new();
    for (r0, r1) in parse_two_row_maf_blocks(cons_maf) {
        let Some(&ref_seq) = ctg_map.get(r0.src.as_str()) else {
            continue;
        };
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
        let r0_size = r0.seq.bytes().filter(|&b| b != b'-').count();
        let r1_size = r1.seq.bytes().filter(|&b| b != b'-').count();
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

#[allow(clippy::too_many_lines)]
fn main() -> Result<()> {
    let raw_args: Vec<String> = std::env::args().skip(1).collect();
    let mut progressive = false;
    let mut respect_softmask = false;
    let mut write_hal = false;
    let mut top_k_opt: Option<usize> = None;
    let mut pos = Vec::new();
    let mut arg_i = 0usize;
    while arg_i < raw_args.len() {
        match raw_args[arg_i].as_str() {
            "--progressive" => {
                progressive = true;
                arg_i += 1;
            }
            "--respect-softmask" => {
                respect_softmask = true;
                arg_i += 1;
            }
            "--hal" => {
                write_hal = true;
                arg_i += 1;
            }
            "--rbest" => {
                top_k_opt = Some(1);
                arg_i += 1;
            }
            "--all-pairs" => {
                top_k_opt = Some(0);
                arg_i += 1;
            }
            "--top-k" if arg_i + 1 < raw_args.len() => {
                if let Ok(v) = raw_args[arg_i + 1].parse() {
                    top_k_opt = Some(v);
                }
                arg_i += 2;
            }
            other => {
                pos.push(other.to_owned());
                arg_i += 1;
            }
        }
    }
    if pos.len() < 4 {
        bail!(
            "usage: horizon_guide [--progressive] [--respect-softmask] [--hal] [--rbest | --top-k K | --all-pairs] <out_prefix> <max_contigs> <sp1=fa1> <sp2=fa2>..."
        );
    }
    let t0 = Instant::now();
    let out_prefix = &pos[0];
    let max_contigs: usize = pos[1].parse().context("parsing max_contigs")?;
    let top_k = top_k_opt.unwrap_or(usize::from(max_contigs != 1));

    let mut specs: Vec<(String, PathBuf)> = Vec::new();
    for spec in &pos[2..] {
        let Some((sp, path)) = spec.split_once('=') else {
            bail!("expected <species>=<fasta>, got {spec}");
        };
        specs.push((sp.to_owned(), PathBuf::from(path)));
    }

    let indexed: Vec<IndexedSpecies> = specs
        .par_iter()
        .map(|(sp, path)| {
            let contigs = read_fasta_contigs(path, max_contigs)?;
            let idx = contigs
                .par_iter()
                .map(|(n, s)| Ok((n.clone(), prepare_contig_index_opts(s, respect_softmask)?)))
                .collect::<Result<Vec<_>>>()?;
            Ok((sp.clone(), contigs, idx))
        })
        .collect::<Result<Vec<_>>>()?;

    let hinge_idx = select_central_hinge(&indexed);
    let (hinge_sp, hinge_contigs, hinge_indices) = &indexed[hinge_idx];
    let cfg = AlignConfig {
        skip_softmasked: respect_softmask,
        ..AlignConfig::default()
    };

    let mut pairwise_mafs: Vec<(String, String)> = Vec::new();
    let mut hal_queries: Vec<(String, Vec<u8>)> = Vec::new();
    let mut hal_all_blocks: Vec<(usize, usize, usize, usize, bool)> = Vec::new();

    let mut sp_to_pw_idx: HashMap<usize, usize> = HashMap::new();
    for (i, (qry_sp, qry_contigs, qry_indices)) in indexed.iter().enumerate() {
        if i == hinge_idx {
            continue;
        }
        let (maf_buf, hblocks) = align_species_pair(
            hinge_contigs,
            hinge_indices,
            qry_contigs,
            qry_indices,
            &cfg,
            top_k,
        )?;
        let q_idx = hal_queries.len();
        hal_queries.push((
            qry_sp.clone(),
            qry_contigs.first().map_or_else(Vec::new, |c| c.1.clone()),
        ));
        for (rs, qs, len, rev) in hblocks {
            hal_all_blocks.push((q_idx, rs, qs, len, rev));
        }
        sp_to_pw_idx.insert(i, pairwise_mafs.len());
        pairwise_mafs.push((qry_sp.clone(), maf_buf));
    }

    if progressive && indexed.len() >= 3 {
        // Align non-hinge sister pairs and project transitive hinge <-> i <-> j alignments.
        let mut aligned_pairs: HashSet<(usize, usize)> = HashSet::new();
        let mut rescues: Vec<(usize, String)> = Vec::new();
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
                    (sc, j)
                })
                .max()
                .map(|(_, j)| j);
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
                    &cfg,
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
                    rescues.push((pw_j, rescued_j));
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
                    rescues.push((pw_i, rescued_i));
                }
            }
        }
        for (pw_idx, rescued_maf) in rescues {
            pairwise_mafs[pw_idx].1.push_str(&rescued_maf);
        }

        // Ancestral consensus rescue for distant outgroups (< 50% hinge coverage).
        let covs: Vec<(usize, f64)> = (0..indexed.len())
            .filter(|&i| i != hinge_idx)
            .map(|i| {
                let pw = sp_to_pw_idx[&i];
                (i, maf_ref_coverage(hinge_contigs, &pairwise_mafs[pw].1))
            })
            .collect();
        let high_cov_mafs: Vec<&str> = covs
            .iter()
            .filter(|&&(_, c)| c >= 0.50)
            .map(|&(i, _)| pairwise_mafs[sp_to_pw_idx[&i]].1.as_str())
            .collect();
        let low_cov_species: Vec<usize> = covs
            .iter()
            .filter(|&&(_, c)| c < 0.50)
            .map(|&(i, _)| i)
            .collect();
        if !high_cov_mafs.is_empty() && !low_cov_species.is_empty() {
            let cons_contigs = build_ancestral_consensus_contigs(hinge_contigs, &high_cov_mafs);
            let cons_indices: Vec<(String, ContigIndex)> = cons_contigs
                .par_iter()
                .map(|(n, s)| Ok((n.clone(), prepare_contig_index_opts(s, respect_softmask)?)))
                .collect::<Result<Vec<_>>>()?;
            for k in low_cov_species {
                let (cons_maf, _) = align_species_pair(
                    &cons_contigs,
                    &cons_indices,
                    &indexed[k].1,
                    &indexed[k].2,
                    &cfg,
                    top_k,
                )?;
                let remapped = remap_consensus_maf_to_hinge(hinge_contigs, &cons_maf);
                if !remapped.is_empty() {
                    pairwise_mafs[sp_to_pw_idx[&k]].1.push_str(&remapped);
                }
            }
        }
    }

    let hal = HalAlignment::from_star_blocks(
        hinge_sp,
        hinge_contigs.first().map_or(&[], |c| c.1.as_slice()),
        &hal_queries,
        &hal_all_blocks,
    );
    let (merged_maf, merged_gfa, total_blocks, multi_blocks) =
        merge_pairwise_mafs(hinge_sp, &pairwise_mafs, 1000);
    fs::write(format!("{out_prefix}.maf"), merged_maf)?;
    fs::write(format!("{out_prefix}.gfa"), merged_gfa)?;
    if write_hal {
        hal.write_hdf5_hal(std::path::Path::new(&format!("{out_prefix}.hal")))?;
    }
    println!(
        "hinge={} species={} hal_tree={} hal_bottom_segs={} total_blocks={} multi_species_blocks={} elapsed_s={:.4}",
        hinge_sp,
        indexed.len(),
        hal.tree_newick,
        hal.genomes[0].bottom_segments.len(),
        total_blocks,
        multi_blocks,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}
