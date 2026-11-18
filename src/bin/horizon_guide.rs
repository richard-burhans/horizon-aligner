//! Nystrom/sketch-selected hinge + multi-species WGA pipeline (`horizon_wga` + `horizon_merge`).

use anyhow::{Context, Result, bail};
use horizon_aligner::{
    AlignConfig, ContigIndex, align_contig_pair_with_config, format_maf_blocks, hal::HalAlignment,
    merge_pairwise_mafs, prepare_contig_index_opts, shared_seed_count,
};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

type ContigList = Vec<(String, Vec<u8>)>;
type IndexedSpecies = (String, ContigList, Vec<(String, ContigIndex)>);

fn read_fasta_contigs(path: &Path, max_contigs: usize) -> Result<ContigList> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut out = Vec::new();
    let mut cur_name: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = cur_name.take() {
                out.push((name, std::mem::take(&mut cur_seq)));
                if out.len() >= max_contigs {
                    return Ok(out);
                }
            }
            let header = rest
                .split_whitespace()
                .next()
                .unwrap_or("contig")
                .to_owned();
            cur_name = Some(header);
        } else if cur_name.is_some() {
            cur_seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some(name) = cur_name {
        out.push((name, cur_seq));
    }
    Ok(out)
}

fn select_central_hinge(species: &[IndexedSpecies]) -> usize {
    let mut best_idx = 0usize;
    let mut best_score = 0usize;
    for (i, (_, _, idx_i)) in species.iter().enumerate() {
        let mut total = 0usize;
        for (j, (_, _, idx_j)) in species.iter().enumerate() {
            if i == j {
                continue;
            }
            for (_, ix) in idx_i {
                for (_, iy) in idx_j {
                    total += shared_seed_count(ix, iy);
                }
            }
        }
        if total > best_score {
            best_score = total;
            best_idx = i;
        }
    }
    best_idx
}

type HalBlockList = Vec<(usize, usize, usize, bool)>;

fn align_species_pair(
    ref_contigs: &ContigList,
    ref_indices: &[(String, ContigIndex)],
    qry_contigs: &ContigList,
    qry_indices: &[(String, ContigIndex)],
    cfg: &AlignConfig,
) -> Result<(String, HalBlockList)> {
    let mut maf_buf = String::from("##maf version=1 scoring=horizon_mdl\n");
    let mut hal_blocks = Vec::new();
    for (xi, (nx, ix)) in ref_indices.iter().enumerate() {
        for (yi, (ny, iy)) in qry_indices.iter().enumerate() {
            if shared_seed_count(ix, iy) < 8 {
                continue;
            }
            let aln = align_contig_pair_with_config(ix, iy, cfg)?;
            if aln.aligned_pairs.is_empty() {
                continue;
            }
            let sx = &ref_contigs[xi].1;
            let sy = &qry_contigs[yi].1;
            maf_buf.push_str(&format_maf_blocks(nx, sx.len(), ny, sy.len(), sx, sy, &aln));
            if let Some(&(x0, y0)) = aln.aligned_pairs.first()
                && let Some(&(x1, y1)) = aln.aligned_pairs.last()
                && x1 > x0
            {
                hal_blocks.push((x0, y0.min(y1), x1 - x0 + 1, y1 < y0));
            }
        }
    }
    Ok((maf_buf, hal_blocks))
}

struct MafRow {
    src: String,
    start: usize,
    strand: char,
    src_size: usize,
    seq: String,
}

fn parse_two_row_maf_blocks(maf_text: &str) -> Vec<(MafRow, MafRow)> {
    let mut blocks = Vec::new();
    let mut rows: Vec<MafRow> = Vec::new();
    for line in maf_text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            if rows.len() >= 2 {
                let r1 = rows.remove(0);
                let r2 = rows.remove(0);
                blocks.push((r1, r2));
            }
            rows.clear();
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix("s ") {
            let parts: Vec<&str> = rest.split_whitespace().collect();
            if parts.len() >= 6
                && let (Ok(start), Ok(src_size)) =
                    (parts[1].parse::<usize>(), parts[4].parse::<usize>())
            {
                let strand = parts[3].chars().next().unwrap_or('+');
                rows.push(MafRow {
                    src: parts[0].to_owned(),
                    start,
                    strand,
                    src_size,
                    seq: parts[5].to_owned(),
                });
            }
        }
    }
    if rows.len() >= 2 {
        let r1 = rows.remove(0);
        let r2 = rows.remove(0);
        blocks.push((r1, r2));
    }
    blocks
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
        if run.len() >= min_run_len && matches * 100 >= run.len() * min_ident_pct {
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

#[allow(clippy::too_many_lines)]
fn main() -> Result<()> {
    let raw_args: Vec<String> = std::env::args().skip(1).collect();
    let mut progressive = false;
    let mut respect_softmask = false;
    let mut pos = Vec::new();
    for a in raw_args {
        if a == "--progressive" {
            progressive = true;
        } else if a == "--respect-softmask" {
            respect_softmask = true;
        } else {
            pos.push(a);
        }
    }
    if pos.len() < 4 {
        bail!(
            "usage: horizon_guide [--progressive] [--respect-softmask] <out_prefix> <max_contigs> <sp1=fa1> <sp2=fa2>..."
        );
    }
    let t0 = Instant::now();
    let out_prefix = &pos[0];
    let max_contigs: usize = pos[1].parse().context("parsing max_contigs")?;

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
                .iter()
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
        let (maf_buf, hblocks) =
            align_species_pair(hinge_contigs, hinge_indices, qry_contigs, qry_indices, &cfg)?;
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
                .filter(|&j| j != i && j != hinge_idx)
                .max_by_key(|&j| {
                    indexed[i]
                        .2
                        .iter()
                        .flat_map(|(_, ix)| {
                            indexed[j]
                                .2
                                .iter()
                                .map(move |(_, iy)| shared_seed_count(ix, iy))
                        })
                        .sum::<usize>()
                });
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
                    60,
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
                    60,
                );
                if !rescued_i.is_empty() {
                    rescues.push((pw_i, rescued_i));
                }
            }
        }
        for (pw_idx, rescued_maf) in rescues {
            pairwise_mafs[pw_idx].1.push_str(&rescued_maf);
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
