//! Canonical HAL segment sweep-line partitioning, gap filling, and paralogy-cycle normalization.

use super::{HalAlignment, HalBottomSegment, HalSequence, HalTopSegment};

fn interval_overlaps(
    occupied: &std::collections::BTreeMap<usize, usize>,
    start: usize,
    end: usize,
) -> bool {
    if let Some((&_prev_s, &prev_e)) = occupied.range(..end).next_back()
        && prev_e > start
    {
        return true;
    }
    occupied.range(start..end).next().is_some()
}

fn seq_contains_interval(seqs: &[HalSequence], start: usize, end: usize) -> bool {
    let idx = seqs.partition_point(|s| s.offset.saturating_add(s.length) <= start);
    seqs.get(idx)
        .is_some_and(|s| start >= s.offset && end <= s.offset.saturating_add(s.length))
}

fn finalize_child_top_segments(
    sequences: &[HalSequence],
    mut segs: Vec<(usize, usize, usize, bool)>,
    c_idx: usize,
    root_bottoms: &mut [HalBottomSegment],
) -> Vec<HalTopSegment> {
    segs.sort_by_key(|&(s, _, _, _)| s);

    let mut final_tops: Vec<HalTopSegment> = Vec::new();
    let mut seg_idx = 0usize;
    for seq in sequences {
        let mut pos = seq.offset;
        let seq_end = seq.offset.saturating_add(seq.length);
        while seg_idx < segs.len() && segs[seg_idx].0 < seq_end {
            let (s, l, p_b_idx, rev) = segs[seg_idx];
            seg_idx += 1;
            if s < pos {
                continue;
            }
            if s > pos {
                final_tops.push(HalTopSegment {
                    start: pos,
                    length: s - pos,
                    parent_bottom_idx: None,
                    parent_reversed: false,
                    next_paralog: None,
                });
            }
            final_tops.push(HalTopSegment {
                start: s,
                length: l,
                parent_bottom_idx: Some(p_b_idx),
                parent_reversed: rev,
                next_paralog: None,
            });
            pos = s.saturating_add(l);
        }
        if pos < seq_end {
            final_tops.push(HalTopSegment {
                start: pos,
                length: seq_end - pos,
                parent_bottom_idx: None,
                parent_reversed: false,
                next_paralog: None,
            });
        }
    }

    let mut chain_head = vec![usize::MAX; root_bottoms.len()];
    let mut chain_tail = vec![usize::MAX; root_bottoms.len()];
    for t_idx in 0..final_tops.len() {
        let Some(pb) = final_tops[t_idx].parent_bottom_idx else {
            continue;
        };
        if chain_head[pb] == usize::MAX {
            chain_head[pb] = t_idx;
            chain_tail[pb] = t_idx;
            root_bottoms[pb].child_tops[c_idx] = Some((t_idx, final_tops[t_idx].parent_reversed));
        } else {
            let prev_tail = chain_tail[pb];
            final_tops[prev_tail].next_paralog = Some(t_idx);
            final_tops[t_idx].next_paralog = Some(chain_head[pb]);
            chain_tail[pb] = t_idx;
        }
    }

    final_tops
}

pub(super) struct NormalizedSegments {
    pub(super) tops_by_genome: Vec<Vec<HalTopSegment>>,
    pub(super) bottoms_by_genome: Vec<Vec<HalBottomSegment>>,
}

#[allow(clippy::too_many_lines)]
pub(super) fn normalize_hal_segments_for_export(hal: &HalAlignment) -> NormalizedSegments {
    use std::collections::{BTreeMap, BTreeSet};
    let n_genomes = hal.genomes.len();
    if n_genomes == 0 {
        return NormalizedSegments {
            tops_by_genome: Vec::new(),
            bottoms_by_genome: Vec::new(),
        };
    }
    let root = &hal.genomes[0];

    let num_children = root.children.len();
    let mut pts = BTreeSet::new();
    pts.insert(0usize);
    pts.insert(root.dna.len());
    for s in &root.sequences {
        pts.insert(s.offset.min(root.dna.len()));
        pts.insert(s.offset.saturating_add(s.length).min(root.dna.len()));
    }
    let mut sorted_orig_b: Vec<&HalBottomSegment> = root
        .bottom_segments
        .iter()
        .filter(|b| b.length > 0 && b.start < root.dna.len())
        .collect();
    for &b in &sorted_orig_b {
        pts.insert(b.start);
        pts.insert(b.start.saturating_add(b.length).min(root.dna.len()));
    }
    sorted_orig_b.sort_by_key(|b| b.start);
    let pts_vec: Vec<usize> = pts.into_iter().collect();

    let mut new_bottoms: Vec<HalBottomSegment> = Vec::with_capacity(pts_vec.len());
    let mut placed_tops: Vec<Vec<(usize, usize, usize, bool)>> = vec![Vec::new(); num_children];
    let mut occupied_child: Vec<BTreeMap<usize, usize>> = vec![BTreeMap::new(); num_children];
    let mut active: Vec<&HalBottomSegment> = Vec::new();
    let mut next_b = 0usize;

    for win in pts_vec.windows(2) {
        let a = win[0];
        let b = win[1];
        if b <= a {
            continue;
        }
        while next_b < sorted_orig_b.len() && sorted_orig_b[next_b].start <= a {
            active.push(sorted_orig_b[next_b]);
            next_b += 1;
        }
        active.retain(|orig_b| orig_b.start.saturating_add(orig_b.length) >= b);

        let sub_len = b - a;
        let b_idx = new_bottoms.len();
        let child_links = vec![None; num_children];

        for &orig_b in &active {
            for c_idx in 0..num_children {
                let Some(child_g_idx) = root.children.get(c_idx).copied() else {
                    continue;
                };
                let Some(child_g) = hal.genomes.get(child_g_idx) else {
                    continue;
                };
                let Some((mut cur_t_idx, mut cur_rev)) =
                    orig_b.child_tops.get(c_idx).and_then(|x| *x)
                else {
                    continue;
                };
                let start_t_idx = cur_t_idx;
                let mut seen_t = 0usize;
                while let Some(orig_t) = child_g.top_segments.get(cur_t_idx) {
                    let off_in_orig = a.saturating_sub(orig_b.start);
                    let t_start = if cur_rev {
                        let tail_off = orig_b
                            .length
                            .saturating_sub(off_in_orig.saturating_add(sub_len));
                        orig_t.start.checked_add(tail_off)
                    } else {
                        orig_t.start.checked_add(off_in_orig)
                    };
                    if let Some((t_start, t_end)) =
                        t_start.and_then(|s| s.checked_add(sub_len).map(|e| (s, e)))
                        && seq_contains_interval(&child_g.sequences, t_start, t_end)
                        && !interval_overlaps(&occupied_child[c_idx], t_start, t_end)
                    {
                        occupied_child[c_idx].insert(t_start, t_end);
                        placed_tops[c_idx].push((t_start, sub_len, b_idx, cur_rev));
                    }
                    seen_t += 1;
                    let Some(next_t) = orig_t.next_paralog else {
                        break;
                    };
                    if next_t == start_t_idx
                        || next_t == cur_t_idx
                        || seen_t >= child_g.top_segments.len()
                    {
                        break;
                    }
                    cur_t_idx = next_t;
                    cur_rev = child_g
                        .top_segments
                        .get(cur_t_idx)
                        .is_some_and(|t| t.parent_reversed);
                }
            }
        }

        new_bottoms.push(HalBottomSegment {
            start: a,
            length: sub_len,
            child_tops: child_links,
        });
    }

    let mut tops_by_genome = vec![Vec::new(); n_genomes];
    let mut bottoms_by_genome = vec![Vec::new(); n_genomes];
    for (c_idx, &child_g_idx) in root.children.iter().enumerate() {
        let Some(child_g) = hal.genomes.get(child_g_idx) else {
            continue;
        };
        if let Some(dst) = tops_by_genome.get_mut(child_g_idx) {
            *dst = finalize_child_top_segments(
                &child_g.sequences,
                std::mem::take(&mut placed_tops[c_idx]),
                c_idx,
                &mut new_bottoms,
            );
        }
    }
    bottoms_by_genome[0] = new_bottoms;

    NormalizedSegments {
        tops_by_genome,
        bottoms_by_genome,
    }
}

pub(super) fn normalize_star_hal_in_place(star: &mut HalAlignment) {
    let norm = normalize_hal_segments_for_export(star);
    for (g, (tops, bottoms)) in star
        .genomes
        .iter_mut()
        .zip(norm.tops_by_genome.into_iter().zip(norm.bottoms_by_genome))
    {
        g.top_segments = tops;
        g.bottom_segments = bottoms;
    }
}

pub(super) struct CompactTop {
    pub(super) start: u64,
    pub(super) parent_bot: u32,
    pub(super) next_paralog: u32,
    pub(super) rev: bool,
    pub(super) is_chain_head: bool,
}

pub(super) const NO_PARALOG: u32 = u32::MAX;

pub(super) struct ProgressivePartition {
    pub(super) bot_starts: Vec<u64>,
    pub(super) bot_lens: Vec<u64>,
    pub(super) q_tops: Vec<Vec<CompactTop>>,
}

fn normalize_query_tops(
    qry_seqs: &[HalSequence],
    mut segs: Vec<(usize, usize, u32, bool)>,
    n_bottoms: usize,
) -> std::io::Result<Vec<CompactTop>> {
    segs.sort_by_key(|&(s, _, _, _)| s);
    let mut final_tops: Vec<CompactTop> = Vec::with_capacity(segs.len() * 2 + 1);
    let mut seg_i = 0usize;
    for seq in qry_seqs {
        let mut pos = seq.offset;
        let seq_end = seq.offset.saturating_add(seq.length);
        while seg_i < segs.len() && segs[seg_i].0 < seq_end {
            let (s, l, p_b_idx, rev) = segs[seg_i];
            seg_i += 1;
            if s < pos {
                continue;
            }
            if s > pos {
                final_tops.push(CompactTop {
                    start: pos as u64,
                    parent_bot: NO_PARALOG,
                    next_paralog: NO_PARALOG,
                    rev: false,
                    is_chain_head: false,
                });
            }
            final_tops.push(CompactTop {
                start: s as u64,
                parent_bot: p_b_idx,
                next_paralog: NO_PARALOG,
                rev,
                is_chain_head: false,
            });
            pos = s + l;
        }
        if pos < seq_end {
            final_tops.push(CompactTop {
                start: pos as u64,
                parent_bot: NO_PARALOG,
                next_paralog: NO_PARALOG,
                rev: false,
                is_chain_head: false,
            });
        }
    }

    let mut chain_head = vec![NO_PARALOG; n_bottoms];
    let mut chain_tail = vec![NO_PARALOG; n_bottoms];
    for t_idx in 0..final_tops.len() {
        let pb = final_tops[t_idx].parent_bot;
        if pb == NO_PARALOG {
            continue;
        }
        let t_u32 = u32::try_from(t_idx)
            .ok()
            .filter(|&v| v < NO_PARALOG)
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "HAL top segment count exceeds u32::MAX - 1",
                )
            })?;
        let pb_idx = pb as usize;
        if chain_head[pb_idx] == NO_PARALOG {
            chain_head[pb_idx] = t_u32;
            chain_tail[pb_idx] = t_u32;
            final_tops[t_idx].is_chain_head = true;
        } else {
            let prev_tail = chain_tail[pb_idx] as usize;
            final_tops[prev_tail].next_paralog = t_u32;
            final_tops[t_idx].next_paralog = chain_head[pb_idx];
            chain_tail[pb_idx] = t_u32;
        }
    }
    Ok(final_tops)
}

type ProgressiveHit = (u64, u64, u32, u64, bool);

fn collect_progressive_hits_and_points(
    ref_dna_len: usize,
    ref_seqs: &[HalSequence],
    qry_seqs: &[Vec<HalSequence>],
    blocks: &[(usize, usize, usize, usize, usize, usize, bool)],
) -> (Vec<ProgressiveHit>, Vec<u64>) {
    let mut hits: Vec<ProgressiveHit> = Vec::with_capacity(blocks.len());
    let mut pts: Vec<u64> = Vec::with_capacity(2 + ref_seqs.len() * 2 + blocks.len() * 2);
    pts.push(0);
    pts.push(ref_dna_len as u64);
    for s in ref_seqs {
        pts.push(s.offset.min(ref_dna_len) as u64);
        pts.push(s.offset.saturating_add(s.length).min(ref_dna_len) as u64);
    }
    for &(q_idx, r_ci, q_ci, r_local, q_local, len, rev) in blocks {
        if len == 0 {
            continue;
        }
        let Some(r_seq) = ref_seqs.get(r_ci) else {
            continue;
        };
        let Some(q_seq) = qry_seqs.get(q_idx).and_then(|s| s.get(q_ci)) else {
            continue;
        };
        if r_local.checked_add(len).is_none_or(|e| e > r_seq.length)
            || q_local.checked_add(len).is_none_or(|e| e > q_seq.length)
        {
            continue;
        }
        let (Some(r_global), Some(q_global)) = (
            r_seq.offset.checked_add(r_local),
            q_seq.offset.checked_add(q_local),
        ) else {
            continue;
        };
        let Ok(q_u32) = u32::try_from(q_idx) else {
            continue;
        };
        let r_g = r_global as u64;
        let len_u64 = len as u64;
        pts.push(r_g);
        pts.push(r_g.saturating_add(len_u64).min(ref_dna_len as u64));
        hits.push((r_g, len_u64, q_u32, q_global as u64, rev));
    }
    pts.sort_unstable();
    pts.dedup();
    hits.sort_by_key(|&(r_g, len, _, _, _)| (r_g, len));
    (hits, pts)
}

pub(super) fn partition_progressive_segments(
    ref_dna_len: usize,
    ref_seqs: &[HalSequence],
    qry_seqs: &[Vec<HalSequence>],
    blocks: &[(usize, usize, usize, usize, usize, usize, bool)],
) -> std::io::Result<ProgressivePartition> {
    let n_q = qry_seqs.len();
    let (hits, pts) = collect_progressive_hits_and_points(ref_dna_len, ref_seqs, qry_seqs, blocks);

    let mut bot_starts: Vec<u64> = Vec::with_capacity(pts.len().saturating_sub(1));
    let mut bot_lens: Vec<u64> = Vec::with_capacity(pts.len().saturating_sub(1));
    let mut placed_tops: Vec<Vec<(usize, usize, u32, bool)>> =
        (0..n_q).map(|_| Vec::new()).collect();
    let mut occupied_q: Vec<std::collections::BTreeMap<usize, usize>> = (0..n_q)
        .map(|_| std::collections::BTreeMap::new())
        .collect();
    let mut active: Vec<ProgressiveHit> = Vec::new();
    let mut next_h = 0usize;

    for win in pts.windows(2) {
        let a = win[0];
        let b = win[1];
        if b <= a {
            continue;
        }
        while next_h < hits.len() && hits[next_h].0 <= a {
            active.push(hits[next_h]);
            next_h += 1;
        }
        active.retain(|h| h.0.saturating_add(h.1) >= b);

        let sub_len = b - a;
        let b_idx = u32::try_from(bot_starts.len())
            .ok()
            .filter(|&v| v < NO_PARALOG)
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "HAL bottom segment count exceeds u32::MAX - 1",
                )
            })?;
        bot_starts.push(a);
        bot_lens.push(sub_len);

        for &(r_g, orig_len, q_u32, q_g, rev) in &active {
            let q_idx = q_u32 as usize;
            let off_in_orig = a - r_g;
            let t_start_u64 = if rev {
                let tail_off = orig_len.saturating_sub(off_in_orig + sub_len);
                q_g.saturating_add(tail_off)
            } else {
                q_g.saturating_add(off_in_orig)
            };
            let (Ok(t_start), Ok(sub_len_usize)) =
                (usize::try_from(t_start_u64), usize::try_from(sub_len))
            else {
                continue;
            };
            let Some(t_end) = t_start.checked_add(sub_len_usize) else {
                continue;
            };
            if seq_contains_interval(&qry_seqs[q_idx], t_start, t_end)
                && !interval_overlaps(&occupied_q[q_idx], t_start, t_end)
            {
                occupied_q[q_idx].insert(t_start, t_end);
                placed_tops[q_idx].push((t_start, sub_len_usize, b_idx, rev));
            }
        }
    }

    let n_bottoms = bot_starts.len();
    let q_tops = placed_tops
        .into_iter()
        .enumerate()
        .map(|(q_idx, segs)| normalize_query_tops(&qry_seqs[q_idx], segs, n_bottoms))
        .collect::<std::io::Result<Vec<_>>>()?;

    Ok(ProgressivePartition {
        bot_starts,
        bot_lens,
        q_tops,
    })
}
