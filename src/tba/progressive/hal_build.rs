//! Hierarchical pure-Rust HAL graph construction from bottom-up progressive parent->child runs.

use super::{GenomeContigs, HalEdgeRecord};
use crate::hal::{HalAlignment, HalBottomSegment, HalGenome, HalSequence, HalTopSegment};
use crate::tba::model::PairwiseTbaBlock;
use std::collections::{BTreeSet, HashMap};

fn concat_genome_sequences(
    contigs: &[(String, Vec<u8>)],
) -> (Vec<HalSequence>, Vec<u8>, HashMap<String, usize>) {
    let mut sequences = Vec::with_capacity(contigs.len());
    let mut dna = Vec::new();
    let mut offsets = HashMap::with_capacity(contigs.len());
    let mut cur = 0usize;
    for (name, seq) in contigs {
        offsets.insert(name.clone(), cur);
        sequences.push(HalSequence {
            name: name.clone(),
            length: seq.len(),
            offset: cur,
        });
        dna.extend(seq.iter().map(u8::to_ascii_uppercase));
        cur += seq.len();
    }
    (sequences, dna, offsets)
}

/// Extract colinear ungapped `(parent_global_start, child_global_start, len, child_reversed)` runs
/// from `parent_to_child` pairwise blocks.
fn pairwise_blocks_to_global_runs(
    parent_name: &str,
    parent_offsets: &HashMap<String, usize>,
    child_offsets: &HashMap<String, usize>,
    blocks: &[PairwiseTbaBlock],
) -> Vec<(usize, usize, usize, bool)> {
    let mut runs = Vec::new();
    for pb in blocks {
        let (p_row, c_row) = if pb.row_a.species == parent_name {
            (&pb.row_a, &pb.row_b)
        } else {
            (&pb.row_b, &pb.row_a)
        };
        let Some(&p_base) = parent_offsets.get(&p_row.contig) else {
            continue;
        };
        let Some(&c_base) = child_offsets.get(&c_row.contig) else {
            continue;
        };
        let p_rev = p_row.strand == '-';
        let c_rev = c_row.strand == '-';
        let rel_rev = p_rev != c_rev;
        let mut p_idx = 0usize;
        let mut c_idx = 0usize;
        let mut cur_run: Option<(usize, usize, usize)> = None;

        for (&pb_ch, &cb_ch) in p_row.seq.iter().zip(&c_row.seq) {
            let p_has = pb_ch != b'-';
            let c_has = cb_ch != b'-';
            if p_has && c_has {
                let Some(p_fwd) = (if p_rev {
                    p_row.src_size.checked_sub(p_row.start + p_idx + 1)
                } else {
                    Some(p_row.start + p_idx)
                }) else {
                    if let Some((ps, cs, len)) = cur_run.take() {
                        runs.push((ps, cs, len, rel_rev));
                    }
                    p_idx += 1;
                    c_idx += 1;
                    continue;
                };
                let Some(c_fwd) = (if c_rev {
                    c_row.src_size.checked_sub(c_row.start + c_idx + 1)
                } else {
                    Some(c_row.start + c_idx)
                }) else {
                    if let Some((ps, cs, len)) = cur_run.take() {
                        runs.push((ps, cs, len, rel_rev));
                    }
                    p_idx += 1;
                    c_idx += 1;
                    continue;
                };
                let gp = p_base.saturating_add(p_fwd);
                let gc = c_base.saturating_add(c_fwd);
                let extends = match cur_run {
                    Some((ps, cs, len)) => {
                        let p_ok = if p_rev {
                            gp.checked_add(1) == Some(ps)
                        } else {
                            ps.checked_add(len) == Some(gp)
                        };
                        let c_ok = if c_rev {
                            gc.checked_add(1) == Some(cs)
                        } else {
                            cs.checked_add(len) == Some(gc)
                        };
                        p_ok && c_ok
                    }
                    None => false,
                };
                if extends {
                    if let Some((ref mut ps, ref mut cs, ref mut len)) = cur_run {
                        if p_rev {
                            *ps = gp;
                        }
                        if c_rev {
                            *cs = gc;
                        }
                        *len += 1;
                    }
                } else {
                    if let Some((ps, cs, len)) = cur_run.take() {
                        runs.push((ps, cs, len, rel_rev));
                    }
                    cur_run = Some((gp, gc, 1));
                }
            } else if let Some((ps, cs, len)) = cur_run.take() {
                runs.push((ps, cs, len, rel_rev));
            }
            if p_has {
                p_idx += 1;
            }
            if c_has {
                c_idx += 1;
            }
        }
        if let Some((ps, cs, len)) = cur_run {
            runs.push((ps, cs, len, rel_rev));
        }
    }
    runs.sort_unstable_by_key(|&(ps, _, _, _)| ps);
    runs
}

fn dfs_hal_genome_order<'a>(
    name: &'a str,
    edge_by_parent: &HashMap<&'a str, &'a HalEdgeRecord>,
    ordered: &mut Vec<String>,
    map: &mut HashMap<String, usize>,
) {
    if map.contains_key(name) {
        return;
    }
    let idx = ordered.len();
    ordered.push(name.to_owned());
    map.insert(name.to_owned(), idx);
    if let Some(edge) = edge_by_parent.get(name) {
        dfs_hal_genome_order(&edge.left_child, edge_by_parent, ordered, map);
        dfs_hal_genome_order(&edge.right_child, edge_by_parent, ordered, map);
    }
}

#[allow(clippy::too_many_lines)]
pub(super) fn build_progressive_hal<S: std::hash::BuildHasher>(
    root_name: &str,
    tree_newick: &str,
    leaf_genomes: &HashMap<String, GenomeContigs, S>,
    ancestors: &[super::AncestralGenomeEntry],
    edges: &[HalEdgeRecord],
) -> HalAlignment {
    // Order genomes top-down starting with root so root is index 0.
    let mut name_to_idx: HashMap<String, usize> = HashMap::new();
    let mut ordered_names: Vec<String> = Vec::new();
    let edge_by_parent: HashMap<&str, &HalEdgeRecord> =
        edges.iter().map(|e| (e.parent.as_str(), e)).collect();
    dfs_hal_genome_order(
        root_name,
        &edge_by_parent,
        &mut ordered_names,
        &mut name_to_idx,
    );

    let anc_by_name: HashMap<&str, &[(String, Vec<u8>)]> = ancestors
        .iter()
        .map(|(n, ctgs, _)| (n.as_str(), ctgs.as_slice()))
        .collect();

    let mut genomes: Vec<HalGenome> = Vec::with_capacity(ordered_names.len());
    let mut offsets_by_idx: Vec<HashMap<String, usize>> = Vec::with_capacity(ordered_names.len());
    for name in &ordered_names {
        let ctgs = leaf_genomes
            .get(name)
            .map(Vec::as_slice)
            .or_else(|| anc_by_name.get(name.as_str()).copied())
            .unwrap_or(&[]);
        let (sequences, dna, offsets) = concat_genome_sequences(ctgs);
        genomes.push(HalGenome {
            name: name.clone(),
            parent: None,
            children: Vec::new(),
            sequences,
            dna,
            top_segments: Vec::new(),
            bottom_segments: Vec::new(),
        });
        offsets_by_idx.push(offsets);
    }

    for edge in edges {
        let p_idx = name_to_idx[&edge.parent];
        let l_idx = name_to_idx[&edge.left_child];
        let r_idx = name_to_idx[&edge.right_child];
        genomes[p_idx].children = vec![l_idx, r_idx];
        genomes[l_idx].parent = Some(p_idx);
        genomes[r_idx].parent = Some(p_idx);

        let p_len = genomes[p_idx].dna.len();
        let l_runs = pairwise_blocks_to_global_runs(
            &edge.parent,
            &offsets_by_idx[p_idx],
            &offsets_by_idx[l_idx],
            &edge.parent_to_left,
        );
        let r_runs = pairwise_blocks_to_global_runs(
            &edge.parent,
            &offsets_by_idx[p_idx],
            &offsets_by_idx[r_idx],
            &edge.parent_to_right,
        );

        // Partition parent coordinate space at all run start/end breakpoints.
        let mut breaks: BTreeSet<usize> = BTreeSet::new();
        breaks.insert(0);
        breaks.insert(p_len);
        for &(ps, _, len, _) in l_runs.iter().chain(&r_runs) {
            if ps <= p_len {
                breaks.insert(ps);
            }
            if ps + len <= p_len {
                breaks.insert(ps + len);
            }
        }
        let pts: Vec<usize> = breaks.into_iter().collect();

        // Map each parent interval [a, b) to a single HalBottomSegment and corresponding child TopSegments,
        // using IntervalSet availability so split inversions and translocations are never dropped.
        let l_len = genomes[l_idx].dna.len();
        let r_len = genomes[r_idx].dna.len();
        let mut l_top_mapped: Vec<(usize, usize, usize, bool)> = Vec::new();
        let mut r_top_mapped: Vec<(usize, usize, usize, bool)> = Vec::new();
        let mut l_avail = crate::netting::IntervalSet::from_single(0, l_len);
        let mut r_avail = crate::netting::IntervalSet::from_single(0, r_len);
        let mut avail_scratch: Vec<usize> = Vec::new();
        let mut bottom_segs = Vec::with_capacity(pts.len().saturating_sub(1));

        // O(log R) lookup into runs sorted by parent start `ps`.
        let find_run = |runs: &[(usize, usize, usize, bool)], a: usize, b: usize| {
            let hi_idx = runs.partition_point(|&(ps, _, _, _)| ps <= a);
            for &(ps, cs, len, rev) in runs[..hi_idx].iter().rev() {
                if ps.saturating_add(len) >= b {
                    let off = a - ps;
                    let sub_len = b - a;
                    let c_start = if rev {
                        cs + len - off - sub_len
                    } else {
                        cs + off
                    };
                    return Some((c_start, sub_len, rev));
                }
                if ps.saturating_add(len) <= a {
                    break;
                }
            }
            None
        };

        for win in pts.windows(2) {
            let (a, b) = (win[0], win[1]);
            if b <= a {
                continue;
            }
            let b_idx = bottom_segs.len();
            if let Some((cs, sl, rev)) = find_run(&l_runs, a, b)
                && cs.saturating_add(sl) <= l_len
                && l_avail.contains_range(cs, cs + sl)
            {
                l_avail.subtract(cs, cs + sl, &mut avail_scratch);
                l_top_mapped.push((cs, sl, b_idx, rev));
            }
            if let Some((cs, sl, rev)) = find_run(&r_runs, a, b)
                && cs.saturating_add(sl) <= r_len
                && r_avail.contains_range(cs, cs + sl)
            {
                r_avail.subtract(cs, cs + sl, &mut avail_scratch);
                r_top_mapped.push((cs, sl, b_idx, rev));
            }
            bottom_segs.push(HalBottomSegment {
                start: a,
                length: b - a,
                child_tops: vec![None, None],
            });
        }

        // Sort mapped child intervals by child coordinate so split reversed runs and
        // rearrangements tile cleanly in ascending child order.
        l_top_mapped.sort_unstable_by_key(|&(cs, _, _, _)| cs);
        r_top_mapped.sort_unstable_by_key(|&(cs, _, _, _)| cs);

        // Convert mapped child intervals into complete tiled TopSegments (filling unaligned insertions).
        let build_child_tops = |c_len: usize,
                                mapped: &[(usize, usize, usize, bool)],
                                b_segs: &mut [HalBottomSegment],
                                child_slot: usize|
         -> Vec<HalTopSegment> {
            let mut tops = Vec::new();
            let mut cur = 0usize;
            for &(cs, sl, b_idx, rev) in mapped {
                if cs < cur || cs.saturating_add(sl) > c_len {
                    continue;
                }
                if cs > cur {
                    tops.push(HalTopSegment {
                        start: cur,
                        length: cs - cur,
                        parent_bottom_idx: None,
                        parent_reversed: false,
                        next_paralog: None,
                    });
                }
                let final_t_idx = tops.len();
                tops.push(HalTopSegment {
                    start: cs,
                    length: sl,
                    parent_bottom_idx: Some(b_idx),
                    parent_reversed: rev,
                    next_paralog: None,
                });
                b_segs[b_idx].child_tops[child_slot] = Some((final_t_idx, rev));
                cur = cs + sl;
            }
            if cur < c_len {
                tops.push(HalTopSegment {
                    start: cur,
                    length: c_len - cur,
                    parent_bottom_idx: None,
                    parent_reversed: false,
                    next_paralog: None,
                });
            }
            tops
        };

        let l_tops = build_child_tops(l_len, &l_top_mapped, &mut bottom_segs, 0);
        let r_tops = build_child_tops(r_len, &r_top_mapped, &mut bottom_segs, 1);
        genomes[l_idx].top_segments = l_tops;
        genomes[r_idx].top_segments = r_tops;
        genomes[p_idx].bottom_segments = bottom_segs;
    }

    HalAlignment {
        tree_newick: format!("{tree_newick};"),
        genomes,
    }
}
