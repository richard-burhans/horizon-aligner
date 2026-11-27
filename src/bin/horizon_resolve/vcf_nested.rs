//! Per-chromosome discovery and Rayon-parallel genotyping of nested reference-anchored
//! local bubbles (`LV=1`) for `horizon_resolve --vcf --nested-snarls`.

use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

use crate::vcf::{
    CapturedWalk, HapId, SnarlKey, TraversalEntry, VcfBuildCtx, VcfRecord, build_vcf_record,
    canonical_snarl_key, insert_traversal,
};

const REF_ABSENT: u32 = u32::MAX;
const REF_DUPLICATED: u32 = u32::MAX - 1;

#[derive(Clone, Copy)]
pub(super) struct RefAnchorLoc {
    walk_idx: u32,
    step_idx: u32,
    ref_fwd: bool,
}

impl RefAnchorLoc {
    const ABSENT: Self = Self {
        walk_idx: REF_ABSENT,
        step_idx: 0,
        ref_fwd: true,
    };

    #[inline]
    fn is_unique(self) -> bool {
        self.walk_idx < REF_DUPLICATED
    }
}

pub(super) struct NestedPlan {
    pub(super) seg_ref_loc: Vec<RefAnchorLoc>,
    pub(super) bubbles_by_chrom: Vec<Vec<(u32, u32, u32)>>,
    pub(super) walks_by_chrom: Vec<Vec<usize>>,
}

/// Pass A: map unique sequence-bearing reference segments and partition candidate `LV=1` nested
/// bubbles and touching walks by reference chromosome so Pass B runs one chromosome at a time.
#[allow(clippy::too_many_arguments)]
pub(super) fn plan_nested_bubbles_by_chrom(
    ctx: &VcfBuildCtx<'_>,
    walks: &[CapturedWalk],
    ref_walk_indices: &[usize],
    ref_walk_chrom_ids: &[usize],
    num_chroms: usize,
    lv0_keys: &HashSet<SnarlKey>,
    max_nested_steps: usize,
) -> NestedPlan {
    let mut seg_ref_loc = vec![RefAnchorLoc::ABSENT; ctx.graph.num_segs()];
    for (rw_i, &w_idx) in ref_walk_indices.iter().enumerate() {
        let Ok(rw_u32) = u32::try_from(rw_i) else {
            break;
        };
        if rw_u32 >= REF_DUPLICATED {
            break;
        }
        for (s_idx, &step) in walks[w_idx].steps.iter().enumerate() {
            if !ctx.seg_seqs.has_seq(step.id()) {
                continue;
            }
            let Ok(step_u32) = u32::try_from(s_idx) else {
                continue;
            };
            let slot = &mut seg_ref_loc[step.id() as usize];
            if slot.walk_idx == REF_ABSENT {
                *slot = RefAnchorLoc {
                    walk_idx: rw_u32,
                    step_idx: step_u32,
                    ref_fwd: step.fwd(),
                };
            } else {
                slot.walk_idx = REF_DUPLICATED;
            }
        }
    }

    let mut bubbles_by_chrom: Vec<Vec<(u32, u32, u32)>> = vec![Vec::new(); num_chroms];
    let mut walks_by_chrom: Vec<Vec<usize>> = vec![Vec::new(); num_chroms];
    let mut seen_bubbles: HashSet<(u32, u32, u32)> = HashSet::new();
    let mut touched_stamp: Vec<u32> = vec![0; num_chroms];

    for (w_idx, w) in walks.iter().enumerate() {
        #[allow(clippy::cast_possible_truncation)]
        let stamp = (w_idx as u32).wrapping_add(1);
        let mut prev: Option<(usize, u32, u32, bool)> = None;
        for (k, &step) in w.steps.iter().enumerate() {
            let loc = seg_ref_loc[step.id() as usize];
            if !loc.is_unique() {
                continue;
            }
            let c_idx = ref_walk_chrom_ids[loc.walk_idx as usize];
            if touched_stamp[c_idx] != stamp {
                touched_stamp[c_idx] = stamp;
                walks_by_chrom[c_idx].push(w_idx);
            }
            let aligned_fwd = step.fwd() == loc.ref_fwd;
            if let Some((pk, prw, pr_step, p_aligned)) = prev
                && prw == loc.walk_idx
                && p_aligned == aligned_fwd
            {
                let (r_start, r_end) = if aligned_fwd {
                    (pr_step, loc.step_idx)
                } else {
                    (loc.step_idx, pr_step)
                };
                if r_end > r_start {
                    let ref_span = (r_end - r_start) as usize;
                    let alt_span = k - pk;
                    if ref_span <= max_nested_steps
                        && alt_span <= max_nested_steps
                        && (ref_span > 1 || alt_span > 1)
                    {
                        let rw = &walks[ref_walk_indices[prw as usize]];
                        let u = rw.steps[r_start as usize];
                        let v = rw.steps[r_end as usize];
                        let (ck, _) = canonical_snarl_key(u, v);
                        if !lv0_keys.contains(&ck) && seen_bubbles.insert((prw, r_start, r_end)) {
                            bubbles_by_chrom[c_idx].push((prw, r_start, r_end));
                        }
                    }
                }
            }
            prev = Some((k, loc.walk_idx, loc.step_idx, aligned_fwd));
        }
    }

    NestedPlan {
        seg_ref_loc,
        bubbles_by_chrom,
        walks_by_chrom,
    }
}

/// Pass B (single chromosome): genotype only the walks touching `chrom_idx` across its `LV=1`
/// bubbles, then build polymorphic VCF records in parallel via Rayon and drop traversal tables.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn genotype_chrom_nested_bubbles(
    ctx: &VcfBuildCtx<'_>,
    walks: &[CapturedWalk],
    walk_hap_ids: &[HapId],
    ref_walk_indices: &[usize],
    ref_walk_chrom_ids: &[usize],
    ref_step_ends: &[Vec<usize>],
    seg_ref_loc: &[RefAnchorLoc],
    chrom_idx: usize,
    bubbles: &[(u32, u32, u32)],
    touching_walks: &[usize],
    max_nested_steps: usize,
) -> Vec<VcfRecord> {
    if bubbles.is_empty() {
        return Vec::new();
    }
    let mut nested_starts: HashSet<(u32, u32)> = HashSet::with_capacity(bubbles.len());
    let mut nested_ends: HashSet<(u32, u32)> = HashSet::with_capacity(bubbles.len());
    let mut seen_bubbles: HashSet<(u32, u32, u32)> = HashSet::with_capacity(bubbles.len());
    for &(rw, rs, re) in bubbles {
        nested_starts.insert((rw, rs));
        nested_ends.insert((rw, re));
        seen_bubbles.insert((rw, rs, re));
    }

    let mut nested_traversals: HashMap<(u32, u32, u32), HashMap<HapId, TraversalEntry>> =
        HashMap::with_capacity(bubbles.len());
    let mut hits: Vec<(usize, u32, u32, bool)> = Vec::new();

    for &w_idx in touching_walks {
        let w = &walks[w_idx];
        let hap_id = walk_hap_ids[w_idx];
        hits.clear();
        for (k, &step) in w.steps.iter().enumerate() {
            let loc = seg_ref_loc[step.id() as usize];
            if loc.is_unique() && ref_walk_chrom_ids[loc.walk_idx as usize] == chrom_idx {
                hits.push((k, loc.walk_idx, loc.step_idx, step.fwd() == loc.ref_fwd));
            }
        }
        for hi in 0..hits.len() {
            let (step_a, rw_a, ref_pos_a, fwd_a) = hits[hi];
            if fwd_a {
                if !nested_starts.contains(&(rw_a, ref_pos_a)) {
                    continue;
                }
                for &(step_b, rw_b, ref_pos_b, fwd_b) in &hits[(hi + 1)..] {
                    if step_b - step_a > max_nested_steps {
                        break;
                    }
                    if rw_b != rw_a || !fwd_b {
                        continue;
                    }
                    if ref_pos_b > ref_pos_a && seen_bubbles.contains(&(rw_a, ref_pos_a, ref_pos_b))
                    {
                        insert_traversal(
                            nested_traversals
                                .entry((rw_a, ref_pos_a, ref_pos_b))
                                .or_default(),
                            hap_id,
                            &w.steps[(step_a + 1)..step_b],
                            false,
                        );
                    }
                }
            } else if nested_ends.contains(&(rw_a, ref_pos_a)) {
                for &(step_b, rw_b, ref_pos_b, fwd_b) in &hits[(hi + 1)..] {
                    if step_b - step_a > max_nested_steps {
                        break;
                    }
                    if rw_b != rw_a || fwd_b {
                        continue;
                    }
                    if ref_pos_b < ref_pos_a && seen_bubbles.contains(&(rw_a, ref_pos_b, ref_pos_a))
                    {
                        insert_traversal(
                            nested_traversals
                                .entry((rw_a, ref_pos_b, ref_pos_a))
                                .or_default(),
                            hap_id,
                            &w.steps[(step_a + 1)..step_b],
                            true,
                        );
                    }
                }
            }
        }
    }

    drop(nested_starts);
    drop(nested_ends);
    drop(seen_bubbles);
    drop(hits);

    bubbles
        .par_iter()
        .map_init(
            || (Vec::new(), Vec::new()),
            |(scratch_int, scratch_seq), &(rw_i, r_start, r_end)| {
                let hap_map = nested_traversals.get(&(rw_i, r_start, r_end))?;
                let rw = &walks[ref_walk_indices[rw_i as usize]];
                let step_ends = &ref_step_ends[rw_i as usize];
                let u = rw.steps[r_start as usize];
                let v = rw.steps[r_end as usize];
                let ref_interior = &rw.steps[(r_start as usize + 1)..(r_end as usize)];
                let ref_interior_start_0 = step_ends[r_start as usize];
                build_vcf_record(
                    ctx,
                    u,
                    v,
                    ref_interior,
                    ref_interior_start_0,
                    false,
                    hap_map,
                    1,
                    scratch_int,
                    scratch_seq,
                )
            },
        )
        .flatten()
        .collect()
}
