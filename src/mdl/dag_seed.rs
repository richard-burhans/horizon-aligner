//! Seeded subgraph selection and component partitioning for [`align_seq_to_dag_mdl`].

use super::dag::{GraphAlignment, GraphNode, align_seq_to_dag_mdl};
use crate::mdl::dp::HmmParams;
use crate::seeding::{FlatSeedIndex, MASK_RY18, MASK_W10, aa5_keys, sdust_mask};
use std::collections::HashMap;

/// Pre-built per-node spaced-seed indices (`MASK_W10`, `MASK_RY18`, and compact 5-of-6 `aa5`
/// translated keys) for repeated queries against the same sequence DAG.
#[derive(Debug, Clone)]
pub struct DagSeedIndex {
    w10: Vec<Option<FlatSeedIndex>>,
    ry18: Vec<Option<FlatSeedIndex>>,
    aa5: Vec<Vec<u16>>,
}

impl DagSeedIndex {
    /// Build `MASK_W10`, `MASK_RY18`, and compact `aa5` seed indices once for all nodes in `nodes`.
    #[must_use]
    pub fn build(nodes: &[GraphNode]) -> Self {
        let w10 = nodes
            .iter()
            .map(|nd| {
                (nd.seq.len() >= MASK_W10.len())
                    .then(|| {
                        FlatSeedIndex::build_with_sdust(&nd.seq, MASK_W10, 8, false, false, None)
                            .ok()
                    })
                    .flatten()
            })
            .collect();
        let ry18 = nodes
            .iter()
            .map(|nd| {
                (nd.seq.len() >= MASK_RY18.len())
                    .then(|| {
                        FlatSeedIndex::build_with_sdust(&nd.seq, MASK_RY18, 8, false, true, None)
                            .ok()
                    })
                    .flatten()
            })
            .collect();
        let aa5 = nodes
            .iter()
            .map(|nd| {
                if nd.seq.len() >= 18 {
                    aa5_keys(&nd.seq)
                } else {
                    Vec::new()
                }
            })
            .collect();
        Self { w10, ry18, aa5 }
    }
}

/// Estimate divergence-adaptive [`HmmParams`] for `seq_y` against a sequence DAG by reusing
/// the pre-built `MASK_W10` indices in `dag_index`.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn estimate_dag_hmm_params(
    nodes: &[GraphNode],
    dag_index: &DagSeedIndex,
    seq_y: &[u8],
) -> HmmParams {
    use crate::seeding::complement_base;
    let Ok(q_idx) = FlatSeedIndex::build_with_sdust(seq_y, MASK_W10, 4, false, false, None) else {
        return HmmParams::default();
    };
    let span = MASK_W10.len();
    let mut matches = 0usize;
    let mut total = 0usize;
    let mut indel_events = 0usize;
    let mut total_indel_len = 0usize;
    let mut total_span = 0usize;
    for (nd, n_idx_opt) in nodes.iter().zip(dag_index.w10.iter()) {
        let Some(n_idx) = n_idx_opt.as_ref() else {
            continue;
        };
        let (fwd_hits, rc_hits) = n_idx.intersect(&q_idx);
        let mut node_hits_fwd: Vec<(usize, usize)> = Vec::new();
        let mut node_hits_rc: Vec<(usize, usize)> = Vec::new();
        // `FlatSeedIndex::intersect` emits paired canonical entries in `fwd_hits` and `rc_hits`.
        // Compare both orientations per canonical hit and retain only the matching strand so
        // reverse-complement queries never saturate `total >= 2000` on bogus forward comparisons.
        for (&(np, qp_fwd), &(_, qp_rc_raw)) in fwd_hits.iter().zip(rc_hits.iter()) {
            if np + span > nd.seq.len() {
                continue;
            }
            let eval_strand = |is_rc: bool, qp: usize| -> (usize, usize) {
                if qp + span > seq_y.len() {
                    return (0, 0);
                }
                let mut m_cnt = 0usize;
                let mut valid = 0usize;
                for off in 0..span {
                    let bx = nd.seq[np + off].to_ascii_uppercase();
                    let by = if is_rc {
                        complement_base(seq_y[qp + span - 1 - off]).to_ascii_uppercase()
                    } else {
                        seq_y[qp + off].to_ascii_uppercase()
                    };
                    if matches!(bx, b'A' | b'C' | b'G' | b'T')
                        && matches!(by, b'A' | b'C' | b'G' | b'T')
                    {
                        valid += 1;
                        if bx == by {
                            m_cnt += 1;
                        }
                    }
                }
                (m_cnt, valid)
            };
            let (m_fwd, v_fwd) = eval_strand(false, qp_fwd);
            let (m_rc, v_rc) = eval_strand(true, qp_fwd);
            if v_fwd == 0 && v_rc == 0 {
                continue;
            }
            if m_rc > m_fwd && v_rc > 0 {
                node_hits_rc.push((np, qp_rc_raw));
                total += v_rc;
                matches += m_rc;
            } else if v_fwd > 0 {
                node_hits_fwd.push((np, qp_fwd));
                total += v_fwd;
                matches += m_fwd;
            }
            if total >= 2000 {
                break;
            }
        }
        for mut strand_hits in [node_hits_fwd, node_hits_rc] {
            strand_hits.sort_unstable();
            for w in strand_hits.windows(2).take(64) {
                let dx = w[1].0.saturating_sub(w[0].0);
                let dy = w[1].1.saturating_sub(w[0].1);
                if dx > 0 && dy > 0 && dx <= 500 && dy <= 500 {
                    total_span += dx.max(dy);
                    if dx != dy {
                        indel_events += 1;
                        total_indel_len += dx.abs_diff(dy);
                    }
                }
            }
        }
        if total >= 2000 {
            break;
        }
    }
    if total < 60 {
        return HmmParams::default();
    }
    let defaults = HmmParams::default();
    #[allow(clippy::cast_precision_loss)]
    let identity_prior = ((matches as f64 + 195.0) / (total as f64 + 300.0)).clamp(0.60, 0.72);
    #[allow(clippy::cast_precision_loss)]
    let gap_open_bits = if total_span >= 200 {
        ((indel_events as f64 + 1.0) / (total_span as f64 + 22.0))
            .log2()
            .clamp(-5.5, -4.0)
    } else {
        defaults.gap_open_bits
    };
    #[allow(clippy::cast_precision_loss)]
    let gap_extend_bits = if indel_events >= 3 && total_indel_len > 0 {
        let mean_len = (total_indel_len as f64 + 2.5) / (indel_events as f64 + 1.0);
        ((mean_len - 1.0) / mean_len)
            .max(0.01)
            .log2()
            .clamp(-1.8, -0.5)
    } else {
        defaults.gap_extend_bits
    };
    HmmParams {
        identity_prior,
        gap_open_bits,
        gap_extend_bits,
    }
}

fn sorted_u16_intersection_count(a: &[u16], b: &[u16]) -> usize {
    let mut i = 0usize;
    let mut j = 0usize;
    let mut hits = 0usize;
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                hits += 1;
                i += 1;
                j += 1;
            }
        }
    }
    hits
}

fn try_small_or_invalid_dag(
    nodes: &[GraphNode],
    seq_y: &[u8],
    params: &HmmParams,
    min_block_bits: f64,
) -> Option<GraphAlignment> {
    if nodes.is_empty()
        || nodes
            .iter()
            .enumerate()
            .any(|(u, nd)| nd.preds.iter().any(|&p| p >= u))
    {
        return Some(GraphAlignment {
            path_nodes: Vec::new(),
            kept_pairs: Vec::new(),
            bubbles: usize::from(seq_y.len() >= 6),
            score_bits: 0.0,
        });
    }
    let total_bp: usize = nodes.iter().map(|n| n.seq.len()).sum();
    if seq_y.is_empty()
        || (nodes.len() <= 8
            && total_bp <= 50_000
            && total_bp.saturating_mul(seq_y.len()) <= 16_000_000)
    {
        return Some(align_seq_to_dag_mdl(nodes, seq_y, params, min_block_bits));
    }
    None
}

/// Cluster forward seed hits `(node_pos, q_pos)` by diagonal `d = node_pos - q_pos` and return
/// the projected genomic bounds `(est_start, est_end)` of the densest collinear diagonal band
/// (immune to sporadic off-diagonal 10-mer collisions on long rGFA backbone segments).
fn densest_hit_projected_bounds(
    hits: &[(usize, usize)],
    query_len: usize,
    diag_band: isize,
) -> Option<(usize, usize)> {
    if hits.is_empty() {
        return None;
    }
    let mut diags: Vec<isize> = hits
        .iter()
        .map(|&(np, qp)| np.cast_signed().saturating_sub(qp.cast_signed()))
        .collect();
    diags.sort_unstable();
    let mut best_count = 0usize;
    let mut best_span = isize::MAX;
    let mut best_lr = (0usize, 0usize);
    let mut left = 0usize;
    for right in 0..diags.len() {
        while diags[right].saturating_sub(diags[left]) > diag_band {
            left += 1;
        }
        let count = right - left + 1;
        let span = diags[right].saturating_sub(diags[left]);
        if count > best_count || (count == best_count && span < best_span) {
            best_count = count;
            best_span = span;
            best_lr = (left, right);
        }
    }
    let d_lo = diags[best_lr.0];
    let d_hi = diags[best_lr.1];
    let est_start = usize::try_from(d_lo.max(0)).unwrap_or(0);
    let est_end = usize::try_from(d_hi.saturating_add(query_len.cast_signed()).max(0)).unwrap_or(0);
    Some((est_start, est_end))
}

/// Compute per-node slice interval `[slice_start, slice_end)` for a long node (`seq.len() > 4_000`)
/// around its densest seed-hit cluster's projected query bounds, keeping the flank margin strictly
/// inside `align_seq_to_dag_mdl`'s diagonal band (`(m/8).clamp(64, 512)`).
fn compute_node_slice(
    nid: usize,
    nd: &GraphNode,
    hit_bounds: &HashMap<usize, (usize, usize)>,
    query_len: usize,
    is_pred_of_seeded: bool,
    is_succ_of_seeded: bool,
) -> (usize, usize) {
    let n_len = nd.seq.len();
    let max_unsliced = (query_len.saturating_mul(2)).max(4_096);
    if n_len <= max_unsliced {
        return (0, n_len);
    }
    let band_w = (query_len / 8).clamp(64, 512);
    let margin = (band_w / 2).clamp(32, 128);
    if let Some(&(est_s, est_e)) = hit_bounds.get(&nid) {
        let mut s = est_s.saturating_sub(margin);
        let mut e = est_e.saturating_add(margin).min(n_len);
        // Snap near-boundary slices to exact node ends so graph edges stay intact across bubbles.
        if s <= margin {
            s = 0;
        }
        if n_len.saturating_sub(e) <= margin {
            e = n_len;
        }
        return (s, e.max(s));
    }
    let flank = margin.min(n_len);
    if is_pred_of_seeded && !is_succ_of_seeded {
        (n_len - flank, n_len)
    } else {
        (0, flank)
    }
}

/// Select the highest-scoring topologically ordered anchor chain in the DAG when multiple
/// disconnected or branching seeded regions exist, returning the set of anchor nodes that lie
/// on the winning chain (within `max_bridge_bp` between consecutive anchors).
fn select_chained_dag_anchors(
    nodes: &[GraphNode],
    seed_scores: &HashMap<usize, usize>,
    anchor_floor: usize,
    max_bridge_bp: usize,
) -> Vec<usize> {
    let n = nodes.len();
    // DP state per node: (best_accumulated_seed_score, bp_since_last_anchor, prev_node_in_chain)
    let mut best_score = vec![0usize; n];
    let mut dist_since_anchor = vec![usize::MAX; n];
    let mut prev_anchor = vec![usize::MAX; n];
    let mut last_anchor_at = vec![usize::MAX; n];

    for (u, nd) in nodes.iter().enumerate() {
        let mut in_score = 0usize;
        let mut in_dist = usize::MAX;
        let mut in_last_anchor = usize::MAX;
        for &p in &nd.preds {
            // dist_since_anchor[p] holds intervening non-anchor bp up to and including p
            // (0 when p itself is an anchor, so long anchor nodes still chain to successors).
            if p < u && dist_since_anchor[p] <= max_bridge_bp {
                let d = dist_since_anchor[p];
                if best_score[p] > in_score || (best_score[p] == in_score && d < in_dist) {
                    in_score = best_score[p];
                    in_dist = d;
                    in_last_anchor = last_anchor_at[p];
                }
            }
        }
        let own_score = seed_scores
            .get(&u)
            .copied()
            .filter(|&sc| sc >= anchor_floor)
            .unwrap_or(0);
        if own_score > 0 {
            best_score[u] = in_score.saturating_add(own_score);
            dist_since_anchor[u] = 0;
            prev_anchor[u] = in_last_anchor;
            last_anchor_at[u] = u;
        } else if in_dist <= max_bridge_bp {
            let next_dist = in_dist.saturating_add(nd.seq.len());
            if next_dist <= max_bridge_bp {
                best_score[u] = in_score;
                dist_since_anchor[u] = next_dist;
                last_anchor_at[u] = in_last_anchor;
            }
        }
    }

    // Backtrack from the highest-scoring anchor node.
    let Some(mut cur) = (0..n)
        .filter(|&u| dist_since_anchor[u] == 0 && best_score[u] > 0)
        .max_by_key(|&u| (best_score[u], std::cmp::Reverse(u)))
    else {
        return Vec::new();
    };
    let mut chain = Vec::new();
    while cur != usize::MAX {
        chain.push(cur);
        cur = prev_anchor[cur];
    }
    chain.reverse();
    chain
}

/// Seeded DAG alignment using a pre-built [`DagSeedIndex`] to isolate the
/// topologically chained DAG subgraph matching `seq_y`, then runs [`align_seq_to_dag_mdl`].
/// Maps returned node indices back to the original `nodes` slice.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn align_seq_to_dag_seeded_indexed(
    nodes: &[GraphNode],
    index: &DagSeedIndex,
    seq_y: &[u8],
    params: &HmmParams,
    min_block_bits: f64,
) -> GraphAlignment {
    if let Some(aln) = try_small_or_invalid_dag(nodes, seq_y, params, min_block_bits) {
        return aln;
    }
    let total_bp: usize = nodes.iter().map(|n| n.seq.len()).sum();

    let y_sdust = sdust_mask(seq_y);
    let mut seed_scores: HashMap<usize, usize> = HashMap::new();
    let mut node_hit_bounds: HashMap<usize, (usize, usize)> = HashMap::new();
    let long_node_thresh = (seq_y.len().saturating_mul(2)).max(4_096);
    let diag_band = (seq_y.len() / 4).clamp(128, 1_024).cast_signed();
    for (mask, ry, node_indices) in [
        (MASK_W10, false, &index.w10),
        (MASK_RY18, true, &index.ry18),
    ] {
        if ry && seed_scores.values().copied().max().unwrap_or(0) >= 2 {
            break;
        }
        if let Ok(q_idx) =
            FlatSeedIndex::build_with_sdust(seq_y, mask, 8, false, ry, Some(&y_sdust))
            && !q_idx.keys.is_empty()
        {
            for (nid, nd) in nodes.iter().enumerate() {
                let Some(Some(n_idx)) = node_indices.get(nid) else {
                    continue;
                };
                let min_hits = if nd.seq.len() > 4_000 { 3 } else { 1 };
                let shared = q_idx.shared_unique_keys(n_idx);
                if shared >= min_hits {
                    *seed_scores.entry(nid).or_insert(0) += shared;
                    if nd.seq.len() > long_node_thresh && !node_hit_bounds.contains_key(&nid) {
                        let (fwd_hits, rc_hits) = n_idx.intersect(&q_idx);
                        let hit_buf: &[(usize, usize)] = if fwd_hits.is_empty() {
                            &rc_hits
                        } else {
                            &fwd_hits
                        };
                        if let Some(bounds) =
                            densest_hit_projected_bounds(hit_buf, seq_y.len(), diag_band)
                        {
                            node_hit_bounds.insert(nid, bounds);
                        }
                    }
                }
            }
        }
    }

    if seed_scores.is_empty() && seq_y.len() >= 18 {
        let q_aa5 = aa5_keys(seq_y);
        if !q_aa5.is_empty() {
            for (nid, n_aa5) in index.aa5.iter().take(nodes.len()).enumerate() {
                if !n_aa5.is_empty() {
                    let shared = sorted_u16_intersection_count(&q_aa5, n_aa5);
                    if shared > 0 {
                        *seed_scores.entry(nid).or_insert(0) += shared;
                    }
                }
            }
        }
    }

    if seed_scores.is_empty() {
        if total_bp <= 50_000 {
            return align_seq_to_dag_mdl(nodes, seq_y, params, min_block_bits);
        }
        return GraphAlignment {
            path_nodes: Vec::new(),
            kept_pairs: Vec::new(),
            bubbles: usize::from(seq_y.len() >= 6),
            score_bits: 0.0,
        };
    }

    // Expand backward and forward reachability from seeded nodes within a local
    // sequence budget so large pangenome DAGs do not pull in whole chromosomes.
    let budget_bp = (seq_y.len() / 2).max(128);
    let mut succs: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
    for (u, nd) in nodes.iter().enumerate() {
        for &p in &nd.preds {
            if p < u {
                succs[p].push(u);
            }
        }
    }

    let max_seed_score = seed_scores.values().copied().max().unwrap_or(0);
    let anchor_floor = if nodes.len() > 64 && max_seed_score >= 4 {
        (max_seed_score / 4).max(2)
    } else {
        1
    };
    let mut dist_back: Vec<usize> = vec![usize::MAX; nodes.len()];
    let mut dist_fwd: Vec<usize> = vec![usize::MAX; nodes.len()];
    let chained_anchors = if nodes.len() > 64 {
        select_chained_dag_anchors(nodes, &seed_scores, anchor_floor, budget_bp)
    } else {
        Vec::new()
    };
    if chained_anchors.len() >= 2 {
        for &s in &chained_anchors {
            dist_back[s] = 0;
            dist_fwd[s] = 0;
        }
    } else {
        for (&s, &sc) in &seed_scores {
            if sc >= anchor_floor {
                dist_back[s] = 0;
                dist_fwd[s] = 0;
            }
        }
    }
    for u in (0..nodes.len()).rev() {
        if dist_back[u] < budget_bp {
            for &p in &nodes[u].preds {
                if p < u {
                    let cand = dist_back[u].saturating_add(nodes[p].seq.len());
                    if cand <= budget_bp && cand < dist_back[p] {
                        dist_back[p] = cand;
                    }
                }
            }
        }
    }
    for u in 0..nodes.len() {
        if dist_fwd[u] < budget_bp {
            for &s in &succs[u] {
                let cand = dist_fwd[u].saturating_add(nodes[s].seq.len());
                if cand <= budget_bp && cand < dist_fwd[s] {
                    dist_fwd[s] = cand;
                }
            }
        }
    }
    let anchor_count = dist_back.iter().filter(|&&d| d == 0).count();
    let use_bridge_only = nodes.len() > 64 && anchor_count >= 2;
    let flank_budget = (budget_bp / 4).max(64);
    let active_nodes: Vec<usize> = (0..nodes.len())
        .filter(|&u| {
            if use_bridge_only {
                dist_back[u] == 0
                    || (dist_back[u] <= budget_bp && dist_fwd[u] <= budget_bp)
                    || dist_back[u] <= flank_budget
                    || dist_fwd[u] <= flank_budget
            } else {
                dist_back[u] <= budget_bp || dist_fwd[u] <= budget_bp
            }
        })
        .collect();

    // Partition active_nodes into weakly connected components so an unrelated
    // decoy hit never distorts `total_max` banding for the true component.
    let mut active_pos: HashMap<usize, usize> = HashMap::with_capacity(active_nodes.len());
    for (idx, &u) in active_nodes.iter().enumerate() {
        active_pos.insert(u, idx);
    }
    let mut undir: Vec<Vec<usize>> = vec![Vec::new(); active_nodes.len()];
    for (idx, &u) in active_nodes.iter().enumerate() {
        for &p in &nodes[u].preds {
            if let Some(&p_idx) = active_pos.get(&p) {
                undir[idx].push(p_idx);
                undir[p_idx].push(idx);
            }
        }
    }
    let mut visited = vec![false; active_nodes.len()];
    let mut components: Vec<(usize, usize, Vec<usize>)> = Vec::new();
    for start in 0..active_nodes.len() {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let mut stack = vec![start];
        let mut comp = Vec::new();
        let mut comp_score = 0usize;
        while let Some(cur) = stack.pop() {
            let orig = active_nodes[cur];
            comp.push(orig);
            comp_score += seed_scores.get(&orig).copied().unwrap_or(0);
            for &nb in &undir[cur] {
                if !visited[nb] {
                    visited[nb] = true;
                    stack.push(nb);
                }
            }
        }
        comp.sort_unstable();
        let sliced_bp: usize = comp
            .iter()
            .map(|&orig| {
                let is_pred = succs[orig].iter().any(|s| seed_scores.contains_key(s));
                let is_succ = nodes[orig]
                    .preds
                    .iter()
                    .any(|p| seed_scores.contains_key(p));
                let (s, e) = compute_node_slice(
                    orig,
                    &nodes[orig],
                    &node_hit_bounds,
                    seq_y.len(),
                    is_pred,
                    is_succ,
                );
                e - s
            })
            .sum();
        if sliced_bp <= 200_000 {
            components.push((comp_score, sliced_bp, comp));
        }
    }
    components.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));

    let mut best_aln = GraphAlignment {
        path_nodes: Vec::new(),
        kept_pairs: Vec::new(),
        bubbles: usize::from(seq_y.len() >= 6),
        score_bits: 0.0,
    };
    for (_, _, comp) in components.into_iter().take(2) {
        let mut old_to_new: HashMap<usize, usize> = HashMap::with_capacity(comp.len());
        let mut new_to_old: Vec<usize> = Vec::with_capacity(comp.len());
        let mut slice_ranges: Vec<(usize, usize)> = Vec::with_capacity(comp.len());
        let mut sub_nodes: Vec<GraphNode> = Vec::with_capacity(comp.len());
        for old_id in comp {
            let nd = &nodes[old_id];
            let new_id = sub_nodes.len();
            let is_pred = succs[old_id].iter().any(|s| seed_scores.contains_key(s));
            let is_succ = nd.preds.iter().any(|p| seed_scores.contains_key(p));
            let (s, e) =
                compute_node_slice(old_id, nd, &node_hit_bounds, seq_y.len(), is_pred, is_succ);
            old_to_new.insert(old_id, new_id);
            new_to_old.push(old_id);
            slice_ranges.push((s, e));
            let mapped_preds = if s == 0 {
                nd.preds
                    .iter()
                    .filter_map(|&p| {
                        let p_new = *old_to_new.get(&p)?;
                        let (_, p_end) = slice_ranges[p_new];
                        (p_end == nodes[p].seq.len()).then_some(p_new)
                    })
                    .collect()
            } else {
                Vec::new()
            };
            sub_nodes.push(GraphNode {
                name: nd.name.clone(),
                seq: nd.seq[s..e].to_vec(),
                preds: mapped_preds,
            });
        }
        let mut aln = align_seq_to_dag_mdl(&sub_nodes, seq_y, params, min_block_bits);
        for nid in &mut aln.path_nodes {
            *nid = new_to_old[*nid];
        }
        for (nid, node_pos, _) in &mut aln.kept_pairs {
            let off = slice_ranges[*nid].0;
            *nid = new_to_old[*nid];
            *node_pos += off;
        }
        if aln.kept_pairs.len() > best_aln.kept_pairs.len()
            || (aln.kept_pairs.len() == best_aln.kept_pairs.len()
                && aln.score_bits > best_aln.score_bits)
        {
            best_aln = aln;
        }
    }
    best_aln
}

#[cfg(test)]
#[path = "dag_seed_tests.rs"]
mod tests;
