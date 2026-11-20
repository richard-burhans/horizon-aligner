//! Seeded subgraph selection and component partitioning for [`align_seq_to_dag_mdl`].

use super::dag::{GraphAlignment, GraphNode, align_seq_to_dag_mdl};
use crate::mdl::dp::HmmParams;
use crate::seeding::{FlatSeedIndex, MASK_RY20, MASK_W10, sdust_mask};
use std::collections::HashMap;

/// Pre-built per-node spaced-seed indices (`MASK_W10` and `MASK_RY20`) for repeated
/// queries against the same sequence DAG.
#[derive(Debug, Clone)]
pub struct DagSeedIndex {
    w10: Vec<Option<FlatSeedIndex>>,
    ry20: Vec<Option<FlatSeedIndex>>,
}

impl DagSeedIndex {
    /// Build `MASK_W10` and `MASK_RY20` seed indices once for all nodes in `nodes`.
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
        let ry20 = nodes
            .iter()
            .map(|nd| {
                (nd.seq.len() >= MASK_RY20.len())
                    .then(|| {
                        FlatSeedIndex::build_with_sdust(&nd.seq, MASK_RY20, 8, false, true, None)
                            .ok()
                    })
                    .flatten()
            })
            .collect();
        Self { w10, ry20 }
    }
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
    (nodes.len() <= 8 && total_bp <= 50_000)
        .then(|| align_seq_to_dag_mdl(nodes, seq_y, params, min_block_bits))
}

/// Seeded DAG alignment: builds a [`DagSeedIndex`] and delegates to [`align_seq_to_dag_seeded_indexed`].
#[must_use]
pub fn align_seq_to_dag_seeded(
    nodes: &[GraphNode],
    seq_y: &[u8],
    params: &HmmParams,
    min_block_bits: f64,
) -> GraphAlignment {
    if let Some(aln) = try_small_or_invalid_dag(nodes, seq_y, params, min_block_bits) {
        return aln;
    }
    let index = DagSeedIndex::build(nodes);
    align_seq_to_dag_seeded_indexed(nodes, &index, seq_y, params, min_block_bits)
}

/// Seeded DAG alignment using a pre-built [`DagSeedIndex`] to isolate the
/// connected/reachable DAG subgraph matching `seq_y`, then runs [`align_seq_to_dag_mdl`].
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
    for (mask, ry, node_indices) in [
        (MASK_W10, false, &index.w10),
        (MASK_RY20, true, &index.ry20),
    ] {
        if ry && !seed_scores.is_empty() {
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
    for (&s, &sc) in &seed_scores {
        if sc >= anchor_floor {
            dist_back[s] = 0;
            dist_fwd[s] = 0;
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
    let active_nodes: Vec<usize> = (0..nodes.len())
        .filter(|&u| {
            if use_bridge_only {
                dist_back[u] == 0 || (dist_back[u] <= budget_bp && dist_fwd[u] <= budget_bp)
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
        let mut comp_bp = 0usize;
        while let Some(cur) = stack.pop() {
            let orig = active_nodes[cur];
            comp.push(orig);
            comp_score += seed_scores.get(&orig).copied().unwrap_or(0);
            comp_bp += nodes[orig].seq.len();
            for &nb in &undir[cur] {
                if !visited[nb] {
                    visited[nb] = true;
                    stack.push(nb);
                }
            }
        }
        comp.sort_unstable();
        if comp_bp <= 200_000 {
            components.push((comp_score, comp_bp, comp));
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
        let mut sub_nodes: Vec<GraphNode> = Vec::with_capacity(comp.len());
        for old_id in comp {
            let nd = &nodes[old_id];
            let new_id = sub_nodes.len();
            old_to_new.insert(old_id, new_id);
            new_to_old.push(old_id);
            let mapped_preds = nd
                .preds
                .iter()
                .filter_map(|p| old_to_new.get(p).copied())
                .collect();
            sub_nodes.push(GraphNode {
                name: nd.name.clone(),
                seq: nd.seq.clone(),
                preds: mapped_preds,
            });
        }
        let mut aln = align_seq_to_dag_mdl(&sub_nodes, seq_y, params, min_block_bits);
        for nid in &mut aln.path_nodes {
            *nid = new_to_old[*nid];
        }
        for (nid, _, _) in &mut aln.kept_pairs {
            *nid = new_to_old[*nid];
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
mod tests {
    use super::super::dag::GraphNode;
    use super::*;

    #[test]
    fn seeded_dag_aligns_ignoring_disconnected_decoy_component() {
        let mut nodes = Vec::new();
        // 6 unrelated decoy nodes
        for idx in 0..6 {
            nodes.push(GraphNode {
                name: format!("decoy_{idx}"),
                seq: b"CCCCCCCCCCCCCCCCCCCCCCCCC".to_vec(),
                preds: if idx == 0 { vec![] } else { vec![idx - 1] },
            });
        }
        // Real 4-node diamond component (indices 6..10)
        nodes.push(GraphNode {
            name: "r0".into(),
            seq: b"ACTGACTGACTGACTGACTGACTGACTG".to_vec(),
            preds: vec![],
        });
        nodes.push(GraphNode {
            name: "r1".into(),
            seq: b"TTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec(),
            preds: vec![6],
        });
        nodes.push(GraphNode {
            name: "r2".into(),
            seq: b"GATTACAGATTACAGATTACAGATTACA".to_vec(),
            preds: vec![6],
        });
        nodes.push(GraphNode {
            name: "r3".into(),
            seq: b"CATGCATGCATGCATGCATGCATGCATG".to_vec(),
            preds: vec![7, 8],
        });
        let query =
            b"ACTGACTGACTGACTGACTGACTGACTGGATTACAGATTACAGATTACAGATTACACATGCATGCATGCATGCATGCATGCATG";
        let aln = align_seq_to_dag_seeded(&nodes, query, &HmmParams::default(), 10.0);
        assert_eq!(aln.path_nodes, vec![6, 8, 9]);
        assert_eq!(aln.kept_pairs.len(), query.len());
    }

    #[test]
    fn dag_handles_query_insertion_spanning_node_boundary() {
        let nodes = vec![
            GraphNode {
                name: "n0".into(),
                seq: b"ACTGACTGACTGACTGACTGACTGACTGACTG".to_vec(),
                preds: vec![],
            },
            GraphNode {
                name: "n1".into(),
                seq: b"GATTACAGATTACAGATTACAGATTACAGATT".to_vec(),
                preds: vec![0],
            },
        ];
        let query = b"ACTGACTGACTGACTGACTGACTGACTGACTGCCCCCC\
                      GATTACAGATTACAGATTACAGATTACAGATT";
        let aln = align_seq_to_dag_seeded(&nodes, query, &HmmParams::default(), 10.0);
        assert_eq!(aln.path_nodes, vec![0, 1]);
        assert!(aln.kept_pairs.iter().any(|&(n, _, _)| n == 0));
        assert!(aln.kept_pairs.iter().any(|&(n, _, _)| n == 1));
    }

    #[test]
    fn seeded_dag_abstains_without_oom_on_large_unseeded_graph() {
        let nodes: Vec<GraphNode> = (0..12)
            .map(|idx| GraphNode {
                name: format!("big_{idx}"),
                seq: vec![b'A'; 5_000],
                preds: if idx == 0 { vec![] } else { vec![idx - 1] },
            })
            .collect();
        let query = vec![b'N'; 400];
        let aln = align_seq_to_dag_seeded(&nodes, &query, &HmmParams::default(), 10.0);
        assert_eq!(aln.path_nodes, [] as [usize; 0]);
        assert_eq!(aln.kept_pairs, [] as [(usize, usize, usize); 0]);
        assert_eq!(aln.bubbles, 1);
    }
}
