//! E-graph tandem duplication resolution and overlap-clustered MDL paralog pruning.

use crate::egraph;
use crate::mdl::RunAlignment;
use crate::netting::partition_colinear_run_groups;

/// Run e-graph duplication/inversion collapse on a sequence of signed block labels.
fn resolve_lcb_duplications(blocks: &[(&str, i8)]) -> Option<egraph::Expr> {
    let mut eg = egraph::EGraph::new();
    let root = egraph::build_seq_expr(&mut eg, blocks)?;
    egraph::saturate_synteny(&mut eg, 4, 800);
    egraph::extract_min_cost(&eg, root).map(|(_, expr)| expr)
}

pub(crate) fn collect_windowed_duplications(
    tokens: &[(String, i8)],
    out: &mut std::collections::BTreeSet<String>,
) {
    // Use 50% overlapping 12-token windows within each contig so duplications
    // straddling a window boundary are still detected without mixing contigs.
    let mut start = 0usize;
    while start < tokens.len() {
        let end = (start + 12).min(tokens.len());
        let win = &tokens[start..end];
        if win.len() >= 2 {
            let refs: Vec<(&str, i8)> = win.iter().map(|(s, st)| (s.as_str(), *st)).collect();
            if let Some(expr) = resolve_lcb_duplications(&refs) {
                for d in egraph::collect_dups(&expr) {
                    out.insert(d);
                }
            }
        }
        if end == tokens.len() {
            break;
        }
        start += 6;
    }
}

/// Score a colinear run-group by strand-aware ACGT match identity (`+2` exact match, `-1` mismatch).
fn score_run_group_identity(
    grp: &[crate::mdl::AlignedRun],
    seq_x: &[u8],
    seq_y: &[u8],
) -> (i64, usize) {
    let mut score: i64 = 0;
    let mut bp: usize = 0;
    for r in grp {
        bp = bp.saturating_add(r.len);
        for k in 0..r.len {
            let x_pos = r.x_start.saturating_add(k);
            let y_pos = if r.is_rc {
                let Some(p) = r.y_start.checked_sub(k) else {
                    break;
                };
                p
            } else {
                r.y_start.saturating_add(k)
            };
            let (Some(&raw_ref), Some(&raw_qry)) = (seq_x.get(x_pos), seq_y.get(y_pos)) else {
                break;
            };
            let ref_base = raw_ref.to_ascii_uppercase();
            let qry_base = if r.is_rc {
                crate::seeding::complement_base(raw_qry).to_ascii_uppercase()
            } else {
                raw_qry.to_ascii_uppercase()
            };
            if matches!(ref_base, b'A' | b'C' | b'G' | b'T')
                && matches!(qry_base, b'A' | b'C' | b'G' | b'T')
            {
                if ref_base == qry_base {
                    score += 2;
                } else {
                    score -= 1;
                }
            }
        }
    }
    (score, bp)
}

/// Sequence-aware e-graph duplication resolution output for a contig-pair alignment.
#[derive(Debug, Clone, Default)]
pub struct ResolvedDuplications {
    /// Human-readable labels for collapsed loci (e.g. `C0@990-1910(2x)`).
    pub dup_labels: Vec<String>,
    /// Cluster indices collapsed by `MdlCostModel` (`Dup` or `Repeat`).
    pub dup_cluster_ids: Vec<usize>,
    /// Reference-overlap cluster index assigned to each colinear run-group.
    pub group_clusters: Vec<usize>,
}

/// Measure strand-aware ACGT mismatch divergence `(mismatches, compared_acgt_bp)` for a run-group.
fn run_group_divergence(
    grp: &[crate::mdl::AlignedRun],
    seq_x: &[u8],
    seq_y: &[u8],
) -> (usize, usize) {
    let mut mismatches = 0usize;
    let mut compared = 0usize;
    for r in grp {
        for k in 0..r.len {
            let x_pos = r.x_start.saturating_add(k);
            let y_pos = if r.is_rc {
                let Some(p) = r.y_start.checked_sub(k) else {
                    break;
                };
                p
            } else {
                r.y_start.saturating_add(k)
            };
            let (Some(&raw_ref), Some(&raw_qry)) = (seq_x.get(x_pos), seq_y.get(y_pos)) else {
                break;
            };
            let ref_base = raw_ref.to_ascii_uppercase();
            let qry_base = if r.is_rc {
                crate::seeding::complement_base(raw_qry).to_ascii_uppercase()
            } else {
                raw_qry.to_ascii_uppercase()
            };
            if matches!(ref_base, b'A' | b'C' | b'G' | b'T')
                && matches!(qry_base, b'A' | b'C' | b'G' | b'T')
            {
                compared += 1;
                if ref_base != qry_base {
                    mismatches += 1;
                }
            }
        }
    }
    (mismatches, compared)
}

/// Cluster colinear run-groups by reference span overlap (`overlap >= 50%` of shorter span),
/// returning `(cluster_spans, group_cluster_ids)`.
///
/// Uses an active-cluster sweep when `groups` are sorted by ascending reference start (the
/// standard output of `partition_colinear_run_groups`), reducing runtime from `O(G^2)` to
/// `O(G * D)` where `D` is local overlap depth while preserving identical cluster assignments.
fn cluster_run_groups_by_ref_overlap(
    groups: &[(Vec<crate::mdl::AlignedRun>, bool)],
) -> (Vec<(usize, usize)>, Vec<usize>) {
    let mut cluster_spans: Vec<(usize, usize)> = Vec::new();
    let mut group_clusters: Vec<usize> = Vec::with_capacity(groups.len());
    let sorted_by_x = groups
        .windows(2)
        .all(|w| w[0].0[0].x_start <= w[1].0[0].x_start);
    let mut active_cids: Vec<usize> = Vec::new();

    for (grp, _) in groups {
        let g_lo = grp[0].x_start;
        let g_hi = grp.last().map_or(g_lo, |r| r.x_start.saturating_add(r.len));
        let g_len = g_hi.saturating_sub(g_lo).max(1);
        let mut best_cid: Option<(usize, usize)> = None;

        if sorted_by_x {
            // Since `g_lo` is non-decreasing, any cluster ending at or before `g_lo` can never
            // overlap this or any later group.
            active_cids.retain(|&cid| cluster_spans[cid].1 > g_lo);
            for &cid in &active_cids {
                let (c_lo, c_hi) = cluster_spans[cid];
                let ov = g_hi.min(c_hi).saturating_sub(g_lo.max(c_lo));
                if ov == 0 {
                    continue;
                }
                let c_len = c_hi.saturating_sub(c_lo).max(1);
                let min_len = g_len.min(c_len);
                if ov.saturating_mul(2) >= min_len
                    && best_cid.is_none_or(|(_, best_ov)| ov > best_ov)
                {
                    best_cid = Some((cid, ov));
                }
            }
        } else {
            for (cid, &(c_lo, c_hi)) in cluster_spans.iter().enumerate() {
                let ov = g_hi.min(c_hi).saturating_sub(g_lo.max(c_lo));
                if ov == 0 {
                    continue;
                }
                let c_len = c_hi.saturating_sub(c_lo).max(1);
                let min_len = g_len.min(c_len);
                if ov.saturating_mul(2) >= min_len
                    && best_cid.is_none_or(|(_, best_ov)| ov > best_ov)
                {
                    best_cid = Some((cid, ov));
                }
            }
        }

        if let Some((cid, _)) = best_cid {
            let span = &mut cluster_spans[cid];
            span.0 = span.0.min(g_lo);
            span.1 = span.1.max(g_hi);
            group_clusters.push(cid);
        } else {
            let cid = cluster_spans.len();
            cluster_spans.push((g_lo, g_hi));
            if sorted_by_x {
                active_cids.push(cid);
            }
            group_clusters.push(cid);
        }
    }
    (cluster_spans, group_clusters)
}

/// Cluster colinear run-groups by reference coordinate overlap (`overlap >= 50%` of shorter span),
/// evaluate real per-cluster `(length, divergence)` under `MdlCostModel` with `k`-copy e-graph
/// saturation, and return collapsed paralog clusters.
#[must_use]
pub fn resolve_alignment_duplications_mdl(
    aln: &RunAlignment,
    seq_x: &[u8],
    seq_y: &[u8],
) -> ResolvedDuplications {
    let groups = partition_colinear_run_groups(&aln.aligned_runs);
    if groups.len() < 2 {
        return ResolvedDuplications {
            dup_labels: Vec::new(),
            dup_cluster_ids: Vec::new(),
            group_clusters: (0..groups.len()).collect(),
        };
    }

    let (cluster_spans, group_clusters) = cluster_run_groups_by_ref_overlap(&groups);

    // Step 2: Per-cluster empirical (mean_len_bp, divergence) for MdlCostModel.
    let n_clusters = cluster_spans.len();
    let mut cluster_bp_sum = vec![0usize; n_clusters];
    let mut cluster_cnt = vec![0usize; n_clusters];
    // Track the two lowest-divergence copies per cluster so one noisy 3rd paralog does not
    // poison a true high-identity 2-copy duplication.
    let mut cluster_div_min1 = vec![1.0_f64; n_clusters];
    let mut cluster_div_min2 = vec![1.0_f64; n_clusters];

    let mut tokens: Vec<(String, i8)> = Vec::with_capacity(groups.len());
    for ((grp, is_rev), &cid) in groups.iter().zip(&group_clusters) {
        let bp: usize = grp.iter().map(|r| r.len).sum();
        let (mm, cmp) = run_group_divergence(grp, seq_x, seq_y);
        #[allow(clippy::cast_precision_loss)]
        let div = if cmp > 0 {
            (mm as f64 / cmp as f64).clamp(1e-6, 0.75)
        } else {
            0.75
        };
        cluster_bp_sum[cid] += bp;
        cluster_cnt[cid] += 1;
        if div < cluster_div_min1[cid] {
            cluster_div_min2[cid] = cluster_div_min1[cid];
            cluster_div_min1[cid] = div;
        } else if div < cluster_div_min2[cid] {
            cluster_div_min2[cid] = div;
        }
        tokens.push((format!("C{cid}"), if *is_rev { -1 } else { 1 }));
    }

    let mut block_stats: std::collections::HashMap<String, (usize, f64)> =
        std::collections::HashMap::with_capacity(n_clusters);
    let mut saturated_clusters: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for cid in 0..n_clusters {
        let mean_bp = (cluster_bp_sum[cid] / cluster_cnt[cid].max(1)).max(1);
        // Pairwise copy divergence bound from the two closest copies clamped to JC69 0.75.
        let est_div = if cluster_cnt[cid] >= 2 {
            (cluster_div_min1[cid] + cluster_div_min2[cid]).clamp(1e-6, 0.75)
        } else {
            cluster_div_min1[cid]
        };
        if est_div >= 0.72 {
            saturated_clusters.insert(cid);
        }
        block_stats.insert(format!("C{cid}"), (mean_bp, est_div));
    }

    // Step 3: Windowed e-graph saturation + sequence-aware MdlCostModel extraction.
    let mut collapsed_mults: std::collections::BTreeMap<usize, usize> =
        std::collections::BTreeMap::new();
    let mut start = 0usize;
    while start < tokens.len() {
        let end = (start + 12).min(tokens.len());
        let win = &tokens[start..end];
        if win.len() >= 2 {
            let refs: Vec<(&str, i8)> = win.iter().map(|(s, st)| (s.as_str(), *st)).collect();
            let mut eg = egraph::EGraph::new();
            if let Some(root) = egraph::build_seq_expr(&mut eg, &refs) {
                let mut model = egraph::MdlCostModel::from_egraph(&eg);
                for (lbl, _) in &refs {
                    if let Some(&stats) = block_stats.get(*lbl) {
                        model.block_stats.insert((*lbl).to_owned(), stats);
                    }
                }
                egraph::saturate_synteny(&mut eg, 6, 1024);
                if let Some((_, expr)) = egraph::extract_min_mdl_cost(&eg, root, &model) {
                    for (lbl, mult) in egraph::collect_dup_multiplicities(&expr) {
                        if let Some(cid) =
                            lbl.strip_prefix('C').and_then(|s| s.parse::<usize>().ok())
                            && !saturated_clusters.contains(&cid)
                        {
                            let entry = collapsed_mults.entry(cid).or_insert(mult);
                            *entry = (*entry).max(mult);
                        }
                    }
                }
            }
        }
        if end == tokens.len() {
            break;
        }
        start += 6;
    }

    let mut dup_cluster_ids = Vec::with_capacity(collapsed_mults.len());
    let mut dup_labels = Vec::with_capacity(collapsed_mults.len());
    for (cid, mult) in collapsed_mults {
        let (lo, hi) = cluster_spans[cid];
        dup_cluster_ids.push(cid);
        dup_labels.push(format!("C{cid}@{lo}-{hi}({mult}x)"));
    }

    ResolvedDuplications {
        dup_labels,
        dup_cluster_ids,
        group_clusters,
    }
}

/// Prune redundant paralog run-groups for every cluster collapsed by
/// [`resolve_alignment_duplications_mdl`], keeping the single highest-identity ortholog per cluster.
#[must_use]
pub fn prune_resolved_dup_clusters(
    aln: &RunAlignment,
    resolved: &ResolvedDuplications,
    seq_x: &[u8],
    seq_y: &[u8],
) -> RunAlignment {
    if resolved.dup_cluster_ids.is_empty() || aln.aligned_runs.is_empty() {
        return aln.clone();
    }
    let dup_set: std::collections::HashSet<usize> =
        resolved.dup_cluster_ids.iter().copied().collect();
    let groups = partition_colinear_run_groups(&aln.aligned_runs);
    if groups.len() != resolved.group_clusters.len() {
        return aln.clone();
    }
    let mut best_for_cluster: std::collections::HashMap<usize, (usize, i64, usize, usize, usize)> =
        std::collections::HashMap::new();

    for (idx, ((grp, _), &cid)) in groups.iter().zip(&resolved.group_clusters).enumerate() {
        if !dup_set.contains(&cid) {
            continue;
        }
        let (id_score, bp) = score_run_group_identity(grp, seq_x, seq_y);
        let ref_start = grp[0].x_start;
        let qry_start = grp[0].y_start;
        let replace = match best_for_cluster.get(&cid) {
            None => true,
            Some(&(_, cur_score, cur_bp, best_ref_start, best_qry_start)) => {
                (
                    id_score,
                    bp,
                    std::cmp::Reverse(ref_start),
                    std::cmp::Reverse(qry_start),
                ) > (
                    cur_score,
                    cur_bp,
                    std::cmp::Reverse(best_ref_start),
                    std::cmp::Reverse(best_qry_start),
                )
            }
        };
        if replace {
            best_for_cluster.insert(cid, (idx, id_score, bp, ref_start, qry_start));
        }
    }

    let mut kept_runs = Vec::with_capacity(aln.aligned_runs.len());
    let mut lcbs_fwd = 0usize;
    let mut lcbs_rc = 0usize;
    for (idx, ((grp, is_rev), &cid)) in groups.into_iter().zip(&resolved.group_clusters).enumerate()
    {
        let keep = if dup_set.contains(&cid) {
            best_for_cluster
                .get(&cid)
                .is_some_and(|&(best_idx, _, _, _, _)| best_idx == idx)
        } else {
            true
        };
        if keep {
            if is_rev {
                lcbs_rc += 1;
            } else {
                lcbs_fwd += 1;
            }
            kept_runs.extend(grp);
        }
    }
    let total_pairs = kept_runs.iter().map(|r| r.len).sum();
    RunAlignment {
        aligned_runs: kept_runs,
        total_pairs,
        lcbs_fwd,
        lcbs_rc,
        macro_bubbles: aln.macro_bubbles,
        internal_bubbles: aln.internal_bubbles,
    }
}
