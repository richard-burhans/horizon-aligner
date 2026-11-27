//! GFA-to-DAG loader and unit tests for `horizon_wga --graph`.

use anyhow::Result;
use horizon_aligner::gfa::{GfaRecord, stream_gfa_records};
#[cfg(test)]
use horizon_aligner::{DagSeedIndex, estimate_dag_hmm_params};
use horizon_aligner::{GraphNode, revcomp};
use std::collections::{HashMap, VecDeque};

type OrientedSegs = Vec<(String, Vec<u8>)>;
type DirectedEdges = Vec<(usize, usize)>;

#[allow(clippy::too_many_lines)]
pub(crate) fn load_gfa_dag(path: &std::path::Path) -> Result<Vec<GraphNode>> {
    let mut segs = Vec::<(String, Vec<u8>)>::new();
    let mut name_to_raw = HashMap::<String, usize>::new();
    let mut raw_edges = Vec::<(String, u8, String, u8)>::new();

    stream_gfa_records(path, |rec| {
        match rec {
            GfaRecord::Segment(s) if s.sequence != b"*" => {
                anyhow::ensure!(
                    !name_to_raw.contains_key(&s.name),
                    "duplicate GFA segment '{}'",
                    s.name
                );
                let idx = segs.len();
                name_to_raw.insert(s.name.clone(), idx);
                segs.push((s.name, s.sequence));
            }
            GfaRecord::Link(l) => {
                raw_edges.push((l.from, l.from_orient, l.to, l.to_orient));
            }
            _ => {}
        }
        Ok(std::ops::ControlFlow::Continue(()))
    })?;

    let has_rev = raw_edges
        .iter()
        .any(|(_, fo, _, to)| *fo == b'-' || *to == b'-');

    let (oriented_segs, directed_edges): (OrientedSegs, DirectedEdges) = if has_rev {
        let mut o_segs = Vec::with_capacity(segs.len() * 2);
        for (name, seq) in &segs {
            o_segs.push((format!("{name}+"), seq.clone()));
            o_segs.push((format!("{name}-"), revcomp(seq)));
        }
        let orient_idx = |raw_i: usize, orient: u8| -> usize {
            if orient == b'-' {
                raw_i * 2 + 1
            } else {
                raw_i * 2
            }
        };
        let flip = |o: u8| -> u8 { if o == b'-' { b'+' } else { b'-' } };
        let mut d_edges = Vec::with_capacity(raw_edges.len() * 2);
        let mut skipped_links = 0usize;
        for (u_name, uo, v_name, vo) in &raw_edges {
            if let (Some(&u_raw), Some(&v_raw)) = (name_to_raw.get(u_name), name_to_raw.get(v_name))
            {
                let u = orient_idx(u_raw, *uo);
                let v = orient_idx(v_raw, *vo);
                if u != v {
                    d_edges.push((u, v));
                }
                let rv = orient_idx(v_raw, flip(*vo));
                let ru = orient_idx(u_raw, flip(*uo));
                if rv != ru {
                    d_edges.push((rv, ru));
                }
            } else {
                skipped_links += 1;
            }
        }
        if skipped_links > 0 {
            eprintln!(
                "horizon_wga: warning: skipped {skipped_links} GFA link(s) referencing missing or `*`-sequence segments"
            );
        }
        (o_segs, d_edges)
    } else {
        let mut d_edges = Vec::with_capacity(raw_edges.len());
        let mut skipped_links = 0usize;
        for (u_name, _, v_name, _) in &raw_edges {
            if let (Some(&u), Some(&v)) = (name_to_raw.get(u_name), name_to_raw.get(v_name)) {
                if u != v {
                    d_edges.push((u, v));
                }
            } else {
                skipped_links += 1;
            }
        }
        if skipped_links > 0 {
            eprintln!(
                "horizon_wga: warning: skipped {skipped_links} GFA link(s) referencing missing or `*`-sequence segments"
            );
        }
        (segs, d_edges)
    };

    let n = oriented_segs.len();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut indeg: Vec<usize> = vec![0; n];

    for (u, v) in directed_edges {
        adj[u].push(v);
    }
    for nbrs in &mut adj {
        nbrs.sort_unstable();
        nbrs.dedup();
        for &v in nbrs.iter() {
            indeg[v] += 1;
        }
    }

    let mut q: VecDeque<usize> = (0..n).filter(|&i| indeg[i] == 0).collect();
    let mut order = Vec::with_capacity(n);
    while let Some(u) = q.pop_front() {
        order.push(u);
        for &v in &adj[u] {
            indeg[v] = indeg[v].saturating_sub(1);
            if indeg[v] == 0 {
                q.push_back(v);
            }
        }
    }
    // Append any remaining nodes (if cycles existed) in input order
    if order.len() < n {
        eprintln!(
            "warning: cyclic GFA detected ({} of {} oriented nodes in cycles); back-edges will be dropped to form a DAG",
            n - order.len(),
            n
        );
        let seen: std::collections::HashSet<usize> = order.iter().copied().collect();
        for i in 0..n {
            if !seen.contains(&i) {
                order.push(i);
            }
        }
    }

    let mut raw_to_topo = vec![0usize; n];
    for (topo_i, &raw_i) in order.iter().enumerate() {
        raw_to_topo[raw_i] = topo_i;
    }

    let mut nodes: Vec<GraphNode> = order
        .iter()
        .map(|&raw_i| GraphNode {
            name: oriented_segs[raw_i].0.clone(),
            seq: oriented_segs[raw_i].1.clone(),
            preds: Vec::new(),
        })
        .collect();

    for (u_raw, vs) in adj.iter().enumerate() {
        let u_topo = raw_to_topo[u_raw];
        for &v_raw in vs {
            let v_topo = raw_to_topo[v_raw];
            if u_topo < v_topo {
                nodes[v_topo].preds.push(u_topo);
            }
        }
    }
    Ok(nodes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn load_gfa_dag_handles_bidirected_minus_links() -> Result<()> {
        let dir = std::env::temp_dir().join(format!("horizon_wga_dag_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir)?;
        let gfa = dir.join("inv.gfa");
        // Seg a -> Seg b(-) -> Seg c(+)
        std::fs::write(
            &gfa,
            "H\tVN:Z:1.0\n\
             S\ta\tACTGACTGACTGACTGACTG\n\
             S\tb\tCCCCCCCCCCCCCCCCAAAA\n\
             S\tc\tGATTACAGATTACAGATTAC\n\
             L\ta\t+\tb\t-\t0M\n\
             L\tb\t-\tc\t+\t0M\n",
        )?;
        let nodes = load_gfa_dag(&gfa)?;
        assert_eq!(nodes.len(), 6);
        let b_rc = revcomp(b"CCCCCCCCCCCCCCCCAAAA");
        let mut query = b"ACTGACTGACTGACTGACTG".to_vec();
        query.extend_from_slice(&b_rc);
        query.extend_from_slice(b"GATTACAGATTACAGATTAC");
        let dag_idx = DagSeedIndex::build(&nodes);
        let aln = horizon_aligner::align_seq_to_dag_seeded_indexed(
            &nodes,
            &dag_idx,
            &query,
            &horizon_aligner::HmmParams::default(),
            10.0,
        );
        let path_names: Vec<&str> = aln
            .path_nodes
            .iter()
            .map(|&i| nodes[i].name.as_str())
            .collect();
        assert_eq!(path_names, vec!["a+", "b-", "c+"]);
        assert_eq!(aln.kept_pairs.len(), query.len());
        let _ = std::fs::remove_dir_all(&dir);
        Ok(())
    }

    #[test]
    fn estimate_dag_hmm_params_adapts_identity_and_gaps() {
        let short_nodes = vec![GraphNode {
            name: "s+".to_string(),
            seq: b"ACGTACGTACGT".to_vec(),
            preds: vec![],
        }];
        let def = horizon_aligner::HmmParams::default();
        let short_idx = DagSeedIndex::build(&short_nodes);
        let fallback = estimate_dag_hmm_params(&short_nodes, &short_idx, b"ACGTACGTACGT");
        assert!((fallback.identity_prior - def.identity_prior).abs() < 1e-9);

        let bases = *b"ACGT";
        let mut ref_seq = Vec::with_capacity(400);
        let mut state: u64 = 0x517c_c1b7_2722_0a95;
        for _ in 0..400 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            ref_seq.push(bases[((state >> 32) & 3) as usize]);
        }
        // Introduce 4 small deletions in query so dx != dy across seed pairs
        let mut query = Vec::new();
        for (idx, &b) in ref_seq.iter().enumerate() {
            if matches!(idx, 70 | 140 | 210 | 280 | 281 | 282) {
                continue;
            }
            query.push(b);
        }
        let nodes = vec![GraphNode {
            name: "n0+".to_string(),
            seq: ref_seq,
            preds: vec![],
        }];
        let dag_idx = DagSeedIndex::build(&nodes);
        let adapted = estimate_dag_hmm_params(&nodes, &dag_idx, &query);
        assert!(adapted.identity_prior >= 0.70);
        assert!(adapted.gap_open_bits >= -5.5 && adapted.gap_open_bits <= -4.0);
        assert!(adapted.gap_extend_bits >= -1.8 && adapted.gap_extend_bits <= -0.5);

        // R8 regression: a >=400 bp 100%-identical reverse-complement query must estimate high
        // identity (~0.72), not clamp to 0.60 from bogus forward canonical-hit comparisons.
        let rc_query = revcomp(&nodes[0].seq);
        let adapted_rc = estimate_dag_hmm_params(&nodes, &dag_idx, &rc_query);
        assert!(
            adapted_rc.identity_prior >= 0.70,
            "RC query identity_prior corrupted to {}",
            adapted_rc.identity_prior
        );
    }

    #[test]
    fn global_one_to_one_prefers_syntenic_ortholog_over_ectopic_paralog() {
        // R2 regression: identity-scored global netting prioritizes high-identity orthologs
        // over shorter or lower-identity ectopic repeats and returns exact unique_x/unique_y bp.
        let seq_x = vec![b'A'; 100];
        let seq_y_paralog = vec![b'A'; 30];
        let seq_y_syntenic = vec![b'A'; 100];
        let mut pairs = vec![
            (
                0,
                0,
                vec![horizon_aligner::AlignedRun {
                    x_start: 10,
                    y_start: 5,
                    len: 12,
                    is_rc: false,
                }],
            ),
            (
                0,
                1,
                vec![horizon_aligner::AlignedRun {
                    x_start: 0,
                    y_start: 0,
                    len: 80,
                    is_rc: false,
                }],
            ),
        ];
        let (ux, uy) = horizon_aligner::net_global_one_to_one_runs_scored(
            &mut pairs,
            &[seq_x.as_slice()],
            &[seq_y_paralog.as_slice(), seq_y_syntenic.as_slice()],
        );
        assert_eq!(ux, 80);
        assert_eq!(uy, 80);
        assert_eq!(pairs[1].2.len(), 1);
        assert_eq!(pairs[1].2[0].len, 80);
        assert_eq!(pairs[0].2.len(), 0);

        // Also test high-identity ortholog beats a longer low-identity ectopic repeat on another contig.
        let x2 = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(); // 32 bp
        let y_low = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec(); // 25% id on 32 bp (score -8)
        let y_high = b"ACGTACGTACGTACGTACGTACGT".to_vec(); // 100% id on 24 bp (score +48)
        let mut pairs2 = vec![
            (
                0,
                0,
                vec![horizon_aligner::AlignedRun {
                    x_start: 0,
                    y_start: 0,
                    len: 32,
                    is_rc: false,
                }],
            ),
            (
                0,
                1,
                vec![horizon_aligner::AlignedRun {
                    x_start: 0,
                    y_start: 0,
                    len: 24,
                    is_rc: false,
                }],
            ),
        ];
        let (ux2, uy2) = horizon_aligner::net_global_one_to_one_runs_scored(
            &mut pairs2,
            &[x2.as_slice()],
            &[y_low.as_slice(), y_high.as_slice()],
        );
        assert_eq!(pairs2[1].2.len(), 1);
        assert_eq!(pairs2[1].2[0].len, 24);
        // Remaining 8 bp tail [24..32) of y_low stays disjoint on X and Y (length >= 6).
        assert_eq!(ux2, 32);
        assert_eq!(uy2, 32);
    }
}
