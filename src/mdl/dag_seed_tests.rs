use super::super::dag::GraphNode;
use super::*;

#[test]
fn select_chained_dag_anchors_chains_long_adjacent_backbone_anchors() {
    // Two directly adjacent 500 bp anchors must chain even when max_bridge_bp = 128 < 500.
    let nodes = vec![
        GraphNode {
            name: "long_a".to_owned(),
            seq: vec![b'A'; 500],
            preds: Vec::new(),
        },
        GraphNode {
            name: "long_b".to_owned(),
            seq: vec![b'C'; 500],
            preds: vec![0],
        },
    ];
    let mut scores = HashMap::new();
    scores.insert(0, 10);
    scores.insert(1, 10);
    let chain = select_chained_dag_anchors(&nodes, &scores, 1, 128);
    assert_eq!(chain, vec![0, 1]);
}

#[test]
fn seeded_dag_chains_longest_topological_anchor_path_over_off_chain_decoy() {
    // Build a >64-node DAG with a shared root feeding two branches:
    // - Branch A (true 3-node homologous chain): n1 -> n2 -> n3
    // - Branch B (long decoy chain with 1 sporadic seed hit at the end): 70 unrelated nodes
    let seg1 = b"ACGTACGTTTGCATGCATGCATGCAACCGGTTAA".to_vec();
    let seg2 = b"GGCCAATTCCGGTTAACGTACGTACGTGGGGTTT".to_vec();
    let seg3 = b"TGCATGCAAACCGGTTCCGGTTAACGTACGTAAA".to_vec();
    let mut query = Vec::new();
    query.extend_from_slice(&seg1);
    query.extend_from_slice(&seg2);
    query.extend_from_slice(&seg3);

    let mut nodes = vec![GraphNode {
        name: "root".to_owned(),
        seq: b"CCCCCCCCCCCCCCCCCCCC".to_vec(),
        preds: Vec::new(),
    }];
    // True branch: nodes 1, 2, 3
    nodes.push(GraphNode {
        name: "true1".to_owned(),
        seq: seg1,
        preds: vec![0],
    });
    nodes.push(GraphNode {
        name: "true2".to_owned(),
        seq: seg2,
        preds: vec![1],
    });
    nodes.push(GraphNode {
        name: "true3".to_owned(),
        seq: seg3.clone(),
        preds: vec![2],
    });
    // Decoy branch: 70 unrelated nodes chained from root, ending with a duplicate of seg3
    for k in 0..70usize {
        let pred = if k == 0 { 0 } else { 3 + k };
        nodes.push(GraphNode {
            name: format!("decoy_{k}"),
            seq: if k == 69 {
                seg3.clone()
            } else {
                b"ATATATATATATATATATATATATATATATAT".to_vec()
            },
            preds: vec![pred],
        });
    }

    let idx = DagSeedIndex::build(&nodes);
    let aln = align_seq_to_dag_seeded_indexed(&nodes, &idx, &query, &HmmParams::default(), 10.0);
    assert!(
        aln.path_nodes.contains(&1) && aln.path_nodes.contains(&2) && aln.path_nodes.contains(&3),
        "expected topological chain [1, 2, 3], got {:?}",
        aln.path_nodes
    );
    assert!(
        !aln.path_nodes.contains(&73),
        "off-chain decoy node 73 should be pruned by anchor chaining"
    );
}

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
    let idx = DagSeedIndex::build(&nodes);
    let aln = align_seq_to_dag_seeded_indexed(&nodes, &idx, query, &HmmParams::default(), 10.0);
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
    let idx = DagSeedIndex::build(&nodes);
    let aln = align_seq_to_dag_seeded_indexed(&nodes, &idx, query, &HmmParams::default(), 10.0);
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
    let idx = DagSeedIndex::build(&nodes);
    let aln = align_seq_to_dag_seeded_indexed(&nodes, &idx, &query, &HmmParams::default(), 10.0);
    assert_eq!(aln.path_nodes, [] as [usize; 0]);
    assert_eq!(aln.kept_pairs, [] as [(usize, usize, usize); 0]);
    assert_eq!(aln.bubbles, 1);
}

#[test]
fn seeded_dag_falls_back_to_aa5_translated_seeds_on_synonymous_divergence() {
    // Synonymous mutations (AGA->CGA, CCA->CCT, GTT->GTA, ACT->ACA) flip purine/pyrimidine
    // every 3-5 bp, breaking every MASK_W10 and MASK_RY18 seed while keeping 100%
    // amino-acid identity (RPVT * 4 = 16 aa) so aa5_keys seeds the target node.
    let target_coding = b"AGACCAGTTACTAGACCAGTTACTAGACCAGTTACTAGACCAGTTACT";
    let query_coding = b"CGACCTGTAACACGACCTGTAACACGACCTGTAACACGACCTGTAACA";
    assert_eq!(target_coding.len(), 48);
    assert_eq!(query_coding.len(), 48);

    let mut nodes: Vec<GraphNode> = (0..11)
        .map(|idx| GraphNode {
            name: format!("decoy_{idx}"),
            seq: vec![b'N'; 5_000],
            preds: if idx == 0 { vec![] } else { vec![idx - 1] },
        })
        .collect();
    let target_id = nodes.len();
    nodes.push(GraphNode {
        name: "coding_target".into(),
        seq: target_coding.to_vec(),
        preds: vec![],
    });

    let index = DagSeedIndex::build(&nodes);
    let q_w10 = FlatSeedIndex::build_with_sdust(query_coding, MASK_W10, 8, false, false, None)
        .expect("build w10 seed index");
    let q_ry18 = FlatSeedIndex::build_with_sdust(query_coding, MASK_RY18, 8, false, true, None)
        .expect("build ry18 seed index");
    assert_eq!(
        q_w10.shared_unique_keys(index.w10[target_id].as_ref().expect("target w10 index")),
        0
    );
    assert_eq!(
        q_ry18.shared_unique_keys(index.ry18[target_id].as_ref().expect("target ry18 index")),
        0
    );

    let aln =
        align_seq_to_dag_seeded_indexed(&nodes, &index, query_coding, &HmmParams::default(), 8.0);
    assert_eq!(aln.path_nodes, vec![target_id]);
    assert!(aln.kept_pairs.len() >= 30);
}

#[test]
fn seeded_dag_aa5_fallback_handles_minus_strand_coding_genes() {
    let target_coding =
        crate::seeding::revcomp(b"AGACCAGTTACTAGACCAGTTACTAGACCAGTTACTAGACCAGTTACT");
    let query_coding = crate::seeding::revcomp(b"CGACCTGTAACACGACCTGTAACACGACCTGTAACACGACCTGTAACA");
    let mut nodes: Vec<GraphNode> = (0..11)
        .map(|idx| GraphNode {
            name: format!("decoy_{idx}"),
            seq: vec![b'N'; 5_000],
            preds: if idx == 0 { vec![] } else { vec![idx - 1] },
        })
        .collect();
    let target_id = nodes.len();
    nodes.push(GraphNode {
        name: "minus_coding_target".into(),
        seq: target_coding,
        preds: vec![],
    });
    let index = DagSeedIndex::build(&nodes);
    let aln =
        align_seq_to_dag_seeded_indexed(&nodes, &index, &query_coding, &HmmParams::default(), 8.0);
    assert_eq!(aln.path_nodes, vec![target_id]);
    assert!(aln.kept_pairs.len() >= 30);
}

#[test]
fn seeded_dag_slices_coarse_rgfa_long_node_and_remaps_offsets() {
    let mut state: u64 = 0xCAFE_BABE_1234_5678;
    let bases = *b"ACGT";
    let big_seq: Vec<u8> = (0..220_000)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            bases[((state >> 32) & 3) as usize]
        })
        .collect();
    let query = big_seq[110_000..111_000].to_vec();
    let nodes = vec![GraphNode {
        name: "s1".to_owned(),
        seq: big_seq,
        preds: Vec::new(),
    }];
    let index = DagSeedIndex::build(&nodes);
    let aln = align_seq_to_dag_seeded_indexed(&nodes, &index, &query, &HmmParams::default(), 10.0);
    assert_eq!(aln.path_nodes, vec![0]);
    assert!(
        aln.kept_pairs.len() >= 950,
        "expected >=950/1000 bp aligned on 220 kb coarse rGFA node, got {}",
        aln.kept_pairs.len()
    );
    let min_node_pos = aln
        .kept_pairs
        .iter()
        .map(|&(_, np, _)| np)
        .min()
        .unwrap_or(0);
    let max_node_pos = aln
        .kept_pairs
        .iter()
        .map(|&(_, np, _)| np)
        .max()
        .unwrap_or(0);
    assert!(
        (109_990..=110_010).contains(&min_node_pos) && (110_990..=111_010).contains(&max_node_pos),
        "remapped node offsets should lie at [110_000..111_000), got [{min_node_pos}..{max_node_pos}]"
    );
}
