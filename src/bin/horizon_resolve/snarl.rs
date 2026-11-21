use anyhow::Result;
use horizon_aligner::gfa::{
    GfaRecord, GfaStreamOptions, parse_path_steps, parse_walk_steps, stream_gfa_records_with_opts,
};
use std::collections::HashMap;
use std::ops::ControlFlow;

/// Compact `u32`-interned biedged pangenome graph in CSR format with 2-edge-connected
/// snarl/bubble boundary annotations.
pub(crate) struct PackedGraph {
    pub(crate) seg_to_id: HashMap<String, u32>,
    pub(crate) id_to_seg: Vec<String>,
    pub(crate) seg_lens: Vec<usize>,
    pub(crate) is_anchor: Vec<bool>,
    pub(crate) num_bridges: usize,
    pub(crate) num_components: usize,
}

impl PackedGraph {
    pub(crate) fn get_id(&self, name: &str) -> Option<u32> {
        self.seg_to_id.get(name).copied()
    }

    pub(crate) fn intern(&mut self, name: &str, len: usize) -> u32 {
        if let Some(&id) = self.seg_to_id.get(name) {
            if len > 0 && self.seg_lens[id as usize] == 0 {
                self.seg_lens[id as usize] = len;
            }
            return id;
        }
        let id = u32::try_from(self.id_to_seg.len()).expect("segment count fits u32");
        let owned = name.to_owned();
        self.seg_to_id.insert(owned.clone(), id);
        self.id_to_seg.push(owned);
        self.seg_lens.push(len);
        id
    }
}

#[inline]
const fn exit_side(seg: u32, fwd: bool) -> u32 {
    2 * seg + if fwd { 1 } else { 0 }
}

#[inline]
const fn entry_side(seg: u32, fwd: bool) -> u32 {
    2 * seg + if fwd { 0 } else { 1 }
}

fn tarjan_bridges(
    num_sides: usize,
    total_edges: usize,
    head: &[usize],
    csr_to: &[u32],
    csr_eid: &[usize],
) -> (Vec<bool>, usize) {
    let mut tin = vec![0u32; num_sides];
    let mut low = vec![0u32; num_sides];
    let mut is_bridge = vec![false; total_edges];
    let mut timer = 0u32;
    let mut num_bridges = 0usize;
    let mut stack: Vec<(usize, usize, usize)> = Vec::new();

    for root in 0..num_sides {
        if tin[root] != 0 {
            continue;
        }
        timer += 1;
        tin[root] = timer;
        low[root] = timer;
        stack.push((root, usize::MAX, head[root]));

        while let Some((u, peid, idx)) = stack.last_mut() {
            let u_node = *u;
            let parent_eid = *peid;
            let end = head[u_node + 1];
            if *idx < end {
                let cur = *idx;
                *idx += 1;
                let v = csr_to[cur] as usize;
                let eid = csr_eid[cur];
                if eid == parent_eid {
                    continue;
                }
                if tin[v] != 0 {
                    low[u_node] = low[u_node].min(tin[v]);
                } else {
                    timer += 1;
                    tin[v] = timer;
                    low[v] = timer;
                    stack.push((v, eid, head[v]));
                }
            } else {
                let finished_u = u_node;
                let finished_peid = parent_eid;
                stack.pop();
                if let Some(&(p, _, _)) = stack.last() {
                    low[p] = low[p].min(low[finished_u]);
                    if low[finished_u] > tin[p] && !is_bridge[finished_peid] {
                        is_bridge[finished_peid] = true;
                        num_bridges += 1;
                    }
                }
            }
        }
    }
    (is_bridge, num_bridges)
}

/// Build CSR adjacency for the biedged graph (black segment edges + gray link edges)
/// and run iterative Tarjan bridge-finding to identify snarl boundary anchors.
fn compute_biedged_snarls(
    num_segs: usize,
    seg_lens: &[usize],
    seg_visits: &[u32],
    gray_edges: &[(u32, u32)],
    min_anchor_bp: usize,
) -> (Vec<bool>, usize, usize) {
    let num_sides = num_segs * 2;
    if num_sides == 0 {
        return (Vec::new(), 0, 0);
    }

    let total_edges = num_segs + gray_edges.len();
    let mut deg = vec![1u32; num_sides];
    for &(u, v) in gray_edges {
        deg[u as usize] += 1;
        deg[v as usize] += 1;
    }

    let mut head = vec![0usize; num_sides + 1];
    for s in 0..num_sides {
        head[s + 1] = head[s] + deg[s] as usize;
    }
    let mut csr_to = vec![0u32; head[num_sides]];
    let mut csr_eid = vec![0usize; head[num_sides]];
    let mut cursor = head[..num_sides].to_vec();

    for v in 0..num_segs {
        let s0 = 2 * v;
        let s1 = s0 + 1;
        let p0 = cursor[s0];
        cursor[s0] += 1;
        csr_to[p0] = u32::try_from(s1).expect("side fits u32");
        csr_eid[p0] = v;

        let p1 = cursor[s1];
        cursor[s1] += 1;
        csr_to[p1] = u32::try_from(s0).expect("side fits u32");
        csr_eid[p1] = v;
    }

    for (idx, &(u, v)) in gray_edges.iter().enumerate() {
        let eid = num_segs + idx;
        let pu = cursor[u as usize];
        cursor[u as usize] += 1;
        csr_to[pu] = v;
        csr_eid[pu] = eid;

        let pv = cursor[v as usize];
        cursor[v as usize] += 1;
        csr_to[pv] = u;
        csr_eid[pv] = eid;
    }

    let (is_bridge, num_bridges) = tarjan_bridges(num_sides, total_edges, &head, &csr_to, &csr_eid);

    // Count 2-edge-connected components by BFS ignoring bridge edges.
    let mut comp = vec![u32::MAX; num_sides];
    let mut num_components = 0usize;
    let mut queue = Vec::new();
    for start in 0..num_sides {
        if comp[start] != u32::MAX {
            continue;
        }
        let cid = u32::try_from(num_components).expect("component fits u32");
        num_components += 1;
        comp[start] = cid;
        queue.push(start);
        let mut q_head = 0usize;
        while q_head < queue.len() {
            let u = queue[q_head];
            q_head += 1;
            for cur in head[u]..head[u + 1] {
                let eid = csr_eid[cur];
                if is_bridge[eid] {
                    continue;
                }
                let v = csr_to[cur] as usize;
                if comp[v] == u32::MAX {
                    comp[v] = cid;
                    queue.push(v);
                }
            }
        }
        queue.clear();
    }

    // A segment is a snarl boundary anchor if:
    // - its black edge is a bridge AND (shared by >=2 paths, or >= min_anchor_bp, or a graph tip), OR
    // - it is shared by >=2 paths AND >= min_anchor_bp.
    let mut is_anchor = vec![false; num_segs];
    for v in 0..num_segs {
        let black_bridge = is_bridge[v];
        let is_tip = deg[2 * v] == 1 || deg[2 * v + 1] == 1;
        let multi_visit = seg_visits[v] >= 2;
        let long_seg = seg_lens[v] >= min_anchor_bp;
        if (black_bridge && (multi_visit || long_seg || is_tip)) || (multi_visit && long_seg) {
            is_anchor[v] = true;
        }
    }

    (is_anchor, num_bridges, num_components)
}

fn record_walk_steps(
    steps: &[(&str, bool)],
    graph: &mut PackedGraph,
    seg_visits: &mut Vec<u32>,
    visit_stamp: &mut Vec<u32>,
    stamp: u32,
    gray_edges: &mut Vec<(u32, u32)>,
) {
    let mut prev: Option<(u32, bool)> = None;
    for &(seg, fwd) in steps {
        let id = graph.intern(seg, 0);
        let idx = id as usize;
        if idx >= seg_visits.len() {
            seg_visits.resize(idx + 1, 0);
            visit_stamp.resize(idx + 1, 0);
        }
        if visit_stamp[idx] != stamp {
            visit_stamp[idx] = stamp;
            seg_visits[idx] += 1;
        }
        if let Some((pu, pf)) = prev {
            let su = exit_side(pu, pf);
            let sv = entry_side(id, fwd);
            gray_edges.push(if su <= sv { (su, sv) } else { (sv, su) });
        }
        prev = Some((id, fwd));
    }
}

pub(crate) fn pass1_index(gfa_path: &str, min_anchor_bp: usize) -> Result<PackedGraph> {
    let mut graph = PackedGraph {
        seg_to_id: HashMap::new(),
        id_to_seg: Vec::new(),
        seg_lens: Vec::new(),
        is_anchor: Vec::new(),
        num_bridges: 0,
        num_components: 0,
    };
    let mut seg_visits: Vec<u32> = Vec::new();
    let mut visit_stamp: Vec<u32> = Vec::new();
    let mut cur_stamp = 0u32;
    let mut gray_edges: Vec<(u32, u32)> = Vec::new();

    stream_gfa_records_with_opts(
        gfa_path,
        GfaStreamOptions {
            skip_segment_seq: true,
        },
        |rec| {
            match rec {
                GfaRecord::Segment(s) => {
                    let len = s.seq_len();
                    let id = graph.intern(&s.name, len);
                    let idx = id as usize;
                    if idx >= seg_visits.len() {
                        seg_visits.resize(idx + 1, 0);
                        visit_stamp.resize(idx + 1, 0);
                    }
                }
                GfaRecord::Link(l) => {
                    let u = graph.intern(&l.from, 0);
                    let v = graph.intern(&l.to, 0);
                    let max_idx = (u.max(v) as usize) + 1;
                    if max_idx > seg_visits.len() {
                        seg_visits.resize(max_idx, 0);
                        visit_stamp.resize(max_idx, 0);
                    }
                    let su = exit_side(u, l.from_orient != b'-');
                    let sv = entry_side(v, l.to_orient != b'-');
                    gray_edges.push(if su <= sv { (su, sv) } else { (sv, su) });
                }
                GfaRecord::Path(p) => {
                    cur_stamp += 1;
                    let steps = parse_path_steps(&p.segment_names);
                    record_walk_steps(
                        &steps,
                        &mut graph,
                        &mut seg_visits,
                        &mut visit_stamp,
                        cur_stamp,
                        &mut gray_edges,
                    );
                }
                GfaRecord::Walk(w) => {
                    cur_stamp += 1;
                    let steps = parse_walk_steps(&w.walk);
                    record_walk_steps(
                        &steps,
                        &mut graph,
                        &mut seg_visits,
                        &mut visit_stamp,
                        cur_stamp,
                        &mut gray_edges,
                    );
                }
            }
            Ok(ControlFlow::Continue(()))
        },
    )?;

    gray_edges.sort_unstable();
    gray_edges.dedup();
    let (is_anchor, num_bridges, num_components) = compute_biedged_snarls(
        graph.id_to_seg.len(),
        &graph.seg_lens,
        &seg_visits,
        &gray_edges,
        min_anchor_bp,
    );
    graph.is_anchor = is_anchor;
    graph.num_bridges = num_bridges;
    graph.num_components = num_components;
    Ok(graph)
}
