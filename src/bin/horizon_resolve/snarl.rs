use anyhow::Result;
use horizon_aligner::gfa::{
    GfaRecordRef, parse_cigar_overlap_bp, parse_path_steps, parse_walk_steps,
    stream_gfa_record_refs,
};
use std::borrow::Cow;
use std::collections::HashMap;
use std::io::Write;
use std::ops::ControlFlow;

/// 4-byte packed oriented segment `(id << 1) | fwd`.
///
/// Ordering matches `(u32, bool)` (`false < true`) so canonical snarl keys stay identical.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub(crate) struct PackedStep(u32);

impl PackedStep {
    #[inline]
    pub(crate) const fn new(id: u32, fwd: bool) -> Self {
        Self((id << 1) | if fwd { 1 } else { 0 })
    }

    #[inline]
    pub(crate) const fn id(self) -> u32 {
        self.0 >> 1
    }

    #[inline]
    pub(crate) const fn fwd(self) -> bool {
        (self.0 & 1) != 0
    }

    #[inline]
    pub(crate) const fn rev(self) -> Self {
        Self(self.0 ^ 1)
    }

    #[inline]
    pub(crate) const fn raw(self) -> u32 {
        self.0
    }
}

/// 8-byte packed directed GFA edge `(u_step << 32) | v_step`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub(crate) struct PackedEdge(u64);

impl PackedEdge {
    #[inline]
    pub(crate) const fn directed(u: PackedStep, v: PackedStep) -> Self {
        Self(((u.raw() as u64) << 32) | (v.raw() as u64))
    }

    #[inline]
    pub(crate) const fn src(self) -> PackedStep {
        #[allow(clippy::cast_possible_truncation)]
        PackedStep((self.0 >> 32) as u32)
    }

    #[inline]
    pub(crate) const fn dst(self) -> PackedStep {
        #[allow(clippy::cast_possible_truncation)]
        PackedStep(self.0 as u32)
    }

    #[inline]
    pub(crate) const fn rev(self) -> Self {
        Self::directed(self.dst().rev(), self.src().rev())
    }

    #[inline]
    pub(crate) const fn canonical(self) -> Self {
        let r = self.rev();
        if self.0 <= r.0 { self } else { r }
    }
}

/// Parse a canonical decimal `u32` segment ID (`"0"`, `"1"`, `"42"`, ...), rejecting
/// leading zeros (`"01"`) and `u32::MAX` (reserved as the non-numeric sentinel).
#[inline]
fn parse_canonical_u32(name: &str) -> Option<u32> {
    let bytes = name.as_bytes();
    if bytes.is_empty() || bytes.len() > 10 || (bytes.len() > 1 && bytes[0] == b'0') {
        return None;
    }
    let mut val = 0u64;
    for &b in bytes {
        if !b.is_ascii_digit() {
            return None;
        }
        val = val * 10 + u64::from(b - b'0');
    }
    if val >= u64::from(u32::MAX) {
        return None;
    }
    #[allow(clippy::cast_possible_truncation)]
    Some(val as u32)
}

/// Append decimal `u32` directly into a `String` without `core::fmt` overhead.
#[inline]
fn push_u32_decimal(mut n: u32, out: &mut String) {
    let mut buf = [0u8; 10];
    let mut i = buf.len();
    loop {
        i -= 1;
        buf[i] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    // SAFETY-free ASCII decimal conversion (forbid(unsafe_code)).
    if let Ok(s) = std::str::from_utf8(&buf[i..]) {
        out.push_str(s);
    }
}

/// Write decimal `u32` directly into a byte `Write` sink without heap allocation.
#[inline]
fn write_u32_decimal<W: Write>(mut n: u32, out: &mut W) -> std::io::Result<()> {
    let mut buf = [0u8; 10];
    let mut i = buf.len();
    loop {
        i -= 1;
        buf[i] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    out.write_all(&buf[i..])
}

/// Compact `u32`-interned biedged pangenome graph in CSR format with 2-edge-connected
/// snarl/bubble boundary annotations.
///
/// Standard 1-based sequential variation GFAs (`1`, `2`, ..., `N`) use zero bytes for
/// segment-name tables; arbitrary numeric or alphanumeric GFA segment names transparently
/// promote to compact `u32` + string-fallback tables.
pub(crate) struct PackedGraph {
    sequential_1based: bool,
    num_segs: usize,
    num_to_id: HashMap<u32, u32>,
    str_to_id: HashMap<String, u32>,
    id_num: Vec<u32>,
    non_num_names: HashMap<u32, String>,
    pub(crate) seg_lens: Vec<u32>,
    pub(crate) is_anchor: Vec<bool>,
    pub(crate) is_bridge_seg: Vec<bool>,
    pub(crate) link_overlaps: HashMap<PackedEdge, u32>,
    pub(crate) num_bridges: usize,
    pub(crate) num_components: usize,
}

impl PackedGraph {
    #[inline]
    pub(crate) fn num_segs(&self) -> usize {
        self.num_segs
    }

    #[inline]
    pub(crate) fn seg_len(&self, id: u32) -> usize {
        self.seg_lens.get(id as usize).copied().unwrap_or(0) as usize
    }

    #[inline]
    pub(crate) fn get_id(&self, name: &str) -> Option<u32> {
        if self.sequential_1based {
            let num = parse_canonical_u32(name)?;
            let id = num.checked_sub(1)?;
            ((id as usize) < self.num_segs).then_some(id)
        } else if let Some(num) = parse_canonical_u32(name) {
            self.num_to_id.get(&num).copied()
        } else {
            self.str_to_id.get(name).copied()
        }
    }

    pub(crate) fn push_name(&self, id: u32, out: &mut String) {
        if self.sequential_1based {
            push_u32_decimal(id + 1, out);
            return;
        }
        let num = self.id_num[id as usize];
        if num != u32::MAX {
            push_u32_decimal(num, out);
        } else if let Some(s) = self.non_num_names.get(&id) {
            out.push_str(s);
        }
    }

    pub(crate) fn write_name<W: Write>(&self, id: u32, out: &mut W) -> std::io::Result<()> {
        if self.sequential_1based {
            return write_u32_decimal(id + 1, out);
        }
        let num = self.id_num[id as usize];
        if num != u32::MAX {
            write_u32_decimal(num, out)
        } else if let Some(s) = self.non_num_names.get(&id) {
            out.write_all(s.as_bytes())
        } else {
            Ok(())
        }
    }

    pub(crate) fn seg_name_cow(&self, id: u32) -> Cow<'_, str> {
        if self.sequential_1based {
            let mut s = String::with_capacity(10);
            push_u32_decimal(id + 1, &mut s);
            return Cow::Owned(s);
        }
        let num = self.id_num[id as usize];
        if num == u32::MAX {
            Cow::Borrowed(self.non_num_names.get(&id).map_or("", String::as_str))
        } else {
            let mut s = String::with_capacity(10);
            push_u32_decimal(num, &mut s);
            Cow::Owned(s)
        }
    }

    fn demote_sequential(&mut self) {
        self.sequential_1based = false;
        self.num_to_id.reserve(self.num_segs + 1);
        self.id_num.reserve(self.num_segs + 1);
        for prev in 0..self.num_segs {
            #[allow(clippy::cast_possible_truncation)]
            let id = prev as u32;
            let num = id + 1;
            self.num_to_id.insert(num, id);
            self.id_num.push(num);
        }
    }

    pub(crate) fn intern(&mut self, name: &str, len: usize) -> Result<u32> {
        let len_u32 = u32::try_from(len).unwrap_or(u32::MAX);
        if let Some(id) = self.get_id(name) {
            if len_u32 > 0 && self.seg_lens[id as usize] == 0 {
                self.seg_lens[id as usize] = len_u32;
            }
            return Ok(id);
        }
        // `exit_side`/`entry_side` and `PackedStep` compute `2 * seg + 1` in `u32`, so cap at `u32::MAX / 2`.
        anyhow::ensure!(
            self.num_segs < (u32::MAX / 2) as usize,
            "GFA segment count exceeds maximum representable biedged side ID (u32::MAX / 2)"
        );
        #[allow(clippy::cast_possible_truncation)]
        let id = self.num_segs as u32;
        if self.sequential_1based {
            if parse_canonical_u32(name) == Some(id + 1) {
                self.num_segs += 1;
                self.seg_lens.push(len_u32);
                return Ok(id);
            }
            self.demote_sequential();
        }
        if let Some(num) = parse_canonical_u32(name) {
            self.num_to_id.insert(num, id);
            self.id_num.push(num);
        } else {
            let owned = name.to_owned();
            self.str_to_id.insert(owned.clone(), id);
            self.id_num.push(u32::MAX);
            self.non_num_names.insert(id, owned);
        }
        self.num_segs += 1;
        self.seg_lens.push(len_u32);
        Ok(id)
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
    csr_eid: &[u32],
) -> (Vec<bool>, usize) {
    let mut tin = vec![0u32; num_sides];
    let mut low = vec![0u32; num_sides];
    let mut is_bridge = vec![false; total_edges];
    let mut timer = 0u32;
    let mut num_bridges = 0usize;
    let mut stack: Vec<(u32, u32, usize)> = Vec::new();

    for root in 0..num_sides {
        if tin[root] != 0 {
            continue;
        }
        timer += 1;
        tin[root] = timer;
        low[root] = timer;
        #[allow(clippy::cast_possible_truncation)]
        stack.push((root as u32, u32::MAX, head[root]));

        while let Some((u, peid, idx)) = stack.last_mut() {
            let u_node = *u as usize;
            let parent_eid = *peid;
            let end = head[u_node + 1];
            if *idx < end {
                let cur = *idx;
                *idx += 1;
                let v_u32 = csr_to[cur];
                let v = v_u32 as usize;
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
                    stack.push((v_u32, eid, head[v]));
                }
            } else {
                let finished_u = u_node;
                let finished_peid = parent_eid;
                stack.pop();
                if let Some(&(p_u32, _, _)) = stack.last() {
                    let p = p_u32 as usize;
                    low[p] = low[p].min(low[finished_u]);
                    let pe = finished_peid as usize;
                    if low[finished_u] > tin[p] && !is_bridge[pe] {
                        is_bridge[pe] = true;
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
    seg_lens: &[u32],
    seg_visits: &[u32],
    gray_edges: &[(u32, u32)],
    min_anchor_bp: usize,
) -> Result<(Vec<bool>, Vec<bool>, usize, usize)> {
    let num_sides = num_segs * 2;
    if num_sides == 0 {
        return Ok((Vec::new(), Vec::new(), 0, 0));
    }

    let total_edges = num_segs + gray_edges.len();
    anyhow::ensure!(
        total_edges < u32::MAX as usize,
        "GFA biedged edge count ({total_edges}) exceeds u32::MAX - 1"
    );
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
    let mut csr_eid = vec![0u32; head[num_sides]];
    let mut cursor = head[..num_sides].to_vec();

    for v in 0..num_segs {
        let s0 = 2 * v;
        let s1 = s0 + 1;
        #[allow(clippy::cast_possible_truncation)]
        let eid = v as u32;
        let p0 = cursor[s0];
        cursor[s0] += 1;
        csr_to[p0] = u32::try_from(s1)
            .map_err(|_| anyhow::anyhow!("GFA graph side count exceeds u32::MAX"))?;
        csr_eid[p0] = eid;

        let p1 = cursor[s1];
        cursor[s1] += 1;
        csr_to[p1] = u32::try_from(s0)
            .map_err(|_| anyhow::anyhow!("GFA graph side count exceeds u32::MAX"))?;
        csr_eid[p1] = eid;
    }

    for (idx, &(u, v)) in gray_edges.iter().enumerate() {
        #[allow(clippy::cast_possible_truncation)]
        let eid = (num_segs + idx) as u32;
        let pu = cursor[u as usize];
        cursor[u as usize] += 1;
        csr_to[pu] = v;
        csr_eid[pu] = eid;

        let pv = cursor[v as usize];
        cursor[v as usize] += 1;
        csr_to[pv] = u;
        csr_eid[pv] = eid;
    }
    drop(cursor);

    let (is_bridge, num_bridges) = tarjan_bridges(num_sides, total_edges, &head, &csr_to, &csr_eid);

    // Count 2-edge-connected components by BFS ignoring bridge edges.
    let mut comp = vec![u32::MAX; num_sides];
    let mut num_components = 0usize;
    let mut queue: Vec<u32> = Vec::new();
    for start in 0..num_sides {
        if comp[start] != u32::MAX {
            continue;
        }
        let cid = u32::try_from(num_components)
            .map_err(|_| anyhow::anyhow!("GFA component count exceeds u32::MAX"))?;
        num_components += 1;
        comp[start] = cid;
        #[allow(clippy::cast_possible_truncation)]
        queue.push(start as u32);
        let mut q_head = 0usize;
        while q_head < queue.len() {
            let u = queue[q_head] as usize;
            q_head += 1;
            for cur in head[u]..head[u + 1] {
                let eid = csr_eid[cur] as usize;
                if is_bridge[eid] {
                    continue;
                }
                let v_u32 = csr_to[cur];
                let v = v_u32 as usize;
                if comp[v] == u32::MAX {
                    comp[v] = cid;
                    queue.push(v_u32);
                }
            }
        }
        queue.clear();
    }

    // A segment is a snarl boundary anchor if:
    // - its black edge is a bridge AND (shared by >=2 paths, or >= min_anchor_bp, or a graph tip), OR
    // - it is shared by >=2 paths AND >= min_anchor_bp.
    let mut is_anchor = vec![false; num_segs];
    let mut is_bridge_seg = vec![false; num_segs];
    for v in 0..num_segs {
        let black_bridge = is_bridge[v];
        is_bridge_seg[v] = black_bridge;
        let is_tip = deg[2 * v] == 1 || deg[2 * v + 1] == 1;
        let multi_visit = seg_visits[v] >= 2;
        let long_seg = (seg_lens[v] as usize) >= min_anchor_bp;
        if (black_bridge && (multi_visit || long_seg || is_tip)) || (multi_visit && long_seg) {
            is_anchor[v] = true;
        }
    }

    Ok((is_anchor, is_bridge_seg, num_bridges, num_components))
}

fn maybe_compact_gray_edges(gray_edges: &mut Vec<(u32, u32)>, dedup_threshold: &mut usize) {
    if gray_edges.len() >= *dedup_threshold {
        gray_edges.sort_unstable();
        gray_edges.dedup();
        *dedup_threshold = (gray_edges.len().saturating_mul(2)).max(1 << 20);
    }
}

fn record_walk_steps(
    steps: &[(&str, bool)],
    graph: &mut PackedGraph,
    seg_visits: &mut Vec<u32>,
    visit_stamp: &mut Vec<u32>,
    stamp: u32,
    gray_edges: &mut Vec<(u32, u32)>,
    dedup_threshold: &mut usize,
) -> Result<()> {
    let mut prev: Option<(u32, bool)> = None;
    for &(seg, fwd) in steps {
        let id = graph.intern(seg, 0)?;
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
    maybe_compact_gray_edges(gray_edges, dedup_threshold);
    Ok(())
}

pub(crate) fn pass1_index(
    gfa_path: impl AsRef<std::path::Path>,
    min_anchor_bp: usize,
) -> Result<PackedGraph> {
    let gfa_path = gfa_path.as_ref();
    let mut graph = PackedGraph {
        sequential_1based: true,
        num_segs: 0,
        num_to_id: HashMap::new(),
        str_to_id: HashMap::new(),
        id_num: Vec::new(),
        non_num_names: HashMap::new(),
        seg_lens: Vec::new(),
        is_anchor: Vec::new(),
        is_bridge_seg: Vec::new(),
        link_overlaps: HashMap::new(),
        num_bridges: 0,
        num_components: 0,
    };
    let mut seg_visits: Vec<u32> = Vec::new();
    let mut visit_stamp: Vec<u32> = Vec::new();
    let mut cur_stamp = 0u32;
    let mut gray_edges: Vec<(u32, u32)> = Vec::new();
    let mut dedup_threshold = 1usize << 20;

    stream_gfa_record_refs(gfa_path, |rec| {
        match rec {
            GfaRecordRef::Segment(s) => {
                let len = s.seq_len();
                let id = graph.intern(s.name, len)?;
                let idx = id as usize;
                if idx >= seg_visits.len() {
                    seg_visits.resize(idx + 1, 0);
                    visit_stamp.resize(idx + 1, 0);
                }
            }
            GfaRecordRef::Link(l) => {
                let u = graph.intern(l.from, 0)?;
                let v = graph.intern(l.to, 0)?;
                let max_idx = (u.max(v) as usize) + 1;
                if max_idx > seg_visits.len() {
                    seg_visits.resize(max_idx, 0);
                    visit_stamp.resize(max_idx, 0);
                }
                let uf = l.from_orient != b'-';
                let vf = l.to_orient != b'-';
                let ov = parse_cigar_overlap_bp(l.overlap);
                if ov > 0 {
                    let ov_u32 = u32::try_from(ov).unwrap_or(u32::MAX);
                    let edge = PackedEdge::directed(PackedStep::new(u, uf), PackedStep::new(v, vf));
                    graph.link_overlaps.insert(edge, ov_u32);
                    graph.link_overlaps.insert(edge.rev(), ov_u32);
                }
                let su = exit_side(u, uf);
                let sv = entry_side(v, vf);
                gray_edges.push(if su <= sv { (su, sv) } else { (sv, su) });
                maybe_compact_gray_edges(&mut gray_edges, &mut dedup_threshold);
            }
            GfaRecordRef::Path(p) => {
                cur_stamp += 1;
                let steps = parse_path_steps(p.segment_names);
                record_walk_steps(
                    &steps,
                    &mut graph,
                    &mut seg_visits,
                    &mut visit_stamp,
                    cur_stamp,
                    &mut gray_edges,
                    &mut dedup_threshold,
                )?;
            }
            GfaRecordRef::Walk(w) => {
                cur_stamp += 1;
                let steps = parse_walk_steps(w.walk);
                record_walk_steps(
                    &steps,
                    &mut graph,
                    &mut seg_visits,
                    &mut visit_stamp,
                    cur_stamp,
                    &mut gray_edges,
                    &mut dedup_threshold,
                )?;
            }
        }
        Ok(ControlFlow::Continue(()))
    })?;

    gray_edges.sort_unstable();
    gray_edges.dedup();
    drop(visit_stamp);
    let (is_anchor, is_bridge_seg, num_bridges, num_components) = compute_biedged_snarls(
        graph.num_segs(),
        &graph.seg_lens,
        &seg_visits,
        &gray_edges,
        min_anchor_bp,
    )?;
    graph.is_anchor = is_anchor;
    graph.is_bridge_seg = is_bridge_seg;
    graph.num_bridges = num_bridges;
    graph.num_components = num_components;
    Ok(graph)
}
