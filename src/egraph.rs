//! Congruence-closed e-graph for resolving tandem duplications and inversions across LCBs.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

/// Immutable e-node referencing canonical child e-class IDs.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ENode {
    /// Operator name (`blk:<label>`, `seq`, `inv`, or `dup`).
    pub op: String,
    /// Child e-class IDs.
    pub children: Vec<usize>,
}

/// Extracted minimum-cost synteny expression tree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Expr {
    /// Leaf synteny block label.
    Block(String),
    /// Sequential concatenation.
    Seq(Box<Expr>, Box<Expr>),
    /// Inverted block or span.
    Inv(Box<Expr>),
    /// Collapsed 2-copy tandem duplication.
    Dup(Box<Expr>),
    /// Collapsed `k`-copy tandem repeat (`k >= 3`).
    Repeat(Box<Expr>, usize),
}

const MAX_REPEAT_ARITY: usize = 64;

fn parse_rep_count(op: &str) -> Option<usize> {
    let k: usize = op.strip_prefix("rep:")?.parse().ok()?;
    (k >= 3).then_some(k)
}

fn repeat_arity_in(eg: &mut EGraph, node: &ENode) -> Option<(usize, usize)> {
    if node.children.len() != 1 {
        return None;
    }
    if node.op == "dup" {
        return Some((eg.find(node.children[0]), 2));
    }
    parse_rep_count(&node.op).map(|k| (eg.find(node.children[0]), k))
}

fn add_repeat_node(eg: &mut EGraph, base: usize, k: usize) -> usize {
    let op = if k == 2 {
        "dup".to_owned()
    } else {
        format!("rep:{k}")
    };
    eg.add(ENode {
        op,
        children: vec![base],
    })
}

/// Minimal congruence-closed e-graph (Willsey et al. 2021).
#[derive(Debug, Default)]
pub struct EGraph {
    parent: Vec<usize>,
    classes: BTreeMap<usize, BTreeSet<ENode>>,
    parents: Vec<Vec<usize>>,
    memo: BTreeMap<ENode, usize>,
    worklist: Vec<usize>,
}

impl EGraph {
    /// Create an empty e-graph.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Find canonical e-class representative with path compression.
    ///
    /// # Panics
    /// Panics if `eid >= self.parent.len()`.
    pub fn find(&mut self, mut eid: usize) -> usize {
        assert!(
            eid < self.parent.len(),
            "e-class id {eid} out of bounds (len={})",
            self.parent.len()
        );
        while self.parent[eid] != eid {
            let p = self.parent[eid];
            self.parent[eid] = self.parent[p];
            eid = self.parent[eid];
        }
        eid
    }

    fn find_imm(&self, mut eid: usize) -> usize {
        debug_assert!(
            eid < self.parent.len(),
            "e-class id {eid} out of bounds (len={})",
            self.parent.len()
        );
        while let Some(&p) = self.parent.get(eid) {
            if p == eid {
                break;
            }
            eid = p;
        }
        eid
    }

    fn canonicalize(&mut self, node: &ENode) -> ENode {
        let children = node.children.iter().map(|&c| self.find(c)).collect();
        ENode {
            op: node.op.clone(),
            children,
        }
    }

    /// Insert an e-node and return its canonical e-class ID.
    #[allow(clippy::needless_pass_by_value)]
    pub fn add(&mut self, node: ENode) -> usize {
        let cnode = self.canonicalize(&node);
        if let Some(&id) = self.memo.get(&cnode) {
            return self.find(id);
        }
        let eid = self.parent.len();
        self.parent.push(eid);
        self.parents.push(Vec::new());
        let mut seen_children = Vec::with_capacity(cnode.children.len());
        for &child in &cnode.children {
            if !seen_children.contains(&child) {
                seen_children.push(child);
                self.parents[child].push(eid);
            }
        }
        let mut set = BTreeSet::new();
        set.insert(cnode.clone());
        self.classes.insert(eid, set);
        self.memo.insert(cnode, eid);
        eid
    }

    /// Merge two e-classes and queue the winner for congruence repair.
    pub fn union(&mut self, a: usize, b: usize) -> usize {
        let mut ra = self.find(a);
        let mut rb = self.find(b);
        if ra == rb {
            return ra;
        }
        let len_a = self.classes.get(&ra).map_or(0, BTreeSet::len);
        let len_b = self.classes.get(&rb).map_or(0, BTreeSet::len);
        if len_a < len_b {
            std::mem::swap(&mut ra, &mut rb);
        }
        self.parent[rb] = ra;
        if let Some(nodes_b) = self.classes.remove(&rb) {
            self.classes.entry(ra).or_default().extend(nodes_b);
        }
        let parents_b = std::mem::take(&mut self.parents[rb]);
        self.parents[ra].extend(parents_b);
        self.worklist.push(ra);
        ra
    }

    /// Restore congruence invariants after pending unions.
    pub fn rebuild(&mut self) {
        while !self.worklist.is_empty() {
            let drained = std::mem::take(&mut self.worklist);
            let todo: BTreeSet<usize> = drained.into_iter().map(|id| self.find(id)).collect();
            for eid in todo {
                self.repair(eid);
            }
        }
    }

    fn repair(&mut self, eid: usize) {
        let eid = self.find(eid);
        let parent_entries = std::mem::take(&mut self.parents[eid]);
        let unique_parents: BTreeSet<usize> =
            parent_entries.into_iter().map(|p| self.find(p)).collect();
        let mut to_union: Vec<(usize, usize)> = Vec::new();
        for peid_raw in &unique_parents {
            let peid = self.find(*peid_raw);
            let Some(old_nodes) = self.classes.remove(&peid) else {
                continue;
            };
            let mut new_nodes = BTreeSet::new();
            for node in old_nodes {
                let cnode = self.canonicalize(&node);
                if cnode != node {
                    self.memo.remove(&node);
                }
                if let Some(&other_raw) = self.memo.get(&cnode) {
                    let other = self.find(other_raw);
                    if other != peid {
                        to_union.push((peid, other));
                    }
                }
                self.memo.insert(cnode.clone(), peid);
                new_nodes.insert(cnode);
            }
            self.classes.insert(peid, new_nodes);
        }
        for (a, b) in to_union {
            let merged = self.union(a, b);
            if let Some(nodes) = self.classes.get(&merged) {
                for node in nodes {
                    self.memo.insert(node.clone(), merged);
                }
            }
        }
        let final_eid = self.find(eid);
        let kept_parents: BTreeSet<usize> =
            unique_parents.into_iter().map(|p| self.find(p)).collect();
        self.parents[final_eid].extend(kept_parents);
    }
}

/// Build a right-associated e-graph expression from signed block labels `(label, strand)`.
pub fn build_seq_expr(eg: &mut EGraph, blocks: &[(&str, i8)]) -> Option<usize> {
    if blocks.is_empty() {
        return None;
    }
    let mut terms = Vec::with_capacity(blocks.len());
    for &(label, strand) in blocks {
        let atom = eg.add(ENode {
            op: format!("blk:{label}"),
            children: Vec::new(),
        });
        if strand >= 0 {
            terms.push(atom);
        } else {
            terms.push(eg.add(ENode {
                op: "inv".to_owned(),
                children: vec![atom],
            }));
        }
    }
    let mut cur = terms[terms.len() - 1];
    for &term in terms[..terms.len() - 1].iter().rev() {
        cur = eg.add(ENode {
            op: "seq".to_owned(),
            children: vec![term, cur],
        });
    }
    Some(cur)
}

fn sorted_class_nodes(eg: &EGraph, id: usize) -> Vec<ENode> {
    let mut nodes: Vec<ENode> = eg
        .classes
        .get(&id)
        .map_or_default(|s| s.iter().cloned().collect());
    nodes.sort_unstable_by(|a, b| a.op.cmp(&b.op).then_with(|| a.children.cmp(&b.children)));
    nodes
}

/// Run bounded equality saturation for inversion cancellation, reversal, tandem collapse, and associativity.
#[allow(clippy::too_many_lines)]
pub fn saturate_synteny(eg: &mut EGraph, max_iters: usize, max_nodes: usize) {
    for _ in 0..max_iters {
        if eg.memo.len() >= max_nodes {
            break;
        }
        let mut merges: Vec<(usize, usize)> = Vec::new();
        let mut class_ids: Vec<usize> = eg.classes.keys().copied().collect();
        class_ids.sort_unstable();
        let snapshot_map: HashMap<usize, Vec<ENode>> = class_ids
            .iter()
            .map(|&id| (id, sorted_class_nodes(eg, id)))
            .collect();

        for &eid in &class_ids {
            let Some(nodes) = snapshot_map.get(&eid) else {
                continue;
            };
            if eg.memo.len() >= max_nodes {
                break;
            }
            for node in nodes {
                if eg.memo.len() >= max_nodes {
                    break;
                }
                if node.op == "inv" && node.children.len() == 1 {
                    let child = eg.find(node.children[0]);
                    let empty = Vec::new();
                    for cnode in snapshot_map.get(&child).unwrap_or(&empty) {
                        if eg.memo.len() >= max_nodes {
                            break;
                        }
                        if cnode.op == "inv" && cnode.children.len() == 1 {
                            merges.push((eid, cnode.children[0]));
                        } else if let Some((base, arity)) = repeat_arity_in(eg, cnode) {
                            let inv_base = eg.add(ENode {
                                op: "inv".to_owned(),
                                children: vec![base],
                            });
                            let rhs = add_repeat_node(eg, inv_base, arity);
                            merges.push((eid, rhs));
                        } else if cnode.op == "seq" && cnode.children.len() == 2 {
                            let (a, b) = (cnode.children[0], cnode.children[1]);
                            let inv_b = eg.add(ENode {
                                op: "inv".to_owned(),
                                children: vec![b],
                            });
                            let inv_a = eg.add(ENode {
                                op: "inv".to_owned(),
                                children: vec![a],
                            });
                            let rhs = eg.add(ENode {
                                op: "seq".to_owned(),
                                children: vec![inv_b, inv_a],
                            });
                            merges.push((eid, rhs));
                        }
                    }
                } else if node.op == "seq" && node.children.len() == 2 {
                    let a = eg.find(node.children[0]);
                    let b = eg.find(node.children[1]);
                    if a == b {
                        let dup_x = eg.add(ENode {
                            op: "dup".to_owned(),
                            children: vec![a],
                        });
                        merges.push((eid, dup_x));
                    }
                    let empty = Vec::new();
                    let a_nodes = snapshot_map.get(&a).unwrap_or(&empty);
                    let b_nodes = snapshot_map.get(&b).unwrap_or(&empty);
                    for anode in a_nodes {
                        if eg.memo.len() >= max_nodes {
                            break;
                        }
                        if let Some((base_a, k_a)) = repeat_arity_in(eg, anode)
                            && base_a == b
                            && k_a < MAX_REPEAT_ARITY
                        {
                            let rhs = add_repeat_node(eg, base_a, k_a + 1);
                            merges.push((eid, rhs));
                        }
                        if anode.op == "seq" && anode.children.len() == 2 {
                            let (x, y) = (anode.children[0], anode.children[1]);
                            let yz = eg.add(ENode {
                                op: "seq".to_owned(),
                                children: vec![y, b],
                            });
                            let rhs = eg.add(ENode {
                                op: "seq".to_owned(),
                                children: vec![x, yz],
                            });
                            merges.push((eid, rhs));
                        }
                    }
                    for bnode in b_nodes {
                        if eg.memo.len() >= max_nodes {
                            break;
                        }
                        if let Some((base_b, k_b)) = repeat_arity_in(eg, bnode)
                            && base_b == a
                            && k_b < MAX_REPEAT_ARITY
                        {
                            let rhs = add_repeat_node(eg, base_b, k_b + 1);
                            merges.push((eid, rhs));
                        }
                        if bnode.op == "seq" && bnode.children.len() == 2 {
                            let (y, z) = (bnode.children[0], bnode.children[1]);
                            let xy = eg.add(ENode {
                                op: "seq".to_owned(),
                                children: vec![a, y],
                            });
                            let rhs = eg.add(ENode {
                                op: "seq".to_owned(),
                                children: vec![xy, z],
                            });
                            merges.push((eid, rhs));
                        }
                    }
                    for anode in a_nodes {
                        if eg.memo.len() >= max_nodes {
                            break;
                        }
                        let Some((base_a, k_a)) = repeat_arity_in(eg, anode) else {
                            continue;
                        };
                        for bnode in b_nodes {
                            if eg.memo.len() >= max_nodes {
                                break;
                            }
                            if let Some((base_b, k_b)) = repeat_arity_in(eg, bnode)
                                && base_a == base_b
                                && k_a + k_b <= MAX_REPEAT_ARITY
                            {
                                let rhs = add_repeat_node(eg, base_a, k_a + k_b);
                                merges.push((eid, rhs));
                            }
                        }
                    }
                }
            }
        }

        if merges.is_empty() {
            break;
        }
        let mut changed = false;
        for (a, b) in merges {
            if eg.find(a) != eg.find(b) {
                eg.union(a, b);
                changed = true;
            }
        }
        eg.rebuild();
        if !changed {
            break;
        }
    }
}

/// MDL bit-cost parameters for e-graph extraction.
#[derive(Debug, Clone)]
pub struct MdlCostModel {
    /// Number of distinct block classes (for `log2(K)` catalog pointer bits).
    pub num_blocks: usize,
    /// Optional per-block `(length_bp, divergence)` for sequence-aware MDL bit costs.
    pub block_stats: HashMap<String, (usize, f64)>,
    /// Bit cost for a `seq` concatenation node.
    pub seq_op_bits: f64,
    /// Bit cost for an `inv` strand-flip node.
    pub inv_op_bits: f64,
    /// Bit overhead for a `dup` copy header (op tag + Elias-gamma copy count).
    pub dup_header_bits: f64,
}

impl Default for MdlCostModel {
    fn default() -> Self {
        Self {
            num_blocks: 4,
            block_stats: HashMap::new(),
            seq_op_bits: 1.0,
            inv_op_bits: 1.0,
            dup_header_bits: 2.5,
        }
    }
}

impl MdlCostModel {
    /// Construct an MDL cost model calibrated to the distinct `blk:` leaves in `eg`.
    #[must_use]
    pub fn from_egraph(eg: &EGraph) -> Self {
        let mut labels = HashSet::new();
        for nodes in eg.classes.values() {
            for n in nodes {
                if let Some(lbl) = n.op.strip_prefix("blk:") {
                    labels.insert(lbl);
                }
            }
        }
        Self {
            num_blocks: labels.len().max(2),
            ..Self::default()
        }
    }

    /// Bit overhead for a `k`-copy tandem repeat (`k >= 2`): `(k - 1) * dup_header_bits`
    /// (one residual-stream framing header per extra copy) plus the Elias-gamma copy-count delta
    /// `2 * floor(log2(k)) + 1 - 3` (so `k = 2` equals `dup_header_bits` identically, and every
    /// `k >= 2` rejects JC69-saturated `d = 0.75` noise at `K = 2`).
    #[must_use]
    pub fn repeat_header_bits(&self, k: usize) -> f64 {
        let k_clamped = k.max(2);
        #[allow(clippy::cast_precision_loss)]
        let extra = (k_clamped - 1) as f64;
        #[allow(clippy::cast_precision_loss)]
        let gamma_k = 2.0 * (k_clamped as f64).log2().floor() + 1.0;
        extra * self.dup_header_bits + (gamma_k - 3.0)
    }

    /// Description length (bits) of a leaf block `label`.
    #[must_use]
    pub fn block_bits(&self, label: &str) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let catalog_bits = (self.num_blocks.max(2) as f64).log2().max(1.0);
        if let Some(&(len_bp, _div)) = self.block_stats.get(label) {
            #[allow(clippy::cast_precision_loss)]
            let seq_bits = 2.0 * (len_bp as f64);
            catalog_bits + seq_bits
        } else {
            // Structural minimum: catalog pointer + 1.0 bit length-tag floor
            catalog_bits + 1.0
        }
    }

    fn block_residual_bits(&self, label: &str) -> f64 {
        if self.block_stats.is_empty() {
            return 0.0;
        }
        if let Some(&(len_bp, div)) = self.block_stats.get(label) {
            let p = if div.is_finite() {
                div.clamp(1e-6, 0.75)
            } else {
                0.75
            };
            let h2 = -p * p.log2() - (1.0 - p) * (1.0 - p).log2() + p * 3.0_f64.log2();
            #[allow(clippy::cast_precision_loss)]
            {
                (len_bp as f64) * h2
            }
        } else {
            0.0
        }
    }
}

/// Extract minimum-MDL-bit-cost synteny expression rooted at `root`.
/// Returns `None` if `root` is out of bounds or has no extractable expression.
#[must_use]
pub fn extract_min_cost(eg: &EGraph, root: usize) -> Option<(f64, Expr)> {
    let model = MdlCostModel::from_egraph(eg);
    extract_min_mdl_cost(eg, root, &model)
}

fn build_expr_from_best(
    eg: &EGraph,
    eid: usize,
    best: &HashMap<usize, (f64, f64, usize, ENode)>,
    visiting: &mut HashSet<usize>,
) -> Option<Expr> {
    if !visiting.insert(eid) {
        return None;
    }
    let (_, _, _, node) = best.get(&eid)?;
    let expr = if let Some(lbl) = node.op.strip_prefix("blk:") {
        Some(Expr::Block(lbl.to_owned()))
    } else if node.op == "inv" && node.children.len() == 1 {
        let cid = eg.find_imm(node.children[0]);
        build_expr_from_best(eg, cid, best, visiting).map(|e| Expr::Inv(Box::new(e)))
    } else if node.op == "dup" && node.children.len() == 1 {
        let cid = eg.find_imm(node.children[0]);
        build_expr_from_best(eg, cid, best, visiting).map(|e| Expr::Dup(Box::new(e)))
    } else if let Some(k) = parse_rep_count(&node.op)
        && node.children.len() == 1
    {
        let cid = eg.find_imm(node.children[0]);
        build_expr_from_best(eg, cid, best, visiting).map(|e| Expr::Repeat(Box::new(e), k))
    } else if node.op == "seq" && node.children.len() == 2 {
        let l_id = eg.find_imm(node.children[0]);
        let r_id = eg.find_imm(node.children[1]);
        let left = build_expr_from_best(eg, l_id, best, visiting)?;
        let right = build_expr_from_best(eg, r_id, best, visiting)?;
        Some(Expr::Seq(Box::new(left), Box::new(right)))
    } else {
        None
    };
    visiting.remove(&eid);
    expr
}

/// Extract minimum-parsimony synteny expression rooted at `root` under an explicit [`MdlCostModel`].
/// Returns `None` if `root` is out of bounds or has no extractable expression.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn extract_min_mdl_cost(eg: &EGraph, root: usize, model: &MdlCostModel) -> Option<(f64, Expr)> {
    if root >= eg.parent.len() {
        return None;
    }
    // Map canonical e-class -> (total_bits, residual_bits_per_copy, leaf_block_count, best_node).
    let mut best: HashMap<usize, (f64, f64, usize, ENode)> = HashMap::new();
    let max_rounds = eg.classes.len() + 2;
    let sorted_classes: Vec<(usize, Vec<&ENode>)> = eg
        .classes
        .keys()
        .copied()
        .filter_map(|raw_id| {
            let eid = eg.find_imm(raw_id);
            let node_set = eg.classes.get(&eid)?;
            let nodes: Vec<&ENode> = node_set.iter().collect();
            Some((eid, nodes))
        })
        .collect();

    for _ in 0..max_rounds {
        let mut updated = false;
        for (eid, nodes) in &sorted_classes {
            let eid = *eid;
            for &node in nodes {
                if let Some(lbl) = node.op.strip_prefix("blk:") {
                    let cost = model.block_bits(lbl);
                    let res = model.block_residual_bits(lbl);
                    if best.get(&eid).is_none_or(|(cur, _, _, _)| cost < *cur) {
                        best.insert(eid, (cost, res, 1, (*node).clone()));
                        updated = true;
                    }
                } else if node.op == "inv" && node.children.len() == 1 {
                    let cid = eg.find_imm(node.children[0]);
                    if let Some(&(c, res, leaves, _)) = best.get(&cid) {
                        let cost = c + model.inv_op_bits;
                        if best.get(&eid).is_none_or(|(cur, _, _, _)| cost < *cur) {
                            best.insert(eid, (cost, res, leaves, (*node).clone()));
                            updated = true;
                        }
                    }
                } else if node.op == "dup" && node.children.len() == 1 {
                    let cid = eg.find_imm(node.children[0]);
                    if let Some(&(c, res, leaves, _)) = best.get(&cid) {
                        // Charge dup framing overhead per duplicated leaf block so Dup(Seq(B,B))
                        // cannot bypass the d = 0.75 JC69 saturation rejection gate.
                        #[allow(clippy::cast_precision_loss)]
                        let leaf_scale = leaves.max(1) as f64;
                        let cost = c + leaf_scale * model.dup_header_bits + res;
                        if best.get(&eid).is_none_or(|(cur, _, _, _)| cost < *cur) {
                            best.insert(
                                eid,
                                (cost, 2.0 * res, leaves.saturating_mul(2), (*node).clone()),
                            );
                            updated = true;
                        }
                    }
                } else if let Some(k) = parse_rep_count(&node.op)
                    && node.children.len() == 1
                {
                    let cid = eg.find_imm(node.children[0]);
                    if let Some(&(c, res, leaves, _)) = best.get(&cid) {
                        #[allow(clippy::cast_precision_loss)]
                        let extra_copies = (k - 1) as f64;
                        #[allow(clippy::cast_precision_loss)]
                        let total_copies = k as f64;
                        #[allow(clippy::cast_precision_loss)]
                        let leaf_scale = leaves.max(1) as f64;
                        #[allow(clippy::cast_precision_loss)]
                        let gamma_k = 2.0 * (k.max(2) as f64).log2().floor() + 1.0;
                        let header =
                            extra_copies * leaf_scale * model.dup_header_bits + (gamma_k - 3.0);
                        let cost = c + header + extra_copies * res;
                        if best.get(&eid).is_none_or(|(cur, _, _, _)| cost < *cur) {
                            best.insert(
                                eid,
                                (
                                    cost,
                                    total_copies * res,
                                    leaves.saturating_mul(k),
                                    (*node).clone(),
                                ),
                            );
                            updated = true;
                        }
                    }
                } else if node.op == "seq" && node.children.len() == 2 {
                    let l_id = eg.find_imm(node.children[0]);
                    let r_id = eg.find_imm(node.children[1]);
                    if let (Some(&(c1, r1, l1, _)), Some(&(c2, r2, l2, _))) =
                        (best.get(&l_id), best.get(&r_id))
                    {
                        let cost = c1 + c2 + model.seq_op_bits;
                        let res = r1 + r2;
                        if best.get(&eid).is_none_or(|(cur, _, _, _)| cost < *cur) {
                            best.insert(eid, (cost, res, l1.saturating_add(l2), (*node).clone()));
                            updated = true;
                        }
                    }
                }
            }
        }
        if !updated {
            break;
        }
    }
    let root_id = eg.find_imm(root);
    let cost = best.get(&root_id)?.0;
    let mut visiting = HashSet::new();
    let expr = build_expr_from_best(eg, root_id, &best, &mut visiting)?;
    Some((cost, expr))
}

/// Collect all block labels wrapped in `Expr::Dup(...)` or `Expr::Repeat(..., k)`.
#[must_use]
pub fn collect_dups(expr: &Expr) -> Vec<String> {
    collect_dup_multiplicities(expr)
        .into_iter()
        .map(|(lbl, _)| lbl)
        .collect()
}

/// Collect `(block_label, copy_count)` for every block collapsed under `Expr::Dup` (`2`)
/// or `Expr::Repeat(_, k)` (`k`), multiplying nested repeat factors.
#[must_use]
pub fn collect_dup_multiplicities(expr: &Expr) -> Vec<(String, usize)> {
    fn walk(e: &Expr, mult: usize, out: &mut Vec<(String, usize)>) {
        match e {
            Expr::Block(lbl) => {
                if mult > 1 {
                    out.push((lbl.clone(), mult));
                }
            }
            Expr::Seq(a, b) => {
                walk(a, mult, out);
                walk(b, mult, out);
            }
            Expr::Inv(a) => walk(a, mult, out),
            Expr::Dup(inner) => walk(inner, mult.saturating_mul(2), out),
            Expr::Repeat(inner, k) => walk(inner, mult.saturating_mul(*k), out),
        }
    }
    let mut out = Vec::new();
    walk(expr, 1, &mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn has_dup_of(expr: &Expr, label: &str) -> bool {
        fn has_leaf(e: &Expr, label: &str) -> bool {
            match e {
                Expr::Block(lbl) => lbl == label,
                Expr::Seq(a, b) => has_leaf(a, label) || has_leaf(b, label),
                Expr::Inv(a) | Expr::Dup(a) | Expr::Repeat(a, _) => has_leaf(a, label),
            }
        }
        match expr {
            Expr::Dup(inner) | Expr::Repeat(inner, _) => has_leaf(inner, label),
            Expr::Seq(a, b) => has_dup_of(a, label) || has_dup_of(b, label),
            Expr::Inv(a) => has_dup_of(a, label),
            Expr::Block(_) => false,
        }
    }

    #[test]
    fn egraph_collapses_inverted_tandem_duplication() {
        // Ancestral A B C D -> tandem dup B -> invert [B, B, C] => A, -C, -B, -B, D
        let haplotype = [("A", 1), ("C", -1), ("B", -1), ("B", -1), ("D", 1)];
        let mut eg = EGraph::new();
        let root = build_seq_expr(&mut eg, &haplotype).expect("nonempty");
        saturate_synteny(&mut eg, 4, 800);
        let (_, expr) = extract_min_cost(&eg, root).expect("extracted");
        assert!(has_dup_of(&expr, "B"), "expected dup(B) in {expr:?}");
    }

    #[test]
    fn mdl_cost_model_gates_dup_by_sequence_divergence() {
        let haplotype = [("B", 1), ("B", 1)];
        let mut eg = EGraph::new();
        let root = build_seq_expr(&mut eg, &haplotype).expect("nonempty");
        saturate_synteny(&mut eg, 3, 200);

        let mut low_div = MdlCostModel::from_egraph(&eg);
        low_div.block_stats.insert("B".to_owned(), (1000, 0.02));
        let (cost_low, expr_low) = extract_min_mdl_cost(&eg, root, &low_div).expect("extracted");
        assert!(
            has_dup_of(&expr_low, "B"),
            "2% div should collapse dup: {expr_low:?}"
        );
        assert!(cost_low < 2500.0, "expected <2500 bits, got {cost_low}");

        let mut high_div = MdlCostModel::from_egraph(&eg);
        high_div.block_stats.insert("B".to_owned(), (1000, 0.75));
        let (_, expr_high) = extract_min_mdl_cost(&eg, root, &high_div).expect("extracted");
        assert!(
            !has_dup_of(&expr_high, "B"),
            "75% (JC-saturated) div should NOT collapse dup: {expr_high:?}"
        );
    }

    #[test]
    fn congruence_closure_propagates_to_parent_enodes() {
        let mut eg = EGraph::new();
        let a = eg.add(ENode {
            op: "blk:A".to_owned(),
            children: Vec::new(),
        });
        let b = eg.add(ENode {
            op: "blk:B".to_owned(),
            children: Vec::new(),
        });
        let inv_a = eg.add(ENode {
            op: "inv".to_owned(),
            children: vec![a],
        });
        let inv_inv_a = eg.add(ENode {
            op: "inv".to_owned(),
            children: vec![inv_a],
        });
        let seq_nested = eg.add(ENode {
            op: "seq".to_owned(),
            children: vec![inv_inv_a, b],
        });
        let seq_plain = eg.add(ENode {
            op: "seq".to_owned(),
            children: vec![a, b],
        });
        assert_ne!(eg.find(seq_nested), eg.find(seq_plain));

        eg.union(inv_inv_a, a);
        eg.rebuild();
        assert_eq!(eg.find(seq_nested), eg.find(seq_plain));
    }

    #[test]
    fn repair_cleans_stale_memo_when_both_children_merge_separately() {
        // Reviewer 2 #2 regression: merging a->a2 then b->b2 in separate rebuilds must not leave
        // stale intermediate `seq(a2, b)` in `eg.memo` or `eg.classes`.
        let mut eg = EGraph::new();
        let a1 = eg.add(ENode {
            op: "blk:A1".to_owned(),
            children: vec![],
        });
        let a2 = eg.add(ENode {
            op: "blk:A2".to_owned(),
            children: vec![],
        });
        let b1 = eg.add(ENode {
            op: "blk:B1".to_owned(),
            children: vec![],
        });
        let b2 = eg.add(ENode {
            op: "blk:B2".to_owned(),
            children: vec![],
        });
        let seq_ab = eg.add(ENode {
            op: "seq".to_owned(),
            children: vec![a1, b1],
        });

        eg.union(a2, a1);
        eg.rebuild();
        let ra = eg.find(a1);
        let stale_mid = ENode {
            op: "seq".to_owned(),
            children: vec![ra, b1],
        };
        assert!(eg.memo.contains_key(&stale_mid));

        eg.union(b2, b1);
        eg.rebuild();
        let rb = eg.find(b1);
        let final_node = ENode {
            op: "seq".to_owned(),
            children: vec![ra, rb],
        };
        let got = eg.memo.get(&final_node).copied();
        let expected = eg.find(seq_ab);
        assert_eq!(got.map(|id| eg.find(id)), Some(expected));
        if rb != b1 {
            assert!(
                !eg.memo.contains_key(&stale_mid),
                "stale intermediate e-node remained in memo after second child merged"
            );
        }
    }

    #[test]
    fn egraph_collapses_three_copy_tandem_array_into_repeat_3() {
        let haplotype = [("A", 1), ("B", 1), ("B", 1), ("B", 1), ("C", 1)];
        let mut eg = EGraph::new();
        let root = build_seq_expr(&mut eg, &haplotype).expect("nonempty");
        saturate_synteny(&mut eg, 6, 2048);

        let mut model = MdlCostModel::from_egraph(&eg);
        model.block_stats.insert("A".to_owned(), (500, 0.01));
        model.block_stats.insert("B".to_owned(), (500, 0.02));
        model.block_stats.insert("C".to_owned(), (500, 0.01));
        let (cost, expr) = extract_min_mdl_cost(&eg, root, &model).expect("extracted");

        let mults = collect_dup_multiplicities(&expr);
        assert!(
            mults.contains(&("B".to_owned(), 3)),
            "expected Repeat(B, 3) in {expr:?}, got mults={mults:?}"
        );
        assert!(
            collect_dups(&expr).contains(&"B".to_owned()),
            "collect_dups must include Repeat leaves: {expr:?}"
        );
        // 3 copies via rep:3 costs strictly less than encoding 3 separate copies
        let uncollapsed_b_only = 3.0 * model.block_bits("B") + 2.0 * model.seq_op_bits;
        assert!(
            cost < model.block_bits("A") + model.block_bits("C") + uncollapsed_b_only,
            "expected rep:3 cost ({cost}) < uncollapsed cost"
        );
    }

    #[test]
    fn egraph_collapses_inverted_three_copy_tandem_array() {
        let haplotype = [
            ("A", 1),
            ("C", -1),
            ("B", -1),
            ("B", -1),
            ("B", -1),
            ("D", 1),
        ];
        let mut eg = EGraph::new();
        let root = build_seq_expr(&mut eg, &haplotype).expect("nonempty");
        saturate_synteny(&mut eg, 6, 2048);

        let (_, expr) = extract_min_cost(&eg, root).expect("extracted");
        let mults = collect_dup_multiplicities(&expr);
        assert!(
            mults.contains(&("B".to_owned(), 3)),
            "expected inverted 3-copy Repeat(B, 3) in {expr:?}, got mults={mults:?}"
        );
    }

    #[test]
    fn mdl_cost_model_gates_repeat_3_by_divergence() {
        let haplotype = [("B", 1), ("B", 1), ("B", 1)];
        let mut eg = EGraph::new();
        let root = build_seq_expr(&mut eg, &haplotype).expect("nonempty");
        saturate_synteny(&mut eg, 5, 1024);

        let mut low_div = MdlCostModel::from_egraph(&eg);
        low_div.block_stats.insert("B".to_owned(), (1000, 0.02));
        let (_, expr_low) = extract_min_mdl_cost(&eg, root, &low_div).expect("low div");
        assert!(
            matches!(expr_low, Expr::Repeat(_, 3)),
            "2% div 3-copy array should collapse to Repeat(_, 3), got {expr_low:?}"
        );

        let mut high_div = MdlCostModel::from_egraph(&eg);
        high_div.block_stats.insert("B".to_owned(), (1000, 0.75));
        let (_, expr_high) = extract_min_mdl_cost(&eg, root, &high_div).expect("high div");
        assert!(
            collect_dup_multiplicities(&expr_high).is_empty(),
            "75% (JC-saturated) 3-copy array must not collapse, got {expr_high:?}"
        );
    }

    #[test]
    fn mdl_cost_model_rejects_dup_of_seq_at_jc69_saturation() {
        // 4-copy JC69-saturated array [B, B, B, B] must NOT collapse into Dup(Seq(B, B)).
        let haplotype = [("B", 1), ("B", 1), ("B", 1), ("B", 1)];
        let mut eg = EGraph::new();
        let root = build_seq_expr(&mut eg, &haplotype).expect("nonempty");
        saturate_synteny(&mut eg, 6, 1024);

        let mut high_div = MdlCostModel::from_egraph(&eg);
        high_div.block_stats.insert("B".to_owned(), (1000, 0.75));
        let (_, expr_high) = extract_min_mdl_cost(&eg, root, &high_div).expect("high div");
        assert!(
            collect_dup_multiplicities(&expr_high).is_empty(),
            "JC69-saturated 4-copy array must not bypass gate via Dup(Seq(B, B)), got {expr_high:?}"
        );
    }
}
