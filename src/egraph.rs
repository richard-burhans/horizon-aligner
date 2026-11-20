//! Congruence-closed e-graph for resolving tandem duplications and inversions across LCBs.

use std::collections::{HashMap, HashSet};

/// Immutable e-node referencing canonical child e-class IDs.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
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
    /// Collapsed tandem duplication.
    Dup(Box<Expr>),
}

/// Minimal congruence-closed e-graph (Willsey et al. 2021).
#[derive(Debug, Default)]
pub struct EGraph {
    parent: Vec<usize>,
    classes: HashMap<usize, HashSet<ENode>>,
    parents: Vec<Vec<(ENode, usize)>>,
    memo: HashMap<ENode, usize>,
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
        for &child in &cnode.children {
            self.parents[child].push((cnode.clone(), eid));
        }
        let mut set = HashSet::new();
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
        let len_a = self.classes.get(&ra).map_or(0, HashSet::len);
        let len_b = self.classes.get(&rb).map_or(0, HashSet::len);
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
            let todo: HashSet<usize> = drained.into_iter().map(|id| self.find(id)).collect();
            for eid in todo {
                self.repair(eid);
            }
        }
    }

    fn repair(&mut self, eid: usize) {
        let eid = self.find(eid);
        let parent_entries = std::mem::take(&mut self.parents[eid]);
        let mut new_parents: HashMap<ENode, usize> = HashMap::new();
        for (pnode, peid) in parent_entries {
            self.memo.remove(&pnode);
            let peid = self.find(peid);
            let cpnode = self.canonicalize(&pnode);
            if let Some(nodes) = self.classes.get_mut(&peid) {
                nodes.remove(&pnode);
                nodes.insert(cpnode.clone());
            }
            if let Some(&other_raw) = self.memo.get(&cpnode) {
                let other = self.find(other_raw);
                if other != peid {
                    self.union(peid, other);
                }
            }
            let peid = self.find(peid);
            self.memo.insert(cpnode.clone(), peid);
            if let Some(&prev_peid) = new_parents.get(&cpnode) {
                let prev = self.find(prev_peid);
                if prev != peid {
                    let merged = self.union(prev, peid);
                    new_parents.insert(cpnode, merged);
                    continue;
                }
            }
            new_parents.insert(cpnode, peid);
        }
        let final_eid = self.find(eid);
        for (cpnode, peid) in new_parents {
            let rep = self.find(peid);
            self.parents[final_eid].push((cpnode, rep));
        }
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

/// Run bounded equality saturation for inversion cancellation, reversal, tandem collapse, and associativity.
#[allow(clippy::too_many_lines)]
pub fn saturate_synteny(eg: &mut EGraph, max_iters: usize, max_nodes: usize) {
    for _ in 0..max_iters {
        if eg.memo.len() >= max_nodes {
            break;
        }
        let mut merges: Vec<(usize, usize)> = Vec::new();
        let mut snapshot: Vec<(usize, Vec<ENode>)> = eg
            .classes
            .iter()
            .map(|(&id, set)| {
                let mut nodes: Vec<ENode> = set.iter().cloned().collect();
                nodes.sort_unstable_by(|a, b| {
                    a.op.cmp(&b.op).then_with(|| a.children.cmp(&b.children))
                });
                (id, nodes)
            })
            .collect();
        snapshot.sort_unstable_by_key(|(id, _)| *id);

        for (eid, nodes) in snapshot {
            if eg.memo.len() >= max_nodes {
                break;
            }
            for node in nodes {
                if eg.memo.len() >= max_nodes {
                    break;
                }
                if node.op == "inv" && node.children.len() == 1 {
                    let child = eg.find(node.children[0]);
                    let mut cnodes: Vec<ENode> = eg
                        .classes
                        .get(&child)
                        .map_or_default(|s| s.iter().cloned().collect());
                    cnodes.sort_unstable_by(|a, b| {
                        a.op.cmp(&b.op).then_with(|| a.children.cmp(&b.children))
                    });
                    for cnode in cnodes {
                        if eg.memo.len() >= max_nodes {
                            break;
                        }
                        if cnode.op == "inv" && cnode.children.len() == 1 {
                            merges.push((eid, cnode.children[0]));
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
                    let mut anodes: Vec<ENode> = eg
                        .classes
                        .get(&a)
                        .map_or_default(|s| s.iter().cloned().collect());
                    anodes.sort_unstable_by(|x, y| {
                        x.op.cmp(&y.op).then_with(|| x.children.cmp(&y.children))
                    });
                    for anode in anodes {
                        if eg.memo.len() >= max_nodes {
                            break;
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
                    let mut bnodes: Vec<ENode> = eg
                        .classes
                        .get(&b)
                        .map_or_default(|s| s.iter().cloned().collect());
                    bnodes.sort_unstable_by(|x, y| {
                        x.op.cmp(&y.op).then_with(|| x.children.cmp(&y.children))
                    });
                    for bnode in bnodes {
                        if eg.memo.len() >= max_nodes {
                            break;
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

    /// Incremental bit cost of collapsing a second copy via `dup(expr)`.
    #[must_use]
    pub fn dup_cost(&self, child_cost: f64, child_expr: &Expr) -> f64 {
        let residual_bits = if self.block_stats.is_empty() {
            0.0
        } else {
            Self::expr_residual_bits(child_expr, &self.block_stats)
        };
        child_cost + self.dup_header_bits + residual_bits
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

    fn expr_residual_bits(expr: &Expr, stats: &HashMap<String, (usize, f64)>) -> f64 {
        match expr {
            Expr::Block(lbl) => {
                if let Some(&(len_bp, div)) = stats.get(lbl) {
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
            Expr::Inv(inner) | Expr::Dup(inner) => Self::expr_residual_bits(inner, stats),
            Expr::Seq(a, b) => {
                Self::expr_residual_bits(a, stats) + Self::expr_residual_bits(b, stats)
            }
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
    best: &HashMap<usize, (f64, f64, ENode)>,
    visiting: &mut HashSet<usize>,
) -> Option<Expr> {
    if !visiting.insert(eid) {
        return None;
    }
    let (_, _, node) = best.get(&eid)?;
    let expr = if let Some(lbl) = node.op.strip_prefix("blk:") {
        Some(Expr::Block(lbl.to_owned()))
    } else if node.op == "inv" && node.children.len() == 1 {
        let cid = eg.find_imm(node.children[0]);
        build_expr_from_best(eg, cid, best, visiting).map(|e| Expr::Inv(Box::new(e)))
    } else if node.op == "dup" && node.children.len() == 1 {
        let cid = eg.find_imm(node.children[0]);
        build_expr_from_best(eg, cid, best, visiting).map(|e| Expr::Dup(Box::new(e)))
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
pub fn extract_min_mdl_cost(
    eg: &EGraph,
    root: usize,
    model: &MdlCostModel,
) -> Option<(f64, Expr)> {
    if root >= eg.parent.len() {
        return None;
    }
    let mut best: HashMap<usize, (f64, f64, ENode)> = HashMap::new();
    let max_rounds = eg.classes.len() + 2;
    let mut class_ids: Vec<usize> = eg.classes.keys().copied().collect();
    class_ids.sort_unstable();

    for _ in 0..max_rounds {
        let mut updated = false;
        for raw_id in &class_ids {
            let eid = eg.find_imm(*raw_id);
            let Some(node_set) = eg.classes.get(&eid) else {
                continue;
            };
            let mut nodes: Vec<&ENode> = node_set.iter().collect();
            nodes
                .sort_unstable_by(|a, b| a.op.cmp(&b.op).then_with(|| a.children.cmp(&b.children)));
            for node in nodes {
                if let Some(lbl) = node.op.strip_prefix("blk:") {
                    let cost = model.block_bits(lbl);
                    let res = model.block_residual_bits(lbl);
                    if best.get(&eid).is_none_or(|(cur, _, _)| cost < *cur) {
                        best.insert(eid, (cost, res, (*node).clone()));
                        updated = true;
                    }
                } else if node.op == "inv" && node.children.len() == 1 {
                    let cid = eg.find_imm(node.children[0]);
                    if let Some(&(c, res, _)) = best.get(&cid) {
                        let cost = c + model.inv_op_bits;
                        if best.get(&eid).is_none_or(|(cur, _, _)| cost < *cur) {
                            best.insert(eid, (cost, res, (*node).clone()));
                            updated = true;
                        }
                    }
                } else if node.op == "dup" && node.children.len() == 1 {
                    let cid = eg.find_imm(node.children[0]);
                    if let Some(&(c, res, _)) = best.get(&cid) {
                        let cost = c + model.dup_header_bits + res;
                        if best.get(&eid).is_none_or(|(cur, _, _)| cost < *cur) {
                            best.insert(eid, (cost, res, (*node).clone()));
                            updated = true;
                        }
                    }
                } else if node.op == "seq" && node.children.len() == 2 {
                    let l_id = eg.find_imm(node.children[0]);
                    let r_id = eg.find_imm(node.children[1]);
                    if let (Some(&(c1, r1, _)), Some(&(c2, r2, _))) =
                        (best.get(&l_id), best.get(&r_id))
                    {
                        let cost = c1 + c2 + model.seq_op_bits;
                        let res = r1 + r2;
                        if best.get(&eid).is_none_or(|(cur, _, _)| cost < *cur) {
                            best.insert(eid, (cost, res, (*node).clone()));
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

/// Collect all block labels wrapped in `Expr::Dup(...)`.
#[must_use]
pub fn collect_dups(expr: &Expr) -> Vec<String> {
    fn collect_leaves(e: &Expr, out: &mut Vec<String>) {
        match e {
            Expr::Block(lbl) => out.push(lbl.clone()),
            Expr::Seq(a, b) => {
                collect_leaves(a, out);
                collect_leaves(b, out);
            }
            Expr::Inv(a) | Expr::Dup(a) => collect_leaves(a, out),
        }
    }
    fn walk(e: &Expr, out: &mut Vec<String>) {
        match e {
            Expr::Dup(inner) => collect_leaves(inner, out),
            Expr::Seq(a, b) => {
                walk(a, out);
                walk(b, out);
            }
            Expr::Inv(a) => walk(a, out),
            Expr::Block(_) => {}
        }
    }
    let mut out = Vec::new();
    walk(expr, &mut out);
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
                Expr::Inv(a) | Expr::Dup(a) => has_leaf(a, label),
            }
        }
        match expr {
            Expr::Dup(inner) => has_leaf(inner, label),
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
        let (cost_low, expr_low) =
            extract_min_mdl_cost(&eg, root, &low_div).expect("extracted");
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
}
