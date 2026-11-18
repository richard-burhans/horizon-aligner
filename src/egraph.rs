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
    pub fn find(&mut self, mut eid: usize) -> usize {
        while self.parent[eid] != eid {
            let p = self.parent[eid];
            self.parent[eid] = self.parent[p];
            eid = self.parent[eid];
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

    fn repair(&mut self, mut eid: usize) {
        eid = self.find(eid);
        let Some(nodes) = self.classes.remove(&eid) else {
            return;
        };
        let mut new_set = HashSet::new();
        for node in nodes {
            self.memo.remove(&node);
            let cnode = self.canonicalize(&node);
            if let Some(&other_raw) = self.memo.get(&cnode) {
                let other = self.find(other_raw);
                if other != eid {
                    // Re-insert partially built set before union.
                    self.classes.entry(eid).or_default().extend(new_set);
                    eid = self.union(eid, other);
                    new_set = HashSet::new();
                }
            }
            self.memo.insert(cnode.clone(), eid);
            new_set.insert(cnode);
        }
        self.classes.entry(eid).or_default().extend(new_set);
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
        let snapshot: Vec<(usize, Vec<ENode>)> = eg
            .classes
            .iter()
            .map(|(&id, set)| (id, set.iter().cloned().collect()))
            .collect();

        for (eid, nodes) in snapshot {
            if eg.memo.len() >= max_nodes {
                break;
            }
            for node in nodes {
                if node.op == "inv" && node.children.len() == 1 {
                    let child = eg.find(node.children[0]);
                    let cnodes: Vec<ENode> = eg
                        .classes
                        .get(&child)
                        .map_or_default(|s| s.iter().cloned().collect());
                    for cnode in cnodes {
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
                    let anodes: Vec<ENode> = eg
                        .classes
                        .get(&a)
                        .map_or_default(|s| s.iter().cloned().collect());
                    for anode in anodes {
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
                    let bnodes: Vec<ENode> = eg
                        .classes
                        .get(&b)
                        .map_or_default(|s| s.iter().cloned().collect());
                    for bnode in bnodes {
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

/// Extract minimum-parsimony synteny expression rooted at `root`.
#[must_use]
pub fn extract_min_cost(eg: &mut EGraph, root: usize) -> Option<(f64, Expr)> {
    let mut best: HashMap<usize, (f64, Expr)> = HashMap::new();
    let max_rounds = eg.classes.len() + 2;
    let class_ids: Vec<usize> = eg.classes.keys().copied().collect();

    for _ in 0..max_rounds {
        let mut updated = false;
        for raw_id in &class_ids {
            let eid = eg.find(*raw_id);
            let nodes: Vec<ENode> = eg
                .classes
                .get(&eid)
                .map_or_default(|s| s.iter().cloned().collect());
            for node in nodes {
                let cand: Option<(f64, Expr)> = if let Some(lbl) = node.op.strip_prefix("blk:") {
                    Some((1.0, Expr::Block(lbl.to_owned())))
                } else if node.op == "inv" && node.children.len() == 1 {
                    best.get(&eg.find(node.children[0]))
                        .map(|(c, e)| (c + 1.2, Expr::Inv(Box::new(e.clone()))))
                } else if node.op == "dup" && node.children.len() == 1 {
                    best.get(&eg.find(node.children[0]))
                        .map(|(c, e)| (c + 1.3, Expr::Dup(Box::new(e.clone()))))
                } else if node.op == "seq" && node.children.len() == 2 {
                    let left = best.get(&eg.find(node.children[0]));
                    let right = best.get(&eg.find(node.children[1]));
                    match (left, right) {
                        (Some((c1, e1)), Some((c2, e2))) => Some((
                            c1 + c2 + 1.0,
                            Expr::Seq(Box::new(e1.clone()), Box::new(e2.clone())),
                        )),
                        _ => None,
                    }
                } else {
                    None
                };

                if let Some((cost, expr)) = cand {
                    let replace = best.get(&eid).is_none_or(|(cur_cost, _)| cost < *cur_cost);
                    if replace {
                        best.insert(eid, (cost, expr));
                        updated = true;
                    }
                }
            }
        }
        if !updated {
            break;
        }
    }
    best.remove(&eg.find(root))
}

/// Collect all block labels wrapped in `Expr::Dup(...)`.
#[must_use]
pub fn collect_dups(expr: &Expr) -> Vec<String> {
    fn walk(e: &Expr, out: &mut Vec<String>) {
        match e {
            Expr::Dup(inner) => match inner.as_ref() {
                Expr::Block(lbl) => out.push(lbl.clone()),
                Expr::Inv(b) => {
                    if let Expr::Block(lbl) = b.as_ref() {
                        out.push(lbl.clone());
                    } else {
                        walk(b, out);
                    }
                }
                _ => walk(inner, out),
            },
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
        match expr {
            Expr::Dup(inner) => match inner.as_ref() {
                Expr::Block(lbl) => lbl == label,
                Expr::Inv(b) => matches!(b.as_ref(), Expr::Block(lbl) if lbl == label),
                _ => has_dup_of(inner, label),
            },
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
        let (_, expr) = extract_min_cost(&mut eg, root).expect("extracted");
        assert!(has_dup_of(&expr, "B"), "expected dup(B) in {expr:?}");
    }
}
