//! Anchor chaining (`lis_chain`, `cluster_strand_lcbs`, `merge_colinear_lcbs`, `IntervalSet`).

use crate::AlignError;

use super::canonical_spaced_seeds;

/// Compute the Longest Increasing Subsequence chain over `(x, y)` anchors.
#[must_use]
pub fn lis_chain(anchors: &[(usize, usize)]) -> Vec<(usize, usize)> {
    if anchors.is_empty() {
        return Vec::new();
    }
    let mut sorted = anchors.to_vec();
    sorted.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| b.1.cmp(&a.1)));

    let mut tails_y: Vec<usize> = Vec::with_capacity(sorted.len());
    let mut tail_idx: Vec<usize> = Vec::with_capacity(sorted.len());
    let mut prev = vec![usize::MAX; sorted.len()];

    for (idx, &(_, y)) in sorted.iter().enumerate() {
        let pos = tails_y.partition_point(|&ty| ty < y);
        if pos == tails_y.len() {
            tails_y.push(y);
            tail_idx.push(idx);
        } else {
            tails_y[pos] = y;
            tail_idx[pos] = idx;
        }
        if pos > 0 {
            prev[idx] = tail_idx[pos - 1];
        }
    }

    let mut chain = Vec::with_capacity(tail_idx.len());
    let mut cur = tail_idx[tail_idx.len() - 1];
    while cur != usize::MAX {
        chain.push(sorted[cur]);
        cur = prev[cur];
    }
    chain.reverse();
    chain
}

/// Cluster `(x, y)` hits on one strand into colinear LCB chains.
///
/// Uses `BTreeSet` indices on `(lx, id)` for `max_gap` expiration and `(ly, id)`
/// for $O(\log |\text{active}|)$ local `y`-window range queries instead of
/// scanning all active chains across the chromosome.
#[must_use]
pub fn cluster_strand_lcbs(
    hits: &[(usize, usize)],
    max_gap: usize,
    max_diag_drift: usize,
    min_anchors: usize,
    min_span_bp: usize,
    max_mean_step: f64,
) -> Vec<Vec<(usize, usize)>> {
    use std::collections::BTreeSet;

    if hits.is_empty() {
        return Vec::new();
    }
    let mut sorted = hits.to_vec();
    sorted.sort_unstable();

    let keep = |ch: &[(usize, usize)]| -> bool {
        let (Some(first), Some(last)) = (ch.first(), ch.last()) else {
            return false;
        };
        let span_x = last.0.saturating_sub(first.0);
        #[allow(clippy::cast_precision_loss)]
        let mean_step = (span_x as f64) / (ch.len() as f64);
        ch.len() >= min_anchors && span_x >= min_span_bp && mean_step <= max_mean_step
    };

    let mut chains: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut tails: Vec<(usize, usize)> = Vec::new();
    let mut by_x: BTreeSet<(usize, usize)> = BTreeSet::new();
    let mut by_y: BTreeSet<(usize, usize)> = BTreeSet::new();
    let mut finished: Vec<Vec<(usize, usize)>> = Vec::new();

    for &(x, y) in &sorted {
        let x_floor = x.saturating_sub(max_gap);
        while let Some(&(lx, id)) = by_x.first() {
            if lx >= x_floor {
                break;
            }
            by_x.remove(&(lx, id));
            let (_, ly) = tails[id];
            by_y.remove(&(ly, id));
            let ch = std::mem::take(&mut chains[id]);
            if keep(&ch) {
                finished.push(ch);
            }
        }

        let y_floor = y.saturating_sub(max_gap);
        let mut best_id: Option<usize> = None;
        let mut best_cost = usize::MAX;
        for &(_, id) in by_y.range((y_floor, 0)..(y, 0)) {
            let (lx, ly) = tails[id];
            if x > lx && y > ly {
                let dx = x - lx;
                let dy = y - ly;
                let drift = dx.abs_diff(dy);
                if dx <= max_gap && dy <= max_gap && drift <= max_diag_drift {
                    let cost = drift + dx;
                    if cost < best_cost {
                        best_cost = cost;
                        best_id = Some(id);
                    }
                }
            }
        }

        if let Some(id) = best_id {
            let (prev_x, prev_y) = tails[id];
            by_x.remove(&(prev_x, id));
            by_y.remove(&(prev_y, id));
            chains[id].push((x, y));
            tails[id] = (x, y);
            by_x.insert((x, id));
            by_y.insert((y, id));
        } else {
            let id = chains.len();
            chains.push(vec![(x, y)]);
            tails.push((x, y));
            by_x.insert((x, id));
            by_y.insert((y, id));
        }
    }

    for (_, id) in by_x {
        let ch = std::mem::take(&mut chains[id]);
        if keep(&ch) {
            finished.push(ch);
        }
    }
    finished
}

/// Find local `w=9` spaced-seed anchors with 2-hit diagonal consistency.
///
/// # Errors
///
/// Returns [`AlignError::InvalidMask`] if `mask` is invalid.
pub fn local_w9_subchain(
    sub_x: &[u8],
    sub_y: &[u8],
    off_x: usize,
    off_y: usize,
    mask: &str,
    max_diag_diff: usize,
    max_step_bp: usize,
) -> Result<Vec<(usize, usize)>, AlignError> {
    let ix = canonical_spaced_seeds(sub_x, mask, 2)?;
    let iy = canonical_spaced_seeds(sub_y, mask, 2)?;
    let mut hits = Vec::new();
    for (k, pxs) in &ix {
        if let Some(pys) = iy.get(k) {
            for &px in pxs {
                for &py in pys {
                    hits.push((px, py));
                }
            }
        }
    }
    let chain = lis_chain(&hits);
    let mut kept = Vec::new();

    let max_dd = max_diag_diff as u128;
    for (idx, &(x, y)) in chain.iter().enumerate() {
        let diag = (x as i128) - (y as i128);
        let has_prev = idx > 0 && {
            let (px, py) = chain[idx - 1];
            let pdiag = (px as i128) - (py as i128);
            (diag - pdiag).unsigned_abs() <= max_dd && (x - px) <= max_step_bp
        };
        let has_next = idx + 1 < chain.len() && {
            let (nx, ny) = chain[idx + 1];
            let ndiag = (nx as i128) - (ny as i128);
            (diag - ndiag).unsigned_abs() <= max_dd && (nx - x) <= max_step_bp
        };
        if has_prev || has_next {
            kept.push((off_x + x, off_y + y));
        }
    }
    Ok(kept)
}

pub(crate) fn merge_colinear_lcbs(
    lcb_list: &[Vec<(usize, usize)>],
    max_local_rescue_bp: usize,
) -> Vec<Vec<(usize, usize)>> {
    let mut ordered: Vec<Vec<(usize, usize)>> =
        lcb_list.iter().filter(|c| !c.is_empty()).cloned().collect();
    if ordered.is_empty() {
        return Vec::new();
    }
    ordered.sort_by_key(|c| c[0].0);
    let mut iter = ordered.into_iter();
    let Some(first) = iter.next() else {
        return Vec::new();
    };
    let mut merged = vec![first];
    for nxt in iter {
        let Some(&(px, py)) = merged.last().and_then(|last| last.last()) else {
            merged.push(nxt);
            continue;
        };
        let Some(&(nx, ny)) = nxt.first() else {
            continue;
        };
        if nx > px && ny > py {
            let dx = nx - px;
            let dy = ny - py;
            if dx <= max_local_rescue_bp
                && dy <= max_local_rescue_bp
                && dx.abs_diff(dy) <= 4000
                && let Some(last_mut) = merged.last_mut()
            {
                last_mut.extend(nxt);
                continue;
            }
        }
        merged.push(nxt);
    }
    merged
}

/// Disjoint closed-interval set `[sx, ex]` backed by `BTreeMap` for `O(log K)` overlap checks.
pub(crate) struct IntervalSet {
    map: std::collections::BTreeMap<usize, usize>,
}

impl IntervalSet {
    pub(crate) fn new() -> Self {
        Self {
            map: std::collections::BTreeMap::new(),
        }
    }

    pub(crate) fn from_lcbs<'a, I>(lcbs: I) -> Self
    where
        I: IntoIterator<Item = &'a Vec<(usize, usize)>>,
    {
        let mut set = Self::new();
        for c in lcbs {
            if let (Some(first), Some(last)) = (c.first(), c.last()) {
                set.insert_merge(first.0, last.0);
            }
        }
        set
    }

    fn insert_merge(&mut self, sx: usize, ex: usize) {
        let (mut sx, mut ex) = (sx.min(ex), sx.max(ex));
        while let Some((&s, &e)) = self.map.range(..=ex).next_back() {
            if e < sx {
                break;
            }
            self.map.remove(&s);
            sx = sx.min(s);
            ex = ex.max(e);
        }
        self.map.insert(sx, ex);
    }

    pub(crate) fn insert_if_disjoint(&mut self, sx: usize, ex: usize) -> bool {
        let (sx, ex) = (sx.min(ex), sx.max(ex));
        if let Some((_, &e)) = self.map.range(..=ex).next_back()
            && e >= sx
        {
            return false;
        }
        self.map.insert(sx, ex);
        true
    }
}

pub(crate) fn push_disjoint_lcbs(
    target: &mut super::LcbList,
    existing: &[Vec<(usize, usize)>],
    mut candidates: super::LcbList,
) {
    candidates.sort_unstable_by_key(|lcb| {
        let span = lcb
            .last()
            .zip(lcb.first())
            .map_or(0, |(last, first)| last.0.saturating_sub(first.0));
        std::cmp::Reverse((lcb.len(), span))
    });
    let mut intervals = IntervalSet::from_lcbs(existing.iter().chain(target.iter()));
    for lcb in candidates {
        if let (Some(first), Some(last)) = (lcb.first(), lcb.last())
            && intervals.insert_if_disjoint(first.0, last.0)
        {
            target.push(lcb);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{IntervalSet, lis_chain, merge_colinear_lcbs};

    #[test]
    fn lis_chain_enforces_strict_x_colinearity_on_duplicate_x() {
        let anchors = vec![(10, 50), (10, 100), (20, 120)];
        let chain = lis_chain(&anchors);
        assert_eq!(chain.len(), 2);
        assert!(chain[0].0 < chain[1].0 && chain[0].1 < chain[1].1);
    }

    #[test]
    fn merge_colinear_lcbs_merges_nearby_diagonal_blocks() {
        let a = vec![(0, 0), (15, 15)];
        let b = vec![(30, 30), (45, 45)];
        let c = vec![(500, 10), (515, 25)];
        let merged = merge_colinear_lcbs(&[b, c, a], 50);
        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0], vec![(0, 0), (15, 15), (30, 30), (45, 45)]);
        assert_eq!(merged[1], vec![(500, 10), (515, 25)]);
    }

    #[test]
    fn interval_set_merges_overlaps_and_rejects_non_disjoint() {
        let lcbs = [vec![(10, 0), (20, 10)], vec![(15, 0), (30, 15)]];
        let mut set = IntervalSet::from_lcbs(lcbs.iter());
        assert!(!set.insert_if_disjoint(25, 35));
        assert!(!set.insert_if_disjoint(5, 10));
        assert!(set.insert_if_disjoint(31, 40));
        assert!(set.insert_if_disjoint(0, 9));
    }
}
