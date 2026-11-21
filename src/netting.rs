//! Colinear block/run-group partitioning and 1-to-1 interval netting.

use crate::mdl::AlignedRun;

/// Partition aligned pairs into colinear runs of at least 12 bases (sorting by `x` if unsorted).
pub(crate) fn partition_colinear_blocks(
    pairs: &[(usize, usize)],
) -> Vec<(Vec<(usize, usize)>, bool)> {
    type ActiveChain = (Vec<(usize, usize)>, Option<bool>);
    let sorted_storage: Vec<(usize, usize)>;
    let ordered: &[(usize, usize)] = if pairs.windows(2).all(|w| w[0].0 <= w[1].0) {
        pairs
    } else {
        let mut tmp = pairs.to_vec();
        tmp.sort_unstable();
        sorted_storage = tmp;
        &sorted_storage
    };
    let mut blocks: Vec<(Vec<(usize, usize)>, bool)> = Vec::new();
    let mut active: Vec<ActiveChain> = Vec::new();

    for &(x, y) in ordered {
        let mut i = 0;
        while i < active.len() {
            let px = active[i].0.last().map_or(x, |&(px, _)| px);
            if x.saturating_sub(px) > 20 {
                let (chain, rev_opt) = active.swap_remove(i);
                if chain.len() >= 12 {
                    blocks.push((chain, rev_opt.unwrap_or(false)));
                }
            } else {
                i += 1;
            }
        }
        let mut placed = false;
        for (chain, rev_opt) in &mut active {
            let (px, py) = chain[chain.len() - 1];
            let dx = x - px;
            let (dy, rev_step) = if y >= py {
                (y - py, false)
            } else {
                (py - y, true)
            };
            if dx > 0
                && dy > 0
                && dx <= 20
                && dy <= 20
                && rev_opt.is_none_or(|r| r == rev_step)
                && dx.abs_diff(dy) <= 4
            {
                chain.push((x, y));
                if rev_opt.is_none() {
                    *rev_opt = Some(rev_step);
                }
                placed = true;
                break;
            }
        }
        if !placed {
            if active.len() >= 32
                && let Some((oldest_idx, _)) = active
                    .iter()
                    .enumerate()
                    .min_by_key(|(_, (c, _))| c.last().map_or(0, |&(px, _)| px))
            {
                let (chain, rev_opt) = active.swap_remove(oldest_idx);
                if chain.len() >= 12 {
                    blocks.push((chain, rev_opt.unwrap_or(false)));
                }
            }
            active.push((vec![(x, y)], None));
        }
    }

    for (chain, rev_opt) in active {
        if chain.len() >= 12 {
            blocks.push((chain, rev_opt.unwrap_or(false)));
        }
    }
    blocks.sort_unstable_by_key(|(b, _)| (b[0].0, b[0].1));
    blocks
}

fn free_subintervals(
    occupied: &std::collections::BTreeMap<usize, usize>,
    lo: usize,
    hi: usize,
    out: &mut Vec<(usize, usize)>,
) {
    out.clear();
    if lo >= hi {
        return;
    }
    let mut cur = lo;
    if let Some((&_s, &e)) = occupied.range(..lo).next_back()
        && e > cur
    {
        cur = e;
    }
    for (&s, &e) in occupied.range(lo..hi) {
        if cur >= hi {
            break;
        }
        if s > cur {
            out.push((cur, s.min(hi)));
        }
        cur = cur.max(e);
    }
    if cur < hi {
        out.push((cur, hi));
    }
}

fn mark_occupied(
    occupied: &mut std::collections::BTreeMap<usize, usize>,
    mut lo: usize,
    mut hi: usize,
) {
    if lo >= hi {
        return;
    }
    let mut to_remove = Vec::new();
    if let Some((&ps, &pe)) = occupied.range(..=lo).next_back()
        && pe >= lo
    {
        lo = lo.min(ps);
        hi = hi.max(pe);
        to_remove.push(ps);
    }
    for (&s, &e) in occupied.range(lo..=hi) {
        lo = lo.min(s);
        hi = hi.max(e);
        to_remove.push(s);
    }
    for k in to_remove {
        occupied.remove(&k);
    }
    occupied.insert(lo, hi);
}

/// Net overlapping [`AlignedRun`]s within a contig pair to strict 1-to-1 coverage
/// by prioritizing the longest colinear run-groups first and clipping overlapping bases.
#[must_use]
pub fn net_one_to_one_runs(runs: &[AlignedRun], len_x: usize, len_y: usize) -> Vec<AlignedRun> {
    use std::collections::BTreeMap;
    if runs.is_empty() || len_x == 0 || len_y == 0 {
        return Vec::new();
    }
    let mut blocks: Vec<(usize, Vec<AlignedRun>)> = partition_colinear_run_groups(runs)
        .into_iter()
        .map(|(grp, _)| {
            let total_bp: usize = grp.iter().map(|r| r.len).sum();
            (total_bp, grp)
        })
        .collect();
    blocks.sort_unstable_by_key(|(bp, grp)| (std::cmp::Reverse(*bp), grp[0].x_start));

    let mut used_x: BTreeMap<usize, usize> = BTreeMap::new();
    let mut used_y: BTreeMap<usize, usize> = BTreeMap::new();
    let mut free_x = Vec::new();
    let mut free_y = Vec::new();
    let mut kept: Vec<AlignedRun> = Vec::with_capacity(runs.len());

    for (_, grp) in blocks {
        for r in grp {
            if r.len == 0 || r.x_start >= len_x {
                continue;
            }
            let Some(y_end_rc) = r.y_start.checked_add(1) else {
                continue;
            };
            let valid_len = if r.is_rc {
                r.len.min(y_end_rc).min(len_x - r.x_start)
            } else {
                r.len
                    .min(len_x - r.x_start)
                    .min(len_y.saturating_sub(r.y_start))
            };
            if valid_len < 6 {
                continue;
            }
            let x_lo = r.x_start;
            let x_hi = r.x_start + valid_len;
            free_subintervals(&used_x, x_lo, x_hi, &mut free_x);
            if free_x.is_empty() {
                continue;
            }
            let (y_lo, y_hi) = if r.is_rc {
                (y_end_rc - valid_len, y_end_rc.min(len_y))
            } else {
                (r.y_start, (r.y_start + valid_len).min(len_y))
            };
            free_subintervals(&used_y, y_lo, y_hi, &mut free_y);
            if free_y.is_empty() {
                continue;
            }

            // Intersect free_x (projected to k-offsets) with free_y (projected to k-offsets).
            let mut ix = 0usize;
            let mut iy = 0usize;
            while ix < free_x.len() && iy < free_y.len() {
                let kx0 = free_x[ix].0 - r.x_start;
                let kx1 = free_x[ix].1 - r.x_start;
                let (ky0, ky1) = if r.is_rc {
                    let rev_idx = free_y.len() - 1 - iy;
                    let (ya, yb) = free_y[rev_idx];
                    (y_end_rc - yb, y_end_rc - ya)
                } else {
                    let (ya, yb) = free_y[iy];
                    (ya - r.y_start, yb - r.y_start)
                };
                let k0 = kx0.max(ky0);
                let k1 = kx1.min(ky1);
                if k1 >= k0 + 6 {
                    let sub_len = k1 - k0;
                    let sx = r.x_start + k0;
                    let sy = if r.is_rc {
                        r.y_start - k0
                    } else {
                        r.y_start + k0
                    };
                    mark_occupied(&mut used_x, sx, sx + sub_len);
                    if r.is_rc {
                        mark_occupied(&mut used_y, sy + 1 - sub_len, sy + 1);
                    } else {
                        mark_occupied(&mut used_y, sy, sy + sub_len);
                    }
                    kept.push(AlignedRun {
                        x_start: sx,
                        y_start: sy,
                        len: sub_len,
                        is_rc: r.is_rc,
                    });
                }
                if kx1 <= ky1 {
                    ix += 1;
                } else {
                    iy += 1;
                }
            }
        }
    }
    kept.sort_unstable_by_key(|r| (r.x_start, r.y_start, r.is_rc));
    kept
}

/// Partition [`AlignedRun`]s into colinear block groups of at least 12 bases.
pub(crate) fn partition_colinear_run_groups(runs: &[AlignedRun]) -> Vec<(Vec<AlignedRun>, bool)> {
    let ordered: std::borrow::Cow<'_, [AlignedRun]> =
        if runs.windows(2).all(|w| w[0].x_start <= w[1].x_start) {
            std::borrow::Cow::Borrowed(runs)
        } else {
            let mut copy = runs.to_vec();
            copy.sort_unstable_by_key(|r| (r.x_start, r.y_start, r.is_rc));
            std::borrow::Cow::Owned(copy)
        };

    let mut blocks: Vec<(Vec<AlignedRun>, bool)> = Vec::new();
    let mut active: Vec<Vec<AlignedRun>> = Vec::new();
    let retire_chain = |ch: Vec<AlignedRun>, out: &mut Vec<(Vec<AlignedRun>, bool)>| {
        let total_bp: usize = ch.iter().map(|x| x.len).sum();
        if total_bp >= 12 {
            let is_rev = ch[0].is_rc;
            out.push((ch, is_rev));
        }
    };
    for &r in ordered.iter() {
        if r.len == 0 || (r.is_rc && r.y_start.checked_add(1).is_none_or(|y1| y1 < r.len)) {
            continue;
        }
        let mut idx = 0usize;
        while idx < active.len() {
            let last = active[idx][active[idx].len() - 1];
            if last.x_start + last.len + 40 < r.x_start {
                let ch = active.swap_remove(idx);
                retire_chain(ch, &mut blocks);
            } else {
                idx += 1;
            }
        }
        let mut placed = false;
        for ch in &mut active {
            let prev = ch[ch.len() - 1];
            if prev.is_rc != r.is_rc {
                continue;
            }
            let prev_x_end = prev.x_start + prev.len;
            if r.x_start < prev_x_end || r.x_start - prev_x_end > 40 {
                continue;
            }
            let ok_y = if r.is_rc {
                let Some(prev_y_lim) = prev
                    .y_start
                    .checked_add(1)
                    .and_then(|y1| y1.checked_sub(prev.len))
                else {
                    continue;
                };
                r.y_start < prev_y_lim && (prev_y_lim - 1 - r.y_start) <= 40
            } else {
                let prev_y_lim = prev.y_start + prev.len;
                r.y_start >= prev_y_lim && (r.y_start - prev_y_lim) <= 40
            };
            if ok_y {
                ch.push(r);
                placed = true;
                break;
            }
        }
        if !placed {
            if active.len() >= 32
                && let Some((oldest_idx, _)) = active.iter().enumerate().min_by_key(|(_, ch)| {
                    let last = ch[ch.len() - 1];
                    last.x_start + last.len
                })
            {
                let ch = active.swap_remove(oldest_idx);
                retire_chain(ch, &mut blocks);
            }
            active.push(vec![r]);
        }
    }
    for ch in active {
        retire_chain(ch, &mut blocks);
    }
    blocks.sort_by_key(|(grp, _)| grp[0].x_start);
    blocks
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partition_colinear_blocks_preserves_overlapping_chains() {
        let mut pairs = Vec::new();
        for i in 0..15 {
            pairs.push((100 + i, 200 + i));
            pairs.push((100 + i, 500 - i));
        }
        let blocks = partition_colinear_blocks(&pairs);
        assert_eq!(blocks.len(), 2);
        assert!(!blocks[0].1);
        assert!(blocks[1].1);
        assert_eq!(blocks[0].0.len(), 15);
        assert_eq!(blocks[1].0.len(), 15);
    }

    #[test]
    fn net_one_to_one_runs_keeps_longest_and_clips_overlaps() {
        let runs = vec![
            AlignedRun {
                x_start: 10,
                y_start: 20,
                len: 30,
                is_rc: false,
            },
            AlignedRun {
                x_start: 25,
                y_start: 80,
                len: 20,
                is_rc: false,
            },
        ];
        let netted = net_one_to_one_runs(&runs, 100, 120);
        assert_eq!(netted.len(), 1);
        assert_eq!(netted[0].x_start, 10);
        assert_eq!(netted[0].len, 30);
    }

    #[test]
    fn partition_colinear_run_groups_handles_unsorted_and_malformed_rc() {
        let runs = vec![
            AlignedRun {
                x_start: 100,
                y_start: 200,
                len: 15,
                is_rc: false,
            },
            AlignedRun {
                x_start: 10,
                y_start: 2,
                len: 15,
                is_rc: true,
            },
            AlignedRun {
                x_start: 20,
                y_start: 50,
                len: 15,
                is_rc: false,
            },
        ];
        let groups = partition_colinear_run_groups(&runs);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].0[0].x_start, 20);
        assert!(!groups[0].1);
        assert_eq!(groups[1].0[0].x_start, 100);
        assert!(!groups[1].1);
    }

    #[test]
    fn partition_colinear_blocks_rejects_zero_step_duplicates() {
        let mut pairs: Vec<(usize, usize)> = (0..14).map(|i| (i, i)).collect();
        pairs.insert(5, (5, 4)); // duplicate x=5 with backward y=4
        pairs.insert(8, (7, 7)); // duplicate y=7
        let blocks = partition_colinear_blocks(&pairs);
        assert_eq!(blocks.len(), 1);
        for w in blocks[0].0.windows(2) {
            assert!(w[1].0 > w[0].0);
            assert!(w[1].1 > w[0].1);
        }
    }
}
