//! Colinear run-group partitioning and 1-to-1 interval netting.

pub(crate) use crate::coord::IntervalSet;
use crate::mdl::AlignedRun;

/// Maximum inter-run gap (bp) merged into a single colinear block and eligible for
/// bounded Gotoh affine-gap micro-realignment (`src/formats.rs`).
pub(crate) const MAX_COLINEAR_GAP_BP: usize = 40;
/// Minimum net sub-run length (bp) retained after 1-to-1 interval subtraction.
const MIN_NETTED_SUBRUN_BP: usize = 6;
/// Minimum colinear run-group length (bp) emitted as a block.
const MIN_COLINEAR_GROUP_BP: usize = 12;

/// Canonical DNA column identity score (`+2` per canonical A/C/G/T match, `-1` per canonical
/// A/C/G/T mismatch, `0` for gaps or ambiguous `N`/non-ACGT bases).
#[inline]
#[must_use]
pub(crate) fn score_acgt_columns(a: &[u8], b: &[u8]) -> i64 {
    let mut score = 0i64;
    for (&ba, &bb) in a.iter().zip(b.iter()) {
        let ux = ba.to_ascii_uppercase();
        let uy = bb.to_ascii_uppercase();
        if matches!(ux, b'A' | b'C' | b'G' | b'T') && matches!(uy, b'A' | b'C' | b'G' | b'T') {
            if ux == uy {
                score += 2;
            } else {
                score -= 1;
            }
        }
    }
    score
}

fn score_run_group_acgt(grp: &[AlignedRun], seq_x: &[u8], seq_y: &[u8]) -> i64 {
    let mut score = 0i64;
    for r in grp {
        let Some(avail_x) = seq_x.len().checked_sub(r.x_start) else {
            continue;
        };
        let len = if r.is_rc {
            if r.y_start >= seq_y.len() {
                continue;
            }
            r.len.min(avail_x).min(r.y_start + 1)
        } else {
            let Some(avail_y) = seq_y.len().checked_sub(r.y_start) else {
                continue;
            };
            r.len.min(avail_x).min(avail_y)
        };
        if len == 0 {
            continue;
        }
        let sx = &seq_x[r.x_start..(r.x_start + len)];
        if r.is_rc {
            let y_lo = r.y_start + 1 - len;
            for (&bx, &by_raw) in sx.iter().zip(seq_y[y_lo..=r.y_start].iter().rev()) {
                let ux = bx.to_ascii_uppercase();
                let uy = crate::seeding::complement_base(by_raw).to_ascii_uppercase();
                if matches!(ux, b'A' | b'C' | b'G' | b'T')
                    && matches!(uy, b'A' | b'C' | b'G' | b'T')
                {
                    score += if ux == uy { 2 } else { -1 };
                }
            }
        } else {
            score += score_acgt_columns(sx, &seq_y[r.y_start..(r.y_start + len)]);
        }
    }
    score
}

#[derive(Default)]
struct ClipScratch {
    free_x: Vec<(usize, usize)>,
    free_y: Vec<(usize, usize)>,
    merge_keys: Vec<usize>,
}

fn clip_run_against_intervals(
    r: AlignedRun,
    (len_x, len_y): (usize, usize),
    used_x: &mut IntervalSet,
    used_y: &mut IntervalSet,
    scratch: &mut ClipScratch,
    kept: &mut Vec<AlignedRun>,
) {
    if r.len == 0 || r.x_start >= len_x {
        return;
    }
    let Some(y_end_rc) = r.y_start.checked_add(1) else {
        return;
    };
    let valid_len = if r.is_rc {
        r.len.min(y_end_rc).min(len_x - r.x_start)
    } else {
        r.len
            .min(len_x - r.x_start)
            .min(len_y.saturating_sub(r.y_start))
    };
    if valid_len < MIN_NETTED_SUBRUN_BP {
        return;
    }
    let x_lo = r.x_start;
    let x_hi = r.x_start + valid_len;
    used_x.free_subintervals(x_lo, x_hi, &mut scratch.free_x);
    if scratch.free_x.is_empty() {
        return;
    }
    let (y_lo, y_hi) = if r.is_rc {
        (y_end_rc - valid_len, y_end_rc.min(len_y))
    } else {
        (r.y_start, (r.y_start + valid_len).min(len_y))
    };
    used_y.free_subintervals(y_lo, y_hi, &mut scratch.free_y);
    if scratch.free_y.is_empty() {
        return;
    }

    let mut ix = 0usize;
    let mut iy = 0usize;
    while ix < scratch.free_x.len() && iy < scratch.free_y.len() {
        let kx0 = scratch.free_x[ix].0 - r.x_start;
        let kx1 = scratch.free_x[ix].1 - r.x_start;
        let (ky0, ky1) = if r.is_rc {
            let rev_idx = scratch.free_y.len() - 1 - iy;
            let (ya, yb) = scratch.free_y[rev_idx];
            (y_end_rc - yb, y_end_rc - ya)
        } else {
            let (ya, yb) = scratch.free_y[iy];
            (ya - r.y_start, yb - r.y_start)
        };
        let k0 = kx0.max(ky0);
        let k1 = kx1.min(ky1);
        if k1 >= k0 + MIN_NETTED_SUBRUN_BP {
            let sub_len = k1 - k0;
            let sx = r.x_start + k0;
            let sy = if r.is_rc {
                r.y_start - k0
            } else {
                r.y_start + k0
            };
            used_x.insert(sx, sx + sub_len, &mut scratch.merge_keys);
            if r.is_rc {
                used_y.insert(sy + 1 - sub_len, sy + 1, &mut scratch.merge_keys);
            } else {
                used_y.insert(sy, sy + sub_len, &mut scratch.merge_keys);
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

/// Net overlapping [`AlignedRun`]s to strict 1-to-1 coverage by prioritizing the highest
/// sequence-identity-scored (`2 * matches - mismatches`) colinear run-groups first, so
/// true orthologous chains beat slightly longer low-identity repeat chains.
#[must_use]
pub fn net_one_to_one_runs_scored(
    runs: &[AlignedRun],
    seq_x: &[u8],
    seq_y: &[u8],
) -> Vec<AlignedRun> {
    net_one_to_one_runs_with_scorer(runs, seq_x.len(), seq_y.len(), |grp| {
        score_run_group_acgt(grp, seq_x, seq_y)
    })
}

struct GlobalGroup {
    pair_idx: usize,
    x_i: usize,
    y_i: usize,
    pair_score: i64,
    group_score: i64,
    group_bp: usize,
    runs: Vec<AlignedRun>,
}

/// Enforce species-wide global 1-to-1 orthology across all aligned `(contig_x, contig_y)` pairs
/// in `O(R log R)` time and `O(R)` memory using disjoint interval sets.
///
/// Ranks colinear run-groups across all contig pairs descending by canonical A/C/G/T identity
/// score (`2 * matches - mismatches`, with tie-breaks by syntenic contig-pair total score and
/// group bp), subtracting already-claimed coordinates in `X` and `Y` and discarding sub-runs
/// shorter than 6 bp. Returns `(unique_x_bp, unique_y_bp)` covered across all contigs.
pub fn net_global_one_to_one_runs_scored<T>(
    aligned_pairs: &mut [(usize, usize, Vec<AlignedRun>)],
    seqs_x: &[T],
    seqs_y: &[T],
) -> (usize, usize)
where
    T: AsRef<[u8]>,
{
    let mut used_x: Vec<IntervalSet> = (0..seqs_x.len()).map(|_| IntervalSet::new()).collect();
    let mut used_y: Vec<IntervalSet> = (0..seqs_y.len()).map(|_| IntervalSet::new()).collect();

    let mut groups: Vec<GlobalGroup> = Vec::new();
    for (pair_idx, (x_i, y_i, runs)) in aligned_pairs.iter().enumerate() {
        if runs.is_empty() {
            continue;
        }
        let Some(sx) = seqs_x.get(*x_i).map(AsRef::as_ref) else {
            continue;
        };
        let Some(sy) = seqs_y.get(*y_i).map(AsRef::as_ref) else {
            continue;
        };
        let mut pair_groups = Vec::new();
        let mut pair_score = 0i64;
        for (grp, _) in partition_colinear_run_groups(runs) {
            let g_bp: usize = grp.iter().map(|r| r.len).sum();
            let g_sc = score_run_group_acgt(&grp, sx, sy);
            pair_score = pair_score.saturating_add(g_sc.max(0));
            pair_groups.push((g_sc, g_bp, grp));
        }
        for (group_score, group_bp, grp) in pair_groups {
            groups.push(GlobalGroup {
                pair_idx,
                x_i: *x_i,
                y_i: *y_i,
                pair_score,
                group_score,
                group_bp,
                runs: grp,
            });
        }
    }

    // Prioritize highest-identity colinear run-group first; when a shorter paralogous group
    // ties or sits inside a syntenic scaffold, break ties by total syntenic contig-pair score
    // and deterministic coordinates.
    groups.sort_unstable_by_key(|g| {
        let r0 = g.runs[0];
        (
            std::cmp::Reverse(g.group_score),
            std::cmp::Reverse(g.pair_score),
            std::cmp::Reverse(g.group_bp),
            g.x_i,
            g.y_i,
            r0.x_start,
            r0.y_start,
            r0.is_rc,
            g.pair_idx,
        )
    });

    for (_, _, runs) in aligned_pairs.iter_mut() {
        runs.clear();
    }

    let mut scratch = ClipScratch::default();
    for g in groups {
        let len_x = seqs_x[g.x_i].as_ref().len();
        let len_y = seqs_y[g.y_i].as_ref().len();
        let out_runs = &mut aligned_pairs[g.pair_idx].2;
        for r in g.runs {
            clip_run_against_intervals(
                r,
                (len_x, len_y),
                &mut used_x[g.x_i],
                &mut used_y[g.y_i],
                &mut scratch,
                out_runs,
            );
        }
    }

    for (_, _, runs) in aligned_pairs.iter_mut() {
        runs.sort_unstable_by_key(|r| (r.x_start, r.y_start, r.is_rc));
    }

    let unique_x: usize = used_x
        .iter()
        .flat_map(IntervalSet::iter)
        .map(|(s, e)| e - s)
        .sum();
    let unique_y: usize = used_y
        .iter()
        .flat_map(IntervalSet::iter)
        .map(|(s, e)| e - s)
        .sum();
    (unique_x, unique_y)
}

/// Compute exact unique covered bases across all contigs in `O(R log R)` time and `O(R)` memory.
#[must_use]
pub fn count_unique_covered_bp(
    aligned_pairs: &[(usize, usize, &[AlignedRun])],
    n_contigs_x: usize,
    n_contigs_y: usize,
) -> (usize, usize) {
    let mut cov_x: Vec<IntervalSet> = (0..n_contigs_x).map(|_| IntervalSet::new()).collect();
    let mut cov_y: Vec<IntervalSet> = (0..n_contigs_y).map(|_| IntervalSet::new()).collect();
    let mut scratch = Vec::new();
    for &(x_i, y_i, runs) in aligned_pairs {
        for r in runs {
            if r.len == 0 {
                continue;
            }
            if let Some(sx) = cov_x.get_mut(x_i) {
                sx.insert(r.x_start, r.x_start.saturating_add(r.len), &mut scratch);
            }
            if let Some(sy) = cov_y.get_mut(y_i) {
                if r.is_rc {
                    let y_hi = r.y_start.saturating_add(1);
                    let y_lo = y_hi.saturating_sub(r.len);
                    sy.insert(y_lo, y_hi, &mut scratch);
                } else {
                    sy.insert(r.y_start, r.y_start.saturating_add(r.len), &mut scratch);
                }
            }
        }
    }
    let ux: usize = cov_x
        .iter()
        .flat_map(IntervalSet::iter)
        .map(|(s, e)| e - s)
        .sum();
    let uy: usize = cov_y
        .iter()
        .flat_map(IntervalSet::iter)
        .map(|(s, e)| e - s)
        .sum();
    (ux, uy)
}

fn net_one_to_one_runs_with_scorer<F>(
    runs: &[AlignedRun],
    len_x: usize,
    len_y: usize,
    mut score_grp: F,
) -> Vec<AlignedRun>
where
    F: FnMut(&[AlignedRun]) -> i64,
{
    if runs.is_empty() || len_x == 0 || len_y == 0 {
        return Vec::new();
    }
    let mut blocks: Vec<(i64, usize, Vec<AlignedRun>)> = partition_colinear_run_groups(runs)
        .into_iter()
        .map(|(grp, _)| {
            let total_bp: usize = grp.iter().map(|r| r.len).sum();
            let sc = score_grp(&grp);
            (sc, total_bp, grp)
        })
        .collect();
    blocks.sort_unstable_by_key(|(sc, bp, grp)| {
        (
            std::cmp::Reverse(*sc),
            std::cmp::Reverse(*bp),
            grp[0].x_start,
            grp[0].y_start,
            grp[0].is_rc,
        )
    });

    let mut used_x = IntervalSet::new();
    let mut used_y = IntervalSet::new();
    let mut scratch = ClipScratch::default();
    let mut kept: Vec<AlignedRun> = Vec::with_capacity(runs.len());

    for (_, _, grp) in blocks {
        for r in grp {
            clip_run_against_intervals(
                r,
                (len_x, len_y),
                &mut used_x,
                &mut used_y,
                &mut scratch,
                &mut kept,
            );
        }
    }
    kept.sort_unstable_by_key(|r| (r.x_start, r.y_start, r.is_rc));
    kept
}

/// Partition [`AlignedRun`]s into colinear block groups of at least 12 bases.
pub(crate) fn partition_colinear_run_groups(runs: &[AlignedRun]) -> Vec<(Vec<AlignedRun>, bool)> {
    partition_colinear_run_groups_avoiding_gaps(runs, None, None)
}

/// Partition [`AlignedRun`]s into colinear block groups of at least 12 bases, refusing to bridge
/// any intervening gap on `x` or `y` that overlaps `forbidden_x` or `forbidden_y`.
///
/// Used by `--one-to-one` MAF/PAF emission so a `<= 40 bp` gap claimed by another contig pair
/// (or another diagonal) is never bridged into an overlapping target or query block span.
pub(crate) fn partition_colinear_run_groups_avoiding_gaps(
    runs: &[AlignedRun],
    forbidden_x: Option<&IntervalSet>,
    forbidden_y: Option<&IntervalSet>,
) -> Vec<(Vec<AlignedRun>, bool)> {
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
        if total_bp >= MIN_COLINEAR_GROUP_BP {
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
            if last
                .x_start
                .saturating_add(last.len)
                .saturating_add(MAX_COLINEAR_GAP_BP)
                < r.x_start
            {
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
            let prev_x_end = prev.x_start.saturating_add(prev.len);
            if r.x_start < prev_x_end || r.x_start - prev_x_end > MAX_COLINEAR_GAP_BP {
                continue;
            }
            if forbidden_x.is_some_and(|fx| fx.overlaps(prev_x_end, r.x_start)) {
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
                r.y_start < prev_y_lim
                    && (prev_y_lim - 1 - r.y_start) <= MAX_COLINEAR_GAP_BP
                    && !forbidden_y.is_some_and(|fy| fy.overlaps(r.y_start + 1, prev_y_lim))
            } else {
                let prev_y_lim = prev.y_start.saturating_add(prev.len);
                r.y_start >= prev_y_lim
                    && (r.y_start - prev_y_lim) <= MAX_COLINEAR_GAP_BP
                    && !forbidden_y.is_some_and(|fy| fy.overlaps(prev_y_lim, r.y_start))
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
                    last.x_start.saturating_add(last.len)
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
    blocks.sort_by_key(|(grp, _)| (grp[0].x_start, grp[0].y_start, grp[0].is_rc));
    blocks
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partition_colinear_run_groups_preserves_overlapping_chains() {
        let runs = vec![
            AlignedRun {
                x_start: 100,
                y_start: 200,
                len: 15,
                is_rc: false,
            },
            AlignedRun {
                x_start: 100,
                y_start: 500,
                len: 15,
                is_rc: true,
            },
        ];
        let groups = partition_colinear_run_groups(&runs);
        assert_eq!(groups.len(), 2);
        assert!(!groups[0].1);
        assert!(groups[1].1);
        assert_eq!(groups[0].0[0].len, 15);
        assert_eq!(groups[1].0[0].len, 15);
    }

    #[test]
    fn net_one_to_one_runs_scored_keeps_longest_and_clips_overlaps() {
        let seq_x = vec![b'A'; 100];
        let seq_y = vec![b'A'; 120];
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
        let netted = net_one_to_one_runs_scored(&runs, &seq_x, &seq_y);
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
    fn partition_colinear_run_groups_keeps_strictly_increasing_x_per_group() {
        let runs = vec![
            AlignedRun {
                x_start: 0,
                y_start: 0,
                len: 8,
                is_rc: false,
            },
            AlignedRun {
                x_start: 5,
                y_start: 4,
                len: 8,
                is_rc: false,
            },
            AlignedRun {
                x_start: 10,
                y_start: 10,
                len: 8,
                is_rc: false,
            },
        ];
        let groups = partition_colinear_run_groups(&runs);
        assert_eq!(groups.len(), 1);
        for w in groups[0].0.windows(2) {
            assert!(w[1].x_start >= w[0].x_start + w[0].len);
            assert!(w[1].y_start >= w[0].y_start + w[0].len);
        }
    }

    #[test]
    fn net_one_to_one_runs_scored_prefers_exact_ortholog_over_longer_divergent_repeat() {
        // seq_x[0..60] matches seq_y[0..50] at 100% identity (score = 100),
        // and overlaps a 60 bp alignment to seq_y[100..160] at 50% identity (score = 30).
        let seq_x = b"ACGT".repeat(15); // 60 bp
        let mut seq_y = Vec::with_capacity(200);
        seq_y.extend_from_slice(&seq_x[..50]); // 0..50: 100% match
        seq_y.resize(100, b'N');
        for (i, &bx) in seq_x.iter().enumerate().take(60) {
            // 100..160: 50% identity
            seq_y.push(if i % 2 == 0 { bx } else { b'T' });
        }
        let runs = vec![
            AlignedRun {
                x_start: 0,
                y_start: 100,
                len: 60,
                is_rc: false,
            },
            AlignedRun {
                x_start: 0,
                y_start: 0,
                len: 50,
                is_rc: false,
            },
        ];
        let netted = net_one_to_one_runs_scored(&runs, &seq_x, &seq_y);
        assert!(
            netted.iter().any(|r| r.y_start == 0 && r.len == 50),
            "expected 50 bp exact ortholog to win over 60 bp divergent repeat: {netted:?}"
        );
    }

    #[test]
    fn interval_set_supports_insert_subtract_queries_and_scale() {
        let mut set = IntervalSet::from_single(0, 100);
        let mut scratch = Vec::new();
        assert!(set.contains_range(10, 90));
        set.subtract(20, 40, &mut scratch);
        set.subtract(60, 80, &mut scratch);
        assert_eq!(
            set.iter().collect::<Vec<_>>(),
            vec![(0, 20), (40, 60), (80, 100)]
        );
        assert!(set.contains_pos(19));
        assert!(!set.contains_pos(20));
        assert!(!set.contains_range(15, 25));
        assert!(set.contains_range(40, 60));

        // Insert back overlapping intervals and coalesce.
        set.insert(15, 45, &mut scratch);
        assert_eq!(set.iter().collect::<Vec<_>>(), vec![(0, 60), (80, 100)]);

        let mut free = Vec::new();
        set.free_subintervals(50, 90, &mut free);
        assert_eq!(free, vec![(60, 80)]);

        // Scale check: 10,000 non-monotonic insertions + subtractions stay O(M log M).
        let mut big = IntervalSet::new();
        for i in (0..10_000usize).rev() {
            big.insert(i * 10, i * 10 + 5, &mut scratch);
        }
        assert_eq!(big.iter().count(), 10_000);
        for i in 0..10_000usize {
            big.subtract(i * 10 + 1, i * 10 + 4, &mut scratch);
        }
        assert_eq!(big.iter().count(), 20_000);
    }

    #[test]
    fn count_unique_covered_bp_handles_fwd_rc_overlaps_and_oob() {
        let r_fwd = [
            AlignedRun {
                x_start: 0,
                y_start: 10,
                len: 20,
                is_rc: false,
            },
            AlignedRun {
                x_start: 10,
                y_start: 20,
                len: 20,
                is_rc: false,
            },
        ];
        // RC run covering y=[5..15) (y_start=14, len=10) overlapping fwd y=[10..40) on y[10..15)
        let r_rc = [AlignedRun {
            x_start: 5,
            y_start: 14,
            len: 10,
            is_rc: true,
        }];
        let pairs: Vec<(usize, usize, &[AlignedRun])> = vec![
            (0, 0, &r_fwd),
            (1, 0, &r_rc),
            (99, 99, &r_fwd), // OOB contig indices ignored cleanly
        ];
        let (ux, uy) = count_unique_covered_bp(&pairs, 2, 1);
        // X: contig 0 has [0..30) = 30, contig 1 has [5..15) = 10 -> 40
        // Y: contig 0 has [5..40) = 35
        assert_eq!((ux, uy), (40, 35));
    }

    #[test]
    fn net_global_one_to_one_handles_rc_clipping_y_overlap_tiebreak_and_oob() {
        // Two X contigs compete for the same Y contig:
        // Pair 0 (x=0, y=0) has an extra syntenic anchor so its pair_score wins a group_score tie.
        // Pair 1 (x=1, y=0) has an RC run partially overlapping y=[20..40) at y=[10..30).
        let x0 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 40 bp
        let x1 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 40 bp (self-complementary repeats!)
        let y0 = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 40 bp
        let mut pairs = vec![
            (
                0usize,
                0usize,
                vec![
                    AlignedRun {
                        x_start: 0,
                        y_start: 0,
                        len: 10,
                        is_rc: false,
                    },
                    AlignedRun {
                        x_start: 20,
                        y_start: 20,
                        len: 20,
                        is_rc: false,
                    },
                ],
            ),
            (
                1usize,
                0usize,
                vec![AlignedRun {
                    x_start: 0,
                    y_start: 29,
                    len: 20,
                    is_rc: true, // covers y=[10..30), overlaps [20..30) with pair 0
                }],
            ),
            (
                9usize,
                0usize,
                vec![AlignedRun {
                    x_start: 0,
                    y_start: 0,
                    len: 10,
                    is_rc: false,
                }],
            ),
        ];
        let (ux, uy) = net_global_one_to_one_runs_scored(&mut pairs, &[x0, x1], &[y0]);
        // Pair 0 keeps both runs ([0..10) + [20..40) = 30 bp).
        assert_eq!(pairs[0].2.len(), 2);
        // Pair 1's RC run y=[10..30) is clipped to y=[10..20) (len 10, x=[10..20), y_start=19).
        assert_eq!(
            pairs[1].2,
            vec![AlignedRun {
                x_start: 10,
                y_start: 19,
                len: 10,
                is_rc: true,
            }]
        );
        // OOB pair 2 is cleared cleanly.
        assert_eq!(pairs[2].2, [] as [AlignedRun; 0]);
        assert_eq!((ux, uy), (40, 40));
    }
}
