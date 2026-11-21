//! Run-length encoded alignment types and diagonal merging helpers.

/// Contiguous ungapped colinear run of aligned base pairs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct AlignedRun {
    /// 0-based start coordinate on reference `x` (`+` strand).
    pub x_start: usize,
    /// 0-based start coordinate on query `y` (`+` strand; for `is_rc`, `y` decreases as `x` increases).
    pub y_start: usize,
    /// Number of aligned base pairs in this ungapped run (`>= 1`).
    pub len: usize,
    /// Whether the query `y` steps in reverse (`y_start - k`) as `x` increases (`x_start + k`).
    pub is_rc: bool,
}

impl AlignedRun {
    /// Expand this run into `(x, y)` pairs.
    pub fn push_pairs(&self, out: &mut Vec<(usize, usize)>) {
        out.reserve(self.len);
        if self.is_rc {
            for k in 0..self.len {
                out.push((self.x_start + k, self.y_start - k));
            }
        } else {
            for k in 0..self.len {
                out.push((self.x_start + k, self.y_start + k));
            }
        }
    }
}

/// Compress a slice of `(x, y)` pairs (ordered along a single gap/extension) into [`AlignedRun`]s.
pub fn push_pairs_as_runs(pairs: &[(usize, usize)], is_rc: bool, out: &mut Vec<AlignedRun>) {
    let mut cur: Option<AlignedRun> = None;
    for &(x, y) in pairs {
        if let Some(mut r) = cur {
            let matches_next = if is_rc {
                x == r.x_start + r.len && r.y_start >= r.len && y == r.y_start - r.len
            } else {
                x == r.x_start + r.len && y == r.y_start + r.len
            };
            if matches_next {
                r.len += 1;
                cur = Some(r);
                continue;
            }
            out.push(r);
        }
        cur = Some(AlignedRun {
            x_start: x,
            y_start: y,
            len: 1,
            is_rc,
        });
    }
    if let Some(r) = cur {
        out.push(r);
    }
}

/// Sort and merge overlapping or adjacent [`AlignedRun`]s that lie on the exact same diagonal.
pub fn merge_colinear_runs(runs: &mut Vec<AlignedRun>) {
    if runs.len() <= 1 {
        return;
    }
    let diag_key = |r: &AlignedRun| -> i128 {
        if r.is_rc {
            (r.x_start as i128) + (r.y_start as i128)
        } else {
            (r.y_start as i128) - (r.x_start as i128)
        }
    };
    runs.sort_unstable_by(|a, b| {
        a.is_rc
            .cmp(&b.is_rc)
            .then_with(|| diag_key(a).cmp(&diag_key(b)))
            .then_with(|| a.x_start.cmp(&b.x_start))
            .then_with(|| a.len.cmp(&b.len))
    });
    let mut write_idx = 0usize;
    for read_idx in 1..runs.len() {
        let next = runs[read_idx];
        let cur = &mut runs[write_idx];
        let cur_x_end = cur.x_start + cur.len;
        if cur.is_rc == next.is_rc
            && diag_key(cur) == diag_key(&next)
            && next.x_start <= cur_x_end
        {
            let next_x_end = next.x_start + next.len;
            if next_x_end > cur_x_end {
                cur.len = next_x_end - cur.x_start;
            }
            continue;
        }
        write_idx += 1;
        runs[write_idx] = next;
    }
    runs.truncate(write_idx + 1);
    runs.sort_unstable();
}

/// Expand [`AlignedRun`]s into sorted, deduplicated `(x, y)` pairs.
#[must_use]
pub fn runs_to_pairs(runs: &[AlignedRun]) -> Vec<(usize, usize)> {
    let total: usize = runs.iter().map(|r| r.len).sum();
    let mut pairs = Vec::with_capacity(total);
    for r in runs {
        r.push_pairs(&mut pairs);
    }
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

/// Result of aligning a full contig pair in compact run-encoded form.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunAlignment {
    /// Sorted, diagonal-merged ungapped alignment runs.
    pub aligned_runs: Vec<AlignedRun>,
    /// Total aligned base-pair count (sum of run lengths).
    pub total_pairs: usize,
    /// Number of forward LCBs.
    pub lcbs_fwd: usize,
    /// Number of reverse-complement LCBs.
    pub lcbs_rc: usize,
    /// Number of macro Horizon Bubbles (`> max_gap_bp`).
    pub macro_bubbles: usize,
    /// Number of internal Horizon Bubbles emitted inside DP gaps.
    pub internal_bubbles: usize,
}

/// Result of aligning a full contig pair.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContigAlignment {
    /// Kept 0-based `+` strand `(x, y)` base pairs.
    pub aligned_pairs: Vec<(usize, usize)>,
    /// Number of forward LCBs.
    pub lcbs_fwd: usize,
    /// Number of reverse-complement LCBs.
    pub lcbs_rc: usize,
    /// Number of macro Horizon Bubbles (`> max_gap_bp`).
    pub macro_bubbles: usize,
    /// Number of internal Horizon Bubbles emitted inside DP gaps.
    pub internal_bubbles: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aligned_run_compression_and_merge_roundtrip() {
        let mut runs = Vec::new();
        push_pairs_as_runs(&[(10, 20), (11, 21), (12, 22)], false, &mut runs);
        push_pairs_as_runs(&[(12, 22), (13, 23), (14, 24)], false, &mut runs);
        merge_colinear_runs(&mut runs);
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].x_start, 10);
        assert_eq!(runs[0].y_start, 20);
        assert_eq!(runs[0].len, 5);
        assert_eq!(
            runs_to_pairs(&runs),
            vec![(10, 20), (11, 21), (12, 22), (13, 23), (14, 24)]
        );
    }

    #[test]
    fn merge_colinear_runs_merges_across_intervening_off_diagonal_run() {
        let mut runs = vec![
            AlignedRun {
                x_start: 10,
                y_start: 10,
                len: 10,
                is_rc: false,
            },
            AlignedRun {
                x_start: 15,
                y_start: 200,
                len: 5,
                is_rc: false,
            },
            AlignedRun {
                x_start: 18,
                y_start: 18,
                len: 7,
                is_rc: false,
            },
        ];
        merge_colinear_runs(&mut runs);
        assert_eq!(runs.len(), 2);
        assert_eq!(runs[0].x_start, 10);
        assert_eq!(runs[0].y_start, 10);
        assert_eq!(runs[0].len, 15);
        assert_eq!(runs[1].x_start, 15);
        assert_eq!(runs[1].y_start, 200);
    }
}
