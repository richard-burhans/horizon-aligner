//! Genomic strand, interval coordinate arithmetic, and half-open `IntervalSet`.

use std::collections::BTreeMap;

/// Orientation of a genomic alignment or graph traversal step.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    /// Forward (`'+'`) strand.
    Plus,
    /// Reverse-complement (`'-'`) strand.
    Minus,
}

impl Strand {
    /// Parse a MAF/BED/GFA strand character (`'+'` or `'-'`).
    #[inline]
    #[must_use]
    pub const fn from_char(ch: char) -> Option<Self> {
        match ch {
            '+' => Some(Self::Plus),
            '-' => Some(Self::Minus),
            _ => None,
        }
    }

    /// Parse a GFA v1.1 walk orientation character (`'>'` -> `Plus`, `'<'` -> `Minus`).
    #[inline]
    #[must_use]
    pub const fn from_gfa_walk_orient(ch: char) -> Option<Self> {
        match ch {
            '>' => Some(Self::Plus),
            '<' => Some(Self::Minus),
            _ => None,
        }
    }

    /// Convert to canonical `'+'` or `'-'` character.
    #[inline]
    #[must_use]
    pub const fn as_char(self) -> char {
        match self {
            Self::Plus => '+',
            Self::Minus => '-',
        }
    }

    /// Flip strand (`Plus <-> Minus`).
    #[inline]
    #[must_use]
    pub const fn flip(self) -> Self {
        match self {
            Self::Plus => Self::Minus,
            Self::Minus => Self::Plus,
        }
    }

    /// Returns `true` for [`Strand::Plus`].
    #[inline]
    #[must_use]
    pub const fn is_plus(self) -> bool {
        matches!(self, Self::Plus)
    }

    /// Returns `true` for [`Strand::Minus`].
    #[inline]
    #[must_use]
    pub const fn is_minus(self) -> bool {
        matches!(self, Self::Minus)
    }
}

/// A genomic interval enforcing strand-aware UCSC coordinate arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval {
    /// 0-based start on the forward strand.
    pub start: usize,
    /// 0-based end (exclusive) on the forward strand.
    pub end: usize,
    /// Strand (`+` or `-`).
    pub strand: char,
    /// Length of the parent sequence (required for `-` strand UCSC translation).
    pub src_size: usize,
}

impl Interval {
    /// Construct and validate a new interval (bounds on forward strand).
    ///
    /// # Errors
    /// Returns an error if `start > end`, `end > src_size`, or `strand` is not `'+'` or `'-'`.
    pub fn try_new(
        start: usize,
        end: usize,
        strand: char,
        src_size: usize,
    ) -> anyhow::Result<Self> {
        anyhow::ensure!(start <= end, "Interval start ({start}) > end ({end})");
        anyhow::ensure!(
            end <= src_size,
            "Interval end ({end}) > src_size ({src_size})"
        );
        anyhow::ensure!(
            Strand::from_char(strand).is_some(),
            "Interval strand must be '+' or '-', got {strand:?}"
        );
        Ok(Self {
            start,
            end,
            strand,
            src_size,
        })
    }

    /// Construct a new interval (bounds on forward strand).
    ///
    /// # Panics
    /// Panics if `start > end`, `end > src_size`, or `strand` is not `'+'` or `'-'`.
    #[must_use]
    pub fn new(start: usize, end: usize, strand: char, src_size: usize) -> Self {
        assert!(start <= end, "Interval start ({start}) > end ({end})");
        assert!(
            end <= src_size,
            "Interval end ({end}) > src_size ({src_size})"
        );
        assert!(
            Strand::from_char(strand).is_some(),
            "Interval strand must be '+' or '-', got {strand:?}"
        );
        Self {
            start,
            end,
            strand,
            src_size,
        }
    }

    /// Construct and validate an interval from UCSC/MAF strand-relative `(ucsc_start, size, strand, src_size)`.
    ///
    /// # Errors
    /// Returns an error if `strand` is not `'+'`/`'-'` or `ucsc_start + size > src_size`.
    pub fn from_ucsc(
        ucsc_start: usize,
        size: usize,
        strand: char,
        src_size: usize,
    ) -> anyhow::Result<Self> {
        let st = Strand::from_char(strand)
            .ok_or_else(|| anyhow::anyhow!("Interval strand must be '+' or '-', got {strand:?}"))?;
        let end_rel = ucsc_start
            .checked_add(size)
            .ok_or_else(|| anyhow::anyhow!("Interval UCSC end overflow"))?;
        anyhow::ensure!(
            end_rel <= src_size,
            "Interval UCSC end ({end_rel}) > src_size ({src_size})"
        );
        let (start, end) = match st {
            Strand::Plus => (ucsc_start, end_rel),
            Strand::Minus => (src_size - end_rel, src_size - ucsc_start),
        };
        Ok(Self {
            start,
            end,
            strand,
            src_size,
        })
    }

    /// Convert a UCSC/MAF strand-relative span to a forward-strand `[fwd_start, fwd_end)` interval
    /// using saturating arithmetic (never panics on malformed external inputs).
    #[inline]
    #[must_use]
    pub const fn ucsc_to_fwd_interval(
        ucsc_start: usize,
        size: usize,
        strand: Strand,
        src_size: usize,
    ) -> (usize, usize) {
        match strand {
            Strand::Plus => (ucsc_start, ucsc_start.saturating_add(size)),
            Strand::Minus => {
                let fwd_end = src_size.saturating_sub(ucsc_start);
                let fwd_start = fwd_end.saturating_sub(size);
                (fwd_start, fwd_end)
            }
        }
    }

    /// Convert a 0-based ungapped offset `k` along a UCSC/MAF row (`ucsc_start`, `strand`, `src_size`)
    /// to its 0-based forward-strand coordinate using saturating arithmetic.
    #[inline]
    #[must_use]
    pub const fn ucsc_offset_to_fwd_pos(
        ucsc_start: usize,
        k: usize,
        strand: Strand,
        src_size: usize,
    ) -> usize {
        match strand {
            Strand::Plus => ucsc_start.saturating_add(k),
            Strand::Minus => src_size
                .saturating_sub(1)
                .saturating_sub(ucsc_start.saturating_add(k)),
        }
    }

    /// Convert a UCSC/MAF `(ucsc_start, size, strand, src_size)` to the opposite strand's UCSC start.
    #[inline]
    #[must_use]
    pub const fn flip_ucsc_start(ucsc_start: usize, size: usize, src_size: usize) -> usize {
        src_size.saturating_sub(ucsc_start.saturating_add(size))
    }

    /// Typed [`Strand`] for this interval (`Plus` for `'+'`, `Minus` for `'-'`).
    #[inline]
    #[must_use]
    pub const fn strand_enum(&self) -> Strand {
        if self.strand == '-' {
            Strand::Minus
        } else {
            Strand::Plus
        }
    }

    /// The span length of the interval.
    #[inline]
    #[must_use]
    pub const fn span(&self) -> usize {
        self.end - self.start
    }

    /// The 0-based start coordinate relative to the recorded strand (UCSC/MAF format).
    #[inline]
    #[must_use]
    pub const fn ucsc_start(&self) -> usize {
        if self.strand == '-' {
            self.src_size - self.end
        } else {
            self.start
        }
    }

    /// Forward-strand 0-based coordinate of the `k`-th base along this interval's strand (`0 <= k < span`).
    #[inline]
    #[must_use]
    pub const fn fwd_pos(&self, k: usize) -> usize {
        Self::ucsc_offset_to_fwd_pos(self.ucsc_start(), k, self.strand_enum(), self.src_size)
    }

    /// Return the same genomic span on the opposite strand.
    #[inline]
    #[must_use]
    pub const fn revcomp(&self) -> Self {
        Self {
            start: self.start,
            end: self.end,
            strand: self.strand_enum().flip().as_char(),
            src_size: self.src_size,
        }
    }
}

/// Disjoint half-open interval set `[start, end)` backed by a `BTreeMap` for `O(log M)`
/// membership, insertion, subtraction, overlap-length, and free-gap queries.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct IntervalSet {
    map: BTreeMap<usize, usize>,
}

impl IntervalSet {
    /// Construct an empty half-open interval set.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Construct an interval set containing the single half-open range `[lo, hi)` (if `lo < hi`).
    #[must_use]
    pub fn from_single(lo: usize, hi: usize) -> Self {
        let mut set = Self::new();
        if lo < hi {
            set.map.insert(lo, hi);
        }
        set
    }

    /// Number of disjoint coalesced intervals currently stored.
    #[must_use]
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Returns `true` when no intervals are stored.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Returns `true` if position `pos` lies inside any stored `[s, e)`.
    #[must_use]
    pub fn contains_pos(&self, pos: usize) -> bool {
        self.map
            .range(..=pos)
            .next_back()
            .is_some_and(|(_, &e)| pos < e)
    }

    /// Returns `true` if any stored interval overlaps the half-open range `[lo, hi)`.
    #[must_use]
    pub fn overlaps(&self, lo: usize, hi: usize) -> bool {
        if lo >= hi {
            return false;
        }
        if self
            .map
            .range(..=lo)
            .next_back()
            .is_some_and(|(_, &e)| e > lo)
        {
            return true;
        }
        self.map.range(lo..hi).next().is_some()
    }

    /// Returns `true` if the entire half-open range `[lo, hi)` is covered by a single stored interval.
    #[must_use]
    pub fn contains_range(&self, lo: usize, hi: usize) -> bool {
        if lo >= hi {
            return false;
        }
        self.map
            .range(..=lo)
            .next_back()
            .is_some_and(|(_, &e)| hi <= e)
    }

    /// Populate `out` with all uncovered sub-intervals of `[lo, hi)` in ascending order.
    pub fn free_subintervals(&self, lo: usize, hi: usize, out: &mut Vec<(usize, usize)>) {
        out.clear();
        if lo >= hi {
            return;
        }
        let mut cur = lo;
        if let Some((&_s, &e)) = self.map.range(..lo).next_back()
            && e > cur
        {
            cur = e;
        }
        for (&s, &e) in self.map.range(lo..hi) {
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

    /// Insert `[lo, hi)` and coalesce any overlapping or adjacent intervals, reusing `scratch`.
    pub fn insert(&mut self, mut lo: usize, mut hi: usize, scratch: &mut Vec<usize>) {
        if lo >= hi {
            return;
        }
        scratch.clear();
        if let Some((&ps, &pe)) = self.map.range(..=lo).next_back()
            && pe >= lo
        {
            lo = lo.min(ps);
            hi = hi.max(pe);
            scratch.push(ps);
        }
        for (&s, &e) in self.map.range(lo..=hi) {
            lo = lo.min(s);
            hi = hi.max(e);
            scratch.push(s);
        }
        for &k in scratch.iter() {
            self.map.remove(&k);
        }
        self.map.insert(lo, hi);
    }

    /// Convenience wrapper around [`Self::insert`] without an external scratch buffer.
    pub fn insert_simple(&mut self, lo: usize, hi: usize) {
        let mut scratch = Vec::new();
        self.insert(lo, hi, &mut scratch);
    }

    /// Insert `[lo, hi)` only if it does not overlap any existing interval (returning `true` if inserted).
    pub fn insert_if_disjoint(&mut self, lo: usize, hi: usize) -> bool {
        if lo >= hi {
            return false;
        }
        if let Some((_, &e)) = self.map.range(..hi).next_back()
            && e > lo
        {
            return false;
        }
        self.map.insert(lo, hi);
        true
    }

    /// Subtract `[cut_s, cut_e)` from the set, splitting any spanning interval in place.
    pub fn subtract(&mut self, cut_s: usize, cut_e: usize, scratch: &mut Vec<usize>) {
        if cut_s >= cut_e {
            return;
        }
        scratch.clear();
        let mut reinsert_left: Option<(usize, usize)> = None;
        let mut reinsert_right: Option<(usize, usize)> = None;
        if let Some((&ps, &pe)) = self.map.range(..cut_s).next_back()
            && pe > cut_s
        {
            scratch.push(ps);
            reinsert_left = Some((ps, cut_s));
            if pe > cut_e {
                reinsert_right = Some((cut_e, pe));
            }
        }
        for (&s, &e) in self.map.range(cut_s..cut_e) {
            scratch.push(s);
            if e > cut_e {
                reinsert_right = Some((cut_e, e));
            }
        }
        for &k in scratch.iter() {
            self.map.remove(&k);
        }
        if let Some((s, e)) = reinsert_left {
            self.map.insert(s, e);
        }
        if let Some((s, e)) = reinsert_right {
            self.map.insert(s, e);
        }
    }

    /// Total bases in `[lo, hi)` already covered by this set (`O(log M + k)`).
    #[must_use]
    pub fn overlap_bp(&self, lo: usize, hi: usize) -> usize {
        if lo >= hi {
            return 0;
        }
        let mut ov = 0usize;
        if let Some((&a, &b)) = self.map.range(..lo).next_back()
            && b > lo
        {
            ov += b.min(hi) - a.max(lo);
        }
        for (&a, &b) in self.map.range(lo..hi) {
            ov += b.min(hi) - a;
        }
        ov
    }

    /// Total bases covered across all disjoint intervals in this set.
    #[must_use]
    pub fn covered_bp(&self) -> usize {
        self.map.iter().map(|(&s, &e)| e.saturating_sub(s)).sum()
    }

    /// Iterate over all disjoint half-open intervals `(start, end)` in ascending order.
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.map.iter().map(|(&s, &e)| (s, e))
    }
}

#[cfg(test)]
mod tests {
    use super::{Interval, IntervalSet, Strand};

    #[test]
    fn forward_and_reverse_ucsc_coords() {
        let fwd = Interval::new(10, 35, '+', 100);
        assert_eq!(fwd.span(), 25);
        assert_eq!(fwd.ucsc_start(), 10);
        assert_eq!(fwd.fwd_pos(0), 10);
        assert_eq!(fwd.fwd_pos(24), 34);
        assert_eq!(fwd.strand_enum(), Strand::Plus);

        let rev = Interval::new(10, 35, '-', 100);
        assert_eq!(rev.span(), 25);
        assert_eq!(rev.ucsc_start(), 65);
        assert_eq!(rev.fwd_pos(0), 34);
        assert_eq!(rev.fwd_pos(24), 10);
        assert_eq!(rev.strand_enum(), Strand::Minus);
        assert_eq!(rev.revcomp(), fwd);
        assert_eq!(Interval::from_ucsc(65, 25, '-', 100).expect("valid"), rev);
    }

    #[test]
    fn try_new_and_from_ucsc_reject_invalid_bounds_and_strand() {
        assert!(Interval::try_new(40, 35, '+', 100).is_err());
        assert!(Interval::try_new(10, 105, '+', 100).is_err());
        assert!(Interval::try_new(10, 35, '?', 100).is_err());
        assert!(Interval::try_new(10, 35, '-', 100).is_ok());
        assert!(Interval::from_ucsc(80, 30, '-', 100).is_err());
        assert!(Interval::from_ucsc(10, 20, '?', 100).is_err());
    }

    #[test]
    fn interval_set_coalesces_subtracts_disjoint_and_overlaps() {
        let mut iv = IntervalSet::new();
        iv.insert_simple(10, 20);
        iv.insert_simple(20, 30);
        iv.insert_simple(12, 18);
        iv.insert_simple(40, 50);
        iv.insert_simple(25, 45);
        assert_eq!(iv.len(), 1);
        assert_eq!(iv.covered_bp(), 40);
        assert_eq!(iv.overlap_bp(0, 10), 0);
        assert_eq!(iv.overlap_bp(5, 15), 5);
        assert_eq!(iv.overlap_bp(45, 60), 5);
        assert_eq!(iv.overlap_bp(0, 100), 40);
        assert!(!iv.insert_if_disjoint(49, 60));
        assert!(iv.insert_if_disjoint(50, 60));
        assert_eq!(iv.covered_bp(), 50);
    }
}
