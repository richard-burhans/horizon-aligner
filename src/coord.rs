//! Genomic interval coordinate arithmetic.

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
    /// Construct a new interval (bounds on forward strand).
    #[must_use]
    pub fn new(start: usize, end: usize, strand: char, src_size: usize) -> Self {
        Self { start, end, strand, src_size }
    }

    /// The span length of the interval.
    #[must_use]
    pub fn span(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    /// The 0-based start coordinate relative to the recorded strand (UCSC/MAF format).
    #[must_use]
    pub fn ucsc_start(&self) -> usize {
        if self.strand == '-' {
            self.src_size.saturating_sub(self.end)
        } else {
            self.start
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Interval;

    #[test]
    fn forward_and_reverse_ucsc_coords() {
        let fwd = Interval::new(10, 35, '+', 100);
        assert_eq!(fwd.span(), 25);
        assert_eq!(fwd.ucsc_start(), 10);

        let rev = Interval::new(10, 35, '-', 100);
        assert_eq!(rev.span(), 25);
        assert_eq!(rev.ucsc_start(), 65);
    }
}
