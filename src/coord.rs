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
            strand == '+' || strand == '-',
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
            strand == '+' || strand == '-',
            "Interval strand must be '+' or '-', got {strand:?}"
        );
        Self {
            start,
            end,
            strand,
            src_size,
        }
    }

    /// The span length of the interval.
    #[must_use]
    pub fn span(&self) -> usize {
        self.end - self.start
    }

    /// The 0-based start coordinate relative to the recorded strand (UCSC/MAF format).
    #[must_use]
    pub fn ucsc_start(&self) -> usize {
        if self.strand == '-' {
            self.src_size - self.end
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

    #[test]
    fn try_new_rejects_invalid_bounds_and_strand() {
        assert!(Interval::try_new(40, 35, '+', 100).is_err());
        assert!(Interval::try_new(10, 105, '+', 100).is_err());
        assert!(Interval::try_new(10, 35, '?', 100).is_err());
        assert!(Interval::try_new(10, 35, '-', 100).is_ok());
    }
}
