//! Core TBA data structures (`TbaRow`, `TbaBlock`, `PairwiseTbaBlock`) and MAF parsing.

use crate::coord::{Interval, Strand};
use crate::dna::complement_base;
use crate::formats::{
    MafRow, for_each_maf_stanza_reader, matches_sample_prefix, strip_sample_prefix,
};
use std::borrow::Cow;
use std::collections::HashSet;

/// One species row inside a reference-free [`TbaBlock`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TbaRow {
    /// Species name.
    pub species: String,
    /// Contig/sequence name (without duplicated `<species>.` prefix).
    pub contig: String,
    /// 0-based start coordinate on `strand` (UCSC MAF convention).
    pub start: usize,
    /// Number of ungapped bases in `seq`.
    pub size: usize,
    /// Alignment strand (`'+'` or `'-'`).
    pub strand: char,
    /// Full contig length in bases.
    pub src_size: usize,
    /// ASCII alignment columns (`A/C/G/T/N/-`).
    pub seq: Vec<u8>,
}

impl TbaRow {
    /// Forward-strand half-open interval `[fwd_start, fwd_end)` covered by this row.
    #[must_use]
    pub fn fwd_interval(&self) -> (usize, usize) {
        let strand = Strand::from_char(self.strand).unwrap_or(Strand::Plus);
        Interval::ucsc_to_fwd_interval(self.start, self.size, strand, self.src_size)
    }

    /// Forward-strand 0-based genomic coordinate of the `k`-th ungapped base (`0 <= k < self.size`).
    #[must_use]
    pub fn fwd_pos(&self, k: usize) -> usize {
        let strand = Strand::from_char(self.strand).unwrap_or(Strand::Plus);
        Interval::ucsc_offset_to_fwd_pos(self.start, k, strand, self.src_size)
    }

    /// Reverse-complement this row while preserving exact UCSC MAF coordinates.
    #[must_use]
    pub fn revcomp(&self) -> Self {
        let strand = Strand::from_char(self.strand).unwrap_or(Strand::Plus);
        let new_strand = strand.flip().as_char();
        let new_start = Interval::flip_ucsc_start(self.start, self.size, self.src_size);
        let new_seq: Vec<u8> = self
            .seq
            .iter()
            .rev()
            .map(|&b| if b == b'-' { b'-' } else { complement_base(b) })
            .collect();
        Self {
            species: self.species.clone(),
            contig: self.contig.clone(),
            start: new_start,
            size: self.size,
            strand: new_strand,
            src_size: self.src_size,
            seq: new_seq,
        }
    }
}

/// One reference-free multi-species alignment block.
#[derive(Debug, Clone, PartialEq)]
pub struct TbaBlock {
    /// Alignment score (higher wins when resolving conflicting cross-clade links).
    pub score: f64,
    /// Rows with identical column count and distinct species names.
    pub rows: Vec<TbaRow>,
}

impl TbaBlock {
    /// Number of alignment columns in this block.
    #[must_use]
    pub fn n_cols(&self) -> usize {
        self.rows.first().map_or(0, |r| r.seq.len())
    }

    /// Reverse-complement every row in this block.
    #[must_use]
    pub fn revcomp(&self) -> Self {
        Self {
            score: self.score,
            rows: self.rows.iter().map(TbaRow::revcomp).collect(),
        }
    }

    /// Slice columns `[lo_col, hi_col)`, stripping all-gap columns and dropping rows with 0 bases.
    #[must_use]
    pub fn slice_cols(&self, lo_col: usize, hi_col: usize) -> Option<Self> {
        let n_cols = self.n_cols();
        let hi = hi_col.min(n_cols);
        if lo_col >= hi {
            return None;
        }
        let mut keep_col = vec![false; hi - lo_col];
        let mut row_slices: Vec<(usize, usize, usize)> = Vec::with_capacity(self.rows.len());
        for (ri, r) in self.rows.iter().enumerate() {
            if r.seq.len() < hi {
                return None;
            }
            let skipped = r.seq[..lo_col].iter().filter(|&&b| b != b'-').count();
            let mut kept_bases = 0usize;
            for (k, &b) in r.seq[lo_col..hi].iter().enumerate() {
                if b != b'-' {
                    kept_bases += 1;
                    keep_col[k] = true;
                }
            }
            if kept_bases > 0 {
                row_slices.push((ri, skipped, kept_bases));
            }
        }
        let out_cols = keep_col.iter().filter(|&&k| k).count();
        if out_cols == 0 || row_slices.is_empty() {
            return None;
        }
        let mut rows = Vec::with_capacity(row_slices.len());
        for (ri, skipped, kept_bases) in row_slices {
            let r = &self.rows[ri];
            let mut seq = Vec::with_capacity(out_cols);
            for (&b, &keep) in r.seq[lo_col..hi].iter().zip(keep_col.iter()) {
                if keep {
                    seq.push(b);
                }
            }
            rows.push(TbaRow {
                species: r.species.clone(),
                contig: r.contig.clone(),
                start: r.start.saturating_add(skipped),
                size: kept_bases,
                strand: r.strand,
                src_size: r.src_size,
                seq,
            });
        }
        #[allow(clippy::cast_precision_loss)]
        let frac = (out_cols as f64) / (n_cols.max(1) as f64);
        Some(Self {
            score: self.score * frac,
            rows,
        })
    }
}

/// One two-row pairwise alignment block between species `a` and species `b`.
#[derive(Debug, Clone)]
pub struct PairwiseTbaBlock {
    /// Pairwise alignment score (defaults to ungapped match count if unspecified).
    pub score: f64,
    /// First species row.
    pub row_a: TbaRow,
    /// Second species row.
    pub row_b: TbaRow,
}

impl PairwiseTbaBlock {
    /// Convert a parsed pair of [`MafRow`]s with explicit species names into a [`PairwiseTbaBlock`].
    #[must_use]
    pub fn from_maf_rows(sp_a: &str, r_a: MafRow, sp_b: &str, r_b: MafRow) -> Option<Self> {
        if r_a.seq.len() != r_b.seq.len() || r_a.seq.is_empty() {
            return None;
        }
        let a_bytes = r_a.seq.into_bytes();
        let b_bytes = r_b.seq.into_bytes();
        let a_size = a_bytes.iter().filter(|&&b| b != b'-').count();
        let b_size = b_bytes.iter().filter(|&&b| b != b'-').count();
        if a_size == 0 || b_size == 0 {
            return None;
        }
        #[allow(clippy::cast_precision_loss)]
        let score = crate::netting::score_acgt_columns(&a_bytes, &b_bytes) as f64;
        Some(Self {
            score,
            row_a: TbaRow {
                species: sp_a.to_owned(),
                contig: strip_sample_prefix(&r_a.src, sp_a).to_owned(),
                start: r_a.start,
                size: a_size,
                strand: r_a.strand,
                src_size: r_a.src_size,
                seq: a_bytes,
            },
            row_b: TbaRow {
                species: sp_b.to_owned(),
                contig: strip_sample_prefix(&r_b.src, sp_b).to_owned(),
                start: r_b.start,
                size: b_size,
                strand: r_b.strand,
                src_size: r_b.src_size,
                seq: b_bytes,
            },
        })
    }

    /// Materialize a [`crate::mdl::RunAlignment`] directly into [`PairwiseTbaBlock`]s without
    /// serializing and re-parsing intermediate MAF text.
    #[must_use]
    pub fn from_run_alignment(
        sp_a: &str,
        ctg_a: &str,
        seq_a: &[u8],
        sp_b: &str,
        ctg_b: &str,
        seq_b: &[u8],
        aln: &crate::mdl::RunAlignment,
    ) -> Vec<Self> {
        let stripped_a = strip_sample_prefix(ctg_a, sp_a);
        let stripped_b = strip_sample_prefix(ctg_b, sp_b);
        let mut out = Vec::new();
        for (grp, is_rev) in crate::netting::partition_colinear_run_groups(&aln.aligned_runs) {
            let Some(m) = crate::formats::materialize_run_group(
                seq_a.len(),
                seq_b.len(),
                seq_a,
                seq_b,
                &grp,
                is_rev,
            ) else {
                continue;
            };
            let a_size = m.ax.iter().filter(|&&b| b != b'-').count();
            let b_size = m.ay.iter().filter(|&&b| b != b'-').count();
            if a_size == 0 || b_size == 0 {
                continue;
            }
            #[allow(clippy::cast_precision_loss)]
            let score = crate::netting::score_acgt_columns(&m.ax, &m.ay) as f64;
            out.push(Self {
                score,
                row_a: TbaRow {
                    species: sp_a.to_owned(),
                    contig: stripped_a.to_owned(),
                    start: m.x0,
                    size: a_size,
                    strand: '+',
                    src_size: seq_a.len(),
                    seq: m.ax,
                },
                row_b: TbaRow {
                    species: sp_b.to_owned(),
                    contig: stripped_b.to_owned(),
                    start: m.y_start,
                    size: b_size,
                    strand: m.y_strand,
                    src_size: seq_b.len(),
                    seq: m.ay,
                },
            });
        }
        out
    }

    pub(super) fn oriented_for_left_right<'a, S: std::hash::BuildHasher>(
        &'a self,
        left_species: &HashSet<&str, S>,
    ) -> Option<Cow<'a, Self>> {
        let a_in_left = left_species.contains(self.row_a.species.as_str());
        let b_in_left = left_species.contains(self.row_b.species.as_str());
        match (a_in_left, b_in_left) {
            (true, false) if self.row_a.strand == '+' => Some(Cow::Borrowed(self)),
            (true, false) => Some(Cow::Owned(Self {
                score: self.score,
                row_a: self.row_a.revcomp(),
                row_b: self.row_b.revcomp(),
            })),
            (false, true) if self.row_b.strand == '+' => Some(Cow::Owned(Self {
                score: self.score,
                row_a: self.row_b.clone(),
                row_b: self.row_a.clone(),
            })),
            (false, true) => Some(Cow::Owned(Self {
                score: self.score,
                row_a: self.row_b.revcomp(),
                row_b: self.row_a.revcomp(),
            })),
            _ => None,
        }
    }
}

/// Parse all 2-row MAF stanzas from a buffered reader into [`PairwiseTbaBlock`]s.
///
/// Automatically detects and swaps row order when `r0.src`/`r1.src` carry explicit species
/// prefixes matching `(sp_b, sp_a)` instead of `(sp_a, sp_b)`.
///
/// # Errors
/// Returns an error if reading fails.
pub fn parse_pairwise_tba_reader<R: std::io::BufRead>(
    sp_a: &str,
    sp_b: &str,
    reader: R,
) -> anyhow::Result<Vec<PairwiseTbaBlock>> {
    let mut out = Vec::new();
    for_each_maf_stanza_reader(reader, |rows| {
        if rows.len() != 2 {
            return;
        }
        let mut it = rows.drain(..);
        let (Some((r0, _)), Some((r1, _))) = (it.next(), it.next()) else {
            return;
        };
        let r0_is_a = matches_sample_prefix(&r0.src, sp_a);
        let r0_is_b = matches_sample_prefix(&r0.src, sp_b);
        let r1_is_a = matches_sample_prefix(&r1.src, sp_a);
        let r1_is_b = matches_sample_prefix(&r1.src, sp_b);
        let (ra, rb) = if r0_is_b && r1_is_a && !r0_is_a && !r1_is_b {
            (r1, r0)
        } else {
            (r0, r1)
        };
        if let Some(blk) = PairwiseTbaBlock::from_maf_rows(sp_a, ra, sp_b, rb) {
            out.push(blk);
        }
    })?;
    Ok(out)
}

#[derive(Default)]
struct MultiMafStanza {
    score: Option<f64>,
    rows: Vec<TbaRow>,
    species: HashSet<String>,
    cols: Option<usize>,
    valid: bool,
    skipped_dup: bool,
}

fn split_maf_src_species_contig(src: &str) -> (&str, &str) {
    // Recognize versioned NCBI assembly accessions (`GCA_000002765.3.<contig>` or `GCF_...`)
    // so dots inside the accession version are not misread as the species/contig separator.
    let b = src.as_bytes();
    if b.len() > 6 && (src.starts_with("GCA_") || src.starts_with("GCF_")) {
        let mut i = 4usize;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
        if i > 4 && i < b.len() && b[i] == b'.' {
            let mut j = i + 1;
            while j < b.len() && b[j].is_ascii_digit() {
                j += 1;
            }
            if j > i + 1 && j + 1 < b.len() && b[j] == b'.' {
                return (&src[..j], &src[j + 1..]);
            }
        }
    }
    src.split_once('.').map_or((src, src), |(s, c)| (s, c))
}

impl MultiMafStanza {
    fn new() -> Self {
        Self {
            valid: true,
            ..Self::default()
        }
    }

    fn flush_with<F: FnMut(TbaBlock)>(&mut self, mut on_block: F) {
        if self.valid && self.rows.len() >= 2 {
            let score = self.score.unwrap_or_else(|| {
                #[allow(clippy::cast_precision_loss)]
                let fallback = self.rows.iter().map(|r| r.size).sum::<usize>() as f64;
                fallback
            });
            let mut rows = std::mem::take(&mut self.rows);
            if self.skipped_dup {
                let n_cols = rows.first().map_or(0, |r| r.seq.len());
                let keep: Vec<bool> = (0..n_cols)
                    .map(|c| rows.iter().any(|r| r.seq[c] != b'-'))
                    .collect();
                if keep.iter().any(|&k| !k) {
                    for r in &mut rows {
                        r.seq = r
                            .seq
                            .iter()
                            .zip(&keep)
                            .filter_map(|(&b, &k)| k.then_some(b))
                            .collect();
                    }
                }
            }
            on_block(TbaBlock { score, rows });
        } else {
            self.rows.clear();
        }
        self.score = None;
        self.species.clear();
        self.cols = None;
        self.valid = true;
        self.skipped_dup = false;
    }

    fn parse_a_header(&mut self, trimmed: &str) {
        for tok in trimmed.split_ascii_whitespace().skip(1) {
            if let Some(raw_sc) = tok.strip_prefix("score=")
                && let Ok(v) = raw_sc.parse::<f64>()
                && v.is_finite()
            {
                self.score = Some(v);
            }
        }
    }

    fn push_s_line(&mut self, trimmed: &str) {
        if !self.valid {
            return;
        }
        let Some((mr, ungapped)) = crate::formats::parse_maf_s_line(trimmed) else {
            self.valid = false;
            return;
        };
        if ungapped == 0 {
            return;
        }
        let seq_bytes = mr.seq.into_bytes();
        if let Some(expected) = self.cols {
            if seq_bytes.len() != expected {
                self.valid = false;
                return;
            }
        } else {
            self.cols = Some(seq_bytes.len());
        }
        let (sp, ctg) = split_maf_src_species_contig(&mr.src);
        if !self.species.insert(sp.to_owned()) {
            // Keep first row per species inside a block (TBA single-occurrence invariant).
            self.skipped_dup = true;
            return;
        }
        self.rows.push(TbaRow {
            species: sp.to_owned(),
            contig: ctg.to_owned(),
            start: mr.start,
            size: ungapped,
            strand: mr.strand,
            src_size: mr.src_size,
            seq: seq_bytes,
        });
    }
}

/// Stream multi-species MAF stanzas (`>= 2` rows) from a buffered reader, invoking `on_block`
/// for each parsed [`TbaBlock`] without holding the entire MAF in memory.
///
/// # Errors
/// Returns an error if reading fails, a single line exceeds 256 MiB, or `on_block` returns an error.
pub fn for_each_multi_tba_block_reader<R, F>(mut reader: R, mut on_block: F) -> anyhow::Result<()>
where
    R: std::io::BufRead,
    F: FnMut(TbaBlock) -> anyhow::Result<()>,
{
    const MAX_LINE_BYTES: usize = 256 * 1024 * 1024;
    let mut stanza = MultiMafStanza::new();
    let mut line = String::new();
    let mut cb_err: Option<anyhow::Error> = None;
    loop {
        line.clear();
        let n = std::io::BufRead::read_line(
            &mut std::io::Read::take(&mut reader, MAX_LINE_BYTES as u64 + 1),
            &mut line,
        )?;
        if n == 0 {
            break;
        }
        if line.len() > MAX_LINE_BYTES {
            anyhow::bail!("MAF line exceeds 256 MiB safety limit");
        }
        let trimmed = line.trim();
        if trimmed.is_empty() {
            stanza.flush_with(|blk| {
                if cb_err.is_none()
                    && let Err(e) = on_block(blk)
                {
                    cb_err = Some(e);
                }
            });
            if let Some(e) = cb_err {
                return Err(e);
            }
        } else if trimmed == "a" || trimmed.starts_with("a ") || trimmed.starts_with("a\t") {
            stanza.flush_with(|blk| {
                if cb_err.is_none()
                    && let Err(e) = on_block(blk)
                {
                    cb_err = Some(e);
                }
            });
            if let Some(e) = cb_err {
                return Err(e);
            }
            stanza.parse_a_header(trimmed);
        } else if trimmed.starts_with("s ") || trimmed.starts_with("s\t") {
            stanza.push_s_line(trimmed);
        } else if !trimmed.starts_with('#')
            && !trimmed
                .strip_prefix(['i', 'e', 'q'])
                .is_some_and(|r| r.starts_with(char::is_whitespace))
        {
            anyhow::bail!("malformed MAF line: unrecognized record prefix");
        }
    }
    stanza.flush_with(|blk| {
        if cb_err.is_none()
            && let Err(e) = on_block(blk)
        {
            cb_err = Some(e);
        }
    });
    if let Some(e) = cb_err {
        return Err(e);
    }
    Ok(())
}

/// Parse all multi-species MAF stanzas (`>= 2` rows) from a buffered reader into [`TbaBlock`]s,
/// preserving `a score=` values and `<species>.<contig>` prefixes (for post-hoc tip projection
/// onto an existing TBA/MULTIZ backbone MAF).
///
/// # Errors
/// Returns an error if reading fails or a single line exceeds 256 MiB.
pub fn parse_multi_tba_reader<R: std::io::BufRead>(reader: R) -> anyhow::Result<Vec<TbaBlock>> {
    let mut out = Vec::new();
    for_each_multi_tba_block_reader(reader, |blk| {
        out.push(blk);
        Ok(())
    })?;
    Ok(out)
}

/// Canonical unordered species-pair key `(min(a, b), max(a, b))`.
#[must_use]
pub fn canonical_pair_key(a: &str, b: &str) -> (String, String) {
    if a <= b {
        (a.to_owned(), b.to_owned())
    } else {
        (b.to_owned(), a.to_owned())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multi_maf_preserves_versioned_ncbi_accession_species_names() {
        let maf = "\
##maf version=1
a score=20.0
s GCA_000002765.3.chr1 0 8 + 100 ACGTACGT
s GCF_002759435.1.PEKT02000001.1 0 8 + 100 ACGTACGT

";
        let blocks = parse_multi_tba_reader(std::io::Cursor::new(maf)).expect("parse");
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].rows[0].species, "GCA_000002765.3");
        assert_eq!(blocks[0].rows[0].contig, "chr1");
        assert_eq!(blocks[0].rows[1].species, "GCF_002759435.1");
        assert_eq!(blocks[0].rows[1].contig, "PEKT02000001.1");
    }
}
