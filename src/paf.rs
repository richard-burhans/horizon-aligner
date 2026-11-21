//! Streaming Pairwise mApping Format (PAF) parser, writer, and CIGAR iterator.
//!
//! Supports `.paf`, multi-member BGZF `.paf.gz`, `.paf.zst`, and `"-"` (`stdin`) without
//! materializing entire alignment files in memory.

use crate::io_util::open_text_reader;
use anyhow::{Context, Result, bail};
use std::io::BufRead;
use std::ops::ControlFlow;
use std::path::Path;

/// A single CIGAR operation decoded from a PAF `cg:Z:` tag.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    /// `M`: alignment match (sequence match or mismatch).
    Match(usize),
    /// `I`: insertion to the target (consumes query only).
    Ins(usize),
    /// `D`: deletion from the target (consumes target only).
    Del(usize),
    /// `N`: skipped region from the target (consumes target only).
    RefSkip(usize),
    /// `S`: soft clipping (consumes query only).
    SoftClip(usize),
    /// `H`: hard clipping (consumes neither).
    HardClip(usize),
    /// `P`: padding (consumes neither).
    Pad(usize),
    /// `=`: sequence match.
    Equal(usize),
    /// `X`: sequence mismatch.
    Diff(usize),
}

/// Zero-allocation iterator over a PAF `cg:Z:` CIGAR string.
pub struct CigarIter<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> CigarIter<'a> {
    /// Create a new CIGAR iterator from the raw `cg:Z:` value (without the `cg:Z:` prefix).
    #[must_use]
    pub const fn new(cigar: &'a str) -> Self {
        Self {
            bytes: cigar.as_bytes(),
            pos: 0,
        }
    }
}

impl Iterator for CigarIter<'_> {
    type Item = Result<CigarOp>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.bytes.len() {
            return None;
        }
        let mut len = 0usize;
        let mut saw_digit = false;
        while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_digit() {
            saw_digit = true;
            len = match len
                .checked_mul(10)
                .and_then(|v| v.checked_add(usize::from(self.bytes[self.pos] - b'0')))
            {
                Some(v) => v,
                None => return Some(Err(anyhow::anyhow!("CIGAR length overflow"))),
            };
            self.pos += 1;
        }
        if !saw_digit || self.pos >= self.bytes.len() {
            self.pos = self.bytes.len();
            return Some(Err(anyhow::anyhow!("malformed CIGAR string")));
        }
        let op = self.bytes[self.pos];
        self.pos += 1;
        let parsed = match op {
            b'M' => CigarOp::Match(len),
            b'I' => CigarOp::Ins(len),
            b'D' => CigarOp::Del(len),
            b'N' => CigarOp::RefSkip(len),
            b'S' => CigarOp::SoftClip(len),
            b'H' => CigarOp::HardClip(len),
            b'P' => CigarOp::Pad(len),
            b'=' => CigarOp::Equal(len),
            b'X' => CigarOp::Diff(len),
            other => {
                self.pos = self.bytes.len();
                return Some(Err(anyhow::anyhow!(
                    "unsupported CIGAR op '{}'",
                    char::from(other)
                )));
            }
        };
        Some(Ok(parsed))
    }
}

/// A parsed 12-column PAF record plus optional SAM-style tags.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PafRecord {
    /// Query sequence name.
    pub qname: String,
    /// Query sequence length.
    pub qlen: usize,
    /// Query start (0-based, BED-like).
    pub qstart: usize,
    /// Query end (0-based, half-open).
    pub qend: usize,
    /// Relative strand (`'+'` or `'-'`).
    pub strand: char,
    /// Target sequence name.
    pub tname: String,
    /// Target sequence length.
    pub tlen: usize,
    /// Target start (0-based, BED-like).
    pub tstart: usize,
    /// Target end (0-based, half-open).
    pub tend: usize,
    /// Number of residue matches.
    pub matches: usize,
    /// Alignment block length.
    pub aln_len: usize,
    /// Mapping quality (`0..=255`).
    pub mapq: u8,
    /// Optional trailing SAM-style tags (e.g., `cg:Z:...`, `dv:f:...`).
    pub tags: Vec<String>,
}

impl PafRecord {
    /// Parse a single non-empty PAF line.
    ///
    /// # Errors
    ///
    /// Returns an error if fewer than 12 tab-delimited columns are present or numeric fields fail to parse.
    pub fn parse_line(line: &str) -> Result<Self> {
        let mut fields = line.split('\t');
        let next_str = |it: &mut std::str::Split<'_, char>, name: &str| -> Result<String> {
            it.next()
                .map(str::to_owned)
                .with_context(|| format!("missing PAF field {name}"))
        };
        let next_usize = |it: &mut std::str::Split<'_, char>, name: &str| -> Result<usize> {
            let s = it
                .next()
                .with_context(|| format!("missing PAF field {name}"))?;
            s.parse::<usize>()
                .with_context(|| format!("invalid integer in PAF field {name}: '{s}'"))
        };

        let qname = next_str(&mut fields, "qname")?;
        let qlen = next_usize(&mut fields, "qlen")?;
        let qstart = next_usize(&mut fields, "qstart")?;
        let qend = next_usize(&mut fields, "qend")?;
        let strand_str = fields.next().context("missing PAF field strand")?;
        let strand = match strand_str {
            "+" => '+',
            "-" => '-',
            other => bail!("invalid PAF strand '{other}'"),
        };
        let tname = next_str(&mut fields, "tname")?;
        let tlen = next_usize(&mut fields, "tlen")?;
        let tstart = next_usize(&mut fields, "tstart")?;
        let tend = next_usize(&mut fields, "tend")?;
        let matches = next_usize(&mut fields, "matches")?;
        let aln_len = next_usize(&mut fields, "aln_len")?;
        let mapq_str = fields.next().context("missing PAF field mapq")?;
        let mapq = mapq_str
            .parse::<u8>()
            .with_context(|| format!("invalid PAF mapq '{mapq_str}'"))?;
        let tags: Vec<String> = fields
            .filter(|t| !t.is_empty())
            .map(str::to_owned)
            .collect();

        Ok(Self {
            qname,
            qlen,
            qstart,
            qend,
            strand,
            tname,
            tlen,
            tstart,
            tend,
            matches,
            aln_len,
            mapq,
            tags,
        })
    }

    /// Return the value of a SAM-style tag with the given prefix (e.g., `"cg:Z:"`), if present.
    #[must_use]
    pub fn tag_value(&self, prefix: &str) -> Option<&str> {
        self.tags
            .iter()
            .find_map(|tag| tag.strip_prefix(prefix))
    }

    /// Return the raw `cg:Z:` CIGAR string, if present.
    #[must_use]
    pub fn cigar(&self) -> Option<&str> {
        self.tag_value("cg:Z:")
    }

    /// Iterate over CIGAR operations in `cg:Z:`, if present.
    #[must_use]
    pub fn cigar_ops(&self) -> Option<CigarIter<'_>> {
        self.cigar().map(CigarIter::new)
    }

    /// Format as a single PAF line (without trailing newline), appending `extra_tags`.
    #[must_use]
    pub fn format_line(&self, extra_tags: &[String]) -> String {
        use std::fmt::Write as _;
        let mut out = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.qname,
            self.qlen,
            self.qstart,
            self.qend,
            self.strand,
            self.tname,
            self.tlen,
            self.tstart,
            self.tend,
            self.matches,
            self.aln_len,
            self.mapq
        );
        for tag in &self.tags {
            let _ = write!(out, "\t{tag}");
        }
        for tag in extra_tags {
            let _ = write!(out, "\t{tag}");
        }
        out
    }
}

fn stream_paf_reader<R: BufRead, F>(mut reader: R, mut on_record: F) -> Result<()>
where
    F: FnMut(PafRecord) -> Result<ControlFlow<()>>,
{
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(&['\r', '\n'][..]);
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let rec = PafRecord::parse_line(trimmed)?;
        if matches!(on_record(rec)?, ControlFlow::Break(())) {
            break;
        }
    }
    Ok(())
}

/// Stream PAF records from a file path (`.paf`, `.paf.gz`, `.paf.zst`, or `"-"` for `stdin`).
///
/// # Errors
///
/// Returns an error if the file cannot be opened or any PAF line is malformed.
pub fn stream_paf_records<F>(path: &str, on_record: F) -> Result<()>
where
    F: FnMut(PafRecord) -> Result<ControlFlow<()>>,
{
    let reader = open_text_reader(Path::new(path))
        .with_context(|| format!("opening PAF '{path}'"))?;
    stream_paf_reader(reader, on_record)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_paf_and_cigar_roundtrip() {
        let line = "qry1\t1000\t10\t50\t-\ttgt1\t2000\t100\t141\t35\t42\t60\tdv:f:0.10\tcg:Z:10=1X10=1D19=";
        let rec = PafRecord::parse_line(line).expect("valid paf");
        assert_eq!(rec.qname, "qry1");
        assert_eq!(rec.strand, '-');
        assert_eq!(rec.matches, 35);
        assert_eq!(rec.aln_len, 42);
        let ops: Vec<CigarOp> = rec
            .cigar_ops()
            .expect("cigar")
            .collect::<Result<Vec<_>>>()
            .expect("valid ops");
        assert_eq!(
            ops,
            vec![
                CigarOp::Equal(10),
                CigarOp::Diff(1),
                CigarOp::Equal(10),
                CigarOp::Del(1),
                CigarOp::Equal(19),
            ]
        );
        assert_eq!(rec.format_line(&[]), line);
    }
}
