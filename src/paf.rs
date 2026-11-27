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
            let Some(next_len) = len
                .checked_mul(10)
                .and_then(|v| v.checked_add(usize::from(self.bytes[self.pos] - b'0')))
            else {
                self.pos = self.bytes.len();
                return Some(Err(anyhow::anyhow!("CIGAR length overflow")));
            };
            len = next_len;
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
        anyhow::ensure!(
            qstart <= qend && qend <= qlen,
            "invalid PAF query coordinates {qstart}..{qend} for qlen {qlen}"
        );
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
        anyhow::ensure!(
            tstart <= tend && tend <= tlen,
            "invalid PAF target coordinates {tstart}..{tend} for tlen {tlen}"
        );
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
        self.tags.iter().find_map(|tag| tag.strip_prefix(prefix))
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

    /// Write a single PAF line directly to `out` (with trailing newline) plus `md:f:` and `hb:i:` tags,
    /// avoiding per-record `String` and `Vec` heap allocations.
    ///
    /// # Errors
    /// Returns an `std::io::Error` if writing to `out` fails.
    pub fn write_line_with_mdl(
        &self,
        out: &mut impl std::io::Write,
        mdl_bits: f64,
        hb: u8,
    ) -> std::io::Result<()> {
        write!(
            out,
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
        )?;
        for tag in &self.tags {
            if tag.starts_with("md:f:") || tag.starts_with("mb:f:") || tag.starts_with("hb:i:") {
                continue;
            }
            write!(out, "\t{tag}")?;
        }
        writeln!(out, "\tmd:f:{mdl_bits:.2}\thb:i:{hb}")
    }
}

const MAX_PAF_LINE_BYTES: usize = 256 * 1024 * 1024;

fn stream_paf_reader<R: BufRead, F>(mut reader: R, mut on_record: F) -> Result<()>
where
    F: FnMut(PafRecord) -> Result<ControlFlow<()>>,
{
    let mut line = String::new();
    loop {
        line.clear();
        let n = std::io::Read::take(&mut reader, (MAX_PAF_LINE_BYTES + 1) as u64)
            .read_line(&mut line)?;
        if n == 0 {
            break;
        }
        if n > MAX_PAF_LINE_BYTES {
            bail!("PAF line exceeds {MAX_PAF_LINE_BYTES} bytes");
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
    let reader =
        open_text_reader(Path::new(path)).with_context(|| format!("opening PAF '{path}'"))?;
    stream_paf_reader(reader, on_record)
}

/// Format a pairwise [`crate::mdl::RunAlignment`] into standard 12-column PAF v1 lines with
/// extended `cg:Z:` CIGAR (`=`/`X`/`I`/`D`) and `tp:A:P` tags.
///
/// Target coordinates correspond to `seq_x` (`tname`, `tlen`) and query coordinates correspond to
/// `seq_y` (`qname`, `qlen`) using 0-based half-open forward-strand PAF conventions.
#[must_use]
pub fn format_paf_records_runs(
    target_name: &str,
    target_len: usize,
    query_name: &str,
    query_len: usize,
    seq_x: &[u8],
    seq_y: &[u8],
    aln: &crate::mdl::RunAlignment,
) -> String {
    format_paf_records_runs_avoiding_gaps(
        target_name,
        target_len,
        query_name,
        query_len,
        seq_x,
        seq_y,
        aln,
        None,
        None,
    )
}

/// Build extended `cg:Z:` CIGAR (`=`/`X`/`I`/`D`) and composition-aware BOCPD+MDL bit-savings
/// for a materialized alignment block `(ax, ay)`.
fn build_cigar_and_mdl_bits(
    ax: &[u8],
    ay: &[u8],
    gc_x: &[f64],
    gc_y: &[f64],
    cigar: &mut String,
) -> (usize, usize, f64) {
    use crate::bocpd::{
        GAP_EXTEND_BITS, GAP_OPEN_BITS, IDENTITY_PRIOR, base_probs, joint_pair_bits,
    };
    use std::fmt::Write as _;

    cigar.clear();
    let mut matches = 0usize;
    let mut aln_len = 0usize;
    let mut cur_op: Option<char> = None;
    let mut cur_len = 0usize;
    let mut ix = 0usize;
    let mut iy = 0usize;
    let mut mdl_bits = 0.0f64;

    for (&bx, &by) in ax.iter().zip(ay.iter()) {
        let op = match (bx, by) {
            (b'-', b'-') => continue,
            (b'-', _) => {
                mdl_bits += if cur_op == Some('I') {
                    GAP_EXTEND_BITS
                } else {
                    GAP_OPEN_BITS
                };
                iy += 1;
                'I'
            }
            (_, b'-') => {
                mdl_bits += if cur_op == Some('D') {
                    GAP_EXTEND_BITS
                } else {
                    GAP_OPEN_BITS
                };
                ix += 1;
                'D'
            }
            (a, b) => {
                let px = base_probs(gc_x.get(ix).copied().unwrap_or(0.5));
                let py = base_probs(gc_y.get(iy).copied().unwrap_or(0.5));
                mdl_bits += joint_pair_bits(a, b, px, py, IDENTITY_PRIOR);
                ix += 1;
                iy += 1;
                if a.eq_ignore_ascii_case(&b)
                    && matches!(a.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
                {
                    matches += 1;
                    '='
                } else {
                    'X'
                }
            }
        };
        aln_len += 1;
        if cur_op == Some(op) {
            cur_len += 1;
        } else {
            if let Some(prev) = cur_op {
                let _ = write!(cigar, "{cur_len}{prev}");
            }
            cur_op = Some(op);
            cur_len = 1;
        }
    }
    if let Some(prev) = cur_op {
        let _ = write!(cigar, "{cur_len}{prev}");
    }
    (matches, aln_len, mdl_bits)
}

/// Format [`crate::mdl::RunAlignment`] into PAF records, refusing to bridge gaps overlapping
/// `forbidden_x` or `forbidden_y` (enforcing strict 0-bp overlap under `--one-to-one`).
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn format_paf_records_runs_avoiding_gaps(
    target_name: &str,
    target_len: usize,
    query_name: &str,
    query_len: usize,
    seq_x: &[u8],
    seq_y: &[u8],
    aln: &crate::mdl::RunAlignment,
    forbidden_x: Option<&crate::coord::IntervalSet>,
    forbidden_y: Option<&crate::coord::IntervalSet>,
) -> String {
    use crate::bocpd::{BocpdScratch, bocpd_gc_profile_into};
    use std::fmt::Write as _;

    let mut out = String::new();
    let mut cigar = String::new();
    let mut bocpd_scratch = BocpdScratch::default();
    let mut ux: Vec<u8> = Vec::new();
    let mut uy: Vec<u8> = Vec::new();
    let mut gc_x: Vec<f64> = Vec::new();
    let mut gc_y: Vec<f64> = Vec::new();

    for (grp, is_rev) in crate::netting::partition_colinear_run_groups_avoiding_gaps(
        &aln.aligned_runs,
        forbidden_x,
        forbidden_y,
    ) {
        let Some(m) = crate::formats::materialize_run_group(
            target_len, query_len, seq_x, seq_y, &grp, is_rev,
        ) else {
            continue;
        };
        ux.clear();
        uy.clear();
        for (&bx, &by) in m.ax.iter().zip(m.ay.iter()) {
            if bx == b'-' && by == b'-' {
                continue;
            }
            if bx != b'-' {
                ux.push(bx);
            }
            if by != b'-' {
                uy.push(by);
            }
        }
        let x_ungapped = ux.len();
        let y_ungapped = uy.len();
        if x_ungapped == 0 || y_ungapped == 0 {
            continue;
        }
        bocpd_gc_profile_into(&ux, 1.0 / 120.0, 80, &mut bocpd_scratch, &mut gc_x);
        bocpd_gc_profile_into(&uy, 1.0 / 120.0, 80, &mut bocpd_scratch, &mut gc_y);

        let tstart = m.x0;
        let tend = tstart.saturating_add(x_ungapped).min(target_len);
        let (qstart, qend) = if m.y_strand == '-' {
            let qe = query_len.saturating_sub(m.y_start);
            let qs = qe.saturating_sub(y_ungapped);
            (qs, qe)
        } else {
            let qs = m.y_start;
            let qe = qs.saturating_add(y_ungapped).min(query_len);
            (qs, qe)
        };

        let (matches, aln_len, mdl_bits) =
            build_cigar_and_mdl_bits(&m.ax, &m.ay, &gc_x, &gc_y, &mut cigar);
        #[allow(clippy::cast_precision_loss)]
        let bits_per_col = mdl_bits / (aln_len.max(1) as f64);
        let _ = writeln!(
            out,
            "{query_name}\t{query_len}\t{qstart}\t{qend}\t{}\t{target_name}\t{target_len}\t{tstart}\t{tend}\t{matches}\t{aln_len}\t60\ttp:A:P\tmd:f:{mdl_bits:.2}\tmb:f:{bits_per_col:.4}\tcg:Z:{cigar}",
            m.y_strand
        );
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::too_many_lines)]
    fn format_paf_records_runs_roundtrips_plus_and_minus_strands() {
        use crate::mdl::{AlignedRun, RunAlignment};

        let seq_x = b"ACGTACGTGGGGTTTTACGTACGT";
        let seq_y = b"ACGTACGTTTGGGGTTTTACGTACGT";
        let aln_fwd = RunAlignment {
            aligned_runs: vec![
                AlignedRun {
                    x_start: 0,
                    y_start: 0,
                    len: 8,
                    is_rc: false,
                },
                AlignedRun {
                    x_start: 8,
                    y_start: 10,
                    len: 16,
                    is_rc: false,
                },
            ],
            total_pairs: 24,
            lcbs_fwd: 1,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let paf = format_paf_records_runs(
            "tgt1",
            seq_x.len(),
            "qry1",
            seq_y.len(),
            seq_x,
            seq_y,
            &aln_fwd,
        );
        let rec = PafRecord::parse_line(paf.trim()).expect("valid fwd paf");
        assert_eq!(rec.strand, '+');
        assert_eq!((rec.tstart, rec.tend), (0, 24));
        assert_eq!((rec.qstart, rec.qend), (0, 26));
        assert_eq!(rec.matches, 24);
        assert_eq!(rec.aln_len, 26);
        let ops: Vec<CigarOp> = rec
            .cigar_ops()
            .expect("cg")
            .collect::<Result<Vec<_>>>()
            .expect("ops");
        assert_eq!(
            ops,
            vec![CigarOp::Equal(8), CigarOp::Ins(2), CigarOp::Equal(16)]
        );
        let score_bits = rec
            .tags
            .iter()
            .find_map(|t| t.strip_prefix("md:f:"))
            .and_then(|v| v.parse::<f64>().ok())
            .expect("horizon_wga PAF must emit continuous md:f: bit-savings tag");
        let density_bits = rec
            .tags
            .iter()
            .find_map(|t| t.strip_prefix("mb:f:"))
            .and_then(|v| v.parse::<f64>().ok())
            .expect("horizon_wga PAF must emit per-column mb:f: bit-savings density tag");
        assert!(
            score_bits > 10.0 && score_bits.is_finite(),
            "md:f:={score_bits}"
        );
        let expected_density = score_bits / f64::from(u32::try_from(rec.aln_len).expect("u32"));
        assert!(
            (density_bits - expected_density).abs() < 1e-3,
            "mb:f:={density_bits} vs md/aln={expected_density}"
        );
        let mut rescore_buf = Vec::new();
        rec.write_line_with_mdl(&mut rescore_buf, 42.5, 0)
            .expect("write_line_with_mdl");
        let rescore_line = String::from_utf8(rescore_buf).expect("utf8");
        assert_eq!(
            rescore_line.matches("\tmd:f:").count(),
            1,
            "write_line_with_mdl must replace existing md:f: instead of duplicating:\n{rescore_line}"
        );
        assert!(rescore_line.contains("\tmd:f:42.50\thb:i:0"));

        let seq_y_rc = crate::seeding::revcomp(seq_x);
        let aln_rev = RunAlignment {
            aligned_runs: vec![AlignedRun {
                x_start: 2,
                y_start: 21,
                len: 18,
                is_rc: true,
            }],
            total_pairs: 18,
            lcbs_fwd: 0,
            lcbs_rc: 1,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let paf_rev = format_paf_records_runs(
            "tgt1",
            seq_x.len(),
            "qry_rc",
            seq_y_rc.len(),
            seq_x,
            &seq_y_rc,
            &aln_rev,
        );
        let rec_rev = PafRecord::parse_line(paf_rev.trim()).expect("valid rev paf");
        assert_eq!(rec_rev.strand, '-');
        assert_eq!((rec_rev.tstart, rec_rev.tend), (2, 20));
        assert_eq!((rec_rev.qstart, rec_rev.qend), (4, 22));
        assert_eq!(rec_rev.matches, 18);

        // Wave B regression: N==N columns must emit 'X' in cg:Z: and never inflate `matches`.
        let ref_n = b"ACGTNNNNACGT";
        let qry_n = b"ACGTNNNNACGT";
        let aln_n = RunAlignment {
            aligned_runs: vec![AlignedRun {
                x_start: 0,
                y_start: 0,
                len: 12,
                is_rc: false,
            }],
            total_pairs: 12,
            lcbs_fwd: 1,
            lcbs_rc: 0,
            macro_bubbles: 0,
            internal_bubbles: 0,
        };
        let paf_n = format_paf_records_runs("t", 12, "q", 12, ref_n, qry_n, &aln_n);
        let parsed_n = PafRecord::parse_line(paf_n.trim()).expect("valid N-tract paf");
        assert_eq!(parsed_n.matches, 8);
        assert_eq!(parsed_n.aln_len, 12);
        assert_eq!(parsed_n.cigar(), Some("4=4X4="));
    }

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
        let mut buf = Vec::new();
        rec.write_line_with_mdl(&mut buf, 12.5, 0).expect("write");
        let rendered = String::from_utf8(buf).expect("utf8");
        assert!(rendered.starts_with(line));
        assert!(rendered.contains("\tmd:f:12.50\thb:i:0\n"));

        assert!(PafRecord::parse_line("q\t100\t60\t50\t+\tt\t100\t0\t10\t10\t10\t60").is_err());
        assert!(PafRecord::parse_line("q\t100\t0\t101\t+\tt\t100\t0\t10\t10\t10\t60").is_err());
        assert!(PafRecord::parse_line("q\t100\t0\t10\t+\tt\t100\t0\t101\t10\t10\t60").is_err());
    }
}
