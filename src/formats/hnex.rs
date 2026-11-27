//! Lossless `.hnex` multi-scale alignment & Horizon Bubble provenance archive I/O.

use anyhow::{Context, Result, bail};
use std::io::{BufRead, Write};

/// Sequence metadata entry in an `.hnex` archive header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HnexSeqMeta {
    pub species: String,
    pub contig: String,
    pub length: usize,
}

/// One aligned row inside an `.hnex` homology block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HnexRow {
    pub species: String,
    pub contig: String,
    pub start: usize,
    pub size: usize,
    pub strand: char,
    pub src_size: usize,
    pub seq: String,
}

/// An aligned homology block with optional per-column posterior homology probabilities.
#[derive(Debug, Clone, PartialEq)]
pub struct HnexBlock {
    pub block_id: usize,
    pub score_bits: f64,
    pub rows: Vec<HnexRow>,
    pub col_posteriors: Vec<f32>,
}

/// An explicit unaligned sequence span enclosed inside a Horizon Bubble.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HnexBubbleSpan {
    pub species: String,
    pub contig: String,
    pub start: usize,
    pub end: usize,
    pub seq: String,
}

/// Explicit epistemic abstention record (Coarse Horizon Bubble) with MDL evidence delta.
#[derive(Debug, Clone, PartialEq)]
pub struct HnexBubble {
    pub bubble_id: usize,
    pub left_block: Option<usize>,
    pub right_block: Option<usize>,
    pub delta_l_bits: f64,
    pub reason: String,
    pub spans: Vec<HnexBubbleSpan>,
}

/// Lossless `.hnex` multi-scale alignment, guide-tree, and Horizon Bubble archive.
#[derive(Debug, Clone, PartialEq)]
pub struct HnexArchive {
    pub version: String,
    pub guide_tree_newick: Option<String>,
    pub sequences: Vec<HnexSeqMeta>,
    pub blocks: Vec<HnexBlock>,
    pub bubbles: Vec<HnexBubble>,
}

impl HnexArchive {
    /// Serialize the archive in lossless line-oriented `.hnex` format.
    ///
    /// # Errors
    /// Returns an error if any field contains unescaped tab/newline control characters or writing fails.
    pub fn write_to<W: Write>(&self, out: &mut W) -> Result<()> {
        let check_field = |val: &str, name: &str| -> Result<()> {
            anyhow::ensure!(
                !val.bytes().any(|b| matches!(b, b'\t' | b'\n' | b'\r')),
                ".hnex field `{name}` must not contain tab or newline characters"
            );
            Ok(())
        };
        check_field(&self.version, "version")?;
        writeln!(out, "HNEX\tVN:{}", self.version)?;
        if let Some(ref tree) = self.guide_tree_newick {
            check_field(tree, "guide_tree_newick")?;
            writeln!(out, "TREE\t{tree}")?;
        }
        for s in &self.sequences {
            check_field(&s.species, "SEQ.species")?;
            check_field(&s.contig, "SEQ.contig")?;
            writeln!(out, "SEQ\t{}\t{}\t{}", s.species, s.contig, s.length)?;
        }
        for blk in &self.blocks {
            anyhow::ensure!(
                blk.score_bits.is_finite(),
                ".hnex BLK {} has non-finite score",
                blk.block_id
            );
            writeln!(out, "BLK\t{}\t{:.6}", blk.block_id, blk.score_bits)?;
            if !blk.col_posteriors.is_empty() {
                write!(out, "POST\t")?;
                for (i, &p) in blk.col_posteriors.iter().enumerate() {
                    anyhow::ensure!(
                        p.is_finite() && (0.0..=1.0).contains(&p),
                        ".hnex POST value out of [0,1]: {p}"
                    );
                    if i > 0 {
                        write!(out, ",")?;
                    }
                    write!(out, "{p:.4}")?;
                }
                writeln!(out)?;
            }
            for r in &blk.rows {
                check_field(&r.species, "ROW.species")?;
                check_field(&r.contig, "ROW.contig")?;
                check_field(&r.seq, "ROW.seq")?;
                writeln!(
                    out,
                    "ROW\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    r.species, r.contig, r.start, r.size, r.strand, r.src_size, r.seq
                )?;
            }
        }
        for bub in &self.bubbles {
            check_field(&bub.reason, "BUB.reason")?;
            anyhow::ensure!(
                bub.delta_l_bits.is_finite(),
                ".hnex BUB {} has non-finite delta_l_bits",
                bub.bubble_id
            );
            let left_s = bub
                .left_block
                .map_or_else(|| "-".to_owned(), |id| id.to_string());
            let right_s = bub
                .right_block
                .map_or_else(|| "-".to_owned(), |id| id.to_string());
            writeln!(
                out,
                "BUB\t{}\t{}\t{}\t{:.6}\t{}",
                bub.bubble_id, left_s, right_s, bub.delta_l_bits, bub.reason
            )?;
            for sp in &bub.spans {
                check_field(&sp.species, "SPAN.species")?;
                check_field(&sp.contig, "SPAN.contig")?;
                check_field(&sp.seq, "SPAN.seq")?;
                writeln!(
                    out,
                    "SPAN\t{}\t{}\t{}\t{}\t{}",
                    sp.species, sp.contig, sp.start, sp.end, sp.seq
                )?;
            }
        }
        Ok(())
    }

    /// Deserialize an `.hnex` archive losslessly from a buffered reader.
    ///
    /// # Errors
    /// Returns an error if the header or any record is malformed, violates block/row invariants,
    /// or exceeds the 256 MiB single-line safety limit.
    #[allow(clippy::too_many_lines)]
    pub fn read_from<R: BufRead>(mut reader: R) -> Result<Self> {
        const MAX_HNEX_LINE_BYTES: usize = 256 * 1024 * 1024;
        let mut version: Option<String> = None;
        let mut guide_tree_newick = None;
        let mut sequences = Vec::new();
        let mut blocks: Vec<HnexBlock> = Vec::new();
        let mut bubbles: Vec<HnexBubble> = Vec::new();
        let mut raw_line = String::new();
        let mut line_no = 0usize;

        loop {
            let n = crate::io_util::read_bounded_line(
                &mut reader,
                &mut raw_line,
                MAX_HNEX_LINE_BYTES,
                ".hnex",
            )?;
            if n == 0 {
                break;
            }
            line_no += 1;
            let line = raw_line.trim_end_matches(['\r', '\n']);
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            match fields[0] {
                "HNEX" => {
                    anyhow::ensure!(fields.len() == 2, "malformed HNEX header on line {line_no}");
                    let vn = fields[1]
                        .strip_prefix("VN:")
                        .with_context(|| format!("invalid HNEX header on line {line_no}"))?;
                    version = Some(vn.to_owned());
                }
                "TREE" => {
                    anyhow::ensure!(fields.len() == 2, "malformed TREE on line {line_no}");
                    guide_tree_newick = Some(fields[1].to_owned());
                }
                "SEQ" => {
                    anyhow::ensure!(fields.len() == 4, "malformed SEQ on line {line_no}");
                    sequences.push(HnexSeqMeta {
                        species: fields[1].to_owned(),
                        contig: fields[2].to_owned(),
                        length: fields[3].parse().context("invalid SEQ length")?,
                    });
                }
                "BLK" => {
                    anyhow::ensure!(fields.len() == 3, "malformed BLK on line {line_no}");
                    let score_bits: f64 = fields[2].parse().context("invalid BLK score")?;
                    anyhow::ensure!(
                        score_bits.is_finite(),
                        "non-finite BLK score on line {line_no}"
                    );
                    blocks.push(HnexBlock {
                        block_id: fields[1].parse().context("invalid BLK id")?,
                        score_bits,
                        rows: Vec::new(),
                        col_posteriors: Vec::new(),
                    });
                }
                "POST" => {
                    let cur = blocks
                        .last_mut()
                        .with_context(|| format!("POST outside BLK on line {line_no}"))?;
                    anyhow::ensure!(fields.len() == 2, "malformed POST on line {line_no}");
                    if !fields[1].is_empty() {
                        cur.col_posteriors = fields[1]
                            .split(',')
                            .map(|tok| {
                                let p = tok.parse::<f32>().context("invalid POST probability")?;
                                anyhow::ensure!(
                                    p.is_finite() && (0.0..=1.0).contains(&p),
                                    "POST probability out of [0,1] on line {line_no}: {p}"
                                );
                                Ok(p)
                            })
                            .collect::<Result<Vec<_>>>()?;
                    }
                }
                "ROW" => {
                    let cur = blocks
                        .last_mut()
                        .with_context(|| format!("ROW outside BLK on line {line_no}"))?;
                    anyhow::ensure!(fields.len() == 8, "malformed ROW on line {line_no}");
                    let strand = match fields[5] {
                        "+" => '+',
                        "-" => '-',
                        other => bail!("invalid ROW strand `{other}` on line {line_no}"),
                    };
                    let start: usize = fields[3].parse().context("invalid ROW start")?;
                    let size: usize = fields[4].parse().context("invalid ROW size")?;
                    let src_size: usize = fields[6].parse().context("invalid ROW src_size")?;
                    anyhow::ensure!(
                        start.checked_add(size).is_some_and(|e| e <= src_size),
                        "ROW interval [{start}, {start}+{size}) exceeds src_size {src_size} on line {line_no}"
                    );
                    let seq_str = fields[7];
                    anyhow::ensure!(
                        seq_str.is_ascii(),
                        "non-ASCII ROW sequence on line {line_no}"
                    );
                    let ungapped = seq_str.bytes().filter(|&b| b != b'-').count();
                    anyhow::ensure!(
                        ungapped == size,
                        "ROW ungapped base count ({ungapped}) != declared size ({size}) on line {line_no}"
                    );
                    if let Some(first) = cur.rows.first() {
                        anyhow::ensure!(
                            seq_str.len() == first.seq.len(),
                            "ROW column width mismatch ({} != {}) on line {line_no}",
                            seq_str.len(),
                            first.seq.len()
                        );
                    }
                    cur.rows.push(HnexRow {
                        species: fields[1].to_owned(),
                        contig: fields[2].to_owned(),
                        start,
                        size,
                        strand,
                        src_size,
                        seq: seq_str.to_owned(),
                    });
                }
                "BUB" => {
                    anyhow::ensure!(fields.len() == 6, "malformed BUB on line {line_no}");
                    let parse_opt_id = |s: &str| -> Result<Option<usize>> {
                        if s == "-" {
                            Ok(None)
                        } else {
                            Ok(Some(s.parse().context("invalid BUB neighbor block id")?))
                        }
                    };
                    let delta_l_bits: f64 =
                        fields[4].parse().context("invalid BUB delta_l_bits")?;
                    anyhow::ensure!(
                        delta_l_bits.is_finite(),
                        "non-finite BUB delta_l_bits on line {line_no}"
                    );
                    bubbles.push(HnexBubble {
                        bubble_id: fields[1].parse().context("invalid BUB id")?,
                        left_block: parse_opt_id(fields[2])?,
                        right_block: parse_opt_id(fields[3])?,
                        delta_l_bits,
                        reason: fields[5].to_owned(),
                        spans: Vec::new(),
                    });
                }
                "SPAN" => {
                    let cur = bubbles
                        .last_mut()
                        .with_context(|| format!("SPAN outside BUB on line {line_no}"))?;
                    anyhow::ensure!(fields.len() == 6, "malformed SPAN on line {line_no}");
                    let start: usize = fields[3].parse().context("invalid SPAN start")?;
                    let end: usize = fields[4].parse().context("invalid SPAN end")?;
                    anyhow::ensure!(
                        start <= end,
                        "invalid SPAN interval [{start}, {end}) on line {line_no}"
                    );
                    cur.spans.push(HnexBubbleSpan {
                        species: fields[1].to_owned(),
                        contig: fields[2].to_owned(),
                        start,
                        end,
                        seq: fields[5].to_owned(),
                    });
                }
                other => bail!("unknown .hnex record `{other}` on line {line_no}"),
            }
        }
        let version = version.context("missing HNEX header record")?;
        for blk in &blocks {
            if !blk.col_posteriors.is_empty()
                && let Some(first) = blk.rows.first()
            {
                anyhow::ensure!(
                    blk.col_posteriors.len() == first.seq.len(),
                    "BLK {} POST length ({}) != column width ({})",
                    blk.block_id,
                    blk.col_posteriors.len(),
                    first.seq.len()
                );
            }
        }
        Ok(Self {
            version,
            guide_tree_newick,
            sequences,
            blocks,
            bubbles,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hnex_rejects_malformed_rows_strands_bounds_and_control_chars() {
        // Mismatched row lengths inside a block.
        let bad_width = "HNEX\tVN:1.0\nBLK\t0\t10.0\nROW\tspA\tchr1\t0\t4\t+\t100\tACGT\nROW\tspB\tchr1\t0\t2\t+\t100\tAC\n";
        assert!(HnexArchive::read_from(bad_width.as_bytes()).is_err());

        // Invalid strand character.
        let bad_strand = "HNEX\tVN:1.0\nBLK\t0\t10.0\nROW\tspA\tchr1\t0\t4\t?\t100\tACGT\n";
        assert!(HnexArchive::read_from(bad_strand.as_bytes()).is_err());

        // start + size > src_size.
        let bad_bounds = "HNEX\tVN:1.0\nBLK\t0\t10.0\nROW\tspA\tchr1\t98\t4\t+\t100\tACGT\n";
        assert!(HnexArchive::read_from(bad_bounds.as_bytes()).is_err());

        // Tab in reason rejected by write_to.
        let bad_archive = HnexArchive {
            version: "1.0".to_owned(),
            guide_tree_newick: None,
            sequences: vec![HnexSeqMeta {
                species: "spC".to_owned(),
                contig: "chr9".to_owned(),
                length: 4,
            }],
            blocks: Vec::new(),
            bubbles: vec![HnexBubble {
                bubble_id: 0,
                left_block: None,
                right_block: None,
                delta_l_bits: -1.0,
                reason: "bad\treason".to_owned(),
                spans: vec![HnexBubbleSpan {
                    species: "spC".to_owned(),
                    contig: "chr9".to_owned(),
                    start: 0,
                    end: 4,
                    seq: "ACGT".to_owned(),
                }],
            }],
        };
        let mut out = Vec::new();
        assert!(bad_archive.write_to(&mut out).is_err());
    }
}
