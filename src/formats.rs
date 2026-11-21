//! GFA v1.1, MAF v1, FASTA I/O, and e-graph paralog resolution helpers.

use crate::egraph;
use crate::mdl::{ContigAlignment, RunAlignment};
use crate::netting::{partition_colinear_blocks, partition_colinear_run_groups};
use crate::seeding::complement_base;
use std::fmt::Write as _;

/// Read peak resident memory (`VmHWM` in KiB) from `/proc/self/status` on Linux.
#[must_use]
pub fn read_vm_hwm_kb() -> Option<usize> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in status.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            return rest.split_whitespace().next()?.parse::<usize>().ok();
        }
    }
    None
}

/// Check whether `name` equals `prefix` or starts with `prefix` followed by `.`, `#`, or `:`.
#[must_use]
pub fn matches_sample_prefix(name: &str, prefix: &str) -> bool {
    name == prefix
        || name.strip_prefix(prefix).is_some_and(|rest| {
            rest.starts_with('.') || rest.starts_with('#') || rest.starts_with(':')
        })
}

/// Strip a leading `prefix.` / `prefix#` / `prefix:` sample qualifier from `name` when present.
#[must_use]
pub fn strip_sample_prefix<'a>(name: &'a str, prefix: &str) -> &'a str {
    if let Some(rest) = name.strip_prefix(prefix)
        && (rest.starts_with('.') || rest.starts_with('#') || rest.starts_with(':'))
    {
        &rest[1..]
    } else {
        name
    }
}

/// Format aligned LCBs and Horizon Bubbles as minimal GFA v1.1 records.
#[must_use]
pub fn format_gfa_summary(name_x: &str, name_y: &str, aln: &ContigAlignment) -> String {
    let mut out = String::from("H\tVN:Z:1.1\n");
    let total_bubbles = aln.macro_bubbles + aln.internal_bubbles;
    let _ = writeln!(
        out,
        "S\t{name_x}__{name_y}\t*\tLN:i:{}\tLF:i:{}\tLR:i:{}\tHB:i:{total_bubbles}",
        aln.aligned_pairs.len(),
        aln.lcbs_fwd,
        aln.lcbs_rc
    );
    out
}

fn emit_gfa_blocks(
    name_x: &str,
    name_y: &str,
    runs: &[(usize, usize, usize, usize, char)],
) -> String {
    let mut out = String::from("H\tVN:Z:1.1\n");
    for (idx, &(bx, by, blen, _x_end, bst)) in runs.iter().enumerate() {
        let seg_id = format!("{name_x}__{name_y}__blk_{idx}");
        let _ = writeln!(
            out,
            "S\t{seg_id}\t*\tLN:i:{blen}\tXS:i:{bx}\tYS:i:{by}\tST:Z:{bst}"
        );
        if idx > 0 {
            let prev_seg = format!("{name_x}__{name_y}__blk_{}", idx - 1);
            let prev_bst = runs[idx - 1].4;
            let gap_x = bx.saturating_sub(runs[idx - 1].3);
            if gap_x >= 6 {
                let bub_id = format!("{name_x}__{name_y}__bub_{idx}");
                let _ = writeln!(out, "S\t{bub_id}\t*\tLN:i:{gap_x}\tHB:i:1");
                let _ = writeln!(out, "L\t{prev_seg}\t{prev_bst}\t{bub_id}\t+\t0M");
                let _ = writeln!(out, "L\t{bub_id}\t+\t{seg_id}\t{bst}\t0M");
            } else {
                let _ = writeln!(out, "L\t{prev_seg}\t{prev_bst}\t{seg_id}\t{bst}\t0M");
            }
        }
    }
    out
}

/// Emit detailed GFA v1.1 block segments, Horizon-Bubble nodes, and links for a contig pair.
#[must_use]
pub fn format_gfa_blocks(name_x: &str, name_y: &str, aln: &ContigAlignment) -> String {
    if aln.aligned_pairs.is_empty() {
        return String::from("H\tVN:Z:1.1\n");
    }

    let runs: Vec<(usize, usize, usize, usize, char)> =
        partition_colinear_blocks(&aln.aligned_pairs)
            .into_iter()
            .map(|(blk, is_rev)| {
                let strand = if is_rev { '-' } else { '+' };
                let x_end = blk.last().map_or(blk[0].0, |p| p.0 + 1);
                (blk[0].0, blk[0].1, blk.len(), x_end, strand)
            })
            .collect();

    emit_gfa_blocks(name_x, name_y, &runs)
}

/// Run e-graph duplication/inversion collapse on a sequence of signed block labels.
#[must_use]
pub fn resolve_lcb_duplications(blocks: &[(&str, i8)]) -> Option<egraph::Expr> {
    let mut eg = egraph::EGraph::new();
    let root = egraph::build_seq_expr(&mut eg, blocks)?;
    egraph::saturate_synteny(&mut eg, 4, 800);
    egraph::extract_min_cost(&eg, root).map(|(_, expr)| expr)
}

pub(crate) fn collect_windowed_duplications(
    tokens: &[(String, i8)],
    out: &mut std::collections::BTreeSet<String>,
) {
    // Use 50% overlapping 12-token windows within each contig so duplications
    // straddling a window boundary are still detected without mixing contigs.
    let mut start = 0usize;
    while start < tokens.len() {
        let end = (start + 12).min(tokens.len());
        let win = &tokens[start..end];
        if win.len() >= 2 {
            let refs: Vec<(&str, i8)> = win.iter().map(|(s, st)| (s.as_str(), *st)).collect();
            if let Some(expr) = resolve_lcb_duplications(&refs) {
                for d in egraph::collect_dups(&expr) {
                    out.insert(d);
                }
            }
        }
        if end == tokens.len() {
            break;
        }
        start += 6;
    }
}

/// Format e-graph-resolved paralog annotations as GFA v1.1 comment/tag lines.
#[must_use]
pub fn format_gfa_paralog_tags(name_x: &str, name_y: &str, dups: &[String]) -> String {
    if dups.is_empty() {
        return String::new();
    }
    format!("# EG:Z:{name_x}__{name_y}\t{}\n", dups.join(","))
}

/// Emit a compact GFA v1.1 summary segment for a [`RunAlignment`].
#[must_use]
pub fn format_gfa_summary_runs(name_x: &str, name_y: &str, aln: &RunAlignment) -> String {
    let mut out = String::from("H\tVN:Z:1.1\n");
    let total_bubbles = aln.macro_bubbles + aln.internal_bubbles;
    let _ = writeln!(
        out,
        "S\t{name_x}__{name_y}\t*\tLN:i:{}\tLF:i:{}\tLR:i:{}\tHB:i:{total_bubbles}",
        aln.total_pairs, aln.lcbs_fwd, aln.lcbs_rc
    );
    out
}

/// Emit detailed GFA v1.1 block segments, Horizon-Bubble nodes, and links for a [`RunAlignment`].
#[must_use]
pub fn format_gfa_blocks_runs(name_x: &str, name_y: &str, aln: &RunAlignment) -> String {
    if aln.aligned_runs.is_empty() {
        return String::from("H\tVN:Z:1.1\n");
    }
    let blocks = partition_colinear_run_groups(&aln.aligned_runs);
    let runs: Vec<(usize, usize, usize, usize, char)> = blocks
        .into_iter()
        .map(|(grp, is_rev)| {
            let strand = if is_rev { '-' } else { '+' };
            let first = grp[0];
            let last = grp[grp.len() - 1];
            let total_bp: usize = grp.iter().map(|r| r.len).sum();
            let x_end = last.x_start + last.len;
            (first.x_start, first.y_start, total_bp, x_end, strand)
        })
        .collect();
    emit_gfa_blocks(name_x, name_y, &runs)
}

/// Bucket colinear alignment run-groups into `bin_bp` reference windows and return
/// deduplicated block labels collapsed by the e-graph duplication/inversion rules.
#[must_use]
pub fn resolve_alignment_duplications_runs(aln: &RunAlignment, bin_bp: usize) -> Vec<String> {
    let bins: Vec<(String, i8)> = partition_colinear_run_groups(&aln.aligned_runs)
        .into_iter()
        .map(|(grp, is_rev)| {
            let b_idx = grp[0].x_start / bin_bp.max(1);
            (format!("B{b_idx}"), if is_rev { -1 } else { 1 })
        })
        .collect();
    let mut dups = std::collections::BTreeSet::new();
    collect_windowed_duplications(&bins, &mut dups);
    dups.into_iter().collect()
}

/// Group run-encoded alignments into colinear MAF blocks and emit standard MAF v1 text.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn format_maf_blocks_runs(
    name_x: &str,
    len_x: usize,
    name_y: &str,
    len_y: usize,
    seq_x: &[u8],
    seq_y: &[u8],
    aln: &RunAlignment,
) -> String {
    let mut out = String::from("##maf version=1 scoring=horizon_mdl\n");
    if aln.aligned_runs.is_empty() {
        return out;
    }
    for (grp, is_rev) in partition_colinear_run_groups(&aln.aligned_runs) {
        let valid = grp.iter().all(|r| {
            r.len > 0
                && r.x_start + r.len <= seq_x.len()
                && if r.is_rc {
                    r.y_start < seq_y.len() && r.y_start + 1 >= r.len
                } else {
                    r.y_start + r.len <= seq_y.len()
                }
        });
        if !valid {
            continue;
        }
        let x0 = grp[0].x_start;
        let last = grp[grp.len() - 1];
        let x_end = last.x_start + last.len;
        let score_bp: usize = grp.iter().map(|r| r.len).sum();
        let mut ax = String::with_capacity(x_end - x0 + 16);
        let mut ay = String::with_capacity(x_end - x0 + 16);

        let mut prev_x_end = x0;
        let mut prev_y_next = grp[0].y_start;
        for (ri, r) in grp.iter().enumerate() {
            if ri > 0 {
                for &bx in &seq_x[prev_x_end..r.x_start] {
                    ax.push(bx.to_ascii_uppercase() as char);
                    ay.push('-');
                }
                if is_rev {
                    for &by in seq_y[(r.y_start + 1)..prev_y_next].iter().rev() {
                        ax.push('-');
                        ay.push(complement_base(by.to_ascii_uppercase()) as char);
                    }
                } else {
                    for &by in &seq_y[prev_y_next..r.y_start] {
                        ax.push('-');
                        ay.push(by.to_ascii_uppercase() as char);
                    }
                }
            }
            for &bx in &seq_x[r.x_start..(r.x_start + r.len)] {
                ax.push(bx.to_ascii_uppercase() as char);
            }
            if is_rev {
                let y_lo = r.y_start + 1 - r.len;
                for &by in seq_y[y_lo..=r.y_start].iter().rev() {
                    ay.push(complement_base(by.to_ascii_uppercase()) as char);
                }
                prev_y_next = y_lo;
            } else {
                let y_hi = r.y_start + r.len;
                for &by in &seq_y[r.y_start..y_hi] {
                    ay.push(by.to_ascii_uppercase() as char);
                }
                prev_y_next = y_hi;
            }
            prev_x_end = r.x_start + r.len;
        }

        let (y_start, y_span, y_strand) = if is_rev {
            let y_lo = last.y_start + 1 - last.len;
            let y_hi = grp[0].y_start + 1;
            (len_y.saturating_sub(y_hi), y_hi - y_lo, '-')
        } else {
            let y_lo = grp[0].y_start;
            let y_hi = last.y_start + last.len;
            (y_lo, y_hi - y_lo, '+')
        };

        let _ = writeln!(out, "a score={score_bp}.0");
        let _ = writeln!(out, "s {name_x} {x0} {} + {len_x} {ax}", x_end - x0);
        let _ = writeln!(
            out,
            "s {name_y} {y_start} {y_span} {y_strand} {len_y} {ay}\n"
        );
    }
    out
}

/// Group aligned `(x, y)` pairs into colinear MAF blocks and emit standard MAF v1 text.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn format_maf_blocks(
    name_x: &str,
    len_x: usize,
    name_y: &str,
    len_y: usize,
    seq_x: &[u8],
    seq_y: &[u8],
    aln: &ContigAlignment,
) -> String {
    let mut out = String::from("##maf version=1 scoring=horizon_mdl\n");
    if aln.aligned_pairs.is_empty() {
        return out;
    }

    for (blk, is_rev) in partition_colinear_blocks(&aln.aligned_pairs) {
        if blk
            .iter()
            .any(|&(x, y)| x >= seq_x.len() || y >= seq_y.len())
        {
            continue;
        }
        let x0 = blk[0].0;
        let x_end = blk[blk.len() - 1].0 + 1;
        let mut ax = String::with_capacity(blk.len() + 8);
        let mut ay = String::with_capacity(blk.len() + 8);

        let mut prev_x = blk[0].0;
        let mut prev_y = blk[0].1;
        ax.push(seq_x[prev_x].to_ascii_uppercase() as char);
        let first_y = seq_y[prev_y].to_ascii_uppercase();
        ay.push(if is_rev {
            complement_base(first_y) as char
        } else {
            first_y as char
        });

        for &(x, y) in blk.iter().skip(1) {
            for &bx in &seq_x[(prev_x + 1)..x] {
                ax.push(bx.to_ascii_uppercase() as char);
                ay.push('-');
            }
            if is_rev {
                for &by in seq_y[(y + 1)..prev_y].iter().rev() {
                    ax.push('-');
                    ay.push(complement_base(by.to_ascii_uppercase()) as char);
                }
            } else {
                for &by in &seq_y[(prev_y + 1)..y] {
                    ax.push('-');
                    ay.push(by.to_ascii_uppercase() as char);
                }
            }
            ax.push(seq_x[x].to_ascii_uppercase() as char);
            let by = seq_y[y].to_ascii_uppercase();
            ay.push(if is_rev {
                complement_base(by) as char
            } else {
                by as char
            });
            prev_x = x;
            prev_y = y;
        }

        let (y_start, y_span, y_strand) = if is_rev {
            let y_lo = blk[blk.len() - 1].1;
            let y_hi = blk[0].1 + 1;
            (len_y.saturating_sub(y_hi), y_hi - y_lo, '-')
        } else {
            let y_lo = blk[0].1;
            let y_hi = blk[blk.len() - 1].1 + 1;
            (y_lo, y_hi - y_lo, '+')
        };

        let _ = writeln!(out, "a score={}.0", blk.len());
        let _ = writeln!(out, "s {name_x} {x0} {} + {len_x} {ax}", x_end - x0);
        let _ = writeln!(
            out,
            "s {name_y} {y_start} {y_span} {y_strand} {len_y} {ay}\n"
        );
    }
    out
}

/// Read up to `max_contigs` FASTA records `(contig_name, sequence_bytes)` from `path`
/// (`max_contigs == 0` reads all records).
///
/// # Errors
/// Returns an `std::io::Error` if `path` cannot be opened or read.
pub fn read_fasta_contigs(
    path: &std::path::Path,
    max_contigs: usize,
) -> std::io::Result<Vec<(String, Vec<u8>)>> {
    use std::io::BufRead as _;
    let f = std::fs::File::open(path)?;
    let mut reader = std::io::BufReader::with_capacity(256 * 1024, f);
    let mut out = Vec::new();
    let mut name = String::new();
    let mut seq = Vec::new();
    let mut line_buf = Vec::with_capacity(4096);
    loop {
        line_buf.clear();
        let n = reader.read_until(b'\n', &mut line_buf)?;
        if n == 0 {
            break;
        }
        let mut s = line_buf.as_slice();
        while s.last().is_some_and(u8::is_ascii_whitespace) {
            s = &s[..s.len() - 1];
        }
        while s.first().is_some_and(u8::is_ascii_whitespace) {
            s = &s[1..];
        }
        if let Some(rest) = s.strip_prefix(b">") {
            if !name.is_empty() {
                out.push((std::mem::take(&mut name), std::mem::take(&mut seq)));
                if max_contigs > 0 && out.len() >= max_contigs {
                    return Ok(out);
                }
            }
            let header = String::from_utf8_lossy(rest);
            header
                .split_whitespace()
                .next()
                .unwrap_or("contig")
                .clone_into(&mut name);
        } else if !name.is_empty() && !s.is_empty() {
            seq.extend_from_slice(s);
        }
    }
    if !name.is_empty() && (max_contigs == 0 || out.len() < max_contigs) {
        out.push((name, seq));
    }
    Ok(out)
}

/// A parsed `s`-line row from a MAF v1 alignment block.
#[derive(Debug, Clone)]
pub struct MafRow {
    /// Sequence/contig identifier.
    pub src: String,
    /// 0-based start offset on `strand`.
    pub start: usize,
    /// Strand (`'+'` or `'-'`).
    pub strand: char,
    /// Total ungapped length of the source sequence.
    pub src_size: usize,
    /// Gapped alignment row string.
    pub seq: String,
}

impl MafRow {
    /// Convert a 0-based ungapped index `nongap_k` along this row to a forward-strand coordinate.
    #[must_use]
    pub fn fwd_pos(&self, nongap_k: usize) -> usize {
        if self.strand == '-' {
            self.src_size
                .saturating_sub(1usize.saturating_add(self.start).saturating_add(nongap_k))
        } else {
            self.start.saturating_add(nongap_k)
        }
    }
}

/// Parse a single trimmed MAF `s` line into `(MafRow, ungapped_size)`.
#[must_use]
pub fn parse_maf_s_line(trimmed: &str) -> Option<(MafRow, usize)> {
    let rest = trimmed
        .strip_prefix('s')
        .filter(|r| r.starts_with(char::is_whitespace))?;
    let mut parts = rest.split_whitespace();
    let (
        Some(src),
        Some(start_str),
        Some(size_str),
        Some(strand_str),
        Some(src_size_str),
        Some(seq),
    ) = (
        parts.next(),
        parts.next(),
        parts.next(),
        parts.next(),
        parts.next(),
        parts.next(),
    )
    else {
        return None;
    };
    let (Ok(start), Ok(size), Ok(src_size)) = (
        start_str.parse::<usize>(),
        size_str.parse::<usize>(),
        src_size_str.parse::<usize>(),
    ) else {
        return None;
    };
    let strand = match strand_str {
        "+" => '+',
        "-" => '-',
        _ => return None,
    };
    Some((
        MafRow {
            src: src.to_owned(),
            start,
            strand,
            src_size,
            seq: seq.to_owned(),
        },
        size,
    ))
}

/// Parse two-row pairwise alignment blocks from a MAF v1 text stream.
#[must_use]
pub fn parse_two_row_maf_blocks(maf_text: &str) -> Vec<(MafRow, MafRow)> {
    let mut blocks = Vec::new();
    let mut rows: Vec<MafRow> = Vec::new();
    for line in maf_text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed == "a"
            || trimmed
                .strip_prefix('a')
                .is_some_and(|r| r.starts_with(char::is_whitespace))
        {
            if rows.len() >= 2 {
                let mut drain = rows.drain(..);
                if let (Some(r1), Some(r2)) = (drain.next(), drain.next()) {
                    blocks.push((r1, r2));
                }
            } else {
                rows.clear();
            }
            continue;
        }
        if let Some((row, _)) = parse_maf_s_line(trimmed) {
            rows.push(row);
        }
    }
    if rows.len() >= 2 {
        let mut drain = rows.drain(..);
        if let (Some(r1), Some(r2)) = (drain.next(), drain.next()) {
            blocks.push((r1, r2));
        }
    }
    blocks
}

#[cfg(test)]
mod tests {
    use super::read_fasta_contigs;

    #[test]
    fn reads_fasta_contigs_with_limits_and_blank_lines() {
        let fa_path =
            std::env::temp_dir().join(format!("horizon_formats_fa_test_{}.fa", std::process::id()));
        std::fs::write(&fa_path, ">chr1 desc\nACGT\n\nTGCA\n>   \nGGGG\n").expect("write fa");
        let all = read_fasta_contigs(&fa_path, 0).expect("read all");
        assert_eq!(
            all,
            vec![
                ("chr1".to_owned(), b"ACGTTGCA".to_vec()),
                ("contig".to_owned(), b"GGGG".to_vec()),
            ]
        );
        let one = read_fasta_contigs(&fa_path, 1).expect("read one");
        assert_eq!(one, vec![("chr1".to_owned(), b"ACGTTGCA".to_vec())]);
        let _ = std::fs::remove_file(&fa_path);
    }

    #[test]
    fn maf_row_fwd_pos_and_parse_s_line_handle_both_strands_and_errors() {
        let (fwd, len_f) = super::parse_maf_s_line("s hg38.chr1 10 4 + 100 ACGT").expect("fwd");
        assert_eq!(len_f, 4);
        assert_eq!(fwd.fwd_pos(0), 10);
        assert_eq!(fwd.fwd_pos(3), 13);

        let (rev, len_r) = super::parse_maf_s_line("s mm10.chr2 10 4 - 100 ACGT").expect("rev");
        assert_eq!(len_r, 4);
        // UCSC '-' start=10, size=4, src_size=100 -> forward [86, 90), so idx 0 -> 89, idx 3 -> 86.
        assert_eq!(rev.fwd_pos(0), 89);
        assert_eq!(rev.fwd_pos(3), 86);

        assert!(super::parse_maf_s_line("a score=1").is_none());
        assert!(super::parse_maf_s_line("s chr1 0 4 ? 100 ACGT").is_none());
        assert!(super::parse_maf_s_line("s chr1 bad 4 + 100 ACGT").is_none());
    }

    #[test]
    fn sample_prefix_helpers_and_two_row_maf_blocks_work() {
        assert!(super::matches_sample_prefix("hg38.chr1", "hg38"));
        assert!(super::matches_sample_prefix("hg38#0#chr1", "hg38"));
        assert!(!super::matches_sample_prefix("hg38_alt.chr1", "hg38"));
        assert_eq!(super::strip_sample_prefix("hg38#0#chr1", "hg38"), "0#chr1");
        assert_eq!(super::strip_sample_prefix("hg38#chr1", "hg38"), "chr1");
        assert_eq!(super::strip_sample_prefix("hg38.chr1", "hg38"), "chr1");

        let maf = "a score=1\ns hg38.chr1 0 4 + 10 ACGT\ns mm10.chr1 2 4 - 20 TGCA\n";
        let blocks = super::parse_two_row_maf_blocks(maf);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].0.src, "hg38.chr1");
        assert_eq!(blocks[0].1.src, "mm10.chr1");
    }
}
