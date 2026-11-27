//! Streaming MAF v1 and GFA v1.1 output writers for reference-free TBA blocksets.

use super::TbaBlock;
use crate::formats::strip_sample_prefix;
use std::collections::{BTreeMap, HashSet};
use std::fmt::Write as _;
use std::io::Write;

/// Compact `Copy` segment identifier replacing per-visit heap `String` allocations:
/// `sub_plus_one == 0` formats as `tba_{blk}`, `sub_plus_one == s + 1` formats as `tba_{blk}_{s}`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct SegId {
    blk: usize,
    sub_plus_one: usize,
}

impl std::fmt::Display for SegId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.sub_plus_one == 0 {
            write!(f, "tba_{}", self.blk)
        } else {
            write!(f, "tba_{}_{}", self.blk, self.sub_plus_one - 1)
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct WalkVisit {
    fwd_start: usize,
    fwd_end: usize,
    seg_id: SegId,
    orient: char,
}

/// Write a single multi-species [`TbaBlock`] stanza to a MAF v1 writer.
pub(super) fn write_maf_block<W: Write>(
    idx: usize,
    blk: &TbaBlock,
    maf_out: &mut W,
) -> std::io::Result<()> {
    writeln!(maf_out, "a score={:.1}", blk.score)?;
    for row in &blk.rows {
        let stripped = strip_sample_prefix(&row.contig, &row.species);
        let seq_str = std::str::from_utf8(&row.seq).map_err(|e| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("non-UTF-8 sequence in TBA block {idx}: {e}"),
            )
        })?;
        writeln!(
            maf_out,
            "s {}.{stripped} {} {} {} {} {seq_str}",
            row.species, row.start, row.size, row.strand, row.src_size
        )?;
    }
    maf_out.write_all(b"\n")
}

/// Write a reference-free TBA blockset to MAF v1 and GFA v1.1 writers.
///
/// Gapped TBA blocks are decomposed into gapless sub-segments (`tba_{idx}_{sub}`) by column
/// presence mask so every `S` record's `LN:i` equals the exact ungapped base count traversed by
/// every visiting row. `W` walk records are split across unaligned genomic gaps so every emitted
/// walk satisfies `sum(LN(s) for s in walk) == walk_end - walk_start` (GFA v1.1 / `horizon_lift` /
/// `odgi` / `vg` invariant), and `L` links are emitted only between genuinely contiguous steps.
///
/// Returns `(total_blocks, multi_species_blocks)` where `multi_species_blocks` counts blocks
/// with `>= 3` species (matching `horizon_merge`'s `>= 2` query rows convention).
///
/// # Errors
/// Returns an `std::io::Error` if writing to `maf_out` or `gfa_out` fails.
#[allow(clippy::too_many_lines)]
pub fn write_tba_maf_and_gfa<W1: Write, W2: Write>(
    blocks: &[TbaBlock],
    maf_out: &mut W1,
    gfa_out: &mut W2,
) -> std::io::Result<(usize, usize)> {
    maf_out.write_all(b"##maf version=1 scoring=horizon_tba\n")?;
    gfa_out.write_all(b"H\tVN:Z:1.1\n")?;

    let mut total_blocks = 0usize;
    let mut multi_species_blocks = 0usize;
    let mut walks: BTreeMap<(String, String), Vec<WalkVisit>> = BTreeMap::new();
    let mut row_offsets: Vec<usize> = Vec::new();
    let mut mask: Vec<bool> = Vec::new();

    for (idx, blk) in blocks.iter().enumerate() {
        let n_cols = blk.n_cols();
        if blk.rows.len() < 2 || n_cols == 0 {
            continue;
        }
        total_blocks += 1;
        if blk.rows.len() >= 3 {
            multi_species_blocks += 1;
        }

        write_maf_block(idx, blk, maf_out)?;

        let has_internal_gaps = blk
            .rows
            .iter()
            .any(|r| r.seq.len() != n_cols || r.seq.contains(&b'-'));
        if !has_internal_gaps {
            let seg_id = SegId {
                blk: idx,
                sub_plus_one: 0,
            };
            let anchor = &blk.rows[0];
            let anchor_ctg = strip_sample_prefix(&anchor.contig, &anchor.species);
            let (fwd_s, _) = anchor.fwd_interval();
            writeln!(
                gfa_out,
                "S\t{seg_id}\t*\tLN:i:{n_cols}\tSP:i:{}\tSN:Z:{}.{anchor_ctg}\tSO:i:{fwd_s}\tSR:i:0",
                blk.rows.len(),
                anchor.species,
            )?;
            for row in &blk.rows {
                let stripped = strip_sample_prefix(&row.contig, &row.species);
                let (fwd_s, fwd_e) = row.fwd_interval();
                let orient = if row.strand == '-' { '-' } else { '+' };
                walks
                    .entry((row.species.clone(), stripped.to_owned()))
                    .or_default()
                    .push(WalkVisit {
                        fwd_start: fwd_s,
                        fwd_end: fwd_e,
                        seg_id,
                        orient,
                    });
            }
            continue;
        }

        // Decompose gapped TBA blocks into maximal gapless sub-segments by per-column presence mask
        // so every S record has exact ungapped length LN:i for all visiting rows.
        row_offsets.clear();
        row_offsets.resize(blk.rows.len(), 0);
        mask.clear();
        mask.resize(blk.rows.len(), false);
        let mut col = 0usize;
        let mut sub_idx = 0usize;
        while col < n_cols {
            for (slot, r) in mask.iter_mut().zip(&blk.rows) {
                *slot = r.seq.get(col).is_some_and(|&b| b != b'-');
            }
            let mut end_col = col + 1;
            while end_col < n_cols
                && mask
                    .iter()
                    .zip(&blk.rows)
                    .all(|(&present, r)| r.seq.get(end_col).is_some_and(|&b| b != b'-') == present)
            {
                end_col += 1;
            }
            let sub_len = end_col - col;
            let active_count = mask.iter().filter(|&&present| present).count();
            if active_count > 0 {
                let seg_id = SegId {
                    blk: idx,
                    sub_plus_one: sub_idx + 1,
                };
                sub_idx += 1;
                let first_active = mask.iter().position(|&p| p).unwrap_or(0);
                let anchor = &blk.rows[first_active];
                let anchor_ctg = strip_sample_prefix(&anchor.contig, &anchor.species);
                let (a_fwd_s, a_fwd_e) = anchor.fwd_interval();
                let a_so = if anchor.strand == '-' {
                    a_fwd_e.saturating_sub(row_offsets[first_active] + sub_len)
                } else {
                    a_fwd_s + row_offsets[first_active]
                };
                writeln!(
                    gfa_out,
                    "S\t{seg_id}\t*\tLN:i:{sub_len}\tSP:i:{active_count}\tSN:Z:{}.{anchor_ctg}\tSO:i:{a_so}\tSR:i:0",
                    anchor.species,
                )?;
                for (ri, row) in blk.rows.iter().enumerate() {
                    if !mask[ri] {
                        continue;
                    }
                    let stripped = strip_sample_prefix(&row.contig, &row.species);
                    let (fwd_s, fwd_e) = row.fwd_interval();
                    let orient = if row.strand == '-' { '-' } else { '+' };
                    let (sub_s, sub_e) = if row.strand == '-' {
                        let e = fwd_e.saturating_sub(row_offsets[ri]);
                        let s = e.saturating_sub(sub_len);
                        (s, e)
                    } else {
                        let s = fwd_s + row_offsets[ri];
                        (s, s + sub_len)
                    };
                    row_offsets[ri] += sub_len;
                    walks
                        .entry((row.species.clone(), stripped.to_owned()))
                        .or_default()
                        .push(WalkVisit {
                            fwd_start: sub_s,
                            fwd_end: sub_e,
                            seg_id,
                            orient,
                        });
                }
            }
            col = end_col;
        }
    }

    // Sort each (species, contig) visit list by forward coordinate and split into maximal
    // contiguous runs (next.fwd_start == prev.fwd_end) so L links and W walk spans remain exact.
    let mut emitted_links: HashSet<(SegId, char, SegId, char)> = HashSet::new();
    let mut step_buf = String::new();
    for ((species, contig), mut visits) in walks {
        visits.sort_unstable_by_key(|v| (v.fwd_start, v.fwd_end, v.seg_id));

        let mut run_start = 0usize;
        while run_start < visits.len() {
            let mut run_end = run_start + 1;
            while run_end < visits.len() && visits[run_end].fwd_start == visits[run_end - 1].fwd_end
            {
                run_end += 1;
            }
            let run = &visits[run_start..run_end];
            let w_start = run[0].fwd_start;
            let w_end = run[run.len() - 1].fwd_end;

            for win in run.windows(2) {
                let s1 = win[0].seg_id;
                let o1 = win[0].orient;
                let s2 = win[1].seg_id;
                let o2 = win[1].orient;
                let flip = |o: char| if o == '-' { '+' } else { '-' };
                let fwd_key = (s1, o1, s2, o2);
                let rev_key = (s2, flip(o2), s1, flip(o1));
                let canon = fwd_key.min(rev_key);
                if emitted_links.insert(canon) {
                    writeln!(gfa_out, "L\t{s1}\t{o1}\t{s2}\t{o2}\t0M")?;
                }
            }

            step_buf.clear();
            for v in run {
                let dir = if v.orient == '-' { '<' } else { '>' };
                let _ = write!(step_buf, "{dir}{}", v.seg_id);
            }
            writeln!(
                gfa_out,
                "W\t{species}\t0\t{contig}\t{w_start}\t{w_end}\t{step_buf}"
            )?;
            run_start = run_end;
        }
    }

    Ok((total_blocks, multi_species_blocks))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tba::TbaRow;

    fn mk_row(species: &str, contig: &str, start: usize, seq: &[u8]) -> TbaRow {
        let size = seq.iter().filter(|&&b| b != b'-').count();
        TbaRow {
            species: species.to_owned(),
            contig: contig.to_owned(),
            start,
            size,
            strand: '+',
            src_size: 1000,
            seq: seq.to_vec(),
        }
    }

    #[test]
    fn write_tba_maf_and_gfa_emits_valid_maf_gfa_and_exact_walk_spans() {
        let blocks = vec![
            // 3-species gapped block on sp1.chr1 -> decomposes into 3 gapless sub-segments (4, 1, 3 bp)
            TbaBlock {
                score: 42.5,
                rows: vec![
                    mk_row("sp1", "sp1.chr1", 10, b"ACGTACGT"),
                    mk_row("sp2", "chr1", 20, b"ACGTACGT"),
                    mk_row("sp3", "chr1", 30, b"ACGT-CGT"),
                ],
            },
            // Adjacent gapless block immediately following sp1:[10..18) -> contiguous W walk
            TbaBlock {
                score: 16.0,
                rows: vec![
                    mk_row("sp1", "chr1", 18, b"TTTTGGGG"),
                    mk_row("sp2", "chr1", 60, b"TTTTGGGG"),
                ],
            },
            // Clade-specific block on a different anchor (sp2, chr2) -> separate W walk, no chimeric link
            TbaBlock {
                score: 12.0,
                rows: vec![
                    mk_row("sp2", "chr2", 5, b"CCCCAAAA"),
                    mk_row("sp3", "chr2", 5, b"CCCCAAAA"),
                ],
            },
            // Degenerate 1-row block -> skipped
            TbaBlock {
                score: 0.0,
                rows: vec![mk_row("sp1", "chr1", 90, b"ACGT")],
            },
        ];

        let mut maf_buf = Vec::new();
        let mut gfa_buf = Vec::new();
        let (total, multi) =
            write_tba_maf_and_gfa(&blocks, &mut maf_buf, &mut gfa_buf).expect("write");
        assert_eq!(total, 3);
        assert_eq!(multi, 1);

        let maf_text = String::from_utf8(maf_buf).expect("utf8 maf");
        assert!(maf_text.starts_with("##maf version=1 scoring=horizon_tba\n"));
        assert!(maf_text.contains("s sp1.chr1 10 8 + 1000 ACGTACGT\n"));
        assert!(maf_text.contains("s sp3.chr1 30 7 + 1000 ACGT-CGT\n"));
        assert!(!maf_text.contains("s sp1.chr1 90"));

        let gfa_text = String::from_utf8(gfa_buf).expect("utf8 gfa");
        assert!(gfa_text.starts_with("H\tVN:Z:1.1\n"));
        // Gapped block 0 splits into gapless sub-segments: 4 bp, 1 bp indel bubble, 3 bp.
        assert!(gfa_text.contains("S\ttba_0_0\t*\tLN:i:4\tSP:i:3\t"));
        assert!(gfa_text.contains("S\ttba_0_1\t*\tLN:i:1\tSP:i:2\t"));
        assert!(gfa_text.contains("S\ttba_0_2\t*\tLN:i:3\tSP:i:3\t"));
        assert!(gfa_text.contains("S\ttba_1\t*\tLN:i:8\tSP:i:2\t"));
        // sp1 is contiguous across [10..18) + [18..26) -> single walk with sum(LN) = 4+1+3+8 = 16 == 26-10
        assert!(gfa_text.contains("L\ttba_0_2\t+\ttba_1\t+\t0M\n"));
        assert!(gfa_text.contains("W\tsp1\t0\tchr1\t10\t26\t>tba_0_0>tba_0_1>tba_0_2>tba_1\n"));
        // sp2 has an unaligned gap between [20..28) and [60..68) -> two separate exact W walks
        assert!(gfa_text.contains("W\tsp2\t0\tchr1\t20\t28\t>tba_0_0>tba_0_1>tba_0_2\n"));
        assert!(gfa_text.contains("W\tsp2\t0\tchr1\t60\t68\t>tba_1\n"));
        // sp3 skips the 1 bp insertion bubble -> sum(LN) = 4+3 = 7 == 37-30
        assert!(gfa_text.contains("W\tsp3\t0\tchr1\t30\t37\t>tba_0_0>tba_0_2\n"));
        assert!(gfa_text.contains("W\tsp2\t0\tchr2\t5\t13\t>tba_2\n"));
        assert!(gfa_text.contains("W\tsp3\t0\tchr2\t5\t13\t>tba_2\n"));
        assert!(
            !gfa_text.contains("L\ttba_1\t+\ttba_2"),
            "must not emit chimeric link across different anchor contigs"
        );
    }
}
