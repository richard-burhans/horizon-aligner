//! Reference-guided multi-species MAF merging and coordinate projection.

use std::fmt::Write as _;

struct PairwiseMafBlock {
    ref_start: usize,
    ref_end: usize,
    ref_src_size: usize,
    qry_species: String,
    qry_contig: String,
    qry_start: usize,
    qry_strand: char,
    qry_src_size: usize,
    ref_seq: Vec<u8>,
    qry_seq: Vec<u8>,
}

struct ProjectedInterval {
    ref_bases: Vec<u8>,
    qry_bases: Vec<u8>,
    ins_after: Vec<Vec<u8>>,
    qry_start: usize,
    qry_size: usize,
    qry_src_size: usize,
}

#[derive(Clone, Copy)]
struct BlockCursor {
    col_idx: usize,
    ref_pos: usize,
    qry_consumed: usize,
}

fn project_block_interval_with_cursor(
    blk: &PairwiseMafBlock,
    cursor: &mut BlockCursor,
    lo: usize,
    hi: usize,
) -> Option<ProjectedInterval> {
    if blk.ref_start > lo || blk.ref_end < hi {
        return None;
    }
    if cursor.ref_pos > lo {
        *cursor = BlockCursor {
            col_idx: 0,
            ref_pos: blk.ref_start,
            qry_consumed: 0,
        };
    }
    let span = hi - lo;
    let mut col_idx = cursor.col_idx;
    let mut ref_pos = cursor.ref_pos;
    let mut qry_consumed = cursor.qry_consumed;
    let mut qry_in = 0usize;
    let mut ref_bases = Vec::with_capacity(span);
    let mut qry_bases = Vec::with_capacity(span);
    let mut ins_after = vec![Vec::new(); span];
    let mut saved_start_cursor: Option<BlockCursor> = None;

    while col_idx < blk.ref_seq.len() {
        let rb = blk.ref_seq[col_idx];
        let qb = blk.qry_seq[col_idx];
        if rb != b'-' {
            if ref_pos == lo && saved_start_cursor.is_none() {
                saved_start_cursor = Some(BlockCursor {
                    col_idx,
                    ref_pos,
                    qry_consumed,
                });
            }
            if (lo..hi).contains(&ref_pos) {
                ref_bases.push(rb);
                qry_bases.push(qb);
                if qb != b'-' {
                    qry_in += 1;
                }
            } else if ref_pos >= hi {
                break;
            }
            if ref_pos < lo && qb != b'-' {
                qry_consumed += 1;
            }
            ref_pos += 1;
        } else if ref_pos <= lo {
            if qb != b'-' {
                qry_consumed += 1;
            }
        } else {
            if ref_pos > hi {
                break;
            }
            let prev_idx = ref_pos - lo - 1;
            if qb != b'-' {
                ins_after[prev_idx].push(qb);
                qry_in += 1;
            }
        }
        col_idx += 1;
    }

    if let Some(sc) = saved_start_cursor {
        *cursor = BlockCursor {
            col_idx,
            ref_pos,
            qry_consumed: sc.qry_consumed + qry_in,
        };
    }

    if ref_bases.len() != span || qry_in == 0 {
        return None;
    }
    let qry_start = saved_start_cursor.map_or(qry_consumed, |sc| sc.qry_consumed);
    Some(ProjectedInterval {
        ref_bases,
        qry_bases,
        ins_after,
        qry_start: blk.qry_start.checked_add(qry_start)?,
        qry_size: qry_in,
        qry_src_size: blk.qry_src_size,
    })
}

#[cfg(test)]
fn project_block_interval(
    blk: &PairwiseMafBlock,
    lo: usize,
    hi: usize,
) -> Option<ProjectedInterval> {
    let mut cursor = BlockCursor {
        col_idx: 0,
        ref_pos: blk.ref_start,
        qry_consumed: 0,
    };
    project_block_interval_with_cursor(blk, &mut cursor, lo, hi)
}

struct ThreadedQueryRow {
    sp: String,
    qry_contig: String,
    qry_start: usize,
    qry_size: usize,
    qry_strand: char,
    qry_src_size: usize,
    seq: String,
}

fn thread_interval_columns(
    ref_bases: &[u8],
    sp_proj: &[(String, String, char, ProjectedInterval)],
) -> (String, Vec<ThreadedQueryRow>) {
    let span = ref_bases.len();
    let max_ins: Vec<usize> = (0..span)
        .map(|i| {
            sp_proj
                .iter()
                .map(|(_, _, _, p)| p.ins_after[i].len())
                .max()
                .unwrap_or(0)
        })
        .collect();
    let total_cols = span + max_ins.iter().sum::<usize>();
    let mut ref_out = String::with_capacity(total_cols);
    let mut q_outs: Vec<String> = (0..sp_proj.len())
        .map(|_| String::with_capacity(total_cols))
        .collect();

    for (i, &rb) in ref_bases.iter().enumerate() {
        ref_out.push(rb as char);
        for (q_idx, (_, _, _, proj)) in sp_proj.iter().enumerate() {
            q_outs[q_idx].push(proj.qry_bases[i] as char);
        }
        for k in 0..max_ins[i] {
            ref_out.push('-');
            for (q_idx, (_, _, _, proj)) in sp_proj.iter().enumerate() {
                let ch = proj.ins_after[i].get(k).copied().unwrap_or(b'-') as char;
                q_outs[q_idx].push(ch);
            }
        }
    }

    let rows = sp_proj
        .iter()
        .zip(q_outs)
        .map(|((sp, qc, st, proj), qseq)| ThreadedQueryRow {
            sp: sp.clone(),
            qry_contig: qc.clone(),
            qry_start: proj.qry_start,
            qry_size: proj.qry_size,
            qry_strand: *st,
            qry_src_size: proj.qry_src_size,
            seq: qseq,
        })
        .collect();
    (ref_out, rows)
}

fn parse_pairwise_maf_blocks(
    query_sp: &str,
    maf_text: &str,
    step_bp: usize,
    by_contig: &mut std::collections::BTreeMap<String, Vec<PairwiseMafBlock>>,
) -> std::collections::BTreeSet<String> {
    let mut cur_rows: Vec<(crate::formats::MafRow, usize)> = Vec::new();
    let mut tokens_by_contig: std::collections::BTreeMap<String, Vec<(usize, String, i8)>> =
        std::collections::BTreeMap::new();

    let flush = |rows: &mut Vec<(crate::formats::MafRow, usize)>,
                 by_contig: &mut std::collections::BTreeMap<String, Vec<PairwiseMafBlock>>,
                 tokens_by_contig: &mut std::collections::BTreeMap<
        String,
        Vec<(usize, String, i8)>,
    >| {
        if rows.len() >= 2 {
            let mut it = rows.drain(..);
            let (r, r_size) = it.next().expect("rows >= 2");
            let (q, _q_size) = it.next().expect("rows >= 2");
            if r.strand != '+' || r.seq.len() != q.seq.len() {
                return;
            }
            let Some(ref_end) = r.start.checked_add(r_size) else {
                return;
            };
            let b_idx = r.start / step_bp;
            tokens_by_contig.entry(r.src.clone()).or_default().push((
                r.start,
                format!("B{b_idx}"),
                if q.strand == '-' { -1 } else { 1 },
            ));
            by_contig.entry(r.src).or_default().push(PairwiseMafBlock {
                ref_start: r.start,
                ref_end,
                ref_src_size: r.src_size,
                qry_species: query_sp.to_owned(),
                qry_contig: q.src,
                qry_start: q.start,
                qry_strand: q.strand,
                qry_src_size: q.src_size,
                ref_seq: r.seq.into_bytes(),
                qry_seq: q.seq.into_bytes(),
            });
        }
        rows.clear();
    };

    for raw_line in maf_text.lines() {
        let line = raw_line.trim();
        if line.starts_with('#') {
            continue;
        }
        if line.is_empty() {
            flush(&mut cur_rows, by_contig, &mut tokens_by_contig);
            continue;
        }
        if line == "a"
            || line
                .strip_prefix('a')
                .is_some_and(|r| r.starts_with(char::is_whitespace))
        {
            flush(&mut cur_rows, by_contig, &mut tokens_by_contig);
        } else if let Some(parsed) = crate::formats::parse_maf_s_line(line) {
            cur_rows.push(parsed);
        }
    }
    flush(&mut cur_rows, by_contig, &mut tokens_by_contig);

    let mut sp_dups = std::collections::BTreeSet::new();
    for mut tokens in tokens_by_contig.into_values() {
        tokens.sort_unstable_by_key(|(start, _, _)| *start);
        let ordered: Vec<(String, i8)> = tokens.into_iter().map(|(_, tok, st)| (tok, st)).collect();
        crate::formats::collect_windowed_duplications(&ordered, &mut sp_dups);
    }
    sp_dups
}

struct MergeOut {
    maf_out: String,
    gfa_out: String,
    total_blocks: usize,
    multi_species_blocks: usize,
}

fn emit_merged_block(
    ref_contig: &str,
    (lo, hi): (usize, usize),
    ref_src_size: usize,
    ref_sub: &str,
    sp_rows: &[ThreadedQueryRow],
    out: &mut MergeOut,
    prev_seg: &mut Option<String>,
) {
    let _ = writeln!(out.maf_out, "a score={}.0", sp_rows.len() + 1);
    let _ = writeln!(
        out.maf_out,
        "s {ref_contig} {lo} {} + {ref_src_size} {ref_sub}",
        hi - lo
    );
    for row in sp_rows {
        let stripped = crate::formats::strip_sample_prefix(&row.qry_contig, &row.sp);
        let _ = writeln!(
            out.maf_out,
            "s {}.{stripped} {} {} {} {} {}",
            row.sp, row.qry_start, row.qry_size, row.qry_strand, row.qry_src_size, row.seq
        );
    }
    out.maf_out.push('\n');

    let seg_id = format!("ms_{ref_contig}_{lo}_{hi}");
    let _ = writeln!(
        out.gfa_out,
        "S\t{seg_id}\t*\tLN:i:{}\tSP:i:{}\tSN:Z:{ref_contig}\tSO:i:{lo}\tSR:i:0",
        hi - lo,
        sp_rows.len() + 1
    );
    if let Some(prev) = prev_seg.as_ref() {
        let _ = writeln!(out.gfa_out, "L\t{prev}\t+\t{seg_id}\t+\t0M");
    }
    *prev_seg = Some(seg_id);
}

fn merge_contig_blocks(ref_contig: &str, mut blocks: Vec<PairwiseMafBlock>, out: &mut MergeOut) {
    let mut prev_seg: Option<String> = None;
    blocks.sort_unstable_by_key(|b| (b.ref_start, b.ref_end));
    let mut bounds: Vec<usize> = blocks
        .iter()
        .flat_map(|b| [b.ref_start, b.ref_end])
        .collect();
    bounds.sort_unstable();
    bounds.dedup();

    let mut active: Vec<(usize, BlockCursor)> = Vec::new();
    let mut next_idx = 0usize;

    for w in bounds.windows(2) {
        let (lo, hi) = (w[0], w[1]);
        if lo >= hi {
            continue;
        }
        while next_idx < blocks.len() && blocks[next_idx].ref_start <= lo {
            if blocks[next_idx].ref_end >= hi {
                active.push((
                    next_idx,
                    BlockCursor {
                        col_idx: 0,
                        ref_pos: blocks[next_idx].ref_start,
                        qry_consumed: 0,
                    },
                ));
            }
            next_idx += 1;
        }
        active.retain(|&(idx, _)| blocks[idx].ref_end >= hi);
        if active.is_empty() {
            continue;
        }

        let mut ref_src_size_opt: Option<usize> = None;
        let mut ref_bases_opt: Option<Vec<u8>> = None;
        // Net overlapping alignments per query species by highest match score in [lo, hi).
        let mut best_by_sp: std::collections::BTreeMap<
            String,
            (isize, String, char, ProjectedInterval),
        > = std::collections::BTreeMap::new();
        for (idx, cursor) in &mut active {
            let blk = &blocks[*idx];
            if let Some(proj) = project_block_interval_with_cursor(blk, cursor, lo, hi) {
                ref_src_size_opt.get_or_insert(blk.ref_src_size);
                if ref_bases_opt.is_none() {
                    ref_bases_opt = Some(proj.ref_bases.clone());
                }
                let score = interval_match_score(&proj);
                match best_by_sp.get(&blk.qry_species) {
                    Some((best_sc, _, _, _)) if *best_sc >= score => {}
                    _ => {
                        best_by_sp.insert(
                            blk.qry_species.clone(),
                            (score, blk.qry_contig.clone(), blk.qry_strand, proj),
                        );
                    }
                }
            }
        }
        let sp_proj: Vec<(String, String, char, ProjectedInterval)> = best_by_sp
            .into_iter()
            .map(|(sp, (_, ctg, strand, proj))| (sp, ctg, strand, proj))
            .collect();
        let (Some(ref_src_size), Some(ref_bases)) = (ref_src_size_opt, ref_bases_opt) else {
            continue;
        };
        let (ref_sub, sp_rows) = thread_interval_columns(&ref_bases, &sp_proj);
        out.total_blocks += 1;
        if sp_rows.len() >= 2 {
            out.multi_species_blocks += 1;
        }
        emit_merged_block(
            ref_contig,
            (lo, hi),
            ref_src_size,
            &ref_sub,
            &sp_rows,
            out,
            &mut prev_seg,
        );
    }
}

/// Merge multiple pairwise MAF texts anchored on `ref_species` using exact reference-interval
/// splitting and base-level threading into multi-species MAF v1 + GFA v1.1 + e-graph paralog tags.
#[must_use]
pub fn merge_pairwise_mafs(
    ref_species: &str,
    pairwise_mafs: &[(String, String)],
    bin_bp: usize,
) -> (String, String, usize, usize) {
    use std::collections::{BTreeMap, BTreeSet};

    let mut by_contig: BTreeMap<String, Vec<PairwiseMafBlock>> = BTreeMap::new();
    let mut dups_by_species: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let step_bp = bin_bp.max(1);

    for (query_sp, maf_text) in pairwise_mafs {
        let sp_dups = parse_pairwise_maf_blocks(query_sp, maf_text, step_bp, &mut by_contig);
        dups_by_species
            .entry(query_sp.clone())
            .or_default()
            .extend(sp_dups);
    }

    let mut out = MergeOut {
        maf_out: String::from("##maf version=1 scoring=horizon_merge\n"),
        gfa_out: String::from("H\tVN:Z:1.1\n"),
        total_blocks: 0,
        multi_species_blocks: 0,
    };

    for (ref_contig, blocks) in by_contig {
        merge_contig_blocks(&ref_contig, blocks, &mut out);
    }

    for (sp, dups) in dups_by_species {
        if !dups.is_empty() {
            let joined: Vec<&str> = dups.iter().map(String::as_str).collect();
            let _ = writeln!(
                out.gfa_out,
                "# EG:Z:{ref_species}__{sp}\t{}",
                joined.join(",")
            );
        }
    }

    (
        out.maf_out,
        out.gfa_out,
        out.total_blocks,
        out.multi_species_blocks,
    )
}

fn interval_match_score(proj: &ProjectedInterval) -> isize {
    let mut score = 0isize;
    for (&rb, &qb) in proj.ref_bases.iter().zip(&proj.qry_bases) {
        if rb != b'-' && qb != b'-' {
            if rb.eq_ignore_ascii_case(&qb) {
                score += 2;
            } else {
                score -= 1;
            }
        }
    }
    score
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn project_block_interval_tracks_insertions_and_minus_strand() {
        let blk = PairwiseMafBlock {
            ref_start: 10,
            ref_end: 13,
            ref_src_size: 100,
            qry_species: "dog".to_owned(),
            qry_contig: "chr2".to_owned(),
            qry_start: 50,
            qry_strand: '-',
            qry_src_size: 200,
            ref_seq: b"A-CG".to_vec(),
            qry_seq: b"TTGC".to_vec(),
        };
        let proj = project_block_interval(&blk, 10, 13).expect("projected");
        assert_eq!(proj.ref_bases, b"ACG");
        assert_eq!(proj.qry_bases, b"TGC");
        assert_eq!(proj.ins_after, vec![vec![b'T'], Vec::new(), Vec::new()]);
        assert_eq!(proj.qry_size, 4);
        assert_eq!(proj.qry_start, 50);
    }

    #[test]
    fn project_block_interval_preserves_insertion_at_split_boundary() {
        let blk = PairwiseMafBlock {
            ref_start: 10,
            ref_end: 14,
            ref_src_size: 100,
            qry_species: "dog".to_owned(),
            qry_contig: "chr2".to_owned(),
            qry_start: 20,
            qry_strand: '+',
            qry_src_size: 200,
            ref_seq: b"AC-GT".to_vec(),
            qry_seq: b"AAAGG".to_vec(),
        };
        let left = project_block_interval(&blk, 10, 12).expect("left");
        assert_eq!(left.ref_bases, b"AC");
        assert_eq!(left.qry_bases, b"AA");
        assert_eq!(left.ins_after, vec![Vec::new(), vec![b'A']]);
        assert_eq!(left.qry_start, 20);
        assert_eq!(left.qry_size, 3);

        let right = project_block_interval(&blk, 12, 14).expect("right");
        assert_eq!(right.ref_bases, b"GT");
        assert_eq!(right.qry_bases, b"GG");
        assert_eq!(right.ins_after, vec![Vec::new(), Vec::new()]);
        assert_eq!(right.qry_start, 23);
        assert_eq!(right.qry_size, 2);
    }

    #[test]
    fn merge_pairwise_mafs_handles_prefixed_contigs_and_unsorted_blocks() {
        let maf_dog = "\
##maf version=1
a score=1
s human.chr1 10 4 + 100 ACGT
s dog.chr2 20 4 + 100 ACGT

a score=1
s human.chr1 0 4 + 100 TTTT
s dog.chr2 0 4 + 100 TTTT
";
        let (maf_out, _gfa_out, total_blocks, _multi) =
            merge_pairwise_mafs("human", &[("dog".to_owned(), maf_dog.to_owned())], 10);
        assert_eq!(total_blocks, 2);
        // Should not double-prefix "dog.dog.chr2"
        assert!(maf_out.contains("s dog.chr2 "));
        assert!(!maf_out.contains("dog.dog.chr2"));
        // Blocks should be emitted in sorted reference order (0 before 10).
        let first_offset = maf_out.find("s human.chr1 0 4").expect("block 0");
        let second_offset = maf_out.find("s human.chr1 10 4").expect("block 10");
        assert!(first_offset < second_offset);
    }
}
