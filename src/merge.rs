//! Reference-guided multi-species MAF merging and coordinate projection.

use rayon::prelude::*;
use std::fmt::Write as _;

struct PairwiseMafBlock {
    ref_start: usize,
    ref_end: usize,
    ref_src_size: usize,
    qry_sp_id: usize,
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
    /// Sparse `(column_offset, inserted_bases)` sorted by `column_offset`.
    ins_after: Vec<(usize, Vec<u8>)>,
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
    let mut ins_after: Vec<(usize, Vec<u8>)> = Vec::new();
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
                if let Some((last_idx, buf)) = ins_after.last_mut()
                    && *last_idx == prev_idx
                {
                    buf.push(qb);
                } else {
                    ins_after.push((prev_idx, vec![qb]));
                }
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

#[derive(Default)]
struct ThreadScratch {
    max_ins: Vec<usize>,
    ins_cursors: Vec<usize>,
    ref_out: String,
    q_outs: Vec<String>,
}

fn thread_interval_columns(
    ref_bases: &[u8],
    sp_proj: &[(&str, &str, char, ProjectedInterval)],
    scratch: &mut ThreadScratch,
) {
    let span = ref_bases.len();
    scratch.max_ins.clear();
    scratch.max_ins.resize(span, 0);
    for (_, _, _, proj) in sp_proj {
        for (col_idx, ins) in &proj.ins_after {
            scratch.max_ins[*col_idx] = scratch.max_ins[*col_idx].max(ins.len());
        }
    }
    let total_cols = span + scratch.max_ins.iter().sum::<usize>();
    scratch.ref_out.clear();
    scratch.ref_out.reserve(total_cols);
    if scratch.q_outs.len() < sp_proj.len() {
        scratch.q_outs.resize_with(sp_proj.len(), String::new);
    }
    for q in &mut scratch.q_outs[..sp_proj.len()] {
        q.clear();
        q.reserve(total_cols);
    }
    scratch.ins_cursors.clear();
    scratch.ins_cursors.resize(sp_proj.len(), 0);

    for (i, &rb) in ref_bases.iter().enumerate() {
        scratch.ref_out.push(rb as char);
        for (q_idx, (_, _, _, proj)) in sp_proj.iter().enumerate() {
            scratch.q_outs[q_idx].push(proj.qry_bases[i] as char);
        }
        let n_ins = scratch.max_ins[i];
        if n_ins > 0 {
            for k in 0..n_ins {
                scratch.ref_out.push('-');
                for (q_idx, (_, _, _, proj)) in sp_proj.iter().enumerate() {
                    let cur = scratch.ins_cursors[q_idx];
                    let ch = proj
                        .ins_after
                        .get(cur)
                        .filter(|(col_idx, _)| *col_idx == i)
                        .and_then(|(_, buf)| buf.get(k).copied())
                        .unwrap_or(b'-') as char;
                    scratch.q_outs[q_idx].push(ch);
                }
            }
            for (q_idx, (_, _, _, proj)) in sp_proj.iter().enumerate() {
                if proj
                    .ins_after
                    .get(scratch.ins_cursors[q_idx])
                    .is_some_and(|(col_idx, _)| *col_idx == i)
                {
                    scratch.ins_cursors[q_idx] += 1;
                }
            }
        }
    }
}

fn finish_parsed_maf_blocks(
    tokens_by_contig: std::collections::BTreeMap<String, Vec<(usize, String, i8)>>,
) -> std::collections::BTreeSet<String> {
    let mut sp_dups = std::collections::BTreeSet::new();
    for mut tokens in tokens_by_contig.into_values() {
        tokens.sort_unstable_by_key(|(start, _, _)| *start);
        let ordered: Vec<(String, i8)> = tokens.into_iter().map(|(_, tok, st)| (tok, st)).collect();
        crate::formats::collect_windowed_duplications(&ordered, &mut sp_dups);
    }
    sp_dups
}

fn push_parsed_stanza(
    rows: &mut Vec<(crate::formats::MafRow, usize)>,
    query_sp_id: usize,
    step_bp: usize,
    by_contig: &mut std::collections::BTreeMap<String, Vec<PairwiseMafBlock>>,
    tokens_by_contig: &mut std::collections::BTreeMap<String, Vec<(usize, String, i8)>>,
) {
    if rows.len() < 2 {
        return;
    }
    let mut it = rows.drain(..);
    let (mut r, r_size) = it.next().expect("rows >= 2");
    let (mut q, q_size) = it.next().expect("rows >= 2");
    if r.seq.len() != q.seq.len() {
        return;
    }
    if r.strand == '-' {
        let Some(r_fwd_start) = r.src_size.checked_sub(r.start.saturating_add(r_size)) else {
            return;
        };
        let Some(q_flipped_start) = q.src_size.checked_sub(q.start.saturating_add(q_size)) else {
            return;
        };
        r.start = r_fwd_start;
        r.strand = '+';
        q.start = q_flipped_start;
        q.strand = if q.strand == '-' { '+' } else { '-' };
        let rc_bytes = |s: &str| -> String {
            s.bytes()
                .rev()
                .map(|b| {
                    if b == b'-' {
                        b'-'
                    } else {
                        crate::seeding::complement_base(b)
                    }
                })
                .map(char::from)
                .collect()
        };
        r.seq = rc_bytes(&r.seq);
        q.seq = rc_bytes(&q.seq);
    } else if r.strand != '+' {
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
        qry_sp_id: query_sp_id,
        qry_contig: q.src,
        qry_start: q.start,
        qry_strand: q.strand,
        qry_src_size: q.src_size,
        ref_seq: r.seq.into_bytes(),
        qry_seq: q.seq.into_bytes(),
    });
}

type ParsedMafContigs = (
    std::collections::BTreeMap<String, Vec<PairwiseMafBlock>>,
    std::collections::BTreeSet<String>,
);

fn parse_pairwise_maf_blocks(
    query_sp_id: usize,
    maf_text: &str,
    step_bp: usize,
) -> ParsedMafContigs {
    let mut by_contig: std::collections::BTreeMap<String, Vec<PairwiseMafBlock>> =
        std::collections::BTreeMap::new();
    let mut tokens_by_contig: std::collections::BTreeMap<String, Vec<(usize, String, i8)>> =
        std::collections::BTreeMap::new();

    crate::formats::for_each_maf_stanza(maf_text, |rows| {
        push_parsed_stanza(
            rows,
            query_sp_id,
            step_bp,
            &mut by_contig,
            &mut tokens_by_contig,
        );
    });

    (by_contig, finish_parsed_maf_blocks(tokens_by_contig))
}

fn parse_pairwise_maf_reader<R: std::io::BufRead>(
    query_sp_id: usize,
    reader: R,
    step_bp: usize,
) -> anyhow::Result<ParsedMafContigs> {
    let mut by_contig: std::collections::BTreeMap<String, Vec<PairwiseMafBlock>> =
        std::collections::BTreeMap::new();
    let mut tokens_by_contig: std::collections::BTreeMap<String, Vec<(usize, String, i8)>> =
        std::collections::BTreeMap::new();

    crate::formats::for_each_maf_stanza_reader(reader, |rows| {
        push_parsed_stanza(
            rows,
            query_sp_id,
            step_bp,
            &mut by_contig,
            &mut tokens_by_contig,
        );
    })?;

    Ok((by_contig, finish_parsed_maf_blocks(tokens_by_contig)))
}

struct MergeOut {
    maf_out: String,
    gfa_out: String,
    total_blocks: usize,
    multi_species_blocks: usize,
}

#[allow(clippy::too_many_arguments)]
fn emit_merged_block(
    ref_species_prefix: Option<&str>,
    ref_contig: &str,
    (lo, hi): (usize, usize),
    ref_src_size: usize,
    ref_sub: &str,
    sp_proj: &[(&str, &str, char, ProjectedInterval)],
    q_outs: &[String],
    out: &mut MergeOut,
    prev_seg: &mut Option<(String, usize)>,
) {
    let _ = writeln!(out.maf_out, "a score={}.0", sp_proj.len() + 1);
    if let Some(sp) = ref_species_prefix {
        let stripped_ref = crate::formats::strip_sample_prefix(ref_contig, sp);
        let _ = writeln!(
            out.maf_out,
            "s {sp}.{stripped_ref} {lo} {} + {ref_src_size} {ref_sub}",
            hi - lo
        );
    } else {
        let _ = writeln!(
            out.maf_out,
            "s {ref_contig} {lo} {} + {ref_src_size} {ref_sub}",
            hi - lo
        );
    }
    for ((sp, qry_contig, qry_strand, proj), seq) in sp_proj.iter().zip(q_outs.iter()) {
        let stripped = crate::formats::strip_sample_prefix(qry_contig, sp);
        let _ = writeln!(
            out.maf_out,
            "s {sp}.{stripped} {} {} {qry_strand} {} {seq}",
            proj.qry_start, proj.qry_size, proj.qry_src_size
        );
    }
    out.maf_out.push('\n');

    let seg_id = format!("ms_{ref_contig}_{lo}_{hi}");
    let _ = writeln!(
        out.gfa_out,
        "S\t{seg_id}\t*\tLN:i:{}\tSP:i:{}\tSN:Z:{ref_contig}\tSO:i:{lo}\tSR:i:0",
        hi - lo,
        sp_proj.len() + 1
    );
    if let Some((prev, prev_hi)) = prev_seg.as_ref()
        && *prev_hi == lo
    {
        let _ = writeln!(out.gfa_out, "L\t{prev}\t+\t{seg_id}\t+\t0M");
    }
    *prev_seg = Some((seg_id, hi));
}

#[allow(clippy::too_many_lines)]
fn merge_contig_blocks(
    ref_species_prefix: Option<&str>,
    ref_contig: &str,
    mut blocks: Vec<PairwiseMafBlock>,
    sp_names: &[String],
) -> MergeOut {
    let mut out = MergeOut {
        maf_out: String::new(),
        gfa_out: String::new(),
        total_blocks: 0,
        multi_species_blocks: 0,
    };
    let mut prev_seg: Option<(String, usize)> = None;
    blocks.sort_unstable_by(|a, b| {
        (a.ref_start, a.ref_end, a.qry_sp_id)
            .cmp(&(b.ref_start, b.ref_end, b.qry_sp_id))
            .then_with(|| a.qry_contig.cmp(&b.qry_contig))
            .then_with(|| (a.qry_start, a.qry_strand).cmp(&(b.qry_start, b.qry_strand)))
    });
    let mut bounds: Vec<usize> = blocks
        .iter()
        .flat_map(|b| [b.ref_start, b.ref_end])
        .collect();
    bounds.sort_unstable();
    bounds.dedup();

    let mut active: Vec<(usize, BlockCursor)> = Vec::new();
    let mut next_idx = 0usize;
    let mut best_by_sp: Vec<Option<(i64, usize, ProjectedInterval)>> =
        (0..sp_names.len()).map(|_| None).collect();
    let mut sp_proj: Vec<(&str, &str, char, ProjectedInterval)> =
        Vec::with_capacity(sp_names.len());
    let mut scratch = ThreadScratch::default();

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
        for slot in &mut best_by_sp {
            *slot = None;
        }
        // Net overlapping alignments per query species by highest match score in [lo, hi).
        for (idx, cursor) in &mut active {
            let blk = &blocks[*idx];
            if let Some(proj) = project_block_interval_with_cursor(blk, cursor, lo, hi) {
                ref_src_size_opt.get_or_insert(blk.ref_src_size);
                if ref_bases_opt.is_none() {
                    ref_bases_opt = Some(proj.ref_bases.clone());
                }
                let score = interval_match_score(&proj);
                match &best_by_sp[blk.qry_sp_id] {
                    Some((best_sc, _, _)) if *best_sc >= score => {}
                    _ => {
                        best_by_sp[blk.qry_sp_id] = Some((score, *idx, proj));
                    }
                }
            }
        }
        sp_proj.clear();
        for (sp_id, slot) in best_by_sp.iter_mut().enumerate() {
            if let Some((_, blk_idx, proj)) = slot.take() {
                let blk = &blocks[blk_idx];
                sp_proj.push((
                    sp_names[sp_id].as_str(),
                    blk.qry_contig.as_str(),
                    blk.qry_strand,
                    proj,
                ));
            }
        }
        let (Some(ref_src_size), Some(ref_bases)) = (ref_src_size_opt, ref_bases_opt) else {
            continue;
        };
        thread_interval_columns(&ref_bases, &sp_proj, &mut scratch);
        out.total_blocks += 1;
        if sp_proj.len() >= 2 {
            out.multi_species_blocks += 1;
        }
        emit_merged_block(
            ref_species_prefix,
            ref_contig,
            (lo, hi),
            ref_src_size,
            &scratch.ref_out,
            &sp_proj,
            &scratch.q_outs[..sp_proj.len()],
            &mut out,
            &mut prev_seg,
        );
    }
    out
}

/// Direct in-memory contig-pair alignment input for [`merge_direct_alignments_prefixed_to_writers`],
/// avoiding intermediate MAF text serialization and re-parsing.
pub struct DirectContigAlignment<'a> {
    /// Reference contig name.
    pub ref_contig: &'a str,
    /// Reference contig bases.
    pub ref_seq: &'a [u8],
    /// Query species name.
    pub qry_species: &'a str,
    /// Query contig name.
    pub qry_contig: &'a str,
    /// Query contig bases.
    pub qry_seq: &'a [u8],
    /// Run-length encoded alignment between `ref_seq` and `qry_seq`.
    pub aln: &'a crate::mdl::RunAlignment,
}

fn build_species_index<'a>(
    names: impl Iterator<Item = &'a str>,
) -> (Vec<String>, std::collections::HashMap<String, usize>) {
    let mut sp_names: Vec<String> = names.map(str::to_owned).collect();
    sp_names.sort();
    sp_names.dedup();
    let sp_to_id = sp_names
        .iter()
        .enumerate()
        .map(|(idx, s)| (s.clone(), idx))
        .collect();
    (sp_names, sp_to_id)
}

type DupToken = (usize, String, i8);
type ContigBlockMap = std::collections::BTreeMap<String, Vec<PairwiseMafBlock>>;

fn direct_alignment_to_blocks(
    item: &DirectContigAlignment<'_>,
    sp_id: usize,
    step_bp: usize,
) -> (ContigBlockMap, Vec<DupToken>) {
    let mut by_contig: std::collections::BTreeMap<String, Vec<PairwiseMafBlock>> =
        std::collections::BTreeMap::new();
    let mut tokens: Vec<(usize, String, i8)> = Vec::new();

    for (grp, is_rev) in crate::netting::partition_colinear_run_groups(&item.aln.aligned_runs) {
        let Some(m) = crate::formats::materialize_run_group(
            item.ref_seq.len(),
            item.qry_seq.len(),
            item.ref_seq,
            item.qry_seq,
            &grp,
            is_rev,
        ) else {
            continue;
        };
        tokens.push((
            m.x0,
            format!("B{}", m.x0 / step_bp),
            if m.y_strand == '-' { -1 } else { 1 },
        ));
        by_contig
            .entry(item.ref_contig.to_owned())
            .or_default()
            .push(PairwiseMafBlock {
                ref_start: m.x0,
                ref_end: m.x_end,
                ref_src_size: item.ref_seq.len(),
                qry_sp_id: sp_id,
                qry_contig: item.qry_contig.to_owned(),
                qry_start: m.y_start,
                qry_strand: m.y_strand,
                qry_src_size: item.qry_seq.len(),
                ref_seq: m.ax,
                qry_seq: m.ay,
            });
    }

    (by_contig, tokens)
}

type ParsedSpeciesBlocks = (
    usize,
    std::collections::BTreeMap<String, Vec<PairwiseMafBlock>>,
    std::collections::BTreeSet<String>,
);

fn merge_parsed_blocks_to_writers<W1: std::io::Write, W2: std::io::Write>(
    ref_species: &str,
    ref_species_prefix: Option<&str>,
    sp_names: &[String],
    parsed: Vec<ParsedSpeciesBlocks>,
    maf_out: &mut W1,
    gfa_out: &mut W2,
) -> anyhow::Result<(usize, usize)> {
    use std::collections::{BTreeMap, BTreeSet};

    let mut by_contig: BTreeMap<String, Vec<PairwiseMafBlock>> = BTreeMap::new();
    let mut dups_by_species: BTreeMap<&str, BTreeSet<String>> = BTreeMap::new();
    for (sp_id, local_contigs, sp_dups) in parsed {
        for (ctg, mut blks) in local_contigs {
            by_contig.entry(ctg).or_default().append(&mut blks);
        }
        dups_by_species
            .entry(sp_names[sp_id].as_str())
            .or_default()
            .extend(sp_dups);
    }

    maf_out.write_all(b"##maf version=1 scoring=horizon_merge\n")?;
    gfa_out.write_all(b"H\tVN:Z:1.1\n")?;
    let mut total_blocks = 0usize;
    let mut multi_species_blocks = 0usize;

    // Process contigs in bounded parallel batches and flush each batch immediately so we do not
    // hold all formatted multi-species MAF/GFA strings in RAM simultaneously.
    let batch_size = rayon::current_num_threads().max(1);
    let mut contig_iter = by_contig.into_iter();
    loop {
        let batch: Vec<(String, Vec<PairwiseMafBlock>)> =
            contig_iter.by_ref().take(batch_size).collect();
        if batch.is_empty() {
            break;
        }
        let contig_outs: Vec<MergeOut> = batch
            .into_par_iter()
            .map(|(ref_contig, blocks)| {
                merge_contig_blocks(ref_species_prefix, &ref_contig, blocks, sp_names)
            })
            .collect();
        for chunk in contig_outs {
            maf_out.write_all(chunk.maf_out.as_bytes())?;
            gfa_out.write_all(chunk.gfa_out.as_bytes())?;
            total_blocks += chunk.total_blocks;
            multi_species_blocks += chunk.multi_species_blocks;
        }
    }
    for (sp, dups) in dups_by_species {
        if !dups.is_empty() {
            let joined: Vec<&str> = dups.iter().map(String::as_str).collect();
            writeln!(gfa_out, "# EG:Z:{ref_species}__{sp}\t{}", joined.join(","))?;
        }
    }
    Ok((total_blocks, multi_species_blocks))
}

/// Merge direct in-memory [`DirectContigAlignment`]s anchored on `ref_species` (prefixing
/// reference rows with `<ref_species>.`) and stream directly to `maf_out` and `gfa_out`.
///
/// # Errors
/// Returns an error if writing to `maf_out` or `gfa_out` fails.
pub fn merge_direct_alignments_prefixed_to_writers<W1: std::io::Write, W2: std::io::Write>(
    ref_species: &str,
    alignments: &[DirectContigAlignment<'_>],
    bin_bp: usize,
    maf_out: &mut W1,
    gfa_out: &mut W2,
) -> anyhow::Result<(usize, usize)> {
    use std::collections::{BTreeMap, BTreeSet};

    let (sp_names, sp_to_id) = build_species_index(alignments.iter().map(|a| a.qry_species));
    let step_bp = bin_bp.max(1);
    let raw_parsed: Vec<(usize, &str, ContigBlockMap, Vec<DupToken>)> = alignments
        .par_iter()
        .map(|item| {
            let sp_id = sp_to_id[item.qry_species];
            let (local_contigs, tokens) = direct_alignment_to_blocks(item, sp_id, step_bp);
            (sp_id, item.ref_contig, local_contigs, tokens)
        })
        .collect();

    let mut tokens_by_sp_contig: BTreeMap<(usize, &str), Vec<DupToken>> = BTreeMap::new();
    let mut by_sp: Vec<BTreeMap<String, Vec<PairwiseMafBlock>>> =
        (0..sp_names.len()).map(|_| BTreeMap::new()).collect();
    for (sp_id, ref_ctg, local_contigs, mut tokens) in raw_parsed {
        for (ctg, mut blks) in local_contigs {
            by_sp[sp_id].entry(ctg).or_default().append(&mut blks);
        }
        tokens_by_sp_contig
            .entry((sp_id, ref_ctg))
            .or_default()
            .append(&mut tokens);
    }

    let mut dups_by_sp: Vec<BTreeSet<String>> = vec![BTreeSet::new(); sp_names.len()];
    for ((sp_id, _), mut tokens) in tokens_by_sp_contig {
        tokens.sort_unstable_by_key(|(start, _, _)| *start);
        let ordered: Vec<(String, i8)> = tokens.into_iter().map(|(_, tok, st)| (tok, st)).collect();
        crate::formats::collect_windowed_duplications(&ordered, &mut dups_by_sp[sp_id]);
    }

    let parsed: Vec<ParsedSpeciesBlocks> = by_sp
        .into_iter()
        .zip(dups_by_sp)
        .enumerate()
        .map(|(sp_id, (contigs, dups))| (sp_id, contigs, dups))
        .collect();

    merge_parsed_blocks_to_writers(
        ref_species,
        Some(ref_species),
        &sp_names,
        parsed,
        maf_out,
        gfa_out,
    )
}

/// Stream pairwise MAF files from disk (`(species, path)`) in parallel, avoiding loading
/// full raw MAF files into memory, and write multi-species MAF v1 + GFA v1.1 directly.
///
/// # Errors
/// Returns an error if opening/reading any input MAF or writing to `maf_out`/`gfa_out` fails.
pub fn merge_pairwise_maf_paths_to_writers<W1: std::io::Write, W2: std::io::Write>(
    ref_species: &str,
    pairwise_paths: &[(String, std::path::PathBuf)],
    bin_bp: usize,
    maf_out: &mut W1,
    gfa_out: &mut W2,
) -> anyhow::Result<(usize, usize)> {
    use anyhow::Context;

    let (sp_names, sp_to_id) =
        build_species_index(pairwise_paths.iter().map(|(sp, _)| sp.as_str()));

    let step_bp = bin_bp.max(1);
    let parsed: Vec<ParsedSpeciesBlocks> = pairwise_paths
        .par_iter()
        .map(|(query_sp, path)| {
            let sp_id = sp_to_id[query_sp.as_str()];
            let reader = crate::io_util::open_text_reader(path)
                .with_context(|| format!("opening {}", path.display()))?;
            let (local_contigs, sp_dups) = parse_pairwise_maf_reader(sp_id, reader, step_bp)
                .with_context(|| format!("reading {}", path.display()))?;
            Ok((sp_id, local_contigs, sp_dups))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;

    merge_parsed_blocks_to_writers(ref_species, None, &sp_names, parsed, maf_out, gfa_out)
}

/// Merge multiple pairwise MAF texts anchored on `ref_species`, prefixing reference `s` rows
/// with `<ref_species>.` and streaming multi-species MAF v1 + GFA v1.1 + e-graph paralog tags
/// directly to `maf_out` and `gfa_out`.
///
/// # Errors
/// Returns an error if writing to `maf_out` or `gfa_out` fails.
pub fn merge_pairwise_mafs_prefixed_to_writers<W1: std::io::Write, W2: std::io::Write>(
    ref_species: &str,
    pairwise_mafs: &[(String, String)],
    bin_bp: usize,
    maf_out: &mut W1,
    gfa_out: &mut W2,
) -> anyhow::Result<(usize, usize)> {
    merge_pairwise_mafs_impl(
        ref_species,
        Some(ref_species),
        pairwise_mafs,
        bin_bp,
        maf_out,
        gfa_out,
    )
}

fn merge_pairwise_mafs_impl<W1: std::io::Write, W2: std::io::Write>(
    ref_species: &str,
    ref_species_prefix: Option<&str>,
    pairwise_mafs: &[(String, String)],
    bin_bp: usize,
    maf_out: &mut W1,
    gfa_out: &mut W2,
) -> anyhow::Result<(usize, usize)> {
    let (sp_names, sp_to_id) = build_species_index(pairwise_mafs.iter().map(|(sp, _)| sp.as_str()));

    let step_bp = bin_bp.max(1);
    let parsed: Vec<ParsedSpeciesBlocks> = pairwise_mafs
        .par_iter()
        .map(|(query_sp, maf_text)| {
            let sp_id = sp_to_id[query_sp.as_str()];
            let (local_contigs, sp_dups) = parse_pairwise_maf_blocks(sp_id, maf_text, step_bp);
            (sp_id, local_contigs, sp_dups)
        })
        .collect();

    merge_parsed_blocks_to_writers(
        ref_species,
        ref_species_prefix,
        &sp_names,
        parsed,
        maf_out,
        gfa_out,
    )
}

/// Merge multiple pairwise MAF texts anchored on `ref_species` using exact reference-interval
/// splitting and base-level threading into multi-species MAF v1 + GFA v1.1 + e-graph paralog tags.
///
/// # Panics
/// Panics if internal UTF-8 formatting of the merged MAF or GFA output is invalid.
#[must_use]
pub fn merge_pairwise_mafs(
    ref_species: &str,
    pairwise_mafs: &[(String, String)],
    bin_bp: usize,
) -> (String, String, usize, usize) {
    let mut maf_buf = Vec::new();
    let mut gfa_buf = Vec::new();
    let (total_blocks, multi_species_blocks) = merge_pairwise_mafs_impl(
        ref_species,
        None,
        pairwise_mafs,
        bin_bp,
        &mut maf_buf,
        &mut gfa_buf,
    )
    .expect("in-memory Vec write cannot fail");
    let maf_out = String::from_utf8(maf_buf).expect("merged MAF is valid UTF-8");
    let gfa_out = String::from_utf8(gfa_buf).expect("merged GFA is valid UTF-8");
    (maf_out, gfa_out, total_blocks, multi_species_blocks)
}

fn interval_match_score(proj: &ProjectedInterval) -> i64 {
    crate::netting::score_acgt_columns(&proj.ref_bases, &proj.qry_bases)
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
            qry_sp_id: 0,
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
        assert_eq!(proj.ins_after, vec![(0, vec![b'T'])]);
        assert_eq!(proj.qry_size, 4);
        assert_eq!(proj.qry_start, 50);
    }

    #[test]
    fn project_block_interval_preserves_insertion_at_split_boundary() {
        let blk = PairwiseMafBlock {
            ref_start: 10,
            ref_end: 14,
            ref_src_size: 100,
            qry_sp_id: 0,
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
        assert_eq!(left.ins_after, vec![(1, vec![b'A'])]);
        assert_eq!(left.qry_start, 20);
        assert_eq!(left.qry_size, 3);

        let right = project_block_interval(&blk, 12, 14).expect("right");
        assert_eq!(right.ref_bases, b"GT");
        assert_eq!(right.qry_bases, b"GG");
        assert_eq!(right.ins_after, []);
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
        let (maf_out, gfa_out, total_blocks, _multi) =
            merge_pairwise_mafs("human", &[("dog".to_owned(), maf_dog.to_owned())], 10);
        assert_eq!(total_blocks, 2);
        // Should not double-prefix "dog.dog.chr2"
        assert!(maf_out.contains("s dog.chr2 "));
        assert!(!maf_out.contains("dog.dog.chr2"));
        // Blocks should be emitted in sorted reference order (0 before 10).
        let first_offset = maf_out.find("s human.chr1 0 4").expect("block 0");
        let second_offset = maf_out.find("s human.chr1 10 4").expect("block 10");
        assert!(first_offset < second_offset);
        // Unaligned gap [4..10) must NOT produce a bogus 0M GFA L-link.
        assert!(
            !gfa_out.contains("\nL\t"),
            "unexpected 0M link across unaligned gap:\n{gfa_out}"
        );

        // Adjacent blocks [0..4) and [4..8) MUST emit a 0M GFA L-link.
        let maf_adj = "\
##maf version=1
a score=1
s human.chr1 0 4 + 100 TTTT
s dog.chr2 0 4 + 100 TTTT

a score=1
s human.chr1 4 4 + 100 ACGT
s dog.chr2 4 4 + 100 ACGT
";
        let (_, gfa_adj, adj_blocks, _) =
            merge_pairwise_mafs("human", &[("dog".to_owned(), maf_adj.to_owned())], 10);
        assert_eq!(adj_blocks, 2);
        assert!(
            gfa_adj.contains("L\tms_human.chr1_0_4\t+\tms_human.chr1_4_8\t+\t0M"),
            "expected contiguous 0M link, got:\n{gfa_adj}"
        );
    }
}
