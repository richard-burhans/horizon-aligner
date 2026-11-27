//! Streaming progressive HDF5 v0 HAL exporter with windowed two-pass Fitch parsimony.

use super::fitch::{
    FITCH_ALL, add_consensus_base, fitch_base_mask, fitch_resolve_child, fitch_resolve_root,
    fitch_step, pack_fitch_masks_into,
};
use super::hdf5_write::{
    HalGenomeDatasets, encode_hal_nibble, encode_hal_seq_names, init_hdf5_root,
    write_bottom_sentinel, write_hal_genome_group, write_top_record, write_top_sentinel,
};
use super::hdf5_writer::Hdf5V0Builder;
use super::normalize::{
    CompactTop, NO_PARALOG, ProgressivePartition, normalize_star_hal_in_place,
    partition_progressive_segments,
};
use super::{HalAlignment, HalSequence, build_hal_sequences};
use std::io::{BufWriter, Seek, Write};

fn write_hal_dna_contigs_stream<W: Write>(
    contigs: &[(String, Vec<u8>)],
    w: &mut W,
) -> std::io::Result<u64> {
    let mut buf = [0u8; 8192];
    let mut written = 0u64;
    let mut pos = 0usize;
    let mut pending_hi: Option<u8> = None;
    for (_, seq) in contigs {
        for &base in seq {
            let nib = encode_hal_nibble(base);
            if let Some(hi) = pending_hi.take() {
                buf[pos] = (hi << 4) | nib;
                pos += 1;
                if pos == buf.len() {
                    w.write_all(&buf)?;
                    written += pos as u64;
                    pos = 0;
                }
            } else {
                pending_hi = Some(nib);
            }
        }
    }
    if let Some(hi) = pending_hi {
        buf[pos] = hi << 4;
        pos += 1;
    }
    if pos > 0 {
        w.write_all(&buf[..pos])?;
        written += pos as u64;
    }
    Ok(written)
}

fn encode_hal_dna_into(dna: &[u8], dst: &mut Vec<u8>) {
    dst.clear();
    dst.reserve(dna.len().div_ceil(2));
    for chunk in dna.chunks(2) {
        let hi = encode_hal_nibble(chunk[0]);
        let lo = chunk.get(1).copied().map_or(0, encode_hal_nibble);
        dst.push((hi << 4) | lo);
    }
}

fn copy_contig_range(
    contigs: &[(String, Vec<u8>)],
    seqs: &[HalSequence],
    start: usize,
    end: usize,
    dst: &mut Vec<u8>,
) {
    dst.clear();
    if start >= end || contigs.is_empty() {
        return;
    }
    dst.reserve(end - start);
    let mut cur = start;
    let mut ci = seqs.partition_point(|s| s.offset.saturating_add(s.length) <= cur);
    while cur < end && ci < contigs.len() {
        let s = &seqs[ci];
        let local_s = cur.saturating_sub(s.offset);
        let avail = contigs[ci].1.len().saturating_sub(local_s);
        let take = avail.min(end - cur);
        if take > 0 {
            dst.extend_from_slice(&contigs[ci].1[local_s..local_s + take]);
            cur += take;
        }
        ci += 1;
    }
}

struct RefStreamMeta<'a> {
    dna_len: usize,
    dna_packed_len: u64,
    rup_idx: u32,
    seqs: &'a [HalSequence],
    seqname_bytes: &'a [u8],
    seqname_stride: usize,
    bot_offsets: &'a [u64],
    bot_starts: &'a [u64],
    bot_lens: &'a [u64],
}

#[allow(
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::type_complexity
)]
fn write_windowed_fitch_ancestors<W: Write + Seek>(
    w: &mut Hdf5V0Builder<W>,
    ref_contigs: &[(String, Vec<u8>)],
    ref_meta: &RefStreamMeta<'_>,
    queries: &[(String, Vec<(String, Vec<u8>)>)],
    qry_seqs: &[Vec<HalSequence>],
    q_tops: &[Vec<CompactTop>],
    window_bp: usize,
) -> std::io::Result<Vec<u64>> {
    let n_q = queries.len();
    let ref_dna_len = ref_meta.dna_len;
    let mut anc_dna_hdrs = vec![0u64; n_q];
    let mut anc_dna_offsets = vec![0u64; n_q];
    for level in (0..n_q).rev() {
        let (hdr, data_off) =
            w.reserve_contiguous_dataset(&[ref_meta.dna_packed_len], 1, ref_meta.dna_packed_len)?;
        anc_dna_hdrs[level] = hdr;
        anc_dna_offsets[level] = data_off;
    }
    if ref_dna_len == 0 {
        return Ok(anc_dna_hdrs);
    }

    // Index each query's chain-head top segments sorted by parent bottom index so every
    // reference window queries only its overlapping segments in O(log T).
    let chain_heads: Vec<Vec<u32>> = q_tops
        .iter()
        .map(|tops| {
            let mut heads: Vec<u32> = tops
                .iter()
                .enumerate()
                .filter_map(|(ti, t)| t.is_chain_head.then(|| u32::try_from(ti).ok()).flatten())
                .collect();
            heads.sort_unstable_by_key(|&ti| tops[ti as usize].parent_bot);
            heads
        })
        .collect();

    // Force an even window >= 2 so every window starts on an exact packed-byte offset win_s/2.
    let win_step = window_bp.max(2) & !1usize;
    let max_win = ref_dna_len.min(win_step);
    let mut ref_win = Vec::with_capacity(max_win);
    let mut leaf_counts: Vec<[u16; 4]> = Vec::with_capacity(max_win);
    let mut prev_masks: Vec<u8> = Vec::with_capacity(max_win);
    let mut q_masks: Vec<u8> = Vec::with_capacity(max_win);
    // Preallocate and reuse `n_q` packed Fitch mask buffers across all windows (zero per-window
    // heap churn on 100-way 3 Gbp alignments).
    let mut level_sets: Vec<Vec<u8>> = (0..n_q)
        .map(|_| Vec::with_capacity(max_win.div_ceil(2)))
        .collect();
    let mut parent_dna = vec![0u8; max_win];
    let mut anc_dna = vec![0u8; max_win];
    let mut packed_buf = Vec::with_capacity(max_win.div_ceil(2));

    let mut win_s = 0usize;
    while win_s < ref_dna_len {
        let win_e = win_s.saturating_add(win_step).min(ref_dna_len);
        let win_len = win_e - win_s;
        copy_contig_range(ref_contigs, ref_meta.seqs, win_s, win_e, &mut ref_win);

        let bi_lo = ref_meta
            .bot_starts
            .partition_point(|&st| st <= win_s as u64)
            .saturating_sub(1);
        let bi_hi = ref_meta.bot_starts.partition_point(|&st| st < win_e as u64);

        leaf_counts.clear();
        prev_masks.clear();
        for &rb in &ref_win {
            let mut c = [0u16; 4];
            add_consensus_base(&mut c, rb);
            leaf_counts.push(c);
            prev_masks.push(fitch_base_mask(rb));
        }
        q_masks.resize(win_len, FITCH_ALL);

        for level in 0..n_q {
            q_masks.fill(FITCH_ALL);
            let tops = &q_tops[level];
            let heads = &chain_heads[level];
            let h_lo = heads.partition_point(|&ti| (tops[ti as usize].parent_bot as usize) < bi_lo);
            let h_hi = heads.partition_point(|&ti| (tops[ti as usize].parent_bot as usize) < bi_hi);
            let q_contigs = &queries[level].1;
            let q_seq_meta = &qry_seqs[level];
            for &ti in &heads[h_lo..h_hi] {
                let top = &tops[ti as usize];
                let bi = top.parent_bot as usize;
                let r_start = usize::try_from(ref_meta.bot_starts[bi]).unwrap_or(0);
                let b_len = usize::try_from(ref_meta.bot_lens[bi]).unwrap_or(0);
                let r_end = r_start.saturating_add(b_len).min(ref_dna_len);
                let ov_s = r_start.max(win_s);
                let ov_e = r_end.min(win_e);
                if ov_s >= ov_e {
                    continue;
                }
                let t_start = usize::try_from(top.start).unwrap_or(0);
                let q_ci =
                    q_seq_meta.partition_point(|s| s.offset.saturating_add(s.length) <= t_start);
                let (Some(q_smeta), Some((_, q_cseq))) =
                    (q_seq_meta.get(q_ci), q_contigs.get(q_ci))
                else {
                    continue;
                };
                let Some(q_local_base) = t_start.checked_sub(q_smeta.offset) else {
                    continue;
                };
                for r_pos in ov_s..ov_e {
                    let off = r_pos - r_start;
                    let q_off = if top.rev { b_len - 1 - off } else { off };
                    if let Some(&qb) = q_cseq.get(q_local_base + q_off) {
                        let oriented = if top.rev {
                            crate::seeding::complement_base(qb)
                        } else {
                            qb
                        };
                        let local = r_pos - win_s;
                        q_masks[local] = fitch_base_mask(oriented);
                        add_consensus_base(&mut leaf_counts[local], oriented);
                    }
                }
            }
            for (prev, &r) in prev_masks.iter_mut().zip(q_masks.iter()) {
                *prev = fitch_step(*prev, r);
            }
            pack_fitch_masks_into(&prev_masks, &mut level_sets[level]);
        }

        let p_slice = &mut parent_dna[..win_len];
        let a_slice = &mut anc_dna[..win_len];
        for level in (0..n_q).rev() {
            if level + 1 == n_q {
                fitch_resolve_root(&ref_win, &level_sets[level], &leaf_counts, a_slice);
            } else {
                fitch_resolve_child(&ref_win, p_slice, &level_sets[level], a_slice);
            }
            encode_hal_dna_into(a_slice, &mut packed_buf);
            w.write_at(anc_dna_offsets[level] + (win_s / 2) as u64, &packed_buf)?;
            p_slice.copy_from_slice(a_slice);
        }
        win_s = win_e;
    }

    Ok(anc_dna_hdrs)
}

impl HalAlignment {
    /// Stream a progressive multi-contig `HalAlignment` directly to an HDF5 `v0` `.hal` file
    /// without materializing all ancestral `HalGenome` segment tables in memory simultaneously.
    ///
    /// Produces the same tree, ancestral consensus sequences, and sorted HAL datasets as
    /// [`Self::from_multi_contig_progressive_blocks`] followed by [`Self::write_hdf5_hal`],
    /// while reducing peak memory by an order of magnitude on multi-genome datasets.
    ///
    /// Returns `(tree_newick, root_bottom_segment_count)`.
    ///
    /// # Errors
    /// Returns an `std::io::Error` if creating or writing the file fails.
    #[allow(clippy::type_complexity)]
    pub fn write_multi_contig_progressive_hdf5_hal(
        path: &std::path::Path,
        ref_name: &str,
        ref_contigs: &[(String, Vec<u8>)],
        queries: &[(String, Vec<(String, Vec<u8>)>)],
        blocks: &[(usize, usize, usize, usize, usize, usize, bool)],
    ) -> std::io::Result<(String, usize)> {
        // 4 MiB bases per Fitch window (even so each window aligns to packed byte offset win_s/2):
        // caps working RAM at ~240 MiB even for 100-way 3 Gbp mammalian alignments.
        const DEFAULT_FITCH_WINDOW_BP: usize = 4 * 1024 * 1024;
        Self::write_multi_contig_progressive_hdf5_hal_with_window(
            path,
            ref_name,
            ref_contigs,
            queries,
            blocks,
            DEFAULT_FITCH_WINDOW_BP,
        )
    }

    #[allow(clippy::type_complexity)]
    pub(super) fn write_multi_contig_progressive_hdf5_hal_with_window(
        path: &std::path::Path,
        ref_name: &str,
        ref_contigs: &[(String, Vec<u8>)],
        queries: &[(String, Vec<(String, Vec<u8>)>)],
        blocks: &[(usize, usize, usize, usize, usize, usize, bool)],
        window_bp: usize,
    ) -> std::io::Result<(String, usize)> {
        if queries.len() < 2 {
            let mut hal =
                Self::from_multi_contig_star_blocks(ref_name, ref_contigs, queries, blocks);
            normalize_star_hal_in_place(&mut hal);
            let desc = hal.tree_newick.clone();
            let nb = hal.genomes[0].bottom_segments.len();
            hal.write_hdf5_hal(path)?;
            return Ok((desc, nb));
        }

        let n_q = queries.len();
        let (ref_dna_len, ref_seqs) = build_hal_sequences(ref_name, ref_contigs);
        let mut qry_lens: Vec<usize> = Vec::with_capacity(n_q);
        let mut qry_seqs: Vec<Vec<HalSequence>> = Vec::with_capacity(n_q);
        for (q_name, q_contigs) in queries {
            let (q_len, seqs) = build_hal_sequences(q_name, q_contigs);
            qry_lens.push(q_len);
            qry_seqs.push(seqs);
        }

        // Partition reference coordinate space [0, ref_dna_len) into contiguous, non-overlapping
        // bottom segments and normalize per-query top segments with paralogy cycles for libHal.
        let ProgressivePartition {
            bot_starts,
            bot_lens,
            mut q_tops,
        } = partition_progressive_segments(ref_dna_len, &ref_seqs, &qry_seqs, blocks)?;

        let n_b = bot_starts.len();
        let mut tree_newick = ref_name.to_owned();
        for (j, (q_name, _)) in queries.iter().enumerate() {
            tree_newick = format!("({tree_newick},{q_name})Anc{j}");
        }
        tree_newick.push(';');

        let file = std::fs::File::create(path)?;
        let mut w = Hdf5V0Builder::new(BufWriter::with_capacity(256 * 1024, file))?;
        let (gcol_addr, mut owned_root_entries) =
            init_hdf5_root(&mut w, &tree_newick, 1 + 2 * n_q)?;

        // Precompute ref SEQNAME_ARRAY once and reuse for all ancestors + ref leaf.
        let (ref_seqname_bytes, ref_seqname_stride) = encode_hal_seq_names(&ref_seqs)?;
        let ref_bot_offsets: Vec<u64> = ref_seqs
            .iter()
            .map(|s| bot_starts.partition_point(|&st| st < s.offset as u64) as u64)
            .collect();
        let ref_meta = RefStreamMeta {
            dna_len: ref_dna_len,
            dna_packed_len: ref_dna_len.div_ceil(2) as u64,
            rup_idx: if ref_dna_len % 2 == 1 { 4 } else { 3 },
            seqs: &ref_seqs,
            seqname_bytes: &ref_seqname_bytes,
            seqname_stride: ref_seqname_stride,
            bot_offsets: &ref_bot_offsets,
            bot_starts: &bot_starts,
            bot_lens: &bot_lens,
        };

        let anc_dna_hdrs = write_windowed_fitch_ancestors(
            &mut w,
            ref_contigs,
            &ref_meta,
            queries,
            &qry_seqs,
            &q_tops,
            window_bp,
        )?;

        let mut right_links: Vec<(i64, u8)> = vec![(-1, 0); n_b];
        for level in (0..n_q).rev() {
            let is_root = level + 1 == n_q;
            let g_hdr = write_progressive_ancestor_genome(
                &mut w,
                &ref_meta,
                &q_tops[level],
                is_root,
                gcol_addr,
                anc_dna_hdrs[level],
                &mut right_links,
            )?;
            owned_root_entries.push((format!("Anc{level}"), g_hdr));
        }
        drop(right_links);

        let ref_hdr = write_ref_leaf_genome(&mut w, ref_contigs, &ref_meta, gcol_addr)?;
        owned_root_entries.push((ref_name.to_owned(), ref_hdr));
        drop(ref_seqs);
        drop(ref_seqname_bytes);
        drop(ref_bot_offsets);
        drop(bot_starts);
        drop(bot_lens);

        // Stream each query leaf genome directly from its contigs and free its top segments.
        for (j, (q_name, q_contigs)) in queries.iter().enumerate() {
            let q_seqs = std::mem::take(&mut qry_seqs[j]);
            let tops = std::mem::take(&mut q_tops[j]);
            let g_hdr =
                write_query_leaf_genome(&mut w, q_contigs, qry_lens[j], &q_seqs, &tops, gcol_addr)?;
            owned_root_entries.push((q_name.clone(), g_hdr));
        }

        let root_refs: Vec<(&str, u64)> = owned_root_entries
            .iter()
            .map(|(name, addr)| (name.as_str(), *addr))
            .collect();
        let (root_hdr, root_bt, root_hp) = w.write_group(&root_refs)?;
        let _ = w.finish(root_hdr, root_bt, root_hp)?;
        Ok((tree_newick, n_b))
    }
}

#[allow(clippy::too_many_arguments)]
fn write_progressive_ancestor_genome<W: Write + Seek>(
    w: &mut Hdf5V0Builder<W>,
    meta: &RefStreamMeta<'_>,
    q_tops: &[CompactTop],
    is_root: bool,
    gcol_addr: u64,
    dna_hdr_addr: u64,
    right_links: &mut [(i64, u8)],
) -> std::io::Result<u64> {
    let n_b = meta.bot_starts.len();
    right_links.fill((-1, 0));
    for (t_idx, top) in q_tops.iter().enumerate() {
        if !top.is_chain_head {
            continue;
        }
        let bi = top.parent_bot as usize;
        right_links[bi] = (i64::try_from(t_idx).unwrap_or(-1), u8::from(top.rev));
    }

    let n_tops = if is_root { 0 } else { n_b };
    let dna = dna_hdr_addr;
    let seqidx = w.write_dataset_stream(&[(meta.seqs.len() + 1) as u64], 24, |out| {
        let mut written = 0u64;
        for (s, &bot_off) in meta.seqs.iter().zip(meta.bot_offsets.iter()) {
            let top_off = if is_root { 0 } else { bot_off };
            out.write_all(&(s.offset as u64).to_le_bytes())?;
            out.write_all(&top_off.to_le_bytes())?;
            out.write_all(&bot_off.to_le_bytes())?;
            written += 24;
        }
        out.write_all(&(meta.dna_len as u64).to_le_bytes())?;
        out.write_all(&(n_tops as u64).to_le_bytes())?;
        out.write_all(&(n_b as u64).to_le_bytes())?;
        written += 24;
        Ok(written)
    })?;
    let seqname = w.write_dataset(
        &[meta.seqs.len() as u64],
        meta.seqname_stride,
        meta.seqname_bytes,
    )?;
    let top = w.write_dataset_stream(&[(n_tops + 1) as u64], 33, |out| {
        let mut written = 0u64;
        if !is_root {
            for (bi, &st) in meta.bot_starts.iter().enumerate() {
                let start = i64::try_from(st).unwrap_or(0);
                let bi_i64 = i64::try_from(bi).unwrap_or(-1);
                write_top_record(out, start, bi_i64, -1, bi_i64, 0)?;
                written += 33;
            }
        }
        write_top_sentinel(out, meta.dna_len)?;
        written += 33;
        Ok(written)
    })?;
    let bot = w.write_dataset_stream(&[(n_b + 1) as u64], 24 + 9 * 2, |out| {
        let mut written = 0u64;
        let mut rec = [0u8; 42];
        for bi in 0..n_b {
            let start = i64::try_from(meta.bot_starts[bi]).unwrap_or(0);
            let left_i = i64::try_from(bi).unwrap_or(-1);
            let top_parse = if is_root { -1i64 } else { left_i };
            let (right_i, right_rev) = right_links[bi];
            rec[0..8].copy_from_slice(&start.to_le_bytes());
            rec[8..16].copy_from_slice(&meta.bot_lens[bi].to_le_bytes());
            rec[16..24].copy_from_slice(&top_parse.to_le_bytes());
            rec[24..32].copy_from_slice(&left_i.to_le_bytes());
            rec[32] = 0;
            rec[33..41].copy_from_slice(&right_i.to_le_bytes());
            rec[41] = right_rev;
            out.write_all(&rec)?;
            written += 42;
        }
        write_bottom_sentinel(out, 2, meta.dna_len)?;
        written += 42;
        Ok(written)
    })?;
    write_hal_genome_group(
        w,
        &HalGenomeDatasets {
            dna,
            seqidx,
            seqname,
            top,
            bot,
        },
        gcol_addr,
        meta.rup_idx,
    )
}

fn write_ref_leaf_genome<W: Write + Seek>(
    w: &mut Hdf5V0Builder<W>,
    ref_contigs: &[(String, Vec<u8>)],
    meta: &RefStreamMeta<'_>,
    gcol_addr: u64,
) -> std::io::Result<u64> {
    let n_b = meta.bot_starts.len();
    let dna = w.write_dataset_stream(&[meta.dna_packed_len], 1, |out| {
        write_hal_dna_contigs_stream(ref_contigs, out)
    })?;
    let seqidx = w.write_dataset_stream(&[(meta.seqs.len() + 1) as u64], 24, |out| {
        let mut written = 0u64;
        for (s, &top_off) in meta.seqs.iter().zip(meta.bot_offsets.iter()) {
            out.write_all(&(s.offset as u64).to_le_bytes())?;
            out.write_all(&top_off.to_le_bytes())?;
            out.write_all(&0u64.to_le_bytes())?;
            written += 24;
        }
        out.write_all(&(meta.dna_len as u64).to_le_bytes())?;
        out.write_all(&(n_b as u64).to_le_bytes())?;
        out.write_all(&0u64.to_le_bytes())?;
        written += 24;
        Ok(written)
    })?;
    let seqname = w.write_dataset(
        &[meta.seqs.len() as u64],
        meta.seqname_stride,
        meta.seqname_bytes,
    )?;
    let top = w.write_dataset_stream(&[(n_b + 1) as u64], 33, |out| {
        let mut written = 0u64;
        for (bi, &st) in meta.bot_starts.iter().enumerate() {
            let start = i64::try_from(st).unwrap_or(0);
            let parent = i64::try_from(bi).unwrap_or(-1);
            write_top_record(out, start, -1, -1, parent, 0)?;
            written += 33;
        }
        write_top_sentinel(out, meta.dna_len)?;
        written += 33;
        Ok(written)
    })?;
    let bot = w.write_dataset_stream(&[1], 24, |out| {
        write_bottom_sentinel(out, 0, meta.dna_len)?;
        Ok(24)
    })?;
    write_hal_genome_group(
        w,
        &HalGenomeDatasets {
            dna,
            seqidx,
            seqname,
            top,
            bot,
        },
        gcol_addr,
        meta.rup_idx,
    )
}

fn write_query_leaf_genome<W: Write + Seek>(
    w: &mut Hdf5V0Builder<W>,
    q_contigs: &[(String, Vec<u8>)],
    q_dna_len: usize,
    q_seqs: &[HalSequence],
    tops: &[CompactTop],
    gcol_addr: u64,
) -> std::io::Result<u64> {
    let dna = w.write_dataset_stream(&[q_dna_len.div_ceil(2) as u64], 1, |out| {
        write_hal_dna_contigs_stream(q_contigs, out)
    })?;
    let seqidx = w.write_dataset_stream(&[(q_seqs.len() + 1) as u64], 24, |out| {
        let mut written = 0u64;
        for s in q_seqs {
            let top_off = tops.partition_point(|t| t.start < s.offset as u64) as u64;
            out.write_all(&(s.offset as u64).to_le_bytes())?;
            out.write_all(&top_off.to_le_bytes())?;
            out.write_all(&0u64.to_le_bytes())?;
            written += 24;
        }
        out.write_all(&(q_dna_len as u64).to_le_bytes())?;
        out.write_all(&(tops.len() as u64).to_le_bytes())?;
        out.write_all(&0u64.to_le_bytes())?;
        written += 24;
        Ok(written)
    })?;
    let (seqname_bytes, stride) = encode_hal_seq_names(q_seqs)?;
    let seqname = w.write_dataset(&[q_seqs.len() as u64], stride, &seqname_bytes)?;
    let top = w.write_dataset_stream(&[(tops.len() + 1) as u64], 33, |out| {
        let mut written = 0u64;
        for t in tops {
            let start = i64::try_from(t.start).unwrap_or(0);
            let paralogy = if t.next_paralog == NO_PARALOG {
                -1i64
            } else {
                i64::from(t.next_paralog)
            };
            let parent = if t.parent_bot == NO_PARALOG {
                -1i64
            } else {
                i64::from(t.parent_bot)
            };
            write_top_record(out, start, -1, paralogy, parent, u8::from(t.rev))?;
            written += 33;
        }
        write_top_sentinel(out, q_dna_len)?;
        written += 33;
        Ok(written)
    })?;
    let bot = w.write_dataset_stream(&[1], 24, |out| {
        write_bottom_sentinel(out, 0, q_dna_len)?;
        Ok(24)
    })?;
    let rup_idx = if q_dna_len % 2 == 1 { 4 } else { 3 };
    write_hal_genome_group(
        w,
        &HalGenomeDatasets {
            dna,
            seqidx,
            seqname,
            top,
            bot,
        },
        gcol_addr,
        rup_idx,
    )
}
