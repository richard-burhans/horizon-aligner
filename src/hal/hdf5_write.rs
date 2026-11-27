//! Pure-Rust HDF5 v0 binary HAL exporter for in-memory `HalAlignment` models.

use super::hdf5_writer::Hdf5V0Builder;
use super::normalize::normalize_hal_segments_for_export;
use super::{HalAlignment, HalBottomSegment, HalGenome, HalSequence, HalTopSegment};
use std::io::{BufWriter, Cursor, Seek, Write};

pub(super) fn encode_hal_nibble(base: u8) -> u8 {
    match base {
        b'a' => 0,
        b'c' => 1,
        b'g' => 2,
        b't' => 3,
        b'n' => 4,
        b'A' => 8,
        b'C' => 9,
        b'G' => 10,
        b'T' => 11,
        _ => 12,
    }
}

pub(super) fn encode_hal_seq_names(seqs: &[HalSequence]) -> std::io::Result<(Vec<u8>, usize)> {
    const MAX_SEQNAME_TABLE_BYTES: usize = 256 * 1024 * 1024;
    let stride = seqs
        .iter()
        .map(|s| s.name.len().saturating_add(1))
        .max()
        .unwrap_or(1);
    let total = seqs
        .len()
        .checked_mul(stride)
        .filter(|&n| n <= MAX_SEQNAME_TABLE_BYTES)
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "HAL SEQNAME_ARRAY size ({} x {stride}) exceeds {MAX_SEQNAME_TABLE_BYTES}-byte limit",
                    seqs.len()
                ),
            )
        })?;
    let mut out = vec![0u8; total];
    for (i, s) in seqs.iter().enumerate() {
        let bytes = s.name.as_bytes();
        let copy_len = bytes.len().min(stride - 1);
        out[i * stride..i * stride + copy_len].copy_from_slice(&bytes[..copy_len]);
    }
    Ok((out, stride))
}

fn write_hal_dna_stream<W: Write>(dna: &[u8], w: &mut W) -> std::io::Result<u64> {
    let mut buf = [0u8; 8192];
    let mut written = 0u64;
    let mut pos = 0usize;
    for chunk in dna.chunks(2) {
        let hi = encode_hal_nibble(chunk[0]);
        let lo = chunk.get(1).copied().map_or(0, encode_hal_nibble);
        buf[pos] = (hi << 4) | lo;
        pos += 1;
        if pos == buf.len() {
            w.write_all(&buf)?;
            written += pos as u64;
            pos = 0;
        }
    }
    if pos > 0 {
        w.write_all(&buf[..pos])?;
        written += pos as u64;
    }
    Ok(written)
}

fn is_non_star_tree(hal: &HalAlignment) -> bool {
    hal.genomes.first().is_some_and(|root| {
        root.parent.is_some()
            || hal.genomes[1..]
                .iter()
                .any(|g| g.parent != Some(0) || !g.children.is_empty())
    })
}

fn build_sort_permutation<T, F: Fn(&T) -> usize>(
    items: &[T],
    key: F,
) -> Option<(Vec<usize>, Vec<usize>)> {
    if items.windows(2).all(|w| key(&w[0]) <= key(&w[1])) {
        return None;
    }
    let mut ord: Vec<usize> = (0..items.len()).collect();
    ord.sort_by_key(|&i| key(&items[i]));
    let mut inv = vec![0usize; items.len()];
    for (new_i, &old_i) in ord.iter().enumerate() {
        inv[old_i] = new_i;
    }
    Some((ord, inv))
}

pub(super) fn write_top_record<W: Write>(
    w: &mut W,
    start: i64,
    bottom_parse: i64,
    paralogy: i64,
    parent: i64,
    rev: u8,
) -> std::io::Result<()> {
    let mut rec = [0u8; 33];
    rec[0..8].copy_from_slice(&start.to_le_bytes());
    rec[8..16].copy_from_slice(&bottom_parse.to_le_bytes());
    rec[16..24].copy_from_slice(&paralogy.to_le_bytes());
    rec[24..32].copy_from_slice(&parent.to_le_bytes());
    rec[32] = rev;
    w.write_all(&rec)
}

pub(super) fn write_top_sentinel<W: Write>(w: &mut W, total_dna_len: usize) -> std::io::Result<()> {
    let end_pos = i64::try_from(total_dna_len).unwrap_or(0);
    write_top_record(w, end_pos, -1, -1, -1, 0)
}

pub(super) fn write_bottom_sentinel<W: Write>(
    w: &mut W,
    num_children: usize,
    total_dna_len: usize,
) -> std::io::Result<()> {
    let end_pos = i64::try_from(total_dna_len).unwrap_or(0);
    w.write_all(&end_pos.to_le_bytes())?;
    w.write_all(&0u64.to_le_bytes())?;
    w.write_all(&(-1i64).to_le_bytes())?;
    for _ in 0..num_children {
        w.write_all(&(-1i64).to_le_bytes())?;
        w.write_all(&[0u8])?;
    }
    Ok(())
}

pub(super) struct HalGenomeDatasets {
    pub(super) dna: u64,
    pub(super) seqidx: u64,
    pub(super) seqname: u64,
    pub(super) top: u64,
    pub(super) bot: u64,
}

pub(super) fn init_hdf5_root<W: Write + Seek>(
    w: &mut Hdf5V0Builder<W>,
    tree_newick: &str,
    n_genomes: usize,
) -> std::io::Result<(u64, Vec<(String, u64)>)> {
    // Store VLS attribute strings in a single Global Heap collection:
    // idx 1 = tree_newick, idx 2 = "2.2", idx 3 = "0", idx 4 = "1"
    let gcol_addr = w.write_gcol(&[tree_newick, "2.2", "0", "1"])?;
    let tree_len = u32::try_from(tree_newick.len()).unwrap_or(u32::MAX);

    let (meta_hdr, _, _) = w.write_group(&[])?;
    let (genomes_hdr, _, _) = w.write_group(&[])?;
    let (phylo_hdr, _, _) =
        w.write_group_with_attrs(&[], &[("Phylogeny", tree_len, gcol_addr, 1)])?;
    // Intentional historic spelling required by UCSC libHal ("Verison", not "Version").
    let (ver_hdr, _, _) = w.write_group_with_attrs(&[], &[("Verison", 3, gcol_addr, 2)])?;

    let mut owned_root_entries: Vec<(String, u64)> = Vec::with_capacity(n_genomes + 4);
    owned_root_entries.push(("Meta".to_owned(), meta_hdr));
    owned_root_entries.push(("Phylogeny".to_owned(), phylo_hdr));
    owned_root_entries.push(("Genomes".to_owned(), genomes_hdr));
    owned_root_entries.push(("Verison".to_owned(), ver_hdr));
    Ok((gcol_addr, owned_root_entries))
}

pub(super) fn write_hal_genome_group<W: Write + Seek>(
    w: &mut Hdf5V0Builder<W>,
    ds: &HalGenomeDatasets,
    gcol_addr: u64,
    rup_idx: u32,
) -> std::io::Result<u64> {
    let (g_meta_hdr, _, _) = w.write_group(&[])?;
    let (g_rup_hdr, _, _) = w.write_group_with_attrs(&[], &[("Rup", 1, gcol_addr, rup_idx)])?;
    let (g_hdr, _, _) = w.write_group(&[
        ("DNA_ARRAY", ds.dna),
        ("SEQIDX_ARRAY", ds.seqidx),
        ("SEQNAME_ARRAY", ds.seqname),
        ("TOP_ARRAY", ds.top),
        ("BOTTOM_ARRAY", ds.bot),
        ("Meta", g_meta_hdr),
        ("Rup", g_rup_hdr),
    ])?;
    Ok(g_hdr)
}

fn find_covering_seg<T, FStart, FLen>(
    segs: &[T],
    ord: Option<&[usize]>,
    pos: usize,
    start_of: FStart,
    len_of: FLen,
) -> i64
where
    FStart: Fn(&T) -> usize,
    FLen: Fn(&T) -> usize,
{
    let idx_opt = if let Some(perm) = ord {
        perm.partition_point(|&i| start_of(&segs[i]) <= pos)
            .checked_sub(1)
    } else {
        segs.partition_point(|s| start_of(s) <= pos).checked_sub(1)
    };
    let Some(sorted_pos) = idx_opt else {
        return -1;
    };
    let raw_idx = ord.map_or(sorted_pos, |perm| perm[sorted_pos]);
    let seg = &segs[raw_idx];
    if start_of(seg).saturating_add(len_of(seg)) > pos {
        i64::try_from(sorted_pos).unwrap_or(-1)
    } else {
        -1
    }
}

struct GenomeExportView<'a> {
    genome: &'a HalGenome,
    tops: &'a [HalTopSegment],
    bottoms: &'a [HalBottomSegment],
    top_ord: Option<&'a [usize]>,
    top_inv: Option<&'a [usize]>,
    bot_ord: Option<&'a [usize]>,
    parent_bot_inv: Option<&'a [usize]>,
    top_perms: &'a [Option<(Vec<usize>, Vec<usize>)>],
}

#[allow(clippy::too_many_lines)]
fn write_hal_genome_datasets<W: Write + Seek>(
    w: &mut Hdf5V0Builder<W>,
    view: &GenomeExportView<'_>,
    gcol_addr: u64,
) -> std::io::Result<u64> {
    let g = view.genome;
    let tops = view.tops;
    let bottoms = view.bottoms;
    let top_ord = view.top_ord;
    let top_inv = view.top_inv;
    let bot_ord = view.bot_ord;

    let dna_packed_len = g.dna.len().div_ceil(2) as u64;
    let dna = w.write_dataset_stream(&[dna_packed_len], 1, |out| {
        write_hal_dna_stream(&g.dna, out)
    })?;

    let seqidx = w.write_dataset_stream(&[(g.sequences.len() + 1) as u64], 24, |out| {
        let mut written = 0u64;
        for s in &g.sequences {
            let top_off = if let Some(ord) = top_ord {
                ord.partition_point(|&old_i| tops[old_i].start < s.offset)
            } else {
                tops.partition_point(|t| t.start < s.offset)
            } as u64;
            let bot_off = if let Some(ord) = bot_ord {
                ord.partition_point(|&old_i| bottoms[old_i].start < s.offset)
            } else {
                bottoms.partition_point(|b| b.start < s.offset)
            } as u64;
            out.write_all(&(s.offset as u64).to_le_bytes())?;
            out.write_all(&top_off.to_le_bytes())?;
            out.write_all(&bot_off.to_le_bytes())?;
            written += 24;
        }
        out.write_all(&(g.dna.len() as u64).to_le_bytes())?;
        out.write_all(&(tops.len() as u64).to_le_bytes())?;
        out.write_all(&(bottoms.len() as u64).to_le_bytes())?;
        written += 24;
        Ok(written)
    })?;

    let (seqname_bytes, stride) = encode_hal_seq_names(&g.sequences)?;
    let seqname = w.write_dataset(&[g.sequences.len() as u64], stride, &seqname_bytes)?;

    let top = w.write_dataset_stream(&[(tops.len() + 1) as u64], 33, |out| {
        let mut written = 0u64;
        for pos in 0..tops.len() {
            let old_i = top_ord.map_or(pos, |ord| ord[pos]);
            let t = &tops[old_i];
            let parent_remapped = t.parent_bottom_idx.and_then(|pb| {
                view.parent_bot_inv
                    .map_or(Some(pb), |inv| inv.get(pb).copied())
            });
            let paralogy_remapped = t
                .next_paralog
                .and_then(|np| top_inv.map_or(Some(np), |inv| inv.get(np).copied()));
            let parent = parent_remapped.map_or(-1i64, |v| i64::try_from(v).unwrap_or(-1));
            let paralogy = paralogy_remapped.map_or(-1i64, |v| i64::try_from(v).unwrap_or(-1));
            let bottom_parse =
                find_covering_seg(bottoms, bot_ord, t.start, |b| b.start, |b| b.length);
            let start = i64::try_from(t.start).unwrap_or(0);
            write_top_record(
                out,
                start,
                bottom_parse,
                paralogy,
                parent,
                u8::from(t.parent_reversed),
            )?;
            written += 33;
        }
        write_top_sentinel(out, g.dna.len())?;
        written += 33;
        Ok(written)
    })?;

    let num_children = g.children.len();
    let rec_size = 24 + 9 * num_children;
    let bot = w.write_dataset_stream(&[(bottoms.len() + 1) as u64], rec_size, |out| {
        let mut written = 0u64;
        for pos in 0..bottoms.len() {
            let old_i = bot_ord.map_or(pos, |ord| ord[pos]);
            let b = &bottoms[old_i];
            let start = i64::try_from(b.start).unwrap_or(0);
            let top_parse = find_covering_seg(tops, top_ord, b.start, |t| t.start, |t| t.length);
            out.write_all(&start.to_le_bytes())?;
            out.write_all(&(b.length as u64).to_le_bytes())?;
            out.write_all(&top_parse.to_le_bytes())?;
            for c_idx in 0..num_children {
                let link = b.child_tops.get(c_idx).and_then(|opt| *opt);
                let remapped = link.and_then(|(old_t, rev)| {
                    let child_inv = g
                        .children
                        .get(c_idx)
                        .and_then(|&cg| view.top_perms.get(cg))
                        .and_then(|p| p.as_ref())
                        .map(|(_, inv)| inv.as_slice());
                    let new_t = child_inv.map_or(Some(old_t), |inv| inv.get(old_t).copied())?;
                    Some((new_t, rev))
                });
                let (child_i, rev) = remapped.map_or((-1i64, 0u8), |(idx, r)| {
                    (i64::try_from(idx).unwrap_or(-1), u8::from(r))
                });
                out.write_all(&child_i.to_le_bytes())?;
                out.write_all(&[rev])?;
            }
            written += rec_size as u64;
        }
        write_bottom_sentinel(out, num_children, g.dna.len())?;
        written += rec_size as u64;
        Ok(written)
    })?;

    let rup_idx = if g.dna.len() % 2 == 1 { 4 } else { 3 };
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

impl HalAlignment {
    fn write_hdf5_to<W: Write + Seek>(&self, writer: W) -> std::io::Result<W> {
        let mut w = Hdf5V0Builder::new(writer)?;
        let (gcol_addr, mut owned_root_entries) =
            init_hdf5_root(&mut w, &self.tree_newick, self.genomes.len())?;

        if is_non_star_tree(self) {
            // Stream non-star trees using only lightweight permutation tables for unsorted
            // genomes, avoiding cloning all HalTopSegment/HalBottomSegment tables in RAM.
            let top_perms: Vec<Option<(Vec<usize>, Vec<usize>)>> = self
                .genomes
                .iter()
                .map(|g| build_sort_permutation(&g.top_segments, |t| t.start))
                .collect();
            let bot_perms: Vec<Option<(Vec<usize>, Vec<usize>)>> = self
                .genomes
                .iter()
                .map(|g| build_sort_permutation(&g.bottom_segments, |b| b.start))
                .collect();
            for (g_idx, g) in self.genomes.iter().enumerate() {
                let view = GenomeExportView {
                    genome: g,
                    tops: &g.top_segments,
                    bottoms: &g.bottom_segments,
                    top_ord: top_perms[g_idx].as_ref().map(|(ord, _)| ord.as_slice()),
                    top_inv: top_perms[g_idx].as_ref().map(|(_, inv)| inv.as_slice()),
                    bot_ord: bot_perms[g_idx].as_ref().map(|(ord, _)| ord.as_slice()),
                    parent_bot_inv: g
                        .parent
                        .and_then(|pg| bot_perms.get(pg))
                        .and_then(|p| p.as_ref())
                        .map(|(_, inv)| inv.as_slice()),
                    top_perms: &top_perms,
                };
                let g_hdr = write_hal_genome_datasets(&mut w, &view, gcol_addr)?;
                owned_root_entries.push((g.name.clone(), g_hdr));
            }
        } else {
            let norm = normalize_hal_segments_for_export(self);
            for (g_idx, g) in self.genomes.iter().enumerate() {
                let view = GenomeExportView {
                    genome: g,
                    tops: &norm.tops_by_genome[g_idx],
                    bottoms: &norm.bottoms_by_genome[g_idx],
                    top_ord: None,
                    top_inv: None,
                    bot_ord: None,
                    parent_bot_inv: None,
                    top_perms: &[],
                };
                let g_hdr = write_hal_genome_datasets(&mut w, &view, gcol_addr)?;
                owned_root_entries.push((g.name.clone(), g_hdr));
            }
        }
        let root_refs: Vec<(&str, u64)> = owned_root_entries
            .iter()
            .map(|(name, addr)| (name.as_str(), *addr))
            .collect();
        let (root_hdr, root_bt, root_hp) = w.write_group(&root_refs)?;
        w.finish(root_hdr, root_bt, root_hp)
    }

    /// Serialize this `HalAlignment` into a valid, self-contained HDF5 `v0` `.hal` byte stream
    /// using pure Rust without linking `libhdf5`.
    #[must_use]
    pub fn to_hdf5_bytes(&self) -> Vec<u8> {
        self.write_hdf5_to(Cursor::new(Vec::new()))
            .map_or_default(Cursor::into_inner)
    }

    /// Write this `HalAlignment` to disk as an HDF5 `v0` `.hal` file, streaming directly
    /// through a buffered file writer without materializing the full HDF5 image in memory.
    ///
    /// # Errors
    /// Returns an `std::io::Error` if creating or writing the file fails.
    pub fn write_hdf5_hal(&self, path: &std::path::Path) -> std::io::Result<()> {
        let file = std::fs::File::create(path)?;
        let writer = BufWriter::with_capacity(256 * 1024, file);
        let _ = self.write_hdf5_to(writer)?;
        Ok(())
    }
}
