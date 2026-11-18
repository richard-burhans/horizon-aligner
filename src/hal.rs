//! HAL-compatible hierarchical alignment model (`TopSegment`/`BottomSegment` + guide tree).
//!
//! Models the HAL graph schema (Hickey et al. 2013) in pure Rust:
//! - Phylogenetic Newick tree linking leaf and ancestral genomes
//! - Per-genome `TopSegment` (child->parent) and `BottomSegment` (parent->child) coordinate tables
//! - Reference-free projection (`to_maf_for_genome`) onto any genome in the tree
//! - HDF5 container magic-byte inspection (`is_hdf5_hal`) for external `.hal` files

use std::collections::BTreeMap;
use std::fmt::Write as _;

/// HDF5 8-byte file signature used by `.hal` containers.
pub const HDF5_SIGNATURE: [u8; 8] = [0x89, b'H', b'D', b'F', 0x0d, 0x0a, 0x1a, 0x0a];

/// One sequence (chromosome/contig) within a HAL genome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalSequence {
    /// Sequence name.
    pub name: String,
    /// Total length in bases.
    pub length: usize,
    /// Global start offset in the genome's concatenated coordinate space.
    pub offset: usize,
}

/// Child-to-parent colinear segment in HAL.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalTopSegment {
    /// Start coordinate in the child genome's global space.
    pub start: usize,
    /// Aligned span in bases.
    pub length: usize,
    /// Index of the matching `BottomSegment` in the parent genome (`None` if unaligned/insertion).
    pub parent_bottom_idx: Option<usize>,
    /// Whether the child aligns on the reverse strand relative to the parent.
    pub parent_reversed: bool,
    /// Optional index of next paralogous `TopSegment` in the same genome.
    pub next_paralog: Option<usize>,
}

/// Parent-to-child colinear segment in HAL.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalBottomSegment {
    /// Start coordinate in the parent genome's global space.
    pub start: usize,
    /// Aligned span in bases.
    pub length: usize,
    /// For each child genome index: `Some((top_seg_idx, reversed))` or `None` if deleted/absent.
    pub child_tops: Vec<Option<(usize, bool)>>,
}

/// One genome (leaf or ancestral) in a HAL alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalGenome {
    /// Species/genome name.
    pub name: String,
    /// Parent genome index (`None` for root).
    pub parent: Option<usize>,
    /// Child genome indices.
    pub children: Vec<usize>,
    /// Sequences belonging to this genome.
    pub sequences: Vec<HalSequence>,
    /// Concatenated DNA sequence (uppercase ASCII).
    pub dna: Vec<u8>,
    /// Child->parent segments.
    pub top_segments: Vec<HalTopSegment>,
    /// Parent->child segments.
    pub bottom_segments: Vec<HalBottomSegment>,
}

/// Complete hierarchical HAL alignment graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalAlignment {
    /// Newick guide tree.
    pub tree_newick: String,
    /// All genomes (root + internal + leaves).
    pub genomes: Vec<HalGenome>,
}

/// Check whether a byte slice begins with the standard HDF5 `.hal` container signature.
#[must_use]
pub fn is_hdf5_hal(bytes: &[u8]) -> bool {
    bytes.len() >= HDF5_SIGNATURE.len() && bytes[..8] == HDF5_SIGNATURE
}

impl HalAlignment {
    /// Construct a star-tree `HalAlignment` from a reference genome and pairwise `(species, seq_name, ref_start, qry_start, len, rev, ref_bases, qry_bases)` blocks.
    #[must_use]
    pub fn from_star_blocks(
        ref_name: &str,
        ref_dna: &[u8],
        queries: &[(String, Vec<u8>)],
        blocks: &[(usize, usize, usize, usize, bool)], // (query_idx, ref_start, qry_start, len, rev)
    ) -> Self {
        type HalIntervalHits = Vec<(usize, usize, bool)>;
        let mut genomes = Vec::with_capacity(1 + queries.len());
        let mut bottom_segments = Vec::new();

        // Group blocks by (ref_start, len) on the root/reference genome.
        let mut by_interval: BTreeMap<(usize, usize), HalIntervalHits> = BTreeMap::new();
        for &(q_idx, r_start, q_start, len, rev) in blocks {
            by_interval
                .entry((r_start, len))
                .or_default()
                .push((q_idx, q_start, rev));
        }

        let mut child_tops_per_q: Vec<Vec<HalTopSegment>> = vec![Vec::new(); queries.len()];
        for ((r_start, len), hits) in by_interval {
            let b_idx = bottom_segments.len();
            let mut child_links = vec![None; queries.len()];
            for (q_idx, q_start, rev) in hits {
                let t_idx = child_tops_per_q[q_idx].len();
                child_tops_per_q[q_idx].push(HalTopSegment {
                    start: q_start,
                    length: len,
                    parent_bottom_idx: Some(b_idx),
                    parent_reversed: rev,
                    next_paralog: None,
                });
                child_links[q_idx] = Some((t_idx, rev));
            }
            bottom_segments.push(HalBottomSegment {
                start: r_start,
                length: len,
                child_tops: child_links,
            });
        }

        genomes.push(HalGenome {
            name: ref_name.to_owned(),
            parent: None,
            children: (1..=queries.len()).collect(),
            sequences: vec![HalSequence {
                name: format!("{ref_name}.chr1"),
                length: ref_dna.len(),
                offset: 0,
            }],
            dna: ref_dna.to_vec(),
            top_segments: Vec::new(),
            bottom_segments,
        });

        let mut leaf_names = Vec::with_capacity(queries.len());
        for (q_idx, (q_name, q_dna)) in queries.iter().enumerate() {
            leaf_names.push(q_name.clone());
            genomes.push(HalGenome {
                name: q_name.clone(),
                parent: Some(0),
                children: Vec::new(),
                sequences: vec![HalSequence {
                    name: format!("{q_name}.chr1"),
                    length: q_dna.len(),
                    offset: 0,
                }],
                dna: q_dna.clone(),
                top_segments: std::mem::take(&mut child_tops_per_q[q_idx]),
                bottom_segments: Vec::new(),
            });
        }

        let tree_newick = format!("({}){ref_name};", leaf_names.join(","));
        Self {
            tree_newick,
            genomes,
        }
    }

    /// Project the HAL alignment onto `target_genome_name` as MAF v1 text.
    #[must_use]
    pub fn to_maf_for_genome(&self, target_genome_name: &str) -> Option<String> {
        let g_idx = self
            .genomes
            .iter()
            .position(|g| g.name == target_genome_name)?;
        let root = &self.genomes[g_idx];
        let mut out = String::from("##maf version=1 scoring=horizon_hal\n");

        for bseg in &root.bottom_segments {
            let r_end = (bseg.start + bseg.length).min(root.dna.len());
            if r_end <= bseg.start {
                continue;
            }
            let r_seq = String::from_utf8_lossy(&root.dna[bseg.start..r_end]);
            let active_children = bseg.child_tops.iter().filter(|c| c.is_some()).count();
            let _ = writeln!(out, "a score={}.0", active_children + 1);
            let _ = writeln!(
                out,
                "s {} {} {} + {} {}",
                root.name,
                bseg.start,
                r_end - bseg.start,
                root.dna.len(),
                r_seq
            );
            for (c_pos, link) in bseg.child_tops.iter().enumerate() {
                let Some((t_idx, rev)) = *link else {
                    continue;
                };
                let child_idx = root.children[c_pos];
                let child = &self.genomes[child_idx];
                let tseg = &child.top_segments[t_idx];
                let c_end = (tseg.start + tseg.length).min(child.dna.len());
                if c_end <= tseg.start {
                    continue;
                }
                let slice = &child.dna[tseg.start..c_end];
                let q_seq = if rev {
                    String::from_utf8_lossy(&crate::revcomp(slice)).into_owned()
                } else {
                    String::from_utf8_lossy(slice).into_owned()
                };
                let strand = if rev { '-' } else { '+' };
                let q_start = if rev {
                    child.dna.len().saturating_sub(c_end)
                } else {
                    tseg.start
                };
                let _ = writeln!(
                    out,
                    "s {} {} {} {} {} {}",
                    child.name,
                    q_start,
                    c_end - tseg.start,
                    strand,
                    child.dna.len(),
                    q_seq
                );
            }
            out.push('\n');
        }
        Some(out)
    }
}

/// Summary metadata parsed from an HDF5 `.hal` file superblock (v0/v1/v2/v3).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hdf5Superblock {
    /// Superblock format version (`0`, `1`, `2`, or `3`).
    pub version: u8,
    /// Size in bytes of file offsets (`sizeof_addr`).
    pub offset_size: u8,
    /// Size in bytes of object lengths (`sizeof_len`).
    pub length_size: u8,
    /// Root group object header address when `offset_size == 8` and present in the header.
    pub root_group_addr: Option<u64>,
}

/// Parse the HDF5 superblock header from the leading bytes of a `.hal` file.
#[must_use]
pub fn read_hdf5_superblock(bytes: &[u8]) -> Option<Hdf5Superblock> {
    if bytes.len() < 16 || !is_hdf5_hal(bytes) {
        return None;
    }
    let version = bytes[8];
    match version {
        0 | 1 => {
            let offset_size = bytes[13];
            let length_size = bytes[14];
            let root_group_addr = if offset_size == 8 && length_size == 8 && bytes.len() >= 72 {
                Some(u64::from_le_bytes(bytes[64..72].try_into().ok()?))
            } else {
                None
            };
            Some(Hdf5Superblock {
                version,
                offset_size,
                length_size,
                root_group_addr,
            })
        }
        2 | 3 => {
            let offset_size = bytes[9];
            let length_size = bytes[10];
            let root_group_addr = if offset_size == 8 && bytes.len() >= 44 {
                Some(u64::from_le_bytes(bytes[36..44].try_into().ok()?))
            } else {
                None
            };
            Some(Hdf5Superblock {
                version,
                offset_size,
                length_size,
                root_group_addr,
            })
        }
        _ => None,
    }
}

/// Inspect a `.hal` file on disk and parse its HDF5 superblock header.
///
/// # Errors
/// Returns an `std::io::Error` if the file cannot be opened or read.
pub fn inspect_hal_file(path: &std::path::Path) -> std::io::Result<Option<Hdf5Superblock>> {
    use std::io::Read as _;
    let mut f = std::fs::File::open(path)?;
    let mut buf = [0u8; 96];
    let n = f.read(&mut buf)?;
    Ok(read_hdf5_superblock(&buf[..n]))
}

/// Summary of an HDF5 object header (v1 or v2) parsed from a `.hal` file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Hdf5ObjectHeader {
    /// Object header format version (`1` or `2`).
    pub version: u8,
    /// Message type IDs found in the header prefix (e.g., `0x0011` `SymbolTable`, `0x000C` Attribute).
    pub message_types: Vec<u16>,
}

fn iter_v1_header_messages<'a, F>(bytes: &'a [u8], obj_addr: usize, mut f: F) -> Option<()>
where
    F: FnMut(u16, &'a [u8]) -> std::ops::ControlFlow<()>,
{
    let slice = bytes.get(obj_addr..)?;
    if slice.len() < 16 || slice[0] != 1 {
        return None;
    }
    let num_msgs = usize::from(u16::from_le_bytes([slice[2], slice[3]]));
    let hdr_size = usize::try_from(u32::from_le_bytes([
        slice[8], slice[9], slice[10], slice[11],
    ]))
    .ok()?;
    let first_end = obj_addr
        .saturating_add(16)
        .saturating_add(hdr_size)
        .min(bytes.len());
    let mut chunks = std::collections::VecDeque::from([(obj_addr + 16, first_end)]);
    let mut visited_addrs = std::collections::HashSet::from([obj_addr]);
    let mut seen = 0usize;

    while let Some((mut pos, end)) = chunks.pop_front() {
        while seen < num_msgs && pos + 8 <= end {
            let mtype = u16::from_le_bytes([bytes[pos], bytes[pos + 1]]);
            let msize = usize::from(u16::from_le_bytes([bytes[pos + 2], bytes[pos + 3]]));
            let body_start = pos + 8;
            let body_end = body_start.saturating_add(msize).min(end);
            let body = &bytes[body_start..body_end];
            seen += 1;
            pos += 8 + msize;

            if mtype == 0x0010 && body.len() >= 16 {
                if let (Some(cont_u64), Some(len_u64)) =
                    (read_u64_le(body, 0), read_u64_le(body, 8))
                    && let (Ok(cont_addr), Ok(cont_len)) =
                        (usize::try_from(cont_u64), usize::try_from(len_u64))
                    && visited_addrs.len() < 1024
                    && visited_addrs.insert(cont_addr)
                    && cont_addr < bytes.len()
                {
                    let cont_end = cont_addr.saturating_add(cont_len).min(bytes.len());
                    chunks.push_back((cont_addr, cont_end));
                }
                continue;
            }

            if f(mtype, body).is_break() {
                return Some(());
            }
        }
    }
    Some(())
}

/// Parse an HDF5 v1 or v2 object header starting at `slice` (positioned at `root_group_addr`).
#[must_use]
pub fn read_hdf5_object_header(slice: &[u8]) -> Option<Hdf5ObjectHeader> {
    if slice.len() < 16 {
        return None;
    }
    if slice.starts_with(b"OHDR") {
        // HDF5 Object Header v2
        let flags = slice[5];
        let mut pos = 6usize;
        if flags & 0x20 != 0 {
            pos += 16; // access, mod, change, birth times (4 x u32)
        }
        if flags & 0x10 != 0 {
            pos += 4; // max_compact, min_dense
        }
        let chunk_sz_bytes = 1usize << (flags & 0x03);
        if pos + chunk_sz_bytes > slice.len() {
            return None;
        }
        let mut chunk_sz = 0usize;
        for b in 0..chunk_sz_bytes {
            chunk_sz |= usize::from(slice[pos + b]) << (8 * b);
        }
        pos += chunk_sz_bytes;
        let end = (pos + chunk_sz).min(slice.len());
        let mut msg_types = Vec::new();
        while pos + 4 <= end {
            let mtype = u16::from(slice[pos]);
            let msize = usize::from(u16::from_le_bytes([slice[pos + 1], slice[pos + 2]]));
            msg_types.push(mtype);
            pos += 4 + msize;
        }
        return Some(Hdf5ObjectHeader {
            version: 2,
            message_types: msg_types,
        });
    }
    if slice[0] == 1 {
        let mut msg_types = Vec::new();
        iter_v1_header_messages(slice, 0, |mtype, _| {
            msg_types.push(mtype);
            std::ops::ControlFlow::Continue(())
        })?;
        return Some(Hdf5ObjectHeader {
            version: 1,
            message_types: msg_types,
        });
    }
    None
}

fn read_u64_le(bytes: &[u8], off: usize) -> Option<u64> {
    Some(u64::from_le_bytes(
        bytes.get(off..off + 8)?.try_into().ok()?,
    ))
}

/// A child entry inside an HDF5 v0 group (`SymbolTable` + `TREE`/`SNOD` + `HEAP`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Hdf5GroupEntry {
    /// Member name in the group.
    pub name: String,
    /// File byte offset of the member's HDF5 object header.
    pub object_header_addr: u64,
}

fn collect_snod_entries(
    bytes: &[u8],
    snod_addr: usize,
    heap_data: &[u8],
    out: &mut Vec<Hdf5GroupEntry>,
) {
    if bytes.get(snod_addr..snod_addr + 4) != Some(b"SNOD") {
        return;
    }
    let Some(n_syms) = bytes
        .get(snod_addr + 6..snod_addr + 8)
        .map(|b| usize::from(u16::from_le_bytes([b[0], b[1]])))
    else {
        return;
    };
    for i in 0..n_syms {
        let entry_off = snod_addr + 8 + i * 40;
        let Some(name_off_u64) = read_u64_le(bytes, entry_off) else {
            break;
        };
        let Some(obj_addr) = read_u64_le(bytes, entry_off + 8) else {
            break;
        };
        let Ok(name_off) = usize::try_from(name_off_u64) else {
            continue;
        };
        if let Some(tail) = heap_data.get(name_off..) {
            let end = tail.iter().position(|&b| b == 0).unwrap_or(tail.len());
            if let Ok(s) = std::str::from_utf8(&tail[..end])
                && !s.is_empty()
            {
                out.push(Hdf5GroupEntry {
                    name: s.to_owned(),
                    object_header_addr: obj_addr,
                });
            }
        }
    }
}

fn walk_btree_v1_group_entries(
    bytes: &[u8],
    node_addr: usize,
    heap_data: &[u8],
    out: &mut Vec<Hdf5GroupEntry>,
) {
    if bytes.get(node_addr..node_addr + 4) != Some(b"TREE") {
        return;
    }
    let node_type = *bytes.get(node_addr + 4).unwrap_or(&255);
    let node_level = *bytes.get(node_addr + 5).unwrap_or(&255);
    let Some(n_entries) = bytes
        .get(node_addr + 6..node_addr + 8)
        .map(|b| usize::from(u16::from_le_bytes([b[0], b[1]])))
    else {
        return;
    };
    if node_type != 0 {
        return;
    }
    // B-tree v1 header with 8-byte offsets: 8 bytes prefix + 2 * 8 bytes sibling addrs = 24 bytes.
    // Entries are alternating (key: u64, child_addr: u64) for `n_entries`.
    for i in 0..n_entries {
        let child_off = node_addr + 24 + i * 16 + 8;
        let Some(child_u64) = read_u64_le(bytes, child_off) else {
            break;
        };
        let Ok(child_addr) = usize::try_from(child_u64) else {
            continue;
        };
        if node_level == 0 {
            collect_snod_entries(bytes, child_addr, heap_data, out);
        } else {
            walk_btree_v1_group_entries(bytes, child_addr, heap_data, out);
        }
    }
}

/// Extract child member entries `(name, object_header_addr)` from any HDF5 v0 group object header.
#[must_use]
pub fn read_hdf5_v0_group_entries(
    bytes: &[u8],
    group_obj_addr: usize,
) -> Option<Vec<Hdf5GroupEntry>> {
    let mut symtab_addrs: Option<(usize, usize)> = None;
    iter_v1_header_messages(bytes, group_obj_addr, |mtype, body| {
        if mtype == 0x0011
            && body.len() >= 16
            && let Some(bt_u64) = read_u64_le(body, 0)
            && let Some(hp_u64) = read_u64_le(body, 8)
            && let (Ok(btree_addr), Ok(heap_addr)) =
                (usize::try_from(bt_u64), usize::try_from(hp_u64))
        {
            symtab_addrs = Some((btree_addr, heap_addr));
            return std::ops::ControlFlow::Break(());
        }
        std::ops::ControlFlow::Continue(())
    })?;
    let (btree_addr, heap_addr) = symtab_addrs?;
    if bytes.get(heap_addr..heap_addr + 4) != Some(b"HEAP") {
        return None;
    }
    let data_seg_size = usize::try_from(read_u64_le(bytes, heap_addr + 8)?).ok()?;
    let data_seg_addr = usize::try_from(read_u64_le(bytes, heap_addr + 24)?).ok()?;
    let heap_data = bytes.get(data_seg_addr..data_seg_addr + data_seg_size)?;
    let mut entries = Vec::new();
    walk_btree_v1_group_entries(bytes, btree_addr, heap_data, &mut entries);
    Some(entries)
}

/// Extract root-group member names from an HDF5 v0 `.hal` file via its `SymbolTable` (`0x0011`),
/// `HEAP` local heap, and `TREE`/`SNOD` B-tree v1 nodes.
#[must_use]
pub fn read_hdf5_v0_root_names(bytes: &[u8]) -> Option<Vec<String>> {
    let sb = read_hdf5_superblock(bytes)?;
    if sb.version > 1 || sb.offset_size != 8 || sb.length_size != 8 {
        return None;
    }
    let root_addr = usize::try_from(sb.root_group_addr?).ok()?;
    let entries = read_hdf5_v0_group_entries(bytes, root_addr)?;
    Some(entries.into_iter().map(|e| e.name).collect())
}

/// Parse HDF5 `Dataspace` (`0x0001`) dimension sizes from an object header (`v1`).
#[must_use]
pub fn read_hdf5_dataspace_dims(bytes: &[u8], obj_addr: usize) -> Option<Vec<u64>> {
    let mut result = None;
    iter_v1_header_messages(bytes, obj_addr, |mtype, body| {
        if mtype == 0x0001 && body.len() >= 4 {
            let ds_ver = body[0];
            let rank = usize::from(body[1]);
            let dims_off = match ds_ver {
                1 => 8usize,
                2 => 4usize,
                _ => return std::ops::ControlFlow::Break(()),
            };
            let mut dims = Vec::with_capacity(rank);
            for d in 0..rank {
                let Some(dim) = read_u64_le(body, dims_off + d * 8) else {
                    return std::ops::ControlFlow::Break(());
                };
                dims.push(dim);
            }
            result = Some(dims);
            return std::ops::ControlFlow::Break(());
        }
        std::ops::ControlFlow::Continue(())
    })?;
    result
}

/// Summary of an HDF5 `.hal` file read natively in pure Rust.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalFileInfo {
    /// Parsed HDF5 superblock metadata.
    pub superblock: Hdf5Superblock,
    /// Root-group member names (genomes and `Meta`/`Phylogeny` groups) for HDF5 v0 files.
    pub root_names: Vec<String>,
}

fn read_at<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    offset: u64,
    len: usize,
) -> std::io::Result<Vec<u8>> {
    r.seek(std::io::SeekFrom::Start(offset))?;
    let mut buf = vec![0u8; len];
    let mut filled = 0usize;
    while filled < len {
        let n = r.read(&mut buf[filled..])?;
        if n == 0 {
            break;
        }
        filled += n;
    }
    buf.truncate(filled);
    Ok(buf)
}

fn read_v1_obj_header_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    obj_addr: u64,
) -> std::io::Result<Option<Vec<u8>>> {
    let prefix = read_at(r, obj_addr, 16)?;
    if prefix.len() < 16 || prefix[0] != 1 {
        return Ok(None);
    }
    let hdr_size = usize::try_from(u32::from_le_bytes([
        prefix[8], prefix[9], prefix[10], prefix[11],
    ]))
    .unwrap_or(0);
    read_at(r, obj_addr, 16 + hdr_size).map(Some)
}

fn collect_snod_names_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    snod_addr: u64,
    heap_data: &[u8],
    out: &mut Vec<String>,
) -> std::io::Result<()> {
    let hdr = read_at(r, snod_addr, 8)?;
    if hdr.len() < 8 || &hdr[0..4] != b"SNOD" {
        return Ok(());
    }
    let n_syms = usize::from(u16::from_le_bytes([hdr[6], hdr[7]]));
    let entries_buf = read_at(r, snod_addr + 8, n_syms * 40)?;
    for i in 0..n_syms {
        let Some(name_off_u64) = read_u64_le(&entries_buf, i * 40) else {
            break;
        };
        let Ok(name_off) = usize::try_from(name_off_u64) else {
            continue;
        };
        if let Some(tail) = heap_data.get(name_off..) {
            let end = tail.iter().position(|&b| b == 0).unwrap_or(tail.len());
            if let Ok(s) = std::str::from_utf8(&tail[..end])
                && !s.is_empty()
            {
                out.push(s.to_owned());
            }
        }
    }
    Ok(())
}

fn walk_btree_v1_group_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    node_addr: u64,
    heap_data: &[u8],
    out: &mut Vec<String>,
) -> std::io::Result<()> {
    let hdr = read_at(r, node_addr, 24)?;
    if hdr.len() < 24 || &hdr[0..4] != b"TREE" || hdr[4] != 0 {
        return Ok(());
    }
    let node_level = hdr[5];
    let n_entries = usize::from(u16::from_le_bytes([hdr[6], hdr[7]]));
    let body = read_at(r, node_addr + 24, n_entries * 16)?;
    for i in 0..n_entries {
        let Some(child_addr) = read_u64_le(&body, i * 16 + 8) else {
            break;
        };
        if node_level == 0 {
            collect_snod_names_seek(r, child_addr, heap_data, out)?;
        } else {
            walk_btree_v1_group_seek(r, child_addr, heap_data, out)?;
        }
    }
    Ok(())
}

/// Extract root-group member names from an HDF5 v0 `.hal` stream via `Seek` without loading the full file.
///
/// # Errors
/// Returns an `std::io::Error` if seeking or reading fails.
pub fn read_hdf5_v0_root_names_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    sb: Hdf5Superblock,
) -> std::io::Result<Option<Vec<String>>> {
    if sb.version > 1 || sb.offset_size != 8 || sb.length_size != 8 {
        return Ok(None);
    }
    let Some(root_addr) = sb.root_group_addr else {
        return Ok(None);
    };
    let Some(root_slice) = read_v1_obj_header_seek(r, root_addr)? else {
        return Ok(None);
    };
    let num_msgs = usize::from(u16::from_le_bytes([root_slice[2], root_slice[3]]));
    let mut pos = 16usize;
    let end = root_slice.len();
    let mut symtab_addrs: Option<(u64, u64)> = None;
    let mut seen = 0usize;
    while seen < num_msgs && pos + 8 <= end {
        let mtype = u16::from_le_bytes([root_slice[pos], root_slice[pos + 1]]);
        let msize = usize::from(u16::from_le_bytes([
            root_slice[pos + 2],
            root_slice[pos + 3],
        ]));
        if mtype == 0x0011 && msize >= 16 && pos + 8 + 16 <= end {
            if let (Some(bt), Some(hp)) = (
                read_u64_le(&root_slice, pos + 8),
                read_u64_le(&root_slice, pos + 16),
            ) {
                symtab_addrs = Some((bt, hp));
            }
            break;
        }
        seen += 1;
        pos += 8 + msize;
    }
    let Some((btree_addr, heap_addr)) = symtab_addrs else {
        return Ok(None);
    };
    let heap_hdr = read_at(r, heap_addr, 32)?;
    if heap_hdr.len() < 32 || &heap_hdr[0..4] != b"HEAP" {
        return Ok(None);
    }
    let Some(data_seg_size) = read_u64_le(&heap_hdr, 8).and_then(|v| usize::try_from(v).ok())
    else {
        return Ok(None);
    };
    let Some(data_seg_addr) = read_u64_le(&heap_hdr, 24) else {
        return Ok(None);
    };
    let heap_data = read_at(r, data_seg_addr, data_seg_size)?;
    let mut names = Vec::new();
    walk_btree_v1_group_seek(r, btree_addr, &heap_data, &mut names)?;
    Ok(Some(names))
}

/// Inspect a `.hal` file on disk via `Seek` and extract its HDF5 superblock and v0 root member names.
///
/// # Errors
/// Returns an `std::io::Error` if the file cannot be opened or read.
pub fn inspect_hal_genomes(path: &std::path::Path) -> std::io::Result<Option<HalFileInfo>> {
    use std::io::Read as _;
    let mut f = std::fs::File::open(path)?;
    let mut buf = [0u8; 96];
    let n = f.read(&mut buf)?;
    let Some(superblock) = read_hdf5_superblock(&buf[..n]) else {
        return Ok(None);
    };
    let root_names = read_hdf5_v0_root_names_seek(&mut f, superblock)?.unwrap_or_default();
    Ok(Some(HalFileInfo {
        superblock,
        root_names,
    }))
}

/// Read contiguous dataset bytes from an HDF5 v1 object header (`DataLayout` `0x0008`, class `1`).
#[must_use]
pub fn read_hdf5_v0_contiguous_data(bytes: &[u8], obj_addr: usize) -> Option<&[u8]> {
    let mut result = None;
    iter_v1_header_messages(bytes, obj_addr, |mtype, body| {
        if mtype == 0x0008
            && body.len() >= 18
            && body[0] == 3
            && body[1] == 1
            && let Some(addr_u64) = read_u64_le(body, 2)
            && let Some(size_u64) = read_u64_le(body, 10)
            && let (Ok(data_addr), Ok(data_size)) =
                (usize::try_from(addr_u64), usize::try_from(size_u64))
        {
            result = bytes.get(data_addr..data_addr + data_size);
            return std::ops::ControlFlow::Break(());
        }
        std::ops::ControlFlow::Continue(())
    })?;
    result
}

/// One raw HDF5 chunk descriptor from a B-tree v1 type-1 chunk index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hdf5Chunk {
    /// Compressed/stored byte length of the chunk in the file.
    pub size: u32,
    /// File byte offset of the chunk payload.
    pub addr: u64,
}

fn walk_chunk_btree_v1(bytes: &[u8], node_addr: usize, ndim: usize, out: &mut Vec<Hdf5Chunk>) {
    if bytes.get(node_addr..node_addr + 4) != Some(b"TREE") {
        return;
    }
    let node_type = *bytes.get(node_addr + 4).unwrap_or(&255);
    let node_level = *bytes.get(node_addr + 5).unwrap_or(&255);
    let Some(n_entries) = bytes
        .get(node_addr + 6..node_addr + 8)
        .map(|b| usize::from(u16::from_le_bytes([b[0], b[1]])))
    else {
        return;
    };
    if node_type != 1 {
        return;
    }
    // Type-1 chunk key: 4 bytes chunk_size + 4 bytes filter_mask + (ndim + 1) * 8 bytes coords.
    let key_bytes = 8 + (ndim + 1) * 8;
    let rec_bytes = key_bytes + 8;
    for i in 0..n_entries {
        let key_off = node_addr + 24 + i * rec_bytes;
        let Some(sz_bytes) = bytes.get(key_off..key_off + 4) else {
            break;
        };
        let chunk_size = u32::from_le_bytes([sz_bytes[0], sz_bytes[1], sz_bytes[2], sz_bytes[3]]);
        let Some(child_addr) = read_u64_le(bytes, key_off + key_bytes) else {
            break;
        };
        if node_level == 0 {
            out.push(Hdf5Chunk {
                size: chunk_size,
                addr: child_addr,
            });
        } else if let Ok(sub) = usize::try_from(child_addr) {
            walk_chunk_btree_v1(bytes, sub, ndim, out);
        }
    }
}

/// Read chunked dataset descriptors (`DataLayout` `0x0008`, class `2`) from an HDF5 v1 object header.
#[must_use]
pub fn read_hdf5_v0_chunk_addrs(bytes: &[u8], obj_addr: usize) -> Option<Vec<Hdf5Chunk>> {
    let mut result = None;
    iter_v1_header_messages(bytes, obj_addr, |mtype, body| {
        if mtype == 0x0008
            && body.len() >= 11
            && body[0] == 3
            && body[1] == 2
            && let Some(bt_u64) = read_u64_le(body, 3)
            && let Ok(btree_addr) = usize::try_from(bt_u64)
        {
            let ndim = usize::from(body[2]);
            let mut chunks = Vec::new();
            walk_chunk_btree_v1(bytes, btree_addr, ndim, &mut chunks);
            result = Some(chunks);
            return std::ops::ControlFlow::Break(());
        }
        std::ops::ControlFlow::Continue(())
    })?;
    result
}

/// Check whether an HDF5 v1 object header enables the deflate (`id == 1`) filter in `FilterPipeline` (`0x000B`).
#[must_use]
pub fn is_hdf5_deflated(bytes: &[u8], obj_addr: usize) -> bool {
    let mut deflated = false;
    let _ = iter_v1_header_messages(bytes, obj_addr, |mtype, body| {
        if mtype == 0x000B && body.len() >= 4 {
            let ver = body[0];
            let n_filters = usize::from(body[1]);
            let mut fpos = match ver {
                1 => 8usize,
                2 => 2usize,
                _ => return std::ops::ControlFlow::Break(()),
            };
            for _ in 0..n_filters {
                if fpos + 2 > body.len() {
                    break;
                }
                let filter_id = u16::from_le_bytes([body[fpos], body[fpos + 1]]);
                if filter_id == 1 {
                    deflated = true;
                    return std::ops::ControlFlow::Break(());
                }
                if ver == 1 && fpos + 8 <= body.len() {
                    let name_len =
                        usize::from(u16::from_le_bytes([body[fpos + 2], body[fpos + 3]]));
                    let n_client =
                        usize::from(u16::from_le_bytes([body[fpos + 6], body[fpos + 7]]));
                    let padded_name = (name_len + 7) & !7;
                    let padded_client = ((n_client * 4) + 7) & !7;
                    fpos += 8 + padded_name + padded_client;
                } else if ver == 2 && fpos + 6 <= body.len() {
                    let flags = u16::from_le_bytes([body[fpos + 2], body[fpos + 3]]);
                    let n_client =
                        usize::from(u16::from_le_bytes([body[fpos + 4], body[fpos + 5]]));
                    let name_len = if filter_id >= 256 && fpos + 8 <= body.len() {
                        usize::from(u16::from_le_bytes([body[fpos + 6], body[fpos + 7]])) + 2
                    } else {
                        0
                    };
                    let _ = flags;
                    fpos += 6 + name_len + n_client * 4;
                } else {
                    break;
                }
            }
        }
        std::ops::ControlFlow::Continue(())
    });
    deflated
}

/// Read element byte size from an HDF5 `Datatype` (`0x0003`) message in an object header (`v1`).
#[must_use]
pub fn read_hdf5_datatype_size(bytes: &[u8], obj_addr: usize) -> Option<usize> {
    let mut size = None;
    iter_v1_header_messages(bytes, obj_addr, |mtype, body| {
        if mtype == 0x0003 && body.len() >= 8 {
            let sz = u32::from_le_bytes([body[4], body[5], body[6], body[7]]);
            size = usize::try_from(sz).ok();
            return std::ops::ControlFlow::Break(());
        }
        std::ops::ControlFlow::Continue(())
    })?;
    size
}

const MAX_CHUNK_DECOMPRESSED_BYTES: u64 = 256 * 1024 * 1024;

/// Read raw dataset bytes from an HDF5 v0 dataset object header, supporting both contiguous
/// (`DataLayout` class `1`) and chunked (`DataLayout` class `2`, optionally zlib-deflated) storage.
#[must_use]
pub fn read_hdf5_v0_dataset_bytes(bytes: &[u8], obj_addr: usize) -> Option<Vec<u8>> {
    if let Some(contig) = read_hdf5_v0_contiguous_data(bytes, obj_addr) {
        return Some(contig.to_vec());
    }
    let chunks = read_hdf5_v0_chunk_addrs(bytes, obj_addr)?;
    let deflated = is_hdf5_deflated(bytes, obj_addr);
    let mut out = Vec::new();
    for ch in chunks {
        let start = usize::try_from(ch.addr).ok()?;
        let len = usize::try_from(ch.size).ok()?;
        let raw_chunk = bytes.get(start..start + len)?;
        if deflated {
            use std::io::Read as _;
            let dec = flate2::read::ZlibDecoder::new(raw_chunk);
            let mut bounded = dec.take(MAX_CHUNK_DECOMPRESSED_BYTES + 1);
            let n = bounded.read_to_end(&mut out).ok()?;
            if u64::try_from(n).unwrap_or(u64::MAX) > MAX_CHUNK_DECOMPRESSED_BYTES {
                return None;
            }
        } else {
            out.extend_from_slice(raw_chunk);
        }
    }
    Some(out)
}

/// Decode HAL `SEQNAME_ARRAY` fixed-width null-terminated strings.
#[must_use]
pub fn decode_hal_seq_names(raw: &[u8], stride: usize) -> Vec<String> {
    if stride == 0 {
        return Vec::new();
    }
    raw.chunks_exact(stride)
        .map(|slot| {
            let end = slot.iter().position(|&b| b == 0).unwrap_or(slot.len());
            String::from_utf8_lossy(&slot[..end]).into_owned()
        })
        .collect()
}

/// Decoded HAL `BOTTOM_ARRAY` record: `(start_idx, length, top_parse_idx, child_links)`.
pub type HalBottomRecord = (i64, u64, i64, Vec<(i64, bool)>);

/// Decode HAL `BOTTOM_ARRAY` packed records (`24 + 9 * num_children` bytes per segment).
#[must_use]
pub fn decode_hal_bottom_array(raw: &[u8], num_children: usize) -> Vec<HalBottomRecord> {
    let rec_size = 24 + 9 * num_children;
    if rec_size == 0 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(raw.len() / rec_size);
    for rec in raw.chunks_exact(rec_size) {
        let start = i64::from_le_bytes(rec[0..8].try_into().unwrap_or([0; 8]));
        let length = u64::from_le_bytes(rec[8..16].try_into().unwrap_or([0; 8]));
        let top_parse = i64::from_le_bytes(rec[16..24].try_into().unwrap_or([0; 8]));
        let mut children = Vec::with_capacity(num_children);
        for c in 0..num_children {
            let off = 24 + c * 9;
            let child_idx = i64::from_le_bytes(rec[off..off + 8].try_into().unwrap_or([0; 8]));
            let rev = rec[off + 8] != 0;
            children.push((child_idx, rev));
        }
        out.push((start, length, top_parse, children));
    }
    out
}

/// Decode HAL `SEQIDX_ARRAY` 24-byte records `(start, top_offset, bottom_offset)`.
#[must_use]
pub fn decode_hal_seq_index(raw: &[u8]) -> Vec<(u64, u64, u64)> {
    let mut out = Vec::with_capacity(raw.len() / 24);
    for chunk in raw.as_chunks::<24>().0 {
        if let (Some(a), Some(b), Some(c)) = (
            read_u64_le(chunk, 0),
            read_u64_le(chunk, 8),
            read_u64_le(chunk, 16),
        ) {
            out.push((a, b, c));
        }
    }
    out
}

/// Decode HAL `TOP_ARRAY` 33-byte packed records into `(genome_idx, bottom_idx, paralogy_idx, parent_idx, reverse_flag)`.
#[must_use]
pub fn decode_hal_top_array(raw: &[u8]) -> Vec<(i64, i64, i64, i64, bool)> {
    let mut out = Vec::with_capacity(raw.len() / 33);
    for chunk in raw.as_chunks::<33>().0 {
        let g = i64::from_le_bytes(chunk[0..8].try_into().unwrap_or([0; 8]));
        let b = i64::from_le_bytes(chunk[8..16].try_into().unwrap_or([0; 8]));
        let p = i64::from_le_bytes(chunk[16..24].try_into().unwrap_or([0; 8]));
        let par = i64::from_le_bytes(chunk[24..32].try_into().unwrap_or([0; 8]));
        let rev = chunk[32] != 0;
        out.push((g, b, p, par, rev));
    }
    out
}

const fn decode_hal_nibble(nibble: u8) -> u8 {
    match nibble & 0x0f {
        0 => b'a',
        1 => b'c',
        2 => b'g',
        3 => b't',
        8 => b'A',
        9 => b'C',
        10 => b'G',
        11 => b'T',
        4 => b'n',
        _ => b'N',
    }
}

/// Decode HAL 4-bit nibble-packed `DNA_ARRAY` bytes (high nibble first) into ASCII DNA.
#[must_use]
pub fn decode_hal_dna(packed: &[u8], length: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(length);
    for &b in packed {
        if out.len() >= length {
            break;
        }
        out.push(decode_hal_nibble(b >> 4));
        if out.len() >= length {
            break;
        }
        out.push(decode_hal_nibble(b));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hal_star_roundtrip_and_hdf5_signature() {
        assert!(is_hdf5_hal(&[
            0x89, b'H', b'D', b'F', 0x0d, 0x0a, 0x1a, 0x0a, 0
        ]));
        assert!(!is_hdf5_hal(b"##maf version=1"));

        let ref_dna = b"ACGTACGTACGT";
        let q1_dna = b"ACGTACGTACGT";
        let q2_dna = b"ACGTACGTACGT";
        let hal = HalAlignment::from_star_blocks(
            "hg38",
            ref_dna,
            &[
                ("mm10".to_owned(), q1_dna.to_vec()),
                ("tupChi1".to_owned(), q2_dna.to_vec()),
            ],
            &[(0, 0, 0, 12, false), (1, 0, 0, 12, true)],
        );
        let maf = hal.to_maf_for_genome("hg38").expect("projected");
        assert!(maf.contains("s hg38 0 12 + 12 ACGTACGTACGT"));
        assert!(maf.contains("s mm10 0 12 + 12 ACGTACGTACGT"));
        assert!(maf.contains("s tupChi1 0 12 - 12 ACGTACGTACGT"));
    }

    #[test]
    fn hal_dna_nibble_decode_and_superblock_parse() {
        let packed = [0x8b_u8, 0xb8, 0x88, 0xaa, 0xbb, 0xb8];
        assert_eq!(decode_hal_dna(&packed, 11), b"ATTAAAGGTTT");

        let mut hdr = [0u8; 16];
        hdr[..8].copy_from_slice(&HDF5_SIGNATURE);
        hdr[8] = 2;
        hdr[9] = 8;
        hdr[10] = 8;
        assert_eq!(
            read_hdf5_superblock(&hdr),
            Some(Hdf5Superblock {
                version: 2,
                offset_size: 8,
                length_size: 8,
                root_group_addr: None,
            })
        );

        let mut ohdr = vec![0u8; 20];
        ohdr[0..4].copy_from_slice(b"OHDR");
        ohdr[4] = 2;
        ohdr[5] = 0; // 1-byte chunk size
        ohdr[6] = 6; // chunk size = 6
        ohdr[7] = 0x0A; // GroupInfo message
        ohdr[8..10].copy_from_slice(&2_u16.to_le_bytes());
        let parsed = read_hdf5_object_header(&ohdr).expect("valid OHDR");
        assert_eq!(parsed.version, 2);
        assert_eq!(parsed.message_types, vec![0x0A]);

        let mut seqidx = Vec::new();
        seqidx.extend_from_slice(&0_u64.to_le_bytes());
        seqidx.extend_from_slice(&0_u64.to_le_bytes());
        seqidx.extend_from_slice(&0_u64.to_le_bytes());
        seqidx.extend_from_slice(&29903_u64.to_le_bytes());
        seqidx.extend_from_slice(&1116_u64.to_le_bytes());
        seqidx.extend_from_slice(&0_u64.to_le_bytes());
        let parsed_idx = decode_hal_seq_index(&seqidx);
        assert_eq!(parsed_idx, vec![(0, 0, 0), (29903, 1116, 0)]);

        let mut top = Vec::new();
        top.extend_from_slice(&10_i64.to_le_bytes());
        top.extend_from_slice(&(-1_i64).to_le_bytes());
        top.extend_from_slice(&(-1_i64).to_le_bytes());
        top.extend_from_slice(&1_i64.to_le_bytes());
        top.push(1);
        assert_eq!(decode_hal_top_array(&top), vec![(10, -1, -1, 1, true)]);
    }

    #[test]
    fn hdf5_v0_group_entries_and_dataspace_parse() {
        // Construct a minimal synthetic HDF5 v0 buffer with:
        // - Root group object header (v1) at offset 0x60 pointing to TREE@0x100, HEAP@0x180
        // - HEAP at 0x180 with data segment at 0x1c0 containing "Gallus_gallus\0"
        // - TREE (type=0, level=0, n_entries=1) at 0x100 pointing to SNOD@0x200
        // - SNOD (n_syms=1) at 0x200 with name_off=0 and obj_hdr_addr=0x280
        // - Dataset object header (v1) at 0x280 with Dataspace message (0x0001, v1, rank=1, dim=42)
        let mut buf = vec![0u8; 0x300];
        buf[..8].copy_from_slice(&HDF5_SIGNATURE);
        buf[8] = 0; // v0
        buf[13] = 8; // offset_size
        buf[14] = 8; // length_size
        buf[64..72].copy_from_slice(&0x60_u64.to_le_bytes());

        // Root group v1 object header at 0x60
        buf[0x60] = 1; // version 1
        buf[0x62..0x64].copy_from_slice(&1_u16.to_le_bytes()); // num_msgs = 1
        buf[0x68..0x6c].copy_from_slice(&24_u32.to_le_bytes()); // hdr_size = 24
        buf[0x70..0x72].copy_from_slice(&0x0011_u16.to_le_bytes()); // SymbolTable
        buf[0x72..0x74].copy_from_slice(&16_u16.to_le_bytes());
        buf[0x78..0x80].copy_from_slice(&0x100_u64.to_le_bytes()); // btree_addr
        buf[0x80..0x88].copy_from_slice(&0x180_u64.to_le_bytes()); // heap_addr

        // TREE at 0x100
        buf[0x100..0x104].copy_from_slice(b"TREE");
        buf[0x104] = 0; // group
        buf[0x105] = 0; // level 0
        buf[0x106..0x108].copy_from_slice(&1_u16.to_le_bytes());
        buf[0x100 + 24 + 8..0x100 + 24 + 16].copy_from_slice(&0x200_u64.to_le_bytes());

        // HEAP at 0x180
        buf[0x180..0x184].copy_from_slice(b"HEAP");
        buf[0x188..0x190].copy_from_slice(&32_u64.to_le_bytes()); // data_seg_size
        buf[0x198..0x1a0].copy_from_slice(&0x1c0_u64.to_le_bytes()); // data_seg_addr
        buf[0x1c0..0x1c0 + 14].copy_from_slice(b"Gallus_gallus\0");

        // SNOD at 0x200
        buf[0x200..0x204].copy_from_slice(b"SNOD");
        buf[0x206..0x208].copy_from_slice(&1_u16.to_le_bytes());
        buf[0x208..0x210].copy_from_slice(&0_u64.to_le_bytes()); // name_off
        buf[0x210..0x218].copy_from_slice(&0x280_u64.to_le_bytes()); // obj_hdr_addr

        // Dataset v1 object header at 0x280 with Dataspace (0x0001)
        buf[0x280] = 1;
        buf[0x282..0x284].copy_from_slice(&1_u16.to_le_bytes());
        buf[0x288..0x28c].copy_from_slice(&24_u32.to_le_bytes());
        buf[0x290..0x292].copy_from_slice(&0x0001_u16.to_le_bytes()); // Dataspace
        buf[0x292..0x294].copy_from_slice(&16_u16.to_le_bytes());
        buf[0x298] = 1; // ds_ver = 1
        buf[0x299] = 1; // rank = 1
        buf[0x2a0..0x2a8].copy_from_slice(&42_u64.to_le_bytes()); // dim[0] = 42

        let entries = read_hdf5_v0_group_entries(&buf, 0x60).expect("group entries");
        assert_eq!(
            entries,
            vec![Hdf5GroupEntry {
                name: "Gallus_gallus".to_owned(),
                object_header_addr: 0x280,
            }]
        );
        assert_eq!(
            read_hdf5_v0_root_names(&buf),
            Some(vec!["Gallus_gallus".to_owned()])
        );
        let sb = read_hdf5_superblock(&buf).expect("superblock");
        let mut cursor = std::io::Cursor::new(&buf);
        assert_eq!(
            read_hdf5_v0_root_names_seek(&mut cursor, sb).expect("seek"),
            Some(vec!["Gallus_gallus".to_owned()])
        );
        assert_eq!(read_hdf5_dataspace_dims(&buf, 0x280), Some(vec![42]));
    }

    #[test]
    fn hdf5_deflated_chunk_and_hal_bottom_seqname_decode() {
        use std::io::Write as _;
        let payload = b"ACGTACGTACGTACGT";
        let mut enc = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(payload).expect("zlib encode");
        let compressed = enc.finish().expect("zlib finish");

        // Synthetic HDF5 v1 object header at 0x20 with:
        // - Datatype (0x0003) size = 16
        // - FilterPipeline (0x000B) v2 with 1 filter (id=1 deflate)
        // - DataLayout (0x0008) v3 class=2 (chunked, ndim=1) pointing to TREE@0x100
        let mut buf = vec![0u8; 0x300];
        buf[0x20] = 1; // v1 obj header
        buf[0x22..0x24].copy_from_slice(&3_u16.to_le_bytes());
        buf[0x28..0x2c].copy_from_slice(&64_u32.to_le_bytes());

        // Msg 1: Datatype (0x0003), size=8
        let mut p = 0x30;
        buf[p..p + 2].copy_from_slice(&0x0003_u16.to_le_bytes());
        buf[p + 2..p + 4].copy_from_slice(&8_u16.to_le_bytes());
        buf[p + 12..p + 16].copy_from_slice(&16_u32.to_le_bytes());
        p += 16;

        // Msg 2: FilterPipeline (0x000B) v2, 1 filter (id=1 deflate), size=8
        buf[p..p + 2].copy_from_slice(&0x000B_u16.to_le_bytes());
        buf[p + 2..p + 4].copy_from_slice(&8_u16.to_le_bytes());
        buf[p + 8] = 2; // version 2
        buf[p + 9] = 1; // n_filters = 1
        buf[p + 10..p + 12].copy_from_slice(&1_u16.to_le_bytes()); // deflate filter_id = 1
        p += 16;

        // Msg 3: DataLayout (0x0008) v3 class=2 ndim=1, size=12
        buf[p..p + 2].copy_from_slice(&0x0008_u16.to_le_bytes());
        buf[p + 2..p + 4].copy_from_slice(&12_u16.to_le_bytes());
        buf[p + 8] = 3; // layout_ver
        buf[p + 9] = 2; // chunked
        buf[p + 10] = 1; // ndim = 1
        buf[p + 11..p + 19].copy_from_slice(&0x100_u64.to_le_bytes());

        // Chunk B-tree v1 type-1 at 0x100 (ndim=1 => key_bytes=24, child_addr at +24)
        buf[0x100..0x104].copy_from_slice(b"TREE");
        buf[0x104] = 1; // chunk btree
        buf[0x105] = 0; // leaf
        buf[0x106..0x108].copy_from_slice(&1_u16.to_le_bytes());
        let ch_sz = u32::try_from(compressed.len()).expect("len");
        buf[0x118..0x11c].copy_from_slice(&ch_sz.to_le_bytes());
        buf[0x130..0x138].copy_from_slice(&0x200_u64.to_le_bytes());
        buf[0x200..0x200 + compressed.len()].copy_from_slice(&compressed);

        assert!(is_hdf5_deflated(&buf, 0x20));
        assert_eq!(read_hdf5_datatype_size(&buf, 0x20), Some(16));
        assert_eq!(
            read_hdf5_v0_dataset_bytes(&buf, 0x20).as_deref(),
            Some(payload.as_slice())
        );

        let names_raw = b"chr1\0\0\0\0chr2\0\0\0\0";
        assert_eq!(decode_hal_seq_names(names_raw, 8), vec!["chr1", "chr2"]);

        // 1 bottom segment with 2 children: 24 + 9*2 = 42 bytes
        let mut bot = Vec::new();
        bot.extend_from_slice(&100_i64.to_le_bytes());
        bot.extend_from_slice(&50_u64.to_le_bytes());
        bot.extend_from_slice(&(-1_i64).to_le_bytes());
        bot.extend_from_slice(&5_i64.to_le_bytes());
        bot.push(0);
        bot.extend_from_slice(&(-1_i64).to_le_bytes());
        bot.push(1);
        assert_eq!(
            decode_hal_bottom_array(&bot, 2),
            vec![(100, 50, -1, vec![(5, false), (-1, true)])]
        );
    }
}
