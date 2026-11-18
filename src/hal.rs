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
        // HDF5 Object Header v1
        let num_msgs = usize::from(u16::from_le_bytes([slice[2], slice[3]]));
        let hdr_size = usize::try_from(u32::from_le_bytes([
            slice[8], slice[9], slice[10], slice[11],
        ]))
        .ok()?;
        let mut pos = 16usize;
        let end = (16 + hdr_size).min(slice.len());
        let mut msg_types = Vec::new();
        while msg_types.len() < num_msgs && pos + 8 <= end {
            let mtype = u16::from_le_bytes([slice[pos], slice[pos + 1]]);
            let msize = usize::from(u16::from_le_bytes([slice[pos + 2], slice[pos + 3]]));
            msg_types.push(mtype);
            pos += 8 + msize;
        }
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

fn collect_snod_names(bytes: &[u8], snod_addr: usize, heap_data: &[u8], out: &mut Vec<String>) {
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
}

fn walk_btree_v1_group(bytes: &[u8], node_addr: usize, heap_data: &[u8], out: &mut Vec<String>) {
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
            collect_snod_names(bytes, child_addr, heap_data, out);
        } else {
            walk_btree_v1_group(bytes, child_addr, heap_data, out);
        }
    }
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
    let root_slice = bytes.get(root_addr..)?;
    if root_slice.len() < 16 || root_slice[0] != 1 {
        return None;
    }
    let num_msgs = usize::from(u16::from_le_bytes([root_slice[2], root_slice[3]]));
    let hdr_size = usize::try_from(u32::from_le_bytes([
        root_slice[8],
        root_slice[9],
        root_slice[10],
        root_slice[11],
    ]))
    .ok()?;
    let mut pos = 16usize;
    let end = (16 + hdr_size).min(root_slice.len());
    let mut symtab_addrs: Option<(usize, usize)> = None;
    let mut seen = 0usize;
    while seen < num_msgs && pos + 8 <= end {
        let mtype = u16::from_le_bytes([root_slice[pos], root_slice[pos + 1]]);
        let msize = usize::from(u16::from_le_bytes([
            root_slice[pos + 2],
            root_slice[pos + 3],
        ]));
        if mtype == 0x0011 && msize >= 16 && pos + 8 + 16 <= end {
            let btree_addr = usize::try_from(read_u64_le(root_slice, pos + 8)?).ok()?;
            let heap_addr = usize::try_from(read_u64_le(root_slice, pos + 16)?).ok()?;
            symtab_addrs = Some((btree_addr, heap_addr));
            break;
        }
        seen += 1;
        pos += 8 + msize;
    }
    let (btree_addr, heap_addr) = symtab_addrs?;
    if bytes.get(heap_addr..heap_addr + 4) != Some(b"HEAP") {
        return None;
    }
    let data_seg_size = usize::try_from(read_u64_le(bytes, heap_addr + 8)?).ok()?;
    let data_seg_addr = usize::try_from(read_u64_le(bytes, heap_addr + 24)?).ok()?;
    let heap_data = bytes.get(data_seg_addr..data_seg_addr + data_seg_size)?;
    let mut names = Vec::new();
    walk_btree_v1_group(bytes, btree_addr, heap_data, &mut names);
    Some(names)
}

/// Summary of an HDF5 `.hal` file read natively in pure Rust.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalFileInfo {
    /// Parsed HDF5 superblock metadata.
    pub superblock: Hdf5Superblock,
    /// Root-group member names (genomes and `Meta`/`Phylogeny` groups) for HDF5 v0 files.
    pub root_names: Vec<String>,
}

/// Inspect a `.hal` file on disk and extract its HDF5 superblock and v0 root member names.
///
/// # Errors
/// Returns an `std::io::Error` if the file cannot be opened or read.
pub fn inspect_hal_genomes(path: &std::path::Path) -> std::io::Result<Option<HalFileInfo>> {
    let bytes = std::fs::read(path)?;
    let Some(superblock) = read_hdf5_superblock(&bytes) else {
        return Ok(None);
    };
    let root_names = read_hdf5_v0_root_names(&bytes).unwrap_or_default();
    Ok(Some(HalFileInfo {
        superblock,
        root_names,
    }))
}

/// Read contiguous dataset bytes from an HDF5 v1 object header (`DataLayout` `0x0008`, class `1`).
#[must_use]
pub fn read_hdf5_v0_contiguous_data(bytes: &[u8], obj_addr: usize) -> Option<&[u8]> {
    let slice = bytes.get(obj_addr..)?;
    if slice.len() < 16 || slice[0] != 1 {
        return None;
    }
    let num_msgs = usize::from(u16::from_le_bytes([slice[2], slice[3]]));
    let hdr_size = usize::try_from(u32::from_le_bytes([
        slice[8], slice[9], slice[10], slice[11],
    ]))
    .ok()?;
    let mut pos = 16usize;
    let end = (16 + hdr_size).min(slice.len());
    let mut seen = 0usize;
    while seen < num_msgs && pos + 8 <= end {
        let mtype = u16::from_le_bytes([slice[pos], slice[pos + 1]]);
        let msize = usize::from(u16::from_le_bytes([slice[pos + 2], slice[pos + 3]]));
        if mtype == 0x0008 && msize >= 18 && pos + 8 + 18 <= end {
            let layout_ver = slice[pos + 8];
            let layout_class = slice[pos + 9];
            if layout_ver == 3 && layout_class == 1 {
                let data_addr = usize::try_from(read_u64_le(slice, pos + 10)?).ok()?;
                let data_size = usize::try_from(read_u64_le(slice, pos + 18)?).ok()?;
                return bytes.get(data_addr..data_addr + data_size);
            }
        }
        seen += 1;
        pos += 8 + msize;
    }
    None
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
    let slice = bytes.get(obj_addr..)?;
    if slice.len() < 16 || slice[0] != 1 {
        return None;
    }
    let num_msgs = usize::from(u16::from_le_bytes([slice[2], slice[3]]));
    let hdr_size = usize::try_from(u32::from_le_bytes([
        slice[8], slice[9], slice[10], slice[11],
    ]))
    .ok()?;
    let mut pos = 16usize;
    let end = (16 + hdr_size).min(slice.len());
    let mut seen = 0usize;
    while seen < num_msgs && pos + 8 <= end {
        let mtype = u16::from_le_bytes([slice[pos], slice[pos + 1]]);
        let msize = usize::from(u16::from_le_bytes([slice[pos + 2], slice[pos + 3]]));
        if mtype == 0x0008 && msize >= 11 && pos + 8 + 11 <= end {
            let layout_ver = slice[pos + 8];
            let layout_class = slice[pos + 9];
            if layout_ver == 3 && layout_class == 2 {
                let ndim = usize::from(slice[pos + 10]);
                let btree_addr = usize::try_from(read_u64_le(slice, pos + 11)?).ok()?;
                let mut chunks = Vec::new();
                walk_chunk_btree_v1(bytes, btree_addr, ndim, &mut chunks);
                return Some(chunks);
            }
        }
        seen += 1;
        pos += 8 + msize;
    }
    None
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
}
