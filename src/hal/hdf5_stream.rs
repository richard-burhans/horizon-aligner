//! Streaming chunk-slice HDF5 v0/v1 dataset reader and GCOL/attribute string decoder.
//!
//! Removes the `256 MiB` full-dataset buffer cap for multi-gigabyte HAL alignments (e.g.
//! `4-t2t-orangs-mc-2023v2.full.hal` and `jarvis_lab_birds.hal`) by reading and caching only
//! the specific contiguous byte range or compressed chunks overlapping a queried element slice.

use super::hdf5_dataset::{FilterSummary, Hdf5Chunk, parse_dataspace_dims, parse_filter_pipeline};
use super::hdf5_read::{iter_v1_header_messages_seek, read_at, read_at_into, read_u64_le};
use std::collections::{HashMap, VecDeque};
use std::io::{Read, Seek};

/// Safety limit on a single decompressed HDF5 chunk (64 MiB), while allowing arbitrary multi-GiB
/// logical datasets to be queried via chunk-slice reads.
const MAX_SINGLE_CHUNK_BYTES: u64 = 64 * 1024 * 1024;
const MAX_GCOL_BYTES: usize = 16 * 1024 * 1024;
/// Default maximum number of decompressed chunks cached per HDF5 dataset reader.
pub const DEFAULT_CHUNK_CACHE_CAP: usize = 64;
/// Soft byte budget per dataset chunk cache (16 MiB, always retaining at least 2 chunks).
const MAX_DATASET_CACHE_BYTES: usize = 16 * 1024 * 1024;
/// Maximum number of raw HDF5 chunk B-tree nodes cached per dataset reader.
const MAX_BTREE_NODE_CACHE: usize = 512;

/// Read a variable-length string object from an HDF5 Global Heap collection (`GCOL`).
fn read_gcol_object<R: Read + Seek>(
    r: &mut R,
    gcol_addr: u64,
    obj_idx: u32,
) -> std::io::Result<Option<Vec<u8>>> {
    if gcol_addr == u64::MAX || obj_idx == 0 {
        return Ok(None);
    }
    let hdr = read_at(r, gcol_addr, 16)?;
    if hdr.len() < 16 || &hdr[0..4] != b"GCOL" || hdr[4] != 1 {
        return Ok(None);
    }
    let Some(col_sz_u64) = read_u64_le(&hdr, 8) else {
        return Ok(None);
    };
    let Ok(col_sz) = usize::try_from(col_sz_u64) else {
        return Ok(None);
    };
    if !(16..=MAX_GCOL_BYTES).contains(&col_sz) {
        return Ok(None);
    }
    let buf = read_at(r, gcol_addr, col_sz)?;
    if buf.len() < col_sz {
        return Ok(None);
    }
    let mut pos = 16usize;
    while pos + 16 <= buf.len() {
        let oid = u16::from_le_bytes([buf[pos], buf[pos + 1]]);
        if oid == 0 {
            break;
        }
        let Some(osz_u64) = read_u64_le(&buf, pos + 8) else {
            break;
        };
        let Ok(osz) = usize::try_from(osz_u64) else {
            break;
        };
        let data_start = pos + 16;
        let Some(data_end) = data_start.checked_add(osz) else {
            break;
        };
        if data_end > buf.len() {
            break;
        }
        if u32::from(oid) == obj_idx {
            return Ok(Some(buf[data_start..data_end].to_vec()));
        }
        let padded = (osz + 7) & !7;
        pos = data_start.saturating_add(padded);
    }
    Ok(None)
}

/// Read a scalar string attribute (fixed-length `class=3` or VLS `class=9` via `GCOL`) from an
/// HDF5 v1 object header.
///
/// # Errors
/// Returns an `std::io::Error` if seeking or reading fails.
pub fn read_hdf5_v0_attr_string_seek<R: Read + Seek>(
    r: &mut R,
    obj_addr: u64,
    attr_name: &str,
) -> std::io::Result<Option<String>> {
    enum AttrValue {
        Fixed(Vec<u8>),
        Vls {
            len: usize,
            gcol_addr: u64,
            obj_idx: u32,
        },
    }

    let mut found: Option<AttrValue> = None;
    let Some(()) = iter_v1_header_messages_seek(r, obj_addr, |mtype, body| {
        if mtype != 0x000C || body.len() < 8 {
            return std::ops::ControlFlow::Continue(());
        }
        let ver = body[0];
        let name_sz = usize::from(u16::from_le_bytes([body[2], body[3]]));
        let type_sz = usize::from(u16::from_le_bytes([body[4], body[5]]));
        let space_sz = usize::from(u16::from_le_bytes([body[6], body[7]]));
        let hdr_len = match ver {
            1 | 2 => 8usize,
            3 => 9usize,
            _ => return std::ops::ControlFlow::Continue(()),
        };
        let pad = |n: usize| if ver == 1 { (n + 7) & !7 } else { n };
        let name_off = hdr_len;
        let type_off = name_off.saturating_add(pad(name_sz));
        let space_off = type_off.saturating_add(pad(type_sz));
        let val_off = space_off.saturating_add(pad(space_sz));
        if val_off > body.len() || name_off + name_sz > body.len() || type_off + 8 > body.len() {
            return std::ops::ControlFlow::Continue(());
        }
        let raw_name = &body[name_off..name_off + name_sz];
        let trimmed = raw_name
            .iter()
            .position(|&b| b == 0)
            .map_or(raw_name, |z| &raw_name[..z]);
        if trimmed != attr_name.as_bytes() {
            return std::ops::ControlFlow::Continue(());
        }
        let dt_class = body[type_off] & 0x0F;
        let val = &body[val_off..];
        if dt_class == 3 {
            let sz = usize::try_from(u32::from_le_bytes([
                body[type_off + 4],
                body[type_off + 5],
                body[type_off + 6],
                body[type_off + 7],
            ]))
            .unwrap_or(val.len())
            .min(val.len());
            found = Some(AttrValue::Fixed(val[..sz].to_vec()));
            return std::ops::ControlFlow::Break(());
        }
        if dt_class == 9 && val.len() >= 16 {
            let vlen =
                usize::try_from(u32::from_le_bytes([val[0], val[1], val[2], val[3]])).unwrap_or(0);
            let gcol_addr = read_u64_le(val, 4).unwrap_or(u64::MAX);
            let obj_idx = u32::from_le_bytes([val[12], val[13], val[14], val[15]]);
            found = Some(AttrValue::Vls {
                len: vlen,
                gcol_addr,
                obj_idx,
            });
            return std::ops::ControlFlow::Break(());
        }
        std::ops::ControlFlow::Continue(())
    })?
    else {
        return Ok(None);
    };

    match found {
        None => Ok(None),
        Some(AttrValue::Fixed(bytes)) => {
            let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
            Ok(Some(String::from_utf8_lossy(&bytes[..end]).into_owned()))
        }
        Some(AttrValue::Vls {
            len,
            gcol_addr,
            obj_idx,
        }) => {
            let Some(raw) = read_gcol_object(r, gcol_addr, obj_idx)? else {
                return Ok(None);
            };
            let take = if len > 0 {
                len.min(raw.len())
            } else {
                raw.len()
            };
            let slice = &raw[..take];
            let end = slice.iter().position(|&b| b == 0).unwrap_or(slice.len());
            Ok(Some(String::from_utf8_lossy(&slice[..end]).into_owned()))
        }
    }
}

enum DatasetStorage {
    Empty,
    Contiguous {
        addr: u64,
        byte_len: u64,
    },
    Chunked {
        btree_addr: u64,
        ndim: usize,
        chunk_elems: u64,
        deflate_layers: usize,
        /// Tunable LRU capacity (maximum cached decompressed chunks).
        cache_cap: usize,
        /// Current total decompressed bytes held in `cache`.
        cache_bytes: usize,
        /// True LRU cache keyed by chunk's `offset_elem` (bounded RAM).
        cache: VecDeque<(u64, Vec<u8>)>,
        /// Cached raw B-tree nodes `(24-byte header + body)` keyed by file address.
        btree_nodes: HashMap<u64, Vec<u8>>,
        /// Reusable temporary buffers for chunk reading & zlib decompression.
        raw_tmp: Vec<u8>,
        dec_tmp: Vec<u8>,
    },
}

const MAX_BTREE_NODES_PER_QUERY: usize = 4_096;
const MAX_CHUNKS_PER_QUERY: usize = 4_096;
const MAX_DECOMPRESSED_PER_CALL: u64 = 256 * 1024 * 1024;

struct RangeQueryBudget {
    visited: std::collections::HashSet<u64>,
    nodes_seen: usize,
}

impl RangeQueryBudget {
    fn new() -> Self {
        Self {
            visited: std::collections::HashSet::new(),
            nodes_seen: 0,
        }
    }
}

/// Descend an HDF5 v1 chunk B-tree (`TREE` type 1) in `O(log C)` to collect ONLY the chunks
/// overlapping `[query_start, query_end)`, failing closed on cycles, non-decreasing levels,
/// truncation, or budget exhaustion.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn query_chunk_btree_range<R: Read + Seek>(
    r: &mut R,
    node_addr: u64,
    ndim: usize,
    chunk_elems: u64,
    query_start: u64,
    query_end: u64,
    parent_level: Option<u8>,
    depth: usize,
    budget: &mut RangeQueryBudget,
    btree_nodes: &mut HashMap<u64, Vec<u8>>,
    out: &mut Vec<Hdf5Chunk>,
) -> std::io::Result<()> {
    if depth > 64 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "HDF5 chunk B-tree exceeds max depth 64",
        ));
    }
    if node_addr == u64::MAX {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "HDF5 chunk B-tree node has HADDR_UNDEF",
        ));
    }
    if !budget.visited.insert(node_addr) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "HDF5 chunk B-tree contains a cycle or shared node DAG",
        ));
    }
    budget.nodes_seen += 1;
    if budget.nodes_seen > MAX_BTREE_NODES_PER_QUERY {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "HDF5 chunk B-tree query exceeded node visit safety cap",
        ));
    }
    let key_bytes = 8usize.saturating_add(ndim.saturating_mul(8));
    let rec_bytes = key_bytes.saturating_add(8);

    // Load full node `(24-byte header + body)` from cache or disk once.
    let node_buf: Vec<u8> = if let Some(cached) = btree_nodes.get(&node_addr) {
        cached.clone()
    } else {
        let hdr = read_at(r, node_addr, 24)?;
        if hdr.len() < 24 || &hdr[0..4] != b"TREE" || hdr[4] != 1 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "invalid HDF5 chunk B-tree node signature or type",
            ));
        }
        let node_level = hdr[5];
        let n_entries = usize::from(u16::from_le_bytes([hdr[6], hdr[7]]));
        let extra = if node_level > 0 { key_bytes } else { 0 };
        let read_len = n_entries
            .checked_mul(rec_bytes)
            .and_then(|n| n.checked_add(extra))
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "B-tree record size overflow",
                )
            })?;
        let mut full = Vec::with_capacity(24 + read_len);
        full.extend_from_slice(&hdr);
        if read_len > 0 {
            let body_off = node_addr.checked_add(24).ok_or_else(|| {
                std::io::Error::new(std::io::ErrorKind::InvalidData, "B-tree offset overflow")
            })?;
            let body = read_at(r, body_off, read_len)?;
            full.extend_from_slice(&body);
        }
        if btree_nodes.len() < MAX_BTREE_NODE_CACHE {
            btree_nodes.insert(node_addr, full.clone());
        }
        full
    };

    if node_buf.len() < 24 || &node_buf[0..4] != b"TREE" || node_buf[4] != 1 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "invalid HDF5 chunk B-tree node signature or type",
        ));
    }
    let node_level = node_buf[5];
    if let Some(p_lvl) = parent_level
        && node_level >= p_lvl
    {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "HDF5 chunk B-tree level failed to strictly decrease",
        ));
    }
    let n_entries = usize::from(u16::from_le_bytes([node_buf[6], node_buf[7]]));
    if n_entries == 0 {
        return Ok(());
    }
    let body = &node_buf[24..];
    if body.len() < n_entries * rec_bytes {
        return Err(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            "truncated HDF5 chunk B-tree node body",
        ));
    }
    if node_level == 0 {
        for i in 0..n_entries {
            let key_off = i * rec_bytes;
            let offset_elem = read_u64_le(body, key_off + 8).ok_or_else(|| {
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "truncated chunk key")
            })?;
            if offset_elem >= query_end {
                break;
            }
            if offset_elem.saturating_add(chunk_elems) <= query_start {
                continue;
            }
            if out.len() >= MAX_CHUNKS_PER_QUERY {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "HDF5 chunk range query exceeded max chunk cap",
                ));
            }
            let chunk_size = u32::from_le_bytes([
                body[key_off],
                body[key_off + 1],
                body[key_off + 2],
                body[key_off + 3],
            ]);
            let filter_mask = u32::from_le_bytes([
                body[key_off + 4],
                body[key_off + 5],
                body[key_off + 6],
                body[key_off + 7],
            ]);
            let child_addr = read_u64_le(body, key_off + key_bytes).ok_or_else(|| {
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "truncated chunk addr")
            })?;
            out.push(Hdf5Chunk {
                size: chunk_size,
                filter_mask,
                addr: child_addr,
                offset_elem,
            });
        }
    } else {
        let mut descend = Vec::new();
        for i in 0..n_entries {
            let key_off = i * rec_bytes;
            let k_lo = read_u64_le(body, key_off + 8).ok_or_else(|| {
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "truncated internal key")
            })?;
            if k_lo >= query_end {
                break;
            }
            let next_key_off = (i + 1) * rec_bytes;
            let k_hi = if next_key_off + 16 <= body.len() {
                read_u64_le(body, next_key_off + 8)
            } else {
                None
            };
            if k_hi.is_some_and(|hi| hi > k_lo && hi.saturating_add(chunk_elems) <= query_start) {
                continue;
            }
            let child_addr = read_u64_le(body, key_off + key_bytes).ok_or_else(|| {
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "truncated child addr")
            })?;
            descend.push(child_addr);
        }
        drop(node_buf);
        for child_addr in descend {
            query_chunk_btree_range(
                r,
                child_addr,
                ndim,
                chunk_elems,
                query_start,
                query_end,
                Some(node_level),
                depth + 1,
                budget,
                btree_nodes,
                out,
            )?;
        }
    }
    Ok(())
}

/// Random-access streaming reader for 1D HDF5 v0/v1 contiguous or chunked (optionally zlib-deflated)
/// datasets of arbitrary multi-GiB size.
pub struct Hdf5DatasetReader {
    total_elems: u64,
    elem_size: usize,
    storage: DatasetStorage,
}

impl Hdf5DatasetReader {
    /// Open a 1D HDF5 dataset at object header `obj_addr` without loading its payload into memory,
    /// using [`DEFAULT_CHUNK_CACHE_CAP`].
    ///
    /// # Errors
    /// Returns an `std::io::Error` if reading metadata or chunk B-tree nodes fails.
    pub fn open<R: Read + Seek>(r: &mut R, obj_addr: u64) -> std::io::Result<Option<Self>> {
        Self::open_with_cache_cap(r, obj_addr, DEFAULT_CHUNK_CACHE_CAP)
    }

    /// Open a 1D HDF5 dataset with an explicit per-dataset LRU chunk cache capacity (`>= 1`).
    ///
    /// # Errors
    /// Returns an `std::io::Error` if reading metadata or chunk B-tree nodes fails.
    pub fn open_with_cache_cap<R: Read + Seek>(
        r: &mut R,
        obj_addr: u64,
        cache_cap: usize,
    ) -> std::io::Result<Option<Self>> {
        let cap = cache_cap.max(1);
        let mut dims: Option<Vec<u64>> = None;
        let mut elem_size: Option<usize> = None;
        let mut contig: Option<(u64, u64)> = None;
        let mut chunked: Option<(u64, usize, u64)> = None;
        let mut filters = FilterSummary::default();

        let Some(()) = iter_v1_header_messages_seek(r, obj_addr, |mtype, body| {
            if mtype == 0x0001 {
                dims = parse_dataspace_dims(body);
            } else if mtype == 0x0003 && body.len() >= 8 {
                let sz = u32::from_le_bytes([body[4], body[5], body[6], body[7]]);
                elem_size = (sz > 0).then(|| usize::try_from(sz).ok()).flatten();
            } else if mtype == 0x0008 && body.len() >= 18 && body[0] == 3 && body[1] == 1 {
                if let (Some(addr), Some(sz_u64)) = (read_u64_le(body, 2), read_u64_le(body, 10)) {
                    contig = Some((addr, sz_u64));
                }
            } else if mtype == 0x0008 && body.len() >= 15 && body[0] == 3 && body[1] == 2 {
                if let Some(bt_addr) = read_u64_le(body, 3) {
                    let ndim = usize::from(body[2]);
                    let chunk_dim0 =
                        u64::from(u32::from_le_bytes([body[11], body[12], body[13], body[14]]));
                    chunked = Some((bt_addr, ndim, chunk_dim0));
                }
            } else if mtype == 0x000B {
                let fs = parse_filter_pipeline(body);
                filters.deflate_layers = filters.deflate_layers.saturating_add(fs.deflate_layers);
                filters.has_unsupported |= fs.has_unsupported;
            }
            std::ops::ControlFlow::Continue(())
        })?
        else {
            return Ok(None);
        };

        if filters.has_unsupported || filters.deflate_layers > 2 {
            return Ok(None);
        }
        let Some(es) = elem_size else {
            return Ok(None);
        };
        let Some(dims_vec) = dims else {
            return Ok(None);
        };
        let Some(total_elems) = dims_vec.iter().try_fold(1u64, |acc, &x| acc.checked_mul(x)) else {
            return Ok(None);
        };

        if total_elems == 0 {
            return Ok(Some(Self {
                total_elems: 0,
                elem_size: es,
                storage: DatasetStorage::Empty,
            }));
        }
        if let Some((addr, byte_len)) = contig {
            if byte_len == 0 || addr == u64::MAX {
                return Ok(Some(Self {
                    total_elems: 0,
                    elem_size: es,
                    storage: DatasetStorage::Empty,
                }));
            }
            return Ok(Some(Self {
                total_elems,
                elem_size: es,
                storage: DatasetStorage::Contiguous { addr, byte_len },
            }));
        }
        let Some((btree_addr, ndim, chunk_dim0)) = chunked else {
            return Ok(None);
        };
        if btree_addr == u64::MAX || ndim == 0 || ndim > 2 {
            return Ok(None);
        }
        let chunk_elems = if chunk_dim0 > 0 {
            chunk_dim0
        } else {
            total_elems.max(1)
        };
        Ok(Some(Self {
            total_elems,
            elem_size: es,
            storage: DatasetStorage::Chunked {
                btree_addr,
                ndim,
                chunk_elems,
                deflate_layers: filters.deflate_layers,
                cache_cap: cap,
                cache_bytes: 0,
                cache: VecDeque::with_capacity(cap.min(64)),
                btree_nodes: HashMap::new(),
                raw_tmp: Vec::new(),
                dec_tmp: Vec::new(),
            },
        }))
    }

    /// Tune the per-dataset LRU chunk cache capacity (`>= 1`), immediately evicting oldest entries
    /// if the new cap is smaller than the current cached count.
    pub fn set_cache_cap(&mut self, cap: usize) {
        let cap = cap.max(1);
        if let DatasetStorage::Chunked {
            cache_cap,
            cache_bytes,
            cache,
            ..
        } = &mut self.storage
        {
            *cache_cap = cap;
            while cache.len() > *cache_cap {
                if let Some((_, evicted)) = cache.pop_front() {
                    *cache_bytes = cache_bytes.saturating_sub(evicted.len());
                }
            }
        }
    }

    /// Return the configured chunk cache capacity (or `0` for non-chunked datasets).
    #[must_use]
    pub const fn cache_cap(&self) -> usize {
        match &self.storage {
            DatasetStorage::Chunked { cache_cap, .. } => *cache_cap,
            DatasetStorage::Empty | DatasetStorage::Contiguous { .. } => 0,
        }
    }

    /// Total number of logical elements in this dataset.
    #[must_use]
    pub const fn total_elems(&self) -> u64 {
        self.total_elems
    }

    /// Byte size of one dataset element.
    #[must_use]
    pub const fn elem_size(&self) -> usize {
        self.elem_size
    }

    /// Total logical size in bytes (`total_elems * elem_size`).
    #[must_use]
    pub const fn total_bytes(&self) -> u64 {
        self.total_elems.saturating_mul(self.elem_size as u64)
    }

    /// Read `count` consecutive elements starting at `start_elem` into a newly allocated byte
    /// buffer (clamped to `total_elems`).
    ///
    /// # Errors
    /// Returns an `std::io::Error` if seeking, reading, or zlib chunk decompression fails.
    pub fn read_elem_range<R: Read + Seek>(
        &mut self,
        r: &mut R,
        start_elem: u64,
        count: usize,
    ) -> std::io::Result<Vec<u8>> {
        let mut out = Vec::new();
        self.read_elem_range_into(r, start_elem, count, &mut out)?;
        Ok(out)
    }

    /// Read `count` consecutive elements starting at `start_elem` into reusable buffer `out`
    /// (clamped to `total_elems`), avoiding per-call heap allocations on hot segment lookups.
    ///
    /// # Errors
    /// Returns an `std::io::Error` if seeking, reading, or zlib chunk decompression fails.
    #[allow(clippy::too_many_lines)]
    pub fn read_elem_range_into<R: Read + Seek>(
        &mut self,
        r: &mut R,
        start_elem: u64,
        count: usize,
        out: &mut Vec<u8>,
    ) -> std::io::Result<()> {
        out.clear();
        if count == 0 || start_elem >= self.total_elems {
            return Ok(());
        }
        let avail = self.total_elems - start_elem;
        let actual_count = usize::try_from(avail.min(count as u64)).unwrap_or(0);
        if actual_count == 0 {
            return Ok(());
        }
        let out_bytes = actual_count.checked_mul(self.elem_size).ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "range overflow")
        })?;
        if u64::try_from(out_bytes).unwrap_or(u64::MAX)
            > super::hdf5_read::MAX_HDF5_READ_BYTES as u64
        {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "HDF5 read_elem_range exceeds 256 MiB per-call safety cap",
            ));
        }
        let es_u64 = self.elem_size as u64;

        match &mut self.storage {
            DatasetStorage::Empty => {
                out.resize(out_bytes, 0);
                Ok(())
            }
            DatasetStorage::Contiguous { addr, byte_len } => {
                let rel_off = start_elem.checked_mul(es_u64).ok_or_else(|| {
                    std::io::Error::new(std::io::ErrorKind::InvalidData, "offset overflow")
                })?;
                if rel_off >= *byte_len {
                    out.resize(out_bytes, 0);
                    return Ok(());
                }
                let abs_off = addr.checked_add(rel_off).ok_or_else(|| {
                    std::io::Error::new(std::io::ErrorKind::InvalidData, "addr overflow")
                })?;
                let readable = usize::try_from((*byte_len - rel_off).min(out_bytes as u64))
                    .unwrap_or(out_bytes);
                read_at_into(r, abs_off, readable, out)?;
                if out.len() != readable {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::UnexpectedEof,
                        "truncated HDF5 contiguous dataset payload",
                    ));
                }
                out.resize(out_bytes, 0);
                Ok(())
            }
            DatasetStorage::Chunked {
                btree_addr,
                ndim,
                chunk_elems,
                deflate_layers,
                cache_cap,
                cache_bytes,
                cache,
                btree_nodes,
                raw_tmp,
                dec_tmp,
            } => {
                let end_elem = start_elem.saturating_add(actual_count as u64);
                let ce = *chunk_elems;
                // Fast path 1: entire requested slice falls inside one cached chunk (zero B-tree/disk I/O).
                // Scan from back (MRU) to front for O(1) sequential/hot hits.
                if let Some(pos) = cache.iter().rposition(|(cs, buf)| {
                    let c_end = cs.saturating_add((buf.len() as u64) / es_u64);
                    *cs <= start_elem && end_elem <= c_end
                }) {
                    if pos + 1 < cache.len()
                        && let Some(entry) = cache.remove(pos)
                    {
                        cache.push_back(entry);
                    }
                    if let Some((c_start, chunk_data)) = cache.back() {
                        let src_off =
                            usize::try_from((start_elem - *c_start) * es_u64).unwrap_or(0);
                        out.extend_from_slice(&chunk_data[src_off..src_off + out_bytes]);
                        return Ok(());
                    }
                }
                // Fast path 2: slice straddles two adjacent chunks that are BOTH already cached.
                if let Some(p1) = cache.iter().rposition(|(cs, buf)| {
                    let c_end = cs.saturating_add((buf.len() as u64) / es_u64);
                    *cs <= start_elem && start_elem < c_end && c_end < end_elem
                }) {
                    let mid_elem = cache[p1]
                        .0
                        .saturating_add((cache[p1].1.len() as u64) / es_u64);
                    if let Some(p2) = cache.iter().rposition(|(cs, buf)| {
                        let c_end = cs.saturating_add((buf.len() as u64) / es_u64);
                        *cs == mid_elem && end_elem <= c_end
                    }) {
                        let c1_start = cache[p1].0;
                        let s1 = usize::try_from((start_elem - c1_start) * es_u64).unwrap_or(0);
                        let n1 = usize::try_from((mid_elem - start_elem) * es_u64).unwrap_or(0);
                        let n2 = usize::try_from((end_elem - mid_elem) * es_u64).unwrap_or(0);
                        out.reserve(out_bytes);
                        out.extend_from_slice(&cache[p1].1[s1..s1 + n1]);
                        out.extend_from_slice(&cache[p2].1[..n2]);
                        return Ok(());
                    }
                }

                out.resize(out_bytes, 0);
                let mut matching_chunks = Vec::new();
                let mut budget = RangeQueryBudget::new();
                query_chunk_btree_range(
                    r,
                    *btree_addr,
                    *ndim,
                    ce,
                    start_elem,
                    end_elem,
                    None,
                    0,
                    &mut budget,
                    btree_nodes,
                    &mut matching_chunks,
                )?;
                let mut decompressed_budget: u64 = 0;
                for ch in matching_chunks {
                    let chunk_data = load_chunk_cached(
                        r,
                        ch,
                        *deflate_layers,
                        *cache_cap,
                        cache_bytes,
                        cache,
                        raw_tmp,
                        dec_tmp,
                        &mut decompressed_budget,
                    )?;
                    let chunk_n_elems = (chunk_data.len() as u64) / es_u64;
                    let c_start = ch.offset_elem;
                    let c_end = c_start.saturating_add(chunk_n_elems);
                    let ov_s = start_elem.max(c_start);
                    let ov_e = end_elem.min(c_end);
                    if ov_s >= ov_e {
                        continue;
                    }
                    let src_off = usize::try_from((ov_s - c_start) * es_u64).unwrap_or(0);
                    let dst_off = usize::try_from((ov_s - start_elem) * es_u64).unwrap_or(0);
                    let copy_len = usize::try_from((ov_e - ov_s) * es_u64).unwrap_or(0);
                    if src_off + copy_len <= chunk_data.len() && dst_off + copy_len <= out.len() {
                        out[dst_off..dst_off + copy_len]
                            .copy_from_slice(&chunk_data[src_off..src_off + copy_len]);
                    }
                }
                Ok(())
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn load_chunk_cached<'a, R: Read + Seek>(
    r: &mut R,
    ch: Hdf5Chunk,
    deflate_layers: usize,
    cache_cap: usize,
    cache_bytes: &mut usize,
    cache: &'a mut VecDeque<(u64, Vec<u8>)>,
    raw_tmp: &mut Vec<u8>,
    dec_tmp: &mut Vec<u8>,
    decompressed_budget: &mut u64,
) -> std::io::Result<&'a [u8]> {
    if let Some(pos) = cache.iter().rposition(|(off, _)| *off == ch.offset_elem) {
        if pos + 1 < cache.len()
            && let Some(entry) = cache.remove(pos)
        {
            cache.push_back(entry);
        }
        return Ok(cache.back().map_or(&[], |(_, buf)| buf.as_slice()));
    }
    let len = usize::try_from(ch.size).map_err(|_| {
        std::io::Error::new(std::io::ErrorKind::InvalidData, "HDF5 chunk size overflow")
    })?;
    if u64::try_from(len).unwrap_or(u64::MAX) > MAX_SINGLE_CHUNK_BYTES {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "HDF5 chunk exceeds 64 MiB single-chunk safety cap",
        ));
    }
    if len > 0 && ch.addr == u64::MAX {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "HDF5 chunk has HADDR_UNDEF",
        ));
    }
    read_at_into(r, ch.addr, len, raw_tmp)?;
    if raw_tmp.len() != len {
        return Err(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            "truncated HDF5 chunk payload",
        ));
    }
    let active_layers: usize = (0..deflate_layers)
        .filter(|&layer| (ch.filter_mask & (1u32 << layer)) == 0)
        .count();
    let payload = if active_layers > 0 {
        for _ in 0..active_layers {
            dec_tmp.clear();
            let dec = flate2::read::ZlibDecoder::new(raw_tmp.as_slice());
            let mut bounded = dec.take(MAX_SINGLE_CHUNK_BYTES.saturating_add(1));
            let n = bounded.read_to_end(dec_tmp)?;
            let n_u64 = u64::try_from(n).unwrap_or(u64::MAX);
            if n_u64 > MAX_SINGLE_CHUNK_BYTES {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "HDF5 decompressed chunk exceeds 64 MiB single-chunk safety cap",
                ));
            }
            *decompressed_budget = decompressed_budget.saturating_add(n_u64);
            if *decompressed_budget > MAX_DECOMPRESSED_PER_CALL {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "HDF5 cumulative decompressed chunk bytes exceeded 256 MiB per-call safety cap",
                ));
            }
            std::mem::swap(raw_tmp, dec_tmp);
        }
        std::mem::take(raw_tmp)
    } else {
        *decompressed_budget = decompressed_budget.saturating_add(len as u64);
        if *decompressed_budget > MAX_DECOMPRESSED_PER_CALL {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "HDF5 cumulative chunk payload exceeded 256 MiB per-call safety cap",
            ));
        }
        std::mem::take(raw_tmp)
    };
    let payload_len = payload.len();
    let min_keep = cache_cap.min(2).saturating_sub(1);
    while cache.len() >= cache_cap
        || (cache.len() > min_keep
            && cache_bytes.saturating_add(payload_len) > MAX_DATASET_CACHE_BYTES)
    {
        if let Some((_, evicted)) = cache.pop_front() {
            *cache_bytes = cache_bytes.saturating_sub(evicted.len());
        } else {
            break;
        }
    }
    *cache_bytes = cache_bytes.saturating_add(payload_len);
    cache.push_back((ch.offset_elem, payload));
    Ok(cache.back().map_or(&[], |(_, buf)| buf.as_slice()))
}
