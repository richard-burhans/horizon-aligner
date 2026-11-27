//! Pure-Rust HDF5 v0/v1 container and dataset readers (unified over `Read + Seek` and `Cursor<&[u8]>`).

use super::is_hdf5_hal;
use std::io::Cursor;

pub(super) const MAX_HDF5_READ_BYTES: usize = 256 * 1024 * 1024;
const MAX_HDF5_META_BYTES: usize = 16 * 1024 * 1024;
const MAX_HDF5_NAME_BYTES: usize = 4096;
pub(super) const MAX_CHUNK_DECOMPRESSED_BYTES: u64 = 256 * 1024 * 1024;

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
        0 => {
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
        1 => {
            // Superblock v1 inserts 4 bytes at 20..24 (Indexed Storage Internal Node K + reserved),
            // shifting the root-group symbol-table entry's object header address to 68..76.
            let offset_size = bytes[13];
            let length_size = bytes[14];
            let root_group_addr = if offset_size == 8 && length_size == 8 && bytes.len() >= 76 {
                Some(u64::from_le_bytes(bytes[68..76].try_into().ok()?))
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

/// Summary of an HDF5 object header (v1 or v2) parsed from a `.hal` file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Hdf5ObjectHeader {
    /// Object header format version (`1` or `2`).
    pub version: u8,
    /// Message type IDs found in the header prefix (e.g., `0x0011` `SymbolTable`, `0x000C` Attribute).
    pub message_types: Vec<u16>,
}

/// A child entry inside an HDF5 v0 group (`SymbolTable` + `TREE`/`SNOD` + `HEAP`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Hdf5GroupEntry {
    /// Member name in the group.
    pub name: String,
    /// File byte offset of the member's HDF5 object header.
    pub object_header_addr: u64,
}

/// Summary of an HDF5 `.hal` file read natively in pure Rust.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalFileInfo {
    /// Parsed HDF5 superblock metadata.
    pub superblock: Hdf5Superblock,
    /// Root-group member names (genomes and `Meta`/`Phylogeny` groups) for HDF5 v0 files.
    pub root_names: Vec<String>,
}

pub(super) fn read_u64_le(bytes: &[u8], off: usize) -> Option<u64> {
    let end = off.checked_add(8)?;
    Some(u64::from_le_bytes(bytes.get(off..end)?.try_into().ok()?))
}

pub(super) fn read_at_into<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    offset: u64,
    len: usize,
    buf: &mut Vec<u8>,
) -> std::io::Result<()> {
    use std::io::Read as _;
    if len > MAX_HDF5_READ_BYTES {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("HDF5 read length {len} exceeds limit {MAX_HDF5_READ_BYTES}"),
        ));
    }
    buf.clear();
    // Reject HDF5 HADDR_UNDEF (u64::MAX) and offsets beyond signed 64-bit file seek limits so
    // OS-backed File handles and in-memory Cursors behave identically instead of raising EINVAL.
    if offset == u64::MAX || offset > i64::MAX as u64 {
        return Ok(());
    }
    r.seek(std::io::SeekFrom::Start(offset))?;
    r.take(len as u64).read_to_end(buf)?;
    Ok(())
}

pub(super) fn read_at<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    offset: u64,
    len: usize,
) -> std::io::Result<Vec<u8>> {
    let mut buf = Vec::with_capacity(len.min(64 * 1024));
    read_at_into(r, offset, len, &mut buf)?;
    Ok(buf)
}

pub(super) fn iter_v1_header_messages_seek<R: std::io::Read + std::io::Seek, F>(
    r: &mut R,
    obj_addr: u64,
    mut f: F,
) -> std::io::Result<Option<()>>
where
    F: FnMut(u16, &[u8]) -> std::ops::ControlFlow<()>,
{
    let prefix = read_at(r, obj_addr, 16)?;
    if prefix.len() < 16 || prefix[0] != 1 {
        return Ok(None);
    }
    let num_msgs = usize::from(u16::from_le_bytes([prefix[2], prefix[3]]));
    let Ok(hdr_size) = usize::try_from(u32::from_le_bytes([
        prefix[8], prefix[9], prefix[10], prefix[11],
    ])) else {
        return Ok(None);
    };
    if hdr_size > MAX_HDF5_META_BYTES {
        return Ok(None);
    }
    let Some(msg_off) = obj_addr.checked_add(16) else {
        return Ok(None);
    };
    let first_buf = read_at(r, msg_off, hdr_size)?;
    let mut chunks = std::collections::VecDeque::from([first_buf]);
    let mut visited_addrs = std::collections::HashSet::from([obj_addr]);
    let mut total_meta_bytes = hdr_size;
    let mut seen = 0usize;

    while let Some(buf) = chunks.pop_front() {
        let mut pos = 0usize;
        let end = buf.len();
        while seen < num_msgs && pos + 8 <= end {
            let mtype = u16::from_le_bytes([buf[pos], buf[pos + 1]]);
            let msize = usize::from(u16::from_le_bytes([buf[pos + 2], buf[pos + 3]]));
            let body_start = pos + 8;
            let body_end = body_start.saturating_add(msize).min(end);
            let body = &buf[body_start..body_end];
            seen += 1;
            pos += 8 + msize;

            if mtype == 0x0010 && body.len() >= 16 {
                if let (Some(cont_addr), Some(len_u64)) =
                    (read_u64_le(body, 0), read_u64_le(body, 8))
                    && let Ok(cont_len) = usize::try_from(len_u64)
                    && total_meta_bytes
                        .checked_add(cont_len)
                        .is_some_and(|t| t <= MAX_HDF5_META_BYTES)
                    && visited_addrs.len() < 1024
                    && visited_addrs.insert(cont_addr)
                {
                    total_meta_bytes += cont_len;
                    chunks.push_back(read_at(r, cont_addr, cont_len)?);
                }
                continue;
            }

            if f(mtype, body).is_break() {
                return Ok(Some(()));
            }
        }
    }
    Ok(Some(()))
}

/// Parse an HDF5 v1 or v2 object header at file-absolute offset `obj_addr` inside `bytes`.
///
/// Unlike [`read_hdf5_object_header`], `bytes` is expected to begin at file offset 0 so v1
/// Object Header Continuation (`0x0010`) messages resolve file-absolute offsets accurately.
#[must_use]
pub fn read_hdf5_object_header_at(bytes: &[u8], obj_addr: u64) -> Option<Hdf5ObjectHeader> {
    let start = usize::try_from(obj_addr).ok()?;
    let slice = bytes.get(start..)?;
    if slice.len() < 16 {
        return None;
    }
    if slice.starts_with(b"OHDR") {
        if slice[4] != 2 {
            return None;
        }
        let flags = slice[5];
        let mut pos = 6usize;
        if flags & 0x20 != 0 {
            pos += 16;
        }
        if flags & 0x10 != 0 {
            pos += 4;
        }
        let chunk_sz_bytes = 1usize << (flags & 0x03);
        if pos + chunk_sz_bytes > slice.len() {
            return None;
        }
        let mut chunk_sz_u64 = 0u64;
        for b in 0..chunk_sz_bytes {
            chunk_sz_u64 |= u64::from(slice[pos + b]) << (8 * b);
        }
        let chunk_sz = usize::try_from(chunk_sz_u64).ok()?;
        pos += chunk_sz_bytes;
        let end = pos.saturating_add(chunk_sz).min(slice.len());
        let msg_hdr_len = if flags & 0x04 != 0 { 6usize } else { 4usize };
        let mut msg_types = Vec::new();
        while pos.saturating_add(msg_hdr_len) <= end {
            let mtype = u16::from(slice[pos]);
            let msize = usize::from(u16::from_le_bytes([slice[pos + 1], slice[pos + 2]]));
            msg_types.push(mtype);
            let Some(next) = pos
                .checked_add(msg_hdr_len)
                .and_then(|p| p.checked_add(msize))
            else {
                break;
            };
            pos = next;
        }
        return Some(Hdf5ObjectHeader {
            version: 2,
            message_types: msg_types,
        });
    }
    if slice[0] == 1 {
        let mut msg_types = Vec::new();
        let mut cursor = Cursor::new(bytes);
        iter_v1_header_messages_seek(&mut cursor, obj_addr, |mtype, _| {
            msg_types.push(mtype);
            std::ops::ControlFlow::Continue(())
        })
        .ok()
        .flatten()?;
        return Some(Hdf5ObjectHeader {
            version: 1,
            message_types: msg_types,
        });
    }
    None
}

/// Parse an HDF5 v1 or v2 object header starting at `slice` (positioned at `root_group_addr`).
///
/// Note: v1 Object Header Continuation (`0x0010`) messages store file-absolute offsets; use
/// [`read_hdf5_object_header_at`] when parsing an object header at a non-zero offset inside a
/// full file buffer.
#[must_use]
pub fn read_hdf5_object_header(slice: &[u8]) -> Option<Hdf5ObjectHeader> {
    read_hdf5_object_header_at(slice, 0)
}

pub(super) const MAX_BTREE_NODES: usize = 65_536;
pub(super) const MAX_BTREE_RECORDS: usize = 1_048_576;

pub(super) struct BtreeBudget {
    pub(super) meta_bytes_remaining: usize,
    pub(super) entries_remaining: usize,
    pub(super) rejected: bool,
}

impl BtreeBudget {
    pub(super) const fn new() -> Self {
        Self {
            meta_bytes_remaining: MAX_HDF5_META_BYTES,
            entries_remaining: MAX_BTREE_RECORDS,
            rejected: false,
        }
    }

    pub(super) const fn reject(&mut self) {
        self.rejected = true;
    }

    pub(super) fn consume_bytes(&mut self, bytes: usize) -> bool {
        if bytes > self.meta_bytes_remaining {
            self.rejected = true;
            return false;
        }
        self.meta_bytes_remaining -= bytes;
        true
    }

    pub(super) fn consume_entry(&mut self) -> bool {
        if self.entries_remaining == 0 {
            self.rejected = true;
            return false;
        }
        self.entries_remaining -= 1;
        true
    }
}

struct GroupWalkCtx<'a> {
    heap_data: &'a [u8],
    out: Vec<Hdf5GroupEntry>,
    visited: std::collections::HashSet<u64>,
    name_bytes: usize,
    budget: BtreeBudget,
}

fn collect_snod_entries_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    snod_addr: u64,
    ctx: &mut GroupWalkCtx<'_>,
) -> std::io::Result<()> {
    if ctx.budget.rejected {
        return Ok(());
    }
    if ctx.visited.len() >= MAX_BTREE_NODES
        || ctx.out.len() >= MAX_BTREE_RECORDS
        || ctx.name_bytes >= MAX_HDF5_META_BYTES
        || !ctx.visited.insert(snod_addr)
        || !ctx.budget.consume_bytes(8)
    {
        ctx.budget.reject();
        return Ok(());
    }
    let hdr = read_at(r, snod_addr, 8)?;
    if hdr.len() < 8 || &hdr[0..4] != b"SNOD" {
        ctx.budget.reject();
        return Ok(());
    }
    let n_syms = usize::from(u16::from_le_bytes([hdr[6], hdr[7]]));
    let Some(entries_off) = snod_addr.checked_add(8) else {
        ctx.budget.reject();
        return Ok(());
    };
    let Some(read_len) = n_syms.checked_mul(40) else {
        ctx.budget.reject();
        return Ok(());
    };
    if !ctx.budget.consume_bytes(read_len) {
        ctx.budget.reject();
        return Ok(());
    }
    let entries_buf = read_at(r, entries_off, read_len)?;
    if entries_buf.len() < read_len {
        ctx.budget.reject();
        return Ok(());
    }
    for i in 0..n_syms {
        if ctx.out.len() >= MAX_BTREE_RECORDS
            || ctx.name_bytes >= MAX_HDF5_META_BYTES
            || !ctx.budget.consume_entry()
        {
            ctx.budget.reject();
            return Ok(());
        }
        let Some(name_off_u64) = read_u64_le(&entries_buf, i * 40) else {
            ctx.budget.reject();
            return Ok(());
        };
        let Some(obj_addr) = read_u64_le(&entries_buf, i * 40 + 8) else {
            ctx.budget.reject();
            return Ok(());
        };
        let Ok(name_off) = usize::try_from(name_off_u64) else {
            continue;
        };
        if let Some(tail) = ctx.heap_data.get(name_off..) {
            let scan = &tail[..tail.len().min(MAX_HDF5_NAME_BYTES)];
            let end = scan.iter().position(|&b| b == 0).unwrap_or(scan.len());
            if let Ok(s) = std::str::from_utf8(&scan[..end])
                && !s.is_empty()
            {
                ctx.name_bytes = ctx.name_bytes.saturating_add(s.len());
                if ctx.name_bytes > MAX_HDF5_META_BYTES {
                    ctx.budget.reject();
                    return Ok(());
                }
                ctx.out.push(Hdf5GroupEntry {
                    name: s.to_owned(),
                    object_header_addr: obj_addr,
                });
            }
        }
    }
    Ok(())
}

fn walk_btree_v1_group_entries_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    node_addr: u64,
    ctx: &mut GroupWalkCtx<'_>,
    parent_level: Option<u8>,
    depth: usize,
) -> std::io::Result<()> {
    if ctx.budget.rejected {
        return Ok(());
    }
    if depth > 64
        || ctx.visited.len() >= MAX_BTREE_NODES
        || ctx.out.len() >= MAX_BTREE_RECORDS
        || ctx.name_bytes >= MAX_HDF5_META_BYTES
        || !ctx.visited.insert(node_addr)
        || !ctx.budget.consume_bytes(24)
    {
        ctx.budget.reject();
        return Ok(());
    }
    let hdr = read_at(r, node_addr, 24)?;
    if hdr.len() < 24 || &hdr[0..4] != b"TREE" || hdr[4] != 0 {
        ctx.budget.reject();
        return Ok(());
    }
    let node_level = hdr[5];
    if parent_level.is_some_and(|pl| node_level >= pl) {
        ctx.budget.reject();
        return Ok(());
    }
    let n_entries = usize::from(u16::from_le_bytes([hdr[6], hdr[7]]));
    let Some(body_off) = node_addr.checked_add(24) else {
        ctx.budget.reject();
        return Ok(());
    };
    let Some(read_len) = n_entries.checked_mul(16) else {
        ctx.budget.reject();
        return Ok(());
    };
    if !ctx.budget.consume_bytes(read_len) {
        ctx.budget.reject();
        return Ok(());
    }
    let body = read_at(r, body_off, read_len)?;
    if body.len() < read_len {
        ctx.budget.reject();
        return Ok(());
    }
    let mut children = Vec::with_capacity(n_entries.min(ctx.budget.entries_remaining));
    for i in 0..n_entries {
        if ctx.out.len() >= MAX_BTREE_RECORDS
            || ctx.name_bytes >= MAX_HDF5_META_BYTES
            || !ctx.budget.consume_entry()
        {
            ctx.budget.reject();
            return Ok(());
        }
        let Some(child_addr) = read_u64_le(&body, i * 16 + 8) else {
            ctx.budget.reject();
            return Ok(());
        };
        children.push(child_addr);
    }
    drop(body);
    for child_addr in children {
        if ctx.budget.rejected {
            return Ok(());
        }
        if ctx.out.len() >= MAX_BTREE_RECORDS || ctx.name_bytes >= MAX_HDF5_META_BYTES {
            ctx.budget.reject();
            return Ok(());
        }
        if node_level == 0 {
            collect_snod_entries_seek(r, child_addr, ctx)?;
        } else {
            walk_btree_v1_group_entries_seek(r, child_addr, ctx, Some(node_level), depth + 1)?;
        }
    }
    Ok(())
}

/// Extract child member entries `(name, object_header_addr)` from an HDF5 v0 group via `Seek`.
///
/// # Errors
/// Returns an `std::io::Error` if seeking or reading fails.
pub fn read_hdf5_v0_group_entries_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    group_obj_addr: u64,
) -> std::io::Result<Option<Vec<Hdf5GroupEntry>>> {
    let mut symtab_addrs: Option<(u64, u64)> = None;
    let Some(()) = iter_v1_header_messages_seek(r, group_obj_addr, |mtype, body| {
        if mtype == 0x0011
            && body.len() >= 16
            && let (Some(bt), Some(hp)) = (read_u64_le(body, 0), read_u64_le(body, 8))
        {
            symtab_addrs = Some((bt, hp));
            return std::ops::ControlFlow::Break(());
        }
        std::ops::ControlFlow::Continue(())
    })?
    else {
        return Ok(None);
    };
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
    if data_seg_size > MAX_HDF5_META_BYTES {
        return Ok(None);
    }
    let Some(data_seg_addr) = read_u64_le(&heap_hdr, 24) else {
        return Ok(None);
    };
    let heap_data = read_at(r, data_seg_addr, data_seg_size)?;
    if heap_data.len() < data_seg_size {
        return Ok(None);
    }
    let mut ctx = GroupWalkCtx {
        heap_data: &heap_data,
        out: Vec::new(),
        visited: std::collections::HashSet::new(),
        name_bytes: 0,
        budget: BtreeBudget::new(),
    };
    walk_btree_v1_group_entries_seek(r, btree_addr, &mut ctx, None, 0)?;
    if ctx.budget.rejected {
        return Ok(None);
    }
    Ok(Some(ctx.out))
}

/// Extract child member entries `(name, object_header_addr)` from any HDF5 v0 group object header.
#[must_use]
pub fn read_hdf5_v0_group_entries(
    bytes: &[u8],
    group_obj_addr: usize,
) -> Option<Vec<Hdf5GroupEntry>> {
    let mut cursor = Cursor::new(bytes);
    read_hdf5_v0_group_entries_seek(&mut cursor, group_obj_addr as u64)
        .ok()
        .flatten()
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
    Ok(read_hdf5_v0_group_entries_seek(r, root_addr)?
        .map(|entries| entries.into_iter().map(|e| e.name).collect()))
}

/// Extract root-group member names from an HDF5 v0 `.hal` file via its `SymbolTable` (`0x0011`),
/// `HEAP` local heap, and `TREE`/`SNOD` B-tree v1 nodes.
#[must_use]
pub fn read_hdf5_v0_root_names(bytes: &[u8]) -> Option<Vec<String>> {
    let sb = read_hdf5_superblock(bytes)?;
    let mut cursor = Cursor::new(bytes);
    read_hdf5_v0_root_names_seek(&mut cursor, sb).ok().flatten()
}

/// Inspect a `.hal` file on disk via `Seek` and extract its HDF5 superblock and v0 root member names.
///
/// # Errors
/// Returns an `std::io::Error` if the file cannot be opened or read.
pub fn inspect_hal_genomes(path: &std::path::Path) -> std::io::Result<Option<HalFileInfo>> {
    let mut f = std::io::BufReader::new(std::fs::File::open(path)?);
    let buf = read_at(&mut f, 0, 96)?;
    let Some(superblock) = read_hdf5_superblock(&buf) else {
        return Ok(None);
    };
    let root_names = read_hdf5_v0_root_names_seek(&mut f, superblock)?.unwrap_or_default();
    Ok(Some(HalFileInfo {
        superblock,
        root_names,
    }))
}
