//! Pure-Rust HDF5 v0/v1 container and dataset readers (unified over `Read + Seek` and `Cursor<&[u8]>`).

use super::is_hdf5_hal;
use std::io::Cursor;

const MAX_HDF5_READ_BYTES: usize = 256 * 1024 * 1024;
const MAX_HDF5_META_BYTES: usize = 16 * 1024 * 1024;
const MAX_HDF5_NAME_BYTES: usize = 4096;
const MAX_CHUNK_DECOMPRESSED_BYTES: u64 = 256 * 1024 * 1024;

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

/// A chunked dataset descriptor (`DataLayout` `0x0008` class `2`, `TREE` type `1`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Hdf5Chunk {
    /// Compressed (or raw) chunk size in bytes.
    size: u32,
    /// HDF5 filter mask (bit `i` set means filter `i` was skipped).
    filter_mask: u32,
    /// File byte offset of the chunk payload.
    addr: u64,
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
    Some(u64::from_le_bytes(
        bytes.get(off..off + 8)?.try_into().ok()?,
    ))
}

fn read_at_into<R: std::io::Read + std::io::Seek>(
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

fn iter_v1_header_messages_seek<R: std::io::Read + std::io::Seek, F>(
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

/// Parse an HDF5 v1 or v2 object header starting at `slice` (positioned at `root_group_addr`).
///
/// Note: v1 Object Header Continuation (`0x0010`) messages store file-absolute offsets, so v1
/// continuations are only resolved accurately when `slice` starts at file offset 0 (for general
/// file-backed reads, prefer `read_hdf5_v0_dataset_bytes_seek`).
#[must_use]
pub fn read_hdf5_object_header(slice: &[u8]) -> Option<Hdf5ObjectHeader> {
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
        let mut cursor = Cursor::new(slice);
        iter_v1_header_messages_seek(&mut cursor, 0, |mtype, _| {
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

fn collect_snod_entries_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    snod_addr: u64,
    heap_data: &[u8],
    out: &mut Vec<Hdf5GroupEntry>,
) -> std::io::Result<()> {
    let hdr = read_at(r, snod_addr, 8)?;
    if hdr.len() < 8 || &hdr[0..4] != b"SNOD" {
        return Ok(());
    }
    let n_syms = usize::from(u16::from_le_bytes([hdr[6], hdr[7]]));
    let Some(entries_off) = snod_addr.checked_add(8) else {
        return Ok(());
    };
    let entries_buf = read_at(r, entries_off, n_syms * 40)?;
    for i in 0..n_syms {
        let Some(name_off_u64) = read_u64_le(&entries_buf, i * 40) else {
            break;
        };
        let Some(obj_addr) = read_u64_le(&entries_buf, i * 40 + 8) else {
            break;
        };
        let Ok(name_off) = usize::try_from(name_off_u64) else {
            continue;
        };
        if let Some(tail) = heap_data.get(name_off..) {
            let scan = &tail[..tail.len().min(MAX_HDF5_NAME_BYTES)];
            let end = scan.iter().position(|&b| b == 0).unwrap_or(scan.len());
            if let Ok(s) = std::str::from_utf8(&scan[..end])
                && !s.is_empty()
            {
                out.push(Hdf5GroupEntry {
                    name: s.to_owned(),
                    object_header_addr: obj_addr,
                });
            }
        }
    }
    Ok(())
}

const MAX_BTREE_NODES: usize = 65_536;
const MAX_BTREE_RECORDS: usize = 1_048_576;

fn walk_btree_v1_group_entries_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    node_addr: u64,
    heap_data: &[u8],
    out: &mut Vec<Hdf5GroupEntry>,
    visited: &mut std::collections::HashSet<u64>,
    depth: usize,
) -> std::io::Result<()> {
    if depth > 64
        || visited.len() >= MAX_BTREE_NODES
        || out.len() >= MAX_BTREE_RECORDS
        || !visited.insert(node_addr)
    {
        return Ok(());
    }
    let hdr = read_at(r, node_addr, 24)?;
    if hdr.len() < 24 || &hdr[0..4] != b"TREE" || hdr[4] != 0 {
        return Ok(());
    }
    let node_level = hdr[5];
    let n_entries = usize::from(u16::from_le_bytes([hdr[6], hdr[7]]));
    let Some(body_off) = node_addr.checked_add(24) else {
        return Ok(());
    };
    let body = read_at(r, body_off, n_entries.saturating_mul(16))?;
    for i in 0..n_entries {
        if out.len() >= MAX_BTREE_RECORDS {
            break;
        }
        let Some(child_addr) = read_u64_le(&body, i * 16 + 8) else {
            break;
        };
        if node_level == 0 {
            collect_snod_entries_seek(r, child_addr, heap_data, out)?;
        } else {
            walk_btree_v1_group_entries_seek(r, child_addr, heap_data, out, visited, depth + 1)?;
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
    let mut entries = Vec::new();
    let mut visited = std::collections::HashSet::new();
    walk_btree_v1_group_entries_seek(r, btree_addr, &heap_data, &mut entries, &mut visited, 0)?;
    Ok(Some(entries))
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

fn parse_dataspace_dims(body: &[u8]) -> Option<Vec<u64>> {
    if body.len() < 4 {
        return None;
    }
    let ds_ver = body[0];
    let rank = usize::from(body[1]);
    let dims_off = match ds_ver {
        1 => 8usize,
        2 => 4usize,
        _ => return None,
    };
    let mut dims = Vec::with_capacity(rank);
    for d in 0..rank {
        dims.push(read_u64_le(body, dims_off + d * 8)?);
    }
    Some(dims)
}

/// Parse HDF5 `Dataspace` (`0x0001`) dimension sizes from an object header (`v1`) via `Seek`.
///
/// # Errors
/// Returns an `std::io::Error` if seeking or reading fails.
pub fn read_hdf5_dataspace_dims_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    obj_addr: u64,
) -> std::io::Result<Option<Vec<u64>>> {
    let mut result = None;
    let Some(()) = iter_v1_header_messages_seek(r, obj_addr, |mtype, body| {
        if mtype == 0x0001 {
            result = parse_dataspace_dims(body);
            return std::ops::ControlFlow::Break(());
        }
        std::ops::ControlFlow::Continue(())
    })?
    else {
        return Ok(None);
    };
    Ok(result)
}

/// Parse HDF5 `Dataspace` (`0x0001`) dimension sizes from an object header (`v1`).
#[must_use]
pub fn read_hdf5_dataspace_dims(bytes: &[u8], obj_addr: usize) -> Option<Vec<u64>> {
    let mut cursor = Cursor::new(bytes);
    read_hdf5_dataspace_dims_seek(&mut cursor, obj_addr as u64)
        .ok()
        .flatten()
}

/// Read element byte size from an HDF5 `Datatype` (`0x0003`) message via `Seek`.
///
/// # Errors
/// Returns an `std::io::Error` if seeking or reading fails.
pub fn read_hdf5_datatype_size_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    obj_addr: u64,
) -> std::io::Result<Option<usize>> {
    let mut size = None;
    let Some(()) = iter_v1_header_messages_seek(r, obj_addr, |mtype, body| {
        if mtype == 0x0003 && body.len() >= 8 {
            let sz = u32::from_le_bytes([body[4], body[5], body[6], body[7]]);
            size = usize::try_from(sz).ok();
            return std::ops::ControlFlow::Break(());
        }
        std::ops::ControlFlow::Continue(())
    })?
    else {
        return Ok(None);
    };
    Ok(size)
}

/// Read element byte size from an HDF5 `Datatype` (`0x0003`) message in an object header (`v1`).
#[must_use]
pub fn read_hdf5_datatype_size(bytes: &[u8], obj_addr: usize) -> Option<usize> {
    let mut cursor = Cursor::new(bytes);
    read_hdf5_datatype_size_seek(&mut cursor, obj_addr as u64)
        .ok()
        .flatten()
}

fn walk_chunk_btree_v1_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    node_addr: u64,
    ndim: usize,
    out: &mut Vec<Hdf5Chunk>,
    visited: &mut std::collections::HashSet<u64>,
    depth: usize,
) -> std::io::Result<()> {
    if depth > 64
        || visited.len() >= MAX_BTREE_NODES
        || out.len() >= MAX_BTREE_RECORDS
        || !visited.insert(node_addr)
    {
        return Ok(());
    }
    let hdr = read_at(r, node_addr, 24)?;
    if hdr.len() < 24 || &hdr[0..4] != b"TREE" || hdr[4] != 1 {
        return Ok(());
    }
    let node_level = hdr[5];
    let n_entries = usize::from(u16::from_le_bytes([hdr[6], hdr[7]]));
    // DataLayout v3 class 2 `ndim` already equals `rank + 1` (includes element-size dimension),
    // and B-tree v1 type-1 chunk keys store `8 + ndim * 8` bytes before `child_addr`.
    let key_bytes = 8usize.saturating_add(ndim.saturating_mul(8));
    let rec_bytes = key_bytes.saturating_add(8);
    let Some(body_off) = node_addr.checked_add(24) else {
        return Ok(());
    };
    let body = read_at(r, body_off, n_entries.saturating_mul(rec_bytes))?;
    for i in 0..n_entries {
        if out.len() >= MAX_BTREE_RECORDS {
            break;
        }
        let key_off = i * rec_bytes;
        let Some(sz_bytes) = body.get(key_off..key_off + 8) else {
            break;
        };
        let chunk_size = u32::from_le_bytes([sz_bytes[0], sz_bytes[1], sz_bytes[2], sz_bytes[3]]);
        let filter_mask = u32::from_le_bytes([sz_bytes[4], sz_bytes[5], sz_bytes[6], sz_bytes[7]]);
        let Some(child_addr) = read_u64_le(&body, key_off + key_bytes) else {
            break;
        };
        if node_level == 0 {
            out.push(Hdf5Chunk {
                size: chunk_size,
                filter_mask,
                addr: child_addr,
            });
        } else {
            walk_chunk_btree_v1_seek(r, child_addr, ndim, out, visited, depth + 1)?;
        }
    }
    Ok(())
}

#[derive(Default)]
struct FilterSummary {
    deflated: bool,
    has_unsupported: bool,
}

fn parse_filter_pipeline(body: &[u8]) -> FilterSummary {
    let mut summary = FilterSummary::default();
    if body.len() < 4 {
        summary.has_unsupported = true;
        return summary;
    }
    let ver = body[0];
    let n_filters = usize::from(body[1]);
    let mut fpos = match ver {
        1 if body.len() >= 8 => 8usize,
        2 => 2usize,
        _ => {
            summary.has_unsupported = true;
            return summary;
        }
    };
    for _ in 0..n_filters {
        if fpos + 2 > body.len() {
            summary.has_unsupported = true;
            break;
        }
        let filter_id = u16::from_le_bytes([body[fpos], body[fpos + 1]]);
        if filter_id == 1 {
            summary.deflated = true;
        } else {
            summary.has_unsupported = true;
        }
        if ver == 1 && fpos + 8 <= body.len() {
            let name_len = usize::from(u16::from_le_bytes([body[fpos + 2], body[fpos + 3]]));
            let n_client = usize::from(u16::from_le_bytes([body[fpos + 6], body[fpos + 7]]));
            fpos += 8 + ((name_len + 7) & !7) + (((n_client * 4) + 7) & !7);
        } else if ver == 2 && fpos + 6 <= body.len() {
            let (name_len, n_client, hdr_len) = if filter_id >= 256 {
                if fpos + 8 > body.len() {
                    summary.has_unsupported = true;
                    break;
                }
                (
                    usize::from(u16::from_le_bytes([body[fpos + 2], body[fpos + 3]])),
                    usize::from(u16::from_le_bytes([body[fpos + 6], body[fpos + 7]])),
                    8usize,
                )
            } else {
                (
                    0usize,
                    usize::from(u16::from_le_bytes([body[fpos + 4], body[fpos + 5]])),
                    6usize,
                )
            };
            fpos += hdr_len + name_len + n_client * 4;
        } else {
            summary.has_unsupported = true;
            break;
        }
    }
    summary
}

/// Check whether an HDF5 v1 object header enables the deflate (`id == 1`) filter in `FilterPipeline` (`0x000B`).
#[must_use]
pub fn is_hdf5_deflated(bytes: &[u8], obj_addr: usize) -> bool {
    let mut deflated = false;
    let mut cursor = Cursor::new(bytes);
    let _ = iter_v1_header_messages_seek(&mut cursor, obj_addr as u64, |mtype, body| {
        if mtype == 0x000B {
            let fs = parse_filter_pipeline(body);
            if fs.deflated && !fs.has_unsupported {
                deflated = true;
                return std::ops::ControlFlow::Break(());
            }
        }
        std::ops::ControlFlow::Continue(())
    });
    deflated
}

/// Read raw dataset bytes from an HDF5 v0 dataset object header via `Seek`, supporting both
/// contiguous (`DataLayout` class `1`) and chunked (`DataLayout` class `2`, optionally zlib-deflated) storage.
///
/// # Errors
/// Returns an `std::io::Error` if seeking, reading, or zlib decompression fails.
#[allow(clippy::too_many_lines)]
pub fn read_hdf5_v0_dataset_bytes_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    obj_addr: u64,
) -> std::io::Result<Option<Vec<u8>>> {
    let mut contig: Option<(u64, usize)> = None;
    let mut chunk_btree: Option<(u64, usize)> = None;
    let mut filters = FilterSummary::default();
    let mut dims: Option<Vec<u64>> = None;
    let mut elem_size: Option<usize> = None;

    let Some(()) = iter_v1_header_messages_seek(r, obj_addr, |mtype, body| {
        if mtype == 0x0001 {
            dims = parse_dataspace_dims(body);
        } else if mtype == 0x0003 && body.len() >= 8 {
            let sz = u32::from_le_bytes([body[4], body[5], body[6], body[7]]);
            elem_size = usize::try_from(sz).ok();
        } else if mtype == 0x0008 && body.len() >= 18 && body[0] == 3 && body[1] == 1 {
            if let (Some(addr), Some(sz_u64)) = (read_u64_le(body, 2), read_u64_le(body, 10))
                && let Ok(sz) = usize::try_from(sz_u64)
            {
                contig = Some((addr, sz));
            }
        } else if mtype == 0x0008 && body.len() >= 11 && body[0] == 3 && body[1] == 2 {
            if let Some(bt_addr) = read_u64_le(body, 3) {
                chunk_btree = Some((bt_addr, usize::from(body[2])));
            }
        } else if mtype == 0x000B {
            let fs = parse_filter_pipeline(body);
            filters.deflated |= fs.deflated;
            filters.has_unsupported |= fs.has_unsupported;
        }
        std::ops::ControlFlow::Continue(())
    })?
    else {
        return Ok(None);
    };

    if filters.has_unsupported {
        return Ok(None);
    }

    if let Some((addr, size)) = contig {
        if size == 0 {
            return Ok(Some(Vec::new()));
        }
        if addr == u64::MAX
            || u64::try_from(size).unwrap_or(u64::MAX) > MAX_CHUNK_DECOMPRESSED_BYTES
        {
            return Ok(None);
        }
        let buf = read_at(r, addr, size)?;
        return if buf.len() == size {
            Ok(Some(buf))
        } else {
            Ok(None)
        };
    }
    let Some((btree_addr, ndim)) = chunk_btree else {
        return Ok(None);
    };
    let mut chunks = Vec::new();
    let mut visited = std::collections::HashSet::new();
    walk_chunk_btree_v1_seek(r, btree_addr, ndim, &mut chunks, &mut visited, 0)?;
    let prealloc = dims
        .as_ref()
        .zip(elem_size)
        .and_then(|(ds_dims, es)| {
            ds_dims.iter().try_fold(es, |acc, &d| {
                usize::try_from(d).ok().and_then(|du| acc.checked_mul(du))
            })
        })
        .unwrap_or(0)
        .min(usize::try_from(MAX_CHUNK_DECOMPRESSED_BYTES).unwrap_or(usize::MAX));
    let mut out = Vec::with_capacity(prealloc);
    let mut raw_chunk = Vec::new();
    for ch in chunks {
        let Ok(len) = usize::try_from(ch.size) else {
            return Ok(None);
        };
        if (len > 0 && ch.addr == u64::MAX)
            || u64::try_from(len).unwrap_or(u64::MAX) > MAX_CHUNK_DECOMPRESSED_BYTES
        {
            return Ok(None);
        }
        read_at_into(r, ch.addr, len, &mut raw_chunk)?;
        if raw_chunk.len() != len {
            return Ok(None);
        }
        if filters.deflated && (ch.filter_mask & 1 == 0) {
            use std::io::Read as _;
            let rem = MAX_CHUNK_DECOMPRESSED_BYTES
                .saturating_sub(u64::try_from(out.len()).unwrap_or(u64::MAX));
            let dec = flate2::read::ZlibDecoder::new(raw_chunk.as_slice());
            let mut bounded = dec.take(rem + 1);
            let n = bounded.read_to_end(&mut out)?;
            if u64::try_from(n).unwrap_or(u64::MAX) > rem {
                return Ok(None);
            }
        } else {
            if u64::try_from(out.len().saturating_add(raw_chunk.len())).unwrap_or(u64::MAX)
                > MAX_CHUNK_DECOMPRESSED_BYTES
            {
                return Ok(None);
            }
            out.extend_from_slice(&raw_chunk);
        }
    }
    if let (Some(ds_dims), Some(es)) = (dims, elem_size)
        && let Some(logical_bytes) = ds_dims.iter().try_fold(es, |acc, &d| {
            usize::try_from(d).ok().and_then(|du| acc.checked_mul(du))
        })
        && out.len() > logical_bytes
    {
        out.truncate(logical_bytes);
    }
    Ok(Some(out))
}

/// Read raw dataset bytes from an HDF5 v0 dataset object header, supporting both contiguous
/// (`DataLayout` class `1`) and chunked (`DataLayout` class `2`, optionally zlib-deflated) storage.
#[must_use]
pub fn read_hdf5_v0_dataset_bytes(bytes: &[u8], obj_addr: usize) -> Option<Vec<u8>> {
    let mut cursor = Cursor::new(bytes);
    read_hdf5_v0_dataset_bytes_seek(&mut cursor, obj_addr as u64)
        .ok()
        .flatten()
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
