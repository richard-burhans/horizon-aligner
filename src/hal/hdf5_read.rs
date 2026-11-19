//! Pure-Rust HDF5 v0/v1 container and dataset readers (slice and Seek).

use super::is_hdf5_hal;

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

pub(super) fn read_u64_le(bytes: &[u8], off: usize) -> Option<u64> {
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
    let Some(magic_end) = snod_addr.checked_add(4) else {
        return;
    };
    if bytes.get(snod_addr..magic_end) != Some(b"SNOD") {
        return;
    }
    let Some(n_syms) = snod_addr
        .checked_add(8)
        .and_then(|end| bytes.get(snod_addr + 6..end))
        .map(|b| usize::from(u16::from_le_bytes([b[0], b[1]])))
    else {
        return;
    };
    for i in 0..n_syms {
        let Some(entry_off) = i.checked_mul(40).and_then(|o| snod_addr.checked_add(8 + o)) else {
            break;
        };
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
    visited: &mut std::collections::HashSet<usize>,
    depth: usize,
) {
    if depth > 64 || !visited.insert(node_addr) {
        return;
    }
    let Some(magic_end) = node_addr.checked_add(8) else {
        return;
    };
    if bytes.get(node_addr..node_addr + 4) != Some(b"TREE") || magic_end > bytes.len() {
        return;
    }
    let node_type = bytes[node_addr + 4];
    let node_level = bytes[node_addr + 5];
    let n_entries = usize::from(u16::from_le_bytes([
        bytes[node_addr + 6],
        bytes[node_addr + 7],
    ]));
    if node_type != 0 {
        return;
    }
    // B-tree v1 header with 8-byte offsets: 8 bytes prefix + 2 * 8 bytes sibling addrs = 24 bytes.
    // Entries are alternating (key: u64, child_addr: u64) for `n_entries`.
    for i in 0..n_entries {
        let Some(child_off) = i
            .checked_mul(16)
            .and_then(|o| node_addr.checked_add(24 + o + 8))
        else {
            break;
        };
        let Some(child_u64) = read_u64_le(bytes, child_off) else {
            break;
        };
        let Ok(child_addr) = usize::try_from(child_u64) else {
            continue;
        };
        if node_level == 0 {
            collect_snod_entries(bytes, child_addr, heap_data, out);
        } else {
            walk_btree_v1_group_entries(bytes, child_addr, heap_data, out, visited, depth + 1);
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
    let heap_end = heap_addr.checked_add(4)?;
    if bytes.get(heap_addr..heap_end) != Some(b"HEAP") {
        return None;
    }
    let data_seg_size = usize::try_from(read_u64_le(bytes, heap_addr.checked_add(8)?)?).ok()?;
    let data_seg_addr = usize::try_from(read_u64_le(bytes, heap_addr.checked_add(24)?)?).ok()?;
    let seg_end = data_seg_addr.checked_add(data_seg_size)?;
    let heap_data = bytes.get(data_seg_addr..seg_end)?;
    let mut entries = Vec::new();
    let mut visited = std::collections::HashSet::new();
    walk_btree_v1_group_entries(bytes, btree_addr, heap_data, &mut entries, &mut visited, 0);
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

const MAX_HDF5_READ_BYTES: usize = 256 * 1024 * 1024;

pub(super) fn read_at<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    offset: u64,
    len: usize,
) -> std::io::Result<Vec<u8>> {
    use std::io::Read;
    if len > MAX_HDF5_READ_BYTES {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("HDF5 read length {len} exceeds limit {MAX_HDF5_READ_BYTES}"),
        ));
    }
    r.seek(std::io::SeekFrom::Start(offset))?;
    let mut buf = Vec::with_capacity(len.min(64 * 1024));
    r.take(len as u64).read_to_end(&mut buf)?;
    Ok(buf)
}

const MAX_HDF5_META_BYTES: usize = 16 * 1024 * 1024;

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
    let hdr_size = usize::try_from(u32::from_le_bytes([
        prefix[8], prefix[9], prefix[10], prefix[11],
    ]))
    .unwrap_or(0);
    if hdr_size > MAX_HDF5_META_BYTES {
        return Ok(None);
    }
    let Some(msg_off) = obj_addr.checked_add(16) else {
        return Ok(None);
    };
    let first_buf = read_at(r, msg_off, hdr_size)?;
    let mut chunks = std::collections::VecDeque::from([first_buf]);
    let mut visited_addrs = std::collections::HashSet::from([obj_addr]);
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
                    && cont_len <= MAX_HDF5_META_BYTES
                    && visited_addrs.len() < 1024
                    && visited_addrs.insert(cont_addr)
                {
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
    Ok(())
}

fn walk_btree_v1_group_entries_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    node_addr: u64,
    heap_data: &[u8],
    out: &mut Vec<Hdf5GroupEntry>,
    visited: &mut std::collections::HashSet<u64>,
    depth: usize,
) -> std::io::Result<()> {
    if depth > 64 || !visited.insert(node_addr) {
        return Ok(());
    }
    let hdr = read_at(r, node_addr, 24)?;
    if hdr.len() < 24 || &hdr[0..4] != b"TREE" || hdr[4] != 0 {
        return Ok(());
    }
    let node_level = hdr[5];
    let n_entries = usize::from(u16::from_le_bytes([hdr[6], hdr[7]]));
    let body = read_at(
        r,
        node_addr.saturating_add(24),
        n_entries.saturating_mul(16),
    )?;
    for i in 0..n_entries {
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
    })?
    else {
        return Ok(None);
    };
    Ok(result)
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

fn walk_chunk_btree_v1_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    node_addr: u64,
    ndim: usize,
    out: &mut Vec<Hdf5Chunk>,
    visited: &mut std::collections::HashSet<u64>,
    depth: usize,
) -> std::io::Result<()> {
    if depth > 64 || !visited.insert(node_addr) {
        return Ok(());
    }
    let hdr = read_at(r, node_addr, 24)?;
    if hdr.len() < 24 || &hdr[0..4] != b"TREE" || hdr[4] != 1 {
        return Ok(());
    }
    let node_level = hdr[5];
    let n_entries = usize::from(u16::from_le_bytes([hdr[6], hdr[7]]));
    let key_bytes = 8usize.saturating_add(ndim.saturating_add(1).saturating_mul(8));
    let rec_bytes = key_bytes.saturating_add(8);
    let body = read_at(
        r,
        node_addr.saturating_add(24),
        n_entries.saturating_mul(rec_bytes),
    )?;
    for i in 0..n_entries {
        let key_off = i * rec_bytes;
        let Some(sz_bytes) = body.get(key_off..key_off + 4) else {
            break;
        };
        let chunk_size = u32::from_le_bytes([sz_bytes[0], sz_bytes[1], sz_bytes[2], sz_bytes[3]]);
        let Some(child_addr) = read_u64_le(&body, key_off + key_bytes) else {
            break;
        };
        if node_level == 0 {
            out.push(Hdf5Chunk {
                size: chunk_size,
                addr: child_addr,
            });
        } else {
            walk_chunk_btree_v1_seek(r, child_addr, ndim, out, visited, depth + 1)?;
        }
    }
    Ok(())
}

/// Read raw dataset bytes from an HDF5 v0 dataset object header via `Seek`, supporting both
/// contiguous (`DataLayout` class `1`) and chunked (`DataLayout` class `2`, optionally zlib-deflated) storage.
///
/// # Errors
/// Returns an `std::io::Error` if seeking, reading, or zlib decompression fails.
pub fn read_hdf5_v0_dataset_bytes_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    obj_addr: u64,
) -> std::io::Result<Option<Vec<u8>>> {
    let mut contig: Option<(u64, usize)> = None;
    let mut chunk_btree: Option<(u64, usize)> = None;
    let mut deflated = false;

    let Some(()) = iter_v1_header_messages_seek(r, obj_addr, |mtype, body| {
        if mtype == 0x0008 && body.len() >= 18 && body[0] == 3 && body[1] == 1 {
            if let (Some(addr), Some(sz_u64)) = (read_u64_le(body, 2), read_u64_le(body, 10))
                && let Ok(sz) = usize::try_from(sz_u64)
            {
                contig = Some((addr, sz));
            }
        } else if mtype == 0x0008 && body.len() >= 11 && body[0] == 3 && body[1] == 2 {
            if let Some(bt_addr) = read_u64_le(body, 3) {
                chunk_btree = Some((bt_addr, usize::from(body[2])));
            }
        } else if mtype == 0x000B && body.len() >= 4 {
            let ver = body[0];
            let n_filters = usize::from(body[1]);
            let mut fpos = match ver {
                1 => 8usize,
                2 => 2usize,
                _ => return std::ops::ControlFlow::Continue(()),
            };
            for _ in 0..n_filters {
                if fpos + 2 > body.len() {
                    break;
                }
                let filter_id = u16::from_le_bytes([body[fpos], body[fpos + 1]]);
                if filter_id == 1 {
                    deflated = true;
                    break;
                }
                if ver == 1 && fpos + 8 <= body.len() {
                    let name_len =
                        usize::from(u16::from_le_bytes([body[fpos + 2], body[fpos + 3]]));
                    let n_client =
                        usize::from(u16::from_le_bytes([body[fpos + 6], body[fpos + 7]]));
                    fpos += 8 + ((name_len + 7) & !7) + (((n_client * 4) + 7) & !7);
                } else if ver == 2 && fpos + 6 <= body.len() {
                    let n_client =
                        usize::from(u16::from_le_bytes([body[fpos + 4], body[fpos + 5]]));
                    let name_len = if filter_id >= 256 && fpos + 8 <= body.len() {
                        usize::from(u16::from_le_bytes([body[fpos + 6], body[fpos + 7]])) + 2
                    } else {
                        0
                    };
                    fpos += 6 + name_len + n_client * 4;
                } else {
                    break;
                }
            }
        }
        std::ops::ControlFlow::Continue(())
    })?
    else {
        return Ok(None);
    };

    if let Some((addr, size)) = contig {
        if u64::try_from(size).unwrap_or(u64::MAX) > MAX_CHUNK_DECOMPRESSED_BYTES {
            return Ok(None);
        }
        return read_at(r, addr, size).map(Some);
    }
    let Some((btree_addr, ndim)) = chunk_btree else {
        return Ok(None);
    };
    let mut chunks = Vec::new();
    let mut visited = std::collections::HashSet::new();
    walk_chunk_btree_v1_seek(r, btree_addr, ndim, &mut chunks, &mut visited, 0)?;
    let mut out = Vec::new();
    for ch in chunks {
        let Ok(len) = usize::try_from(ch.size) else {
            return Ok(None);
        };
        if u64::try_from(len).unwrap_or(u64::MAX) > MAX_CHUNK_DECOMPRESSED_BYTES {
            return Ok(None);
        }
        let raw_chunk = read_at(r, ch.addr, len)?;
        if deflated {
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
    Ok(Some(out))
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
            && let Some(data_end) = data_addr.checked_add(data_size)
        {
            result = bytes.get(data_addr..data_end);
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

fn walk_chunk_btree_v1(
    bytes: &[u8],
    node_addr: usize,
    ndim: usize,
    out: &mut Vec<Hdf5Chunk>,
    visited: &mut std::collections::HashSet<usize>,
    depth: usize,
) {
    if depth > 64 || !visited.insert(node_addr) {
        return;
    }
    let Some(magic_end) = node_addr.checked_add(8) else {
        return;
    };
    if bytes.get(node_addr..node_addr + 4) != Some(b"TREE") || magic_end > bytes.len() {
        return;
    }
    let node_type = bytes[node_addr + 4];
    let node_level = bytes[node_addr + 5];
    let n_entries = usize::from(u16::from_le_bytes([
        bytes[node_addr + 6],
        bytes[node_addr + 7],
    ]));
    if node_type != 1 {
        return;
    }
    // Type-1 chunk key: 4 bytes chunk_size + 4 bytes filter_mask + (ndim + 1) * 8 bytes coords.
    let key_bytes = 8usize.saturating_add(ndim.saturating_add(1).saturating_mul(8));
    let rec_bytes = key_bytes.saturating_add(8);
    for i in 0..n_entries {
        let Some(key_off) = i
            .checked_mul(rec_bytes)
            .and_then(|o| node_addr.checked_add(24 + o))
        else {
            break;
        };
        let Some(sz_end) = key_off.checked_add(4) else {
            break;
        };
        let Some(sz_bytes) = bytes.get(key_off..sz_end) else {
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
            walk_chunk_btree_v1(bytes, sub, ndim, out, visited, depth + 1);
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
            let mut visited = std::collections::HashSet::new();
            walk_chunk_btree_v1(bytes, btree_addr, ndim, &mut chunks, &mut visited, 0);
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
        let end = start.checked_add(len)?;
        let raw_chunk = bytes.get(start..end)?;
        if deflated {
            use std::io::Read as _;
            let rem = MAX_CHUNK_DECOMPRESSED_BYTES
                .saturating_sub(u64::try_from(out.len()).unwrap_or(u64::MAX));
            let dec = flate2::read::ZlibDecoder::new(raw_chunk);
            let mut bounded = dec.take(rem + 1);
            let n = bounded.read_to_end(&mut out).ok()?;
            if u64::try_from(n).unwrap_or(u64::MAX) > rem {
                return None;
            }
        } else {
            if u64::try_from(out.len().saturating_add(raw_chunk.len())).unwrap_or(u64::MAX)
                > MAX_CHUNK_DECOMPRESSED_BYTES
            {
                return None;
            }
            out.extend_from_slice(raw_chunk);
        }
    }
    Some(out)
}
