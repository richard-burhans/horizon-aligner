//! Pure-Rust HDF5 v0/v1 dataset, chunk B-tree, and filter-pipeline readers.

use super::hdf5_read::{
    BtreeBudget, MAX_BTREE_NODES, MAX_BTREE_RECORDS, MAX_CHUNK_DECOMPRESSED_BYTES,
    iter_v1_header_messages_seek, read_at, read_at_into, read_u64_le,
};
use std::io::Cursor;

/// A chunked dataset descriptor (`DataLayout` `0x0008` class `2`, `TREE` type `1`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct Hdf5Chunk {
    pub(super) size: u32,
    pub(super) filter_mask: u32,
    pub(super) addr: u64,
    pub(super) offset_elem: u64,
}

pub(super) fn parse_dataspace_dims(body: &[u8]) -> Option<Vec<u64>> {
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

struct ChunkWalkCtx {
    ndim: usize,
    out: Vec<Hdf5Chunk>,
    visited: std::collections::HashSet<u64>,
    budget: BtreeBudget,
}

fn walk_chunk_btree_v1_seek<R: std::io::Read + std::io::Seek>(
    r: &mut R,
    node_addr: u64,
    ctx: &mut ChunkWalkCtx,
    parent_level: Option<u8>,
    depth: usize,
) -> std::io::Result<()> {
    if ctx.budget.rejected {
        return Ok(());
    }
    if depth > 64
        || ctx.visited.len() >= MAX_BTREE_NODES
        || ctx.out.len() >= MAX_BTREE_RECORDS
        || !ctx.visited.insert(node_addr)
        || !ctx.budget.consume_bytes(24)
    {
        ctx.budget.reject();
        return Ok(());
    }
    let hdr = read_at(r, node_addr, 24)?;
    if hdr.len() < 24 || &hdr[0..4] != b"TREE" || hdr[4] != 1 {
        ctx.budget.reject();
        return Ok(());
    }
    let node_level = hdr[5];
    if parent_level.is_some_and(|pl| node_level >= pl) {
        ctx.budget.reject();
        return Ok(());
    }
    let n_entries = usize::from(u16::from_le_bytes([hdr[6], hdr[7]]));
    // DataLayout v3 class 2 `ndim` already equals `rank + 1` (includes element-size dimension),
    // and B-tree v1 type-1 chunk keys store `8 + ndim * 8` bytes before `child_addr`.
    let key_bytes = 8usize.saturating_add(ctx.ndim.saturating_mul(8));
    let rec_bytes = key_bytes.saturating_add(8);
    let Some(body_off) = node_addr.checked_add(24) else {
        ctx.budget.reject();
        return Ok(());
    };
    let Some(read_len) = n_entries.checked_mul(rec_bytes) else {
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
    if node_level == 0 {
        for i in 0..n_entries {
            if ctx.out.len() >= MAX_BTREE_RECORDS || !ctx.budget.consume_entry() {
                ctx.budget.reject();
                return Ok(());
            }
            let key_off = i * rec_bytes;
            let Some(sz_bytes) = body.get(key_off..key_off + 8) else {
                ctx.budget.reject();
                return Ok(());
            };
            let chunk_size =
                u32::from_le_bytes([sz_bytes[0], sz_bytes[1], sz_bytes[2], sz_bytes[3]]);
            let filter_mask =
                u32::from_le_bytes([sz_bytes[4], sz_bytes[5], sz_bytes[6], sz_bytes[7]]);
            let offset_elem = read_u64_le(&body, key_off + 8).unwrap_or(0);
            let Some(child_addr) = read_u64_le(&body, key_off + key_bytes) else {
                ctx.budget.reject();
                return Ok(());
            };
            ctx.out.push(Hdf5Chunk {
                size: chunk_size,
                filter_mask,
                addr: child_addr,
                offset_elem,
            });
        }
    } else {
        let mut child_addrs = Vec::with_capacity(n_entries.min(ctx.budget.entries_remaining));
        for i in 0..n_entries {
            if ctx.out.len() >= MAX_BTREE_RECORDS || !ctx.budget.consume_entry() {
                ctx.budget.reject();
                return Ok(());
            }
            let key_off = i * rec_bytes;
            let Some(child_addr) = read_u64_le(&body, key_off + key_bytes) else {
                ctx.budget.reject();
                return Ok(());
            };
            child_addrs.push(child_addr);
        }
        drop(body);
        for child_addr in child_addrs {
            if ctx.budget.rejected {
                return Ok(());
            }
            if ctx.out.len() >= MAX_BTREE_RECORDS {
                ctx.budget.reject();
                return Ok(());
            }
            walk_chunk_btree_v1_seek(r, child_addr, ctx, Some(node_level), depth + 1)?;
        }
    }
    Ok(())
}

#[derive(Default)]
pub(super) struct FilterSummary {
    pub(super) deflate_layers: usize,
    pub(super) has_unsupported: bool,
}

pub(super) fn parse_filter_pipeline(body: &[u8]) -> FilterSummary {
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
            summary.deflate_layers += 1;
        } else {
            summary.has_unsupported = true;
        }
        if ver == 1 && fpos + 8 <= body.len() {
            let name_len = usize::from(u16::from_le_bytes([body[fpos + 2], body[fpos + 3]]));
            let n_client = usize::from(u16::from_le_bytes([body[fpos + 6], body[fpos + 7]]));
            fpos += 8 + ((name_len + 7) & !7) + (((n_client * 4) + 7) & !7);
            if fpos > body.len() {
                summary.has_unsupported = true;
                break;
            }
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
            if fpos > body.len() {
                summary.has_unsupported = true;
                break;
            }
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
    let mut filters = FilterSummary::default();
    let mut cursor = Cursor::new(bytes);
    let _ = iter_v1_header_messages_seek(&mut cursor, obj_addr as u64, |mtype, body| {
        if mtype == 0x000B {
            let fs = parse_filter_pipeline(body);
            filters.deflate_layers = filters.deflate_layers.saturating_add(fs.deflate_layers);
            filters.has_unsupported |= fs.has_unsupported;
        }
        std::ops::ControlFlow::Continue(())
    });
    !filters.has_unsupported && filters.deflate_layers > 0
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
            elem_size = (sz > 0).then(|| usize::try_from(sz).ok()).flatten();
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
    // Only 1D chunked datasets (ndim == 1 or rank 1 + element-size dim => ndim == 2) are supported;
    // reject multi-dimensional chunked datasets rather than returning scrambled tiles.
    if btree_addr == u64::MAX || ndim == 0 || ndim > 2 {
        return Ok(None);
    }
    let mut ctx = ChunkWalkCtx {
        ndim,
        out: Vec::new(),
        visited: std::collections::HashSet::new(),
        budget: BtreeBudget::new(),
    };
    walk_chunk_btree_v1_seek(r, btree_addr, &mut ctx, None, 0)?;
    if ctx.budget.rejected {
        return Ok(None);
    }
    let chunks = ctx.out;
    if elem_size.is_none() && chunks.iter().any(|c| c.offset_elem > 0) {
        return Ok(None);
    }
    // Compute logical byte size upfront so sparse chunks beyond the declared Dataspace cannot
    // trigger 256 MiB zero-fill amplification before truncation.
    let logical_bytes_u64: Option<u64> = dims.as_ref().zip(elem_size).and_then(|(ds_dims, es)| {
        let es_u64 = u64::try_from(es).ok()?;
        ds_dims
            .iter()
            .try_fold(es_u64, |acc, &d| acc.checked_mul(d))
    });
    // Cap initial preallocation to 1 MiB so untrusted Dataspace dims cannot trigger a 256 MiB OOM
    // before any chunk payload is verified.
    let prealloc = logical_bytes_u64
        .and_then(|b| usize::try_from(b).ok())
        .unwrap_or(0)
        .min(1024 * 1024);
    let mut out = Vec::with_capacity(prealloc);
    let mut raw_chunk = Vec::new();
    let mut tmp_chunk = Vec::new();
    let mut total_inflated: u64 = 0;
    let mut written_spans: std::collections::BTreeMap<usize, usize> =
        std::collections::BTreeMap::new();
    let use_coord_offsets = elem_size.is_some();
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
        let active_layers: usize = (0..filters.deflate_layers)
            .filter(|&layer| (ch.filter_mask & (1u32 << layer)) == 0)
            .count();
        if active_layers > 0 {
            use std::io::Read as _;
            for _ in 0..active_layers {
                let remaining = MAX_CHUNK_DECOMPRESSED_BYTES.saturating_sub(total_inflated);
                tmp_chunk.clear();
                let dec = flate2::read::ZlibDecoder::new(raw_chunk.as_slice());
                let mut bounded = dec.take(remaining.saturating_add(1));
                let n = bounded.read_to_end(&mut tmp_chunk)?;
                let n_u64 = u64::try_from(n).unwrap_or(u64::MAX);
                total_inflated = total_inflated.saturating_add(n_u64);
                if total_inflated > MAX_CHUNK_DECOMPRESSED_BYTES {
                    return Ok(None);
                }
                std::mem::swap(&mut raw_chunk, &mut tmp_chunk);
            }
        } else {
            let len_u64 = u64::try_from(raw_chunk.len()).unwrap_or(u64::MAX);
            total_inflated = total_inflated.saturating_add(len_u64);
            if total_inflated > MAX_CHUNK_DECOMPRESSED_BYTES {
                return Ok(None);
            }
        }
        if use_coord_offsets {
            let es = elem_size.unwrap_or(1) as u64;
            let Some(byte_off_u64) = ch.offset_elem.checked_mul(es) else {
                return Ok(None);
            };
            if let Some(logical_limit) = logical_bytes_u64
                && !raw_chunk.is_empty()
                && byte_off_u64 >= logical_limit
            {
                return Ok(None);
            }
            let Some(end_u64) = byte_off_u64.checked_add(raw_chunk.len() as u64) else {
                return Ok(None);
            };
            if end_u64 > MAX_CHUNK_DECOMPRESSED_BYTES {
                return Ok(None);
            }
            let Ok(byte_off) = usize::try_from(byte_off_u64) else {
                return Ok(None);
            };
            let Ok(end) = usize::try_from(end_u64) else {
                return Ok(None);
            };
            if byte_off < end {
                let overlaps_prev = written_spans
                    .range(..=byte_off)
                    .next_back()
                    .is_some_and(|(_, &prev_end)| prev_end > byte_off);
                let overlaps_next = written_spans.range(byte_off..end).next().is_some();
                if overlaps_prev || overlaps_next {
                    return Ok(None);
                }
                written_spans.insert(byte_off, end);
            }
            if end > out.len() {
                out.resize(end, 0);
            }
            out[byte_off..end].copy_from_slice(&raw_chunk);
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
