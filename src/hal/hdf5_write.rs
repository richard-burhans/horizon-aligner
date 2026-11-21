//! Pure-Rust HDF5 v0 binary HAL writer.

use super::{HalAlignment, HalBottomSegment, HalSequence, HalTopSegment};

fn encode_hal_nibble(base: u8) -> u8 {
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

/// Encode ASCII DNA into HAL 4-bit nibble-packed `DNA_ARRAY` bytes (high nibble first).
#[must_use]
pub fn encode_hal_dna(dna: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(dna.len().div_ceil(2));
    for chunk in dna.chunks(2) {
        let hi = encode_hal_nibble(chunk[0]);
        let lo = chunk.get(1).copied().map_or(0, encode_hal_nibble);
        out.push((hi << 4) | lo);
    }
    out
}

fn encode_hal_seq_index(
    seqs: &[HalSequence],
    tops: &[HalTopSegment],
    bottoms: &[HalBottomSegment],
    total_dna_len: usize,
) -> Vec<u8> {
    let mut out = Vec::with_capacity((seqs.len() + 1) * 24);
    for s in seqs {
        let top_off = tops.partition_point(|t| t.start < s.offset) as u64;
        let bot_off = bottoms.partition_point(|b| b.start < s.offset) as u64;
        out.extend_from_slice(&(s.offset as u64).to_le_bytes());
        out.extend_from_slice(&top_off.to_le_bytes());
        out.extend_from_slice(&bot_off.to_le_bytes());
    }
    // Sentinel record required by libHal (`SEQIDX_ARRAY` has `num_seqs + 1` entries).
    out.extend_from_slice(&(total_dna_len as u64).to_le_bytes());
    out.extend_from_slice(&(tops.len() as u64).to_le_bytes());
    out.extend_from_slice(&(bottoms.len() as u64).to_le_bytes());
    out
}

fn encode_hal_seq_names(seqs: &[HalSequence]) -> (Vec<u8>, usize) {
    let stride = seqs.iter().map(|s| s.name.len() + 1).max().unwrap_or(1);
    let mut out = vec![0u8; seqs.len() * stride];
    for (i, s) in seqs.iter().enumerate() {
        let bytes = s.name.as_bytes();
        let copy_len = bytes.len().min(stride - 1);
        out[i * stride..i * stride + copy_len].copy_from_slice(&bytes[..copy_len]);
    }
    (out, stride)
}

fn encode_hal_top_array(tops: &[HalTopSegment], total_dna_len: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity((tops.len() + 1) * 33);
    for t in tops {
        let parent = t
            .parent_bottom_idx
            .map_or(-1i64, |v| i64::try_from(v).unwrap_or(-1));
        let paralogy = t
            .next_paralog
            .map_or(-1i64, |v| i64::try_from(v).unwrap_or(-1));
        let start = i64::try_from(t.start).unwrap_or(0);
        out.extend_from_slice(&start.to_le_bytes()); // start coordinate in child genome
        out.extend_from_slice(&(-1i64).to_le_bytes()); // bottom_idx
        out.extend_from_slice(&paralogy.to_le_bytes());
        out.extend_from_slice(&parent.to_le_bytes());
        out.push(u8::from(t.parent_reversed));
    }
    // Sentinel record required by libHal (`TOP_ARRAY` has `num_top + 1` entries; length = next.start - cur.start).
    let end_pos = i64::try_from(total_dna_len).unwrap_or(0);
    out.extend_from_slice(&end_pos.to_le_bytes());
    out.extend_from_slice(&(-1i64).to_le_bytes());
    out.extend_from_slice(&(-1i64).to_le_bytes());
    out.extend_from_slice(&(-1i64).to_le_bytes());
    out.push(0);
    out
}

fn encode_hal_bottom_array(
    bots: &[HalBottomSegment],
    num_children: usize,
    total_dna_len: usize,
) -> Vec<u8> {
    let mut out = Vec::with_capacity((bots.len() + 1) * (24 + 9 * num_children));
    for b in bots {
        let start = i64::try_from(b.start).unwrap_or(0);
        let len = b.length as u64;
        out.extend_from_slice(&start.to_le_bytes());
        out.extend_from_slice(&len.to_le_bytes());
        out.extend_from_slice(&(-1i64).to_le_bytes()); // top_parse
        for c_idx in 0..num_children {
            let link = b.child_tops.get(c_idx).and_then(|opt| *opt);
            let (child_i, rev) = link.map_or((-1i64, 0u8), |(idx, r)| {
                (i64::try_from(idx).unwrap_or(-1), u8::from(r))
            });
            out.extend_from_slice(&child_i.to_le_bytes());
            out.push(rev);
        }
    }
    // Sentinel record required by libHal (`BOTTOM_ARRAY` has `num_bottom + 1` entries).
    let end_pos = i64::try_from(total_dna_len).unwrap_or(0);
    out.extend_from_slice(&end_pos.to_le_bytes());
    out.extend_from_slice(&0u64.to_le_bytes());
    out.extend_from_slice(&(-1i64).to_le_bytes());
    for _ in 0..num_children {
        out.extend_from_slice(&(-1i64).to_le_bytes());
        out.push(0);
    }
    out
}

struct Hdf5V0Builder {
    buf: Vec<u8>,
}

impl Hdf5V0Builder {
    fn new() -> Self {
        Self { buf: vec![0u8; 96] }
    }

    fn align8(&mut self) {
        while !self.buf.len().is_multiple_of(8) {
            self.buf.push(0);
        }
    }

    fn write_v1_obj_header(&mut self, messages: &[(u16, &[u8])]) -> u64 {
        self.align8();
        let hdr_addr = self.buf.len() as u64;
        let mut msg_bytes = Vec::new();
        for &(mtype, body) in messages {
            let padded_len = (body.len() + 7) & !7;
            let padded_u16 = u16::try_from(padded_len).unwrap_or_else(|_| {
                panic!("HDF5 v1 object header message exceeds 65535 bytes: {padded_len}")
            });
            msg_bytes.extend_from_slice(&mtype.to_le_bytes());
            msg_bytes.extend_from_slice(&padded_u16.to_le_bytes());
            msg_bytes.push(0); // flags
            msg_bytes.extend_from_slice(&[0u8; 3]); // reserved
            msg_bytes.extend_from_slice(body);
            msg_bytes.resize(msg_bytes.len() + (padded_len - body.len()), 0);
        }
        let n_msgs_u16 = u16::try_from(messages.len()).unwrap_or_else(|_| {
            panic!(
                "HDF5 v1 object header message count exceeds 65535: {}",
                messages.len()
            )
        });
        let msg_bytes_u32 = u32::try_from(msg_bytes.len()).unwrap_or_else(|_| {
            panic!(
                "HDF5 v1 object header size exceeds u32::MAX: {}",
                msg_bytes.len()
            )
        });
        self.buf.push(1); // version 1
        self.buf.push(0);
        self.buf.extend_from_slice(&n_msgs_u16.to_le_bytes());
        self.buf.extend_from_slice(&1u32.to_le_bytes()); // ref_count
        self.buf.extend_from_slice(&msg_bytes_u32.to_le_bytes());
        self.buf.extend_from_slice(&[0u8; 4]); // padding to 16-byte header
        self.buf.extend_from_slice(&msg_bytes);
        hdr_addr
    }

    fn write_dataset(&mut self, dims: &[u64], elem_size: usize, data: &[u8]) -> u64 {
        self.align8();
        let data_addr = self.buf.len() as u64;
        self.buf.extend_from_slice(data);

        // Dataspace (0x0001), version 1
        let mut ds = vec![1u8, u8::try_from(dims.len()).unwrap_or(1), 0, 0, 0, 0, 0, 0];
        for &d in dims {
            ds.extend_from_slice(&d.to_le_bytes());
        }

        // Datatype (0x0003), Opaque class (5), version 1, element size at [4..8]
        let mut dt = vec![0u8; 8];
        dt[0] = 0x15;
        dt[4..8].copy_from_slice(&u32::try_from(elem_size).unwrap_or(1).to_le_bytes());

        // DataLayout (0x0008), version 3, class 1 (contiguous)
        let layout_addr = if data.is_empty() { u64::MAX } else { data_addr };
        let mut dl = Vec::with_capacity(18);
        dl.push(3);
        dl.push(1);
        dl.extend_from_slice(&layout_addr.to_le_bytes());
        dl.extend_from_slice(&(data.len() as u64).to_le_bytes());

        self.write_v1_obj_header(&[(0x0001, &ds), (0x0003, &dt), (0x0008, &dl)])
    }

    fn write_group_btree(&mut self, snods: Vec<(u64, u64, u64)>) -> u64 {
        // Group Internal Node K = 16 => max 32 child entries per TREE node.
        const MAX_TREE_CHILDREN: usize = 32;
        if snods.is_empty() {
            self.align8();
            let empty_addr = self.buf.len() as u64;
            self.buf.extend_from_slice(b"TREE");
            self.buf.push(0);
            self.buf.push(0);
            self.buf.extend_from_slice(&0u16.to_le_bytes());
            self.buf.extend_from_slice(&u64::MAX.to_le_bytes());
            self.buf.extend_from_slice(&u64::MAX.to_le_bytes());
            self.buf.extend_from_slice(&0u64.to_le_bytes());
            for _ in 0..MAX_TREE_CHILDREN {
                self.buf.extend_from_slice(&0u64.to_le_bytes());
                self.buf.extend_from_slice(&0u64.to_le_bytes());
            }
            return empty_addr;
        }

        let mut level_nodes: Vec<(u64, u64, u64)> = snods;
        let mut level: u8 = 0;

        loop {
            let num_chunks = level_nodes.len().div_ceil(MAX_TREE_CHILDREN);
            let mut next_level: Vec<(u64, u64, u64)> = Vec::with_capacity(num_chunks);
            let mut prev_last_off = 0u64;
            for chunk in level_nodes.chunks(MAX_TREE_CHILDREN) {
                self.align8();
                let node_addr = self.buf.len() as u64;
                self.buf.extend_from_slice(b"TREE");
                self.buf.push(0); // type 0 (group)
                self.buf.push(level);
                self.buf.extend_from_slice(
                    &u16::try_from(chunk.len()).unwrap_or(u16::MAX).to_le_bytes(),
                );
                self.buf.extend_from_slice(&u64::MAX.to_le_bytes());
                self.buf.extend_from_slice(&u64::MAX.to_le_bytes());
                self.buf.extend_from_slice(&prev_last_off.to_le_bytes());
                for &(_, child_addr, last_off) in chunk {
                    self.buf.extend_from_slice(&child_addr.to_le_bytes());
                    self.buf.extend_from_slice(&last_off.to_le_bytes());
                }
                let first_key = chunk.first().map_or(0, |&(first_off, _, _)| first_off);
                let end_key = chunk.last().map_or(0, |&(_, _, last_off)| last_off);
                for _ in chunk.len()..MAX_TREE_CHILDREN {
                    self.buf.extend_from_slice(&0u64.to_le_bytes());
                    self.buf.extend_from_slice(&0u64.to_le_bytes());
                }
                next_level.push((first_key, node_addr, end_key));
                prev_last_off = end_key;
            }
            for idx in 0..next_level.len() {
                let addr = usize::try_from(next_level[idx].1).unwrap_or(0);
                if idx > 0 {
                    let left = next_level[idx - 1].1;
                    self.buf[addr + 8..addr + 16].copy_from_slice(&left.to_le_bytes());
                }
                if idx + 1 < next_level.len() {
                    let right = next_level[idx + 1].1;
                    self.buf[addr + 16..addr + 24].copy_from_slice(&right.to_le_bytes());
                }
            }
            if next_level.len() == 1 {
                return next_level[0].1;
            }
            level_nodes = next_level;
            level = level.saturating_add(1);
        }
    }

    fn write_gcol(&mut self, strings: &[&str]) -> u64 {
        self.align8();
        let gcol_addr = self.buf.len() as u64;
        let mut body = Vec::new();
        for (i, s) in strings.iter().enumerate() {
            let idx = u16::try_from(i + 1).unwrap_or(u16::MAX);
            let bytes = s.as_bytes();
            let padded = (bytes.len() + 7) & !7;
            body.extend_from_slice(&idx.to_le_bytes());
            body.extend_from_slice(&0u16.to_le_bytes()); // ref_count
            body.extend_from_slice(&[0u8; 4]); // reserved
            body.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
            body.extend_from_slice(bytes);
            body.resize(body.len() + (padded - bytes.len()), 0);
        }
        // Minimum GCOL size is 4096 bytes (including 16-byte header + 16-byte free sentinel).
        let total_size = (16 + body.len() + 16).max(4096);
        let total_size = (total_size + 7) & !7;
        let free_space = total_size - 16 - body.len();

        self.buf.extend_from_slice(b"GCOL");
        self.buf.push(1); // version 1
        self.buf.extend_from_slice(&[0u8; 3]); // reserved
        self.buf.extend_from_slice(&(total_size as u64).to_le_bytes());
        self.buf.extend_from_slice(&body);
        // End-of-collection free object (index = 0)
        self.buf.extend_from_slice(&0u16.to_le_bytes());
        self.buf.extend_from_slice(&0u16.to_le_bytes());
        self.buf.extend_from_slice(&[0u8; 4]);
        self.buf.extend_from_slice(&(free_space as u64).to_le_bytes());
        self.buf.resize(self.buf.len() + free_space - 16, 0);
        gcol_addr
    }

    fn write_group_with_attrs(
        &mut self,
        entries: &[(&str, u64)],
        attrs: &[(&str, u32, u64, u32)],
    ) -> (u64, u64, u64) {
        // Group Leaf Node K = 4 => max 8 entries per SNOD
        const MAX_SNOD: usize = 8;
        self.align8();
        let mut sorted: Vec<(&str, u64)> = entries.to_vec();
        sorted.sort_by_key(|&(name, _)| name);

        let mut heap_data = vec![0u8]; // offset 0 reserved for empty string
        let mut name_offsets = Vec::with_capacity(sorted.len());
        for &(name, _) in &sorted {
            name_offsets.push(heap_data.len() as u64);
            heap_data.extend_from_slice(name.as_bytes());
            heap_data.push(0);
        }
        while !heap_data.len().is_multiple_of(8) {
            heap_data.push(0);
        }

        let heap_data_addr = self.buf.len() as u64;
        self.buf.extend_from_slice(&heap_data);

        self.align8();
        let heap_addr = self.buf.len() as u64;
        self.buf.extend_from_slice(b"HEAP");
        self.buf.extend_from_slice(&[0u8; 4]); // version 0 + reserved
        self.buf
            .extend_from_slice(&(heap_data.len() as u64).to_le_bytes());
        self.buf.extend_from_slice(&1u64.to_le_bytes()); // H5HL_FREE_NULL (1): no free blocks
        self.buf.extend_from_slice(&heap_data_addr.to_le_bytes());

        let mut snods: Vec<(u64, u64, u64)> = Vec::new();
        for (chunk_idx, chunk) in sorted.chunks(MAX_SNOD).enumerate() {
            self.align8();
            let snod_addr = self.buf.len() as u64;
            self.buf.extend_from_slice(b"SNOD");
            self.buf.push(1); // version 1
            self.buf.push(0);
            self.buf
                .extend_from_slice(&u16::try_from(chunk.len()).unwrap_or(u16::MAX).to_le_bytes());
            for (local_idx, &(_, obj_addr)) in chunk.iter().enumerate() {
                let idx = chunk_idx * MAX_SNOD + local_idx;
                let mut rec = [0u8; 40];
                rec[0..8].copy_from_slice(&name_offsets[idx].to_le_bytes());
                rec[8..16].copy_from_slice(&obj_addr.to_le_bytes());
                self.buf.extend_from_slice(&rec);
            }
            // Pad SNOD to 2 * leaf_k = 8 slots so libhdf5 reads a full node cleanly.
            for _ in chunk.len()..MAX_SNOD {
                self.buf.extend_from_slice(&[0u8; 40]);
            }
            let base_idx = chunk_idx * MAX_SNOD;
            let first_off = name_offsets[base_idx];
            let last_off = name_offsets[base_idx + chunk.len() - 1];
            snods.push((first_off, snod_addr, last_off));
        }

        let btree_addr = self.write_group_btree(snods);

        let mut symtab = [0u8; 16];
        symtab[0..8].copy_from_slice(&btree_addr.to_le_bytes());
        symtab[8..16].copy_from_slice(&heap_addr.to_le_bytes());

        let attr_msgs: Vec<Vec<u8>> = attrs
            .iter()
            .map(|&(k, len, gcol_addr, obj_idx)| encode_vls_attr_msg(k, len, gcol_addr, obj_idx))
            .collect();
        let mut msgs: Vec<(u16, &[u8])> = Vec::with_capacity(1 + attr_msgs.len());
        msgs.push((0x0011, &symtab));
        for a in &attr_msgs {
            msgs.push((0x000C, a.as_slice()));
        }
        let hdr_addr = self.write_v1_obj_header(&msgs);
        (hdr_addr, btree_addr, heap_addr)
    }

    fn write_group(&mut self, entries: &[(&str, u64)]) -> (u64, u64, u64) {
        self.write_group_with_attrs(entries, &[])
    }

    fn finish(mut self, root_hdr: u64, root_bt: u64, root_hp: u64) -> Vec<u8> {
        let eof_addr = self.buf.len() as u64;
        self.buf[0..8].copy_from_slice(b"\x89HDF\r\n\x1a\n");
        self.buf[8] = 0; // superblock version 0
        self.buf[13] = 8; // offset_size
        self.buf[14] = 8; // length_size
        self.buf[16..18].copy_from_slice(&4u16.to_le_bytes()); // group leaf node K
        self.buf[18..20].copy_from_slice(&16u16.to_le_bytes()); // group internal node K
        self.buf[24..32].copy_from_slice(&0u64.to_le_bytes());
        self.buf[32..40].copy_from_slice(&u64::MAX.to_le_bytes());
        self.buf[40..48].copy_from_slice(&eof_addr.to_le_bytes());
        self.buf[48..56].copy_from_slice(&u64::MAX.to_le_bytes());
        // Root group symbol-table entry at [56..96]
        self.buf[56..64].copy_from_slice(&0u64.to_le_bytes());
        self.buf[64..72].copy_from_slice(&root_hdr.to_le_bytes());
        self.buf[72..76].copy_from_slice(&1u32.to_le_bytes());
        self.buf[80..88].copy_from_slice(&root_bt.to_le_bytes());
        self.buf[88..96].copy_from_slice(&root_hp.to_le_bytes());
        self.buf
    }
}

fn encode_vls_attr_msg(name: &str, str_len: u32, gcol_addr: u64, obj_idx: u32) -> Vec<u8> {
    let name_bytes = name.as_bytes();
    let name_size = u16::try_from(name_bytes.len() + 1).unwrap_or_else(|_| {
        panic!("HDF5 attribute name exceeds 65535 bytes: {}", name.len())
    });

    let mut msg = Vec::new();
    msg.push(1); // Attribute message version 1
    msg.push(0); // reserved
    msg.extend_from_slice(&name_size.to_le_bytes());
    msg.extend_from_slice(&20u16.to_le_bytes()); // datatype size (20 bytes for VLS of char)
    msg.extend_from_slice(&8u16.to_le_bytes()); // dataspace size (8 bytes for scalar)

    // Name (null-terminated, padded to 8-byte boundary)
    msg.extend_from_slice(name_bytes);
    msg.push(0);
    while !msg.len().is_multiple_of(8) {
        msg.push(0);
    }

    // Datatype: VarLen class (9), version 1, string type (1), null-terminated ASCII (20 bytes, padded to 24)
    msg.extend_from_slice(&[0x19, 0x01, 0x00, 0x00]);
    msg.extend_from_slice(&16u32.to_le_bytes());
    // Parent base type: Integer class (0), version 1, size 1, precision 8 bits (12 bytes)
    msg.extend_from_slice(&[
        0x10, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08, 0x00,
    ]);
    // Pad datatype to 8-byte boundary (20 -> 24 bytes)
    msg.extend_from_slice(&[0u8; 4]);

    // Dataspace: Scalar (version 1, rank 0)
    msg.extend_from_slice(&[1, 0, 0, 0, 0, 0, 0, 0]);

    // Data: hvl_t on disk = (len: u32, gcol_addr: u64, obj_idx: u32)
    msg.extend_from_slice(&str_len.to_le_bytes());
    msg.extend_from_slice(&gcol_addr.to_le_bytes());
    msg.extend_from_slice(&obj_idx.to_le_bytes());
    msg
}

fn interval_overlaps(
    occupied: &std::collections::BTreeMap<usize, usize>,
    start: usize,
    end: usize,
) -> bool {
    if let Some((&_prev_s, &prev_e)) = occupied.range(..end).next_back()
        && prev_e > start
    {
        return true;
    }
    occupied.range(start..end).next().is_some()
}

fn seq_contains_interval(seqs: &[HalSequence], start: usize, end: usize) -> bool {
    let idx = seqs.partition_point(|s| s.offset.saturating_add(s.length) <= start);
    seqs.get(idx)
        .is_some_and(|s| start >= s.offset && end <= s.offset.saturating_add(s.length))
}

fn finalize_child_top_segments(
    sequences: &[HalSequence],
    mut segs: Vec<(usize, usize, usize, bool)>,
    c_idx: usize,
    root_bottoms: &mut [HalBottomSegment],
) -> Vec<HalTopSegment> {
    use std::collections::HashMap;
    segs.sort_by_key(|&(s, _, _, _)| s);

    let mut final_tops: Vec<HalTopSegment> = Vec::new();
    let mut seg_idx = 0usize;
    for seq in sequences {
        let mut pos = seq.offset;
        let seq_end = seq.offset + seq.length;
        while seg_idx < segs.len() && segs[seg_idx].0 < seq_end {
            let (s, l, p_b_idx, rev) = segs[seg_idx];
            seg_idx += 1;
            if s < pos {
                continue;
            }
            if s > pos {
                final_tops.push(HalTopSegment {
                    start: pos,
                    length: s - pos,
                    parent_bottom_idx: None,
                    parent_reversed: false,
                    next_paralog: None,
                });
            }
            final_tops.push(HalTopSegment {
                start: s,
                length: l,
                parent_bottom_idx: Some(p_b_idx),
                parent_reversed: rev,
                next_paralog: None,
            });
            pos = s + l;
        }
        if pos < seq_end {
            final_tops.push(HalTopSegment {
                start: pos,
                length: seq_end - pos,
                parent_bottom_idx: None,
                parent_reversed: false,
                next_paralog: None,
            });
        }
    }

    let mut by_parent: HashMap<usize, Vec<(usize, bool)>> = HashMap::new();
    for (t_idx, t) in final_tops.iter().enumerate() {
        if let Some(pb) = t.parent_bottom_idx {
            by_parent
                .entry(pb)
                .or_default()
                .push((t_idx, t.parent_reversed));
        }
    }
    for (pb, list) in by_parent {
        if let Some(&(first_t, first_rev)) = list.first() {
            root_bottoms[pb].child_tops[c_idx] = Some((first_t, first_rev));
        }
        if list.len() > 1 {
            for w in 0..list.len() {
                let cur_t = list[w].0;
                let next_t = list[(w + 1) % list.len()].0;
                final_tops[cur_t].next_paralog = Some(next_t);
            }
        }
    }

    final_tops
}

struct NormalizedSegments {
    tops_by_genome: Vec<Vec<HalTopSegment>>,
    bottoms_by_genome: Vec<Vec<HalBottomSegment>>,
}

#[allow(clippy::too_many_lines)]
fn normalize_hal_segments_for_export(hal: &HalAlignment) -> NormalizedSegments {
    use std::collections::{BTreeMap, BTreeSet};
    let n_genomes = hal.genomes.len();
    if n_genomes == 0 {
        return NormalizedSegments {
            tops_by_genome: Vec::new(),
            bottoms_by_genome: Vec::new(),
        };
    }
    let root = &hal.genomes[0];
    if root.parent.is_some()
        || hal.genomes[1..]
            .iter()
            .any(|g| g.parent != Some(0) || !g.children.is_empty())
    {
        // Non-star tree: sort segments per genome and remap all cross-genome indices.
        let mut top_orders: Vec<Vec<usize>> = Vec::with_capacity(n_genomes);
        let mut top_invs: Vec<Vec<usize>> = Vec::with_capacity(n_genomes);
        let mut bot_orders: Vec<Vec<usize>> = Vec::with_capacity(n_genomes);
        let mut bot_invs: Vec<Vec<usize>> = Vec::with_capacity(n_genomes);
        for g in &hal.genomes {
            let mut t_ord: Vec<usize> = (0..g.top_segments.len()).collect();
            t_ord.sort_by_key(|&i| g.top_segments[i].start);
            let mut t_inv = vec![0usize; g.top_segments.len()];
            for (new_i, &old_i) in t_ord.iter().enumerate() {
                t_inv[old_i] = new_i;
            }
            top_orders.push(t_ord);
            top_invs.push(t_inv);

            let mut b_ord: Vec<usize> = (0..g.bottom_segments.len()).collect();
            b_ord.sort_by_key(|&i| g.bottom_segments[i].start);
            let mut b_inv = vec![0usize; g.bottom_segments.len()];
            for (new_i, &old_i) in b_ord.iter().enumerate() {
                b_inv[old_i] = new_i;
            }
            bot_orders.push(b_ord);
            bot_invs.push(b_inv);
        }

        let tops_by_genome: Vec<Vec<HalTopSegment>> = hal
            .genomes
            .iter()
            .enumerate()
            .map(|(gi, g)| {
                top_orders[gi]
                    .iter()
                    .map(|&old_i| {
                        let mut seg = g.top_segments[old_i].clone();
                        seg.parent_bottom_idx = seg.parent_bottom_idx.map(|pb| {
                            g.parent
                                .and_then(|pg| bot_invs.get(pg))
                                .and_then(|inv| inv.get(pb))
                                .copied()
                                .unwrap_or(pb)
                        });
                        seg.next_paralog = seg
                            .next_paralog
                            .and_then(|np| top_invs[gi].get(np).copied());
                        seg
                    })
                    .collect()
            })
            .collect();
        let bottoms_by_genome: Vec<Vec<HalBottomSegment>> = hal
            .genomes
            .iter()
            .enumerate()
            .map(|(gi, g)| {
                bot_orders[gi]
                    .iter()
                    .map(|&old_i| {
                        let mut seg = g.bottom_segments[old_i].clone();
                        for (c_idx, link) in seg.child_tops.iter_mut().enumerate() {
                            if let Some((old_t, rev)) = *link {
                                let new_t = g
                                    .children
                                    .get(c_idx)
                                    .and_then(|&cg| top_invs.get(cg))
                                    .and_then(|inv| inv.get(old_t))
                                    .copied()
                                    .unwrap_or(old_t);
                                *link = Some((new_t, rev));
                            }
                        }
                        seg
                    })
                    .collect()
            })
            .collect();
        return NormalizedSegments {
            tops_by_genome,
            bottoms_by_genome,
        };
    }

    let num_children = root.children.len();
    let mut pts = BTreeSet::new();
    pts.insert(0usize);
    pts.insert(root.dna.len());
    for s in &root.sequences {
        pts.insert(s.offset.min(root.dna.len()));
        pts.insert((s.offset + s.length).min(root.dna.len()));
    }
    let mut sorted_orig_b: Vec<&HalBottomSegment> = root
        .bottom_segments
        .iter()
        .filter(|b| b.length > 0 && b.start < root.dna.len())
        .collect();
    for &b in &sorted_orig_b {
        pts.insert(b.start);
        pts.insert((b.start + b.length).min(root.dna.len()));
    }
    sorted_orig_b.sort_by_key(|b| b.start);
    let pts_vec: Vec<usize> = pts.into_iter().collect();

    let mut new_bottoms: Vec<HalBottomSegment> = Vec::with_capacity(pts_vec.len());
    let mut placed_tops: Vec<Vec<(usize, usize, usize, bool)>> = vec![Vec::new(); num_children];
    let mut occupied_child: Vec<BTreeMap<usize, usize>> = vec![BTreeMap::new(); num_children];
    let mut active: Vec<&HalBottomSegment> = Vec::new();
    let mut next_b = 0usize;

    for win in pts_vec.windows(2) {
        let a = win[0];
        let b = win[1];
        if b <= a {
            continue;
        }
        while next_b < sorted_orig_b.len() && sorted_orig_b[next_b].start <= a {
            active.push(sorted_orig_b[next_b]);
            next_b += 1;
        }
        active.retain(|orig_b| orig_b.start + orig_b.length >= b);

        let sub_len = b - a;
        let b_idx = new_bottoms.len();
        let child_links = vec![None; num_children];

        for &orig_b in &active {
            for c_idx in 0..num_children {
                let Some(child_g_idx) = root.children.get(c_idx).copied() else {
                    continue;
                };
                let Some(child_g) = hal.genomes.get(child_g_idx) else {
                    continue;
                };
                let Some((mut cur_t_idx, mut cur_rev)) =
                    orig_b.child_tops.get(c_idx).and_then(|x| *x)
                else {
                    continue;
                };
                let mut seen_t = 0usize;
                while let Some(orig_t) = child_g.top_segments.get(cur_t_idx) {
                    let off_in_orig = a - orig_b.start;
                    let t_start = if cur_rev {
                        let tail_off = orig_b.length.saturating_sub(off_in_orig + sub_len);
                        orig_t.start + tail_off
                    } else {
                        orig_t.start + off_in_orig
                    };
                    let t_end = t_start + sub_len;
                    if seq_contains_interval(&child_g.sequences, t_start, t_end)
                        && !interval_overlaps(&occupied_child[c_idx], t_start, t_end)
                    {
                        occupied_child[c_idx].insert(t_start, t_end);
                        placed_tops[c_idx].push((t_start, sub_len, b_idx, cur_rev));
                    }
                    seen_t += 1;
                    let Some(next_t) = orig_t.next_paralog else {
                        break;
                    };
                    if next_t == cur_t_idx || seen_t >= child_g.top_segments.len() {
                        break;
                    }
                    cur_t_idx = next_t;
                    cur_rev = child_g
                        .top_segments
                        .get(cur_t_idx)
                        .is_some_and(|t| t.parent_reversed);
                }
            }
        }

        new_bottoms.push(HalBottomSegment {
            start: a,
            length: sub_len,
            child_tops: child_links,
        });
    }

    let mut tops_by_genome = vec![Vec::new(); n_genomes];
    let mut bottoms_by_genome = vec![Vec::new(); n_genomes];
    for (c_idx, &child_g_idx) in root.children.iter().enumerate() {
        let Some(child_g) = hal.genomes.get(child_g_idx) else {
            continue;
        };
        if let Some(dst) = tops_by_genome.get_mut(child_g_idx) {
            *dst = finalize_child_top_segments(
                &child_g.sequences,
                std::mem::take(&mut placed_tops[c_idx]),
                c_idx,
                &mut new_bottoms,
            );
        }
    }
    bottoms_by_genome[0] = new_bottoms;

    NormalizedSegments {
        tops_by_genome,
        bottoms_by_genome,
    }
}

impl HalAlignment {
    /// Serialize this `HalAlignment` into a valid, self-contained HDF5 `v0` `.hal` byte stream
    /// using pure Rust without linking `libhdf5`.
    #[must_use]
    pub fn to_hdf5_bytes(&self) -> Vec<u8> {
        let norm = normalize_hal_segments_for_export(self);
        let mut w = Hdf5V0Builder::new();
        // Store VLS attribute strings in a single Global Heap collection:
        // idx 1 = tree_newick, idx 2 = "2.2", idx 3 = "0", idx 4 = "1"
        let gcol_addr = w.write_gcol(&[self.tree_newick.as_str(), "2.2", "0", "1"]);
        let tree_len = u32::try_from(self.tree_newick.len()).unwrap_or(0);

        let (meta_hdr, _, _) = w.write_group(&[]);
        let (genomes_hdr, _, _) = w.write_group(&[]);
        let (phylo_hdr, _, _) =
            w.write_group_with_attrs(&[], &[("Phylogeny", tree_len, gcol_addr, 1)]);
        let (ver_hdr, _, _) = w.write_group_with_attrs(&[], &[("Verison", 3, gcol_addr, 2)]);

        let mut root_entries: Vec<(&str, u64)> = Vec::with_capacity(self.genomes.len() + 4);
        root_entries.push(("Meta", meta_hdr));
        root_entries.push(("Phylogeny", phylo_hdr));
        root_entries.push(("Genomes", genomes_hdr));
        root_entries.push(("Verison", ver_hdr));

        for (g_idx, g) in self.genomes.iter().enumerate() {
            let tops = &norm.tops_by_genome[g_idx];
            let bottoms = &norm.bottoms_by_genome[g_idx];

            let dna_packed = encode_hal_dna(&g.dna);
            let dna_addr = w.write_dataset(&[dna_packed.len() as u64], 1, &dna_packed);

            let seqidx_bytes =
                encode_hal_seq_index(&g.sequences, tops, bottoms, g.dna.len());
            let seqidx_addr =
                w.write_dataset(&[(g.sequences.len() + 1) as u64], 24, &seqidx_bytes);

            let (seqname_bytes, stride) = encode_hal_seq_names(&g.sequences);
            let seqname_addr = w.write_dataset(&[g.sequences.len() as u64], stride, &seqname_bytes);

            let top_bytes = encode_hal_top_array(tops, g.dna.len());
            let top_addr = w.write_dataset(&[(tops.len() + 1) as u64], 33, &top_bytes);

            let num_children = g.children.len();
            let bot_bytes = encode_hal_bottom_array(bottoms, num_children, g.dna.len());
            let bot_addr = w.write_dataset(
                &[(bottoms.len() + 1) as u64],
                24 + 9 * num_children,
                &bot_bytes,
            );

            let (g_meta_hdr, _, _) = w.write_group(&[]);
            let rup_idx = if g.dna.len() % 2 == 1 { 4 } else { 3 };
            let (g_rup_hdr, _, _) = w.write_group_with_attrs(&[], &[("Rup", 1, gcol_addr, rup_idx)]);

            let (g_hdr, _, _) = w.write_group(&[
                ("DNA_ARRAY", dna_addr),
                ("SEQIDX_ARRAY", seqidx_addr),
                ("SEQNAME_ARRAY", seqname_addr),
                ("TOP_ARRAY", top_addr),
                ("BOTTOM_ARRAY", bot_addr),
                ("Meta", g_meta_hdr),
                ("Rup", g_rup_hdr),
            ]);
            root_entries.push((g.name.as_str(), g_hdr));
        }
        let (root_hdr, root_bt, root_hp) = w.write_group(&root_entries);
        w.finish(root_hdr, root_bt, root_hp)
    }

    /// Write this `HalAlignment` to disk as an HDF5 `v0` `.hal` file.
    ///
    /// # Errors
    /// Returns an `std::io::Error` if writing the file fails.
    pub fn write_hdf5_hal(&self, path: &std::path::Path) -> std::io::Result<()> {
        std::fs::write(path, self.to_hdf5_bytes())
    }
}
