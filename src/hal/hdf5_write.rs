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

fn encode_hal_seq_index(seqs: &[HalSequence]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seqs.len() * 24);
    for s in seqs {
        out.extend_from_slice(&(s.offset as u64).to_le_bytes());
        out.extend_from_slice(&0u64.to_le_bytes());
        out.extend_from_slice(&0u64.to_le_bytes());
    }
    out
}

fn encode_hal_seq_names(seqs: &[HalSequence]) -> (Vec<u8>, usize) {
    let stride = seqs
        .iter()
        .map(|s| s.name.len() + 1)
        .max()
        .unwrap_or(1)
        .max(1);
    let mut out = vec![0u8; seqs.len() * stride];
    for (i, s) in seqs.iter().enumerate() {
        let bytes = s.name.as_bytes();
        let copy_len = bytes.len().min(stride - 1);
        out[i * stride..i * stride + copy_len].copy_from_slice(&bytes[..copy_len]);
    }
    (out, stride)
}

fn encode_hal_top_array(tops: &[HalTopSegment]) -> Vec<u8> {
    let mut out = Vec::with_capacity(tops.len() * 33);
    for t in tops {
        let parent = t
            .parent_bottom_idx
            .map_or(-1i64, |v| i64::try_from(v).unwrap_or(-1));
        let paralogy = t
            .next_paralog
            .map_or(-1i64, |v| i64::try_from(v).unwrap_or(-1));
        out.extend_from_slice(&0i64.to_le_bytes()); // genome_idx
        out.extend_from_slice(&(-1i64).to_le_bytes()); // bottom_idx
        out.extend_from_slice(&paralogy.to_le_bytes());
        out.extend_from_slice(&parent.to_le_bytes());
        out.push(u8::from(t.parent_reversed));
    }
    out
}

fn encode_hal_bottom_array(bots: &[HalBottomSegment], num_children: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(bots.len() * (24 + 9 * num_children));
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
            msg_bytes.extend_from_slice(&mtype.to_le_bytes());
            msg_bytes
                .extend_from_slice(&u16::try_from(padded_len).unwrap_or(u16::MAX).to_le_bytes());
            msg_bytes.push(0); // flags
            msg_bytes.extend_from_slice(&[0u8; 3]); // reserved
            msg_bytes.extend_from_slice(body);
            msg_bytes.resize(msg_bytes.len() + (padded_len - body.len()), 0);
        }
        self.buf.push(1); // version 1
        self.buf.push(0);
        self.buf.extend_from_slice(
            &u16::try_from(messages.len())
                .unwrap_or(u16::MAX)
                .to_le_bytes(),
        );
        self.buf.extend_from_slice(&1u32.to_le_bytes()); // ref_count
        self.buf.extend_from_slice(
            &u32::try_from(msg_bytes.len())
                .unwrap_or(u32::MAX)
                .to_le_bytes(),
        );
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

        // Datatype (0x0003), minimal 8-byte header with element size at [4..8]
        let mut dt = vec![0u8; 8];
        dt[4..8].copy_from_slice(&u32::try_from(elem_size).unwrap_or(1).to_le_bytes());

        // DataLayout (0x0008), version 3, class 1 (contiguous)
        let mut dl = Vec::with_capacity(18);
        dl.push(3);
        dl.push(1);
        dl.extend_from_slice(&data_addr.to_le_bytes());
        dl.extend_from_slice(&(data.len() as u64).to_le_bytes());

        self.write_v1_obj_header(&[(0x0001, &ds), (0x0003, &dt), (0x0008, &dl)])
    }

    fn write_group(&mut self, entries: &[(&str, u64)]) -> (u64, u64, u64) {
        self.align8();
        let mut heap_data = vec![0u8]; // offset 0 reserved for empty string
        let mut name_offsets = Vec::with_capacity(entries.len());
        for &(name, _) in entries {
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
        self.buf.extend_from_slice(&0u64.to_le_bytes()); // free offset
        self.buf.extend_from_slice(&heap_data_addr.to_le_bytes());

        self.align8();
        let snod_addr = self.buf.len() as u64;
        self.buf.extend_from_slice(b"SNOD");
        self.buf.push(1); // version 1
        self.buf.push(0);
        self.buf.extend_from_slice(
            &u16::try_from(entries.len())
                .unwrap_or(u16::MAX)
                .to_le_bytes(),
        );
        for (idx, &(_, obj_addr)) in entries.iter().enumerate() {
            let mut rec = [0u8; 40];
            rec[0..8].copy_from_slice(&name_offsets[idx].to_le_bytes());
            rec[8..16].copy_from_slice(&obj_addr.to_le_bytes());
            self.buf.extend_from_slice(&rec);
        }

        self.align8();
        let btree_addr = self.buf.len() as u64;
        self.buf.extend_from_slice(b"TREE");
        self.buf.push(0); // type 0 (group)
        self.buf.push(0); // level 0 (leaf)
        self.buf.extend_from_slice(&1u16.to_le_bytes()); // n_entries
        self.buf.extend_from_slice(&u64::MAX.to_le_bytes()); // left sibling
        self.buf.extend_from_slice(&u64::MAX.to_le_bytes()); // right sibling
        self.buf.extend_from_slice(&0u64.to_le_bytes()); // key 0
        self.buf.extend_from_slice(&snod_addr.to_le_bytes()); // child SNOD
        self.buf.extend_from_slice(&u64::MAX.to_le_bytes()); // end key

        let mut symtab = [0u8; 16];
        symtab[0..8].copy_from_slice(&btree_addr.to_le_bytes());
        symtab[8..16].copy_from_slice(&heap_addr.to_le_bytes());
        let hdr_addr = self.write_v1_obj_header(&[(0x0011, &symtab)]);
        (hdr_addr, btree_addr, heap_addr)
    }

    fn finish(mut self, root_hdr: u64, root_bt: u64, root_hp: u64) -> Vec<u8> {
        let eof_addr = self.buf.len() as u64;
        self.buf[0..8].copy_from_slice(b"\x89HDF\r\n\x1a\n");
        self.buf[8] = 0; // superblock version 0
        self.buf[13] = 8; // offset_size
        self.buf[14] = 8; // length_size
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

impl HalAlignment {
    /// Serialize this `HalAlignment` into a valid, self-contained HDF5 `v0` `.hal` byte stream
    /// using pure Rust without linking `libhdf5`.
    #[must_use]
    pub fn to_hdf5_bytes(&self) -> Vec<u8> {
        let mut w = Hdf5V0Builder::new();
        let mut root_entries: Vec<(&str, u64)> = Vec::with_capacity(self.genomes.len());
        for g in &self.genomes {
            let dna_packed = encode_hal_dna(&g.dna);
            let dna_addr = w.write_dataset(&[dna_packed.len() as u64], 1, &dna_packed);

            let seqidx_bytes = encode_hal_seq_index(&g.sequences);
            let seqidx_addr = w.write_dataset(&[g.sequences.len() as u64], 24, &seqidx_bytes);

            let (seqname_bytes, stride) = encode_hal_seq_names(&g.sequences);
            let seqname_addr = w.write_dataset(&[g.sequences.len() as u64], stride, &seqname_bytes);

            let top_bytes = encode_hal_top_array(&g.top_segments);
            let top_addr = w.write_dataset(&[g.top_segments.len() as u64], 33, &top_bytes);

            let num_children = g.children.len();
            let bot_bytes = encode_hal_bottom_array(&g.bottom_segments, num_children);
            let bot_addr = w.write_dataset(
                &[g.bottom_segments.len() as u64],
                24 + 9 * num_children,
                &bot_bytes,
            );

            let (g_hdr, _, _) = w.write_group(&[
                ("DNA_ARRAY", dna_addr),
                ("SEQIDX_ARRAY", seqidx_addr),
                ("SEQNAME_ARRAY", seqname_addr),
                ("TOP_ARRAY", top_addr),
                ("BOTTOM_ARRAY", bot_addr),
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
