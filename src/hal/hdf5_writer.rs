//! Low-level HDF5 v0 container writer (superblock, B-tree, HEAP, GCOL, contiguous datasets).

use std::io::{Seek, SeekFrom, Write};

pub(super) struct Hdf5V0Builder<W: Write + Seek> {
    w: W,
    pos: u64,
}

impl<W: Write + Seek> Hdf5V0Builder<W> {
    pub(super) fn new(mut w: W) -> std::io::Result<Self> {
        w.write_all(&[0u8; 96])?;
        Ok(Self { w, pos: 96 })
    }

    fn write_bytes(&mut self, data: &[u8]) -> std::io::Result<()> {
        self.w.write_all(data)?;
        self.pos += data.len() as u64;
        Ok(())
    }

    fn align8(&mut self) -> std::io::Result<()> {
        let rem = usize::try_from(self.pos & 7).unwrap_or(0);
        if rem != 0 {
            const ZEROS: [u8; 8] = [0u8; 8];
            self.write_bytes(&ZEROS[..8 - rem])?;
        }
        Ok(())
    }

    fn write_v1_obj_header(&mut self, messages: &[(u16, &[u8])]) -> std::io::Result<u64> {
        self.align8()?;
        let hdr_addr = self.pos;
        let mut msg_bytes = Vec::new();
        for &(mtype, body) in messages {
            let padded_len = (body.len() + 7) & !7;
            let padded_u16 = u16::try_from(padded_len).map_err(|_| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    format!("HDF5 v1 object header message exceeds 65535 bytes: {padded_len}"),
                )
            })?;
            msg_bytes.extend_from_slice(&mtype.to_le_bytes());
            msg_bytes.extend_from_slice(&padded_u16.to_le_bytes());
            msg_bytes.push(0); // flags
            msg_bytes.extend_from_slice(&[0u8; 3]); // reserved
            msg_bytes.extend_from_slice(body);
            msg_bytes.resize(msg_bytes.len() + (padded_len - body.len()), 0);
        }
        let n_msgs_u16 = u16::try_from(messages.len()).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "HDF5 v1 object header message count exceeds 65535: {}",
                    messages.len()
                ),
            )
        })?;
        let msg_bytes_u32 = u32::try_from(msg_bytes.len()).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "HDF5 v1 object header size exceeds u32::MAX: {}",
                    msg_bytes.len()
                ),
            )
        })?;
        let mut hdr = [0u8; 16];
        hdr[0] = 1; // version 1
        hdr[2..4].copy_from_slice(&n_msgs_u16.to_le_bytes());
        hdr[4..8].copy_from_slice(&1u32.to_le_bytes()); // ref_count
        hdr[8..12].copy_from_slice(&msg_bytes_u32.to_le_bytes());
        self.write_bytes(&hdr)?;
        self.write_bytes(&msg_bytes)?;
        Ok(hdr_addr)
    }

    fn write_dataset_header(
        &mut self,
        dims: &[u64],
        elem_size: usize,
        data_addr: u64,
        data_len: u64,
    ) -> std::io::Result<u64> {
        // Dataspace (0x0001), version 1
        let rank = u8::try_from(dims.len()).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("HDF5 dataspace rank {} exceeds u8::MAX", dims.len()),
            )
        })?;
        let mut ds = vec![1u8, rank, 0, 0, 0, 0, 0, 0];
        for &d in dims {
            ds.extend_from_slice(&d.to_le_bytes());
        }

        // Datatype (0x0003), Opaque class (5), version 1, element size at [4..8]
        let elem_u32 = u32::try_from(elem_size).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("HDF5 datatype element size {elem_size} exceeds u32::MAX"),
            )
        })?;
        let mut dt = vec![0u8; 8];
        dt[0] = 0x15;
        dt[4..8].copy_from_slice(&elem_u32.to_le_bytes());

        // DataLayout (0x0008), version 3, class 1 (contiguous)
        let layout_addr = if data_len == 0 { u64::MAX } else { data_addr };
        let mut dl = Vec::with_capacity(18);
        dl.push(3);
        dl.push(1);
        dl.extend_from_slice(&layout_addr.to_le_bytes());
        dl.extend_from_slice(&data_len.to_le_bytes());

        self.write_v1_obj_header(&[(0x0001, &ds), (0x0003, &dt), (0x0008, &dl)])
    }

    pub(super) fn write_dataset_stream<F>(
        &mut self,
        dims: &[u64],
        elem_size: usize,
        writer: F,
    ) -> std::io::Result<u64>
    where
        F: FnOnce(&mut W) -> std::io::Result<u64>,
    {
        self.align8()?;
        let data_addr = self.pos;
        let data_len = writer(&mut self.w)?;
        self.pos += data_len;
        self.write_dataset_header(dims, elem_size, data_addr, data_len)
    }

    pub(super) fn reserve_contiguous_dataset(
        &mut self,
        dims: &[u64],
        elem_size: usize,
        data_len: u64,
    ) -> std::io::Result<(u64, u64)> {
        self.align8()?;
        let data_addr = self.pos;
        if data_len > 0 {
            let end_pos = data_addr.checked_add(data_len).ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "HDF5 dataset extent overflows u64",
                )
            })?;
            self.w.seek(SeekFrom::Start(end_pos - 1))?;
            self.w.write_all(&[0u8])?;
            self.pos = end_pos;
        }
        let hdr_addr = self.write_dataset_header(dims, elem_size, data_addr, data_len)?;
        Ok((hdr_addr, data_addr))
    }

    pub(super) fn write_at(&mut self, offset: u64, bytes: &[u8]) -> std::io::Result<()> {
        if bytes.is_empty() {
            return Ok(());
        }
        let end = offset
            .checked_add(bytes.len() as u64)
            .filter(|&e| e <= self.pos)
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "HDF5 write_at exceeds reserved file extent",
                )
            })?;
        debug_assert!(end <= self.pos);
        self.w.seek(SeekFrom::Start(offset))?;
        self.w.write_all(bytes)?;
        self.w.seek(SeekFrom::Start(self.pos))?;
        Ok(())
    }

    pub(super) fn write_dataset(
        &mut self,
        dims: &[u64],
        elem_size: usize,
        data: &[u8],
    ) -> std::io::Result<u64> {
        self.write_dataset_stream(dims, elem_size, |w| {
            w.write_all(data)?;
            Ok(data.len() as u64)
        })
    }

    fn write_group_btree(&mut self, snods: Vec<(u64, u64)>) -> std::io::Result<u64> {
        // Group Internal Node K = 16 => max 32 child entries per TREE node.
        const MAX_TREE_CHILDREN: usize = 32;
        if snods.is_empty() {
            self.align8()?;
            let empty_addr = self.pos;
            let mut node = Vec::with_capacity(32 + MAX_TREE_CHILDREN * 16);
            node.extend_from_slice(b"TREE");
            node.push(0);
            node.push(0);
            node.extend_from_slice(&0u16.to_le_bytes());
            node.extend_from_slice(&u64::MAX.to_le_bytes());
            node.extend_from_slice(&u64::MAX.to_le_bytes());
            node.extend_from_slice(&0u64.to_le_bytes());
            node.resize(node.len() + MAX_TREE_CHILDREN * 16, 0);
            self.write_bytes(&node)?;
            return Ok(empty_addr);
        }

        let mut level_nodes: Vec<(u64, u64)> = snods;
        let mut level: u8 = 0;

        loop {
            const TREE_NODE_LEN: usize = 32 + MAX_TREE_CHILDREN * 16; // 544 bytes (8-byte aligned)
            const TREE_NODE_BYTES: u64 = TREE_NODE_LEN as u64;
            let num_chunks = level_nodes.len().div_ceil(MAX_TREE_CHILDREN);
            self.align8()?;
            let level_base_addr = self.pos;
            let mut next_level: Vec<(u64, u64)> = Vec::with_capacity(num_chunks);
            let mut prev_last_off = 0u64;
            for (chunk_idx, chunk) in level_nodes.chunks(MAX_TREE_CHILDREN).enumerate() {
                let node_addr = level_base_addr + (chunk_idx as u64) * TREE_NODE_BYTES;
                debug_assert_eq!(self.pos, node_addr);
                let left_sib = if chunk_idx > 0 {
                    level_base_addr + ((chunk_idx - 1) as u64) * TREE_NODE_BYTES
                } else {
                    u64::MAX
                };
                let right_sib = if chunk_idx + 1 < num_chunks {
                    level_base_addr + ((chunk_idx + 1) as u64) * TREE_NODE_BYTES
                } else {
                    u64::MAX
                };
                let mut node = Vec::with_capacity(TREE_NODE_LEN);
                node.extend_from_slice(b"TREE");
                node.push(0); // type 0 (group)
                node.push(level);
                node.extend_from_slice(
                    &u16::try_from(chunk.len()).unwrap_or(u16::MAX).to_le_bytes(),
                );
                node.extend_from_slice(&left_sib.to_le_bytes());
                node.extend_from_slice(&right_sib.to_le_bytes());
                node.extend_from_slice(&prev_last_off.to_le_bytes());
                for &(child_addr, last_off) in chunk {
                    node.extend_from_slice(&child_addr.to_le_bytes());
                    node.extend_from_slice(&last_off.to_le_bytes());
                }
                let end_key = chunk.last().map_or(0, |&(_, last_off)| last_off);
                node.resize(TREE_NODE_LEN, 0);
                self.write_bytes(&node)?;
                next_level.push((node_addr, end_key));
                prev_last_off = end_key;
            }
            if next_level.len() == 1 {
                return Ok(next_level[0].0);
            }
            level_nodes = next_level;
            level = level.saturating_add(1);
        }
    }

    pub(super) fn write_gcol(&mut self, strings: &[&str]) -> std::io::Result<u64> {
        self.align8()?;
        let gcol_addr = self.pos;
        let mut body = Vec::new();
        for (i, s) in strings.iter().enumerate() {
            let idx = u16::try_from(i + 1).map_err(|_| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    format!("HDF5 global heap object index {} exceeds u16::MAX", i + 1),
                )
            })?;
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

        let mut gcol = Vec::with_capacity(total_size);
        gcol.extend_from_slice(b"GCOL");
        gcol.push(1); // version 1
        gcol.extend_from_slice(&[0u8; 3]); // reserved
        gcol.extend_from_slice(&(total_size as u64).to_le_bytes());
        gcol.extend_from_slice(&body);
        // End-of-collection free object (index = 0)
        gcol.extend_from_slice(&0u16.to_le_bytes());
        gcol.extend_from_slice(&0u16.to_le_bytes());
        gcol.extend_from_slice(&[0u8; 4]);
        gcol.extend_from_slice(&(free_space as u64).to_le_bytes());
        gcol.resize(gcol.len() + free_space - 16, 0);
        self.write_bytes(&gcol)?;
        Ok(gcol_addr)
    }

    pub(super) fn write_group_with_attrs(
        &mut self,
        entries: &[(&str, u64)],
        attrs: &[(&str, u32, u64, u32)],
    ) -> std::io::Result<(u64, u64, u64)> {
        // Group Leaf Node K = 4 => max 8 entries per SNOD
        const MAX_SNOD: usize = 8;
        self.align8()?;
        let mut sorted: Vec<(&str, u64)> = entries.to_vec();
        sorted.sort_by_key(|&(name, _)| name);
        if let Some(dup) = sorted.windows(2).find(|w| w[0].0 == w[1].0) {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("duplicate HDF5 group entry name `{}`", dup[0].0),
            ));
        }

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

        let heap_data_addr = self.pos;
        self.write_bytes(&heap_data)?;

        self.align8()?;
        let heap_addr = self.pos;
        let mut heap_hdr = Vec::with_capacity(32);
        heap_hdr.extend_from_slice(b"HEAP");
        heap_hdr.extend_from_slice(&[0u8; 4]); // version 0 + reserved
        heap_hdr.extend_from_slice(&(heap_data.len() as u64).to_le_bytes());
        heap_hdr.extend_from_slice(&1u64.to_le_bytes()); // H5HL_FREE_NULL (1): no free blocks
        heap_hdr.extend_from_slice(&heap_data_addr.to_le_bytes());
        self.write_bytes(&heap_hdr)?;

        let mut snods: Vec<(u64, u64)> = Vec::new();
        for (chunk_idx, chunk) in sorted.chunks(MAX_SNOD).enumerate() {
            self.align8()?;
            let snod_addr = self.pos;
            let mut snod = Vec::with_capacity(8 + MAX_SNOD * 40);
            snod.extend_from_slice(b"SNOD");
            snod.push(1); // version 1
            snod.push(0);
            snod.extend_from_slice(&u16::try_from(chunk.len()).unwrap_or(u16::MAX).to_le_bytes());
            for (local_idx, &(_, obj_addr)) in chunk.iter().enumerate() {
                let idx = chunk_idx * MAX_SNOD + local_idx;
                let mut rec = [0u8; 40];
                rec[0..8].copy_from_slice(&name_offsets[idx].to_le_bytes());
                rec[8..16].copy_from_slice(&obj_addr.to_le_bytes());
                snod.extend_from_slice(&rec);
            }
            // Pad SNOD to 2 * leaf_k = 8 slots so libhdf5 reads a full node cleanly.
            snod.resize(snod.len() + (MAX_SNOD - chunk.len()) * 40, 0);
            self.write_bytes(&snod)?;
            let base_idx = chunk_idx * MAX_SNOD;
            let last_off = name_offsets[base_idx + chunk.len() - 1];
            snods.push((snod_addr, last_off));
        }

        let btree_addr = self.write_group_btree(snods)?;

        let mut symtab = [0u8; 16];
        symtab[0..8].copy_from_slice(&btree_addr.to_le_bytes());
        symtab[8..16].copy_from_slice(&heap_addr.to_le_bytes());

        let attr_msgs: Vec<Vec<u8>> = attrs
            .iter()
            .map(|&(k, len, gcol_addr, obj_idx)| encode_vls_attr_msg(k, len, gcol_addr, obj_idx))
            .collect::<std::io::Result<Vec<_>>>()?;
        let mut msgs: Vec<(u16, &[u8])> = Vec::with_capacity(1 + attr_msgs.len());
        msgs.push((0x0011, &symtab));
        for a in &attr_msgs {
            msgs.push((0x000C, a.as_slice()));
        }
        let hdr_addr = self.write_v1_obj_header(&msgs)?;
        Ok((hdr_addr, btree_addr, heap_addr))
    }

    pub(super) fn write_group(
        &mut self,
        entries: &[(&str, u64)],
    ) -> std::io::Result<(u64, u64, u64)> {
        self.write_group_with_attrs(entries, &[])
    }

    pub(super) fn finish(
        mut self,
        root_hdr: u64,
        root_bt: u64,
        root_hp: u64,
    ) -> std::io::Result<W> {
        let eof_addr = self.pos;
        let mut sb = [0u8; 96];
        sb[0..8].copy_from_slice(b"\x89HDF\r\n\x1a\n");
        sb[8] = 0; // superblock version 0
        sb[13] = 8; // offset_size
        sb[14] = 8; // length_size
        sb[16..18].copy_from_slice(&4u16.to_le_bytes()); // group leaf node K
        sb[18..20].copy_from_slice(&16u16.to_le_bytes()); // group internal node K
        sb[24..32].copy_from_slice(&0u64.to_le_bytes());
        sb[32..40].copy_from_slice(&u64::MAX.to_le_bytes());
        sb[40..48].copy_from_slice(&eof_addr.to_le_bytes());
        sb[48..56].copy_from_slice(&u64::MAX.to_le_bytes());
        // Root group symbol-table entry at [56..96]
        sb[56..64].copy_from_slice(&0u64.to_le_bytes());
        sb[64..72].copy_from_slice(&root_hdr.to_le_bytes());
        sb[72..76].copy_from_slice(&1u32.to_le_bytes());
        sb[80..88].copy_from_slice(&root_bt.to_le_bytes());
        sb[88..96].copy_from_slice(&root_hp.to_le_bytes());
        self.w.seek(SeekFrom::Start(0))?;
        self.w.write_all(&sb)?;
        self.w.flush()?;
        Ok(self.w)
    }
}

fn encode_vls_attr_msg(
    name: &str,
    str_len: u32,
    gcol_addr: u64,
    obj_idx: u32,
) -> std::io::Result<Vec<u8>> {
    let name_bytes = name.as_bytes();
    let name_size = u16::try_from(name_bytes.len() + 1).map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("HDF5 attribute name exceeds 65535 bytes: {}", name.len()),
        )
    })?;

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
    Ok(msg)
}
