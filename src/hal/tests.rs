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

    assert_eq!(read_hdf5_superblock(&hdr[..15]), None);
    hdr[8] = 3;
    assert_eq!(
        read_hdf5_superblock(&hdr),
        Some(Hdf5Superblock {
            version: 3,
            offset_size: 8,
            length_size: 8,
            root_group_addr: None,
        })
    );
    assert_eq!(decode_hal_seq_names(&[0u8; 8], 0), Vec::<String>::new());

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

    // Msg 3: DataLayout (0x0008) v3 class=2 ndim=2 (1D + element-size dim), size=12
    buf[p..p + 2].copy_from_slice(&0x0008_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&12_u16.to_le_bytes());
    buf[p + 8] = 3; // layout_ver
    buf[p + 9] = 2; // chunked
    buf[p + 10] = 2; // ndim = 2 (rank 1 + element-size dimension)
    buf[p + 11..p + 19].copy_from_slice(&0x100_u64.to_le_bytes());

    // Chunk B-tree v1 type-1 at 0x100 (ndim=2 => key_bytes=24, child_addr at +24)
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

#[test]
#[ignore = "requires HORIZON_TEST_HAL=/path/to/file.hal"]
fn real_hal_decode_from_env() {
    let path = std::env::var("HORIZON_TEST_HAL").expect("set HORIZON_TEST_HAL");
    let mut f = std::io::BufReader::new(std::fs::File::open(&path).expect("open hal file"));
    let hdr = read_at(&mut f, 0, 512).expect("read header");
    let sb = read_hdf5_superblock(&hdr).expect("superblock");
    let root_addr = sb.root_group_addr.expect("root addr");
    let root = read_hdf5_v0_group_entries_seek(&mut f, root_addr)
        .expect("seek root")
        .expect("root entries");
    assert_ne!(root.len(), 0);
    let first_genome = root
        .iter()
        .find(|e| !matches!(e.name.as_str(), "Meta" | "Phylogeny" | "Genomes" | "Verison"))
        .expect("genome group");
    let g_entries = read_hdf5_v0_group_entries_seek(&mut f, first_genome.object_header_addr)
        .expect("seek genome")
        .expect("genome entries");
    let seqname_obj = g_entries
        .iter()
        .find(|e| e.name == "SEQNAME_ARRAY")
        .expect("SEQNAME_ARRAY");
    let sn_dims = read_hdf5_dataspace_dims_seek(&mut f, seqname_obj.object_header_addr)
        .expect("seek dims")
        .expect("SEQNAME dims");
    let sn_stride = read_hdf5_datatype_size_seek(&mut f, seqname_obj.object_header_addr)
        .expect("seek stride")
        .expect("SEQNAME stride");
    let sn_raw = read_hdf5_v0_dataset_bytes_seek(&mut f, seqname_obj.object_header_addr)
        .expect("seek bytes")
        .expect("SEQNAME bytes");
    let names = decode_hal_seq_names(&sn_raw, sn_stride);
    assert_eq!(names.len(), usize::try_from(sn_dims[0]).expect("len"));
    assert_ne!(names[0], "");
}

#[test]
fn hal_hdf5_write_and_seek_roundtrip() {
    let ref_dna = b"ACGTACGTACGT";
    let q1_dna = b"ACGTACGTACGT";
    let hal = HalAlignment::from_star_blocks(
        "hg38",
        ref_dna,
        &[("mm10".to_owned(), q1_dna.to_vec())],
        &[(0, 0, 0, 12, false)],
    );
    let bytes = hal.to_hdf5_bytes();
    assert_eq!(u16::from_le_bytes([bytes[16], bytes[17]]), 4);
    assert_eq!(u16::from_le_bytes([bytes[18], bytes[19]]), 16);
    let sb = read_hdf5_superblock(&bytes).expect("superblock");
    let root_addr = sb.root_group_addr.expect("root addr");
    let mut cursor = std::io::Cursor::new(&bytes);
    let root_names = read_hdf5_v0_root_names_seek(&mut cursor, sb)
        .expect("root names")
        .expect("some names");
    assert_eq!(
        root_names,
        vec![
            "Genomes".to_owned(),
            "Meta".to_owned(),
            "Phylogeny".to_owned(),
            "Verison".to_owned(),
            "hg38".to_owned(),
            "mm10".to_owned(),
        ]
    );

    let root_entries = read_hdf5_v0_group_entries_seek(&mut cursor, root_addr)
        .expect("root entries")
        .expect("some entries");
    let hg38 = root_entries
        .iter()
        .find(|e| e.name == "hg38")
        .expect("hg38");
    let hg38_ds = read_hdf5_v0_group_entries_seek(&mut cursor, hg38.object_header_addr)
        .expect("hg38 ds")
        .expect("some ds");

    let dna_obj = hg38_ds
        .iter()
        .find(|e| e.name == "DNA_ARRAY")
        .expect("DNA_ARRAY");
    let dna_raw = read_hdf5_v0_dataset_bytes_seek(&mut cursor, dna_obj.object_header_addr)
        .expect("dna read")
        .expect("dna bytes");
    assert_eq!(decode_hal_dna(&dna_raw, ref_dna.len()), ref_dna);

    let sn_obj = hg38_ds
        .iter()
        .find(|e| e.name == "SEQNAME_ARRAY")
        .expect("SEQNAME_ARRAY");
    let sn_stride = read_hdf5_datatype_size_seek(&mut cursor, sn_obj.object_header_addr)
        .expect("stride")
        .expect("some stride");
    let sn_raw = read_hdf5_v0_dataset_bytes_seek(&mut cursor, sn_obj.object_header_addr)
        .expect("sn read")
        .expect("sn bytes");
    assert_eq!(decode_hal_seq_names(&sn_raw, sn_stride), vec!["hg38.chr1"]);

    let bot_obj = hg38_ds
        .iter()
        .find(|e| e.name == "BOTTOM_ARRAY")
        .expect("BOTTOM_ARRAY");
    let bot_raw = read_hdf5_v0_dataset_bytes_seek(&mut cursor, bot_obj.object_header_addr)
        .expect("bot read")
        .expect("bot bytes");
    assert_eq!(
        decode_hal_bottom_array(&bot_raw, 1),
        vec![(0, 12, -1, vec![(0, false)]), (12, 0, -1, vec![(-1, false)])]
    );
}

#[test]
fn hal_hdf5_write_group_exceeding_256_entries_roundtrips() {
    let ref_dna = b"ACGTACGT";
    let queries: Vec<(String, Vec<u8>)> = (0..300)
        .map(|i| (format!("q_{i:04}"), ref_dna.to_vec()))
        .collect();
    let hal = HalAlignment::from_star_blocks("ref0", ref_dna, &queries, &[]);
    let bytes = hal.to_hdf5_bytes();
    let sb = read_hdf5_superblock(&bytes).expect("superblock");
    let root_addr = usize::try_from(sb.root_group_addr.expect("root addr")).expect("usize");
    let entries = read_hdf5_v0_group_entries(&bytes, root_addr).expect("root entries");
    assert_eq!(entries.len(), 305);
    assert!(entries.iter().any(|e| e.name == "ref0"));
    assert!(entries.iter().any(|e| e.name == "q_0000"));
    assert!(entries.iter().any(|e| e.name == "q_0299"));

    // Verify HDF5 B-tree v1 group keys satisfy Key[0] == 0 on leftmost node and
    // Key[i+1] == last_off(Child[i]) so libhdf5 lookups work across > 8 entries.
    let root_btree_addr = usize::try_from(u64::from_le_bytes(
        bytes[80..88].try_into().expect("btree addr"),
    ))
    .expect("usize");
    assert_eq!(&bytes[root_btree_addr..root_btree_addr + 4], b"TREE");
    let root_level = bytes[root_btree_addr + 5];
    assert_eq!(root_level, 1);
    let root_used = usize::from(u16::from_le_bytes(
        bytes[root_btree_addr + 6..root_btree_addr + 8]
            .try_into()
            .expect("used"),
    ));
    let key0 = u64::from_le_bytes(
        bytes[root_btree_addr + 24..root_btree_addr + 32]
            .try_into()
            .expect("key0"),
    );
    assert_eq!(key0, 0);
    let mut prev_end = 0u64;
    for i in 0..root_used {
        let off = root_btree_addr + 32 + i * 16;
        let child_addr = usize::try_from(u64::from_le_bytes(
            bytes[off..off + 8].try_into().expect("child"),
        ))
        .expect("usize");
        let key_after = u64::from_le_bytes(bytes[off + 8..off + 16].try_into().expect("key"));
        assert!(key_after > prev_end);
        prev_end = key_after;

        assert_eq!(&bytes[child_addr..child_addr + 4], b"TREE");
        let child_used = usize::from(u16::from_le_bytes(
            bytes[child_addr + 6..child_addr + 8]
                .try_into()
                .expect("c_used"),
        ));
        let child_key0 = u64::from_le_bytes(
            bytes[child_addr + 24..child_addr + 32]
                .try_into()
                .expect("c_key0"),
        );
        if i == 0 {
            assert_eq!(child_key0, 0);
        } else {
            assert!(child_key0 > 0);
        }
        let child_last_key = u64::from_le_bytes(
            bytes[child_addr + 24 + child_used * 16..child_addr + 32 + child_used * 16]
                .try_into()
                .expect("c_last"),
        );
        assert_eq!(child_last_key, key_after);
    }
}

#[test]
fn hal_leaf_to_maf_paralog_chaining_and_ohdr_v2_flags() {
    let ref_dna = b"ACGTACGTACGT";
    let q1_dna = b"ACGTACGTACGTACGT";
    // Two hits for q_idx=0 on the same ref interval (0, 8) exercises next_paralog chaining,
    // plus one out-of-range q_idx=99 that must be safely ignored.
    let hal = HalAlignment::from_star_blocks(
        "hg38",
        ref_dna,
        &[("mm10".to_owned(), q1_dna.to_vec())],
        &[
            (0, 0, 0, 8, true),
            (0, 0, 8, 8, false),
            (99, 0, 0, 8, false),
        ],
    );
    assert_eq!(hal.genomes[1].top_segments.len(), 2);
    assert_eq!(hal.genomes[1].top_segments[1].next_paralog, Some(0));

    let leaf_maf = hal.to_maf_for_genome("mm10").expect("leaf MAF");
    assert!(leaf_maf.contains("s mm10 0 8 + 16 ACGTACGT"));
    assert!(leaf_maf.contains("s hg38 4 8 - 12 ACGTACGT"));

    let root_maf = hal.to_maf_for_genome("hg38").expect("root MAF");
    assert!(root_maf.contains("s mm10 8 8 - 16 ACGTACGT"));
    assert!(root_maf.contains("s mm10 8 8 + 16 ACGTACGT"));

    // Exercise HDF5 v2 OHDR with flags & 0x04 != 0 (6-byte message headers with creation order).
    let ohdr_v2: Vec<u8> = vec![
        b'O', b'H', b'D', b'R', // magic
        2,    // version 2
        0x04, // flags: 1-byte chunk size (0), creation order tracked (0x04)
        12,   // chunk size = 12 bytes (two 6-byte empty messages)
        0x01, 0, 0, 0, 0, 0, // msg type 1, size 0, flags 0, creation_order 0
        0x08, 0, 0, 0, 1, 0, // msg type 8, size 0, flags 0, creation_order 1
        0, 0, 0, 0, // padding to >= 16 bytes
    ];
    let parsed = read_hdf5_object_header(&ohdr_v2).expect("OHDR v2");
    assert_eq!(parsed.version, 2);
    assert_eq!(parsed.message_types, vec![0x01, 0x08]);
}

#[test]
fn hal_superblock_v1_and_inspect_and_chunk_padding() {
    // Superblock v1 places root_group_addr at offset 68..76 (shifted +4 vs v0).
    let mut sb_v1 = vec![0u8; 80];
    sb_v1[..8].copy_from_slice(&HDF5_SIGNATURE);
    sb_v1[8] = 1;
    sb_v1[13] = 8;
    sb_v1[14] = 8;
    sb_v1[68..76].copy_from_slice(&0x1234_u64.to_le_bytes());
    let parsed_v1 = read_hdf5_superblock(&sb_v1).expect("superblock v1");
    assert_eq!(parsed_v1.version, 1);
    assert_eq!(parsed_v1.root_group_addr, Some(0x1234));

    // Exercise inspect_hal_genomes on a temp file.
    let tmp_path = std::env::temp_dir().join(format!(
        "horizon_hal_inspect_{}_{:?}.hal",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::write(&tmp_path, &sb_v1).expect("write temp hal");
    let inspected = inspect_hal_genomes(&tmp_path).expect("inspect_hal_genomes");
    let _ = std::fs::remove_file(&tmp_path);
    assert_eq!(inspected.map(|info| info.superblock), Some(parsed_v1));

    // Chunked dataset with Dataspace (3 x 4 bytes = 12 bytes) and 16-byte chunk trims padding.
    let mut buf = vec![0u8; 0x300];
    buf[0x20] = 1;
    buf[0x22..0x24].copy_from_slice(&4_u16.to_le_bytes());
    buf[0x28..0x2c].copy_from_slice(&80_u32.to_le_bytes());

    let mut p = 0x30;
    // Msg 1: Dataspace (0x0001) v1 rank=1 dim=3
    buf[p..p + 2].copy_from_slice(&0x0001_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&16_u16.to_le_bytes());
    buf[p + 8] = 1;
    buf[p + 9] = 1;
    buf[p + 16..p + 24].copy_from_slice(&3_u64.to_le_bytes());
    p += 24;

    // Msg 2: Datatype (0x0003) size=4
    buf[p..p + 2].copy_from_slice(&0x0003_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&8_u16.to_le_bytes());
    buf[p + 12..p + 16].copy_from_slice(&4_u32.to_le_bytes());
    p += 16;

    // Msg 3: DataLayout (0x0008) v3 class=2 ndim=2 -> TREE@0x100
    buf[p..p + 2].copy_from_slice(&0x0008_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&12_u16.to_le_bytes());
    buf[p + 8] = 3;
    buf[p + 9] = 2;
    buf[p + 10] = 2;
    buf[p + 11..p + 19].copy_from_slice(&0x100_u64.to_le_bytes());

    buf[0x100..0x104].copy_from_slice(b"TREE");
    buf[0x104] = 1;
    buf[0x105] = 0;
    buf[0x106..0x108].copy_from_slice(&1_u16.to_le_bytes());
    buf[0x118..0x11c].copy_from_slice(&16_u32.to_le_bytes());
    buf[0x130..0x138].copy_from_slice(&0x200_u64.to_le_bytes());
    buf[0x200..0x210].copy_from_slice(b"ABCDEFGHIJKLMNOP");

    assert_eq!(
        read_hdf5_v0_dataset_bytes(&buf, 0x20).as_deref(),
        Some(b"ABCDEFGHIJKL".as_slice())
    );

    // Unsupported filter (filter_id=2 shuffle) must return None instead of corrupt bytes.
    buf[0x22..0x24].copy_from_slice(&4_u16.to_le_bytes());
    buf[0x28..0x2c].copy_from_slice(&88_u32.to_le_bytes());
    let fp = p + 20;
    buf[fp..fp + 2].copy_from_slice(&0x000B_u16.to_le_bytes());
    buf[fp + 2..fp + 4].copy_from_slice(&8_u16.to_le_bytes());
    buf[fp + 8] = 2;
    buf[fp + 9] = 1;
    buf[fp + 10..fp + 12].copy_from_slice(&2_u16.to_le_bytes());
    assert_eq!(read_hdf5_v0_dataset_bytes(&buf, 0x20), None);

    // Unpadded 18-byte DataLayout v3 class 1 (contiguous) and size=0 empty dataset.
    let mut cbuf = vec![0u8; 0x80];
    cbuf[0x20] = 1;
    cbuf[0x22..0x24].copy_from_slice(&1_u16.to_le_bytes());
    cbuf[0x28..0x2c].copy_from_slice(&26_u32.to_le_bytes());
    let mp = 0x30;
    cbuf[mp..mp + 2].copy_from_slice(&0x0008_u16.to_le_bytes());
    cbuf[mp + 2..mp + 4].copy_from_slice(&18_u16.to_le_bytes());
    cbuf[mp + 8] = 3;
    cbuf[mp + 9] = 1;
    cbuf[mp + 10..mp + 18].copy_from_slice(&0x60_u64.to_le_bytes());
    cbuf[mp + 18..mp + 26].copy_from_slice(&5_u64.to_le_bytes());
    cbuf[0x60..0x65].copy_from_slice(b"HELLO");
    assert_eq!(
        read_hdf5_v0_dataset_bytes(&cbuf, 0x20).as_deref(),
        Some(b"HELLO".as_slice())
    );
    cbuf[mp + 10..mp + 18].copy_from_slice(&u64::MAX.to_le_bytes());
    cbuf[mp + 18..mp + 26].copy_from_slice(&0_u64.to_le_bytes());
    assert_eq!(
        read_hdf5_v0_dataset_bytes(&cbuf, 0x20).as_deref(),
        Some([].as_slice())
    );
}

#[test]
fn hdf5_v1_continuation_and_multilevel_chunk_btree() {
    // Test 0x0010 v1 object header continuation and node_level == 1 chunk B-tree.
    let mut buf = vec![0u8; 0x300];
    // v1 object header at 0x20 with 2 messages: 0x0010 continuation -> 0x80 (len 28)
    buf[0x20] = 1;
    buf[0x22..0x24].copy_from_slice(&2_u16.to_le_bytes());
    buf[0x28..0x2c].copy_from_slice(&24_u32.to_le_bytes());
    // Message 1 at 0x30: type 0x0010, size 16 -> cont_addr 0x80, cont_len 28
    buf[0x30..0x32].copy_from_slice(&0x0010_u16.to_le_bytes());
    buf[0x32..0x34].copy_from_slice(&16_u16.to_le_bytes());
    buf[0x38..0x40].copy_from_slice(&0x80_u64.to_le_bytes());
    buf[0x40..0x48].copy_from_slice(&28_u64.to_le_bytes());

    // Continuation chunk at 0x80: DataLayout v3 class 2 (chunked, ndim=2, btree at 0x100)
    buf[0x80..0x82].copy_from_slice(&0x0008_u16.to_le_bytes());
    buf[0x82..0x84].copy_from_slice(&20_u16.to_le_bytes());
    buf[0x88] = 3; // version 3
    buf[0x89] = 2; // chunked
    buf[0x8a] = 2; // ndim = 2 -> key_bytes = 8 + 16 = 24, rec_bytes = 32
    buf[0x8b..0x93].copy_from_slice(&0x100_u64.to_le_bytes());

    // Root chunk B-tree at 0x100: node_level = 1, n_entries = 1 -> leaf at 0x180
    buf[0x100..0x104].copy_from_slice(b"TREE");
    buf[0x104] = 1; // chunk B-tree
    buf[0x105] = 1; // node_level = 1
    buf[0x106..0x108].copy_from_slice(&1_u16.to_le_bytes());
    buf[0x118 + 24..0x118 + 32].copy_from_slice(&0x180_u64.to_le_bytes());

    // Leaf chunk B-tree at 0x180: node_level = 0, n_entries = 1 -> 4 bytes at 0x220
    buf[0x180..0x184].copy_from_slice(b"TREE");
    buf[0x184] = 1;
    buf[0x185] = 0;
    buf[0x186..0x188].copy_from_slice(&1_u16.to_le_bytes());
    buf[0x198..0x19c].copy_from_slice(&4_u32.to_le_bytes());
    buf[0x198 + 24..0x198 + 32].copy_from_slice(&0x220_u64.to_le_bytes());
    buf[0x220..0x224].copy_from_slice(b"ACGT");

    assert_eq!(
        read_hdf5_v0_dataset_bytes(&buf, 0x20).as_deref(),
        Some(b"ACGT".as_slice())
    );

    // Also test read_hdf5_object_header when object header starts at slice[0].
    buf[0x38..0x40].copy_from_slice(&0x60_u64.to_le_bytes());
    buf[0x80..0x9c].copy_within(0..28, 0); // already at 0x80 = 0x20 + 0x60!
    let hdr = read_hdf5_object_header(&buf[0x20..]).expect("parse v1 header");
    assert_eq!(hdr.version, 1);
    assert_eq!(hdr.message_types, vec![0x0008]);
}

#[test]
fn hdf5_dataspace_v2_and_filter_pipeline_v1_chunked_dataset() {
    use flate2::{Compression, write::ZlibEncoder};
    use std::io::Write as _;

    // Build a minimal v0 dataset at 0x20 using:
    // - Dataspace v2 (0x0001): ndims=1, flags=0, type=1, dim=4
    // - DataLayout v3 chunked (0x0008): chunk_dim=4, btree_addr=0x100
    // - FilterPipeline v1 (0x000b): deflate (id=1) + custom filter (id=256)
    let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(b"TGCA").expect("compress");
    let compressed = encoder.finish().expect("finish zlib");

    let mut buf = vec![0u8; 0x280];
    // Object header v1 at 0x20: 3 messages
    buf[0x20] = 1;
    buf[0x22..0x24].copy_from_slice(&3_u16.to_le_bytes());

    // Msg 1 at 0x30: Dataspace v2 (type 0x0001, size 12)
    // [ver=2, ndims=1, flags=0, type=1, dim0=4 (u64)]
    let mut p = 0x30;
    buf[p..p + 2].copy_from_slice(&0x0001_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&12_u16.to_le_bytes());
    p += 8;
    buf[p] = 2;
    buf[p + 1] = 1;
    buf[p + 2] = 0;
    buf[p + 3] = 1;
    buf[p + 4..p + 12].copy_from_slice(&4_u64.to_le_bytes());
    p += 12;

    // Msg 2: DataLayout v3 chunked (type 0x0008, size 19)
    buf[p..p + 2].copy_from_slice(&0x0008_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&19_u16.to_le_bytes());
    p += 8;
    buf[p] = 3; // layout v3
    buf[p + 1] = 2; // chunked
    buf[p + 2] = 2; // ndims + 1
    buf[p + 3..p + 11].copy_from_slice(&0x100_u64.to_le_bytes());
    buf[p + 11..p + 15].copy_from_slice(&4_u32.to_le_bytes());
    buf[p + 15..p + 19].copy_from_slice(&1_u32.to_le_bytes());
    p += 19;

    // Msg 3: FilterPipeline v1 (type 0x000b)
    // Header: 8 bytes [ver=1, nfilters=1, 6 reserved bytes]
    // Filter 0 (deflate id=1): 8 + 0 + 8 = 16 bytes (n_client_vals=1, padded to 8)
    // Total msg_size = 8 + 16 = 24
    buf[p..p + 2].copy_from_slice(&0x000b_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&24_u16.to_le_bytes());
    p += 8;
    buf[p] = 1; // ver=1
    buf[p + 1] = 1; // nfilters=1
    p += 8;
    // Filter 0: id=1 (deflate), name_len=0, flags=0, n_client_vals=1, client_val=6
    buf[p..p + 2].copy_from_slice(&1_u16.to_le_bytes());
    buf[p + 6..p + 8].copy_from_slice(&1_u16.to_le_bytes());
    buf[p + 8..p + 12].copy_from_slice(&6_u32.to_le_bytes());

    let obj_hdr_size = u32::try_from((p + 16) - 0x30).expect("u32");
    buf[0x28..0x2c].copy_from_slice(&obj_hdr_size.to_le_bytes());

    // Leaf chunk B-tree at 0x100: node_level=0, n_entries=1
    buf[0x100..0x104].copy_from_slice(b"TREE");
    buf[0x104] = 1;
    buf[0x105] = 0;
    buf[0x106..0x108].copy_from_slice(&1_u16.to_le_bytes());
    let comp_len = u32::try_from(compressed.len()).expect("u32");
    buf[0x118..0x11c].copy_from_slice(&comp_len.to_le_bytes());
    buf[0x118 + 24..0x118 + 32].copy_from_slice(&0x180_u64.to_le_bytes());
    buf[0x180..0x180 + compressed.len()].copy_from_slice(&compressed);

    assert_eq!(
        read_hdf5_v0_dataset_bytes(&buf, 0x20).as_deref(),
        Some(b"TGCA".as_slice())
    );
}

#[test]
fn multi_contig_hal_hdf5_write_and_seek_roundtrip() {
    let ref_contigs = vec![
        ("chr1".to_owned(), b"ACGTACGT".to_vec()),
        ("chr2".to_owned(), b"GGGGTTTT".to_vec()),
    ];
    let q_contigs = vec![
        ("q_c1".to_owned(), b"ACGTACGT".to_vec()),
        ("q_c2".to_owned(), b"GGGGTTTT".to_vec()),
    ];
    let hal = HalAlignment::from_multi_contig_star_blocks(
        "hg38",
        &ref_contigs,
        &[("mm10".to_owned(), q_contigs)],
        &[(0, 0, 0, 0, 0, 8, false), (0, 1, 1, 0, 0, 8, false)],
    );
    assert_eq!(hal.genomes[0].sequences.len(), 2);
    assert_eq!(hal.genomes[0].sequences[1].offset, 8);
    assert_eq!(hal.genomes[0].bottom_segments.len(), 2);
    assert_eq!(hal.genomes[0].bottom_segments[1].start, 8);

    let bytes = hal.to_hdf5_bytes();
    let sb = read_hdf5_superblock(&bytes).expect("superblock");
    let root_addr = sb.root_group_addr.expect("root addr");
    let mut cursor = std::io::Cursor::new(&bytes);
    let root_entries = read_hdf5_v0_group_entries_seek(&mut cursor, root_addr)
        .expect("root entries")
        .expect("some entries");
    let hg38 = root_entries.iter().find(|e| e.name == "hg38").expect("hg38");
    let hg38_ds = read_hdf5_v0_group_entries_seek(&mut cursor, hg38.object_header_addr)
        .expect("hg38 ds")
        .expect("some ds");
    let sn_obj = hg38_ds
        .iter()
        .find(|e| e.name == "SEQNAME_ARRAY")
        .expect("SEQNAME_ARRAY");
    let sn_stride = read_hdf5_datatype_size_seek(&mut cursor, sn_obj.object_header_addr)
        .expect("stride")
        .expect("stride");
    let sn_raw = read_hdf5_v0_dataset_bytes_seek(&mut cursor, sn_obj.object_header_addr)
        .expect("bytes")
        .expect("bytes");
    assert_eq!(decode_hal_seq_names(&sn_raw, sn_stride), vec!["chr1", "chr2"]);
}

#[test]
fn hal_hdf5_export_normalizes_gaps_overlaps_and_paralogy_cycles() {
    // Odd-length ref (15bp) with two overlapping query blocks ([0..10) and [4..12))
    // plus unaligned tail [12..15), and a paralogous block on q1.
    let ref_dna = b"ACGTACGTACGTACG";
    let q1_dna = b"ACGTACGTACNNACGTACGTAC";
    let q2_dna = b"ACGTACGT";
    let hal = HalAlignment::from_star_blocks(
        "hg38",
        ref_dna,
        &[
            ("q1".to_owned(), q1_dna.to_vec()),
            ("q2".to_owned(), q2_dna.to_vec()),
        ],
        &[
            (0, 0, 0, 10, false),
            (0, 0, 12, 10, false),
            (1, 4, 0, 8, false),
        ],
    );
    let bytes = hal.to_hdf5_bytes();
    assert_eq!(&bytes[96..100], b"GCOL");

    let sb = read_hdf5_superblock(&bytes).expect("superblock");
    let root_addr = sb.root_group_addr.expect("root addr");
    let mut cursor = std::io::Cursor::new(&bytes);
    let root_entries = read_hdf5_v0_group_entries_seek(&mut cursor, root_addr)
        .expect("root entries")
        .expect("some entries");
    for req in ["Meta", "Phylogeny", "Genomes", "Verison", "hg38", "q1", "q2"] {
        assert!(
            root_entries.iter().any(|e| e.name == req),
            "missing root entry {req}"
        );
    }

    let q1 = root_entries.iter().find(|e| e.name == "q1").expect("q1");
    let q1_ds = read_hdf5_v0_group_entries_seek(&mut cursor, q1.object_header_addr)
        .expect("q1 ds")
        .expect("some ds");
    assert!(q1_ds.iter().any(|e| e.name == "Meta"));
    assert!(q1_ds.iter().any(|e| e.name == "Rup"));

    let top_obj = q1_ds
        .iter()
        .find(|e| e.name == "TOP_ARRAY")
        .expect("TOP_ARRAY");
    let top_raw = read_hdf5_v0_dataset_bytes_seek(&mut cursor, top_obj.object_header_addr)
        .expect("top read")
        .expect("top bytes");
    let tops = decode_hal_top_array(&top_raw);
    // Last entry is sentinel with start == q1_dna.len()
    assert_eq!(tops.last().expect("sentinel").0, 22);
    // First two segments ([0..4), [4..10)) and their paralogs ([12..16), [16..22)) form 2-cycles
    let linked: Vec<_> = tops[..tops.len() - 1]
        .iter()
        .enumerate()
        .filter(|(_, t)| t.3 >= 0)
        .collect();
    assert_eq!(linked.len(), 4);
    for &(idx, t) in &linked {
        assert!(t.2 >= 0, "paralogous segment {idx} must have next_paralog");
    }
}

#[test]
fn hal_to_maf_handles_cyclic_next_paralog_without_duplicate_rows() {
    let mut hal = HalAlignment::from_star_blocks(
        "hg38",
        b"ACGTACGT",
        &[("q1".to_owned(), b"ACGTACGTACGTACGT".to_vec())],
        &[(0, 0, 0, 8, false), (0, 0, 8, 8, false)],
    );
    // Close the paralog chain into a cycle (0 -> 1 -> 0) as normalized HALs do.
    hal.genomes[1].top_segments[0].next_paralog = Some(1);
    let maf = hal.to_maf_for_genome("hg38").expect("maf");
    let q1_rows = maf.lines().filter(|l| l.starts_with("s q1 ")).count();
    assert_eq!(q1_rows, 2, "cyclic next_paralog must emit each paralog once:\n{maf}");
}
