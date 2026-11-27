use super::tests_hdf5::read_hal_array_bytes;
use super::*;

#[test]
fn progressive_multi_contig_hal_builds_ancestral_tree_and_roundtrips() {
    let ref_contigs = vec![
        ("c1".to_owned(), b"ACGTACGT".to_vec()),
        ("c2".to_owned(), b"TTTTCCCC".to_vec()),
    ];
    let q1_contigs = vec![("q1_c1".to_owned(), b"ACGTACGT".to_vec())];
    let q2_contigs = vec![("q2_c2".to_owned(), b"GGGGAAAA".to_vec())];
    let hal = HalAlignment::from_multi_contig_progressive_blocks(
        "hg38",
        &ref_contigs,
        &[
            ("panTro4".to_owned(), q1_contigs),
            ("mm10".to_owned(), q2_contigs),
        ],
        &[(0, 0, 0, 0, 0, 8, false), (1, 1, 0, 0, 0, 8, true)],
    );
    assert_eq!(hal.tree_newick, "((hg38,panTro4)Anc0,mm10)Anc1;");
    assert_eq!(hal.genomes.len(), 5);
    assert_eq!(hal.genomes[0].name, "Anc1");
    assert_eq!(hal.genomes[0].bottom_segments.len(), 2);
    // Verify two-pass Fitch parsimony ancestral inference:
    // 1) If panTro4 and mm10 share 'G' at pos 0 while hg38 has 'A', BOTH Anc1 (root) AND Anc0
    //    must infer 'G' (1 substitution on hg38 instead of 2).
    // 2) If (hg38, sp1) share 'C' while outgroups (sp2, sp3) share 'T', Anc0 must infer 'C' while
    //    Anc1 and Anc2 infer 'T' (exact phylogenetic clade split).
    let consensus_hal = HalAlignment::from_multi_contig_progressive_blocks(
        "hg38",
        &[("c1".to_owned(), b"ACAAAAAA".to_vec())],
        &[
            (
                "panTro4".to_owned(),
                vec![("q1".to_owned(), b"GCAAAAAA".to_vec())],
            ),
            (
                "mm10".to_owned(),
                vec![("q2".to_owned(), b"GTAAAAAA".to_vec())],
            ),
            (
                "canFam3".to_owned(),
                vec![("q3".to_owned(), b"GTAAAAAA".to_vec())],
            ),
        ],
        &[
            (0, 0, 0, 0, 0, 8, false),
            (1, 0, 0, 0, 0, 8, false),
            (2, 0, 0, 0, 0, 8, false),
        ],
    );
    // Genomes layout: [Anc2 (root), Anc1, Anc0, hg38, panTro4, mm10, canFam3]
    assert_eq!(consensus_hal.genomes[0].name, "Anc2");
    assert_eq!(consensus_hal.genomes[1].name, "Anc1");
    assert_eq!(consensus_hal.genomes[2].name, "Anc0");
    // Site 0: hg38=A, panTro4=G, mm10=G, canFam3=G -> Anc2=G, Anc1=G, Anc0=G
    assert_eq!(consensus_hal.genomes[0].dna[0], b'G');
    assert_eq!(consensus_hal.genomes[1].dna[0], b'G');
    assert_eq!(consensus_hal.genomes[2].dna[0], b'G');
    // Site 1: hg38=C, panTro4=C, mm10=T, canFam3=T -> Anc2=T, Anc1=T, Anc0=C
    assert_eq!(consensus_hal.genomes[0].dna[1], b'T');
    assert_eq!(consensus_hal.genomes[1].dna[1], b'T');
    assert_eq!(consensus_hal.genomes[2].dna[1], b'C');

    let bytes = hal.to_hdf5_bytes();
    let sb = read_hdf5_superblock(&bytes).expect("superblock");
    let root_addr = sb.root_group_addr.expect("root addr");
    let mut cursor = std::io::Cursor::new(&bytes);
    let root_entries = read_hdf5_v0_group_entries_seek(&mut cursor, root_addr)
        .expect("root entries")
        .expect("some entries");
    assert!(root_entries.iter().any(|e| e.name == "Anc0"));
    assert!(root_entries.iter().any(|e| e.name == "Anc1"));
    assert!(root_entries.iter().any(|e| e.name == "hg38"));
    assert!(root_entries.iter().any(|e| e.name == "panTro4"));
    assert!(root_entries.iter().any(|e| e.name == "mm10"));
}

#[test]
#[allow(clippy::too_many_lines)]
fn streaming_progressive_hdf5_enforces_libhal_invariants_and_consensus() {
    let ref_contigs = vec![
        ("c1".to_owned(), b"AAAAAAAACCCC".to_vec()),
        ("c2".to_owned(), b"TTTTGGGG".to_vec()),
    ];
    let queries = vec![
        (
            "sp1".to_owned(),
            vec![
                ("q1_c1".to_owned(), b"GAAAAAAACCCCTTTT".to_vec()),
                ("q1_c2".to_owned(), b"GAAAAAAA".to_vec()),
            ],
        ),
        (
            "sp2".to_owned(),
            vec![("q2_c1".to_owned(), b"CCCCAAAAGAAAAAAA".to_vec())],
        ),
        (
            "sp3".to_owned(),
            vec![("q3_c1".to_owned(), b"GAAAAAAACCCC".to_vec())],
        ),
    ];
    // Includes unsorted query coordinates, paralogs on sp1, overlapping ref blocks, and RC hit on sp2.
    let blocks = vec![
        (0, 0, 1, 0, 0, 8, false),
        (0, 0, 0, 0, 0, 8, false),
        (1, 0, 0, 0, 8, 8, false),
        (1, 1, 0, 0, 0, 8, true),
        (2, 0, 0, 0, 0, 12, false),
    ];

    let tmp_path = std::env::temp_dir().join(format!(
        "horizon_streaming_prog_{}_{:?}.hal",
        std::process::id(),
        std::thread::current().id()
    ));
    let (tree_newick, n_b) = HalAlignment::write_multi_contig_progressive_hdf5_hal(
        &tmp_path,
        "hg38",
        &ref_contigs,
        &queries,
        &blocks,
    )
    .expect("streaming progressive hal write");
    let bytes = std::fs::read(&tmp_path).expect("read streamed hal");
    let _ = std::fs::remove_file(&tmp_path);

    assert_eq!(tree_newick, "(((hg38,sp1)Anc0,sp2)Anc1,sp3)Anc2;");
    // Partitioned ref [0..20) with breakpoints {0, 8, 12, 20} -> 3 contiguous bottom segments.
    assert_eq!(n_b, 3);

    // Anc2 (root) must infer 'G' at pos 0 via majority vote across sp1/sp2/sp3.
    let anc2_dna = decode_hal_dna(&read_hal_array_bytes(&bytes, "Anc2", "DNA_ARRAY"), 20);
    assert_eq!(anc2_dna[0], b'G');

    // Verify libHal TOP_ARRAY / BOTTOM_ARRAY invariants for every child->parent edge:
    // 1) tops partition [0, dna_len) strictly monotonically with sentinel at dna_len
    // 2) every aligned top segment's implicit length equals its parent bottom segment's length
    // 3) paralogy pointers form closed cycles
    let edges = [
        ("Anc0", "Anc1", 20usize),
        ("Anc1", "Anc2", 20usize),
        ("hg38", "Anc0", 20usize),
        ("sp1", "Anc0", 24usize),
        ("sp2", "Anc1", 16usize),
        ("sp3", "Anc2", 12usize),
    ];
    for (child, parent, child_len) in edges {
        let tops = decode_hal_top_array(&read_hal_array_bytes(&bytes, child, "TOP_ARRAY"));
        let bots =
            decode_hal_bottom_array(&read_hal_array_bytes(&bytes, parent, "BOTTOM_ARRAY"), 2);
        assert!(tops.len() >= 2, "{child} TOP_ARRAY missing segments");
        assert_eq!(tops[0].0, 0, "{child} TOP_ARRAY must start at 0");
        assert_eq!(
            tops.last().expect("sentinel").0,
            i64::try_from(child_len).expect("i64"),
            "{child} TOP_ARRAY sentinel must equal dna_len"
        );
        for i in 0..tops.len() - 1 {
            let len = tops[i + 1].0 - tops[i].0;
            assert!(len > 0, "{child} top[{i}] has non-positive length {len}");
            let parent_idx = tops[i].3;
            if parent_idx >= 0 {
                let b = &bots[usize::try_from(parent_idx).expect("usize")];
                assert_eq!(
                    u64::try_from(len).expect("u64"),
                    b.1,
                    "{child} top[{i}] length {len} != {parent} bot[{parent_idx}] length {}",
                    b.1
                );
            }
            if tops[i].2 >= 0 {
                let mut cur = usize::try_from(tops[i].2).expect("usize");
                let mut steps = 0usize;
                while cur != i && steps < tops.len() {
                    assert!(tops[cur].2 >= 0, "{child} open paralogy chain at {cur}");
                    cur = usize::try_from(tops[cur].2).expect("usize");
                    steps += 1;
                }
                assert_eq!(cur, i, "{child} paralogy chain did not close at {i}");
            }
        }
    }
    let sp1_tops = decode_hal_top_array(&read_hal_array_bytes(&bytes, "sp1", "TOP_ARRAY"));
    let paralog_count = sp1_tops[..sp1_tops.len() - 1]
        .iter()
        .filter(|t| t.2 >= 0)
        .count();
    assert_eq!(paralog_count, 2, "expected 2-cycle of paralogs in sp1");

    // C3 regression: internal progressive ancestors must link top <-> bottom parse edges
    // (`bottomIdx == i`, `topIdx == i`); root bottom `topIdx == -1`; leaves `bottomIdx == -1`.
    for internal in ["Anc0", "Anc1"] {
        let tops = decode_hal_top_array(&read_hal_array_bytes(&bytes, internal, "TOP_ARRAY"));
        let bots =
            decode_hal_bottom_array(&read_hal_array_bytes(&bytes, internal, "BOTTOM_ARRAY"), 2);
        assert_eq!(tops.len(), bots.len());
        for i in 0..tops.len() - 1 {
            let i_i64 = i64::try_from(i).expect("i64");
            assert_eq!(
                tops[i].1, i_i64,
                "{internal} top[{i}] bottomIdx must link to {i_i64}"
            );
            assert_eq!(
                bots[i].2, i_i64,
                "{internal} bot[{i}] topIdx must link to {i_i64}"
            );
        }
    }
    let root_bots =
        decode_hal_bottom_array(&read_hal_array_bytes(&bytes, "Anc2", "BOTTOM_ARRAY"), 2);
    for (i, b) in root_bots[..root_bots.len() - 1].iter().enumerate() {
        assert_eq!(b.2, -1, "root Anc2 bot[{i}] topIdx must be -1");
    }
    for leaf in ["hg38", "sp1", "sp2", "sp3"] {
        let tops = decode_hal_top_array(&read_hal_array_bytes(&bytes, leaf, "TOP_ARRAY"));
        for (i, t) in tops[..tops.len() - 1].iter().enumerate() {
            assert_eq!(t.1, -1, "leaf {leaf} top[{i}] bottomIdx must be -1");
        }
    }

    // In-memory `from_multi_contig_progressive_blocks` must now produce identical normalized
    // HDF5 datasets to the streaming progressive writer across every genome and array.
    let in_mem =
        HalAlignment::from_multi_contig_progressive_blocks("hg38", &ref_contigs, &queries, &blocks);
    let in_mem_bytes = in_mem.to_hdf5_bytes();
    assert_eq!(in_mem.tree_newick, tree_newick);
    assert_eq!(in_mem.genomes[0].bottom_segments.len(), n_b);
    for genome in ["Anc0", "Anc1", "Anc2", "hg38", "sp1", "sp2", "sp3"] {
        for ds in [
            "DNA_ARRAY",
            "SEQIDX_ARRAY",
            "SEQNAME_ARRAY",
            "TOP_ARRAY",
            "BOTTOM_ARRAY",
        ] {
            assert_eq!(
                read_hal_array_bytes(&in_mem_bytes, genome, ds),
                read_hal_array_bytes(&bytes, genome, ds),
                "mismatch in {genome}/{ds} between in-memory and streaming progressive HAL"
            );
        }
    }
}

#[test]
fn hdf5_chunked_reader_rejects_multidim_and_places_sparse_1d_chunks() {
    let mut buf = vec![0u8; 0x300];
    buf[0x20] = 1;
    buf[0x22..0x24].copy_from_slice(&3_u16.to_le_bytes());
    buf[0x28..0x2c].copy_from_slice(&60_u32.to_le_bytes());
    let mut p = 0x30;
    // Dataspace v1: rank=1, dim=12
    buf[p..p + 2].copy_from_slice(&0x0001_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&16_u16.to_le_bytes());
    buf[p + 8] = 1;
    buf[p + 9] = 1;
    buf[p + 16..p + 24].copy_from_slice(&12_u64.to_le_bytes());
    p += 24;
    // Datatype: size=1
    buf[p..p + 2].copy_from_slice(&0x0003_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&8_u16.to_le_bytes());
    buf[p + 12..p + 16].copy_from_slice(&1_u32.to_le_bytes());
    p += 16;
    // DataLayout v3 class=2 ndim=2 -> TREE@0x100
    let dl_off = p;
    buf[p..p + 2].copy_from_slice(&0x0008_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&12_u16.to_le_bytes());
    buf[p + 8] = 3;
    buf[p + 9] = 2;
    buf[p + 10] = 2;
    buf[p + 11..p + 19].copy_from_slice(&0x100_u64.to_le_bytes());

    // Sparse 1D chunk B-tree at 0x100: single chunk at coord[0]=8 with payload "ACGT".
    buf[0x100..0x104].copy_from_slice(b"TREE");
    buf[0x104] = 1;
    buf[0x105] = 0;
    buf[0x106..0x108].copy_from_slice(&1_u16.to_le_bytes());
    buf[0x118..0x11c].copy_from_slice(&4_u32.to_le_bytes());
    buf[0x120..0x128].copy_from_slice(&8_u64.to_le_bytes());
    buf[0x130..0x138].copy_from_slice(&0x200_u64.to_le_bytes());
    buf[0x200..0x204].copy_from_slice(b"ACGT");

    assert_eq!(
        read_hdf5_v0_dataset_bytes(&buf, 0x20).as_deref(),
        Some(b"\0\0\0\0\0\0\0\0ACGT".as_slice())
    );

    // Setting ndim=3 (2D chunked dataset) must be rejected with None.
    buf[dl_off + 10] = 3;
    assert_eq!(read_hdf5_v0_dataset_bytes(&buf, 0x20), None);
}

#[test]
fn hdf5_v1_stacked_deflate_and_filter_mask() {
    use flate2::{Compression, write::ZlibEncoder};
    use std::io::Write as _;

    let zlib_once = |data: &[u8]| -> Vec<u8> {
        let mut enc = ZlibEncoder::new(Vec::new(), Compression::default());
        enc.write_all(data).expect("zlib write");
        enc.finish().expect("zlib finish")
    };

    // Chunk 0: compressed twice (filter_mask = 0 -> both layers active).
    let c0_twice = zlib_once(&zlib_once(b"HELLO_"));
    // Chunk 1: compressed once (filter_mask = 0x2 -> layer 1 skipped, only 1 layer active).
    let c1_once = zlib_once(b"STACKED");

    let mut buf = vec![0u8; 0x400];
    // v1 object header at 0x20 with 2 messages: DataLayout (0x0008) + FilterPipeline (0x000B).
    buf[0x20] = 1;
    buf[0x22..0x24].copy_from_slice(&2_u16.to_le_bytes());
    buf[0x28..0x2c].copy_from_slice(&42_u32.to_le_bytes());

    // Msg 1 at 0x30: DataLayout v3 class 2 (chunked, ndim=1, btree at 0x100), msg_size=12 (total 20).
    let mut p = 0x30;
    buf[p..p + 2].copy_from_slice(&0x0008_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&12_u16.to_le_bytes());
    buf[p + 8] = 3;
    buf[p + 9] = 2;
    buf[p + 10] = 1;
    buf[p + 11..p + 19].copy_from_slice(&0x100_u64.to_le_bytes());

    // Msg 2 at 0x44: FilterPipeline v2 with 2 deflate (id=1) filters, msg_size=14 (total 22).
    p += 20;
    buf[p..p + 2].copy_from_slice(&0x000B_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&14_u16.to_le_bytes());
    buf[p + 8] = 2; // version 2
    buf[p + 9] = 2; // n_filters = 2
    buf[p + 10..p + 12].copy_from_slice(&1_u16.to_le_bytes()); // filter 0: id=1 (deflate)
    buf[p + 16..p + 18].copy_from_slice(&1_u16.to_le_bytes()); // filter 1: id=1 (deflate)

    // Chunk B-tree leaf at 0x100: ndim=1, 2 entries.
    buf[0x100..0x104].copy_from_slice(b"TREE");
    buf[0x104] = 1;
    buf[0x105] = 0;
    buf[0x106..0x108].copy_from_slice(&2_u16.to_le_bytes());
    buf[0x108..0x110].copy_from_slice(&u64::MAX.to_le_bytes());
    buf[0x110..0x118].copy_from_slice(&u64::MAX.to_le_bytes());

    // Entry 0 at 0x118 (ndim=1 -> key_bytes=16, rec_bytes=24): filter_mask=0, payload at 0x200.
    let c0_addr = 0x200usize;
    buf[0x118..0x11c].copy_from_slice(
        &u32::try_from(c0_twice.len())
            .expect("chunk len")
            .to_le_bytes(),
    );
    buf[0x11c..0x120].copy_from_slice(&0_u32.to_le_bytes());
    buf[0x128..0x130].copy_from_slice(&(c0_addr as u64).to_le_bytes());
    buf[c0_addr..c0_addr + c0_twice.len()].copy_from_slice(&c0_twice);

    // Entry 1 at 0x130: filter_mask=0x2 (bit 1 set), payload at 0x300.
    let c1_addr = 0x300usize;
    buf[0x130..0x134].copy_from_slice(
        &u32::try_from(c1_once.len())
            .expect("chunk len")
            .to_le_bytes(),
    );
    buf[0x134..0x138].copy_from_slice(&0x2_u32.to_le_bytes());
    buf[0x140..0x148].copy_from_slice(&(c1_addr as u64).to_le_bytes());
    buf[c1_addr..c1_addr + c1_once.len()].copy_from_slice(&c1_once);

    assert!(is_hdf5_deflated(&buf, 0x20));
    let mut cursor = std::io::Cursor::new(&buf);
    let out = read_hdf5_v0_dataset_bytes_seek(&mut cursor, 0x20)
        .expect("decode stacked deflate")
        .expect("some bytes");
    assert_eq!(out, b"HELLO_STACKED");

    // Adding a second 0x000B message with an unsupported filter must make both
    // is_hdf5_deflated and read_hdf5_v0_dataset_bytes_seek reject the dataset.
    buf[0x22..0x24].copy_from_slice(&3_u16.to_le_bytes());
    buf[0x28..0x2c].copy_from_slice(&58_u32.to_le_bytes());
    let fp2 = p + 22;
    buf[fp2..fp2 + 2].copy_from_slice(&0x000B_u16.to_le_bytes());
    buf[fp2 + 2..fp2 + 4].copy_from_slice(&8_u16.to_le_bytes());
    buf[fp2 + 8] = 2;
    buf[fp2 + 9] = 1;
    buf[fp2 + 10..fp2 + 12].copy_from_slice(&2_u16.to_le_bytes()); // id=2 (shuffle, unsupported)

    assert!(!is_hdf5_deflated(&buf, 0x20));
    let mut cursor2 = std::io::Cursor::new(&buf);
    assert_eq!(
        read_hdf5_v0_dataset_bytes_seek(&mut cursor2, 0x20).expect("no io error"),
        None
    );
}

#[test]
fn progressive_in_memory_and_streaming_agree_on_overlapping_same_interval_hits() {
    // C4 regression: overlapping query hits mapped to the same reference interval must follow
    // identical FIFO priority in both in-memory and streaming progressive HAL builders.
    let ref_contigs = vec![("c1".to_owned(), b"AAAAAAAACCCC".to_vec())];
    let queries = vec![
        (
            "sp1".to_owned(),
            vec![("q1".to_owned(), b"GGGGTTTTCCCC".to_vec())],
        ),
        (
            "sp2".to_owned(),
            vec![("q2".to_owned(), b"GGGGTTTTCCCC".to_vec())],
        ),
    ];
    let blocks = vec![
        (0, 0, 0, 0, 0, 8, false),
        (0, 0, 0, 0, 4, 8, false),
        (1, 0, 0, 0, 0, 8, false),
        (1, 0, 0, 0, 0, 0, false), // zero-length block must be ignored identically
    ];

    let tmp_path = std::env::temp_dir().join(format!(
        "horizon_c4_parity_{}_{:?}.hal",
        std::process::id(),
        std::thread::current().id()
    ));
    let (stream_tree, stream_nb) = HalAlignment::write_multi_contig_progressive_hdf5_hal(
        &tmp_path,
        "hg38",
        &ref_contigs,
        &queries,
        &blocks,
    )
    .expect("streaming write");
    let stream_bytes = std::fs::read(&tmp_path).expect("read streamed hal");
    let _ = std::fs::remove_file(&tmp_path);

    let in_mem =
        HalAlignment::from_multi_contig_progressive_blocks("hg38", &ref_contigs, &queries, &blocks);
    let in_mem_bytes = in_mem.to_hdf5_bytes();
    assert_eq!(in_mem.tree_newick, stream_tree);
    assert_eq!(in_mem.genomes[0].bottom_segments.len(), stream_nb);

    for genome in ["Anc0", "Anc1", "hg38", "sp1", "sp2"] {
        for ds in [
            "DNA_ARRAY",
            "SEQIDX_ARRAY",
            "SEQNAME_ARRAY",
            "TOP_ARRAY",
            "BOTTOM_ARRAY",
        ] {
            assert_eq!(
                read_hal_array_bytes(&in_mem_bytes, genome, ds),
                read_hal_array_bytes(&stream_bytes, genome, ds),
                "C4 parity mismatch in {genome}/{ds}"
            );
        }
    }
}

#[test]
fn chunked_dataset_rejects_overlapping_duplicate_chunk_zip_bomb() {
    // C5 regression: two chunk B-tree entries writing to the same coordinate offset must be rejected.
    use flate2::Compression;
    use flate2::write::ZlibEncoder;
    use std::io::Write as _;

    let mut enc = ZlibEncoder::new(Vec::new(), Compression::default());
    enc.write_all(b"ABCDEFGH").expect("deflate");
    let c_bytes = enc.finish().expect("finish");

    let mut buf = vec![0u8; 0x400];
    // Object header at 0x20: 3 messages (Datatype, DataLayout chunked, FilterPipeline).
    buf[0x20] = 1;
    buf[0x22..0x24].copy_from_slice(&3_u16.to_le_bytes());
    buf[0x28..0x2c].copy_from_slice(&56_u32.to_le_bytes());

    // Msg 1: Datatype (0x0003), elem_size = 1 byte.
    let mut p = 0x30;
    buf[p..p + 2].copy_from_slice(&0x0003_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&8_u16.to_le_bytes());
    buf[p + 12..p + 16].copy_from_slice(&1_u32.to_le_bytes());
    p += 16;

    // Msg 2: DataLayout v3 class=2 (chunked), ndim=1, btree_addr=0x100.
    buf[p..p + 2].copy_from_slice(&0x0008_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&16_u16.to_le_bytes());
    buf[p + 8] = 3;
    buf[p + 9] = 2;
    buf[p + 10] = 1;
    buf[p + 11..p + 19].copy_from_slice(&0x100_u64.to_le_bytes());
    p += 24;

    // Msg 3: FilterPipeline v2, 1 deflate layer.
    buf[p..p + 2].copy_from_slice(&0x000B_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&8_u16.to_le_bytes());
    buf[p + 8] = 2;
    buf[p + 9] = 1;
    buf[p + 10..p + 12].copy_from_slice(&1_u16.to_le_bytes());

    // Chunk B-tree leaf at 0x100: ndim=1, 2 entries both at offset_elem=4 (overlapping!).
    buf[0x100..0x104].copy_from_slice(b"TREE");
    buf[0x104] = 1;
    buf[0x105] = 0;
    buf[0x106..0x108].copy_from_slice(&2_u16.to_le_bytes());
    buf[0x108..0x110].copy_from_slice(&u64::MAX.to_le_bytes());
    buf[0x110..0x118].copy_from_slice(&u64::MAX.to_le_bytes());

    let c_addr = 0x200usize;
    buf[c_addr..c_addr + c_bytes.len()].copy_from_slice(&c_bytes);
    for (entry_off, coord) in [(0x118usize, 4_u64), (0x130usize, 4_u64)] {
        buf[entry_off..entry_off + 4]
            .copy_from_slice(&u32::try_from(c_bytes.len()).expect("len").to_le_bytes());
        buf[entry_off + 8..entry_off + 16].copy_from_slice(&coord.to_le_bytes());
        buf[entry_off + 16..entry_off + 24].copy_from_slice(&(c_addr as u64).to_le_bytes());
    }

    let mut cursor = std::io::Cursor::new(&buf);
    assert_eq!(
        read_hdf5_v0_dataset_bytes_seek(&mut cursor, 0x20).expect("no io error"),
        None,
        "overlapping chunk offsets at coord=4 must be rejected"
    );

    // R16 regression: duplicate chunks at offset_elem = 0 MUST also be rejected (not bypass overlap check).
    for entry_off in [0x118usize, 0x130usize] {
        buf[entry_off + 8..entry_off + 16].copy_from_slice(&0_u64.to_le_bytes());
    }
    let mut cursor_zero = std::io::Cursor::new(&buf);
    assert_eq!(
        read_hdf5_v0_dataset_bytes_seek(&mut cursor_zero, 0x20).expect("no io error"),
        None,
        "overlapping chunk offsets at coord=0 must be rejected"
    );
}
