use super::tests_hdf5::read_hal_array_bytes;
use super::*;

#[test]
fn fitch_progressive_preserves_n_scaffold_gaps_in_all_ancestors() {
    // C4 regression: N / scaffold-gap positions must remain 'N' in EVERY ancestor (not turn into 'A').
    let ref_contigs = vec![("chr1".to_owned(), b"ACGTNNNNACGT".to_vec())];
    let queries = vec![
        (
            "q0".to_owned(),
            vec![("c".to_owned(), b"ACGTNNNNACGT".to_vec())],
        ),
        (
            "q1".to_owned(),
            vec![("c".to_owned(), b"ACGTNNNNACGT".to_vec())],
        ),
        (
            "q2".to_owned(),
            vec![("c".to_owned(), b"ACGTNNNNACGT".to_vec())],
        ),
    ];
    let blocks = vec![
        (0, 0, 0, 0, 0, 12, false),
        (1, 0, 0, 0, 0, 12, false),
        (2, 0, 0, 0, 0, 12, false),
    ];
    let in_mem =
        HalAlignment::from_multi_contig_progressive_blocks("hg38", &ref_contigs, &queries, &blocks);
    for anc in ["Anc0", "Anc1", "Anc2"] {
        let g = in_mem
            .genomes
            .iter()
            .find(|g| g.name == anc)
            .expect("ancestor genome");
        assert_eq!(
            &g.dna,
            b"ACGTNNNNACGT",
            "in-memory {anc} corrupted N scaffold gap to {:?}",
            String::from_utf8_lossy(&g.dna)
        );
    }

    let tmp_path = std::env::temp_dir().join(format!(
        "horizon_c4_fitch_n_{}_{:?}.hal",
        std::process::id(),
        std::thread::current().id()
    ));
    HalAlignment::write_multi_contig_progressive_hdf5_hal(
        &tmp_path,
        "hg38",
        &ref_contigs,
        &queries,
        &blocks,
    )
    .expect("streaming write");
    let stream_bytes = std::fs::read(&tmp_path).expect("read hal");
    let _ = std::fs::remove_file(&tmp_path);
    let in_mem_bytes = in_mem.to_hdf5_bytes();
    for anc in ["Anc0", "Anc1", "Anc2"] {
        assert_eq!(
            read_hal_array_bytes(&stream_bytes, anc, "DNA_ARRAY"),
            read_hal_array_bytes(&in_mem_bytes, anc, "DNA_ARRAY"),
            "streaming {anc} DNA_ARRAY diverged on N scaffold gap"
        );
    }
}

#[test]
fn btree_reader_rejects_non_decreasing_level_and_bounds_overlapping_node_reads() {
    // C5 regression: a crafted HDF5 chunk B-tree where an internal level-1 node points to
    // another level-1 node (non-decreasing level) and overlapping nodes claim 65,535 entries
    // must terminate cleanly without amplification.
    let mut buf = vec![0u8; 0x2000];
    // Object header at 0x20: DataLayout v3 class=2 (chunked, ndim=1, btree at 0x100).
    buf[0x20] = 1;
    buf[0x22..0x24].copy_from_slice(&1_u16.to_le_bytes());
    buf[0x28..0x2c].copy_from_slice(&20_u32.to_le_bytes());
    buf[0x30..0x32].copy_from_slice(&0x0008_u16.to_le_bytes());
    buf[0x32..0x34].copy_from_slice(&12_u16.to_le_bytes());
    buf[0x38] = 3;
    buf[0x39] = 2;
    buf[0x3a] = 1;
    buf[0x3b..0x43].copy_from_slice(&0x100_u64.to_le_bytes());

    // Root TREE at 0x100: type=1, level=1, n_entries=2 -> points to 0x200 (also level=1, illegal!).
    buf[0x100..0x104].copy_from_slice(b"TREE");
    buf[0x104] = 1;
    buf[0x105] = 1;
    buf[0x106..0x108].copy_from_slice(&2_u16.to_le_bytes());
    buf[0x128..0x130].copy_from_slice(&0x200_u64.to_le_bytes());
    buf[0x140..0x148].copy_from_slice(&0x300_u64.to_le_bytes());

    // Node at 0x200: level=1 (non-decreasing -> must be ignored even though it points to a leaf).
    buf[0x200..0x204].copy_from_slice(b"TREE");
    buf[0x204] = 1;
    buf[0x205] = 1;
    buf[0x206..0x208].copy_from_slice(&1_u16.to_le_bytes());
    buf[0x228..0x230].copy_from_slice(&0x400_u64.to_le_bytes());

    // Node at 0x300: level=0, n_entries=65535 (claims ~1.5 MB read -> bounded cleanly).
    buf[0x300..0x304].copy_from_slice(b"TREE");
    buf[0x304] = 1;
    buf[0x305] = 0;
    buf[0x306..0x308].copy_from_slice(&u16::MAX.to_le_bytes());

    // Valid leaf at 0x400 that must NOT be reached through illegal level-1 -> level-1 descent.
    buf[0x400..0x404].copy_from_slice(b"TREE");
    buf[0x404] = 1;
    buf[0x405] = 0;
    buf[0x406..0x408].copy_from_slice(&1_u16.to_le_bytes());
    buf[0x418..0x41c].copy_from_slice(&4_u32.to_le_bytes());
    buf[0x428..0x430].copy_from_slice(&0x500_u64.to_le_bytes());
    buf[0x500..0x504].copy_from_slice(b"NOPE");

    let mut cursor = std::io::Cursor::new(&buf);
    let out = read_hdf5_v0_dataset_bytes_seek(&mut cursor, 0x20).expect("no io error");
    assert_eq!(
        out, None,
        "illegal non-decreasing B-tree level or truncated node must fail closed with None"
    );

    // Valid sibling leaf + malformed sibling must ALSO fail closed with None (not partial bytes).
    buf[0x128..0x130].copy_from_slice(&0x400_u64.to_le_bytes()); // valid sibling first
    buf[0x140..0x148].copy_from_slice(&0x200_u64.to_le_bytes()); // illegal level-1 sibling second
    let mut cursor2 = std::io::Cursor::new(&buf);
    let out2 = read_hdf5_v0_dataset_bytes_seek(&mut cursor2, 0x20).expect("no io error");
    assert_eq!(
        out2, None,
        "valid sibling followed by malformed sibling must fail closed, not return partial bytes"
    );
}

#[test]
fn streaming_progressive_hal_single_query_returns_normalized_bottom_count() {
    // Reviewer 3 #6 regression: queries.len() < 2 must return the post-normalization bottom count
    // actually written to disk (3 segments for interior [2..6) in 8 bp ref, not raw 1).
    let ref_contigs = vec![("c1".to_owned(), b"ACGTACGT".to_vec())];
    let queries = vec![(
        "q0".to_owned(),
        vec![("qc1".to_owned(), b"ACGTACGT".to_vec())],
    )];
    let blocks = vec![(0, 0, 0, 2, 2, 4, false)];
    let tmp_path = std::env::temp_dir().join(format!(
        "horizon_single_q_nb_{}_{:?}.hal",
        std::process::id(),
        std::thread::current().id()
    ));
    let (_, nb) = HalAlignment::write_multi_contig_progressive_hdf5_hal(
        &tmp_path,
        "hg38",
        &ref_contigs,
        &queries,
        &blocks,
    )
    .expect("write single-query hal");
    let bytes = std::fs::read(&tmp_path).expect("read");
    let _ = std::fs::remove_file(&tmp_path);
    let bots = decode_hal_bottom_array(&read_hal_array_bytes(&bytes, "hg38", "BOTTOM_ARRAY"), 1);
    // bots includes +1 libHal sentinel -> active bottom segments = bots.len() - 1 == 3.
    assert_eq!(nb, bots.len() - 1);
    assert_eq!(nb, 3);
}

#[test]
fn hdf5_read_handles_haddr_undef_on_file_consistently_with_cursor() {
    // Build a minimal v1 dataset header whose DataLayout (0x0008, v3 contiguous) has address = u64::MAX (HADDR_UNDEF).
    let mut buf = vec![0u8; 0x80];
    buf[0] = 1; // v1
    buf[2..4].copy_from_slice(&2_u16.to_le_bytes()); // 2 messages
    buf[8..12].copy_from_slice(&48_u32.to_le_bytes());
    let mut p = 16usize;
    // Dataspace (1D, len=4)
    buf[p..p + 2].copy_from_slice(&0x0001_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&16_u16.to_le_bytes());
    buf[p + 8] = 1;
    buf[p + 9] = 1;
    buf[p + 16..p + 24].copy_from_slice(&4_u64.to_le_bytes());
    p += 24;
    // DataLayout v3 contiguous at u64::MAX, len=16
    buf[p..p + 2].copy_from_slice(&0x0008_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&24_u16.to_le_bytes());
    buf[p + 8] = 3;
    buf[p + 9] = 1;
    buf[p + 10..p + 18].copy_from_slice(&u64::MAX.to_le_bytes());
    buf[p + 18..p + 26].copy_from_slice(&16_u64.to_le_bytes());

    let tmp_path = std::env::temp_dir().join(format!(
        "horizon_haddr_undef_{}_{:?}.hal",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::write(&tmp_path, &buf).expect("write");
    let mut file = std::fs::File::open(&tmp_path).expect("open");
    let mut cursor = std::io::Cursor::new(&buf);
    let from_cursor = read_hdf5_v0_dataset_bytes_seek(&mut cursor, 0).expect("cursor ok");
    let from_file = read_hdf5_v0_dataset_bytes_seek(&mut file, 0).expect("file must not EINVAL");
    let _ = std::fs::remove_file(&tmp_path);
    assert_eq!(from_file, from_cursor);
}

#[test]
fn hdf5_write_rejects_duplicate_genome_names() {
    let hal = HalAlignment {
        tree_newick: "(dup)dup;".to_owned(),
        genomes: vec![
            HalGenome {
                name: "dup".to_owned(),
                dna: b"ACGT".to_vec(),
                sequences: vec![HalSequence {
                    name: "c1".to_owned(),
                    offset: 0,
                    length: 4,
                }],
                top_segments: Vec::new(),
                bottom_segments: Vec::new(),
                parent: None,
                children: vec![1],
            },
            HalGenome {
                name: "dup".to_owned(),
                dna: b"ACGT".to_vec(),
                sequences: vec![HalSequence {
                    name: "c1".to_owned(),
                    offset: 0,
                    length: 4,
                }],
                top_segments: Vec::new(),
                bottom_segments: Vec::new(),
                parent: Some(0),
                children: Vec::new(),
            },
        ],
    };
    let tmp_path = std::env::temp_dir().join(format!(
        "horizon_dup_genome_{}_{:?}.hal",
        std::process::id(),
        std::thread::current().id()
    ));
    let err = hal
        .write_hdf5_hal(&tmp_path)
        .expect_err("duplicate genome names must be rejected");
    let _ = std::fs::remove_file(&tmp_path);
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
}

#[test]
fn read_hdf5_object_header_at_resolves_v1_continuation_at_nonzero_offset() {
    // v1 object header placed at offset 0x20 with a 0x0010 continuation pointing to file offset 0x80.
    let mut buf = vec![0u8; 0x100];
    let hdr_off = 0x20usize;
    buf[hdr_off] = 1; // v1
    buf[hdr_off + 2..hdr_off + 4].copy_from_slice(&2_u16.to_le_bytes()); // 2 total msgs (1 cont + 1 target)
    buf[hdr_off + 8..hdr_off + 12].copy_from_slice(&24_u32.to_le_bytes());
    let p = hdr_off + 16;
    buf[p..p + 2].copy_from_slice(&0x0010_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&16_u16.to_le_bytes());
    buf[p + 8..p + 16].copy_from_slice(&0x80_u64.to_le_bytes());
    buf[p + 16..p + 24].copy_from_slice(&16_u64.to_le_bytes());
    // Continuation chunk at 0x80 contains message type 0x0008.
    buf[0x80..0x82].copy_from_slice(&0x0008_u16.to_le_bytes());
    buf[0x82..0x84].copy_from_slice(&8_u16.to_le_bytes());

    let parsed = read_hdf5_object_header_at(&buf, hdr_off as u64).expect("parse at 0x20");
    assert_eq!(parsed.version, 1);
    assert_eq!(parsed.message_types, vec![0x0008]);

    // Also confirm a v1 continuation pointing to HADDR_UNDEF (u64::MAX) does not raise EINVAL on File.
    buf[p + 8..p + 16].copy_from_slice(&u64::MAX.to_le_bytes());
    let tmp_path = std::env::temp_dir().join(format!(
        "horizon_v1_cont_undef_{}_{:?}.hal",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::write(&tmp_path, &buf).expect("write");
    let mut file = std::fs::File::open(&tmp_path).expect("open");
    let mut cursor = std::io::Cursor::new(&buf);
    let from_file = read_hdf5_v0_dataset_bytes_seek(&mut file, hdr_off as u64).expect("no EINVAL");
    let from_cursor =
        read_hdf5_v0_dataset_bytes_seek(&mut cursor, hdr_off as u64).expect("cursor ok");
    let _ = std::fs::remove_file(&tmp_path);
    assert_eq!(from_file, from_cursor);
}

#[test]
fn windowed_progressive_fitch_matches_in_memory_across_straddling_and_odd_windows() {
    // Odd total length (17 bp across two contigs: 9 + 8) with window_bp = 4 forces 5 windows
    // that split forward and reverse-complement segments AND contig boundaries mid-segment.
    let ref_contigs = vec![
        ("c1".to_owned(), b"ACGTACGTA".to_vec()),
        ("c2".to_owned(), b"TTTTCCCC".to_vec()),
    ];
    let queries = vec![
        (
            "sp1".to_owned(),
            vec![
                ("q1a".to_owned(), b"GCGTACGTG".to_vec()),
                ("q1b".to_owned(), b"GGGGAAAA".to_vec()),
            ],
        ),
        (
            "sp2".to_owned(),
            vec![
                ("q2a".to_owned(), b"GCGTACGTT".to_vec()),
                ("q2b".to_owned(), b"GGGGAAAA".to_vec()),
            ],
        ),
        (
            "sp3".to_owned(),
            vec![
                ("q3a".to_owned(), b"TCGTACGTT".to_vec()),
                ("q3b".to_owned(), b"ATTTCCCC".to_vec()),
            ],
        ),
    ];
    let blocks = vec![
        (0, 0, 0, 0, 0, 9, false),
        (0, 1, 1, 0, 0, 8, true),
        (1, 0, 0, 0, 0, 9, false),
        (1, 1, 1, 0, 0, 8, true),
        (2, 0, 0, 0, 0, 9, false),
        (2, 1, 1, 0, 0, 8, false),
    ];
    let in_mem =
        HalAlignment::from_multi_contig_progressive_blocks("hg38", &ref_contigs, &queries, &blocks);
    let tmp_path = std::env::temp_dir().join(format!(
        "horizon_windowed_fitch_{}_{:?}.hal",
        std::process::id(),
        std::thread::current().id()
    ));
    let (tree, nb) = HalAlignment::write_multi_contig_progressive_hdf5_hal_with_window(
        &tmp_path,
        "hg38",
        &ref_contigs,
        &queries,
        &blocks,
        4,
    )
    .expect("windowed streaming HAL");
    let loaded_bytes = std::fs::read(&tmp_path).expect("read hal");
    let _ = std::fs::remove_file(&tmp_path);

    assert_eq!(tree, in_mem.tree_newick);
    assert_eq!(nb, in_mem.genomes[0].bottom_segments.len());

    // Decode every ancestor and leaf DNA_ARRAY and compare byte-for-byte with in-memory Fitch.
    let in_mem_bytes = in_mem.to_hdf5_bytes();
    for g in &in_mem.genomes {
        let dna_raw = read_hal_array_bytes(&loaded_bytes, &g.name, "DNA_ARRAY");
        let in_mem_raw = read_hal_array_bytes(&in_mem_bytes, &g.name, "DNA_ARRAY");
        assert_eq!(
            dna_raw, in_mem_raw,
            "raw packed DNA_ARRAY bytes mismatch for genome {}",
            g.name
        );
        let decoded = decode_hal_dna(&dna_raw, g.dna.len());
        assert_eq!(
            decoded, g.dna,
            "windowed Fitch DNA mismatch for genome {}",
            g.name
        );
    }
}

#[test]
fn windowed_progressive_fitch_handles_empty_contigs_empty_ref_and_odd_window() {
    // Interleaved empty + odd-length contigs with odd window_bp = 3 (clamped to even 2).
    let ref_contigs = vec![
        ("empty0".to_owned(), Vec::new()),
        ("c1".to_owned(), b"ACGTA".to_vec()),
        ("empty1".to_owned(), Vec::new()),
        ("c2".to_owned(), b"TTCC".to_vec()),
    ];
    let queries = vec![
        (
            "q1".to_owned(),
            vec![
                ("e0".to_owned(), Vec::new()),
                ("qc1".to_owned(), b"GCGTA".to_vec()),
                ("qc2".to_owned(), b"GGAA".to_vec()),
            ],
        ),
        (
            "q2".to_owned(),
            vec![
                ("qc1".to_owned(), b"GCGTT".to_vec()),
                ("e1".to_owned(), Vec::new()),
                ("qc2".to_owned(), b"GGAA".to_vec()),
            ],
        ),
    ];
    let blocks = vec![
        (0, 1, 1, 0, 0, 5, false),
        (0, 3, 2, 0, 0, 4, true),
        (1, 1, 0, 0, 0, 5, false),
        (1, 3, 2, 0, 0, 4, true),
    ];
    let in_mem =
        HalAlignment::from_multi_contig_progressive_blocks("ref", &ref_contigs, &queries, &blocks);
    let in_mem_bytes = in_mem.to_hdf5_bytes();
    let tmp_path = std::env::temp_dir().join(format!(
        "horizon_windowed_odd_empty_{}_{:?}.hal",
        std::process::id(),
        std::thread::current().id()
    ));
    let (tree, nb) = HalAlignment::write_multi_contig_progressive_hdf5_hal_with_window(
        &tmp_path,
        "ref",
        &ref_contigs,
        &queries,
        &blocks,
        3,
    )
    .expect("odd-window empty-contig HAL");
    let loaded = std::fs::read(&tmp_path).expect("read");
    let _ = std::fs::remove_file(&tmp_path);
    assert_eq!(tree, in_mem.tree_newick);
    assert_eq!(nb, in_mem.genomes[0].bottom_segments.len());
    for g in &in_mem.genomes {
        assert_eq!(
            read_hal_array_bytes(&loaded, &g.name, "DNA_ARRAY"),
            read_hal_array_bytes(&in_mem_bytes, &g.name, "DNA_ARRAY"),
            "raw DNA_ARRAY mismatch for {}",
            g.name
        );
    }

    // Zero-length reference with >= 2 queries must succeed cleanly.
    let tmp_empty = std::env::temp_dir().join(format!(
        "horizon_windowed_zero_ref_{}_{:?}.hal",
        std::process::id(),
        std::thread::current().id()
    ));
    let empty_ref = vec![("empty".to_owned(), Vec::new())];
    let (_, nb_empty) = HalAlignment::write_multi_contig_progressive_hdf5_hal_with_window(
        &tmp_empty,
        "ref",
        &empty_ref,
        &queries,
        &[],
        4,
    )
    .expect("zero-length ref HAL");
    let _ = std::fs::remove_file(&tmp_empty);
    assert_eq!(nb_empty, 0);
}
