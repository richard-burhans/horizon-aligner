use super::*;

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
    // - FilterPipeline v1 (0x000b): deflate (id=1)
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
    let hg38 = root_entries
        .iter()
        .find(|e| e.name == "hg38")
        .expect("hg38");
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
    assert_eq!(
        decode_hal_seq_names(&sn_raw, sn_stride),
        vec!["chr1", "chr2"]
    );
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
    for req in [
        "Meta",
        "Phylogeny",
        "Genomes",
        "Verison",
        "hg38",
        "q1",
        "q2",
    ] {
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
    assert_eq!(hal.genomes[1].top_segments[0].next_paralog, Some(1));
    hal.genomes[1].top_segments[1].next_paralog = Some(0);
    let maf = hal.to_maf_for_genome("hg38").expect("maf");
    let q1_rows = maf.lines().filter(|l| l.starts_with("s q1 ")).count();
    assert_eq!(
        q1_rows, 2,
        "cyclic next_paralog must emit each paralog once:\n{maf}"
    );
}

pub(super) fn read_hal_array_bytes(bytes: &[u8], genome: &str, dataset: &str) -> Vec<u8> {
    let sb = read_hdf5_superblock(bytes).expect("superblock");
    let root_addr = sb.root_group_addr.expect("root addr");
    let mut cursor = std::io::Cursor::new(bytes);
    let root_entries = read_hdf5_v0_group_entries_seek(&mut cursor, root_addr)
        .expect("root entries")
        .expect("some entries");
    let grp = root_entries
        .iter()
        .find(|e| e.name == genome)
        .expect("genome group");
    let ds_entries = read_hdf5_v0_group_entries_seek(&mut cursor, grp.object_header_addr)
        .expect("genome ds")
        .expect("some ds");
    let ds_obj = ds_entries
        .iter()
        .find(|e| e.name == dataset)
        .expect("dataset");
    read_hdf5_v0_dataset_bytes_seek(&mut cursor, ds_obj.object_header_addr)
        .expect("ds read")
        .expect("ds bytes")
}

#[test]
#[allow(clippy::too_many_lines)]
fn non_star_unsorted_segments_remap_cross_references_on_hdf5_export() {
    // Genuine 3-level non-star tree: ((leaf1)anc,leaf2)root;
    // All segment arrays are deliberately unsorted AND leaf1 carries an unsorted paralogy cycle
    // so `is_non_star_tree` is true and `build_sort_permutation` + cross-reference remapping run.
    let hal = HalAlignment {
        tree_newick: "((leaf1)anc,leaf2)root;".to_owned(),
        genomes: vec![
            HalGenome {
                name: "root".to_owned(),
                dna: b"ACGTACGT".to_vec(),
                sequences: vec![HalSequence {
                    name: "chr1".to_owned(),
                    offset: 0,
                    length: 8,
                }],
                top_segments: Vec::new(),
                bottom_segments: vec![
                    HalBottomSegment {
                        start: 4,
                        length: 4,
                        child_tops: vec![Some((1, false)), None],
                    },
                    HalBottomSegment {
                        start: 0,
                        length: 4,
                        child_tops: vec![Some((0, false)), None],
                    },
                ],
                parent: None,
                children: vec![1, 3],
            },
            HalGenome {
                name: "anc".to_owned(),
                dna: b"ACGTACGT".to_vec(),
                sequences: vec![HalSequence {
                    name: "chr1".to_owned(),
                    offset: 0,
                    length: 8,
                }],
                top_segments: vec![
                    HalTopSegment {
                        start: 4,
                        length: 4,
                        parent_bottom_idx: Some(1),
                        parent_reversed: false,
                        next_paralog: None,
                    },
                    HalTopSegment {
                        start: 0,
                        length: 4,
                        parent_bottom_idx: Some(0),
                        parent_reversed: false,
                        next_paralog: None,
                    },
                ],
                bottom_segments: vec![
                    HalBottomSegment {
                        start: 4,
                        length: 4,
                        child_tops: vec![Some((1, false))],
                    },
                    HalBottomSegment {
                        start: 0,
                        length: 4,
                        child_tops: vec![Some((0, false))],
                    },
                ],
                parent: Some(0),
                children: vec![2],
            },
            HalGenome {
                name: "leaf1".to_owned(),
                dna: b"ACGTACGTACGT".to_vec(),
                sequences: vec![HalSequence {
                    name: "chr1".to_owned(),
                    offset: 0,
                    length: 12,
                }],
                top_segments: vec![
                    // Old idx 0 (start=8) <-> old idx 2 (start=0) form a 2-node paralogy cycle under anc old bot[1] (start=0).
                    HalTopSegment {
                        start: 8,
                        length: 4,
                        parent_bottom_idx: Some(1),
                        parent_reversed: false,
                        next_paralog: Some(2),
                    },
                    // Old idx 1 (start=4) -> anc old bot[0] (start=4).
                    HalTopSegment {
                        start: 4,
                        length: 4,
                        parent_bottom_idx: Some(0),
                        parent_reversed: false,
                        next_paralog: None,
                    },
                    HalTopSegment {
                        start: 0,
                        length: 4,
                        parent_bottom_idx: Some(1),
                        parent_reversed: false,
                        next_paralog: Some(0),
                    },
                ],
                bottom_segments: Vec::new(),
                parent: Some(1),
                children: Vec::new(),
            },
            HalGenome {
                name: "leaf2".to_owned(),
                dna: b"AAAA".to_vec(),
                sequences: vec![HalSequence {
                    name: "chr1".to_owned(),
                    offset: 0,
                    length: 4,
                }],
                top_segments: vec![HalTopSegment {
                    start: 0,
                    length: 4,
                    parent_bottom_idx: None,
                    parent_reversed: false,
                    next_paralog: None,
                }],
                bottom_segments: Vec::new(),
                parent: Some(0),
                children: Vec::new(),
            },
        ],
    };
    let bytes = hal.to_hdf5_bytes();
    let root_bots =
        decode_hal_bottom_array(&read_hal_array_bytes(&bytes, "root", "BOTTOM_ARRAY"), 2);
    let anc_tops = decode_hal_top_array(&read_hal_array_bytes(&bytes, "anc", "TOP_ARRAY"));
    let anc_bots = decode_hal_bottom_array(&read_hal_array_bytes(&bytes, "anc", "BOTTOM_ARRAY"), 1);
    let leaf1_tops = decode_hal_top_array(&read_hal_array_bytes(&bytes, "leaf1", "TOP_ARRAY"));

    // root bot[0] (start=0, old bot[1]) -> anc old top[0] (start=4) -> new anc top[1].
    assert_eq!(root_bots[0].0, 0);
    assert_eq!(root_bots[0].2, -1, "root bottom topIdx must be NULL_INDEX");
    assert_eq!(root_bots[0].3[0], (1, false));
    // anc top[0] (start=0, old top[1]) -> root old bot[0] (start=4) -> new root bot[1].
    assert_eq!(anc_tops[0].0, 0);
    assert_eq!(
        anc_tops[0].1, 0,
        "internal anc top[0] must link to bottomIdx=0"
    );
    assert_eq!(anc_tops[0].3, 1);
    assert_eq!(
        anc_tops[1].1, 1,
        "internal anc top[1] must link to bottomIdx=1"
    );
    // anc bot[0] (start=0, old bot[1]) -> leaf1 old top[0] (start=8) -> new leaf1 top[2].
    assert_eq!(anc_bots[0].0, 0);
    assert_eq!(
        anc_bots[0].2, 0,
        "internal anc bot[0] must link to topIdx=0"
    );
    assert_eq!(anc_bots[0].3[0], (2, false));
    assert_eq!(
        anc_bots[1].2, 1,
        "internal anc bot[1] must link to topIdx=1"
    );
    // leaf1 sorted by start: [0 (old 2), 4 (old 1), 8 (old 0)].
    // Paralogy cycle old 0 <-> old 2 remaps to new 2 <-> new 0; parent old bot[1] remaps to new bot[0].
    assert_eq!(leaf1_tops[0].0, 0);
    assert_eq!(leaf1_tops[0].1, -1, "leaf bottomIdx must be NULL_INDEX");
    assert_eq!(leaf1_tops[0].3, 0);
    assert_eq!(leaf1_tops[0].2, 2);
    assert_eq!(leaf1_tops[1].0, 4);
    assert_eq!(leaf1_tops[1].3, 1);
    assert_eq!(leaf1_tops[1].2, -1);
    assert_eq!(leaf1_tops[2].0, 8);
    assert_eq!(leaf1_tops[2].3, 0);
    assert_eq!(leaf1_tops[2].2, 0);
}
