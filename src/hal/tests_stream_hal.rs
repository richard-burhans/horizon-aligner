use super::*;
use std::io::Cursor;

#[test]
#[allow(clippy::too_many_lines)]
fn streaming_hal_multilevel_progressive_stats_liftover_and_maf() {
    // Build a genuine 3-level progressive HAL: ((hg38,mm10)Anc0,rn6)Anc1;
    // with a non-palindromic, odd-length chr2 (15 bp) so reverse-complement DNA, partial-slice
    // UCSC `-`-strand MAF coordinates, odd global_start DNA decoding, and total_bp clamping are
    // all strictly verified (none can be masked by palindromic symmetry!).
    let hg_c1 = b"ACGTACGTACGTACGT".to_vec(); // 16 bp (even)
    let hg_c2 = b"AAACCCGGGTTTACG".to_vec(); // 15 bp (odd, non-palindrome; revcomp = CGTAAACCCGGGTTT)
    let mm_c1 = b"ACGTACGTACGTACGT".to_vec();
    let mm_c2 = crate::seeding::revcomp(&hg_c2);
    let rn_c1 = b"ACGTACGTACGTACGT".to_vec();
    let rn_c2 = hg_c2.clone();

    let dir = std::env::temp_dir().join(format!(
        "horizon_stream_hal_test_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).expect("mkdir");
    let hal_path = dir.join("prog.hal");
    let (tree, _nb) = HalAlignment::write_multi_contig_progressive_hdf5_hal(
        &hal_path,
        "hg38",
        &[("chr1".to_owned(), hg_c1), ("chr2".to_owned(), hg_c2)],
        &[
            (
                "mm10".to_owned(),
                vec![("m_chr1".to_owned(), mm_c1), ("m_chr2".to_owned(), mm_c2)],
            ),
            (
                "rn6".to_owned(),
                vec![("r_chr1".to_owned(), rn_c1), ("r_chr2".to_owned(), rn_c2)],
            ),
        ],
        &[
            (0, 0, 0, 0, 0, 16, false),
            (0, 1, 1, 0, 0, 15, true),
            (1, 0, 0, 0, 0, 16, false),
            (1, 1, 1, 0, 0, 15, false),
        ],
    )
    .expect("write progressive HAL");
    assert_eq!(tree, "((hg38,mm10)Anc0,rn6)Anc1;");

    let mut hal = StreamingHalFile::open(&hal_path).expect("open streaming HAL");
    assert_eq!(hal.tree_newick(), "((hg38,mm10)Anc0,rn6)Anc1;");
    // Phase 5A: open() must be lazy — no genome SEQIDX/SEQNAME loaded yet.
    assert!(!hal.is_genome_sequences_loaded("hg38"));
    assert!(!hal.is_genome_offsets_loaded("Anc0"));

    let stats = hal.stats().expect("stats");
    assert_eq!(stats.len(), 5, "Anc1, Anc0, hg38, mm10, rn6");
    // stats() resolves total_bp from SEQIDX sentinel only, without loading SEQNAME tables.
    assert!(!hal.is_genome_sequences_loaded("hg38"));
    assert!(!hal.is_genome_sequences_loaded("mm10"));
    let hg_st = stats.iter().find(|s| s.name == "hg38").expect("hg38");
    assert_eq!(hg_st.total_bp, 31);
    assert_eq!(hg_st.num_sequences, 2);
    assert_eq!(hg_st.parent.as_deref(), Some("Anc0"));

    // Multi-level liftover: hg38 -> Anc0 -> Anc1 -> rn6 (3 tree hops!)
    let hits_rn = hal
        .liftover_interval("hg38", "chr2", 2, 12, "rn6")
        .expect("liftover hg38 -> rn6");
    // Phase 5A: untouched leaf `mm10` never loads SEQIDX/SEQNAME; intermediate ancestors load
    // compact SEQIDX offsets only without SEQNAME strings!
    assert!(hal.is_genome_sequences_loaded("hg38"));
    assert!(hal.is_genome_sequences_loaded("rn6"));
    assert!(!hal.is_genome_sequences_loaded("mm10"));
    assert!(!hal.is_genome_offsets_loaded("mm10"));
    assert!(hal.is_genome_offsets_loaded("Anc0"));
    assert!(!hal.is_genome_sequences_loaded("Anc0"));
    assert_eq!(hits_rn.len(), 1);
    assert_eq!(hits_rn[0].seq_name, "r_chr2");
    assert_eq!(hits_rn[0].start, 2);
    assert_eq!(hits_rn[0].end, 12);
    assert_eq!(hits_rn[0].strand, '+');

    // Multi-level reverse-complement liftover: rn6 -> Anc1 -> Anc0 -> mm10 (3 hops + rev!)
    // On 15-bp chr2, [2..12) on + maps to [15-12 .. 15-2) = [3..13) on -.
    let hits_mm = hal
        .liftover_interval("rn6", "r_chr2", 2, 12, "mm10")
        .expect("liftover rn6 -> mm10");
    assert_eq!(hits_mm.len(), 1);
    assert_eq!(hits_mm[0].seq_name, "m_chr2");
    assert_eq!(hits_mm[0].start, 3);
    assert_eq!(hits_mm[0].end, 13);
    assert_eq!(hits_mm[0].strand, '-');

    // Partial-slice MAF projection on odd-offset non-palindromic chr2[1..10) (9 bp: "AACCCGGGT"):
    // On mm10 (- strand, len 15), local + coords are [5..14), so UCSC `-` start is 15 - 14 = 1!
    let mut slice_maf = Vec::new();
    let slice_stats = hal
        .project_to_maf("hg38", Some(("chr2", 1, Some(9))), None, &mut slice_maf)
        .expect("partial slice project_to_maf");
    assert_eq!(slice_stats.blocks_emitted, 1);
    assert_eq!(slice_stats.ref_bp_emitted, 9);
    let slice_text = String::from_utf8(slice_maf).expect("utf8");
    assert!(
        slice_text.contains("s hg38.chr2 1 9 + 15 AACCCGGGT"),
        "unexpected ref row:\n{slice_text}"
    );
    assert!(
        slice_text.contains("s mm10.m_chr2 1 9 - 15 AACCCGGGT"),
        "expected exact non-palindromic `-` MAF coordinate & revcomp DNA:\n{slice_text}"
    );
    assert!(
        slice_text.contains("s rn6.r_chr2 1 9 + 15 AACCCGGGT"),
        "unexpected rn6 row:\n{slice_text}"
    );

    // Odd total_bp clamping + odd global_start DNA read (must NOT leak padding nibble 'a').
    let hg_idx = hal
        .stats()
        .expect("stats")
        .iter()
        .position(|s| s.name == "hg38")
        .expect("hg38 idx");
    let tail_dna = hal.read_dna_slice(hg_idx, 29, 100).expect("clamped tail");
    assert_eq!(
        tail_dna, b"CG",
        "must clamp to total_bp=31 without padding nibble"
    );
    assert!(hal.read_dna_slice(999, 0, 4).is_err());

    // Verify in-memory Cursor reader works identically and tunable cache cap works cleanly.
    let bytes = std::fs::read(&hal_path).expect("read");
    let mut hal_mem = StreamingHalFile::open_reader(Cursor::new(bytes)).expect("open cursor");
    hal_mem.set_chunk_cache_cap(16);
    assert_eq!(hal_mem.chunk_cache_cap(), 16);
    assert_eq!(hal_mem.stats().expect("stats mem"), stats);
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn streaming_hal_newick_parser_rejects_deep_nesting_and_trailing_garbage() {
    let deep = format!("{}a{}", "(".repeat(5_000), ")".repeat(5_000));
    assert!(super::stream_hal_maf::parse_hal_newick(&deep).is_err());
    assert!(super::stream_hal_maf::parse_hal_newick("(a,b)root; GARBAGE").is_err());
    assert!(super::stream_hal_maf::parse_hal_newick("(a,b)root;").is_ok());
}

#[test]
#[allow(clippy::too_many_lines)]
fn hdf5_stream_chunked_deflated_multilevel_btree_and_lru_eviction() {
    use flate2::{Compression, write::ZlibEncoder};
    use std::io::Write as _;

    // Construct a synthetic HDF5 v0 chunked + zlib-deflated 2-level B-tree dataset:
    // 12 chunks * 4 bytes/chunk = 48 elements (elem_size=1), split across 2 leaf B-tree nodes
    // (6 chunks each) under a level-1 root B-tree node. Exceeds CHUNK_CACHE_CAP=8 to exercise LRU.
    const NUM_CHUNKS: usize = 12;
    const CHUNK_LEN: usize = 4;
    let payload: Vec<u8> = (0..(NUM_CHUNKS * CHUNK_LEN))
        .map(|i| b"ACGT"[i % 4] + u8::try_from((i / 4) % 3).expect("u8"))
        .collect();

    let mut compressed_chunks: Vec<Vec<u8>> = Vec::with_capacity(NUM_CHUNKS);
    for c in 0..NUM_CHUNKS {
        let mut enc = ZlibEncoder::new(Vec::new(), Compression::default());
        enc.write_all(&payload[c * CHUNK_LEN..(c + 1) * CHUNK_LEN])
            .expect("write");
        compressed_chunks.push(enc.finish().expect("zlib"));
    }

    let mut buf = vec![0u8; 0x2000];
    // Object header v1 at 0x20: 4 messages (Dataspace v2, Datatype 1B, DataLayout v3 chunked, FilterPipeline v1 deflate)
    buf[0x20] = 1;
    buf[0x22..0x24].copy_from_slice(&4_u16.to_le_bytes());

    let mut p = 0x30usize;
    // Msg 1: Dataspace v2 (0x0001, size 12): ndims=1, dim0=48
    buf[p..p + 2].copy_from_slice(&0x0001_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&12_u16.to_le_bytes());
    p += 8;
    buf[p] = 2;
    buf[p + 1] = 1;
    buf[p + 2] = 0;
    buf[p + 3] = 1;
    buf[p + 4..p + 12].copy_from_slice(&((NUM_CHUNKS * CHUNK_LEN) as u64).to_le_bytes());
    p += 12;

    // Msg 2: Datatype (0x0003, size 8): fixed-point 1 byte
    buf[p..p + 2].copy_from_slice(&0x0003_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&8_u16.to_le_bytes());
    p += 8;
    buf[p] = 0x10; // class 0, version 1
    buf[p + 4..p + 8].copy_from_slice(&1_u32.to_le_bytes());
    p += 8;

    // Msg 3: DataLayout v3 chunked (0x0008, size 19): ndims+1 = 2, btree=0x200, chunk_dim=4, elem=1
    let root_btree: u64 = 0x200;
    buf[p..p + 2].copy_from_slice(&0x0008_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&19_u16.to_le_bytes());
    p += 8;
    buf[p] = 3;
    buf[p + 1] = 2;
    buf[p + 2] = 2;
    buf[p + 3..p + 11].copy_from_slice(&root_btree.to_le_bytes());
    buf[p + 11..p + 15].copy_from_slice(&u32::try_from(CHUNK_LEN).expect("u32").to_le_bytes());
    buf[p + 15..p + 19].copy_from_slice(&1_u32.to_le_bytes());
    p += 19;

    // Msg 4: FilterPipeline v1 (0x000b, size 24): deflate id=1
    buf[p..p + 2].copy_from_slice(&0x000b_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&24_u16.to_le_bytes());
    p += 8;
    buf[p] = 1;
    buf[p + 1] = 1;
    p += 8;
    buf[p..p + 2].copy_from_slice(&1_u16.to_le_bytes());
    buf[p + 6..p + 8].copy_from_slice(&1_u16.to_le_bytes());
    buf[p + 8..p + 12].copy_from_slice(&6_u32.to_le_bytes());
    p += 16;

    let obj_hdr_size = u32::try_from(p - 0x30).expect("u32");
    buf[0x28..0x2c].copy_from_slice(&obj_hdr_size.to_le_bytes());

    // Place chunk payloads starting at 0x1000
    let mut chunk_addrs = Vec::with_capacity(NUM_CHUNKS);
    let mut data_pos = 0x1000usize;
    for comp in &compressed_chunks {
        chunk_addrs.push(data_pos as u64);
        buf[data_pos..data_pos + comp.len()].copy_from_slice(comp);
        data_pos += (comp.len() + 7) & !7;
    }

    // ndim=2 -> key_bytes = 8 + 16 = 24, rec_bytes = 32
    let key_bytes = 24usize;
    let rec_bytes = 32usize;
    let leaf0: u64 = 0x400;
    let leaf1: u64 = 0x700;

    // Level-1 root B-tree at 0x200 with 2 children: leaf0 [0..24), leaf1 [24..48)
    let rb = 0x200usize;
    buf[rb..rb + 4].copy_from_slice(b"TREE");
    buf[rb + 4] = 1; // chunk B-tree
    buf[rb + 5] = 1; // level 1
    buf[rb + 6..rb + 8].copy_from_slice(&2_u16.to_le_bytes());
    // Entry 0: key offset_elem = 0, child = leaf0
    buf[rb + 24 + 8..rb + 24 + 16].copy_from_slice(&0_u64.to_le_bytes());
    buf[rb + 24 + key_bytes..rb + 24 + rec_bytes].copy_from_slice(&leaf0.to_le_bytes());
    // Entry 1: key offset_elem = 24, child = leaf1
    buf[rb + 24 + rec_bytes + 8..rb + 24 + rec_bytes + 16].copy_from_slice(&24_u64.to_le_bytes());
    buf[rb + 24 + rec_bytes + key_bytes..rb + 24 + 2 * rec_bytes]
        .copy_from_slice(&leaf1.to_le_bytes());
    // Trailing upper-bound key: offset_elem = 48
    buf[rb + 24 + 2 * rec_bytes + 8..rb + 24 + 2 * rec_bytes + 16]
        .copy_from_slice(&48_u64.to_le_bytes());

    // Populate leaf0 (chunks 0..6) and leaf1 (chunks 6..12)
    for (leaf_addr, range) in [(0x400usize, 0..6usize), (0x700usize, 6..12usize)] {
        buf[leaf_addr..leaf_addr + 4].copy_from_slice(b"TREE");
        buf[leaf_addr + 4] = 1;
        buf[leaf_addr + 5] = 0; // leaf
        buf[leaf_addr + 6..leaf_addr + 8].copy_from_slice(&6_u16.to_le_bytes());
        for (slot, c) in range.enumerate() {
            let off = leaf_addr + 24 + slot * rec_bytes;
            let comp_len = u32::try_from(compressed_chunks[c].len()).expect("u32");
            buf[off..off + 4].copy_from_slice(&comp_len.to_le_bytes());
            buf[off + 8..off + 16].copy_from_slice(&((c * CHUNK_LEN) as u64).to_le_bytes());
            buf[off + key_bytes..off + rec_bytes].copy_from_slice(&chunk_addrs[c].to_le_bytes());
        }
    }

    let mut cur = Cursor::new(&buf);
    let mut ds = Hdf5DatasetReader::open_with_cache_cap(&mut cur, 0x20, 8)
        .expect("open chunked reader")
        .expect("dataset present");
    assert_eq!(ds.total_elems(), 48);
    assert_eq!(ds.cache_cap(), 8);

    // Cross-chunk slice spanning leaf0 chunks [0..3) (start=2, count=8 -> [2..10))
    assert_eq!(
        ds.read_elem_range(&mut cur, 2, 8).expect("r1"),
        payload[2..10]
    );
    // Cross-leaf slice spanning all 12 chunks (start=10, count=36 -> [10..46), forces LRU eviction beyond cap=8)
    assert_eq!(
        ds.read_elem_range(&mut cur, 10, 36).expect("r2"),
        payload[10..46]
    );
    // Shrink cache cap dynamically to 2 and verify exact recovery from disk + two-chunk fast path
    ds.set_cache_cap(2);
    assert_eq!(ds.cache_cap(), 2);
    assert_eq!(
        ds.read_elem_range(&mut cur, 0, 6).expect("r3"),
        payload[0..6]
    );
    // Second read across the same two cached chunks [0..4) + [4..8) hits two-chunk cache fast path
    assert_eq!(
        ds.read_elem_range(&mut cur, 2, 5).expect("r4"),
        payload[2..7]
    );
}

#[test]
fn hdf5_stream_rejects_btree_cycle_nondecreasing_level_and_truncated_contiguous() {
    // 1. Self-referential B-tree cycle at 0x200 -> 0x200 must fail closed with InvalidData.
    let mut buf = vec![0u8; 0x400];
    buf[0x20] = 1;
    buf[0x22..0x24].copy_from_slice(&3_u16.to_le_bytes());
    let mut p = 0x30usize;
    buf[p..p + 2].copy_from_slice(&0x0001_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&12_u16.to_le_bytes());
    p += 8;
    buf[p] = 2;
    buf[p + 1] = 1;
    buf[p + 3] = 1;
    buf[p + 4..p + 12].copy_from_slice(&16_u64.to_le_bytes());
    p += 12;
    buf[p..p + 2].copy_from_slice(&0x0003_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&8_u16.to_le_bytes());
    p += 8;
    buf[p] = 0x10;
    buf[p + 4..p + 8].copy_from_slice(&1_u32.to_le_bytes());
    p += 8;
    buf[p..p + 2].copy_from_slice(&0x0008_u16.to_le_bytes());
    buf[p + 2..p + 4].copy_from_slice(&19_u16.to_le_bytes());
    p += 8;
    buf[p] = 3;
    buf[p + 1] = 2;
    buf[p + 2] = 2;
    buf[p + 3..p + 11].copy_from_slice(&0x200_u64.to_le_bytes());
    buf[p + 11..p + 15].copy_from_slice(&4_u32.to_le_bytes());
    buf[p + 15..p + 19].copy_from_slice(&1_u32.to_le_bytes());
    p += 19;
    buf[0x28..0x2c].copy_from_slice(&u32::try_from(p - 0x30).expect("u32").to_le_bytes());

    // B-tree node at 0x200 pointing back to itself (0x200)
    buf[0x200..0x204].copy_from_slice(b"TREE");
    buf[0x204] = 1;
    buf[0x205] = 1;
    buf[0x206..0x208].copy_from_slice(&1_u16.to_le_bytes());
    buf[0x218 + 24..0x218 + 32].copy_from_slice(&0x200_u64.to_le_bytes());
    buf[0x218 + 32 + 8..0x218 + 32 + 16].copy_from_slice(&16_u64.to_le_bytes());

    let mut cur = Cursor::new(&buf);
    let mut ds = Hdf5DatasetReader::open(&mut cur, 0x20)
        .expect("open")
        .expect("dataset");
    let err = ds
        .read_elem_range(&mut cur, 0, 8)
        .expect_err("self-cycle B-tree must error");
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);

    // 2. Non-decreasing level (level 1 -> child at 0x280 also level 1) must fail closed.
    buf[0x218 + 24..0x218 + 32].copy_from_slice(&0x280_u64.to_le_bytes());
    buf[0x280..0x284].copy_from_slice(b"TREE");
    buf[0x284] = 1;
    buf[0x285] = 1; // NOT strictly less than parent_level=1!
    buf[0x286..0x288].copy_from_slice(&1_u16.to_le_bytes());
    let mut cur2 = Cursor::new(&buf);
    let mut ds2 = Hdf5DatasetReader::open(&mut cur2, 0x20)
        .expect("open")
        .expect("dataset");
    let err2 = ds2
        .read_elem_range(&mut cur2, 0, 8)
        .expect_err("non-decreasing B-tree level must error");
    assert_eq!(err2.kind(), std::io::ErrorKind::InvalidData);

    // 3. Truncated contiguous dataset pointing beyond EOF must fail with UnexpectedEof.
    let mut cbuf = vec![0u8; 0x100];
    cbuf[0x20] = 1;
    cbuf[0x22..0x24].copy_from_slice(&3_u16.to_le_bytes());
    let mut cp = 0x30usize;
    cbuf[cp..cp + 2].copy_from_slice(&0x0001_u16.to_le_bytes());
    cbuf[cp + 2..cp + 4].copy_from_slice(&12_u16.to_le_bytes());
    cp += 8;
    cbuf[cp] = 2;
    cbuf[cp + 1] = 1;
    cbuf[cp + 3] = 1;
    cbuf[cp + 4..cp + 12].copy_from_slice(&64_u64.to_le_bytes());
    cp += 12;
    cbuf[cp..cp + 2].copy_from_slice(&0x0003_u16.to_le_bytes());
    cbuf[cp + 2..cp + 4].copy_from_slice(&8_u16.to_le_bytes());
    cp += 8;
    cbuf[cp] = 0x10;
    cbuf[cp + 4..cp + 8].copy_from_slice(&1_u32.to_le_bytes());
    cp += 8;
    cbuf[cp..cp + 2].copy_from_slice(&0x0008_u16.to_le_bytes());
    cbuf[cp + 2..cp + 4].copy_from_slice(&18_u16.to_le_bytes());
    cp += 8;
    cbuf[cp] = 3;
    cbuf[cp + 1] = 1; // contiguous
    cbuf[cp + 2..cp + 10].copy_from_slice(&0x0E0_u64.to_le_bytes()); // only 32 bytes left in 0x100 file!
    cbuf[cp + 10..cp + 18].copy_from_slice(&64_u64.to_le_bytes());
    cp += 18;
    cbuf[0x28..0x2c].copy_from_slice(&u32::try_from(cp - 0x30).expect("u32").to_le_bytes());
    let mut cur3 = Cursor::new(&cbuf);
    let mut ds3 = Hdf5DatasetReader::open(&mut cur3, 0x20)
        .expect("open")
        .expect("dataset");
    let err3 = ds3
        .read_elem_range(&mut cur3, 0, 64)
        .expect_err("truncated contiguous read must error");
    assert_eq!(err3.kind(), std::io::ErrorKind::UnexpectedEof);
}

#[test]
fn streaming_hal_multilevel_maf_splits_target_breakpoints() {
    // 3-level progressive HAL `((hg38,mm10)Anc0,rn6)Anc1;` where hg38 chr1 is 20 bp
    // aligned contiguously to rn6 [0..20), while mm10 splits into two halves on opposite strands
    // (chrA [0..10) '+' and chrB [0..10) '-'). Projecting hg38 -> MAF must split the ref segment
    // at bp 10 and emit BOTH elementary sub-blocks with exact coordinates and DNA!
    let hg_c1 = b"AAACCCGGGTTTTTACGACG".to_vec(); // 20 bp
    let mm_a = b"AAACCCGGGT".to_vec(); // matches hg[0..10) on '+'
    let mm_b = crate::seeding::revcomp(&hg_c1[10..20]); // matches hg[10..20) on '-'
    let rn_c1 = hg_c1.clone();

    let dir = std::env::temp_dir().join(format!(
        "horizon_stream_hal_split_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).expect("mkdir");
    let hal_path = dir.join("split.hal");
    let _ = HalAlignment::write_multi_contig_progressive_hdf5_hal(
        &hal_path,
        "hg38",
        &[("chr1".to_owned(), hg_c1)],
        &[
            (
                "mm10".to_owned(),
                vec![("chrA".to_owned(), mm_a), ("chrB".to_owned(), mm_b)],
            ),
            ("rn6".to_owned(), vec![("r_chr1".to_owned(), rn_c1)]),
        ],
        &[
            (0, 0, 0, 0, 0, 10, false),
            (0, 0, 1, 10, 0, 10, true),
            (1, 0, 0, 0, 0, 20, false),
        ],
    )
    .expect("write split HAL");

    let mut hal = StreamingHalFile::open(&hal_path).expect("open");
    let mut maf = Vec::new();
    let stats = hal
        .project_to_maf("hg38", Some(("chr1", 0, Some(20))), None, &mut maf)
        .expect("project_to_maf");
    let _ = std::fs::remove_dir_all(&dir);
    let text = String::from_utf8(maf).expect("utf8");
    assert_eq!(
        stats.blocks_emitted, 2,
        "expected 2 split MAF blocks:\n{text}"
    );
    assert_eq!(stats.ref_bp_emitted, 20);
    assert!(text.contains("s hg38.chr1 0 10 + 20 AAACCCGGGT"), "{text}");
    assert!(text.contains("s mm10.chrA 0 10 + 10 AAACCCGGGT"), "{text}");
    assert!(text.contains("s rn6.r_chr1 0 10 + 20 AAACCCGGGT"), "{text}");
    assert!(
        text.contains("s hg38.chr1 10 10 + 20 TTTACGACG")
            || text.contains("s hg38.chr1 10 10 + 20 TTTTACGACG"),
        "{text}"
    );
    assert!(text.contains("s mm10.chrB 0 10 - 10 TTTTACGACG"), "{text}");
    assert!(
        text.contains("s rn6.r_chr1 10 10 + 20 TTTTACGACG"),
        "{text}"
    );
}
