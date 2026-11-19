//! HAL-compatible hierarchical alignment model (`TopSegment`/`BottomSegment` + guide tree).
//!
//! Models the HAL graph schema (Hickey et al. 2013) in pure Rust:
//! - Phylogenetic Newick tree linking leaf and ancestral genomes
//! - Per-genome `TopSegment` (child->parent) and `BottomSegment` (parent->child) coordinate tables
//! - Reference-free projection (`to_maf_for_genome`) onto any genome in the tree
//! - HDF5 container magic-byte inspection (`is_hdf5_hal`) for external `.hal` files

use std::collections::BTreeMap;
use std::fmt::Write as _;

/// HDF5 8-byte file signature used by `.hal` containers.
pub const HDF5_SIGNATURE: [u8; 8] = [0x89, b'H', b'D', b'F', 0x0d, 0x0a, 0x1a, 0x0a];

/// One sequence (chromosome/contig) within a HAL genome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalSequence {
    /// Sequence name.
    pub name: String,
    /// Total length in bases.
    pub length: usize,
    /// Global start offset in the genome's concatenated coordinate space.
    pub offset: usize,
}

/// Child-to-parent colinear segment in HAL.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalTopSegment {
    /// Start coordinate in the child genome's global space.
    pub start: usize,
    /// Aligned span in bases.
    pub length: usize,
    /// Index of the matching `BottomSegment` in the parent genome (`None` if unaligned/insertion).
    pub parent_bottom_idx: Option<usize>,
    /// Whether the child aligns on the reverse strand relative to the parent.
    pub parent_reversed: bool,
    /// Optional index of next paralogous `TopSegment` in the same genome.
    pub next_paralog: Option<usize>,
}

/// Parent-to-child colinear segment in HAL.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalBottomSegment {
    /// Start coordinate in the parent genome's global space.
    pub start: usize,
    /// Aligned span in bases.
    pub length: usize,
    /// For each child genome index: `Some((top_seg_idx, reversed))` or `None` if deleted/absent.
    pub child_tops: Vec<Option<(usize, bool)>>,
}

/// One genome (leaf or ancestral) in a HAL alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalGenome {
    /// Species/genome name.
    pub name: String,
    /// Parent genome index (`None` for root).
    pub parent: Option<usize>,
    /// Child genome indices.
    pub children: Vec<usize>,
    /// Sequences belonging to this genome.
    pub sequences: Vec<HalSequence>,
    /// Concatenated DNA sequence (uppercase ASCII).
    pub dna: Vec<u8>,
    /// Child->parent segments.
    pub top_segments: Vec<HalTopSegment>,
    /// Parent->child segments.
    pub bottom_segments: Vec<HalBottomSegment>,
}

/// Complete hierarchical HAL alignment graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalAlignment {
    /// Newick guide tree.
    pub tree_newick: String,
    /// All genomes (root + internal + leaves).
    pub genomes: Vec<HalGenome>,
}

/// Check whether a byte slice begins with the standard HDF5 `.hal` container signature.
#[must_use]
pub fn is_hdf5_hal(bytes: &[u8]) -> bool {
    bytes.len() >= HDF5_SIGNATURE.len() && bytes[..8] == HDF5_SIGNATURE
}

impl HalAlignment {
    /// Construct a star-tree `HalAlignment` from a reference genome and pairwise `(species, seq_name, ref_start, qry_start, len, rev, ref_bases, qry_bases)` blocks.
    #[must_use]
    pub fn from_star_blocks(
        ref_name: &str,
        ref_dna: &[u8],
        queries: &[(String, Vec<u8>)],
        blocks: &[(usize, usize, usize, usize, bool)], // (query_idx, ref_start, qry_start, len, rev)
    ) -> Self {
        type HalIntervalHits = Vec<(usize, usize, bool)>;
        let mut genomes = Vec::with_capacity(1 + queries.len());
        let mut bottom_segments = Vec::new();

        // Group blocks by (ref_start, len) on the root/reference genome.
        let mut by_interval: BTreeMap<(usize, usize), HalIntervalHits> = BTreeMap::new();
        for &(q_idx, r_start, q_start, len, rev) in blocks {
            by_interval
                .entry((r_start, len))
                .or_default()
                .push((q_idx, q_start, rev));
        }

        let mut child_tops_per_q: Vec<Vec<HalTopSegment>> = vec![Vec::new(); queries.len()];
        for ((r_start, len), hits) in by_interval {
            let b_idx = bottom_segments.len();
            let mut child_links = vec![None; queries.len()];
            for (q_idx, q_start, rev) in hits {
                if q_idx >= queries.len() {
                    continue;
                }
                let t_idx = child_tops_per_q[q_idx].len();
                let prev_paralog = child_links[q_idx].map(|(prev_t, _)| prev_t);
                child_tops_per_q[q_idx].push(HalTopSegment {
                    start: q_start,
                    length: len,
                    parent_bottom_idx: Some(b_idx),
                    parent_reversed: rev,
                    next_paralog: prev_paralog,
                });
                child_links[q_idx] = Some((t_idx, rev));
            }
            bottom_segments.push(HalBottomSegment {
                start: r_start,
                length: len,
                child_tops: child_links,
            });
        }

        genomes.push(HalGenome {
            name: ref_name.to_owned(),
            parent: None,
            children: (1..=queries.len()).collect(),
            sequences: vec![HalSequence {
                name: format!("{ref_name}.chr1"),
                length: ref_dna.len(),
                offset: 0,
            }],
            dna: ref_dna.to_vec(),
            top_segments: Vec::new(),
            bottom_segments,
        });

        let mut leaf_names = Vec::with_capacity(queries.len());
        for (q_idx, (q_name, q_dna)) in queries.iter().enumerate() {
            leaf_names.push(q_name.clone());
            genomes.push(HalGenome {
                name: q_name.clone(),
                parent: Some(0),
                children: Vec::new(),
                sequences: vec![HalSequence {
                    name: format!("{q_name}.chr1"),
                    length: q_dna.len(),
                    offset: 0,
                }],
                dna: q_dna.clone(),
                top_segments: std::mem::take(&mut child_tops_per_q[q_idx]),
                bottom_segments: Vec::new(),
            });
        }

        let tree_newick = format!("({}){ref_name};", leaf_names.join(","));
        Self {
            tree_newick,
            genomes,
        }
    }

    fn write_leaf_maf(leaf: &HalGenome, parent: &HalGenome, out: &mut String) {
        for tseg in &leaf.top_segments {
            let Some(b_idx) = tseg.parent_bottom_idx else {
                continue;
            };
            let Some(bseg) = parent.bottom_segments.get(b_idx) else {
                continue;
            };
            let c_end = tseg.start.saturating_add(tseg.length).min(leaf.dna.len());
            let p_end = bseg.start.saturating_add(bseg.length).min(parent.dna.len());
            let aln_len = c_end
                .saturating_sub(tseg.start)
                .min(p_end.saturating_sub(bseg.start));
            if aln_len == 0 {
                continue;
            }
            let c_end = tseg.start + aln_len;
            let p_end = bseg.start + aln_len;
            let q_seq = String::from_utf8_lossy(&leaf.dna[tseg.start..c_end]);
            let p_slice = &parent.dna[bseg.start..p_end];
            let p_seq = if tseg.parent_reversed {
                String::from_utf8_lossy(&crate::revcomp(p_slice)).into_owned()
            } else {
                String::from_utf8_lossy(p_slice).into_owned()
            };
            let p_strand = if tseg.parent_reversed { '-' } else { '+' };
            let p_start = if tseg.parent_reversed {
                parent.dna.len().saturating_sub(p_end)
            } else {
                bseg.start
            };
            let _ = writeln!(out, "a score=2.0");
            let _ = writeln!(
                out,
                "s {} {} {} + {} {}",
                leaf.name,
                tseg.start,
                aln_len,
                leaf.dna.len(),
                q_seq
            );
            let _ = writeln!(
                out,
                "s {} {} {} {} {} {}\n",
                parent.name,
                p_start,
                aln_len,
                p_strand,
                parent.dna.len(),
                p_seq
            );
        }
    }

    /// Project the HAL alignment onto `target_genome_name` as MAF v1 text.
    #[must_use]
    pub fn to_maf_for_genome(&self, target_genome_name: &str) -> Option<String> {
        let g_idx = self
            .genomes
            .iter()
            .position(|g| g.name == target_genome_name)?;
        let root = &self.genomes[g_idx];
        let mut out = String::from("##maf version=1 scoring=horizon_hal\n");

        if root.bottom_segments.is_empty()
            && let Some(p_idx) = root.parent
            && let Some(parent) = self.genomes.get(p_idx)
        {
            Self::write_leaf_maf(root, parent, &mut out);
            return Some(out);
        }

        for bseg in &root.bottom_segments {
            let r_end = bseg.start.saturating_add(bseg.length).min(root.dna.len());
            if r_end <= bseg.start {
                continue;
            }
            let r_seq = String::from_utf8_lossy(&root.dna[bseg.start..r_end]);
            let active_children = bseg.child_tops.iter().filter(|c| c.is_some()).count();
            let _ = writeln!(out, "a score={}.0", active_children + 1);
            let _ = writeln!(
                out,
                "s {} {} {} + {} {}",
                root.name,
                bseg.start,
                r_end - bseg.start,
                root.dna.len(),
                r_seq
            );
            for (c_pos, link) in bseg.child_tops.iter().enumerate() {
                let Some((t_idx, rev)) = *link else {
                    continue;
                };
                let Some(&child_idx) = root.children.get(c_pos) else {
                    continue;
                };
                let Some(child) = self.genomes.get(child_idx) else {
                    continue;
                };
                let mut cur_t = Some(t_idx);
                let mut steps = 0usize;
                while let Some(idx) = cur_t {
                    if steps >= child.top_segments.len() {
                        break;
                    }
                    steps += 1;
                    let Some(tseg) = child.top_segments.get(idx) else {
                        break;
                    };
                    let c_end = tseg.start.saturating_add(tseg.length).min(child.dna.len());
                    if c_end > tseg.start {
                        let slice = &child.dna[tseg.start..c_end];
                        let q_seq = if rev {
                            String::from_utf8_lossy(&crate::revcomp(slice)).into_owned()
                        } else {
                            String::from_utf8_lossy(slice).into_owned()
                        };
                        let strand = if rev { '-' } else { '+' };
                        let q_start = if rev {
                            child.dna.len().saturating_sub(c_end)
                        } else {
                            tseg.start
                        };
                        let _ = writeln!(
                            out,
                            "s {} {} {} {} {} {}",
                            child.name,
                            q_start,
                            c_end - tseg.start,
                            strand,
                            child.dna.len(),
                            q_seq
                        );
                    }
                    cur_t = tseg.next_paralog;
                }
            }
            out.push('\n');
        }
        Some(out)
    }
}

mod hdf5_read;
mod hdf5_write;
pub use hdf5_read::*;
pub use hdf5_write::*;

/// Decode HAL `SEQNAME_ARRAY` fixed-width null-terminated strings.
#[must_use]
pub fn decode_hal_seq_names(raw: &[u8], stride: usize) -> Vec<String> {
    if stride == 0 {
        return Vec::new();
    }
    raw.chunks_exact(stride)
        .map(|slot| {
            let end = slot.iter().position(|&b| b == 0).unwrap_or(slot.len());
            String::from_utf8_lossy(&slot[..end]).into_owned()
        })
        .collect()
}

/// Decoded HAL `BOTTOM_ARRAY` record: `(start_idx, length, top_parse_idx, child_links)`.
pub type HalBottomRecord = (i64, u64, i64, Vec<(i64, bool)>);

/// Decode HAL `BOTTOM_ARRAY` packed records (`24 + 9 * num_children` bytes per segment).
#[must_use]
pub fn decode_hal_bottom_array(raw: &[u8], num_children: usize) -> Vec<HalBottomRecord> {
    let Some(rec_size) = num_children.checked_mul(9).and_then(|n| n.checked_add(24)) else {
        return Vec::new();
    };
    let mut out = Vec::with_capacity(raw.len() / rec_size);
    for rec in raw.chunks_exact(rec_size) {
        let start = i64::from_le_bytes(rec[0..8].try_into().unwrap_or([0; 8]));
        let length = u64::from_le_bytes(rec[8..16].try_into().unwrap_or([0; 8]));
        let top_parse = i64::from_le_bytes(rec[16..24].try_into().unwrap_or([0; 8]));
        let mut children = Vec::with_capacity(num_children);
        for c in 0..num_children {
            let off = 24 + c * 9;
            let child_idx = i64::from_le_bytes(rec[off..off + 8].try_into().unwrap_or([0; 8]));
            let rev = rec[off + 8] != 0;
            children.push((child_idx, rev));
        }
        out.push((start, length, top_parse, children));
    }
    out
}

/// Decode HAL `SEQIDX_ARRAY` 24-byte records `(start, top_offset, bottom_offset)`.
#[must_use]
pub fn decode_hal_seq_index(raw: &[u8]) -> Vec<(u64, u64, u64)> {
    let mut out = Vec::with_capacity(raw.len() / 24);
    for chunk in raw.as_chunks::<24>().0 {
        if let (Some(a), Some(b), Some(c)) = (
            read_u64_le(chunk, 0),
            read_u64_le(chunk, 8),
            read_u64_le(chunk, 16),
        ) {
            out.push((a, b, c));
        }
    }
    out
}

/// Decode HAL `TOP_ARRAY` 33-byte packed records into `(genome_idx, bottom_idx, paralogy_idx, parent_idx, reverse_flag)`.
#[must_use]
pub fn decode_hal_top_array(raw: &[u8]) -> Vec<(i64, i64, i64, i64, bool)> {
    let mut out = Vec::with_capacity(raw.len() / 33);
    for chunk in raw.as_chunks::<33>().0 {
        let g = i64::from_le_bytes(chunk[0..8].try_into().unwrap_or([0; 8]));
        let b = i64::from_le_bytes(chunk[8..16].try_into().unwrap_or([0; 8]));
        let p = i64::from_le_bytes(chunk[16..24].try_into().unwrap_or([0; 8]));
        let par = i64::from_le_bytes(chunk[24..32].try_into().unwrap_or([0; 8]));
        let rev = chunk[32] != 0;
        out.push((g, b, p, par, rev));
    }
    out
}

const fn decode_hal_nibble(nibble: u8) -> u8 {
    match nibble & 0x0f {
        0 => b'a',
        1 => b'c',
        2 => b'g',
        3 => b't',
        8 => b'A',
        9 => b'C',
        10 => b'G',
        11 => b'T',
        4 => b'n',
        _ => b'N',
    }
}

/// Decode HAL 4-bit nibble-packed `DNA_ARRAY` bytes (high nibble first) into ASCII DNA.
#[must_use]
pub fn decode_hal_dna(packed: &[u8], length: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(length.min(packed.len().saturating_mul(2)));
    for &b in packed {
        if out.len() >= length {
            break;
        }
        out.push(decode_hal_nibble(b >> 4));
        if out.len() >= length {
            break;
        }
        out.push(decode_hal_nibble(b));
    }
    out
}

#[cfg(test)]
mod tests {
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

        // Msg 3: DataLayout (0x0008) v3 class=2 ndim=1, size=12
        buf[p..p + 2].copy_from_slice(&0x0008_u16.to_le_bytes());
        buf[p + 2..p + 4].copy_from_slice(&12_u16.to_le_bytes());
        buf[p + 8] = 3; // layout_ver
        buf[p + 9] = 2; // chunked
        buf[p + 10] = 1; // ndim = 1
        buf[p + 11..p + 19].copy_from_slice(&0x100_u64.to_le_bytes());

        // Chunk B-tree v1 type-1 at 0x100 (ndim=1 => key_bytes=24, child_addr at +24)
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
            .find(|e| e.name != "Meta" && e.name != "Phylogeny")
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
        assert_eq!(root_names, vec!["hg38".to_owned(), "mm10".to_owned()]);

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
            vec![(0, 12, -1, vec![(0, false)])]
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
        assert_eq!(entries.len(), 301);
        assert!(entries.iter().any(|e| e.name == "ref0"));
        assert!(entries.iter().any(|e| e.name == "q_0000"));
        assert!(entries.iter().any(|e| e.name == "q_0299"));
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
}
