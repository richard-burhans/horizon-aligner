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
    bytes.starts_with(&HDF5_SIGNATURE)
}

fn concat_hal_contigs(
    species_name: &str,
    contigs: &[(String, Vec<u8>)],
) -> (Vec<u8>, Vec<HalSequence>) {
    let mut dna = Vec::new();
    let mut seqs = Vec::with_capacity(contigs.len().max(1));
    for (cname, cseq) in contigs {
        let offset = dna.len();
        dna.extend_from_slice(cseq);
        seqs.push(HalSequence {
            name: cname.clone(),
            length: cseq.len(),
            offset,
        });
    }
    if seqs.is_empty() {
        seqs.push(HalSequence {
            name: format!("{species_name}.chr1"),
            length: 0,
            offset: 0,
        });
    }
    (dna, seqs)
}

impl HalAlignment {
    /// Construct a star-tree `HalAlignment` from a reference genome, query `(name, dna)` pairs,
    /// and `(query_idx, ref_start, qry_start, len, rev)` blocks.
    #[must_use]
    pub fn from_star_blocks(
        ref_name: &str,
        ref_dna: &[u8],
        queries: &[(String, Vec<u8>)],
        blocks: &[(usize, usize, usize, usize, bool)], // (query_idx, ref_start, qry_start, len, rev)
    ) -> Self {
        type HalSpeciesContigs = (String, Vec<(String, Vec<u8>)>);
        let ref_contigs = [(format!("{ref_name}.chr1"), ref_dna.to_vec())];
        let qry_contigs: Vec<HalSpeciesContigs> = queries
            .iter()
            .map(|(q_name, q_dna)| {
                (
                    q_name.clone(),
                    vec![(format!("{q_name}.chr1"), q_dna.clone())],
                )
            })
            .collect();
        let mc_blocks: Vec<(usize, usize, usize, usize, usize, usize, bool)> = blocks
            .iter()
            .map(|&(q_idx, r_start, q_start, len, rev)| (q_idx, 0, 0, r_start, q_start, len, rev))
            .collect();
        Self::from_multi_contig_star_blocks(ref_name, &ref_contigs, &qry_contigs, &mc_blocks)
    }

    /// Construct a multi-contig star-tree `HalAlignment` from named reference contigs,
    /// query `(species_name, contigs)` pairs, and `(query_idx, ref_contig_idx, qry_contig_idx,
    /// ref_local_start, qry_local_start, len, rev)` blocks.
    #[must_use]
    #[allow(clippy::type_complexity)]
    pub fn from_multi_contig_star_blocks(
        ref_name: &str,
        ref_contigs: &[(String, Vec<u8>)],
        queries: &[(String, Vec<(String, Vec<u8>)>)],
        blocks: &[(usize, usize, usize, usize, usize, usize, bool)],
    ) -> Self {
        type HalIntervalHits = Vec<(usize, usize, bool)>;
        let (ref_dna, ref_seqs) = concat_hal_contigs(ref_name, ref_contigs);

        let mut qry_dnas: Vec<Vec<u8>> = Vec::with_capacity(queries.len());
        let mut qry_seqs: Vec<Vec<HalSequence>> = Vec::with_capacity(queries.len());
        for (q_name, q_contigs) in queries {
            let (dna, seqs) = concat_hal_contigs(q_name, q_contigs);
            qry_dnas.push(dna);
            qry_seqs.push(seqs);
        }

        let mut by_interval: BTreeMap<(usize, usize), HalIntervalHits> = BTreeMap::new();
        for &(q_idx, r_ci, q_ci, r_local, q_local, len, rev) in blocks {
            let Some(r_seq) = ref_seqs.get(r_ci) else {
                continue;
            };
            let Some(q_seq) = qry_seqs.get(q_idx).and_then(|s| s.get(q_ci)) else {
                continue;
            };
            let r_global = r_seq.offset + r_local;
            let q_global = q_seq.offset + q_local;
            by_interval
                .entry((r_global, len))
                .or_default()
                .push((q_idx, q_global, rev));
        }

        let mut bottom_segments = Vec::new();
        let mut child_tops_per_q: Vec<Vec<HalTopSegment>> = vec![Vec::new(); queries.len()];
        for ((r_start, len), hits) in by_interval {
            let b_idx = bottom_segments.len();
            let mut child_links = vec![None; queries.len()];
            let mut any_valid = false;
            for (q_idx, q_start, rev) in hits {
                if q_idx >= queries.len() {
                    continue;
                }
                any_valid = true;
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
            if any_valid {
                bottom_segments.push(HalBottomSegment {
                    start: r_start,
                    length: len,
                    child_tops: child_links,
                });
            }
        }

        let mut genomes = Vec::with_capacity(1 + queries.len());
        genomes.push(HalGenome {
            name: ref_name.to_owned(),
            parent: None,
            children: (1..=queries.len()).collect(),
            sequences: ref_seqs,
            dna: ref_dna,
            top_segments: Vec::new(),
            bottom_segments,
        });
        let mut leaf_names = Vec::with_capacity(queries.len());
        for (q_idx, (q_name, _)) in queries.iter().enumerate() {
            leaf_names.push(q_name.clone());
            genomes.push(HalGenome {
                name: q_name.clone(),
                parent: Some(0),
                children: Vec::new(),
                sequences: std::mem::take(&mut qry_seqs[q_idx]),
                dna: std::mem::take(&mut qry_dnas[q_idx]),
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
            let (p_seq, p_strand, p_start) =
                format_oriented_maf_slice(&parent.dna, bseg.start, p_end, tseg.parent_reversed);
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
                let Some((t_idx, _)) = *link else {
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
                    if (steps > 0 && idx == t_idx) || steps >= child.top_segments.len() {
                        break;
                    }
                    steps += 1;
                    let Some(tseg) = child.top_segments.get(idx) else {
                        break;
                    };
                    let rev = tseg.parent_reversed;
                    let r_len = r_end - bseg.start;
                    let Some(c_end) = tseg.start.checked_add(r_len) else {
                        cur_t = tseg.next_paralog;
                        continue;
                    };
                    if tseg.length >= r_len && c_end <= child.dna.len() {
                        let (q_seq, strand, q_start) =
                            format_oriented_maf_slice(&child.dna, tseg.start, c_end, rev);
                        let _ = writeln!(
                            out,
                            "s {} {} {} {} {} {}",
                            child.name,
                            q_start,
                            r_len,
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

fn format_oriented_maf_slice(
    dna: &[u8],
    start: usize,
    end: usize,
    reversed: bool,
) -> (std::borrow::Cow<'_, str>, char, usize) {
    let slice = &dna[start..end];
    if reversed {
        let rev = crate::revcomp(slice);
        (
            std::borrow::Cow::Owned(String::from_utf8_lossy(&rev).into_owned()),
            '-',
            dna.len().saturating_sub(end),
        )
    } else {
        (String::from_utf8_lossy(slice), '+', start)
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

/// Decode HAL `TOP_ARRAY` 33-byte packed records into `(start, bottom_parse_idx, paralogy_idx, parent_idx, reverse_flag)`.
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
mod tests;
