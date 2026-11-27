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

fn build_hal_sequences(
    species_name: &str,
    contigs: &[(String, Vec<u8>)],
) -> (usize, Vec<HalSequence>) {
    let mut offset = 0usize;
    let mut seqs = Vec::with_capacity(contigs.len().max(1));
    for (cname, cseq) in contigs {
        seqs.push(HalSequence {
            name: cname.clone(),
            length: cseq.len(),
            offset,
        });
        offset = offset.saturating_add(cseq.len());
    }
    if seqs.is_empty() {
        seqs.push(HalSequence {
            name: format!("{species_name}.chr1"),
            length: 0,
            offset: 0,
        });
    }
    (offset, seqs)
}

fn concat_hal_contigs(
    species_name: &str,
    contigs: &[(String, Vec<u8>)],
) -> (Vec<u8>, Vec<HalSequence>) {
    let (total_len, seqs) = build_hal_sequences(species_name, contigs);
    let mut dna = Vec::with_capacity(total_len);
    for (_, cseq) in contigs {
        dna.extend_from_slice(cseq);
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
            if len == 0 {
                continue;
            }
            let Some(r_seq) = ref_seqs.get(r_ci) else {
                continue;
            };
            let Some(q_seq) = qry_seqs.get(q_idx).and_then(|s| s.get(q_ci)) else {
                continue;
            };
            if r_local.checked_add(len).is_none_or(|e| e > r_seq.length)
                || q_local.checked_add(len).is_none_or(|e| e > q_seq.length)
            {
                continue;
            }
            let (Some(r_global), Some(q_global)) = (
                r_seq.offset.checked_add(r_local),
                q_seq.offset.checked_add(q_local),
            ) else {
                continue;
            };
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
            let mut child_tails: Vec<Option<usize>> = vec![None; queries.len()];
            for (q_idx, q_start, rev) in hits {
                let t_idx = child_tops_per_q[q_idx].len();
                child_tops_per_q[q_idx].push(HalTopSegment {
                    start: q_start,
                    length: len,
                    parent_bottom_idx: Some(b_idx),
                    parent_reversed: rev,
                    next_paralog: None,
                });
                if let Some(prev_tail) = child_tails[q_idx] {
                    child_tops_per_q[q_idx][prev_tail].next_paralog = Some(t_idx);
                } else {
                    child_links[q_idx] = Some((t_idx, rev));
                }
                child_tails[q_idx] = Some(t_idx);
            }
            bottom_segments.push(HalBottomSegment {
                start: r_start,
                length: len,
                child_tops: child_links,
            });
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

    /// Construct a progressive (non-star) multi-contig `HalAlignment` by threading
    /// `queries` along a pectinate guide tree rooted at `Anc_{k-1}`:
    /// `((...((ref_name,q0)Anc0,q1)Anc1...),q_{k-1})Anc_{k-1};`.
    /// Falls back to [`Self::from_multi_contig_star_blocks`] when `queries.len() < 2`.
    #[must_use]
    #[allow(clippy::type_complexity, clippy::too_many_lines)]
    pub fn from_multi_contig_progressive_blocks(
        ref_name: &str,
        ref_contigs: &[(String, Vec<u8>)],
        queries: &[(String, Vec<(String, Vec<u8>)>)],
        blocks: &[(usize, usize, usize, usize, usize, usize, bool)],
    ) -> Self {
        if queries.len() < 2 {
            return Self::from_multi_contig_star_blocks(ref_name, ref_contigs, queries, blocks);
        }
        let mut star = Self::from_multi_contig_star_blocks(ref_name, ref_contigs, queries, blocks);
        normalize::normalize_star_hal_in_place(&mut star);
        let n_q = queries.len();
        let n_b = star.genomes[0].bottom_segments.len();

        // Build internal ancestral nodes Anc0..Anc_{n_q-1} sharing ref_genome's coordinate frame,
        // inferring ancestral bases via exact two-pass Fitch parsimony (Fitch 1971) with
        // leaf majority-vote root tie-breaking across aligned leaves.
        // Genome layout:
        //   0..n_q         : Anc_{n_q-1} (root, idx 0) .. Anc0 (idx n_q - 1)
        //   n_q            : ref_name (leaf, child of Anc0)
        //   n_q + 1 ..     : queries[0..n_q] (query j is child of Anc_j)
        let ref_leaf_idx = n_q;
        let q_leaf_base = n_q + 1;
        let mut anc_dnas = infer_ancestral_dnas(&star, n_q);

        let mut star_iter = star.genomes.into_iter();
        let Some(ref_genome) = star_iter.next() else {
            return Self {
                tree_newick: star.tree_newick,
                genomes: Vec::new(),
            };
        };
        let query_genomes: Vec<HalGenome> = star_iter.collect();

        let mut genomes = Vec::with_capacity(1 + 2 * n_q);

        for level_rev in 0..n_q {
            let level = n_q - 1 - level_rev; // n_q-1 down to 0
            let anc_name = format!("Anc{level}");
            let parent = if level_rev == 0 {
                None
            } else {
                Some(level_rev - 1)
            };
            let left_child = if level == 0 {
                ref_leaf_idx
            } else {
                level_rev + 1
            };
            let right_child = q_leaf_base + level;

            let top_segments = if parent.is_some() {
                (0..n_b)
                    .map(|bi| HalTopSegment {
                        start: ref_genome.bottom_segments[bi].start,
                        length: ref_genome.bottom_segments[bi].length,
                        parent_bottom_idx: Some(bi),
                        parent_reversed: false,
                        next_paralog: None,
                    })
                    .collect()
            } else {
                Vec::new()
            };

            let bottom_segments = (0..n_b)
                .map(|bi| {
                    let bseg = &ref_genome.bottom_segments[bi];
                    let left_link = Some((bi, false));
                    let right_link = bseg.child_tops[level];
                    HalBottomSegment {
                        start: bseg.start,
                        length: bseg.length,
                        child_tops: vec![left_link, right_link],
                    }
                })
                .collect();

            genomes.push(HalGenome {
                name: anc_name,
                parent,
                children: vec![left_child, right_child],
                sequences: ref_genome.sequences.clone(),
                dna: std::mem::take(&mut anc_dnas[level]),
                top_segments,
                bottom_segments,
            });
        }

        // ref_name leaf (child of Anc0 at index n_q - 1)
        let anc0_idx = n_q - 1;
        let ref_tops = (0..n_b)
            .map(|bi| HalTopSegment {
                start: ref_genome.bottom_segments[bi].start,
                length: ref_genome.bottom_segments[bi].length,
                parent_bottom_idx: Some(bi),
                parent_reversed: false,
                next_paralog: None,
            })
            .collect();
        genomes.push(HalGenome {
            name: ref_genome.name,
            parent: Some(anc0_idx),
            children: Vec::new(),
            sequences: ref_genome.sequences,
            dna: ref_genome.dna,
            top_segments: ref_tops,
            bottom_segments: Vec::new(),
        });

        // Query leaves (query j is child of Anc_j at index n_q - 1 - j)
        for (j, mut g) in query_genomes.into_iter().enumerate() {
            let parent_idx = n_q - 1 - j;
            g.parent = Some(parent_idx);
            genomes.push(g);
        }

        let mut tree_newick = ref_name.to_owned();
        for (j, (q_name, _)) in queries.iter().enumerate() {
            tree_newick = format!("({tree_newick},{q_name})Anc{j}");
        }
        tree_newick.push(';');

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

use fitch::{
    FITCH_ALL, add_consensus_base, fitch_base_mask, fitch_resolve_child, fitch_resolve_root,
    fitch_step, pack_fitch_masks,
};

fn infer_ancestral_dnas(star: &HalAlignment, n_q: usize) -> Vec<Vec<u8>> {
    let ref_genome = &star.genomes[0];
    let dna_len = ref_genome.dna.len();
    let mut leaf_counts: Vec<[u16; 4]> = ref_genome
        .dna
        .iter()
        .map(|&rb| {
            let mut c = [0u16; 4];
            add_consensus_base(&mut c, rb);
            c
        })
        .collect();

    // Pass 1 (bottom-up): compute Fitch state sets for Anc0..Anc_{n_q-1}.
    let mut prev_sets: Vec<u8> = ref_genome
        .dna
        .iter()
        .map(|&rb| fitch_base_mask(rb))
        .collect();
    let mut level_sets: Vec<Vec<u8>> = Vec::with_capacity(n_q);
    let mut q_masks = vec![FITCH_ALL; dna_len];

    for level in 0..n_q {
        q_masks.fill(FITCH_ALL);
        if let Some(q_gen) = star.genomes.get(1 + level) {
            for bseg in &ref_genome.bottom_segments {
                let Some((t_idx, _)) = bseg.child_tops[level] else {
                    continue;
                };
                let Some(tseg) = q_gen.top_segments.get(t_idx) else {
                    continue;
                };
                let r_end = bseg.start.saturating_add(bseg.length).min(dna_len);
                for (off, (qm, counts)) in q_masks[bseg.start..r_end]
                    .iter_mut()
                    .zip(leaf_counts[bseg.start..r_end].iter_mut())
                    .enumerate()
                {
                    let q_pos = if tseg.parent_reversed {
                        tseg.start + bseg.length - 1 - off
                    } else {
                        tseg.start + off
                    };
                    if let Some(&qb) = q_gen.dna.get(q_pos) {
                        let oriented = if tseg.parent_reversed {
                            crate::seeding::complement_base(qb)
                        } else {
                            qb
                        };
                        *qm = fitch_base_mask(oriented);
                        add_consensus_base(counts, oriented);
                    }
                }
            }
        }
        for (prev, &r) in prev_sets.iter_mut().zip(q_masks.iter()) {
            *prev = fitch_step(*prev, r);
        }
        level_sets.push(pack_fitch_masks(&prev_sets));
    }

    // Pass 2 (top-down): resolve root Anc_{n_q-1} down to Anc0.
    let mut anc_dnas: Vec<Vec<u8>> = vec![Vec::new(); n_q];
    if n_q == 0 {
        return anc_dnas;
    }
    let mut root_dna = vec![0u8; dna_len];
    fitch_resolve_root(
        &ref_genome.dna,
        &level_sets[n_q - 1],
        &leaf_counts,
        &mut root_dna,
    );
    anc_dnas[n_q - 1] = root_dna;
    for level in (0..n_q - 1).rev() {
        let mut child_dna = vec![0u8; dna_len];
        fitch_resolve_child(
            &ref_genome.dna,
            &anc_dnas[level + 1],
            &level_sets[level],
            &mut child_dna,
        );
        anc_dnas[level] = child_dna;
    }
    anc_dnas
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
mod fitch;
mod hdf5_dataset;
mod hdf5_progressive;
mod hdf5_read;
mod hdf5_stream;
mod hdf5_write;
mod hdf5_writer;
mod normalize;
mod stream_hal;
mod stream_hal_genome;
mod stream_hal_maf;
pub use hdf5_dataset::*;
pub use hdf5_read::*;
pub use hdf5_stream::*;
pub use stream_hal::*;
pub use stream_hal_maf::{HalProjectedBlock, HalProjectedRow};

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
mod tests_hdf5;
#[cfg(test)]
mod tests_model;
#[cfg(test)]
mod tests_progressive;
#[cfg(test)]
mod tests_progressive_hdf5;
#[cfg(test)]
mod tests_stream_hal;
