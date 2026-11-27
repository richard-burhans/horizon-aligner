//! Pure-Rust streaming multi-level HAL reader (`stats`, `liftover`, `to-maf`).
//!
//! Uses [`Hdf5DatasetReader`] chunk-slice reads so multi-gigabyte HAL alignments (e.g.
//! `5.6 GiB` orangutan and `57 GiB` 363-way bird HALs) can be inspected, lifted over across
//! arbitrary multi-level phylogenetic tree paths, and projected to MAF/GFA in bounded RAM
//! with zero `libhdf5` or C++ `halTools` dependencies.

use super::decode_hal_dna;
use super::hdf5_read::{read_at, read_hdf5_superblock, read_hdf5_v0_group_entries_seek};
use super::hdf5_stream::{DEFAULT_CHUNK_CACHE_CAP, read_hdf5_v0_attr_string_seek};
use super::stream_hal_genome::StreamGenomeHandle;
pub use super::stream_hal_genome::{HalGenomeStats, HalLiftoverHit, HalMafStats, HalStreamSeq};
use super::stream_hal_maf::parse_hal_newick;
use std::collections::HashMap;
use std::io::{BufReader, Read, Seek};
use std::path::Path;

/// Directed tree hop from `from` to `to` (either up to parent or down into `child_slot`).
#[derive(Clone, Copy)]
pub(super) struct TreeStep {
    pub(super) from: usize,
    pub(super) to: usize,
    pub(super) going_up: bool,
    pub(super) child_slot: usize,
}

/// Source-tracked lifted alignment fragment `[src_s, src_e) -> [dst_s, dst_e)` with cumulative
/// relative orientation `rev`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct LiftedFrag {
    pub(super) src_s: u64,
    pub(super) src_e: u64,
    pub(super) dst_s: u64,
    pub(super) dst_e: u64,
    pub(super) rev: bool,
}

/// Decoded single top segment `(start, length, paralogy_idx, parent_bot_idx, rev)`.
#[derive(Clone, Copy)]
pub(super) struct TopSegRec {
    pub(super) start: u64,
    pub(super) length: u64,
    paralogy_idx: i64,
    parent_bot_idx: i64,
    rev: bool,
}

/// Decoded single bottom segment `(start, length, child_top_idx)` for one child slot.
#[derive(Clone, Copy)]
pub(super) struct BotSegChildRec {
    pub(super) start: u64,
    pub(super) length: u64,
    child_top_idx: i64,
}

/// Streaming pure-Rust reader over any HDF5 v0/v1 `.hal` file.
pub struct StreamingHalFile<R: Read + Seek> {
    pub(super) reader: R,
    tree_newick: String,
    pub(super) genomes: Vec<StreamGenomeHandle>,
    pub(super) name_to_idx: HashMap<String, usize>,
    chunk_cache_cap: usize,
    seg_scratch: Vec<u8>,
}

impl StreamingHalFile<BufReader<std::fs::File>> {
    /// Open an HDF5 `.hal` file on disk for streaming queries.
    ///
    /// # Errors
    /// Returns an `std::io::Error` if the file cannot be opened or is not a valid HDF5 v0/v1 HAL.
    pub fn open(path: &Path) -> std::io::Result<Self> {
        let f = std::fs::File::open(path)?;
        Self::open_reader(BufReader::with_capacity(64 * 1024, f))
    }
}

impl<R: Read + Seek> StreamingHalFile<R> {
    /// Open any `Read + Seek` HDF5 `.hal` stream lazily (parsing only the root group and
    /// `Phylogeny` tree; per-genome metadata and `SEQIDX`/`SEQNAME` tables load on demand).
    ///
    /// # Errors
    /// Returns an `std::io::Error` if the stream is not a valid HDF5 v0/v1 HAL container.
    pub fn open_reader(mut reader: R) -> std::io::Result<Self> {
        let lead = read_at(&mut reader, 0, 96)?;
        let sb = read_hdf5_superblock(&lead).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "not a valid HDF5 .hal file (invalid superblock)",
            )
        })?;
        let root_addr = sb.root_group_addr.ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "unsupported HDF5 superblock (missing 64-bit root group address)",
            )
        })?;
        let root_entries =
            read_hdf5_v0_group_entries_seek(&mut reader, root_addr)?.ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "failed to read HDF5 v0 root group",
                )
            })?;
        let root_map: HashMap<String, u64> = root_entries
            .into_iter()
            .map(|e| (e.name, e.object_header_addr))
            .collect();
        let phylo_addr = *root_map.get("Phylogeny").ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "HAL file missing root `Phylogeny` group",
            )
        })?;
        let tree_newick = read_hdf5_v0_attr_string_seek(&mut reader, phylo_addr, "Phylogeny")?
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "HAL `Phylogeny` group missing readable Newick string attribute",
                )
            })?;
        let tree_nodes = parse_hal_newick(&tree_newick)?;
        let mut genomes = Vec::with_capacity(tree_nodes.len());
        let mut name_to_idx = HashMap::with_capacity(tree_nodes.len());
        for (idx, tn) in tree_nodes.into_iter().enumerate() {
            if name_to_idx.contains_key(&tn.name) {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("duplicate genome name `{}` in HAL Phylogeny tree", tn.name),
                ));
            }
            let group_addr = root_map.get(&tn.name).copied().unwrap_or(u64::MAX);
            name_to_idx.insert(tn.name.clone(), idx);
            genomes.push(StreamGenomeHandle::new_lazy(
                tn.name,
                tn.parent,
                tn.children,
                group_addr,
            ));
        }
        Ok(Self {
            reader,
            tree_newick,
            genomes,
            name_to_idx,
            chunk_cache_cap: DEFAULT_CHUNK_CACHE_CAP,
            seg_scratch: Vec::with_capacity(256),
        })
    }

    /// Tune the per-dataset LRU chunk cache capacity (`>= 1` chunks) for all open and future
    /// dataset readers in this HAL file.
    pub fn set_chunk_cache_cap(&mut self, cap: usize) {
        let cap = cap.max(1);
        self.chunk_cache_cap = cap;
        for g in &mut self.genomes {
            g.set_cache_cap(cap);
        }
    }

    /// Return the configured per-dataset LRU chunk cache capacity.
    #[must_use]
    pub const fn chunk_cache_cap(&self) -> usize {
        self.chunk_cache_cap
    }

    /// Return the parsed Newick guide tree string.
    #[must_use]
    pub fn tree_newick(&self) -> &str {
        &self.tree_newick
    }

    /// Return genome-level summary statistics matching `halStats` without loading full
    /// `SEQNAME_ARRAY` tables.
    ///
    /// # Errors
    /// Returns an `std::io::Error` if reading any genome's HDF5 metadata fails.
    pub fn stats(&mut self) -> std::io::Result<Vec<HalGenomeStats>> {
        let n = self.genomes.len();
        let mut out = Vec::with_capacity(n);
        for i in 0..n {
            let total_bp = self.genomes[i].ensure_total_bp(&mut self.reader)?;
            let g = &self.genomes[i];
            let parent = g.parent.map(|p| self.genomes[p].name.clone());
            out.push(HalGenomeStats {
                name: g.name.clone(),
                parent,
                num_children: g.children.len(),
                num_sequences: g.num_sequences,
                total_bp,
                num_top_segments: g.num_top_segments,
                num_bottom_segments: g.num_bottom_segments,
            });
        }
        Ok(out)
    }

    /// Return sequences for `genome_name`, lazily loading its `SEQIDX`/`SEQNAME` tables on first
    /// access.
    ///
    /// # Errors
    /// Returns an `std::io::Error` if reading HDF5 sequence metadata fails.
    pub fn genome_sequences(
        &mut self,
        genome_name: &str,
    ) -> std::io::Result<Option<&[HalStreamSeq]>> {
        let Some(&idx) = self.name_to_idx.get(genome_name) else {
            return Ok(None);
        };
        self.genomes[idx].ensure_sequences(&mut self.reader)?;
        Ok(Some(&self.genomes[idx].sequences))
    }

    pub(super) fn ensure_genome_sequences(&mut self, genome_idx: usize) -> std::io::Result<()> {
        self.genomes[genome_idx].ensure_sequences(&mut self.reader)
    }

    pub(super) fn ensure_genome_header(&mut self, genome_idx: usize) -> std::io::Result<()> {
        self.genomes[genome_idx].ensure_header(&mut self.reader)
    }

    #[cfg(test)]
    pub(super) fn is_genome_sequences_loaded(&self, genome_name: &str) -> bool {
        self.name_to_idx
            .get(genome_name)
            .is_some_and(|&i| self.genomes[i].sequences_loaded)
    }

    #[cfg(test)]
    pub(super) fn is_genome_offsets_loaded(&self, genome_name: &str) -> bool {
        self.name_to_idx
            .get(genome_name)
            .is_some_and(|&i| self.genomes[i].seq_offsets_loaded)
    }

    /// Compute the directed edge path from `src` to `dst` through their LCA:
    /// each step is `(from_idx, to_idx, going_up, child_slot)`.
    pub(super) fn find_tree_path(&self, src: usize, dst: usize) -> Vec<TreeStep> {
        if src == dst {
            return Vec::new();
        }
        let mut src_anc = Vec::new();
        let mut cur = Some(src);
        while let Some(u) = cur {
            src_anc.push(u);
            cur = self.genomes[u].parent;
        }
        let src_pos: HashMap<usize, usize> =
            src_anc.iter().enumerate().map(|(i, &u)| (u, i)).collect();

        let mut dst_anc = Vec::new();
        let mut cur = Some(dst);
        let mut lca_src_pos = None;
        while let Some(u) = cur {
            if let Some(&sp) = src_pos.get(&u) {
                lca_src_pos = Some(sp);
                break;
            }
            dst_anc.push(u);
            cur = self.genomes[u].parent;
        }
        let Some(lca_pos) = lca_src_pos else {
            return Vec::new();
        };
        let mut steps = Vec::new();
        for i in 0..lca_pos {
            steps.push(TreeStep {
                from: src_anc[i],
                to: src_anc[i + 1],
                going_up: true,
                child_slot: 0,
            });
        }
        let mut parent = src_anc[lca_pos];
        for &child in dst_anc.iter().rev() {
            let slot = self.genomes[parent]
                .children
                .iter()
                .position(|&c| c == child)
                .unwrap_or(0);
            steps.push(TreeStep {
                from: parent,
                to: child,
                going_up: false,
                child_slot: slot,
            });
            parent = child;
        }
        steps
    }

    pub(super) fn read_top_segment(
        &mut self,
        genome_idx: usize,
        seg_idx: u64,
    ) -> std::io::Result<Option<TopSegRec>> {
        self.genomes[genome_idx].ensure_header(&mut self.reader)?;
        let n_tops = self.genomes[genome_idx].num_top_segments;
        if seg_idx >= n_tops {
            return Ok(None);
        }
        let cache_cap = self.chunk_cache_cap;
        {
            let g = &mut self.genomes[genome_idx];
            let rdr = g.ensure_top_reader(&mut self.reader, cache_cap)?;
            rdr.read_elem_range_into(&mut self.reader, seg_idx, 2, &mut self.seg_scratch)?;
        }
        let raw = &self.seg_scratch;
        if raw.len() < 66 {
            return Ok(None);
        }
        let start = i64::from_le_bytes(raw[0..8].try_into().unwrap_or([0; 8]));
        let paralogy_idx = i64::from_le_bytes(raw[16..24].try_into().unwrap_or([0xFF; 8]));
        let parent_bot_idx = i64::from_le_bytes(raw[24..32].try_into().unwrap_or([0xFF; 8]));
        let rev = raw[32] != 0;
        let next_start = i64::from_le_bytes(raw[33..41].try_into().unwrap_or([0; 8]));
        if start < 0 || next_start <= start {
            return Ok(None);
        }
        Ok(Some(TopSegRec {
            start: start.cast_unsigned(),
            length: (next_start - start).cast_unsigned(),
            paralogy_idx,
            parent_bot_idx,
            rev,
        }))
    }

    pub(super) fn read_bot_segment_child(
        &mut self,
        genome_idx: usize,
        seg_idx: u64,
        child_slot: Option<usize>,
    ) -> std::io::Result<Option<BotSegChildRec>> {
        self.genomes[genome_idx].ensure_header(&mut self.reader)?;
        let n_bots = self.genomes[genome_idx].num_bottom_segments;
        if seg_idx >= n_bots {
            return Ok(None);
        }
        let cache_cap = self.chunk_cache_cap;
        let es = {
            let g = &mut self.genomes[genome_idx];
            let rdr = g.ensure_bot_reader(&mut self.reader, cache_cap)?;
            let es = rdr.elem_size();
            rdr.read_elem_range_into(&mut self.reader, seg_idx, 2, &mut self.seg_scratch)?;
            es
        };
        let raw = &self.seg_scratch;
        if es < 24 || raw.len() < es * 2 {
            return Ok(None);
        }
        let start = i64::from_le_bytes(raw[0..8].try_into().unwrap_or([0; 8]));
        let next_start = i64::from_le_bytes(raw[es..es + 8].try_into().unwrap_or([0; 8]));
        if start < 0 || next_start <= start {
            return Ok(None);
        }
        let child_top_idx = if let Some(slot) = child_slot {
            let off = 24 + slot * 9;
            if off + 9 <= es {
                i64::from_le_bytes(raw[off..off + 8].try_into().unwrap_or([0xFF; 8]))
            } else {
                -1
            }
        } else {
            -1
        };
        Ok(Some(BotSegChildRec {
            start: start.cast_unsigned(),
            length: (next_start - start).cast_unsigned(),
            child_top_idx,
        }))
    }

    /// Bounded binary search over top-segment indices `[lo_init, hi_init)` returning the first
    /// segment whose end `start + length > pos`.
    fn binary_search_top_range(
        &mut self,
        genome_idx: usize,
        pos: u64,
        mut lo: u64,
        mut hi: u64,
    ) -> std::io::Result<u64> {
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let Some(seg) = self.read_top_segment(genome_idx, mid)? else {
                break;
            };
            if seg.start.saturating_add(seg.length) <= pos {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        Ok(lo)
    }

    /// Bounded binary search over bottom-segment indices `[lo_init, hi_init)` returning the first
    /// segment whose end `start + length > pos`.
    fn binary_search_bot_range(
        &mut self,
        genome_idx: usize,
        pos: u64,
        mut lo: u64,
        mut hi: u64,
    ) -> std::io::Result<u64> {
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let Some(seg) = self.read_bot_segment_child(genome_idx, mid, None)? else {
                break;
            };
            if seg.start.saturating_add(seg.length) <= pos {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        Ok(lo)
    }

    pub(super) fn find_first_top_covering(
        &mut self,
        genome_idx: usize,
        pos: u64,
    ) -> std::io::Result<u64> {
        self.genomes[genome_idx].ensure_seq_offsets(&mut self.reader)?;
        let n = self.genomes[genome_idx].num_top_segments;
        if n == 0 {
            return Ok(0);
        }
        let (seq_lo, seq_hi) = self.genomes[genome_idx]
            .seq_seg_bounds(pos, true, n)
            .unwrap_or((0, n));

        // Fast path: check sequential / local galloping hint inside `[seq_lo, seq_hi)`.
        if let Some(hint) = self.genomes[genome_idx].last_top_hint
            && hint >= seq_lo
            && hint < seq_hi
            && let Some(hseg) = self.read_top_segment(genome_idx, hint)?
        {
            let h_end = hseg.start.saturating_add(hseg.length);
            if pos >= hseg.start && pos < h_end {
                return Ok(hint);
            }
            if pos >= h_end {
                // Gallop forward (1, 2, 4, 8, 16, ...) to bracket local fragments in the same chunk.
                let mut prev = hint + 1;
                let mut step = 1u64;
                let mut found_bound = None;
                while prev < seq_hi && step <= 64 {
                    let probe = (hint + step).min(seq_hi - 1);
                    let Some(pseg) = self.read_top_segment(genome_idx, probe)? else {
                        break;
                    };
                    if pseg.start.saturating_add(pseg.length) > pos {
                        found_bound = Some((prev, probe + 1));
                        break;
                    }
                    prev = probe + 1;
                    if probe + 1 >= seq_hi {
                        break;
                    }
                    step <<= 1;
                }
                if let Some((g_lo, g_hi)) = found_bound {
                    let idx = self.binary_search_top_range(genome_idx, pos, g_lo, g_hi)?;
                    self.genomes[genome_idx].last_top_hint = Some(idx);
                    return Ok(idx);
                }
            }
        }

        // Per-sequence offset-bounded binary search over `[seq_lo, seq_hi)`.
        let mut idx = self.binary_search_top_range(genome_idx, pos, seq_lo, seq_hi)?;
        // Safety fallback for malformed/non-monotonic SEQIDX tables: verify bracketed result.
        if (seq_lo > 0 || seq_hi < n)
            && (idx == seq_hi
                || (idx == seq_lo
                    && self
                        .read_top_segment(genome_idx, idx)?
                        .is_some_and(|s| s.start > pos)))
        {
            idx = self.binary_search_top_range(genome_idx, pos, 0, n)?;
        }
        self.genomes[genome_idx].last_top_hint = Some(idx);
        Ok(idx)
    }

    pub(super) fn find_first_bot_covering(
        &mut self,
        genome_idx: usize,
        pos: u64,
    ) -> std::io::Result<u64> {
        self.genomes[genome_idx].ensure_seq_offsets(&mut self.reader)?;
        let n = self.genomes[genome_idx].num_bottom_segments;
        if n == 0 {
            return Ok(0);
        }
        let (seq_lo, seq_hi) = self.genomes[genome_idx]
            .seq_seg_bounds(pos, false, n)
            .unwrap_or((0, n));

        // Fast path: check sequential / local galloping hint inside `[seq_lo, seq_hi)`.
        if let Some(hint) = self.genomes[genome_idx].last_bot_hint
            && hint >= seq_lo
            && hint < seq_hi
            && let Some(hseg) = self.read_bot_segment_child(genome_idx, hint, None)?
        {
            let h_end = hseg.start.saturating_add(hseg.length);
            if pos >= hseg.start && pos < h_end {
                return Ok(hint);
            }
            if pos >= h_end {
                let mut prev = hint + 1;
                let mut step = 1u64;
                let mut found_bound = None;
                while prev < seq_hi && step <= 64 {
                    let probe = (hint + step).min(seq_hi - 1);
                    let Some(pseg) = self.read_bot_segment_child(genome_idx, probe, None)? else {
                        break;
                    };
                    if pseg.start.saturating_add(pseg.length) > pos {
                        found_bound = Some((prev, probe + 1));
                        break;
                    }
                    prev = probe + 1;
                    if probe + 1 >= seq_hi {
                        break;
                    }
                    step <<= 1;
                }
                if let Some((g_lo, g_hi)) = found_bound {
                    let idx = self.binary_search_bot_range(genome_idx, pos, g_lo, g_hi)?;
                    self.genomes[genome_idx].last_bot_hint = Some(idx);
                    return Ok(idx);
                }
            }
        }

        // Per-sequence offset-bounded binary search over `[seq_lo, seq_hi)`.
        let mut idx = self.binary_search_bot_range(genome_idx, pos, seq_lo, seq_hi)?;
        if (seq_lo > 0 || seq_hi < n)
            && (idx == seq_hi
                || (idx == seq_lo
                    && self
                        .read_bot_segment_child(genome_idx, idx, None)?
                        .is_some_and(|s| s.start > pos)))
        {
            idx = self.binary_search_bot_range(genome_idx, pos, 0, n)?;
        }
        self.genomes[genome_idx].last_bot_hint = Some(idx);
        Ok(idx)
    }

    fn step_up(
        &mut self,
        child_idx: usize,
        parent_idx: usize,
        frag: LiftedFrag,
    ) -> std::io::Result<Vec<LiftedFrag>> {
        let mut out = Vec::new();
        let mut idx = self.find_first_top_covering(child_idx, frag.dst_s)?;
        let n_tops = self.genomes[child_idx].num_top_segments;
        while idx < n_tops {
            let Some(tseg) = self.read_top_segment(child_idx, idx)? else {
                break;
            };
            if tseg.start >= frag.dst_e {
                break;
            }
            if tseg.parent_bot_idx >= 0 {
                let ov_s = frag.dst_s.max(tseg.start);
                let ov_e = frag.dst_e.min(tseg.start.saturating_add(tseg.length));
                if ov_e > ov_s
                    && let Some(bseg) = self.read_bot_segment_child(
                        parent_idx,
                        tseg.parent_bot_idx.cast_unsigned(),
                        None,
                    )?
                {
                    if bseg.length != tseg.length {
                        return Err(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            "HAL top/bottom segment length mismatch on step_up",
                        ));
                    }
                    let off_s = ov_s - tseg.start;
                    let off_e = ov_e - tseg.start;
                    let d_lo = ov_s - frag.dst_s;
                    let d_hi = ov_e - frag.dst_s;
                    let (sub_src_s, sub_src_e) = if frag.rev {
                        (frag.src_e - d_hi, frag.src_e - d_lo)
                    } else {
                        (frag.src_s + d_lo, frag.src_s + d_hi)
                    };
                    let (p_s, p_e) = if tseg.rev {
                        (
                            bseg.start.saturating_add(bseg.length - off_e),
                            bseg.start.saturating_add(bseg.length - off_s),
                        )
                    } else {
                        (
                            bseg.start.saturating_add(off_s),
                            bseg.start.saturating_add(off_e),
                        )
                    };
                    if p_e > p_s {
                        out.push(LiftedFrag {
                            src_s: sub_src_s,
                            src_e: sub_src_e,
                            dst_s: p_s,
                            dst_e: p_e,
                            rev: frag.rev ^ tseg.rev,
                        });
                    }
                }
            }
            idx += 1;
        }
        self.genomes[child_idx].last_top_hint = Some(idx.saturating_sub(1));
        Ok(out)
    }

    fn step_down(
        &mut self,
        parent_idx: usize,
        child_slot: usize,
        child_idx: usize,
        frag: LiftedFrag,
    ) -> std::io::Result<Vec<LiftedFrag>> {
        let mut out = Vec::new();
        let mut idx = self.find_first_bot_covering(parent_idx, frag.dst_s)?;
        let n_bots = self.genomes[parent_idx].num_bottom_segments;
        while idx < n_bots {
            let Some(bseg) = self.read_bot_segment_child(parent_idx, idx, Some(child_slot))? else {
                break;
            };
            if bseg.start >= frag.dst_e {
                break;
            }
            if bseg.child_top_idx >= 0 {
                let ov_s = frag.dst_s.max(bseg.start);
                let ov_e = frag.dst_e.min(bseg.start.saturating_add(bseg.length));
                if ov_e > ov_s {
                    let off_s = ov_s - bseg.start;
                    let off_e = ov_e - bseg.start;
                    let d_lo = ov_s - frag.dst_s;
                    let d_hi = ov_e - frag.dst_s;
                    let (sub_src_s, sub_src_e) = if frag.rev {
                        (frag.src_e - d_hi, frag.src_e - d_lo)
                    } else {
                        (frag.src_s + d_lo, frag.src_s + d_hi)
                    };
                    let mut cur_t = bseg.child_top_idx;
                    let mut seen_t: [i64; 64] = [-1; 64];
                    let mut seen_len = 0usize;
                    while cur_t >= 0 && seen_len < 64 && !seen_t[..seen_len].contains(&cur_t) {
                        seen_t[seen_len] = cur_t;
                        seen_len += 1;
                        let Some(tseg) = self.read_top_segment(child_idx, cur_t.cast_unsigned())?
                        else {
                            break;
                        };
                        if tseg.length != bseg.length {
                            return Err(std::io::Error::new(
                                std::io::ErrorKind::InvalidData,
                                "HAL bottom/top segment length mismatch on step_down",
                            ));
                        }
                        let (c_s, c_e) = if tseg.rev {
                            (
                                tseg.start.saturating_add(tseg.length - off_e),
                                tseg.start.saturating_add(tseg.length - off_s),
                            )
                        } else {
                            (
                                tseg.start.saturating_add(off_s),
                                tseg.start.saturating_add(off_e),
                            )
                        };
                        if c_e > c_s {
                            out.push(LiftedFrag {
                                src_s: sub_src_s,
                                src_e: sub_src_e,
                                dst_s: c_s,
                                dst_e: c_e,
                                rev: frag.rev ^ tseg.rev,
                            });
                        }
                        cur_t = tseg.paralogy_idx;
                    }
                }
            }
            idx += 1;
        }
        self.genomes[parent_idx].last_bot_hint = Some(idx.saturating_sub(1));
        Ok(out)
    }

    pub(super) fn liftover_global(
        &mut self,
        path: &[TreeStep],
        start: u64,
        end: u64,
    ) -> std::io::Result<Vec<LiftedFrag>> {
        const MAX_ACTIVE_FRAGMENTS: usize = 65_536;
        let mut cur = vec![LiftedFrag {
            src_s: start,
            src_e: end,
            dst_s: start,
            dst_e: end,
            rev: false,
        }];
        for &step in path {
            let mut next = Vec::new();
            for frag in cur {
                let projected = if step.going_up {
                    self.step_up(step.from, step.to, frag)?
                } else {
                    self.step_down(step.from, step.child_slot, step.to, frag)?
                };
                next.extend(projected);
                if next.len() > MAX_ACTIVE_FRAGMENTS {
                    break;
                }
            }
            if next.is_empty() {
                return Ok(Vec::new());
            }
            // Coalesce ONLY strictly adjacent source-order, destination-collinear fragments on the
            // same destination contig — preventing cross-contig merging, inversion/transposition
            // scrambling, and paralog interval inflation. Uses compact `seq_offsets` without
            // loading `SEQNAME_ARRAY` on intermediate ancestors.
            self.genomes[step.to].ensure_seq_offsets(&mut self.reader)?;
            next.sort_unstable_by_key(|f| (f.src_s, f.dst_s, f.rev));
            next.dedup();
            let dst_g = &self.genomes[step.to];
            let mut merged: Vec<LiftedFrag> = Vec::with_capacity(next.len());
            for f in next {
                if let Some(last) = merged.last_mut()
                    && last.rev == f.rev
                    && last.src_e == f.src_s
                    && ((!f.rev && last.dst_e == f.dst_s) || (f.rev && f.dst_e == last.dst_s))
                    && dst_g.seq_idx_for_pos_offsets(last.dst_s)
                        == dst_g.seq_idx_for_pos_offsets(f.dst_s)
                {
                    last.src_e = f.src_e;
                    if f.rev {
                        last.dst_s = f.dst_s;
                    } else {
                        last.dst_e = f.dst_e;
                    }
                } else {
                    merged.push(f);
                }
            }
            if merged.len() > MAX_ACTIVE_FRAGMENTS {
                merged.truncate(MAX_ACTIVE_FRAGMENTS);
            }
            cur = merged;
        }
        Ok(cur)
    }

    /// Read and unpack an exact ASCII DNA slice `[global_start .. global_start + len)` from
    /// `genome_idx`'s `DNA_ARRAY` via chunk-slice I/O (clamped to `genome.total_bp`).
    ///
    /// # Errors
    /// Returns an `std::io::Error` if `genome_idx` is invalid or HDF5 reading fails.
    pub fn read_dna_slice(
        &mut self,
        genome_idx: usize,
        global_start: u64,
        len: usize,
    ) -> std::io::Result<Vec<u8>> {
        if genome_idx >= self.genomes.len() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("invalid genome index {genome_idx}"),
            ));
        }
        let total_bp = self.genomes[genome_idx].ensure_total_bp(&mut self.reader)?;
        if len == 0 || global_start >= total_bp {
            return Ok(Vec::new());
        }
        let clamped_len =
            usize::try_from(total_bp.saturating_sub(global_start).min(len as u64)).unwrap_or(0);
        if clamped_len == 0 {
            return Ok(Vec::new());
        }
        let byte_start = global_start / 2;
        let byte_end = global_start.saturating_add(clamped_len as u64).div_ceil(2);
        let n_bytes = usize::try_from(byte_end.saturating_sub(byte_start)).unwrap_or(0);
        let cache_cap = self.chunk_cache_cap;
        let packed = {
            let g = &mut self.genomes[genome_idx];
            let rdr = g.ensure_dna_reader(&mut self.reader, cache_cap)?;
            rdr.read_elem_range(&mut self.reader, byte_start, n_bytes)?
        };
        let odd_offset = usize::try_from(global_start & 1).unwrap_or(0);
        let dec_len = odd_offset.checked_add(clamped_len).ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "DNA slice length overflow")
        })?;
        let mut decoded = decode_hal_dna(&packed, dec_len);
        if odd_offset == 0 {
            decoded.truncate(clamped_len);
            Ok(decoded)
        } else if odd_offset < decoded.len() {
            Ok(decoded[odd_offset..decoded.len().min(dec_len)].to_vec())
        } else {
            Ok(Vec::new())
        }
    }
}
