//! Per-genome lazy HDF5 handle and public summary types for [`super::stream_hal::StreamingHalFile`].

use super::hdf5_read::read_hdf5_v0_group_entries_seek;
use super::hdf5_stream::Hdf5DatasetReader;
use super::{decode_hal_seq_index, decode_hal_seq_names};
use std::collections::HashMap;
use std::io::{Read, Seek};

/// Parsed sequence coordinate entry inside a HAL genome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalStreamSeq {
    /// Contig/chromosome name.
    pub name: String,
    /// Global 0-based start offset in the genome's concatenated DNA space.
    pub offset: u64,
    /// Sequence length in base pairs.
    pub length: u64,
    /// First top-segment index belonging to this sequence.
    pub top_offset: u64,
    /// First bottom-segment index belonging to this sequence.
    pub bottom_offset: u64,
}

/// Summary statistics for one genome in a `.hal` file (matching `halStats`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalGenomeStats {
    /// Genome name.
    pub name: String,
    /// Parent genome name (`None` for tree root).
    pub parent: Option<String>,
    /// Number of direct children in the guide tree.
    pub num_children: usize,
    /// Number of chromosomes/scaffolds.
    pub num_sequences: usize,
    /// Total genome length in base pairs.
    pub total_bp: u64,
    /// Number of top segments (excluding sentinel).
    pub num_top_segments: u64,
    /// Number of bottom segments (excluding sentinel).
    pub num_bottom_segments: u64,
}

/// One projected BED interval emitted by [`super::stream_hal::StreamingHalFile::liftover_interval`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalLiftoverHit {
    /// Destination contig/chromosome name.
    pub seq_name: String,
    /// 0-based start on `seq_name`.
    pub start: u64,
    /// 0-based exclusive end on `seq_name`.
    pub end: u64,
    /// Total length of `seq_name`.
    pub seq_len: u64,
    /// Relative orientation (`'+'` or `'-'`).
    pub strand: char,
}

/// Telemetry returned by [`super::stream_hal::StreamingHalFile::project_to_maf`].
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct HalMafStats {
    /// Number of MAF blocks emitted.
    pub blocks_emitted: usize,
    /// Number of reference bases emitted across MAF blocks.
    pub ref_bp_emitted: u64,
    /// Number of target alignment rows emitted across MAF blocks.
    pub target_rows_emitted: usize,
}

const MAX_HAL_META_ARRAY_BYTES: u64 = 64 * 1024 * 1024;

pub(super) struct StreamGenomeHandle {
    pub(super) name: String,
    pub(super) parent: Option<usize>,
    pub(super) children: Vec<usize>,
    pub(super) group_addr: u64,
    header_loaded: bool,
    dna_addr: Option<u64>,
    top_addr: Option<u64>,
    bot_addr: Option<u64>,
    seqidx_addr: Option<u64>,
    seqname_addr: Option<u64>,
    dna_reader: Option<Hdf5DatasetReader>,
    top_reader: Option<Hdf5DatasetReader>,
    bot_reader: Option<Hdf5DatasetReader>,
    /// Lazy compact `SEQIDX_ARRAY` triples `(dna_offset, top_offset, bottom_offset)` (`N + 1` with sentinel).
    pub(super) seq_offsets_loaded: bool,
    pub(super) seq_offsets: Vec<(u64, u64, u64)>,
    /// Lazy full sequence metadata (`SEQNAME_ARRAY` + `HalStreamSeq` + name lookup map).
    pub(super) sequences_loaded: bool,
    pub(super) sequences: Vec<HalStreamSeq>,
    pub(super) seq_name_to_idx: HashMap<String, usize>,
    pub(super) total_bp: Option<u64>,
    pub(super) num_sequences: usize,
    pub(super) num_top_segments: u64,
    pub(super) num_bottom_segments: u64,
    /// Per-genome sequential/galloping segment hints for fast multi-fragment tree hops.
    pub(super) last_top_hint: Option<u64>,
    pub(super) last_bot_hint: Option<u64>,
}

impl StreamGenomeHandle {
    pub(super) fn new_lazy(
        name: String,
        parent: Option<usize>,
        children: Vec<usize>,
        group_addr: u64,
    ) -> Self {
        let missing = group_addr == u64::MAX;
        Self {
            name,
            parent,
            children,
            group_addr,
            header_loaded: missing,
            dna_addr: None,
            top_addr: None,
            bot_addr: None,
            seqidx_addr: None,
            seqname_addr: None,
            dna_reader: None,
            top_reader: None,
            bot_reader: None,
            seq_offsets_loaded: missing,
            seq_offsets: Vec::new(),
            sequences_loaded: missing,
            sequences: Vec::new(),
            seq_name_to_idx: HashMap::new(),
            total_bp: missing.then_some(0),
            num_sequences: 0,
            num_top_segments: 0,
            num_bottom_segments: 0,
            last_top_hint: None,
            last_bot_hint: None,
        }
    }

    /// Lazily load only the genome's HDF5 group directory and dataset dataspace dimensions
    /// (zero payload/chunk reads).
    pub(super) fn ensure_header<R: Read + Seek>(&mut self, r: &mut R) -> std::io::Result<()> {
        if self.header_loaded {
            return Ok(());
        }
        let entries = read_hdf5_v0_group_entries_seek(r, self.group_addr)?.unwrap_or_default();
        let find = |k: &str| {
            entries
                .iter()
                .find(|e| e.name == k)
                .map(|e| e.object_header_addr)
        };
        self.dna_addr = find("DNA_ARRAY");
        self.top_addr = find("TOP_ARRAY");
        self.bot_addr = find("BOTTOM_ARRAY");
        self.seqidx_addr = find("SEQIDX_ARRAY");
        self.seqname_addr = find("SEQNAME_ARRAY");

        self.num_top_segments = if let Some(addr) = self.top_addr {
            super::hdf5_dataset::read_hdf5_dataspace_dims_seek(r, addr)?
                .and_then(|d| d.first().copied())
                .unwrap_or(0)
                .saturating_sub(1)
        } else {
            0
        };
        self.num_bottom_segments = if let Some(addr) = self.bot_addr {
            super::hdf5_dataset::read_hdf5_dataspace_dims_seek(r, addr)?
                .and_then(|d| d.first().copied())
                .unwrap_or(0)
                .saturating_sub(1)
        } else {
            0
        };
        self.num_sequences =
            if let (Some(si_addr), Some(sn_addr)) = (self.seqidx_addr, self.seqname_addr) {
                let idx_count = super::hdf5_dataset::read_hdf5_dataspace_dims_seek(r, si_addr)?
                    .and_then(|d| d.first().copied())
                    .unwrap_or(0)
                    .saturating_sub(1);
                let name_count = super::hdf5_dataset::read_hdf5_dataspace_dims_seek(r, sn_addr)?
                    .and_then(|d| d.first().copied())
                    .unwrap_or(0);
                usize::try_from(idx_count.min(name_count)).unwrap_or(0)
            } else {
                0
            };
        self.header_loaded = true;
        Ok(())
    }

    /// Lazily resolve `total_bp` by reading ONLY the final sentinel element of `SEQIDX_ARRAY`
    /// (zero `SEQNAME_ARRAY` reads and zero full-table decompression for `stats`).
    pub(super) fn ensure_total_bp<R: Read + Seek>(&mut self, r: &mut R) -> std::io::Result<u64> {
        if let Some(bp) = self.total_bp {
            return Ok(bp);
        }
        self.ensure_header(r)?;
        let mut bp = 0u64;
        if let Some(si_addr) = self.seqidx_addr
            && let Some(mut si_rdr) = Hdf5DatasetReader::open_with_cache_cap(r, si_addr, 2)?
        {
            if si_rdr.elem_size() != 24 || si_rdr.total_bytes() > MAX_HAL_META_ARRAY_BYTES {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "HAL genome `{}` SEQIDX metadata exceeds safety limits or invalid elem_size",
                        self.name
                    ),
                ));
            }
            let total_elems = si_rdr.total_elems();
            if total_elems > 0 {
                let last_bytes = si_rdr.read_elem_range(r, total_elems - 1, 1)?;
                if last_bytes.len() >= 8 {
                    bp = u64::from_le_bytes(last_bytes[0..8].try_into().unwrap_or([0; 8]));
                }
            }
        } else if let Some(addr) = self.dna_addr
            && let Some(dims) = super::hdf5_dataset::read_hdf5_dataspace_dims_seek(r, addr)?
        {
            bp = dims.first().copied().unwrap_or(0).saturating_mul(2);
        }
        self.total_bp = Some(bp);
        Ok(bp)
    }

    /// Lazily load compact `SEQIDX_ARRAY` triples `(dna_offset, top_offset, bot_offset)` without
    /// touching `SEQNAME_ARRAY` (used by intermediate ancestors during multi-hop liftover).
    pub(super) fn ensure_seq_offsets<R: Read + Seek>(&mut self, r: &mut R) -> std::io::Result<()> {
        if self.seq_offsets_loaded {
            return Ok(());
        }
        self.ensure_header(r)?;
        if let Some(si_addr) = self.seqidx_addr
            && let Some(mut si_rdr) = Hdf5DatasetReader::open_with_cache_cap(r, si_addr, 4)?
        {
            if si_rdr.elem_size() != 24 || si_rdr.total_bytes() > MAX_HAL_META_ARRAY_BYTES {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "HAL genome `{}` SEQIDX metadata exceeds safety limits or invalid elem_size",
                        self.name
                    ),
                ));
            }
            let idx_n = usize::try_from(si_rdr.total_elems()).unwrap_or(0);
            let idx_bytes = si_rdr.read_elem_range(r, 0, idx_n)?;
            let seqidx = decode_hal_seq_index(&idx_bytes);
            if let Some(&(end_bp, _, _)) = seqidx.last() {
                self.total_bp = Some(end_bp);
            }
            self.seq_offsets = seqidx;
        } else {
            let _ = self.ensure_total_bp(r)?;
        }
        self.seq_offsets_loaded = true;
        Ok(())
    }

    /// Lazily load full `HalStreamSeq` records + `seq_name_to_idx` map (only for genomes whose
    /// sequence names are queried: `src`, `dst`, `ref`, `targets`, or `--chroms`).
    pub(super) fn ensure_sequences<R: Read + Seek>(&mut self, r: &mut R) -> std::io::Result<()> {
        if self.sequences_loaded {
            return Ok(());
        }
        self.ensure_seq_offsets(r)?;
        if let Some(sn_addr) = self.seqname_addr
            && !self.seq_offsets.is_empty()
            && let Some(mut sn_rdr) = Hdf5DatasetReader::open_with_cache_cap(r, sn_addr, 4)?
        {
            if sn_rdr.total_bytes() > MAX_HAL_META_ARRAY_BYTES {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "HAL genome `{}` SEQNAME metadata exceeds safety limits",
                        self.name
                    ),
                ));
            }
            let name_n = usize::try_from(sn_rdr.total_elems()).unwrap_or(0);
            let name_bytes = sn_rdr.read_elem_range(r, 0, name_n)?;
            let seqnames = decode_hal_seq_names(&name_bytes, sn_rdr.elem_size());
            let n = seqnames.len().min(self.seq_offsets.len().saturating_sub(1));
            self.sequences.reserve(n);
            self.seq_name_to_idx.reserve(n);
            for i in 0..n {
                let (off, top_off, bot_off) = self.seq_offsets[i];
                let next_off = self.seq_offsets[i + 1].0;
                let len = next_off.saturating_sub(off);
                self.seq_name_to_idx.insert(seqnames[i].clone(), i);
                self.sequences.push(HalStreamSeq {
                    name: seqnames[i].clone(),
                    offset: off,
                    length: len,
                    top_offset: top_off,
                    bottom_offset: bot_off,
                });
            }
            self.num_sequences = self.sequences.len();
        }
        self.sequences_loaded = true;
        Ok(())
    }

    pub(super) fn set_cache_cap(&mut self, cap: usize) {
        if let Some(rdr) = &mut self.dna_reader {
            rdr.set_cache_cap(cap);
        }
        if let Some(rdr) = &mut self.top_reader {
            rdr.set_cache_cap(cap);
        }
        if let Some(rdr) = &mut self.bot_reader {
            rdr.set_cache_cap(cap);
        }
    }

    fn open_typed_reader<'a, R: Read + Seek>(
        r: &mut R,
        genome_name: &str,
        dset_name: &str,
        addr: Option<u64>,
        expected_es: usize,
        cache_cap: usize,
        slot: &'a mut Option<Hdf5DatasetReader>,
    ) -> std::io::Result<&'a mut Hdf5DatasetReader> {
        if slot.is_none() {
            let a = addr.ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    format!("genome `{genome_name}` missing {dset_name}"),
                )
            })?;
            let rdr =
                Hdf5DatasetReader::open_with_cache_cap(r, a, cache_cap)?.ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("unable to open {dset_name} for `{genome_name}`"),
                    )
                })?;
            if rdr.elem_size() != expected_es {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "genome `{genome_name}` {dset_name} elem_size={} (expected {expected_es})",
                        rdr.elem_size()
                    ),
                ));
            }
            *slot = Some(rdr);
        }
        slot.as_mut().ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("unable to open {dset_name} for `{genome_name}`"),
            )
        })
    }

    pub(super) fn ensure_dna_reader<R: Read + Seek>(
        &mut self,
        r: &mut R,
        cache_cap: usize,
    ) -> std::io::Result<&mut Hdf5DatasetReader> {
        self.ensure_header(r)?;
        Self::open_typed_reader(
            r,
            &self.name,
            "DNA_ARRAY",
            self.dna_addr,
            1,
            cache_cap,
            &mut self.dna_reader,
        )
    }

    pub(super) fn ensure_top_reader<R: Read + Seek>(
        &mut self,
        r: &mut R,
        cache_cap: usize,
    ) -> std::io::Result<&mut Hdf5DatasetReader> {
        self.ensure_header(r)?;
        Self::open_typed_reader(
            r,
            &self.name,
            "TOP_ARRAY",
            self.top_addr,
            33,
            cache_cap,
            &mut self.top_reader,
        )
    }

    pub(super) fn ensure_bot_reader<R: Read + Seek>(
        &mut self,
        r: &mut R,
        cache_cap: usize,
    ) -> std::io::Result<&mut Hdf5DatasetReader> {
        self.ensure_header(r)?;
        let expected_es = 24usize
            .checked_add(self.children.len().saturating_mul(9))
            .unwrap_or(24);
        Self::open_typed_reader(
            r,
            &self.name,
            "BOTTOM_ARRAY",
            self.bot_addr,
            expected_es,
            cache_cap,
            &mut self.bot_reader,
        )
    }

    /// Locate the sequence index containing `global_pos` using compact `seq_offsets` alone.
    pub(super) fn seq_idx_for_pos_offsets(&self, global_pos: u64) -> Option<usize> {
        let n = self.seq_offsets.len().saturating_sub(1);
        if n == 0 {
            return None;
        }
        let idx = self.seq_offsets[..n]
            .partition_point(|&(off, _, _)| off <= global_pos)
            .checked_sub(1)?;
        let start = self.seq_offsets[idx].0;
        let end = self.seq_offsets[idx + 1].0;
        (global_pos >= start && global_pos < end).then_some(idx)
    }

    /// Return `[lo_seg, hi_seg)` segment index bounds for the sequence containing `global_pos`,
    /// using `SEQIDX_ARRAY` (`top_offset` if `is_top`, else `bottom_offset`).
    pub(super) fn seq_seg_bounds(
        &self,
        global_pos: u64,
        is_top: bool,
        total_segs: u64,
    ) -> Option<(u64, u64)> {
        let idx = self.seq_idx_for_pos_offsets(global_pos)?;
        let (lo, hi) = if is_top {
            (self.seq_offsets[idx].1, self.seq_offsets[idx + 1].1)
        } else {
            (self.seq_offsets[idx].2, self.seq_offsets[idx + 1].2)
        };
        (lo < hi && hi <= total_segs).then_some((lo, hi))
    }

    pub(super) fn map_global_to_seq(
        &self,
        global_start: u64,
        global_end: u64,
    ) -> Option<(&HalStreamSeq, u64, u64)> {
        if self.sequences.is_empty() || global_end <= global_start {
            return None;
        }
        let idx = self
            .sequences
            .partition_point(|s| s.offset.saturating_add(s.length) <= global_start);
        let seq = self.sequences.get(idx)?;
        if global_start < seq.offset || global_start >= seq.offset.saturating_add(seq.length) {
            return None;
        }
        let local_s = global_start - seq.offset;
        let local_e = global_end.min(seq.offset.saturating_add(seq.length)) - seq.offset;
        (local_e > local_s).then_some((seq, local_s, local_e))
    }
}
