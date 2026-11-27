# ADR-002: Streaming Progressive HDF5 HAL Export

## Status
Accepted (updated 2026-11-24 with `4 MiB`-windowed streaming Fitch parsimony and internal parse edges)

## Date
2026-11-22

## Context
During `horizon_guide --progressive --hal` on 20 *Plasmodium* genomes (GitHub Issue #1), assembling all ancestral `HalGenome` segment tables and consensus DNA arrays simultaneously in memory caused peak RSS to exceed 30.7 GiB and trigger OOM kills. Scaling further to 100-way 3 Gbp mammalian alignments would require ~480 GiB if full-genome ancestral DNA or Fitch state-set tables were materialized in RAM.

## Decision
Implement `HalAlignment::write_multi_contig_progressive_hdf5_hal` in `src/hal/hdf5_progressive.rs` (sharing `init_hdf5_root`, `write_hal_genome_group`, and binary record writers with `src/hal/hdf5_write.rs`) to stream progressive HAL v2.2 HDF5 `v0` archives directly to disk:
1. **Single-pass coordinate partition & flat paralogy cycles:** Partition the reference coordinate space `[0, ref_dna_len)` into contiguous, non-overlapping bottom segments once (`partition_progressive_segments` in `src/hal/normalize.rs`), and normalize each query's top segments with unaligned gap fills and flat `chain_head`/`chain_tail` circular paralogy arrays to satisfy `libHal` invariants.
2. **Zero-concatenation DNA streaming:** Stream leaf `DNA_ARRAY` nibbles directly from per-contig slices (`write_hal_dna_contigs_stream`) without concatenating genome-wide DNA vectors.
3. **`4 MiB`-windowed streaming two-pass Fitch parsimony:** Pre-reserve contiguous `DNA_ARRAY` extents for `Anc0..Anc_{n-1}` (`Hdf5V0Builder::reserve_contiguous_dataset`), then slide an even `4 MiB`-base window across the reference. Each window reuses preallocated `n_q` nibble-packed Fitch state-set buffers (`pack_fitch_masks_into`), preserves `N`-scaffold gaps (`src/hal/fitch.rs`), resolves top-down ancestral bases, and writes each window in place via `Hdf5V0Builder::write_at`. Working RAM stays bounded at ~240 MiB even for 100-way 3 Gbp mammalian alignments.
4. **Exact internal parse edges:** Emit `bottom_parse` and `top_parse` indices on all internal ancestral nodes (`Anc0..Anc_{n-2}`) so `libHal` tree traversal works across multi-level progressive hierarchies, and consume `pairwise_mafs` by value before HDF5 streaming begins.

## Alternatives Considered

### Keep in-memory `from_multi_contig_progressive_blocks` + `write_hdf5_hal`
- **Pros:** Single builder for both in-memory inspection and disk export.
- **Cons:** Memory scales as $O(N_{\text{species}} \times (\text{DNA} + \text{segments}))$, exceeding 30 GiB on 20-genome eukaryote cohorts.
- **Rejected:** Retained only for small in-memory round-trip tests (`windowed_progressive_fitch_matches_in_memory_*` asserts byte-for-byte raw `DNA_ARRAY` parity against it).

### Full-genome single-buffer ancestral streaming (`running_counts` + full `anc_dna`)
- **Pros:** Simpler single pass without `reserve_contiguous_dataset` / `write_at`.
- **Cons:** Still allocates $O(N_{\text{species}} \times L_{\text{ref}} / 2)$ for two-pass Fitch state sets across the entire reference (~150–480 GiB on 100-way 3 Gbp mammals).
- **Rejected:** Replaced in Wave C-opt-5 by `4 MiB`-windowed streaming Fitch parsimony with reused `level_sets` buffers.

### Link against C++ `libhdf5` / `libHal`
- **Pros:** Delegates HDF5 chunking and B-tree management to upstream C libraries.
- **Cons:** Breaks `horizon-aligner`'s pure-Rust, zero-C-FFI build portability and complicates static/HPC deployment.
- **Rejected:** Pure-Rust HDF5 `v0` streaming keeps the entire crate free of C/C++ FFI dependencies.

## Consequences
- Peak memory during progressive HAL export stays bounded at `~240 MiB` even at 100-way mammalian genome scale while producing raw `DNA_ARRAY` byte-identical output to in-memory two-pass Fitch parsimony.
- Multi-contig progressive HALs round-trip cleanly through both Horizon's fail-closed pure-Rust HDF5 reader (`hdf5_read.rs`, `hdf5_dataset.rs`) and UCSC `libHal` / `h5ls`.
- For reference-free de novo bottom-up ancestral contigs and multi-level HAL construction without a single hinge frame, see [ADR-004](0004-progressive-ancestral-tba-and-ml-transducers.md) (`horizon_tba --progressive-ancestors --hal`).
