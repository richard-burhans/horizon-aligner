# Horizon Suite Architecture

This document establishes the conceptual map of the major abstractions in `horizon-aligner` and the contracts between them.

## Tier Map

The architecture of `horizon-aligner` follows a unidirectional dependency flow from core algorithms and streaming I/O up to output formatters and the 9 CLI drivers.

| Layer / Package | Abstraction / Responsibility | Key Public API |
|-----------------|------------------------------|----------------|
| **CLI Drivers**<br>`src/bin/` | 9 suite executables: `horizon_wga` (pairwise & `--graph` DAG/POA & `--hal-info`), `horizon_guide` (progressive WGA + HAL export), `horizon_merge` (multi-species MAF/rGFA merge), `horizon_gfa` (streaming stats/extract), `horizon_lift` (BED/GFF3 liftover), `horizon_clean` (BOCPD+MDL junction cleaner), `horizon_resolve` (GFA walk e-graph resolver), `horizon_eval` (MAF evaluator), `bench_chr1`. | Executable entry points |
| **I/O, GFA & HAL**<br>`src/formats.rs`<br>`src/gfa/mod.rs`<br>`src/hal/` | Serialization and streaming deserialization. `formats.rs` emits MAF v1 and GFA v1.1. `gfa/mod.rs` streams gzipped/plain/stdin GFA v1.0/v1.1 (`MultiGzDecoder` + bounded `SpoolableGfaInput`). `hal/` provides `HalAlignment`, streaming `Read+Seek` HDF5 `v0`/`v1` decoding (`hdf5_read.rs`), and `libhdf5`-compatible HDF5 `v0` writing (`hdf5_write.rs`). | `format_maf_blocks`, `format_gfa_blocks`, `read_fasta_contigs`, `stream_gfa`, `HalAlignment`, `to_hdf5_bytes` |
| **Multi-species Merger**<br>`src/merge.rs` | `O(N log N)` sweep-line exact-coordinate multi-species MAF merger and query-insertion gap-column threading into MAF and `rGFA`. | `merge_pairwise_mafs` |
| **E-Graph Resolution**<br>`src/egraph.rs` | Congruence-closed e-graph with `MdlCostModel` for identifying and collapsing tandem duplications and nested inversions in both pairwise LCBs (`# EG:Z:`) and GFA walks. | `resolve_lcb_duplications`, `resolve_alignment_duplications`, `MdlCostModel` |
| **Alignment (MDL+BOCPD & DAG)**<br>`src/mdl/` | Pairwise banded Pair-HMM DP + X-drop (`dp.rs`), adaptive `HmmParams` shrinkage (`mod.rs`), and seeded topological DAG/POA alignment (`dag.rs`) with Horizon Bubble abstention. | `align_contig_pair`, `banded_mdl_gap`, `mdl_xdrop_extend`, `align_seq_to_dag_mdl`, `align_seq_to_dag_seeded` |
| **Seeding & Chaining**<br>`src/seeding/` | SDUST repeat masking, `FlatSeedIndex` (`flat_index.rs`, `12 B/entry` + 16-bit prefix buckets), Rayon-parallel `LASTZ12` transition seeding + 6-frame biochemical `AA5` (`translated.rs`), LIS chaining, and `BTreeSet`-indexed LCB clustering. | `FlatSeedIndex`, `prepare_contig_index`, `lis_chain`, `local_w9_subchain` |
| **Core Primitives**<br>`src/bocpd.rs`<br>`src/coord.rs` | Beta-Bernoulli BOCPD compositional segmentation (`bocpd_gc_profile`) and 0-based half-open / UCSC reverse-strand interval arithmetic (`Interval`). | `bocpd_gc_profile`, `Interval` |

## Cross-Cutting Contracts

- **Dependency Direction:** `src/mdl/` imports `src/seeding/` and `src/bocpd.rs`, never vice versa. `src/gfa/` and `src/coord.rs` sit at the I/O and coordinate base; `src/merge.rs` and `src/formats.rs` consume lower layers without circular dependencies.
- **Coordinate Systems:** All internal intervals use 0-based, half-open `[start, end)` coordinates (`coord::Interval`). Reverse-strand coordinates use UCSC forward-length-relative conventions unless explicitly translated for DP.
- **Source of Truth for Alignments:** `ContigAlignment` in `mdl/mod.rs` is the definitive pairwise representation; `formats.rs` projects it to MAF or GFA without mutating topology.

## Known Asymmetries

- `horizon_clean` currently severs only unsupported (`walk_support == 0`) discordant edges into `HB:i:1` Horizon Bubbles while preserving walk-supported links unchanged.
- `horizon_resolve` rewrites walk sequences (`P`/`W`) for tandem collapses/double-inversion cancellation without regenerating `W` coordinate bounds or `L` links.
- `hal/hdf5_write.rs` emits `libhdf5`-compatible HDF5 `v0` archives with multi-level B-trees (`>256` entries), whereas `hal/hdf5_read.rs` supports both `v0` and `v1` object headers.
