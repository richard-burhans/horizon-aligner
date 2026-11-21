# Horizon Suite Architecture

This document establishes the conceptual map of the major abstractions in `horizon-aligner` and the contracts between them.

## Tier Map

The architecture of `horizon-aligner` follows a unidirectional dependency flow from core algorithms and streaming I/O up to output formatters and the 12 CLI drivers (`11` suite tools + `bench_chr1`).

| Layer / Package | Abstraction / Responsibility | Key Public API |
|-----------------|------------------------------|----------------|
| **CLI Drivers**<br>`src/bin/` | `11` suite executables + `bench_chr1`: `horizon` (subcommand dispatcher), `horizon_wga` (pairwise & `--graph` DAG/POA & `--hal-info`), `horizon_guide` (MinHash-hinge & parallel UPGMA progressive WGA + HAL export in `main.rs`, `rescue.rs`), `horizon_merge` (streaming multi-species MAF/rGFA merge), `horizon_gfa` (streaming stats/extract), `horizon_lift` (BED/GFF3 liftover in `main.rs`, `lift.rs`), `horizon_clean` (BOCPD+MDL junction cleaner), `horizon_resolve` (modular `main.rs`, `snarl.rs`, `resolve.rs` GFA `PackedGraph` snarl + e-graph resolver), `horizon_filter` (streaming SDUST+MDL PAF filter & netter), `horizon_locus` (regional multi-FASTA locus aligner & similarity matrix), `horizon_eval` (MAF evaluator), `bench_chr1`. | Executable entry points |
| **I/O, PAF, GFA & HAL**<br>`src/io_util.rs`<br>`src/paf.rs`<br>`src/formats.rs`<br>`src/gfa/mod.rs`<br>`src/hal/` | Serialization and transparent streaming I/O (`io_util.rs` handles plain, `.gz`/`.bgz`, `.zst`/`.zstd`, and `stdin`/`stdout`). `paf.rs` streams PAF + CIGAR (`cg:Z:`). `formats.rs` emits MAF v1 and GFA v1.1 (`format_maf_blocks_to_writer`, `format_gfa_blocks_to_writer`, `materialize_run_group`). `gfa/mod.rs` streams GFA v1.0/v1.1 (`SpoolableGfaInput`). `hal/` provides `HalAlignment`, streaming `Read+Seek` HDF5 `v0`/`v1` decoding (`hdf5_read.rs`), `libhdf5`-compatible HDF5 `v0` writing (`hdf5_write.rs`), and unit/integration tests (`tests_hdf5.rs`, `tests_model.rs`). | `open_text_reader`, `open_text_writer`, `stream_paf_records`, `format_maf_blocks`, `format_gfa_blocks`, `read_fasta_contigs`, `stream_gfa`, `HalAlignment`, `to_hdf5_bytes` |
| **Multi-species Merger & Netting**<br>`src/merge.rs`<br>`src/netting.rs` | `O(N log N)` sweep-line exact-coordinate multi-species MAF merger, interval netting, and streaming query-insertion gap-column threading into MAF and `rGFA` (`merge_pairwise_maf_paths_to_writers`). | `merge_pairwise_mafs`, `merge_pairwise_maf_paths_to_writers` |
| **E-Graph Resolution**<br>`src/egraph.rs` | Congruence-closed e-graph with `MdlCostModel` for identifying and collapsing tandem duplications and nested inversions in both pairwise LCBs (`# EG:Z:`) and GFA snarl/walk resolution. | `resolve_lcb_duplications`, `resolve_alignment_duplications`, `MdlCostModel` |
| **Alignment (MDL+BOCPD & DAG)**<br>`src/mdl/` | Pairwise banded Pair-HMM DP + X-drop (`dp.rs`, `dp_fast.rs` diagonal/narrow-Viterbi fast paths, `xdrop.rs`, `runs.rs`, `fallback.rs`), Conformal Risk Control calibration (`crc.rs`), adaptive `HmmParams` shrinkage (`mod.rs`), and seeded topological DAG/POA alignment (`dag.rs`, `dag_seed.rs` with `MASK_W10`, `MASK_RY20`, and `aa5_keys`) with Horizon Bubble abstention. | `align_contig_pair`, `banded_mdl_gap`, `mdl_xdrop_extend`, `calibrate_crc_threshold`, `align_seq_to_dag_mdl`, `align_seq_to_dag_seeded` |
| **Seeding & Chaining**<br>`src/seeding/` | SDUST repeat masking, `FlatSeedIndex` (`flat_index.rs`, `12 B/entry` + 16-bit prefix buckets), `ContigIndex::minhash_sketch` / `shared_sketch_count` for fast hinge selection, Rayon-parallel `LASTZ12` transition seeding + 6-frame biochemical `AA5` (`translated.rs`), LIS chaining, and `BTreeSet`-indexed LCB clustering. | `FlatSeedIndex`, `prepare_contig_index`, `shared_sketch_count`, `lis_chain`, `local_w9_subchain` |
| **Core Primitives**<br>`src/bocpd.rs`<br>`src/coord.rs` | Beta-Bernoulli BOCPD compositional segmentation (`bocpd_gc_profile`) and 0-based half-open / UCSC reverse-strand interval arithmetic (`Interval`). | `bocpd_gc_profile`, `Interval` |

## Cross-Cutting Contracts

- **Dependency Direction:** `src/mdl/` imports `src/seeding/` and `src/bocpd.rs`, never vice versa. `src/gfa/` and `src/coord.rs` sit at the I/O and coordinate base; `src/merge.rs` and `src/formats.rs` consume lower layers without circular dependencies.
- **Coordinate Systems:** All internal intervals use 0-based, half-open `[start, end)` coordinates (`coord::Interval`). Reverse-strand coordinates use UCSC forward-length-relative conventions unless explicitly translated for DP.
- **Source of Truth for Alignments:** `ContigAlignment` in `mdl/mod.rs` is the definitive pairwise representation; `formats.rs` projects it to MAF or GFA without mutating topology.

## Known Asymmetries

- `horizon_clean` severs low-support discordant edges (`walk_support <= 1` by default, or `0` with `--prune-unsupported`) into `HB:i:1` Horizon Bubbles while preserving higher-support links unchanged.
- `horizon_resolve` (`src/bin/horizon_resolve/{main.rs,snarl.rs,resolve.rs}`) decomposes graphs via `u32`-interned `PackedGraph` + $O(V + E)$ Tarjan biedged snarls, estimates window divergence via a topological `[0.02, 0.15]` heuristic (unless `--divergence` is passed), streams `S`-line sequences intact on rewrite, rewrites `P`/`W` walks, updates `W`-line `end` bounds, and prunes/synthesizes `L` links (assuming `0M` overlaps).
- `horizon` dispatcher (`src/bin/horizon.rs`) forwards `wga`, `merge`, `guide`, `gfa`, `lift`, `clean`, `resolve`, `filter`, `locus`, `eval`, and `bench-chr1`.
- `hal/hdf5_write.rs` emits `libhdf5`-compatible HDF5 `v0` archives with multi-level B-trees (`>256` entries), whereas `hal/hdf5_read.rs` supports both `v0` and `v1` object headers.
