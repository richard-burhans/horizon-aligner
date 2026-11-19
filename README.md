# horizon-aligner

A pure-Rust **information-theoretic genome alignment and pangenome suite** combining **strand-symmetric spaced-seed indexing (`FlatSeedIndex`)**, **Bayesian Online Change-Point Detection (BOCPD)** compositional segmentation, **Minimum Description Length (MDL)** banded & DAG dynamic programming with explicit epistemic abstention (**Horizon Bubbles**), **MDL-calibrated e-graph paralog resolution**, and **streaming `GFA v1.0/v1.1` + pure-Rust `HAL`/HDF5** interoperability.

## Key Capabilities

- **Pairwise & DAG/POA MDL alignment (`horizon_wga`):** Strand-symmetric `KP-14` (`11-of-15`) indexing in `FlatSeedIndex` (`12 bytes/entry` + 16-bit prefix buckets) with SDUST repeat masking, Rayon-parallel `LASTZ12`/`LASTZ11` transition seeding, `w=10`/`w=9`, 1-bit `R/Y`, and 6-frame `AA5` fallback seeding, driving flat `BandGrid` Pair-HMM forward-backward alignment (`src/mdl/dp.rs`) and seeded DAG/POA alignment (`horizon_wga --graph`, `src/mdl/dag.rs`).
- **Calibrated epistemic abstention (Horizon Bubbles):** Intervals falling below the search-space-scaled MDL floor (`0.55 log2(N*M)` + conformal `10.0`-bit floor) or `min_post_prob` abstain explicitly as `# HB:Z:` / `HB:i:1` Horizon Bubbles rather than forcing spurious base pairs.
- **Streaming `12 GB+` GFA pangenome toolkit (`src/gfa/mod.rs`, `src/coord.rs`):**
  - `horizon_gfa`: Single-pass BGZF/gzip (`MultiGzDecoder`) `stats` (`99.97 s`, `97.0 MB` RSS across `81.4M` segments / `3.29 Gbp`) and two-pass chromosome `extract`.
  - `horizon_lift`: Graph-aware BED/GFF3 coordinate liftover with `--bridge-bp` micro-bubble coalescing and UCSC reverse-strand arithmetic.
  - `horizon_clean`: Two-pass `S`-flank-indexed BOCPD+MDL junction scoring that severs unsupported (`sup=0`) discordant links into Horizon Bubbles.
  - `horizon_resolve`: Standalone MDL-calibrated e-graph resolution (`MdlCostModel`) over GFA `P`/`W` walks (`--collapse-dups`, `--tag-walks`).
- **Multi-species merging & progressive guide-tree WGA (`horizon_merge`, `horizon_guide`):** `O(N log N)` coordinate sweep-line MAF/rGFA merging with best-score interval netting (`src/merge.rs`) and iterative ancestral consensus rescue.
- **Pure-Rust HDF5 / HAL reader & writer (`src/hal/`):** Streaming `Read + Seek` HDF5 `v0`/`v1` parser (`hdf5_read.rs`) and `libhdf5`-compatible HDF5 `v0` writer (`hdf5_write.rs`) with zero `libhdf5` dependency.

## Repository Structure

```text
src/
├── lib.rs              # Public API re-exports and core unit tests
├── bocpd.rs            # BOCPD compositional segmentation & MDL pair scoring
├── coord.rs            # 0-based half-open & UCSC reverse-strand coordinate arithmetic
├── egraph.rs           # Congruence-closed e-graph & MdlCostModel paralog resolution
├── formats.rs          # MAF v1, GFA v1.1, FASTA I/O, and e-graph paralog formatting
├── merge.rs            # Multi-species reference-guided MAF/rGFA sweep-line merger
├── gfa/
│   └── mod.rs          # Streaming GFA v1.0/v1.1 reader (MultiGzDecoder + SpoolableGfaInput)
├── seeding/
│   ├── mod.rs          # SDUST masking, multi-shape seeding, LIS & LCB clustering
│   ├── flat_index.rs   # Compact FlatSeedIndex with 16-bit prefix buckets
│   └── translated.rs   # 6-frame biochemical AA5 seeding & HSP-gated LASTZ12
├── mdl/
│   ├── mod.rs          # Contig-pair alignment pipeline & adaptive HmmParams shrinkage
│   ├── dp.rs           # Flat BandGrid Pair-HMM Forward-Backward & X-drop extension
│   └── dag.rs          # Seeded topological DAG/POA MDL alignment (--graph)
├── hal/
│   ├── mod.rs          # In-memory HalAlignment model, paralog chaining, MAF export
│   ├── hdf5_read.rs    # Streaming Read+Seek HDF5 v0/v1 & HAL binary array decoder
│   └── hdf5_write.rs   # Pure-Rust HDF5 v0 writer (multi-level B-trees)
└── bin/
    ├── horizon_wga.rs      # Pairwise WGA, DAG/POA (--graph), and --hal-info CLI
    ├── horizon_merge.rs    # Multi-species pairwise MAF -> MAF + rGFA merger
    ├── horizon_guide.rs    # Progressive guide-tree / central-hinge WGA + HAL export
    ├── horizon_gfa.rs      # Streaming GFA stats & two-pass chromosome extract
    ├── horizon_lift.rs     # Graph-aware BED/GFF3 coordinate liftover
    ├── horizon_clean.rs    # BOCPD+MDL GFA junction scorer & bubble cutter
    ├── horizon_resolve.rs  # Standalone GFA walk e-graph duplication resolver
    ├── horizon_eval.rs     # Pairwise MAF recall/precision evaluator
    └── bench_chr1.rs       # Pairwise chromosome benchmark harness
```

## Build & Test

Requires Rust (`edition = "2024"`):

```bash
cargo build --release --bins
cargo test --all-features
cargo clippy --all-targets --all-features -- -D warnings
cargo fmt --check
```

To run the ignored real-`.hal` integration test against a downloaded HAL file (e.g., Ensembl `sars-cov-2.hal`):

```bash
HORIZON_TEST_HAL=/path/to/sars-cov-2.hal \
  cargo test real_hal_decode_from_env -- --ignored --nocapture
```

## CLI Quick Reference

```bash
# Pairwise WGA (MAF + GFA + e-graph paralog resolution)
target/release/horizon_wga target.fa query.fa 0 --maf out.maf --gfa out.gfa --resolve-dups

# Sequence-to-graph DAG/POA alignment
target/release/horizon_wga --graph target.gfa.gz query.fa --maf out.maf

# Streaming 12 GB+ GFA stats & chromosome extraction
target/release/horizon_gfa stats hprc-v1.0-mc-grch38.gfa.gz
target/release/horizon_gfa extract hprc-v1.0-mc-grch38.gfa.gz --chrom GRCh38#0#chr6 --out chr6.gfa

# Graph-aware BED/GFF3 liftover, junction cleaning, and e-graph walk resolution
target/release/horizon_lift hprc.gfa.gz regions.bed --from GRCh38#0#chr22 --to CHM13#0#chr22 --bridge-bp 10
target/release/horizon_clean --gfa in.gfa.gz --min-cp-bits 12.0 > clean.gfa
target/release/horizon_resolve in.gfa.gz --collapse-dups --tag-walks > resolved.gfa

# Multi-species sweep-line MAF/rGFA merge & progressive guide-tree WGA
target/release/horizon_merge ref_species out.maf out.gfa pair1.maf pair2.maf
target/release/horizon_guide --progressive --rbest --hal out.hal tree.nwk genomes.tsv
```

## Performance Highlights

- **Pairwise WGA benchmarks (`scripts/benchmarks/run_all_evals.py`):**
  - `multiz100way` (`hg38` vs `panTro4` chrM): **94.4%** recall / **97.4%** precision (`0.18s`).
  - `multiz100way` (`hg38` vs `mm10` chrM): **76.4%** recall / **93.6%** precision (`0.16s`).
  - `cactus447way` (`hg38` vs `Macaca_mulatta`): **90.8%** recall / **86.9%** precision (`0.39s`).
  - `multiz470way` (`hg38` vs `HLcalPym1`): **85.7%** recall / **87.7%** precision (`0.41s`).
  - `cactus447way` (`hg38` vs `Equus_caballus`): **54.0%** recall / **72.4%** precision (`0.24s`).
- **12 GB HPRC v1.0 Minigraph-Cactus Graph & Whole-Chromosome Scale:**
  - Full-graph `horizon_gfa stats` (`81.4M` segments, `3.29 Gbp`, `113.0M` links): **99.97 s**, **97.0 MB** peak RSS.
  - `horizon_lift --bridge-bp 10` (`500` `GRCh38.chr22` 1-kb windows → `CHM13.chr22`): **182.67 s**, **619.7 MB** peak RSS.
  - `horizon_wga --resolve-dups` (`GRCh38.chr22` vs `CHM13.chr22`): **116.58 s**, **3,336.4 MB** peak RSS (`95.6%` of `CHM13.chr22` q-arm).
  - `horizon_wga --resolve-dups` (`GRCh38.chr6` `170.8 Mbp` vs `CHM13.chr6` `110.7 Mbp`): **451.35 s** (`7m 31s`, **2.35×** faster than `1,062.45 s` baseline via 16-bit prefix-bucketed `FlatSeedIndex` + Rayon `LASTZ12` transition seeding), **10.27 GB** peak RSS.

## Documentation

- [Technical Manual](docs/manual.md) — Processing stages, MDL/BOCPD theory, and CLI usage for all suite binaries.
- [Architecture Map](ARCHITECTURE.md) — Tier map, module contracts, and dependency flow.
- [Pangenome Suite](docs/PANGENOME_SUITE.md) — Status and design of `horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_resolve`, and `horizon_wga --graph`.
- [Benchmark Datasets](docs/BENCHMARK_DATASETS.md) & [Benchmark Results](docs/BENCHMARK_RESULTS.md) — Full UCSC, Ensembl, and HPRC v1.0 evaluation tables.
- [Manuscript](manuscript/horizon_paper.md) & [Technical Companion](manuscript/horizon_companion.md) — Main paper and mathematical derivations.
- [Guides](manuscript/guide_making_alignments.md): [MDL vs HOXD70](manuscript/guide_mdl_vs_hoxd70.md) | [Pure-Rust HAL](manuscript/guide_pure_rust_hal.md) | [Horizon Bubbles](manuscript/guide_horizon_bubbles.md)
