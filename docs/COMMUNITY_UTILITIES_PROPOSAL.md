# Community Utilities Proposal: Large-Scale Alignment & Pangenome Tooling

*Grounded in the project's alignment and pangenome literature corpus (`references.bib`, `genome-alignment-synteny-duplications`, and `bioinformatics-tooling-benchmarking`), mapped onto existing primitives in `horizon-aligner`, and **adversarially triaged in Round 8** (`reports/rounds/round-8/`).*

---

## Executive Summary & Round 8 Triage Verdict

Across recent large-scale whole-genome alignment (WGA) and pangenome graph papers—including HPRC (`liao_2023_draft`), Minigraph-Cactus (`hickey_2023_pangenome`), PGGB (`garrison_2024_building`), ODGI (`guarracino_2022_odgi`), comparative pangenome benchmarks (`andreace_2023_comparing`, `kopalli_2025_toward`), inversion & snarl theory (`romain_2026_inversions-pangenome-graphs`, `sena_2026_bubble-subgraphs-spqr`), and maximal-match diagnostics (`shivakumar_2025_mumemto`)—researchers repeatedly report dependency walls (`libhdf5` / C++ `halTools`), multi-tool glue pipelines, and lack of lightweight QC.

**Round 8 Adversarial Triage** evaluated six candidate utilities across `correctness`, `prior-art`, `provenance`, and `reproduction` lenses:

| Utility | Round 8 Verdict | Decision | Summary of Adjudicated Findings |
|---|---|---|---|
| `horizon hal` | **NARROWED (Worth Building — P1)** | **Build** | Only non-duplicative value vs `halTools`/`hal2vg` is zero-`libhdf5` static-binary packaging. Requires removing `hdf5_read.rs`'s `256 MiB` dataset cap via streaming chunk-slice reads and extending `HalAlignment::to_maf_for_genome` beyond 1 tree level. Must clarify relationship to `horizon lift` (GFA) vs `horizon hal liftover` (HAL). |
| `horizon gfa-qc` | **NARROWED (Conditional — P2)** | **Build only if scoped to low-RSS gaps** | `horizon_gfa stats` (`96.9 MB` RSS) is a stateless 5-counter pass skipping `C`/`J` lines; whole-genome (`81.4M` segments) components + per-node depth/AF require ~6–8 GB RAM, and `odgi`/`panacus`/`gfatools` already cover depth/AF/components. FASTA lengths cannot replace `kopalli_2025_toward`'s `minimap2`+`Liftoff` gene checks. |
| `horizon gfa-prune` | **NARROWED (Low Priority — P3)** | **Defer** | `vg clip` (`-d`, `-D`, `-m`) already prunes GFAs without `XG` (`liao_2023_draft` ran `vcfbub` on VCFs, not GFAs). `horizon_clean.rs` counts raw traversals rather than haplotype AF, splits walks rather than bridging pruned segments, projects to ~35+ GB RSS at `81.4M` segments, and cannot detect `--max-del-bp` shortcuts in 2 streaming passes. |
| `horizon inversions` | **REFUTED** | **Do not build** | `PackedGraph` (`horizon_resolve.rs`) stores no adjacency list, snarl tree, or segment sequences (only counts Tarjan bridges and slices 1D walk windows, projecting to `~38.9 GB` RSS at `81.4M` segments); `FlatSeedIndex::intersect` conflates fwd/rc hits and lacks banded extension; `INVPG_annot` already solves this; `sena_2026` is an unreviewed preprint (linear-time ultrabubbles, not all-snarl enumeration). |
| `horizon eval-null` | **REFUTED** | **Do not build as CLI** | `calibrate_null_crc_from_seqs` returns only a scalar `64`-bp reversed-window threshold, `horizon_eval.rs` materializes all aligned pairs in `HashSet<(usize, usize)>` (unscalable to whole genomes), `quinlan-hall_2010_bedtools` is not a WGA ground-truth critique, and `tantan`/`LASTEX` + `Alignathon`/`mafComparator` + `bedtools` already cover this. Keep biological CRC calibration internal (`R7-CRC-BIO`). |
| `horizon synteny` | **REFUTED** | **Do not build** | Dense 1-bp `FlatSeedIndex` has a hard `4 Gbp` `u32` limit and peaks at `~93 GB` RAM on a single `3.1 Gbp` human genome; `lis_chain` is single-contig monotonic; `sdust_mask` misses interspersed LTRs; Mumemto's `707–800 GB` figure was 8-thread total (`139 GB` serial); pairwise seeding cannot replace N-way multi-MUMs or `minimap2`+`SyRI`/`RagTag`/`Panagram`. |

---

## Utilities Worth Building (Narrowed Scope)

### 1. `horizon hal` (`P1`) — Pure-Rust Streaming HAL Toolkit
- **Community Value:** Eliminates C++ `halTools` + `libhdf5` build/deployment friction by shipping a single static binary for `stats`, `to-maf`, `liftover` (distinct from GFA-based `horizon lift`), and `to-gfa`.
- **Host HAL Fixtures Available:** `/home/rcb112/project/claude/bench/data/horizon_hal_bench/4-t2t-orangs-mc-2023v2.full.hal` (`5.53 GiB`) and `jarvis_lab_birds.hal` (`56.43 GiB`) pass `horizon_wga --hal-info`.
- **Required Engineering Before Shipping:**
  1. Remove `MAX_CHUNK_DECOMPRESSED_BYTES` (`256 MiB`) full-dataset buffering in `src/hal/hdf5_read.rs` in favor of streaming chunk-slice reads.
  2. Implement an HDF5-to-HAL segment/tree streaming reader (beyond the current in-memory single-level `HalAlignment::to_maf_for_genome` in `src/hal/mod.rs`).
  3. Define explicit GFA v1.1 `W`-line / `0M` link conventions for `to-gfa` (aligning with `vg`/`odgi` dialect assumptions).
- **Acceptance Criteria:** Interval-identical `stats`, `to-maf`, and `liftover` output compared against `halTools` on `4-t2t-orangs-mc-2023v2.full.hal` and `jarvis_lab_birds.hal` with zero `libhdf5` runtime dependencies.

### 2. `horizon gfa-qc` (`P2`, Conditional) — Lightweight GFA Path & Component Summary
- **Scope Adjustment:** Rather than duplicating `odgi` (`stats`, `depth`, `paths`), `panacus`, `gfatools`, or `minimap2`+`Liftoff` gene-loss checks, extend `horizon_gfa` (`src/bin/horizon_gfa.rs`) only where a single streaming Rust binary adds value:
  - Per-path retained bp and `W`/`P` dialect consistency (avoiding mixed `P`+`W` double counting).
  - Explicitly budgeted (~6–8 GB RAM at `81.4M` HPRC segments) connected-component / inter-chromosomal bridge reporting, clearly distinguished from `horizon_gfa stats`'s `96.9 MB` stateless 5-counter pass.

---

## Utilities Rejected / Deferred After Round 8

- **`horizon gfa-prune` (`P3`, Deferred):** `vg clip` (`-d`, `-D`, `-m`) already operates on GFA without `XG`, and `horizon_clean.rs` lacks haplotype-AF tracking, segment bridging, and multi-pass deletion-shortcut detection.
- **`horizon inversions` (Rejected):** Existing `horizon_resolve.rs` `PackedGraph` does not implement snarl-tree enumeration or sequence storage and OOMs at `81.4M` segments (`~38.9 GB` projected); `INVPG_annot` already fills the workflow.
- **`horizon eval-null` (Rejected as CLI):** Reversed-sequence negative-control calibration belongs inside `horizon_filter` / `R7-CRC-BIO` rather than a separate CLI duplicating `LASTEX` / `mafComparator`.
- **`horizon synteny` (Rejected):** `FlatSeedIndex`'s `u32` `4 Gbp` cap, `~93 GB` single-genome build peak, and pairwise-only LIS cannot compete with `Mumemto` (`139 GB` serial N-way multi-MUMs), `RagTag`, `SyRI`, or `Panagram`.
