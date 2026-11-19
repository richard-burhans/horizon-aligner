# Horizon Pangenome Suite

> **Implementation Status (2026-11-19):** `horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_resolve`, `horizon_wga --graph` (seeded DAG/POA MDL DP in `src/mdl/dag.rs` with bidirected `-` GFA links), `src/gfa/mod.rs`, `src/coord.rs`, and MDL-calibrated `src/egraph.rs` (`MdlCostModel`) are all implemented and unit-tested.

The `horizon` suite extends beyond pairwise whole-genome alignment (WGA) by building on our buffered GFA v1.0/v1.1 streaming reader (`src/gfa/mod.rs`, `MultiGzDecoder` + `stdin`) to provide a composable toolkit grounded in Minimum Description Length (MDL) and E-graphs.

## 1. Pangenome De-tangling and Polishing (`horizon_clean` — Implemented)
**The Problem:** Current graph builders (like Minigraph-Cactus) rely on affine gap penalties, which can over-align repetitive or highly divergent regions, creating spurious edges and tangled "hairballs."
**The Implementation (`src/bin/horizon_clean.rs`):**
* Streams an external `.gfa[.gz]` (or `-` stdin) using a two-pass `S`-flank-indexed reader (bounded by `--flank`, default `128` bp per segment end).
* Scores candidate `L`-line junctions with `bocpd_gc_profile` and an MDL two-regime changepoint bit gain (`--min-cp-bits`).
* Severs unsupported discordant edges into **Horizon Bubbles** (`*`-sequence `HB:i:1` markers with `CP:f:` scores) and splits affected `P`/`W` walks with updated `start`/`end` coordinates to preserve path validity.

## 2. Sequence-to-Graph Alignment (`horizon_wga --graph` — Implemented)
**The Problem:** Aligning novel, divergent genomes to an existing pangenome (e.g., HPRC) is challenging across deep evolutionary splits.
**The Implementation (`src/mdl/dag.rs` + `src/bin/horizon_wga.rs`):**
* `align_seq_to_dag_mdl` generalizes `BandGrid` (`src/mdl/dp.rs`) to topological DAG / partial-order Pair-HMM MDL + BOCPD dynamic programming across graph bifurcations, with Forward-Backward posterior gating and Horizon Bubble abstention.
* `align_seq_to_dag_seeded` uses `MASK_W10` + `MASK_RY20` spaced seeds and forward/backward reachability closure to isolate relevant DAG components before running `align_seq_to_dag_mdl`.
* `horizon_wga --graph <target.gfa[.gz]> <query.fa>` loads GFA `S`/`L` records (expanding `-`-orientation links into bidirected `+`/`-` nodes) via Kahn topological sort and emits the optimal node path, kept base count, MDL bit score, and Horizon Bubble count.

## 3. Global Structural Paralog Resolution (`src/egraph.rs` + `horizon_resolve` — Implemented)
**The Problem:** Complex structural variants (nested inversions, tandem duplications) produce cyclical or tangled graph topologies that complicate downstream variant calling and visualization.
**The Implementation (`src/egraph.rs` + `src/bin/horizon_resolve.rs`):**
* `src/egraph.rs` replaces the original unitless weights (`1.0`/`1.2`/`1.3`) with `MdlCostModel` (`extract_min_mdl_cost`), calibrating leaf catalog pointer bits ($\log_2 K$), `seq` (`1.0` bit), `inv` (`1.0` bit), and `dup` copy-header + sequence residual bits ($H_2(d) + d \log_2 3$).
* Active both inside `horizon_wga --resolve-dups` (benchmarked on the 5 Mbp MHC/HLA locus and 20 Mbp `chr6` q-arm) and as standalone `horizon_resolve` (`src/bin/horizon_resolve.rs`), which streams GFA v1.0/v1.1 `P`/`W` lines in two passes, decomposes walks into bounded snarl/bubble windows (`--max-bubble-steps`), cancels double inversions, and detects/collapses tandem duplications (`--collapse-dups`, `--tag-walks`).

## 4. Graph-Aware Exact Liftover (`horizon_lift` — Implemented)
**The Problem:** Lifting over genomic annotations between haplotypes in a pangenome is prone to coordinate-translation bugs and fragmentation.
**The Implementation (`src/bin/horizon_lift.rs`):**
* Uses `horizon_aligner::coord::Interval` (0-based half-open + UCSC reverse-strand arithmetic) and two-pass GFA `P`/`W` streaming (`wanted_chroms` filtered so peak RSS on the 12 GB HPRC `chr22` graph is `619.2 MB`).
* Supports BED and GFF3 inputs, forward and inverted segments, `--bridge-bp` micro-bubble coalescing, `--min-frac` filtering, and explicit `EXACT` / `PARTIAL` / `SPLIT_k/N` / `DUP` / `--unmapped` reporting.
