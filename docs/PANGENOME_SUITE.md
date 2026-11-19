# Horizon Pangenome Suite: Proposed Roadmap

> **Status Note (2026-11-19):** None of the four pangenome binaries/flags below (`horizon_clean`, `horizon_wga --graph`, `horizon_resolve`, `horizon_lift`) are implemented yet in `src/bin/`. `src/gfa/mod.rs` (`stream_gfa`) is currently an uncalled, untested GFA v1.0/v1.1 reader, and `src/coord.rs` is not yet compiled into `src/lib.rs`. This document outlines the proposed engineering roadmap.

The target vision for the `horizon` suite extends beyond pairwise whole-genome alignment (WGA) by building on our buffered GFA v1.0/v1.1 streaming reader (`src/gfa/mod.rs`) to provide a unified toolkit grounded in Minimum Description Length (MDL) and E-graphs.

## 1. Pangenome De-tangling and Polishing (`horizon_clean` — Proposed)
**The Problem:** Current graph builders (like Minigraph-Cactus) rely on affine gap penalties, which can over-align repetitive or highly divergent regions, creating spurious edges and tangled "hairballs."
**The Proposed Solution:**
* Stream an external `.gfa` file using a two-pass or `S`-line-indexed reader (since `L` lines carry no DNA sequences and `S`/`L` lines can appear in arbitrary order).
* Evaluate candidate junction flanks across `L` lines using the Horizon BOCPD/MDL scoring engine.
* If a junction does not yield sufficient bit-savings, explicitly sever the edge, emit a **Horizon Bubble** (`*`-sequence `HB:i:1` marker), and update affected `W`/`P` walks to preserve path validity.
* **Target Impact:** Relaxes over-aligned pangenomes without arbitrary graph pruning heuristics.

## 2. Sequence-to-Graph Alignment (`horizon_wga --graph` — Proposed)
**The Problem:** Aligning novel, divergent genomes to an existing pangenome (e.g., HPRC) is challenging across deep evolutionary splits.
**The Proposed Solution:**
* Ingest the sequence (`S`) lines of a GFA into a graph-aware spaced-seed index, utilizing transition-tolerant (`MASK_RY20`, 25-bp) and translated Amino Acid (`AA5`, requiring >=18-bp / 6-codon spans) masks across node boundaries.
* Extend `BandGrid` (`src/mdl/dp.rs`, currently strictly 2D linear-diagonal) to banded partial-order / DAG MDL dynamic programming across graph bifurcations (building on prior art such as `GraphAligner` and `minigraph`).
* **Target Impact:** Enables mapping distant relatives into established pangenomes while abstaining at low-confidence divergence thresholds (noting that on our `cactus241way` 1MB nuclear test, abstention reduces recall from 45.3% to 9.9% with 494 Horizon Bubbles).

## 3. Global Structural Paralog Resolution (`horizon_resolve` — Proposed)
**The Problem:** Complex structural variants (nested inversions, tandem duplications) produce cyclical or tangled graph topologies that complicate downstream variant calling and visualization.
**The Proposed Solution:**
* Decompose complex GFA subgraphs ("hairballs") into localized bubbles and map them into a congruence-closed **E-graph** (`src/egraph.rs`).
* Replace the current fixed unitless weights (`1.0`/`1.2`/`1.3`) with an MDL-calibrated structural cost model and extract minimum-cost paths.
* **Target Impact:** Disentangles complex loci (such as HLA or centromeric regions, once benchmarked) while preserving structural variant history.

## 4. Graph-Aware Exact Liftover (`horizon_lift` — Proposed)
**The Problem:** Lifting over genomic annotations between haplotypes in a pangenome is prone to coordinate-translation bugs and fragmentation.
**The Proposed Solution:**
* Wire `src/coord.rs` (`coord::Interval`, zero-based half-open arithmetic) into `src/lib.rs`, build an `S`-line segment-length index, and traverse GFA Walk (`W`) and Path (`P`) lines.
* Project functional annotations from a reference path (e.g., `GRCh38`) to alternate haplotypes via high-confidence structural paths, explicitly flagging annotations dropped or fragmented across non-syntenic breaks.
* **Target Impact:** Provides a boundary-safe liftover tool directly on GFA walks/paths.
