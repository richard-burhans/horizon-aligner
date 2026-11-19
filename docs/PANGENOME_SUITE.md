# Horizon Pangenome Suite: Re-imagining Large-Scale Genomic Tooling

The vision for the `horizon` suite extends beyond pairwise whole-genome alignment (WGA). By implementing a native Graphical Fragment Assembly (GFA v1.0/v1.1) streaming engine, `horizon` aims to replace the fragmented ecosystem of current pangenome tools (e.g., `vg`, `odgi`, `hal`, `liftoff`) with a unified, theoretically grounded toolkit based on Minimum Description Length (MDL) and E-graphs.

## 1. Pangenome De-tangling and Polishing (`horizon_clean`)
**The Problem:** Current graph builders (like Minigraph-Cactus) rely on affine gap penalties, which force alignments in repetitive or highly divergent regions, creating spurious edges and tangled "hairballs."
**The Horizon Solution:**
* Stream an external `.gfa` file and evaluate every Link (`L` line) using the Horizon BOCPD/MDL scoring engine.
* If a path through the graph does not yield sufficient bit-savings, explicitly sever the edge and emit a **Horizon Bubble**.
* **Impact:** Automatically relaxes and cleans over-aligned pangenomes, producing sparse, high-confidence biological topologies without relying on arbitrary graph pruning heuristics.

## 2. Sequence-to-Graph Alignment (`horizon_wga --graph`)
**The Problem:** Aligning novel, highly divergent genomes to an existing human pangenome (e.g., HPRC) is computationally prohibitive and fails across deep evolutionary splits.
**The Horizon Solution:**
* Ingest the sequence (`S`) lines of a GFA and build a graph-aware spaced-seed index, utilizing transition-tolerant (`MASK_RY20`) and Amino Acid masks.
* Thread the query sequence through the graph, performing banded MDL dynamic programming across graph bifurcations.
* **Impact:** Enables the direct mapping of distant relatives (e.g., archaic hominids, distant primates) into established pangenomes while gracefully abstaining at low-confidence divergence thresholds.

## 3. Global Structural Paralog Resolution (`horizon_resolve`)
**The Problem:** Complex structural variants (nested inversions, tandem duplications) cause standard tools to generate cyclical graph topologies that break downstream variant calling and visualization.
**The Horizon Solution:**
* Ingest complex GFA subgraphs (hairballs) and map them directly into a congruence-closed **E-graph**.
* Apply e-graph saturation rules to model the evolutionary topology and extract the minimum-cost structural path.
* **Impact:** Automatically disentangles complex loci (like the HLA region or centromeres), resolving paralogs without losing structural variant history.

## 4. Graph-Aware Exact Liftover (`horizon_lift`)
**The Problem:** Lifting over genomic annotations (genes, regulatory elements) between different haplotypes in a pangenome is lossy and prone to off-by-one coordinate errors.
**The Horizon Solution:**
* Traverse the GFA Walk (`W`) lines and Path (`P`) lines using Horizon's strict `coord::Interval` zero-based half-open arithmetic.
* Project functional annotations from a reference path (e.g., `GRCh38`) to any alternate haplotype strictly via high-confidence structural paths.
* **Impact:** Provides a mathematically proven, boundary-safe liftover tool directly on the graph, replacing fragile chain/net workflows.

---
*This document outlines the strategic roadmap for expanding `horizon-aligner` into a comprehensive graph processing suite.*
