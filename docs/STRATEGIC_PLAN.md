# Horizon Pangenome Suite: Strategic Plan and Risk Mitigation

> **Status Note (2026-11-19):** This document defines the target architecture and roadmap for expanding `horizon-aligner` into a composable pangenome suite. Currently implemented binaries in `src/bin/` are `horizon_wga`, `horizon_merge`, `horizon_guide`, `horizon_eval`, `bench_chr1`, and `horizon_gfa` (streaming GFA `stats`/`list-walks`/`extract`). The remaining pangenome tools (`horizon_clean`, `horizon_resolve`, `horizon_lift`, `horizon_wga --graph`) and `stdin`/`stdout` pipe streaming across all stages are **planned**, not yet implemented.

## 1. The Vision: Composition over Monolith
The history of bioinformatics alignment is a history of composition: `lastz` (pairwise) feeding into `multiz` (MSA), or `cactus` (graph) feeding into `hal` (storage). The vision for `horizon-aligner` is to unify these capabilities under a single framework based on Minimum Description Length (MDL) and E-graphs.

Rather than building a monolithic "run-all" executable, the target architecture is a suite of specialized tools (`horizon_clean`, `horizon_wga`, `horizon_resolve`, `horizon_lift`) that compose cleanly—following the Unix-pipe precedent of `gfatools` and `odgi`—wherever single-pass or indexed streaming is feasible.

## 2. Current Foundation and Known Gaps
* **MDL Bit-Savings (`src/bocpd.rs`, `src/mdl/`):** We replace classical substitution matrices with BOCPD joint-pair bit-savings (`joint_pair_bits`) and use normalized Pair-HMM transitions (`src/mdl/dp.rs`, `M` and `Ix`/`Iy` outgoing probabilities sum to `1.0`). However, the current implementation still retains fixed gap-bit penalties (`GAP_OPEN_BITS = -4.5`, `GAP_EXTEND_BITS = -0.8`), `IDENTITY_PRIOR = 0.65`, and heuristic chaining thresholds (`max_gap_bp`, `hsp_threshold`, `min_block_bits`). Furthermore, MDL abstention suppresses false homology on *divergent* sequences; it does not by itself prevent high-identity paralogous or repetitive "hairballs."
* **E-Graphs for Structural Variants (`src/egraph.rs`):** Congruence-closure rewriting currently operates on 1D linear signed-block expressions (`seq`/`inv`/`dup` in 12-block chunks, `max_nodes = 800`) using fixed unitless parsimony weights (`1.0`/`1.2`/`1.3`). It does not yet model evolutionary trajectories, eliminate GFA cycles, or have HLA/centromere benchmarks.
* **Buffered GFA Streaming Parser (`src/gfa/mod.rs`, `src/bin/horizon_gfa.rs`):** `stream_gfa` provides a buffered line-streaming GFA v1.0/v1.1 reader (supporting `.gfa`, multi-member BGZF `.gfa.gz` via `MultiGzDecoder`, and `"-"` `stdin`). It is **not** zero-copy (it allocates owned `String`/`Vec<u8>` fields per record), but on the 12 GB HPRC v1.0 Minigraph-Cactus graph it streams `81.4M` segments (`3.29 Gbp`) and `113.0M` links in `88.84s` at **97.1 MB peak RSS**, and extracts `GRCh38.chr22` + `CHM13.chr22` in `166.73s` at **456.1 MB peak RSS**. While modern indexed/disk-backed representations (`GBZ` / `GBZ-base` by Sirén et al., `gfaidx`) also address subgraph queries without full-graph RAM materialization, native streaming gives Horizon direct GFA v1.0/v1.1 interoperability without external `vg`/`odgi` dependencies.

## 3. Strategic Adjustments & Risk Mitigations

### A. Avoid a Proprietary "Walled Garden" Format
**Risk:** Binary-only formats with heavy external C dependencies (such as HDF5 `.hal` in `Cactus`) complicate downstream parsing and ecosystem interoperability. (For backwards-compatibility, `horizon-aligner` currently includes a pure-Rust `.hal` reader/writer in `src/hal/` via `horizon_wga --hal-info`, which will remain optional.)
**Mitigation:** We will **not invent a proprietary binary format** for the pangenome suite. Target tools will exchange standard GFA (v1.0/v1.1) and MAF (plus VCF where variant output is added), using only valid GFA optional tags/comments (e.g., `HB:i:`, `XS:i:`, `YS:i:`, `ST:Z:`, `# EG:Z:`). Any required random-access indexing (for `S`-line sequences or `W`/`P` coordinates) will use lightweight, optional sidecar files (analogous to `.fai`/`.gzi`), not a locked database.

### B. Emulate Unix Composition Where Topology Allows
**Risk:** Building a monolithic CLI that couples all stages internally.
**Mitigation:** Support `stdin`/`stdout` GFA/MAF streaming on linear filter/transform tools (`horizon_clean`, `horizon_resolve`, `horizon_lift`), while acknowledging that multi-input or random-access stages (`horizon_wga --graph --query`, or unsorted GFA `S`/`L` resolution) require file/sidecar arguments or a two-pass index rather than pure single-pass `stdin` piping.

### C. Localize the E-Graphs (Prevent Combinatorial Explosion)
**Risk:** E-graph saturation under associative `seq` rules grows combinatorially (hitting `max_nodes = 800` at small block counts) and cannot scale to whole chromosomes.
**Mitigation:** Strict **bubble extraction** via a topological decomposition pass (analogous to `vg` snarls / `smoothxg`) rather than single-pass unsorted GFA streaming. Once tangled subgraphs ("hairballs") are isolated into bounded bubbles, the E-graph resolves each bubble locally and writes the resolved path back into the graph.

### D. Target Graph Tooling Gaps First
**Risk:** Competing head-on against `lastz` or `minimap2` on routine pairwise alignment offers limited differentiation.
**Mitigation:** Focus our next implementation and benchmarking milestones on pangenome graph tooling (`horizon_clean`, `horizon_resolve`, `horizon_lift`, and HPRC `GRCh38` vs. `CHM13` evaluation, which is currently unrun), where MDL abstention and bounded structural resolution address concrete pain points in complex loci.
