# Horizon Pangenome Suite

> **Implementation Status (2026-11-20, `9dfec4a`):** `horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_resolve` (`u32` `PackedGraph` + $O(V+E)$ Tarjan biedged snarl decomposition + per-window divergence), `horizon_filter` (streaming PAF SDUST+BOCPD+MDL filter & netter), `horizon_locus` (regional `impg`-style multi-FASTA locus aligner), `horizon_wga --graph` (seeded DAG/POA MDL DP in `src/mdl/dag.rs` with `MASK_W10`, `MASK_RY20`, and `aa5_keys` 6-frame translated fallback), `src/io_util.rs` (plain, `.gz`/`.bgz`, `.zst`/`.zstd`, `stdin`/`stdout`), `src/paf.rs`, `src/gfa/mod.rs` (`SpoolableGfaInput`), `src/coord.rs`, `src/seeding/flat_index.rs`, and MDL-calibrated `src/egraph.rs` (`MdlCostModel`) are implemented (**118** passing unit tests: `86` lib + `32` bin, `1` ignored; `11` suite binaries + `bench_chr1` in `src/bin/`). Note that `stream_gfa` parses `S`/`L`/`P`/`W` records while skipping `C`/`J` records and defaulting empty `L` orientations to `+`.

The `horizon` suite extends beyond pairwise whole-genome alignment (WGA) by building on our buffered GFA v1.0/v1.1 streaming reader (`src/gfa/mod.rs`, `MultiGzDecoder` + single-pass `stdin`) to provide a composable toolkit grounded in Minimum Description Length (MDL) and E-graphs.

## 1. Pangenome De-tangling and Polishing (`horizon_clean` — Implemented)
**The Problem:** Current graph builders (like Minigraph-Cactus) rely on affine gap penalties, which can over-align repetitive or highly divergent regions, creating spurious edges and tangled "hairballs."
**The Implementation (`src/bin/horizon_clean.rs`):**
* Streams an external `.gfa[.gz|.zst]` using a two-pass `S`-flank-indexed reader (bounded by `--flank`, default `128` bp per segment end; strips `S`-line optional tags on output).
* Scores candidate `L`-line junctions with `bocpd_gc_profile` and an MDL two-regime changepoint bit gain (`--min-cp-bits`).
* Severs low-support discordant edges (`walk_support <= 1` by default, or `0` with `--prune-unsupported`) into **Horizon Bubbles** (`*`-sequence `HB:i:1` markers with `CP:f:` scores) and accounts for `L`-line overlap CIGARs when computing split `W`-line `start`/`end` coordinates.

## 2. Sequence-to-Graph Alignment (`horizon_wga --graph` — Implemented)
**The Problem:** Aligning novel, divergent genomes to an existing pangenome (e.g., HPRC) is challenging across deep evolutionary splits.
**The Implementation (`src/mdl/dag.rs` + `src/bin/horizon_wga.rs`):**
* `align_seq_to_dag_mdl` generalizes `BandGrid` (`src/mdl/dp.rs`) to topological DAG / partial-order Pair-HMM MDL + BOCPD dynamic programming across graph bifurcations (omitting `-log2(out_degree)` branching penalties), with Forward-Backward posterior gating and Horizon Bubble abstention.
* `align_seq_to_dag_seeded` uses `MASK_W10` + `MASK_RY20` intra-node spaced seeds (requiring `>=14` bp segments) plus compact `aa5_keys` 6-frame translated seeding fallback (`src/mdl/dag_seed.rs`), and ancestor-then-descendant reachability closure to isolate relevant DAG components before running `align_seq_to_dag_mdl` (bypassed when `nodes.len() <= 8` or `seed_hits.is_empty()`).
* `horizon_wga --graph <target.gfa[.gz|.zst]> <query.fa[.gz|.zst]>` loads the GFA into a topological DAG (expanding `-`-orientation links into bidirected `+`/`-` nodes) via Kahn topological sort and passes `&HmmParams::default()`.

## 3. Global Structural Paralog Resolution (`src/egraph.rs` + `horizon_resolve` — Implemented)
**The Problem:** Complex structural variants (nested inversions, tandem duplications) produce cyclical or tangled graph topologies that complicate downstream variant calling and visualization.
**The Implementation (`src/egraph.rs` + `src/bin/horizon_resolve.rs`):**
* `src/egraph.rs` replaces the original unitless weights (`1.0`/`1.2`/`1.3`) with `MdlCostModel` (`extract_min_mdl_cost`), calibrating leaf catalog pointer bits ($\log_2 K$), `seq` (`1.0` bit), `inv` (`1.0` bit), and `dup` copy-header + sequence residual bits ($H_2(d) + d \log_2 3$).
* Active both inside `horizon_wga --resolve-dups` (emitting `# EG:Z:` tags) and as standalone `horizon_resolve` (`src/bin/horizon_resolve.rs`), which streams GFA v1.0/v1.1 `P`/`W` lines in two passes into a `u32`-interned `PackedGraph`, decomposes graphs via $O(V + E)$ Tarjan biedged snarl decomposition, estimates per-window divergence (with optional `--divergence` override), cancels double inversions, detects/collapses adjacent tandem duplications (`--collapse-dups`, `--tag-walks`), updates `W`-line `end` bounds, and prunes/synthesizes `L` links.

## 5. Streaming PAF Filtering & Regional Locus Alignment (`horizon_filter`, `horizon_locus` — Implemented)
* **`horizon_filter` (`src/bin/horizon_filter.rs`, `src/paf.rs`):** Streams `.paf[.gz|.zst]` or `stdin`, evaluates alignments with SDUST repeat masking (`--min-nonrepeat-bp`), BOCPD compositional profiles, and Pair-HMM MDL bit-savings (`md:f:` / `hb:i:`), and optionally nets overlapping alignments (`--net` / `--max-occ`).
* **`horizon_locus` (`src/bin/horizon_locus.rs`):** Aligns all pairs in a regional `impg`-style multi-FASTA (`.fa[.gz|.zst]`) using `SDUST + BOCPD/MDL + Aa5Index`, emitting `<out>.similarity.tsv` (`jaccard_cov`, `p_dist`, `jc69_dist`, Horizon Bubbles), `<out>.maf`, and `<out>.gfa`.

## 4. Graph-Aware Exact Liftover (`horizon_lift` — Implemented)
**The Problem:** Lifting over genomic annotations between haplotypes in a pangenome is prone to coordinate-translation bugs and fragmentation.
**The Implementation (`src/bin/horizon_lift.rs`):**
* Uses `horizon_aligner::coord::Interval` (0-based half-open + UCSC reverse-strand arithmetic) and two-pass GFA `P`/`W` streaming (`wanted_chroms` filtered so peak RSS on the 12 GB HPRC `chr22` graph is `619.2 MB`).
* Supports BED and GFF3 inputs, forward and inverted segments, `--bridge-bp` micro-bubble coalescing, `--min-frac` filtering, and explicit `EXACT` / `PARTIAL` / `SPLIT_k/N` / `DUP` / `--unmapped` reporting.
