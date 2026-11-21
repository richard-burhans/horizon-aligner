# Horizon Aligner: Complete Technical Manual

**Version:** 0.1.0  
**Author:** Richard Burhans  

---

## Table of Contents
1. [Introduction and Philosophy](#1-introduction-and-philosophy)
2. [Availability and Compilation](#2-availability-and-compilation)
3. [Overview of Processing Stages](#3-overview-of-processing-stages)
4. [Command Line Syntax](#4-command-line-syntax)
   - [horizon_wga (Pairwise Alignment)](#horizon_wga)
   - [horizon_merge (Multi-Species Merging)](#horizon_merge)
   - [horizon_guide (Progressive WGA)](#horizon_guide)
5. [Processing Stages in Detail](#5-processing-stages-in-detail)
   - [Sequence Input](#sequence-input)
   - [Indexing Target Seed Words](#indexing-target-seed-words)
   - [Seeding](#seeding)
   - [Gap-free Extension (HSP Gating)](#gap-free-extension-hsp-gating)
   - [HSP Chaining and LCB Generation](#hsp-chaining-and-lcb-generation)
   - [Gapped Extension (MDL+BOCPD)](#gapped-extension-mdlbocpd)
   - [Paralog Resolution via E-Graphs](#paralog-resolution-via-e-graphs)
6. [Output Formats](#6-output-formats)
7. [Advanced Topics](#7-advanced-topics)

---

## 1. Introduction and Philosophy

`horizon-aligner` is a fast, reference-free, whole-genome alignment tool in the lineage of BLASTZ and LASTZ. While traditional aligners rely on static, empirically derived substitution matrices (e.g., HOXD70) and affine gap penalties, these static parameters fundamentally miscalibrate when applied to compositionally extreme genomes (e.g., *Plasmodium falciparum* at ~19% GC) or across wide evolutionary spans where background AT/GC drift overwhelms homologous signal. Under Viterbi traceback, random AT-rich segments accumulate positive scores, leading to spurious alignments.

`horizon-aligner` solves this by replacing the static affine penalty model with a **Minimum Description Length (MDL)** framework, governed by **Bayesian Online Changepoint Detection (BOCPD)**. 

### Horizon Bubbles

![Horizon Bubbles Concept](images/bubble_v4.svg)

*Figure 2. Explicit abstention produces unaligned Horizon Bubbles preserving graph topology.* Orange blocks represent chained locally colinear blocks (LCBs). The solid green path depicts an aligned region passing the minimum bit-savings threshold. Dashed paths (blue for target, pink for query) represent a "Horizon Bubble"—a region where the algorithm explicitly abstains from alignment, preserving the independent topologies of the divergent sequences without forcing spurious matches.

Instead of forcing low-confidence sequences to align or arbitrarily breaking the alignment, the aligner explicitly abstains. Regions lacking sufficient bit-savings are emitted as "Horizon Bubbles"—unaligned graphical paths preserving the structural topology.

---

## 2. Availability and Compilation

`horizon-aligner` is a pure Rust project. It requires `cargo` and `rustc` (Edition 2024).

```bash
cargo build --release
```

Binaries will be placed in `target/release/`:
* `horizon`: Unified subcommand dispatcher (`wga`, `merge`, `guide`, `gfa`, `lift`, `clean`, `resolve`, `eval`, `bench-chr1`).
* `horizon_wga`: Pairwise WGA, DAG/POA (`--graph`), and HAL inspection (`--hal-info`).
* `horizon_merge`: Exact-coordinate multi-species MAF/rGFA merger.
* `horizon_guide`: Progressive guide-tree aligner with HAL export.
* `horizon_gfa`: Streaming GFA v1.0/v1.1 `stats` and chromosome `extract`.
* `horizon_lift`: Graph-aware BED/GFF3 coordinate liftover.
* `horizon_clean`: BOCPD+MDL GFA junction scorer and Horizon Bubble cutter.
* `horizon_resolve`: Standalone GFA `PackedGraph` snarl + e-graph walk resolver.
* `horizon_filter`: Streaming SDUST+BOCPD+MDL PAF filter and `--net` interval netter.
* `horizon_locus`: Regional `impg`-style multi-FASTA locus aligner and similarity matrix emitter.
* `horizon_eval`: Pairwise MAF recall/precision evaluator.

---

## 3. Overview of Processing Stages

![Horizon Aligner Pipeline](images/pipeline_v4.svg)

*Figure 1. Flowchart of the Horizon Aligner processing stages.* The algorithm proceeds through target indexing, spaced-seed querying, ungapped extension, banded gapped alignment, and e-graph paralog resolution. Arrows indicate the flow of sequence data, highlighting the explicit abstention pathway that bypasses forced alignment in low-confidence regions. Colors denote distinct algorithmic phases.

The general flow of `horizon-aligner` mimics the traditional BLASTZ/LASTZ pipeline but introduces novel theoretical frameworks at key bottlenecks:

1. **Target Indexing:** The target sequence is mapped using exact spaced-seed hashing.
2. **Seeding:** The query sequence is hashed and intersected with the target to find seed hits.
3. **Ungapped Extension:** Seeds are extended without gaps into High-Scoring Segment Pairs (HSPs).
4. **Chaining:** HSPs are chained into Locally Colinear Blocks (LCBs).
5. **Gapped Extension:** The gaps between LCBs are evaluated using banded DP guided by MDL bit-savings.
6. **Resolution:** Overlapping structures are saturated in a congruence-closed e-graph.

---

## 4. Command Line Syntax

### `horizon_wga`
Performs pairwise whole-genome alignment.

```bash
horizon_wga <target.fa> <query.fa> [max_contigs] [OPTIONS]
```
**Options:**
* `--maf <path>`: Write the resulting colinear alignment blocks to a MAF v1 file.
* `--gfa <path>`: Write the resulting alignment graph, including Horizon Bubbles, to a GFA v1.1 file.
* `--hal-info <path>`: Decode an HDF5-encoded HAL file directly.
* `--graph <target.gfa[.gz]>`: Align `query.fa` against a GFA pangenome graph using seeded topological DAG/POA MDL dynamic programming (`src/mdl/dag.rs`).
* `--top-k <K>`: Retain only the highest scoring $K$ alignments per query contig.
* `--rbest`: Syntactic sugar for `--top-k 1` (Reciprocal Best Hit).
* `--resolve-dups`: Enable the E-Graph structural paralog resolver.
* `--respect-softmask`: Treat lowercase nucleotides as masked.
* `--ry-seeds` / `--no-ry-seeds`: Toggle the RY (purine/pyrimidine) fallback seeding for deep divergence.
* `--min-shared-seeds <N>`: Set the minimum number of shared seeds (default: 8).
* `--min-top-lcb <N>`: Set the minimum length of the top LCB block.

### `horizon_merge`
Executes an $O(N \log N)$ sweep-line exact-coordinate algorithm to merge multiple pairwise MAFs against a common reference into a single multi-species MAF/rGFA.

```bash
horizon_merge <ref_species> <bin_bp> <maf_out> <gfa_out> <sp1=path1.maf> [<sp2=path2.maf> ...]
```
* `bin_bp`: Window size for the exact-coordinate projection (e.g. 1000).

### `horizon_guide`
Orchestrates a progressive alignment across a guide tree using Nyström low-rank Gram approximations.

```bash
horizon_guide [--progressive] [--respect-softmask] [--hal] [--rbest | --top-k K | --all-pairs] <out_prefix> <max_contigs> <sp1=fa1> <sp2=fa2>...
```

### `horizon_gfa`
Streams plain, BGZF/gzipped (`.gfa.gz`), or `stdin` (`-`) GFA v1.0/v1.1 graphs for `O(1)`-memory structural `stats` or two-pass chromosome subgraph `extract`.

```bash
horizon_gfa stats <in.gfa[.gz]>
horizon_gfa extract <in.gfa[.gz]> --chrom <prefix> [--chrom <prefix2>] --out <out.gfa>
```

### `horizon_lift`
Lifts BED or GFF3 coordinates between haplotypes in a GFA pangenome graph using two-pass `P`/`W` streaming, `--bridge-bp` micro-bubble coalescing, and UCSC reverse-strand arithmetic.

```bash
horizon_lift <in.gfa[.gz]> <intervals.bed|gff3> --from <src_prefix> --to <tgt_prefix> [--bridge-bp 10] [--min-frac 0.5] [--unmapped]
```

### `horizon_clean`
Scores candidate GFA `L`-line junctions with `bocpd_gc_profile` and MDL two-regime changepoint bit gain (`--min-cp-bits`), severing low-support (`walk_support <= 1` by default, or `0` with `--prune-unsupported`) discordant edges into `HB:i:1` Horizon Bubbles while accounting for `L`-line overlap CIGARs when splitting `W` walks.

```bash
horizon_clean --gfa <in.gfa[.gz|.zst]|-> [--flank 128] [--min-cp-bits 12.0] [--prune-unsupported] > clean.gfa
```

### `horizon_resolve`
Runs standalone MDL-calibrated e-graph resolution (`MdlCostModel`) over a `u32`-interned `PackedGraph` with $O(V + E)$ Tarjan biedged snarl decomposition and per-window divergence estimation (optional `--divergence` override), canceling double inversions, collapsing adjacent tandem duplications, updating `W`-line `end` bounds, and pruning/synthesizing `L` links.

```bash
horizon_resolve <in.gfa[.gz|.zst]|-> [--max-bubble-steps 16] [--min-anchor-bp 200] [--divergence 0.02] [--collapse-dups] [--tag-walks] > resolved.gfa
```

### `horizon_filter`
Streams `.paf[.gz|.zst]` or `stdin`, evaluates `FastGA`/`wfmash` PAF alignments with SDUST repeat masking, BOCPD compositional tracking, and Pair-HMM MDL bit-savings (`md:f:` and `hb:i:` tags), and optionally nets overlapping alignments (`--net` / `--max-occ`).

```bash
horizon_filter --paf <in.paf[.gz|.zst]|-> [--query <q.fa>] [--target <t.fa>] [--seqs <all.fa>] [--out <out.paf[.gz|.zst]>] [--min-bits <F>] [--min-nonrepeat-bp <N>] [--min-len <N>] [--net | --max-occ <K>] [--keep-all]
```

### `horizon_locus`
Aligns all sequence pairs in a regional `impg`-style multi-FASTA (`.fa[.gz|.zst]`) using `SDUST + BOCPD/MDL + Aa5Index`, emitting `<prefix>.similarity.tsv` (`jaccard_cov`, `p_dist`, `jc69_dist`, Horizon Bubbles), `<prefix>.maf`, and `<prefix>.gfa`.

```bash
horizon_locus --fasta <locus.fa[.gz|.zst]> --out-prefix <prefix> [--ref <seq_name>] [--threads <N>] [--respect-softmask]
```

### `horizon_eval`
Evaluates pairwise alignment recall and precision against a reference MAF file.

```bash
horizon_eval <ref_name> <tgt_name> <ref.fa> <tgt.fa> <truth.maf> [extra_horizon_wga_flags...]
```

---

## 5. Processing Stages in Detail

### Sequence Input
Sequences are read from plain, `.gz`/`.bgz`, or `.zst`/`.zstd` FASTA files via `src/io_util.rs`. If `--respect-softmask` is supplied, lowercase bases (e.g. `a`, `c`, `g`, `t`) are ignored during indexing, dramatically reducing seed collision depth in repetitive regions.

### Indexing Target Seed Words
`horizon-aligner` utilizes a custom `ContigIndex` struct mapping spaced-seed hashes to chromosomal offsets.

### Seeding
Instead of contiguous k-mers, `horizon-aligner` leverages **spaced seeds**. 
* **Primary Seeding:** Uses `MASK_PRIMARY_11_15` (`111011011011011`). Matches strictly require 11 bases to match across a 15bp window.
* **Fallback Seeding (Distant Divergence):** If the primary seed coverage drops below `1.0` (as tuned dynamically), `horizon-aligner` engages aggressive fallback seeds:
  * **Transition-Tolerant Seeds (`MASK_LASTZ12`, `MASK_LASTZ11`):** Captures transition mutations.
  * **R/Y Seeds (`MASK_RY20`):** Translates the genome into a 1-bit purine (R) / pyrimidine (Y) alphabet before seeding. Highly effective for deep evolutionary splits.
  * **Amino Acid Translated Seeds:** Executes a 6-frame translation. Seeds are matched on a 5-of-6 amino acid shape.

### Gap-free Extension (HSP Gating)
Candidate seeds are extended in both directions without gaps to form High-Scoring Segment Pairs (HSPs). HSPs must exceed dynamic thresholds (e.g., `180 + hsp_bonus` for LASTZ seeds, `220` for RY seeds, `90` for Amino Acid seeds) to pass.

### HSP Chaining and LCB Generation
Valid HSPs are chained using a **Longest Increasing Subsequence (LIS)** algorithm to form Locally Colinear Blocks (LCBs). For gaps between LCBs $\le 2.5$ kb, the aligner lowers its stringency and sweeps the gap with a dense `w=9` (`MASK_W9`) or `w=10` (`MASK_W10`) seed mask.

### Gapped Extension (MDL+BOCPD)
The core DP algorithm (analogous to LASTZ's Y-drop) is performed on the `band_w` unaligned gaps between LCB anchors:
1. **X-Drop:** Flanks of LCBs are extended into the gap using an X-drop algorithm (max `1500bp`).
2. **BOCPD:** The base-by-base match/mismatch probabilities are inferred locally via a Beta-Bernoulli model to detect GC shifts.
3. **MDL Evaluation:** The dynamic programming calculates the Minimum Description Length (MDL) score. 
4. **Horizon Bubbles:** If the DP bit-savings fall below `min_delta_bits`, or if the gap exceeds `max_gap_bp` (2400bp), the aligner explicitly abstains, yielding a macro or internal bubble in the output `GapResult`.

### Paralog Resolution via E-Graphs
To handle nested inversions and tandem duplications, overlapping `MafRow` blocks are modeled as an `ENode` inside a congruence-closed `EGraph`. The minimum-cost structural path is extracted, and paralogs are deduplicated without losing the evolutionary topology.

---

## 6. Output Formats

### Multiple Alignment Format (MAF)
A standard, text-based format for colinear alignments. 
* Coordinates are strictly **0-based, half-open** (managed by `coord::Interval`).
* For alignments on the reverse strand (`-`), the start coordinate is `src_size - end_pos`.
```text
##maf version=1 scoring=horizon_mdl
a score=50.0
s hg38.chr22 100 50 + 50000000 ACTGACTG...
s mm10.chr11 200 50 - 120000000 ACTGACTG...
```

### Graphical Fragment Assembly (GFA v1.1)
Outputs the alignment as a sequence graph. Homologous blocks are `S` (Sequence) lines. Unaligned Horizon Bubbles are represented as parallel, disjoint sequences connecting the same bounding LCBs.
```text
H VN:Z:1.1
S hg38__mm10__blk_0 * LN:i:50 XS:i:100 YS:i:200 ST:Z:-
# E-graph paralog metadata is appended to S lines when `--resolve-dups` is used.
```

---

## 7. Advanced Topics

### Coordinate Arithmetic (`coord.rs`)
Coordinate boundaries can easily drift by off-by-one errors during strand inversion. `horizon-aligner` handles all UCSC-style coordinate translation explicitly inside the `coord::Interval` struct, ensuring that span sizes and endpoints remain invariant when mapped to the complementary sequence.

### HDF5 HAL Parsing & Writing (`hal/`)
While MAF and GFA are the primary emission targets, `horizon` natively decodes and writes Hierarchical Alignment (HAL) graphs stored in HDF5 (`v0` or `v1` reading in `src/hal/hdf5_read.rs`; `libhdf5`-compatible `v0` writing with multi-level B-trees in `src/hal/hdf5_write.rs`). Streaming `Read + Seek` superblock, object-header, and zlib chunk decoders extract sequence and alignment arrays (`DNA_ARRAY`, `TOP_ARRAY`, `BOTTOM_ARRAY`) without loading the entire archive into memory.
