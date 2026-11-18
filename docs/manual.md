# Horizon Aligner: Complete Technical Manual

**Version:** 0.1.0  
**Author:** Richard Burhans  

---

## Table of Contents
1. [Introduction and Philosophy](#1-introduction-and-philosophy)
2. [The Alignment Pipeline](#2-the-alignment-pipeline)
   - [Phase 1: Multi-Shape Spaced Seeding](#phase-1-multi-shape-spaced-seeding)
   - [Phase 2: Ungapped Extension & HSP Gating](#phase-2-ungapped-extension--hsp-gating)
   - [Phase 3: Chaining & Sparse Exon Fallback](#phase-3-chaining--sparse-exon-fallback)
   - [Phase 4: MDL + BOCPD Gapped Alignment](#phase-4-mdl--bocpd-gapped-alignment)
   - [Phase 5: Paralog Resolution via E-Graphs](#phase-5-paralog-resolution-via-e-graphs)
3. [Command Line Reference](#3-command-line-reference)
   - [horizon_wga (Pairwise Alignment)](#horizon_wga)
   - [horizon_merge (Multi-Species Merging)](#horizon_merge)
   - [horizon_guide (Progressive WGA)](#horizon_guide)
4. [Output Formats](#4-output-formats)
   - [Multiple Alignment Format (MAF)](#multiple-alignment-format-maf)
   - [Graphical Fragment Assembly (GFA v1.1)](#graphical-fragment-assembly-gfa-v11)
5. [Performance & Tuning](#5-performance--tuning)

---

## 1. Introduction and Philosophy

Standard whole-genome aligners (e.g., LASTZ, Minimap2) rely on static, empirically derived substitution matrices (like HOXD70) and affine gap penalties. While highly optimized for vertebrate genomes with ~47.5% GC content, these static parameters fundamentally miscalibrate when applied to compositionally extreme genomes (e.g., *Plasmodium falciparum* at ~19% GC) or across wide evolutionary spans where background AT/GC drift overwhelms homologous signal. Under Viterbi traceback, random AT-rich segments accumulate positive scores, leading to spurious alignments.

`horizon-aligner` solves this by replacing the static affine penalty model with a **Minimum Description Length (MDL)** framework, governed by **Bayesian Online Changepoint Detection (BOCPD)**. 

### Core Concepts
* **BOCPD:** The genome is treated as a sequence of regimes with distinct background nucleotide frequencies (e.g., GC-rich exons vs AT-rich intergenic regions). The aligner dynamically infers these regimes using a Beta-Bernoulli conjugate model.
* **MDL Bit-Savings:** An alignment is only accepted if encoding the target sequence *given* the reference sequence requires strictly fewer bits than encoding the target sequence against the independent BOCPD null model.
* **Horizon Bubbles:** Instead of forcing low-confidence sequences to align or arbitrarily breaking the alignment, the aligner explicitly abstains. Regions lacking sufficient bit-savings are emitted as "Horizon Bubbles"—unaligned graphical paths preserving the structural topology.

---

## 2. The Alignment Pipeline

`horizon-aligner` executes a multi-tiered pipeline that filters massive search spaces down to exact coordinate-level base pairings.

### Phase 1: Multi-Shape Spaced Seeding
The aligner utilizes exact hashing of spaced seeds to identify candidate homologies, avoiding the sensitivity drops associated with contiguous k-mers.

* **Primary Seeding:** Uses `MASK_PRIMARY_11_15` (`111011011011011`). Matches are strictly required to share this sequence of 11 matching bases spanning a 15bp window.
* **Fallback Seeding (Distant Divergence):** If the primary seed coverage drops below 35%, `horizon-aligner` assumes a highly divergent homology and engages aggressive fallback seeds:
  * **Transition-Tolerant Seeds (`MASK_LASTZ12`, `MASK_LASTZ11`):** Uses the classic 12-of-19 and 11-of-18 shapes, tuned to allow 1 or 2 transition mutations (A<->G, C<->T) within the seed window.
  * **R/Y Seeds (`MASK_RY20`):** Translates the genome into a 1-bit purine (R) / pyrimidine (Y) alphabet before seeding. Highly effective for deep evolutionary splits where transitions have saturated but the structural purine/pyrimidine backbone is conserved.
  * **Amino Acid Translated Seeds:** Executes a 6-frame translation. Seeds are matched on a 5-of-6 amino acid shape, enabling recovery of synonymous divergent coding genes.

### Phase 2: Ungapped Extension & HSP Gating
Candidate seeds are extended in both directions without gaps to form High-Scoring Segment Pairs (HSPs).
* **Gating:** HSPs must exceed predefined dynamic thresholds (e.g., `180 + hsp_bonus` for LASTZ seeds, `220` for RY seeds, `90` for Amino Acid seeds). The `hsp_bonus` scales inversely with the size of the search space to maintain a constant false-discovery rate across different contig lengths.

### Phase 3: Chaining & Sparse Exon Fallback
Valid HSPs are chained using a **Longest Increasing Subsequence (LIS)** algorithm to form Locally Colinear Blocks (LCBs).
* **Inter-LCB Rescue:** For gaps between LCBs that are $\le 2.5$ kb (often indicating missing sparse exons), the aligner lowers its stringency and sweeps the gap with a dense `w=9` (`MASK_W9`) or `w=10` (`MASK_W10`) seed mask, capturing micro-exons that failed the initial HSP gate.

### Phase 4: MDL + BOCPD Gapped Alignment
Gaps between chained anchors are evaluated at base-pair resolution.
1. **X-Drop Extension:** The flanks of LCBs are extended into the gap using an X-drop algorithm, continuing until the alignment score drops by `max_ext` from its maximum.
2. **Banded DP:** For the remaining unaligned gap, a banded dynamic programming (DP) matrix is allocated (`band_w` half-width).
3. **Evaluation:** The DP computes the MDL score. The base-by-base match/mismatch probabilities are derived from the local BOCPD-inferred GC frequencies. 
4. **Abstention:** If the total bit-savings of the DP path does not exceed `min_delta_bits`, the gap is designated an *internal bubble*. If a gap exceeds `max_gap_bp` (default 2400bp), it bypasses DP and becomes a *macro bubble*.

### Phase 5: Paralog Resolution via E-Graphs
Genomes are replete with tandem duplications and nested inversions. Instead of arbitrarily picking one copy (which breaks downstream synteny) or crashing, `horizon-aligner` models overlapping colinear blocks as equivalence classes in a **congruence-closed E-Graph**.
* The e-graph applies structural rewrite rules (e.g., tandem duplication collapse, inversion flattening) until saturation.
* The minimum-cost expression is extracted. Duplicated blocks are tagged with `# EG:Z:` metadata in the output, completely preserving the paralogous history for downstream pangenome graphs.

---

## 3. Command Line Reference

### `horizon_wga`
Performs pairwise whole-genome alignment between a target and a query.

```bash
horizon_wga <target.fa> <query.fa> [max_contigs] [OPTIONS]
```
**Options:**
* `--maf <path>`: Write the resulting colinear alignment blocks to a MAF v1 file.
* `--gfa <path>`: Write the resulting alignment graph, including Horizon Bubbles and E-Graph nodes, to a GFA v1.1 file.
* `--hal-info <path>`: Read an HDF5-encoded HAL file, extracting the tree structure and underlying FASTA sequences directly.
* `--top-k <K>`: Retain only the highest scoring $K$ alignments per query contig. Use `--rbest` for exactly 1 reciprocal best hit.
* `--resolve-dups`: Enable the Phase 5 E-Graph structural paralog resolver.
* `--respect-softmask`: Treat lowercase nucleotides as masked; they will be skipped during indexing, vastly speeding up alignment of repeat-heavy genomes.
* `--ry-seeds` / `--no-ry-seeds`: Explicitly enable or disable the RY (purine/pyrimidine) fallback seeding for deep divergence.

### `horizon_merge`
Executes an $O(N \log N)$ sweep-line exact-coordinate algorithm to merge multiple pairwise MAFs against a common reference into a single multi-species MAF/rGFA.

```bash
horizon_merge <ref_species> <bin_bp> <maf_out> <gfa_out> <sp1=path1.maf> [<sp2=path2.maf> ...]
```
* `bin_bp`: Controls the spatial partitioning of the sweep-line (typically `1000` or `10000`).

### `horizon_guide`
Orchestrates a progressive alignment across a guide tree. It uses Nyström low-rank Gram approximations and Determinantal Point Processes (DPP) to automatically select optimal "hinge" genomes, aligns to them, and projects the remaining genomes transitively.

```bash
horizon_guide [--progressive] [--respect-softmask] [--hal] [--rbest | --top-k K | --all-pairs] <out_prefix> <max_contigs> <sp1=fa1> <sp2=fa2>...
```

---

## 4. Output Formats

### Multiple Alignment Format (MAF)
A standard, text-based format for colinear alignments. 
* Coordinates are strictly **0-based, half-open** `[start, end)`.
* For alignments on the reverse strand (`-`), the start coordinate is mapped relative to the *forward* strand (i.e., `src_size - end_pos`).
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

## 5. Performance & Tuning
`horizon-aligner` processes pairwise alignments extremely efficiently. 
* **Parallelization:** Gap dynamic programming (the computationally heaviest phase) is parallelized across available CPU cores using Rayon.
* **Complexity Limits:** Banded DP guarantees bounded memory consumption. 
* **Tuning for Distant Divergence:** For pairs separated by >50 million years (e.g., Human vs. Horse), ensure `--ry-seeds` is enabled or rely on the dynamic `primary_cov < 1.0` fallback triggers, as primary KP-15 seeds will silently drop true homologies due to accumulated transition mutations.
