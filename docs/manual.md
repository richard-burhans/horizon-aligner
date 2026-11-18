# Horizon Aligner: Comprehensive Manual

**Version:** 0.1.0  
**Author:** Richard Burhans  

## 1. Introduction

`horizon-aligner` is a highly scalable, reference-free, whole-genome alignment tool designed to address the unique challenges of aligning evolutionarily distant sequences and compositionally extreme genomes (e.g., *Plasmodium falciparum*). 

Unlike traditional aligners (e.g., LASTZ) that use a static affine gap scoring matrix, `horizon-aligner` utilizes **Minimum Description Length (MDL)** paired with **Bayesian Online Changepoint Detection (BOCPD)** to dynamically adapt to local GC-content shifts and reliably distinguish homologous sequence from random compositional similarities. Where the aligner lacks statistical confidence to bridge a gap, it explicitly abstains, emitting a **Horizon Bubble**.

## 2. Algorithm Overview

The alignment process proceeds through five distinct stages, from initial heuristic seeding to exact coordinate merging and e-graph-based structural deduplication.

### Stage 1: Multi-Shape Seeding
Instead of canonical k-mers, `horizon-aligner` relies on **spaced seeds**. 
* **Primary Seeding:** Uses a weight-11, span-15 mask (`MASK_PRIMARY_11_15`: `111011011011011`). Matches are found using exact hashing. 
* **Fallback Seeding:** For divergent clades (where primary coverage drops below thresholds), the aligner engages:
  * **HSP-gated Multi-shape Seeds:** Complementary 11-of-18 and 12-of-19 shapes (`MASK_LASTZ11`, `MASK_LASTZ12`) to capture transition mutations.
  * **Amino-Acid (AA) Seeds:** 6-frame translated 5-of-6 amino acid seeding for highly divergent coding homologs.
  * **RY (Purine/Pyrimidine) Seeds:** Transition-tolerant 1-bit alphabet seeding (`MASK_RY20`) to capture structurally conserved but nucleotide-divergent regions.

### Stage 2: HSP Gating & LCB Chaining
Matches from the seeding stage form High-Scoring Segment Pairs (HSPs). 
* The aligner computes the **Longest Increasing Subsequence (LIS)** of these HSPs to form Locally Colinear Blocks (LCBs).
* Short overlapping or weak off-diagonal pairs are pruned. A fallback to `w=10` or `w=9` seeds is engaged inside gaps of `<= 2.5 kb` to bridge sparse exons.

### Stage 3: Banded MDL + BOCPD Gapped Alignment
This is the core evaluation stage. Large unaligned gaps between LCBs (or gaps at the ends of contigs) are passed to the dynamic programming (DP) module.
* **X-Drop Extension:** Flanks are extended using an X-drop algorithm to naturally define the end of homologous signal.
* **Banded DP:** A banded DP gap alignment evaluates the edit path using **MDL bit-savings**. The underlying null model dynamically adjusts to local GC content shifts using a Beta-Bernoulli BOCPD model.
* **Horizon Bubbles:** If the DP bit-savings fall below a predefined conformal risk threshold (`min_delta_bits`), the aligner explicitly abstains. Large gaps become *macro bubbles*, and DP dropouts become *internal bubbles*.

### Stage 4: Paralog Resolution (E-Graph)
Structural variation (tandem duplications, nested inversions) results in overlapping colinear blocks. `horizon-aligner` models these using a **congruence-closed E-Graph** to extract the minimum-cost synteny expression, natively collapsing duplicated paralogs into clean path expressions rather than dropping them.

### Stage 5: Formatting and Merging
The abstract `ContigAlignment` state is losslessly projected into standard output formats:
* **MAF v1:** A colinear Multiple Alignment Format.
* **rGFA v1.1:** A graphical format capturing the exact structure of Horizon Bubbles.
* **Multi-Species Merge:** `horizon_merge` performs an $O(N \log N)$ exact-coordinate sweep-line merge of pairwise MAFs into a unified multi-species alignment, threading gap insertions cleanly.

---

## 3. Command-Line Options

### 3.1 `horizon_wga`
The primary pairwise whole-genome aligner.

**Usage:**
```bash
horizon_wga <target.fa> <query.fa> [max_contigs] [OPTIONS]
```

#### Output Options
* `--maf <file.maf>`: Output the alignment in MAF format.
* `--gfa <file.gfa>`: Output the alignment in GFA v1.1 format (includes Horizon Bubbles).
* `--hal-info <file.hal>`: Parse and decode HDF5 Hierarchical Alignment Format (HAL) metadata.

#### Filtering and Scoring Options
* `--top-k <K>`: Keep only the top `K` scoring contig-pair alignments. Prunes off-diagonal weak homologies.
* `--rbest`: Syntactic sugar for `--top-k 1`. Computes Reciprocal Best Hits.
* `--min-shared-seeds <N>`: Minimum number of shared seeds required to evaluate a contig pair (default: `8`).
* `--min-top-lcb <N>`: Minimum size of the top Locally Colinear Block (LCB) to retain the pair.

#### Seeding Strategies
* `--respect-softmask`: Skip lowercase softmasked bases during indexing/seeding. Prevents repeats from dominating the seed tables.
* `--ry-seeds`: Opt-in to transition-tolerant purine/pyrimidine (RY) fallback seeds for highly divergent pairs.
* `--no-ry-seeds`: Explicitly disable RY seeding.

#### Structural Resolution
* `--resolve-dups`: Enable the E-Graph duplication/inversion resolver. Automatically tags paralogs (`# EG:Z:`) in the GFA output.

---

### 3.2 `horizon_merge`
Merges multiple pairwise MAFs into a single reference-guided multi-species MAF/GFA.

**Usage:**
```bash
horizon_merge <ref_species> <bin_bp> <maf_out> <gfa_out> <sp1=path1.maf> [<sp2=path2.maf> ...]
```
* `<ref_species>`: The name of the reference species sequence.
* `<bin_bp>`: Bin size for the exact-coordinate projection.
* `<spX=pathX.maf>`: Pairs of species identifiers and their pairwise MAF paths.

---

### 3.3 `horizon_guide`
A progressive multi-species WGA pipeline that utilizes Nyström/sketch-selected hinges.

**Usage:**
```bash
horizon_guide [OPTIONS] <out_prefix> <max_contigs> <sp1=fa1> <sp2=fa2>...
```
#### Options
* `--progressive`: Enables identity-gated transitive sister-clade rescue across the guide tree.
* `--hal`: Enables hierarchical alignment output generation.
* `--all-pairs`: Generate all pairwise alignments (`--top-k 0`).

---

## 4. Output Formats

### 4.1 Multiple Alignment Format (MAF)
`horizon-aligner` emits standard MAF v1. A block begins with an `a` line denoting the score, followed by `s` lines for the sequences. Coordinates are strictly 0-based, half-open, with reverse strand (`-`) coordinates given relative to the full length of the forward strand.
```text
##maf version=1 scoring=horizon_mdl
a score=432.0
s hg38.chr22 1000 432 + 50818468 ACGT...
s mm10.chr11 5000 432 - 122082543 ACGT...
```

### 4.2 Graphical Fragment Assembly (GFA v1.1)
GFA output captures homologous blocks as Sequence (`S`) lines and represents structural abstentions (Horizon Bubbles) as divergence in the graph topology, linked by `L` lines. Paralogs resolved by the E-graph are tagged using `EG:Z:` tags on sequence nodes.

```text
H VN:Z:1.1
S hg38__mm10__blk_0 * LN:i:432 XS:i:1000 YS:i:5000 ST:Z:-
```