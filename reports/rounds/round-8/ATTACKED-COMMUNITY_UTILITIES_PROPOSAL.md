# Community Utilities Proposal: Large-Scale Alignment & Pangenome Tooling

*Grounded in the project's alignment and pangenome literature corpus (`references.bib`, `genome-alignment-synteny-duplications`, and `bioinformatics-tooling-benchmarking`) and mapped onto existing primitives in `horizon-aligner`.*

---

## Executive Summary

Across recent large-scale whole-genome alignment (WGA) and pangenome graph papers—including HPRC (`liao_2023_draft`), Minigraph-Cactus (`hickey_2023_pangenome`), PGGB (`garrison_2024_building`), ODGI (`guarracino_2022_odgi`), comparative pangenome benchmarks (`andreace_2023_comparing`, `kopalli_2025_toward`), inversion & snarl theory (`romain_2026_inversions-pangenome-graphs`, `sena_2026_bubble-subgraphs-spqr`), and maximal-match diagnostics (`shivakumar_2025_mumemto`)—researchers repeatedly report the same operational gaps:

1. **Dependency & memory walls** (`libhdf5` / C++ `halTools`, `700–800 GB` suffix-array RAM, `70–154 GB` graph builders/indexers).
2. **Ad-hoc multi-tool glue** (chaining `odgi paths` + `minimap2` + `samtools` + `Liftoff`, or `vg deconstruct` + `minimap2` + `bedtools`, just to audit clipping or annotate inversions).
3. **Lack of ground-truth-free specificity & retention diagnostics** (silent sequence clipping in Minigraph-Cactus, giant-component collapses, and aligner benchmarks that rely solely on self-referential coverage or cross-aligner agreement).

Because `horizon-aligner` already provides a zero-`libhdf5` streaming HAL/HDF5 reader (`src/hal/`), a streaming `12 GB+` GFA parser (`src/gfa/mod.rs`), a `u32` `PackedGraph` + `O(V+E)` Tarjan snarl decomposer (`src/bin/horizon_resolve.rs`), strand-symmetric `FlatSeedIndex` (`src/seeding/flat_index.rs`), and reversed-sequence Conformal Risk Control calibration (`src/mdl/dp.rs`), six high-impact community utilities can be built directly on existing internal primitives.

---

## Prioritized Roadmap

| Priority | Proposed Subcommand | Primary Pain Point Addressed | Key Literature Motivation | Existing `horizon-aligner` Primitives |
|---|---|---|---|---|
| **P1** | `horizon hal` | Zero-`libhdf5` HAL inspection, MAF export, and coordinate liftover | `hickey_2013_hierarchical`, `liao_2023_draft`, `armstrong_2020_progressive` | `src/hal/hdf5_read.rs`, `src/hal/mod.rs`, `src/coord.rs` |
| **P1** | `horizon gfa-qc` | Single-pass graph retention, clipping, component & node-depth audit | `kopalli_2025_toward`, `liao_2023_draft`, `andreace_2023_comparing`, `guarracino_2022_odgi` | `src/gfa/mod.rs`, `src/bin/horizon_gfa.rs`, `src/io_util.rs` |
| **P2** | `horizon inversions` | Native GFA path-explicit & alignment-rescued inversion / snarl caller | `romain_2026_inversions-pangenome-graphs`, `sena_2026_bubble-subgraphs-spqr` | `PackedGraph` + Tarjan snarls (`horizon_resolve.rs`), `FlatSeedIndex` |
| **P2** | `horizon gfa-prune` | Streaming AF-filtering & deletion-bubble capping without `vg`/`XG` | `liao_2023_draft`, `hickey_2023_pangenome`, `kopalli_2025_toward` | Two-pass `W`/`P` rewriting (`horizon_clean.rs`), `src/gfa/mod.rs` |
| **P3** | `horizon eval-null` | Reversed/permuted-sequence empirical FPR + cross-aligner concordance | `earl_2014_alignathon`, `quinlan-hall_2010_bedtools`, `frith_2011_tantan` | `calibrate_null_crc_from_seqs` (`src/mdl/dp.rs`), `horizon_eval.rs` |
| **P3** | `horizon synteny` | Low-RAM chimeric-join, orientation-flip & repeat-spike assembly scan | `shivakumar_2025_mumemto`, `liu_2018_inferring` | `FlatSeedIndex`, LIS chaining (`src/seeding/chain.rs`), SDUST |

---

## Detailed Tool Proposals

### 1. `horizon hal` — Pure-Rust Streaming HAL Toolkit
- **Problem:** Multi-genome WGA ecosystems (Progressive Cactus, Minigraph-Cactus, UCSC/VGP pipelines) store alignments in `.hal` (HDF5), requiring C++ `halTools` (`halStats`, `hal2maf`, `halLiftover`, `hal2vg`) linked against `libhdf5`. Building and deploying `halTools` in containers or HPC environments is brittle, and `horizon_wga --hal-info` currently exposes only root superblock/member inspection.
- **Proposed Subcommands:**
  - `horizon hal stats <file.hal>` — genome sizes, top/bottom segment counts, contig counts, and tree topology.
  - `horizon hal to-maf <file.hal> --ref <genome> [--target-genomes ...] [--region chr:start-end] [--out out.maf.gz]` — streaming clade/region MAF projection.
  - `horizon hal liftover <file.hal> --from <genome> --to <genome> <regions.bed>` — coordinate liftover reusing `src/coord.rs` UCSC strand arithmetic.
  - `horizon hal to-gfa <file.hal> [--no-ancestors] --out out.gfa.zst` — direct HAL-to-GFA v1.1 conversion without `hal2vg | vg convert`.
- **Reused Internals:** `src/hal/hdf5_read.rs` (already verified against 5.5 GiB `4-t2t-orangs-mc-2023v2.full.hal` and 56.4 GiB `jarvis_lab_birds.hal`, including stacked multi-layer zlib deflate filters) and `src/hal/mod.rs` (`HalAlignment::to_maf_for_genome`).
- **Acceptance Criteria:** Byte-accurate or interval-identical `stats`, `to-maf`, and `liftover` output compared against `halTools` on `4-t2t-orangs-mc-2023v2.full.hal` and `jarvis_lab_birds.hal` with zero `libhdf5` runtime dependencies.

---

### 2. `horizon gfa-qc` — Streaming Graph Retention, Clipping & Topology Auditor
- **Problem:**
  - `liao_2023_draft` and `garrison_2024_building` (Supplementary E) document that Minigraph-Cactus clips an average of `188.6 Mb` per non-reference human haplotype and `15–22%` of sequence in *Arabidopsis* / *E. coli* / gibbon, sometimes removing non-centromeric genes.
  - `kopalli_2025_toward` had to combine `odgi paths`, `minimap2`, `samtools`, and `Liftoff` to compute per-assembly retention in crop pangenomes.
  - `andreace_2023_comparing` and `guarracino_2022_odgi` highlight giant connected components linking unrelated chromosomes and `10^5–10^6`-fold node-depth spikes in collapsed repeats.
- **Proposed CLI:**
  ```bash
  horizon gfa-qc <graph.gfa.gz> [--fasta-manifest assemblies.tsv] [--out-prefix qc]
  ```
- **Outputs:**
  - `<out>.paths.tsv`: per-path/haplotype retained bp, unrepresented gap count/bp (when input FASTA lengths or `.fai` are supplied), and chromosome-component assignments.
  - `<out>.components.tsv`: connected components, node/bp totals, and flags for multi-chromosome bridges.
  - `<out>.depth_hist.tsv`: node-depth & haplotype-support (`AF`) histograms (allowing users to preview `vg clip -d` thresholds before pruning) plus top-N extreme-depth repeat hotspots.
- **Reused Internals:** `src/gfa/mod.rs` streaming reader, `src/io_util.rs` transparent `.gz`/`.bgz`/`.zst` I/O, and the low-memory counter architecture of `src/bin/horizon_gfa.rs` (`~97 MB` peak RSS on the `12 GB` HPRC graph).

---

### 3. `horizon inversions` — GFA-Native Inversion & Snarl/Ultrabubble Caller
- **Problem:** `romain_2026_inversions-pangenome-graphs` (`INVPG_annot`) demonstrates that inversions manifest in pangenome graphs in two distinct topologies:
  1. **Path-explicit:** alleles traverse shared nodes in opposite orientations.
  2. **Alignment-rescued:** alleles traverse distinct nodes whose sequences are reverse complements of each other.
  Current workflows require `vg deconstruct` + external `minimap2` alignment of bubble alleles + `bedtools`, and still miss or fragment inversions into dense small-bubble clusters. Meanwhile `sena_2026_bubble-subgraphs-spqr` shows SPQR-tree traversal avoids Paten (2018)'s quadratic worst case when enumerating all snarls/ultrabubbles.
- **Proposed CLI:**
  ```bash
  horizon inversions <graph.gfa.gz> [--min-len 50] [--min-cov 0.5] [--out inversions.bedpe]
  ```
- **Design:**
  - Build `PackedGraph` + Tarjan biedged snarl decomposition (reusing `horizon_resolve.rs`).
  - Classify each snarl/bubble as `path_explicit` (opposite-strand shared node coverage $\ge \text{min-cov}$) or `alignment_rescued` (using `FlatSeedIndex` strand-symmetric `11-of-15` + banded extension directly on candidate allele sequences, eliminating external `minimap2` calls).
  - Emit standardized BEDPE/TSV with breakpoint coordinates, topology class, allele lengths, and haplotype traversal counts.

---

### 4. `horizon gfa-prune` — Streaming Allele-Frequency & Bubble-Size Graph Pruner
- **Problem:** `liao_2023_draft` required `vg clip -d 9 -m 10000` (minimum 10% haplotype support) to make HPRC graphs usable for short-read mapping and `vcfbub -l 0 -r 100000` to drop >100 kb / >10 Mb spurious deletion edges. `kopalli_2025_toward` reported that PGGB graphs failed `Giraffe` indexing due to memory/complexity.
- **Proposed CLI:**
  ```bash
  horizon gfa-prune <in.gfa.gz> --min-haps 2 [--min-af 0.10] [--max-del-bp 100000] --out pruned.gfa.zst
  ```
- **Design:**
  - Pass 1: stream `W`/`P` paths to count per-segment/per-link haplotype support and span.
  - Pass 2: prune unsupported segments/links and oversized deletion shortcuts while rewriting affected `W`/`P` walks (splitting or bridging cleanly using the walk-splitting and CIGAR-overlap logic already implemented in `horizon_clean.rs`).

---

### 5. `horizon eval-null` — Empirical Negative-Control & Concordance Benchmarker
- **Problem:** `earl_2014_alignathon`, `quinlan-hall_2010_bedtools`, `frith_2011_tantan`, `garrison_2024_building`, and `kopalli_2025_toward` all emphasize that real-genome WGA lacks ground truth. Most papers report only self-referential coverage or agreement with another aligner (`NUCMER`/`MULTIZ`), rarely measuring empirical false-positive rates on reversed or shuffled sequences.
- **Proposed CLI:**
  ```bash
  horizon eval-null --query q.fa --target t.fa --aln out.maf [--decoy-aln reversed.maf] [--n-shuffles 20]
  ```
- **Design:**
  - Generate reversed / dinucleotide-preserving shuffled decoy alignments or score existing decoy MAF/PAF outputs.
  - Report empirical FPR (aligned bp/kb on reversed negative controls), GC-stratified bit-score distributions (reusing `calibrate_null_crc_from_seqs` from `src/mdl/dp.rs`), and pairwise-vs-graph alignment concordance alongside `horizon_eval.rs` recall/precision.

---

### 6. `horizon synteny` — Low-Memory Chimeric-Join & Orientation-Flip Assembly Scanner
- **Problem:** `shivakumar_2025_mumemto` showed that multi-MUM collinearity catches real HPRC misassemblies (such as a `chr17–chr19` chimeric join), RagTag orientation flips, and pericentromeric LTR repeat spikes—but full suffix-array + LCP construction required `707–800 GB` RAM on human pangenomes. `liu_2018_inferring` also highlights how anchor-based synteny methods degrade on fragmented assemblies.
- **Proposed CLI:**
  ```bash
  horizon synteny <ref.fa> <assemblies...> [--min-block-bp 10000] [--out synteny_flags.tsv]
  ```
- **Design:**
  - Index reference/assemblies with `FlatSeedIndex` (`12 bytes/entry`) + SDUST repeat masking (`src/seeding/mod.rs`) and run LIS chaining (`src/seeding/chain.rs`).
  - Flag contigs with split multi-chromosome collinear chains (chimeric joins), strand-inverted scaffolds relative to the consensus, and repeat-density spikes—using a few GB of RAM rather than hundreds of GB.

---

## Recommended Implementation Order

1. **`horizon hal` (`stats`, `to-maf`, `liftover`)** — highest immediate leverage because `src/hal/` is already built and verified on multi-gigabyte real HAL files (`jarvis_lab_birds.hal`, `4-t2t-orangs-mc-2023v2.full.hal`).
2. **`horizon gfa-qc`** — natural extension of `horizon_gfa` (`~97 MB` RSS on `12 GB` HPRC GFA) that directly automates the retention/depth/component diagnostics pangenome papers currently script by hand.
3. **`horizon inversions` & `horizon gfa-prune`** — leverage `PackedGraph` + Tarjan snarls (`horizon_resolve.rs`) and walk-rewriting (`horizon_clean.rs`) to replace `vg` + `minimap2` + `bedtools` glue pipelines.
4. **`horizon eval-null` & `horizon synteny`** — round out benchmarking and assembly-QC workflows using existing CRC null calibration and `FlatSeedIndex` chaining.
