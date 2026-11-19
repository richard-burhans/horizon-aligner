# Benchmark Results

This document tracks the precision, recall, runtime, and memory usage of `horizon-aligner` across deep-divergence WGA datasets and the 12 GB HPRC pangenome graph.

All results were generated out-of-the-box using the default Minimum Description Length (MDL) thresholds and priors in `src/bocpd.rs`, `src/mdl/mod.rs`, and normalized Pair-HMM transitions in `src/mdl/dp.rs` (note that fixed gap-bit penalties `GAP_OPEN_BITS = -4.5`, `GAP_EXTEND_BITS = -0.8`, and `IDENTITY_PRIOR = 0.65` remain hardcoded, and MAF-extracted FASTAs sample homologous blocks).

## Executive Summary: Performance Highlights

* **Multi-species merging (`horizon_merge`):** Fast and exact at 400+ species scale, with best-score interval netting.
* **Close-to-medium divergence (`horizon_wga`):**
  * UCSC 470-way primates: **89.4%** recall / **87.4%** precision (Marmoset).
  * UCSC cactus447way primates: **92.5%** recall / **85.3%** precision (Macaque).
  * UCSC 100-way chrM: **95.4%** recall / **100.0%** precision (Chimp).
* **Distant divergence (`horizon_wga` & `horizon_guide`):**
  * UCSC cactus447way Laurasiatheria: **50.3%** recall / **63.7%** precision (Horse).
  * UCSC cactus241way `chr22 1MB`: **65.6%** recall / **53.7%** precision (Dog, across 129 target scaffolds).
  * UCSC 100-way chrM: **81.7%** recall / **94.1%** precision (Mouse).
* **12 GB HPRC v1.0 Minigraph-Cactus Graph (`horizon_gfa` + `horizon_wga`):**
  * Full-graph streaming (`horizon_gfa stats`): **88.84s**, **97.1 MB** peak RSS across `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks.
  * Two-pass chromosome extraction (`horizon_gfa extract` for `GRCh38.chr22` 50.8 Mbp + `CHM13.chr22` 35.6 Mbp): **166.73s**, **456.1 MB** peak RSS (`921,252` cached segments).
  * Pairwise alignment (`horizon_wga --rbest --resolve-dups` on `GRCh38.chr22` vs `CHM13.chr22`): **154.92s**, **24,365,611** aligned pairs (`68.5%` of `CHM13.chr22` q-arm), `60` fwd / `4` rev LCBs, `9,690` Horizon Bubbles.

---

## Tuning Sets
These sets were used during development to verify that the MDL formulations behave correctly.

| Dataset | Reference | Target | Clade | Recall vs Truth MAF | Precision | Horizon Bubbles | Alignment Time |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `multiz100way` (chrM) | Human (`hg38`) | Chimp (`panTro4`) | Very Close Primate | **95.4%** | 100.0% | 63 | ~0.05s |
| `multiz100way` (chrM) | Human (`hg38`) | Mouse (`mm10`) | Distant Mammal | **81.7%** | 94.1% | 173 | ~0.23s |

## Held-Out Test Sets
These sets were completely locked away during algorithm development to prevent overfitting.

| Dataset | Reference | Target | Clade | Recall vs Truth MAF | Precision | Horizon Bubbles | Alignment Time |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `cactus447way` | Human (`hg38`) | Horse (`Equus_caballus`) | Distant Mammal | **50.3%** | 63.7% | 153 | ~0.32s |
| `cactus447way` | Human (`hg38`) | Macaque (`Macaca_mulatta`) | Close Primate | **92.5%** | 85.3% | 279 | ~0.46s |
| `multiz470way` | Human (`hg38`) | Marmoset (`HLcalPym1`) | Medium Primate | **89.4%** | 87.4% | 539 | ~0.36s |
| `cactus241way` (chr22_random) | Human (`hg38`) | Dog (`Canis_lupus`) | Distant Mammal | **57.9%** | 75.9% | 266 | ~0.31s |
| `cactus241way` (chr22 1MB) | Human (`hg38`) | Dog (`Canis_lupus`) | Distant Mammal | **65.6%** | 53.7% | 4,285 | ~52.5s |
| B10K 363-way | Chicken (`galGal6`) | Zebra Finch (`taeGut2`) | Deep Avian | *(Pending Data)* | - | - | - |

## Pangenome Graph Benchmarks (12 GB `hprc-v1.0-mc-grch38.gfa.gz`)

| Stage | Binary | Scope | Result | Peak RSS | Wall Time |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Full BGZF Graph Stream | `horizon_gfa stats` | `81.4M` segments (`3.29 Gbp`), `113.0M` links, `2,580` `P`, `24,456` `W` | Full pass | **97.1 MB** | 88.84s |
| Two-Pass Path/Walk Extract | `horizon_gfa extract` | `GRCh38.chr22` (`50.8 Mbp`) + `CHM13.chr22` (`35.6 Mbp`) | `921,252` cached segs | **456.1 MB** | 166.73s |
| Chr22 Haplotype Alignment | `horizon_wga --rbest --resolve-dups` | `GRCh38.chr22` vs `CHM13.chr22:15729460-51324926` | `24,365,611` pairs (`9,690` bubbles) | - | 154.92s |
