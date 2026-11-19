# Benchmark Results

This document tracks the precision, recall, runtime, and memory usage of `horizon-aligner` across deep-divergence WGA datasets and the 12 GB HPRC pangenome graph.

All results were generated out-of-the-box using the divergence-adaptive Pair-HMM (`HmmParams` estimated from LCB anchors with Beta prior shrinkage in `src/mdl/mod.rs`), BOCPD compositional profiles (`src/bocpd.rs`), and normalized Pair-HMM transitions (`src/mdl/dp.rs`). Note that `GAP_EXTEND_BITS = -0.8` and the Karlin-Altschul `KARLIN_MDL_COEFF = 0.55` multiplier remain fixed defaults, and MAF-extracted FASTAs sample homologous blocks.

## Executive Summary: Performance Highlights

* **Multi-species merging (`horizon_merge`):** Fast and exact at 400+ species scale, with best-score interval netting.
* **Close-to-medium divergence (`horizon_wga`):**
  * UCSC 470-way primates: **89.6%** recall / **88.0%** precision (Marmoset).
  * UCSC cactus447way primates: **92.6%** recall / **85.4%** precision (Macaque).
  * UCSC 100-way chrM: **95.3%** recall / **100.0%** precision (Chimp).
* **Distant divergence (`horizon_wga` & `horizon_guide`):**
  * UCSC cactus447way Laurasiatheria: **50.5%** recall / **64.1%** precision (Horse).
  * UCSC cactus241way `chr22 1MB`: **64.3%** recall / **61.0%** precision (Dog, across 129 target scaffolds; **62.6% F1** in **15.0s**).
  * UCSC 100-way chrM: **81.8%** recall / **94.5%** precision (Mouse).
* **12 GB HPRC v1.0 Minigraph-Cactus Graph (`horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_wga`):**
  * Full-graph streaming (`horizon_gfa stats`): **88.84s**, **97.1 MB** peak RSS across `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks.
  * Two-pass chromosome extraction (`horizon_gfa extract` for `GRCh38.chr22` 50.8 Mbp + `CHM13.chr22` 35.6 Mbp): **166.73s**, **456.1 MB** peak RSS (`921,252` cached segments).
  * Two-pass graph liftover (`horizon_lift --bridge-bp 10` on 500 `GRCh38.chr22` 1-kb windows -> `CHM13.chr22`): **256.60s**, **619.2 MB** peak RSS (`469` exact, `1` partial, `19` split, `11` unmapped).
  * Pairwise alignment (`horizon_wga --rbest --resolve-dups` on `GRCh38.chr22` vs `CHM13.chr22`): **154.92s**, **24,365,611** aligned pairs (`68.5%` of `CHM13.chr22` q-arm), `60` fwd / `4` rev LCBs, `9,690` Horizon Bubbles.
  * MHC / HLA 5 Mbp locus (`GRCh38.chr6:28.5M-33.5M` vs synthetic C4 tandem-dup haplotype): **31.76s**, **6,609,627** aligned pairs, `543` fwd / `355` rev LCBs, `24,984` Horizon Bubbles, `# EG:Z:` MDL e-graph paralog resolution.
  * `chr6` q-arm 20 Mbp (`GRCh38.chr6:65M-85M` vs `CHM13.chr6:61.4M-81.4M`): **80.21s**, **15,084,227** aligned pairs (`75.4%`), `36` fwd / `41` rev LCBs, `17,317` Horizon Bubbles, `56` e-graph paralog bins. (Full 170.8 Mbp `chr6` hit ~24.3 GB RAM in `pair_set` before completion and was cancelled on the 31.4 GB host, highlighting `HashSet<(usize, usize)>` scaling at >100 Mbp.)

---

## Tuning Sets
These sets were used during development to verify that the MDL formulations behave correctly.

| Dataset | Reference | Target | Clade | Recall vs Truth MAF | Precision | Horizon Bubbles | Alignment Time |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `multiz100way` (chrM) | Human (`hg38`) | Chimp (`panTro4`) | Very Close Primate | **95.3%** | 100.0% | 64 | ~0.04s |
| `multiz100way` (chrM) | Human (`hg38`) | Mouse (`mm10`) | Distant Mammal | **81.8%** | 94.5% | 177 | ~0.23s |

## Held-Out Test Sets
These sets were completely locked away during algorithm development to prevent overfitting.

| Dataset | Reference | Target | Clade | Recall vs Truth MAF | Precision | Horizon Bubbles | Alignment Time |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `cactus447way` | Human (`hg38`) | Horse (`Equus_caballus`) | Distant Mammal | **50.5%** | 64.1% | 157 | ~0.29s |
| `cactus447way` | Human (`hg38`) | Macaque (`Macaca_mulatta`) | Close Primate | **92.6%** | 85.4% | 280 | ~0.43s |
| `multiz470way` | Human (`hg38`) | Marmoset (`HLcalPym1`) | Medium Primate | **89.6%** | 88.0% | 536 | ~0.40s |
| `cactus241way` (chr22_random) | Human (`hg38`) | Dog (`Canis_lupus`) | Distant Mammal | **57.9%** | 75.9% | 266 | ~0.31s |
| `cactus241way` (chr22 1MB) | Human (`hg38`) | Dog (`Canis_lupus`) | Distant Mammal | **64.3%** | **61.0%** | 3,362 | 15.02s |
| B10K 363-way | Chicken (`galGal6`) | Zebra Finch (`taeGut2`) | Deep Avian | *(Pending Data)* | - | - | - |

## Pangenome Graph Benchmarks (12 GB `hprc-v1.0-mc-grch38.gfa.gz` & Subgraphs)

| Stage | Binary | Scope | Result | Peak RSS | Wall Time |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Full BGZF Graph Stream | `horizon_gfa stats` | `81.4M` segments (`3.29 Gbp`), `113.0M` links, `2,580` `P`, `24,456` `W` | Full pass | **97.1 MB** | 88.84s |
| Two-Pass Path/Walk Extract | `horizon_gfa extract` | `GRCh38.chr22` (`50.8 Mbp`) + `CHM13.chr22` (`35.6 Mbp`) | `921,252` cached segs | **456.1 MB** | 166.73s |
| Two-Pass Path/Walk Extract | `horizon_gfa extract` | `GRCh38.chr6` (`170.8 Mbp`) + `CHM13.chr6` (`110.7 Mbp`) | `3,057,602` cached segs | **1,255.1 MB** | 244.93s |
| Two-Pass Graph Liftover | `horizon_lift --bridge-bp 10` | 500 `GRCh38.chr22` 1-kb windows -> `CHM13.chr22` | `469` exact, `1` partial, `19` split, `11` unmapped | **619.2 MB** | 256.60s |
| Chr22 Haplotype Alignment | `horizon_wga --rbest --resolve-dups` | `GRCh38.chr22` vs `CHM13.chr22:15729460-51324926` | `24,365,611` pairs (`9,690` bubbles) | - | 154.92s |
| MHC/HLA 5 Mbp Paralog Test | `horizon_wga --rbest --resolve-dups` | `GRCh38.chr6:28.5M-33.5M` vs synthetic C4 dup | `6,609,627` pairs (`24,984` bubbles) | - | 31.76s |
| Chr6 q-arm 20 Mbp Alignment | `horizon_wga --rbest --resolve-dups` | `GRCh38.chr6:65M-85M` vs `CHM13.chr6:61.4M-81.4M` | `15,084,227` pairs (`17,317` bubbles) | - | 80.21s |
| GFA Junction Cleaner | `horizon_clean` | `/tmp/hprc_chr22.gfa` (`2,405` links) | `2,405` kept, `0` severed | **5.0 MB** | <0.01s |
| DAG/POA MDL Alignment | `horizon_wga --graph` | 4-node diamond GFA + 96-bp query | `>n0>n2>n3`, `96/96` bp (`124.4` bits) | - | 0.0024s |
