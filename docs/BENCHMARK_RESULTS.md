# Benchmark Results

This document tracks the precision, recall, runtime, and memory usage of `horizon-aligner` across deep-divergence WGA datasets and the 12 GB HPRC pangenome graph.

Pairwise WGA results use the divergence-adaptive Pair-HMM (`HmmParams`: `identity_prior`, `gap_open_bits`, and `gap_extend_bits` estimated from LCB anchors with Beta prior shrinkage in `src/mdl/mod.rs`), BOCPD compositional profiles (`src/bocpd.rs`), and normalized Pair-HMM transitions (`src/mdl/dp.rs`, `src/mdl/dag.rs`). Note that `mdl_xdrop_extend` (`src/mdl/dp.rs`) still uses static defaults (`IDENTITY_PRIOR`, `GAP_OPEN_BITS`, `GAP_EXTEND_BITS`), `horizon_wga --graph` passes `&HmmParams::default()`, `KARLIN_MDL_COEFF = 0.55` remains a fixed constant, and MAF-extracted FASTAs sample homologous blocks. Earlier WGA rows (`multiz100way`, `cactus447way`, `multiz470way`, `cactus241way chr22_random`, `GRCh38.chr22 vs CHM13.chr22`) were recorded prior to `0e70693`/`b18f25b` and have not yet been re-run on `b18f25b` (`cactus241way chr22 1MB` was re-run on `b18f25b` and shifted from `64.3%` to `63.6%` recall).

## Executive Summary: Performance Highlights

* **Multi-species merging (`horizon_merge`):** Implemented in `src/bin/horizon_merge.rs` (best-score interval netting; no multi-species benchmark row or binary unit tests yet).
* **Close-to-medium divergence (`horizon_wga`):**
  * UCSC 470-way primates: **89.6%** recall / **88.0%** precision (Marmoset).
  * UCSC cactus447way primates: **92.6%** recall / **85.4%** precision (Macaque).
  * UCSC 100-way chrM: **95.3%** recall / **100.0%** precision (Chimp).
* **Distant divergence (`horizon_wga` & `horizon_guide`):**
  * UCSC cactus447way Laurasiatheria: **50.5%** recall / **64.1%** precision (Horse).
  * UCSC cactus241way `chr22 1MB`: **63.6%** recall / **61.0%** precision (Dog, across 129 target scaffolds; **62.3% F1** in **15.1s**).
  * UCSC 100-way chrM: **81.8%** recall / **94.5%** precision (Mouse).
* **12 GB HPRC v1.0 Minigraph-Cactus Graph (`horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_resolve`, `horizon_wga`):**
  * Full-graph streaming (`horizon_gfa stats`): **88.84s**, **97.1 MB** peak RSS across `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks.
  * Two-pass chromosome extraction (`horizon_gfa extract` for `GRCh38.chr22` 50.8 Mbp + `CHM13.chr22` 35.6 Mbp): **166.73s**, **456.1 MB** peak RSS (`921,252` cached segments).
  * Two-pass graph liftover (`horizon_lift --bridge-bp 10` on 500 `+`-strand `GRCh38.chr22:16.0M-49.9M` 1-kb windows -> `CHM13.chr22`, excluding `0-16.0M` p-arm/centromere): **256.60s**, **619.2 MB** peak RSS (`469` exact including `<=10 bp` bridged gaps, `1` partial, `19` split, `11` unmapped; no external coordinate truth oracle).
  * Pairwise alignment (`horizon_wga --resolve-dups` on `GRCh38.chr22` vs `CHM13.chr22`, pre-`0e70693`): **154.92s**, **24,365,611** many-to-many `(x, y)` aligned pairs (non-deduplicated ratio `68.5%` vs `35.6 Mbp` `CHM13.chr22` q-arm), `60` fwd / `4` rev LCBs, `9,690` Horizon Bubbles.
  * Real HPRC `chr22` Pericentromere 5 Mbp (`GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M`): **23.45s**, **3,084.5 MB** peak RSS, **5,894,191** many-to-many `(x, y)` aligned pairs, `246` fwd / `184` rev LCBs, `17,789` Horizon Bubbles (`horizon_clean` on `*`-sequence `sup>0` pairwise GFA: `6,030` links kept, `0` severed, **0.01s**, **7.8 MB** RSS).
  * Real HPRC `chr6` Pericentromere 5 Mbp (`GRCh38.chr6:61.0M-66.0M` vs `CHM13.chr6:61.4M-66.4M`): **17.60s**, **3,084.5 MB** peak RSS, **4,785,880** many-to-many `(x, y)` aligned pairs, `264` fwd / `326` rev LCBs, `24,916` Horizon Bubbles (`horizon_clean` on `*`-sequence `sup>0` pairwise GFA: `2,776` links kept, `0` severed, **0.01s**, **5.0 MB** RSS).
  * MHC / HLA 5 Mbp locus (`GRCh38.chr6:28.5M-33.5M` vs synthetic C4 tandem-dup haplotype): **25.08s**, **3,345.9 MB** peak RSS, **6,679,345** aligned pairs, `543` fwd / `355` rev LCBs, `24,827` Horizon Bubbles, `# EG:Z:` MDL e-graph paralog resolution.
  * `chr6` q-arm 20 Mbp (`GRCh38.chr6:65M-85M` vs `CHM13.chr6:61.4M-81.4M`): **81.33s**, **5,863.3 MB** peak RSS (`0.293 GB/Mbp`), **15,021,911** many-to-many `(x, y)` aligned pairs (`75.1%` non-deduplicated ratio), `36` fwd / `41` rev LCBs, `21,475` Horizon Bubbles. (After replacing `HashSet<(usize, usize)>` with sorted/deduplicated `Vec<(usize, usize)>`, full `GRCh38.chr6` `170.8 Mbp` vs `CHM13.chr6` `110.7 Mbp` reached **26.5 GB** mid-run RSS — `0.094 GB/Mbp` combined, `0.155 GB/Mbp` ref-only — at `47:33` CPU before cancellation on the 30 GiB host due to `ContigIndex` `HashMap<u64, Vec<usize>>` seed maps.)

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
| `cactus241way` (chr22 1MB) | Human (`hg38`) | Dog (`Canis_lupus`) | Distant Mammal | **63.6%** | **61.0%** | 3,411 | 15.11s |
| B10K 363-way | Chicken (`galGal6`) | Zebra Finch (`taeGut2`) | Deep Avian | *(Pending Data)* | - | - | - |

## Pangenome Graph Benchmarks (12 GB `hprc-v1.0-mc-grch38.gfa.gz` & Subgraphs)

| Stage | Binary | Scope | Result | Peak RSS | Wall Time |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Full BGZF Graph Stream | `horizon_gfa stats` | `81.4M` segments (`3.29 Gbp`), `113.0M` links, `2,580` `P`, `24,456` `W` | Full pass | **97.1 MB** | 88.84s |
| Two-Pass Path/Walk Extract | `horizon_gfa extract` | `GRCh38.chr22` (`50.8 Mbp`) + `CHM13.chr22` (`35.6 Mbp`) | `921,252` cached segs | **456.1 MB** | 166.73s |
| Two-Pass Path/Walk Extract | `horizon_gfa extract` | `GRCh38.chr6` (`170.8 Mbp`) + `CHM13.chr6` (`110.7 Mbp`) | `3,057,602` cached segs | **1,255.1 MB** | 244.93s |
| Two-Pass Graph Liftover | `horizon_lift --bridge-bp 10` | 500 `GRCh38.chr22` 1-kb windows -> `CHM13.chr22` | `469` exact, `1` partial, `19` split, `11` unmapped | **619.2 MB** | 256.60s |
| Chr22 Haplotype Alignment | `horizon_wga --rbest --resolve-dups` | `GRCh38.chr22` vs `CHM13.chr22:15729460-51324926` | `24,365,611` pairs (`9,690` bubbles) | - | 154.92s |
| Chr22 Pericentromere 5 Mbp | `horizon_wga --resolve-dups` | `GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M` | `5,894,191` pairs (`17,789` bubbles) | **3,084.5 MB** | 23.45s |
| Chr6 Pericentromere 5 Mbp | `horizon_wga --resolve-dups` | `GRCh38.chr6:61.0M-66.0M` vs `CHM13.chr6:61.4M-66.4M` | `4,785,880` pairs (`24,916` bubbles) | **3,084.5 MB** | 17.60s |
| MHC/HLA 5 Mbp Paralog Test | `horizon_wga --resolve-dups` | `GRCh38.chr6:28.5M-33.5M` vs synthetic C4 dup | `6,679,345` pairs (`24,827` bubbles) | **3,345.9 MB** | 25.08s |
| Chr6 q-arm 20 Mbp Alignment | `horizon_wga --resolve-dups` | `GRCh38.chr6:65M-85M` vs `CHM13.chr6:61.4M-81.4M` | `15,021,911` pairs (`21,475` bubbles) | **5,863.3 MB** | 81.33s |
| GFA Junction Cleaner | `horizon_clean` | `*`-seq pairwise `chr22`/`chr6` peri GFAs (`6,030` / `2,776` `sup>0` links) | `6,030` / `2,776` kept, `0` severed (no-op) | **7.8 MB** / **5.0 MB** | 0.01s |
| Synthetic E-Graph Resolver | `horizon_resolve --collapse-dups` | Synthetic `200`-seg `ACGT` repeat, `20` identical walks (`--divergence 0.02`) | `280` windows, `840` dups collapsed, `308,760.9` analytical `2*len` bits | **24.0 MB** | 4.48s |
| Seeded DAG/POA Alignment | `horizon_wga --graph` | 4-node diamond GFA (`<=8` node bypass) + 96-bp query | `>n0>n2>n3`, `96/96` bp (`124.4` bits) | - | 0.0031s |
