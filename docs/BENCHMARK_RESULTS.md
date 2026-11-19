# Benchmark Results

This document tracks the precision, recall, runtime, and memory usage of `horizon-aligner` across deep-divergence WGA datasets and the 12 GB HPRC pangenome graph.

Pairwise WGA results use the divergence-adaptive Pair-HMM (`HmmParams`: `identity_prior`, `gap_open_bits`, and `gap_extend_bits` estimated from LCB anchors with Beta prior shrinkage in `src/mdl/mod.rs`), BOCPD compositional profiles (`src/bocpd.rs`, ignoring non-ACGT `'N'` bases), and normalized Pair-HMM transitions (`src/mdl/dp.rs`, `src/mdl/dag.rs`), including `mdl_xdrop_extend_with_params` and `horizon_wga --graph` (`estimate_dag_hmm_params` now also calibrates `gap_open_bits`/`gap_extend_bits` from adjacent seed pairs). Note that `KARLIN_MDL_COEFF = 0.55` remains a fixed constant, and MAF-extracted FASTAs sample homologous blocks (`horizon_eval` evaluates reference-MAF `(ref_id, tgt_id)` contig pairs individually rather than all-vs-all across scaffolds). All 7 WGA rows and `horizon_gfa stats`/`extract` (`chr22`), `horizon_lift`, `horizon_resolve`, `horizon_wga --graph`, and `horizon_clean --gfa -` below were re-measured on the current release build (head-to-head `minimap2`/`odgi`/`vg`/`LASTZ` baselines remain pending external tool availability on this host).

## Executive Summary: Performance Highlights

* **Multi-species merging (`horizon_merge`):** Implemented in `src/bin/horizon_merge.rs` (best-score interval netting; unit-tested; no multi-species benchmark row yet).
* **Close-to-medium divergence (`horizon_wga`):**
  * UCSC 470-way primates: **88.8%** recall / **86.1%** precision (Marmoset, `0.411s`).
  * UCSC cactus447way primates: **94.4%** recall / **84.6%** precision (Macaque, `0.437s`).
  * UCSC 100-way chrM: **95.5%** recall / **98.8%** precision (Chimp, `0.139s`).
* **Distant divergence (`horizon_wga` & `horizon_guide`):**
  * UCSC cactus447way Laurasiatheria: **53.0%** recall / **65.4%** precision (Horse, `0.296s`).
  * UCSC cactus241way `chr22 1MB`: **62.3%** recall / **60.3%** precision (Dog, across 31 target contigs evaluated per-pair vs reference MAF; **61.3% F1** in **14.15s**).
  * UCSC 100-way chrM: **78.8%** recall / **93.1%** precision (Mouse, `0.179s`).
* **12 GB HPRC v1.0 Minigraph-Cactus Graph (`horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_resolve`, `horizon_wga`):**
  * Full-graph streaming (`horizon_gfa stats`, counter-only pass on current build): **99.97s**, **97.0 MB** peak RSS across `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks.
  * Two-pass chromosome extraction (`horizon_gfa extract` for `GRCh38.chr22` 50.8 Mbp + `CHM13.chr22` 35.6 Mbp): **183.29s**, **456.4 MB** peak RSS (`921,252` cached segments).
  * Two-pass graph liftover (`horizon_lift --bridge-bp 10` on 500 `+`-strand `GRCh38.chr22:16.0M-49.9M` 1-kb windows -> `CHM13.chr22`, excluding `0-16.0M` p-arm/centromere): **182.67s**, **619.7 MB** peak RSS (`472` exact including `<=10 bp` bridged gaps, `5` partial, `15` split, `8` unmapped; no external coordinate truth oracle).
  * Pairwise alignment (`horizon_wga --resolve-dups` on `GRCh38.chr22` vs `CHM13.chr22`, pre-`0e70693`): **154.92s**, **24,365,611** many-to-many `(x, y)` aligned pairs (non-deduplicated ratio `68.5%` vs `35.6 Mbp` `CHM13.chr22` q-arm), `60` fwd / `4` rev LCBs, `9,690` Horizon Bubbles.
  * Real HPRC `chr22` Pericentromere 5 Mbp (`GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M`): **14.57s**, **1,861.9 MB** peak RSS, **5,106,541** many-to-many `(x, y)` aligned pairs (`1.164x` multiplicity; `unique_x_bp=4,387,655`, `unique_y_bp=4,549,858`; no external truth baseline), `16,083` Horizon Bubbles.
  * Real HPRC `chr6` Pericentromere 5 Mbp (`GRCh38.chr6:61.0M-66.0M` vs `CHM13.chr6:61.4M-66.4M`): **9.82s**, **1,048.2 MB** peak RSS, **4,247,511** many-to-many `(x, y)` aligned pairs (`unique_x_bp=4,012,961`, `unique_y_bp=3,998,723`), `23,299` Horizon Bubbles.
  * MHC / HLA 5 Mbp locus (`GRCh38.chr6:28.5M-33.5M` vs synthetic C4 tandem-dup haplotype): **25.08s**, **3,345.9 MB** peak RSS, **6,679,345** aligned pairs, `543` fwd / `355` rev LCBs, `24,827` Horizon Bubbles, `# EG:Z:` MDL e-graph paralog resolution.
  * `chr6` q-arm 20 Mbp (`GRCh38.chr6:65M-85M` vs `CHM13.chr6:61M-81M`): **38.23s** (with `BTreeSet`-indexed `cluster_strand_lcbs`, `seed_end=0.00s`, down from `47.51s`), **2,409.7 MB** peak RSS (`0.120 GB/Mbp`, down from `5,863.3 MB` / `0.293 GB/Mbp` via `FlatSeedIndex` [`12 bytes/entry` post-build, `~28 bytes/entry` build peak] + `rev_offsets`/`bocpd`/`dp` fixes), **15,025,501** many-to-many `(x, y)` aligned pairs (`unique_x_bp=14,999,823`, `unique_y_bp=14,998,488`, `75.0%` deduplicated coverage), `16,877` Horizon Bubbles. (Prior to `FlatSeedIndex`, full `GRCh38.chr6` `170.8 Mbp` vs `CHM13.chr6` `110.7 Mbp` hit `26.5 GB` RSS; a completed full-`chr6` run remains in progress.)

---

## Tuning Sets
These sets were used during development to verify that the MDL formulations behave correctly.

| Dataset | Reference | Target | Clade | Recall vs Reference MAF | Precision | Horizon Bubbles | Alignment Time |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `multiz100way` (chrM) | Human (`hg38`) | Chimp (`panTro4`) | Very Close Primate | **95.5%** | 98.8% | 83 | 0.139s |
| `multiz100way` (chrM) | Human (`hg38`) | Mouse (`mm10`) | Distant Mammal | **78.8%** | 93.1% | 201 | 0.179s |

## Held-Out Test Sets
These sets were completely locked away during algorithm development to prevent overfitting.

| Dataset | Reference | Target | Clade | Recall vs Reference MAF | Precision | Horizon Bubbles | Alignment Time |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `cactus447way` | Human (`hg38`) | Horse (`Equus_caballus`) | Distant Mammal | **53.0%** | 65.4% | 187 | 0.296s |
| `cactus447way` | Human (`hg38`) | Macaque (`Macaca_mulatta`) | Close Primate | **94.4%** | 84.6% | 279 | 0.437s |
| `multiz470way` | Human (`hg38`) | Marmoset (`HLcalPym1`) | Medium Primate | **88.8%** | 86.1% | 561 | 0.411s |
| `cactus241way` (chr22_random) | Human (`hg38`) | Dog (`Canis_lupus_familiaris`) | Distant Mammal | **54.7%** | 82.6% | 212 | 0.221s |
| `cactus241way` (chr22 1MB) | Human (`hg38`) | Dog (`Canis_lupus_familiaris`, 31 contigs) | Distant Mammal | **62.3%** | **60.3%** | 3,490 | 14.15s |
| B10K 363-way | Chicken (`galGal6`) | Zebra Finch (`taeGut2`) | Deep Avian | *(Pending Data)* | - | - | - |

## Pangenome Graph Benchmarks (12 GB `hprc-v1.0-mc-grch38.gfa.gz` & Subgraphs)

| Stage | Binary | Scope | Result | Peak RSS | Wall Time |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Full BGZF Graph Stream | `horizon_gfa stats` | `81.4M` segments (`3.29 Gbp`), `113.0M` links, `2,580` `P`, `24,456` `W` | Counter-only pass | **97.0 MB** | 99.97s |
| Two-Pass Path/Walk Extract | `horizon_gfa extract` | `GRCh38.chr22` (`50.8 Mbp`) + `CHM13.chr22` (`35.6 Mbp`) | `921,252` cached segs | **456.4 MB** | 183.29s |
| Two-Pass Path/Walk Extract (`b857377`) | `horizon_gfa extract` | `GRCh38.chr6` (`170.8 Mbp`) + `CHM13.chr6` (`110.7 Mbp`) | `3,057,602` cached segs | **1,255.1 MB** | 244.93s |
| Two-Pass Graph Liftover | `horizon_lift --bridge-bp 10` | 500 `GRCh38.chr22` 1-kb windows -> `CHM13.chr22` | `472` exact, `5` partial, `15` split, `8` unmapped | **619.7 MB** | 182.67s |
| Chr22 Haplotype Alignment | `horizon_wga --resolve-dups` | `GRCh38.chr22` vs `CHM13.chr22:15729460-51324926` (pre-`0e70693`) | `24,365,611` pairs (`9,690` bubbles) | - | 154.92s |
| Chr22 Pericentromere 5 Mbp | `horizon_wga --resolve-dups` | `GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M` | `5,106,541` pairs (`4.39M`/`4.55M` unique bp, `16,083` bubbles) | **1,861.9 MB** | 14.57s |
| Chr6 Pericentromere 5 Mbp | `horizon_wga --resolve-dups` | `GRCh38.chr6:61.0M-66.0M` vs `CHM13.chr6:61.4M-66.4M` | `4,247,511` pairs (`4.01M`/`4.00M` unique bp, `23,299` bubbles) | **1,048.2 MB** | 9.82s |
| MHC/HLA 5 Mbp Paralog Test | `horizon_wga --resolve-dups` | `GRCh38.chr6:28.5M-33.5M` vs synthetic C4 dup | `6,679,345` pairs (`24,827` bubbles) | **3,345.9 MB** | 25.08s |
| Chr6 q-arm 20 Mbp Alignment | `horizon_wga --resolve-dups` | `GRCh38.chr6:65M-85M` vs `CHM13.chr6:61M-81M` | `15,025,501` pairs (`15.00M`/`15.00M` unique bp, `16,877` bubbles) | **2,409.7 MB** | 38.23s |
| GFA Junction Cleaner (`stdin` spool) | `horizon_clean --gfa -` | `*`-seq pairwise `chr6` peri GFA (`2,776` `sup>0` links) | `2,776` kept, `0` severed (no-op on `*`-seq) | **4.9 MB** | 0.15s |
| Synthetic E-Graph Resolver | `horizon_resolve --collapse-dups` | Synthetic `200`-seg `ACGT` repeat, `20` identical walks (`--divergence 0.02`, `L`-links rewritten) | `280` windows, `840` dups collapsed, `308,760.9` analytical `2*len` bits | **24.5 MB** | 4.58s |
| Seeded DAG/POA Alignment | `horizon_wga --graph` | 4-node diamond GFA (`<=8` node bypass) + 96-bp query | `>n0>n2>n3`, `96/96` bp (`133.2` bits) | - | 0.0021s |
