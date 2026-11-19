# Benchmark Results

This document tracks the precision, recall, runtime, and memory usage of `horizon-aligner` across deep-divergence WGA datasets, multi-species `horizon_merge` / `horizon_guide` runs, and the 12 GB HPRC pangenome graph.

Pairwise WGA results use the divergence-adaptive Pair-HMM (`HmmParams`: `identity_prior`, `gap_open_bits`, and `gap_extend_bits` estimated from LCB anchors with Beta prior shrinkage in `src/mdl/mod.rs`), BOCPD compositional profiles (`src/bocpd.rs`, ignoring non-ACGT `'N'` bases), composition-aware BOCPD LCB gating (`joint_pair_bits - GAP_OPEN_BITS` across all LCB lengths in `src/mdl/mod.rs`), and normalized Pair-HMM transitions (`src/mdl/dp.rs`, `src/mdl/dag.rs`), including `mdl_xdrop_extend_with_params` and `horizon_wga --graph` (`estimate_dag_hmm_params` calibrates `gap_open_bits`/`gap_extend_bits` from adjacent seed pairs, and `align_seq_to_dag_seeded` uses bp-budgeted component-isolated subgraph extraction). Note that `KARLIN_MDL_COEFF = 0.55` remains a fixed constant, and MAF-extracted FASTAs sample homologous blocks (`horizon_eval` evaluates reference-MAF `(ref_id, tgt_id)` contig pairs individually rather than all-vs-all across scaffolds). All 5 open workstreams (`R1-IMPL`, `R1-HMM`, `R2-HMM-XDROP`, `R2-RESOLVE-LIFT`, `R2-INDEX-RAM`) have been verified on `HEAD` (`bb43dd9`), including a head-to-head `vg v1.60.0` comparison on a real 10,000-segment HPRC `chr22` subgraph (`minimap2`/`odgi`/`LASTZ` remain unavailable on this host).

## Executive Summary: Performance Highlights

* **Multi-species merging & progressive guide tree (`horizon_merge`, `horizon_guide`):**
  * `horizon_merge` across `92` vertebrate species (`hg38` + `91` targets from `chrM.maf`): **4.451s**, **26.8 MB** peak RSS, `460` total / `460` multi-species MAF blocks.
  * `horizon_guide --progressive --hal --rbest` across `5` mammals (`hg38`, `Macaca_mulatta`, `HLcalPym1`, `Equus_caballus`, `Canis_lupus`): **2.527s**, **610.1 MB** peak RSS, `1,225` MAF blocks (`951` multi-species), `4` HAL bottom segments verified via `horizon_wga --hal-info` (`hdf5_v=0`, `root_members=5`).
* **BOCPD composition-aware LCB gate (`R1-HMM`):**
  * `10,000`-bp random `18%`-GC DNA: **0** false-positive aligned pairs (`0.0068s`, down from `18,077` pairs / `74.9%` coverage before gating); `30%`-GC: `247` pairs (`0.1079s`); `45%`-GC: `149` pairs (`0.1182s`).
* **Close-to-medium divergence (`horizon_wga`):**
  * UCSC 470-way primates: **85.7%** recall / **87.7%** precision (Marmoset, `0.333s`, `569` bubbles).
  * UCSC cactus447way primates: **90.8%** recall / **86.9%** precision (Macaque, `0.324s`, `267` bubbles).
  * UCSC 100-way chrM: **94.4%** recall / **97.4%** precision (Chimp, `0.129s`, `103` bubbles).
* **Distant divergence (`horizon_wga` & `horizon_guide`):**
  * UCSC cactus447way Laurasiatheria: **54.0%** recall / **72.3%** precision (Horse, `0.220s`, `145` bubbles).
  * UCSC cactus241way `chr22 1MB`: **62.3%** recall / **60.3%** precision (Dog, across 31 target contigs evaluated per-pair vs reference MAF; **61.3% F1** in **14.15s**).
  * UCSC 100-way chrM: **76.4%** recall / **93.6%** precision (Mouse, `0.141s`, `200` bubbles).
* **12 GB HPRC v1.0 Minigraph-Cactus Graph (`horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_resolve`, `horizon_wga`):**
  * Full-graph streaming (`horizon_gfa stats`, counter-only pass): **99.97s**, **97.0 MB** peak RSS across `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks (`104.98s` / `97.1 MB` on `33bede8`).
  * Real `10,000`-segment sequence-bearing `chr22` subgraph (`10,937,960` bp, `10,300` links, `5` walks): `horizon_gfa stats` runs in **0.0256s** at **22.7 MB** RSS vs `vg stats -z -l -N` (`4.2984s`, `77.0 MB` RSS — **168x faster**, **3.4x lower RSS**) and `vg convert -g -p` (`1.8292s`, `76.6 MB` RSS); `horizon_clean` runs in **0.0525s** (`36.9 MB` RSS, `10,300` links kept); `horizon_resolve --collapse-dups --tag-walks` runs in **8.106s** (`36.9 MB` RSS, `223` anchors, `755` snarl/bubble windows); `horizon_wga --graph` (`20,000` bidirected DAG nodes) runs in **0.821s–1.102s** (`403.6–404.4 MB` RSS, cleanly abstaining on `400`-bp all-`N` and recovering `122/258` bp and `197/356` bp on real HPRC walk queries).
  * Two-pass chromosome extraction (`horizon_gfa extract` on `HEAD`): `chr22` (`50.8 Mbp` + `35.6 Mbp`, `921,252` segs) in **295.62s**, **455.8 MB** RSS; `chr6` (`170.8 Mbp` + `110.7 Mbp`, `3,057,602` segs) in **399.48s**, **1,254.9 MB** RSS.
  * Two-pass graph liftover (`horizon_lift --from GRCh38 --to CHM13 --to-chrom chr22 --bridge-bp 10` on `500` `chr22` 1-kb windows): **178.93s**, **708.6 MB** peak RSS (`247` exact, `230` partial, `15` split, `8` unmapped).
  * Full `chr6` (`GRCh38.chr6` `170.8 Mbp` vs `CHM13.chr6` `110.7 Mbp` q-arm `61401689-172126628`): **451.35s** (`53130b8`, **10,515.0 MB** RSS, **2.35x** faster vs `1e47fe0` `1,062.45s`); re-verified on `HEAD` at **592.68s** (`12,533.0 MB` RSS, `233,384,527` aligned pairs, `122,218,291` `unique_x_bp` [`71.6%`], `105,962,950` `unique_y_bp` [`95.7%`], `175,249` Horizon Bubbles, `797,103` DP jobs).

---

## Tuning Sets
These sets were used during development to verify that the MDL formulations behave correctly.

| Dataset | Reference | Target | Clade | Recall vs Reference MAF | Precision | Horizon Bubbles | Alignment Time |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `multiz100way` (chrM) | Human (`hg38`) | Chimp (`panTro4`) | Very Close Primate | **94.4%** | 97.4% | 103 | 0.129s |
| `multiz100way` (chrM) | Human (`hg38`) | Mouse (`mm10`) | Distant Mammal | **76.4%** | 93.6% | 200 | 0.141s |

## Held-Out Test Sets
These sets were completely locked away during algorithm development to prevent overfitting.

| Dataset | Reference | Target | Clade | Recall vs Reference MAF | Precision | Horizon Bubbles | Alignment Time |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `cactus447way` | Human (`hg38`) | Horse (`Equus_caballus`) | Distant Mammal | **54.0%** | 72.3% | 145 | 0.220s |
| `cactus447way` | Human (`hg38`) | Macaque (`Macaca_mulatta`) | Close Primate | **90.8%** | 86.9% | 267 | 0.324s |
| `multiz470way` | Human (`hg38`) | Marmoset (`HLcalPym1`) | Medium Primate | **85.7%** | 87.7% | 569 | 0.333s |
| `cactus241way` (chr22_random) | Human (`hg38`) | Dog (`Canis_lupus_familiaris`) | Distant Mammal | **54.7%** | 82.6% | 212 | 0.221s |
| `cactus241way` (chr22 1MB) | Human (`hg38`) | Dog (`Canis_lupus_familiaris`, 31 contigs) | Distant Mammal | **62.3%** | **60.3%** | 3,490 | 14.15s |
| B10K 363-way | Chicken (`galGal6`) | Zebra Finch (`taeGut2`) | Deep Avian | *(Pending Data)* | - | - | - |

## Multi-Species & Progressive HAL Benchmarks

| Stage | Binary | Scope | Result | Peak RSS | Wall Time |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 92-Species Netting Merge | `horizon_merge` | `hg38` + `91` vertebrate targets (`chrM.maf`) | `460` total / `460` multi-species blocks | **26.8 MB** | 4.451s |
| 5-Species Progressive + HAL | `horizon_guide --progressive --hal --rbest` | `hg38`, `Macaca`, `HLcalPym1`, `Equus`, `Canis` | `1,225` MAF blocks (`951` multi-sp), `4` HAL bottom segs | **610.1 MB** | 2.527s |

## Pangenome Graph Benchmarks (12 GB `hprc-v1.0-mc-grch38.gfa.gz` & Subgraphs)

| Stage | Binary | Scope | Result | Peak RSS | Wall Time |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Full BGZF Graph Stream | `horizon_gfa stats` | `81.4M` segments (`3.29 Gbp`), `113.0M` links, `2,580` `P`, `24,456` `W` | Counter-only pass | **97.0 MB** | 99.97s |
| Two-Pass Path/Walk Extract (`HEAD`) | `horizon_gfa extract` | `GRCh38.chr22` (`50.8 Mbp`) + `CHM13.chr22` (`35.6 Mbp`) | `921,252` cached segs | **455.8 MB** | 295.62s |
| Two-Pass Path/Walk Extract (`HEAD`) | `horizon_gfa extract` | `GRCh38.chr6` (`170.8 Mbp`) + `CHM13.chr6` (`110.7 Mbp`) | `3,057,602` cached segs | **1,254.9 MB** | 399.48s |
| Two-Pass Graph Liftover (`--to-chrom`) | `horizon_lift --bridge-bp 10` | 500 `GRCh38.chr22` 1-kb windows -> `CHM13.chr22` | `247` exact, `230` partial, `15` split, `8` unmapped | **708.6 MB** | 178.93s |
| Chr22 Haplotype Alignment | `horizon_wga --resolve-dups` | `GRCh38.chr22` vs `CHM13.chr22:15729460-51324926` | `34,151,373` pairs (`34.03M`/`34.02M` unique bp, `21,759` bubbles) | **3,336.4 MB** | 116.58s |
| Chr22 Pericentromere 5 Mbp | `horizon_wga --resolve-dups` | `GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M` | `5,106,541` pairs (`4.39M`/`4.55M` unique bp, `16,083` bubbles) | **1,861.9 MB** | 14.57s |
| Chr6 Pericentromere 5 Mbp | `horizon_wga --resolve-dups` | `GRCh38.chr6:61.0M-66.0M` vs `CHM13.chr6:61.4M-66.4M` | `4,247,511` pairs (`4.01M`/`4.00M` unique bp, `23,299` bubbles) | **1,048.2 MB** | 9.82s |
| MHC/HLA 5 Mbp Paralog Test | `horizon_wga --resolve-dups` | `GRCh38.chr6:28.5M-33.5M` vs synthetic C4 dup | `6,679,345` pairs (`24,827` bubbles) | **3,345.9 MB** | 25.08s |
| Chr6 q-arm 20 Mbp Alignment | `horizon_wga --resolve-dups` | `GRCh38.chr6:65M-85M` vs `CHM13.chr6:61M-81M` | `15,025,501` pairs (`15.00M`/`15.00M` unique bp, `16,877` bubbles) | **2,409.7 MB** | 38.23s |
| Full Chr6 Alignment (`53130b8`) | `horizon_wga --resolve-dups` | `GRCh38.chr6` (`170.8 Mbp`) vs `CHM13.chr6` (`110.7 Mbp`) | `301,156` fwd / `298,178` rev LCBs (`778,471` DP jobs) | **10,515.0 MB** | **451.35s** |
| Full Chr6 Alignment (`HEAD`) | `horizon_wga --resolve-dups` | `GRCh38.chr6` (`170.8 Mbp`) vs `CHM13.chr6` (`110.7 Mbp`) | `233,384,527` pairs (`71.6%` X, `95.7%` Y, `175,249` bubbles) | **12,533.0 MB** | 592.68s |
| Real 10k HPRC Subgraph Stats | `horizon_gfa stats` vs `vg stats` | `10,000` segs (`10.94 Mbp`), `10,300` links, `5` walks | `horizon_gfa` **168x faster**, **3.4x lower RSS** than `vg v1.60.0` | **22.7 MB** vs `77.0 MB` | **0.0256s** vs `4.2984s` |
| Real 10k HPRC Junction Cleaner | `horizon_clean --min-cp-bits 10.0` | `10,000` segs (`10.94 Mbp`), `10,300` links, `5` walks | `10,300` kept, `0` severed on valid HPRC walks | **36.9 MB** | 0.0525s |
| Real 10k HPRC E-Graph Resolver | `horizon_resolve --collapse-dups --tag-walks` | `10,000` segs (`10.94 Mbp`), `10,300` links, `5` walks | `223` anchors, `755` snarl/bubble windows | **36.9 MB** | 8.106s |
| Real 20k-Node HPRC DAG Alignment | `horizon_wga --graph` | `20,000` bidirected DAG nodes (`10.94 Mbp`) + `258`/`356`/`400`-bp queries | `122/258` bp (`6` bubbles), `197/356` bp (`8` bubbles), `0/400` all-`N` | **404.4 MB** | 0.821s–1.102s |
