# Benchmark Results

This document tracks the precision, recall, runtime, and memory usage of `horizon-aligner` across deep-divergence WGA datasets, multi-species `horizon_merge` / `horizon_guide` runs, and the 12 GB HPRC pangenome graph.

Pairwise WGA results use the divergence-adaptive Pair-HMM (`HmmParams`: `identity_prior`, `gap_open_bits`, and `gap_extend_bits` estimated from LCB anchors with Beta prior shrinkage in `src/mdl/mod.rs`), BOCPD compositional profiles (`src/bocpd.rs`, ignoring non-ACGT `'N'` bases), composition-aware BOCPD LCB gating (`joint_pair_bits - GAP_OPEN_BITS` across all LCB lengths in `src/mdl/mod.rs`), and normalized Pair-HMM transitions (`src/mdl/dp.rs`, `src/mdl/dag.rs`), including `mdl_xdrop_extend_with_params` and `horizon_wga --graph` (`estimate_dag_hmm_params` calibrates `gap_open_bits`/`gap_extend_bits` from adjacent seed pairs, and `align_seq_to_dag_seeded` uses bp-budgeted component-isolated subgraph extraction). Note that `KARLIN_MDL_COEFF = 0.55` remains a fixed constant, and MAF-extracted FASTAs sample homologous blocks (`horizon_eval` evaluates reference-MAF `(ref_id, tgt_id)` contig pairs individually rather than all-vs-all across scaffolds). All benchmarks below were re-verified on `HEAD`, including a head-to-head `vg v1.60.0` comparison on a real 10,000-segment HPRC `chr22` subgraph, the **B10K 363-way** avian slice (`Gallus_gallus` vs `Taeniopygia_guttata`), and a base-level coordinate concordance evaluation of full `GRCh38.chr6` × `CHM13.chr6` against both the HPRC Minigraph-Cactus graph and `minimap2 -x asm5` (`v2.28-r1209`; `odgi`/`LASTZ` remain unavailable on this host).

## Executive Summary: Performance Highlights

* **Multi-species merging & progressive guide tree (`horizon_merge`, `horizon_guide`):**
  * `horizon_merge` across `92` vertebrate species (`hg38` + `91` targets from `chrM.maf`): **21.78s** (`21.82s` wall), **26.6 MB** peak RSS, `460` total / `460` multi-species MAF blocks.
  * `horizon_guide --progressive --hal --rbest` across `5` mammals (`hg38`, `Macaca_mulatta`, `HLcalPym1`, `Equus_caballus`, `Canis_lupus`): **0.859s**, `1,109` MAF blocks (`967` multi-species), `4` HAL bottom segments verified via `horizon_wga --hal-info` (`hdf5_v=0`, `root_members=5`).
* **BOCPD composition-aware LCB gate (`R1-HMM`):**
  * `10,000`-bp random `18%`-GC DNA: **0** false-positive aligned pairs (`0.0069s`); `30%`-GC: `288` pairs (`0.4150s`); `45%`-GC: `149` pairs (`0.3829s`).
* **Close-to-medium divergence (`horizon_wga`):**
  * UCSC 470-way primates: **86.3%** recall / **87.7%** precision (Marmoset, `0.189s`, `567` bubbles).
  * UCSC cactus447way primates: **91.4%** recall / **87.0%** precision (Macaque, `0.102s`, `267` bubbles).
  * UCSC 100-way chrM: **95.6%** recall / **97.4%** precision (Chimp, `0.058s`, `103` bubbles).
* **Distant divergence (`horizon_wga` & `horizon_guide`):**
  * UCSC cactus447way Laurasiatheria: **54.0%** recall / **67.6%** precision (Horse, `0.068s`, `160` bubbles).
  * UCSC cactus241way `chr22 1MB`: **59.8%** recall / **78.5%** precision (Dog, across 31 target contigs evaluated per-pair vs reference MAF; **67.9% F1** in **0.852s**, down from `14.15s`).
  * UCSC 100-way chrM: **76.5%** recall / **93.7%** precision (Mouse, `0.064s`, `200` bubbles).
  * **B10K 363-way avian** (`chr1:5900000-6400000`, `Gallus_gallus` vs `Taeniopygia_guttata`): **71.1%** recall / **95.1%** precision (`0.642s`, `3,814` bubbles).
* **12 GB HPRC v1.0 Minigraph-Cactus Graph (`horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_resolve`, `horizon_wga`):**
  * Full-graph streaming (`horizon_gfa stats`, counter-only pass): **196.93s**, **96.9 MB** peak RSS across `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks.
  * Real `10,000`-segment sequence-bearing `chr22` subgraph (`10,937,960` bp, `10,300` links, `1` path, `4` walks): `horizon_gfa stats` runs in **0.08s** at **22.6 MB** RSS; `horizon_clean --min-cp-bits 10.0` runs in **0.13s** (`33.4 MB` RSS, `10,300` kept, `0` severed); `horizon_resolve --collapse-dups --tag-walks` runs in **5.571s** (`23.6 MB` RSS, `223` anchors, `755` snarl/bubble windows); `horizon_wga --graph` (`20,000` bidirected DAG nodes) runs in **2.98s–6.09s** (`153.9–154.6 MB` RSS, cleanly abstaining on `400`-bp all-`N` [`0/400` bp] and recovering `122/258` bp and `197/356` bp on real HPRC walk queries).
  * Two-pass chromosome extraction (`horizon_gfa extract` on `8cd5db0`): `chr22` (`50.8 Mbp` + `35.6 Mbp`, `855,439` segs) in **193.89s**, **474.2 MB** RSS; `chr6` (`170.8 Mbp` + `110.7 Mbp`, `2,958,762` segs) in **317.00s**, **1,526.7 MB** RSS.
  * Two-pass graph liftover (`horizon_lift --from GRCh38 --to CHM13 --to-chrom chr22 --bridge-bp 10` on `500` `chr22` 1-kb windows): **177.55s**, **708.8 MB** peak RSS (`472` exact, `5` partial, `15` split, `8` unmapped).
  * Full `chr6` (`GRCh38.chr6` `170.8 Mbp` vs `CHM13.chr6` `110.7 Mbp` q-arm `61401689-172126628`): **258.88s–267.86s** on `8cd5db0` (`2.21×` faster than `bb43dd9`'s `592.68s`), `282,193,564` aligned pairs, `124,160,423` `unique_x_bp` (`72.7%`), `109,464,447` `unique_y_bp` (`98.9%`), `178,737` Horizon Bubbles.
  * **Full `chr6` vs HPRC Minigraph-Cactus graph & `minimap2 -x asm5`:** Against `1,628,018` `1-to-1` shared segments (`110,082,990` bp) between `GRCh38.chr6` and `CHM13.chr6:61401689-172126628`, Horizon's many-to-many MAF (`370,937` blocks, `282,163,046` pairs after fixing `partition_colinear_blocks`) achieves **99.54%** base-coverage recall, **97.96%** exact base-coordinate recall (`107,836,669 / 110,082,990` bp; **98.14%** within `1 kb`), and `48.63%` exact pair precision (`~2.6×` repeat/paralog multiplicity), compared to `minimap2 -x asm5` (`v2.28-r1209`, `46.95s`, `2,947.5 MB` RSS, `1-to-1` primary: **99.71%** exact recall, **99.75%** exact precision; `manuscript/figures/fig7_chr6_hprc_comparison.pdf`).

---

## Tuning Sets
These sets were used during development to verify that the MDL formulations behave correctly.

| Dataset | Reference | Target | Clade | Recall vs Reference MAF | Precision | Horizon Bubbles | Alignment Time |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `multiz100way` (chrM) | Human (`hg38`) | Chimp (`panTro4`) | Very Close Primate | **95.6%** | 97.4% | 103 | 0.058s |
| `multiz100way` (chrM) | Human (`hg38`) | Mouse (`mm10`) | Distant Mammal | **76.5%** | 93.7% | 200 | 0.064s |

## Held-Out Test Sets
These sets were completely locked away during algorithm development to prevent overfitting.

| Dataset | Reference | Target | Clade | Recall vs Reference MAF | Precision | Horizon Bubbles | Alignment Time |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `cactus447way` | Human (`hg38`) | Horse (`Equus_caballus`) | Distant Mammal | **54.0%** | 67.6% | 160 | 0.068s |
| `cactus447way` | Human (`hg38`) | Macaque (`Macaca_mulatta`) | Close Primate | **91.4%** | 87.0% | 267 | 0.102s |
| `multiz470way` | Human (`hg38`) | Marmoset (`HLcalPym1`) | Medium Primate | **86.3%** | 87.7% | 567 | 0.189s |
| `cactus241way` (chr22_random) | Human (`hg38`) | Dog (`Canis_lupus_familiaris`) | Distant Mammal | **48.0%** | 83.3% | 196 | 0.050s |
| `cactus241way` (chr22 1MB) | Human (`hg38`) | Dog (`Canis_lupus_familiaris`, 31 contigs) | Distant Mammal | **59.8%** | **78.5%** | 2,076 | 0.852s |
| B10K 363-way (`chr1:5.9M-6.4M`) | Chicken (`Gallus_gallus`) | Zebra Finch (`Taeniopygia_guttata`) | Deep Avian | **71.1%** | **95.1%** | 3,814 | 0.642s |

## Multi-Species & Progressive HAL Benchmarks

| Stage | Binary | Scope | Result | Peak RSS | Wall Time |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 92-Species Netting Merge | `horizon_merge` | `hg38` + `91` vertebrate targets (`chrM.maf`) | `460` total / `460` multi-species blocks | **26.6 MB** | 21.78s |
| 5-Species Progressive + HAL | `horizon_guide --progressive --hal --rbest` | `hg38`, `Macaca`, `HLcalPym1`, `Equus`, `Canis` | `1,109` MAF blocks (`967` multi-sp), `4` HAL bottom segs | — | 0.859s |

## Pangenome Graph Benchmarks (12 GB `hprc-v1.0-mc-grch38.gfa.gz` & Subgraphs)

| Stage | Binary | Scope | Result | Peak RSS | Wall Time |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Full BGZF Graph Stream | `horizon_gfa stats` | `81.4M` segments (`3.29 Gbp`), `113.0M` links, `2,580` `P`, `24,456` `W` | Counter-only pass | **96.9 MB** | 196.93s |
| Two-Pass Path/Walk Extract | `horizon_gfa extract` | `GRCh38.chr22` (`50.8 Mbp`) + `CHM13.chr22` (`35.6 Mbp`) | `855,439` cached segs | **474.2 MB** | 193.89s |
| Two-Pass Path/Walk Extract | `horizon_gfa extract` | `GRCh38.chr6` (`170.8 Mbp`) + `CHM13.chr6` (`110.7 Mbp`) | `2,958,762` cached segs | **1,526.7 MB** | 317.00s |
| Two-Pass Graph Liftover (`--to-chrom`) | `horizon_lift --bridge-bp 10` | 500 `GRCh38.chr22` 1-kb windows -> `CHM13.chr22` | `472` exact, `5` partial, `15` split, `8` unmapped | **708.8 MB** | 177.55s |
| Chr22 Haplotype Alignment | `horizon_wga --resolve-dups` | `GRCh38.chr22` vs `CHM13.chr22:15729460-51324926` | `34,165,619` pairs (`34.04M`/`34.02M` unique bp, `21,787` bubbles) | **3,971.6 MB** | 69.10s |
| Chr22 Pericentromere 5 Mbp | `horizon_wga --resolve-dups` | `GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M` | `5,158,327` pairs (`4.37M`/`4.57M` unique bp, `14,596` bubbles) | — | 4.682s |
| Chr6 Pericentromere 5 Mbp | `horizon_wga --resolve-dups` | `GRCh38.chr6:61.0M-66.0M` vs `CHM13.chr6:61.4M-66.4M` | `4,248,262` pairs (`3.99M`/`3.98M` unique bp, `22,560` bubbles) | — | 2.213s |
| MHC/HLA 5 Mbp Paralog Test | `horizon_wga --resolve-dups` | `GRCh38.chr6:28.5M-33.5M` vs synthetic C4 dup | `6,148,942` pairs (`4.85M`/`4.95M` unique bp, `21,522` bubbles) | **1,845.4 MB** | 24.15s |
| Chr6 q-arm 20 Mbp Alignment | `horizon_wga --resolve-dups` | `GRCh38.chr6:65M-85M` vs `CHM13.chr6:61M-81M` | `15,019,451` pairs (`14.99M`/`14.99M` unique bp, `16,417` bubbles) | — | 12.84s |
| Full Chr6 Alignment (`8cd5db0`) | `horizon_wga --resolve-dups` | `GRCh38.chr6` (`170.8 Mbp`) vs `CHM13.chr6` (`110.7 Mbp`) | `282,193,564` pairs (`72.7%` X, `98.9%` Y, `178,737` bubbles) | — | **258.88s–267.86s** |
| Full Chr6 vs HPRC Graph Truth | `horizon_wga` vs `1.63M` HPRC segs | `GRCh38.chr6` × `CHM13.chr6:61.4M-172.1M` (`110.08 Mbp` 1-to-1) | **97.96%** exact recall (**99.54%** cov), `48.63%` many-to-many prec. | — | 304.53s |
| Full Chr6 (`minimap2 -x asm5`) | `minimap2 v2.28` vs `1.63M` HPRC segs | `GRCh38.chr6` × `CHM13.chr6:61.4M-172.1M` (`1-to-1` primary) | **99.71%** exact recall, **99.75%** exact precision (`110.20M` pairs) | **2,947.5 MB** | 46.95s |
| Real 10k HPRC Subgraph Stats | `horizon_gfa stats` | `10,000` segs (`10.94 Mbp`), `10,300` links, `1` `P`, `4` `W` | Streaming graph statistics | **22.6 MB** | **0.08s** |
| Real 10k HPRC Junction Cleaner | `horizon_clean --min-cp-bits 10.0` | `10,000` segs (`10.94 Mbp`), `10,300` links, `1` `P`, `4` `W` | `10,300` kept, `0` severed on valid HPRC walks | **33.4 MB** | 0.13s |
| Real 10k HPRC E-Graph Resolver | `horizon_resolve --collapse-dups --tag-walks` | `10,000` segs (`10.94 Mbp`), `10,300` links, `1` `P`, `4` `W` | `223` anchors, `755` snarl/bubble windows | **23.6 MB** | 5.571s |
| Real 20k-Node HPRC DAG Alignment | `horizon_wga --graph` | `20,000` bidirected DAG nodes (`10.94 Mbp`) + `258`/`356`/`400`-bp queries | `122/258` bp (`6` bubbles), `197/356` bp (`8` bubbles), `0/400` all-`N` | **154.6 MB** | 2.98s–6.09s |
