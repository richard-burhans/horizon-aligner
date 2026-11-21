# Benchmark Results

This document tracks the precision, recall, runtime, and memory usage of `horizon-aligner` across deep-divergence WGA datasets, multi-species `horizon_merge` / `horizon_guide` runs, and the 12 GB HPRC pangenome graph.

Pairwise WGA results use the divergence-adaptive Pair-HMM (`HmmParams`: `identity_prior`, `gap_open_bits`, and `gap_extend_bits` estimated from LCB anchors with Beta prior shrinkage in `src/mdl/mod.rs`), BOCPD compositional profiles (`src/bocpd.rs`, ignoring non-ACGT `'N'` bases), composition-aware BOCPD LCB gating (`joint_pair_bits - GAP_OPEN_BITS` across all LCB lengths in `src/mdl/mod.rs`), and normalized Pair-HMM transitions (`src/mdl/dp.rs`, `src/mdl/dag.rs`), including `mdl_xdrop_extend_with_params` and `horizon_wga --graph` (`estimate_dag_hmm_params` calibrates `gap_open_bits`/`gap_extend_bits` from adjacent seed pairs, and `align_seq_to_dag_seeded` uses bp-budgeted component-isolated subgraph extraction). Note that `KARLIN_MDL_COEFF = 0.55` remains a fixed constant, and MAF-extracted FASTAs sample homologous blocks (`horizon_eval` evaluates reference-MAF `(ref_id, tgt_id)` contig pairs individually rather than all-vs-all across scaffolds). All benchmarks below are pinned to their respective commits (`8cd5db0` / `2896229` for core WGA & HPRC, `42700e1` / `9dfec4a` for VGP `horizon_filter` & `horizon_locus`), including a head-to-head `vg v1.60.0` / `odgi` / `gfatools` comparison on `10k`–`200k` HPRC subgraphs (`bench_external_graph_baselines.py`), a head-to-head `FastGA` / `minimap2` comparison on `5 Mbp` `chr6` pericentromere (`bench_external_wga_baselines.py`; `LASTZ` timed out at `15 s`), the **B10K 363-way** avian slice (`Gallus_gallus` vs `Taeniopygia_guttata`), and a base-level coordinate concordance evaluation of full `GRCh38.chr6` × `CHM13.chr6` against both the HPRC Minigraph-Cactus graph and `minimap2 -x asm5` (`v2.28-r1209`).

## Executive Summary: Performance Highlights

* **Multi-species merging & progressive guide tree (`horizon_merge`, `horizon_guide`):**
  * `horizon_merge` across `92` vertebrate species (`hg38` + `91` targets from `chrM.maf`): **21.78s** (`21.82s` wall), **26.6 MB** peak RSS, `460` total / `460` multi-species MAF blocks.
  * `horizon_guide --progressive --hal --rbest` across `5` mammals (`hg38`, `Macaca_mulatta`, `HLcalPym1`, `Equus_caballus`, `Canis_lupus`): **0.859s**, `1,109` MAF blocks (`967` multi-species) at `2896229` (`657–707` blocks at `5cb83fb`); on 5-mammal UCSC `chrM.maf` (`verify_merge_guide_maf.py`, `5cb83fb`), **0.258s**, progressive ancestral-node HDF5 `v0` HAL (`((((hg38,panTro4)Anc0,rheMac3)Anc1,mm10)Anc2,canFam3)Anc3;`, `2,103` bottom segments).
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
| 5-Species Progressive + HAL | `horizon_guide --progressive --hal --rbest` | `hg38`, `Macaca`, `HLcalPym1`, `Equus`, `Canis`; UCSC `chrM.maf` | `1,109` MAF blocks (`2896229`; `657–707` at `5cb83fb`); `2,103` HAL bottom segs on `chrM.maf` (`5cb83fb`, `0.258s`) | — | 0.859s / 0.258s |

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

## VGP / GenomeArk PAF Filtering & Regional Locus Benchmarks (`42700e1` / `9dfec4a`)

| Stage | Binary | Scope | Result | Peak RSS | Wall Time |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `chr6` FastGA PAF Netting (seq-free) | `horizon_filter --net` | `93 MB` (`56,857` `FastGA` records, `GRCh38.chr6` × `CHM13.chr6`) | `565` records (`110.00 Mbp` `aln_len`); **99.12%** exact recall (`109,119,624` bp), **99.65%** cond. / **99.32%** uncond. exact precision (`99.92%` `≤1 kb` precision) vs 2D HPRC | — | **0.45s–0.48s** (`0.50s` `.paf.zst`) |
| `chr6` FastGA PAF SDUST+MDL Filter | `horizon_filter --net --query --target` | `56,857` `FastGA` records with `GRCh38.chr6` & `CHM13.chr6` FASTAs | `553` post-netting records (`436` pre-netting filtered by MDL; removes `12` p-arm `< 46 Mbp` records outside q-arm `1-to-1` truth while preserving `100%` of exact `109,119,624` bp coordinates: **99.12%** exact recall, **99.65%** cond. / **99.32%** uncond. exact precision) vs. `sweepga v0.2.0 --num-mappings 1:1 -t 8` (`0.376s`, `1,085` records, `99.20%` exact recall, `98.86%` cond. / `98.62%` uncond. exact precision, `99.16%` `≤1 kb` precision) | — | 4.68s–4.95s (`0.376s` `sweepga`) |
| `10k`–`200k` HPRC Subgraph Baselines | `bench_external_graph_baselines.py` | `10,000` / `50,000` / `200,000` HPRC `chr22` segments | `200k`: `horizon_gfa stats` `0.058s` / `12.0 MB`, `horizon_resolve` `0.258s` / `95.5 MB`, `horizon_clean` `0.346s` / `130.6 MB`, `horizon_wga --graph` `1.561s` / `501.4 MB` vs `gfatools stat` (`0.129s`, `91.3 MB`), `odgi build` (`0.279s`, `78.2 MB`), `odgi stats -S` (`0.097s`, `60.3 MB`), `vg stats` (`0.564s`, `43.9 MB`) | `12.0–501.4 MB` | 0.058s–1.561s |
| `5 Mbp` `chr6` Pericentromere WGA | `bench_external_wga_baselines.py` | `GRCh38.chr6` × `CHM13.chr6` (`5 Mbp`) | `horizon_wga` (`1.520s`, `876.6 MB`) / `--one-to-one` (`1.519s`, `897.5 MB`) vs `FastGA -T8` (`2.894s`, `321.7 MB`) and `minimap2 -x asm5 -t8` (`2.062s`, `437.1 MB`); `LASTZ` timed out at `15s` | `876.6 MB` | **1.520s** |
| GenomeArk VGP Koala × Amphioxus PAF | `horizon_filter --net` | `GCA_003287225.2` × `GCA_015852565.1` (`36,601` raw `FastGA` records, pre-`9dfec4a`) | `1,976` records (`71.6%` exact overlap with VGP `FastGA_11`'s `1,535` records) | — | **0.10s** |
| Koala × Amphioxus 25 kb Locus Audit | `horizon_filter --query --target` | VGP `FastGA_11` #1 longest hit (`3,548 bp`, `2,646` matches) | Rejected (`md:f:-210.51 hb:i:1`; caught `GA`/`AT`/`CA` simple tandem repeat) | — | <0.05s |
| Koala × Amphioxus 25 kb Locus Align | `horizon_locus` | Extracted 25 kb Koala × Amphioxus locus (`42700e1`) | `9,697 / 25,000` bp (`38.8%`) aligned; emits `similarity.tsv`, `.maf`, `.gfa` | **40.1 MB** | **0.26s** |
| Avian 2-Way Regional Locus | `horizon_locus` (`4` threads) | `galGal` (`386.7 kb`) × `taeGut` (`349.4 kb`, `42700e1`) | `231,444` aligned bp, `p_dist = 0.192552`, `jc69_dist = 0.222528` (aligned core), `1,364` MAF blocks | **236.7–239.3 MB** | **0.83s–0.86s** |
| Human 4-Way 25 kb Regional Locus | `horizon_locus` (`4` threads) | `grch38_mhc`, `hla_dup`, `grch38_peri`, `chm13_peri` (`6` pairs, `42700e1` / `6c1a546`) | `102` MAF blocks (`92` multi-sp); `grch38_mhc` vs `hla_dup` `p_dist = 0.041539` (all `6` pairs) vs. `impg v0.5.0 similarity --gfa-engine poa` (`11.593s` `impg` + `0.59s` `minimap2` = `12.183s`, `mhc`×`hla_dup` only) | **171.8 MB** | **0.41s–0.46s** (`12.183s` `impg`) |
