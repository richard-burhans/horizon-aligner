# Benchmark Results

This document tracks the precision and recall of `horizon-aligner` across a suite of deep-divergence datasets.

All results were generated out-of-the-box using the default Minimum Description Length (MDL) thresholds and priors in `src/bocpd.rs` and `src/mdl/mod.rs` (note that fixed gap-bit penalties `GAP_OPEN_BITS = -4.5`, `GAP_EXTEND_BITS = -0.8`, and `IDENTITY_PRIOR = 0.65` remain hardcoded, and MAF-extracted FASTAs sample homologous blocks).

## Executive Summary: Performance Highlights

* **Multi-species merging (`horizon_merge`):** Fast and exact at 400+ species scale, with best-score interval netting.
* **Close-to-medium divergence (`horizon_wga`):**
  * UCSC 470-way primates: 82.2% recall / 88.6% precision (Marmoset).
  * UCSC cactus447way primates: 92.5% recall / 85.7% precision (Macaque).
  * UCSC 100-way chrM: 95.3% recall / 100.0% precision (Chimp).
* **Distant divergence (`horizon_wga` & `horizon_guide`):**
  * UCSC cactus447way Laurasiatheria: 47.3% recall / 64.2% precision (Horse).
  * UCSC 100-way chrM: 79.1% recall / 94.2% precision (Mouse).

---

## Tuning Sets
These sets were used during development to verify that the MDL formulations behave correctly.

| Dataset | Reference | Target | Clade | Recall vs Truth MAF | Horizon Bubbles | Alignment Time |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `multiz100way` (chrM) | Human (`hg38`) | Chimp (`panTro4`) | Very Close Primate | **95.3%** | 64 | ~0.04s |
| `multiz100way` (chrM) | Human (`hg38`) | Mouse (`mm10`) | Distant Mammal | **79.1%** | 187 | ~0.18s |

## Held-Out Test Sets
These sets were completely locked away during algorithm development to prevent overfitting.

| Dataset | Reference | Target | Clade | Recall vs Truth MAF | Horizon Bubbles | Alignment Time |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `cactus447way` | Human (`hg38`) | Horse (`Equus_caballus`) | Distant Mammal | **47.3%** | 146 | ~0.19s |
| `cactus447way` | Human (`hg38`) | Macaque (`Macaca_mulatta`) | Close Primate | **92.5%** | 275 | ~0.26s |
| `multiz470way` | Human (`hg38`) | Marmoset (`HLcalPym1`) | Medium Primate | **82.2%** | 497 | ~0.29s |
| `cactus241way` (chr22_random) | Human (`hg38`) | Dog (`Canis_lupus`) | Distant Mammal | **45.3%** | 233 | ~0.66s |
| `cactus241way` (chr22 1MB) | Human (`hg38`) | Dog (`Canis_lupus`) | Distant Mammal | **9.9%** | 494 | ~0.68s |
| B10K 363-way | Chicken (`galGal6`) | Zebra Finch (`taeGut2`) | Deep Avian | *(Pending Data)* | - | - |
| HPRC (v1.1) | `GRCh38` | `CHM13` | Human Pangenome | *(Graph Eval)* | - | - |

> **Note on `cactus241way` (chr22 1MB):** Recall drops 4.5x from `45.3%` (`chr22_random`, 38 kb) to `9.9%` (`chr22 1MB`, 494 Horizon Bubbles) due to aggressive MDL abstention across divergent nuclear regions.

