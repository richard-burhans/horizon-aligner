# Benchmark Results

This document tracks the precision and recall of `horizon-aligner` across a suite of deep-divergence datasets.

All results were generated out-of-the-box using the parameter-free, Information-Theoretic Minimum Description Length (MDL) threshold defaults.

## Executive Summary: Performance Highlights

* **Multi-species merging (`horizon_merge`):** Fast and exact at 400+ species scale (3.84 s for 406 species in cactus447way, 2.53 s for 317 species in 470-way, 3.07 s for 92 species in 100-way chrM, 1.38 s for 92 species in Ensembl 92_mammals), with best-score interval netting.
* **Close-to-medium divergence (`horizon_wga`):**
  * UCSC 470-way primates: 72.38–93.49% recall / 96.69–99.46% precision.
  * UCSC cactus447way primates: 71.90–77.90% recall / 91.80–93.87% precision.
  * Ensembl sars-cov-2.hal (.hal-only) SARS-CoV-1: 76.87% recall / 98.36% precision.
  * Ensembl 92_mammals (0.72–1.43 Mb): 65.22–71.54% recall / 74.10–96.24% precision.
  * UCSC 100-way chrM: 65.61–88.94% recall / 88.58–99.21% precision.
  * Plasmodium (Pf × Pr): 85.58% recall / 98.33% precision (50.2× faster than Galaxy Growler UCSC pair on Pm × Pb).
* **Distant divergence (`horizon_wga` & `horizon_guide`):**
  * UCSC cactus447way Laurasiatheria: 33.73–41.45% recall / 74.08–78.99% precision.
  * Ensembl sars-cov-2.hal distant CoVs (MERS, OC43, Rousettus_bat_CoV): 36.23–44.63% recall / 58.80–72.93% precision.
  * UCSC hg38 5-way (1 kb promoters): 18.72–26.53% recall / 40.35–57.31% precision.

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

