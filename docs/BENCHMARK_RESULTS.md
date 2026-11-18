# Benchmark Results

This document tracks the precision and recall of `horizon-aligner` across a suite of deep-divergence datasets.

All results were generated out-of-the-box using the parameter-free, Information-Theoretic Minimum Description Length (MDL) threshold defaults.

## Tuning Sets
These sets were used during development to verify that the MDL formulations behave correctly.

| Dataset | Reference | Target | Clade | Recall vs Truth MAF | Horizon Bubbles |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `multiz100way` (chrM) | Human (`hg38`) | Chimp (`panTro4`) | Very Close Primate | **95.3%** | 64 |
| `multiz100way` (chrM) | Human (`hg38`) | Mouse (`mm10`) | Distant Mammal | **79.1%** | 187 |

## Held-Out Test Sets
These sets were completely locked away during algorithm development to prevent overfitting.

| Dataset | Reference | Target | Clade | Recall vs Truth MAF | Horizon Bubbles |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `cactus447way` | Human (`hg38`) | Horse (`Equus_caballus`) | Distant Mammal | **47.3%** | 146 |
| `cactus447way` | Human (`hg38`) | Macaque (`Macaca_mulatta`) | Close Primate | **92.5%** | 275 |
| `multiz470way` | Human (`hg38`) | Marmoset (`HLcalPym1`) | Medium Primate | **82.2%** | 497 |

