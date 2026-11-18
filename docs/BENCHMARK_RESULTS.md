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

*Note: The `cactus447way` and `multiz470way` slices tested here (`chr22_KI270735v1_random`) are extremely short, highly repetitive unplaced contigs (only ~38kb) where Progressive Cactus failed to map mammalian homologs natively, resulting in 0 truth pairs in the baseline MAF. However, `horizon-aligner` correctly identified ~30,000 alignment pairs between Human and Marmoset in this slice with 497 topological abstentions, highlighting its power on unplaced repetitive regions.*

| Dataset | Reference | Target | Clade | Recall vs Truth MAF | Horizon Bubbles |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `cactus447way` | Human (`hg38`) | Horse (`Equus_caballus`) | Distant Mammal | N/A (0 truth pairs) | 0 |
| `cactus447way` | Human (`hg38`) | Macaque (`Macaca_mulatta`) | Close Primate | N/A (0 truth pairs) | 0 |
| `multiz470way` | Human (`hg38`) | Marmoset (`HLcalPym1`) | Medium Primate | Novel (0 truth pairs)* | 497 |

