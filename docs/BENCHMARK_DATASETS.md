# Benchmark Datasets

This document tracks the large-scale multiple alignment datasets used to benchmark `horizon-aligner`'s recall and precision, especially focusing on distant divergence testing (e.g., Laurasiatheria, distant CoVs).

## UCSC GoldenPath Datasets
These are specific slice files downloaded to evaluate mammalian and vertebrate alignment capabilities.

- **UCSC `cactus447way`**
  - **File:** `chr22_KI270735v1_random.maf.gz`
  - **Source:** `http://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus447way/maf/`
  - **Purpose:** Testing Progressive Cactus's latest 447-mammal alignment graph, specifically challenging the aligner with Laurasiatheria and other distant clades.

- **UCSC `hg38` 470-way**
  - **File:** `chr22_KI270735v1_random.maf`
  - **Source:** `http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/maf/`
  - **Purpose:** Comparing against the classical MULTIZ/TBA approach for deep mammalian phylogeny.

- **UCSC `hg38` 100-way**
  - **File:** `chrM.maf.gz`
  - **Source:** `http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/maf/`
  - **Purpose:** Testing 92-vertebrate mitochondrial alignments (highly divergent, circular genomes).

## Ensembl Datasets
- **Ensembl `92_mammals.epo_extended`**
  - **File:** `92_mammals.epo_extended.1_14.maf.gz`
  - **Release:** Ensembl release-115
  - **Purpose:** Evaluating against the Ensembl EPO (Enredo, Pecan, Ortheus) pipeline which computes posterior match probabilities.

- **Ensembl `sars-cov-2.hal`**
  - **File:** `sars-cov-2.hal`
  - **Purpose:** Testing the pure-Rust HDF5 `.hal` parser and evaluating distant viral divergence (SARS-CoV-1, MERS-CoV, HCoV-OC43, Rousettus bat CoV).

## Future Evaluation Targets
- Zoonomia `241-way` (`cactus241way`)
- B10K `363-way` avian (`b10k.genomics.cn`)
- HPRC pangenome `.gfa` graphs

## Recovery Script
To quickly re-download the core UCSC slice files after a sandbox reset, run:
```bash
bash scripts/fetch_benchmarks.sh
```
