# Benchmark Datasets

This document tracks the large-scale multiple alignment datasets used to benchmark `horizon-aligner`'s recall and precision, especially focusing on distant divergence testing (e.g., Laurasiatheria, distant CoVs).

## UCSC GoldenPath Datasets
These are specific slice files downloaded to evaluate mammalian and vertebrate alignment capabilities.

- **UCSC `cactus447way` (447 Mammals)**
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

- **Zoonomia 241-way (`cactus241way`)**
  - **File:** `chrM.maf.gz` (and other subsets)
  - **Source:** `http://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus241way/`
  - **Purpose:** Evaluating deep mammalian conservation scores and topological abstention using the Zoonomia Project's 241 placental mammals tree.

## Ensembl Datasets
- **Ensembl `92_mammals.epo_extended`**
  - **File:** `92_mammals.epo_extended.1_14.maf.gz`
  - **Source:** `https://ftp.ensembl.org/pub/release-115/maf/ensembl-compara/multiple_alignments/92_mammals.epo_extended/`
  - **Purpose:** Evaluating against the Ensembl EPO (Enredo, Pecan, Ortheus) pipeline which computes posterior match probabilities.

- **Ensembl `sars-cov-2.hal`**
  - **File:** `sars-cov-2.hal`
  - **Source:** `https://ftp.ensemblgenomes.ebi.ac.uk/pub/viruses/fasta/sars_cov_2/`
  - **Purpose:** Testing the pure-Rust HDF5 `.hal` parser and evaluating distant viral divergence (SARS-CoV-1, MERS-CoV, HCoV-OC43, Rousettus bat CoV).

## Avian and Pangenome Datasets
- **B10K 363-way avian**
  - **File:** `b10k_363way.maf`
  - **Source:** `http://b10k.genomics.cn/` (Bird 10,000 Genomes Project)
  - **Purpose:** Testing deep avian and reptilian genomic splits, assessing scaling to hundreds of compact genomes.

- **HPRC Pangenome (Minigraph-Cactus)**
  - **File:** `hprc-v1.0-mc-grch38.gfa.gz`
  - **Source:** `https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/`
  - **Purpose:** Validating `H`/`S`/`L`/`P`/`W` GFA v1.0/v1.1 graph streaming (`horizon_gfa`, `horizon_lift`), 10k–200k subgraph scaling (`bench_hprc_subgraph_scaling.py`, `bench_external_graph_baselines.py`), and whole-chromosome `GRCh38.chr6` × `CHM13.chr6` 1-to-1 concordance (`extract_hprc_chr6_truth.py`, `compare_chr6_horizon_vs_hprc.py`, `score_paf_vs_hprc.py`).

- **VGP / GenomeArk Pairwise PAF & Regional Locus Suites**
  - **Files:** `chr6` `FastGA` PAF (`chr6_fastga.paf` / `.paf.zst`), GenomeArk Koala (`GCA_003287225.2`) × Amphioxus (`GCA_015852565.1`) `FastGA` PAF, and regional multi-FASTA loci (`galGal` × `taeGut`, 4-way human `25 kb` locus).
  - **Purpose:** Evaluating streaming PAF SDUST+MDL filtering/netting (`horizon_filter`) against 2D HPRC 1-to-1 truth and multi-FASTA `impg`-style locus extraction/distance estimation (`horizon_locus`; `bench_vgp_filter_locus.py`).

## Recovery Script
To quickly re-download the core UCSC slice files after a sandbox reset, run:
```bash
bash scripts/benchmarks/fetch_benchmarks.sh
```
