# horizon-aligner

A whole-genome aligner and multi-species alignment toolkit written in Rust, combining **KP-14 spaced-seed indexing**, **Bayesian Online Change-Point Detection (BOCPD)** compositional segmentation, and **Minimum Description Length (MDL)** banded dynamic programming with explicit abstention (**Horizon Bubbles**).

## Key Features

- **Strand-symmetric spaced-seed & multi-shape seeding** (`src/seeding/`): KP-14 (`11-of-15`) primary indexing with low-complexity repeat filtering, plus transition-tolerant `12-of-19` / `11-of-18`, `w=10`, 1-bit `R/Y`, and 6-frame translated `AA5` fallback seeding in unanchored gaps.
- **BOCPD + MDL banded alignment** (`src/bocpd.rs`, `src/mdl/`): Online GC/composition segmentation coupled with log-odds Viterbi DP and X-drop flank extension. Regions whose MDL bit gain falls below the search-space-scaled threshold abstain as explicit **Horizon Bubbles** rather than emitting spurious alignments.
- **Reference-guided multi-species merging** (`src/merge.rs`): Parallel projection and threading of pairwise MAFs onto a reference coordinate frame with net-style overlap resolution (`horizon_merge`).
- **E-graph synteny & paralog collapse** (`src/egraph.rs`): Equality-saturation rules over signed block expressions to identify and tag tandem/inverted duplications (`EG:Z:` in GFA v1.1).
- **Pure-Rust HDF5 / HAL inspection & decoding** (`src/hal/`): Native streaming superblock, `v1` object-header (including `0x0010` continuation messages), `TREE`/`SNOD`/`HEAP` group traversal, contiguous and zlib-deflated chunked dataset reading, and `SEQNAME_ARRAY` / `BOTTOM_ARRAY` decoding—zero `libhdf5` dependency.

## Repository Structure

```text
src/
├── lib.rs              # Public API re-exports and core unit tests
├── bocpd.rs            # BOCPD compositional segmentation & MDL pair scoring
├── seeding/            # Spaced-seed indexing, transition/AA/R-Y seeding, LIS & LCB chaining
├── mdl/                # Banded MDL+BOCPD DP, X-drop extension, contig-pair pipeline
├── merge.rs            # Multi-species reference-guided MAF merging
├── formats.rs          # MAF v1, GFA v1.1, FASTA I/O, and e-graph paralog formatting
├── egraph.rs           # E-graph synteny / duplication collapse
├── hal/                # Pure-Rust HDF5/HAL reader and HDF5 v0 writer
└── bin/
    ├── horizon_wga.rs    # Pairwise whole-genome aligner & `--hal-info` CLI
    ├── horizon_merge.rs  # Multi-species pairwise MAF -> MAF + GFA merger
    ├── horizon_guide.rs  # Guide-tree / central-hinge multi-species pipeline
    ├── bench_chr1.rs     # Pairwise chromosome benchmark harness
    └── eval_pf_pv.rs     # Pairwise MAF recall/precision evaluator
manuscript/
├── horizon_paper.md / .pdf          # Main manuscript
├── horizon_companion.md / .pdf      # Mathematical & data-structure companion
├── guide_mdl_vs_hoxd70.md / .pdf    # Guide 1: Why ΔL beats HOXD70 in AT/GC extremes
├── guide_pure_rust_hal.md / .pdf    # Guide 2: Pure-Rust Seek HDF5/HAL reader & writer
├── guide_horizon_bubbles.md / .pdf  # Guide 3: What is a Horizon Bubble?
├── guide_making_alignments.md / .pdf # Guide 4: Build genome alignments with Horizon
├── references.bib                   # BibTeX references
├── figures/                         # Publication figures (SVG & PNG)
└── scripts/generate_figures.py      # Reproducible figure generator
```

## Build & Test

Requires Rust (`edition = "2024"`, tested with `rustc 1.100.0`):

```bash
cargo build --release --bins
cargo test --all-features
cargo clippy --all-targets --all-features -- -D warnings
cargo fmt --check
```

To run the ignored real-`.hal` integration test against a downloaded HAL file (e.g., Ensembl `sars-cov-2.hal`):

```bash
HORIZON_TEST_HAL=/path/to/sars-cov-2.hal \
  cargo test real_hal_decode_from_env -- --ignored --nocapture
```

## CLI Usage

### Pairwise Alignment (`horizon_wga`)

```bash
target/release/horizon_wga target.fa query.fa 0 \
  --maf out.maf \
  --gfa out.gfa \
  --resolve-dups
```

Common options:
- `max_contigs`: third positional argument (`0` = all contigs, default `1`).
- `--top-k K` / `--rbest`: restrict contig-pair alignment to top-K (or reciprocal-best) LCB pairs.
- `--respect-softmask`: skip lowercase softmasked bases during fallback seeding.
- `--ry-seeds`: enable un-gated `R/Y` fallback clustering.
- `--hal-info file.hal`: inspect an HDF5 `.hal` file's superblock and root-group genomes without FASTA inputs.

### Multi-Species MAF Merging (`horizon_merge`)

```bash
target/release/horizon_merge <ref_species> <out.maf> <out.gfa> <pair1.maf> [pair2.maf ...]
```

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



## Documentation

- [Manuscript (Markdown)](manuscript/horizon_paper.md) | [PDF](manuscript/horizon_paper.pdf)
- [Technical Companion: Data Formats, Derivations & Worked Trace (Markdown)](manuscript/horizon_companion.md) | [PDF](manuscript/horizon_companion.pdf)
- [Guide 1: Why Bit-Savings (ΔL) Outperforms HOXD70 in Compositionally Extreme Genomes (Markdown)](manuscript/guide_mdl_vs_hoxd70.md) | [PDF](manuscript/guide_mdl_vs_hoxd70.pdf)
- [Guide 2: Reading and Streaming HAL Alignments in Pure Rust Without libhdf5 (Markdown)](manuscript/guide_pure_rust_hal.md) | [PDF](manuscript/guide_pure_rust_hal.pdf)
- [Guide 3: What Is a Horizon Bubble? Calibrated Epistemic Abstention (Markdown)](manuscript/guide_horizon_bubbles.md) | [PDF](manuscript/guide_horizon_bubbles.pdf)

## Architecture
See [ARCHITECTURE.md](ARCHITECTURE.md) for a map of the core abstractions and dependency flow.


## Manual
For a detailed explanation of the algorithmic stages (MDL, BOCPD, multi-shape seeding) and command-line usage, see the [Comprehensive Manual](docs/manual.md).



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



## Documentation

* [Technical Manual](docs/manual.md): Comprehensive documentation covering theory, processing stages, CLI usage, and outputs.
* [Benchmark Datasets](docs/BENCHMARK_DATASETS.md): A complete listing of the UCSC and Ensembl datasets used to evaluate aligner performance across different evolutionary divergences.
* [Benchmark Results](docs/BENCHMARK_RESULTS.md): Zero-shot evaluation metrics on deep divergence datasets using the default Information-Theoretic MDL thresholds. of the UCSC and Ensembl datasets used to evaluate aligner performance across different evolutionary divergences.
