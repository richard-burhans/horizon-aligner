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

