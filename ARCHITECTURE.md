# Horizon Aligner Architecture

This document establishes the conceptual map of the major abstractions in `horizon-aligner` and the contracts between them.

## Tier Map

The architecture of `horizon-aligner` follows a unidirectional dependency flow from the core algorithms up to the output formatters and CLI drivers. 

| Layer / Package | Abstraction / Responsibility | Key Public API |
|-----------------|------------------------------|----------------|
| **CLI Drivers**<br>`src/bin/` | Command-line interfaces parsing arguments and driving the pipeline. `horizon_wga` for pairwise, `horizon_guide` for progressive multi-species WGA, and `horizon_merge` for merging pairwise outputs. | N/A (Executable entry points) |
| **I/O & Formats**<br>`src/formats.rs`<br>`src/hal.rs` | Serialization and deserialization of alignments. Emits MAF v1 and GFA v1.1. Parses FASTA sequences. `hal.rs` decodes Hierarchical Alignment (HAL) via HDF5. | `format_maf_blocks`, `format_gfa_blocks`, `read_fasta_contigs` |
| **Multi-species Merger**<br>`src/merge.rs` | `O(N log N)` sweep-line exact-coordinate multi-species MAF merger and query-insertion gap-column threading. Relies on the `egraph` tier for deduplication. | `merge_pairwise_mafs` |
| **E-Graph Resolution**<br>`src/egraph.rs` | Congruence-closed e-graph for identifying and collapsing tandem duplications and nested inversions in synteny blocks. | `resolve_lcb_duplications`, `resolve_alignment_duplications` |
| **Alignment (MDL+BOCPD)**<br>`src/mdl/` | Dynamic programming gap alignment using Minimum Description Length (MDL) and Bayesian Online Changepoint Detection (BOCPD) for alignment evaluation and X-drop extensions. | `align_contig_pair`, `banded_mdl_gap`, `mdl_xdrop_extend` |
| **Seeding & Chaining**<br>`src/seeding/` | Multi-shape spaced-seed hashing (KP-14, W9, LASTZ, RY20), Longest Increasing Subsequence (LIS) chaining, and colinear block (LCB) clustering. | `canonical_spaced_seeds`, `lis_chain`, `local_w9_subchain` |
| **BOCPD Core**<br>`src/bocpd.rs` | Core Bayesian Online Changepoint Detection for base frequencies (GC shifts). Base layer utility used by MDL scoring. | `bocpd_gc_profile` |

## Cross-Cutting Contracts

- **Dependency Direction:** Code in `src/mdl/` imports `src/seeding/` and `src/bocpd.rs`, but not vice versa. `src/merge.rs` and `src/formats.rs` consume the lower layers but do not depend on each other's internals (except formats using egraph tags).
- **Coordinate Systems:** The pipeline strictly uses 0-based, half-open intervals `[start, end)`. Reverse-strand coordinates must be relative to the forward strand's length (i.e. typical UCSC/MAF convention) unless explicitly translated for DP.
- **Source of Truth for Alignments:** The `ContigAlignment` struct in `mdl/mod.rs` is the definitive representation of a completed pairwise alignment. Formatters (`formats.rs`) project it to MAF or GFA without modifying the structural topology.

## Known Asymmetries

- `hal.rs` is primarily a reader/decoder for HDF5 HAL files and heavily relies on bit-twiddling and chunk unzipping. It currently lacks full structural emission equivalent to what `formats.rs` provides for MAF/GFA. 
- While `seeding` contains complex logic for clustering and chaining (`local_w9_subchain`, `lis_chain`), the resolution of paralogs and structural duplications defers to `egraph.rs`, meaning overlapping blocks survive the seeding phase and are untangled downstream.