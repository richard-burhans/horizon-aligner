# Horizon Suite Architecture

This document establishes the conceptual baseline of the major abstractions in `horizon-aligner`, the strictly unidirectional module dependency DAG between them, and the contracts verified mechanically by `scripts/check_docs.py` and the [Architecture Audit](docs/architecture_audit.md). Every `.rs` file in `src/` (`107` files) is maintained strictly under `1,000` LOC (current max `942` LOC).

## Dependency-Direction Diagram (Enforced Strict DAG)

```text
  src/bin/ (13 CLI binaries: horizon + 12 subcommands)
      │
      ├──► Tier 8: src/tba/ · src/merge.rs · src/rgfa.rs · src/paf_filter/
      │       │        │          │              │
      │       │        │          │              ▼
      │       │        │          │      Tier 7: src/paf.rs
      │       │        │          │              │
      │       │        ▼          │              ▼
      │       ├──► Tier 6: src/formats/ ◄────────┘
      │       │       │
      │       ├──► Tier 5: src/hal/
      │       │       │
      │       └──► Tier 4: src/netting.rs
      │               │
      │               ▼
      │            Tier 3: src/mdl/
      │               │
      │               ▼
      │            Tier 2: src/seeding/
      │               │
      │               ▼
      │            Tier 1: src/bocpd.rs · src/gfa/mod.rs
      │               │                       │
      │               ▼                       ▼
      └──────────► Tier 0: src/dna.rs · src/coord.rs · src/io_util.rs · src/egraph.rs
```

## Tier Map

| Layer / Package | Abstraction / Responsibility | Key Public API |
|-----------------|------------------------------|----------------|
| **CLI Drivers**<br>`src/bin/` | `13` suite executables with stdout/stdin hygiene (`-t`/`--threads` across all 9 multithreaded binaries `horizon_wga`, `horizon_guide`, `horizon_tba`, `horizon_merge`, `horizon_locus`, `horizon_filter`, `horizon_eval`, `horizon_clean`, `horizon_resolve` + `horizon` forwarding; the 3 streaming I/O binaries `horizon_gfa`, `horizon_hal`, `horizon_lift` are single-threaded): `horizon` (subcommand dispatcher), `horizon_wga` (`main.rs`, `graph.rs`: pairwise, `--maf`, `--gfa`, `--paf`, native UCSC `--chain`/`--net`, `--one-to-one`, `--tantan`, `--train-hmm`, `--graph` DAG/POA, `--hal-info`), `horizon_guide` (`main.rs`, `store.rs`, `rescue.rs`, `sister_align.rs`, `consensus.rs`, `hal_export.rs`: adaptive `>512 MiB` streaming `SpeciesStore`, `auN`-weighted MinHash hinge, BFS UPGMA rescue waves, single-pass MAF-bucketed `16 Mb` consensus, streaming HAL), `horizon_tba` (`main.rs`, `cli.rs`, `pipeline.rs`, `stream_extend.rs`: all-vs-all, constant-RAM streaming two-tier `--hinges`/`--auto-hinges`/`--extend-maf`, bottom-up `--progressive-ancestors` Fitch/ML + `--hal`, `--variation-gfa`, `--out-gz`, `--tantan`, `--pairwise-gz`, `--from-mafs`), `horizon_merge` (streaming multi-species MAF/rGFA merge), `horizon_gfa` (streaming `stats`/`qc`/`list-walks`/`extract`/`from-maf` + incremental `rgfa --respect-softmask`), `horizon_hal` (pure-Rust zero-`libhdf5` streaming multi-level HAL `stats [--chroms]`/`liftover`/`to-maf`/`to-gfa`), `horizon_lift` (`main.rs`, `lift.rs`, `chain.rs`, `features.rs`: GFA walk & UCSC `.chain` BED/GFF3 liftover), `horizon_clean` (`main.rs`, `engine.rs`: symmetric BOCPD+MDL junction cleaner, `--min-haplotypes`, `--clip-bubbles-bp`), `horizon_resolve` (`main.rs`, `snarl.rs`, `resolve.rs`, `vcf.rs`, `vcf_nested.rs`: `PackedGraph` Tarjan bridges, `k`-copy e-graph resolver, phased VCFv4.2 export + `--nested-snarls` `LV=1` local bubble decomposition), `horizon_filter` (thin CLI over `src/paf_filter/` with `--net-score mdl\|log-length-ani`), `horizon_locus` (regional multi-FASTA locus aligner, BioNJ, `--hnex`, `--consistency`, `--repair-frameshifts`), `horizon_eval` ($O(R \log R)$ `DiagRun` MAF evaluator). | Executable entry points, `PairSchedule`, `parse_species_fasta_spec` |
| **Multi-Species Pipelines, rGFA & PAF Filter**<br>`src/tba/`<br>`src/merge.rs`<br>`src/rgfa.rs`<br>`src/paf_filter/` | Reference-free Threaded Blockset Aligner (`tba/{mod,model,merge,poa_refine,hnex,tree,two_tier,variation_gfa,writer}.rs` and `tba/progressive/{mod,ancestral,hal_build}.rs` for bottom-up Fitch/ML ancestral reconstruction with non-root child-union retention and hierarchical HAL construction, plus streaming tip-by-tip projection and packed-RAM `24 B`/`16 B` native base-level sequence-bearing variation GFA v1.1 emission), reference-guided $O(N \log N)$ sweep-line exact-coordinate multi-species MAF/rGFA merger (`merge.rs`), Minigraph-style incremental rGFA SV-grafting builder (`rgfa.rs`), and SDUST+MDL PAF scoring + nested `CrcCache` + `max_occ` interval & symmetric reciprocal-best `log-length-ani` netting (`paf_filter/{mod,net}.rs`). | `merge_pairwise_mafs`, `merge_pairwise_maf_paths_to_writers`, `run_tba_tree_merge`, `tba_project_tips_streaming`, `run_progressive_ancestral_tba`, `refine_shared_tip_insertions_poa`, `parse_newick`, `build_upgma_tree`, `build_bionj_tree`, `write_tba_maf_and_gfa`, `write_tba_maf_and_variation_gfa`, `write_variation_gfa_from_maf_reader`, `RgfaGraph`, `score_paf_record_with_floor`, `net_top_k_records_with_mode` |
| **Formats, PAF & HAL Serialization**<br>`src/paf.rs`<br>`src/formats/`<br>`src/hal/` | `paf.rs` streams PAF + CIGAR (`cg:Z:`). `formats/{mod,maf,dup_resolve,hnex,realign,ucsc,ucsc_net}.rs` emits MAF v1, GFA v1.1, native UCSC `.chain` and hierarchical `.net` (`top`/`syn`/`inv`/`dup`/`nonSyn`), lossless `.hnex` archives (with `TbaBlock` adapters in `tba/hnex.rs`), overlap-clustered empirical MDL paralog pruning, and zero-alloc stack HOXD70 affine micro-gap realignment (`MAX_COLINEAR_GAP_BP = 40`, neutral `N`). `hal/` provides `HalAlignment` (`mod.rs`), fail-closed `Read+Seek` HDF5 `v0`/`v1` container & dataset readers (`hdf5_read.rs`, `hdf5_dataset.rs`), $O(\log C)$ chunk B-tree range-descent streaming dataset reader (`hdf5_stream.rs`, removing the `256 MiB` dataset cap), arbitrary-depth LCA multi-level streaming HAL traversal & MAF/GFA projection (`stream_hal.rs`, `stream_hal_genome.rs`, `stream_hal_maf.rs`), nibble-packed Fitch parsimony with `N`-scaffold preservation (`fitch.rs`), flat-array paralogy normalization (`normalize.rs`), HDF5 `v0` builder primitives (`hdf5_writer.rs`), unified in-memory exporter (`hdf5_write.rs`), and `4 MiB`-windowed streaming Fitch progressive exporter (`hdf5_progressive.rs`). | `stream_paf_records`, `format_paf_records_runs`, `format_maf_blocks_runs`, `format_gfa_blocks_runs`, `build_ucsc_chains`, `build_ucsc_nets`, `parse_ucsc_chains`, `HnexArchive`, `resolve_alignment_duplications_mdl`, `prune_resolved_dup_clusters`, `read_fasta_contigs`, `HalAlignment`, `StreamingHalFile`, `write_multi_contig_progressive_hdf5_hal` |
| **Netting**<br>`src/netting.rs` | Species-wide $O(R \log R)$ global 1-to-1 orthology netting built over `coord::IntervalSet` + shared canonical `score_acgt_columns`. | `net_one_to_one_runs_scored`, `net_global_one_to_one_runs_scored`, `count_unique_covered_bp` |
| **Alignment (MDL+BOCPD, DAG & Consistency)**<br>`src/mdl/` | Shared DP parameter & `BandLayout` types (`dp_params.rs`), pairwise banded Pair-HMM Forward-Backward MEA DP + local-BOCPD flank trimming (`dp.rs`), composition-aware diagonal & narrow-Viterbi fast paths with per-block indel gating (`dp_fast.rs`), banded X-drop (`xdrop.rs`), run-length encoding (`runs.rs`), uncovered-gap fallback (`fallback.rs`), biological + GC-matched Conformal Risk Control (`crc.rs`), adaptive `HmmParams` shrinkage (`mod.rs`), seeded topological DAG/POA alignment (`dag.rs` capped at `MAX_DAG_DP_CELLS = 32_000_000`, `dag_seed.rs` + `dag_seed_tests.rs` with `MASK_W10`, `MASK_RY18`, 6-frame `aa5_keys`, and strand-aware `estimate_dag_hmm_params`), rolling-row banded ProbCons probabilistic consistency MSA with 2-bit packed MEA DP (`consistency.rs`, `MAX_PAIR_HMM_CELLS = 16_000_000`), directional Baum–Welch Pair-HMM training (`train.rs`, `--train-hmm`), and fixed-share null / codon frameshift repair (`squint.rs`, `--repair-frameshifts`). | `align_contig_pair_runs`, `banded_mdl_gap`, `banded_mdl_gap_posterior_trace`, `mdl_xdrop_extend_with_budget`, `calibrate_null_crc_matched`, `align_seq_to_dag_mdl`, `align_seq_to_dag_seeded_indexed`, `estimate_dag_hmm_params`, `align_locus_probcons_msa`, `train_asymmetric_pair_params`, `repair_coding_frameshifts` |
| **Seeding & Chaining**<br>`src/seeding/` | SDUST repeat masking, `FlatSeedIndex` (`flat_index.rs`, `12 B/entry` + 16-bit prefix buckets), zero-copy owned `ContigIndex` + `IndexedGenome` + `effective_contig_count` + `adaptive_top_k_for_fragmentation` (`sketch.rs`), bottom-$k$ MinHash sketching & Broder `sketch_jaccard`, streaming sixth-root Eisenstein hashing, auto-hinge sketching, & power-sum set reconciliation (`eisenstein/{mod,sketch,reconcile}.rs`), constant-$E$ Karlin–Altschul HSP gating (`HSP_SCALE = 7.25`) & Rayon-parallel `LASTZ12` transition seeding (`hsp.rs`), 6-frame biochemical `AA5` (`translated.rs`), LIS chaining & `BTreeSet`-indexed LCB clustering (`chain.rs`), palindrome-aware colored compacted de Bruijn graphs (`cdbg.rs`), and 64 kb-chunked `tantan` HMM tandem-repeat masking (`tantan.rs`, `--tantan`). | `FlatSeedIndex`, `IndexedGenome`, `effective_contig_count`, `adaptive_top_k_for_fragmentation`, `prepare_contig_index`, `prepare_contig_index_owned`, `sketch_jaccard`, `merge_contig_sketches`, `EisensteinSamplerSketch`, `select_auto_hinges_eisenstein`, `PowerSumSketch`, `build_colored_cdbg`, `tantan_mask`, `lis_chain`, `local_w9_subchain` |
| **Core Models & Streaming GFA I/O**<br>`src/bocpd.rs`<br>`src/gfa/mod.rs` | Exact Adams & MacKay `(1.0, 1.0)` Beta-Bernoulli BOCPD compositional segmentation (`bocpd_gc_profile`) + finite-guarded `joint_pair_bits` (`bocpd.rs`); leaf streaming GFA v1.0/v1.1 reader, `SpoolableGfaInput`, and shared CIGAR/PanSN helpers (`gfa/mod.rs`, importing only `io_util`). | `bocpd_gc_profile`, `stream_gfa_record_refs`, `parse_cigar_overlap_bp`, `parse_path_sample_fields` |
| **Dependency-Free Base Primitives**<br>`src/dna.rs`<br>`src/coord.rs`<br>`src/io_util.rs`<br>`src/egraph.rs` | Alphabet & Watson–Crick complement primitives (`dna.rs`: `base_to_idx`, `complement_base`, `revcomp`); typed `Strand`, UCSC `Interval`, and unified half-open `IntervalSet` (`coord.rs`); transparent plain / `.gz` / `.bgz` / `.zst` / `stdin` / `stdout` streaming, `read_bounded_line`, `AtomicTempFile`, and `ProgressLogger` (`io_util.rs`); deterministic (`BTreeMap`/`BTreeSet`) congruence-closed `k`-copy e-graph & leaf-scaled `MdlCostModel` (`egraph.rs`). | `base_to_idx`, `complement_base`, `revcomp`, `Strand`, `Interval`, `IntervalSet`, `open_text_reader`, `open_text_writer`, `read_bounded_line`, `AtomicTempFile`, `ProgressLogger`, `EGraph`, `MdlCostModel` |

## Cross-Cutting Contracts

- **Memory Safety:** `#![forbid(unsafe_code)]` is enforced crate-wide in `Cargo.toml`; zero `unsafe` blocks exist anywhere in `horizon-aligner`.
- **Module Size Budget:** Every `.rs` file in `src/` must remain strictly `< 1,000` LOC (enforced mechanically by `scripts/check_docs.py`).
- **Strict Unidirectional Dependency DAG:** Every production import `crate::<dst>` across all `16` library modules must flow strictly downward across `MODULE_TIERS` (`Tier 8 -> ... -> Tier 0`), with zero cycles and zero upward imports (enforced mechanically by `scripts/check_docs.py::check_module_boundaries`).
- **Coordinate Systems:** All internal intervals and `coord::IntervalSet` ranges use 0-based, half-open `[start, end)` coordinates. Reverse-strand UCSC translation (`Interval::from_ucsc`, `ucsc_to_fwd_interval`, `ucsc_offset_to_fwd_pos`, `flip_ucsc_start`) and typed `Strand` live in `src/coord.rs`.
- **Source of Truth for Alignments:** `RunAlignment` (`AlignedRun` coordinate runs in `mdl/runs.rs`) is the compact pairwise representation consumed directly by `netting.rs`, `formats/`, `paf.rs`, `merge.rs`, `rgfa.rs`, `tba/`, `horizon_eval`, and `horizon_guide` without per-base pair expansion.
- **Fail-Closed Binary & Text Readers:** HDF5 B-tree walkers (`hdf5_read.rs`, `hdf5_dataset.rs`) enforce strict metadata byte/entry/depth/cycle/level budgets and fail closed (`Ok(None)`) on any truncation or structural anomaly; text parsers (`hnex.rs`, `rgfa.rs`, `two_tier.rs`, `tree.rs`) enforce explicit bounded line/recursion limits via `read_bounded_line` and `MAX_NEWICK_DEPTH`.

## Known Asymmetries

- `horizon_clean` preserves `S`-line optional tags, severs low-support discordant edges (`walk_support <= 1` by default, or `0` with `--prune-unsupported`) into `HB:i:1` Horizon Bubbles using an orientation-symmetric Bernoulli/BIC + BOCPD junction test, and supports `--min-haplotypes` (PanSN `W`/`P` parsing) and `--clip-bubbles-bp` pruning.
- `horizon_resolve` decomposes graphs via `u32`-interned `PackedGraph` + $O(V + E)$ iterative Tarjan bridges, estimates window divergence via a topological `[0.02, 0.15]` heuristic (unless `--divergence` is passed), streams `S`-line sequences intact on rewrite, rewrites `P`/`W` walks with `k`-copy repeat saturation, updates `W`-line `end` bounds, prunes/synthesizes `L` links (assuming `0M` overlaps for walk rewriting), and projects Tarjan snarls to phased VCFv4.2 (`--vcf`, trimming closing link overlaps and marking multi-visit conflicts as `.`).
- `horizon_merge` emits flat rank-0 rGFA (`SR:i:0`), whereas `horizon_gfa rgfa` emits multi-rank incremental SV rGFA (`SR:i:0..R`).
- `horizon` dispatcher (`src/bin/horizon.rs`) forwards `wga`, `merge`, `guide`, `tba`, `gfa`, `hal`, `lift`, `clean`, `resolve`, `filter`, `locus`, `eval`, and `help <subcommand>`.
- `hal/hdf5_write.rs` and `hal/hdf5_progressive.rs` emit `libhdf5`-compatible HDF5 `v0` archives with multi-level B-trees (`>256` entries), whereas `hal/hdf5_read.rs` supports both `v0` and `v1` object headers.

## Reference Index

| Abstraction | Defining File(s) | Decision Record |
| :--- | :--- | :--- |
| Reference-free TBA & bounded-RAM scheduling | `src/tba/{mod,model,merge,writer}.rs`, `src/bin/horizon_tba/` | [ADR-001](docs/decisions/0001-standalone-horizon-tba.md) |
| Streaming `4 MiB`-windowed progressive HAL export | `src/hal/{hdf5_progressive,fitch,normalize}.rs` | [ADR-002](docs/decisions/0002-streaming-progressive-hal-export.md) |
| HOXD70 micro-gap realignment, `auN` hinge, global 1-to-1 netting | `src/formats/realign.rs`, `src/seeding/sketch.rs`, `src/netting.rs` | [ADR-003](docs/decisions/0003-hoxd70-micro-gap-realignment.md) |
| Bottom-up progressive ancestral TBA, ML transducers, POA | `src/tba/progressive/`, `src/tba/poa_refine.rs` | [ADR-004](docs/decisions/0004-progressive-ancestral-tba-and-ml-transducers.md) |
| Eisenstein hashing, `--auto-hinges`, two-tier TBA, `--pairwise-gz` | `src/seeding/eisenstein/`, `src/tba/two_tier.rs` | [ADR-005](docs/decisions/0005-eisenstein-hashing-auto-hinges-and-two-tier-tba.md) |
| `k`-copy e-graph saturation & empirical MDL dup pruning | `src/egraph.rs`, `src/formats/dup_resolve.rs` | [ADR-006](docs/decisions/0006-k-copy-egraph-repeat-saturation-and-mdl-dup-pruning.md) |
| Incremental rGFA, Tarjan-snarl VCFv4.2, lossless `.hnex` | `src/rgfa.rs`, `src/bin/horizon_resolve/vcf.rs`, `src/formats/hnex.rs`, `src/tba/hnex.rs` | [ADR-007](docs/decisions/0007-pangenome-interoperability-rgfa-snarl-vcf-and-hnex.md) |
| Locus ProbCons consistency, BioNJ, library-only primitives | `src/mdl/{consistency,squint,train}.rs`, `src/seeding/{cdbg,tantan}.rs`, `src/tba/tree.rs` | [ADR-008](docs/decisions/0008-locus-probcons-consistency-bionj-and-library-primitives.md) |

## Architecture Audit

See [docs/architecture_audit.md](docs/architecture_audit.md) for the 7-dimension audit report and verification status (`A1`–`A14`).
