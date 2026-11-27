# Architecture Audit — `horizon-aligner`

- **Audit Date:** 2026-11-25
- **Initial Audited Ref:** `0ab93be` → **All Proposals Implemented (`A1`–`A8`):** `f9b938b` (A3–A5), `a4d375e` (A1, A2, A6), and Commit 3 (A7, A8, docs)
- **Baseline:** [`ARCHITECTURE.md`](../ARCHITECTURE.md)

## Headline Verdict

**The suite architecture is clean, strictly acyclic, and mechanically verified.** Every source file under `src/` (`82` `.rs` files) stays strictly under `1,000` LOC (max `969`), `#![forbid(unsafe_code)]` holds crate-wide, external dependencies stay minimal (`anyhow`, `flate2`, `rayon`, `thiserror` — all four actively used), dead public surface has been purged, and all 8 structural/abstraction proposals (`A1`–`A8`) identified in the audit have been implemented and enforced via `scripts/check_docs.py`.

---

## Findings Summary

| ID | Dimension | Severity | Status | Summary |
| :--- | :--- | :---: | :---: | :--- |
| **A1** | Boundary integrity | Medium | `[fixed]` | Broke the `formats <-> tba` cycle by moving `HnexArchive::{from_tba_blocks, to_tba_blocks}` into `src/tba/hnex.rs` (`a4d375e`). |
| **A2** | Boundary integrity | Medium | `[fixed]` | Moved aligner-driven rGFA construction out of `src/gfa/` to top-level `src/rgfa.rs`, restoring `src/gfa/mod.rs` as a leaf I/O parser (`a4d375e`). |
| **A3** | Boundary integrity | Low | `[fixed]` | Extracted dependency-free `src/dna.rs` (`base_to_idx`, `complement_base`, `revcomp`), deleted duplicate `base_to_idx` in `horizon_guide/consensus.rs`, and moved `seeding/translated.rs`'s `align_contig_pair_runs` test to `src/lib.rs` (`f9b938b`). |
| **A4** | Abstraction consistency | Medium | `[fixed]` | Added typed `coord::Strand` enum and UCSC helpers (`from_ucsc`, `ucsc_to_fwd_interval`, `ucsc_offset_to_fwd_pos`, `flip_ucsc_start`, `fwd_pos`, `revcomp`) to `coord::Interval`; delegated `MafRow` and `TbaRow` coordinate math to `coord` (`f9b938b`). |
| **A5** | Abstraction consistency | Medium | `[fixed]` | Promoted half-open `coord::IntervalSet` into `src/coord.rs` and replaced all 3 duplicate interval-set implementations across `netting.rs`, `seeding/chain.rs` (eliminating the closed/half-open name collision), `horizon_guide/consensus.rs`, and `horizon_locus.rs` (`f9b938b`). |
| **A6** | Duplication / missed reuse | Medium | `[fixed]` | Promoted `AtomicTempFile` to `src/io_util.rs`; added `IndexedGenome` + `effective_contig_count` in `src/seeding/sketch.rs`; added shared `PairSchedule` + `parse_species_fasta_spec` in `src/lib.rs` adopted by `horizon_guide` and `horizon_tba` (`a4d375e`). |
| **A7** | Contract enforcement | Medium | `[fixed]` | Added `check_module_boundaries()` in `scripts/check_docs.py` enforcing a strict 9-tier unidirectional DAG across all `16` library modules (0 cycles, 0 upward imports). |
| **A8** | Duplication / missed reuse | Low | `[fixed]` | Added shared `tests/common/mod.rs` (`TestTempDir`, `synthetic_dna`, `mutate_dna`) adopted across `tests/{suite_cli,tba_cli,wga_paf}_e2e.rs`. |
| **A9** | Abstraction consistency | Low | `[accepted]` | `paf_filter::BinnedIntervals` (depth-`k` multiplicity) vs `coord::IntervalSet` (boolean occupancy) serve distinct invariants and remain separate. |
| **A10** | Abstraction consistency | Low | `[accepted]` | `horizon_resolve::snarl::PackedGraph` vs `horizon_clean::engine` flank-indexed state remain separate to avoid bloating `horizon_clean`'s streaming memory footprint; both share `src/gfa/mod.rs` helpers. |
| **A11** | Naming / vocabulary | Low | `[accepted]` | `horizon_resolve/resolve.rs:358::parse_path_sample_fields` is a 3-line owning adapter delegating directly to `horizon_aligner::gfa::parse_path_sample_fields`. |
| **A12** | Dead / reserved surface | Low | `[accepted]` | Library-only primitives (`cdbg.rs`, `tantan.rs`, `mdl/train.rs`, `mdl/squint.rs`, `PowerSumSketch`) are intentionally exported without dedicated CLIs (documented in `ADR-005` & `ADR-008`). |
| **A13** | Doc / code agreement | Low | `[fixed]` | Updated `manuscript/guide_egraph_synteny.md` (+ rebuilt PDF) to fix invalid `--gfa` flag syntax on `horizon_resolve`, `src/formats/dup_resolve.rs` module path, leaf-scaled `MdlCostModel` framing note, and default divergence heuristic (`51a7f81`). |
| **A14** | Doc / code agreement | Low | `[fixed]` | Upgraded `ARCHITECTURE.md` with the enforced strict DAG diagram and ADR reference index, and hardened `scripts/check_docs.py` against passing `--gfa` to `horizon_resolve` and stale `rgfa` paths (`51a7f81` + Commit 3). |

---

## Resolution Notes (`A1`–`A8`)

1. **Strict Module DAG (`A1`, `A2`, `A3`, `A7`):** All 16 library modules now sit in strictly ordered tiers (`dna`/`coord`/`io_util`/`egraph` = 0 → `bocpd`/`gfa` = 1 → `seeding` = 2 → `mdl` = 3 → `netting` = 4 → `hal` = 5 → `formats` = 6 → `paf` = 7 → `paf_filter`/`merge`/`rgfa`/`tba` = 8), verified on every `python3 scripts/check_docs.py` run.
2. **Unified Base Abstractions (`A4`, `A5`, `A6`):**
   - `src/dna.rs` is the single dependency-free source of truth for ASCII base lookup and complement/revcomp.
   - `src/coord.rs` is the single source of truth for `Strand`, UCSC `Interval` arithmetic, and half-open `IntervalSet`.
   - `src/io_util.rs::AtomicTempFile`, `src/seeding/sketch.rs::{IndexedGenome, effective_contig_count}`, and `src/lib.rs::{PairSchedule, parse_species_fasta_spec}` eliminate duplicated scaffolding across multi-species CLI binaries.
3. **Integration Test Hygiene (`A8`):** `tests/common/mod.rs` centralizes `TestTempDir` and deterministic `synthetic_dna`/`mutate_dna` across the three CLI E2E test suites.
