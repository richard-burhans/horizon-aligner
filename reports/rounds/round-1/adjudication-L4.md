# Round 1 Adjudication — L4 (Provenance)

**Target attack:** `reports/rounds/round-1/attack-L4-provenance.md` (SHA-256 `d8c4fb1194dcd7f71579eaeb6902ee377256a317f65fe60cab44d72565ac1e91`, matches sealed digest `d8c4fb1194dc`).
**Prior rounds & method cards:** `./scripts/round history "STRATEGIC_PLAN" --kind adjudication` confirms 0 closed rounds. No `refs/methods/` cards exist in the repository, and no statistical dispute is raised here.

---

## F1 -- C1: Tool inventory and Unix-pipe composition
- **Diagnosis:** **NARROWED.** The attacker's codebase check is verified: `src/bin/` contains `bench_chr1.rs`, `horizon_eval.rs`, `horizon_guide.rs`, `horizon_merge.rs`, and `horizon_wga.rs`; `horizon_clean`, `horizon_resolve`, and `horizon_lift` do not exist in `src/` or `Cargo.toml`, and `stdin` has 0 matches in `src/`. However, C1's first sentence (*"We are not building a monolithic 'run-all' executable"*) is true of the 5-binary layout, and *"We are building..."* is aspirational rather than a claim that all four tools are already shipped.
- **What survives:** The named-tool inventory (`horizon_clean`, `horizon_resolve`, `horizon_lift`) and present-tense Unix-pipe composition claim in `docs/STRATEGIC_PLAN.md:6` lack codebase provenance and omit existing binaries.
- **Correction:** **ENDORSED** (also note `bench_chr1` in `src/bin/` for completeness). **The fix is free** (text-only edit to `docs/STRATEGIC_PLAN.md:6`).

## F2 -- C2: MDL vs. heuristic thresholds and pangenome hairball provenance
- **Diagnosis:** **UPHELD.** Verified in `src/mdl/mod.rs` (`max_gap_bp = 2400` at line 123, `hsp_threshold` at line 47, clamped `min_fallback_cov` at lines 187–192) and `src/mdl/dp.rs:81` (`min_block_bits`), as well as `docs/BENCHMARK_RESULTS.md:39` (`HPRC (v1.1)` unrun as `*(Graph Eval)*`).
- **What survives:** MDL replaces substitution/gap scoring in pairwise DP (`src/mdl/dp.rs`), while heuristic thresholds remain in `src/mdl/mod.rs` and pangenome hairball prevention is untested.
- **Correction:** **ENDORSED.** **The fix is free** (text-only edit to `docs/STRATEGIC_PLAN.md:9`).

## F3 -- C3: E-graph HLA and pangenome graph provenance
- **Diagnosis:** **UPHELD.** Verified `src/egraph.rs:1-385` is invoked on pairwise/merged LCBs (`src/bin/horizon_wga.rs:254`, `src/merge.rs`), with zero HLA or pangenome-graph tests or benchmarks in `src/`, `docs/BENCHMARK_DATASETS.md`, or `docs/BENCHMARK_RESULTS.md:39`.
- **What survives:** `src/egraph.rs` resolves tandem duplications and inversions on pairwise/merged LCBs; HLA and cyclical pangenome-graph claims remain unverified.
- **Correction:** **ENDORSED.** **The fix is free** (text-only edit to `docs/STRATEGIC_PLAN.md:10`).

## F4 -- C4: Uncited `vg`/`odgi` crash claim and `src/gfa/mod.rs` zero-copy/dead-code status
- **Diagnosis:** **UPHELD.** Verified: (1) `docs/` and `manuscript/references.bib` contain zero citations or benchmarks for `vg`/`odgi` crashing on 12GB pangenomes; (2) `src/gfa/mod.rs:75-120` uses `reader.read_line(&mut line)` and `.to_string()`/`.to_vec()` allocations for every record field (streaming, not zero-copy); (3) `stream_gfa` has 0 callers in `src/` and 0 unit tests.
- **What survives:** C4's uncited 12GB crash claim, "zero-copy" characterization, and "solves the I/O bottleneck entirely" claim are all refuted.
- **Correction:** **ENDORSED.** **The fix is free** (text-only edit to `docs/STRATEGIC_PLAN.md:11`).

## F5 -- C5: `stdin`/`stdout` GFA/MAF/VCF communication and undisclosed `src/hal/`
- **Diagnosis:** **NARROWED.** Verified: `vcf` and `stdin` have 0 matches in `src/`, and `src/hal/` (`mod.rs`, `hdf5_read.rs`, `hdf5_write.rs`) implements a pure-Rust HDF5 `.hal` reader/writer exposed via `horizon_wga --hal-info`. However, C5 is written in future tense (*"We will never invent a proprietary binary format. Our tools will communicate..."*), and `.hal` is an existing external format (Hickey et al. 2013) rather than a newly invented Horizon format, so the future-tense commitment is not itself refuted by current file-based I/O.
- **What survives:** `docs/STRATEGIC_PLAN.md:15-17` omits `horizon-aligner`'s existing `src/hal/` reader/writer and leaves unstated that VCF and `stdin`/`stdout` streaming are not yet implemented.
- **Correction:** **ENDORSED** as a current-state disclosure (file-based FASTA/MAF/GFA I/O plus `src/hal/`). **The fix is free** (text-only edit to `docs/STRATEGIC_PLAN.md:15-17`).

## F6 -- C6: Normative `stdin`/`stdout` requirement
- **Diagnosis:** **NARROWED.** Verified that 0 of 5 binaries in `src/bin/` read from `stdin`. However, C6 (*"Every tool must read GFA/MAF from `stdin` and write to `stdout`"*) is a normative design rule under `## 3. Strategic Adjustments & Risk Mitigations`, not an assertion that current binaries already do so. The attacker should have distinguished normative future requirements from present-tense implementation claims.
- **What survives:** No current binary implements `stdin` GFA/MAF input, so C6 requires an explicit unimplemented/future-requirement label.
- **Correction:** **ENDORSED.** **The fix is free** (text-only edit to `docs/STRATEGIC_PLAN.md:23`).

## F7 -- C7: Linear GFA-stream-to-E-graph bubble stitching
- **Diagnosis:** **NARROWED.** Verified that `stream_gfa` has 0 callers and `src/egraph.rs` operates only on `ContigAlignment`/`MafBlock` with no GFA bubble extraction or path stitching. However, C7 is explicitly future-tense (*"We will stream the pangenome linearly..."*) under `## 3. Strategic Adjustments & Risk Mitigations`, not a claim of existing code.
- **What survives:** GFA bubble extraction and path stitching have no implementation in `src/` and must be labeled as a future design proposal.
- **Correction:** **ENDORSED.** **The fix is free** (text-only edit to `docs/STRATEGIC_PLAN.md:27`).

## F8 -- C8: Uncited HPRC crash premise and non-existent release targets
- **Diagnosis:** **UPHELD.** Verified: the present-tense premise (*"The HPRC produces massive graphs that researchers cannot use because standard tooling crashes or produces tangled paths"*) is uncited across `docs/` and `manuscript/references.bib`; `horizon_clean`, `horizon_resolve`, and `horizon_lift` do not exist; and `docs/BENCHMARK_RESULTS.md:23-39` contains 7 pairwise/MSA benchmarks and 0 HPRC graph results.
- **What survives:** Both the uncited HPRC tooling-crash premise and the claimed "first" benchmarking/release focus on three non-existent tools lack provenance.
- **Correction:** **ENDORSED.** **The fix is free** (text-only edit to `docs/STRATEGIC_PLAN.md:31`).

## F9 -- C9: `horizon_clean` roadmap item and `GfaLink` sequence limitation
- **Diagnosis:** **NARROWED.** Verified: `horizon_clean` does not exist, `stream_gfa` has 0 callers, and `GfaLink` (`src/gfa/mod.rs:21-27`) stores only segment IDs, orientations, and `overlap` — making stateless `L`-line BOCPD/MDL sequence scoring impossible without retaining `S` sequences. However, `docs/PANGENOME_SUITE.md:34` explicitly frames the document as a *strategic roadmap*, so non-existence of `horizon_clean` alone does not refute a roadmap item; the technical `GfaLink` limitation does.
- **What survives:** `horizon_clean` is an unimplemented roadmap proposal, and its stateless `L`-line BOCPD/MDL scoring design is technically infeasible without `S`-sequence lookup.
- **Correction:** **ENDORSED.** **The fix is free** (text-only edit to `docs/PANGENOME_SUITE.md:5-10`).

## F10 -- C10: `horizon_wga --graph` roadmap item
- **Diagnosis:** **NARROWED.** Verified: `src/bin/horizon_wga.rs:32-145` has no `--graph` flag and `src/mdl/dp.rs:81` (`banded_mdl_gap`) accepts only linear `&[u8]` slices. However, `docs/PANGENOME_SUITE.md:34` explicitly frames the document as a *strategic roadmap* rather than a current-release CLI reference.
- **What survives:** `horizon_wga --graph` and graph-bifurcation MDL DP do not exist in `src/` and must be clearly marked as unimplemented proposals so readers do not mistake them for flags on the existing `horizon_wga` binary.
- **Correction:** **ENDORSED.** **The fix is free** (text-only edit to `docs/PANGENOME_SUITE.md:12-17`).

## F11 -- C11: `horizon_lift` roadmap item and dead `coord::Interval`
- **Diagnosis:** **NARROWED.** Verified: `horizon_lift` does not exist, and `src/coord.rs:5-38` (`Interval`) is not declared in `src/lib.rs:3-10` (`pub mod coord;` is absent), making `coord::Interval` uncompiled dead code. As with F9/F10, `docs/PANGENOME_SUITE.md:34` frames `horizon_lift` as a roadmap item, so non-existence of `horizon_lift` is expected, while citing `coord::Interval` as an active Horizon primitive is false.
- **What survives:** `horizon_lift` is an unimplemented roadmap proposal and `coord::Interval` is currently dead, uncompiled code omitted from `src/lib.rs`.
- **Correction:** **ENDORSED.** **The fix is free** (text-only edit to `docs/PANGENOME_SUITE.md:26-31`).

---

## Verdict Summary

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | 3/4 named tools and `stdin` piping are absent, though non-monolithic layout holds and "we are building" is aspirational; free fix endorsed |
| F2 | UPHELD | `src/mdl/mod.rs` retains heuristic thresholds and `HPRC (v1.1)` hairball prevention is unrun (`docs/BENCHMARK_RESULTS.md:39`); free fix endorsed |
| F3 | UPHELD | `src/egraph.rs` runs only on pairwise/merged LCBs with zero HLA or pangenome-graph provenance; free fix endorsed |
| F4 | UPHELD | `vg`/`odgi` 12GB crash is uncited, `src/gfa/mod.rs` allocates `String`/`Vec<u8>` per field (not zero-copy), and `stream_gfa` is dead code; free fix endorsed |
| F5 | NARROWED | VCF and `stdin` are absent and `src/hal/` is undisclosed, but C5 is a future-tense mitigation ("We will...") rather than a refuted present fact; free fix endorsed |
| F6 | NARROWED | 0 binaries read `stdin`, but C6 is a normative design rule ("must") in a mitigation section rather than a current-state claim; free fix endorsed |
| F7 | NARROWED | GFA-to-E-graph bubble stitching does not exist in `src/`, but C7 is an explicit future-tense mitigation ("We will..."); free fix endorsed |
| F8 | UPHELD | HPRC tooling-crash premise is uncited, all three release targets are absent, and `docs/BENCHMARK_RESULTS.md` has 0 HPRC results; free fix endorsed |
| F9 | NARROWED | `docs/PANGENOME_SUITE.md:34` is a roadmap, but `horizon_clean` is unimplemented and `GfaLink` lacks sequences needed for stateless `L`-line BOCPD/MDL; free fix endorsed |
| F10 | NARROWED | `docs/PANGENOME_SUITE.md:34` is a roadmap, but `--graph` and graph-bifurcation DP are absent on existing `horizon_wga` and need an unimplemented label; free fix endorsed |
| F11 | NARROWED | `docs/PANGENOME_SUITE.md:34` is a roadmap, but `horizon_lift` is unimplemented and `src/coord.rs` (`Interval`) is uncompiled dead code omitted from `src/lib.rs`; free fix endorsed |
