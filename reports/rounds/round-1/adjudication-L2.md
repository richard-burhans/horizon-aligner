# Round 1 Adjudication -- L2 (Sampling)

**Attack evaluated:** `reports/rounds/round-1/attack-L2-sampling.md` (sealed digest `cdcc127c9fdb`, verified via `sha256sum`).
**Prior rounds / method cards:** `scripts/round history "" --kind adjudication` confirmed no prior closed rounds exist. No `refs/methods/` cards exist in the repository; rulings below rely on direct verification of `docs/`, `scripts/benchmarks/`, and `src/`.

---

## F1 -- C8 (`docs/STRATEGIC_PLAN.md:31`): **NARROWED**

* **What survives:** Zero HPRC or pangenome-tool benchmarks exist (`docs/BENCHMARK_RESULTS.md:39` is an unrun `*(Graph Eval)* | - | -` placeholder), `horizon_clean`/`horizon_resolve`/`horizon_lift` do not exist in `src/bin/`, no `vg`/`odgi` failure samples are recorded, and actual benchmarking (`scripts/benchmarks/run_all_evals.py`) targets pairwise WGA first.
* **Diagnosis vs. attacker's scope:** **NARROWED on the strategic-intent leg; UPHELD on the empirical/status leg.** The attacker is right on every factual absence (and on the `v1.1` vs `v1.0` HPRC version mismatch between `docs/BENCHMARK_DATASETS.md:46` and `scripts/benchmarks/fetch_benchmarks.sh:27`), but `REFUTED` overreaches by treating a forward-looking strategic-plan sentence (*"We will focus our benchmarking and release strategy..."* under `## 3. Strategic Adjustments & Risk Mitigations`) as a falsified empirical report. Note also that `data/benchmarks/` being absent from a clean checkout is expected given `docs/BENCHMARK_DATASETS.md:51` ("re-download ... after a sandbox reset").
* **Repair:** **ENDORSED** (clarify in `docs/STRATEGIC_PLAN.md:31` that all current benchmarks sample pairwise UCSC MAF slices via `horizon_wga`/`horizon_eval`, and that `horizon_clean`, `horizon_resolve`, `horizon_lift`, and HPRC graph evaluations are unimplemented/unbenchmarked). **The fix is free** (wording edit requiring no re-measurement).

---

## F2 -- C4 (`docs/STRATEGIC_PLAN.md:11`): **UPHELD**

* **What survives:** `src/gfa/mod.rs:50-124` defines `stream_gfa`, has zero callers or tests anywhere in `src/`, has never been run on a 12GB pangenome, and has no comparative `vg`/`odgi` RAM-exhaustion sample.
* **Diagnosis:** **UPHELD.** Unlike C8, C4 is stated in the present tense under `## 2. What We Are Doing Right` (*"solves the I/O bottleneck entirely"*) and cites an uncalled module (`src/gfa/mod.rs`) with zero benchmarked GFA runs.
* **Repair:** **ENDORSED** (remove or qualify the present-tense claim in `docs/STRATEGIC_PLAN.md:11` until `stream_gfa` is wired into a CLI binary and benchmarked on a multi-gigabyte GFA alongside `vg`/`odgi`). **The fix is free** as a documentation qualification requiring no re-measurement (actually validating the claim requires running the multi-GB benchmark).

---

## F3 -- C2 (`docs/STRATEGIC_PLAN.md:9`): **UPHELD**

* **What survives:** `scripts/benchmarks/extract_fasta_from_maf.py:28-45` (and `59-76`) filters input FASTAs to homologous MAF blocks (`if ref_row and tgt_row:`), stripping unaligned/repetitive inter-block regions and concatenating blocks from coordinate `0` before `horizon_eval` runs; `scripts/benchmarks/run_all_evals.py` only evaluates tiny `chrM` (~16.5 kb) and `chr22_KI270735v1_random` (~38 kb) slices; and `docs/BENCHMARK_RESULTS.md:23-39` reports no graph or hairball topology metrics.
* **Diagnosis:** **UPHELD.** The attacker appropriately chose `NARROWED` on C2 rather than `REFUTED`. Whole-set check: `docs/BENCHMARK_RESULTS.md:37` also lists `cactus241way (chr22 1MB)` (not driven by `run_all_evals.py`), where recall drops to `9.9%`, which further underscores that MDL hairball prevention is unverified at chromosome/graph scale.
* **Repair:** **ENDORSED** (narrow `docs/STRATEGIC_PLAN.md:9` to state that MDL bit-savings has only been evaluated on pairwise UCSC MAF slices—primarily pre-filtered homologous blocks from `chrM` and `chr22_KI270735v1_random`, plus `chr22 1MB`—and not on pangenome graphs or hairball metrics). **The fix is free** (wording edit requiring no re-measurement).

---

## F4 -- C3 (`docs/STRATEGIC_PLAN.md:10`): **UPHELD**

* **What survives:** Neither `chr6` (HLA) nor any centromeric region is sampled in `docs/BENCHMARK_DATASETS.md`, `scripts/benchmarks/run_all_evals.py`, or `docs/BENCHMARK_RESULTS.md`, and `src/` contains zero references to `hla` or `centromer`.
* **Diagnosis:** **UPHELD.** `src/egraph.rs` exists, so the attacker's `NARROWED` verdict on C3 (removing the unverified HLA/centromere empirical assertion rather than denying e-graph code exists) is properly scoped. Whole-set check: `docs/PANGENOME_SUITE.md:24` repeats the same unbenchmarked HLA/centromere claim, though `PANGENOME_SUITE.md:34` frames that document as a strategic roadmap.
* **Repair:** **ENDORSED** (narrow `docs/STRATEGIC_PLAN.md:10` to remove or qualify the unverified HLA paralog assertion until an HLA/`chr6` or centromere dataset is sampled and benchmarked, and ensure `docs/PANGENOME_SUITE.md:24` is clearly marked as an unvalidated target). **The fix is free** (wording edit requiring no re-measurement).

---

## F5 -- Claim-list steering: **UPHELD**

* **What survives:** The runner's `C1`–`C11` claim list omitted `horizon_wga --graph` distant-relative mapping (`docs/PANGENOME_SUITE.md:12-17`), `docs/BENCHMARK_RESULTS.md:14-16` omitted the `cactus241way` `chr22 1MB` recall collapse (`9.9%` at line 37 vs `45.3%` on 38 kb `chr22_random` at line 36), and `scripts/benchmarks/run_all_evals.py:59-75` tunes exclusively on circular mitochondrial `chrM` (`multiz100way`) while testing on nuclear `chr22_KI270735v1_random`.
* **Diagnosis:** **UPHELD** (with the minor clarification that `C1`–`C11` was drawn from `STRATEGIC_PLAN.md` and `PANGENOME_SUITE.md` rather than `BENCHMARK_RESULTS.md`, making the omission of `PANGENOME_SUITE.md:12-17` the direct claim-list steering gap, compounded by the `BENCHMARK_RESULTS.md` Executive Summary omission and the `chrM`/`chr22_random` compartment mismatch).
* **Repair:** **ENDORSED** (include `PANGENOME_SUITE.md:12-17` / distant-relative recall claims in scope, surface the `chr22 1MB` `9.9%` result in `docs/BENCHMARK_RESULTS.md` Executive Summary, and flag the `chrM` vs `chr22_random` tuning/test mismatch). **The fix is free** for the documentation/claim-list updates; fixing the tuning/test compartment mismatch requires re-running calibration on nuclear data.

---

## Verdict Summary

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | Factual absences (no HPRC benchmarks, no `horizon_clean`/`resolve`/`lift` binaries, no `vg`/`odgi` samples) are verified, but `REFUTED` overreaches by treating a forward-looking strategic plan ("We will focus...") as a falsified empirical report; free wording fix endorsed. |
| F2 | UPHELD | Present-tense claim in `STRATEGIC_PLAN.md:11` that `src/gfa/mod.rs` solves 12GB pangenome I/O is unsubstantiated: `stream_gfa` has zero callers/tests and zero GFA or `vg`/`odgi` benchmarks exist; free doc-qualification fix endorsed. |
| F3 | UPHELD | `extract_fasta_from_maf.py` pre-filters to homologous MAF blocks, `run_all_evals.py` only runs tiny `chrM`/`chr22_random` slices, and no graph/hairball metrics exist; `NARROWED` verdict and free wording fix endorsed. |
| F4 | UPHELD | Zero HLA/`chr6` or centromere datasets or tests exist in `docs/`, `scripts/benchmarks/`, or `src/`; `NARROWED` verdict and free wording fix endorsed. |
| F5 | UPHELD | `PANGENOME_SUITE.md:12-17` distant-relative mapping was omitted from `C1`–`C11`, `BENCHMARK_RESULTS.md:14-16` omits the `chr22 1MB` `9.9%` recall collapse, and `run_all_evals.py` tunes on `chrM` while testing on nuclear `chr22_random`. |
