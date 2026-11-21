# Round 7 Cold Post-Round Audit

**Repository:** `horizon-aligner` (`/home/rcb112/project/gemini/bench/horizon-aligner`)
**Pinned baseline:** `5cb83fb6a1346a8193f9e4048055638d0a422b0a` (`5cb83fb`)
**Audit date:** 2026-11-21

---

## 1. Gate & Tree Status

| Check | Command | Result |
| :--- | :--- | :--- |
| Formatting | `cargo fmt --check` | **PASS** (exit 0) |
| Lints | `cargo clippy --all-targets --all-features -- -D warnings` | **PASS** (exit 0, 0 warnings) |
| Unit tests | `cargo test --all-features` | **PASS** (`122` active: `88` lib + `34` binary across 12 bin targets; `1` `#[ignore]`d `hal::tests::real_hal_decode_from_env`; `123` total) |
| Corrections | `./scripts/round corrections 7` | **PASS** (`6` applied: `COR-1`..`COR-6`; `0` not applied/malformed) |
| Drift | `./scripts/round drift 7` | **EXIT 1** (`DRIFTED`: `docs/BENCHMARK_RESULTS.md`, `manuscript/horizon_paper.md`; `UNEXAMINED`: `manuscript/references.bib`, `scripts/benchmarks/bench_vgp_filter_locus.py`; `MOVED`: `README.md`) — `drift_reviewed` entries are recorded in `reports/rounds/round-7/round.json:119-145` |
| Round ledger | `./scripts/round check 7 --full` & `./scripts/round balance` | **FAIL (exit 1)** — `r7: expected 25, derived 1 (-24)` |

---

## 2. Findings

### F1 (Blocking): Round 7 ledger does not reconcile (`-24` items lost)
- **Artifact / Command:** `./scripts/round check 7 --full`, `./scripts/round balance`, `reports/rounds/round-7/round.json:436-488`
- **Output:**
  ```
  FAIL the ledger does not reconcile at round 7: expected 25 outstanding, derived 1. An item was LOST or INVENTED -- most often a discharge or supersede that `outstanding()` does not honour. `round balance` shows the series.
  ```
  `round balance` reports: `r6` ended with `29` outstanding (`carried=25, filed=7, closed=3`), `r7` has `filed=3, closed=7, net=-4, carried=29, dropped=1, OUTSTANDING=1`, giving `expected 25, derived 1 (-24)`.
- **Cause:** `reports/rounds/round-7/round.json:436-449` (`open`) only retains the 3 narrowed residuals (`R2-RESOLVE-LIFT`, `R1-HMM`, `R2-HMM-XDROP`) from the 7 discharged Round 6 items, dropping/omitting 24 carried-over open items from earlier rounds.

### F2 (Minor): Dead `--truth` flag in `bench_vgp_filter_locus.py`
- **Artifact:** `scripts/benchmarks/bench_vgp_filter_locus.py:19,42-46,63-67` vs. `scripts/benchmarks/score_paf_vs_hprc.py:16`
- **Detail:** `COR-1` removed `--truth str(args.truth)` from the `score_paf_vs_hprc.py` subprocess calls because `score_paf_vs_hprc.py:16` hardcodes `TRUTH_TSV = "/tmp/hprc_chr6_truth_segments.tsv"` and does not accept `--truth`. However, `bench_vgp_filter_locus.py:19` still defines `--truth` as a CLI option and uses `args.truth.exists()` as the guard (`lines 42, 63`), so passing a custom `--truth /path/to/other.tsv` will silently score against `/tmp/hprc_chr6_truth_segments.tsv`.

### F3 (Minor): Stale baseline statement in `manuscript/horizon_paper.md` Results
- **Artifact:** `manuscript/horizon_paper.md:61` vs. `manuscript/horizon_paper.md:84` and `docs/BENCHMARK_RESULTS.md:88`
- **Detail:** `manuscript/horizon_paper.md:61` states: `head-to-head external pairwise aligner re-runs (FastGA, LASTZ, KegAlign) remain pending on the benchmark host (minimap2 -x asm5 is compared on full chr6 below).` However, `manuscript/horizon_paper.md:84` and `docs/BENCHMARK_RESULTS.md:88` report head-to-head runs of `FastGA -pafx -T8` (`2.894 s`, `321.7 MB`) and `LASTZ` (timed out at `15 s`) on the `5 Mbp` `chr6` pericentromere slice (`bench_external_wga_baselines.py`). Clarifying that `FastGA`/`LASTZ`/`KegAlign` remain pending *on the eight reference-MAF slices* (with `FastGA`/`LASTZ` already run on `5 Mbp` `chr6`) would remove the contradiction.

### F4 (Informational): Corrected files verified against code
- **`README.md:8` & `manuscript/horizon_paper.md:22,40` (`COR-2`, `COR-4`):** Online GC-conditioned CRC (`calibrate_null_crc_for_gc`, `n_cal = 15`, `α = 0.07`, floored by `max(10.0, 0.55 * log2(N * M))` bits, and `horizon_filter` using `0.55 * log2(N * M)` / `10.0`-bit floor without decoy calibration) matches `src/mdl/mod.rs:151,428` and `src/bin/horizon_filter.rs:199`.
- **`docs/BENCHMARK_RESULTS.md:11,58,86` & `manuscript/horizon_paper.md:10,84` (`COR-3`, `COR-5`):** Progressive ancestral-node HDF5 `v0` HAL (`2,103` bottom segments on 5-mammal UCSC `chrM.maf`, `0.258 s`, `Anc0..Anc3` cloning `hg38` on CLI-order pectinate tree), `chr6` `436` pre-netting vs. `12` post-netting p-arm `< 46 Mbp` record counts, `SweepGA`/`KegAlign` unbenchmarked note, and `122` active / `123` total unit test counts are consistent across both documents and `cargo test --all-features`.
- **`manuscript/references.bib:5,25,57` (`COR-6`):** Added DOIs for `adams2007bayesian`, `angelopoulos2024conformal`, and `chiaromonte2002scoring`.

---

## 3. Coverage & What Was Not Checked

- **Checked:** `git diff 5cb83fb6a134` across all 5 corrected text files; full reads/inspections of `README.md` (lines 1–60), `docs/BENCHMARK_RESULTS.md` (full), `manuscript/horizon_paper.md` (full), `manuscript/references.bib` (lines 1–80), `scripts/benchmarks/bench_vgp_filter_locus.py` (full), `scripts/benchmarks/score_paf_vs_hprc.py` (full), `reports/rounds/round-7/round.json`; `cargo fmt --check`, `cargo clippy`, `cargo test --all-features`, `./scripts/round corrections 7`, `./scripts/round drift 7`, `./scripts/round check 7 --full`, `./scripts/round balance`.
- **Not checked:** Did not inspect `manuscript/horizon_paper.pdf` binary contents; did not execute Python benchmark scripts against `/tmp/` datasets; did not run the `#[ignore]`d `real_hal_decode_from_env` test (`HORIZON_TEST_HAL` unset); did not read `attack-*.md`, `adjudication-*.md`, or `completeness-critic.md` (per cold-audit rules).
