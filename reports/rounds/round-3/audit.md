# Round 3 Cold Post-Round Audit

**Date:** 2026-11-19  
**Target repo:** `/home/rcb112/project/gemini/bench/horizon-aligner`  
**Baseline pinned commit:** `488b020ccc0a73998915e975bad8e25599b16716` (`488b020`)  
**Scope read first:** `git diff 488b020` across 4 modified working-tree files (`docs/BENCHMARK_DATASETS.md`, `docs/BENCHMARK_RESULTS.md`, `docs/STRATEGIC_PLAN.md`, `scripts/benchmarks/run_all_evals.py`) before reading any round artifact.

---

## 1. Tree Greenness & Harness Gates

| Check | Command | Result | Notes |
|---|---|---|---|
| Clippy & Unit Tests | `export PATH="$HOME/.cargo/bin:$PATH" && cargo clippy --all-targets -- -D warnings && cargo test` | **PASS** | `0` clippy warnings; `46` passed (`37` in `src/lib.rs` + `1` in each of the `9` binaries in `src/bin/`), `1` ignored (`hal::tests::real_hal_decode_from_env`), `0` doc-tests. (`cargo` is not on the default shell `PATH`.) |
| Corrections | `scripts/round corrections 3` | **PASS** | `8 correction(s); 0 not applied or malformed` (`COR-1`..`COR-8`). |
| Drift | `scripts/round drift 3` | **FAIL (exit 1)** | `docs/STRATEGIC_PLAN.md` (`+4/-4`, `C1/C3/C4` quotes changed; `drift_reviewed` present in `round.json`), `docs/BENCHMARK_RESULTS.md` (`+16/-16`, `C2/C5/C6` quotes changed; `drift_reviewed` present in `round.json`), `scripts/benchmarks/run_all_evals.py` (`UNEXAMINED`: corrected but never attacked). |
| Full Round Check | `scripts/round check 3 --full` | **FAIL (exit 1)** | Ledger reconciliation failure (`expected -23 outstanding, derived 16`; `scripts/round balance` shows 16 filed vs 39 closed across rounds 1–3) and 39 `discharges` errors (`discharges names an item no round ever opened: '3-h3f39d5'`, etc.). |

---

## 2. Findings

### F1 — Unattacked & Unlisted Working-Tree Edit (`docs/BENCHMARK_DATASETS.md`)
- **Artifact:** `docs/BENCHMARK_DATASETS.md:46-48` (`git diff 488b020 -- docs/BENCHMARK_DATASETS.md`)
- **Finding:** `docs/BENCHMARK_DATASETS.md` was modified in the working tree (`hprc-v1.1-mc-grch38.gfa.gz` → `hprc-v1.0-mc-grch38.gfa.gz`; `GFA v1.1 graph generation` → `` `H`/`S`/`L`/`P`/`W` GFA v1.0/v1.1 graph streaming ``), but it is **not** in `round.json` `attacked`, **not** in `scripts/round corrections 3`, and **not** listed in `reports/rounds/round-3/AUDIT-BRIEF.md` under "Files this round says it corrected".

### F2 — Unattacked Script Edit & Sibling Path Drift (`scripts/benchmarks/run_all_evals.py`)
- **Artifact:** `scripts/benchmarks/run_all_evals.py:10, 46, 85`; `scripts/benchmarks/fetch_benchmarks.sh:4-5`
- **Finding:** `scripts/round drift 3` flags `scripts/benchmarks/run_all_evals.py` as `UNEXAMINED` ("CORRECTED BUT NEVER ATTACKED").
  - Line 10 (`BENCH_DIR = ROOT.parent / "pg" / "data" / "benchmarks"`) resolves to an existing directory on this host (`/home/rcb112/project/gemini/bench/pg/data/benchmarks`), **but** `scripts/benchmarks/fetch_benchmarks.sh:4-5` still uses relative `mkdir -p data/benchmarks && cd data/benchmarks` inside `horizon-aligner`, so the recovery script and evaluation script disagree on where benchmark files live.
  - Lines 46 and 85 hardcode `"/home/agent/.cargo/bin/cargo"`, and the script only evaluates `multiz100way`, `cactus447way`, and `multiz470way` (omitting `cactus241way`).

### F3 — Stale Prose Left Behind in `docs/STRATEGIC_PLAN.md` (Contradicts Round 3 Corrections)
- **Artifact:** `docs/STRATEGIC_PLAN.md:23` and `docs/STRATEGIC_PLAN.md:30` vs `docs/STRATEGIC_PLAN.md:3` and `docs/BENCHMARK_RESULTS.md:18-26, 54-65`
- **Finding:**
  1. `docs/STRATEGIC_PLAN.md:30` (Section 3.D Mitigation) still states: *"Focus our next implementation and benchmarking milestones on pangenome graph tooling (`horizon_clean`, `horizon_resolve`, `horizon_lift`, and HPRC `GRCh38` vs. `CHM13` evaluation, which is currently unrun)"*. This is stale and contradicted by `docs/BENCHMARK_RESULTS.md:18-26, 54-65`, which already reports runs for `horizon_clean`, `horizon_resolve`, `horizon_lift`, and HPRC `GRCh38` vs `CHM13` (`chr22`, `chr22` peri, `chr6` peri, `chr6` q-arm).
  2. `docs/STRATEGIC_PLAN.md:23` (Section 3.B Mitigation) still states: *"Support `stdin`/`stdout` GFA/MAF streaming on linear filter/transform tools (`horizon_clean`, `horizon_resolve`, `horizon_lift`)"*, whereas the Round 3 status note at `docs/STRATEGIC_PLAN.md:3` states: *"Only `horizon_gfa stats` supports single-pass `"-"` `stdin` streaming; two-pass tools (`horizon_clean`, `horizon_resolve`, `horizon_lift`) and `horizon_wga` require seekable file paths."*

### F4 — Derived Quantities & Code-to-Doc Verification
- **Artifacts:** `docs/BENCHMARK_RESULTS.md:16, 22-23, 26, 58`; `docs/STRATEGIC_PLAN.md:3, 11-13`; `src/bin/horizon_wga.rs:367-436`; `src/seeding/mod.rs:332-390`; `src/bin/horizon_resolve.rs:416-478`; `src/gfa/mod.rs:97-120`
- **Verified:**
  - `docs/BENCHMARK_RESULTS.md` derived quantities match their inputs: `61.3% F1` from `62.3%` recall / `60.3%` precision (`61.28%`); `68.5%` non-deduplicated ratio (`24,365,611 / 35,595,466 = 68.45%`); `1.164x` multiplicity (`5,106,541 / 4,387,655 = 1.1638x`, computed against `unique_x_bp`; vs `unique_y_bp` it is `1.122x`); `0.120 GB/Mbp` (`2,409.7 MB / 20 Mbp` in decimal GB) and `75.0%` deduplicated coverage (`14,999,823 / 20,000,000`).
  - `docs/STRATEGIC_PLAN.md:3, 11-13` code descriptions match `src/bin/` (`9` binaries, `46` passing + `1` ignored tests), `estimate_dag_hmm_params` (`src/bin/horizon_wga.rs:367-436`), `FlatSeedIndex` (`src/seeding/mod.rs:332-390`), `horizon_resolve` (`src/bin/horizon_resolve.rs:416-478`), and `stream_gfa` (`src/gfa/mod.rs:97-120`).

---

## 3. What Was Not Reached

- **Benchmark re-execution:** Did not re-run the benchmarks reported in `docs/BENCHMARK_RESULTS.md` (`cactus241way`, `chr22`/`chr6` pericentromere, `chr6` q-arm, 12 GB HPRC `horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_resolve`, `horizon_wga --graph`).
- **Benchmark script execution:** Did not execute `scripts/benchmarks/run_all_evals.py`.
- **Unmodified Rust modules:** Did not inspect `src/mdl/mod.rs`, `src/mdl/dp.rs`, `src/bocpd.rs`, `src/egraph.rs`, or `src/bin/horizon_clean.rs` line-by-line (no Rust source files were modified since `488b020`).
