# Round 4 Attack — L4 (Provenance)

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `33bede8`, 2026-11-19)  
**Lens:** `provenance`  
**Prior history, literature, & method cards:**
- Checked `scripts/round history` for `provenance`, `C1`, `C2`, `C3` and reviewed `reports/rounds/round-3/adjudication-L4.md` (which left `R1-IMPL`, `R1-HMM`, `R2-HMM-XDROP`, `R2-RESOLVE-LIFT`, and `R2-INDEX-RAM` open). Where a Round 4 point overlaps a prior adjudication (`round 3`'s `F3`, `F6`, `F7`), this report states what changed on `33bede8`.
- `../library` does not exist on disk (`ls -d ../library` failed with `No such file or directory`); `refs/LIBRARY_VERSION` is pinned to `d4df9122eddd84694b2c3215e5478c46f4314044`. No prior-art assertions from memory are made.
- Per `BRIEF-provenance.md`, this project keeps **no method cards**.

---

## F1 -- C1: `src/bin/` contains 10 binaries (not 9), and 4 of them do not use the MDL/BOCPD core

- **Claim:** `C1` (`manuscript/horizon_paper.md:18`) -- "Horizon addresses these coupled bottlenecks in a unified, pure-Rust suite of nine composable executables sharing one information-theoretic core."
- **Verdict:** **NARROWED**
- **Flaw:**
  1. `src/bin/` contains **10** executables, not `nine`. Commit `f777935` added `src/bin/horizon.rs` as a 10th CLI binary, and `manuscript/horizon_paper.md:101` itself states `all 10 CLI binaries` (`cargo test --no-fail-fast` compiles and runs unit tests for all 10 binaries).
  2. Four of the 10 binaries do not use the information-theoretic (MDL/BOCPD) core (`src/mdl/`, `src/bocpd.rs`): `src/bin/horizon.rs` has zero `horizon_aligner::` imports (it is a `std::process::Command` dispatcher); `src/bin/horizon_gfa.rs` imports only `horizon_aligner::gfa::*`; `src/bin/horizon_lift.rs` imports only `horizon_aligner::coord::Interval` and `horizon_aligner::gfa::*`; and `src/bin/horizon_merge.rs` imports only `horizon_aligner::merge_pairwise_mafs` (interval sweep-line netting).
- **Evidence:**
  - `ls -la src/bin/`: `bench_chr1.rs`, `horizon.rs`, `horizon_clean.rs`, `horizon_eval.rs`, `horizon_gfa.rs`, `horizon_guide.rs`, `horizon_lift.rs`, `horizon_merge.rs`, `horizon_resolve.rs`, `horizon_wga.rs` (10 files).
  - `grep -n "horizon_aligner::" src/bin/*.rs` confirms zero MDL/BOCPD imports in `horizon.rs`, `horizon_gfa.rs`, `horizon_lift.rs`, and `horizon_merge.rs`.
  - `Cargo.toml` and `Cargo.lock` confirm pure-Rust dependencies (`anyhow`, `flate2` with `miniz_oxide`, `rayon`, `thiserror`).
- **Correction:** In `manuscript/horizon_paper.md:18`, replace `nine composable executables sharing one information-theoretic core` with `ten composable pure-Rust executables (a unified \`horizon\` dispatcher and nine specialized alignment, merging, and pangenome utilities) built over a shared library`.

---

## F2 -- C2: Full `chr6` (`451.35 s`) is pinned to `53130b8` (5 commits behind `HEAD`) and lacks saved runner scripts or profile/RSS logs

- **Claim:** `C2` (`manuscript/horizon_paper.md:101`) -- "full `GRCh38.chr6` (`170.8 Mbp`) × `CHM13.chr6` (`110.7 Mbp`) finished in `451.35 s` (`7m 31s`, `10,515.0 MB` RSS)—a `2.35×` total speedup over our `1,062.45 s` baseline (`12.3×` speedup in `LASTZ12` transition seeding from `588.67 s` to `47.87 s` via 16-bit prefix-bucketed `FlatSeedIndex` and Rayon parallelization; Figure 5)."
- **Verdict:** **NARROWED**
- **Flaw:**
  1. `docs/BENCHMARK_RESULTS.md:27,65` pins the `451.35 s` (`10,515.0 MB`) measurement to commit `53130b8` (and the `1,062.45 s` baseline to `1e47fe0`), whereas `HEAD` is `33bede8`—five commits later (`2695375`, `8c600a0`, `a4e1546`, `f777935`, `33bede8`). `git diff --stat 53130b8..HEAD -- src/` shows **21 modified Rust source files** (`951` insertions, `523` deletions), including `src/seeding/flat_index.rs`, `src/seeding/mod.rs`, `src/seeding/translated.rs`, `src/mdl/dp.rs`, `src/mdl/mod.rs`, and `src/egraph.rs`. Neither `manuscript/horizon_paper.md:10` nor `manuscript/horizon_paper.md:101` discloses that `451.35 s` predates these five commits.
  2. No committed script in `scripts/` or saved `HORIZON_PROFILE=1` / `/usr/bin/time -v` log in `reports/` records either the `1e47fe0` (`1,062.45 s`) or `53130b8` (`451.35 s`) run. `manuscript/scripts/generate_figures.py:669-685,712` embeds these timings and RSS values as hardcoded literals.
- **Evidence:**
  - `docs/BENCHMARK_RESULTS.md:27,64-65`
  - `git diff --stat 53130b8..HEAD -- src/` (`21 files changed, 951 insertions(+), 523 deletions(-)`)
  - `manuscript/scripts/generate_figures.py:669-685,712`
- **Correction:** Disclose in `manuscript/horizon_paper.md:101` that the `451.35 s` (`10,515.0 MB`) and `1,062.45 s` (`10,022.0 MB`) runs were measured on `53130b8` and `1e47fe0` respectively, or re-run full `GRCh38.chr6` × `CHM13.chr6` on `HEAD` (`33bede8`) with a committed script and archived `/usr/bin/time -v` + `HORIZON_PROFILE=1` log.

---

## F3 -- C3: HPRC `v1.0` dataset is now on disk, but `C3` benchmarks predate subsequent `horizon_gfa.rs` / `gfa/mod.rs` edits and lack saved logs

- **Claim:** `C3` (`manuscript/horizon_paper.md:101`) -- "On the `12 GB` HPRC `v1.0` Minigraph-Cactus pangenome graph (`hprc-v1.0-mc-grch38.gfa.gz`), `horizon_gfa stats` streamed all `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks in `99.97 s` using `97.0 MB` peak RSS, while two-pass `horizon_gfa extract` isolated `GRCh38.chr22` + `CHM13.chr22` in `183.29 s` (`456.4 MB` RSS) and `GRCh38.chr6` + `CHM13.chr6` in `244.93 s` (`1,255.1 MB` RSS)."
- **Verdict:** **NARROWED**
- **Prior adjudication & what changed:** Extends `round 3`'s `F6` (`reports/rounds/round-3/adjudication-L4.md:49-53`). Since Round 3, `docs/BENCHMARK_DATASETS.md:46` was fixed from `v1.1` to `v1.0`, `/home/rcb112/project/gemini/bench/pg/data/benchmarks/hprc-v1.0-mc-grch38.gfa.gz` (`12G`) is present on disk, and `horizon_gfa stats` (`88.84s` -> `99.97s`) and `chr22` `extract` were re-measured at `435b456`.
- **Flaw:**
  1. `GRCh38.chr6` + `CHM13.chr6` `horizon_gfa extract` (`244.93 s`, `1,255.1 MB` RSS) remains explicitly tagged `b857377` in `docs/BENCHMARK_RESULTS.md:57`. `git diff --stat b857377..HEAD -- src/bin/horizon_gfa.rs src/gfa/mod.rs` shows `219` insertions and `91` deletions across `src/bin/horizon_gfa.rs` (`147` lines) and `src/gfa/mod.rs` (`163` lines).
  2. `horizon_gfa stats` (`99.97 s`, `97.0 MB`) and `chr22` `extract` (`183.29 s`, `456.4 MB`) were introduced at `435b456` (`git log -S "99.97" --oneline`), which still predates `111` insertions and `16` deletions in `src/bin/horizon_gfa.rs` and `src/gfa/mod.rs` across `435b456..HEAD` (including `abce515` and `33bede8`).
  3. No committed runner script or `/usr/bin/time -v` log in `scripts/` or `reports/` records any of the three `C3` runs.
- **Evidence:**
  - `ls -lh /home/rcb112/project/gemini/bench/pg/data/benchmarks/hprc-v1.0-mc-grch38.gfa.gz` (`12G`)
  - `docs/BENCHMARK_RESULTS.md:55-57`
  - `git diff --stat b857377..HEAD -- src/bin/horizon_gfa.rs src/gfa/mod.rs` and `git diff --stat 435b456..HEAD -- src/bin/horizon_gfa.rs src/gfa/mod.rs`
- **Correction:** Tag `C3` in `manuscript/horizon_paper.md:101` with the generating commits (`435b456` for `stats` and `chr22` `extract`, `b857377` for `chr6` `extract`), or re-run all three on `HEAD` (`33bede8`) with committed scripts and `/usr/bin/time -v` logs.

---

## F4 -- C4: `horizon_clean` severs `sup <= 1` discordant edges by default (not `walk_support == 0`) and is a no-op on `*`-sequence GFAs

- **Claim:** `C4` (`manuscript/horizon_paper.md:49`) -- "3. `horizon_clean` (`src/bin/horizon_clean.rs`) indexes `S`-segment flanks and scores candidate `L`-line junctions with `bocpd_gc_profile` and an MDL two-regime changepoint test (`--min-cp-bits`), severing unsupported (`walk_support == 0`) discordant edges into `HB:i:1` Horizon Bubbles."
- **Verdict:** **NARROWED**
- **Flaw:**
  1. `src/bin/horizon_clean.rs:260` defines the severing predicate as `let sever = (cfg.prune_unsupported && sup == 0) || (sup <= 1 && cp >= cfg.min_cp_bits);`. Under the default configuration (`prune_unsupported = false`, line 476), `horizon_clean` severs discordant edges (`cp >= min_cp_bits`) with `sup <= 1`—including edges with `walk_support == 1`, not solely `walk_support == 0`. Conversely, when `--prune-unsupported` is passed, `sup == 0` edges are severed unconditionally without requiring a changepoint (`cp >= min_cp_bits`).
  2. For `*`-sequence GFA records (`s.sequence == b"*"`), `prefix` and `suffix` are empty vectors (lines 154-175), so `score_junction_cp_bits` returns `0.0` at line 101 (`if flank_u.len() < 12 || flank_v.len() < 12 { return 0.0; }`), making changepoint severing a no-op on `*`-sequence GFAs (as recorded in `docs/BENCHMARK_RESULTS.md:66`).
- **Evidence:**
  - `src/bin/horizon_clean.rs:100-127,154-175,260,305,476`
  - `docs/BENCHMARK_RESULTS.md:66` (`2,776 kept, 0 severed (no-op on *-seq)`)
- **Correction:** In `manuscript/horizon_paper.md:49`, replace `severing unsupported (walk_support == 0) discordant edges` with `severing low-support (walk_support <= 1) discordant edges on sequence-bearing S records (and optionally all walk_support == 0 edges under --prune-unsupported)`.

---

## F5 -- C5: `87.7%`, `54.4%`, `66.8%`, and the *Cannabis sativa* (`cs10` × `Cherry_Pie`) benchmark have zero code, script, dataset, or log provenance in the repository

- **Claim:** `C5` (`manuscript/horizon_paper.md:65`) -- "End-to-end `horizon_wga` abstained on `87.7%` of those MDL-rejected columns while recovering `54.4%` of MDL-supported pairs (`66.8%` with `w = 10` sparse-exon fallback and `band_w = 180`). Applied with zero parameter changes to intraspecific *Cannabis sativa* (`cs10` × `Cherry_Pie`, `45.5%` GC, `418,663` columns), Horizon retained `99.95%` (`418,465`) of reference `MAF` columns and `1.000` recall on scratch realignments."
- **Verdict:** **REFUTED**
- **Flaw:** None of the `C5` measurements, datasets, or evaluation scripts exist anywhere in `src/`, `scripts/`, `docs/`, `reports/`, or `/home/rcb112/project/gemini/bench/pg/data/benchmarks/`. Every occurrence of `87.7`, `54.4`, `66.8`, `418,663`, `418,465`, `99.95`, `cs10`, and `Cherry_Pie` in the repository is prose or hardcoded plotting literals introduced in manuscript commit `ca01d24`.
- **Evidence:**
  - `grep -rn "87.7\|54.4\|66.8\|418,663\|418663\|418,465\|418465\|99.95\|Cherry_Pie\|cs10" src/ scripts/ docs/` returns zero matches.
  - `git log -S "Cherry_Pie" --oneline` returns only `ca01d24 docs(manuscript): add Horizon academic paper, technical companion, and BioVis figures`.
  - `manuscript/scripts/generate_figures.py:217-222` hardcodes `kept = [99.95, 54.4, 12.3]`, `abstained = [0.05, 45.6, 87.7]`, and `n=418,663` as Python literals.
  - Neither `docs/BENCHMARK_DATASETS.md` nor `/home/rcb112/project/gemini/bench/pg/data/benchmarks/` contains *P. falciparum* × *P. vivax* Chr1 or *Cannabis sativa* (`cs10` × `Cherry_Pie`).
- **Correction:** Remove the unsubstantiated *P. falciparum* × *P. vivax* abstention percentages and *Cannabis sativa* (`cs10` × `Cherry_Pie`) claim from `manuscript/horizon_paper.md:10,65,69` and `manuscript/scripts/generate_figures.py:217-222`, or commit the datasets, evaluation scripts, and `HEAD`-pinned logs that produce them.

---

## F6 -- C6: The `30`-trial randomized e-graph vs. greedy synteny benchmark has zero provenance in the repository

- **Claim:** `C6` (`manuscript/horizon_paper.md:99`) -- "Across `30` randomized tandem-duplication and overlapping-inversion trials, e-graph saturation recovered `1.00` ancestral block sets and `1.00` exact `dup(X)` duplications compared to `0.00` under greedy one-to-one synteny filtering."
- **Verdict:** **REFUTED**
- **Flaw:** There is no `30`-trial randomized simulation harness, no greedy one-to-one synteny filter baseline implementation, and no benchmark output log anywhere in `src/`, `scripts/`, `docs/`, or `reports/`. `src/egraph.rs` contains deterministic unit tests only.
- **Evidence:**
  - `grep -rn "greedy\|30.*trial\|egraph.*trial\|randomized" src/ scripts/ docs/` returns zero matches.
  - `grep -rn "30" src/egraph.rs` returns zero matches.
  - `docs/BENCHMARK_RESULTS.md` contains no `30`-trial e-graph vs. greedy synteny row (only line 67's synthetic `horizon_resolve --collapse-dups` 200-segment repeat run).
- **Correction:** Remove the `30`-trial randomized comparison against greedy one-to-one synteny filtering from `manuscript/horizon_paper.md:99`, or commit the simulation script, greedy baseline, and saved trial results.

---

## F7 -- Meta-finding: Claim-list steering omits unlogged `horizon_merge` / Table 1 manuscript claims and undischarged `R2-HMM-XDROP` / `R2-RESOLVE-LIFT` requirements

- **Flaw:** The runner's Round 4 claim list (`C1`–`C6`) steers away from two major provenance gaps in `manuscript/horizon_paper.md`:
  1. `manuscript/horizon_paper.md:10,63,75` and Table 1 (`lines 81-95`) report `horizon_merge` scaling (`406` mammals in `3.84 s`, `317` `470-way` species in `2.53 s`, `92` `100-way chrM` in `3.07 s`, `92_mammals` in `1.38 s`, `65_fish` in `3.66 s`, `60_amniotes` in `16.61 s`, `27_sauropsids` in `0.62 s`), `horizon_guide` (`272.03 s`), and `~3 Mb` `60_amniotes` ablation speedups (`223.56 s` to `25.49 s`), whereas `docs/BENCHMARK_RESULTS.md:9` explicitly states `horizon_merge: Implemented in src/bin/horizon_merge.rs (best-score interval netting; unit-tested; no multi-species benchmark row yet)`.
  2. Extending `round 3`'s `F7` (`reports/rounds/round-3/adjudication-L4.md:57-61`), `R2-HMM-XDROP` (all WGA benchmarks on a single pinned commit) and `R2-RESOLVE-LIFT` (head-to-head `minimap2`/`odgi`/`vg` pangenome baselines) remain unfulfilled on `33bede8` (`docs/BENCHMARK_RESULTS.md:5`: `head-to-head minimap2/odgi/vg/LASTZ baselines remain pending external tool availability on this host`).

---

## Claim Verdict Summary

| Claim | Verdict | Provenance Summary |
| :--- | :--- | :--- |
| C1 | NARROWED | `src/bin/` has 10 binaries (not 9, matching `horizon_paper.md:101`), and 4 (`horizon`, `horizon_gfa`, `horizon_lift`, `horizon_merge`) do not import MDL/BOCPD |
| C2 | NARROWED | `451.35 s` (`53130b8`) and `1,062.45 s` (`1e47fe0`) predate 21 `src/` file edits across 5 commits to `HEAD` (`33bede8`) and lack saved scripts/logs |
| C3 | NARROWED | `hprc-v1.0-mc-grch38.gfa.gz` is on disk, but `chr6` `extract` (`b857377`) and `stats`/`chr22` `extract` (`435b456`) predate `horizon_gfa.rs`/`gfa/mod.rs` edits and lack logs |
| C4 | NARROWED | `src/bin/horizon_clean.rs:260` severs `sup <= 1` discordant edges by default (not `walk_support == 0`) and returns `0.0` cp bits on `*`-sequence GFAs |
| C5 | REFUTED | `87.7%`, `54.4%`, `66.8%`, and *Cannabis sativa* (`cs10` × `Cherry_Pie`, `418,663` cols) have zero code, script, dataset, or log provenance in the repo |
| C6 | REFUTED | `30`-trial randomized e-graph vs. greedy synteny benchmark has zero simulation script, greedy baseline, or log provenance in the repo |
