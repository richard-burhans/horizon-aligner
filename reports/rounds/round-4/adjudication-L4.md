# Round 4 ADJUDICATION of L4 -- Provenance Lens

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `33bede8`, 2026-11-19)  
**Attack file:** `reports/rounds/round-4/attack-L4-provenance.md` (SHA-256 `ed610d5364acd43d72dbc23a477323ed070d584f2fc76c1072be4c34951dcd66`, matches sealed digest `ed610d5364ac`)  
**Prior history & method cards:** Checked `scripts/round history "provenance" --kind adjudication` and reviewed `reports/rounds/round-3/adjudication-L4.md` (which left `R2-HMM-XDROP`, `R2-RESOLVE-LIFT`, and `R2-INDEX-RAM` open). No method cards exist in `refs/methods/` in this repository; F1–F7 are codebase, dataset, and benchmark-provenance verification checks.

---

## F1 -- C1: `src/bin/` contains 10 binaries (not 9), and 4 of them do not use the MDL/BOCPD core

- **Diagnosis:** **UPHELD** (attacker's `NARROWED` ruling on `C1` stands). Independently verified on `33bede8`: `ls -1 src/bin/` lists **10** Rust binaries (`bench_chr1.rs`, `horizon.rs`, `horizon_clean.rs`, `horizon_eval.rs`, `horizon_gfa.rs`, `horizon_guide.rs`, `horizon_lift.rs`, `horizon_merge.rs`, `horizon_resolve.rs`, `horizon_wga.rs`), matching `manuscript/horizon_paper.md:101` (`all 10 CLI binaries`) and contradicting `line 18` (`nine`). `grep -n "horizon_aligner::" src/bin/*.rs` confirms that 4 binaries do not import MDL/BOCPD (`horizon.rs` has zero `horizon_aligner::` imports; `horizon_gfa.rs` imports only `gfa::*`; `horizon_lift.rs` imports only `coord::Interval` and `gfa::*`; `horizon_merge.rs` imports only `merge_pairwise_mafs`). `Cargo.toml` confirms pure-Rust dependencies (`anyhow`, `flate2`, `rayon`, `thiserror`).
- **Repair:** **ENDORSED.** Replacing `nine composable executables sharing one information-theoretic core` in `manuscript/horizon_paper.md:18` with `ten composable pure-Rust executables (a unified \`horizon\` dispatcher and nine specialized alignment, merging, and pangenome utilities) built over a shared library` requires no re-measurement; **the fix is free**.
- **What survives:** `Horizon` is a pure-Rust suite, but `src/bin/` contains 10 binaries rather than 9, and 4 of them (`horizon`, `horizon_gfa`, `horizon_lift`, `horizon_merge`) do not use the MDL/BOCPD core.

---

## F2 -- C2: Full `chr6` (`451.35 s`) is pinned to `53130b8` (5 commits behind `HEAD`) and lacks saved runner scripts or profile/RSS logs

- **Diagnosis:** **UPHELD** (attacker's `NARROWED` ruling on `C2` stands). Independently verified: (1) `docs/BENCHMARK_RESULTS.md:27,64-65` pins `451.35s` (`10,515.0 MB`) to `53130b8` and `1,062.45s` (`10,022.0 MB`) to `1e47fe0`; (2) `git log --oneline 53130b8..HEAD` shows 5 subsequent commits (`2695375`, `8c600a0`, `a4e1546`, `f777935`, `33bede8`) and `git diff --stat 53130b8..HEAD -- src/` confirms **21 modified Rust source files** (`951 insertions(+), 523 deletions(-)`), including `src/seeding/flat_index.rs`, `src/seeding/mod.rs`, `src/seeding/translated.rs`, `src/mdl/dp.rs`, `src/mdl/mod.rs`, and `src/egraph.rs`; (3) `manuscript/horizon_paper.md:10,101` omits the `53130b8` / `1e47fe0` commit tags (note: Figure 5's caption at `line 105` *does* mention `1e47fe0` and `53130b8`, so the omission is specific to the abstract and main text); (4) no runner script in `scripts/` or `/usr/bin/time -v` / `HORIZON_PROFILE=1` log in `reports/` records either run, and `manuscript/scripts/generate_figures.py:669-685,712` embeds the numbers as literals.
- **Repair:** **NARROWED (split by cost).** Disclosing in `manuscript/horizon_paper.md:101` (matching `line 105`) that `451.35 s` (`10,515.0 MB`) and `1,062.45 s` (`10,022.0 MB`) were measured on `53130b8` and `1e47fe0` requires no re-measurement (**the fix is free**); re-running full `GRCh38.chr6` × `CHM13.chr6` on `33bede8` with a committed script and archived logs is non-free (~451s compute).
- **What survives:** `C2`'s `451.35 s` (`53130b8`) and `1,062.45 s` (`1e47fe0`) benchmarks predate 21 `src/` file edits across 5 commits to `33bede8`, are untagged at `manuscript/horizon_paper.md:10,101`, and lack saved runner scripts or profile/RSS logs.

---

## F3 -- C3: HPRC `v1.0` dataset is now on disk, but `C3` benchmarks predate subsequent `horizon_gfa.rs` / `gfa/mod.rs` edits and lack saved logs

- **Diagnosis:** **UPHELD** (attacker's `NARROWED` ruling on `C3` stands). Independently verified: (1) `/home/rcb112/project/gemini/bench/pg/data/benchmarks/hprc-v1.0-mc-grch38.gfa.gz` (`12G`) is present on disk and `docs/BENCHMARK_DATASETS.md:46` now reads `v1.0` (resolving Round 3's `F6` dataset mismatch); (2) `docs/BENCHMARK_RESULTS.md:57` still pins `chr6` `extract` (`244.93s`, `1,255.1 MB`) to `b857377`, and `git diff --stat b857377..HEAD -- src/bin/horizon_gfa.rs src/gfa/mod.rs` confirms `219 insertions(+), 91 deletions(-)` (`147` and `163` lines changed); (3) `horizon_gfa stats` (`99.97s`, `97.0 MB`) and `chr22` `extract` (`183.29s`, `456.4 MB`) were introduced at `435b456`, and `git diff --stat 435b456..HEAD -- src/bin/horizon_gfa.rs src/gfa/mod.rs` confirms `111 insertions(+), 16 deletions(-)` since then; (4) no runner script or `/usr/bin/time -v` log exists in `scripts/` or `reports/`.
- **Repair:** **NARROWED (split by cost).** Tagging `C3` in `manuscript/horizon_paper.md:101` with its generating commits (`435b456` for `stats` and `chr22` `extract`, `b857377` for `chr6` `extract`) requires no re-measurement (**the fix is free**); re-running all three `C3` benchmarks on `33bede8` with committed scripts and `/usr/bin/time -v` logs is non-free (~528s compute).
- **What survives:** The `12 GB` HPRC `v1.0` GFA is on disk, but all three `C3` runs predate subsequent edits to `src/bin/horizon_gfa.rs` and `src/gfa/mod.rs` and lack saved runner scripts or `/usr/bin/time -v` logs.

---

## F4 -- C4: `horizon_clean` severs `sup <= 1` discordant edges by default (not `walk_support == 0`) and is a no-op on `*`-sequence GFAs

- **Diagnosis:** **UPHELD** (attacker's `NARROWED` ruling on `C4` stands). Independently verified in `src/bin/horizon_clean.rs`: (1) line `260` defines `let sever = (cfg.prune_unsupported && sup == 0) || (sup <= 1 && cp >= cfg.min_cp_bits);` with `prune_unsupported = false` by default (line `476`), so discordant edges (`cp >= min_cp_bits`) with `walk_support == 1` are also severed by default, while `--prune-unsupported` severs `sup == 0` edges unconditionally without requiring `cp >= min_cp_bits`; (2) for `*`-sequence GFA records (`s.sequence == b"*"`), `prefix` and `suffix` are empty vectors (lines `154-175`), causing `score_junction_cp_bits` (line `101`) to return `0.0` and making changepoint severing a no-op on `*`-sequence GFAs (matching `docs/BENCHMARK_RESULTS.md:66`).
- **Repair:** **ENDORSED.** Updating `manuscript/horizon_paper.md:49` to `severing low-support (walk_support <= 1) discordant edges on sequence-bearing S records (and optionally all walk_support == 0 edges under --prune-unsupported)` requires no re-measurement; **the fix is free**.
- **What survives:** `horizon_clean` severs `walk_support <= 1` discordant edges by default rather than strictly `walk_support == 0`, and changepoint severing is a no-op on `*`-sequence GFAs.

---

## F5 -- C5: `87.7%`, `54.4%`, `66.8%`, and the *Cannabis sativa* (`cs10` × `Cherry_Pie`) benchmark have zero code, script, dataset, or log provenance in the repository

- **Diagnosis:** **UPHELD** (attacker's `REFUTED` ruling on `C5` stands). Independently verified: `grep -rnE "87\.7|54\.4|66\.8|418,?663|418,?465|99\.95|Cherry_Pie|cs10" src/ scripts/ docs/ reports/` returns zero matches outside `reports/rounds/round-4/` round files; `/home/rcb112/project/gemini/bench/pg/data/benchmarks/` contains neither *P. falciparum* × *P. vivax* Chr1 nor *Cannabis sativa* (`cs10` × `Cherry_Pie`); and `manuscript/scripts/generate_figures.py:217-222` hardcodes `kept = [99.95, 54.4, 12.3]`, `abstained = [0.05, 45.6, 87.7]`, and `n=418,663` as Python literals.
- **Repair:** **NARROWED (split by cost).** Removing the unsubstantiated *P. falciparum* × *P. vivax* abstention percentages and *Cannabis sativa* (`cs10` × `Cherry_Pie`) claim from `manuscript/horizon_paper.md:10,65,69` and `manuscript/scripts/generate_figures.py:217-222` requires no re-measurement (**the fix is free**); acquiring the datasets, committing evaluation scripts, and running `33bede8`-pinned benchmarks is non-free.
- **What survives:** All `C5` percentages (`87.7%`, `54.4%`, `66.8%`, `99.95%`, `1.000`) and the *Cannabis sativa* (`cs10` × `Cherry_Pie`, `418,663` columns) benchmark lack any dataset, script, code, or log provenance in the repository.

---

## F6 -- C6: The `30`-trial randomized e-graph vs. greedy synteny benchmark has zero provenance in the repository

- **Diagnosis:** **UPHELD** (attacker's `REFUTED` ruling on `C6` stands). Independently verified: `grep -rnE "greedy|30.*trial|egraph.*trial|randomized" src/ scripts/ docs/` and `grep -rn "30" src/egraph.rs` return zero matches, and `docs/BENCHMARK_RESULTS.md` contains no `30`-trial e-graph vs. greedy synteny benchmark row (only line `67`'s synthetic `200`-segment `horizon_resolve --collapse-dups` run).
- **Repair:** **NARROWED (split by cost).** Removing the `30`-trial randomized comparison against greedy one-to-one synteny filtering from `manuscript/horizon_paper.md:99` requires no re-measurement (**the fix is free**); writing and running a simulation harness and greedy baseline is non-free.
- **What survives:** The `30`-trial randomized e-graph vs. greedy one-to-one synteny benchmark in `C6` has zero simulation script, greedy baseline, or log provenance in the repository.

---

## F7 -- Meta-finding: Claim-list steering omits unlogged `horizon_merge` / Table 1 manuscript claims and undischarged `R2-HMM-XDROP` / `R2-RESOLVE-LIFT` requirements

- **Diagnosis:** **UPHELD.** Independently verified: (1) `grep -rnE "3\.84|2\.53|3\.07|1\.38|3\.66|16\.61|0\.62|272\.03|223\.56|25\.49" src/ scripts/ docs/` returns zero matches, confirming that the `horizon_merge` multi-species runtimes (`406` mammals in `3.84 s`, etc.), `horizon_guide` (`272.03 s`), and `~3 Mb` `60_amniotes` ablation speedups in `manuscript/horizon_paper.md:10,63,75,81-95` have zero provenance in `docs/BENCHMARK_RESULTS.md`, `scripts/`, or `reports/`, and `docs/BENCHMARK_RESULTS.md:9` explicitly states `no multi-species benchmark row yet`; (2) `R2-HMM-XDROP` (all WGA benchmarks on a single pinned commit) and `R2-RESOLVE-LIFT` (`docs/BENCHMARK_RESULTS.md:5`: `head-to-head minimap2/odgi/vg/LASTZ baselines remain pending external tool availability on this host`) remain unfulfilled on `33bede8`.
- **Repair:** **NARROWED (split by cost).** Keeping `R2-HMM-XDROP` and `R2-RESOLVE-LIFT` open and removing/flagging the unlogged `horizon_merge`, `horizon_guide`, and `~3 Mb` ablation numbers in `manuscript/horizon_paper.md:10,63,75,81-95` is **free**; running the multi-species benchmarks and external tool baselines is non-free.
- **What survives:** Major manuscript benchmark claims (`horizon_merge`, `horizon_guide`, `60_amniotes` ablation, Table 1) lack repository provenance, and carried items `R2-HMM-XDROP` and `R2-RESOLVE-LIFT` remain undischarged on `33bede8`.

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `src/bin/` has 10 binaries (not 9), and 4 (`horizon`, `horizon_gfa`, `horizon_lift`, `horizon_merge`) do not import MDL/BOCPD; updating `manuscript/horizon_paper.md:18` is free |
| F2 | UPHELD | `451.35 s` (`53130b8`) and `1,062.45 s` (`1e47fe0`) predate 21 `src/` file edits across 5 commits to `33bede8`, are untagged at `manuscript/horizon_paper.md:10,101`, and lack saved scripts/logs; commit-tagging `line 101` is free |
| F3 | UPHELD | `hprc-v1.0-mc-grch38.gfa.gz` is on disk, but `chr6` `extract` (`b857377`) and `stats`/`chr22` `extract` (`435b456`) predate `horizon_gfa.rs`/`gfa/mod.rs` edits and lack saved logs; commit-tagging `manuscript/horizon_paper.md:101` is free |
| F4 | UPHELD | `src/bin/horizon_clean.rs:260` severs `sup <= 1` discordant edges by default (not `walk_support == 0`) and returns `0.0` cp bits on `*`-sequence GFAs; updating `manuscript/horizon_paper.md:49` is free |
| F5 | UPHELD | `87.7%`, `54.4%`, `66.8%`, and *Cannabis sativa* (`cs10` × `Cherry_Pie`, `418,663` cols) have zero dataset, script, code, or log provenance in the repo; redacting `horizon_paper.md:10,65,69` and `generate_figures.py:217-222` is free |
| F6 | UPHELD | `30`-trial randomized e-graph vs. greedy synteny benchmark has zero simulation script, greedy baseline, or log provenance in the repo; redacting `manuscript/horizon_paper.md:99` is free |
| F7 | UPHELD | `horizon_merge`, `horizon_guide`, and `~3 Mb` ablation numbers in `horizon_paper.md:10,63,75,81-95` have zero provenance (`docs/BENCHMARK_RESULTS.md:9` confirms `no multi-species benchmark row yet`), and `R2-HMM-XDROP` / `R2-RESOLVE-LIFT` remain open; redacting/flagging unlogged numbers and keeping carried items open is free |
