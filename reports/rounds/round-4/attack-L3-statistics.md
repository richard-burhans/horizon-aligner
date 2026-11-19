# Round 4 Adversarial Attack -- L3 (Statistics)

- **Target:** `reports/rounds/round-4/ATTACKED-horizon_paper.md` (claims `C1`–`C6` from `reports/rounds/round-4/BRIEF-statistics.md`), verified against `HEAD` (`33bede8`), `docs/BENCHMARK_RESULTS.md`, and `src/`.
- **Literature & Method Cards (Rules 2–3):** `../library` is not present on disk (`ls: cannot access '../library': No such file or directory`), so literature checks rely on `refs/catalog_index.json` and `manuscript/references.bib` without asserting held-paper contents from memory. As declared in `BRIEF-statistics.md`, this project keeps **no method cards** (`refs/methods/` absent).
- **Prior Round History (Rule 4):** Checked via `scripts/round history "" --kind adjudication` and `reports/rounds/round-2/adjudication-L3.md`, `reports/rounds/round-3/adjudication-L3.md`. Prior statistics-lens adjudications relevant here are Round 2 L3-F2 (`NARROWED`, `horizon_clean` `0`-severed `*-seq` no-op), Round 2 L3-F3 (`UPHELD`, `MdlCostModel::block_bits` `2.0 * len_bp` charging), Round 3 L3-F5 (`UPHELD`, `horizon_gfa stats` `n = 1` counting pass with no variance/baseline), Round 3 L3-F6 (`UPHELD`, `C1` `UNCHALLENGED` under L3), and Round 3 L3-F7 (`UPHELD`, claim-list steering).

---

## F1 -- C1: No statistical content under the statistics lens

- **Claim:** `C1` (`manuscript/horizon_paper.md:18`) — *"Horizon addresses these coupled bottlenecks in a unified, pure-Rust suite of nine composable executables sharing one information-theoretic core."*
- **Verdict:** **UNCHALLENGED**
- **Flaw / Evidence:** Like Round 3's `C1` (adjudicated `UPHELD` as `UNCHALLENGED` under L3 in `reports/rounds/round-3/adjudication-L3.md:54-57`), this sentence states an architectural summary with no statistical estimator, test statistic, parameter, sample size, or benchmark metric. (The executable count `nine` vs. `10` `.rs` files in `src/bin/`—`bench_chr1.rs`, `horizon.rs`, `horizon_clean.rs`, `horizon_eval.rs`, `horizon_gfa.rs`, `horizon_guide.rs`, `horizon_lift.rs`, `horizon_merge.rs`, `horizon_resolve.rs`, `horizon_wga.rs`—is a non-statistical count left to the correctness lens.)
- **Correction:** None under L3.

---

## F2 -- C2: Full `chr6` `2.35×` speedup is `n = 1`, confounded by `+12.3%` more DP jobs and `+493 MB` RSS across 15 modified `src/` files, and predates `HEAD`

- **Claim:** `C2` (`manuscript/horizon_paper.md:101`) — *"full `GRCh38.chr6` (`170.8 Mbp`) × `CHM13.chr6` (`110.7 Mbp`) finished in `451.35 s` (`7m 31s`, `10,515.0 MB` RSS)—a `2.35×` total speedup over our `1,062.45 s` baseline (`12.3×` speedup in `LASTZ12` transition seeding from `588.67 s` to `47.87 s` via 16-bit prefix-bucketed `FlatSeedIndex` and Rayon parallelization; Figure 5)."*
- **Verdict:** **NARROWED**
- **Flaw / Evidence:**
  1. **Single-run (`n = 1`) timing with profiling overhead and no variance:** Both `1,062.45 s` (`1e47fe0`) and `451.35 s` (`53130b8`) in `docs/BENCHMARK_RESULTS.md:27, 64-65` are `n = 1` runs with `HORIZON_PROFILE=1` enabled, reporting no repeat-run variance (mean ± SD or IQR), thread count, or hardware specification, and no external aligner baseline (`minimap2`, `LASTZ`, `wfmash`).
  2. **Confounded comparison with non-equivalent algorithmic outputs:** `git diff --stat 1e47fe0..53130b8 -- src/` shows **15 modified files** (`948` insertions, `331` deletions), including `src/mdl/dp.rs`, `src/mdl/mod.rs`, `src/seeding/mod.rs`, and `src/bin/horizon_wga.rs` across commits `8c7d929`, `5d341a7`, `b068896`, and `53130b8`. As a result, `docs/BENCHMARK_RESULTS.md:64-65` shows the two runs did **not** compute the same alignment: `53130b8` produced `301,156` fwd / `298,178` rev LCBs and **`778,471` DP jobs** (`+12.3%` over `1e47fe0`'s `693,242` DP jobs) at `+493.0 MB` higher peak RSS (`10,515.0 MB` vs. `10,022.0 MB`), while reporting **zero** accuracy, precision, recall, or aligned-bp parity check for either run on full `chr6`.
  3. **Stale relative to `HEAD` (`33bede8`):** Commits `44d281d`, `2695375`, and `33bede8` subsequently modified `src/mdl/dp.rs` (`70` lines changed), `src/mdl/mod.rs` (`86` lines changed), `src/seeding/flat_index.rs`, `src/seeding/mod.rs`, and `src/seeding/translated.rs`.
- **Correction:** Narrow `manuscript/horizon_paper.md:101` (and `ATTACKED-horizon_paper.md:101`) to state that `451.35 s` (`10,515.0 MB` RSS) vs. `1,062.45 s` (`10,022.0 MB` RSS) is an `n = 1` profiled comparison between commits `53130b8` and `1e47fe0` (predating `HEAD` `33bede8`) across 15 modified `src/` files with `+12.3%` more DP jobs (`778,471` vs. `693,242`) and no aligned-bp or accuracy parity check, rather than a pure `FlatSeedIndex` + Rayon ablation.

---

## F3 -- C3: HPRC `horizon_gfa stats` and `extract` remain `n = 1` measurements without variance or external baselines (`stats` slowed `+12.5%` since Round 3)

- **Claim:** `C3` (`manuscript/horizon_paper.md:101`) — *"On the `12 GB` HPRC `v1.0` Minigraph-Cactus pangenome graph (`hprc-v1.0-mc-grch38.gfa.gz`), `horizon_gfa stats` streamed all `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks in `99.97 s` using `97.0 MB` peak RSS, while two-pass `horizon_gfa extract` isolated `GRCh38.chr22` + `CHM13.chr22` in `183.29 s` (`456.4 MB` RSS) and `GRCh38.chr6` + `CHM13.chr6` in `244.93 s` (`1,255.1 MB` RSS)."*
- **Verdict:** **NARROWED**
- **Flaw / Evidence:**
  1. **Prior adjudication & what changed:** Round 3 L3-F5 (`UPHELD` in `reports/rounds/round-3/adjudication-L3.md:45-50`) already established that `horizon_gfa stats` (`src/bin/horizon_gfa.rs:39-77`) is an `n = 1` single-run measurement for a record-discarding counter-only pass (`s_count`, `s_bases`, `l_count`, `p_count`, `w_count`) with no variance or `gfatools`/`odgi`/`vg` baseline. In Round 4 (`docs/BENCHMARK_RESULTS.md:19, 55`), `horizon_gfa stats` shifted from `88.84 s` / `97.1 MB` to `99.97 s` / `97.0 MB` (`+11.13 s`, `+12.5%`), and `horizon_gfa extract` added `GRCh38.chr6` + `CHM13.chr6` (`244.93 s`, `1,255.1 MB` RSS at `b857377`, predating `HEAD` `33bede8`).
  2. **Remaining statistical defect:** All three HPRC `horizon_gfa` figures (`99.97 s`, `183.29 s`, `244.93 s`) are `n = 1` single-run timings with no repeat-run variance (SD/IQR), no external `gfatools`/`odgi`/`vg` comparison (`docs/BENCHMARK_RESULTS.md:5` confirms baselines remain pending), and `chr6` `extract` is pinned to `b857377` rather than `HEAD` `33bede8`.
- **Correction:** Qualify `manuscript/horizon_paper.md:101` to note that `horizon_gfa stats` and `extract` timings are `n = 1` single-run measurements (`stats` is counter-only; `chr6` `extract` was measured at `b857377`) without repeat-run dispersion or external `gfatools`/`odgi`/`vg` baselines.

---

## F4 -- C4: `horizon_clean` severs `walk_support <= 1` (not `== 0`), uses a `1.5 log2(n)` penalty for a 1-df known-split test, and has zero real-data evaluation

- **Claim:** `C4` (`manuscript/horizon_paper.md:49`) — *"3. `horizon_clean` (`src/bin/horizon_clean.rs`) indexes `S`-segment flanks and scores candidate `L`-line junctions with `bocpd_gc_profile` and an MDL two-regime changepoint test (`--min-cp-bits`), severing unsupported (`walk_support == 0`) discordant edges into `HB:i:1` Horizon Bubbles."*
- **Verdict:** **REFUTED**
- **Flaw / Evidence:**
  1. **Contradicts `walk_support == 0` in code:** In `src/bin/horizon_clean.rs:260`, the severing predicate is:
     `let sever = (cfg.prune_unsupported && sup == 0) || (sup <= 1 && cp >= cfg.min_cp_bits);`
     Thus, when `cp >= cfg.min_cp_bits`, `horizon_clean` severs **`sup <= 1`** (including singleton-supported edges with `walk_support == 1`), directly contradicting `C4`'s claim that it severs `unsupported (walk_support == 0)` discordant edges.
  2. **Miscalibrated MDL penalty and causal BOCPD gate:** In `score_junction_cp_bits` (`src/bin/horizon_clean.rs:100-127`), line `112` sets `let mdl_penalty = 1.5 * n_tot.log2();` for comparing a 1-parameter null Bernoulli model against a 2-parameter two-regime Bernoulli model at a **fixed/known** junction `split = flank_u.len()`. A 1-df known-split BIC/MDL penalty is `0.5 * n_tot.log2()`, not `1.5 * n_tot.log2()`. Lines `121-126` then gate `cp_bits` on a one-sided causal `bocpd_gc_profile` jump between `split - 1` and `(split + 15).min(gc_prof.len() - 1)` (`bocpd_jump < 0.10`), an ad hoc threshold with no statistical calibration.
  3. **Zero real-data evaluation (disclosing Round 2 L3-F2 & Round 3 L3-F7):** As prior rounds noted and `docs/BENCHMARK_RESULTS.md:66` still records, the only real-data `horizon_clean` benchmark (`*-seq pairwise chr6 peri GFA`) kept `2,776` links and severed **`0`** (`no-op on *-seq`, because `has_seq = false` at `src/bin/horizon_clean.rs:154-176` makes `score_junction_cp_bits` return `0.0`). No empirical TPR/FPR or calibration of `--min-cp-bits` exists.
- **Correction:** Fix `manuscript/horizon_paper.md:49` to state `walk_support <= 1` (or fix `src/bin/horizon_clean.rs:260` to `sup == 0` and `0.5 * n_tot.log2()`), and disclose that the BOCPD + MDL junction test has not yet been evaluated on sequence-bearing real GFA data.

---

## F5 -- C5: *Plasmodium* abstention (`87.7%` / `54.4%`) is circularly defined on Horizon's own MDL threshold, *Cannabis* reports recall without precision, and neither has repo data/scripts

- **Claim:** `C5` (`manuscript/horizon_paper.md:65`) — *"End-to-end `horizon_wga` abstained on `87.7%` of those MDL-rejected columns while recovering `54.4%` of MDL-supported pairs (`66.8%` with `w = 10` sparse-exon fallback and `band_w = 180`). Applied with zero parameter changes to intraspecific *Cannabis sativa* (`cs10` × `Cherry_Pie`, `45.5%` GC, `418,663` columns), Horizon retained `99.95%` (`418,465`) of reference `MAF` columns and `1.000` recall on scratch realignments."*
- **Verdict:** **REFUTED**
- **Flaw / Evidence:**
  1. **Circular validation on *Plasmodium*:** `ATTACKED-horizon_paper.md:65` defines the `24,748 / 119,582` (`20.7%`) "MDL-rejected" vs. "MDL-supported" partition using **Horizon's own `10`-bit MDL+BOCPD threshold** on `LASTZ` `HoxD55` blocks, not an independent orthology or synteny ground truth. Reporting that `horizon_wga` abstained on `87.7%` of columns Horizon's own MDL scorer rejected and recovered `54.4%` (`66.8%` with `w = 10` fallback) of columns Horizon's own MDL scorer accepted is circular self-consistency—and still shows `horizon_wga` missed `45.6%` (`33.2%` with fallback) of its own MDL-supported pairs while failing to abstain on `12.3%` of MDL-rejected columns.
  2. **One-sided metric on *Cannabis sativa*:** `418,465 / 418,663 = 99.9527%` (`99.95%`) and `1.000` recall measure only retention/recall against a reference MAF on an `n = 1` aggregate without **precision**, false-positive rate, bubble count, or per-block/per-chromosome variance; high retention alone cannot demonstrate calibrated specificity.
  3. **Zero reproducible artifacts in repository:** Searching the entire repository (`rg -n "Cannabis|cs10|Cherry_Pie|24,748|119,582|418,465|418,663" .`) returns only manuscript text, `manuscript/scripts/generate_figures.py:217-222` (hardcoded plotting arrays), and round briefs. Neither `docs/BENCHMARK_RESULTS.md` nor `scripts/benchmarks/run_all_evals.py` contains *Plasmodium* or *Cannabis* benchmarks or evaluation scripts.
- **Correction:** Remove `C5`'s unverified *Plasmodium* and *Cannabis* figures from `manuscript/horizon_paper.md:10, 65, 69` until reproducible scripts/logs are committed, replace the circular MDL-partition metric with an independent orthology truth set, and report both precision and recall with per-block dispersion on *Cannabis*.

---

## F6 -- C6: `30` randomized e-graph trials have zero code, script, seed, or greedy baseline in the repository, and `saturate_synteny` lacks >2-copy rules

- **Claim:** `C6` (`manuscript/horizon_paper.md:99`) — *"Across `30` randomized tandem-duplication and overlapping-inversion trials, e-graph saturation recovered `1.00` ancestral block sets and `1.00` exact `dup(X)` duplications compared to `0.00` under greedy one-to-one synteny filtering."*
- **Verdict:** **REFUTED**
- **Flaw / Evidence:**
  1. **Completely absent from repository:** Searching the codebase (`rg -n "ancestral block|greedy|30.*trial|tandem.*trial" .`) finds zero Rust tests, zero Python/shell scripts, zero random seeds, and zero greedy-filtering implementations corresponding to these `30` trials. `src/egraph.rs:462-518` contains only **2** deterministic unit tests (`egraph_collapses_inverted_tandem_duplication` on one 5-block vector `[("A", 1), ("C", -1), ("B", -1), ("B", -1), ("D", 1)]` and `mdl_cost_model_gates_dup_by_sequence_divergence` on `[("B", 1), ("B", 1)]`).
  2. **Unstated statistical design & known `MdlCostModel` bias:** No generative model (block count, duplication copy number, inversion nesting depth, divergence), no CI/variance, and no definition of the `0.00` greedy baseline are provided. Moreover, `saturate_synteny` (`src/egraph.rs:165-273`) implements only 4 rewrite rules (`inv(inv(X))`, `inv(seq(A, B))`, associativity, `seq(X, X) -> dup(X)`), with no rule for >2 tandem copies or non-adjacent duplications, and Round 2 L3-F3 (`UPHELD`) already established that `MdlCostModel::block_bits` (`src/egraph.rs:325-328`) charges `2.0 * len_bp` per walk step.
- **Correction:** Remove the `30` randomized trials claim from `manuscript/horizon_paper.md:99` (or commit the exact simulation script, seed, greedy baseline, and trial outputs to the repository and report exact counts with confidence intervals).

---

## F7 -- Claim-list steering: uncalibrated CRC (`n_cal = 15`), Table 1 vs. `BENCHMARK_RESULTS.md` divergence, and missing external baselines

- **Flaw / Evidence:** As noted in `BRIEF-statistics.md:195` and Round 3 L3-F7, `C1`–`C6` steer away from three live statistical weaknesses in `ATTACKED-horizon_paper.md` and `docs/BENCHMARK_RESULTS.md`:
  1. **Uncalibrated Conformal Risk Control (`n_cal = 15`):** `ATTACKED-horizon_paper.md:40` claims CRC at $\alpha = 0.07$ with `n_cal = 15` yields $\hat{\lambda} = 10.0$ bits. For bounded loss $L \in [0, 1]$, finite-sample CRC (Angelopoulos et al. 2024) only guarantees $\mathbb{E}[L] \le \frac{n_{\text{cal}}}{n_{\text{cal}}+1}\alpha + \frac{1}{n_{\text{cal}}+1} = \frac{15}{16}(0.07) + \frac{1}{16} = 0.1281$ (nearly double $\alpha = 0.07$, with step resolution $1/16 = 0.0625$), while `min_delta_bits = 10.0` (`src/mdl/mod.rs:69`) and `KARLIN_MDL_COEFF = 0.55` (`src/mdl/mod.rs:26`) remain hardcoded constants in `AlignConfig::default()`.
  2. **Table 1 vs. `docs/BENCHMARK_RESULTS.md` divergence:** `ATTACKED-horizon_paper.md:81-95` (Table 1) reports `cactus447way` Laurasiatheria `33.73–41.45%` recall / `74.08–78.99%` precision and `470-way` primates `72.38–93.49%` recall / `96.69–99.46%` precision without per-pair sample sizes or variance, diverging from `docs/BENCHMARK_RESULTS.md:44-46` (`cactus447way` Horse `53.0%` / `65.4%`, Macaque `94.4%` / `84.6%`; `multiz470way` Marmoset `88.8%` / `86.1%`) with no reconciliation note.
  3. **No external baselines:** `docs/BENCHMARK_RESULTS.md:5` explicitly notes that head-to-head `minimap2`, `odgi`, `vg`, and `LASTZ` baselines remain pending.

---

## Verdict Table

| Claim | Verdict | Summary |
| :--- | :--- | :--- |
| C1 | UNCHALLENGED | No statistical content under L3 (architectural sentence only; matches Round 3 L3-F6) |
| C2 | NARROWED | `451.35 s` (`2.35×` / `12.3×`) is `n = 1` profiled, confounded by 15-file changes (`+12.3%` DP jobs, `+493 MB` RSS, no accuracy parity check), and predates `HEAD` `33bede8` |
| C3 | NARROWED | `horizon_gfa stats` (`99.97 s`, `+12.5%` vs. Round 3) and `extract` (`chr22`, `chr6` at `b857377`) are `n = 1` runs without variance or `gfatools`/`odgi`/`vg` baselines |
| C4 | REFUTED | `src/bin/horizon_clean.rs:260` severs `sup <= 1` (not `== 0`), line `112` uses `1.5 log2(n)` for a 1-df known-split test, and `docs/BENCHMARK_RESULTS.md:66` shows `0` severed (`*-seq` no-op) |
| C5 | REFUTED | *Plasmodium* (`87.7%`/`54.4%`) is circularly defined on Horizon's own `10`-bit MDL threshold, *Cannabis* (`99.95%`) omits precision/FPR, and neither has data/scripts in the repo |
| C6 | REFUTED | Zero code, script, seed, or greedy baseline exists in the repo for the `30` randomized trials (`src/egraph.rs` has only 2 fixed unit tests and lacks >2-copy `dup` rules) |
