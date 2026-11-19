# Round 4 ADJUDICATION -- Statistics Lens (L3)

- **Attack reviewed:** `reports/rounds/round-4/attack-L3-statistics.md` (sealed digest `c2e5dd17663f`, verified via `sha256sum`: `c2e5dd17663fa830bb49c61920f1a8f078fce8edbdfd2e7addd0b98a1d3d07f9`).
- **Method cards (`refs/methods/`):** Confirmed absent in `/home/rcb112/project/gemini/bench/horizon-aligner/refs/methods/` and `/home/rcb112/project/gemini/bench/refs/methods/`; no project method card governs these statistical disputes.
- **Prior round history (`scripts/round history "" --kind adjudication`):** Round 2 L3-F2 (`NARROWED`, `horizon_clean` `0`-severed `*-seq` no-op), Round 2 L3-F3 (`UPHELD`, `MdlCostModel::block_bits` `2.0 * len_bp` per-step charging), Round 3 L3-F5 (`UPHELD`, `horizon_gfa stats` `n = 1` counting pass), Round 3 L3-F6 (`UPHELD`, `C1` `UNCHALLENGED` under L3), and Round 3 L3-F7 (`UPHELD`, claim-list steering) are all accurately cited by the attacker.

---

## Finding-by-Finding Rulings

### F1 -- C1: No statistical content under the statistics lens
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; no L3 repair required.
  - **Verification:** `reports/rounds/round-4/ATTACKED-horizon_paper.md:18` states an architectural summary with no statistical estimator, parameter, sample size, or benchmark metric, matching Round 3 L3-F6. (`src/bin/` contains `10` `.rs` files—`bench_chr1.rs` plus `9` `horizon*.rs` binaries—left to L1 correctness.)
- **What survives:** `C1` is `UNCHALLENGED` under L3.

---

### F2 -- C2: Full `chr6` `2.35×` speedup is `n = 1`, confounded by `+12.3%` more DP jobs and `+493 MB` RSS across 15 modified `src/` files, and predates `HEAD`
- **Verdict:** **UPHELD**
- **Ruling:** Both the attacker's `NARROWED` diagnosis and proposed qualification repair are **UPHELD**; **the fix is free**.
  - **Verification:** Verified in `docs/BENCHMARK_RESULTS.md:27, 64-65` that both `1,062.45 s` (`1e47fe0`) and `451.35 s` (`53130b8`) are `n = 1` runs with `HORIZON_PROFILE=1` and no repeat-run variance or external aligner baseline. Verified `git diff --stat 1e47fe0..53130b8 -- src/` modifies `15` files (`948` insertions, `331` deletions), producing non-equivalent outputs (`301,156` fwd / `298,178` rev LCBs and `778,471` DP jobs vs. `693,242`, `+12.3%`, at `10,515.0 MB` vs. `10,022.0 MB` RSS, `+493.0 MB`) with no accuracy or aligned-bp parity check. Verified `git diff --stat 53130b8..HEAD -- src/` modifies `21` files (including `src/mdl/dp.rs` `70` lines, `src/mdl/mod.rs` `86` lines, `src/seeding/flat_index.rs`, `src/seeding/mod.rs`, `src/seeding/translated.rs`), so `53130b8` predates `HEAD` (`33bede8`).
- **What survives:** The recorded `451.35 s` (`2.35×` / `12.3×`) `n = 1` timings stand only when qualified as a profiled, pre-`HEAD`, 15-file comparison with non-equivalent DP job counts and no accuracy parity check.

---

### F3 -- C3: HPRC `horizon_gfa stats` and `extract` remain `n = 1` measurements without variance or external baselines (`stats` slowed `+12.5%` since Round 3)
- **Verdict:** **UPHELD**
- **Ruling:** Both the attacker's `NARROWED` diagnosis and proposed qualification repair are **UPHELD**; **the fix is free**.
  - **Verification:** Verified in `src/bin/horizon_gfa.rs:39-77` that `cmd_stats` is a record-discarding counter-only pass, and in `docs/BENCHMARK_RESULTS.md:5, 19, 55-57` that `horizon_gfa stats` (`99.97 s`, `97.0 MB` RSS, `+12.5%` vs. Round 3's `88.84 s`), `chr22` `extract` (`183.29 s`, `456.4 MB`), and `chr6` `extract` (`244.93 s`, `1,255.1 MB` at `b857377`) are `n = 1` measurements without repeat-run dispersion or external `gfatools`/`odgi`/`vg` baselines.
- **What survives:** `C3`'s counts and `n = 1` timings stand when qualified as single-run measurements (`stats` counter-only; `chr6` `extract` at `b857377`) without variance or external tool baselines.

---

### F4 -- C4: `horizon_clean` severs `walk_support <= 1` (not `== 0`), uses a `1.5 log2(n)` penalty for a 1-df known-split test, and has zero real-data evaluation
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; repair **NARROWED** by cost.
  - **Verification:** Verified in `src/bin/horizon_clean.rs:260` that `let sever = (cfg.prune_unsupported && sup == 0) || (sup <= 1 && cp >= cfg.min_cp_bits);` severs singleton-supported (`sup <= 1`) edges rather than only `walk_support == 0`. Verified in `score_junction_cp_bits` (`src/bin/horizon_clean.rs:100-127`) that line `112` sets `mdl_penalty = 1.5 * n_tot.log2()` for a 1-df known-split Bernoulli comparison (where BIC/MDL is `0.5 * n_tot.log2()`), gated by an uncalibrated `bocpd_jump < 0.10` check (lines `121-126`). Verified in `docs/BENCHMARK_RESULTS.md:66` that the only real-data benchmark (`*-seq pairwise chr6 peri GFA`) severed `0` of `2,776` links (leaving only two synthetic unit tests in `src/bin/horizon_clean.rs:558-626`).
  - **Repair:** Updating `manuscript/horizon_paper.md:49` to state `walk_support <= 1` and disclosing the absence of sequence-bearing real-data evaluation is **free**; modifying `src/bin/horizon_clean.rs:112, 260` is **not free** (code change requiring re-testing and benchmark runs).
- **What survives:** `C4` is refuted on its `walk_support == 0` predicate and lacks any sequence-bearing real-data evaluation or calibrated changepoint penalty.

---

### F5 -- C5: *Plasmodium* abstention (`87.7%` / `54.4%`) is circularly defined on Horizon's own MDL threshold, *Cannabis* reports recall without precision, and neither has repo data/scripts
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; repair **NARROWED** by cost.
  - **Verification:** Verified in `ATTACKED-horizon_paper.md:65` that the `24,748 / 119,582` (`20.7%`) *Plasmodium* `LASTZ` partition is defined by Horizon's own `10`-bit MDL+BOCPD threshold rather than an independent orthology truth set, and that *Cannabis* (`418,465 / 418,663 = 99.95%`, `1.000` recall) omits precision, FPR, and per-block dispersion. Searching the repository confirms zero data, scripts, or logs for either benchmark; only hardcoded arrays exist in `manuscript/scripts/generate_figures.py:217-222` (which also mismatch the manuscript text, summing `88,597 + 20,968 = 109,565` rather than `119,582`).
  - **Repair:** Removing the unverified *Plasmodium* and *Cannabis* figures from `manuscript/horizon_paper.md:10, 65, 69` is **free**; running new evaluations against an independent orthology truth set with precision and dispersion on *Cannabis* is **not free**.
- **What survives:** `C5` is refuted as an unverified, non-reproducible claim with circular *Plasmodium* validation and one-sided *Cannabis* metrics.

---

### F6 -- C6: `30` randomized e-graph trials have zero code, script, seed, or greedy baseline in the repository, and `saturate_synteny` lacks >2-copy rules
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; repair **NARROWED** by cost.
  - **Verification:** Verified zero scripts, seeds, or greedy-filtering baselines exist in the repository for the `30` randomized trials. `src/egraph.rs:462-518` contains only `2` deterministic unit tests (`egraph_collapses_inverted_tandem_duplication` and `mdl_cost_model_gates_dup_by_sequence_divergence`), `saturate_synteny` (`src/egraph.rs:165-273`) implements only 4 rewrite rules without >2-copy or non-adjacent duplication rules, and Round 2 L3-F3 already established `MdlCostModel::block_bits` (`src/egraph.rs:325-328`) charges `2.0 * len_bp` per walk step.
  - **Repair:** Removing the `30` randomized trials claim from `manuscript/horizon_paper.md:99` is **free**; implementing and executing the 30-trial simulation with seeds, greedy baseline, and confidence intervals is **not free**.
- **What survives:** `C6` is refuted because the `30` randomized trials are absent from the repository and unverified.

---

### F7 -- Claim-list steering: uncalibrated CRC (`n_cal = 15`), Table 1 vs. `BENCHMARK_RESULTS.md` divergence, and missing external baselines
- **Verdict:** **NARROWED**
- **Ruling:** Steering observation **UPHELD**; mathematical CRC detail **NARROWED**.
  - **Verification:** Verified `ATTACKED-horizon_paper.md:40` (`n_cal = 15`, $\alpha = 0.07$, $\hat{\lambda} = 10.0$ bits), `src/mdl/mod.rs:26, 69` (`KARLIN_MDL_COEFF = 0.55`, `min_delta_bits = 10.0`), Table 1 (`ATTACKED-horizon_paper.md:81-95`) diverging from `docs/BENCHMARK_RESULTS.md:44-46` (`cactus447way` Horse `53.0%` / `65.4%` vs. Table 1 `33.73–41.45%` / `74.08–78.99%`; `470-way` Marmoset `86.1%` precision vs. Table 1 `96.69–99.46%`), and pending external baselines (`docs/BENCHMARK_RESULTS.md:5`).
  - **Narrowing:** Under exact finite-sample CRC (Angelopoulos et al. 2024), choosing $\hat{\lambda}$ via $\frac{n}{n+1}\hat{R}_n(\lambda) + \frac{1}{n+1} \le \alpha$ still guarantees $\mathbb{E}[L] \le \alpha = 0.07$ when $\hat{R}_{15}(\lambda) \le 0.008$ (forcing zero empirical errors on all `15` calibration points); the `0.1281` bound only applies if $\hat{\lambda}$ were selected using unadjusted empirical risk $\hat{R}_{15}(\lambda) \le 0.07$. Either way, `n_cal = 15` is severely underpowered and `min_delta_bits = 10.0` is hardcoded.
- **What survives:** `C1`–`C6` omit underpowered `n_cal = 15` CRC, unreconciled Table 1 vs. `BENCHMARK_RESULTS.md` divergences, and missing external baselines.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `C1` (`manuscript/horizon_paper.md:18`) has no statistical content (`UNCHALLENGED` under L3, matching Round 3 L3-F6) |
| F2 | UPHELD | `C2` `451.35 s` (`2.35×` / `12.3×`) is `n = 1` profiled, confounded by 15-file changes (`+12.3%` DP jobs, `+493.0 MB` RSS, no accuracy parity check), and predates `HEAD` `33bede8`; doc qualification is free |
| F3 | UPHELD | `C3` `horizon_gfa stats` (`99.97 s`, `+12.5%` vs. Round 3) and `extract` (`chr22`, `chr6` at `b857377`) are `n = 1` runs without variance or `gfatools`/`odgi`/`vg` baselines; doc qualification is free |
| F4 | UPHELD | `horizon_clean.rs:260` severs `sup <= 1` (not `== 0`), line `112` uses `1.5 log2(n)` for a 1-df known-split test, and `BENCHMARK_RESULTS.md:66` shows `0` severed (`*-seq` no-op); doc fix is free, code fix is not free |
| F5 | UPHELD | *Plasmodium* (`87.7%`/`54.4%`) is circularly defined on Horizon's `10`-bit MDL threshold, *Cannabis* (`99.95%`) omits precision/FPR, and zero repo scripts/data exist (`generate_figures.py` also mismatches `119,582`); removing unverified claims is free, re-running is not free |
| F6 | UPHELD | Zero code, script, seed, or greedy baseline exists in the repo for the `30` randomized trials (`src/egraph.rs` has 2 unit tests and lacks >2-copy rules); removing the claim is free, running simulations is not free |
| F7 | NARROWED | Steering observation upheld (Table 1 vs. `BENCHMARK_RESULTS.md` divergence, missing baselines, hardcoded `10.0`/`0.55`), but exact CRC guarantees `E[L] <= 0.07` if empirical risk `<= 0.008` (`0/15` errors); `0.1281` applies only to unadjusted empirical risk |
