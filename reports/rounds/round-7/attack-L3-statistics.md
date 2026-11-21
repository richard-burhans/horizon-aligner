# Round 7 Adversarial Attack — L3 (Statistics)

**Environment & Provenance Disclosures (Rules 2–4):**
- **Literature library (Rule 2):** `../library` is not present (`NO_LIBRARY`). All literature checks were performed from [`refs/catalog_index.json`](file:///home/rcb112/project/gemini/bench/horizon-aligner/refs/catalog_index.json), [`manuscript/references.bib`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/references.bib), and the project's own source/scripts rather than from `../library` full texts or summaries.
- **Method cards (Rule 3):** This project keeps **no method cards** yet. Statistical claims are checked against the cited methodology (`Angelopoulos et al. 2024` in `manuscript/references.bib:19-24`) and the repository's own implementation and scoring scripts.
- **Prior-round adjudications (Rule 4):** Checked `scripts/round history` for `Angelopoulos`, `n_cal`, `109,119,624`, `chrM`, `Wittler`, and `horizon_filter`. Round 4 L3 F1 ruled `C1` (unit-test count) had no statistical content (`UNCHALLENGED`); Round 4 L3 F7 (`NARROWED`) ruled that `n_cal = 15, alpha = 0.07` under Angelopoulos et al. (2024) forces `0/15` errors and that `min_delta_bits = 10.0` was hardcoded; Round 5 L2 F3 and L7 F5 (`UPHELD`) flagged `chrM` `n = 1` sampling in `horizon_guide`; Round 6 `R2-RESOLVE-LIFT` carried missing `numpy` for `score_paf_vs_hprc.py`. Each finding below states what changed on `HEAD` (`5cb83fb`).

---

## F1 -- C3: `manuscript/horizon_paper.md:40` falsely claims offline-calibrated constant `λ̂ = 10.0` CRC, while `src/mdl/dp.rs:120-165` runs online `n_cal = 15` decoy calibration that violates exchangeability

- **Claim:** `C3` (`manuscript/horizon_paper.md:40`)
- **Verdict:** `REFUTED`
- **Prior-round context:** Round 4 L3 F7 (`NARROWED`) established that `n_cal = 15, alpha = 0.07` forces `0/15` empirical errors (`1/16 = 0.0625 <= 0.07`) and noted `min_delta_bits = 10.0` was hardcoded in `src/mdl/mod.rs`. Since Round 4, `calibrate_null_crc_for_gc` was added in [`src/mdl/dp.rs:82-165`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/dp.rs#L82-L165) and wired into [`src/mdl/mod.rs:401-429`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/mod.rs#L401-L429), but `manuscript/horizon_paper.md:40` was not updated.
- **Concrete flaw & evidence:**
  1. **Direct contradiction between text and code:** `C3` states Horizon uses *"fixed Conformal Risk Control defaults calibrated offline at risk level `α = 0.07` (`n_cal = 15`; Angelopoulos et al. 2024)"* and *"enforces a constant bit-savings floor `λ̂ = 10.0` bits (`min_delta_bits = 10.0`)"*. In reality, [`src/mdl/mod.rs:428`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/mod.rs#L428) calls `dp::calibrate_null_crc_for_gc(observed_gc, 15, 0.07, ctx.effective_min_bits, &hmm_params)` **online** on every alignment with non-empty `dp_jobs`, producing a data-dependent `λ̂ > 10.0` (the unit test at [`src/mdl/dp.rs:820`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/dp.rs#L820) itself verifies `lambda_hat > 12.6`). Conversely, [`src/bin/horizon_filter.rs:199`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_filter.rs#L199) does not call `calibrate_null_crc_for_gc` at all and uses `(KARLIN_MDL_COEFF * search_space_bits).max(10.0)`.
  2. **Exchangeability violation & degenerate sample-maximum estimator:** `calibrate_null_crc_for_gc` ([`src/mdl/dp.rs:120-165`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/dp.rs#L120-L165)) draws `n_cal = 15` synthetic i.i.d. `64`-bp random windows (`win_len = 64`, `band_w = 16`) at `observed_gc` (estimated from only the first `4096` bases of `seq_x` and `seq_y`), then `calibrate_crc_threshold` ([`src/mdl/dp.rs:86-114`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/dp.rs#L86-L114)) selects `scores[14]`—the sample maximum—because `(n_ge + 1) / 16 <= 0.07` holds only for `n_ge = 0` (`1/16 = 0.0625`). Under Angelopoulos et al. (2024), finite-sample risk control requires **exchangeability** between calibration and test items. Synthetic i.i.d. `64`-bp Bernoulli windows are not exchangeable with real genomic inter-anchor gaps (up to `2,400` bp, `band_w = 180`, with repeats, tandem tracts, and compositional shifts), nor even with `10,000`-bp i.i.d. random sequences (`scripts/benchmarks/bench_random_gc_crc.py` still emits `80–448` false-positive pairs on `4/5` `18%`-GC seeds and `288`/`149` pairs on `30%`/`45%`-GC draws).
- **Correction:** Replace `C3` in [`manuscript/horizon_paper.md:40`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-7/ATTACKED-horizon_paper.md#L40) to state accurately that `horizon_wga` runs an **online** `n_cal = 15` synthetic `64`-bp i.i.d. decoy maximum heuristic (`calibrate_null_crc_for_gc`, `alpha = 0.07`, `λ̂ >= max(10.0, 0.55 log2(NM))`) rather than an offline-calibrated constant `10.0`-bit CRC guarantee, and disclose that synthetic `64`-bp decoys are not exchangeable with genomic gaps (while `horizon_filter` uses the fixed `max(10.0, 0.55 log2(NM))` floor without decoy calibration).

---

## F2 -- C2: `SDUST+MDL` removes `12` p-arm records (`0.0052%` of `aln_len`), leaving `99.65%` cond. / `99.32%` uncond. precision unchanged, and `score_paf_vs_hprc.py` lacks unconditional precision and `numpy`

- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:86`)
- **Verdict:** `NARROWED`
- **Prior-round context:** Round 6 `R2-RESOLVE-LIFT` carried *"Install numpy to re-score 42700e1 horizon_filter 565/553 chr6 outputs against 2D HPRC truth"*. On `HEAD`, `numpy` is still absent in the environment.
- **Concrete flaw & evidence:**
  1. **Negligible statistical effect size (`0.0052%`) confined to off-target p-arm:** Re-running `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf` on `HEAD` with `--net` (`565` records, `elapsed=0.55s`) vs. `--net --query /tmp/hprc_grch38_chr6.fa --target /tmp/hprc_chm13_chr6.fa` (`553` records, `filtered_mdl=436`, `elapsed=4.92s`) shows total `aln_len` drops from `109,995,785` bp to `109,990,091` bp—a net difference of **`5,694` bp (`0.0052%`)**, consisting of `15` removed records and `3` newly selected records. Every one of those `15` removed and `3` added records lies on the `GRCh38.chr6` p-arm (`< 46,000,000` bp), outside the `CHM13.chr6` q-arm (`61401689–172126628`) 1-to-1 truth (`110,082,990` bp). Consequently, `SDUST+MDL` changes conditioned precision (`exact_pairs / pairs_on_shared`) by **0.0000%** and unconditional precision by `< 0.0052` percentage points (invisible at 2 decimal places: `99.65%` / `99.32%` before and after). Almost all filtering (`56,857 -> 565`, `99.0%` of records) is performed by `--net`, not MDL.
  2. **Scoring script gap & `n = 1` scope:** [`scripts/benchmarks/score_paf_vs_hprc.py:165-168`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/score_paf_vs_hprc.py#L165-L168) computes `exact_coord_precision_on_shared` (`exact_pairs / pairs_on_shared`) but **does not compute or emit unconditional precision** (`exact_pairs / total_pairs`, `99.32%`), and fails with `ModuleNotFoundError: No module named 'numpy'` on this host. Moreover, [`scripts/benchmarks/bench_vgp_filter_locus.py:44,65`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/bench_vgp_filter_locus.py#L44-L65) passes `--truth` to `score_paf_vs_hprc.py`, which does not define `--truth` in its `argparse` (`score_paf_vs_hprc.py:9-14`). Finally, the benchmark is `n = 1` chromosome (`chr6`) with no uncertainty bounds.
- **Correction:** Narrow `C2` in [`docs/BENCHMARK_RESULTS.md:86`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-7/ATTACKED-BENCHMARK_RESULTS.md#L86) to note that `--net --query --target` alters only `5,694` bp (`0.0052%` of `109.99 Mbp` `aln_len`, all on the p-arm `< 46 Mbp`) relative to `--net` alone, leaving both conditioned (`99.65%`) and unconditional (`99.32%`) precision unchanged to 2 decimal places, and add unconditional precision + `--truth` support to `score_paf_vs_hprc.py`.

---

## F3 -- C4: `verify_merge_guide_maf.py` only evaluates pairwise `hg38`-projected recall/precision on tuning-set `chrM`, leaving nuclear `horizon_guide` (`0.859 s`) and HAL ancestral nodes unverified

- **Claim:** `C4` (`manuscript/horizon_paper.md:84`)
- **Verdict:** `NARROWED`
- **Prior-round context:** Round 5 L2 F3 and L7 F5 (`UPHELD`) flagged `chrM` `n = 1` sampling. On `HEAD`, [`scripts/benchmarks/verify_merge_guide_maf.py`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/verify_merge_guide_maf.py) was added to verify 5-mammal `chrM.maf`.
- **Concrete flaw & evidence:**
  1. **Pairwise reference projection only, on a tuning-set organelle:** Re-running `python3 scripts/benchmarks/verify_merge_guide_maf.py` reproduces the four accuracy pairs (`hg38 x panTro4`: `0.955/0.998`, `rheMac3`: `0.892/0.971`, `mm10`: `0.718/0.975`, `canFam3`: `0.741/0.982`, `2,103` HAL bottom segments), with `elapsed_s=0.3150` (`+22%` vs. `0.258 s`, reflecting `n = 1` timing variance). However, [`scripts/benchmarks/verify_merge_guide_maf.py:146-155`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/verify_merge_guide_maf.py#L146-L155) evaluates **only pairwise `(hg38, tgt_sp)` projections** from `guide_out.maf` on `~16.5 kb` haploid `chrM` (a **Tuning Set** in [`docs/BENCHMARK_RESULTS.md:38-39`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-7/ATTACKED-BENCHMARK_RESULTS.md#L38-L39)). It does **not** score non-reference pairs (e.g., `mm10 x canFam3`), simultaneous 5-way multi-species columns, or the contents of the emitted `.hal` ancestral nodes (`Anc0`–`Anc3`; lines `139-143` only print `horizon_wga --hal-info` header metadata).
  2. **Zero accuracy check on the 5-mammal nuclear run:** The preceding nuclear claim in `C4` (`hg38`, `Macaca_mulatta`, `HLcalPym1`, `Equus_caballus`, `Canis_lupus` in `0.859 s`, `1,109` MAF blocks) reports only block counts and `n = 1` runtime, with zero precision/recall evaluation against reference MAF.
- **Correction:** Qualify `C4` in [`manuscript/horizon_paper.md:84`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-7/ATTACKED-horizon_paper.md#L84) to state that the `chrM.maf` verification (`n = 1`, tuning set) measures pairwise `hg38`-projected MAF recall/precision only (not 5-way column accuracy or HAL ancestral-node sequence fidelity), and that the 5-mammal nuclear `horizon_guide` run (`0.859 s`) has no MAF accuracy verification.

---

## F4 -- C1 and C5: Non-statistical claims outside the statistics lens

- **Claims:** `C1` (`manuscript/horizon_paper.md:84`) and `C5` (`manuscript/horizon_paper.md:150`)
- **Verdict:** `UNCHALLENGED` (for both `C1` and `C5` under L3)
- **Evidence:**
  - `C1`: Running `/home/agent/.cargo/bin/cargo test --no-fail-fast` on `HEAD` (`5cb83fb`) passes `88` library tests (`1` `.hal` test ignored) + `34` binary tests (`122` active, `123` total) and `/home/agent/.cargo/bin/cargo clippy --all-targets` exits `0` with zero warnings. Consistent with Round 4 L3 F1, a unit-test pass count carries no statistical inference claim.
  - `C5`: `Wittler R (2023)` is present in [`manuscript/references.bib:481-486`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/references.bib#L481-L486) and `refs/catalog_index.json` and contains no statistical claim.

---

## Verdict Table

| claim_id | verdict | summary |
| :--- | :--- | :--- |
| C1 | UNCHALLENGED | Non-statistical unit-test/clippy count (`88` lib + `34` bin = `122` active + `1` ignored `.hal` test verified on `HEAD`) |
| C2 | NARROWED | `SDUST+MDL` removes net `12` p-arm records (`5,694` bp = `0.0052%` of `109.99 Mbp`), leaving `99.65%` cond. / `99.32%` uncond. precision unchanged; `score_paf_vs_hprc.py` lacks uncond. precision and `numpy` |
| C3 | REFUTED | `horizon_wga` runs **online** `n_cal = 15` i.i.d. `64`-bp decoy calibration (`λ̂ > 10.0`), contradicting offline constant `λ̂ = 10.0` text and violating Angelopoulos et al. (2024) exchangeability |
| C4 | NARROWED | `verify_merge_guide_maf.py` scores only pairwise `hg38`-projected accuracy on `~16.5 kb` tuning-set `chrM` (`0.315 s` vs `0.258 s` `n=1`); nuclear `0.859 s` run and HAL `Anc0–Anc3` have no accuracy check |
| C5 | UNCHALLENGED | Bibliography entry (`Wittler 2023`) in `references.bib:481-486` and `catalog_index.json`; no statistical claim |
