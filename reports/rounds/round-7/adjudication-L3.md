# Round 7 ADJUDICATION — L3 (Statistics)

**Attack seal verified:** `sha256sum reports/rounds/round-7/attack-L3-statistics.md` = `3f3b7ac96ac3f2d8a191e4bcc553fd07ada1a3ff6110981b07feedc00a2cd15e` (matches `3f3b7ac96ac3`).
**Method cards (Rule 6):** `refs/methods/` does not exist in this repository; no project method cards govern these disputes. Statistical claims are adjudicated against the cited methodology (`Angelopoulos et al. 2024` in `manuscript/references.bib`), prior round rulings (`scripts/round history`), and direct verification against the live code and benchmark scripts on `HEAD` (`5cb83fb`).

---

## F1 -- C3: Offline constant `λ̂ = 10.0` CRC claim vs. online `n_cal = 15` decoy calibration in `src/mdl/dp.rs`

- **Ruling:** **UPHELD** (diagnosis and prescribed repair both upheld).
- **What survives:** `C3` at [`manuscript/horizon_paper.md:40`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L40) (and the Figure 1 caption at [`manuscript/horizon_paper.md:22`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L22)) is refuted: `horizon_wga` runs an **online** `n_cal = 15` synthetic `64`-bp i.i.d. decoy maximum heuristic (`calibrate_null_crc_for_gc`, `alpha = 0.07`, `λ̂ >= ctx.effective_min_bits`) that violates Angelopoulos et al. (2024) exchangeability, while `horizon_filter` uses `(KARLIN_MDL_COEFF * search_space_bits).max(10.0)` without decoy calibration.
- **Verification & prior rounds:**
  - Verified [`src/mdl/mod.rs:428`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/mod.rs#L428) calls `dp::calibrate_null_crc_for_gc(observed_gc, 15, 0.07, ctx.effective_min_bits, &hmm_params)` online whenever `dp_jobs` is non-empty, estimating `observed_gc` from the first `4,096` bases of `seq_x` and `seq_y`.
  - Verified [`src/mdl/dp.rs:86-165`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/dp.rs#L86-L165) generates `15` i.i.d. `64`-bp decoy pairs (`band_w = 16`) and selects `scores[14]` (when `>= base_floor_bits`) because `(0 + 1)/16 = 0.0625 <= 0.07`; [`src/mdl/dp.rs:820`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/dp.rs#L820) tests `lambda_hat > 12.6`.
  - Verified [`src/bin/horizon_filter.rs:199`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_filter.rs#L199) does not call `calibrate_null_crc_for_gc` and uses `(KARLIN_MDL_COEFF * search_space_bits).max(10.0)`.
  - Round 4 L3 F7 (`NARROWED`) previously established that `n_cal = 15, alpha = 0.07` forces `0/15` empirical errors; the attacker properly disclosed that ruling and identified the post-Round-4 addition of `calibrate_null_crc_for_gc` that left `manuscript/horizon_paper.md:40` (and line `22`) out of sync.
- **Repair & cost:** Endorsed. Updating `manuscript/horizon_paper.md:40` (and the matching `λ̂ = 10.0` mention in Figure 1 at line `22`) to describe the online `n_cal = 15` `64`-bp decoy maximum heuristic, disclose the exchangeability gap, and note `horizon_filter`'s fixed floor requires no re-measurement — **the fix is free**.

---

## F2 -- C2: `SDUST+MDL` `0.0052%` p-arm effect size and `score_paf_vs_hprc.py` gaps

- **Ruling:** **UPHELD** (the attacker's narrowing of `C2` and diagnosis of `score_paf_vs_hprc.py` are upheld).
- **What survives:** On `chr6` (`n = 1`), `--net --query --target` alters only `5,694` bp (`0.0052%` of `109,995,785` bp `aln_len`), all on the `GRCh38.chr6` p-arm (`< 46,000,000` bp) outside the `1-to-1` HPRC truth, leaving conditioned (`99.65%`) and unconditional (`99.32%`) precision unchanged to 2 decimal places; `score_paf_vs_hprc.py` omits unconditional precision, lacks `--truth`, and cannot run on this host without `numpy`.
- **Verification & prior rounds:**
  - Re-ran `./target/release/horizon_filter` on `/tmp/chr6_fastga.paf` (`56,857` records): `--net` kept `565` records (`109,995,785` bp `aln_len`, `elapsed=0.50s`); `--net --query /tmp/hprc_grch38_chr6.fa --target /tmp/hprc_chm13_chr6.fa` kept `553` records (`filtered_mdl=436`, `109,990,091` bp `aln_len`, `elapsed=4.92s`). The exact net difference is `5,694` bp (`0.005177%`), consisting of `15` removed records (`max_qend = 45,643,201`) and `3` added records (`max_qend = 41,436,646`). Zero `1-to-1` segments in `/tmp/hprc_chr6_truth_segments.tsv` have `xs < 46,000,000`.
  - Verified [`scripts/benchmarks/score_paf_vs_hprc.py:165-168`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/score_paf_vs_hprc.py#L165-L168) omits `exact_pairs / total_pairs`, lacks `--truth` in `argparse` (while [`scripts/benchmarks/bench_vgp_filter_locus.py:44,65`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/bench_vgp_filter_locus.py#L44-L65) passes `--truth`), and `python3 -c "import numpy"` fails with `ModuleNotFoundError: No module named 'numpy'` (carrying over Round 6 `R2-RESOLVE-LIFT`).
- **Repair & cost:** Endorsed. Clarifying the `5,694` bp (`0.0052%`) p-arm-only effect in [`docs/BENCHMARK_RESULTS.md:86`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L86) (and [`manuscript/horizon_paper.md:84`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L84)) is a text-only edit (**the fix is free**); adding unconditional precision and `--truth` to `score_paf_vs_hprc.py` is a small script edit whose live re-run remains gated on installing `numpy`.

---

## F3 -- C4: `verify_merge_guide_maf.py` pairwise `hg38`-projected scope on tuning-set `chrM`

- **Ruling:** **UPHELD** (the attacker's narrowing of `C4` is upheld).
- **What survives:** `verify_merge_guide_maf.py` reproduces the four pairwise `hg38`-projected `chrM` accuracy numbers (`n = 1`, tuning set) and HAL header check, but does not verify 5-way column accuracy, HAL ancestral-node (`Anc0`–`Anc3`) sequence fidelity, or the 5-mammal nuclear `horizon_guide` run (`0.859 s`).
- **Verification & prior rounds:**
  - Re-ran `python3 scripts/benchmarks/verify_merge_guide_maf.py` on `HEAD`: `elapsed_s=0.3004` (vs. `0.258 s` in `C4`, `0.3150 s` in attack), `hal_bottom_segs=2103`, `hg38 x panTro4` `0.955/0.998`, `rheMac3` `0.892/0.971`, `mm10` `0.718/0.975`, `canFam3` `0.741/0.982`.
  - Verified [`scripts/benchmarks/verify_merge_guide_maf.py:139-155`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/verify_merge_guide_maf.py#L139-L155) only checks `horizon_wga --hal-info` header metadata and pairwise `(ref_sp, tgt_sp)` projections, and [`docs/BENCHMARK_RESULTS.md:38-39`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L38-L39) classifies `chrM` as a Tuning Set.
  - Round 5 L2 F3 and L7 F5 previously flagged `chrM` `n = 1` sampling; the attacker properly cited those rulings and evaluated the newly added `verify_merge_guide_maf.py` script.
- **Repair & cost:** Endorsed. Qualifying `C4` at [`manuscript/horizon_paper.md:84`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L84) to state the pairwise `hg38`-projected `chrM` scope and unverified nuclear/ancestral accuracy requires no re-measurement — **the fix is free**.

---

## F4 -- C1 and C5: Non-statistical claims outside the statistics lens

- **Ruling:** **UPHELD** (`UNCHALLENGED` under L3 confirmed for both `C1` and `C5`).
- **What survives:** `C1` (Rust unit test and `clippy` counts) and `C5` (`Wittler 2023` bibliography entry) carry no statistical inference claim; both are factually verified on `HEAD`.
- **Verification & prior rounds:**
  - Re-ran `/home/agent/.cargo/bin/cargo test --no-fail-fast --quiet` (`88` library passed + `1` ignored `.hal` test + `34` binary passed = `122` active, `123` total) and `/home/agent/.cargo/bin/cargo clippy --all-targets -- -D warnings` (exited `0`). Consistent with Round 4 L3 F1.
  - Verified `wittler_2023_canonical-kmer-encoding` at [`manuscript/references.bib:481-489`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/references.bib#L481-L489) and [`manuscript/horizon_paper.md:150`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L150).
- **Repair & cost:** None required.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `src/mdl/mod.rs:428` runs online `n_cal = 15` `64`-bp i.i.d. decoy maximum calibration (`λ̂ >= effective_min_bits`), refuting `C3`'s offline constant `λ̂ = 10.0` CRC claim and violating exchangeability; text fix is free |
| F2 | UPHELD | Verified `--net --query --target` alters `5,694` bp (`0.0052%`) on the p-arm (`< 46 Mbp`), leaving `99.65%` cond. / `99.32%` uncond. precision unchanged, while `score_paf_vs_hprc.py` lacks uncond. precision, `--truth`, and `numpy` |
| F3 | UPHELD | Verified `verify_merge_guide_maf.py` only scores pairwise `hg38`-projected accuracy on `~16.5 kb` tuning-set `chrM`, leaving 5-way columns, HAL `Anc0`–`Anc3` sequences, and the `0.859 s` nuclear run unverified; text qualification is free |
| F4 | UPHELD | Verified `C1` (`88` lib + `34` bin = `122` active + `1` ignored `.hal` test, `0` clippy warnings) and `C5` (`Wittler 2023` in `references.bib:481-489`) are accurate non-statistical claims outside L3 |
