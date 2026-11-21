# Round 7 Adjudication — L7 (Reproduction)

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` at `HEAD` (`5cb83fb6a1346a8193f9e4048055638d0a422b0a`)  
**Attack file verified:** `reports/rounds/round-7/attack-L7-reproduction.md` (`sha256: 9ce90679b8a3...`, matches sealed digest)  
**Method cards (`refs/methods/`):** None exist in this project; no statistical method-card dispute applies.

---

## F1 — C1 (`122` active tests, `0` clippy warnings, `halTest.hal` external & `blanchette.hal` panic): **UPHELD**

- **What survives:** `122` active Rust unit tests (`88` library + `34` binary) and `0` `clippy` warnings reproduce on `HEAD`, while the `.hal` integration test (`real_hal_decode_from_env`) relies on an external `halTest.hal` fixture not bundled in the repository and panics on `blanchette.hal`.
- **Verification:**
  - `/home/agent/.cargo/bin/cargo test --no-fail-fast` passes `88` library tests (`1` ignored) and `34` binary tests (`1+2+6+2+1+2+4+11+1+1+1+2 = 34`, `122` active, `123` total); `/home/agent/.cargo/bin/cargo clippy --all-targets -- -D warnings` exits `0`.
  - `find . -name "halTest.hal"` returns nothing; `/home/rcb112/project/claude/bench/data/iuc-port/tools-iuc/tools/haltools/test-data/halTest.hal` exists externally and passes `HORIZON_TEST_HAL=... cargo test --lib real_hal_decode_from_env -- --ignored`.
  - `HORIZON_TEST_HAL=/tmp/ext_tools/src/hal/phyloP/test/blanchette.hal cargo test --lib real_hal_decode_from_env -- --ignored` panics at `src/hal/tests.rs:267:10` (`SEQNAME_ARRAY`).
- **Ruling on repair:** Endorsed. Clarifying in `manuscript/horizon_paper.md:84` that `halTest.hal` is an external `tools-iuc` fixture and that `real_hal_decode_from_env` targets `halTest.hal`-layout HDF5 v0 HAL files is accurate; **the fix is free** (prose-only, no re-measurement).

---

## F2 — C2 (`553` records / `99.12%` recall / `99.65%` cond. / `99.32%` uncond. precision reproduce directly, but `bench_vgp_filter_locus.py` crashes): **UPHELD**

- **What survives:** The `553`-record (`436` MDL-filtered) `chr6` metrics reproduce directly via `horizon_filter` + `score_paf_vs_hprc.py`, but the tracked driver `scripts/benchmarks/bench_vgp_filter_locus.py:44,65` crashes with exit code `2` because it passes `--truth` to `score_paf_vs_hprc.py`, which does not accept `--truth`.
- **Prior rounds:** Round 6 (`adjudication-L4.md:18`, `adjudication-completeness.md:33-38`) verified `553` records and updated `848` -> `553` on `42700e1`, but did not catch the broken `--truth` invocation in `bench_vgp_filter_locus.py` (added in `4fd2240`).
- **Verification:**
  - `/tmp/ext_tools/pdf_venv/bin/python scripts/benchmarks/bench_vgp_filter_locus.py` fails with `CalledProcessError` (exit status `2`): `score_paf_vs_hprc.py: error: unrecognized arguments: --truth /tmp/hprc_chr6_truth_segments.tsv`. (Under system `python3`, it also lacks `numpy`.)
  - Direct execution of `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --net --query /tmp/hprc_grch38_chr6.fa --target /tmp/hprc_chm13_chr6.fa` yields `total=56857 kept=553 filtered_mdl=436 filtered_len=0 elapsed=4.82s`, and `score_paf_vs_hprc.py --all-records` yields `primary_records: 553`, `total_pairs: 109866176`, `pairs_on_hprc_1to1_x: 109503252`, `exact_coord_pairs: 109119624`, `exact_coord_recall_on_hprc_1to1: 0.9912487...` (`99.12%`), `exact_coord_precision_on_shared: 0.9964966...` (`99.65%`), and unconditional exact precision `109119624 / 109866176 = 0.9932049...` (`99.32%`).
- **Ruling on repair:** Endorsed. Adding `--truth` (defaulting to `"/tmp/hprc_chr6_truth_segments.tsv"`) to `scripts/benchmarks/score_paf_vs_hprc.py:9-16` or removing `--truth` from `scripts/benchmarks/bench_vgp_filter_locus.py:44,65` resolves the crash; **the fix is free** (no re-measurement).

---

## F3 — C3 (Online per-job CRC calibration above `max(10.0, 0.55 * log2(N * M))` refutes offline constant `λ̂ = 10.0` bits): **UPHELD**

- **What survives:** `C3` (`manuscript/horizon_paper.md:40`) is false as stated: Horizon calibrates CRC **online** per job at `observed_gc` (`n_cal = 15`, `α = 0.07`) above a search-space-scaled base floor `max(10.0, 0.55 * log2(N * M))`, rather than calibrating offline to a constant `λ̂ = 10.0` bits.
- **Verification:**
  - `src/mdl/mod.rs:151`: `let effective_min_bits = cfg.min_delta_bits.max(KARLIN_MDL_COEFF * search_space_bits);` (with `cfg.min_delta_bits = 10.0` at `:73`).
  - `src/mdl/mod.rs:428`: calls `dp::calibrate_null_crc_for_gc(observed_gc, 15, 0.07, ctx.effective_min_bits, &hmm_params)`, which generates 15 64-bp decoy windows online at `observed_gc` and computes `λ̂` via `calibrate_crc_threshold` (`src/mdl/dp.rs:86-165`).
  - `src/bocpd.rs:10`: `MIN_POST_PROB = 0.15`.
  - **Whole-set check (Rule 3):** `manuscript/horizon_paper.md:22` (Figure 1 caption) also writes `(`λ̂ = 10.0` bits)` and must be updated alongside line `40`.
- **Ruling on repair:** Endorsed, extended to cover `manuscript/horizon_paper.md:22` as well as `:40`. **The fix is free** (prose-only, no re-measurement).

---

## F4 — C4 (`verify_merge_guide_maf.py` reproduces `chrM` accuracy and `2,103` bottom segs, but 5-mammal `horizon_guide` `1,109`/`967` blocks does not reproduce on `HEAD`): **UPHELD**

- **What survives:** `verify_merge_guide_maf.py` reproduces `2,103` bottom segments and all `chrM` recall/precision figures on `HEAD` (`0.3656 s` here vs. `0.258 s` reported), whereas the 5-mammal nuclear `1,109`-block / `967`-multi-species / `0.859 s` run from `2896229` does not reproduce on `HEAD` (`5cb83fb`).
- **Verification:**
  - `git show --stat 5cb83fb` confirms `src/bin/horizon_guide/main.rs` (`+120/-17`) and `src/hal/mod.rs` (`+123`) changed HAL emission from the star-tree stub (`4` bottom segments) to progressive ancestral-node HAL (`root_members=13`).
  - `python3 scripts/benchmarks/verify_merge_guide_maf.py` on `HEAD` (`5cb83fb`) outputs `hal_tree=((((hg38,panTro4)Anc0,rheMac3)Anc1,mm10)Anc2,canFam3)Anc3;`, `hal_bottom_segs=2103`, `total_blocks=266`, `multi_species_blocks=255`, `elapsed_s=0.3656`, and exact recall/precision `0.955/0.998` (`panTro4`), `0.892/0.971` (`rheMac3`), `0.718/0.975` (`mm10`), `0.741/0.982` (`canFam3`).
  - Running `./target/release/horizon_guide --progressive --hal --rbest` on the 5-mammal FASTAs on `HEAD` (`5cb83fb`) with `1` thread yields `hal_bottom_segs=2505 total_blocks=657 multi_species_blocks=586 elapsed_s=0.8546` (vs. attacker's `2650`/`707`/`633`/`1.0780`; neither matches `1,109`/`967` or `4` bottom segments).
- **Ruling on repair:** Endorsed. Removing the `2896229 / HEAD` conflation in `manuscript/horizon_paper.md:84` and pinning `2896229` (or updating `docs/BENCHMARK_RESULTS.md:11,58` and `manuscript/horizon_paper.md:84` to a tracked `HEAD` run) plus updating `0.258 s` -> `0.36 s` for `verify_merge_guide_maf.py` is appropriate; **the text-pinning fix is free** (adding a new tracked 5-mammal nuclear benchmark script is a small script addition).

---

## F5 — C5 (`Wittler R (2023)` bibliographic entry): **UPHELD**

- **What survives:** `Wittler R (2023)` matches `refs/catalog_index.json:1430-1432` and `manuscript/references.bib:481-489`; `volume = {3}` and `pages = {e87}` cannot be independently cross-checked against `../library` because `../library` is absent in this workspace, leaving C5 unchallenged.
- **Prior rounds:** Rounds 5–6 (`adjudication-L6.md`) endorsed adding `wittler_2023_canonical-kmer-encoding`.
- **Verification:** Confirmed in `refs/catalog_index.json:1430-1432` (`doi: 10.24072/pcjournal.323`) and `manuscript/references.bib:481-489` (`volume = {3}`, `pages = {e87}`).

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `122` active tests and `0` clippy warnings reproduce on `HEAD`; `halTest.hal` is external (`tools-iuc`) and `real_hal_decode_from_env` panics on `blanchette.hal`; prose clarification is free. |
| F2 | UPHELD | `553` records and `99.12%` / `99.65%` / `99.32%` metrics reproduce directly, but `scripts/benchmarks/bench_vgp_filter_locus.py` crashes (exit `2`) passing `--truth` to `score_paf_vs_hprc.py`; CLI fix is free. |
| F3 | UPHELD | `src/mdl/mod.rs:151,428` and `src/mdl/dp.rs:86-165` calibrate CRC online per job above `max(10.0, 0.55 * log2(N * M))`, refuting offline constant `λ̂ = 10.0` bits (update both `manuscript/horizon_paper.md:22` and `:40`; fix is free). |
| F4 | UPHELD | `verify_merge_guide_maf.py` reproduces `chrM` accuracy and `2,103` bottom segs (`~0.36 s`), while the 5-mammal `1,109`/`967`-block run from `2896229` does not reproduce on `HEAD` (`5cb83fb`); pinning `2896229` in text is free. |
| F5 | UPHELD | `Wittler R (2023)` matches `refs/catalog_index.json:1430-1432` and `manuscript/references.bib:481-489`; `3:e87` is unverified without `../library`, so C5 remains unchallenged. |
