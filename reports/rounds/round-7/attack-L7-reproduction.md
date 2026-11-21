# Round 7 Adversarial Attack — L7 (Reproduction)

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` at `HEAD` (`5cb83fb6a1346a8193f9e4048055638d0a422b0a`)  
**Lens:** Reproduction (`BRIEF-reproduction.md`)  
**Environment notes (Rules 2, 3, 4):**
- **Held literature (Rule 2):** `../library` is absent (`ls: cannot access '../library': No such file or directory`), so literature verification worked from `refs/catalog_index.json` and `manuscript/references.bib` alone.
- **Method cards (Rule 3):** This project keeps no method cards (`BRIEF-reproduction.md:142`).
- **Prior rounds (Rule 4):** Checked `scripts/round history` for `Wittler`, `122`, `553`, and `halTest`. Round 6 adjudicated `553` records (`adjudication-L4.md:18`, `adjudication-completeness.md:33-38`, updating `848` -> `553` on `42700e1`), and Rounds 5–6 adjudicated `wittler_2023_canonical-kmer-encoding` addition (`adjudication-L6.md`). No prior round mentioned `halTest`.

---

## F1 -- C1: `cargo test` and `cargo clippy` reproduce `122` active tests and zero warnings, but `halTest.hal` is external and `real_hal_decode_from_env` panics on `blanchette.hal`

- **Claim:** `C1` (`manuscript/horizon_paper.md:84`) — *"All `122` active Rust unit tests (`88` library + `34` binary, plus `1` `.hal` integration test verified against `halTest.hal` via `HORIZON_TEST_HAL`; `123` total) pass on `HEAD` with zero `clippy` warnings."*
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. Running `cargo test --no-fail-fast` on `HEAD` (`5cb83fb`) reproduces `88` passing library tests (`src/lib.rs`, `1 ignored`) and `34` passing binary tests (`1 + 2 + 6 + 2 + 1 + 2 + 4 + 11 + 1 + 1 + 1 + 2 = 34`, total `122` active), and `cargo clippy --all-targets -- -D warnings` exits `0` with zero warnings (`cargo fmt --check` also exits `0`).
  2. However, `halTest.hal` is **not** tracked in `/home/rcb112/project/gemini/bench/horizon-aligner` (`find . /tmp -name "halTest.hal"` returns nothing); it exists only in sibling external checkouts such as `/home/rcb112/project/claude/bench/data/iuc-port/tools-iuc/tools/haltools/test-data/halTest.hal`. Running `HORIZON_TEST_HAL=/home/rcb112/project/claude/bench/data/iuc-port/tools-iuc/tools/haltools/test-data/halTest.hal cargo test real_hal_decode_from_env -- --ignored` passes (`1 passed`).
  3. More importantly, `real_hal_decode_from_env` (`src/hal/tests.rs:240-280`) **panics** on another real external HAL fixture present in `/tmp/ext_tools/src/hal/phyloP/test/blanchette.hal`:
     ```
     thread 'hal::tests::real_hal_decode_from_env' panicked at src/hal/tests.rs:267:10:
     SEQNAME_ARRAY
     ```
     (It passes on `/tmp/ext_tools/src/hal/paf/tests/input/mr.hal`, `/tmp/external_chrM.hal`, and `/tmp/guide_5sp.hal`.) Thus, the `.hal` integration test is neither self-contained in the repository nor general across real HAL fixtures in `/tmp/ext_tools/`.
- **Suggested correction:** Note in `manuscript/horizon_paper.md:84` that `halTest.hal` is an external `tools-iuc` fixture (not bundled in the repository) and that `real_hal_decode_from_env` targets `halTest.hal`-layout HDF5 v0 HAL files.

---

## F2 -- C2: Underlying `horizon_filter` metrics (`553` records, `99.12%` recall, `99.65%`/`99.32%` precision) reproduce directly, but the tracked driver `scripts/benchmarks/bench_vgp_filter_locus.py` crashes with exit code 2

- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:86`) — *"`553` records (`436` filtered by MDL; removes `12` noisy/non-1-to-1 records while preserving `100%` of exact `109,119,624` bp coordinates: **99.12%** exact recall, **99.65%** cond. / **99.32%** uncond. exact precision)"*
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. Executing the tracked benchmark driver `scripts/benchmarks/bench_vgp_filter_locus.py` on `HEAD` **fails immediately** after `horizon_filter --net` with `CalledProcessError` (exit status 2):
     ```
     usage: score_paf_vs_hprc.py [-h] [--paf PAF] [--out-npz OUT_NPZ]
                                 [--out-json OUT_JSON] [--all-records]
     score_paf_vs_hprc.py: error: unrecognized arguments: --truth /tmp/hprc_chr6_truth_segments.tsv
     ```
     At `scripts/benchmarks/bench_vgp_filter_locus.py:44` and `:65`, the driver passes `--truth str(args.truth)` to `scripts/benchmarks/score_paf_vs_hprc.py`, which has no `--truth` CLI option (`score_paf_vs_hprc.py:9-16` hardcodes `TRUTH_TSV = "/tmp/hprc_chr6_truth_segments.tsv"`).
  2. When invoked directly without `--truth` using the untracked `/tmp/` inputs (`/tmp/chr6_fastga.paf`, `/tmp/hprc_grch38_chr6.fa`, `/tmp/hprc_chm13_chr6.fa`, `/tmp/hprc_chr6_truth_segments.tsv`), the numbers **do** reproduce on `HEAD`:
     - `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --net --query /tmp/hprc_grch38_chr6.fa --target /tmp/hprc_chm13_chr6.fa --out /tmp/r7_chr6_mdl_net.paf` -> `[horizon_filter] total=56857 kept=553 filtered_mdl=436 filtered_len=0 elapsed=4.86s`.
     - `/tmp/ext_tools/pdf_venv/bin/python scripts/benchmarks/score_paf_vs_hprc.py --paf /tmp/r7_chr6_mdl_net.paf --all-records` -> `primary_records: 553`, `total_pairs: 109866176`, `pairs_on_hprc_1to1_x: 109503252`, `exact_coord_pairs: 109119624`, `exact_coord_recall_on_hprc_1to1: 0.9912487297083773` (`99.12%`), `exact_coord_precision_on_shared: 0.9964966519898423` (`99.65%`), and unconditional exact precision `109119624 / 109866176 = 0.9932049...` (`99.32%`).
     (Note: Round 6 `adjudication-completeness.md:33-38` verified the `553` count via `horizon_filter`, but did not catch the broken `--truth` call in `bench_vgp_filter_locus.py`, which was added in `4fd2240`.)
- **Suggested correction:** Fix `scripts/benchmarks/bench_vgp_filter_locus.py:44,65` (or add `--truth` to `scripts/benchmarks/score_paf_vs_hprc.py:9-16`) so the tracked reproduction script runs without crashing.

---

## F3 -- C3: Horizon does NOT enforce a constant `λ̂ = 10.0` bits calibrated offline; CRC runs online per job atop `max(10.0, 0.55 * log2(N * M))`

- **Claim:** `C3` (`manuscript/horizon_paper.md:40`) — *"Using fixed Conformal Risk Control defaults calibrated offline at risk level `α = 0.07` (`n_cal = 15`; Angelopoulos et al. 2024), Horizon enforces a constant bit-savings floor `λ̂ = 10.0` bits (`min_delta_bits = 10.0`) and a minimum column posterior floor (`0.15`)."*
- **Verdict:** **REFUTED**
- **Concrete flaw:** The actual implementation in `src/mdl/mod.rs` and `src/mdl/dp.rs` contradicts both *"calibrated offline"* and *"enforces a constant bit-savings floor `λ̂ = 10.0` bits"*:
  1. **Search-space scaling (`src/mdl/mod.rs:151`):** `effective_min_bits` is computed as:
     ```rust
     let effective_min_bits = cfg.min_delta_bits.max(KARLIN_MDL_COEFF * search_space_bits);
     ```
     where `KARLIN_MDL_COEFF = 0.55` and `search_space_bits = log2(N * M)`. For any sequence pair with `N * M > 2^(10.0 / 0.55) ≈ 297,500` (e.g., `1000 × 1000` bp gives `10.96` bits; `chr6` `170.8M × 110.7M` gives `29.74` bits), the base floor is strictly greater than `10.0` bits, not constant `10.0`.
  2. **Online per-job CRC calibration (`src/mdl/mod.rs:428`, `src/mdl/dp.rs:120-165`):** When `dp_jobs` is non-empty, `src/mdl/mod.rs:428` calls:
     ```rust
     dp::calibrate_null_crc_for_gc(observed_gc, 15, 0.07, ctx.effective_min_bits, &hmm_params)
     ```
     which generates `n_cal = 15` 64-bp decoy windows **online** at `observed_gc` under `hmm_params` and sets `λ̂ = calibrate_crc_threshold(&scores, 0.07, ctx.effective_min_bits)` (`src/mdl/dp.rs:86-114`). Thus `λ̂` is neither calibrated offline nor fixed at a constant `10.0` bits. (Only `MIN_POST_PROB = 0.15` in `src/bocpd.rs:10` and `cfg.min_delta_bits = 10.0` in `src/mdl/mod.rs:73` are fixed defaults.)
- **Suggested correction:** Rewrite `manuscript/horizon_paper.md:40` to state accurately: *"Using fixed Conformal Risk Control hyperparameters (`α = 0.07`, `n_cal = 15`; Angelopoulos et al. 2024) calibrated online per alignment job at the observed GC composition above a base floor `max(min_delta_bits = 10.0, 0.55 * log2(N * M))` bits, Horizon enforces a minimum column posterior floor (`0.15`)."*

---

## F4 -- C4: 5-mammal `horizon_guide` (`0.859 s`, `1,109` blocks, `967` multi-sp) does NOT reproduce on `HEAD` (`5cb83fb`); `verify_merge_guide_maf.py` reproduces accuracy (`0.358 s`)

- **Claim:** `C4` (`manuscript/horizon_paper.md:84`) — *"and `horizon_guide --progressive --hal --rbest` (`2896229` / `HEAD`) aligned `5` mammals (`hg38`, `Macaca_mulatta`, `HLcalPym1`, `Equus_caballus`, `Canis_lupus`) in `0.859 s` (`1,109` MAF blocks, `967` multi-species) and on 5-mammal UCSC `chrM.maf` (`verify_merge_guide_maf.py`, `0.258 s`) emitted a progressive ancestral-node HDF5 `v0` HAL (`((((hg38,panTro4)Anc0,rheMac3)Anc1,mm10)Anc2,canFam3)Anc3;`, `2,103` bottom segments) with `95.5%/99.8%` (`panTro4`), `89.2%/97.1%` (`rheMac3`), `71.8%/97.5%` (`mm10`), and `74.1%/98.2%` (`canFam3`) MAF recall/precision."*
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. Commit `5cb83fb` (`HEAD`) modified `src/bin/horizon_guide/main.rs` (`+120/-17` lines) and `src/hal/mod.rs` (`+123` lines) to replace the star-tree HAL stub (`4` bottom segments, `root_members=5`) with progressive ancestral-node HAL emission (`root_members=13`). Consequently, `2896229` and `HEAD` do not execute the same `horizon_guide --progressive --hal` code.
  2. Re-running the 5-mammal (`hg38`, `Macaca_mulatta`, `HLcalPym1`, `Equus_caballus`, `Canis_lupus`) command on `HEAD` (`5cb83fb`) using `../pg/data/benchmarks/fastas/cactus447way_hg38.fa`, `cactus447way_Macaca_mulatta.fa`, `multiz470way_HLcalPym1.fa`, `cactus447way_Equus_caballus.fa`, and `/tmp/test_bench_all/d_t.fa` produces:
     ```
     hinge=hg38 species=5 hal_tree=((((hg38,Macaca_mulatta)Anc0,HLcalPym1)Anc1,Equus_caballus)Anc2,Canis_lupus)Anc3; hal_bottom_segs=2650 total_blocks=707 multi_species_blocks=633 elapsed_s=1.0780
     ```
     Neither the `0.859 s` runtime nor the `1,109` MAF blocks / `967` multi-species block counts reproduce on `HEAD` (and no tracked script exists for this 5-mammal run).
  3. By contrast, `python3 scripts/benchmarks/verify_merge_guide_maf.py` **does** reproduce on `HEAD` (`5cb83fb`), yielding `hal_tree=((((hg38,panTro4)Anc0,rheMac3)Anc1,mm10)Anc2,canFam3)Anc3;`, `hal_bottom_segs=2103`, `total_blocks=266`, `multi_species_blocks=255`, `elapsed_s=0.3582` (vs. `0.258 s` reported), and exact MAF recall/precision `0.955/0.998` (`panTro4`), `0.892/0.971` (`rheMac3`), `0.718/0.975` (`mm10`), and `0.741/0.982` (`canFam3`). Note also that `halValidate` is absent from `/tmp/ext_tools/env/bin/`, so external HAL schema validity of the emitted `.hal` file remains unverified by `halValidate`.
- **Suggested correction:** Remove `2896229 / HEAD` conflation in `manuscript/horizon_paper.md:84` and `docs/BENCHMARK_RESULTS.md:11,58`: either re-run and track the 5-mammal nuclear benchmark on `HEAD` (`707` blocks / `633` multi-sp / `2,650` HAL bottom segments / `1.078 s` with the exact FASTA paths) or pin the `0.859 s` / `1,109`-block figure strictly to `2896229` (star-tree stub), and update `0.258 s` -> `0.358 s` for `verify_merge_guide_maf.py` if desired.

---

## F5 -- C5: `Wittler R (2023)` matches `refs/catalog_index.json`, though `3:e87` cannot be independently verified without `../library`

- **Claim:** `C5` (`manuscript/horizon_paper.md:150`) — *"- **Wittler R** (2023). General encoding of canonical *k*-mers. *Peer Community Journal*, 3:e87."*
- **Verdict:** **UNCHALLENGED**
- **Evidence checked:** `refs/catalog_index.json:1430-1432` lists slug `wittler_2023_canonical-kmer-encoding`, `"cite": "Wittler, R. (2023). *General encoding of canonical k-mers*. Peer Community Journal."`, and `"doi": "10.24072/pcjournal.323"`. `manuscript/references.bib:481-489` contains the matching `@article{wittler_2023_canonical-kmer-encoding}` entry with `volume = {3}` and `pages = {e87}`. Because `../library` is unavailable in this workspace, `volume = {3}` and `pages = {e87}` cannot be independently cross-checked against the held PDF/summary, and no concrete defect was found.

---

## Verdict Table

| claim_id | verdict | summary |
| :--- | :--- | :--- |
| C1 | NARROWED | `122` active tests (`88` lib + `34` bin) and `0` clippy warnings reproduce on `HEAD`; `halTest.hal` is external (`tools-iuc`) and `real_hal_decode_from_env` panics on `blanchette.hal`. |
| C2 | NARROWED | Direct `horizon_filter` + `score_paf_vs_hprc.py` reproduces `553` records (`436` MDL), `99.12%` recall, `99.65%` cond. / `99.32%` uncond. precision, but `scripts/benchmarks/bench_vgp_filter_locus.py` crashes (exit 2) passing `--truth`. |
| C3 | REFUTED | `src/mdl/mod.rs:151,428` and `src/mdl/dp.rs:120-165` calibrate CRC **online** per job above `max(10.0, 0.55 * log2(N * M))`, not offline to a constant `10.0`-bit floor. |
| C4 | NARROWED | `verify_merge_guide_maf.py` reproduces `2,103` bottom segs and all `chrM` recall/prec (`0.358 s` vs `0.258 s`), but 5-mammal `horizon_guide` (`0.859 s`, `1,109`/`967` blocks) does not reproduce on `HEAD` (`5cb83fb` gives `1.078 s`, `707`/`633` blocks, `2,650` bottom segs). |
| C5 | UNCHALLENGED | Author, year, title, and journal match `refs/catalog_index.json:1430-1432`; `3:e87` cannot be verified without `../library`. |
