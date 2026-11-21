# Round 7 Adversarial Attack — L4 (Provenance)

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `5cb83fb6a1346a8193f9e4048055638d0a422b0a`)
**Literature & Method-Card Scope (Rules 2–3):**
- `../library` does not exist on disk (`ls: cannot access '../library': No such file or directory`), so all literature checks were performed against `refs/catalog_index.json` (pinned to `refs/LIBRARY_VERSION` = `d4df9122eddd84694b2c3215e5478c46f4314044`) and `manuscript/references.bib` alone.
- Per `BRIEF-provenance.md`, this project keeps **no method cards** yet.

---

## F1 -- C1: `halTest.hal` integration test passes on `HEAD`, but the fixture has zero repository provenance

- **Claim:** `C1` (`manuscript/horizon_paper.md:84`)
- **Verdict:** `NARROWED`
- **Prior-round check:** `scripts/round history "halTest"` returned no prior mentions.
- **Concrete flaw:** The `122` active unit tests (`88` library + `34` binary) and `0` `clippy` warnings hold on `HEAD` (`5cb83fb`), and the `1` ignored `.hal` test (`real_hal_decode_from_env`, `src/hal/tests.rs:240-280`) passes when pointed to an external `halTest.hal` via `HORIZON_TEST_HAL`. However, `halTest.hal` is not shipped in `horizon-aligner`, is not fetched by `scripts/benchmarks/fetch_benchmarks.sh`, and its origin, URL, version, and checksum are undocumented anywhere in `src/hal/tests.rs:240-243`, `docs/BENCHMARK_DATASETS.md`, or `manuscript/horizon_paper.md:84` (it had to be located externally at `/home/rcb112/project/claude/bench/data/iuc-port/tools-iuc/tools/haltools/test-data/halTest.hal`).
- **Evidence:**
  - `cargo test --no-fail-fast`: `88 passed; 0 failed; 1 ignored` (lib) + `34` binary tests (`1+2+6+2+1+2+4+11+1+1+1+2+0 = 34`), `0` failed.
  - `cargo clippy --all-targets -- -D warnings`: exited `0` with zero warnings.
  - `HORIZON_TEST_HAL=/home/rcb112/project/claude/bench/data/iuc-port/tools-iuc/tools/haltools/test-data/halTest.hal cargo test --lib hal::tests::real_hal_decode_from_env -- --ignored`: `1 passed; 0 failed`.
  - `grep -rn "halTest\.hal" README.md docs/ src/ tests/ scripts/ manuscript/`: only `manuscript/horizon_paper.md:84` mentions `halTest.hal`.
- **Correction:** Document the exact upstream source/URL or SHA-256 of `halTest.hal` (e.g., from `tools-iuc/tools/haltools/test-data/halTest.hal`) in `src/hal/tests.rs:241` and `docs/BENCHMARK_DATASETS.md`, or vendor/fetch it in `scripts/benchmarks/fetch_benchmarks.sh`.

---

## F2 -- C2: `bench_vgp_filter_locus.py` crashes on `--truth` and relies on untracked `/tmp/` inputs with no `FastGA` provenance

- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:86`)
- **Verdict:** `NARROWED`
- **Prior-round check:** `round 6`'s `R2-RESOLVE-LIFT` asked to install `numpy` and re-score the `565`/`553` `chr6` outputs against 2D HPRC truth.
- **Concrete flaw:** Running `score_paf_vs_hprc.py --paf /tmp/chr6_fastga_mdl_net.paf --all-records` (using a `venv` with `numpy`, since system `/usr/bin/python3` lacks `numpy`) reproduces `primary_records: 553`, `exact_coord_pairs: 109119624`, `exact_coord_recall_on_hprc_1to1: 0.9912487...` (`99.12%`), `exact_coord_precision_on_shared: 0.9964966...` (`99.65%`), and `109119624 / 109866176 = 0.9932049...` (`99.32%`). However, the tracked provenance chain in `scripts/benchmarks/` is broken in three ways:
  1. `scripts/benchmarks/bench_vgp_filter_locus.py:44,65` invokes `score_paf_vs_hprc.py` with `--truth str(args.truth)`, but `scripts/benchmarks/score_paf_vs_hprc.py:9-14` has **no `--truth` argument** (it hardcodes `TRUTH_TSV = "/tmp/hprc_chr6_truth_segments.tsv"` at line 16), so the tracked driver crashes with `unrecognized arguments: --truth` whenever `args.truth` exists.
  2. `score_paf_vs_hprc.py:151-169` never computes or prints unconditional exact precision (`99.32%`), and `/usr/bin/python3` fails on `import numpy as np` (`ModuleNotFoundError: No module named 'numpy'`).
  3. All inputs (`chr6_fastga.paf`, `hprc_grch38_chr6.fa`, `hprc_chm13_chr6.fa`, `hprc_chr6_truth_segments.tsv`) are hardcoded to `/tmp/` (`bench_vgp_filter_locus.py:16-19`, `score_paf_vs_hprc.py:16`), and neither `docs/BENCHMARK_DATASETS.md:50-52` nor `scripts/benchmarks/fetch_benchmarks.sh` records the `FastGA` command, version, or checksum used to generate `chr6_fastga.paf` (`56,857` records).
- **Evidence:**
  - `scripts/benchmarks/bench_vgp_filter_locus.py:44,65` vs. `scripts/benchmarks/score_paf_vs_hprc.py:9-16`.
  - `/home/rcb112/project/claude/bench/.venv-galaxy/bin/python3 scripts/benchmarks/score_paf_vs_hprc.py --paf /tmp/chr6_fastga_mdl_net.paf --all-records`: `primary_records=553`, `total_pairs=109866176`, `pairs_on_hprc_1to1_x=109503252`, `exact_coord_pairs=109119624`.
- **Correction:** Add `--truth` support and unconditional precision output to `scripts/benchmarks/score_paf_vs_hprc.py`, and document the exact `FastGA` command/version/checksum for `chr6_fastga.paf` in `docs/BENCHMARK_DATASETS.md`.

---

## F3 -- C3: No offline calibration provenance exists for `n_cal = 15`, `α = 0.07`, `λ̂ = 10.0` bits, and `calibrate_null_crc_for_gc` runs online

- **Claim:** `C3` (`manuscript/horizon_paper.md:40`)
- **Verdict:** `REFUTED`
- **Prior-round check:** `round 4`'s `L3-F7` (`reports/rounds/round-4/adjudication-L3.md:62-67`) was adjudicated `NARROWED` on the mathematical CRC bound (`(0 + 1) / 16 <= 0.07` when `0/15` errors occur, while noting `n_cal = 15` is underpowered and `min_delta_bits = 10.0` is hardcoded). `round 5` (`L3-F6`, `L5-F7`) and `round 6` (`L3-F6`, `L5-F6`) also flagged this as an uncalibrated constant.
- **Concrete flaw:** From a provenance standpoint, the claim that *"Using fixed Conformal Risk Control defaults calibrated offline at risk level `α = 0.07` (`n_cal = 15`; Angelopoulos et al. 2024), Horizon enforces a constant bit-savings floor `λ̂ = 10.0` bits (`min_delta_bits = 10.0`) and a minimum column posterior floor (`0.15`)"* is false:
  1. There is **no offline calibration dataset, script, or log** anywhere in the repository showing that `λ̂ = 10.0` bits was ever calibrated offline at `α = 0.07` with `n_cal = 15`.
  2. `min_delta_bits: 10.0` (`src/mdl/mod.rs:73`) and `MIN_POST_PROB: 0.15` (`src/bocpd.rs:10`) are hardcoded constants.
  3. The only function using `n_cal = 15` and `alpha = 0.07` (`dp::calibrate_null_crc_for_gc` called at `src/mdl/mod.rs:428`, defined at `src/mdl/dp.rs:120-165`) runs **online** at alignment time over 15 synthetic 64-bp xorshift decoy windows at `observed_gc`, and it clamps its return value `>= ctx.effective_min_bits` (`max(10.0, 0.55 * log2(N*M))`), so it never produces or validates the `10.0`-bit floor offline.
- **Evidence:**
  - `src/mdl/mod.rs:73` (`min_delta_bits: 10.0`), `src/bocpd.rs:10` (`MIN_POST_PROB: f64 = 0.15`), `src/mdl/mod.rs:401,428` (`calibrate_null_crc_for_gc(observed_gc, 15, 0.07, ctx.effective_min_bits, &hmm_params)`), `src/mdl/dp.rs:120-165` (`win_len = 64`, xorshift PRNG).
- **Correction:** Replace the "calibrated offline" claim in `manuscript/horizon_paper.md:40` with an accurate description of the hardcoded defaults (`min_delta_bits = 10.0`, `MIN_POST_PROB = 0.15`) and the online 15-draw synthetic decoy clamp (`calibrate_null_crc_for_gc`), or commit an actual offline CRC calibration script and dataset.

---

## F4 -- C4: Commit `2896229` lacks progressive ancestral-node HAL emission (`HEAD`-only), and the `0.859 s` 5-mammal run lacks tracked inputs/scripts

- **Claim:** `C4` (`manuscript/horizon_paper.md:84`)
- **Verdict:** `NARROWED`
- **Prior-round check:** `round 5` and `round 6` (`R2-HMM-XDROP`) carried the requirement for multi-contig progressive `.hal` emission; `5cb83fb` (`HEAD`) implemented it.
- **Concrete flaw:** Running `python3 scripts/benchmarks/verify_merge_guide_maf.py` on `HEAD` (`5cb83fb`) reproduces the exact progressive tree `((((hg38,panTro4)Anc0,rheMac3)Anc1,mm10)Anc2,canFam3)Anc3;`, `hal_bottom_segs=2103`, and the exact MAF recall/precision values (`95.5%/99.8%`, `89.2%/97.1%`, `71.8%/97.5%`, `74.1%/98.2%`) in `0.2701 s`. However, C4 has three provenance defects:
  1. `manuscript/horizon_paper.md:84` pins `horizon_guide --progressive --hal --rbest` to `(`2896229` / `HEAD`)`, conflating `2896229` (which only had the 4-bottom-segment star-tree HAL stub, as still stated in the Abstract at `manuscript/horizon_paper.md:10`) with `HEAD` (`5cb83fb`, which first implemented ancestral-node progressive HAL emission).
  2. `scripts/benchmarks/verify_merge_guide_maf.py:109` defaults `--maf` to `ROOT.parent / "pg" / "data" / "benchmarks" / "chrM.maf"`, an external sibling-repo path outside `horizon-aligner` and not fetched by `scripts/benchmarks/fetch_benchmarks.sh`.
  3. The `0.859 s` (`1,109` MAF blocks, `967` multi-species) 5-mammal run (`hg38`, `Macaca_mulatta`, `HLcalPym1`, `Equus_caballus`, `Canis_lupus`) has no tracked script or input FASTAs in `scripts/benchmarks/` or `docs/BENCHMARK_DATASETS.md`.
- **Evidence:**
  - `git show --stat 5cb83fb` vs. `git show 2896229:src/bin/horizon_guide/main.rs`.
  - `manuscript/horizon_paper.md:10` (still describes `2896229` star-tree HAL stub with `4` bottom segments) vs. `manuscript/horizon_paper.md:84`.
  - `python3 scripts/benchmarks/verify_merge_guide_maf.py` output on `HEAD` (`5cb83fb`).
- **Correction:** Pin the progressive ancestral-node HAL result to `HEAD` (`5cb83fb`) rather than `2896229 / HEAD`, update `manuscript/horizon_paper.md:10` to match `5cb83fb`, and document/track the `0.859 s` 5-mammal inputs and `chrM.maf` fetch path.

---

## F5 -- C5: `Wittler R (2023)` bibliography entry matches `refs/catalog_index.json` and `manuscript/references.bib`

- **Claim:** `C5` (`manuscript/horizon_paper.md:150`)
- **Verdict:** `SURVIVED`
- **Prior-round check:** `round 5` (`L6-F2`) and `round 6` (`L6-F4`, `L6-F5`) required adding `wittler_2023_canonical-kmer-encoding` to `manuscript/horizon_paper.md` and `manuscript/references.bib` (`R5-BIB`).
- **Verification:** `manuscript/horizon_paper.md:150` (`- **Wittler R** (2023). General encoding of canonical *k*-mers. *Peer Community Journal*, 3:e87.`) matches `refs/catalog_index.json` (`slug: wittler_2023_canonical-kmer-encoding`, `doi: 10.24072/pcjournal.323`) and `manuscript/references.bib:481-489` (`author = {Wittler, Roland}`, `volume = {3}`, `pages = {e87}`, `year = {2023}`, `doi = {10.24072/pcjournal.323}`), and is cited in-text at `manuscript/horizon_paper.md:28`. (Full-text PDF verification was not possible because `../library` is absent on disk.)

---

## Verdict Table

| claim_id | verdict | summary |
| :--- | :--- | :--- |
| C1 | NARROWED | `122` active tests (`88` lib + `34` bin), `1` `.hal` test, and `0` clippy warnings pass on `HEAD` (`5cb83fb`), but `halTest.hal` has no repo-tracked source, URL, or checksum. |
| C2 | NARROWED | `553` records / `109,119,624` bp / `99.12%` recall / `99.65%` cond. / `99.32%` uncond. precision reproduce on `/tmp/chr6_fastga_mdl_net.paf`, but `bench_vgp_filter_locus.py` crashes passing `--truth` to `score_paf_vs_hprc.py` and `chr6_fastga.paf` lacks generation provenance. |
| C3 | REFUTED | No offline calibration dataset, script, or log exists for `λ̂ = 10.0` bits / `n_cal = 15` / `α = 0.07`; `10.0` and `0.15` are hardcoded and `calibrate_null_crc_for_gc` runs online (disclosing `round 4`'s `L3-F7`). |
| C4 | NARROWED | `verify_merge_guide_maf.py` reproduces the `2,103`-seg progressive HAL and recall/precision on `HEAD` (`5cb83fb`), but C4 conflates `2896229` (star-tree stub) with `HEAD`, uses external `../pg/` `chrM.maf`, and lacks a script for the `0.859 s` 5-mammal run. |
| C5 | SURVIVED | Matches `refs/catalog_index.json` (`10.24072/pcjournal.323`) and `manuscript/references.bib:481-489`, with in-text citation at `manuscript/horizon_paper.md:28` (`../library` absent on disk). |
