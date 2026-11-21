# Round 7 Adjudication — L4 (Provenance)

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `5cb83fb6a1346a8193f9e4048055638d0a422b0a`)
**Attack Sealed Digest:** `89bba6f7b1a8` (verified via `sha256sum reports/rounds/round-7/attack-L4-provenance.md`: `89bba6f7b1a82993c5b88bd0490e805777256cb2a085d90f421566392f96c733`)
**Method Cards:** None exist in this repository (`refs/methods/` is absent), so no method card governs any dispute here.

---

## F1 -- C1: `halTest.hal` integration test passes on `HEAD`, but the fixture has zero repository provenance

- **Adjudication:** **UPHELD** (diagnosis upheld; documentation repair endorsed — **the fix is free**)
- **Prior-round check:** `scripts/round history "halTest" --kind adjudication` confirms no prior mentions across 7 round directories.
- **Verification:** Re-running `/home/agent/.cargo/bin/cargo test --no-fail-fast` on `HEAD` (`5cb83fb`) gives `88 passed; 0 failed; 1 ignored` (lib) plus `34` binary tests (`2+6+2+1+2+4+11+1+1+1+2+0 = 34`, `0` failed), and `cargo clippy --all-targets -- -D warnings` exits `0` with zero warnings. Running `HORIZON_TEST_HAL=/home/rcb112/project/claude/bench/data/iuc-port/tools-iuc/tools/haltools/test-data/halTest.hal /home/agent/.cargo/bin/cargo test --lib hal::tests::real_hal_decode_from_env -- --ignored` passes (`1 passed; 0 failed`). `grep -rn "halTest\.hal"` confirms `manuscript/horizon_paper.md:84` is the sole mention in the repository; `src/hal/tests.rs:240-243` and `docs/BENCHMARK_DATASETS.md` record no upstream URL, version, or checksum, and `scripts/benchmarks/fetch_benchmarks.sh` does not fetch it.
- **What survives:** `C1`'s test counts and zero-warning `clippy` status hold on `HEAD`, but `halTest.hal` has zero repository provenance.
- **Repair:** Endorsed to document the upstream source/URL and SHA-256 of `halTest.hal` in `src/hal/tests.rs:241` and `docs/BENCHMARK_DATASETS.md`; **the fix is free** (text-only, no re-measurement).

---

## F2 -- C2: `bench_vgp_filter_locus.py` crashes on `--truth` and relies on untracked `/tmp/` inputs with no `FastGA` provenance

- **Adjudication:** **UPHELD** (diagnosis upheld; script/documentation repairs endorsed — **the fix is free** for `score_paf_vs_hprc.py`)
- **Prior-round check:** `round 6` (`adjudication-L7.md`, `adjudication-completeness.md`) previously noted missing `numpy` and untracked `/tmp/` inputs on `chr6`; `F2` newly documents the broken `--truth` interface between `bench_vgp_filter_locus.py` and `score_paf_vs_hprc.py` and the missing unconditional precision output.
- **Verification:** Re-running `/home/rcb112/project/claude/bench/.venv-galaxy/bin/python3 scripts/benchmarks/score_paf_vs_hprc.py --paf /tmp/chr6_fastga_mdl_net.paf --all-records` reproduces `primary_records=553`, `total_pairs=109866176`, `pairs_on_hprc_1to1_x=109503252`, `exact_coord_pairs=109119624`, `exact_coord_recall_on_hprc_1to1=0.9912487...` (`99.12%`), `exact_coord_precision_on_shared=0.9964966...` (`99.65%`), and `109119624 / 109866176 = 0.9932049...` (`99.32%`). Inspection confirms:
  1. `scripts/benchmarks/bench_vgp_filter_locus.py:44,65` passes `--truth str(args.truth)` to `score_paf_vs_hprc.py`, whereas `scripts/benchmarks/score_paf_vs_hprc.py:9-16` defines no `--truth` flag and hardcodes `TRUTH_TSV = "/tmp/hprc_chr6_truth_segments.tsv"`.
  2. `score_paf_vs_hprc.py:151-169` does not compute or print unconditional exact precision (`99.32%`), and system `python3 -c "import numpy"` fails with `ModuleNotFoundError: No module named 'numpy'`.
  3. Inputs default to `/tmp/` (`bench_vgp_filter_locus.py:16-19`, `score_paf_vs_hprc.py:16`), and `docs/BENCHMARK_DATASETS.md:50-52` and `scripts/benchmarks/fetch_benchmarks.sh` omit the `FastGA` command, version, and checksum for `chr6_fastga.paf`.
- **What survives:** The `C2` numbers reproduce on `/tmp/chr6_fastga_mdl_net.paf`, but the tracked driver `bench_vgp_filter_locus.py` crashes on `--truth`, `score_paf_vs_hprc.py` omits unconditional precision, and `chr6_fastga.paf` lacks generation provenance.
- **Repair:** Endorsed to add `--truth` and unconditional exact precision output to `scripts/benchmarks/score_paf_vs_hprc.py` (**the fix is free** — no re-measurement needed) and to document the `FastGA` command/version/checksum in `docs/BENCHMARK_DATASETS.md` (free if already recorded externally, otherwise requires re-running `FastGA` to log exact version/flags).

---

## F3 -- C3: No offline calibration provenance exists for `n_cal = 15`, `α = 0.07`, `λ̂ = 10.0` bits, and `calibrate_null_crc_for_gc` runs online

- **Adjudication:** **UPHELD** (diagnosis upheld; wording repair endorsed — **the fix is free**)
- **Prior-round check:** `round 4` (`L3-F7`, `adjudication-L3.md:62-67`), `round 5` (`L3-F6`, `L5-F7`), and `round 6` (`L3-F6`, `L5-F6`) already flagged `min_delta_bits = 10.0` as hardcoded and `n_cal = 15` as underpowered. The attacker properly disclosed `round 4`'s `L3-F7` and `round 5`/`round 6`. What is new in `C3` is the manuscript's explicit provenance claim (`manuscript/horizon_paper.md:40`) that these defaults were *"calibrated offline at risk level `α = 0.07` (`n_cal = 15`; Angelopoulos et al. 2024)"*.
- **Verification:** `src/mdl/mod.rs:73` hardcodes `min_delta_bits: 10.0` and `src/bocpd.rs:10` hardcodes `MIN_POST_PROB: f64 = 0.15`. `dp::calibrate_null_crc_for_gc` (`src/mdl/dp.rs:120-165`, called at `src/mdl/mod.rs:428` with `n_cal = 15, alpha = 0.07, base_floor_bits = ctx.effective_min_bits`) runs **online** at alignment time over 15 synthetic 64-bp xorshift decoy windows at `observed_gc` and clamps `>= ctx.effective_min_bits`. No offline calibration dataset, script, or log exists in the repository.
- **What survives:** `C3`'s claim that `λ̂ = 10.0` bits (`min_delta_bits = 10.0`) was calibrated offline at `α = 0.07` (`n_cal = 15`) is false; `10.0` and `0.15` are hardcoded defaults and `calibrate_null_crc_for_gc` runs online on synthetic decoy draws.
- **Repair:** Endorsed to replace the "calibrated offline" claim in `manuscript/horizon_paper.md:40` with an accurate description of the hardcoded defaults (`min_delta_bits = 10.0`, `MIN_POST_PROB = 0.15`) and the online 15-draw synthetic decoy clamp (`calibrate_null_crc_for_gc`); **the fix is free**.

---

## F4 -- C4: Commit `2896229` lacks progressive ancestral-node HAL emission (`HEAD`-only), and the `0.859 s` 5-mammal run lacks tracked inputs/scripts

- **Adjudication:** **UPHELD** (diagnosis upheld; manuscript pinning repair endorsed — **the fix is free**)
- **Prior-round check:** `round 5` and `round 6` (`R2-HMM-XDROP`) carried the requirement for progressive multi-contig `.hal` emission, implemented in `HEAD` (`5cb83fb`).
- **Verification:** Re-running `python3 scripts/benchmarks/verify_merge_guide_maf.py` on `HEAD` (`5cb83fb`) reproduces `hal_tree=((((hg38,panTro4)Anc0,rheMac3)Anc1,mm10)Anc2,canFam3)Anc3;`, `hal_bottom_segs=2103`, and the exact MAF recall/precision values (`0.955/0.998`, `0.892/0.971`, `0.718/0.975`, `0.741/0.982`) in `0.2997 s`. Inspection confirms:
  1. `git show --stat 5cb83fb` first added ancestral-node progressive HAL emission to `src/bin/horizon_guide/main.rs` and `src/hal/mod.rs`; `2896229` did not have it, and `manuscript/horizon_paper.md:10` still describes `2896229`'s star-tree HAL stub (`4` bottom segments) while `manuscript/horizon_paper.md:84` conflates `2896229 / HEAD`.
  2. `scripts/benchmarks/verify_merge_guide_maf.py:109` defaults `--maf` to `ROOT.parent / "pg" / "data" / "benchmarks" / "chrM.maf"`, an external sibling-repo path not fetched by `scripts/benchmarks/fetch_benchmarks.sh`.
  3. `grep -rn "0\.859" scripts/ docs/ manuscript/` shows `0.859 s` appears only in `manuscript/horizon_paper.md:10,84`, with no tracked script or input FASTAs in `scripts/benchmarks/` or `docs/BENCHMARK_DATASETS.md`.
- **What survives:** The `verify_merge_guide_maf.py` progressive HAL result reproduces on `HEAD` (`5cb83fb`), but `C4` conflates `2896229` with `HEAD`, `verify_merge_guide_maf.py` defaults to an external `../pg/` path, and the `0.859 s` 5-mammal run has no tracked script/inputs.
- **Repair:** Endorsed to pin the progressive ancestral-node HAL result in `manuscript/horizon_paper.md:84` to `HEAD` (`5cb83fb`) rather than `2896229 / HEAD` and update `manuscript/horizon_paper.md:10` consistently (**the fix is free**), and to document the `chrM.maf` fetch path and `0.859 s` 5-mammal inputs in `docs/BENCHMARK_DATASETS.md`.

---

## F5 -- C5: `Wittler R (2023)` bibliography entry matches `refs/catalog_index.json` and `manuscript/references.bib`

- **Adjudication:** **UPHELD** (verification confirmed; no repair needed)
- **Prior-round check:** `round 5` (`L6-F2`) and `round 6` (`L6-F4`, `L6-F5`) required adding `wittler_2023_canonical-kmer-encoding` (`R5-BIB`).
- **Verification:** `manuscript/horizon_paper.md:150` matches `refs/catalog_index.json:1430-1433` (`slug: wittler_2023_canonical-kmer-encoding`, `doi: 10.24072/pcjournal.323`) and `manuscript/references.bib:481-489`, and is cited in-text at `manuscript/horizon_paper.md:28`. `../library` is absent on disk, so full-text PDF verification was not possible.
- **What survives:** `C5` is verified against `refs/catalog_index.json` and `manuscript/references.bib`.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `122` active tests + `1` `.hal` test + `0` clippy warnings reproduce on `HEAD`, and `halTest.hal` has zero repository provenance (free fix to document URL/SHA-256). |
| F2 | UPHELD | `C2` numbers reproduce on `/tmp/chr6_fastga_mdl_net.paf`, but `bench_vgp_filter_locus.py` crashes passing `--truth` to `score_paf_vs_hprc.py`, `score_paf_vs_hprc.py` omits unconditional precision, and `chr6_fastga.paf` lacks `FastGA` provenance. |
| F3 | UPHELD | No offline calibration dataset, script, or log exists for `λ̂ = 10.0` / `n_cal = 15` / `α = 0.07`; `10.0` and `0.15` are hardcoded and `calibrate_null_crc_for_gc` runs online (free wording fix). |
| F4 | UPHELD | `verify_merge_guide_maf.py` reproduces the `2,103`-seg progressive HAL on `HEAD` (`5cb83fb`), but `C4` conflates `2896229` with `HEAD` (free fix), defaults `--maf` to `../pg/`, and lacks a script/inputs for the `0.859 s` 5-mammal run. |
| F5 | UPHELD | `Wittler R (2023)` matches `refs/catalog_index.json` (`10.24072/pcjournal.323`) and `manuscript/references.bib:481-489` and is cited at `manuscript/horizon_paper.md:28`. |
