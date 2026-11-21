# Round 7 Adversarial Adjudication -- L2 (Sampling)

**Sealed attack digest verified:** `sha256sum reports/rounds/round-7/attack-L2-sampling.md` = `e73d057ce8ba912f6da3afa9c1e32757b826e7e1fa907200d5b7eb7f76a533bf` (matches `e73d057ce8ba`).

**Method cards & prior rounds (Rules 5–6):**
- This repository has no `refs/methods/` directory, so no method card governs these disputes.
- Queried `scripts/round history` (`--kind adjudication`) for `n_cal`, `553`, `chrM`, and `halTest`:
  - `n_cal`: Round 4 (`adjudication-L3.md`, `F7`) and Round 5 (`adjudication-L3.md`, `F6`) already flagged `n_cal = 15` (`α = 0.07`, `λ̂ = 10.0` bits) as underpowered; `4fd2240` / `5cb83fb` subsequently wired `calibrate_null_crc_for_gc` into `src/mdl/mod.rs:401-429`, leaving `manuscript/horizon_paper.md:40` stale.
  - `chrM` / `horizon_guide`: Round 5 (`adjudication-L2.md`, `F3`/`F4`; `adjudication-L7.md`, `F5`) already adjudicated MAF concatenation and `chrM` sampling; `5cb83fb` added `scripts/benchmarks/verify_merge_guide_maf.py` and the `5`-mammal `chrM.maf` progressive HAL run.
  - `halTest`: No prior round adjudicated `halTest.hal`.

---

## F1 -- `C1`: `123` total tests rely on a single `55 KB` external `.hal` fixture (`HORIZON_TEST_HAL`) and synthetic/toy unit fixtures

- **Verdict:** **UPHELD**
- **What survives:** `C1`'s test and lint counts are accurate (`88` library + `34` binary = `122` active passing tests, `1` `#[ignore]`d `.hal` integration test), but `C1` requires qualification because the `1` `.hal` integration test (`src/hal/tests.rs:240-280`, `real_hal_decode_from_env`) is skipped under default `cargo test` and exercises only `n = 1` external `55 KB` 4-genome synthetic fixture (`halTest.hal`), while none of the `123` tests validate `horizon_guide`'s emitted `.hal` output or a chromosome-scale real `.hal` file.
- **Premise check:**
  - Verified `/home/agent/.cargo/bin/cargo test --quiet` passes `88` library + `34` binary (`122` active, `1` ignored).
  - Verified `/home/rcb112/project/claude/bench/data/iuc-port/tools-iuc/tools/haltools/test-data/halTest.hal` is `55K` and `src/hal/tests.rs:240-280` is `#[ignore = "requires HORIZON_TEST_HAL=/path/to/file.hal"]`.
  - Checked `src/hal/tests.rs:961-994` (`progressive_multi_contig_hal_builds_ancestral_tree_and_roundtrips`): it tests `HalAlignment::from_multi_contig_progressive_blocks` on synthetic `8`-bp contigs, not the `horizon_guide` binary or real chromosome-scale `.hal` files.
- **Repair:** Endorsed (`NARROWED` qualification on `C1` in `manuscript/horizon_paper.md:84`). **The fix is free** (text qualification only; no re-measurement).

## F2 -- `C2`: `chr6` `horizon_filter` scoring (`100%` preservation of `109,119,624` bp) evaluates `n = 1` q-arm-only `1-to-1` truth (`60.23M–170.74M`), where all `15` removed records lie in the unmeasured p-arm (`0–45.7 Mbp`)

- **Verdict:** **UPHELD**
- **What survives:** `C2` must be narrowed to state that `100%` coordinate preservation applies only to the `1-to-1` `GRCh38.chr6` q-arm truth subset (`60.23M–170.74M`, `n = 1` chromosome pair, `n = 1` aligner `FastGA`), because all `15` MDL-removed netted intervals lie in the unmeasured p-arm (`< 45.7 Mbp`).
- **Premise check:**
  - Verified `/tmp/hprc_chr6_truth_segments.tsv` spans `GRCh38.chr6:60,229,947..170,743,762` (`1,628,018` `mult == "1"` segments = `110,082,990` bp; `38` `mult != "1"` segments = `310` bp), and `scripts/benchmarks/score_paf_vs_hprc.py:34` filters to `mult == "1"`.
  - Verified coordinate diff between `/tmp/chr6_fastga_net_head.paf` (`565` records) and `/tmp/chr6_fastga_mdl_net_head.paf` (`553` records): `15` removed intervals (`6,608` bp total) and `3` added/modified intervals (`914` bp total), with all `15` removed intervals having `qstart` in `147,635..45,643,002` (`< 60,229,947`).
  - Verified `python3 scripts/benchmarks/score_paf_vs_hprc.py --help` fails on this host with `ModuleNotFoundError: No module named 'numpy'`.
- **Repair:** Endorsed (narrow `C2` in `docs/BENCHMARK_RESULTS.md:86` and `manuscript/horizon_paper.md:84`). **The fix is free** (text qualification only; no re-measurement).

## F3 -- `C3`: `horizon_paper.md:40` claims an offline constant `λ̂ = 10.0` bits, whereas `src/mdl/mod.rs:401–428` runs online `calibrate_null_crc_for_gc` with `n_cal = 15` 64-bp i.i.d. draws from a `4,096`-bp prefix sample that is `100%` `'N'` on `GRCh38.chr6`

- **Verdict:** **UPHELD**
- **What survives:** `C3` is refuted against `HEAD` because `src/mdl/mod.rs:401-429` runs online `dp::calibrate_null_crc_for_gc(observed_gc, 15, 0.07, ctx.effective_min_bits, &hmm_params)` (with `n_cal = 15` `64`-bp i.i.d. draws forcing `n_ge = 0` / sample maximum, and a `4,096`-bp prefix GC sampler that reads `4,096` `'N'`s and `0` ACGT bases on `GRCh38.chr6`) rather than enforcing a constant offline `λ̂ = 10.0` bits.
- **Premise check:**
  - Verified `manuscript/horizon_paper.md:40` states Horizon enforces a constant bit-savings floor `λ̂ = 10.0` bits (`min_delta_bits = 10.0`) calibrated offline.
  - Verified `src/mdl/mod.rs:407-429` samples `ctx.seq_x.iter().take(4096).chain(ctx.seq_y.iter().take(4096))` and calls `dp::calibrate_null_crc_for_gc(observed_gc, 15, 0.07, ctx.effective_min_bits, &hmm_params)`, where `effective_min_bits = cfg.min_delta_bits.max(KARLIN_MDL_COEFF * search_space_bits)` (`src/mdl/mod.rs:151`).
  - Verified `/tmp/hprc_grch38_chr6.fa` first `4,096` bp are `4,096` `'N'`s and `0` ACGT bases, while `/tmp/hprc_chm13_chr6.fa` first `4,096` bp are `0` `'N'`s and `4,096` ACGT bases.
  - Verified `src/mdl/dp.rs:86-165` uses `win_len = 64` and `(n_ge + 1) / (15 + 1) <= 0.07`, forcing `n_ge = 0`.
- **Repair:** Split cost:
  1. Updating `manuscript/horizon_paper.md:40` to accurately describe `HEAD`'s online `calibrate_null_crc_for_gc` (`n_cal = 15`, `win_len = 64`, `4,096`-bp prefix GC sample, `effective_min_bits` floor) — **endorsed; the fix is free**.
  2. Replacing the `4,096`-bp prefix GC sampler with chromosome-wide ACGT sampling and increasing `n_cal` beyond the `n_ge = 0` boundary — **endorsed as a code fix, NOT free** (requires modifying `src/mdl/mod.rs` / `src/mdl/dp.rs` and re-running affected `horizon_wga` / `horizon_filter` benchmarks).

## F4 -- `C4`: Progressive HAL accuracy is evaluated only on `n = 1` `~16.5 kb` mitochondrial `chrM.maf` (`15–16 kb` per pair) with `hg38` recall dropping to `71.8%` (`mm10`) and `74.1%` (`canFam3`), while the `5`-mammal nuclear run (`0.859 s`) has no accuracy verification

- **Verdict:** **UPHELD**
- **What survives:** `C4` must be narrowed to state explicitly that progressive HAL accuracy is verified only on `n = 1` single-contig `~16.5 kb` `chrM.maf` (where `hg38` recall drops to `71.8%` on `mm10` and `74.1%` on `canFam3`), whereas the `0.859 s` `5`-mammal nuclear run (`1,109` MAF blocks) has no ground-truth recall/precision verification.
- **Premise check:**
  - Verified by running `python3 scripts/benchmarks/verify_merge_guide_maf.py`: `hg38 x panTro4` `recall=0.955 prec=0.998` (`truth=15980`), `hg38 x rheMac3` `recall=0.892 prec=0.971` (`truth=15025`), `hg38 x mm10` `recall=0.718 prec=0.975` (`truth=15700`), `hg38 x canFam3` `recall=0.741 prec=0.982` (`truth=15575`), `hal_bottom_segs=2103`.
  - Verified no script in `scripts/` measures MAF recall/precision for the `0.859 s` `5`-mammal nuclear run (`hg38`, `Macaca_mulatta`, `HLcalPym1`, `Equus_caballus`, `Canis_lupus`).
- **Repair:** Endorsed (narrow `C4` in `manuscript/horizon_paper.md:84`). **The fix is free** (text qualification only; no re-measurement).

## F5 -- `C5`: Pure bibliographic citation entry (`Wittler R (2023)`) verified in `references.bib` and `catalog_index.json`

- **Verdict:** **UPHELD**
- **What survives:** `C5` remains `UNCHALLENGED` under the sampling lens because it is a verified bibliographic citation (`manuscript/references.bib:481-489`, `refs/catalog_index.json:1430-1433`, `doi:10.24072/pcjournal.323`) making no sampling claim.

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `122+1` test counts hold, but the `1` `.hal` integration test is `#[ignore]`d and exercises only `n = 1` external `55 KB` `halTest.hal`, not `horizon_guide` HAL output or chromosome-scale `.hal` files; free text qualification endorsed. |
| F2 | UPHELD | `100%` coordinate preservation holds only on the `1-to-1` `GRCh38.chr6` q-arm truth (`60.23M–170.74M`); all `15` MDL-removed intervals (`6,608` bp) lie in the unmeasured p-arm (`< 45.7 Mbp`); free text qualification endorsed. |
| F3 | UPHELD | `horizon_paper.md:40` falsely claims offline constant `λ̂ = 10.0` bits, whereas `HEAD` (`src/mdl/mod.rs:407-429`) runs online `n_cal = 15` 64-bp i.i.d. CRC with a `4,096`-bp prefix GC sampler (`100%` `'N'` on `GRCh38.chr6`); manuscript sync is free, code/n_cal fix requires re-benchmarking. |
| F4 | UPHELD | Progressive HAL accuracy is verified only on `n = 1` `~16.5 kb` `chrM.maf` (recall drops to `71.8%`/`74.1%` on `mm10`/`canFam3`), while the `0.859 s` `5`-mammal nuclear run has no accuracy check; free text qualification endorsed. |
| F5 | UPHELD | `C5` (`Wittler R 2023`, `doi:10.24072/pcjournal.323`) is verified in `references.bib` and `catalog_index.json` and makes no sampling claim (`UNCHALLENGED`). |
