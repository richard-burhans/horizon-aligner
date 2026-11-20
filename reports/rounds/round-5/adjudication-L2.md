# Round 5 ADJUDICATION of sampling -- Round 5: Full chr6 (8cd5db0/92e345c), HPRC & minimap2 Verification, B10K 363-way, and Technical Companion

**Attack evaluated:** [`reports/rounds/round-5/attack-sampling.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-5/attack-sampling.md) (sealed digest `cbdb4f838186` verified via `sha256sum`).
**Prior rounds & method cards:** Checked `scripts/round history "sampling" --kind adjudication`. Prior rounds already adjudicated `CHM13.chr6` q-arm-only sampling (`r4 F1`), MAF-block concatenation in `extract_fasta_from_maf.py` (`r2 F8`, `r3 F1`, `r4 F4/F6`), and uncalibrated `n_cal = 15` (`r4 R1-HMM`). No method cards exist in `refs/methods/` (`NO_LIBRARY`).

---

## F1 -- C1: `10,000`-segment HPRC `chr22` subgraph is `96.09%` `N`-bases (`10.51 Mbp` in one segment) and mixes in `chr14` / `chrUn` contigs

- **Verdict:** **UPHELD**
- **What survives:** `C1` is narrowed to an `n = 1` prefix-sampled `10,000`-segment subgraph (`427,960` bp non-`N` ACGT + `10,510,000` bp `'N'`-gap in segment `76067108`) that mixes `GRCh38.chr22` (`3,991` segments) with `GRCh38.chrUn_GL000224v1` (`1,978`), `GRCh38.chr22_KI270736v1_random` (`1,919`), and `GRCh38.chr14_KI270722v1_random` (`1,230`).
- **Premise check:** Direct inspection of [`/tmp/hprc_chr22_sub_real.gfa`](file:///tmp/hprc_chr22_sub_real.gfa) confirms every count: `10,000` `S`-lines, `10,937,960` total bp, `10,510,000` `'N'` bases (`96.09%`) in segment `76067108`, `427,960` non-`N` bp, `1` `P`-line (`2,000` steps), and `4` `W`-lines (`18,000` bases each).
- **Repair:** Endorsed on both [`horizon_paper.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md) and [`docs/BENCHMARK_RESULTS.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md); **the fix is free** (text-only scope disclosure, no re-measurement needed).

---

## F2 -- C2: HPRC `chr6` precision (`48.63%`) conditions on `pairs_on_hprc_1to1_x`, excluding `60.43M` (`21.4%`) off-target pairs

- **Verdict:** **UPHELD**
- **What survives:** `C2` is narrowed because `48.63%` exact pair precision (`107,836,669 / 221,734,653`) conditions on `pairs_on_hprc_1to1_x`, omitting `60,428,393` (`21.4%`) off-target pairs (including `14,584,569` unique `GRCh38.chr6` bases outside `hprc_1to1_x` aligned onto `CHM13.chr6` q-arm), so unconditional exact pair precision across all `282,163,046` emitted pairs is `38.22%`.
- **Premise check:** Verified directly in [`manuscript/data/chr6_hprc_comparison.json`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/data/chr6_hprc_comparison.json), [`/tmp/compare_chr6_horizon_vs_hprc.py`](file:///tmp/compare_chr6_horizon_vs_hprc.py#L136-L148), and [`manuscript/data/chr6_minimap2_summary.json`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/data/chr6_minimap2_summary.json) (`minimap2 -x asm5` has `171,202` [`0.16%`] off-target pairs and `99.60%` unconditional exact precision).
- **Repair:** Endorsed to report `38.22%` unconditional exact pair precision alongside `48.63%` conditioned on `pairs_on_hprc_1to1_x`, plus the `60.43M` (`21.4%`) off-target pair count on `n = 1` chromosome pair; **the fix is free** (all numbers are already in `chr6_hprc_comparison.json`).

---

## F3 -- C3: `horizon_merge` samples `~16.5 kb` `chrM` (`92/100` species) and `horizon_guide` samples `11–33 KB` MAF-concatenated `chr22_random` slices

- **Verdict:** **UPHELD**
- **What survives:** `C3` is narrowed because `horizon_merge` operates on `~16.5 kb` haploid `chrM` (`92/100` species) and `horizon_guide` aligns `11–33 KB` MAF-concatenated `chr22_KI270735v1_random` extracts (`4` HAL bottom segments, `n = 1` run) rather than full nuclear chromosomes.
- **Premise check:** Verified in [`scripts/benchmarks/run_all_evals.py`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/run_all_evals.py#L68-L75), [`scripts/benchmarks/extract_fasta_from_maf.py`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/extract_fasta_from_maf.py#L41-L45), and the extracted FASTAs in `../pg/data/benchmarks/fastas/` (`cactus447way_Equus_caballus.fa` `11K`, `cactus447way_Macaca_mulatta.fa` `22K`, `multiz470way_HLcalPym1.fa` `33K`, `chrM.fa` `17K`).
- **Repair:** Endorsed to disclose the `~16.5 kb` `chrM` (`92/100` species) and `11–33 KB` MAF-concatenated `chr22_KI270735v1_random` (`4` HAL bottom segments) scope in [`horizon_paper.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md); **the fix is free** (text-only clarification).

---

## F4 -- C4: `C4` omits MAF-block concatenation (`11–387 KB` slices incl. `B10K`), cherry-picks `18%`-GC (`0` pairs) over `30%`/`45%`-GC (`288`/`149` pairs) with no script, and cites unimplemented `n_cal = 15`

- **Verdict:** **UPHELD**
- **What survives:** `C4` is narrowed: (1) `horizon_paper.md` omits MAF-block concatenation (`B10K` extracts `/tmp/b10k_galGal.fa` `386,668` bp and `/tmp/b10k_taeGut.fa` `349,383` bp from `500 kb`, plus `11–33 KB` `chr22_KI270735v1_random` extracts); (2) `C4` (`horizon_paper.md:10`) cites only `18%`-GC (`0` pairs) while omitting `30%`-GC (`288` pairs) and `45%`-GC (`149` pairs), and no script/test in `src/` or `scripts/` generates these `10,000`-bp random GC nulls; (3) `n_cal = 15` (`horizon_paper.md:40`, `horizon_companion.md:153`) is unimplemented (`min_delta_bits = 10.0` is hardcoded in [`src/mdl/mod.rs:71`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/mod.rs#L71)).
- **Premise check:** Verified `/tmp/b10k_galGal.fa` (`386,668` bp) and `/tmp/b10k_taeGut.fa` (`349,383` bp), confirmed `grep` across `src/` and `scripts/` returns zero matches for `n_cal` and `18%`, and confirmed `min_delta_bits: 10.0` in `src/mdl/mod.rs:71`.
- **Repair:** Split cost: disclosing MAF-block concatenation (`11–387 KB`), reporting all three GC null bins (`18%`: `0`, `30%`: `288`, `45%`: `149` pairs) in `C4`, and removing `n_cal = 15` / `α = 0.07` (documenting hardcoded `min_delta_bits = 10.0`) are **free** text edits; adding a tracked random-GC null script/replicates or implementing Conformal Risk Control (`n_cal = 15`) requires code and re-measurement (not free).

---

## F5 -- C5: `BandGrid` flat buffer indexing is a deterministic layout identity with no sampling component

- **Verdict:** **UPHELD**
- **What survives:** `C5` survives under the sampling lens because `BandLayout::slot` ([`src/mdl/dp.rs:160-172`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/dp.rs#L160-L172)) is a deterministic array-indexing identity (`i * self.stride + col` with `stride = 2 * bw + 1`, `col = (j - center(i)) + bw`) with no statistical sampling component.
- **Premise check:** Inspected `src/mdl/dp.rs:140-172`; confirmed exact formula and absence of sampling.

---

## F6 -- C6: `# EG:Z:` tag emission is a format claim outside the sampling lens

- **Verdict:** **UPHELD**
- **What survives:** `C6` remains unchallenged under the sampling lens because `# EG:Z:` tag emission ([`src/formats.rs:153-158`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/formats.rs#L153-L158), [`src/bin/horizon_resolve.rs:476`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_resolve.rs#L476)) is a structural format specification, while the empirical `10k`-subgraph sampling defect in `horizon_resolve` is already captured by `F1`.
- **Premise check:** Inspected `src/formats.rs:153-158` and `src/bin/horizon_resolve.rs:476`; confirmed format-only behavior.

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `/tmp/hprc_chr22_sub_real.gfa` is `96.09%` `N`-bases (`10,510,000` bp in 1 seg, `427,960` bp ACGT) and mixes `chrUn`/`chr14` contigs (`n = 1`); text fix is free |
| F2 | UPHELD | `48.63%` conditions on `pairs_on_hprc_1to1_x`, dropping `60,428,393` (`21.4%`) off-target pairs (`38.22%` unconditional across `282.16M` pairs); text fix is free |
| F3 | UPHELD | `horizon_merge` uses `~16.5 kb` `chrM` (`92/100` species) and `horizon_guide` uses `11–33 KB` MAF-concatenated `chr22_KI270735v1_random` slices (`4` HAL bottom segs); text fix is free |
| F4 | UPHELD | `horizon_paper.md` omits MAF concatenation (`11–387 KB` incl. `B10K`), omits `30%`/`45%`-GC nulls in `C4` with no tracked script, and cites unimplemented `n_cal = 15`; text fixes are free, script/implementation is not |
| F5 | UPHELD | `BandLayout::slot` (`src/mdl/dp.rs:160-172`) is a deterministic layout identity with no sampling component |
| F6 | UPHELD | `# EG:Z:` emission (`src/formats.rs:153-158`, `src/bin/horizon_resolve.rs:476`) is a format claim; `10k`-subgraph sampling is covered in `F1` |
