# Round 7 Adversarial Attack — L5 (Inference)

**Environment & Provenance Notes (Rules 2, 3, and 4):**
- **Literature library (Rule 2):** `../library` does not exist in this workspace (`ls: cannot access '../library': No such file or directory`), while `refs/LIBRARY_VERSION` pins `d4df9122eddd84694b2c3215e5478c46f4314044`. Accordingly, all prior-work checks below are grounded strictly in `refs/catalog_index.json`, `manuscript/references.bib`, and the repository's own source/docs; no claims are made from memory about held full texts.
- **Method cards (Rule 3):** This project keeps **no method cards** yet (per `BRIEF-inference.md`). Statistical and inferential claims are checked directly against the codebase, benchmark scripts, and cited definitions.
- **Prior rounds (Rule 4):** Checked `scripts/round history` (`--kind adjudication`) for `"conformal"`, `"horizon_filter"`, `"horizon_guide"`, and `"unit tests"`. Where a prior round adjudicated a related point (e.g., R6 L1/L2/L7 on the `118` vs `120` unit-test count on `42700e1` and `horizon_filter` `565`/`553` records; R5 L2 F4 on random-GC null reporting), the prior ruling and what changed on `HEAD` (`5cb83fb6a1346a8193f9e4048055638d0a422b0a`) are explicitly noted.

---

## F1 — C1 (`manuscript/horizon_paper.md:84`): `1` of `123` tests (`real_hal_decode_from_env`) is `#[ignore]`d on `HEAD` without an untracked external `HORIZON_TEST_HAL` fixture

- **Claim:** C1 (*"All `122` active Rust unit tests (`88` library + `34` binary, plus `1` `.hal` integration test verified against `halTest.hal` via `HORIZON_TEST_HAL`; `123` total) pass on `HEAD` with zero `clippy` warnings."*)
- **Verdict:** **NARROWED**
- **Prior-round context:** Round 6 (L1, L2, L7) adjudicated the `118` vs `120` active unit-test count on `42700e1`. On `HEAD` (`5cb83fb6a1346a8193f9e4048055638d0a422b0a`), `2` additional library tests were added, bringing the active count to `122`.
- **Evidence:**
  - Running `cargo test --all-features --quiet` on `HEAD` (`5cb83fb6`) executes `89` library tests (`88 passed; 0 failed; 1 ignored`) plus `34` binary tests (`1+2+6+2+1+2+4+11+1+1+1+2 = 34`, all passed), totaling `122` active passing tests and `1` ignored test.
  - Running `cargo clippy --all-targets --all-features -- -D warnings` on `HEAD` passes with `0` warnings.
  - However, `src/hal/tests.rs:240-243` marks `real_hal_decode_from_env` with `#[ignore = "requires HORIZON_TEST_HAL=/path/to/file.hal"]`, and `halTest.hal` is neither tracked in the repository nor present at `/tmp/halTest.hal` (`ls: cannot access '/tmp/halTest.hal': No such file or directory`).
- **Inferential flaw & correction:** Inferring from `cargo test` on `HEAD` that all `123` tests pass on `HEAD` overstates repository-self-contained verification: only `122` tests run and pass on `HEAD` without an external untracked `halTest.hal` file. Narrow C1 to state that `122` active unit tests pass on `HEAD` while the `1` `.hal` integration test is `#[ignore]`d by default unless an external `HORIZON_TEST_HAL` fixture is provided.

---

## F2 — C2 (`docs/BENCHMARK_RESULTS.md:86`): SDUST+MDL removes only `0.0047%` of pairs over `--net` alone at `~10×` runtime, and `score_paf_vs_hprc.py` neither runs (`numpy` missing) nor computes unconditional precision

- **Claim:** C2 (*"`553` records (`436` filtered by MDL; removes `12` noisy/non-1-to-1 records while preserving `100%` of exact `109,119,624` bp coordinates: **99.12%** exact recall, **99.65%** cond. / **99.32%** uncond. exact precision)"*)
- **Verdict:** **NARROWED**
- **Prior-round context:** Round 6 (L1 F2, L2) narrowed `horizon_filter --net` from `1,004` to `565` records on `42700e1` and left `CARRIED.md` item `R2-RESOLVE-LIFT` (*"Install numpy to re-score 42700e1 horizon_filter 565/553 chr6 outputs against 2D HPRC truth"*).
- **Evidence:**
  - Re-running `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --net --out /tmp/r7_net.paf` vs. `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --net --query /tmp/hprc_grch38_chr6.fa --target /tmp/hprc_chm13_chr6.fa --out /tmp/r7_mdl_net.paf` on `HEAD` (`5cb83fb6`) produces:
    - `--net` alone: `total=56857 kept=565 filtered_mdl=1 filtered_len=0 elapsed=0.51s` (`109,995,785` `aln_len`, `109,871,329` CIGAR `M/=/X` pairs).
    - `--net --query --target` (SDUST+MDL): `total=56857 kept=553 filtered_mdl=436 filtered_len=0 elapsed=4.89s` (`109,990,091` `aln_len`, `109,866,176` CIGAR `M/=/X` pairs).
  - Comparing the `565` and `553` outputs: `421` of the `436` MDL-filtered records were *already* discarded by `--net` alone. Between `--net` and `--net --query --target`, `15` records (`6,021` pairs) are dropped and `3` alternative records (`868` pairs) are netted in, yielding a net reduction of `5,153` pairs (`0.0047%` of `109,871,329`). Assuming `109,119,624` exact bp in both, unconditional exact precision changes from `99.3158%` (`109,119,624 / 109,871,329`) to `99.3205%` (`109,119,624 / 109,866,176`)—identical at `99.32%`—while runtime increases `9.6×` (`0.51 s` → `4.89 s`).
  - Running `python3 scripts/benchmarks/bench_vgp_filter_locus.py` fails in `scripts/benchmarks/score_paf_vs_hprc.py:7` with `ModuleNotFoundError: No module named 'numpy'` (and `/tmp/ext_tools/env/bin/python3` also lacks `numpy`). Furthermore, `score_paf_vs_hprc.py:151-169` does not even compute or emit unconditional exact precision (`exact_pairs / total_pairs`); it only emits `exact_coord_precision_on_shared`.
- **Inferential flaw & correction:** C2 attributes `99.12%` recall and `99.65%`/`99.32%` precision to the `553`-record SDUST+MDL run (`436` filtered by MDL), obscuring that `--net` alone already achieves the identical rounded metrics (`99.12%` / `99.65%` / `99.32%`) and that SDUST+MDL only alters `0.0047%` of aligned pairs on `FastGA` `chr6` at `9.6×` runtime. Narrow C2 to disclose that `--net` alone accounts for >99.99% of the filtering effect and that `score_paf_vs_hprc.py` cannot currently re-verify the `109,119,624` bp / `99.32%` figures without `numpy` and an unconditional-precision field.

---

## F3 — C3 (`manuscript/horizon_paper.md:40`): Conformal Risk Control is neither fixed/offline nor valid at `α = 0.07` for full alignments

- **Claim:** C3 (*"Using fixed Conformal Risk Control defaults calibrated offline at risk level `α = 0.07` (`n_cal = 15`; Angelopoulos et al. 2024), Horizon enforces a constant bit-savings floor `λ̂ = 10.0` bits (`min_delta_bits = 10.0`) and a minimum column posterior floor (`0.15`)."*)
- **Verdict:** **REFUTED**
- **Prior-round context:** Round 5 (L2 F4) adjudicated random-GC null reporting in `C4`; the specific statistical/calibration mechanics of `C3` were not previously adjudicated.
- **Evidence:**
  1. **Direct contradiction with `src/mdl/mod.rs:401-428`:** The aligner does **not** use a fixed `λ̂ = 10.0` bits calibrated offline. At `src/mdl/mod.rs:401-428`, whenever `!dp_jobs.is_empty()`, it computes `observed_gc` from the first `4,096` bp of `seq_x` and `seq_y` and calls `dp::calibrate_null_crc_for_gc(observed_gc, 15, 0.07, ctx.effective_min_bits, &hmm_params)` **online at runtime**.
  2. **Exchangeability and scale-mismatch violation (`src/mdl/dp.rs:120-165`):** `calibrate_null_crc_for_gc` generates `n_cal = 15` synthetic i.i.d. Bernoulli decoy pairs of fixed length **`win_len = 64` bp** (`band_w = 16`). These `64`-bp synthetic draws are **not exchangeable** with real inter-anchor gaps (up to `max_gap_bp = 2400` bp, `band_w = 180`), X-drop extensions (`max_ext = 1500` bp), or chained LCBs on `10 kb+` / chromosome-scale sequences. Consequently, the finite-sample Conformal Risk Control guarantee (`(count + 1) / (n_cal + 1) <= α`) does **not** hold for the alignment pipeline: as `manuscript/horizon_paper.md:61` and `docs/LARGE_SCALE_DIAGNOSTICS_PLAN.md:66` acknowledge, `10,000`-bp random DNA still emits `80–448` false-positive pairs on `4/5` other `18%`-GC seeds, `288` pairs on `30%`-GC, and `149` pairs on `45%`-GC.
  3. **Degenerate `n_cal = 15, α = 0.07` & `base_floor_bits` masking (`src/mdl/dp.rs:102-113`):** For `n_cal = 15` and `α = 0.07`, `(n_ge + 1) / 16 <= 0.07` forces `n_ge = 0` (always selecting the sample maximum of 15 draws). Moreover, `src/mdl/mod.rs:428` passes `ctx.effective_min_bits = min_delta_bits.max(0.55 * log2(N * M))` as `base_floor_bits` into `banded_mdl_gap_with_min_post` on `64`-bp windows, where `effective_min_bits >= 10.0` already causes almost all `64`-bp null draws to abstain (`0.0` bits), making `calibrate_null_crc_for_gc` a near-no-op wrapper around `effective_min_bits`.
- **Correction:** Remove the claim that Horizon enforces a constant offline-calibrated CRC floor `λ̂ = 10.0` bits at `α = 0.07` (`n_cal = 15`); instead state accurately that `effective_min_bits = max(10.0, 0.55 * log2(N * M))` serves as a heuristic MDL bit-savings floor and that `calibrate_null_crc_for_gc` only checks `15` synthetic `64`-bp Bernoulli windows online without controlling false-positive risk at `α = 0.07` on full-length alignments.

---

## F4 — C4 (`manuscript/horizon_paper.md:84`): `horizon_guide --progressive --hal` emits a reference-cloning pectinate HAL wrapper, `verify_merge_guide_maf.py` never evaluates `.hal` accuracy, and Abstract contradicts Results

- **Claim:** C4 (*"and `horizon_guide --progressive --hal --rbest` (`2896229` / `HEAD`) aligned `5` mammals (`hg38`, `Macaca_mulatta`, `HLcalPym1`, `Equus_caballus`, `Canis_lupus`) in `0.859 s` (`1,109` MAF blocks, `967` multi-species) and on 5-mammal UCSC `chrM.maf` (`verify_merge_guide_maf.py`, `0.258 s`) emitted a progressive ancestral-node HDF5 `v0` HAL (`((((hg38,panTro4)Anc0,rheMac3)Anc1,mm10)Anc2,canFam3)Anc3;`, `2,103` bottom segments) with `95.5%/99.8%` (`panTro4`), `89.2%/97.1%` (`rheMac3`), `71.8%/97.5%` (`mm10`), and `74.1%/98.2%` (`canFam3`) MAF recall/precision."*)
- **Verdict:** **NARROWED**
- **Evidence:**
  1. **No ancestral reconstruction or tree inference in `.hal` (`src/hal/mod.rs:245-359`):** `HalAlignment::from_multi_contig_progressive_blocks` wraps `from_multi_contig_star_blocks` in a hardcoded left-deep pectinate tree in CLI argument order (`((...((ref_name,q0)Anc0,q1)Anc1...),q_{k-1})Anc_{k-1};`, lines 355–359). Every internal ancestral node (`Anc0..Anc3`) literally clones `ref_genome.sequences`, `ref_genome.dna`, and `ref_genome.bottom_segments` (`hg38`, lines 302–323) with 1:1 identity `top_segments` to `ref_genome`. Calling this a "progressive ancestral-node" HAL implies ancestral sequence/coordinate inference that does not occur.
  2. **`.hal` alignment accuracy is never tested (`scripts/benchmarks/verify_merge_guide_maf.py:139-155`):** `verify_merge_guide_maf.py` only runs `horizon_wga --hal-info` (header/segment-count check) on `{out_prefix}.hal` and computes all recall/precision numbers solely from `{out_prefix}.maf`. Furthermore, `verify_merge_guide_maf.py:46` concatenates all UCSC `chrM.maf` blocks per species into a single synthetic `f"{sp}.chrM"` FASTA before aligning.
  3. **Internal manuscript contradiction:** `manuscript/horizon_paper.md:10` (Abstract) still states that `horizon_guide --progressive --hal` (`2896229`) builds *"a single-contig star-tree HDF5 `v0` HAL stub (`4` bottom segments) in `0.859 s`"*, contradicting `manuscript/horizon_paper.md:84` (C4).
- **Correction:** Narrow C4 to clarify that the reported recall/precision metrics are evaluated on `guide_out.maf` (on MAF-concatenated `chrM` FASTAs), that `.hal` output is a reference-cloning pectinate hierarchy verified only for HDF5 `v0` structure (`--hal-info`) rather than HAL alignment accuracy, and reconcile the Abstract (`line 10`) with Results (`line 84`).

---

## F5 — C5 (`manuscript/horizon_paper.md:150`): Citation matches `refs/catalog_index.json` and `manuscript/references.bib`

- **Claim:** C5 (*"- **Wittler R** (2023). General encoding of canonical *k*-mers. *Peer Community Journal*, 3:e87."*)
- **Verdict:** **SURVIVED**
- **Evidence:** Verified in `refs/catalog_index.json:1430-1433` (`slug: wittler_2023_canonical-kmer-encoding`, `doi: 10.24072/pcjournal.323`), `manuscript/references.bib:481-489`, and `manuscript/horizon_paper.md:28, 150`. Because `../library` is unavailable in this workspace, no claim is made about the full text of Wittler (2023).

---

## F6 — Claim-list steering check: C1–C5 omit the manuscript's primary inferential weaknesses

- **Evidence:** Per `reports/rounds/round-7/round.json`, `runner_is_author: true` and `claim_list_author: "runner"`. Selecting `1` test-count claim (C1), `1` PAF-filter claim (C2), `1` CRC paragraph (C3), `1` `chrM` guide-tree claim (C4), and `1` bibliography entry (C5) steers scrutiny away from three broader inferential defects in `manuscript/horizon_paper.md` and `docs/BENCHMARK_RESULTS.md`:
  1. **`n = 1` benchmark reporting:** Runtime, RSS, and random-GC (`30%`, `45%` GC) results are single-run (`n = 1`) measurements with no variance or confidence intervals (`manuscript/horizon_paper.md:78, 84`).
  2. **MAF-concatenated evaluation (`horizon_eval`):** Table 1 (`manuscript/horizon_paper.md:61-72`) evaluates pre-extracted, MAF-concatenated homologous `(ref_id, tgt_id)` pairs rather than whole-genome all-vs-all scaffolds.
  3. **Unfavorable head-to-head `chr6` baseline (`minimap2 -x asm5`):** On full `GRCh38.chr6` × `CHM13.chr6` (`manuscript/horizon_paper.md:78`), `minimap2 -x asm5` is `6.5×` faster (`46.95 s` vs `304.53 s`), has higher exact recall (`99.71%` vs `97.96%`), and has vastly higher exact precision (`99.75%` cond. / `99.60%` uncond. vs `48.63%` cond. / `38.22%` uncond.), while `FastGA` / `LASTZ` / `KegAlign` / `SweepGA` nuclear accuracy baselines remain pending (`manuscript/horizon_paper.md:61, 98`).

---

## Verdict Summary

| claim_id | verdict | summary |
| :--- | :--- | :--- |
| C1 | NARROWED | `122` active tests pass and `clippy` has `0` warnings on `HEAD`, but `1` `.hal` test (`real_hal_decode_from_env`) is `#[ignore]`d without an external `HORIZON_TEST_HAL` file. |
| C2 | NARROWED | `--net` alone (`565` records, `0.51 s`) already yields `99.12%`/`99.65%`/`99.32%`; SDUST+MDL (`553` records, `4.89 s`) changes only `0.0047%` of pairs (`5,153` bp), and `score_paf_vs_hprc.py` lacks `numpy` and unconditional precision. |
| C3 | REFUTED | `src/mdl/mod.rs:401-428` runs online (not fixed offline `λ̂ = 10.0`), calibrates on `15` synthetic `64`-bp Bernoulli draws non-exchangeable with full alignments, and fails to control FP risk at `α = 0.07` on `10`-kb random DNA. |
| C4 | NARROWED | `from_multi_contig_progressive_blocks` clones `hg38` into `Anc0..Anc3` on a hardcoded pectinate tree, `verify_merge_guide_maf.py` never checks `.hal` alignment accuracy, and Abstract (`line 10`) contradicts Results (`line 84`). |
| C5 | SURVIVED | Citation matches `refs/catalog_index.json:1430-1433` and `manuscript/references.bib:481-489` (full text unavailable without `../library`). |
