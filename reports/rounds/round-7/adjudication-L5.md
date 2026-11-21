# Round 7 Adjudication — L5 (Inference)

**Seal verification:** `sha256sum reports/rounds/round-7/attack-L5-inference.md` = `cf804475ca8a` (matches `ADJ-BRIEF-L5.md`).
**Method cards (`refs/methods/`):** Absent in this workspace (`ls: cannot access 'refs/methods/': No such file or directory`); no method card governs these disputes. All rulings below are verified directly against `HEAD` (`5cb83fb6a1346a8193f9e4048055638d0a422b0a`), live test/binary runs, and `scripts/round history --kind adjudication`.

---

## F1 — C1 (`manuscript/horizon_paper.md:84`): `1` of `123` tests (`real_hal_decode_from_env`) is `#[ignore]`d on `HEAD` without an untracked external `HORIZON_TEST_HAL` fixture

- **Prior-round context:** Round 6 (`L1`, `L2`, `L7`) adjudicated the stale `118` vs. `120` active unit-test count on `42700e1`. On `HEAD` (`5cb83fb6`), `2` additional library tests were added, bringing the active count to `122`.
- **Verification:** Re-running `cargo test --all-features --quiet` on `HEAD` executes `88` passed + `1` ignored library tests and `34` passed binary tests (`1+2+6+2+1+2+4+11+1+1+1+2 = 34`), totaling `122` active passing tests and `1` ignored test (`src/hal/tests.rs:240-243`, `#[ignore = "requires HORIZON_TEST_HAL=/path/to/file.hal"]`). `cargo clippy --all-targets --all-features -- -D warnings` exits `0` with zero warnings, and `/tmp/halTest.hal` is absent.
- **Ruling:** **NARROWED.**
  - **What survives:** All `122` active Rust unit tests (`88` library + `34` binary) pass on `HEAD` with zero `clippy` warnings; only the `"123 total"` / `halTest.hal` clause depends on an untracked external fixture and is `#[ignore]`d by default.
  - **Repair:** **ENDORSED (text clarification); the fix is free.** Clarify in `manuscript/horizon_paper.md:84` that `real_hal_decode_from_env` is `#[ignore]`d by default unless an external `HORIZON_TEST_HAL` fixture is provided.

---

## F2 — C2 (`docs/BENCHMARK_RESULTS.md:86`): SDUST+MDL removes only `0.0047%` of pairs over `--net` alone at `~10×` runtime, and `score_paf_vs_hprc.py` neither runs (`numpy` missing) nor computes unconditional precision

- **Prior-round context:** Round 6 (`L1 F2`, `L2`) narrowed `horizon_filter --net` from `1,004` to `565` records on `42700e1` and recorded the missing `numpy` dependency (`R2-RESOLVE-LIFT`).
- **Verification:** Re-running `./target/release/horizon_filter` on `/tmp/chr6_fastga.paf` (`56,857` records) on `HEAD`:
  - `--net`: `kept=565 filtered_mdl=1 elapsed=0.48s` (`109,871,329` CIGAR `M/=/X` pairs).
  - `--net --query /tmp/hprc_grch38_chr6.fa --target /tmp/hprc_chm13_chr6.fa`: `kept=553 filtered_mdl=436 elapsed=4.83s` (`109,866,176` CIGAR `M/=/X` pairs).
  - Set diff confirms `15` records dropped, `3` added, net `-5,153` pairs (`0.0047%` of `109,871,329`) at `10.1×` runtime (`0.48 s` → `4.83 s`).
  - `numpy` is absent in both `python3` and `/tmp/ext_tools/env/bin/python3`, and `scripts/benchmarks/score_paf_vs_hprc.py:151-169` emits `exact_coord_precision_on_shared` (`exact_pairs / pairs_on_shared`) without computing unconditional precision (`exact_pairs / total_pairs`).
- **Ruling:** **NARROWED.**
  - **What survives:** `docs/BENCHMARK_RESULTS.md:85` already reports the identical `99.12%` / `99.65%` / `99.32%` figures for `--net` alone, so C2's `553`-record numbers are consistent with `--net`, but C2's framing (`436 filtered by MDL; removes 12 noisy/non-1-to-1 records`) obscures that `--net` alone accounts for >99.99% of the filtering and that SDUST+MDL alters only `0.0047%` of aligned pairs at `~10×` runtime, while `score_paf_vs_hprc.py` cannot currently re-verify unconditional precision.
  - **Repair:** **ENDORSED** for disclosing the `0.0047%` incremental delta over `--net` in `docs/BENCHMARK_RESULTS.md:86` (**the text fix is free**); re-running `score_paf_vs_hprc.py` requires installing `numpy` and adding an unconditional-precision field.

---

## F3 — C3 (`manuscript/horizon_paper.md:40`): Conformal Risk Control is neither fixed/offline nor valid at `α = 0.07` for full alignments

- **Prior-round context:** Round 5 (`L2 F4`) adjudicated random-GC null reporting; the CRC calibration mechanics in `C3` were not previously adjudicated.
- **Verification:**
  - `src/mdl/mod.rs:401-428` calls `dp::calibrate_null_crc_for_gc(observed_gc, 15, 0.07, ctx.effective_min_bits, &hmm_params)` **online** from the first `4,096` bp of `seq_x` and `seq_y` whenever `!dp_jobs.is_empty()`—contradicting C3's claim of fixed offline calibration at `λ̂ = 10.0` bits.
  - `src/mdl/dp.rs:120-165` generates `n_cal = 15` synthetic i.i.d. Bernoulli pairs of fixed length `win_len = 64` bp (`band_w = 16`), which are non-exchangeable with inter-anchor gaps (`max_gap_bp = 2400`, `band_w = 180`), X-drop extensions (`max_ext = 1500`), and chromosome-scale LCBs; `(n_ge + 1) / 16 <= 0.07` forces `n_ge = 0` (sample maximum). Empirically, `manuscript/horizon_paper.md:61` confirms `80–448` false-positive pairs on `4/5` other `18%`-GC seeds, `288` on `30%`-GC, and `149` on `45%`-GC.
- **Ruling:** **UPHELD.**
  - **What survives:** Nothing in C3's claim of a constant, offline-calibrated CRC floor `λ̂ = 10.0` bits enforcing risk level `α = 0.07` survives.
  - **Repair:** **ENDORSED; the text fix is free.** Replace C3's CRC claim (and matching wording in `manuscript/horizon_paper.md:10, 22, 40, 96`) to state that `effective_min_bits = max(10.0, 0.55 * log2(N * M))` serves as a heuristic MDL bit-savings floor and that `calibrate_null_crc_for_gc` checks `15` synthetic `64`-bp Bernoulli windows online without controlling false-positive risk at `α = 0.07` on full-length alignments.

---

## F4 — C4 (`manuscript/horizon_paper.md:84`): `horizon_guide --progressive --hal` emits a reference-cloning pectinate HAL wrapper, `verify_merge_guide_maf.py` never evaluates `.hal` accuracy, and Abstract contradicts Results

- **Verification:**
  - `src/hal/mod.rs:244-365` (`HalAlignment::from_multi_contig_progressive_blocks`) constructs a hardcoded left-deep pectinate tree in CLI argument order (`((...((ref_name,q0)Anc0,q1)Anc1...),q_{k-1})Anc_{k-1};`) and clones `ref_genome.sequences`, `ref_genome.dna`, and `ref_genome.bottom_segments` into every `Anc0..Anc_{n_q-1}` with 1:1 identity `top_segments` to `ref_genome`.
  - `scripts/benchmarks/verify_merge_guide_maf.py:46, 139-155` concatenates blocks into `f"{sp}.chrM"` (default `multi_contig=False`), runs only `horizon_wga --hal-info` on `{out_prefix}.hal`, and computes all recall/precision metrics from `{out_prefix}.maf`.
  - `manuscript/horizon_paper.md:10` (Abstract) still describes a `"single-contig star-tree HDF5 v0 HAL stub (4 bottom segments) in 0.859 s"`, contradicting `line 84` (`2,103` bottom segments).
- **Ruling:** **NARROWED.**
  - **What survives:** The `0.859 s` (`1,109` MAF blocks) and `0.258 s` `chrM` MAF recall/precision measurements survive, while the "progressive ancestral-node" HAL characterization, the implication that `.hal` alignment accuracy was verified, and the Abstract/Results contradiction are refuted.
  - **Repair:** **ENDORSED; the text fix is free.** Clarify in `manuscript/horizon_paper.md:84` that recall/precision are evaluated on `guide_out.maf` (on MAF-concatenated `chrM` FASTAs), that `.hal` is a reference-cloning pectinate wrapper verified only for HDF5 `v0` structure (`--hal-info`), and reconcile `line 10` with `line 84`.

---

## F5 — C5 (`manuscript/horizon_paper.md:150`): Citation matches `refs/catalog_index.json` and `manuscript/references.bib`

- **Verification:** Confirmed in `refs/catalog_index.json:1430-1433` (`wittler_2023_canonical-kmer-encoding`, `doi: 10.24072/pcjournal.323`), `manuscript/references.bib:481-489`, and `manuscript/horizon_paper.md:28, 150`.
- **Ruling:** **UPHELD** (the attacker's check is accurate; C5 survives intact).
  - **What survives:** The entire C5 citation survives as published (full-text verification remains out of scope without `../library`).

---

## F6 — Claim-list steering check: C1–C5 omit the manuscript's primary inferential weaknesses

- **Verification:** `reports/rounds/round-7/round.json` confirms `runner_is_author: true` and `claim_list_author: "runner"`, and C1–C5 include `C1` (test count) and `C5` (single bibliography entry) while omitting Table 1's MAF-concatenated evaluation, `n = 1` variance-free benchmark reporting, and the `chr6` `minimap2 -x asm5` precision gap. However, `manuscript/horizon_paper.md` already discloses all three limitations in the text (`lines 10, 61, 78, 84, 98`).
- **Ruling:** **NARROWED.**
  - **What survives:** The claim-list steering observation is accurate regarding what C1–C5 selected, though the manuscript itself already discloses `n = 1`, MAF concatenation, and the `minimap2` `chr6` baseline.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | `122` active tests and `0` clippy warnings pass on `HEAD`, but `real_hal_decode_from_env` (`123` total) is `#[ignore]`d without an untracked `HORIZON_TEST_HAL` file; text fix is free. |
| F2 | NARROWED | `--net` alone (`565` records, `0.48 s`) already achieves `99.12%`/`99.65%`/`99.32%`; SDUST+MDL (`553` records, `4.83 s`) changes `0.0047%` of pairs, and `score_paf_vs_hprc.py` lacks `numpy` and unconditional precision; text disclosure is free. |
| F3 | UPHELD | `src/mdl/mod.rs:401-428` runs online (not fixed offline `λ̂ = 10.0`) on `15` non-exchangeable `64`-bp Bernoulli draws (`src/mdl/dp.rs:120-165`) and fails to control FP risk at `α = 0.07`; text fix is free. |
| F4 | NARROWED | `0.859 s` / `chrM` MAF metrics survive, but `.hal` clones `hg38` into `Anc0..Anc3` on a hardcoded pectinate tree, `verify_merge_guide_maf.py` never checks `.hal` accuracy, and Abstract (`line 10`) contradicts Results (`line 84`); text fix is free. |
| F5 | UPHELD | Citation matches `refs/catalog_index.json:1430-1433` and `manuscript/references.bib:481-489`; C5 survives intact. |
| F6 | NARROWED | Claim-list steering check accurately notes `C1`–`C5` omit `n = 1`, MAF-concatenation, and the `chr6` `minimap2` gap, though `manuscript/horizon_paper.md` already discloses them. |
