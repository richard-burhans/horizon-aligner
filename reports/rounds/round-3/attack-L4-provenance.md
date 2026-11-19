# Round 3 Adversarial Attack -- Provenance Lens

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `488b020ccc0a73998915e975bad8e25599b16716`, 2026-11-19 13:11:30 -0500)  
**Lens:** `provenance`  
**Method cards / prior history:** No method cards exist in this repository (checked `find . -maxdepth 3 -not -path '*/target/*' -not -path '*/.git/*'`). Checked `scripts/round history` for `stream_gfa`, `cactus241way`, `FlatSeedIndex`, and `horizon_gfa`, and read `reports/adversarial_log.md` and `reports/rounds/round-2/adjudication-L4.md` / `adjudication-L5.md` / `adjudication-L7.md`: Round 2 discharged `R1-GFA` and left `R2-HMM-XDROP`, `R2-RESOLVE-LIFT`, and `R2-INDEX-RAM` open.

---

## F1 -- C1: `stream_gfa` implementation and unit tests match the quoted claim

- **Claim:** `C1` (`docs/STRATEGIC_PLAN.md:13`) -- "`stream_gfa` provides a buffered line-streaming GFA v1.0/v1.1 reader (supporting `.gfa`, multi-member BGZF `.gfa.gz` via `MultiGzDecoder`, and `"-"` `stdin`, with unit tests for empty `L`-line orientations and `.gfa.gz`)."
- **Verdict:** `SURVIVED`
- **Evidence:** `src/gfa/mod.rs:97-127` implements `stream_gfa` with `"-"` `stdin` (`BufReader::new(stdin.lock())`) and `.gz` (`MultiGzDecoder::new(file)`), `src/gfa/mod.rs:175-194` defaults empty `L`-line orientation fields to `b'+'`, `src/gfa/mod.rs:229-287` (`streams_all_record_types_and_handles_empty_link_orients`) tests empty `L`-line orientations, and `src/gfa/mod.rs:290-329` (`streams_gzipped_gfa_file`) tests two concatenated gzip members. Both tests pass under `cargo test --no-fail-fast`.

---

## F2 -- C2: `cactus241way chr22 1MB` lacks a committed runner script/log, `run_all_evals.py` points to a nonexistent path, and `R2-HMM-XDROP` remains undischarged

- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:16,47`) -- "UCSC cactus241way `chr22 1MB`: **62.3%** recall / **60.3%** precision (Dog, across 129 target scaffolds; **61.3% F1** in **17.16s**)."
- **Verdict:** `NARROWED`
- **Flaw:** While `cactus241way_chr22_1MB.maf` exists on disk, its evaluation has no committed runner script or saved timing/eval log, the repository's evaluation script points to a nonexistent directory and omits `cactus241way`, and `488b020` did not re-run all WGA benchmarks on a single pinned commit as required by `R2-HMM-XDROP`.
- **Evidence:**
  1. `scripts/benchmarks/run_all_evals.py:10` sets `BENCH_DIR = ROOT.parent / "data" / "benchmarks"` (`/home/rcb112/project/gemini/bench/data/benchmarks`), which does not exist (`ls -la /home/rcb112/project/gemini/bench/data/` fails with `No such file or directory`). The actual datasets sit in `/home/rcb112/project/gemini/bench/pg/data/benchmarks/`, and `run_all_evals.py:56-76` has no `cactus241way` evaluation entry.
  2. No script in `scripts/` or log in `reports/` records the `cactus241way chr22 1MB` command or `62.3% / 60.3% / 17.16s` output; `/tmp/test_dog_1m.maf` has mtime `Nov 19 11:26` (predating `488b020` at `13:11:30`).
  3. Round 2 (`CARRIED.md`, `R2-HMM-XDROP`) required *"re-run all WGA benchmarks on a single pinned commit"*. Commit `488b020` claims to *"Discharge R2-HMM-XDROP"*, yet `docs/BENCHMARK_RESULTS.md:5` explicitly states that five earlier WGA rows (`multiz100way`, `cactus447way`, `multiz470way`, `cactus241way chr22_random`, `GRCh38.chr22 vs CHM13.chr22`) were recorded prior to `FlatSeedIndex` / `rev_offsets` fixes.
- **Correction:** Fix `BENCH_DIR` in `scripts/benchmarks/run_all_evals.py:10`, add `cactus241way_chr22_1MB.maf` to `run_all_evals.py`, archive the `488b020` evaluation log, and keep `R2-HMM-XDROP` open until all WGA rows are re-run on one pinned commit.

---

## F3 -- C3: Unit test count is `46` passing (not `45`) and `src/bin/` has 9 binaries (not 10)

- **Claim:** `C3` (`docs/STRATEGIC_PLAN.md:3`) -- "All target pangenome binaries (`horizon_wga`, `horizon_wga --graph`, `horizon_merge`, `horizon_guide`, `horizon_eval`, `horizon_gfa`, `horizon_clean`, `horizon_resolve`, `horizon_lift`, `bench_chr1`) are now implemented in `src/bin/` with unit tests in every binary (`45` passing tests + `1` ignored; `cargo clippy --all-targets -- -D warnings` clean)."
- **Verdict:** `NARROWED`
- **Flaw:** The documented passing test count (`45`) is off by one on `HEAD` (`488b020`), repeating the same count drift adjudicated in Round 2 (`adjudication-L4.md:19`, `adjudication-L5.md:17`, `adjudication-L7.md:28`, where `STRATEGIC_PLAN.md` reported `39` instead of `40`). Additionally, `round.json` summarizes C3 as *"All 10 binaries in src/bin/"*, whereas `src/bin/` contains 9 binary files (`horizon_wga --graph` is a CLI flag, not a separate binary).
- **Evidence:** Running `/home/agent/.cargo/bin/cargo test --no-fail-fast` on `488b020` outputs `37 passed; 0 failed; 1 ignored` in `src/lib.rs` plus `1 passed` in each of the 9 binaries (`bench_chr1`, `horizon_clean`, `horizon_eval`, `horizon_gfa`, `horizon_guide`, `horizon_lift`, `horizon_merge`, `horizon_resolve`, `horizon_wga`) = **46 passing tests + 1 ignored**. `/home/agent/.cargo/bin/cargo clippy --all-targets -- -D warnings` exits `0`.
- **Correction:** Update `docs/STRATEGIC_PLAN.md:3` from `45` passing tests to `46` passing tests (`37` library + `9` binary tests; `1` ignored).

---

## F4 -- C4: `FlatSeedIndex` allocates `28 bytes/entry` at build time, RSS reductions lack saved logs and confound `488b020` changes, and `R2-INDEX-RAM` is uncompleted

- **Claim:** `C4` (`docs/STRATEGIC_PLAN.md:11`) -- "`ContigIndex` (`src/seeding/mod.rs`) now uses a compact `FlatSeedIndex` (`12 bytes/entry` sorted `Vec<u64>` + `Vec<u32>` with on-the-fly reverse-complement coordinates, cutting `20 Mbp chr6 q-arm` peak RSS from `5,863.3 MB` to `2,409.7 MB` and `chr6` 5 Mbp pericentromere RSS from `3,084.5 MB` to `1,048.2 MB`)."
- **Verdict:** `NARROWED`
- **Flaw:** The `12 bytes/entry` figure describes only the post-build steady state, omitting the `28 bytes/entry` peak allocation during `FlatSeedIndex::build`; the RSS before/after figures lack saved `/usr/bin/time -v` logs and confound `FlatSeedIndex` with other `488b020` changes; and `R2-INDEX-RAM`'s >100 Mbp chromosome target remains uncompleted.
- **Evidence:**
  1. `FlatSeedIndex::build` (`src/seeding/mod.rs:386,424-443`) allocates `pairs: Vec<(u64, u32)>` (16 bytes/entry due to 8-byte alignment padding) simultaneously with `keys: Vec<u64>` + `pos: Vec<u32>` (12 bytes/entry) before `drop(pairs)` at line 441 (**28 bytes/entry** peak).
  2. No script or `/usr/bin/time -v` log in `scripts/` or `reports/` records `5,863.3 MB`, `2,409.7 MB`, `3,084.5 MB`, or `1,048.2 MB`; `/tmp/chr6_q20m.gfa` (`Nov 19 11:14`) and `/tmp/chr6_peri_5m.gfa` (`Nov 19 11:36`) predate `488b020` (`13:11:30`). In pre-`488b020` `docs/BENCHMARK_RESULTS.md:59-60`, both `chr22` and `chr6` 5 Mbp pericentromeres reported the exact same `3,084.5 MB` RSS. Moreover, `488b020` also modified `rev_offsets`, `mdl_xdrop_extend_with_params`, `bocpd_gc_profile`, and `horizon_wga.rs`, so the RSS change is not isolated to `FlatSeedIndex`.
  3. `docs/BENCHMARK_RESULTS.md:26` admits full `GRCh38.chr6` (`170.8 Mbp`) vs `CHM13.chr6` (`110.7 Mbp`) has only an uncompleted `"mid-run check at 13:59 used 10.09 GB RSS"`, leaving `R2-INDEX-RAM` undischarged.
- **Correction:** Report `12 bytes/entry` stored (`28 bytes/entry` peak during `FlatSeedIndex::build`), archive `/usr/bin/time -v` logs for the `chr6` benchmarks, and keep `R2-INDEX-RAM` open until a full >100 Mbp chromosome alignment completes.

---

## F5 -- C5: `chr22` 5 Mbp pericentromere `488b020` run has no saved script, output artifact, or timing/RSS log

- **Claim:** `C5` (`docs/BENCHMARK_RESULTS.md:23,59`) -- "Real HPRC `chr22` Pericentromere 5 Mbp (`GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M`): **14.57s**, **1,861.9 MB** peak RSS, **5,106,541** many-to-many `(x, y)` aligned pairs (`unique_x_bp=4,387,655`, `unique_y_bp=4,549,858`), `16,083` Horizon Bubbles."
- **Verdict:** `NARROWED`
- **Flaw:** As previously flagged for the Round 2 run in `reports/rounds/round-2/attack-L4-provenance.md` (F4) and `adjudication-L4.md`, no script in `scripts/` documents how the 5 Mbp slices were extracted or executes `horizon_wga` on them, and no `/usr/bin/time -v` or stdout log in `reports/` records the `488b020` numbers (`14.57s`, `1,861.9 MB`, `5,106,541`, `4,387,655`, `4,549,858`, `16,083`).
- **Evidence:** Input slices `/tmp/grch38_chr22_peri_5m.fa` and `/tmp/chm13_chr22_peri_5m.fa` exist (`Nov 19 11:35`), but the only `/tmp/` output GFA (`/tmp/chr22_peri_5m.gfa`, `Nov 19 11:36`) predates `488b020` (`Nov 19 13:11:30`).
- **Correction:** Commit the pericentromere extraction/benchmark command and save the `488b020` `horizon_wga` + `/usr/bin/time -v` log under `reports/`.

---

## F6 -- C6: `horizon_gfa stats` (`88.84s`, `97.1 MB`) is stale from `b857377` across three subsequent edits to `horizon_gfa.rs` / `gfa/mod.rs`

- **Claim:** `C6` (`docs/BENCHMARK_RESULTS.md:19,54`) -- "Full-graph streaming (`horizon_gfa stats`): **88.84s**, **97.1 MB** peak RSS across `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks."
- **Verdict:** `NARROWED`
- **Flaw:** The `88.84s` / `97.1 MB` measurement was recorded in commit `b857377` (`10:47:45`) and never re-run after `0e70693`, `b18f25b`, and `488b020` modified `src/bin/horizon_gfa.rs` and `src/gfa/mod.rs`; it also has no saved runner script or `/usr/bin/time -v` log, and `docs/BENCHMARK_DATASETS.md:46` cites `hprc-v1.1-mc-grch38.gfa.gz` rather than `hprc-v1.0-mc-grch38.gfa.gz`. All four other benchmarks bundled into `round.json`'s C6 summary (`horizon_gfa extract`, `horizon_lift`, `horizon_resolve`, `horizon_wga --graph`) likewise predate `488b020` edits to their respective source files.
- **Evidence:** `git log -S "88.84s" --oneline` traces the figure to `b857377`; `git show --stat 488b020` shows 68 lines changed in `src/bin/horizon_gfa.rs`, 35 in `src/gfa/mod.rs`, 8 in `src/bin/horizon_lift.rs`, 53 in `src/bin/horizon_resolve.rs`, and 109 in `src/bin/horizon_wga.rs`. The 12 GB file exists at `/home/rcb112/project/gemini/bench/pg/data/benchmarks/hprc-v1.0-mc-grch38.gfa.gz`, but no log exists in `reports/`.
- **Correction:** Tag the `horizon_gfa stats` row with its generating commit (`b857377`) or re-run `horizon_gfa stats`, `extract`, `lift`, `resolve`, and `wga --graph` on `488b020` with saved logs.

---

## F7 -- Claim-list steering: `round.json` and `488b020` claim `R2-HMM-XDROP`, `R2-RESOLVE-LIFT`, and `R2-INDEX-RAM` are discharged while omitting their unfulfilled requirements

- **Claim:** Framing of Round 3 (`round.json`, `runner_is_author: true`, `claim_list_author: "runner"`, and commit `488b020` *"Discharge R2-HMM-XDROP, R2-RESOLVE-LIFT, R2-INDEX-RAM, and 5-axis review fixes"*).
- **Verdict:** `NARROWED`
- **Flaw:** The runner-authored claim list steers away from three explicit carried requirements from Round 2 (`CARRIED.md`) that remain unfulfilled on `488b020`:
  1. `R2-HMM-XDROP` required re-running **all** WGA benchmarks on a single pinned commit; `docs/BENCHMARK_RESULTS.md:5` admits 5 WGA rows remain from pre-`488b020` commits.
  2. `R2-RESOLVE-LIFT` required **head-to-head pangenome baselines (`minimap2`/`odgi`/`vg`)**; `docs/BENCHMARK_RESULTS.md` contains zero `minimap2`, `odgi`, or `vg` rows.
  3. `R2-INDEX-RAM` required bringing full `>100 Mbp` chromosome alignment well below `26.5 GB` RSS; `docs/BENCHMARK_RESULTS.md:26` reports only an uncompleted mid-run check (`10.09 GB` at `13:59`).
- **Evidence:** `reports/rounds/round-3/CARRIED.md`, `docs/BENCHMARK_RESULTS.md:5,26,50-65`, and `git show --stat 488b020`.
- **Correction:** Keep `R2-HMM-XDROP`, `R2-RESOLVE-LIFT`, and `R2-INDEX-RAM` open until all WGA rows are re-run on one commit, head-to-head `minimap2`/`odgi`/`vg` baselines are recorded, and a >100 Mbp chromosome alignment completes.

---

## Verdict Summary

| Claim | Verdict |
| :--- | :--- |
| C1 | SURVIVED |
| C2 | NARROWED |
| C3 | NARROWED |
| C4 | NARROWED |
| C5 | NARROWED |
| C6 | NARROWED |
