# Round 3 ADJUDICATION of L4 -- Provenance Lens

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `488b020ccc0a73998915e975bad8e25599b16716`, 2026-11-19 13:11:30 -0500)  
**Attack file:** `reports/rounds/round-3/attack-L4-provenance.md` (SHA-256 prefix verified: `e1e85bbbf108`)  
**Prior history & method cards:** Checked `scripts/round history` for `stream_gfa`, `FlatSeedIndex`, `cactus241way`, and `horizon_gfa`, and reviewed `reports/rounds/round-2/adjudication-L4.md` (which left `R2-HMM-XDROP`, `R2-RESOLVE-LIFT`, and `R2-INDEX-RAM` open). No method cards exist in `refs/methods/` in this repository.

---

## F1 -- C1: `stream_gfa` implementation and unit tests match the quoted claim

- **Diagnosis:** **UPHELD** (attacker's `SURVIVED` assessment stands). Verified directly in `src/gfa/mod.rs:97-127` (`stream_gfa` handling `"-"`, `.gz` via `MultiGzDecoder`, and plain files), `src/gfa/mod.rs:175-194` (defaulting empty `L`-line orientation fields to `b'+'`), and `src/gfa/mod.rs:229-329` (`streams_all_record_types_and_handles_empty_link_orients` and `streams_gzipped_gfa_file`), both of which pass under `/home/agent/.cargo/bin/cargo test --no-fail-fast`.
- **Repair:** None needed.
- **What survives:** `C1` accurately describes the `stream_gfa` implementation and its passing unit tests on `488b020`.

---

## F2 -- C2: `cactus241way chr22 1MB` lacks a committed runner script/log, `run_all_evals.py` points to a nonexistent path, and `R2-HMM-XDROP` remains undischarged

- **Diagnosis:** **UPHELD** (attacker's `NARROWED` ruling stands). Verified: (1) `scripts/benchmarks/run_all_evals.py:10` sets `BENCH_DIR = ROOT.parent / "data" / "benchmarks"` (`/home/rcb112/project/gemini/bench/data/benchmarks`), which does not exist (`ls` fails with `No such file or directory`), whereas datasets reside in `/home/rcb112/project/gemini/bench/pg/data/benchmarks/`, and `run_all_evals.py:56-76` omits `cactus241way`; (2) no script in `scripts/` or log in `reports/` records the `cactus241way chr22 1MB` run (`62.3%` / `60.3%` / `17.16s`), and `/tmp/test_dog_1m.maf` (`Nov 19 11:26`) predates `488b020` (`13:11:30`); (3) `docs/BENCHMARK_RESULTS.md:5` explicitly notes that five earlier WGA rows (`multiz100way`, `cactus447way`, `multiz470way`, `cactus241way chr22_random`, `GRCh38.chr22 vs CHM13.chr22`) were recorded prior to `FlatSeedIndex` / `rev_offsets` fixes, so `R2-HMM-XDROP`'s requirement to re-run all WGA benchmarks on a single pinned commit is not met.
- **Repair:** **NARROWED (split by cost).** Fixing `BENCH_DIR` and adding the `cactus241way_chr22_1MB.maf` target in `scripts/benchmarks/run_all_evals.py` and keeping `R2-HMM-XDROP` open in the round record require no re-measurement (**the fix is free**); re-running all WGA benchmarks on a single pinned commit and archiving the evaluation logs is non-free (compute/re-measurement).
- **What survives:** `C2`'s dataset exists, but its evaluation lacks a committed runner script or saved log, `run_all_evals.py` points to a nonexistent path and omits `cactus241way`, and `R2-HMM-XDROP` remains undischarged because five WGA rows predate `488b020`.

---

## F3 -- C3: Unit test count is `46` passing (not `45`) and `src/bin/` has 9 binaries (not 10)

- **Diagnosis:** **UPHELD** (attacker's `NARROWED` ruling stands). Verified by running `/home/agent/.cargo/bin/cargo test --no-fail-fast` on `488b020`: `src/lib.rs` reports `37 passed; 0 failed; 1 ignored`, and each of the **9** binary files in `src/bin/` (`bench_chr1`, `horizon_clean`, `horizon_eval`, `horizon_gfa`, `horizon_guide`, `horizon_lift`, `horizon_merge`, `horizon_resolve`, `horizon_wga`) reports `1 passed`, totaling **46 passing tests + 1 ignored** (not `45`), with `cargo clippy --all-targets -- -D warnings` exiting `0`. `horizon_wga --graph` is a CLI flag rather than a 10th binary in `src/bin/`. This repeats the same test-count drift adjudicated in Round 2 (`adjudication-L4.md` F2, `adjudication-L5.md` F2, `adjudication-L7.md` F2).
- **Repair:** **ENDORSED.** Updating `docs/STRATEGIC_PLAN.md:3` from `45` to `46` passing tests (`37` library + `9` binary tests; `1` ignored) and clarifying 9 binaries in `src/bin/` (`horizon_wga --graph` is a CLI mode) requires no re-measurement; **the fix is free**.
- **What survives:** `docs/STRATEGIC_PLAN.md:3` undercounts passing tests (`45` vs. `46` + `1` ignored), and `round.json` / `STRATEGIC_PLAN.md` conflate `horizon_wga --graph` with a 10th binary in `src/bin/`.

---

## F4 -- C4: `FlatSeedIndex` allocates `28 bytes/entry` at build time, RSS reductions lack saved logs and confound `488b020` changes, and `R2-INDEX-RAM` is uncompleted

- **Diagnosis:** **UPHELD** (attacker's `NARROWED` ruling stands). Verified: (1) in `src/seeding/mod.rs:386,424-443`, `FlatSeedIndex::build` allocates `pairs: Vec<(u64, u32)>` (16 bytes/entry with 8-byte alignment padding) concurrently with `keys: Vec<u64>` + `pos: Vec<u32>` (12 bytes/entry) before `drop(pairs)` at line 441, giving **28 bytes/entry** peak build-time allocation versus `12 bytes/entry` stored steady state; (2) no script or `/usr/bin/time -v` log in `scripts/` or `reports/` records `5,863.3 MB`, `2,409.7 MB`, `3,084.5 MB`, or `1,048.2 MB`, `/tmp/chr6_q20m.gfa` (`Nov 19 11:14`) and `/tmp/chr6_peri_5m.gfa` (`Nov 19 11:36`) predate `488b020` (`13:11:30`), `488b020^:docs/BENCHMARK_RESULTS.md` indeed reported `3,084.5 MB` for both `chr22` and `chr6` 5 Mbp pericentromeres, and `git show --stat 488b020` shows simultaneous edits to `src/mdl/dp.rs`, `src/mdl/mod.rs`, `src/bocpd.rs`, and `src/bin/horizon_wga.rs`; (3) `docs/BENCHMARK_RESULTS.md:26` reports only an uncompleted mid-run check (`10.09 GB` RSS at `13:59`) for full `chr6` (`170.8 Mbp` vs `110.7 Mbp`), leaving `R2-INDEX-RAM` undischarged.
- **Repair:** **NARROWED (split by cost).** Clarifying in `docs/STRATEGIC_PLAN.md:11` that `12 bytes/entry` is the stored steady-state size (`28 bytes/entry` peak during `FlatSeedIndex::build`) and keeping `R2-INDEX-RAM` open are **free**; archiving `/usr/bin/time -v` logs and completing a full `>100 Mbp` chromosome alignment are non-free.
- **What survives:** `FlatSeedIndex` stores `12 bytes/entry` post-build (`28 bytes/entry` peak during build), the claimed RSS reductions lack saved `/usr/bin/time -v` logs and are confounded with other `488b020` edits, and `R2-INDEX-RAM` remains incomplete.

---

## F5 -- C5: `chr22` 5 Mbp pericentromere `488b020` run has no saved script, output artifact, or timing/RSS log

- **Diagnosis:** **UPHELD** (attacker's `NARROWED` ruling stands). Verified: no script in `scripts/` extracts the 5 Mbp slices or runs `horizon_wga` on them, and no `/usr/bin/time -v` or stdout log in `reports/` records `14.57s`, `1,861.9 MB`, `5,106,541`, `4,387,655`, `4,549,858`, or `16,083`. In `/tmp/`, the input slices (`/tmp/grch38_chr22_peri_5m.fa`, `/tmp/chm13_chr22_peri_5m.fa`, `Nov 19 11:35`) and the only output GFA (`/tmp/chr22_peri_5m.gfa`, `Nov 19 11:36`) predate `488b020` (`Nov 19 13:11:30`). Round 2 (`adjudication-L4.md` F4) previously flagged the missing extraction script and timing/RSS log for this locus.
- **Repair:** **NARROWED (split by cost).** Marking the `C5` row in `docs/BENCHMARK_RESULTS.md` as unlogged is **free**; committing the extraction/benchmark script and re-running `horizon_wga` with `/usr/bin/time -v` to archive logs under `reports/` is non-free.
- **What survives:** `C5`'s `488b020` `chr22` pericentromere numbers have no committed runner script, `488b020`-dated output artifact, or saved timing/RSS log.

---

## F6 -- C6: `horizon_gfa stats` (`88.84s`, `97.1 MB`) is stale from `b857377` across three subsequent edits to `horizon_gfa.rs` / `gfa/mod.rs`

- **Diagnosis:** **UPHELD** (attacker's `NARROWED` ruling stands). Verified: `git log -S "88.84s" --oneline` traces the measurement to `b857377` (`10:47:45`), while `git show --stat 488b020` confirms `488b020` modified `src/bin/horizon_gfa.rs` (68 lines), `src/gfa/mod.rs` (35 lines), `src/bin/horizon_lift.rs` (8 lines), `src/bin/horizon_resolve.rs` (53 lines), and `src/bin/horizon_wga.rs` (109 lines), and `docs/BENCHMARK_RESULTS.md` lines 54-57 and 64-65 were left unchanged from earlier commits. No runner script or `/usr/bin/time -v` log exists in `scripts/` or `reports/`, and `docs/BENCHMARK_DATASETS.md:46` cites `hprc-v1.1-mc-grch38.gfa.gz` whereas `docs/BENCHMARK_RESULTS.md:18,50` and `/home/rcb112/project/gemini/bench/pg/data/benchmarks/hprc-v1.0-mc-grch38.gfa.gz` are `v1.0`.
- **Repair:** **NARROWED (split by cost).** Tagging `horizon_gfa stats` (and the other pre-`488b020` pangenome rows) in `docs/BENCHMARK_RESULTS.md` with their generating commits and fixing `v1.1` -> `v1.0` in `docs/BENCHMARK_DATASETS.md:46` require no re-measurement (**the fix is free**); re-running `horizon_gfa stats`, `extract`, `lift`, `resolve`, and `wga --graph` on `488b020` with saved logs is non-free.
- **What survives:** `horizon_gfa stats` (`88.84s`, `97.1 MB`) and the other C6 pangenome rows predate `488b020` edits to their binaries, lack saved scripts/logs, and `docs/BENCHMARK_DATASETS.md:46` misstates `v1.1` for `v1.0`.

---

## F7 -- Claim-list steering: `round.json` and `488b020` claim `R2-HMM-XDROP`, `R2-RESOLVE-LIFT`, and `R2-INDEX-RAM` are discharged while omitting their unfulfilled requirements

- **Diagnosis:** **UPHELD** (attacker's `NARROWED` ruling stands). Verified against `reports/rounds/round-3/CARRIED.md` and `docs/BENCHMARK_RESULTS.md`: (1) `R2-HMM-XDROP` required re-running all WGA benchmarks on a single pinned commit, yet `docs/BENCHMARK_RESULTS.md:5` admits five WGA rows remain pre-`488b020`; (2) `R2-RESOLVE-LIFT` required head-to-head pangenome baselines (`minimap2`/`odgi`/`vg`), yet `docs/BENCHMARK_RESULTS.md` contains zero `minimap2`, `odgi`, or `vg` benchmark rows; (3) `R2-INDEX-RAM` required completing a `>100 Mbp` chromosome alignment well below `26.5 GB` RSS, yet `docs/BENCHMARK_RESULTS.md:26` records only an uncompleted mid-run check (`10.09 GB` at `13:59`).
- **Repair:** **NARROWED (split by cost).** Keeping `R2-HMM-XDROP`, `R2-RESOLVE-LIFT`, and `R2-INDEX-RAM` open in the round record is **free**; executing the full pinned WGA suite, head-to-head `minimap2`/`odgi`/`vg` baselines, and full `>100 Mbp` chromosome alignment is non-free.
- **What survives:** All three Round 2 carried items (`R2-HMM-XDROP`, `R2-RESOLVE-LIFT`, `R2-INDEX-RAM`) have unfulfilled requirements on `488b020` despite the commit message claiming discharge.

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `src/gfa/mod.rs` implements `stream_gfa` for `.gfa`, multi-member `.gfa.gz`, and `"-"`, and both cited unit tests pass under `cargo test` |
| F2 | UPHELD | `run_all_evals.py:10` points to nonexistent `../data/benchmarks` and omits `cactus241way`, `cactus241way chr22 1MB` has no saved script/log, and 5 WGA rows predate `488b020`; fixing `BENCH_DIR` and keeping `R2-HMM-XDROP` open is free |
| F3 | UPHELD | `cargo test --no-fail-fast` passes `46` tests (`37` lib + `9` bin, `1` ignored) rather than `45`, and `src/bin/` has `9` binaries; updating `docs/STRATEGIC_PLAN.md:3` is free |
| F4 | UPHELD | `FlatSeedIndex::build` peaks at `28 bytes/entry` (`12 bytes/entry` stored), RSS claims lack `/usr/bin/time -v` logs and confound `488b020` edits, and full `chr6` only has a `10.09 GB` mid-run check; clarifying stored vs. peak and keeping `R2-INDEX-RAM` open is free |
| F5 | UPHELD | `chr22` 5 Mbp pericentromere `488b020` numbers have no committed extraction/runner script, `488b020`-dated output, or timing/RSS log; marking the row unlogged is free |
| F6 | UPHELD | `horizon_gfa stats` (`88.84s`, `97.1 MB`) is from `b857377` prior to `488b020` edits, C6 rows lack saved logs, and `BENCHMARK_DATASETS.md:46` cites `v1.1` instead of `v1.0`; commit-tagging and fixing `v1.0` is free |
| F7 | UPHELD | `R2-HMM-XDROP` (5 stale WGA rows), `R2-RESOLVE-LIFT` (0 `minimap2`/`odgi`/`vg` baselines), and `R2-INDEX-RAM` (uncompleted `chr6` run) remain unfulfilled on `488b020`; keeping all three open is free |
