# Round 2 Adversarial Attack -- Provenance Lens

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `b18f25b444985a5b7f1db38fd9dd312240230839`, 2026-11-19 12:27:59 -0500)  
**Lens:** `provenance`  
**Method cards / prior history:** No method cards exist in this repository (checked via `find . -maxdepth 3 -iname "*card*"`). Checked `./scripts/round history "stream_gfa" --kind adjudication` and `reports/adversarial_log.md`: Round 1 adjudicated `stream_gfa`, `horizon_clean`/`horizon_resolve`/`horizon_lift`, `src/coord.rs`, and `src/mdl/dp.rs` transitions (leaving `R1-IMPL`, `R1-GFA`, `R1-HMM` open); Round 2 evaluates the post-Round-1 implementations and benchmark claims.

---

## F1 -- C2: Stale multi-commit provenance, missing benchmark datasets, and hardcoded `HmmParams::default()` in `--graph`

- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:5`) -- "All results were generated out-of-the-box using the divergence-adaptive Pair-HMM (`HmmParams`: `identity_prior`, `gap_open_bits`, and `gap_extend_bits` estimated from LCB anchors with Beta prior shrinkage in `src/mdl/mod.rs`), BOCPD compositional profiles (`src/bocpd.rs`), and normalized Pair-HMM transitions (`src/mdl/dp.rs`, `src/mdl/dag.rs`)."
- **Verdict:** `REFUTED`
- **Flaw:** `docs/BENCHMARK_RESULTS.md` aggregates results from at least three different git commits (`b857377`, `7bf688c`, `b18f25b`), several of which predate the adaptive `HmmParams`, BOCPD jump-gate activation, and `src/mdl/dag.rs` code claimed to have generated all results. Additionally, `horizon_wga --graph` does not use adaptive `HmmParams` at all, and the WGA/HPRC source datasets are absent from the workspace.
- **Evidence:**
  1. `git show b857377 -- docs/BENCHMARK_RESULTS.md` (10:47:45 -0500) introduced `Chr22 Haplotype Alignment` (`154.92s`, `24,365,611` pairs, `9,690` bubbles), and `/tmp/hprc_chr22.gfa` / `/tmp/hprc_chr22.maf` have mtime `Nov 19 10:47`. That run predates adaptive `HmmParams` and `MdlCostModel` (`7bf688c`, 11:17:25), BOCPD jump-gate activation (`0e70693`, 11:20:04), and adaptive `gap_extend_bits` (`b18f25b`, 12:27:59), and was never re-run in `b18f25b`.
  2. `git show b18f25b -- docs/BENCHMARK_RESULTS.md` re-ran only `cactus241way (chr22 1MB)` (`64.3% -> 63.6%`), leaving `multiz100way` (Chimp/Mouse), `cactus447way` (Horse/Macaque), `multiz470way` (Marmoset), `cactus241way (chr22_random)`, and `horizon_lift` unchanged from `7bf688c` (11:17:25) -- before `0e70693` (11:20:04) activated the BOCPD jump gate and before `b18f25b` (12:27:59) made `gap_extend_bits` adaptive.
  3. `src/bin/horizon_wga.rs:195` and `src/bin/horizon_wga.rs:504` pass `&HmmParams::default()` to `align_seq_to_dag_seeded`, not LCB-estimated `HmmParams`.
  4. `scripts/benchmarks/run_all_evals.py:10` points to `../data/benchmarks`, which does not exist on disk (`ls -lh /home/rcb112/project/gemini/bench/data/` returned empty), and `hprc-v1.0-mc-grch38.gfa.gz` is not in the workspace.
- **Correction:** Restrict the provenance header in `docs/BENCHMARK_RESULTS.md:5` to the specific commit (`b18f25b`) and re-run all WGA, `chr22`, and `horizon_lift` benchmarks on `HEAD`, or tag each table row with its exact generating commit and note that `--graph` uses `HmmParams::default()`.

---

## F2 -- C3: Unit test count is 40 (not 39) and 4 of 9 named binaries have 0 unit tests

- **Claim:** `C3` (`docs/STRATEGIC_PLAN.md:3`) -- "All target pangenome binaries (`horizon_wga`, `horizon_wga --graph`, `horizon_merge`, `horizon_guide`, `horizon_eval`, `horizon_gfa`, `horizon_clean`, `horizon_resolve`, `horizon_lift`, `bench_chr1`) are now implemented in `src/bin/` and unit-tested (`39` passing tests)."
- **Verdict:** `REFUTED`
- **Flaw:** The test count in `docs/STRATEGIC_PLAN.md:3` is stale (`39` instead of `40` + `1` ignored), and 4 of the 9 binaries in `src/bin/` have zero unit tests.
- **Evidence:** Running `cargo test` on `HEAD` (`b18f25b`) outputs:
  - `src/lib.rs`: `35 passed; 0 failed; 1 ignored`
  - `src/bin/horizon_clean.rs`: `1 passed`
  - `src/bin/horizon_gfa.rs`: `1 passed`
  - `src/bin/horizon_lift.rs`: `1 passed`
  - `src/bin/horizon_resolve.rs`: `1 passed` (added in `b18f25b`)
  - `src/bin/horizon_wga.rs`: `1 passed`
  - `src/bin/bench_chr1.rs`, `src/bin/horizon_eval.rs`, `src/bin/horizon_guide.rs`, `src/bin/horizon_merge.rs`: **`0 passed`** (0 unit tests each).
  Total passing unit tests = `40` (not `39`).
- **Correction:** Update `docs/STRATEGIC_PLAN.md:3` to report `40 passing tests (35 library + 5 binary tests across horizon_clean, horizon_gfa, horizon_lift, horizon_resolve, horizon_wga; 1 ignored)` and explicitly note that `bench_chr1`, `horizon_eval`, `horizon_guide`, and `horizon_merge` currently have `0` unit tests.

---

## F3 -- C4: `horizon_lift` `/tmp/` output predates `0e70693` and `b18f25b` edits to `horizon_lift.rs`

- **Claim:** `C4` (`docs/BENCHMARK_RESULTS.md:21,57`) -- "Two-pass graph liftover (`horizon_lift --bridge-bp 10` on 500 `GRCh38.chr22` 1-kb windows -> `CHM13.chr22`): **256.60s**, **619.2 MB** peak RSS (`469` exact, `1` partial, `19` split, `11` unmapped)."
- **Verdict:** `NARROWED`
- **Flaw:** While the `/tmp/` output files match the reported counts, they were generated before commit `0e70693` modified `src/bin/horizon_lift.rs` and `src/gfa/mod.rs`, and neither the source graph, window-generation script, nor RSS/timing log is preserved in the repository.
- **Evidence:**
  - `/tmp/grch38_chr22_windows.bed` (`500` lines, mtime `Nov 19 10:51`), `/tmp/chm13_chr22_lifted_b10.bed` (`514` lines, mtime `Nov 19 11:01`), and `/tmp/chm13_chr22_unmapped_b10.bed` (`11` lines, mtime `Nov 19 11:01`) yield `469 EXACT`, `1 PARTIAL`, `19 SPLIT` (44 fragments), and `11` unmapped.
  - Commit `0e70693` (`2026-11-19 11:20:04 -0500`) subsequently modified `src/bin/horizon_lift.rs` (16 lines) and `src/gfa/mod.rs` (6 lines), and `b18f25b` (`12:27:59 -0500`) never re-ran `horizon_lift`.
  - No script in `scripts/` generates `/tmp/grch38_chr22_windows.bed` or records `256.60s` / `619.2 MB` RSS.
- **Correction:** Re-run `horizon_lift --bridge-bp 10` on `HEAD` (`b18f25b`) and record the exact command, window-generation script, and `/usr/bin/time -v` log.

---

## F4 -- C5: `horizon_clean` benchmarked on `/tmp/` pairwise `horizon_wga` GFA (not HPRC graph) with unlogged runtime/RSS

- **Claim:** `C5` (`docs/BENCHMARK_RESULTS.md:23,59,63`) -- "Real HPRC `chr22` Pericentromere 5 Mbp (`GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M`): **23.45s**, **3,084.5 MB** peak RSS, **5,894,191** aligned pairs, `246` fwd / `184` rev LCBs, `17,789` Horizon Bubbles (`horizon_clean`: `6,030` links kept, **0.01s**, **7.8 MB** RSS)."
- **Verdict:** `NARROWED`
- **Flaw:** The `/tmp/` artifacts match the `horizon_clean` link counts (`6,030` kept), but `horizon_clean` was run on a pairwise `horizon_wga` GFA (`/tmp/chr22_peri_5m.gfa`) rather than the HPRC Minigraph-Cactus graph, and no extraction script or `time`/RSS log exists in the repository.
- **Evidence:**
  - `/tmp/chr22_peri_5m.gfa` (mtime `Nov 19 11:36`) has `6,031` `S` lines, `6,030` `L` lines, and `2,982` `HB:i:1` lines; `/tmp/chr22_peri_5m_clean.gfa` (mtime `Nov 19 11:37`) has `6,030` `L` lines and `0` `HB:i:1` lines.
  - No script in `scripts/` or log in `reports/` documents how `/tmp/grch38_chr22_peri_5m.fa` and `/tmp/chm13_chr22_peri_5m.fa` were sliced or captures the `23.45s` / `3,084.5 MB` / `7.8 MB` measurements.
- **Correction:** Clarify in `docs/BENCHMARK_RESULTS.md:23,63` that `horizon_clean` was benchmarked on the `horizon_wga` pairwise output `/tmp/chr22_peri_5m.gfa` (where `0` links were severed), and archive the benchmark runner/logs in `scripts/benchmarks/`.

---

## F5 -- C6: `horizon_resolve` benchmark is an undisclosed synthetic toy GFA placed under HPRC Pangenome Benchmarks

- **Claim:** `C6` (`docs/BENCHMARK_RESULTS.md:64`) -- "| GFA E-Graph Resolver | `horizon_resolve --collapse-dups` | 200 segs, 20 walks x 225 steps (`4,500` steps) | `280` windows, `840` dups collapsed, `308,760.9` bits saved | **24.0 MB** | 4.48s |"
- **Verdict:** `NARROWED`
- **Flaw:** The input file is an untracked synthetic toy graph in `/tmp/` with identical `ACGT`x50 segments and hardcoded tandem-dup walks, yet it is presented inside `## Pangenome Graph Benchmarks (12 GB hprc-v1.0-mc-grch38.gfa.gz & Subgraphs)` (`docs/BENCHMARK_RESULTS.md:50`) without disclosing that it is synthetic.
- **Evidence:**
  - `/tmp/bench_resolve_10k.gfa` (`Nov 19 11:37`, `61,741` bytes) contains `200` `S` lines (`s0`..`s199`), every one of which is `ACGT` repeated 50 times (`200 bp`), and `20` identical `W` lines (`HG000`..`HG019`) repeating `>s0>s1>s1<s2...`.
  - Meanwhile, on real HPRC `chr22` (`/tmp/hprc_chr22_clean.gfa`, `414,479` bytes), `/tmp/hprc_chr22_resolved.gfa` (`Nov 19 11:27`) has the exact same size (`414,479` bytes), indicating zero collapses on the real graph output.
  - No generator script for `/tmp/bench_resolve_10k.gfa` exists in `scripts/`.
- **Correction:** Move `C6` out of the HPRC pangenome table (or explicitly label it `Synthetic 200-seg/20-walk benchmark`), commit its generator script to `scripts/benchmarks/`, and benchmark `horizon_resolve` on real HPRC `W` lines.

---

## F6 -- C8: Cancelled `chr6` 26.5 GB / 47:33 run has no saved log and its `/tmp/` trace predates `b18f25b`

- **Claim:** `C8` (`docs/BENCHMARK_RESULTS.md:26`) -- "(After replacing `HashSet<(usize, usize)>` with sorted/deduplicated `Vec<(usize, usize)>`, full 170.8 Mbp `GRCh38.chr6` vs 110.7 Mbp `CHM13.chr6` still reached **26.5 GB** RSS and `47:33` CPU before cancellation on the 30 GiB host due to `ContigIndex`'s per-position `HashMap<u64, Vec<usize>>` seed maps, ~0.29 GB/Mbp.)"
- **Verdict:** `NARROWED`
- **Flaw:** While `ContigIndex` (`src/seeding/mod.rs:332-333`) does use `HashMap<u64, Vec<usize>>` on a `31,409 MiB` (`30.67 GiB`) host and `/tmp/hprc_grch38_chr6.fa` / `/tmp/hprc_chm13_chr6.fa` exist, the only `/tmp/` artifacts of the cancelled `chr6` run (`/tmp/hprc_chr6.gfa` and `/tmp/hprc_chr6.maf`, both `0` bytes) have mtime `Nov 19 11:05` -- predating `b18f25b` (`12:27:59`) -- and no RSS or CPU timing log is preserved.
- **Evidence:** `ls -la /tmp/hprc_chr6.*` shows `0` bytes at `Nov 19 11:05`; `free -m` confirms `31409` MiB RAM; no log file in `reports/` or `scripts/` records `26.5 GB` RSS or `47:33` CPU.
- **Correction:** Mark `26.5 GB` RSS / `47:33` CPU as an unlogged observation or re-run with `/usr/bin/time -v` saved to `reports/`.

---

## F7 -- Carried R1-IMPL: `horizon_resolve` and `horizon_wga --graph` still lack HPRC/HLA/centromere benchmarks

- **Claim:** `R1-IMPL` (carried from Round 1) -- "Implement and benchmark horizon_clean, horizon_resolve, horizon_lift, and horizon_wga --graph on HPRC v1.0/v1.1 (GRCh38 vs CHM13) and HLA/centromere loci."
- **Verdict:** `NARROWED`
- **Flaw:** All four tools are now implemented in `src/bin/`, but `R1-IMPL` is not fully discharged on real HPRC/HLA/centromere data: `horizon_resolve` is benchmarked only on synthetic `/tmp/bench_resolve_10k.gfa`, `horizon_wga --graph` is benchmarked only on a 4-node synthetic diamond GFA (`/tmp/test_dag.gfa`), and the 5 Mbp MHC/HLA locus (`docs/BENCHMARK_RESULTS.md:25,61`) uses a synthetic C4 tandem-dup haplotype (`/tmp/hla_dup_5m.fa`) rather than a real HPRC haplotype.
- **Evidence:** `docs/BENCHMARK_RESULTS.md:61,64,65` and `/tmp/bench_resolve_10k.gfa`, `/tmp/test_dag.gfa`, `/tmp/hla_dup_5m.fa`.
- **Correction:** Keep `R1-IMPL` open for real HPRC/HLA/centromere benchmarking of `horizon_resolve` and `horizon_wga --graph` (and `horizon_clean` on real Minigraph-Cactus snarls).

---

## Summary Table

| Claim | Verdict |
| :--- | :--- |
| C1 | SURVIVED |
| C2 | REFUTED |
| C3 | REFUTED |
| C4 | NARROWED |
| C5 | NARROWED |
| C6 | NARROWED |
| C7 | SURVIVED |
| C8 | NARROWED |
| R1-IMPL | NARROWED |
| R1-GFA | SURVIVED |
| R1-HMM | SURVIVED |
