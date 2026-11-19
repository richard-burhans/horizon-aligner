# Round 2 ADJUDICATION of L4 -- Provenance Lens

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `b18f25b444985a5b7f1db38fd9dd312240230839`)  
**Attack file:** `reports/rounds/round-2/attack-L4-provenance.md` (SHA-256 prefix verified: `035443131fdb`)  
**Prior history & method cards:** Checked `./scripts/round history "stream_gfa" --kind adjudication` and `reports/adversarial_log.md` (Round 1 left `R1-IMPL`, `R1-GFA`, `R1-HMM` open). No method cards exist in `refs/methods/`.

---

## F1 -- C2: Stale multi-commit provenance, missing benchmark datasets, and hardcoded `HmmParams::default()` in `--graph`

- **Diagnosis:** **UPHELD.** Verified against `git show b857377`, `git show 7bf688c`, `git show 0e70693`, and `git show b18f25b`: `docs/BENCHMARK_RESULTS.md:5` claims all results were generated with adaptive `identity_prior`, `gap_open_bits`, `gap_extend_bits`, BOCPD profiles, and `src/mdl/dag.rs`, whereas `Chr22 Haplotype Alignment` (`154.92s`, `24,365,611` pairs, `9,690` bubbles) is unchanged from `b857377` (`10:47:45`), `multiz100way`/`cactus447way`/`multiz470way`/`cactus241way (chr22_random)`/`horizon_lift` are unchanged from `7bf688c` (`11:17:25`, before `0e70693` activated the BOCPD jump gate and `b18f25b` made `gap_extend_bits` adaptive), `src/bin/horizon_wga.rs:195,504` passes `&HmmParams::default()` to `align_seq_to_dag_seeded`, and `/home/rcb112/project/gemini/bench/data/` does not exist on disk.
- **Repair:** **NARROWED (split by cost).** The first option—re-running all WGA, `chr22`, and `horizon_lift` benchmarks on `HEAD` (`b18f25b`)—is non-free and currently blocked by the absent `../data/benchmarks` and `hprc-v1.0-mc-grch38.gfa.gz` inputs. The alternative option—tagging each table row/bullet in `docs/BENCHMARK_RESULTS.md` with its exact generating commit and noting at line 5 that `--graph` uses `HmmParams::default()`—is endorsed and **the fix is free**.
- **What survives:** `docs/BENCHMARK_RESULTS.md:5` overstates uniform provenance across results generated on three different commits (`b857377`, `7bf688c`, `b18f25b`), and `horizon_wga --graph` uses `HmmParams::default()` rather than LCB-estimated parameters.

---

## F2 -- C3: Unit test count is 40 (not 39) and 4 of 9 named binaries have 0 unit tests

- **Diagnosis:** **UPHELD.** Verified by running `cargo test --no-fail-fast` on `HEAD` (`b18f25b`): `src/lib.rs` reports `35 passed; 0 failed; 1 ignored`, and the binary crates report `1 passed` each for `horizon_clean`, `horizon_gfa`, `horizon_lift`, `horizon_resolve`, and `horizon_wga`, and `0 passed` for `bench_chr1`, `horizon_eval`, `horizon_guide`, and `horizon_merge` (`40` passing tests total, not `39`, with `4` of `9` binaries having `0` unit tests).
- **Repair:** **ENDORSED.** Updating `docs/STRATEGIC_PLAN.md:3` to report `40 passing tests (35 library + 5 binary tests ...; 1 ignored)` and naming the four binaries with `0` unit tests requires no re-measurement; **the fix is free**.
- **What survives:** `docs/STRATEGIC_PLAN.md:3` undercounts passing tests (`39` vs. `40` + `1` ignored) and implies all named binaries are unit-tested when four (`bench_chr1`, `horizon_eval`, `horizon_guide`, `horizon_merge`) have zero unit tests.

---

## F3 -- C4: `horizon_lift` `/tmp/` output predates `0e70693` and `b18f25b` edits to `horizon_lift.rs`

- **Diagnosis:** **UPHELD (attacker's `NARROWED` ruling stands).** Verified `/tmp/grch38_chr22_windows.bed` (`500` lines, `Nov 19 10:51`), `/tmp/chm13_chr22_lifted_b10.bed` (`514` lines, `Nov 19 11:01`), and `/tmp/chm13_chr22_unmapped_b10.bed` (`11` lines, `Nov 19 11:01`), which match the reported counts, while `git show 0e70693 --stat` confirms `src/bin/horizon_lift.rs` (16 lines) and `src/gfa/mod.rs` (6 lines) were subsequently modified at `11:20:04` and no window-generation script or `time`/RSS log is in the repository.
- **Repair:** **NARROWED by cost.** Re-running `horizon_lift --bridge-bp 10` on `HEAD` (`b18f25b`) with `/usr/bin/time -v` is non-free (~257s runtime and requires the absent HPRC GFA), whereas annotating the `docs/BENCHMARK_RESULTS.md` entry with its pre-`0e70693` provenance and unlogged RSS/time status is free.
- **What survives:** The `horizon_lift` liftover counts match the `/tmp/` BED outputs, but those outputs predate `0e70693` edits to `horizon_lift.rs` and lack a committed generator script or RSS/timing log.

---

## F4 -- C5: `horizon_clean` benchmarked on `/tmp/` pairwise `horizon_wga` GFA (not HPRC graph) with unlogged runtime/RSS

- **Diagnosis:** **UPHELD (attacker's `NARROWED` ruling stands).** Verified `/tmp/chr22_peri_5m.gfa` (`6,031` `S`, `6,030` `L`, `2,982` `HB:i:1`) and `/tmp/chr22_peri_5m_clean.gfa` (`6,030` `L`, `0` `HB:i:1`). While `docs/BENCHMARK_RESULTS.md:63` already names `/tmp/chr22_peri_5m.gfa`, line 23 (quoted in `C5`) nests `horizon_clean` inside the `Real HPRC chr22 Pericentromere 5 Mbp` bullet without clarifying that its input was the pairwise `horizon_wga` output rather than the HPRC graph, and no extraction script or `time`/RSS log exists in the repo.
- **Repair:** **NARROWED (split by cost).** Clarifying in `docs/BENCHMARK_RESULTS.md:23` that `horizon_clean` ran on the pairwise `horizon_wga` GFA `/tmp/chr22_peri_5m.gfa` (`0` links severed) is endorsed and **the fix is free**; archiving benchmark runner scripts and logs requires re-running and is non-free.
- **What survives:** The `horizon_clean` link counts match `/tmp/chr22_peri_5m_clean.gfa`, but the benchmarked input was a pairwise `horizon_wga` GFA rather than the HPRC Minigraph-Cactus graph, and its timing/RSS are unlogged.

---

## F5 -- C6: `horizon_resolve` benchmark is an undisclosed synthetic toy GFA placed under HPRC Pangenome Benchmarks

- **Diagnosis:** **UPHELD (attacker's `NARROWED` ruling stands).** Verified `/tmp/bench_resolve_10k.gfa` contains `200` `S` lines (`s0`..`s199`) all equal to `ACGT` repeated 50 times and `20` identical `W` lines (`HG000`..`HG019`), while `/tmp/hprc_chr22_clean.gfa` and `/tmp/hprc_chr22_resolved.gfa` are identical in byte size (`414,479` bytes). Presenting `C6` inside `## Pangenome Graph Benchmarks (12 GB hprc-v1.0-mc-grch38.gfa.gz & Subgraphs)` without a synthetic label misleads provenance.
- **Repair:** **NARROWED (split by cost).** Labeling `C6` explicitly as a synthetic 200-segment/20-walk benchmark (or moving it out of the HPRC table) in `docs/BENCHMARK_RESULTS.md:64` is endorsed and **the fix is free**; committing a generator script and benchmarking `horizon_resolve` on real HPRC `W` lines is non-free.
- **What survives:** `horizon_resolve --collapse-dups` was benchmarked only on an undisclosed synthetic toy GFA (`/tmp/bench_resolve_10k.gfa`), not on real HPRC walks.

---

## F6 -- C8: Cancelled `chr6` 26.5 GB / 47:33 run has no saved log and its `/tmp/` trace predates `b18f25b`

- **Diagnosis:** **UPHELD (attacker's `NARROWED` ruling stands).** Verified `src/seeding/mod.rs:332-333` uses `HashMap<u64, Vec<usize>>`, and `/tmp/hprc_chr6.gfa` / `/tmp/hprc_chr6.maf` are `0` bytes with mtime `Nov 19 11:05` (predating `b18f25b` at `12:27:59`), with no RSS or CPU timing log preserved in `reports/` or `scripts/`.
- **Repair:** **NARROWED (split by cost).** Marking `26.5 GB` RSS / `47:33` CPU in `docs/BENCHMARK_RESULTS.md:26` as an unlogged observation is endorsed and **the fix is free**; re-running full `chr6` with `/usr/bin/time -v` is non-free (~48 min CPU).
- **What survives:** The `ContigIndex` `HashMap` cause is real in `src/seeding/mod.rs:332-333`, but the reported `26.5 GB` RSS and `47:33` CPU figures have no saved log and only a `0`-byte `/tmp/` trace predating `b18f25b`.

---

## F7 -- Carried R1-IMPL: `horizon_resolve` and `horizon_wga --graph` still lack HPRC/HLA/centromere benchmarks

- **Diagnosis:** **UPHELD (attacker's `NARROWED` ruling stands).** All four tools are now implemented in `src/bin/`, and `horizon_clean` was run on pairwise pericentromere GFAs, but `horizon_resolve` (`bench_resolve_10k.gfa`), `horizon_wga --graph` (`test_dag.gfa`, 4-node diamond), and the 5 Mbp MHC/HLA run (`hla_dup_5m.fa`, synthetic C4 tandem-dup) remain synthetic rather than real HPRC/HLA/centromere graph benchmarks.
- **Repair:** **ENDORSED.** Keeping `R1-IMPL` open in the round record for real HPRC/HLA/centromere benchmarking of `horizon_resolve`, `horizon_wga --graph`, and `horizon_clean` on real Minigraph-Cactus snarls is **free** as a status action (executing those benchmarks is non-free).
- **What survives:** `R1-IMPL` is partially discharged by the binary implementations, but remains open for real HPRC/HLA/centromere benchmarks of `horizon_resolve`, `horizon_wga --graph`, and `horizon_clean` on real Minigraph-Cactus snarls.

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `docs/BENCHMARK_RESULTS.md:5` claims uniform provenance across three commits (`b857377`, `7bf688c`, `b18f25b`), `--graph` passes `HmmParams::default()`, and `../data/benchmarks` is absent; tagging rows with commits and noting `--graph` defaults is free |
| F2 | UPHELD | `cargo test` passes `40` tests (`35` lib + `5` bin, `1` ignored) rather than `39`, and `4` of `9` binaries have `0` unit tests; updating `docs/STRATEGIC_PLAN.md:3` is free |
| F3 | UPHELD | `/tmp/` BED counts match `C4`, but outputs (`11:01`) predate `0e70693` (`11:20`) edits to `horizon_lift.rs` and lack a committed script or RSS/time log |
| F4 | UPHELD | `/tmp/chr22_peri_5m.gfa` matches `horizon_clean` counts, but it is a pairwise `horizon_wga` GFA rather than the HPRC graph and lacks saved logs; clarifying line 23 is free |
| F5 | UPHELD | `/tmp/bench_resolve_10k.gfa` is an undisclosed synthetic toy GFA (`ACGT`x50, 20 identical walks) placed under the HPRC benchmark table; labeling it synthetic is free |
| F6 | UPHELD | `ContigIndex` uses `HashMap<u64, Vec<usize>>`, but `26.5 GB` RSS / `47:33` CPU has no saved log and its 0-byte `/tmp/` trace (`11:05`) predates `b18f25b`; marking it unlogged is free |
| F7 | UPHELD | All four binaries exist, but `horizon_resolve`, `horizon_wga --graph`, and MHC/HLA rely on synthetic inputs, so `R1-IMPL` remains open for real HPRC/HLA/centromere benchmarks |
