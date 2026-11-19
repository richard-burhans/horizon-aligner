# Round 2 ADJUDICATION of L7 (Reproduction) -- Horizon Pangenome Suite, Adaptive Pair-HMM, E-Graph Resolver, and HPRC Benchmarks

Adjudicating `reports/rounds/round-2/attack-L7-reproduction.md` (sha256 prefix `e28e4b5da36f`, verified) against `docs/STRATEGIC_PLAN.md`, `docs/BENCHMARK_RESULTS.md`, `docs/PANGENOME_SUITE.md`, `src/`, `scripts/benchmarks/`, and `/tmp/` artifacts at commit `b18f25b`.

Independent verification executed:
- `/home/agent/.cargo/bin/cargo test`: exited 0 (`35 passed; 0 failed; 1 ignored` in `src/lib.rs`; `1 passed` each in `horizon_clean`, `horizon_gfa`, `horizon_lift`, `horizon_resolve`, `horizon_wga`; `0` tests in `bench_chr1`, `horizon_eval`, `horizon_guide`, `horizon_merge`; total = **40 passing**, `1` ignored).
- `git diff 7bf688c..b18f25b -- docs/BENCHMARK_RESULTS.md`: confirmed `cactus241way (chr22 1MB)` changed (`64.3%` recall / `3,362` bubbles -> `63.6%` recall / `3,411` bubbles), while 7 earlier alignment rows remained unchanged.
- `git show 0e70693 -- src/bin/horizon_lift.rs`: confirmed `0e70693` only removed unused `src_total_len` / `_src_iv` and replaced `.unwrap()` with `.expect()` in tests—it did **not** modify GFA parsing in `horizon_lift.rs`.
- `./target/release/horizon_clean --gfa /tmp/chr22_peri_5m.gfa -o /tmp/adj_test_clean.gfa`: reproduced `total_links=6030 kept=6030 severed=0 bubbles=0 split_walks=0 elapsed=0.02s peak_rss_mb=8.1`, and verified all `6,031` `S`-lines in `/tmp/chr22_peri_5m.gfa` have `*` sequences.
- `./target/release/horizon_resolve --collapse-dups /tmp/bench_resolve_10k.gfa -o /tmp/adj_test_resolve.gfa`: reproduced `segments=200 anchors=180 paths=0 walks=20 windows=280 simplified=280 dups_detected=840 dups_collapsed=840 invs_cancelled=0 mdl_bits_saved=308760.9 elapsed_s=10.886 peak_rss_mb=24.3`.
- `/tmp/chm13_chr22_lifted_b10.bed` (`514` rows, `44` `SPLIT` rows across `19` unique `win_*` IDs) and `/tmp/chm13_chr22_unmapped_b10.bed` (`11` rows) verified against `/tmp/grch38_chr22_windows.bed` (`500` rows).

Prior rounds & method cards:
- Prior round history: `scripts/round history "" --kind adjudication` shows Round 1 adjudicated `stream_gfa`, `horizon_clean`, `horizon_resolve`, and `horizon_lift` before they were implemented in `b857377`–`b18f25b` (disclosed by the attacker).
- Method cards: `refs/methods/` does not exist in this repository; F1–F6 are codebase, artifact, and benchmark-provenance reproduction checks rather than statistical disputes, so no method card governs them.

---

## Rulings

### F1 -- C2 ("All results were generated out-of-the-box...")
- **Verdict:** **UPHELD**
- **What survives:** `docs/BENCHMARK_RESULTS.md:5` claims *"All results were generated out-of-the-box"* using adaptive `gap_extend_bits` (`b18f25b`) and normalized/jump-gated Pair-HMM transitions (`0e70693`), yet `git diff 7bf688c..b18f25b` proves only `cactus241way (chr22 1MB)` and the 5 Mbp/20 Mbp HPRC and synthetic `horizon_resolve` rows were re-run after those behavioral changes, leaving at least 7 earlier alignment rows stale, while `/home/rcb112/project/gemini/bench/data/benchmarks/` is absent from the workspace and `scripts/benchmarks/run_all_evals.py` omits `cactus241way`, HPRC, and synthetic benchmarks. (Minor qualification on leg 2: `scripts/benchmarks/fetch_benchmarks.sh` *does* exist to download the UCSC and HPRC source files, though it has not been populated into `data/benchmarks/`.)
- **Repair ruling:** Endorsed. Qualifying `docs/BENCHMARK_RESULTS.md:5` to identify which rows predate `0e70693`/`b18f25b` is **free** (doc-only edit); re-running all benchmarks on `HEAD` and expanding `run_all_evals.py` costs compute/download time and script authoring.

### F2 -- C3 (`39` passing tests; unit-tested binaries)
- **Verdict:** **UPHELD**
- **What survives:** `cargo test` on `b18f25b` runs **40 passing tests** (`35` in `src/lib.rs` with `1` ignored + `1` each in `horizon_clean`, `horizon_gfa`, `horizon_lift`, `horizon_resolve`, `horizon_wga`), not `39`, and `4` of the `9` binaries in `src/bin/` (`bench_chr1`, `horizon_eval`, `horizon_guide`, `horizon_merge`) have `0` unit tests.
- **Repair ruling:** Endorsed. Updating `docs/STRATEGIC_PLAN.md:3` to `40 passing tests (1 ignored)` and noting the `4` untested binaries is **free** (doc-only edit).

### F3 -- C4 (`horizon_lift --bridge-bp 10` benchmark)
- **Verdict:** **NARROWED**
- **What survives:** The attacker's `NARROWED` reproducibility gap survives on missing repository inputs (`hprc-v1.0-mc-grch38.gfa.gz` is absent on disk, though `fetch_benchmarks.sh` has its URL, and no script generates `/tmp/grch38_chr22_windows.bed`), whereas the attacker's code-staleness premise that `0e70693` modified `src/bin/horizon_lift.rs` to deduplicate GFA parsing via `stream_gfa` is **overturned** (`git show 0e70693 -- src/bin/horizon_lift.rs` shows only unused `src_total_len`/`_src_iv` removal and `.unwrap()` -> `.expect()` in tests).
- **Repair ruling:** Endorsed in part. Documenting the provenance/generation command for the 500 `GRCh38.chr22` windows is **free** if known; adding a generator script and re-running `horizon_lift` on `HEAD` requires downloading the 12 GB GFA and ~257s compute.

### F4 -- C5 (`chr22` pericentromere `horizon_wga` and `horizon_clean`)
- **Verdict:** **UPHELD**
- **What survives:** The attacker's `NARROWED` verdict holds: `horizon_clean` reproduces `6,030` kept links on `/tmp/chr22_peri_5m.gfa` (and all `6,031` `S`-lines have `*` sequences, so flank scoring cannot sever any link), while the 5 Mbp pericentromere FASTAs and `/tmp/chr22_peri_5m.gfa` lack any generator script in the repository.
- **Repair ruling:** Endorsed. Noting in `docs/BENCHMARK_RESULTS.md` that `/tmp/chr22_peri_5m.gfa` has `*`-sequence segments is **free** (doc-only edit); adding a slice-extraction script to `scripts/benchmarks/` requires script authoring.

### F5 -- C6 (`horizon_resolve --collapse-dups`)
- **Verdict:** **UPHELD**
- **What survives:** The attacker's `NARROWED` verdict holds: `./target/release/horizon_resolve --collapse-dups /tmp/bench_resolve_10k.gfa` reproduces the exact counts (`200` segments, `20` walks, `280` windows, `840` dups collapsed, `308,760.9` bits saved), but `/tmp/bench_resolve_10k.gfa` has no generator script or fixture in the repository.
- **Repair ruling:** Endorsed. Committing the generator script or fixture into the repository requires adding the file/script (not doc-only).

### F6 -- C8 (`ContigIndex` RAM scaling / cancelled `chr6` OOM run)
- **Verdict:** **UPHELD**
- **What survives:** The attacker's `NARROWED` verdict holds: `ContigIndex` (`src/seeding/mod.rs:329-333`) and `Vec<(usize, usize)>` (`src/mdl/mod.rs:456`) match the code description and `/tmp/hprc_chr6.gfa`/`.maf` are `0`-byte files, but no log, script, or RSS trace for the `26.5 GB` / `47:33` CPU run exists in the repository, and the `~0.29 GB/Mbp` rate lacks a documented denominator.
- **Repair ruling:** Endorsed. Recording the command, RSS measurement method, and Mbp denominator (or removing the ungrounded `~0.29 GB/Mbp` figure) in `docs/BENCHMARK_RESULTS.md` is **free** (doc-only edit).

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `0e70693`/`b18f25b` altered alignment results (`cactus241way chr22 1MB`), leaving 7+ earlier rows stale under the "All results" header, and `data/benchmarks/` + runner coverage for `cactus241way`/HPRC/synthetic are absent; doc qualification is free, full re-run costs compute |
| F2 | UPHELD | `cargo test` runs 40 passing tests (1 ignored), not 39, and 4 of 9 binaries in `src/bin/` have 0 unit tests; doc fix is free |
| F3 | NARROWED | Missing 12 GB GFA on disk and missing window-generation script hold, but the premise that `0e70693` changed `horizon_lift.rs` GFA parsing is overturned (only dead-arg/test cleanup); doc provenance note is free |
| F4 | UPHELD | Attacker's NARROWED verdict holds: `horizon_clean` reproduces `6,030` kept links on `*`-sequence `/tmp/chr22_peri_5m.gfa` (trivially keeping all links), and no slice-extraction script exists in the repo; `*`-sequence doc note is free |
| F5 | UPHELD | Attacker's NARROWED verdict holds: `horizon_resolve` reproduces exact counts on `/tmp/bench_resolve_10k.gfa`, but no generator script or fixture exists in the repo; adding fixture/script is non-free |
| F6 | UPHELD | Attacker's NARROWED verdict holds: `ContigIndex` structure matches code, but no log/trace exists for the cancelled `26.5 GB` / `47:33` `chr6` run and `~0.29 GB/Mbp` lacks a denominator; doc fix is free |
