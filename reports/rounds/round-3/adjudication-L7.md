# Round 3 ADJUDICATION of L7 (Reproduction) -- FlatSeedIndex, Adaptive X-Drop HMM, GFA Resolver/Cleaner & Review Fixes

Adjudicating `reports/rounds/round-3/attack-L7-reproduction.md` (sha256 prefix `13037802171d`, verified) against `docs/STRATEGIC_PLAN.md`, `docs/BENCHMARK_RESULTS.md`, `src/`, `scripts/benchmarks/`, and `/tmp/` artifacts at commit `488b020` (`HEAD`).

Independent verification executed:
- `cargo test`: exited `0` (`37 passed; 0 failed; 1 ignored` in `src/lib.rs`, plus `1 passed` in each of the `9` binaries in `src/bin/` = `46 passed; 0 failed; 1 ignored`; `cargo test -- --list | grep -c ": test$"` = `47`).
- `cargo clippy --all-targets -- -D warnings`: exited `0` with no warnings.
- `python3 scripts/benchmarks/extract_fasta_from_maf.py /home/rcb112/project/gemini/bench/pg/data/benchmarks/cactus241way_chr22_1MB.maf hg38 Canis_lupus_familiaris ...` + `./target/release/horizon_eval`: extracted `1 ref contigs and 31 target contigs` and reproduced `elapsed=15.636s proto_pairs=136348 maf_pairs=132089 overlap=82228 raw_maf_recall=62.3% raw_maf_precision=60.3% bubbles=3490`.
- `grep -c "^>" /tmp/dog_zoo_1m.fa /tmp/repro_dog_zoo_1m.fa`: both contain `31` contigs.
- `head -n 1 /tmp/{grch38,chm13}_{chr6_q20m,chr6_peri_5m,chr22_peri_5m}.fa`: verified all 6 headers; `/tmp/chm13_chr6_q20m.fa` is `>CHM13.chr6_61M_81M` while the other 5 match `docs/BENCHMARK_RESULTS.md:23-26,59-62`.
- `printf "H\tVN:Z:1.0\nS\ts1\tACGT\nS\ts2\tTGCA\nL\ts1\t+\ts2\t+\t0M\n" | ./target/release/horizon_gfa stats -`: exited `0` (`segments=2 segment_bp=8 links=1`).

Prior rounds & method cards:
- Prior round history: `scripts/round history` confirms Round 2 F1 already flagged missing `data/benchmarks/` and `run_all_evals.py` omitting `cactus241way` (leading to the Round 3 staleness disclosure in `docs/BENCHMARK_RESULTS.md:5`), Round 2 F2 flagged `39` vs `40` tests, and Round 2 F4 flagged missing slice-extraction scripts for `/tmp/chr22_peri_5m.*`. The `129 target scaffolds` vs `31` contigs error and `CHM13.chr6:61.4M-81.4M` vs `61M-81M` mismatch are new in Round 3.
- Method cards: `refs/methods/` does not exist in this repository; F1–F5 are reproduction, artifact, and codebase-consistency checks rather than statistical disputes, so no method card governs them.

---

## Rulings

### F1 -- C2 claims `129 target scaffolds` when the extracted Dog FASTA contains `31` contigs, and `cactus241way_chr22_1MB.maf` is absent from benchmark scripts
- **Verdict:** **UPHELD**
- **What survives:** The attacker's `NARROWED` ruling on `C2` holds: `horizon_eval` reproduces `62.3%` recall, `60.3%` precision, `61.3%` F1, and `3,490` bubbles (`15.636s` here vs `17.16s` claimed), but both `/tmp/dog_zoo_1m.fa` and a fresh extraction from `cactus241way_chr22_1MB.maf` contain `31` Dog contigs rather than `129 target scaffolds`, and neither `scripts/benchmarks/fetch_benchmarks.sh` nor `scripts/benchmarks/run_all_evals.py` fetches or runs `cactus241way_chr22_1MB.maf` (only `/home/rcb112/project/gemini/bench/run_zoo_1mb.py` references it with a stale path).
- **Repair ruling:** Endorsed in part. Changing `across 129 target scaffolds` to `across 31 target contigs` in `docs/BENCHMARK_RESULTS.md:16` is **free** (doc-only edit); adding `cactus241way_chr22_1MB.maf` extraction/evaluation to `fetch_benchmarks.sh` and `run_all_evals.py` costs script authoring.

### F2 -- C3 miscounts the passing unit test suite as `45` instead of `46`
- **Verdict:** **UPHELD**
- **What survives:** `cargo clippy --all-targets -- -D warnings` is clean and all `9` binaries in `src/bin/` have `1` passing unit test each, but `cargo test` runs `37` passing + `1` ignored in `src/lib.rs` plus `9` binary tests = `46` passing + `1` ignored (`47` total), not `45` passing + `1` ignored.
- **Repair ruling:** Endorsed. Updating `docs/STRATEGIC_PLAN.md:3` from `` `45` passing tests + `1` ignored `` to `` `46` passing tests + `1` ignored `` is **free** (doc-only edit).

### F3 -- C4 and C5 rely on untracked `/tmp/*.fa` slices with no extraction script, and C4 misstates the `CHM13.chr6` coordinate window
- **Verdict:** **NARROWED**
- **What survives:** The attacker's reproducibility gap holds—`grep -rn "chr6_q20m\|chr22_peri_5m\|chr6_peri_5m" src/ scripts/` returns no matches, so all 6 `/tmp/*.fa` inputs lack an in-repo extraction script—and `/tmp/chm13_chr6_q20m.fa` has header `>CHM13.chr6_61M_81M` whereas `docs/BENCHMARK_RESULTS.md:26,62` reports `CHM13.chr6:61.4M-81.4M` (all 5 other `/tmp` slice headers match `docs/BENCHMARK_RESULTS.md:23-26,59-62`). Narrowed on two points: (1) `C4` (`docs/STRATEGIC_PLAN.md:11`) does not state the `CHM13.chr6` coordinate window (the misstatement is in `docs/BENCHMARK_RESULTS.md:26,62`), and (2) `C4`'s `12 bytes/entry` accurately describes the stored `FlatSeedIndex` representation (`Vec<u64>` + `Vec<u32>`), even though `FlatSeedIndex::build` (`src/seeding/mod.rs:386-441`) temporarily holds `pairs: Vec<(u64, u32)>` (~28 bytes/entry peak) before `drop(pairs)`.
- **Repair ruling:** Endorsed in part. Fixing `CHM13.chr6:61.4M-81.4M` to `CHM13.chr6:61M-81M` in `docs/BENCHMARK_RESULTS.md:26,62` is **free** (doc-only edit); adding a slice-extraction script under `scripts/benchmarks/` costs script authoring.

### F4 -- `C1` and `C6` survive direct reproduction
- **Verdict:** **UPHELD**
- **What survives:** `C1` (`src/gfa/mod.rs:97-127,230-329`, both `stream_gfa` unit tests, and `horizon_gfa stats -` `stdin` streaming) and `C6` (`horizon_gfa stats` exact counts `81,415,956` segments / `3.29 Gbp` / `112,955,105` links / `2,580` paths / `24,456` walks at `96.7 MB` vs `97.1 MB` RSS) survive direct reproduction.

### F5 -- Framing & steering: `run_all_evals.py` points to a nonexistent data directory, and `C2` omits the 5 stale pre-`488b020` WGA benchmarks
- **Verdict:** **NARROWED**
- **What survives:** `scripts/benchmarks/run_all_evals.py:10` sets `BENCH_DIR = ROOT.parent / "data" / "benchmarks"` (`/home/rcb112/project/gemini/bench/data/benchmarks`, nonexistent), while the downloaded benchmarks reside in `/home/rcb112/project/gemini/bench/pg/data/benchmarks`. Narrowed on the steering leg: although the Round 3 `C2` claim selected only the re-run `cactus241way chr22 1MB` row, `docs/BENCHMARK_RESULTS.md:5` already explicitly discloses the 5 earlier pre-`488b020` WGA rows following Round 2 F1.
- **Repair ruling:** Endorsed on the `run_all_evals.py` path defect; fixing `BENCH_DIR` (or making it configurable via CLI/env) costs a small script edit (not doc-only).

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `C2` metrics reproduce, but extracted Dog FASTA has `31` contigs not `129 target scaffolds` and benchmark scripts omit `cactus241way_chr22_1MB.maf`; doc fix is free, script addition costs authoring |
| F2 | UPHELD | `cargo test` runs `46` passing + `1` ignored (`37` lib + `9` bin), not `45` + `1`; doc fix in `docs/STRATEGIC_PLAN.md:3` is free |
| F3 | NARROWED | Missing `/tmp/*.fa` extraction scripts and `CHM13.chr6:61.4M-81.4M` -> `61M-81M` error in `docs/BENCHMARK_RESULTS.md:26,62` hold (other 5 slices match), but `C4` itself omits the coordinate window and `12 bytes/entry` accurately describes the stored struct despite ~28 bytes/entry build peak; doc fix is free |
| F4 | UPHELD | `C1` (`stream_gfa` `.gfa`/`.gfa.gz`/`-` `stdin` + unit tests) and `C6` (`horizon_gfa stats` counts and ~97 MB RSS) survive direct reproduction |
| F5 | NARROWED | `run_all_evals.py:10` points to nonexistent `data/benchmarks` instead of `pg/data/benchmarks` (non-free script fix), while `docs/BENCHMARK_RESULTS.md:5` already discloses the 5 pre-`488b020` WGA rows |
