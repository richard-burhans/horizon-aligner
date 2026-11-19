# Round 4 Adversarial Attack -- L1 (Correctness)

**Target:** `reports/rounds/round-4/ATTACKED-horizon_paper.md` (`manuscript/horizon_paper.md`) against live codebase at `HEAD` (`33bede8`) and `docs/BENCHMARK_RESULTS.md`.

**Rule 2 / Rule 3 disclosures:**
- **Literature (Rule 2):** No `../library` checkout was accessed; this correctness attack evaluates claims against the live repository (`src/`, `docs/BENCHMARK_RESULTS.md`, `Cargo.toml`, `scripts/`, `data/`) without asserting prior-art claims from memory.
- **Method cards (Rule 3):** Per `reports/rounds/round-4/BRIEF-correctness.md`, this project keeps no method cards yet.

---

## F1 -- C1: Binary count is 10 (not 9), `horizon` shares no core, and 3 tools use no information-theoretic core

- **Claim:** `C1` (`manuscript/horizon_paper.md:18`)
- **Verdict:** `NARROWED`
- **Flaw:** `C1` describes Horizon as *"a unified, pure-Rust suite of nine composable executables sharing one information-theoretic core."* Three concrete mismatches exist in the live codebase:
  1. `src/bin/` contains **10** executables (`bench_chr1.rs`, `horizon.rs`, `horizon_clean.rs`, `horizon_eval.rs`, `horizon_gfa.rs`, `horizon_guide.rs`, `horizon_lift.rs`, `horizon_merge.rs`, `horizon_resolve.rs`, `horizon_wga.rs`), and `manuscript/horizon_paper.md:101` itself states *"all `10` CLI binaries"*.
  2. `src/bin/horizon.rs:7-53` is a `std::process::Command` wrapper dispatching 9 subcommands and imports nothing from `horizon_aligner`.
  3. `src/bin/horizon_gfa.rs:4` imports only `horizon_aligner::gfa::*`; `src/bin/horizon_lift.rs` imports only `coord::Interval` and `gfa::*`; `src/bin/horizon_merge.rs` imports only `merge_pairwise_mafs`. None of these three executables invoke BOCPD, MDL DP, Conformal Risk Control, or e-graphs.
- **Evidence:**
  - `ls -la src/bin/` (10 `.rs` files) and `/home/agent/.cargo/bin/cargo test --no-fail-fast` (10 binary test harnesses, 59 active tests passing: 45 lib + 14 across 10 binaries, 1 ignored).
  - `head -n 20 src/bin/*.rs | grep -E "==>|use horizon_aligner"`
- **Prior-round adjudication:** Round 1 (`C1`, `NARROWED`) and Round 3 (`F3`, `NARROWED`) previously noted binary-count mismatches. Commit `f777935` subsequently added `src/bin/horizon.rs` (bringing `src/bin/` to 10 binaries), creating an internal contradiction between `manuscript/horizon_paper.md:18` ("nine") and `manuscript/horizon_paper.md:101` ("10 CLI binaries").
- **Correction:** Update `manuscript/horizon_paper.md:18` to: *"a unified, pure-Rust suite of ten executables (a `horizon` dispatcher over nine specialized tools), with alignment, cleaning, resolution, and progressive modules sharing the `bocpd`/`mdl`/`egraph` core and GFA/liftover/merge modules sharing the streaming `gfa`/`coord` core."*

---

## F2 -- C2: Full `chr6` `451.35 s` (`2.35×` speedup) is pinned to `53130b8` (6 commits behind `HEAD`) and compares non-identical workloads

- **Claim:** `C2` (`manuscript/horizon_paper.md:101, 105`)
- **Verdict:** `NARROWED`
- **Flaw:** The `451.35 s` (`10,515.0 MB` RSS, `2.35×` total speedup, `12.3×` `LASTZ12` speedup) measurement was recorded on commit `53130b8` with `HORIZON_PROFILE=1`, which is 6 commits behind `HEAD` (`33bede8`). Between `53130b8` and `33bede8`, 21 files in `src/` were modified (`951` insertions, `523` deletions), including `src/seeding/flat_index.rs`, `src/seeding/mod.rs`, `src/seeding/translated.rs`, `src/mdl/dp.rs`, and `src/mdl/mod.rs`. Moreover, `1e47fe0` (`1,062.45 s`) and `53130b8` (`451.35 s`) did not execute the same workload (`292,534` fwd / `291,150` rev LCBs and `693,242` DP jobs vs. `301,156` fwd / `298,178` rev LCBs and `778,471` DP jobs), and `seed_end` also dropped `2.03×` (`184.65 s` → `91.04 s`).
- **Evidence:**
  - `docs/BENCHMARK_RESULTS.md:27, 64-65`
  - `git diff --stat 53130b8..HEAD -- src/` (21 files changed)
- **Correction:** Pin the `451.35 s` / `10,515.0 MB` measurement in `manuscript/horizon_paper.md:101` to commit `53130b8` (`HORIZON_PROFILE=1`) or re-run full `GRCh38.chr6` × `CHM13.chr6` on `HEAD` (`33bede8`).

---

## F3 -- C3: HPRC `chr6` `horizon_gfa extract` (`244.93 s`) is pinned to stale commit `b857377`, and `hprc-v1.0-mc-grch38.gfa.gz` is absent on disk

- **Claim:** `C3` (`manuscript/horizon_paper.md:101`)
- **Verdict:** `NARROWED`
- **Flaw:** `docs/BENCHMARK_RESULTS.md:57` pins `horizon_gfa extract` for `GRCh38.chr6` + `CHM13.chr6` (`244.93 s`, `1,255.1 MB` RSS) to commit `b857377`—11 commits behind `HEAD` (`33bede8`). Six subsequent commits (`7bf688c`, `0e70693`, `488b020`, `7855881`, `8c7d929`, `b068896`, `33bede8`) modified `src/bin/horizon_gfa.rs` and `src/gfa/mod.rs`. Additionally, `hprc-v1.0-mc-grch38.gfa.gz` is not present in `data/` on this host, so none of the `12 GB` HPRC timings (`99.97 s`, `183.29 s`, `244.93 s`) can be re-verified locally. (`horizon_gfa stats` in `src/bin/horizon_gfa.rs:39-77` is verified `O(1)`-memory via 5 `u64` counters.)
- **Evidence:**
  - `docs/BENCHMARK_RESULTS.md:57`
  - `git log --oneline -n 10 -- src/bin/horizon_gfa.rs src/gfa/mod.rs`
  - `ls -lh data/` (empty)
- **Prior-round adjudication:** Round 3 (`C6`, correctness) marked the HPRC graph benchmark `UNCHALLENGED` due to `hprc-v1.0-mc-grch38.gfa.gz` being absent on disk; `docs/BENCHMARK_RESULTS.md:57` remains pinned to `b857377`.
- **Correction:** Pin the `244.93 s` / `1,255.1 MB` `chr6` `extract` result in `manuscript/horizon_paper.md:101` to `b857377` or re-run on `33bede8` where the HPRC archive is available.

---

## F4 -- C4: `horizon_clean` severs `walk_support <= 1` (not only `walk_support == 0`), is a no-op on `*`-sequence GFAs, and sums `W` lengths without overlap adjustment

- **Claim:** `C4` (`manuscript/horizon_paper.md:49`)
- **Verdict:** `NARROWED`
- **Flaw:**
  1. `src/bin/horizon_clean.rs:260` computes `let sever = (cfg.prune_unsupported && sup == 0) || (sup <= 1 && cp >= cfg.min_cp_bits);`. It therefore severs `walk_support == 1` edges (single-walk junctions) when `cp >= cfg.min_cp_bits`, not solely `walk_support == 0` edges as `C4` states.
  2. On `*`-sequence GFAs (including `horizon_wga`'s own `S\t...\t*` output), `has_seq` is `false` (`src/bin/horizon_clean.rs:154-176`), so `score_junction_cp_bits` returns `0.0` for all links (`docs/BENCHMARK_RESULTS.md:66`: `2,776 kept, 0 severed (no-op on *-seq)`).
  3. When splitting `W` lines across severed links (`src/bin/horizon_clean.rs:418-420`), `cur_end` sums `seg_flanks` lengths without subtracting any non-zero `L`-line overlap CIGARs.
- **Evidence:**
  - `src/bin/horizon_clean.rs:154-176, 260, 418-420`
  - `docs/BENCHMARK_RESULTS.md:66`
- **Prior-round adjudication:** Round 2 (`R2-RESOLVE-LIFT`) and Round 3 (`adjudication-L5.md`, `adjudication-completeness.md`) previously noted the `*`-seq no-op and `GfaLink` tag limitation; this finding adds the `sup <= 1` vs `walk_support == 0` discrepancy in `src/bin/horizon_clean.rs:260` against `manuscript/horizon_paper.md:49`.
- **Correction:** Update `manuscript/horizon_paper.md:49` to state `walk_support <= 1` (and note that sequence-less `*` segments yield `0.0` changepoint bits unless `--prune-unsupported` is enabled for `walk_support == 0`).

---

## F5 -- C5: *Plasmodium* (`87.7%`, `54.4%`, `66.8%`) and *Cannabis sativa* (`418,663` / `418,465`, `99.95%`, `1.000`) benchmarks are absent from `docs/BENCHMARK_RESULTS.md`, `scripts/`, and `data/`

- **Claim:** `C5` (`manuscript/horizon_paper.md:65`)
- **Verdict:** `REFUTED`
- **Flaw:** Neither the *P. falciparum* × *P. vivax* (`87.7%` abstention, `54.4%` / `66.8%` recovery) nor the *Cannabis sativa* (`cs10` × `Cherry_Pie`, `418,663` columns, `418,465` retained, `99.95%`, `1.000` recall) benchmark exists in `docs/BENCHMARK_RESULTS.md`, `scripts/`, `src/`, or `data/`. Furthermore, `manuscript/horizon_paper.md:65` contrasts `54.4%` against `66.8%` *"with `w = 10` sparse-exon fallback and `band_w = 180`"*, even though `band_w: 180` (`src/mdl/mod.rs:67`) and `w = 10` fallback (`fallback_min_anchors: 4`, `fallback_min_span: 50`) are already the defaults in `AlignConfig::default()`. Given major algorithmic changes up to `33bede8` (`HmmParams`, normalized transitions, `0.55 log2(NM)` scaling, `FlatSeedIndex`), these unrecorded numbers cannot be verified as current `horizon_wga` behavior.
- **Evidence:**
  - `grep -rn "Cannabis\|418,663\|66\.8" docs/ scripts/ src/` (exit code 1, 0 matches)
  - `src/mdl/mod.rs:57-74` (`AlignConfig::default()`)
- **Correction:** Add reproducible benchmark scripts, dataset citations, and `33bede8` measurements for *P. falciparum* × *P. vivax* and *C. sativa* to `docs/BENCHMARK_RESULTS.md` and `scripts/`, or remove these unverified figures from `manuscript/horizon_paper.md:65`.

---

## F6 -- C6: `30` randomized e-graph trials do not exist in the repository, and `saturate_synteny` cannot recover `1.00` exact `dup(X)` on general tandem duplications or overlapping inversions

- **Claim:** `C6` (`manuscript/horizon_paper.md:99`)
- **Verdict:** `REFUTED`
- **Flaw:**
  1. No script, test, or benchmark in `src/`, `scripts/`, or `docs/BENCHMARK_RESULTS.md` implements or records the `30` randomized tandem-duplication and overlapping-inversion trials.
  2. `saturate_synteny` (`src/egraph.rs:165-273`) only creates `dup(a)` when `a == b` in `seq(a, b)` (`src/egraph.rs:211-217`); it lacks rules for `seq(X, dup(X)) -> dup(X)`, `seq(dup(X), X) -> dup(X)`, `dup(dup(X)) -> dup(X)`, and `inv(dup(X)) <-> dup(inv(X))`, so `>= 3`-copy tandem duplications and inverted duplications cannot collapse to a single `dup(X)`.
  3. `inv(seq(a, b)) -> seq(inv(b), inv(a))` (`src/egraph.rs:191-206`) is unidirectional (no `seq(inv(b), inv(a)) -> inv(seq(a, b))` rule), preventing re-grouping across overlapping inversions.
- **Evidence:**
  - `grep -rn "30 randomized\|tandem-duplication\|ancestral block" . --exclude-dir=.git --exclude-dir=target` (0 matches outside `manuscript/` and `reports/rounds/round-4/`)
  - `src/egraph.rs:165-273`
- **Prior-round adjudication:** Round 1 (`adjudication-L1.md:18`, `F3`) already upheld that `src/egraph.rs:165-273` cannot collapse `>2` copies or arbitrary nested/overlapping inversions. Those rewrite rules remain unchanged at `HEAD` (`33bede8`).
- **Correction:** Replace the unbacked `30`-trial `1.00` claim in `manuscript/horizon_paper.md:99` with the specific 2-copy tandem-duplication / double-inversion cases verified by `horizon_resolve` (`src/bin/horizon_resolve.rs` and `docs/BENCHMARK_RESULTS.md:67`), or implement `n`-ary `dup` / bidirectional `inv` rules and commit the 30-trial benchmark.

---

## F7 -- Framing / Claim-list steering: `C1`–`C6` omit `horizon_merge` benchmark absence and missing external baselines

- **Flaw:** `C1`–`C6` omit major claims in `manuscript/horizon_paper.md:10, 63, 75, 81-95` that conflict with `docs/BENCHMARK_RESULTS.md`:
  - `docs/BENCHMARK_RESULTS.md:9` explicitly states `horizon_merge` has *"no multi-species benchmark row yet"*, whereas `manuscript/horizon_paper.md:10, 63, 75, 81-95` reports `horizon_merge` runtimes across 7 multi-species suites (`406` mammals in `3.84 s`, `317` primates in `2.53 s`, `60_amniotes` in `16.61 s`, etc.).
  - `docs/BENCHMARK_RESULTS.md:5` notes *"head-to-head `minimap2`/`odgi`/`vg`/`LASTZ` baselines remain pending external tool availability on this host"*.
- **Prior-round adjudication:** Round 3 (`F7`, `UPHELD`) flagged claim-list steering for omitting missing baselines and stale rows; `horizon_merge` and Table 1 remain omitted from `C1`–`C6` in Round 4.

---

## Verdict Table

| Claim | Verdict | Summary |
| :--- | :--- | :--- |
| `C1` | `NARROWED` | `src/bin/` has 10 binaries (not 9), `horizon` shares no core, and 3 binaries do not use `bocpd`/`mdl`/`egraph` |
| `C2` | `NARROWED` | `451.35 s` (`2.35×`) measured on `53130b8` (6 commits behind `HEAD` `33bede8`, 21 `src/` files changed) with different LCB/DP counts vs `1e47fe0` |
| `C3` | `NARROWED` | `chr6` `horizon_gfa extract` (`244.93 s`) pinned to `b857377` (11 commits behind `HEAD`); `hprc-v1.0-mc-grch38.gfa.gz` absent on disk |
| `C4` | `NARROWED` | `horizon_clean` severs `walk_support <= 1` (not only `0`), is a no-op on `*`-seq GFAs, and ignores `L`-line overlaps when splitting `W` coordinates |
| `C5` | `REFUTED` | *Plasmodium* (`87.7%`/`54.4%`/`66.8%`) and *Cannabis* (`418,663`/`99.95%`) benchmarks are absent from `docs/BENCHMARK_RESULTS.md`, `scripts/`, and `data/` |
| `C6` | `REFUTED` | `30` randomized trials do not exist in repo, and `src/egraph.rs:165-273` lacks `>2`-copy `dup` and bidirectional `inv` rules needed for `1.00` recovery |
