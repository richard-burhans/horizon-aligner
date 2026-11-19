# Round 4 Adjudication -- L1 (Correctness)

**Target:** `reports/rounds/round-4/attack-L1-correctness.md` (digest `ca60398f47b646ac94ff02bad7dac92d4cb6ccf280a9c7670d179ff6f09d6dd7`, verified against `HEAD` `33bede8`).
**Method cards (`refs/methods/`):** None exist in this repository (`refs/methods/` is absent); rulings rest on direct inspection of `src/`, `docs/BENCHMARK_RESULTS.md`, `git` history, and `reports/rounds/round-4/ATTACKED-horizon_paper.md`.

---

## F1 -- C1: Binary count is 10 (not 9), `horizon` shares no core, and 3 tools use no information-theoretic core

- **Verdict:** `NARROWED`
- **What survives:** `src/bin/` contains 10 executables (`horizon.rs` is a `std::process::Command` dispatcher over 9 binaries and imports nothing from `horizon_aligner`, contradicting `manuscript/horizon_paper.md:18`'s `"nine"` and conflicting with `line 101`'s `"10 CLI binaries"`), and `horizon_gfa` / `horizon_lift` do not use BOCPD/MDL/CRC/e-graphs—but the attacker is **wrong** about `horizon_merge`.
- **Diagnosis ruling:**
  - **UPHELD** on binary count and dispatcher: `src/bin/` has 10 `.rs` files; `src/bin/horizon.rs:7-53` dispatches 9 subcommands via `std::process::Command` with zero `horizon_aligner` imports; `src/bin/horizon_gfa.rs:4` imports only `gfa::*`; `src/bin/horizon_lift.rs` imports only `coord::Interval` and `gfa::*`.
  - **OVERTURNED** on `horizon_merge`: `src/bin/horizon_merge.rs` calls `horizon_aligner::merge_pairwise_mafs` (`src/merge.rs:266`), which explicitly imports `use crate::{egraph, resolve_lcb_duplications};` (`src/merge.rs:3`), invokes e-graph duplication resolution (`src/merge.rs:218-219`), and emits `# EG:Z:` paralog tags (`src/merge.rs:372`). It does **not** import `src/gfa` or `src/coord`.
- **Correction ruling & cost:** **OVERTURNED on the attacker's proposed replacement text**, which would install a new falsehood (*"GFA/liftover/merge modules sharing the streaming `gfa`/`coord` core"*) because `horizon_merge` uses `egraph` and does not use `src/gfa` or `src/coord`. Updating `manuscript/horizon_paper.md:18` to *"ten executables (a `horizon` dispatcher over nine specialized tools)"* and scoping the core-sharing claim to exclude `horizon`, `horizon_gfa`, and `horizon_lift` is **free** (character edit, no re-measurement).
- **Prior rounds:** Extends Round 1 (`C1`, `NARROWED`) and Round 3 (`F3`, `NARROWED`) after `f777935` added `src/bin/horizon.rs`.

---

## F2 -- C2: Full `chr6` `451.35 s` (`2.35×` speedup) is pinned to `53130b8` (6 commits behind `HEAD`) and compares non-identical workloads

- **Verdict:** `UPHELD`
- **What survives:** Every factual premise is verified—`451.35 s` (`10,515.0 MB` RSS) was measured on `53130b8` (`docs/BENCHMARK_RESULTS.md:27, 65`), 6 commits and 21 `src/` files (`951` insertions, `523` deletions, including `src/seeding/flat_index.rs`, `src/seeding/mod.rs`, `src/seeding/translated.rs`, `src/mdl/dp.rs`, and `src/mdl/mod.rs`) behind `HEAD` (`33bede8`), and compares non-identical LCB/DP counts against `1e47fe0` (`292,534` fwd / `291,150` rev LCBs and `693,242` DP jobs vs. `301,156` fwd / `298,178` rev LCBs and `778,471` DP jobs).
- **Diagnosis ruling:** **UPHELD.** Note for completeness that `manuscript/horizon_paper.md:105` (Figure 5 caption) already cites `1e47fe0`, `53130b8`, and `seed_end` (`184.65 s` → `91.04 s`, `2.03×`), whereas `manuscript/horizon_paper.md:101` (body text) omits the `53130b8` commit pin and workload difference.
- **Correction ruling & cost:** Pinning `451.35 s` / `10,515.0 MB` in `manuscript/horizon_paper.md:101` to commit `53130b8` (`HORIZON_PROFILE=1`) and noting the non-identical LCB/DP counts is **free**; re-running full `GRCh38.chr6` × `CHM13.chr6` on `HEAD` (`33bede8`) is **not free**.

---

## F3 -- C3: HPRC `chr6` `horizon_gfa extract` (`244.93 s`) is pinned to stale commit `b857377`, and `hprc-v1.0-mc-grch38.gfa.gz` is absent on disk

- **Verdict:** `NARROWED`
- **What survives:** `docs/BENCHMARK_RESULTS.md:57` pins `chr6` `horizon_gfa extract` (`244.93 s`, `1,255.1 MB` RSS) to `b857377` while `manuscript/horizon_paper.md:101` omits the pin, and `hprc-v1.0-mc-grch38.gfa.gz` is absent locally—but the attacker miscounted both the commit distance and the number of touching commits.
- **Diagnosis ruling:**
  - **UPHELD** on the core defect: `docs/BENCHMARK_RESULTS.md:57` tags `chr6` `extract` with `(`b857377`)`, `manuscript/horizon_paper.md:101` presents it unpinned alongside `HEAD` results, and `data/` does not exist on this host so none of the `12 GB` HPRC timings can be re-verified here (`cmd_stats` in `src/bin/horizon_gfa.rs:39-77` is verified `O(1)`-memory).
  - **NARROWED** on the git counts: `git rev-list --count b857377..HEAD` is **24 commits** behind `HEAD` (not 11), and **7 commits** (`7bf688c`, `0e70693`, `488b020`, `7855881`, `8c7d929`, `b068896`, `33bede8`) modified `src/bin/horizon_gfa.rs` and `src/gfa/mod.rs` (not "Six", even though the attack itself listed 7 hashes).
- **Correction ruling & cost:** Pinning `244.93 s` / `1,255.1 MB` in `manuscript/horizon_paper.md:101` to `b857377` is **free**; re-running on `33bede8` where the HPRC archive is available is **not free**.
- **Prior rounds:** Updates Round 3 (`C6`, `UNCHALLENGED` / completeness Gap 3).

---

## F4 -- C4: `horizon_clean` severs `walk_support <= 1` (not only `walk_support == 0`), is a no-op on `*`-sequence GFAs, and sums `W` lengths without overlap adjustment

- **Verdict:** `UPHELD`
- **What survives:** All three code observations are verified at `src/bin/horizon_clean.rs:154-176` (`has_seq = false` on `*`-sequence segments), `line 260` (`(cfg.prune_unsupported && sup == 0) || (sup <= 1 && cp >= cfg.min_cp_bits)`), and `lines 418-420` (`cur_end += step_len` without subtracting `L`-line overlap CIGARs), plus `docs/BENCHMARK_RESULTS.md:66` (`0 severed (no-op on *-seq)`).
- **Diagnosis ruling:** **UPHELD.** `manuscript/horizon_paper.md:49` states `severing unsupported (walk_support == 0) discordant edges`, which omits that `walk_support == 1` edges are also severed when `cp >= cfg.min_cp_bits` and that `*`-sequence GFAs yield `0.0` changepoint bits.
- **Correction ruling & cost:** Updating `manuscript/horizon_paper.md:49` to state `walk_support <= 1` and noting that sequence-less `*` segments yield `0.0` changepoint bits unless `--prune-unsupported` is enabled for `walk_support == 0` is **free** (character edit, no re-measurement).
- **Prior rounds:** Extends Round 2 (`R2-RESOLVE-LIFT`) and Round 3 (`adjudication-L5.md`, `adjudication-completeness.md`).

---

## F5 -- C5: *Plasmodium* (`87.7%`, `54.4%`, `66.8%`) and *Cannabis sativa* (`418,663` / `418,465`, `99.95%`, `1.000`) benchmarks are absent from `docs/BENCHMARK_RESULTS.md`, `scripts/`, and `data/`

- **Verdict:** `UPHELD`
- **What survives:** Neither benchmark exists anywhere in `docs/BENCHMARK_RESULTS.md`, `scripts/`, `src/`, or `data/` (which does not exist), and `AlignConfig::default()` (`src/mdl/mod.rs:57-75`) already sets `band_w: 180` and `w = 10` fallback defaults (`fallback_min_anchors: 4`, `fallback_min_span: 50`, `fallback_max_occ: 2`), leaving the `54.4%` vs `66.8%` contrast and all *Plasmodium*/*Cannabis* figures ungrounded in the repository.
- **Diagnosis ruling:** **UPHELD.** Note that the same unverified figures also appear in `manuscript/horizon_paper.md:10` (Abstract) and `line 69` (Figure 3 caption), not only `line 65`.
- **Correction ruling & cost:** Removing these unverified figures from `manuscript/horizon_paper.md:10, 65, 69` is **free**; adding reproducible benchmark scripts, dataset citations, and `33bede8` measurements is **not free**.

---

## F6 -- C6: `30` randomized e-graph trials do not exist in the repository, and `saturate_synteny` cannot recover `1.00` exact `dup(X)` on general tandem duplications or overlapping inversions

- **Verdict:** `UPHELD`
- **What survives:** No `30` randomized trials exist in `src/`, `scripts/`, or `docs/BENCHMARK_RESULTS.md`, and `saturate_synteny` (`src/egraph.rs:165-273`) still only creates `dup(a)` for `seq(a, b)` with `a == b` (`lines 211-217`) and only implements unidirectional `inv(seq(a, b)) -> seq(inv(b), inv(a))` (`lines 191-206`, whereas `manuscript/horizon_paper.md:44` writes `<->`).
- **Diagnosis ruling:** **UPHELD.** Round 1 (`adjudication-L1.md:18`, `F3`) already upheld this exact limitation in `src/egraph.rs:165-273`, and the rewrite rules remain unchanged at `HEAD` (`33bede8`).
- **Correction ruling & cost:** Replacing the unbacked `30`-trial `1.00` claim in `manuscript/horizon_paper.md:99` with the verified synthetic `horizon_resolve` benchmark (`docs/BENCHMARK_RESULTS.md:67`) is **free**; implementing `n`-ary `dup` / bidirectional `inv` rules and committing a 30-trial benchmark is **not free**.

---

## F7 -- Framing / Claim-list steering: `C1`–`C6` omit `horizon_merge` benchmark absence and missing external baselines

- **Verdict:** `UPHELD`
- **What survives:** `docs/BENCHMARK_RESULTS.md:9` explicitly states `horizon_merge` has *"no multi-species benchmark row yet"* and `line 5` notes *"head-to-head `minimap2`/`odgi`/`vg`/`LASTZ` baselines remain pending external tool availability on this host"*, yet `manuscript/horizon_paper.md:10, 63, 75, 81-95` reports `horizon_merge` runtimes across 7 multi-species suites that are omitted from `C1`–`C6`.
- **Diagnosis ruling:** **UPHELD.** Continues Round 3 (`F7`, `UPHELD`).
- **Correction ruling & cost:** Including `horizon_merge` / Table 1 claims in the next round's claim list or aligning `manuscript/horizon_paper.md` and `docs/BENCHMARK_RESULTS.md` status notes is **free**; running external baselines is **not free**.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | `src/bin/` has 10 binaries (`horizon.rs` shares no core) and `horizon_gfa`/`horizon_lift` use no BOCPD/MDL/CRC/e-graphs, but `horizon_merge` DOES invoke `egraph`/`resolve_lcb_duplications` (`src/merge.rs:3, 218-219, 372`) and the proposed replacement text falsely claims `horizon_merge` shares `gfa`/`coord`; changing `"nine"` to `"ten"` in `manuscript/horizon_paper.md:18` is free |
| F2 | UPHELD | `451.35 s` (`10,515.0 MB`) is pinned to `53130b8` (6 commits / 21 `src/` files behind `HEAD`) with non-identical LCB/DP counts vs `1e47fe0`; pinning `53130b8` in `manuscript/horizon_paper.md:101` (matching `line 105`) is free |
| F3 | NARROWED | `chr6` `horizon_gfa extract` (`244.93 s`) is pinned to `b857377` in `docs/BENCHMARK_RESULTS.md:57` and `hprc-v1.0-mc-grch38.gfa.gz` is absent locally, but `b857377` is 24 commits behind `HEAD` (not 11) with 7 touching commits (not 6); pinning `b857377` in `manuscript/horizon_paper.md:101` is free |
| F4 | UPHELD | `src/bin/horizon_clean.rs:260` severs `sup <= 1`, `lines 154-176` yield `0.0` bits on `*`-seq GFAs, and `lines 418-420` ignore `L`-line overlaps; updating `manuscript/horizon_paper.md:49` is free |
| F5 | UPHELD | *Plasmodium* (`87.7%`/`54.4%`/`66.8%`) and *Cannabis* (`418,663`/`418,465`/`99.95%`) benchmarks are absent from `docs/BENCHMARK_RESULTS.md`, `scripts/`, `src/`, and `data/`, and `AlignConfig::default()` already uses `band_w: 180` and `w = 10` fallback; removing unverified figures from `manuscript/horizon_paper.md:10, 65, 69` is free |
| F6 | UPHELD | `30` randomized trials do not exist in the repo and `src/egraph.rs:165-273` lacks `>2`-copy `dup` and bidirectional `inv` rules (re-affirming Round 1 `F3`); replacing the claim in `manuscript/horizon_paper.md:99` with `docs/BENCHMARK_RESULTS.md:67` is free |
| F7 | UPHELD | `C1`–`C6` omit `horizon_merge` multi-species runtimes reported in `manuscript/horizon_paper.md:10, 63, 75, 81-95` despite `docs/BENCHMARK_RESULTS.md:9` stating no multi-species benchmark row exists yet (continuing Round 3 `F7`) |
