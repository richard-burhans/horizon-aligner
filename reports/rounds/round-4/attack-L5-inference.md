# Round 4 Attack — Lens L5 (Inference)

**Target:** `horizon-aligner` (`HEAD` `33bede8`, `manuscript/horizon_paper.md` via `reports/rounds/round-4/ATTACKED-horizon_paper.md`)
**Method Cards & Literature Note:** Per `reports/rounds/round-4/BRIEF-inference.md`, this repository keeps **no method cards** (`refs/methods/` does not exist). `../library` is not present on disk (`NO_LIBRARY`), so no held paper summaries or full texts were consulted; all findings below are grounded in direct inspection of the live codebase, `docs/BENCHMARK_RESULTS.md`, git history, and `scripts/round history`.

---

## F1 -- C1: `10` binaries (not `9`), disjoint cores, and non-composable `*`-seq GFA pipe (`NARROWED`)

- **Claim:** `C1` (`manuscript/horizon_paper.md:18`) — *"Horizon addresses these coupled bottlenecks in a unified, pure-Rust suite of nine composable executables sharing one information-theoretic core."*
- **Prior-Round Context:** `round 3` `adjudication-L5.md` (`F1`, `F3`, both `UPHELD`/`NARROWED`) previously established that `stream_gfa` silently drops `C`/`J`/`U`/`O` records and that `stdin` (`-`) lacks `MultiGzDecoder` and CLI `main()` tests. Since `round 3`, `horizon.rs` was added as a unified dispatcher (`f777935`), bringing `src/bin/` to 10 binaries (`ATTACKED-horizon_paper.md:101` itself states *"all `10` CLI binaries"*).
- **Concrete Flaw:** The inference that all executables share *"one information-theoretic core"* and are *"composable"* overstates both code coupling and pipeline interoperability:
  1. `src/bin/` contains **10** binaries (`bench_chr1.rs`, `horizon.rs`, `horizon_clean.rs`, `horizon_eval.rs`, `horizon_gfa.rs`, `horizon_guide.rs`, `horizon_lift.rs`, `horizon_merge.rs`, `horizon_resolve.rs`, `horizon_wga.rs`), not 9.
  2. Importing `horizon_aligner::*` across `src/bin/` shows that `horizon_gfa` imports only `gfa`, `horizon_lift` imports only `coord` and `gfa`, and `horizon_merge` imports only `merge_pairwise_mafs` (`src/merge.rs`); none of these three use `mdl` or `bocpd`. Meanwhile, `horizon_clean` (`src/bin/horizon_clean.rs:90-127`) implements its own ad-hoc Bernoulli GC log-likelihood (`bernoulli_ll_bits`, `1.5 * n_tot.log2()`) rather than `src/mdl/`.
  3. Composability fails on `horizon_wga`'s own GFA output: `docs/BENCHMARK_RESULTS.md:66` records that piping `horizon_wga`'s `*`-sequence pairwise GFA into `horizon_clean --gfa -` is a complete no-op (`2,776 kept, 0 severed (no-op on *-seq)`).
- **Evidence:**
  - `ls -la src/bin/` (10 `.rs` files) and `ATTACKED-horizon_paper.md:101` (`all 10 CLI binaries`).
  - `grep -rn "horizon_aligner::" src/bin/` showing no `mdl`/`bocpd` imports in `horizon_gfa.rs`, `horizon_lift.rs`, or `horizon_merge.rs`.
  - `docs/BENCHMARK_RESULTS.md:66`: `horizon_clean --gfa -` | `*`-seq pairwise `chr6` peri GFA | `2,776 kept, 0 severed (no-op on *-seq)`.
- **Correction:** Narrow `manuscript/horizon_paper.md:18` to state that Horizon provides `10` Rust binaries (`9` `horizon_*` tools + `bench_chr1`) whose alignment, cleaning, and resolution modules share MDL/BOCPD/e-graph primitives where sequence or walk structure is present, while `horizon_gfa`, `horizon_lift`, and `horizon_merge` are format/coordinate/sweep-line utilities and `horizon_clean` requires populated `S` sequences (not `*`-seq GFA).

---

## F2 -- C2: Full `chr6` `2.35×` speedup is confounded across 4 commits, stale vs. `HEAD` (`33bede8`), and unvalidated for accuracy (`NARROWED`)

- **Claim:** `C2` (`manuscript/horizon_paper.md:101`) — *"full `GRCh38.chr6` (`170.8 Mbp`) × `CHM13.chr6` (`110.7 Mbp`) finished in `451.35 s` (`7m 31s`, `10,515.0 MB` RSS)—a `2.35×` total speedup over our `1,062.45 s` baseline (`12.3×` speedup in `LASTZ12` transition seeding from `588.67 s` to `47.87 s` via 16-bit prefix-bucketed `FlatSeedIndex` and Rayon parallelization; Figure 5)."*
- **Prior-Round Context:** `round 3` `adjudication-L5.md` (`F4`, `UPHELD`) flagged that `488b020` bundled `FlatSeedIndex` with alignment-changing edits without an ablation. Since `round 3`, `FlatSeedIndex` was extracted (`2695375`) and `53130b8` parallelized transition seeding, but the full `chr6` comparison still bundles multiple commits and has not been re-run on `HEAD` (`33bede8`).
- **Concrete Flaw:**
  1. **Confounded ablation:** `git log 1e47fe0..53130b8 --oneline` shows three intervening commits (`8c7d929`, `5d341a7`, `b068896`) fixing critical/required bugs across DP, DAG, seeding, and CLIs before `53130b8`. Consequently, `docs/BENCHMARK_RESULTS.md:64-65` shows the LCB counts changed (`292,534` fwd / `291,150` rev at `1e47fe0` vs. `301,156` fwd / `298,178` rev at `53130b8`), DP jobs increased `+12.3%` (`693,242` -> `778,471`), and peak RSS increased (`10,022.0 MB` -> `10,515.0 MB`). The `2.35×` ratio is not a pure runtime ablation of `FlatSeedIndex` + Rayon on identical alignments.
  2. **Stale relative to `HEAD` (`33bede8`):** Six commits landed after `53130b8`, culminating in `33bede8` (*"Fix all Critical and Required code review findings across suite"*), which modified 20 files (`459 insertions(+), 205 deletions(-)`), including `src/mdl/dp.rs`, `src/mdl/mod.rs`, `src/seeding/flat_index.rs`, `src/seeding/mod.rs`, and `src/seeding/translated.rs`. The `451.35 s` / `10,515.0 MB` numbers are unmeasured on `HEAD`.
  3. **No accuracy or external baseline:** `docs/BENCHMARK_RESULTS.md:27, 65` reports only runtime, RSS, LCBs, and DP jobs on full `chr6`—no precision, recall, truth-oracle coverage, or external aligner (`LASTZ`, `minimap2`) comparison.
- **Evidence:**
  - `git log 1e47fe0..53130b8 --oneline` and `git log 53130b8..HEAD --oneline`.
  - `git show --stat 33bede8` (20 files changed).
  - `docs/BENCHMARK_RESULTS.md:27, 64-65`.
- **Correction:** Re-run full `chr6` on `33bede8` (and against `1e47fe0` or with a toggle for `FlatSeedIndex`/Rayon alone), or narrow `manuscript/horizon_paper.md:101` to state that `451.35 s` (`10,515.0 MB` RSS) was measured at commit `53130b8` relative to `1e47fe0` across intervening seeding/DP fixes (increasing DP jobs from `693,242` to `778,471`), with no full-chromosome accuracy oracle or external aligner baseline.

---

## F3 -- C3: `horizon_gfa stats` (`97.0 MB`) is counter-only record discarding; `chr6` extract is pinned to `b857377` (`NARROWED`)

- **Claim:** `C3` (`manuscript/horizon_paper.md:101`) — *"On the `12 GB` HPRC `v1.0` Minigraph-Cactus pangenome graph (`hprc-v1.0-mc-grch38.gfa.gz`), `horizon_gfa stats` streamed all `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks in `99.97 s` using `97.0 MB` peak RSS, while two-pass `horizon_gfa extract` isolated `GRCh38.chr22` + `CHM13.chr22` in `183.29 s` (`456.4 MB` RSS) and `GRCh38.chr6` + `CHM13.chr6` in `244.93 s` (`1,255.1 MB` RSS)."*
- **Prior-Round Context:** `round 3` `adjudication-L5.md` (`F6`, `UPHELD`) already ruled that `horizon_gfa stats` (~97 MB RSS) only increments five `u64` counters and discards records, so it does not license pangenome graph-processing scalability. `ATTACKED-horizon_paper.md:101, 105` (`Figure 5B`) continues to pair `stats` (`97.0 MB`) with `extract` under the framing *"while streaming GFA tools keep `12 GB` HPRC memory bounded"*.
- **Concrete Flaw:**
  1. `horizon_gfa stats` (`src/bin/horizon_gfa.rs:28-76`) increments five `u64` counters (`segments`, `total_bp`, `links`, `paths`, `walks`) and immediately drops each parsed line, while `stream_gfa` (`src/gfa/mod.rs`) silently ignores `C`/`J`/`U`/`O` records and performs no graph construction, indexing, or link/walk referential-integrity validation. `97.0 MB` RSS reflects `flate2` gzip streaming overhead, not pangenome graph memory scaling.
  2. Actual subgraph extraction (`horizon_gfa extract`) caches retained segments in memory and scales linearly to `456.4 MB` (`chr22`, `921,252` segments) and `1,255.1 MB` (`chr6`, `3,057,602` segments) for just **2** haplotypes out of `24,456` walks—and `docs/BENCHMARK_RESULTS.md:57` still pins `chr6` extract to `b857377`, not `HEAD` (`33bede8`).
  3. No external baseline (`odgi stats`, `vg stats`, `gfatools stat`) is compared.
- **Evidence:**
  - `src/bin/horizon_gfa.rs:28-76` (`cmd_stats`).
  - `docs/BENCHMARK_RESULTS.md:5, 55-57` (noting `b857377` for `chr6` extract and missing external baselines).
- **Correction:** Explicitly qualify `manuscript/horizon_paper.md:101, 105` that `horizon_gfa stats` (`97.0 MB`) is a counter-only streaming pass without graph indexing or validation, note that `GRCh38.chr6` extract (`1,255.1 MB`) was measured on `b857377` (or re-run on `33bede8`), and remove the implication that `97.0 MB` bounds graph-processing memory.

---

## F4 -- C4: `horizon_clean` default behavior contradicts `walk_support == 0` severing, tests GC shifts rather than homology, and is a real-data no-op (`REFUTED`)

- **Claim:** `C4` (`manuscript/horizon_paper.md:49`) — *"3. `horizon_clean` (`src/bin/horizon_clean.rs`) indexes `S`-segment flanks and scores candidate `L`-line junctions with `bocpd_gc_profile` and an MDL two-regime changepoint test (`--min-cp-bits`), severing unsupported (`walk_support == 0`) discordant edges into `HB:i:1` Horizon Bubbles."*
- **Concrete Flaw:**
  1. **Code contradicts the severing rule:** In `src/bin/horizon_clean.rs:260`, the severing predicate is:
     `let sever = (cfg.prune_unsupported && sup == 0) || (sup <= 1 && cp >= cfg.min_cp_bits);`
     and `prune_unsupported` defaults to `false` (`src/bin/horizon_clean.rs:476`). Thus: (a) under default CLI flags, `walk_support == 0` edges are **not** severed unless `cp >= min_cp_bits`; (b) when `--prune-unsupported` is passed, **all** `walk_support == 0` edges are severed unconditionally, bypassing `bocpd_gc_profile` and `--min-cp-bits`; and (c) `sup == 1` edges are also severed whenever `cp >= min_cp_bits`, contradicting the claim that only `walk_support == 0` edges are severed.
  2. **Inferential mismatch (GC discontinuity ≠ discordant edge):** `score_junction_cp_bits` (`src/bin/horizon_clean.rs:100-127`) computes a Bernoulli GC-count log-likelihood difference (`ll_alt - ll_null - 1.5 * n_tot.log2()`) gated by `|gc_prof[split + 15] - gc_prof[split - 1]| >= 0.10`. This tests only whether two segment flanks differ in GC% (e.g., an exon/intron or isochore boundary), not whether an `L`-line is chimeric or non-homologous; conversely, any discordant edge between two segments of similar GC% yields `cp_bits <= 0` and is retained.
  3. **Zero real-data validation:** `docs/BENCHMARK_RESULTS.md:66` shows the only real-data run (`horizon_clean --gfa -` on `chr6` peri GFA) severed `0` of `2,776` links (`no-op on *-seq`), because `S` lines with `*` sequences give `flank_u`/`flank_v` length `0` (`cp_bits = 0.0`).
- **Evidence:**
  - `src/bin/horizon_clean.rs:100-127, 260, 476`.
  - `docs/BENCHMARK_RESULTS.md:66`.
- **Correction:** Rewrite `manuscript/horizon_paper.md:49` to state accurately that `horizon_clean` severs `walk_support <= 1` edges exceeding a GC-composition Bernoulli MDL changepoint (`--min-cp-bits`, gated by a `>= 0.10` `bocpd_gc_profile` jump) and optionally prunes all `walk_support == 0` edges with `--prune-unsupported`, and disclose that it has only synthetic unit-test validation (it is a no-op on `*`-sequence GFA output).

---

## F5 -- C5: Circular MDL-on-MDL evaluation, `33.2%–45.6%` loss on Horizon's own supported pairs, and ungrounded *Plasmodium*/*Cannabis* provenance (`REFUTED`)

- **Claim:** `C5` (`manuscript/horizon_paper.md:65`) — *"End-to-end `horizon_wga` abstained on `87.7%` of those MDL-rejected columns while recovering `54.4%` of MDL-supported pairs (`66.8%` with `w = 10` sparse-exon fallback and `band_w = 180`). Applied with zero parameter changes to intraspecific *Cannabis sativa* (`cs10` × `Cherry_Pie`, `45.5%` GC, `418,663` columns), Horizon retained `99.95%` (`418,465`) of reference `MAF` columns and `1.000` recall on scratch realignments."*
- **Concrete Flaw:**
  1. **Circular ground truth:** In `manuscript/horizon_paper.md:65`, `LASTZ` `HoxD55` columns on *P. falciparum* × *P. vivax* Chr 1 are partitioned into "MDL-rejected" (`24,748 / 119,582`) and "MDL-supported" (`94,834`) using **Horizon's own `10`-bit MDL+BOCPD threshold**, and then `horizon_wga` (using the same threshold) is scored on abstaining on the MDL-rejected subset (`87.7%`). Without an external homology oracle (e.g., syntenic orthologs or simulation), abstaining on columns rejected by Horizon's own scorer does not prove *"genuine signal loss"* rather than false-negative rejection of divergent orthologs.
  2. **Severe sensitivity failure on MDL-supported pairs:** Even on the `94,834` columns that Horizon's own MDL scorer classifies as homologous, end-to-end `horizon_wga` recovered only **`54.4%`** (or **`66.8%`** with `w = 10` fallback and `band_w = 180`), missing **`33.2%–45.6%`** of its own MDL-supported pairs.
  3. **Parameter inconsistency and weak control:** C5 reports two different parameterizations (`54.4%` vs. `66.8%` with `w = 10`, `band_w = 180`) immediately before claiming *"Applied with zero parameter changes"* to intraspecific *Cannabis sativa* (`cs10` × `Cherry_Pie`). Retaining `99.95%` of columns on a near-identical intraspecific `45.5%`-GC pair is a trivial floor check, not evidence of inter-specific AT-rich calibration.
  4. **Ungrounded in repo artifacts:** `grep -rn "418,663\|Cherry_Pie\|24,748\|119,582"` finds these figures only in `manuscript/` (`horizon_paper.md`, `guide_*.md`, `generate_figures.py`), with no corresponding script or table row in `scripts/` or `docs/BENCHMARK_RESULTS.md`.
- **Evidence:**
  - `ATTACKED-horizon_paper.md:65` (`manuscript/horizon_paper.md:65`).
  - `grep -rn "418,663\|Cherry_Pie\|24,748\|119,582" . --exclude-dir=.git --exclude-dir=target --exclude-dir=reports` (only `manuscript/` hits).
- **Correction:** Remove the circular inference that `87.7%` abstention on MDL-rejected `LASTZ` columns proves genuine signal loss without an external truth oracle; explicitly report the `33.2%–45.6%` chaining/sensitivity gap on MDL-supported pairs; specify which parameter set (`54.4%` default vs. `66.8%` `w = 10`) was applied to *Cannabis*; and add the reproducible benchmark scripts/logs to `scripts/` and `docs/BENCHMARK_RESULTS.md`.

---

## F6 -- C6: Unbacked `30` randomized trials and tautological `1.00` vs. `0.00` `dup(X)` comparison (`REFUTED`)

- **Claim:** `C6` (`manuscript/horizon_paper.md:99`) — *"Across `30` randomized tandem-duplication and overlapping-inversion trials, e-graph saturation recovered `1.00` ancestral block sets and `1.00` exact `dup(X)` duplications compared to `0.00` under greedy one-to-one synteny filtering."*
- **Prior-Round Context:** `round 1` `adjudication-L6.md` (`F2`, `UPHELD`) established that `src/egraph.rs` implements only syntactic equality saturation (`Block`, `Seq`, `Inv`, `Dup` with `seq(X, X) -> dup(X)`), not general evolutionary rearrangement models.
- **Concrete Flaw:**
  1. **No benchmark artifact in the repository:** Searching `src/`, `tests/`, `scripts/`, and `docs/` for `30` randomized trials, `dup(X)` trial evaluation, or greedy 1-to-1 synteny comparison yields zero code or data artifacts—only a single deterministic unit test (`egraph_collapses_inverted_tandem_duplication` in `src/egraph.rs:483`).
  2. **Tautological strawman baseline:** By definition, a *"greedy one-to-one synteny filter"* enforces 1-to-1 orthology and discards paralogous copies, so it **cannot** emit `dup(X)` (`0.00` by construction). Reporting `1.00` vs. `0.00` against a baseline that structurally lacks a duplication operator licenses no inference about Horizon's e-graph accuracy relative to duplication-aware rearrangement or pangenome tools.
- **Evidence:**
  - `grep -rn "dup(X)\|ancestral block\|one-to-one synteny\|30.*trial" . --exclude-dir=.git --exclude-dir=target` (zero hits in `src/`, `tests/`, `scripts/`, `docs/`).
  - `src/egraph.rs:275-350, 483` (`MdlCostModel`, single unit test).
- **Correction:** Either commit the `30`-trial simulation script/seed and compare against a duplication-aware baseline, or remove the `30`-trial `1.00` vs. `0.00` claim from `manuscript/horizon_paper.md:99` and restrict the claim to the unit-tested algebraic rewrite rules in `src/egraph.rs` and the synthetic `200`-segment `horizon_resolve` benchmark (`docs/BENCHMARK_RESULTS.md:67`).

---

## F7 -- Claim-list steering: `C1`–`C6` omit unbenchmarked `horizon_merge`, no-op/synthetic pangenome rows, oracle `horizon_eval`, and missing external baselines

- **Context:** Per `BRIEF-inference.md:194-195`, the claim list was authored by the runner. Comparing `C1`–`C6` against `docs/BENCHMARK_RESULTS.md` and `round 3` `adjudication-L5.md` (`F2`, `F7`, both `UPHELD`) reveals systematic steering away from the weakest inferential links in `ATTACKED-horizon_paper.md`:
  1. **`horizon_merge` multi-species claims unbacked in `docs/BENCHMARK_RESULTS.md`:** `ATTACKED-horizon_paper.md:10, 63, 75, 83-95` reports `406` mammals in `3.84 s`, `317` species in `2.53 s`, `60_amniotes` in `16.61 s`, etc., yet `docs/BENCHMARK_RESULTS.md:9` explicitly states: *"Multi-species merging (`horizon_merge`): Implemented in `src/bin/horizon_merge.rs` (best-score interval netting; unit-tested; **no multi-species benchmark row yet**)."*
  2. **Oracle `(ref_id, tgt_id)` pairing in `horizon_eval` (`round 3` `adjudication-L5.md` `F2`):** `docs/BENCHMARK_RESULTS.md:5` now discloses that `horizon_eval` evaluates reference-MAF `(ref_id, tgt_id)` pairs individually rather than all-vs-all across scaffolds, but `ATTACKED-horizon_paper.md:81-95` (Table 1) still reports precision without this oracle-pairing caveat.
  3. **No-op/synthetic pangenome rows & missing external baselines (`round 3` `adjudication-L5.md` `F7`):** `C1`–`C6` omit the `horizon_clean` real-data no-op (`docs/BENCHMARK_RESULTS.md:66`), synthetic-only `horizon_resolve` (`line 67`), 4-node synthetic `horizon_wga --graph` (`line 68`), and the complete absence of `minimap2`/`odgi`/`vg`/`LASTZ` head-to-head baselines (`line 5`).

---

## Verdict Table

| Claim | Verdict | Summary of Inference Defect |
| :--- | :--- | :--- |
| **C1** | **NARROWED** | `src/bin/` has 10 binaries (not 9); `horizon_gfa`/`horizon_lift`/`horizon_merge` do not use `mdl`/`bocpd`; `horizon_clean` uses an ad-hoc GC Bernoulli test; `horizon_wga` `*`-seq GFA piped to `horizon_clean` is a no-op (`docs/BENCHMARK_RESULTS.md:66`). |
| **C2** | **NARROWED** | `2.35×` speedup (`1e47fe0` -> `53130b8`) is confounded by 3 intervening bug-fix commits shifting LCBs and DP jobs (`+12.3%`), is stale vs. `HEAD` (`33bede8`, 20 files changed), and has no full-`chr6` accuracy oracle or external baseline. |
| **C3** | **NARROWED** | `horizon_gfa stats` (`97.0 MB`) is a counter-only record-discarding pass (`src/bin/horizon_gfa.rs:28-76`), not graph-processing scalability; `chr6` extract (`1,255.1 MB`) is pinned to `b857377`; no `odgi`/`vg`/`gfatools` baseline. |
| **C4** | **REFUTED** | `prune_unsupported` defaults to `false` (`src/bin/horizon_clean.rs:260, 476`), `sup == 1` edges are also severed, `score_junction_cp_bits` tests GC% shifts rather than homology/discordance, and real `*`-seq GFA cleaning is a no-op (`docs/BENCHMARK_RESULTS.md:66`). |
| **C5** | **REFUTED** | Circular evaluation of `horizon_wga` against Horizon's own `10`-bit MDL+BOCPD partition; `33.2%–45.6%` sensitivity loss on MDL-supported pairs (`54.4%`/`66.8%`); two-parameter ambiguity vs. "zero parameter changes"; no repo script/log for *Plasmodium*/*Cannabis* numbers. |
| **C6** | **REFUTED** | No `30`-trial script, test, or data artifact exists in the repo; `1.00` vs. `0.00` `dup(X)` against a greedy 1-to-1 filter that cannot emit `dup(X)` is a tautological strawman. |
