# Round 4 Adjudication — Lens L5 (Inference)

**Target:** `horizon-aligner` (`HEAD` `33bede8`) — `reports/rounds/round-4/attack-L5-inference.md` (seal digest `44da9655c95a` verified: `44da9655c95a637016abb32df28eedd4618435b83ae2c94014b021bedf86359f`).
**Prior-Round & Method-Card Context:** `refs/methods/` does not exist in this repository, so no method card governs these disputes. `scripts/round history` confirms Round 1 (`F1`–`F10`) and Round 3 (`F1`–`F7` in `reports/rounds/round-3/adjudication-L5.md`) previously adjudicated `stream_gfa`, `horizon_eval` oracle `(ref_id, tgt_id)` pairing, `FlatSeedIndex` confounded ablations, `horizon_gfa stats` counter-only memory, and claim-list steering.

---

## F1 — C1: `10` binaries (not `9`), disjoint cores, and non-composable `*`-seq GFA pipe

* **Ruling:** **NARROWED** (upheld on diagnosis; narrowed on the attacker's proposed count phrasing).
* **What survives:** `src/bin/` contains `10` `.rs` files (`bench_chr1.rs`, `horizon.rs`, `horizon_clean.rs`, `horizon_eval.rs`, `horizon_gfa.rs`, `horizon_guide.rs`, `horizon_lift.rs`, `horizon_merge.rs`, `horizon_resolve.rs`, `horizon_wga.rs`), matching `ATTACKED-horizon_paper.md:101` (`"all 10 CLI binaries"`). `horizon_gfa`, `horizon_lift`, `horizon_merge`, and `horizon` (a process-spawning dispatcher for `9` subcommands) do not import `mdl` or `bocpd`, `horizon_clean` (`src/bin/horizon_clean.rs:90-127`) uses an ad-hoc `bernoulli_ll_bits` test rather than `src/mdl/`, and `docs/BENCHMARK_RESULTS.md:66` confirms `horizon_clean --gfa -` is a `0`-severed no-op on `horizon_wga`'s `*`-sequence GFA. However, the attacker's proposed correction miscounts `horizon.rs` as a `horizon_*` tool (`"9 horizon_* tools + bench_chr1"`): there are **`8` `horizon_*` tools + `horizon` dispatcher + `bench_chr1`** (`9` subcommands dispatched by `horizon`).
* **Repair:** **ENDORSED (with count fix)** — narrowing `manuscript/horizon_paper.md:18` to `10` Rust binaries (`8` `horizon_*` tools, the `horizon` dispatcher, and `bench_chr1`), noting which modules share MDL/BOCPD/e-graph primitives and that `horizon_clean` requires populated `S` sequences, is doc-only; **the fix is free**.

## F2 — C2: Full `chr6` `2.35×` speedup is confounded across 4 commits, stale vs. `HEAD` (`33bede8`), and unvalidated for accuracy

* **Ruling:** **UPHELD**
* **What survives:** `git log 1e47fe0..53130b8 --oneline` confirms three intervening bug-fix commits (`8c7d929`, `5d341a7`, `b068896`) before `53130b8`, shifting LCBs (`292,534`/`291,150` to `301,156`/`298,178`), increasing DP jobs by `+12.3%` (`693,242` to `778,471`), and raising peak RSS (`10,022.0 MB` to `10,515.0 MB`; `docs/BENCHMARK_RESULTS.md:64-65`). Six subsequent commits—culminating in `33bede8` (`20` files changed, including `src/mdl/dp.rs`, `src/mdl/mod.rs`, `src/seeding/flat_index.rs`, `src/seeding/mod.rs`, and `src/seeding/translated.rs`)—make the `451.35 s` figure stale vs. `HEAD`, and `docs/BENCHMARK_RESULTS.md:27, 65` includes no accuracy oracle or external aligner baseline. (Note that `ATTACKED-horizon_paper.md:105` Figure 5 caption already names `1e47fe0` and `53130b8`, whereas line `101` omits the commit pins and `+12.3%` DP job shift.)
* **Repair:** **SPLIT** — narrowing `manuscript/horizon_paper.md:101` to attribute `451.35 s` (`10,515.0 MB` RSS) to `53130b8` vs. `1e47fe0` across intervening seeding/DP fixes (`+12.3%` DP jobs) with no full-chromosome accuracy oracle or external baseline is **free**; re-running full `chr6` on `33bede8` or with an isolated `FlatSeedIndex`/Rayon toggle is **non-free**.

## F3 — C3: `horizon_gfa stats` (`97.0 MB`) is counter-only record discarding; `chr6` extract is pinned to `b857377`

* **Ruling:** **UPHELD** (repeats Round 3 `adjudication-L5.md` `F6`, `UPHELD`).
* **What survives:** `src/bin/horizon_gfa.rs:39-77` (`cmd_stats`) increments five `u64` counters (`s_count`, `s_bases`, `l_count`, `p_count`, `w_count`) and immediately drops each parsed line while `stream_gfa` ignores `C`/`J`/`U`/`O` records without graph indexing or link/walk validation; `97.0 MB` RSS reflects `flate2` streaming overhead rather than graph-processing memory. Meanwhile, `horizon_gfa extract` scales linearly with retained segments (`456.4 MB` for `chr22`, `1,255.1 MB` for `chr6`), `docs/BENCHMARK_RESULTS.md:57` pins `chr6` extract to `b857377` rather than `HEAD` (`33bede8`), and no `odgi`/`vg`/`gfatools` baseline is compared. Although `docs/BENCHMARK_RESULTS.md:19` added `"counter-only pass"` after Round 3, `ATTACKED-horizon_paper.md:101, 105` still pairs `stats` (`97.0 MB`) with `extract` under *"keep `12 GB` HPRC memory bounded"*.
* **Repair:** **SPLIT** — qualifying `manuscript/horizon_paper.md:101, 105` (and line `47`) to state that `horizon_gfa stats` (`97.0 MB`) is a counter-only streaming pass without graph indexing/validation and pinning `chr6` extract to `b857377` is **free**; re-running `chr6` extract on `33bede8` or adding external baselines is **non-free**.

## F4 — C4: `horizon_clean` default behavior contradicts `walk_support == 0` severing, tests GC shifts rather than homology, and is a real-data no-op

* **Ruling:** **UPHELD**
* **What survives:** `src/bin/horizon_clean.rs:260, 476` confirms `prune_unsupported` defaults to `false` and `sever = (cfg.prune_unsupported && sup == 0) || (sup <= 1 && cp >= cfg.min_cp_bits)`, so: (a) under defaults, `walk_support == 0` edges are not severed unless `cp >= min_cp_bits`; (b) `--prune-unsupported` severs all `sup == 0` edges unconditionally, bypassing BOCPD/MDL; and (c) `sup == 1` edges are also severed when `cp >= min_cp_bits`. `score_junction_cp_bits` (`src/bin/horizon_clean.rs:100-127`) tests only a Bernoulli GC-count log-likelihood difference gated by `bocpd_jump >= 0.10`, not sequence homology or chimerism, and `docs/BENCHMARK_RESULTS.md:66` confirms the only real-data run severed `0` of `2,776` links (`no-op on *-seq`).
* **Repair:** **ENDORSED** — rewriting `manuscript/horizon_paper.md:49` to match `horizon_clean.rs:100-127, 260, 476` and disclosing that it has only synthetic unit-test validation (`0`-severed no-op on `*`-seq GFA) is doc-only; **the fix is free**.

## F5 — C5: Circular MDL-on-MDL evaluation, `33.2%–45.6%` loss on Horizon's own supported pairs, and ungrounded *Plasmodium*/*Cannabis* provenance

* **Ruling:** **UPHELD**
* **What survives:** `ATTACKED-horizon_paper.md:65` partitions `LASTZ` `HoxD55` columns on *P. falciparum* × *P. vivax* Chr 1 using Horizon's own `10`-bit MDL+BOCPD threshold (`24,748 / 119,582`) and scores `horizon_wga` on abstaining on the MDL-rejected subset (`87.7%`) without an external homology oracle. On Horizon's own MDL-supported pairs, `horizon_wga` recovered only `54.4%` (`66.8%` with `w = 10`, `band_w = 180`), missing `33.2%–45.6%`. The claim leaves ambiguous which parameterization was applied unchanged to `45.5%`-GC intraspecific *Cannabis* (`cs10` × `Cherry_Pie`), and `418,663`/`Cherry_Pie`/`24,748`/`119,582` appear only in `manuscript/` (with `manuscript/scripts/generate_figures.py:218-219` even hardcoding different `n=88,597` and `n=20,968` counts than `horizon_paper.md:65`).
* **Repair:** **SPLIT** — removing the circular *"genuine signal loss"* inference, disclosing the `33.2%–45.6%` sensitivity gap on MDL-supported pairs, and specifying the *Cannabis* parameter set in `manuscript/horizon_paper.md:65` is **free**; adding external homology oracles and reproducible benchmark scripts/logs to `scripts/` and `docs/BENCHMARK_RESULTS.md` is **non-free**.

## F6 — C6: Unbacked `30` randomized trials and tautological `1.00` vs. `0.00` `dup(X)` comparison

* **Ruling:** **UPHELD**
* **What survives:** Searching `src/`, `tests/`, `scripts/`, and `docs/` finds zero code, test, or data artifacts for `30` randomized trials, ancestral block recovery, or greedy 1-to-1 synteny comparison (only unit tests in `src/egraph.rs:483, 494` and the synthetic `200`-segment `horizon_resolve` row in `docs/BENCHMARK_RESULTS.md:67`). Moreover, a greedy 1-to-1 synteny filter structurally discards paralogs and cannot emit `dup(X)` (`0.00` by construction), making `1.00` vs. `0.00` a tautological baseline.
* **Repair:** **SPLIT** — removing the unbacked `30`-trial `1.00` vs. `0.00` claim from `manuscript/horizon_paper.md:99` and restricting it to the unit-tested algebraic rewrites in `src/egraph.rs` and `docs/BENCHMARK_RESULTS.md:67` is **free**; committing and running a `30`-trial simulation against a duplication-aware baseline is **non-free**.

## F7 — Claim-list steering: `C1`–`C6` omit unbenchmarked `horizon_merge`, no-op/synthetic pangenome rows, oracle `horizon_eval`, and missing external baselines

* **Ruling:** **UPHELD** (extends Round 3 `adjudication-L5.md` `F2` and `F7`, both `UPHELD`).
* **What survives:** `ATTACKED-horizon_paper.md:10, 63, 75, 83-95` reports multi-species `horizon_merge` timings (`406` mammals in `3.84 s`, etc.) while `docs/BENCHMARK_RESULTS.md:9` states `"no multi-species benchmark row yet"`. Table 1 (`ATTACKED-horizon_paper.md:81-95`) omits the oracle `(ref_id, tgt_id)` pairing caveat disclosed in `docs/BENCHMARK_RESULTS.md:5`, and `C1`–`C6` omit the `horizon_clean` real-data no-op (`docs/BENCHMARK_RESULTS.md:66`), synthetic-only `horizon_resolve` (`line 67`) and 4-node `horizon_wga --graph` (`line 68`), and absent `minimap2`/`odgi`/`vg`/`LASTZ` baselines (`line 5`).
* **Repair:** **SPLIT** — surfacing these caveats in `manuscript/horizon_paper.md` (and future claim lists) is **free**; adding multi-species `horizon_merge` benchmark rows and external tool baselines is **non-free**.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | `src/bin/` has 10 binaries (`8` `horizon_*` tools + `horizon` dispatcher + `bench_chr1`), 4 skip `mdl`/`bocpd`, and `horizon_clean` is a `*`-seq no-op; attacker's fix miscounted `9 horizon_* tools`, and with that fixed the doc repair is free |
| F2 | UPHELD | `1e47fe0..53130b8` includes 3 bug-fix commits shifting LCBs and DP jobs (`+12.3%`), 6 later commits (`33bede8`) make it stale vs. `HEAD`, and no `chr6` accuracy oracle/baseline exists; doc narrowing is free |
| F3 | UPHELD | Repeats Round 3 `F6`: `horizon_gfa stats` (`97.0 MB`) is counter-only record discarding, `chr6` extract is pinned to `b857377`, and no `odgi`/`vg`/`gfatools` baseline exists; doc qualification is free |
| F4 | UPHELD | `prune_unsupported` defaults to `false` (`src/bin/horizon_clean.rs:260, 476`), `sup == 1` edges are severed, `score_junction_cp_bits` tests GC% shifts rather than homology, and real `*`-seq GFA is a no-op; doc rewrite is free |
| F5 | UPHELD | Circular evaluation on Horizon's own `10`-bit MDL partition, `33.2%–45.6%` loss on MDL-supported pairs, ambiguous *Cannabis* parameter set, and no repo script/log (`generate_figures.py` hardcodes different counts); doc qualification is free |
| F6 | UPHELD | No `30`-trial script/test/data exists in the repo, and `1.00` vs. `0.00` `dup(X)` against a greedy 1-to-1 filter that cannot emit `dup(X)` is tautological; removing/restricting the claim is free |
| F7 | UPHELD | `ATTACKED-horizon_paper.md` reports `horizon_merge` timings absent from `docs/BENCHMARK_RESULTS.md:9`, Table 1 omits the `horizon_eval` oracle-pairing caveat, and `C1`–`C6` omit no-op/synthetic rows and missing baselines |
