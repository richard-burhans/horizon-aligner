# Round 7 ADJUDICATION -- Completeness Critic

**Critic filing reviewed:** [`reports/rounds/round-7/completeness-critic.md`](completeness-critic.md) (unsealed; key claims quoted below per harness instructions).  
**Target commit:** `5cb83fb6a1346a8193f9e4048055638d0a422b0a` (`HEAD`).  
**Prior rounds & method cards:** Checked `scripts/round history "" --kind adjudication` (Rules 5–6). This repository keeps no method cards (`refs/methods/` does not exist).

---

## Finding-by-Finding Rulings

### G1 -- Refuted `λ̂ = 10.0` / fixed CRC floor persists unchallenged in Figure 1 and `README.md`

- **Critic's claim (quoted):** *"Refuted `λ̂ = 10.0` / fixed CRC floor persists unchallenged in Figure 1 and `README.md` (`ATTACKED-horizon_paper.md:22`, `ATTACKED-README.md:8`). All 7 attackers refuted `C3` at `ATTACKED-horizon_paper.md:40`, but none checked `horizon_paper.md:22` ('gated by Conformal Risk Control (`λ̂ = 10.0` bits)') or `README.md:8` ('`0.55 log2(N*M)` + fixed `10.0`-bit conformal floor', which also writes `+` instead of `max'). ... (Similarly, `ATTACKED-ARCHITECTURE.md` and `ATTACKED-VGP_INTEGRATION_PLAN.md` were not inspected by any attacker.)"*
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; prescribed repair **ENDORSED** (the fix is free). Verification against `ATTACKED-horizon_paper.md:22`, `ATTACKED-README.md:8`, and `src/mdl/mod.rs:151,428` confirms that both unchallenged lines still assert a fixed `10.0`-bit CRC floor (and `README.md:8` writes `+` instead of `max`), while `ATTACKED-ARCHITECTURE.md` and `ATTACKED-VGP_INTEGRATION_PLAN.md` do not contain the `λ̂ = 10.0` defect. Updating `horizon_paper.md:22` and `README.md:8` alongside `horizon_paper.md:40` is a text-only edit requiring no re-measurement (**the fix is free**).
- **What survives:** Both `manuscript/horizon_paper.md:22` and `README.md:8` carry the same refuted `λ̂ = 10.0` / fixed-floor claim as `horizon_paper.md:40` and went unflagged by all 7 attackers.

---

### G2 -- External baseline scripts in `C2` (`bench_external_graph_baselines.py`, `bench_external_wga_baselines.py`) were never executed or audited on `5cb83fb`

- **Critic's claim (quoted):** *"External baseline scripts in `C2` (`scripts/benchmarks/bench_external_graph_baselines.py`, `bench_external_wga_baselines.py`) were never executed or audited on `5cb83fb` (`ATTACKED-horizon_paper.md:84`, `ATTACKED-BENCHMARK_RESULTS.md:87-88`). ... no attacker (not even L7 Reproduction) ran either script on `HEAD` to verify that the numbers reproduce or that `horizon_wga` / `horizon_resolve` / `horizon_clean` produce valid outputs on those inputs."*
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; prescribed repair **ENDORSED** (cost: ~15 s runtime + output validation / replicate logging, not free). Across `attack-L1` through `attack-L7`, no attacker executed either script (`L6` only cited them for missing `KegAlign` and `vg snarls` baselines). Live execution on `5cb83fb` here confirms the runtime/RSS numbers reproduce (`bench_external_graph_baselines.py` at `200k`: `horizon_gfa stats` `0.055 s` / `12.3 MB`, `horizon_clean` `0.353 s` / `130.6 MB`, `horizon_resolve` `0.261 s` / `95.5 MB`, `horizon_wga --graph` `1.582 s` / `501.1 MB`; `bench_external_wga_baselines.py --skip-lastz`: `horizon_wga` `1.594 s` / `842.7 MB`, `--one-to-one` `1.513 s` / `870.4 MB` vs. `FastGA` `2.981 s` and `minimap2` `2.060 s`), while neither script validates output correctness or records `n` replicates.
- **What survives:** Both external baseline scripts in `C2` were unrun and unaudited for output validity and replicate count in Round 7, though their timing/RSS numbers reproduce on `5cb83fb`.

---

### G3 -- Random-GC `median=0` sweep (`bench_random_gc_crc.py`) and `estimate_hmm_params` in `C3` were never run or audited

- **Critic's claim (quoted):** *"Random-GC `median=0` sweep (`scripts/benchmarks/bench_random_gc_crc.py`) and `estimate_hmm_params` in `C3` were never run or audited (`round.json:39`, `ATTACKED-horizon_paper.md:61`, `ATTACKED-BENCHMARK_RESULTS.md:13`). `C3.summary` claims `bench_random_gc_crc.py` (`18%/30%/45%` GC, `10 kb`) yields `median=0` false-positive pairs across 6 seeds on `5cb83fb`, whereas `ATTACKED-horizon_paper.md:61` and `ATTACKED-BENCHMARK_RESULTS.md:13` report `288` pairs (`30%` GC) and `149` pairs (`45%` GC) on single draws."*
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; prescribed repair **ENDORSED** (cost: <1 s runtime for `bench_random_gc_crc.py` + code audit of `estimate_hmm_params`; reconciling the text with the measured sweep is free once reported). Neither `L3` nor `L7` ran `scripts/benchmarks/bench_random_gc_crc.py` or audited `estimate_hmm_params` in `src/mdl/mod.rs`. Live execution of `python3 scripts/benchmarks/bench_random_gc_crc.py` on `5cb83fb` confirms `median=0` holds across all 6 seeds at `18%` (`max=95`), `30%` (`max=216`), and `45%` (`max=307`) GC, leaving the single-draw `288` / `149` figures in `horizon_paper.md:61` and `BENCHMARK_RESULTS.md:13` unreconciled with the 6-seed sweep and `estimate_hmm_params` unaudited.
- **What survives:** `bench_random_gc_crc.py` (`median=0` holds on `5cb83fb`, `max=307` pairs) and `estimate_hmm_params` went unrun/unaudited in Round 7, and the single-draw `30%`/`45%`-GC numbers in the draft need reconciliation.

---

### G4 -- Unverified `C1`/`C5` sub-claims: `horizon` 11-subcommand CLI dispatch, `aa5_keys` test, `horizon_clean` `10M` L-overlap W-split test, and `manuscript/horizon_paper.pdf`

- **Critic's claim (quoted):** *"Unverified `C1`/`C5` sub-claims: `horizon` 11-subcommand CLI dispatch (`src/bin/horizon.rs`), `aa5_keys` test, `horizon_clean` `10M` L-overlap W-split test, and `manuscript/horizon_paper.pdf`. While `cargo test` passed overall, no attacker verified that `horizon <subcommand>` dispatches all 11 subcommands (including `filter` and `locus`) or inspected `aa5_keys` / `10M` L-overlap test coverage, and no attacker checked `manuscript/horizon_paper.pdf` against `horizon_paper.md`."*
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; prescribed repair **ENDORSED** (smoke-testing and citing exact test names is free; rebuilding `manuscript/horizon_paper.pdf` after patching `horizon_paper.md` is a low-cost build step). None of the 7 attackers verified the 11-subcommand dispatcher, cited the exact `aa5_keys` or `10M` unit tests, or inspected `manuscript/horizon_paper.pdf`. Live checks on `5cb83fb` confirm `./target/release/horizon <subcommand> --help` exits `0` for all 11 subcommands (`wga`, `merge`, `guide`, `gfa`, `lift`, `clean`, `resolve`, `filter`, `locus`, `eval`, `bench-chr1`), `seeded_dag_falls_back_to_aa5_translated_seeds_on_synonymous_divergence` (`src/mdl/dag_seed.rs:414`) tests `aa5_keys`, and `walk_split_accounts_for_link_overlap_cigar` (`src/bin/horizon_clean.rs:799`) tests `10M` L-overlap W-splits; `manuscript/horizon_paper.pdf` (`Nov 20 20:10`) still mirrors the unpatched `horizon_paper.md` and must be rebuilt after `C3`/`C4` edits.
- **What survives:** All four `C1`/`C5` sub-items were unverified in Round 7; the CLI and unit tests pass on `5cb83fb`, while `manuscript/horizon_paper.pdf` awaits rebuild after manuscript fixes.

---

### G5 -- Omitted Table 1, HPRC/full-`chr6`/`horizon_merge`/`horizon_locus` benchmarks, and carried `R2-INDEX-RAM` / `estimate_dag_hmm_params` items

- **Critic's claim (quoted):** *"Entire Table 1 (`horizon_paper.md:63-72`), `horizon_gfa` / `horizon_lift` / `horizon_merge` / full `chr6` WGA (`horizon_paper.md:76-78`), VGP Koala × Amphioxus / `horizon_locus` (`BENCHMARK_RESULTS.md:89-93`), and carried `R2-INDEX-RAM` / `estimate_dag_hmm_params` items omitted from `C1`–`C5`. ... Either re-run/pin these benchmarks on `5cb83fb` in a follow-up round or explicitly scope `horizon_paper.md` / `BENCHMARK_RESULTS.md` to their pinned commits (`8cd5db0`, `2896229`, `42700e1`) and carry `R2-INDEX-RAM`."*
- **Verdict:** **NARROWED**
- **Ruling:** Diagnosis **NARROWED**; prescribed repair **NARROWED** (inline commit tagging for Table 1 and `horizon_merge` in `horizon_paper.md` is **free**; re-running full `chr6`/HPRC/Table 1 on `5cb83fb` is an expensive follow-up-round benchmark cost). The coverage gap is real: none of these benchmarks were re-verified on `5cb83fb` in Round 7, and `C2` bundled `R2-INDEX-RAM` without testing it (the tracked `chr6` scripts in `scripts/benchmarks/` still default to `/tmp/`). However, the critic's framing overstates the pinning gap: `BENCHMARK_RESULTS.md:5,26,28,73,81,91-93` and `horizon_paper.md:76,78,84` **already** explicitly pin the HPRC, full-`chr6`, and `horizon_locus` runs to `8cd5db0`, `2896229`, and `42700e1`/`9dfec4a` (and `README.md:110` pins Table 1 to `8cd5db0`), leaving only Table 1 (`horizon_paper.md:61-72`) and `horizon_merge` (`horizon_paper.md:84`) without inline commit tags in `horizon_paper.md`.
- **What survives:** Table 1, HPRC/full-`chr6`/`horizon_merge`/`horizon_locus` benchmarks were not re-run on `5cb83fb` and `R2-INDEX-RAM` / `estimate_dag_hmm_params` remain open, while most of these rows are already pinned to `8cd5db0`/`2896229`/`42700e1` in `BENCHMARK_RESULTS.md` and `horizon_paper.md`.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| G1 | UPHELD | `horizon_paper.md:22` and `README.md:8` still carry the refuted `λ̂ = 10.0` / fixed-floor claim unflagged by all 7 attackers; text fix is free |
| G2 | UPHELD | neither `bench_external_graph_baselines.py` nor `bench_external_wga_baselines.py` was run in Round 7; live runs on `5cb83fb` reproduce timing/RSS, while output validation and `n` replicates remain unrecorded |
| G3 | UPHELD | `bench_random_gc_crc.py` and `estimate_hmm_params` went unrun/unaudited in Round 7; live run confirms `median=0` at 18%/30%/45% GC (`max=307`), needing reconciliation with single-draw text |
| G4 | UPHELD | 11-subcommand `horizon` dispatch, exact `aa5_keys`/`10M` test names, and `horizon_paper.pdf` were unchecked in Round 7; CLI/tests pass on `5cb83fb`, PDF needs rebuild after edits |
| G5 | NARROWED | benchmarks were omitted from Round 7 and `R2-INDEX-RAM` remains open, but `BENCHMARK_RESULTS.md:5` and `horizon_paper.md:76,78,84` already pin most rows to `8cd5db0`/`2896229`/`42700e1` |
