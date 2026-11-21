# Round 6 Completeness Critique

- **Target:** Round 6 (`42700e142b70996c23cdf391ac29e1e5bcf4640d`)
- **Files inspected:** `reports/rounds/round-6/CRITIC-BRIEF.md`, `round.json`, `CARRIED.md`, all 7 `attack-L*.md` files (`L1`–`L7`), `ATTACKED-*.md` snapshots, `docs/manual.md`, and `src/bin/horizon.rs`, plus read-only `target/release/horizon_filter` and `target/release/horizon_locus` checks on `42700e1`.

---

## 1. Lens Execution & Filing Accounting (UNRUN vs. UNCHALLENGED)

- **All 7 declared lenses ran and filed:** `attack-L1-correctness.md`, `attack-L2-sampling.md`, `attack-L3-statistics.md`, `attack-L4-provenance.md`, `attack-L5-inference.md`, `attack-L6-prior-art.md`, `attack-L7-reproduction.md`. In `round.json`, `aborted: []`, `unbriefed: []`, `lens_drops: {}`. **Zero lenses are UNRUN.**
- **Adjudication status:** Per `CRITIC-BRIEF.md:47-53`, **no adjudications (`adjudication-L1.md`..`L7.md`) have been filed yet**.
- **Per-lens filing counts:**
  - **L1 (`correctness`):** 5 findings (`F1`–`F5` covering `C1`–`C5`: `1` `SURVIVED` [`C4`], `3` `NARROWED`, `1` `REFUTED` [`C2`]).
  - **L2 (`sampling`):** 5 findings (`F1`–`F5` covering `C1`–`C5`: `4` `NARROWED`, `1` `REFUTED` [`C2`]).
  - **L3 (`statistics`):** 6 findings (`F1`–`F5` covering `C1`–`C5`, all `NARROWED`, plus `F6` claim-list steering).
  - **L4 (`provenance`):** 5 findings (`F1`–`F5` covering `C1`–`C5`: `3` `NARROWED`, `2` `REFUTED` [`C2`, `C3`]).
  - **L5 (`inference`):** 6 findings (`F1`–`F5` covering `C1`–`C5` plus `F6` claim-list steering: `4` `NARROWED`, `2` `REFUTED` [`C2`, `C5`]).
  - **L6 (`prior-art`):** 6 findings (`F1`–`F5` covering `C1`–`C5` plus `F6` meta-finding: all `NARROWED`).
  - **L7 (`reproduction`):** 5 findings (`F1`–`F5` covering `C1`–`C5`: `1` `SURVIVED` [`C5`], `3` `NARROWED`, `1` `REFUTED` [`C2`]).
- **Thinness assessment:** None of the 7 lenses is empty or perfunctory, though `L1`, `L2`, `L4`, and `L7` filed strictly 1-to-1 per declared claim (`F1`–`F5`) without meta-findings on omitted claims or framing (`L3`/`L5`/`L6` each added `F6`).

---

## 2. Claim-List & Framing Critique

1. **`C5` is over-bundled and references an unsnapshotted document (`docs/manual.md`):**
   - `round.json:68` bundles 6 full documents (`ARCHITECTURE.md`, `README.md`, `docs/manual.md`, `docs/VGP_INTEGRATION_PLAN.md`, `manuscript/horizon_paper.md`, `manuscript/horizon_companion.md`) into a single correctness claim (`C5`), while `round.json:84-134` (`attacked`) snapshots only 7 files and **omits `docs/manual.md`**.
   - Bundling 6 documents into one claim encouraged attackers (`L1`, `L2`, `L7`) to spot-check `horizon_resolve` (`ARCHITECTURE.md:28`) while leaving `docs/manual.md` and `horizon_clean` unverified.
2. **Silent omission of carried Round 5 `.hal` and `horizon_merge`/`horizon_guide` obligations:**
   - `L3 F6`, `L5 F6`, and `L6 F6` caught claim-list steering away from `n_cal = 15` CRC, `KARLIN_MDL_COEFF = 0.55`, the `18%`-GC multi-seed sweep, `horizon_clean` positive controls, and `SweepGA`/`impg`/`odgi`/`gfatools` baselines.
   - However, **no lens** flagged that `C1`–`C5` silently omit Round 5's carried `R2-HMM-XDROP` item (*multi-contig progressive `.hal` emission in `horizon_guide` and un-ignoring `real_hal_decode_from_env`*) and `R2-RESOLVE-LIFT` item (*`horizon_merge`/`horizon_guide` MAF accuracy verification*), even though `ATTACKED-BENCHMARK_RESULTS.md:10-11, 57-58` highlights both tools.

---

## 3. Unattacked Material Gaps

1. **Sequence-aware `horizon_filter --net --query --target` (`848` records, `4.73s`) was never re-run or audited, and fails on `HEAD` (`553` records):**
   - **Artifact:** `ATTACKED-BENCHMARK_RESULTS.md:86` (`848 records (109.76 Mbp), 99.7% of 110.08 Mbp 1-to-1 HPRC coverage retained`, `4.73s`) and `ATTACKED-VGP_INTEGRATION_PLAN.md:72`.
   - **Coverage checked:** `L1 F2`, `L2 F2`, `L4 F2`, and `L7 F2` only re-ran/audited the sequence-free `0.47s` run (`1,004` vs `565` records). `L5 F2` cited line 86 only in passing to contrast its runtime with `0.47s`; no lens executed `--net --query --target`.
   - **Why it matters:** Running `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --net --query /tmp/hprc_grch38_chr6.fa --target /tmp/hprc_chm13_chr6.fa --out /dev/null` on `42700e1` emits **`kept=553`** (`filtered_mdl=436`, `elapsed=4.89s`), **not `848`**. Both `chr6` `horizon_filter` rows in `docs/BENCHMARK_RESULTS.md:85-86` are stale on `HEAD`, and row 86 would survive unchallenged.
   - **What would close it:** Track a `horizon_filter` benchmark script in `scripts/benchmarks/`, update `docs/BENCHMARK_RESULTS.md:86` to `553` records, and re-score its HPRC 1-to-1 coverage and exact pair precision.

2. **GenomeArk Koala × Amphioxus `horizon_filter` (`C2`) and 25 kb `horizon_locus` (`C3`) quantitative claims were never reproduced or statistically audited:**
   - **Artifact:** `round.json` `C2` & `C3`; `ATTACKED-BENCHMARK_RESULTS.md:87-89`; `ATTACKED-VGP_INTEGRATION_PLAN.md:72, 75`.
   - **Coverage checked:** `L4 F2`/`F3` checked `/tmp/vgp_test/` timestamps and missing scripts, and `L6 F2` noted the missing `SweepGA` baseline. No lens (`L1`–`L3`, `L5`, `L7`) re-ran Koala × Amphioxus `horizon_filter` or `horizon_locus` on `HEAD`, audited the definition/denominator of `"71.6% exact overlap with VGP FastGA_11"`, verified the `3,548 bp` simple tandem repeat rejection (`md:f:-210.51 hb:i:1`), or noted that `horizon_locus` on `/tmp/vgp_test/locus_koala_amphioxus.fa` aligns only `9,697 / 25,000` bp (`38.8%`) on `42700e1` (`0.27s`, `51.5 MB` RSS).
   - **Why it matters:** `C2` and `C3` present Koala × Amphioxus as VGP's ~500 Myr proof-of-concept, yet its overlap metric (`71.6%`) and locus alignment yield (`38.8%`) went unverified.
   - **What would close it:** Track Koala × Amphioxus benchmark inputs/scripts in `scripts/benchmarks/`, specify whether `71.6%` is record-level or base-level overlap, and verify `md:f:-210.51 hb:i:1` and `horizon_locus` outputs on `42700e1`.

3. **`docs/manual.md` was unsnapshotted and un-audited against `C5`, and `src/bin/horizon.rs` dispatcher omits `filter` and `locus`:**
   - **Artifact:** `round.json:68, 84-134`; `docs/manual.md:54, 151, 158`; `src/bin/horizon.rs:7-53`.
   - **Coverage checked:** `L3 F5` cited `docs/manual.md:144` and `L4 F5` cited `docs/manual.md:150-162` as present, but no lens noted that `docs/manual.md` was omitted from `round.json`'s `attacked` list, that `docs/manual.md:151, 158` still claim BOCPD in `horizon_filter` and `horizon_locus`, or that `src/bin/horizon.rs:7-53` (the `horizon` dispatcher explicitly cited in `C5`) omits `filter` and `locus`.
   - **Why it matters:** `C5` claims `docs/manual.md` and `horizon` dispatcher subcommands match `42700e1`, yet `docs/manual.md` was unpinned and `horizon filter` / `horizon locus` fail if invoked via the unified dispatcher.
   - **What would close it:** Snapshot `docs/manual.md` in `round.json` when claimed, remove the false BOCPD claims in `docs/manual.md:151, 158`, and either wire `filter`/`locus` into `src/bin/horizon.rs` or document them as standalone-only binaries.

4. **`P1` (`horizon_clean` `L`-overlap CIGAR `W`-split & `walk_support <= 1`) and `P2` (`horizon_wga --graph` `src/mdl/dag_seed.rs` `aa5_keys` fallback) went unverified:**
   - **Artifact:** `round.json` `C1`, `C5`; `ATTACKED-PANGENOME_SUITE.md:12, 18`; `ATTACKED-VGP_INTEGRATION_PLAN.md:73-74`.
   - **Coverage checked:** `L5 F1` audited `align_seq_to_dag_mdl` (`-log2(out_degree)` omission) and `L5 F6` noted `horizon_clean`'s `score_junction_cp_bits` and `0`-severed negative control. No lens (`L1`, `L2`, `L4`, `L7`) inspected `src/mdl/dag_seed.rs` (`aa5_keys`), reproduced `horizon_wga --graph`, or tested `horizon_clean`'s `L`-overlap CIGAR `W`-split arithmetic (`src/bin/horizon_clean.rs:444`) or `walk_support <= 1` default.
   - **Why it matters:** Both are headline `P1`/`P2` deliverables in `C1` and `C5` whose runtime correctness on non-`0M` GFAs and translated DAG seeding remain unverified on `HEAD`.
   - **What would close it:** Audit and test `src/mdl/dag_seed.rs` (`aa5_keys`) and `horizon_clean` `L`-overlap `W`-split + `walk_support <= 1` severing on a positive-control non-`0M` GFA fixture.
