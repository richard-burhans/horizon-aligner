# Round 6 Adjudication of Completeness Critic

- **Target:** `reports/rounds/round-6/completeness-critic.md` (Round 6, `42700e142b70996c23cdf391ac29e1e5bcf4640d`)
- **Method cards:** None (`refs/methods/` does not exist per `docs/round_rubric.md:36`).
- **Verification performed:** Checked `round.json`, `CARRIED.md`, `CRITIC-BRIEF.md`, all 7 `attack-L*.md` filings, `docs/manual.md`, `src/bin/horizon.rs`, `src/bin/horizon_clean.rs`, `src/bin/horizon_locus.rs`, `src/mdl/dag_seed.rs`, and `src/mdl/mod.rs`, and re-ran `target/release/horizon_filter --net --query ... --target ...` and `target/release/horizon_locus` on `/tmp/vgp_test/locus_koala_amphioxus.fa` at `42700e1`.

---

## Accounting Check (Section 1)

The critic's Section 1 accounting is verified: all 7 declared lenses (`L1`–`L7`) ran and filed (`aborted: []`, `unbriefed: []`, `lens_drops: {}` in `round.json`), zero lenses are `UNRUN`, `CRITIC-BRIEF.md:47-53` recorded no adjudications filed at the time of the critic's run, and the per-lens finding counts (`L1`: 5, `L2`: 5, `L3`: 6, `L4`: 5, `L5`: 6, `L6`: 6, `L7`: 5) match the 7 `attack-L*.md` files.

---

## Rulings on Claim-List & Framing Critique (Section 2)

### `2.1` — `C5` is over-bundled and references an unsnapshotted document (`docs/manual.md`)

- **Critic's filing:** "`round.json:68` bundles 6 full documents (`ARCHITECTURE.md`, `README.md`, `docs/manual.md`, `docs/VGP_INTEGRATION_PLAN.md`, `manuscript/horizon_paper.md`, `manuscript/horizon_companion.md`) into a single correctness claim (`C5`), while `round.json:84-134` (`attacked`) snapshots only 7 files and **omits `docs/manual.md`**." "Bundling 6 documents into one claim encouraged attackers (`L1`, `L2`, `L7`) to spot-check `horizon_resolve` (`ARCHITECTURE.md:28`) while leaving `docs/manual.md` and `horizon_clean` unverified."
- **Verification:** Confirmed in `round.json:68, 84-134`. `L1 F5`, `L2 F5`, and `L7 F5` each examined only `horizon_resolve` (`ARCHITECTURE.md:28`) for `C5`, leaving `docs/manual.md` and `horizon_clean` unverified.
- **Verdict:** **UPHELD**. Both the 6-document bundling in `C5` / omission of `docs/manual.md` from `round.json`'s `attacked` list and the resulting spot-check narrowing by `L1`/`L2`/`L7` are verified. (Unbundling `C5` and snapshotting `docs/manual.md` in future rounds is a **free** harness-input repair.)

### `2.2` — Silent omission of carried Round 5 `.hal` and `horizon_merge`/`horizon_guide` obligations

- **Critic's filing:** "However, **no lens** flagged that `C1`–`C5` silently omit Round 5's carried `R2-HMM-XDROP` item (*multi-contig progressive `.hal` emission in `horizon_guide` and un-ignoring `real_hal_decode_from_env`*) and `R2-RESOLVE-LIFT` item (*`horizon_merge`/`horizon_guide` MAF accuracy verification*), even though `ATTACKED-BENCHMARK_RESULTS.md:10-11, 57-58` highlights both tools."
- **Verification:** Confirmed in `round.json:265-274` (Round 5 carried `R2-HMM-XDROP` and `R2-RESOLVE-LIFT`), `ATTACKED-BENCHMARK_RESULTS.md:10-11, 57-58` (highlighting `horizon_merge` on 92 species and `horizon_guide --progressive --hal --rbest` on 5 mammals), and all 7 `attack-L*.md` files (`L1 F1` and `L2 F1` mention the 1 ignored `.hal` test only as a test-count fact; none of `L1`–`L7`, including `L3 F6`, `L5 F6`, and `L6 F6`, flagged the silent omission of Round 5's carried `.hal` and `horizon_merge`/`horizon_guide` MAF-accuracy obligations).
- **Verdict:** **UPHELD**. All 7 lenses missed the silent drop of Round 5's carried `.hal` and `horizon_merge`/`horizon_guide` MAF-accuracy obligations from `C1`–`C5`.

---

## Rulings on Unattacked Material Gaps (Section 3)

### `3.1` — Sequence-aware `horizon_filter --net --query --target` (`848` records, `4.73s`) was never re-run or audited, and fails on `HEAD` (`553` records)

- **Critic's filing:** "No lens executed `--net --query --target`." "Running `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --net --query /tmp/hprc_grch38_chr6.fa --target /tmp/hprc_chm13_chr6.fa --out /dev/null` on `42700e1` emits **`kept=553`** (`filtered_mdl=436`, `elapsed=4.89s`), **not `848`**." Proposed closure: "Track a `horizon_filter` benchmark script in `scripts/benchmarks/`, update `docs/BENCHMARK_RESULTS.md:86` to `553` records, and re-score its HPRC 1-to-1 coverage and exact pair precision."
- **Verification:** Confirmed across `L1`–`L7` that no lens executed `--net --query --target` (`L1 F2`, `L2 F2`, `L7 F2` ran sequence-free `--net` only; `L4 F2` only `stat`-ed `/tmp/chr6_fastga_mdl_net.paf`; `L5 F2` only cited line 86 in passing). Re-running the exact command on `42700e1` returned:
  `[horizon_filter] total=56857 kept=553 filtered_mdl=436 filtered_len=0 elapsed=4.86s`.
- **Verdict:** **UPHELD**. Both the coverage gap and the stale `848` vs `553` record count in `docs/BENCHMARK_RESULTS.md:86` are verified; the proposed repair is valid and **non-free** (requires tracking a benchmark script and re-scoring HPRC 1-to-1 coverage and exact pair precision on the `553`-record output).

### `3.2` — GenomeArk Koala × Amphioxus `horizon_filter` (`C2`) and 25 kb `horizon_locus` (`C3`) quantitative claims were never reproduced or statistically audited

- **Critic's filing:** "No lens (`L1`–`L3`, `L5`, `L7`) re-ran Koala × Amphioxus `horizon_filter` or `horizon_locus` on `HEAD`, audited the definition/denominator of `'71.6% exact overlap with VGP FastGA_11'`, verified the `3,548 bp` simple tandem repeat rejection (`md:f:-210.51 hb:i:1`), or noted that `horizon_locus` on `/tmp/vgp_test/locus_koala_amphioxus.fa` aligns only `9,697 / 25,000` bp (`38.8%`) on `42700e1` (`0.27s`, `51.5 MB` RSS)."
- **Verification:** Confirmed across all 7 `attack-L*.md` files (`L4 F2`/`F3` checked only `/tmp/vgp_test/` timestamps and missing scripts; `L6 F2` noted only the missing `SweepGA` baseline). Re-running `./target/release/horizon_locus --fasta /tmp/vgp_test/locus_koala_amphioxus.fa --out-prefix /tmp/adj_test_koala --threads 4` on `42700e1` yielded `aligned_bp=9697`, `len_a=25000`, `len_b=25000` (`38.8%` aligned, `p_dist=0.222027`, `jc69_dist=0.263272`, `elapsed=0.26s`, `peak_rss_mb=40.1`).
- **Verdict:** **UPHELD**. The GenomeArk Koala × Amphioxus quantitative claims (`C2` and `C3`) went un-reproduced and un-audited across all 7 lenses; the proposed repair is valid and **non-free** (requires tracking inputs/scripts and auditing/re-verifying the `71.6%` overlap definition and `horizon_filter`/`horizon_locus` outputs on `HEAD`).

### `3.3` — `docs/manual.md` was unsnapshotted and un-audited against `C5`, and `src/bin/horizon.rs` dispatcher omits `filter` and `locus`

- **Critic's filing:** "No lens noted that `docs/manual.md` was omitted from `round.json`'s `attacked` list, that `docs/manual.md:151, 158` still claim BOCPD in `horizon_filter` and `horizon_locus`, or that `src/bin/horizon.rs:7-53` (the `horizon` dispatcher explicitly cited in `C5`) omits `filter` and `locus`." Proposed closure: "Snapshot `docs/manual.md` in `round.json` when claimed, remove the false BOCPD claims in `docs/manual.md:151, 158`, and either wire `filter`/`locus` into `src/bin/horizon.rs` or document them as standalone-only binaries."
- **Verification:**
  - `docs/manual.md` is indeed omitted from `round.json:84-134`, and `docs/manual.md:151` falsely claims BOCPD compositional tracking in `horizon_filter` (which `L1`/`L3`/`L5` caught in other docs, missing `docs/manual.md:151`).
  - However, `L4 F5` (`attack-L4-provenance.md:56`) *did* check `src/bin/horizon.rs:7-53` against `ARCHITECTURE.md:29` and `docs/manual.md:150-162`, confirming `ARCHITECTURE.md:29` (and `docs/manual.md:54`) accurately list the 9 subcommands `horizon.rs` actually dispatches; thus `src/bin/horizon.rs` omitting `filter`/`locus` is a CLI feature gap rather than a `C5` doc-to-code mismatch.
  - Additionally, `horizon_locus.rs` delegates pairwise alignment to `src/mdl/mod.rs` (which uses `crate::bocpd`), and the critic did not verify whether `horizon_locus`'s alignment pipeline calls `bocpd_gc_profile` before calling `docs/manual.md:158` false.
- **Verdict:** **NARROWED**. The unsnapshotted `docs/manual.md`, the missed BOCPD claim at `docs/manual.md:151` (`horizon_filter`), and the absence of `filter`/`locus` in `src/bin/horizon.rs:7-53` survive, but `L4 F5` did audit `horizon.rs` against `ARCHITECTURE.md:29`, and `docs/manual.md:158` (`horizon_locus`) was not verified against `src/mdl/mod.rs`. Removing `BOCPD` from `docs/manual.md:151` and documenting `horizon_filter`/`horizon_locus` as standalone binaries (if not wired into `horizon.rs`) are **free** text repairs.

### `3.4` — `P1` (`horizon_clean` `L`-overlap CIGAR `W`-split & `walk_support <= 1`) and `P2` (`horizon_wga --graph` `src/mdl/dag_seed.rs` `aa5_keys` fallback) went unverified

- **Critic's filing:** "No lens (`L1`, `L2`, `L4`, `L7`) inspected `src/mdl/dag_seed.rs` (`aa5_keys`), reproduced `horizon_wga --graph`, or tested `horizon_clean`'s `L`-overlap CIGAR `W`-split arithmetic (`src/bin/horizon_clean.rs:444`) or `walk_support <= 1` default." Proposed closure: "Audit and test `src/mdl/dag_seed.rs` (`aa5_keys`) and `horizon_clean` `L`-overlap `W`-split + `walk_support <= 1` severing on a positive-control non-`0M` GFA fixture."
- **Verification:** Confirmed across all 7 `attack-L*.md` files: `L1 F5` cited `src/bin/horizon_clean.rs:444` only in passing to contrast it with `horizon_resolve`, `L5 F1` audited `align_seq_to_dag_mdl` (`src/mdl/dag.rs:260-265`) rather than `src/mdl/dag_seed.rs`, and `L5 F6` noted only `horizon_clean`'s `0`-severed negative control. Neither `src/mdl/dag_seed.rs` (`aa5_keys` fallback) nor `horizon_clean` on a non-`0M` positive-control GFA was audited or reproduced.
- **Verdict:** **UPHELD**. Both headline `P1`/`P2` capabilities in `C1`/`C5` went unverified across all 7 lenses; the proposed repair is valid and **non-free** (requires code audit + positive-control fixture testing).

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| 2.1 | UPHELD | `round.json:68` bundles 6 docs into `C5` while omitting `docs/manual.md` from `attacked` (`round.json:84-134`), and `L1`/`L2`/`L7` spot-checked only `horizon_resolve`. |
| 2.2 | UPHELD | All 7 lenses missed that `C1`–`C5` silently drop Round 5's carried `.hal` (`R2-HMM-XDROP`) and `horizon_merge`/`horizon_guide` MAF-accuracy (`R2-RESOLVE-LIFT`) obligations. |
| 3.1 | UPHELD | No lens ran `horizon_filter --net --query --target`; on `42700e1` it outputs `kept=553` (`4.86s`), not `848` (`docs/BENCHMARK_RESULTS.md:86`), and re-scoring is non-free. |
| 3.2 | UPHELD | All 7 lenses left Koala × Amphioxus `horizon_filter` (`71.6%` overlap, `md:f:-210.51`) and `horizon_locus` (`9,697/25,000` bp = `38.8%` on `42700e1`) un-reproduced; closure is non-free. |
| 3.3 | NARROWED | `docs/manual.md` was unsnapshotted and `docs/manual.md:151` falsely claims BOCPD in `horizon_filter` (free fix), but `L4 F5` did check `horizon.rs` vs `ARCHITECTURE.md:29` and `docs/manual.md:158` was not verified. |
| 3.4 | UPHELD | No lens audited `src/mdl/dag_seed.rs` (`aa5_keys`), reproduced `horizon_wga --graph`, or tested `horizon_clean` `L`-overlap `W`-split / `walk_support <= 1` on a non-`0M` positive control (non-free). |
