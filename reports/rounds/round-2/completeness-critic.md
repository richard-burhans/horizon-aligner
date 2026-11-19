# Round 2 Completeness Critic — Horizon Pangenome Suite, Adaptive Pair-HMM, E-Graph Resolver, and HPRC Benchmarks

**Target round:** `reports/rounds/round-2/` (`round.json`, `CRITIC-BRIEF.md`, `ATTACKED-BENCHMARK_RESULTS.md`)
**Artifacts read:** `CRITIC-BRIEF.md`, `round.json`, `ATTACKED-BENCHMARK_RESULTS.md`, all 7 attack files (`attack-L1-correctness.md` through `attack-L7-reproduction.md`), all 7 adjudication files (`adjudication-L1.md` through `adjudication-L7.md`), `reports/rounds/round-1/completeness-critic.md`, `docs/BENCHMARK_RESULTS.md`, `docs/STRATEGIC_PLAN.md`, `docs/PANGENOME_SUITE.md`, `src/bin/horizon_clean.rs`, and `src/bin/horizon_merge.rs`.

---

## 1. Execution Status: UNRUN vs. UNCHALLENGED (and Brief Discrepancy)

* **Zero lenses were UNRUN:** All 7 declared lenses (`correctness`, `sampling`, `statistics`, `provenance`, `inference`, `prior-art`, `reproduction`) filed both an attack (`attack-L*.md`) and an adjudication (`adjudication-L*.md`).
* **Metadata discrepancy in `CRITIC-BRIEF.md:49`:** `CRITIC-BRIEF.md:49` states `attack-L3-statistics.md (statistics) — NO ADJUDICATION FILED`, but `reports/rounds/round-2/adjudication-L3.md` (72 lines, SHA-256 prefix `3952ea463b86` verified) exists on disk in the round directory. `L3` was **not** unrun or unadjudicated; `CRITIC-BRIEF.md` reflects stale state prior to `adjudication-L3.md` being registered (the same status-generation lag seen for `L5` in Round 1).
* **Zero claims were UNCHALLENGED across the round:** Every claim `C1`–`C8` was challenged by at least 4 of the 7 lenses. Every claim except `C1` (`SURVIVED` in L2, L4, L5, L6, L7; `NARROWED` in L1; `UNCHALLENGED` in L3) was narrowed or refuted across multiple lenses.

---

## 2. Per-Lens Coverage and Thin Lenses

| Lens | Attack File | Findings Filed | Claims Challenged | Claims UNCHALLENGED by This Lens | Adjudication Rulings |
| :--- | :--- | :---: | :--- | :--- | :--- |
| L1 (`correctness`) | `attack-L1-correctness.md` | 8 | C1–C8 | None | 8 UPHELD |
| L2 (`sampling`) | `attack-L2-sampling.md` | 8 | C1–C8 (+ framing note F8) | None (`C1`, `C8` checked in F7) | 4 UPHELD, 4 NARROWED |
| L3 (`statistics`) | `attack-L3-statistics.md` | 6 | C2, C5, C6, C7, C8 (+ steering note F6) | C1, C3, C4 | 3 UPHELD, 3 NARROWED |
| L4 (`provenance`) | `attack-L4-provenance.md` | 7 | C2–C6, C8, `R1-IMPL` (`C1`, `C7`, `R1-GFA`, `R1-HMM` in table) | None | 7 UPHELD (5 repair-narrowed by cost) |
| L5 (`inference`) | `attack-L5-inference.md` | 7 | C2–C8 (`C1` checked in table) | None | 7 UPHELD |
| L6 (`prior-art`) | `attack-L6-prior-art.md` | 9 | C1, C2, C4, C7, C8 (+ meta-finding F9) | C3, C5, C6 | 8 UPHELD, 1 NARROWED |
| L7 (`reproduction`) | `attack-L7-reproduction.md` | 6 | C2–C6, C8 (`C1`, `C7` checked in prose) | None | 5 UPHELD, 1 NARROWED |

* **Thin lenses:**
  * **L3 (`statistics`)** is the thinnest lens (6 findings; 3 of 8 claims `UNCHALLENGED`: `C1`, `C3`, `C4`). While `C1` (GFA streaming) and `C3` (binary test count) are non-statistical, **`C4` (`horizon_lift`: `469` exact, `1` partial, `19` split, `11` unmapped) has a quantitative classification claim** (`EXACT` inflation from `--bridge-bp 10`, caught by L1/L5) that L3 left unexamined.
  * **L6 (`prior-art`)** left 3 claims `UNCHALLENGED` (`C3`, `C5`, `C6`), though appropriately so because those claims make no priority or novelty assertions, and L6 compensated with meta-finding `F9` documenting missing citations (`GraphAligner`, `vg`, `odgi`, `minigraph`, `smoothxg`, `minimap2`, `last-train`, `Liftoff`).

---

## 3. Claims Nobody Could Attack (Unfalsifiable / Over-Bundled Framing)

Unlike Round 1 (which included future-tense promises from `STRATEGIC_PLAN.md`), Round 2's `C1`–`C8` are all **present-tense capability or empirical assertions** and were therefore falsifiable in principle. However, two claims suffer from **compound bundling** that blurred evaluation:

* **`C2` (`round.json:20-29`, `docs/BENCHMARK_RESULTS.md:5`):** Bundles *"All results were generated out-of-the-box"* across pairwise WGA, `--graph`, and non-HMM GFA utilities (`horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_resolve`), forcing all 7 lenses to spend their `F1`/`F2` finding dissecting which subset of tools actually uses `HmmParams`, `HmmParams::default()`, or no Pair-HMM at all.
* **`C3` (`round.json:32-41`, `docs/STRATEGIC_PLAN.md:3`):** Bundles `horizon_wga --graph` (a flag, not a binary) with 9 `src/bin/*.rs` files and cites a crate-wide `39` test count as evidence that all listed binaries are unit-tested, masking the 4 binaries with 0 binary tests.

---

## 4. Single-Party Framing Bias (`runner_is_author=true`, `claim_list_author="runner"`) & What Every Lens Missed

`round.json:137-139` records `runner_is_author: true`, `claim_list_author: "runner"`, and `human_claim: ""`. Beyond the held-out WGA regressions (`Horse 50.5%`, `Dog 63.6%`), `B10K 363-way` (`Pending Data`), and missing prior-art baselines already flagged by `L2 (F8)`, `L3 (F6)`, and `L6 (F9)`, the runner's `C1`–`C8` list omitted four concrete technical areas that **no Round 2 lens examined**:

1. **Unbenchmarked `horizon_merge` "400+ species scale" claim (`docs/BENCHMARK_RESULTS.md:9`, `src/bin/horizon_merge.rs:25`):**
   * `docs/BENCHMARK_RESULTS.md:9` claims: *"Multi-species merging (`horizon_merge`): Fast and exact at 400+ species scale, with best-score interval netting."*
   * L1, L2, L4, L5, and L7 noted that `src/bin/horizon_merge.rs` has `0` binary unit tests, **but no lens audited the `400+ species scale` claim itself**. `docs/BENCHMARK_RESULTS.md` contains zero `horizon_merge` benchmark rows or runtime/RSS measurements at 400+ species scale, and `src/bin/horizon_merge.rs:25` calls `fs::read_to_string(path)` to load every species' entire MAF file into memory simultaneously.
2. **Coverage percentages computed from many-to-many pair counts on `chr22` and `chr6` (`docs/BENCHMARK_RESULTS.md:22, 26`):**
   * `L3 (F2)` noted that `5,894,191` aligned pairs on `5 Mbp` (`C5`, `117.9%`) is a many-to-many `(x, y)` pair count and parenthetically referenced lines 22 and 26. **However, no lens directly audited the `68.5%` (`24,365,611` pairs on `35.6 Mbp` `CHM13.chr22` q-arm at line 22) and `75.1%` (`15,021,911` pairs on `20 Mbp` `chr6` q-arm at line 26) coverage percentages.** Both percentages divide many-to-many `(x, y)` pair counts by target locus length rather than deduplicating unique aligned bases, inflating apparent base coverage.
3. **`horizon_clean` `S`-tag stripping and untestable walk-splitting (`docs/PANGENOME_SUITE.md:12`, `src/bin/horizon_clean.rs:228-252`):**
   * L1–L5 and L7 noted that `horizon_clean` severed `0` links on `/tmp/chr22_peri_5m.gfa` because all links had `sup > 0` (`*`-sequence GFA). **No lens examined `src/bin/horizon_clean.rs:250`**, where Pass 2 emits `"S\t{}\t{}"` and silently strips all `S`-line optional tags (`LN:i:`, `SN:Z:`, `SO:i:`, etc.) from every surviving segment—nor noted that `walk_support` (`sup > 0` for every consecutive step pair in `P`/`W` lines) makes the advertised `P`/`W` walk-splitting path (`docs/PANGENOME_SUITE.md:12`) unreachable for any link that actually appears in a `P`/`W` walk.
4. **Incomplete tracking of carried Round 1 items (`round.json:140-159`) and cross-lens contradiction on `R1-HMM`:**
   * `round.json:140-159` carried `R1-IMPL`, `R1-GFA`, and `R1-HMM` from Round 1, yet only **L4 (`F7` + summary table)** issued verdicts on them; L1, L2, L3, L5, L6, and L7 did not adjudicate `R1-GFA` or `R1-HMM` explicitly.
   * Worse, `attack-L4-provenance.md:113` marked `R1-HMM` **`SURVIVED`** (as if resolved), directly contradicting `adjudication-L1.md` (`F2`: `mdl_xdrop_extend` still uses unnormalized hardcoded constants) and `adjudication-L3.md` (`F1` + preamble: `R1-HMM` gap-bit calibration and `mdl_xdrop_extend` remain defective).
   * Additionally, Round 1 Completeness Critic's **Gap 4** (scale-dependent recall collapse in `cactus241way` `chr22 1MB`) and **Gap 5** (`horizon_lift` strict structural-path projection failure on non-syntenic haplotypes) were not carried into `round.json`.

---

## 5. Open Gaps Summary

1. **Gap: Unverified `horizon_merge` "400+ species scale" claim (`docs/BENCHMARK_RESULTS.md:9`, `src/bin/horizon_merge.rs:25`).**
   * *Why it matters:* `docs/BENCHMARK_RESULTS.md:9` asserts fast, exact 400+ species merging in the Executive Summary with zero supporting benchmark rows, zero binary unit tests, and full-file `fs::read_to_string` MAF loading.
   * *What would close it:* Either add a reproducible 400+ species benchmark row/script for `horizon_merge` in `docs/BENCHMARK_RESULTS.md` or remove/qualify the `400+ species scale` claim in `docs/BENCHMARK_RESULTS.md:9`.
2. **Gap: Many-to-many `(x, y)` pair counts reported as `%` base coverage (`docs/BENCHMARK_RESULTS.md:22, 26`).**
   * *Why it matters:* `68.5%` on `CHM13.chr22` q-arm (`24,365,611` pairs) and `75.1%` on `20 Mbp` `chr6` q-arm (`15,021,911` pairs) are ratios of many-to-many aligned pairs to locus length, overstating true 1-to-1 base coverage.
   * *What would close it:* Compute and report unique aligned base counts alongside `(x, y)` pair counts, or remove the `%` coverage parentheses on lines 22 and 26 of `docs/BENCHMARK_RESULTS.md`.
3. **Gap: `horizon_clean` strips `S`-line optional tags and cannot split walk-supported edges (`src/bin/horizon_clean.rs:228-252`, `docs/PANGENOME_SUITE.md:12`).**
   * *Why it matters:* `horizon_clean` drops `LN:i:`/`SN:Z:`/`SO:i:` tags on all emitted `S` lines, and its `sup == 0` severing guard prevents it from ever severing or splitting a walk-supported edge.
   * *What would close it:* Preserve `S`-line tags in `src/bin/horizon_clean.rs:250` and clarify in `docs/PANGENOME_SUITE.md:12` when (if ever) walk splitting can trigger.
4. **Gap: Contradictory `R1-HMM` status (`attack-L4-provenance.md:113` vs. `adjudication-L1.md` / `adjudication-L3.md`) and stale `CRITIC-BRIEF.md:49`.**
   * *Why it matters:* Marking `R1-HMM` `SURVIVED` in L4 hides the upheld defects in `mdl_xdrop_extend` and `gap_open_bits`/`gap_extend_bits` calibration, and `CRITIC-BRIEF.md:49` falsely states `L3` has no adjudication.
   * *What would close it:* Keep `R1-HMM` open (alongside `R1-IMPL`) in `round.json` when closing Round 2, and re-render `CRITIC-BRIEF.md` / `round.json` so `adjudication-L3.md` is registered.
