# Round 2 ADJUDICATION -- Completeness Critic

**Target filing:** `reports/rounds/round-2/completeness-critic.md` (unsealed; key passages quoted below per `ADJ-BRIEF-completeness.md`).
**Prior rounds & method cards:** `scripts/round history "R1-HMM" --kind adjudication` confirms Round 1 is the only prior closed round (`round 2` withheld as open). No `refs/methods/` directory exists in this repository; rulings are grounded in direct verification of `src/`, `docs/`, `reports/rounds/round-2/round.json`, `CRITIC-BRIEF.md`, `attack-L1-correctness.md` through `attack-L7-reproduction.md`, and `adjudication-L1.md` through `adjudication-L7.md`.

---

## Finding-by-Finding Rulings

### S1-S2 -- Execution Status, `CRITIC-BRIEF.md:49` Discrepancy, and Per-Lens Coverage
> *"Zero lenses were UNRUN: All 7 declared lenses ... filed both an attack (`attack-L*.md`) and an adjudication (`adjudication-L*.md`). `CRITIC-BRIEF.md:49` states `attack-L3-statistics.md (statistics) — NO ADJUDICATION FILED`, but `reports/rounds/round-2/adjudication-L3.md` (72 lines, SHA-256 prefix `3952ea463b86` verified) exists on disk... L3 (`statistics`) is the thinnest lens (6 findings; 3 of 8 claims `UNCHALLENGED`: `C1`, `C3`, `C4`)."*

- **Verdict:** **NARROWED**
- **What survives:** All 7 declared lenses indeed filed both an attack and an adjudication, `CRITIC-BRIEF.md:49` is stale (`NO ADJUDICATION FILED` despite `adjudication-L3.md` being 72 lines on disk), and every per-lens finding count, challenged/unchallenged claim list, and verdict tally in Section 2 matches `adjudication-L1.md` through `adjudication-L7.md`.
- **Overturned detail:** The critic cites SHA-256 prefix `3952ea463b86` as `adjudication-L3.md`'s digest, whereas `sha256sum` shows `3952ea463b86...` is the digest of `attack-L3-statistics.md`; `adjudication-L3.md`'s actual SHA-256 prefix is `f5ac3a11f1cf`.
- **Repair:** Refreshing `round.json` / `CRITIC-BRIEF.md` registration for `adjudication-L3.md` (and correcting the cited SHA-256 prefix if referenced) is endorsed; **the fix is free** (metadata/text-only edit).

---

### S3 -- Compound Bundling in `C2` and `C3`
> *"`C2` (`round.json:20-29`, `docs/BENCHMARK_RESULTS.md:5`) bundles 'All results were generated out-of-the-box' across pairwise WGA, `--graph`, and non-HMM GFA utilities... `C3` (`round.json:32-41`, `docs/STRATEGIC_PLAN.md:3`) bundles `horizon_wga --graph` (a flag, not a binary) with 9 `src/bin/*.rs` files and cites a crate-wide `39` test count..."*

- **Verdict:** **UPHELD**
- **What survives:** Verified directly in `round.json:20-41`, `docs/BENCHMARK_RESULTS.md:5`, and `docs/STRATEGIC_PLAN.md:3`: both claims bundle heterogeneous components under a single assertion, forcing every lens to spend `F1`/`F2` unpacking which sub-tools use `HmmParams`, `HmmParams::default()`, or no Pair-HMM (`C2`), and which 4 of 9 binaries have `0` binary unit tests (`C3`).
- **Repair:** Endorsed (split compound capability/provenance claims in future claim lists and qualify `docs/BENCHMARK_RESULTS.md:5` / `docs/STRATEGIC_PLAN.md:3`). **The fix is free** (text-only edit).

---

### Gap 1 -- Unverified `horizon_merge` "400+ species scale" Claim (`docs/BENCHMARK_RESULTS.md:9`, `src/bin/horizon_merge.rs:25`)
> *"`docs/BENCHMARK_RESULTS.md:9` claims: 'Multi-species merging (`horizon_merge`): Fast and exact at 400+ species scale, with best-score interval netting.' ... no lens audited the `400+ species scale` claim itself. `docs/BENCHMARK_RESULTS.md` contains zero `horizon_merge` benchmark rows or runtime/RSS measurements at 400+ species scale, and `src/bin/horizon_merge.rs:25` calls `fs::read_to_string(path)` to load every species' entire MAF file into memory simultaneously."*

- **Verdict:** **UPHELD**
- **What survives:** Verified in `docs/BENCHMARK_RESULTS.md:9` (claim present, zero `horizon_merge` rows in any benchmark table), `src/bin/horizon_merge.rs:25` (`fs::read_to_string(path)` loads each species' full MAF into `inputs: Vec<(String, String)>`), and across all `attack-L*.md` and `adjudication-L*.md` files (`400+` appears nowhere; L1, L2, L4, L5, and L7 only noted `0` binary unit tests).
- **Repair:** Endorsed on diagnosis; split by cost. Removing or qualifying the unbenchmarked `400+ species scale` claim in `docs/BENCHMARK_RESULTS.md:9` **is free** (text-only edit), whereas adding a reproducible 400+ species benchmark row/script for `horizon_merge` is **not free**.

---

### Gap 2 -- Many-to-Many `(x, y)` Pair Counts Reported as `%` Base Coverage (`docs/BENCHMARK_RESULTS.md:22, 26`)
> *"no lens directly audited the `68.5%` (`24,365,611` pairs on `35.6 Mbp` `CHM13.chr22` q-arm at line 22) and `75.1%` (`15,021,911` pairs on `20 Mbp` `chr6` q-arm at line 26) coverage percentages. Both percentages divide many-to-many `(x, y)` pair counts by target locus length rather than deduplicating unique aligned bases, inflating apparent base coverage."*

- **Verdict:** **UPHELD**
- **What survives:** Verified in `docs/BENCHMARK_RESULTS.md:22, 26` (`24,365,611 / 35.6M = 68.44%`; `15,021,911 / 20M = 75.11%`) and across all `attack-L*.md` and `adjudication-L*.md` files (`68.5` and `75.1` appear nowhere; `L3 (F2)` only audited line 23's `5,894,191` pairs on `5 Mbp`).
- **Repair:** Endorsed on diagnosis; split by cost. Removing or qualifying the `%` coverage parentheses on `docs/BENCHMARK_RESULTS.md:22, 26` as pair-count-to-locus-length ratios **is free** (text-only edit), whereas computing and reporting unique aligned base counts requires re-running `horizon_wga` and is **not free**.

---

### Gap 3 -- `horizon_clean` `S`-Tag Stripping and Unreachable Walk-Splitting (`src/bin/horizon_clean.rs:228-252`, `docs/PANGENOME_SUITE.md:12`)
> *"No lens examined `src/bin/horizon_clean.rs:250`, where Pass 2 emits `\"S\\t{}\\t{}\"` and silently strips all `S`-line optional tags (`LN:i:`, `SN:Z:`, `SO:i:`, etc.) from every surviving segment—nor noted that `walk_support` (`sup > 0` for every consecutive step pair in `P`/`W` lines) makes the advertised `P`/`W` walk-splitting path (`docs/PANGENOME_SUITE.md:12`) unreachable for any link that actually appears in a `P`/`W` walk."*

- **Verdict:** **UPHELD**
- **What survives:** Verified in `src/bin/horizon_clean.rs:180-200, 228-237, 248-252`: Pass 2 emits `"S\t{}\t{}"` without optional tags, and Pass 1 increments `walk_support` for every consecutive step pair in `P`/`W` lines while `should_sever` requires `sup == 0`, so no `L` link appearing in a `P`/`W` walk can ever be severed or trigger walk splitting. Zero `attack-L*.md` or `adjudication-L*.md` files examined `S`-tag stripping or the unreachable `P`/`W` walk-splitting claim in `docs/PANGENOME_SUITE.md:12`.
- **Repair:** Endorsed on diagnosis; split by cost. Clarifying in `docs/PANGENOME_SUITE.md:12` when (if ever) walk splitting can trigger **is free** (text-only edit), whereas preserving `S`-line optional tags in `src/bin/horizon_clean.rs:250` requires a code change and unit test (**not free**).

---

### Gap 4 -- Contradictory `R1-HMM` Status (`attack-L4-provenance.md:113` vs. `adjudication-L1.md`/`L3.md`), Stale `CRITIC-BRIEF.md:49`, and Uncarried Round 1 Gaps
> *"`attack-L4-provenance.md:113` marked `R1-HMM` `SURVIVED` (as if resolved), directly contradicting `adjudication-L1.md` (`F2`: `mdl_xdrop_extend` still uses unnormalized hardcoded constants) and `adjudication-L3.md` (`F1` + preamble: `R1-HMM` gap-bit calibration and `mdl_xdrop_extend` remain defective). Additionally, Round 1 Completeness Critic's Gap 4 ... and Gap 5 ... were not carried into `round.json`."*

- **Verdict:** **UPHELD**
- **What survives:** Verified in `attack-L4-provenance.md:113` (`| R1-HMM | SURVIVED |`), `adjudication-L1.md:18-22`, `adjudication-L3.md:5, 14`, `CRITIC-BRIEF.md:49`, and `round.json:140-159` (only `R1-IMPL`, `R1-GFA`, and `R1-HMM` were carried; Round 1 Completeness Critic Gaps 4 and 5 were not carried).
- **Repair:** Endorsed (keep `R1-HMM` open alongside `R1-IMPL` in `round.json` when closing Round 2, and re-render `CRITIC-BRIEF.md` / `round.json` so `adjudication-L3.md` is registered). **The fix is free** (metadata/protocol edit).

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| S1-S2 | NARROWED | All 7 lenses filed attacks/adjudications, `CRITIC-BRIEF.md:49` is stale on `adjudication-L3.md`, and Section 2 counts match disk, but `3952ea463b86` is `attack-L3-statistics.md`'s digest (not `adjudication-L3.md`'s `f5ac3a11f1cf`); metadata fix is free |
| S3 | UPHELD | `C2` and `C3` bundle heterogeneous tools/flags and crate-wide test counts, forcing every lens to unpack them in `F1`/`F2`; doc/claim scoping is free |
| Gap 1 | UPHELD | `docs/BENCHMARK_RESULTS.md:9` claims `horizon_merge` at `400+ species scale` with zero benchmark rows and `fs::read_to_string` full-MAF loading (`src/bin/horizon_merge.rs:25`), unexamined by all 7 lenses; doc qualification is free, benchmarking is not |
| Gap 2 | UPHELD | `68.5%` (`chr22` q-arm, line 22) and `75.1%` (`chr6` q-arm, line 26) divide many-to-many `(x, y)` pair counts by locus length and went unexamined by all 7 lenses; doc qualification is free, unique-base re-measurement is not |
| Gap 3 | UPHELD | `horizon_clean.rs:250` strips all `S`-line optional tags and `sup > 0` (`lines 228-237`) makes `P`/`W` walk splitting (`docs/PANGENOME_SUITE.md:12`) unreachable for walk-supported `L` links, unexamined by all 7 lenses; doc fix is free, code fix is not |
| Gap 4 | UPHELD | `attack-L4-provenance.md:113` marks `R1-HMM` `SURVIVED` against `L1`/`L3` upheld defects, `CRITIC-BRIEF.md:49` omits `adjudication-L3.md`, and Round 1 Gaps 4–5 were not carried; metadata/closeout fix is free |
