# Round 3 ADJUDICATION -- Completeness Critic

**Target filing:** `reports/rounds/round-3/completeness-critic.md` (unsealed; ruled passages quoted below per `ADJ-BRIEF-completeness.md`).  
**Verified against:** `reports/rounds/round-3/round.json`, `reports/rounds/round-3/CARRIED.md`, `reports/rounds/round-3/ATTACKED-STRATEGIC_PLAN.md`, `reports/rounds/round-3/ATTACKED-BENCHMARK_RESULTS.md`, all 7 `attack-L*.md` and `adjudication-L*.md` reports, `git show --stat 488b020`, and source files `src/bin/horizon_wga.rs`, `src/bin/horizon_resolve.rs`, `src/bin/horizon_clean.rs`, `src/gfa/mod.rs`.

**Prior round & method-card check:**
- `./scripts/round history "completeness" --kind adjudication` returned Round 1 and Round 2 completeness adjudications (plus Round 2 `adjudication-L3.md` F6 on claim-list steering); none of Round 3's three specific coverage gaps (`estimate_dag_hmm_params`, `horizon_resolve` `L`-link rewriting, or 4 of 5 `C6` benchmarks) were adjudicated in Rounds 1–2.
- `refs/methods/` does not exist in this repository, and `Gap 1`–`Gap 3` are coverage, code-inspection, and reproduction-scope findings rather than statistical disputes, so no method card governs them.

---

## Coverage & Framing Check (Sections 1–2 of `completeness-critic.md`)

- **Lenses:** Confirmed all 7 declared lenses (`L1`–`L7`) ran, filed, and were adjudicated; zero lenses were `UNRUN`, `aborted`, or `dropped`.
- **Compound-claim bundling & runner-authored claim list:** Confirmed `C1`–`C4` and `C6` bundle multiple sub-claims under a single `round.json` verdict, and `C1`–`C6` omit the 5 stale pre-`488b020` WGA rows (`R2-HMM-XDROP`), absent `minimap2`/`odgi`/`vg` baselines (`R2-RESOLVE-LIFT`), and unfinished full `chr6` alignment (`R2-INDEX-RAM`)—already captured by `L1`–`L7` framing findings.

---

## Rulings on Substantive Gaps

### Gap 1 -- `estimate_dag_hmm_params` (`C2`, `src/bin/horizon_wga.rs:367-437`) was never inspected by any lens

> *"Direct inspection of `src/bin/horizon_wga.rs:367-437` shows that: (1) Unlike `estimate_hmm_params`, `estimate_dag_hmm_params` only estimates `identity_prior` (lines 432–436) and leaves `gap_open_bits` and `gap_extend_bits` at `HmmParams::default()`. (2) It samples up to `total >= 2000` bases (lines 411–426) in sorted hash-key order (`keys.sort_unstable()`, line 382) from exact `MASK_W10` seed windows (`canonical_spaced_seeds_with_opts`, lines 370, 377), which still allocates `HashMap`s and biases `identity_prior` upward by conditioning on seed-matching windows. (3) It has no unit test in `src/bin/horizon_wga.rs`."*

- **Verdict:** **UPHELD**
- **Verification:**
  - Coverage: `grep` across all `attack-L*.md` and `adjudication-L*.md` confirms only `attack-L2-sampling.md:12` mentioned `estimate_dag_hmm_params` in passing; no lens inspected `src/bin/horizon_wga.rs:367-437`.
  - Code: `src/bin/horizon_wga.rs:433-436` returns `HmmParams { identity_prior, ..HmmParams::default() }`, leaving `gap_open_bits` and `gap_extend_bits` at defaults; lines `370, 377, 382, 411-426` sample `<= 2000` bases (actually stopping when `total >= 2000`) from `MASK_W10` seed hits (`canonical_spaced_seeds_with_opts`) in sorted key order; and `src/bin/horizon_wga.rs:578-614` contains only `load_gfa_dag_handles_bidirected_minus_links`, with no test for `estimate_dag_hmm_params`.
  - Claim surface: `docs/STRATEGIC_PLAN.md:11` and `docs/BENCHMARK_RESULTS.md:5` both list `horizon_wga --graph` (`estimate_dag_hmm_params`) under the full adaptive `HmmParams` (`identity_prior`, `gap_open_bits`, `gap_extend_bits`) description without qualification.
- **What survives:** `estimate_dag_hmm_params` went uninspected by all 7 lenses and adapts only `identity_prior` (from `<= 2000` bp of `MASK_W10` seed hits, with no unit test), while `docs/STRATEGIC_PLAN.md:11` and `docs/BENCHMARK_RESULTS.md:5` imply full adaptive `HmmParams`. (Lone critic catch across the round.)
- **Repair ruling & cost:** Endorsed in two parts: (1) qualifying `docs/STRATEGIC_PLAN.md:11` and `docs/BENCHMARK_RESULTS.md:5` to state that `estimate_dag_hmm_params` only adapts `identity_prior` from `MASK_W10` seed hits (leaving `gap_open_bits`/`gap_extend_bits` at defaults) is **free** (doc-only edit); (2) adding a unit test in `src/bin/horizon_wga.rs` is **not free** (costs code authoring and `cargo test`).

---

### Gap 2 -- `horizon_resolve` `W`-line/`L`-line rewriting (`src/bin/horizon_resolve.rs:416-480`) and `horizon_clean` non-`S` tag preservation (`src/bin/horizon_clean.rs:276-330`) in `C3` were never inspected by any lens

> *"1. Carried item `R2-RESOLVE-LIFT` explicitly required `'Update W-line start/end coordinates and L links in horizon_resolve --collapse-dups'`. In `src/bin/horizon_resolve.rs:478`, `L` lines are copied verbatim (`writeln!(out, "{trimmed}")?`) without updating or adding links when `--collapse-dups` collapses segments in `P`/`W` lines. 2. In `src/bin/horizon_clean.rs:276-330`, `S`-line tags are preserved (`s.tags`), but `L`-line, `P`-line, and `W`-line optional tags are dropped when rewriting."*

- **Verdict:** **NARROWED**
- **Verification:**
  - Coverage: `grep` across all `attack-L*.md` and `adjudication-L*.md` confirms zero lenses inspected `horizon_resolve.rs:416-480` or `horizon_clean.rs:276-330`.
  - Diagnosis leg 1 (`horizon_resolve` `L` links): **UPHELD**. `src/bin/horizon_resolve.rs:416-477` rewrites `W` lines (updating `end_str` while keeping `start` `cols[4]`, and preserving `cols[7..]` tags) and `P` lines (`cols[4..]`), while `L` lines fall through to line `478` (`writeln!(out, "{trimmed}")?`) verbatim. This leaves `R2-RESOLVE-LIFT`'s explicit `"and L links in horizon_resolve --collapse-dups"` requirement unfulfilled.
  - Diagnosis leg 2 (`horizon_clean` `L`/`P`/`W` optional tags): **OVERTURNED as a `C3`/`R2-RESOLVE-LIFT` defect**. Checking `src/gfa/mod.rs:14-47` shows `GfaLink`, `GfaPath`, and `GfaWalk` do not parse optional tags at all (only `GfaSegment` has `tags: Option<String>`), and both `CARRIED.md` (`R2-RESOLVE-LIFT`) and `round.json` (`C3`) only required/claimed `S-tag preservation`. The critic should have checked `GfaLink`/`GfaPath`/`GfaWalk` in `src/gfa/mod.rs` and the exact `R2-RESOLVE-LIFT` wording before treating `L`/`P`/`W` tag dropping as an unfulfilled `C3` claim.
- **What survives:** All 7 lenses missed inspecting `horizon_resolve.rs:416-480`, where `L` lines are copied verbatim (`line 478`) despite `R2-RESOLVE-LIFT` requiring `L`-link updates in `horizon_resolve --collapse-dups`; `horizon_clean` dropping `L`/`P`/`W` tags is a `stream_gfa` struct limitation rather than a `C3` contradiction. (Lone critic catch on `horizon_resolve` `L` links.)
- **Repair ruling & cost:** Keeping `R2-RESOLVE-LIFT` open until `horizon_resolve --collapse-dups` updates `L` links (or qualifying `C3` in docs) is endorsed: qualifying `C3` in docs is **free** (doc-only edit), whereas implementing `L`-link updates in `src/bin/horizon_resolve.rs` is **not free** (costs code + test authoring).

---

### Gap 3 -- 4 of 5 benchmarks bundled into `C6` (`horizon_gfa extract`, `horizon_lift`, `horizon_resolve`, `horizon_wga --graph`) went unreproduced by L7 and unverified by L1

> *"Neither L7 nor L1 re-executed or verified `horizon_gfa extract` (`456.1/1,255.1 MB`), `horizon_lift` (`256.60s`, `619.2 MB`), `horizon_resolve --collapse-dups` (`4.48s`, `24.0 MB`), or `horizon_wga --graph` (`96/96 bp`, `0.0031s`). Commit `488b020` modified `src/bin/horizon_gfa.rs` (68 lines), `src/bin/horizon_lift.rs` (8 lines), `src/bin/horizon_resolve.rs` (53 lines), and `src/bin/horizon_wga.rs` (109 lines). Marking `C6` as `SURVIVED` in L7 based solely on `horizon_gfa stats` leaves 4 of the 5 bundled benchmarks unverified on `HEAD` (`488b020`)."*

- **Verdict:** **UPHELD**
- **Verification:**
  - `attack-L7-reproduction.md` (`F4`) and `adjudication-L7.md` (`F4`) reproduced only `horizon_gfa stats` (`88.84s`, `97.1 MB`) and recorded `C6` as `SURVIVED`.
  - `attack-L1-correctness.md` (`F6`) inspected only `cmd_stats` (`src/bin/horizon_gfa.rs:28-76`) and marked `C6` `UNCHALLENGED`.
  - `git show --stat 488b020` confirms `488b020` modified `src/bin/horizon_gfa.rs` (68 lines), `src/bin/horizon_lift.rs` (8 lines), `src/bin/horizon_resolve.rs` (53 lines), and `src/bin/horizon_wga.rs` (109 lines), while `adjudication-L4.md` (`F5`) already noted `docs/BENCHMARK_RESULTS.md:54-57, 64-65` were not updated after `488b020`.
- **What survives:** `L7`'s `SURVIVED` and `L1`'s `UNCHALLENGED` on `C6` cover only `horizon_gfa stats`, leaving `horizon_gfa extract`, `horizon_lift`, `horizon_resolve --collapse-dups`, and `horizon_wga --graph` unverified on `488b020`.
- **Repair ruling & cost:** Endorsed in two parts: (1) splitting `C6` into per-benchmark claims in future rounds is **free** (claim-list structuring); (2) re-running `horizon_gfa extract`, `horizon_lift`, `horizon_resolve --collapse-dups`, and `horizon_wga --graph` on `488b020` with saved logs is **not free** (costs benchmark execution on the 12 GB HPRC graph and subgraphs).

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| Gap 1 | UPHELD | `estimate_dag_hmm_params` (`src/bin/horizon_wga.rs:367-437`) was uninspected by all 7 lenses, adapts only `identity_prior` from `<=2000` bp of `MASK_W10` seed hits (leaving gap penalties at defaults), and has no unit test; doc qualification is free, unit test costs authoring |
| Gap 2 | NARROWED | `horizon_resolve.rs:416-480` went uninspected and copies `L` lines verbatim (`line 478`) despite `R2-RESOLVE-LIFT` requiring `L`-link updates (doc qualification free, code fix non-free); `horizon_clean` dropping `L`/`P`/`W` tags is overturned as a `C3` defect because `C3`/`R2-RESOLVE-LIFT` only claimed `S`-tag preservation and `stream_gfa` (`src/gfa/mod.rs:21-47`) lacks `L`/`P`/`W` tag fields |
| Gap 3 | UPHELD | `L7` (`F4`) and `L1` (`F6`) checked only `horizon_gfa stats`/`cmd_stats` before marking `C6` `SURVIVED`/`UNCHALLENGED`, leaving `extract`, `lift`, `resolve`, and `--graph` unreproduced on `488b020` despite `488b020` editing all 4 binaries; splitting `C6` is free, re-running benchmarks costs execution |
