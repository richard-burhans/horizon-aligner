# Round 4 ADJUDICATION -- Completeness Critic

**Target filing:** `reports/rounds/round-4/completeness-critic.md` (unsealed; key passages quoted below per `ADJ-BRIEF-completeness.md`).
**Method cards (`refs/methods/`):** None exist in this repository; rulings rest on direct verification against `reports/rounds/round-4/round.json`, `reports/rounds/round-4/CARRIED.md`, `reports/rounds/round-4/ATTACKED-horizon_paper.md`, `reports/rounds/round-4/attack-L1-*.md`..`attack-L7-*.md`, `reports/rounds/round-4/adjudication-L1.md`..`adjudication-L7.md`, `docs/BENCHMARK_RESULTS.md`, `src/`, and `manuscript/`.
**Prior rounds checked (`scripts/round history "completeness" --kind adjudication`):** Rounds 1–3 completeness adjudications established the distinction between `UNCHALLENGED` and `UNRUN`, flagged unexecuted `horizon_gfa extract` in Round 3 (`adjudication-completeness.md:49`), and tracked carried items `R2-HMM-XDROP` and `R2-RESOLVE-LIFT`.

---

## Sections 1–3: Coverage, `C1` Vagueness, and `round.json` Framing

- **Section 1 (Per-Lens Counts, `UNCHALLENGED` vs. `UNRUN`, and Thin Lenses):** **UPHELD.**
  - *Quoted assertion:* `"All seven declared lenses ... ran and filed sealed reports that were adjudicated (L1–L7). Zero lenses were UNRUN."` / `"L7 (reproduction) is the thinnest lens on execution coverage."` / `"L2 (sampling) and L3 (statistics) marked C1 UNCHALLENGED"`
  - *Verification:* Every finding count, verdict tally, and `UNCHALLENGED` status in Section 1's table matches `adjudication-L1.md` through `adjudication-L7.md` (`L1`: 5 `UPHELD`, 2 `NARROWED`; `L2`: 6 `UPHELD`, `C1` `UNCHALLENGED`; `L3`: 6 `UPHELD`, 1 `NARROWED`, `C1` `UNCHALLENGED`; `L4`: 7 `UPHELD`; `L5`: 6 `UPHELD`, 1 `NARROWED`; `L6`: 7 `UPHELD`; `L7`: 4 `UPHELD`, 2 `NARROWED`). `L7` indeed left `horizon_gfa extract` (`chr22`, `chr6`), `horizon_lift`, `horizon_resolve`, `horizon_wga --graph`, `horizon_guide`, and full `chr6` unexecuted on `HEAD` (`33bede8`).
- **Section 2 (`C1` Underspecified Claim Design):** **UPHELD.**
  - *Quoted assertion:* `"'Unified' and 'composable' have no operational definition in C1 or the manuscript. L2 and L3 had to mark C1 UNCHALLENGED ... while L1, L4, L5, L6, and L7 could only attack the literal executable count (9 vs. 10), library imports, and stdin (-) support."`
  - *Verification:* Confirmed in `round.json:15-19` and `ATTACKED-horizon_paper.md:18`.
- **Section 3 (`round.json` Framing Defects):** **UPHELD.**
  - *Quoted assertion:* `"Overloaded compound claims + single-clause source.quote anchoring"` (`C3`, `C4`, `C5`, `C6`) and `"The Abstract (ATTACKED-horizon_paper.md:10) and Results (line 59) highlight the 4-regime synthetic Plasmodium benchmark ... yet C1–C6 omit it entirely."`
  - *Verification:* Confirmed in `round.json:40-79`: `C3`–`C6` summaries bundle `horizon_lift`, `horizon_resolve`, `horizon_wga --graph`, *P. malariae* `50.2×`, and `horizon_guide --progressive --hal`, while each `source.quote` quotes only the opening clause, and all 7 attackers anchored on the quoted clause.

---

## Section 4: Specific Unattacked Gaps

### Gap 1 -- Synthetic 4-regime *Plasmodium* benchmark (`0.978` recall, `0.949` precision, `0.000` vs. `0.872` `HOXD70`) in Abstract and Results

- **Quoted finding:** `"Artifact: ATTACKED-horizon_paper.md:10, 59; absent from all attack-L1..attack-L7 reports and absent from src/, scripts/, and docs/BENCHMARK_RESULTS.md (verified via grep -rn \"0\\.978\\|0\\.936\\|0\\.872\\|HOXD70\" returning only docs/manual.md:31)."`
- **Verdict:** **UPHELD**
- **What survives:** The 4-regime synthetic *Plasmodium* benchmark (`0.978` recall, `0.949` precision, `0.129` vs. `0.936`, `0.000` vs. `0.872` `HOXD70`) at `ATTACKED-horizon_paper.md:10, 59` was omitted from `C1`–`C6`, went completely unattacked across `L1`–`L7`, and has zero code, script, or benchmark-table provenance in `src/`, `scripts/`, or `docs/BENCHMARK_RESULTS.md`.
- **Verification & minor note:** Verified across `attack-L1`..`attack-L7`, `src/`, `scripts/`, and `docs/BENCHMARK_RESULTS.md`. Note that the critic's `grep` statement omitted `manuscript/` matches (`manuscript/scripts/generate_figures.py:340, 488, 517, 524`, `manuscript/horizon_paper.md:10, 59`, and companion guides), where `HOXD70 false = 0.872` is hardcoded as a plot literal rather than computed.
- **Repair & cost:** Removing the unbacked `0.978` / `0.949` / `0.129` vs. `0.936` / `0.000` vs. `0.872` figures from `manuscript/horizon_paper.md:10, 59` (and `manuscript/scripts/generate_figures.py:488`) **is free**; committing a reproducible benchmark script/seeds and `docs/BENCHMARK_RESULTS.md` row (plus auditing it in the next round) is **not free**.

---

### Gap 2 -- *P. malariae* × *P. brasilianum* `50.2×` speedup (`17.12 min` vs. `14.3 h` / `859.9 min`, `536,041` bubbles) and `Pf` × `Pr` (`68.12 s`)

- **Quoted finding:** `"Artifact: ATTACKED-horizon_paper.md:10, 69, 73 and round.json:66 (C5 summary); zero mentions in attack-L1..attack-L7 or docs/BENCHMARK_RESULTS.md (verified via grep -rn \"50\\.2\\|malariae\\|859\\.9\\|68\\.12\" returning 0 matches)."`
- **Verdict:** **UPHELD**
- **What survives:** Despite appearing in `C5`'s `round.json:66` summary, the Abstract (`line 10`), Figure 3B (`line 69`), and Results (`line 73`), Horizon's headline whole-genome speedup (`50.2×`, `17.12 min` vs. `14.3 h` / `859.9 min`, `536,041` bubbles) and `Pf` × `Pr` (`68.12 s`) received zero scrutiny across all 7 lenses and have zero provenance in `src/`, `scripts/`, or `docs/BENCHMARK_RESULTS.md` (appearing only as hardcoded plot literals in `manuscript/scripts/generate_figures.py:243-255`).
- **Repair & cost:** Removing or qualifying the unbacked `50.2×` (`17.12 min` vs. `14.3 h`) and `68.12 s` figures in `manuscript/horizon_paper.md:10, 69, 73` **is free**; splitting them into a dedicated claim and recording reproducible `docs/BENCHMARK_RESULTS.md` / `scripts/` provenance with Galaxy `Growler UCSC pair` parity is **not free**.

---

### Gap 3 -- `horizon_guide --progressive --hal` (`194,878` MAF blocks, `+105.0%` `mus_musculus` coverage, `272.03 s`, `libhdf5`-validated `.hal`) and `src/hal/` (`sars-cov-2.hal` validation)

- **Quoted finding:** `"Artifact: ATTACKED-horizon_paper.md:75, 99, 111 and round.json:76 (C6 second clause); 194,878, 1,012,274, 493,677, 99,884, and libhdf5 validation were not examined in any attack-L1..attack-L7 report (only attack-L4-provenance.md:103 noted 272.03 s in framing F7), and zero matches exist in src/, scripts/, or docs/BENCHMARK_RESULTS.md."`
- **Verdict:** **UPHELD**
- **What survives:** All `horizon_guide --progressive --hal` empirical numbers (`194,878`, `162,205`, `1,012,274`, `493,677`, `+105.0%`, `99,884`) and `libhdf5` / `sars-cov-2.hal` validation claims went completely unexamined across `L1`–`L7` (with only `attack-L4-provenance.md:103` noting `272.03 s` under `F7`), have zero matches in `src/`, `scripts/`, or `docs/BENCHMARK_RESULTS.md`, and `L7`'s `cargo test` run left `real_hal_decode_from_env` `ignored` while `33bede8` modified `src/bin/horizon_guide.rs` and `src/hal/*`.
- **Repair & cost:** Removing/qualifying the unlogged `horizon_guide --progressive --hal` numbers in `manuscript/horizon_paper.md:75, 99` **is free**; splitting `horizon_guide --progressive --hal` and `src/hal/` into a dedicated claim, logging the `60_amniotes` run in `docs/BENCHMARK_RESULTS.md` / `scripts/`, and running `L1`/`L4`/`L7` checks on `33bede8` is **not free**.

---

### Gap 4 -- `horizon_lift --bridge-bp 10` (`182.67 s`, `619.7 MB`, `472` exact / `5` partial / `15` split / `8` unmapped) in `C3`

- **Quoted finding:** `"Artifact: ATTACKED-horizon_paper.md:48, 101, round.json:40, docs/BENCHMARK_RESULTS.md:21, 58; horizon_lift was only mentioned in C1 import-list checks (attack-L1, L4, L5, L7) and never audited under C3 in attack-L1..attack-L7."`
- **Verdict:** **UPHELD**
- **What survives:** `horizon_lift --bridge-bp 10` (`C3`) was never audited by `L1`–`L7`, leaving unflagged that `ATTACKED-horizon_paper.md:101` omits both the `+`-strand `GRCh38.chr22:16.0M-49.9M` q-arm restriction and *"no external coordinate truth oracle"* caveat disclosed in `docs/BENCHMARK_RESULTS.md:21`, and `L7` did not re-run it after `33bede8` modified `70` lines of `src/bin/horizon_lift.rs`.
- **Repair & cost:** Qualifying `manuscript/horizon_paper.md:101` with the `+`-strand `16.0M–49.9M` q-arm scope and absence of an external coordinate truth oracle **is free**; auditing `src/bin/horizon_lift.rs` and re-running the 500-window `chr22` liftover on `33bede8` is **not free**.

---

### Gap 5 -- `horizon_wga --graph` (`src/mdl/dag.rs`) and `horizon_resolve` (`src/bin/horizon_resolve.rs`) in `C4`, plus carried items `R2-HMM-XDROP` and `R2-RESOLVE-LIFT`

- **Quoted finding:** `"Artifact: round.json:53 (C4), CARRIED.md:5-11 (R2-HMM-XDROP, R2-RESOLVE-LIFT), ATTACKED-horizon_paper.md:50-51; attack-L1..attack-L7 on C4 examined only horizon_clean, leaving horizon_wga --graph (src/mdl/dag.rs), horizon_resolve (W-line coordinate and L-link rewriting), N-base handling in bocpd_gc_profile, mdl_xdrop_extend, and per-binary test coverage unverified on 33bede8."`
- **Verdict:** **NARROWED**
- **What survives:** `L1`–`L7` on `C4` examined only `horizon_clean` and did **not** verify the `33bede8` code fixes in `src/mdl/dag.rs`, `src/bin/horizon_resolve.rs` (`W`/`L` rewriting), `src/bocpd.rs` (`N`-base handling), or `mdl_xdrop_extend_with_params` (`src/mdl/dp.rs`), nor per-binary test depth; narrowed only because `L4` (`F7`) and `L5` (`F7`) already flagged the synthetic-only `horizon_resolve`/`--graph` benchmark status (`docs/BENCHMARK_RESULTS.md:67-68`) and open `R2-HMM-XDROP`/`R2-RESOLVE-LIFT` baseline obligations, and `L7` confirmed `cargo test` passes `14` tests across the `10` binaries.
- **Repair & cost:** Running targeted `L1`/`L7` code and test verification on `src/mdl/dag.rs`, `src/mdl/dp.rs`, `src/bocpd.rs`, and `src/bin/horizon_resolve.rs` on `33bede8` is **not free** (next-round audit/execution).

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| Gap 1 | UPHELD | Synthetic 4-regime *Plasmodium* (`0.978`/`0.949`, `0.000` vs `0.872` `HOXD70`) in `horizon_paper.md:10, 59` was omitted from `C1`–`C6`, unattacked by `L1`–`L7`, and absent from `src/`/`scripts/`/`docs/BENCHMARK_RESULTS.md` (free to remove unbacked figures; non-free to reproduce). |
| Gap 2 | UPHELD | *P. malariae* `50.2×` (`17.12 min` vs `14.3 h`, `536,041` bubbles) and `Pf` × `Pr` (`68.12 s`) in `horizon_paper.md:10, 69, 73` / `C5` summary received zero `L1`–`L7` scrutiny and lack repo provenance (free to remove/qualify; non-free to reproduce). |
| Gap 3 | UPHELD | `horizon_guide --progressive --hal` (`194,878` blocks, `+105.0%` coverage) and `src/hal/` (`sars-cov-2.hal` ignored in `cargo test`) went unaudited by `L1`–`L7` and lack `docs/BENCHMARK_RESULTS.md`/`scripts/` provenance (free to remove/qualify; non-free to audit/re-run). |
| Gap 4 | UPHELD | `horizon_lift --bridge-bp 10` (`C3`) was never audited under `C3`, omitting `docs/BENCHMARK_RESULTS.md:21`'s `+`-strand `16.0M–49.9M` q-arm scope and no-truth-oracle caveats from `horizon_paper.md:101` and unrun after `33bede8` edits (qualifying text is free; re-running is non-free). |
| Gap 5 | NARROWED | `33bede8` code fixes in `src/mdl/dag.rs`, `horizon_resolve.rs`, `bocpd.rs` (`N`-bases), and `mdl_xdrop_extend_with_params` went unverified by `L1`–`L7`, though `L4`/`L5` `F7` already flagged synthetic-only `horizon_resolve`/`--graph` benchmarks and open `R2` baseline items and `L7` ran `cargo test` (14 binary tests). |
