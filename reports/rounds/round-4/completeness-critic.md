# Round 4 Completeness Critic Report

**Artifacts read before asserting any gap:**
- `reports/rounds/round-4/CRITIC-BRIEF.md`
- `reports/rounds/round-4/round.json`
- `reports/rounds/round-4/CARRIED.md`
- `reports/rounds/round-4/ATTACKED-horizon_paper.md` (identical to `manuscript/horizon_paper.md`)
- All seven attack reports: `attack-L1-correctness.md`, `attack-L2-sampling.md`, `attack-L3-statistics.md`, `attack-L4-provenance.md`, `attack-L5-inference.md`, `attack-L6-prior-art.md`, `attack-L7-reproduction.md`
- All seven adjudication reports: `adjudication-L1.md`, `adjudication-L2.md`, `adjudication-L3.md`, `adjudication-L4.md`, `adjudication-L5.md`, `adjudication-L6.md`, `adjudication-L7.md`
- `docs/BENCHMARK_RESULTS.md`

---

## 1. Per-Lens Filing Counts, `UNCHALLENGED` vs. `UNRUN`, and Thin Lenses

All seven declared lenses (`correctness`, `sampling`, `statistics`, `provenance`, `inference`, `prior-art`, `reproduction`) **ran** and filed sealed reports that were adjudicated (`L1`–`L7`). **Zero lenses were `UNRUN`.**

| Lens | Report | Findings Filed | Adjudicated Verdicts | `UNCHALLENGED` Claims | `UNRUN` Checks |
| :--- | :--- | :---: | :--- | :--- | :--- |
| `L1` (`correctness`) | `attack-L1-correctness.md` / `adjudication-L1.md` | 7 (`F1`–`F6` on `C1`–`C6`, `F7` framing) | 5 `UPHELD`, 2 `NARROWED` (`F1`, `F3`) | None | Did not inspect `horizon_lift`, `horizon_resolve`, `horizon_wga --graph`, `horizon_guide`, `src/hal/`, or `R2-HMM-XDROP` fixes |
| `L2` (`sampling`) | `attack-L2-sampling.md` / `adjudication-L2.md` | 6 (`F1`–`F5` on `C2`–`C6`, `F6` framing) | 6 `UPHELD` | `C1` (architectural summary; no empirical sampling assertion) | Did not sample/audit `horizon_lift` (`C3`), `horizon_guide` (`C6`), or synthetic *Plasmodium* (`line 59`) |
| `L3` (`statistics`) | `attack-L3-statistics.md` / `adjudication-L3.md` | 7 (`F1` explicit `UNCHALLENGED` on `C1`, `F2`–`F6` on `C2`–`C6`, `F7` framing) | 6 `UPHELD`, 1 `NARROWED` (`F7` CRC math) | `C1` (no statistical estimator/metric) | Did not audit `horizon_lift` or `horizon_guide` statistics |
| `L4` (`provenance`) | `attack-L4-provenance.md` / `adjudication-L4.md` | 7 (`F1`–`F6` on `C1`–`C6`, `F7` framing) | 7 `UPHELD` | None | Did not audit `horizon_guide` `194,878`/`+105.0%` (`C6`) or *P. malariae* `50.2×` (`C5`) provenance |
| `L5` (`inference`) | `attack-L5-inference.md` / `adjudication-L5.md` | 7 (`F1`–`F6` on `C1`–`C6`, `F7` framing) | 6 `UPHELD`, 1 `NARROWED` (`F1` count phrasing) | None | Did not audit `horizon_wga --graph` (`C4`) or `horizon_guide` (`C6`) inference |
| `L6` (`prior-art`) | `attack-L6-prior-art.md` / `adjudication-L6.md` | 7 (`F1`–`F6` on `C1`–`C6`, `F7` bibliography hygiene) | 7 `UPHELD` | None | None on `C1`–`C6` core; `horizon_lift`/`horizon_guide` prior art not separately checked |
| `L7` (`reproduction`) | `attack-L7-reproduction.md` / `adjudication-L7.md` | 6 (`F1`–`F6` on `C1`–`C6`; no framing finding) | 4 `UPHELD`, 2 `NARROWED` (`F1`, `F3`) | None (`C3` `horizon_gfa stats` `SURVIVED`) | **Did not re-run** `horizon_gfa extract` (`chr22`, `chr6`), `horizon_lift`, `horizon_resolve`, `horizon_wga --graph`, `horizon_guide`, or full `chr6` on `HEAD` (`33bede8`) |

**Which lenses are thin:**
- **`L7` (`reproduction`) is the thinnest lens on execution coverage.** While `L7` re-ran `cargo test`/`clippy`/`build` (`C1`), a 20 Mbp `chr6` q-arm slice (`C2`), `horizon_gfa stats` (`C3`), and the sibling Python prototype `../pg/prototypes/egraph_synteny_resolver.py` (`C6`), it left **unrun** on `HEAD` (`33bede8`): both `horizon_gfa extract` jobs and `horizon_lift --bridge-bp 10` (`C3`), `horizon_wga --graph` (`C4`), full `GRCh38.chr6` × `CHM13.chr6` (`C2`), and `horizon_guide --progressive --hal` (`C6`). As `adjudication-L7.md:34-35` noted, `L7`'s `SURVIVED` verdict on `C3` applied only to `horizon_gfa stats` while `horizon_gfa extract` remained unexecuted.
- **`L2` (`sampling`) and `L3` (`statistics`) marked `C1` `UNCHALLENGED`** (`attack-L2-sampling.md:101`, `attack-L3-statistics.md:11-14`), which is appropriate for an architectural sentence, though `C1` also asserts `"59 active Rust unit tests (45 library + 14 across all 10 CLI binaries...)"` without any lens auditing per-binary test depth (carried item `R2-RESOLVE-LIFT`).

---

## 2. Claims Nobody Could Attack (Underspecified Claim Design)

- **`C1`** (`round.json:15-19`, `ATTACKED-horizon_paper.md:18`): `"Horizon is a unified pure-Rust suite of nine composable executables sharing one information-theoretic core..."`
  - **Why it resisted clean attack:** `"Unified"` and `"composable"` have no operational definition in `C1` or the manuscript. `L2` and `L3` had to mark `C1` `UNCHALLENGED` (`attack-L2-sampling.md:101`, `attack-L3-statistics.md:11-14`), while `L1`, `L4`, `L5`, `L6`, and `L7` could only attack the literal executable count (`9` vs. `10`), library imports, and `stdin` (`-`) support. Without a declared composability contract (e.g., `stdin`/`stdout` streaming vs. shared CLI dispatcher vs. shared core algorithms), attackers had to supply their own criterion.

---

## 3. Framing Defects in `round.json` (`C1`–`C6`)

1. **Overloaded compound claims + single-clause `source.quote` anchoring:**
   - `C3` (`round.json:40-43`) bundles `horizon_gfa stats`, `horizon_gfa extract` (`chr22`, `chr6`), and `horizon_lift --bridge-bp 10`, but `source.quote` quotes only `horizon_gfa`. All 7 attackers stopped after `horizon_gfa` and ignored `horizon_lift`.
   - `C4` (`round.json:53-56`) bundles `horizon_clean`, `horizon_resolve`, and `horizon_wga --graph`, but `source.quote` quotes only `horizon_clean` (`ATTACKED-horizon_paper.md:49`). All 7 attackers attacked `horizon_clean` and ignored `horizon_resolve` and `horizon_wga --graph`.
   - `C5` (`round.json:66-69`) bundles *P. falciparum* × *P. vivax*, *Cannabis sativa*, Table 1 (9 external suites), `horizon_merge` (406 mammals in `3.84 s`), and *P. malariae* × *P. brasilianum* (`50.2×`, `17.12 min` vs. `14.3 h`), but `source.quote` quotes only *Pf* × *Pv* and *Cannabis* (`ATTACKED-horizon_paper.md:65`). All 7 attackers focused on *Pf* × *Pv* and *Cannabis* (and `F7` framing caught `horizon_merge`/Table 1), while **zero** attackers examined the `50.2×` *P. malariae* speedup.
   - `C6` (`round.json:76-79`) bundles the `30` randomized e-graph trials and `horizon_guide --progressive --hal` (`194,878` MAF blocks, `+105.0%` `mus_musculus` coverage, `libhdf5`-validated `.hal`), but `source.quote` quotes only the `30` trials (`ATTACKED-horizon_paper.md:99`). All 7 attackers refuted the `30` trials and ignored `horizon_guide --progressive --hal`.

2. **Abstract and Results claims omitted from `C1`–`C6`:**
   - The Abstract (`ATTACKED-horizon_paper.md:10`) and Results (`line 59`) highlight the 4-regime synthetic *Plasmodium* benchmark (`0.978` recall, `0.949` precision, `0.129` vs. `0.936` hypervariable loop, `0.000` vs. `0.872` under `HOXD70`), yet `C1`–`C6` omit it entirely.

---

## 4. Specific Unattacked Gaps (Verified Against `attack-*.md`, `adjudication-*.md`, and `docs/BENCHMARK_RESULTS.md`)

### Gap 1: Synthetic 4-regime *Plasmodium* benchmark (`0.978` recall, `0.949` precision, `0.000` vs. `0.872` `HOXD70`) in Abstract and Results
- **Artifact:** `ATTACKED-horizon_paper.md:10, 59`; absent from all `attack-L1`..`attack-L7` reports and absent from `src/`, `scripts/`, and `docs/BENCHMARK_RESULTS.md` (verified via `grep -rn "0\.978\|0\.936\|0\.872\|HOXD70"` returning only `docs/manual.md:31`).
- **Why it matters:** This is a headline Abstract claim (*"eliminates random AT-rich false alignments (`0.000` vs. `0.872` under `HOXD70`)"*) that, like `C5` and `C6`, has zero repository provenance or reproduction on `33bede8` and may be an unarchived Python prototype result.
- **What would close it:** Add a dedicated claim for `ATTACKED-horizon_paper.md:59` in the next round, or commit the Rust benchmark script/seeds and `docs/BENCHMARK_RESULTS.md` row (or remove the unbacked figures from the Abstract and Results).

### Gap 2: *P. malariae* × *P. brasilianum* `50.2×` speedup (`17.12 min` vs. `14.3 h` / `859.9 min`, `536,041` bubbles) and `Pf` × `Pr` (`68.12 s`)
- **Artifact:** `ATTACKED-horizon_paper.md:10, 69, 73` and `round.json:66` (`C5` summary); zero mentions in `attack-L1`..`attack-L7` or `docs/BENCHMARK_RESULTS.md` (verified via `grep -rn "50\.2\|malariae\|859\.9\|68\.12"` returning 0 matches).
- **Why it matters:** Featured in the Abstract, Figure 3B, and Results as Horizon's largest whole-genome speedup (`50.2×`), yet it has no provenance in `docs/BENCHMARK_RESULTS.md` or `scripts/`, no commit pin, and no Galaxy `Growler UCSC pair` configuration/hardware parity check.
- **What would close it:** Separate whole-genome *Plasmodium* timing claims from `C5`, record the exact commit/script/dataset provenance in `docs/BENCHMARK_RESULTS.md`, or remove/qualify the unbacked `50.2×` claim.

### Gap 3: `horizon_guide --progressive --hal` (`194,878` MAF blocks, `+105.0%` `mus_musculus` coverage, `272.03 s`, `libhdf5`-validated `.hal`) and `src/hal/` (`sars-cov-2.hal` validation)
- **Artifact:** `ATTACKED-horizon_paper.md:75, 99, 111` and `round.json:76` (`C6` second clause); `194,878`, `1,012,274`, `493,677`, `99,884`, and `libhdf5` validation were not examined in any `attack-L1`..`attack-L7` report (only `attack-L4-provenance.md:103` noted `272.03 s` in framing `F7`), and zero matches exist in `src/`, `scripts/`, or `docs/BENCHMARK_RESULTS.md`.
- **Why it matters:** All 7 attackers refuted `C6` on the `30` e-graph trials alone, leaving the progressive guide-tree builder (`horizon_guide`) and pure-Rust HDF5 v0 `.hal` writer (`src/hal/hdf5_write.rs`) unverified for correctness, provenance, and reproduction.
- **What would close it:** Split `horizon_guide --progressive --hal` and `src/hal/` validation into their own claim, log the 5-species × 40-contig `60_amniotes` run in `docs/BENCHMARK_RESULTS.md` and `scripts/`, and run `L1`/`L4`/`L7` checks on `33bede8`.

### Gap 4: `horizon_lift --bridge-bp 10` (`182.67 s`, `619.7 MB`, `472` exact / `5` partial / `15` split / `8` unmapped) in `C3`
- **Artifact:** `ATTACKED-horizon_paper.md:48, 101`, `round.json:40`, `docs/BENCHMARK_RESULTS.md:21, 58`; `horizon_lift` was only mentioned in `C1` import-list checks (`attack-L1`, `L4`, `L5`, `L7`) and never audited under `C3` in `attack-L1`..`attack-L7`.
- **Why it matters:** `docs/BENCHMARK_RESULTS.md:21` discloses that `horizon_lift` sampled only `+`-strand `GRCh38.chr22:16.0M-49.9M` windows (excluding `0-16.0M` p-arm/centromere) with *"no external coordinate truth oracle"*, yet `ATTACKED-horizon_paper.md:101` omits both caveats and `L7` did not re-run `horizon_lift` on `HEAD` (`33bede8`).
- **What would close it:** Audit `src/bin/horizon_lift.rs` under `L1`/`L2`/`L5`/`L7`, re-run the 500-window `chr22` liftover on `33bede8`, and qualify `manuscript/horizon_paper.md:101` with the `+`-strand q-arm scope and absence of an external truth oracle.

### Gap 5: `horizon_wga --graph` (`src/mdl/dag.rs`) and `horizon_resolve` (`src/bin/horizon_resolve.rs`) in `C4`, plus carried items `R2-HMM-XDROP` and `R2-RESOLVE-LIFT`
- **Artifact:** `round.json:53` (`C4`), `CARRIED.md:5-11` (`R2-HMM-XDROP`, `R2-RESOLVE-LIFT`), `ATTACKED-horizon_paper.md:50-51`; `attack-L1`..`attack-L7` on `C4` examined only `horizon_clean`, leaving `horizon_wga --graph` (`src/mdl/dag.rs`), `horizon_resolve` (`W`-line coordinate and `L`-link rewriting), `N`-base handling in `bocpd_gc_profile`, `mdl_xdrop_extend`, and per-binary test coverage unverified on `33bede8`.
- **Why it matters:** `R2-HMM-XDROP` and `R2-RESOLVE-LIFT` were carried from Round 2 specifically to verify these code fixes on `HEAD`; because `C4` was attacked only on `horizon_clean`, the carried code-level obligations remain unverified.
- **What would close it:** Dedicated `L1`/`L7` verification of `src/mdl/dag.rs`, `mdl_xdrop_extend_with_params`, `src/bocpd.rs` (`N`-base handling), and `src/bin/horizon_resolve.rs` (`W`/`L` rewriting) on `33bede8`, plus a per-binary test-coverage count for `C1`.

---

## Summary Table: Unattacked Gaps

| Gap | Where in Manuscript / Round | Why It Matters | What Would Close It |
| :--- | :--- | :--- | :--- |
| 1. Synthetic 4-regime *Plasmodium* (`0.978` / `0.949`, `0.000` vs `0.872` `HOXD70`) | `ATTACKED-horizon_paper.md:10, 59`; omitted from `C1`–`C6` | Headline Abstract claim has zero repo provenance or `33bede8` verification | Add dedicated claim, commit script/logs to `docs/BENCHMARK_RESULTS.md`, or remove |
| 2. *P. malariae* × *P. brasilianum* `50.2×` speedup (`17.12 min` vs `14.3 h`) | `ATTACKED-horizon_paper.md:10, 69, 73`; `C5` summary | Largest WG speedup in Abstract/Fig 3B has zero `BENCHMARK_RESULTS.md` or script provenance | Split into its own claim; log commit/script/baseline parity or remove |
| 3. `horizon_guide --progressive --hal` (`194,878` blocks, `+105.0%` coverage) & `src/hal/` | `ATTACKED-horizon_paper.md:75, 99, 111`; `C6` second clause | Unverified by all 7 lenses because `C6` attacks stopped at the `30` e-graph trials | Split `horizon_guide`/`.hal` into a separate claim; log run in `BENCHMARK_RESULTS.md` and run `L1`/`L4`/`L7` |
| 4. `horizon_lift --bridge-bp 10` (`182.67 s`, `472` exact) | `ATTACKED-horizon_paper.md:101`; `C3` third clause | `+`-strand q-arm bias and lack of truth oracle (`BENCHMARK_RESULTS.md:21`) omitted from paper; unrun in `L7` | Audit `horizon_lift.rs`, re-run on `33bede8`, and qualify scope in `horizon_paper.md:101` |
| 5. `horizon_wga --graph`, `horizon_resolve`, & carried `R2-HMM-XDROP`/`R2-RESOLVE-LIFT` | `ATTACKED-horizon_paper.md:50-51`; `C4`, `CARRIED.md` | Carried code fixes (`N`-bases, `dag.rs`, `xdrop`, `W`/`L` rewriting) went unverified on `33bede8` | Run targeted `L1`/`L7` code/test checks on `dag.rs`, `bocpd.rs`, `dp.rs`, `horizon_resolve.rs` |
