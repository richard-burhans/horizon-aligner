# Round 3 Completeness Critic Report

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `488b020ccc0a73998915e975bad8e25599b16716`)  
**Artifacts read before asserting any gap:**
- `reports/rounds/round-3/CRITIC-BRIEF.md`, `reports/rounds/round-3/round.json`, `reports/rounds/round-3/CARRIED.md`
- `reports/rounds/round-3/ATTACKED-STRATEGIC_PLAN.md`, `reports/rounds/round-3/ATTACKED-BENCHMARK_RESULTS.md`
- All 7 attack reports: `attack-L1-correctness.md`, `attack-L2-sampling.md`, `attack-L3-statistics.md`, `attack-L4-provenance.md`, `attack-L5-inference.md`, `attack-L6-prior-art.md`, `attack-L7-reproduction.md`
- All 7 adjudication reports: `adjudication-L1.md`, `adjudication-L2.md`, `adjudication-L3.md`, `adjudication-L4.md`, `adjudication-L5.md`, `adjudication-L6.md`, `adjudication-L7.md`
- Source verification: `src/bin/horizon_wga.rs`, `src/bin/horizon_resolve.rs`, `src/bin/horizon_clean.rs`, `src/gfa/mod.rs`

---

## 1. Lens Execution, Filing Counts, and UNCHALLENGED vs. UNRUN

All **7** declared lenses (`correctness`, `sampling`, `statistics`, `provenance`, `inference`, `prior-art`, `reproduction`) ran, filed, and were adjudicated. **Zero** lenses were `UNRUN`, `aborted`, or `dropped`.

| Lens | Attack File | Findings Filed | Adjudication | Claims Challenged | Claims `UNCHALLENGED` (Checked, No Defect Under Lens) |
| :--- | :--- | :---: | :--- | :--- | :--- |
| L1 (`correctness`) | `attack-L1-correctness.md` | 7 (`F1`–`F7`) | `adjudication-L1.md` | `C1`–`C5`, Framing | `C6` (`hprc-v1.0-mc-grch38.gfa.gz` not found in L1's `/tmp/` check; `cmd_stats` code clean) |
| L2 (`sampling`) | `attack-L2-sampling.md` | 7 (`F1`–`F7`) | `adjudication-L2.md` | `C2`–`C6`, Framing | None (`C1` `SURVIVED`) |
| L3 (`statistics`) | `attack-L3-statistics.md` | 7 (`F1`–`F7`) | `adjudication-L3.md` | `C2`–`C6`, Framing | `C1` (API/reader capability with no statistical content) |
| L4 (`provenance`) | `attack-L4-provenance.md` | 7 (`F1`–`F7`) | `adjudication-L4.md` | `C2`–`C6`, Framing | None (`C1` `SURVIVED`) |
| L5 (`inference`) | `attack-L5-inference.md` | 7 (`F1`–`F7`) | `adjudication-L5.md` | `C1`–`C6`, Framing | None |
| L6 (`prior-art`) | `attack-L6-prior-art.md` | 7 (`F1`–`F7`) | `adjudication-L6.md` | `C2`, `C4`, Framing | `C3`, `C5`, `C6` (factual counts/inventory with no novelty assertion; `C1` `SURVIVED`) |
| L7 (`reproduction`) | `attack-L7-reproduction.md` | 5 (`F1`–`F5`) | `adjudication-L7.md` | `C2`–`C5`, Framing | None (`C1` and `C6` `horizon_gfa stats` `SURVIVED`) |

**Thinness assessment:** No lens is thin in total filing volume (5–7 findings each, covering all 6 claims plus framing). However, **L7 (`reproduction`)** and **L1 (`correctness`)** were thin on `C6` and trailing sub-claims of `C2`/`C3`: L7 reproduced only `horizon_gfa stats` out of the 5 pangenome benchmarks bundled into `C6`, and no lens inspected the `horizon_wga --graph` (`estimate_dag_hmm_params`) or `horizon_resolve` (`W`/`L` rewriting) implementations.

---

## 2. Framing & Claim-List Structure Critique

1. **Compound-claim bundling masks unexamined sub-claims (`round.json` `C1`–`C4`, `C6`):**  
   Five of the six claims pack 3–6 distinct assertions into a single summary (`C1`, `C2`, `C3`, `C4`, `C6`), while `round.json` records a single verdict per claim. Attackers consistently focused on the headline clause (`stream_gfa`, `cactus241way`, `45` vs. `46` tests, `FlatSeedIndex`, `horizon_gfa stats`) and left trailing sub-claims unexamined—specifically `estimate_dag_hmm_params` in `C2`, `horizon_resolve` `W`-line/`L`-line updates and `horizon_clean` tag preservation in `C3`, and 4 of the 5 pangenome benchmarks in `C6`.
2. **Runner-authored claim selection (`runner_is_author: true`, `claim_list_author: "runner"`):**  
   As all 7 lenses noted, `C1`–`C6` omitted the 5 stale pre-`488b020` WGA rows (`R2-HMM-XDROP`), the complete absence of `minimap2`/`odgi`/`vg` baselines (`R2-RESOLVE-LIFT`), the unfinished full `chr6` alignment (`R2-INDEX-RAM`), and the no-op/synthetic pangenome rows (`horizon_clean`, `horizon_resolve`, `horizon_wga --graph`), despite commit `488b020` claiming to discharge all three carried items.

---

## 3. Substantive Gaps Missed by Every Lens

### Gap 1: `estimate_dag_hmm_params` (`C2`, `src/bin/horizon_wga.rs:367-437`) was never inspected by any lens
- **Coverage check:** `attack-L2-sampling.md:12` noted in passing that `horizon_wga --graph` was wired to `estimate_dag_hmm_params`, and no other attack or adjudication file (`L1`–`L7`) inspected `src/bin/horizon_wga.rs:367-437`.
- **Why it matters:** `C2` claims adaptive `HmmParams` are wired into `horizon_wga --graph (estimate_dag_hmm_params)`. Direct inspection of `src/bin/horizon_wga.rs:367-437` shows that:
  1. Unlike `estimate_hmm_params`, `estimate_dag_hmm_params` **only estimates `identity_prior`** (lines 432–436) and leaves `gap_open_bits` and `gap_extend_bits` at `HmmParams::default()`.
  2. It samples up to `total >= 2000` bases (lines 411–426) in sorted hash-key order (`keys.sort_unstable()`, line 382) from exact `MASK_W10` seed windows (`canonical_spaced_seeds_with_opts`, lines 370, 377), which still allocates `HashMap`s and biases `identity_prior` upward by conditioning on seed-matching windows.
  3. It has no unit test in `src/bin/horizon_wga.rs`.

### Gap 2: `horizon_resolve` `W`-line/`L`-line rewriting (`src/bin/horizon_resolve.rs:416-480`) and `horizon_clean` non-`S` tag preservation (`src/bin/horizon_clean.rs:276-330`) in `C3` were never inspected by any lens
- **Coverage check:** Zero attack or adjudication files (`L1`–`L7`) inspected `horizon_resolve.rs`'s `W`-line coordinate calculation or `L`-line handling, or `horizon_clean.rs`'s `L`/`P`/`W` tag handling.
- **Why it matters:**
  1. Carried item `R2-RESOLVE-LIFT` explicitly required `"Update W-line start/end coordinates and L links in horizon_resolve --collapse-dups"`. In `src/bin/horizon_resolve.rs:478`, `L` lines are **copied verbatim** (`writeln!(out, "{trimmed}")?`) without updating or adding links when `--collapse-dups` collapses segments in `P`/`W` lines.
  2. In `src/bin/horizon_clean.rs:276-330`, `S`-line tags are preserved (`s.tags`), but `L`-line, `P`-line, and `W`-line optional tags are dropped when rewriting.
  3. Because no lens checked these lines, `C3` was treated as only a test-count (`45` vs. `46`) and `horizon_clean` `*`-seq no-op issue, leaving an unfulfilled `R2-RESOLVE-LIFT` code requirement (`L`-link updates in `horizon_resolve`) unreported.

### Gap 3: 4 of 5 benchmarks bundled into `C6` (`horizon_gfa extract`, `horizon_lift`, `horizon_resolve`, `horizon_wga --graph`) went unreproduced by L7 and unverified by L1
- **Coverage check:** `attack-L7-reproduction.md` (`F4`) and `adjudication-L7.md` (`F4`) reproduced only `horizon_gfa stats` (`88.84s`, `97.1 MB`), and `attack-L1-correctness.md` (`F6`) inspected only `cmd_stats`. Neither L7 nor L1 re-executed or verified `horizon_gfa extract` (`456.1/1,255.1 MB`), `horizon_lift` (`256.60s`, `619.2 MB`), `horizon_resolve --collapse-dups` (`4.48s`, `24.0 MB`), or `horizon_wga --graph` (`96/96 bp`, `0.0031s`).
- **Why it matters:** Commit `488b020` modified `src/bin/horizon_gfa.rs` (68 lines), `src/bin/horizon_lift.rs` (8 lines), `src/bin/horizon_resolve.rs` (53 lines), and `src/bin/horizon_wga.rs` (109 lines). Marking `C6` as `SURVIVED` in L7 based solely on `horizon_gfa stats` leaves 4 of the 5 bundled benchmarks unverified on `HEAD` (`488b020`).

---

## 4. Closing Summary: Gaps, Impact, and Resolution

| Gap | Artifact | Why It Matters | What Would Close It |
| :--- | :--- | :--- | :--- |
| **`estimate_dag_hmm_params` uninspected** | `src/bin/horizon_wga.rs:367-437` (`C2`) | Leaves `gap_open_bits`/`gap_extend_bits` at defaults, biases `identity_prior` via exact `MASK_W10` seed windows (`<=2000` bp), and uses `HashMap`s—contradicting `C2`'s adaptive graph HMM framing. | Qualify `C2` (`docs/STRATEGIC_PLAN.md:11`, `docs/BENCHMARK_RESULTS.md:5`) to note `estimate_dag_hmm_params` only adapts `identity_prior` from `MASK_W10` seed hits (leaving gap penalties at defaults), and add a unit test. |
| **`horizon_resolve` `L`-link & `horizon_clean` tag handling uninspected** | `src/bin/horizon_resolve.rs:416-480`, `src/bin/horizon_clean.rs:276-330` (`C3`, `R2-RESOLVE-LIFT`) | `horizon_resolve --collapse-dups` copies `L` lines verbatim instead of updating `L` links as required by `R2-RESOLVE-LIFT`, and `horizon_clean` drops `L`/`P`/`W` optional tags. | Keep `R2-RESOLVE-LIFT` open until `horizon_resolve --collapse-dups` updates `L` links (or explicitly narrow `C3` to `W`-line `end` + `P`/`W` tag preservation). |
| **4 of 5 `C6` benchmarks unreproduced on `488b020`** | `docs/BENCHMARK_RESULTS.md:55-57, 64-65` (`C6`) | `488b020` edited `horizon_gfa.rs`, `horizon_lift.rs`, `horizon_resolve.rs`, and `horizon_wga.rs`, yet L7 only reproduced `horizon_gfa stats`, leaving `extract`, `lift`, `resolve`, and `--graph` unverified on `HEAD`. | Split `C6` into per-benchmark claims in future rounds and re-run `horizon_gfa extract`, `horizon_lift`, `horizon_resolve`, and `horizon_wga --graph` on `488b020` with saved logs. |
