# Round 2 Adversarial Attack — Inference Lens

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`docs/BENCHMARK_RESULTS.md`, `docs/STRATEGIC_PLAN.md`, `docs/PANGENOME_SUITE.md`)  
**Lens:** `inference`  
**Context checks (Rules 2–4):**
- No method cards or held literature PDFs/summaries directory are declared in `round.toml` or the brief appendix; checked workspace (`manuscript/references.bib` and `manuscript/*.md` exist, no method cards exist).
- Checked `reports/adversarial_log.md` and ran `scripts/round history "HmmParams"`: Round 1 left `R1-IMPL`, `R1-GFA`, and `R1-HMM` open; the inference defects below target post-Round-1 claims (`C1`–`C8`) and do not repeat overturned Round 1 points.

---

## F1 — `C2` falsely attributes all benchmark results to divergence-adaptive `HmmParams`, and `estimate_hmm_params` infers divergence from conserved anchor seeds clamped to `[0.60, 0.72]`

- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:5`)
- **Verdict:** **REFUTED**
- **Concrete flaw:**
  1. `"All results were generated out-of-the-box using the divergence-adaptive Pair-HMM (HmmParams...)"` overgeneralizes across `docs/BENCHMARK_RESULTS.md`: `horizon_gfa`, `horizon_lift`, `horizon_resolve`, and `horizon_clean` do not use Pair-HMM (`src/mdl/dp.rs` or `src/mdl/dag.rs`) at all, and `horizon_wga --graph` explicitly passes `&HmmParams::default()` rather than adaptive `HmmParams`.
  2. Even in pairwise `horizon_wga`, `estimate_hmm_params` infers `identity_prior` from LCB anchor k-mers (`lcb.iter().take(64)`), which are ascertainment-biased toward conserved seed hits, and then clamps the result to `[0.60, 0.72]` (`4000` bases max). It also infers `gap_open_bits` and `gap_extend_bits` from net coordinate differences `|dx - dy|` between adjacent LCB anchors, which misses balanced indels (`dx == dy`) and collapses multiple indels into a single event.
- **Evidence:**
  - `src/bin/horizon_wga.rs:195`: `align_seq_to_dag_seeded(&dag_nodes, qseq, &HmmParams::default(), 10.0)`
  - `src/mdl/mod.rs:759-843`: `estimate_hmm_params` samples `lcb.iter().take(64)` up to `4000` bases (lines 775, 790), clamps `identity_prior` to `[0.60, 0.72]` (line 818), and counts `indel_events` only when `dx != dy` across `lcb.windows(2).take(64)` (lines 794–804).
- **Correction:** Replace `"All results were generated out-of-the-box using the divergence-adaptive Pair-HMM..."` in `docs/BENCHMARK_RESULTS.md:5` with: `"Pairwise WGA results were generated with Pair-HMM transitions normalized in src/mdl/dp.rs and heuristic anchor-derived HmmParams (identity_prior clamped to [0.60, 0.72] in src/mdl/mod.rs); horizon_wga --graph uses HmmParams::default(), and GFA utility benchmarks (horizon_gfa, horizon_lift, horizon_clean, horizon_resolve) do not use the Pair-HMM."`

---

## F2 — `C3` infers from a crate-wide test count that all listed binaries are unit-tested

- **Claim:** `C3` (`docs/STRATEGIC_PLAN.md:3`)
- **Verdict:** **NARROWED**
- **Concrete flaw:** Aggregating the crate-wide test count (`39` in the claim; `40 passed, 1 ignored` on `cargo test --quiet`) to infer that `"All target pangenome binaries ... are now implemented in src/bin/ and unit-tested"` is a composition fallacy: 4 of the 9 binary source files in `src/bin/` have zero unit tests.
- **Evidence:**
  - `cargo test --quiet` reports `35 passed; 0 failed; 1 ignored` in the library crate and `1` passing test each in 5 binaries (`horizon_clean`, `horizon_gfa`, `horizon_lift`, `horizon_resolve`, `horizon_wga`), with `0` tests in the remaining 4 binaries.
  - `grep -c "#\[test\]" src/bin/*.rs`: `src/bin/bench_chr1.rs:0`, `src/bin/horizon_eval.rs:0`, `src/bin/horizon_guide.rs:0`, `src/bin/horizon_merge.rs:0`.
- **Correction:** Narrow `docs/STRATEGIC_PLAN.md:3` to state that all 9 binary files exist in `src/bin/` and 5 of them (`horizon_clean`, `horizon_gfa`, `horizon_lift`, `horizon_resolve`, `horizon_wga`) have binary-level unit tests alongside the library test suite, while `horizon_merge`, `horizon_guide`, `horizon_eval`, and `bench_chr1` have 0 binary unit tests.

---

## F3 — `C4` infers `469 exact` liftover windows under `--bridge-bp 10` without accounting for bridged gap inflation or ground-truth validation

- **Claim:** `C4` (`docs/BENCHMARK_RESULTS.md:21,57`)
- **Verdict:** **NARROWED**
- **Concrete flaw:** `horizon_lift --bridge-bp 10` merges adjacent lifted pieces separated by up to `10` bp on both source and target walks by overwriting `prev.src_end = piece.src_end`. Because `mapped_bp` sums `p.src_end - p.src_start`, up to 10 bp of unshared/indel sequence per gap are counted as mapped bases, allowing intervals with `< 10` bp indels or SNV bubbles to be classified as `EXACT` (`mapped_bp == total_bp`). Furthermore, no external ground-truth liftover (e.g., `halLiftover` or sequence identity check) is reported.
- **Evidence:**
  - `src/bin/horizon_lift.rs:256-269`: when `piece.start - prev.end <= bridge_bp` and `piece.src_start - prev.src_end <= bridge_bp`, `prev.end = piece.end;` and `prev.src_end = piece.src_end;`.
  - `src/bin/horizon_lift.rs:490,507-509`: `mapped_bp` sums `p.src_end - p.src_start`, and `pieces.len() == 1 && mapped_bp == total_bp` increments `summary.lifted_exact`.
- **Correction:** Qualify `469 exact` in `docs/BENCHMARK_RESULTS.md:21,57` as `"469 single-interval projections (including windows coalesced across <=10 bp bubbles by --bridge-bp 10)"`, rather than unconditional exact liftover.

---

## F4 — `C5` cites a zero-severed `horizon_clean` run on `chr22` pericentromere as cleaning evidence and lacks accuracy validation for `horizon_wga`

- **Claim:** `C5` (`docs/BENCHMARK_RESULTS.md:23,59,63`)
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. On `chr22` pericentromere (`6,030` links) and `chr6` pericentromere (`2,776` links), `horizon_clean` severed **0** links (`6,030 / 6,030` and `2,776 / 2,776` kept). Because `horizon_clean` only severs links with `sup == 0` (zero walk support), any link traversed by a haplotype walk is immune regardless of compositional changepoint score. Citing a `0`-severed no-op as evidence of pericentromere graph cleaning is vacuous.
  2. Reporting `5,894,191` aligned pairs and `17,789` Horizon Bubbles for `horizon_wga --resolve-dups` on `chr22` pericentromere without precision/recall or base-identity validation does not establish alignment accuracy.
- **Evidence:**
  - `src/bin/horizon_clean.rs:228-229`: `let should_sever = (cfg.prune_unsupported && sup == 0) || (sup == 0 && cp_bits >= cfg.min_cp_bits);`
  - `docs/BENCHMARK_RESULTS.md:63`: `6,030 / 2,776 kept, 0 severed`.
- **Correction:** State explicitly in `docs/BENCHMARK_RESULTS.md:23,63` that `horizon_clean` severed `0` links on both pericentromere subgraphs (all links had walk support `sup > 0`), and label the `horizon_wga` pericentromere metrics as throughput/yield counts without accuracy ground truth.

---

## F5 — `C6` reports `308,760.9 bits saved` from an assumed `--divergence 0.02` on synthetic walks without sequence comparison or topological snarl decomposition

- **Claim:** `C6` (`docs/BENCHMARK_RESULTS.md:64`)
- **Verdict:** **REFUTED**
- **Concrete flaw:** `horizon_resolve` never compares segment sequences—`pass1_index` records only segment lengths (`s.sequence.len()`) and ignores `L`-lines (`|_| Ok(())`). It assigns a fixed, user-supplied `--divergence` constant (default `0.02`) to every segment in `MdlCostModel`, and defines "anchors" as `visits >= 2 || len >= min_anchor_bp` rather than computing graph snarls/bubbles. Thus `308,760.9 bits saved` on the synthetic `200 segs, 20 walks x 225 steps` dataset is an analytical cost-model evaluation under an assumed `0.02` divergence, not an empirical sequence-compression or snarl-resolution measurement.
- **Evidence:**
  - `src/bin/horizon_resolve.rs:68,201-207`: `divergence` defaults to `0.02` and populates `model.block_stats` for every segment without sequence alignment.
  - `src/bin/horizon_resolve.rs:301-342`: `pass1_index` ignores `L` lines (`|_| Ok(())` at line 312) and sets `anchors` by `visits >= 2 || len >= min_anchor_bp` (line 338).
- **Correction:** Label `C6` in `docs/BENCHMARK_RESULTS.md:64` as a **synthetic** benchmark (`200 segs, 20 walks x 225 steps`) and clarify that `308,760.9` bits saved is computed analytically under fixed `--divergence 0.02` without sequence comparison or topological snarl extraction.

---

## F6 — `C7` claims reachability closure isolates relevant DAG components, but forward-after-backward closure activates the entire downstream connected DAG

- **Claim:** `C7` (`docs/PANGENOME_SUITE.md:18`)
- **Verdict:** **REFUTED**
- **Concrete flaw:** `align_seq_to_dag_seeded` first runs a backward pass from `seed_hits` marking all ancestors (back to the DAG root(s)) as `active`, and then runs a forward pass marking all descendants of `active` nodes as `active`. In any connected DAG sharing a root or upstream anchor, activating all ancestors and then all descendants of those ancestors marks the **entire connected DAG downstream of the root** as `active`, failing to isolate relevant local components. (Additionally, the 4-node diamond benchmark in `docs/BENCHMARK_RESULTS.md:65` bypasses seeding via `if nodes.len() <= 8`.)
- **Evidence:**
  - `src/mdl/dag.rs:503-518`: backward loop (`for u in (0..nodes.len()).rev()`) activates all predecessors of `seed_hits`; forward loop (`for u in 0..nodes.len()`) activates all successors of all `active` nodes.
  - `src/mdl/dag.rs:466`: `if nodes.len() <= 8 { return align_seq_to_dag_mdl(...); }`.
- **Correction:** Remove the claim in `docs/PANGENOME_SUITE.md:18` that forward/backward reachability closure isolates relevant DAG components (or fix the algorithm to restrict to nodes that are *both* ancestors of downstream hits and descendants of upstream hits, i.e., interval/s-t reachability).

---

## F7 — `C8` attaches an inconsistent `~0.29 GB/Mbp` scaling rate to a cancelled `26.5 GB` / `170.8 Mbp` run without heap profiling

- **Claim:** `C8` (`docs/BENCHMARK_RESULTS.md:26`)
- **Verdict:** **NARROWED**
- **Concrete flaw:** `26.5 GB` RSS on `170.8 Mbp` `GRCh38.chr6` vs `110.7 Mbp` `CHM13.chr6` corresponds to `0.155 GB/Mbp` (relative to `170.8 Mbp`) or `0.094 GB/Mbp` (relative to `281.5 Mbp` combined), not `~0.29 GB/Mbp` (which matches the preceding `20 Mbp` q-arm run: `5,863.3 MB / 20 Mbp = 0.293 GB/Mbp`). Moreover, `26.5 GB` was observed *before cancellation* (not a completed peak RSS), and no heap profile is cited to isolate `ContigIndex`'s `HashMap<u64, Vec<usize>>` from `BandGrid` DP or `intersect_seed_maps`.
- **Evidence:**
  - `docs/BENCHMARK_RESULTS.md:26`: `5,863.3 MB` on `20 Mbp` (`0.293 GB/Mbp`) vs `26.5 GB` before cancellation on `170.8 Mbp` (`0.155 GB/Mbp`).
  - `src/seeding/mod.rs:329-334`: `ContigIndex` stores `seeds_fwd` and `seeds_rc` as `HashMap<u64, Vec<usize>>` (qualitatively large, but unprofiled in the cancelled run).
- **Correction:** Remove `~0.29 GB/Mbp` from the `chr6` cancelled-run parenthetical in `docs/BENCHMARK_RESULTS.md:26` and note that `26.5 GB` was a mid-run RSS at cancellation without component-level heap profiling.

---

## Verdict Table

| Claim | Verdict | Summary |
| :--- | :--- | :--- |
| **C1** | `SURVIVED` | Descriptive claim about `stream_gfa` verified in `src/gfa/mod.rs:98-120,250-318`. |
| **C2** | `REFUTED` | Non-HMM tools and `--graph` (`HmmParams::default()`) excluded; `estimate_hmm_params` uses anchor-biased seeds clamped to `[0.60, 0.72]`. |
| **C3** | `NARROWED` | 4 of 9 binaries (`bench_chr1`, `horizon_eval`, `horizon_guide`, `horizon_merge`) have 0 unit tests. |
| **C4** | `NARROWED` | `--bridge-bp 10` absorbs <=10 bp unmapped gaps into `mapped_bp == total_bp` (`EXACT`) with no ground-truth check. |
| **C5** | `NARROWED` | `horizon_clean` severed 0 of 6,030 links (no-op); `horizon_wga` pericentromere counts lack accuracy ground truth. |
| **C6** | `REFUTED` | Synthetic benchmark using fixed `--divergence 0.02` without sequence comparison or topological snarl extraction. |
| **C7** | `REFUTED` | Forward-after-backward closure activates the entire downstream connected DAG rather than isolating components. |
| **C8** | `NARROWED` | `~0.29 GB/Mbp` conflicts with `26.5 GB / 170.8 Mbp` (`0.155 GB/Mbp`), and `26.5 GB` is a cancelled mid-run RSS without heap profiling. |
