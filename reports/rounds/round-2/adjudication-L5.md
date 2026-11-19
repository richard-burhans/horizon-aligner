# Round 2 Adjudication of L5 (Inference Lens)

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner`  
**Attack file:** `reports/rounds/round-2/attack-L5-inference.md` (sha256 prefix verified: `42a67fc30cdc`)  
**Prior history & method cards:** `scripts/round history "HmmParams" --kind adjudication` returned no prior closed-round adjudications. No `refs/methods/` method cards exist in this repository.

---

## F1 — `C2` (`docs/BENCHMARK_RESULTS.md:5`): Pair-HMM attribution and `estimate_hmm_params`

- **Diagnosis:** **UPHELD.** Re-verifying `src/` confirms that `horizon_gfa`, `horizon_lift`, `horizon_clean`, and `horizon_resolve` do not invoke the Pair-HMM at all, and `horizon_wga --graph` explicitly passes `&HmmParams::default()` (`src/bin/horizon_wga.rs:195,504`). In pairwise `horizon_wga`, `estimate_hmm_params` (`src/mdl/mod.rs:759-843`) samples `lcb.iter().take(64)` up to `4000` bases, clamps `identity_prior` to `[0.60, 0.72]`, and estimates `indel_events` from net anchor differences `dx != dy` across `lcb.windows(2).take(64)`.
- **Proposed correction:** **UPHELD.** Replacing the blanket `"All results..."` sentence with an explicit breakdown of which binaries use `estimate_hmm_params` (with its `[0.60, 0.72]` clamp), `HmmParams::default()`, or no Pair-HMM introduces no false claims. **The fix is free** (text-only edit to `docs/BENCHMARK_RESULTS.md:5`, requiring no re-measurement).
- **What survives:** Pairwise `horizon_wga` uses normalized Pair-HMM transitions (`src/mdl/dp.rs`) and heuristic anchor-derived `HmmParams` clamped to `[0.60, 0.72]`, while `horizon_wga --graph` uses `HmmParams::default()` and the GFA utility tools do not use the Pair-HMM.

## F2 — `C3` (`docs/STRATEGIC_PLAN.md:3`): Crate-wide test count vs. per-binary unit tests

- **Diagnosis:** **UPHELD.** Re-running `cargo test --quiet` confirms `35 passed; 0 failed; 1 ignored` in the library crate plus `1` test each in 5 binaries (`horizon_clean`, `horizon_gfa`, `horizon_lift`, `horizon_resolve`, `horizon_wga`), totaling `40` passing tests (`1` ignored), while `grep -c "#\[test\]" src/bin/*.rs` confirms `0` binary unit tests in `bench_chr1.rs`, `horizon_eval.rs`, `horizon_guide.rs`, and `horizon_merge.rs`.
- **Proposed correction:** **UPHELD.** Narrowing `docs/STRATEGIC_PLAN.md:3` to state that all 9 binary files exist in `src/bin/` and 5 have binary-level unit tests alongside the library suite (while 4 have 0 binary unit tests) is accurate. **The fix is free** (text-only edit).
- **What survives:** All 9 target binary source files exist in `src/bin/` and compile, with 5 binaries containing binary-level unit tests and 4 (`bench_chr1`, `horizon_eval`, `horizon_guide`, `horizon_merge`) having 0 binary unit tests.

## F3 — `C4` (`docs/BENCHMARK_RESULTS.md:21,57`): `--bridge-bp 10` inflation of `EXACT` liftover counts

- **Diagnosis:** **UPHELD.** In `src/bin/horizon_lift.rs:256-282,490,507-509`, coalescing adjacent pieces within `bridge_bp` overwrites `prev.end = piece.end` and `prev.src_end = piece.src_end` (or `prev.start = piece.start` on opposite relative direction), and `mapped_bp` sums `p.src_end - p.src_start`. Consequently, windows containing `<= 10` bp gaps/indels or SNV bubbles that coalesce into a single piece satisfy `pieces.len() == 1 && mapped_bp == total_bp` and are counted as `EXACT`, with no external ground-truth liftover or sequence-identity check reported.
- **Proposed correction:** **UPHELD.** Qualifying `469 exact` as `"469 single-interval projections (including windows coalesced across <=10 bp bubbles by --bridge-bp 10)"` without claiming unconditional base-exactness is accurate. **The fix is free** (text-only edit to `docs/BENCHMARK_RESULTS.md:21,57`).
- **What survives:** `horizon_lift --bridge-bp 10` produced 469 single-interval projections (including windows bridged across `<= 10` bp bubbles), 1 partial, 19 split, and 11 unmapped windows in 256.60s at 619.2 MB peak RSS, without external ground-truth validation.

## F4 — `C5` (`docs/BENCHMARK_RESULTS.md:23,59,63`): Zero-severed `horizon_clean` and unvalidated pericentromere yield

- **Diagnosis:** **UPHELD.** In `src/bin/horizon_clean.rs:228-229`, `should_sever` requires `sup == 0` in both disjuncts (`(cfg.prune_unsupported && sup == 0) || (sup == 0 && cp_bits >= cfg.min_cp_bits)`), so any walk-supported link (`sup > 0`) cannot be severed regardless of `cp_bits`. In `docs/BENCHMARK_RESULTS.md:63`, both `chr22` (`6,030` links) and `chr6` (`2,776` links) had `0` severed links, which line 23 omits by citing only `6,030 links kept`. Additionally, the `horizon_wga` pericentromere pair and bubble counts have no precision/recall or base-identity ground truth.
- **Proposed correction:** **UPHELD.** Explicitly stating `0 severed` (all links had `sup > 0`) in `docs/BENCHMARK_RESULTS.md:23` (and `24`) and labeling the `horizon_wga` pericentromere numbers as throughput/yield counts without accuracy ground truth is accurate. **The fix is free** (text-only edit).
- **What survives:** `horizon_wga` produced `5,894,191` aligned pairs and `17,789` Horizon Bubbles on the 5 Mbp `chr22` pericentromere in 23.45s (3,084.5 MB RSS) as unvalidated yield counts, and `horizon_clean` ran in 0.01s (7.8 MB RSS) while severing 0 of 6,030 walk-supported links.

## F5 — `C6` (`docs/BENCHMARK_RESULTS.md:64`): Synthetic `horizon_resolve` benchmark and assumed `--divergence 0.02`

- **Diagnosis:** **UPHELD.** In `src/bin/horizon_resolve.rs:68,201-208,301-342`, `pass1_index` records only segment lengths (`s.sequence.len()`), ignores `L`-lines (`|_| Ok(())`), defines anchors by `visits >= 2 || len >= min_anchor_bp`, and populates `MdlCostModel` using the fixed `--divergence` parameter (default `0.02`) without comparing sequences. Placing the `200 segs, 20 walks x 225 steps` result under the HPRC `12 GB` benchmark table without labeling it synthetic or analytical misrepresents the measurement.
- **Proposed correction:** **UPHELD.** Labeling `C6` in `docs/BENCHMARK_RESULTS.md:64` as a synthetic benchmark and clarifying that `308,760.9` bits saved is an analytical `MdlCostModel` evaluation under fixed `--divergence 0.02` (without sequence comparison or topological snarl extraction) is accurate. **The fix is free** (text-only edit).
- **What survives:** On a synthetic `200 segs, 20 walks x 225 steps` graph, `horizon_resolve --collapse-dups` collapsed `840` dups across `280` windows in 4.48s (24.0 MB RSS), saving `308,760.9` analytical `MdlCostModel` bits under fixed `--divergence 0.02`.

## F6 — `C7` (`docs/PANGENOME_SUITE.md:18`): `align_seq_to_dag_seeded` forward-after-backward reachability closure

- **Diagnosis:** **UPHELD.** In `src/mdl/dag.rs:503-518`, the backward pass from `seed_hits` marks all ancestors back to the DAG root(s) as `active`, and the subsequent forward pass marks all descendants of every `active` node as `active`. In any connected DAG sharing an upstream root, this marks the entire connected DAG downstream of the root as `active` rather than isolating local components. (Additionally, `if nodes.len() <= 8` at line 466 bypasses seeding on the 4-node diamond benchmark in `docs/BENCHMARK_RESULTS.md:65`.)
- **Proposed correction:** **UPHELD on the documentation removal; algorithm rewrite is non-free.** Removing the claim in `docs/PANGENOME_SUITE.md:18` that forward/backward reachability closure isolates relevant DAG components is accurate and **the doc fix is free**; implementing true `s-t`/interval reachability in `src/mdl/dag.rs` would be a separate code change.
- **What survives:** `align_seq_to_dag_seeded` computes `MASK_W10` + `MASK_RY20` spaced-seed hits on DAGs with `> 8` nodes, but its backward-then-forward reachability closure activates the entire downstream connected DAG from shared ancestors rather than isolating local components.

## F7 — `C8` (`docs/BENCHMARK_RESULTS.md:26`): Inconsistent `~0.29 GB/Mbp` rate on cancelled `chr6` run

- **Diagnosis:** **UPHELD.** Dividing `26.5 GB` RSS by `170.8 Mbp` (`GRCh38.chr6`) yields `0.155 GB/Mbp` (or `0.094 GB/Mbp` against `281.5 Mbp` combined), not `~0.29 GB/Mbp` (which matches the preceding `20 Mbp` q-arm run: `5,863.3 MB / 20 Mbp = 0.293 GB/Mbp`). Moreover, `26.5 GB` was observed mid-run before cancellation, with no heap profile isolating `ContigIndex`'s `HashMap<u64, Vec<usize>>` from `BandGrid` DP or `intersect_seed_maps`. (Note: the same `~0.29 GB/Mbp` figure also appears in `docs/STRATEGIC_PLAN.md:11`.)
- **Proposed correction:** **UPHELD.** Removing `~0.29 GB/Mbp` from the `chr6` cancelled-run parenthetical in `docs/BENCHMARK_RESULTS.md:26` (and `docs/STRATEGIC_PLAN.md:11`) and noting that `26.5 GB` was a mid-run RSS at cancellation without component-level heap profiling is accurate. **The fix is free** (text-only edit).
- **What survives:** Full `GRCh38.chr6` vs `CHM13.chr6` reached `26.5 GB` RSS and `47:33` CPU before cancellation on the 30 GiB host, without component-level heap profiling or a `~0.29 GB/Mbp` scaling rate.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `horizon_gfa`/`lift`/`clean`/`resolve` do not use Pair-HMM, `horizon_wga --graph` uses `HmmParams::default()`, and `estimate_hmm_params` clamps `identity_prior` to `[0.60, 0.72]`; text fix is free. |
| F2 | UPHELD | `cargo test` and `grep` confirm 5 binaries have 1 test each and 4 (`bench_chr1`, `horizon_eval`, `horizon_guide`, `horizon_merge`) have 0 binary unit tests; text fix is free. |
| F3 | UPHELD | `horizon_lift --bridge-bp 10` overwrites `prev.src_end` across `<= 10` bp gaps, counting bridged gaps in `mapped_bp == total_bp` (`EXACT`) without ground-truth validation; text fix is free. |
| F4 | UPHELD | `horizon_clean` requires `sup == 0` so it severed 0 of 6,030 links (omitted in `C5`), and `horizon_wga` pericentromere counts lack accuracy ground truth; text fix is free. |
| F5 | UPHELD | `horizon_resolve` ignores `L`-lines and sequence comparisons, computing `308,760.9` bits saved analytically under fixed `--divergence 0.02` on synthetic walks; text fix is free. |
| F6 | UPHELD | Backward-then-forward closure in `align_seq_to_dag_seeded` activates the entire connected DAG downstream of shared roots rather than isolating components; doc removal is free. |
| F7 | UPHELD | `26.5 GB / 170.8 Mbp` is `0.155 GB/Mbp` (not `~0.29 GB/Mbp`), and `26.5 GB` was a mid-run RSS at cancellation without heap profiling; text fix is free. |
