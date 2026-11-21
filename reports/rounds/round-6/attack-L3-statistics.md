# Round 6 Attack -- L3 (Statistics)

- **Target round:** Round 6 (`42700e1`)
- **Lens:** `statistics`
- **Method cards & literature check (Rules 2–3):** Per [`docs/round_rubric.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/round_rubric.md#L34-L37), this project keeps **no method cards** yet. A checkout at `../library` is unavailable in this environment, so literature references are checked against [`refs/catalog_index.json`](file:///home/rcb112/project/gemini/bench/horizon-aligner/refs/catalog_index.json) (pinned at [`refs/LIBRARY_VERSION`](file:///home/rcb112/project/gemini/bench/horizon-aligner/refs/LIBRARY_VERSION) `d4df9122eddd84694b2c3215e5478c46f4314044`) alone, without asserting paper contents from memory.
- **Prior rounds checked (Rule 4):** Checked `scripts/round history` for `statistics`, `258.88`, `horizon_filter`, `horizon_locus`, and `horizon_resolve`, and read [`reports/rounds/round-5/adjudication-L3.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-5/adjudication-L3.md). Prior rulings are disclosed in each relevant finding below.

---

## F1 -- C1: `horizon_filter` omits BOCPD and `horizon_resolve` uses a sequence-blind divergence heuristic

- **Claim:** `C1` ([`docs/PANGENOME_SUITE.md:3`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/PANGENOME_SUITE.md#L3))
- **Verdict:** **NARROWED**
- **Prior adjudication:** Round 3 (`r3 adjudication-L3.md:54 [F6]`) and Round 4 (`r4 adjudication-L3.md:11 [F1]`) ruled earlier `C1` inventory claims had no statistical content. Round 6's `C1` materially changes this by explicitly claiming statistical/MDL methods: `horizon_filter (streaming PAF SDUST+BOCPD+MDL filter & netter)` and `horizon_resolve (... + per-window divergence)`.
- **Concrete flaw:**
  1. `horizon_filter` does **not** run Bayesian Online Changepoint Detection (BOCPD). [`src/bin/horizon_filter.rs:4-6`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_filter.rs#L4-L6) imports only `GAP_EXTEND_BITS`, `GAP_OPEN_BITS`, `IDENTITY_PRIOR`, `base_probs`, and `joint_pair_bits` from `horizon_aligner::bocpd` and never calls `bocpd_gc_profile` or `bocpd_gc_profile_into`. Instead, `IndexedSeq::from_dna` ([`src/bin/horizon_filter.rs:38-60`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_filter.rs#L38-L60)) computes a single static contig-wide GC fraction `(gc_cnt / acgt).clamp(0.05, 0.95)`, and sequence-free mode (`score_from_counts_only`, [`src/bin/horizon_filter.rs:91-136`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_filter.rs#L91-L136)) fixes `base_probs(0.5)`.
  2. `horizon_resolve`'s `"per-window divergence"` is a sequence-blind topological step-fraction formula (`F5`), not a sequence-divergence estimator.
- **Correction:** Replace `"streaming PAF SDUST+BOCPD+MDL filter & netter"` in [`docs/PANGENOME_SUITE.md:3`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/PANGENOME_SUITE.md#L3) with `"streaming PAF SDUST + global-GC MDL filter & netter"` and replace `"per-window divergence"` with `"topological per-window divergence heuristic"`.

---

## F2 -- C2: Unweighted record-level truth hit rate (`1.6% -> 81.2%`) mixed with bp-level coverage (`99.7%`), sequence-free scoring, and `n = 1` timing

- **Claim:** `C2` ([`docs/BENCHMARK_RESULTS.md:85`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L85))
- **Verdict:** **NARROWED**
- **Prior adjudication:** `scripts/round history "horizon_filter"` confirms zero prior mentions across closed rounds.
- **Concrete flaw:**
  1. **Unit-of-analysis mismatch:** `"81.2% 2D HPRC truth hit rate (from 1.6%)"` is an unweighted record-level proportion (`1,004` kept records out of `56,857` raw `FastGA` records, where raw records are dominated by short repeat hits), whereas `"99.7% 1-to-1 cov"` is a base-pair-weighted recall metric (`110.08 Mbp`). Base-pair precision on the `110.22 Mbp` retained by `--net` is omitted.
  2. **Sequence-free `0.47s` run is not BOCPD/sequence-MDL:** Running `horizon_filter --net` without `--query`/`--target` executes `score_from_counts_only` ([`src/bin/horizon_filter.rs:123-130`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_filter.rs#L123-L130)) with fixed `base_probs(0.5)` (`[0.25; 4]`) and `IDENTITY_PRIOR = 0.65`, followed by greedy 50%-reciprocal-overlap netting (`net_top_k_records`, [`src/bin/horizon_filter.rs:253-301`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_filter.rs#L253-L301)).
  3. **Unreplicated `n = 1` timing and missing benchmark script:** `0.47s` (`0.49s` `.paf.zst`) reports `n = 1` wall-clock time with no dispersion, and no script in `scripts/benchmarks/` reproduces `horizon_filter` or the `1.6% -> 81.2%` `FastGA` evaluation.
- **Correction:** Qualify [`docs/BENCHMARK_RESULTS.md:85`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L85) (and [`docs/BENCHMARK_RESULTS.md:29`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L29), [`manuscript/horizon_paper.md:84`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L84)) to state that `1.6% -> 81.2%` is an unweighted *record-level* hit rate (`n = 1` run under sequence-free `base_probs(0.5)` count scoring + greedy `--net`), and report base-pair precision on the `110.22 Mbp` retained.

---

## F3 -- C3: `jc69_dist = 0.222528` is a selection-biased, uncertainty-free point estimate over post-abstention aligned bases with `n = 1` timings

- **Claim:** `C3` ([`README.md:126`](file:///home/rcb112/project/gemini/bench/horizon-aligner/README.md#L126))
- **Verdict:** **NARROWED**
- **Prior adjudication:** `scripts/round history "horizon_locus"` confirms zero prior mentions across closed rounds.
- **Concrete flaw:**
  1. **Selection / survivorship truncation:** While the arithmetic $-\frac{3}{4}\ln(1 - \frac{4}{3}\cdot 0.192552) = 0.222528$ matches `jc69_distance` ([`src/bin/horizon_locus.rs:67-74`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_locus.rs#L67-L74)), `p_dist` is computed exclusively over the `231,444` aligned bp that survived MDL/Horizon Bubble abstention (~60–66% of the `386.7 kb` × `349.4 kb` sequences). Conditioning on post-abstention bases truncates divergent intervals and biases `p_dist` and `jc69_dist` downward relative to true locus divergence.
  2. **Zero uncertainty quantification:** Reporting `p_dist = 0.192552` and `jc69_dist = 0.222528` to 6 decimal places omits the Jukes-Cantor sampling variance / block-bootstrap CI across MAF blocks.
  3. **Run-level double-counting risk and JC69 misspecification:** `compute_pair_metrics` ([`src/bin/horizon_locus.rs:102-138`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_locus.rs#L102-L138)) deduplicates intervals for `jaccard_cov` (`unique_covered_bp`) but sums `aligned_bp = aln.total_pairs` and `matches` over all `aligned_runs` without interval deduplication. JC69 also assumes uniform `[0.25; 4]` base frequencies and equal transition/transversion rates, contradicting Horizon's compositional-skew premise.
  4. **`n = 1` timings with no accuracy verification:** `0.86s` and `0.32s` are single-run timings with no variance, no external truth/concordance metric (for either avian 2-way or human 4-way), and no reproduction script in `scripts/benchmarks/`.
- **Correction:** Qualify [`README.md:126`](file:///home/rcb112/project/gemini/bench/horizon-aligner/README.md#L126) and [`docs/BENCHMARK_RESULTS.md:90-91`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L90-L91) that `p_dist`/`jc69_dist` are conditional on post-abstention aligned bases without SE/CI, and mark `0.86s` / `0.32s` as `n = 1` unverified runtime/RSS measurements.

---

## F4 -- C4: Pre-HEAD `chr6` timing (`n = 1..2`) omits `2.58x` many-to-many multiplicity and `38.22%` unconditional exact pair precision

- **Claim:** `C4` ([`manuscript/horizon_paper.md:78`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L78))
- **Verdict:** **NARROWED**
- **Prior adjudication:** `n = 1` chromosome/q-arm/run was upheld in `r4 adjudication-L2.md:15 [F1]`; the `8cd5db0` vs `2896229` split (`304.53 s`) was upheld in `r5 adjudication-L4.md:23 [F2]`; and the `38.22%` unconditional (`48.63%` conditioned) exact pair precision against HPRC `1-to-1` segments was upheld in `r5 adjudication-L3.md:20-28 [F2]` (`NARROWED` on `C2` there).
- **Concrete flaw:**
  1. `C4` isolates the `8cd5db0` sentence, which reports `282,193,564` aligned pairs covering `124,160,423` bp (`72.7%`) of `GRCh38.chr6` and `109,464,447` bp (`98.9%`) of `CHM13.chr6` q-arm without stating that `282.19M` pairs over `109.46 Mbp` is **`2.58x` many-to-many multiplicity** and achieves **`38.22%` unconditional (`48.63%` conditioned)** exact pair precision against HPRC `1-to-1` segments ([`manuscript/data/chr6_hprc_comparison.json:10-28`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/data/chr6_hprc_comparison.json#L10-L28)).
  2. Speedup ratios (`2.21x` vs `592.68 s` at `bb43dd9`; `3.97x` vs `1,062.45 s` at `1e47fe0`) divide `n = 1` historical baselines by the upper bound (`267.86 s`) of an `n = 2` range (`258.88–267.86 s`) with no mean/SD or replicate count, on a pre-HEAD commit (`8cd5db0`) whose runtime already regressed `+13.7%` to `304.53 s` at `2896229`.
- **Correction:** Qualify `C4` in [`manuscript/horizon_paper.md:78`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L78) with sample sizes (`n = 2` on `8cd5db0`, `n = 1` on `bb43dd9`/`1e47fe0`), the `+13.7%` `2896229` runtime (`304.53 s`), and the `2.58x` multiplicity / `38.22%` unconditional exact pair precision against HPRC `1-to-1` segments.

---

## F5 -- C5: `horizon_resolve` "per-window divergence" is a sequence-blind, orientation-dependent topological heuristic clamped to `[0.01, 0.25]`

- **Claim:** `C5` ([`ARCHITECTURE.md:28`](file:///home/rcb112/project/gemini/bench/horizon-aligner/ARCHITECTURE.md#L28))
- **Verdict:** **NARROWED**
- **Prior adjudication:** `scripts/round history "horizon_resolve"` shows prior rounds audited `horizon_resolve` existence, `W`-line `end` bounds, `L`-link rewriting, and `EG:Z:` tags, but never audited `estimate_window_divergence`.
- **Concrete flaw:** `estimate_window_divergence` ([`src/bin/horizon_resolve.rs:456-477`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_resolve.rs#L456-L477)) never reads DNA sequences or substitution counts (`PackedGraph` does not even store segment sequence strings). Instead it computes:
  $$\hat{d} = \text{clamp}\!\left(0.02 + 0.08 \cdot \frac{\text{non\_anchors}}{n} + 0.05 \cdot \frac{\text{rev\_steps}}{n},\; 0.01,\; 0.25\right)$$
  1. **Orientation-dependent:** Reverse-complementing a walk changes `rev_steps / n` from `0` to `1`, adding `+0.05` to "divergence" for the exact same sequence.
  2. **Biased clamping for `MdlCostModel::block_residual_bits`:** Clamping to `[0.01, 0.25]` prevents estimating intra-species human divergence (`~0.001` in HPRC), inflating `H_2(d) + d log_2 3` ([`src/egraph.rs:416-419`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/egraph.rs#L416-L419)) by `7.4x` (`0.0966` bits/bp at `d = 0.01` vs `0.0130` bits/bp at `d = 0.001`), while capping deep vertebrate divergence (`d > 0.25`) downward.
- **Correction:** Replace `"estimates per-window divergence (unless --divergence is passed)"` in [`ARCHITECTURE.md:28`](file:///home/rcb112/project/gemini/bench/horizon-aligner/ARCHITECTURE.md#L28) (and [`docs/PANGENOME_SUITE.md:3,25`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/PANGENOME_SUITE.md#L3), [`docs/manual.md:144`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/manual.md#L144)) with `"assigns a topological per-window divergence proxy in [0.01, 0.25] from non-anchor and reverse-step fractions (unless --divergence is passed)"`.

---

## F6 -- Claim-list steering: `C1`–`C5` omit uncalibrated CRC (`n_cal = 15`), `KARLIN_MDL_COEFF = 0.55`, and the `18%`-GC null sweep

- **Claim:** Scope critique across `C1`–`C5`
- **Prior adjudication:** Directly extends `r4 adjudication-L3.md:62 [F7]` and `r5 adjudication-L3.md:60 [F6]`.
- **Concrete flaw:** As in Round 5, the runner-authored claim list (`C1`–`C5`) omits the project's core statistical calibration claims in [`manuscript/horizon_paper.md:10, 40, 61`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L10): uncalibrated Conformal Risk Control (`n_cal = 15`, `α = 0.07`, `λ̂ = 10.0` bits; cf. `angelopoulos_2022_conformal` in `refs/catalog_index.json`), fixed `KARLIN_MDL_COEFF = 0.55`, and the `18%`-GC random-DNA false-positive result (`0` on `seed=999`, `80–448` on `4/5` other seeds, `288`/`149` on single `30%`/`45%`-GC draws).

---

## Verdict Summary

| Claim | Verdict | Summary |
|---|---|---|
| C1 | NARROWED | `horizon_filter` never calls BOCPD (`bocpd_gc_profile`), using static contig GC or `0.5`; `horizon_resolve` uses a sequence-blind topological heuristic |
| C2 | NARROWED | Mixes unweighted record-level `81.2%` (`1.6%` raw) with bp-level `99.7%` cov, omits bp precision, uses fixed `0.5` GC counts in `--net`, `n = 1` timing |
| C3 | NARROWED | `jc69_dist = 0.222528` is conditioned on post-abstention aligned bases without SE/CI, un-deduplicated runs, JC69 misspecification, `n = 1` unverified |
| C4 | NARROWED | Pre-HEAD `8cd5db0` (`n = 2`) vs `n = 1` baselines (`+13.7%` regression at `2896229`), omits `2.58x` multiplicity and `38.22%` unconditional exact pair precision |
| C5 | NARROWED | `estimate_window_divergence` never reads sequences, depends on strand orientation (`inv_frac`), and clamps to `[0.01, 0.25]` (`7.4x` residual bias at `d=0.001`) |
