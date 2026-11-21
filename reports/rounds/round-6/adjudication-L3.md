# Round 6 Adjudication -- Statistics Lens (L3)

- **Target attack:** [`reports/rounds/round-6/attack-L3-statistics.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-6/attack-L3-statistics.md) (sealed SHA-256 `2f4570a068b71f4b2899010c5d5920a8744987697bb5d0fc1c71c262d057a68c`, verified via `sha256sum`).
- **Method cards & prior rounds (Rules 5–6):** This repository has no `refs/methods/` directory ([`docs/round_rubric.md:34-37`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/round_rubric.md#L34-L37)), so no method card governs these disputes. Checked `scripts/round history` for `statistics`, `horizon_filter`, `horizon_locus`, and `horizon_resolve`; prior rulings (`r3 adjudication-L3.md [F6]`, `r4 adjudication-L2.md [F1]`, `r4 adjudication-L3.md [F1, F7]`, `r5 adjudication-L3.md [F2, F6]`, `r5 adjudication-L4.md [F2]`) are properly disclosed by the attack.

---

## F1 -- C1: `horizon_filter` omits BOCPD and `horizon_resolve` uses a sequence-blind divergence heuristic

- **Ruling:** **UPHELD** (attacker's `NARROWED` verdict on `C1` stands).
- **What survives:** `C1` ([`docs/PANGENOME_SUITE.md:3`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/PANGENOME_SUITE.md#L3)) claims `horizon_filter` is a `streaming PAF SDUST+BOCPD+MDL filter & netter` and `horizon_resolve` performs `per-window divergence` estimation, neither of which matches the implementation.
- **Verification:**
  1. [`src/bin/horizon_filter.rs:4-6`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_filter.rs#L4-L6) imports only `GAP_EXTEND_BITS`, `GAP_OPEN_BITS`, `IDENTITY_PRIOR`, `base_probs`, and `joint_pair_bits` from `horizon_aligner::bocpd` and never calls `bocpd_gc_profile` or `bocpd_gc_profile_into`. `IndexedSeq::from_dna` ([`src/bin/horizon_filter.rs:38-60`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_filter.rs#L38-L60)) computes a single static contig-wide GC fraction `(gc_cnt / acgt).clamp(0.05, 0.95)`, and sequence-free mode (`score_from_counts_only`, [`src/bin/horizon_filter.rs:91-136`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_filter.rs#L91-L136)) fixes `base_probs(0.5)`.
  2. [`src/bin/horizon_resolve.rs:456-477`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_resolve.rs#L456-L477) (`estimate_window_divergence`) is a sequence-blind topological step-fraction heuristic (`F5`).
- **Repair:** Replacing `"streaming PAF SDUST+BOCPD+MDL filter & netter"` in [`docs/PANGENOME_SUITE.md:3`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/PANGENOME_SUITE.md#L3) (and [`docs/PANGENOME_SUITE.md:28`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/PANGENOME_SUITE.md#L28)) with `"streaming PAF SDUST + global-GC (or 0.5-GC count-only) MDL filter & netter"` and `"per-window divergence"` with `"topological per-window divergence proxy"` is endorsed (**the fix is free**).

---

## F2 -- C2: Unweighted record-level truth hit rate (`1.6% -> 81.2%`) mixed with bp-level coverage (`99.7%`), sequence-free scoring, and `n = 1` timing

- **Ruling:** **UPHELD** on the diagnosis (attacker's `NARROWED` verdict on `C2` stands); **NARROWED** on the repair (text qualification is free; computing bp-level precision is non-free).
- **What survives:** `C2` ([`docs/BENCHMARK_RESULTS.md:85`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L85)) mixes an unweighted record-level hit rate (`1.6% -> 81.2%` across `56,857` raw / `1,004` kept `FastGA` records) with a base-pair-weighted recall metric (`99.7%` of `110.08 Mbp` 1-to-1 coverage), omits base-pair precision on the `110.22 Mbp` retained, reports `n = 1` wall-clock timings (`0.47s` / `0.49s`), and lacks a reproduction script in `scripts/benchmarks/`.
- **Verification:**
  1. [`src/bin/horizon_filter.rs:91-136`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_filter.rs#L91-L136) confirms sequence-free `--net` runs `score_from_counts_only` with fixed `base_probs(0.5)` (`[0.25; 4]`) and `IDENTITY_PRIOR = 0.65`, followed by greedy `--net` (`net_top_k_records`, [`src/bin/horizon_filter.rs:253-301`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_filter.rs#L253-L301), filtering candidates with `>= 50%` of the candidate's own length overlapping already-kept intervals).
  2. [`scripts/benchmarks/`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks) contains [`score_paf_vs_hprc.py`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/score_paf_vs_hprc.py) for scoring PAFs against HPRC `1-to-1` segments, but no script reproducing `horizon_filter` or the `1.6% -> 81.2%` `FastGA` evaluation.
  3. Note that [`docs/BENCHMARK_RESULTS.md:85`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L85) already labels the row `(seq-free)` in column 1, whereas [`README.md:124`](file:///home/rcb112/project/gemini/bench/horizon-aligner/README.md#L124) and [`manuscript/horizon_paper.md:84`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L84) omit both the record-level unit and `n = 1` qualification (and the attacker's citation of `docs/BENCHMARK_RESULTS.md:29` refers to `horizon_wga`, not `horizon_filter`).
- **Repair:** Qualifying [`docs/BENCHMARK_RESULTS.md:85`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L85), [`README.md:124`](file:///home/rcb112/project/gemini/bench/horizon-aligner/README.md#L124), and [`manuscript/horizon_paper.md:84`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L84) to state that `1.6% -> 81.2%` is an unweighted *record-level* hit rate (`n = 1`, sequence-free `base_probs(0.5)` count scoring + greedy `--net`) is endorsed (**the fix is free**). Computing and reporting base-pair precision on the `110.22 Mbp` retained is endorsed as **non-free** (requires re-scoring the `FastGA` PAF).

---

## F3 -- C3: `jc69_dist = 0.222528` is a selection-biased, uncertainty-free point estimate over post-abstention aligned bases with `n = 1` timings

- **Ruling:** **UPHELD** (attacker's `NARROWED` verdict on `C3` stands).
- **What survives:** `p_dist = 0.192552` and `jc69_dist = 0.222528` in [`README.md:126`](file:///home/rcb112/project/gemini/bench/horizon-aligner/README.md#L126) are conditional on the `231,444` post-abstention aligned bases (~60–66% of the avian sequences), omit any SE/CI, sum `aligned_bp = aln.total_pairs` and `matches` across `aligned_runs` without interval deduplication, assume JC69 uniform `[0.25; 4]` base frequencies, and accompany `n = 1` unverified timings (`0.86s`, `0.32s`) with no reproduction script in `scripts/benchmarks/`.
- **Verification:**
  1. Verified `-0.75 * ln(1 - (4/3) * 0.192552) = 0.222528` in [`src/bin/horizon_locus.rs:67-74`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_locus.rs#L67-L74).
  2. Verified [`src/bin/horizon_locus.rs:102-145`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_locus.rs#L102-L145) (`compute_pair_metrics`) deduplicates intervals via `unique_covered_bp` for `jaccard_cov` while summing `aligned_bp = aln.total_pairs` and `matches` over all `aligned_runs` without interval deduplication.
  3. Verified no script in `scripts/benchmarks/` reproduces `horizon_locus`.
- **Repair:** Qualifying [`README.md:126`](file:///home/rcb112/project/gemini/bench/horizon-aligner/README.md#L126) and [`docs/BENCHMARK_RESULTS.md:90-91`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L90-L91) that `p_dist`/`jc69_dist` are conditional on post-abstention aligned bases without SE/CI (and summed over un-deduplicated runs) and marking `0.86s` / `0.32s` as `n = 1` unverified runtime/RSS measurements is endorsed (**the fix is free**).

---

## F4 -- C4: Pre-HEAD `chr6` timing (`n = 1..2`) omits `2.58x` many-to-many multiplicity and `38.22%` unconditional exact pair precision

- **Ruling:** **NARROWED** (**UPHELD** on the `n = 2` vs `n = 1` speedup-ratio reporting and missing mean/SD; **OVERTURNED** on the claim that [`manuscript/horizon_paper.md:78`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L78) omits `304.53 s`, `2.58×` multiplicity, and `38.22%` unconditional exact pair precision).
- **What survives:** In [`manuscript/horizon_paper.md:78`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L78), the speedup ratios (`2.21×` vs `592.68 s` at `bb43dd9`; `3.97×` vs `1,062.45 s` at `1e47fe0`) divide `n = 1` historical baselines by the upper bound (`267.86 s`) of an `n = 2` range (`258.88–267.86 s`) with no mean/SD, while the paragraph's opening parenthetical inaccurately labels all runs `n = 1 per run`.
- **Verification:**
  1. Verified `592.68 / 267.86 = 2.2126` and `1062.45 / 267.86 = 3.9664` in [`manuscript/horizon_paper.md:78`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L78), and confirmed the opening parenthetical states `n = 1 per run` despite reporting a two-value range `258.88–267.86 s`.
  2. **Why the omission leg is OVERTURNED:** The immediately following sentence in [`manuscript/horizon_paper.md:78`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L78) already states `2896229, 304.53 s`, `38.22% unconditional exact pair precision across all 282,163,046 emitted pairs`, and `~2.58× overall multiplicity`. The attacker should have checked the remainder of line 78 rather than isolating a single sentence; the omission applies only to [`docs/BENCHMARK_RESULTS.md:28,73`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L28).
- **Repair:** Clarifying `n = 2` (`258.88–267.86 s`, mean `263.37 s`) vs `n = 1` (`592.68 s`, `1,062.45 s`) in [`manuscript/horizon_paper.md:78`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L78) and cross-referencing the `304.53 s` (`2896229`) / `38.22%` unconditional precision result in [`docs/BENCHMARK_RESULTS.md:28,73`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L28) is endorsed (**the fix is free**).

---

## F5 -- C5: `horizon_resolve` "per-window divergence" is a sequence-blind, orientation-dependent topological heuristic clamped to `[0.01, 0.25]`

- **Ruling:** **UPHELD** (attacker's `NARROWED` verdict on `C5` stands).
- **What survives:** `estimate_window_divergence` ([`src/bin/horizon_resolve.rs:456-477`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_resolve.rs#L456-L477)) never reads DNA sequences (`PackedGraph` stores no segment sequences), instead computing `(0.02 + 0.08 * non_anchor_frac + 0.05 * inv_frac).clamp(0.01, 0.25)`. Because `non_anchor_frac` and `inv_frac` lie in `[0, 1]`, the actual output range is `[0.02, 0.15]` (neither clamp bound is even reachable), it depends on walk orientation (`+0.05` when reverse-complemented), and it inflates `MdlCostModel::block_residual_bits` ([`src/egraph.rs:410-420`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/egraph.rs#L410-L420)) at human intra-species divergence (`0.0966` bits/bp at `d = 0.01` and `0.1731` bits/bp at the actual `0.02` floor vs `0.0130` bits/bp at `d = 0.001`, i.e., `7.4×` to `13.3×`) while capping deep vertebrate divergence (`d > 0.15`).
- **Verification:** Verified [`src/bin/horizon_resolve.rs:399-407, 456-477`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_resolve.rs#L399-L407) and [`src/egraph.rs:410-420`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/egraph.rs#L410-L420).
- **Repair:** Replacing `"estimates per-window divergence (unless --divergence is passed)"` in [`ARCHITECTURE.md:28`](file:///home/rcb112/project/gemini/bench/horizon-aligner/ARCHITECTURE.md#L28) (and [`docs/PANGENOME_SUITE.md:3,25`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/PANGENOME_SUITE.md#L3), [`docs/manual.md:144`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/manual.md#L144)) with `"assigns a topological per-window divergence proxy in [0.02, 0.15] (clamped to [0.01, 0.25]) from non-anchor and reverse-step fractions (unless --divergence is passed)"` is endorsed (**the fix is free**).

---

## F6 -- Claim-list steering: `C1`–`C5` omit uncalibrated CRC (`n_cal = 15`), `KARLIN_MDL_COEFF = 0.55`, and the `18%`-GC null sweep

- **Ruling:** **UPHELD**.
- **What survives:** Extending `r4 adjudication-L3.md [F7]` and `r5 adjudication-L3.md [F6]`, the runner-authored claim list `C1`–`C5` omits the manuscript's core statistical calibration claims in [`manuscript/horizon_paper.md:10,40,61`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L10): uncalibrated CRC (`n_cal = 15`, `α = 0.07`, `λ̂ = 10.0` bits), fixed `KARLIN_MDL_COEFF = 0.55`, and the `18%`-GC random-DNA false-positive result (`0` on `seed=999`, `80–448` on `4/5` other seeds).
- **Repair:** Including CRC calibration provenance and the `18%`-GC null sweep in future claim scopes is noted.

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `horizon_filter` never calls `bocpd_gc_profile` (uses static contig GC or `0.5`), and `horizon_resolve` uses a topological heuristic; free text fix endorsed |
| F2 | UPHELD | Unweighted record-level `1.6% -> 81.2%` mixed with bp-level `99.7%` cov under seq-free `base_probs(0.5)` and `n = 1`; text qualification free, bp precision non-free |
| F3 | UPHELD | `jc69_dist = 0.222528` is conditioned on post-abstention bases without SE/CI, sums un-deduplicated runs, and accompanies `n = 1` unverified timings; free text fix endorsed |
| F4 | NARROWED | `n = 2` (`258.88–267.86 s`) vs `n = 1` speedup-ratio reporting and missing mean/SD upheld, but `manuscript/horizon_paper.md:78` already states `304.53 s`, `2.58×`, and `38.22%` in the next sentence; free text fix endorsed |
| F5 | UPHELD | `estimate_window_divergence` is sequence-blind, orientation-dependent, and bounded to `[0.02, 0.15]` (clamped `[0.01, 0.25]`), inflating `block_residual_bits` at human divergence; free text fix endorsed |
| F6 | UPHELD | `C1`–`C5` omit uncalibrated `n_cal = 15` CRC, `KARLIN_MDL_COEFF = 0.55`, and the `18%`-GC null sweep (`r4 [F7]`, `r5 [F6]`) |
