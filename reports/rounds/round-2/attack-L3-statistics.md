# Round 2 Attack -- Statistics Lens (L3)

**Method cards check (Rule 3):** Searched `/home/rcb112/project/gemini/bench/horizon-aligner` and parent `/home/rcb112/project/gemini/bench`; no statistical method cards exist in this project (`../pg/manuscript/wga_methods.md` belongs to a separate manuscript).
**Prior round history check (Rule 4):** Checked `reports/adversarial_log.md` and `./scripts/round history "" --kind adjudication`. Round 1 left `R1-HMM` open ("Normalize HMM transition weights in `src/mdl/dp.rs:105-167`... and derive/calibrate gap-bit penalties"). Transition normalization was implemented in `banded_mdl_gap_with_params` (`src/mdl/dp.rs:143-150`) and `align_seq_to_dag_mdl` (`src/mdl/dag.rs:88-95`), but gap-bit calibration and `mdl_xdrop_extend` remain defective (see F1).

---

## F1 -- C2: `estimate_hmm_params` is biased and truncated, `mdl_xdrop_extend` ignores `HmmParams`, and `bocpd_gc_profile` misclassifies `N` bases

- **Claim:** **C2** (`docs/BENCHMARK_RESULTS.md:5`) -- *"All results were generated out-of-the-box using the divergence-adaptive Pair-HMM (`HmmParams`: `identity_prior`, `gap_open_bits`, and `gap_extend_bits` estimated from LCB anchors with Beta prior shrinkage in `src/mdl/mod.rs`), BOCPD compositional profiles (`src/bocpd.rs`), and normalized Pair-HMM transitions (`src/mdl/dp.rs`, `src/mdl/dag.rs`)."*
- **Verdict:** **REFUTED**
- **Concrete flaws & artifacts:**
  1. **`mdl_xdrop_extend` ignores `HmmParams` and uses unnormalized transitions (`src/mdl/dp.rs:369-482`, `src/mdl/mod.rs:705`):** Every LCB flank extension and macro-bubble flank (`xdrop_jobs`) calls `mdl_xdrop_extend`, which hardcodes the static constants `IDENTITY_PRIOR` (`0.65`), `GAP_OPEN_BITS` (`-4.5`), and `GAP_EXTEND_BITS` (`-0.8`) rather than `hmm_params`, and omits `mm_bits`/`im_bits` match-transition costs. Thus, "all results" do not use `HmmParams` or normalized Pair-HMM transitions.
  2. **Selection bias, positional truncation, and broken indel estimator in `estimate_hmm_params` (`src/mdl/mod.rs:759-843`):**
     - `identity_prior` is estimated from `lcb.iter().take(64)` over `0..span` (`span = 15`), which are spaced-seed anchor hits (`MASK_PRIMARY_11_15` requires $\ge 11/15 = 73.3\%$ identity by construction). To offset this upward selection bias, line 818 applies an ad-hoc clamp `.clamp(0.60, 0.72)` around a `Beta(195, 105)` prior, capping `identity_prior` at `0.72` even on Human-Chimp (`>98%` identity).
     - Lines 790 and 805 `break` as soon as `total_bases >= 4000`, so on large chromosomes only the first few positionally ordered LCBs are sampled.
     - Lines 794-804 count `if dx != dy { indel_events += 1; total_indel_len += dx.abs_diff(dy); }` across consecutive anchors up to `500` bp apart. Equal-length compensating indels (`dx == dy`) register `0` events, multiple indels in a 500-bp span collapse into `1` event with net difference `|dx - dy|`, and `gap_extend_bits` is clamped to `[-1.8, -0.5]` (mean indel length `[1.40, 3.41]` bp) while `gap_open_bits` is clamped to `[-5.5, -4.0]`.
  3. **`bocpd_gc_profile` treats `N`/ambiguous bases as A/T (`src/bocpd.rs:72`):** `let is_gc = matches!(b.to_ascii_uppercase(), b'G' | b'C');` classifies `'N'` as non-GC, biasing GC estimates downward in masked or centromeric sequence.
- **Correction:** Pass `hmm_params` into `mdl_xdrop_extend` (or replace it with normalized `HmmParams` transitions); sample inter-anchor columns rather than seed-hit positions without positional truncation; estimate indel rates from aligned gap runs rather than `dx != dy`; and skip non-ACGT bases in `bocpd_gc_profile`. Rewrite `docs/BENCHMARK_RESULTS.md:5` to state that `mdl_xdrop_extend` uses fixed defaults (`0.65`, `-4.5`, `-0.8`) and `estimate_hmm_params` estimates parameters from capped LCB-anchor windows.

---

## F2 -- C5: `5,894,191` aligned pairs on a `5 Mbp` locus is 117.9% multi-mapping inflation, and `horizon_clean` never exercised its statistical gate

- **Claim:** **C5** (`docs/BENCHMARK_RESULTS.md:23`) -- *"Real HPRC `chr22` Pericentromere 5 Mbp (`GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M`): **23.45s**, **3,084.5 MB** peak RSS, **5,894,191** aligned pairs, `246` fwd / `184` rev LCBs, `17,789` Horizon Bubbles (`horizon_clean`: `6,030` links kept, **0.01s**, **7.8 MB** RSS)."*
- **Verdict:** **NARROWED**
- **Concrete flaws & artifacts:**
  1. **Many-to-many pair count exceeds locus length:** `horizon_wga --resolve-dups` (run without 1-to-1 base netting) deduplicates `(x, y)` coordinate pairs (`src/mdl/mod.rs:745-746`), allowing each repetitive pericentromeric base to pair with multiple target coordinates across `246` forward and `184` reverse LCBs. Reporting `5,894,191` pairs on a `5,000,000` bp locus (`117.9%` of locus length) without unique-base coverage or multiplicity conflates redundant multi-mapping with alignment coverage (as also occurs in `docs/BENCHMARK_RESULTS.md:22,26`).
  2. **Null-regime `horizon_clean` test:** `horizon_clean` (`src/bin/horizon_clean.rs:228-230`) only severs links with `walk_support == 0`. Because `chr22_peri_5m.gfa` was extracted from `P`/`W` walks (`sup > 0` for all `6,030` links), `should_sever` is identically `false` and the BOCPD/MDL changepoint test (`score_junction_cp_bits`) was never triggered (`0` severed).
- **Correction:** Report unique aligned `GRCh38`/`CHM13` bases alongside total `(x, y)` pairs in `docs/BENCHMARK_RESULTS.md:23`, and note that `horizon_clean` kept all `6,030` links because all were walk-supported (`sup > 0`).

---

## F3 -- C6: `308,760.9` MDL bits saved is an artifact of double-charging segment sequences and hardcoding `divergence = 0.02`

- **Claim:** **C6** (`docs/BENCHMARK_RESULTS.md:64`) -- *"| GFA E-Graph Resolver | `horizon_resolve --collapse-dups` | 200 segs, 20 walks x 225 steps (`4,500` steps) | `280` windows, `840` dups collapsed, `308,760.9` bits saved | **24.0 MB** | 4.48s |"*
- **Verdict:** **REFUTED**
- **Concrete flaws & artifacts:**
  1. **Double-charging segment sequences in `MdlCostModel::block_bits` (`src/egraph.rs:325-328`):** Every walk step referencing segment `B` is charged `catalog_bits + 2.0 * len_bp`, as if the full nucleotide sequence of `B` were transmitted anew on every step rather than referenced from the GFA `S`-line catalog. Collapsing `dup(B, B)` then credits `~ (2.0 - H_2(d) - d log2 3) * len_bp` bits saved per duplicate.
  2. **Hardcoded divergence (`src/bin/horizon_resolve.rs:68, 202-208`):** `horizon_resolve` never estimates sequence divergence from the GFA; it assigns the CLI default `--divergence 0.02` to every segment. At `d = 0.02`, `H_2(0.02) + 0.02 * log2(3) = 0.1731` bits/bp, yielding `1.8269` bits/bp saved (~`367.6` bits per 200-bp dup $\times 840 = 308,784$ bits, matching the reported `308,760.9` bits).
- **Correction:** Charge `catalog_bits` (not `2.0 * len_bp`) per walk reference in `horizon_resolve`, or measure actual sequence divergence between tandem copies. Remove `308,760.9 bits saved` as an empirical MDL compression claim in `docs/BENCHMARK_RESULTS.md:64` or label it as a synthetic model cost under fixed `d = 0.02`.

---

## F4 -- C7: `align_seq_to_dag_mdl` omits graph-branching penalties and resets BOCPD at every node

- **Claim:** **C7** (`docs/PANGENOME_SUITE.md:18`) -- *"`align_seq_to_dag_seeded` uses `MASK_W10` + `MASK_RY20` spaced seeds and forward/backward reachability closure to isolate relevant DAG components before running `align_seq_to_dag_mdl`."*
- **Verdict:** **NARROWED**
- **Concrete flaws & artifacts:**
  1. **Unnormalized graph-branching transitions (`src/mdl/dag.rs:151-161, 260-271`):** Forward and backward passes sum/maximize over predecessor/successor nodes with `0.0` transition bits (omitting `-log2(out_degree(p))`), so branching DAGs inflate `z_fwd` and Viterbi scores relative to linear sequences.
  2. **Per-node BOCPD reset (`src/mdl/dag.rs:73-81`):** `bocpd_gc_profile(&nd.seq, 1.0 / 120.0, 80)` is initialized fresh from `Beta(1, 1)` (`GC = 0.50`) on each node, discarding predecessor GC context on short variant nodes.
- **Correction:** Add `-log2(out_degree(p))` branching transition costs (or document that DAG node transitions are unnormalized) and propagate BOCPD sufficient statistics across DAG edges.

---

## F5 -- C8: `~0.29 GB/Mbp` misapplies the `20 Mbp` q-arm rate to the cancelled `281.5 Mbp` full-`chr6` run

- **Claim:** **C8** (`docs/BENCHMARK_RESULTS.md:26`) -- *"(After replacing `HashSet<(usize, usize)>` with sorted/deduplicated `Vec<(usize, usize)>`, full 170.8 Mbp `GRCh38.chr6` vs 110.7 Mbp `CHM13.chr6` still reached **26.5 GB** RSS and `47:33` CPU before cancellation on the 30 GiB host due to `ContigIndex`'s per-position `HashMap<u64, Vec<usize>>` seed maps, ~0.29 GB/Mbp.)"*
- **Verdict:** **NARROWED**
- **Concrete flaws & artifacts:**
  - `26.5 GB / (170.8 + 110.7) Mbp = 0.094 GB/Mbp` (and `26.5 / 170.8 = 0.155 GB/Mbp`, `26.5 / 110.7 = 0.239 GB/Mbp`). The `~0.29 GB/Mbp` rate is actually from the `20 Mbp` `chr6` q-arm run (`5,863.3 MB / 20 Mbp = 0.293 GB/Mbp` on the same line), not the cancelled full-`chr6` run (whose true per-Mbp rate is also unknown because it was cancelled before completion).
- **Correction:** Replace `~0.29 GB/Mbp` in `docs/BENCHMARK_RESULTS.md:26` with `0.094 GB/Mbp` (across `281.5 Mbp` combined) or move `~0.29 GB/Mbp` to the `20 Mbp` q-arm clause.

---

## F6 -- Claim-list steering: C1–C8 omit held-out WGA test regressions and uncalibrated `KARLIN_MDL_COEFF`

- **Observation:** Per brief instructions, C1–C8 steer away from the Held-Out Test Sets table in `docs/BENCHMARK_RESULTS.md:38-48` (where recall drops to `50.5%` on Horse and `63.6%` recall / `61.0%` precision on Dog `chr22 1MB`, and B10K 363-way is `Pending Data`) and from the uncalibrated `KARLIN_MDL_COEFF = 0.55` and `min_delta_bits = 10.0` constants in `src/mdl/mod.rs:25,68`.

---

## Claim Verdicts

| Claim | Verdict |
| :--- | :--- |
| C1 | UNCHALLENGED |
| C2 | REFUTED |
| C3 | UNCHALLENGED |
| C4 | UNCHALLENGED |
| C5 | NARROWED |
| C6 | REFUTED |
| C7 | NARROWED |
| C8 | NARROWED |
