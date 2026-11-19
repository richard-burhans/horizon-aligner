# Round 2 ADJUDICATION -- Statistics Lens (L3)

- **Attack reviewed:** `reports/rounds/round-2/attack-L3-statistics.md` (sealed digest `3952ea463b86`, verified via `sha256sum`: `3952ea463b863f756d559c03e5d0d400a15321cd0a433b6f4fd0148a6deecd1a`).
- **Method cards (`refs/methods/`):** Neither `/home/rcb112/project/gemini/bench/horizon-aligner/refs/methods/` nor `/home/rcb112/project/gemini/bench/refs/methods/` exists; no project method card governs these statistical disputes.
- **Prior round history:** Checked `reports/adversarial_log.md` and `./scripts/round history "" --kind adjudication`. Round 1 left `R1-HMM` open ("Normalize HMM transition weights in `src/mdl/dp.rs:105-167` ... and derive/calibrate gap-bit penalties"). Transition normalization (`mm_bits`, `im_bits`) was added to `banded_mdl_gap_with_params` (`src/mdl/dp.rs:143-150`) and `align_seq_to_dag_mdl` (`src/mdl/dag.rs:88-95`), but `mdl_xdrop_extend` (`src/mdl/dp.rs:369-482`) and gap-bit calibration were not fixed.

---

## Finding-by-Finding Rulings

### F1 -- C2 (`estimate_hmm_params` bias/truncation, `mdl_xdrop_extend` ignoring `HmmParams`, and `bocpd_gc_profile` misclassifying `N` bases)
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; repair **NARROWED** by cost.
  - **Diagnosis:** Verified against the codebase. `mdl_xdrop_extend` (`src/mdl/dp.rs:369-482`, called at `src/mdl/mod.rs:705` for every LCB flank and macro-bubble flank) takes no `HmmParams`, hardcodes `IDENTITY_PRIOR` (`0.65`), `GAP_OPEN_BITS` (`-4.5`), and `GAP_EXTEND_BITS` (`-0.8`), and omits `mm_bits`/`im_bits` transition costs, directly falsifying C2's claim that *"All results were generated ... using the divergence-adaptive Pair-HMM (`HmmParams`...) ... and normalized Pair-HMM transitions."* In `estimate_hmm_params` (`src/mdl/mod.rs:759-843`), `identity_prior` is indeed sampled from `lcb.iter().take(64)` over `0..span` (`span = 15`) spaced-seed hits, truncated at `total_bases >= 4000` (lines 790, 805), and clamped to `[0.60, 0.72]` around `Beta(195, 105)` (line 818); indel estimation (`dx != dy`, lines 794-804) misses equal-length compensating indels (`dx == dy`), collapses multi-indel intervals up to `500` bp into one event, and clamps `gap_extend_bits` to `[-1.8, -0.5]` (mean indel length `1.40-3.41` bp) and `gap_open_bits` to `[-5.5, -4.0]`. In `src/bocpd.rs:72`, `matches!(b.to_ascii_uppercase(), b'G' | b'C')` treats `'N'` and other ambiguous bytes as non-GC, biasing GC estimates downward in masked or centromeric sequence.
  - **Repair:** Endorsed on the documentation fix to `docs/BENCHMARK_RESULTS.md:5` (clarify that `mdl_xdrop_extend` uses fixed defaults `0.65`/`-4.5`/`-0.8` without normalized transitions and that `estimate_hmm_params` samples capped LCB-anchor windows with `identity_prior` clamped to `[0.60, 0.72]`) -- **the doc fix is free**. The code fixes (passing `hmm_params` into `mdl_xdrop_extend`, sampling inter-anchor columns, replacing `dx != dy` indel counting, and skipping non-ACGT bases in `bocpd_gc_profile`) are valid but **not free**: they alter alignment scores and require re-running the benchmark suite in `docs/BENCHMARK_RESULTS.md`.
- **What survives:** C2 is refuted because `mdl_xdrop_extend` uses fixed unnormalized constants on every `horizon_wga` run and `estimate_hmm_params` / `bocpd_gc_profile` carry the stated estimation biases.

---

### F2 -- C5 (`5,894,191` aligned pairs on a `5 Mbp` locus and null-regime `horizon_clean` test)
- **Verdict:** **NARROWED**
- **Ruling:** Diagnosis **UPHELD**; repair **NARROWED** by cost.
  - **Diagnosis:** Verified against `src/mdl/mod.rs:745-746` and `src/bin/horizon_clean.rs:228-230`. `horizon_wga --resolve-dups` deduplicates `(x, y)` coordinate pairs rather than unique reference/query bases, so `5,894,191` pairs on a `5,000,000` bp pericentromeric locus (`117.9%`) reflects many-to-many multi-mapping across `246` fwd and `184` rev LCBs. In `src/bin/horizon_clean.rs:228-230`, `should_sever` is `(cfg.prune_unsupported && sup == 0) || (sup == 0 && cp_bits >= cfg.min_cp_bits)`, which is identically `false` whenever every link in the `P`/`W`-extracted subgraph has `sup > 0` (leaving `0` links severed and the BOCPD/MDL changepoint gate unexercised). Narrowed because C5 literally reports `5,894,191 aligned pairs` (not `5,894,191 bp` or `117.9%` coverage) and `6,030 links kept`.
  - **Repair:** Endorsed on the documentation clarification in `docs/BENCHMARK_RESULTS.md:23` (label `5,894,191` as many-to-many `(x, y)` pairs rather than unique-base coverage, and note that `horizon_clean` kept all `6,030` links because all were walk-supported `sup > 0` so `should_sever` never triggered) -- **the doc fix is free**. Reporting unique aligned `GRCh38`/`CHM13` bases alongside pair counts requires re-running `horizon_wga` on that locus and is **not free**.
- **What survives:** C5's raw numbers stand, but `5,894,191` pairs is a many-to-many multi-mapping count and `horizon_clean` (`0` severed) did not exercise its statistical changepoint gate.

---

### F3 -- C6 (`308,760.9` MDL bits saved artifact from double-charging segments and hardcoded `divergence = 0.02`)
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; repair **NARROWED** by cost.
  - **Diagnosis:** Verified against `src/egraph.rs:322-365` and `src/bin/horizon_resolve.rs:68, 202-208`. `horizon_resolve` populates `model.block_stats` for every segment in the window using the CLI default `divergence = 0.02` (never estimated from GFA sequences). `MdlCostModel::block_bits` (`src/egraph.rs:325-328`) then charges `catalog_bits + 2.0 * len_bp` on every walk step referencing a segment, and `dup_cost` credits `~(2.0 - H_2(0.02) - 0.02 * log2(3)) * len_bp = 1.8269` bits/bp (`~367.6` bits per 200-bp duplicate $\times 840 \approx 308,784$ bits, matching `308,760.9` bits). Thus `308,760.9` bits saved is a synthetic model cost under fixed `d = 0.02` and `2.0 * len_bp` per-step sequence charging, not an empirical GFA compression measurement.
  - **Repair:** Endorsed on the documentation fix to `docs/BENCHMARK_RESULTS.md:64` (remove `308,760.9 bits saved` as an empirical MDL compression claim, or label it explicitly as synthetic model cost under fixed `--divergence 0.02` and `2.0 * len_bp` per-step sequence charging) -- **the doc fix is free**. Changing `block_bits` to charge `catalog_bits` only or estimating real sequence divergence in `horizon_resolve` is a code change requiring benchmark re-execution (**not free**).
- **What survives:** C6's `308,760.9 bits saved` is refuted as an empirical MDL compression result and survives only as a synthetic cost-model output under hardcoded `d = 0.02`.

---

### F4 -- C7 (`align_seq_to_dag_mdl` omits graph-branching penalties and resets BOCPD at every node)
- **Verdict:** **NARROWED**
- **Ruling:** Diagnosis **UPHELD**; repair **NARROWED** by scope and cost.
  - **Diagnosis:** Verified against `src/mdl/dag.rs:73-81, 151-161, 260-271`. `align_seq_to_dag_mdl` runs `bocpd_gc_profile(&nd.seq, 1.0 / 120.0, 80)` independently per DAG node (resetting to `Beta(1, 1)` at every node boundary) and sums/maximizes over predecessor/successor nodes with `0.0` branching transition bits (omitting `-log2(out_degree(p))`). Narrowed because C7 (`docs/PANGENOME_SUITE.md:18`) only states that `align_seq_to_dag_seeded` uses `MASK_W10` + `MASK_RY20` spaced seeds and reachability closure before calling `align_seq_to_dag_mdl`; it does not claim normalized DAG branching or cross-node BOCPD propagation.
  - **Repair:** Endorsed on documenting in `docs/PANGENOME_SUITE.md:18` that `align_seq_to_dag_mdl` uses unnormalized node-to-node transitions (`0.0` branch bits) and per-node BOCPD initialization -- **the doc fix is free**. Modifying `src/mdl/dag.rs` to add `-log2(out_degree(p))` and propagate BOCPD state across edges is **not free** and would change the `124.4` bits reported in `docs/BENCHMARK_RESULTS.md:65`.
- **What survives:** C7's description of `align_seq_to_dag_seeded` stands, narrowed by the undocumented statistical limitations inside `align_seq_to_dag_mdl`.

---

### F5 -- C8 (`~0.29 GB/Mbp` misapplies the `20 Mbp` q-arm rate to the cancelled `281.5 Mbp` full-`chr6` run)
- **Verdict:** **NARROWED**
- **Ruling:** Both diagnosis and doc repair **UPHELD**.
  - **Diagnosis:** Verified by arithmetic on `docs/BENCHMARK_RESULTS.md:26`. Dividing `26.5 GB` RSS by `170.8 + 110.7 = 281.5 Mbp` combined yields `0.094 GB/Mbp` (a lower bound since the run was cancelled before completion), whereas `~0.29 GB/Mbp` matches the `20 Mbp` `chr6` q-arm run on the same line (`5,863.3 MB / 20 Mbp = 0.293 GB/Mbp`). Narrowed because the `26.5 GB` RSS and `47:33` CPU observation itself is not refuted.
  - **Repair:** Endorsed on `docs/BENCHMARK_RESULTS.md:26`: replace `~0.29 GB/Mbp` on the cancelled full-`chr6` parenthetical with `>= 0.094 GB/Mbp` at cancellation across `281.5 Mbp` combined (or move `~0.29 GB/Mbp` to the `20 Mbp` q-arm clause). **The fix is free.**
- **What survives:** C8's `26.5 GB` / `47:33` cancellation report stands, with the misattributed `~0.29 GB/Mbp` rate corrected to `>= 0.094 GB/Mbp` at cancellation.

---

### F6 -- Claim-list steering (C1-C8 omit held-out WGA test regressions and uncalibrated `KARLIN_MDL_COEFF`)
- **Verdict:** **UPHELD**
- **Ruling:** Structural observation **UPHELD**. C1-C8 omit the Held-Out Test Sets table (`docs/BENCHMARK_RESULTS.md:38-48`, where recall drops to `50.5%` on Horse and `63.6%` recall / `61.0%` precision on Dog `chr22 1MB`, and B10K 363-way is `Pending Data`) as well as the uncalibrated `KARLIN_MDL_COEFF = 0.55` and `min_delta_bits = 10.0` constants (`src/mdl/mod.rs:25, 68`), even though `docs/BENCHMARK_RESULTS.md:5` explicitly notes `KARLIN_MDL_COEFF = 0.55` immediately after the C2 excerpt. No direct claim in C1-C8 is refuted by this observation; it is recorded for the completeness critic and future rounds.
- **What survives:** The claim-list omission is recorded for completeness coverage.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `mdl_xdrop_extend` uses fixed unnormalized constants (`0.65`, `-4.5`, `-0.8`), `estimate_hmm_params` truncates at `4000` bases and clamps `[0.60, 0.72]`, and `bocpd_gc_profile` treats `'N'` as non-GC; doc fix to `docs/BENCHMARK_RESULTS.md:5` is free, code fixes require benchmark re-runs |
| F2 | NARROWED | `5,894,191` pairs is many-to-many multi-mapping and `horizon_clean` (`sup > 0`) never triggered its changepoint gate, though C5 states "aligned pairs"; doc clarification is free |
| F3 | UPHELD | `308,760.9` bits saved is a synthetic artifact of `2.0 * len_bp` per walk step and hardcoded `--divergence 0.02`; labeling/removing the empirical claim in `docs/BENCHMARK_RESULTS.md:64` is free |
| F4 | NARROWED | `align_seq_to_dag_mdl` omits `-log2(out_degree)` branching penalties and resets BOCPD per node, though C7 only describes `align_seq_to_dag_seeded`; doc clarification is free |
| F5 | NARROWED | `26.5 GB / 281.5 Mbp = 0.094 GB/Mbp` (lower bound at cancellation), not `~0.29 GB/Mbp` (from the `20 Mbp` q-arm run); doc fix in `docs/BENCHMARK_RESULTS.md:26` is free |
| F6 | UPHELD | C1-C8 omit the Held-Out Test Sets regressions (`docs/BENCHMARK_RESULTS.md:38-48`) and `KARLIN_MDL_COEFF = 0.55` / `min_delta_bits = 10.0` (`src/mdl/mod.rs:25, 68`) |
