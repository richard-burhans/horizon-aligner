# Round 6 Adjudication — L5 (Inference)

**Sealed attack verified:** `reports/rounds/round-6/attack-L5-inference.md` (`sha256` prefix `f3b365963a0d`, matching brief).
**Method cards (`refs/methods/`):** None exist in this project (`docs/round_rubric.md:36`); disputes are adjudicated against live Rust source (`src/`), published docs/manuscript surfaces, and prior round adjudications (`scripts/round history --kind adjudication`).

---

## F1 -- C1: `horizon_filter` never runs BOCPD, and `horizon_wga --graph` omits DAG branching priors

- **Ruling:** **UPHELD** (attacker's **NARROWED** verdict on `C1` is upheld).
- **Verification:**
  1. In `src/bin/horizon_filter.rs:4-6`, `horizon_filter` imports only `GAP_EXTEND_BITS`, `GAP_OPEN_BITS`, `IDENTITY_PRIOR`, `base_probs`, and `joint_pair_bits` from `horizon_aligner::bocpd` and never calls `bocpd_gc_profile` or `bocpd_gc_profile_into` (`grep` confirms `0` matches). In sequence-aware mode (`IndexedSeq::from_dna`, `src/bin/horizon_filter.rs:38-59`), it computes a single global contig-wide `gc` fraction; in sequence-free mode (`score_from_counts_only`, `src/bin/horizon_filter.rs:91-136`), it uses fixed `base_probs(0.5)`.
  2. In `src/mdl/dag.rs:260-265`, `align_seq_to_dag_mdl` sums Forward log-probabilities across `unique_preds[u]` without `-log2(out_degree)` branching normalization. Note that `docs/PANGENOME_SUITE.md:17` already discloses `(omitting -log2(out_degree) branching penalties)` parenthetically, though the `C1` summary line (`docs/PANGENOME_SUITE.md:3`) omits it.
- **What survives:** `C1` cannot claim `BOCPD` for `horizon_filter` or unnormalized DAG DP as full branching-normalized Pair-HMM inference.
- **Repair:** **ENDORSED — the fix is free.** Removing `BOCPD` from `horizon_filter`'s description in `docs/PANGENOME_SUITE.md:3`, `README.md:15`, `manuscript/horizon_paper.md:52`, and `src/bin/horizon_filter.rs:1` (stating global contig-wide GC in sequence-aware mode and fixed `50%` GC in sequence-free mode) and noting the omitted `-log2(out_degree)` branching penalty in `C1` requires no re-measurement.

---

## F2 -- C2: `horizon_filter --net` 0.47s benchmark is sequence-free greedy length-netting, and `1.6% -> 81.2%` is a record-count denominator artifact

- **Ruling:** **NARROWED** (diagnosis **UPHELD**; attacker's `REFUTED` verdict narrowed because `C2` already labels the row `(seq-free)` and its raw counts are not falsified).
- **Verification:**
  1. Verified in `src/bin/horizon_filter.rs:91-136, 151-153, 253-301`: without `--query` and `--target`, `score_paf_record` falls back to `score_from_counts_only` (zero SDUST, zero BOCPD, fixed `base_probs(0.5)`), and `net_top_k_records` greedily retains `max_occ = 1` record per `>=50%` overlapping interval sorted by linear CIGAR/match-count `mdl_bits`.
  2. Verified arithmetic: `0.016 * 56,857 ≈ 910` truth-hitting records before netting vs. `0.812 * 1,004 ≈ 815` after netting; the `1.6% -> 81.2%` jump is a `56.6×` record-count denominator reduction from 1-to-1 greedy interval netting, not an MDL/BOCPD/SDUST inference.
  3. However, `docs/BENCHMARK_RESULTS.md:85` (`C2`) explicitly labels the row `` `chr6` FastGA PAF Netting (seq-free) `` and separates it from line 86 (`SDUST+MDL Filter`), and the `0.47s`, `1,004` records, `110.22 Mbp`, and `99.7%` 1-to-1 coverage numbers are not refuted. The conflation is sharpest in `README.md:124` and `manuscript/horizon_paper.md:84`, which combine the `0.47s` sequence-free and `4.73s` SDUST+BOCPD+MDL runs into one clause ending with `raising 2D HPRC truth hit rate from 1.6% to 81.2%`.
- **What survives:** `C2`'s `0.47s` / `1,004` records / `110.22 Mbp` / `99.7%` 1-to-1 coverage numbers survive as sequence-free 1-to-1 greedy interval netting; `1.6% -> 81.2%` survives only as a record-count hit rate, not base-level precision or an MDL/BOCPD/SDUST gain.
- **Repair:** **SPLIT.** Clarifying in `docs/BENCHMARK_RESULTS.md:85`, `README.md:124`, and `manuscript/horizon_paper.md:84` that `0.47s` / `1.6% -> 81.2%` is purely sequence-free 1-to-1 greedy interval netting on CIGAR/match counts and is a record-count rate (not base-level precision) is **free**. Computing and reporting base-level precision for `horizon_filter --net` (if not already recorded) is non-free.

---

## F3 -- C3: `horizon_locus` evolutionary distances (`p_dist`, `jc69_dist`) suffer from Horizon Bubble survivorship bias, silent clamping, and mislabeled `jaccard_cov`

- **Ruling:** **UPHELD** (attacker's **NARROWED** verdict on `C3` is upheld).
- **Verification:**
  1. Verified in `src/bin/horizon_locus.rs:102-145`: `p_dist` and `jc69_dist` are computed solely over `aln.total_pairs` (`aligned_bp = 231,444`, covering `~60–66%` of `galGal` `386.7 kb` and `taeGut` `349.4 kb`), omitting `~34–40%` of divergent/unanchored sequence excluded by seeding and Horizon Bubble abstention.
  2. Verified in `src/bin/horizon_locus.rs:67-74`: `jc69_distance` clamps `p_dist` to `[0.0, 0.74]`, silently emitting `~3.24` for saturated `p_dist >= 0.75` (though for the avian pair `p_dist = 0.192552 < 0.74`, so the clamp did not trigger on `C3`).
  3. Verified in `src/bin/horizon_locus.rs:128-134`: `jaccard_cov` is `(cov_x + cov_y) / (seq_x.len() + seq_y.len())` (Sørensen–Dice / total-length coverage fraction), not Jaccard (`|A ∩ B| / |A ∪ B|`); note `jaccard_cov` appears in `docs/BENCHMARK_RESULTS.md:91` rather than `C3` (`README.md:126`).
- **What survives:** `p_dist = 0.192552` and `jc69_dist = 0.222528` survive only as **aligned-core** estimates conditional on the surviving `231,444` aligned bp (`~63%` of the avian locus).
- **Repair:** **ENDORSED.** Qualifying `p_dist` and `jc69_dist` in `README.md:126` and `docs/BENCHMARK_RESULTS.md:90-91` as aligned-core estimates over `231,444` aligned bp is **free**. Returning `NA` / `f64::INFINITY` for `p_dist >= 0.75` in `jc69_distance` and renaming `jaccard_cov` in `src/bin/horizon_locus.rs` are minor code/test edits (non-free).

---

## F4 -- C4: `horizon_wga --resolve-dups` does not resolve emitted pairs, and `124.16 Mbp` on `GRCh38.chr6` vs `110.7 Mbp` `CHM13.chr6` q-arm reflects uncollapsed many-to-many alignments

- **Ruling:** **UPHELD** (attacker's **NARROWED** verdict on `C4` is upheld; properly cites Round 2 `L1-F4`, `L3-F2` and Round 5 `L4-F2`, `L5-F2`).
- **Verification:**
  1. Verified in `src/bin/horizon_wga.rs:444-453`: `--resolve-dups` executes after `maf_writer` (lines 430-438), `gfa_writer` (lines 439-443), and coverage/pair counting, and only appends `# EG:Z:` comment lines to `summary_out` and `gfa_writer` without altering emitted pairs or `unique_x_bp`.
  2. Verified arithmetically: `124,160,423` bp (`72.7%` of `GRCh38.chr6`) exceeds the entire `110,724,939` bp (`110.7 Mbp`) length of the `CHM13.chr6` q-arm target by `~13.4 Mbp`, and `282,193,564` pairs is `2.55×` the q-arm length, confirming uncollapsed many-to-many alignments.
  3. Verified in `manuscript/horizon_paper.md:78`: the `2.21×` / `3.97×` speedup (`258.88–267.86 s`) is computed on pre-fix `8cd5db0`, whereas post-fix `2896229` (`partition_colinear_blocks`) took `304.53 s` (already disclosed in the next sentence, alongside `n = 1` per run).
- **What survives:** `C4`'s raw counts on `8cd5db0` survive only as uncollapsed many-to-many alignment metrics (`# EG:Z:`-annotated, not duplicate-resolved), and the `2.21×`/`3.97×` speedup is relative to pre-fix `8cd5db0` rather than `2896229`.
- **Repair:** **ENDORSED — the wording fix is free.** Clarifying in `manuscript/horizon_paper.md:78` that `--resolve-dups` only appends `# EG:Z:` comments without modifying emitted pairs/`unique_x_bp` and that `124.16 Mbp` > `110.7 Mbp` reflects uncollapsed many-to-many alignments is free; updating the speedup ratio to `304.53 s` (`2896229`) from already-recorded timings (`1.95×` vs `bb43dd9`, `3.49×` vs `1e47fe0`) is also free (re-running `bb43dd9`/`1e47fe0` or `HEAD` would be non-free).

---

## F5 -- C5: `horizon_resolve`'s `estimate_window_divergence` is a sequence-free topological proxy that makes `MdlCostModel` divergence gating unreachable

- **Ruling:** **UPHELD** (core refutation of `C5`'s divergence-estimation and Tarjan snarl-partitioning claims is upheld; properly cites Round 2 `L3-F3`, `L5-F5`).
- **Verification:**
  1. Verified in `src/bin/horizon_resolve.rs:456-477`: `estimate_window_divergence` computes `(0.02 + 0.08 * non_anchor_frac + 0.05 * inv_frac).clamp(0.01, 0.25)` from walk step counts without inspecting DNA sequences or alignments.
  2. Verified in `src/egraph.rs:398, 410-420`: at the upper clamp `p = 0.25`, `H_2(0.25) + 0.25 log_2(3) = 1.207519 < 2.0` bits/bp (`seq_bits = 2.0 * len_bp`, with `len_bp = seg_len.max(100)` at `src/bin/horizon_resolve.rs:402`). Thus `dup(X)` is always cheaper than `seq(X, X)` under `estimate_window_divergence`, making the divergence gating tested in `mdl_cost_model_gates_dup_by_sequence_divergence` (`src/egraph.rs:606-629`, at `d = 0.75`) unreachable unless `--divergence` is passed.
  3. Verified in `src/bin/horizon_resolve.rs:329-375, 479-523`: `comp` (2-edge-connected components) is discarded after counting `num_components`, and `decompose_and_resolve` greedily slices walks up to `max_bubble_steps = 16` between `is_anchor` segments. (`PackedGraph` interning, `P`/`W` rewriting, `W`-line `end` updates, and `L`-link pruning/synthesis remain intact.)
- **What survives:** `horizon_resolve` does not estimate sequence divergence or partition walk windows by Tarjan 2-edge-connected components; `PackedGraph`, `P`/`W` rewriting, `W`-line `end` updates, and `L`-link pruning/synthesis survive.
- **Repair:** **ENDORSED — the fix is free.** Updating `ARCHITECTURE.md:28`, `docs/PANGENOME_SUITE.md:25`, and `manuscript/horizon_paper.md:50` to state that `horizon_resolve` uses a topological non-anchor/inversion proxy clamped to `[0.01, 0.25]` (always favoring `dup` over `seq(X, X)` unless `--divergence` is passed) and greedy `is_anchor` windowing (`comp` is only logged) requires no re-measurement.

---

## F6 -- Claim-list steering: `C1`–`C5` omit uncalibrated CRC/Karlin constants, `18%`-GC null false positives across other seeds, and `horizon_clean`'s ad-hoc score

- **Ruling:** **UPHELD** (meta-finding on Round 6 claim-list selection; properly cites Round 4 `L5-F7` and Round 5 `L3-F6`).
- **Verification:** Verified against `reports/rounds/round-6/ADJ-BRIEF-L5.md`, `src/bin/horizon_filter.rs:17` (`KARLIN_MDL_COEFF = 0.55`), `manuscript/horizon_paper.md:10, 18, 40, 49, 61` (`n_cal = 15`, `80–448` pairs on `4/5` other `18%`-GC seeds, `288`/`149` pairs on `30%`/`45%`-GC draws, ad-hoc `score_junction_cp_bits`), and `docs/BENCHMARK_RESULTS.md:77` (`0`-severed negative control only).
- **What survives:** The Round 6 claim list `C1`–`C5` steers away from these unresolved inferential weaknesses.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `horizon_filter` never calls `bocpd_gc_profile` (uses global contig GC or fixed `0.5`); `C1` omits `-log2(out_degree)` DAG branching omission already noted in `docs/PANGENOME_SUITE.md:17`; wording fix is free. |
| F2 | NARROWED | `0.47s` `horizon_filter --net` is sequence-free greedy netting and `1.6% -> 81.2%` is a record-count denominator artifact, but `C2` (`docs/BENCHMARK_RESULTS.md:85`) already labels the row `(seq-free)` and its raw counts are not falsified; wording clarification is free. |
| F3 | UPHELD | `p_dist = 0.192552` and `jc69_dist = 0.222528` are aligned-core estimates over `~63%` of the avian locus (`src/bin/horizon_locus.rs:102-145`), `jc69_distance` clamps at `0.74`, and `jaccard_cov` is total-length Dice coverage; doc qualification is free. |
| F4 | UPHELD | `--resolve-dups` only appends `# EG:Z:` tags (`src/bin/horizon_wga.rs:444-453`), `124.16 Mbp` > `110.7 Mbp` reflects uncollapsed many-to-many alignments, and `2.21×`/`3.97×` speedup uses pre-fix `8cd5db0` vs `2896229`; wording fix is free. |
| F5 | UPHELD | `estimate_window_divergence` is a topological proxy clamped to `[0.01, 0.25]` making `MdlCostModel` divergence gating unreachable (`1.2075 < 2.0` bits/bp), and `comp` from Tarjan bridge-finding is discarded; wording fix is free. |
| F6 | UPHELD | Round 6 claim list `C1`–`C5` omits uncalibrated `n_cal = 15` / `KARLIN_MDL_COEFF = 0.55`, `18%`-GC null false positives on `4/5` other seeds, and `horizon_clean`'s ad-hoc score / negative-control-only benchmark. |
