# Round 6 Attack — L5 (Inference)

**Scope & Environment Notes (Rubric Rules 2–4):**
- **Literature (`../library`):** No checkout exists at `../library` (`NO_LIBRARY`), so no held summaries or full texts were consulted beyond `manuscript/references.bib` and `refs/catalog_index.json`. No prior-work assertions are made from memory.
- **Method cards:** Per `docs/round_rubric.md:36`, this project keeps **no method cards** yet; inferential claims are evaluated against the frozen snapshots, live Rust source (`src/`), benchmark scripts (`scripts/benchmarks/`), and `manuscript/references.bib`.
- **Prior rounds checked:** `scripts/round history` (`--kind adjudication`) was queried for `horizon_filter`, `horizon_locus`, `horizon_resolve`, `chr6`, `MdlCostModel`, `resolve-dups`, and `bocpd`. Prior adjudications bearing on these claims (Round 2 `L1-F4`, `L3-F2`, `L3-F3`, `L5-F5`; Round 4 `L5-F7`; Round 5 `L4-F2`, `L5-F2`) are cited where applicable below.

---

## F1 -- C1: `horizon_filter` never runs BOCPD, and `horizon_wga --graph` omits DAG branching priors

- **Claim:** `C1` (`docs/PANGENOME_SUITE.md:3`)
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. `C1` describes `horizon_filter` as a `"streaming PAF SDUST+BOCPD+MDL filter & netter"` (echoed in `src/bin/horizon_filter.rs:1`, `README.md:15`, and `manuscript/horizon_paper.md:52`). In `src/bin/horizon_filter.rs:4-6`, `horizon_filter` imports only `GAP_EXTEND_BITS`, `GAP_OPEN_BITS`, `IDENTITY_PRIOR`, `base_probs`, and `joint_pair_bits` from `horizon_aligner::bocpd` — it **never imports or calls `bocpd_gc_profile` or `bocpd_gc_profile_into`**. In sequence-aware mode (`--query` / `--target`), `IndexedSeq::from_dna` (`src/bin/horizon_filter.rs:38-59`) computes a **single global contig-wide GC fraction** (`gc`), and in sequence-free mode (`score_from_counts_only`, `src/bin/horizon_filter.rs:91-136`), it uses fixed `base_probs(0.5)` (50% GC) and fixed `IDENTITY_PRIOR` without SDUST or BOCPD. Inferring local compositional changepoint tracking from `horizon_filter` is therefore unsupported by the implementation.
  2. `C1` claims `horizon_wga --graph` implements `"seeded DAG/POA MDL DP in src/mdl/dag.rs"`. As `docs/PANGENOME_SUITE.md:17` notes parenthetically (but `C1` omits), `align_seq_to_dag_mdl` (`src/mdl/dag.rs:260-265`) sums Forward log-probabilities across incoming DAG predecessors without subtracting `-log2(out_degree)` branching transition costs. In a Pair-HMM / MDL generative model over a DAG, omitting branching normalization causes total transition probability mass to exceed `1.0` at bifurcations, biasing path inference toward higher-degree subgraphs.
- **Evidence:** `src/bin/horizon_filter.rs:4-6, 38-59, 91-136, 151-153`; `grep -rn "bocpd_gc_profile" src/bin/horizon_filter.rs` returns `0` matches; `src/mdl/dag.rs:260-265`; `docs/PANGENOME_SUITE.md:3, 17`.
- **Correction:** Remove `BOCPD` from `horizon_filter`'s description in `docs/PANGENOME_SUITE.md:3` (and `README.md:15`, `manuscript/horizon_paper.md:52`, `src/bin/horizon_filter.rs:1`), stating instead that `horizon_filter` uses a global contig-wide GC fraction in sequence-aware mode and fixed `50%` GC in sequence-free mode; and state in `C1` that `horizon_wga --graph` omits DAG branching transition normalization (`-log2(out_degree)`).

---

## F2 -- C2: `horizon_filter --net` 0.47s benchmark is sequence-free greedy length-netting, and `1.6% -> 81.2%` is a record-count denominator artifact

- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:85`)
- **Verdict:** **REFUTED**
- **Concrete flaw:** `C2` reports that `horizon_filter --net` on `93 MB` (`56,857` `FastGA` records, `GRCh38.chr6` × `CHM13.chr6`) produces `1,004` records (`110.22 Mbp`), raising `"2D HPRC truth hit rate"` from `1.6%` to `81.2%` (`99.7%` 1-to-1 cov) in `0.47s`, within a section/tool framed as `"Streaming MDL + BOCPD + SDUST"` (`src/bin/horizon_filter.rs:1`, `README.md:15`, `docs/VGP_INTEGRATION_PLAN.md:72`).
  1. When run **sequence-free** (`0.47s`, without `--query` and `--target`), `score_paf_record` (`src/bin/horizon_filter.rs:151-153`) immediately falls back to `score_from_counts_only` (`src/bin/horizon_filter.rs:91-136`), which performs **zero SDUST repeat masking**, **zero BOCPD**, and **uses fixed 50% GC** (`base_probs(0.5)`). `net_top_k_records` (`src/bin/horizon_filter.rs:253-301`) then sorts records by `mdl_bits` (a linear function of match/mismatch/gap counts, i.e., net alignment length) and greedily retains `max_occ = 1` record per `>=50%` overlapping query/target interval.
  2. Thus, the `1.6% -> 81.2%` jump is an **arithmetic record-count denominator effect** of standard 1-to-1 greedy interval netting (dropping `55,853` short secondary/repeat `FastGA` records out of `56,857` while retaining `1,004` long primary records covering `110.22 Mbp`), **not** an inference from MDL, BOCPD, or SDUST. By comparison, the sequence-aware SDUST+MDL run (`docs/BENCHMARK_RESULTS.md:86`) takes `4.73s` (`10×` slower), retains `848` records (`109.76 Mbp`, slightly *less* than `110.22 Mbp`), and reports no 2D HPRC truth hit rate.
- **Evidence:** `src/bin/horizon_filter.rs:91-136, 151-153, 253-301`; `docs/BENCHMARK_RESULTS.md:85-86`; `manuscript/horizon_paper.md:84`.
- **Correction:** In `docs/BENCHMARK_RESULTS.md:85`, `README.md:124`, and `manuscript/horizon_paper.md:84`, state explicitly that the `0.47s` `1.6% -> 81.2%` result is **purely 1-to-1 greedy interval netting** on CIGAR/match counts (without SDUST, BOCPD, or sequence-level MDL), and report base-level precision rather than record-count hit rate alone.

---

## F3 -- C3: `horizon_locus` evolutionary distances (`p_dist`, `jc69_dist`) suffer from Horizon Bubble survivorship bias, silent clamping, and mislabeled `jaccard_cov`

- **Claim:** `C3` (`README.md:126`)
- **Verdict:** **NARROWED**
- **Concrete flaw:** `C3` reports `p_dist = 0.192552` and `jc69_dist = 0.222528` on the avian 2-way locus (`galGal` `386.7 kb` × `taeGut` `349.4 kb`, `231,444` aligned bp).
  1. **Survivorship / abstention bias:** In `src/bin/horizon_locus.rs:102-145`, `p_dist = 1.0 - (matches as f64 / aligned_bp as f64)` and `jc69_dist = jc69_distance(p_dist)` are computed **exclusively over `aln.total_pairs` (`aligned_bp`)**. On the avian 2-way locus, `aligned_bp = 231,444` covers only `~60–66%` of the two sequences (`349.4–386.7 kb`); the remaining `~34–40%` of divergent or unanchored sequence was excluded by seeding and Horizon Bubble abstention. Inferring a locus-wide evolutionary distance (`jc69_dist = 0.222528`) from only the conserved, alignable ~63% of the locus **systematically underestimates** true divergence.
  2. **Silent saturation clamping:** `jc69_distance` (`src/bin/horizon_locus.rs:67-74`) clamps `p_dist` at `0.74` (`p_dist.clamp(0.0, 0.74)`), so saturated alignments with $p \ge 0.75$ silently emit a finite `jc69_dist ≈ 3.24` instead of `NA` or `inf`.
  3. **Mislabeled `jaccard_cov`:** `compute_pair_metrics` (`src/bin/horizon_locus.rs:128-134`) defines `jaccard_cov = (cov_x + cov_y) as f64 / (seq_x.len() + seq_y.len()) as f64`, which is **Sørensen–Dice / total-length coverage fraction**, not Jaccard similarity ($|A \cap B| / |A \cup B|$).
- **Evidence:** `src/bin/horizon_locus.rs:67-74, 102-145`; `README.md:126`; `docs/BENCHMARK_RESULTS.md:90-91`.
- **Correction:** Qualify `p_dist` and `jc69_dist` in `README.md:126` and `docs/BENCHMARK_RESULTS.md:90-91` as **aligned-core estimates** conditional on the surviving `231,444` aligned bp (`~63%` of the locus, subject to Horizon Bubble abstention bias), return `NA` when `p_dist >= 0.75` in `jc69_distance`, and rename `jaccard_cov` (or compute true Jaccard).

---

## F4 -- C4: `horizon_wga --resolve-dups` does not resolve emitted pairs, and `124.16 Mbp` on `GRCh38.chr6` vs `110.7 Mbp` `CHM13.chr6` q-arm reflects uncollapsed many-to-many alignments

- **Claim:** `C4` (`manuscript/horizon_paper.md:78`)
- **Verdict:** **NARROWED**
- **Prior round context:** Round 2 (`L1-F4`, `L3-F2`) established that `horizon_wga` emits many-to-many `(x, y)` pairs, and Round 5 (`L4-F2`, `L5-F2`) established that `8cd5db0` (`258.88–267.86 s`) was superseded by `2896229` (`304.53 s`) after fixing `partition_colinear_blocks` (`48.63%` conditioned / `38.22%` unconditional exact precision on `chr6`).
- **Concrete flaw (new in C4):**
  1. `C4` attributes the `chr6` alignment output (`282,193,564` aligned pairs, `124,160,423` bp on `GRCh38.chr6`, `109,464,447` bp on `CHM13.chr6` q-arm) to `horizon_wga --resolve-dups`. In `src/bin/horizon_wga.rs:444-453`, `--resolve-dups` **does not collapse, filter, or alter `aln.aligned_runs`, `total_pairs`, `unique_x_bp`, or `maf_out`**; it only appends `# EG:Z:` comment tags to `summary_out` / `gfa_out` after alignment.
  2. Because `CHM13.chr6` in this benchmark is **only the `110.7 Mbp` q-arm** (`61401689–172126628`), covering `124.16 Mbp` (`72.7%`) of `GRCh38.chr6` (`170.8 Mbp`) is mathematically impossible under 1-to-1 orthology: at least `13.4 Mbp` (and `282.2M` pairs = `2.55×` the q-arm length) represents **many-to-many repeat/p-arm-to-q-arm alignments** onto the `110.7 Mbp` q-arm.
  3. The `2.21×` / `3.97×` speedup (`258.88–267.86 s` at `8cd5db0` vs `592.68 s` at `bb43dd9` and `1,062.45 s` at `1e47fe0`) compares across commits with different alignment outputs (`8cd5db0` predates the `2896229` `partition_colinear_blocks` fix, which raised runtime to `304.53 s`; `n = 1` per run).
- **Evidence:** `src/bin/horizon_wga.rs:444-453`; `manuscript/horizon_paper.md:78`; `docs/BENCHMARK_RESULTS.md:28, 73-74`.
- **Correction:** State in `manuscript/horizon_paper.md:78` that `--resolve-dups` only appends `# EG:Z:` annotations without modifying the `282.2M` emitted pairs or `unique_x_bp`, note that `124.16 Mbp` on `GRCh38.chr6` exceeds the `110.7 Mbp` `CHM13.chr6` q-arm target due to uncollapsed many-to-many repeat alignments, and report speedup on a single pinned post-fix commit (`2896229` / `HEAD`).

---

## F5 -- C5: `horizon_resolve`'s `estimate_window_divergence` is a sequence-free topological proxy that makes `MdlCostModel` divergence gating unreachable

- **Claim:** `C5` (`ARCHITECTURE.md:28`)
- **Verdict:** **REFUTED**
- **Prior round context:** Round 2 (`L3-F3`, `L5-F5`) previously upheld that `horizon_resolve` hardcoded `divergence = 0.02`. Commit `6e5d80c` added `estimate_window_divergence` (`src/bin/horizon_resolve.rs:456-477`) in response.
- **Concrete flaw (in `HEAD`):**
  1. `estimate_window_divergence` (`src/bin/horizon_resolve.rs:456-477`) **does not inspect DNA sequences, base substitutions, or alignments**. It computes a purely topological linear proxy from walk step counts: `(0.02 + 0.08 * non_anchor_frac + 0.05 * inv_frac).clamp(0.01, 0.25)`. Describing this in `ARCHITECTURE.md:28` as `"estimates per-window divergence"` conflates graph topology with evolutionary sequence divergence.
  2. Because `estimate_window_divergence` is **clamped to `[0.01, 0.25]`**, in `MdlCostModel::block_residual_bits` (`src/egraph.rs:410-420`) the per-base residual cost is bounded by $H_2(0.25) + 0.25 \log_2 3 \approx 1.207 < 2.0$ bits/bp (`seq_bits = 2.0 * len_bp`, `src/egraph.rs:398`). Thus, for any segment of length $\ge 4$ bp, `dup(X)` is **always** cheaper than `seq(X, X)` under `estimate_window_divergence`; the divergence-gating mechanism tested in `mdl_cost_model_gates_dup_by_sequence_divergence` (`src/egraph.rs:606-629`, at `d = 0.75`) is **mathematically unreachable** unless `--divergence` is manually supplied.
  3. In `compute_biedged_snarls` (`src/bin/horizon_resolve.rs:329-375`), the 2-edge-connected components (`comp`) are computed only to log `num_components` to `stderr` and are **discarded** before returning; `decompose_and_resolve` (`src/bin/horizon_resolve.rs:479-523`) instead slices walks greedily up to `max_bubble_steps = 16` between `is_anchor` segments.
- **Evidence:** `src/bin/horizon_resolve.rs:329-375, 456-477, 479-523`; `src/egraph.rs:398, 410-420, 606-629`; `ARCHITECTURE.md:28`.
- **Correction:** Update `ARCHITECTURE.md:28` (and `docs/PANGENOME_SUITE.md:25`, `manuscript/horizon_paper.md:50`) to state that `horizon_resolve` uses a **topological non-anchor/inversion proxy** clamped to `[0.01, 0.25]` (which always favors `dup` collapse over `seq(X, X)` for `len >= 4` bp unless `--divergence` is passed), and that `comp` from Tarjan bridge-finding is only logged, not used to partition walk windows.

---

## F6 -- Claim-list steering: `C1`–`C5` omit uncalibrated CRC/Karlin constants, `18%`-GC null false positives across other seeds, and `horizon_clean`'s ad-hoc score

- **Claim:** Round 6 claim list (`C1`–`C5`, `reports/rounds/round-6/round.json`)
- **Verdict:** **NARROWED** (meta-finding per brief)
- **Prior round context:** Round 4 (`L5-F7`) and Round 5 (`L3-F6`) previously flagged claim-list steering away from `n_cal = 15` CRC and `KARLIN_MDL_COEFF = 0.55`.
- **Concrete flaw:** The Round 6 claim list focuses on `P0`–`P4` feature inventory and runtime/coverage highlights while steering away from three inferential weaknesses still live in `manuscript/horizon_paper.md` and `docs/BENCHMARK_RESULTS.md`:
  1. `n_cal = 15` Conformal Risk Control and `KARLIN_MDL_COEFF = 0.55` remain uncalibrated fixed constants (`src/bin/horizon_filter.rs:17`, `manuscript/horizon_paper.md:40`).
  2. `manuscript/horizon_paper.md:10, 61` highlights `0` false-positive pairs on `seed=999` at `18%` GC, while `4/5` other `18%`-GC seeds still emit `80–448` pairs and `30%`/`45%`-GC single draws emit `288`/`149` pairs.
  3. `horizon_clean` (`manuscript/horizon_paper.md:18, 49`) uses an ad-hoc Bernoulli/BIC changepoint score (`score_junction_cp_bits`) rather than `src/mdl/dp.rs`, and `docs/BENCHMARK_RESULTS.md:77` reports only a `0`-severed negative-control run on the `10k`-segment `chr22` subgraph with no positive-control discordant-edge test on real HPRC data.
- **Evidence:** `manuscript/horizon_paper.md:10, 18, 40, 49, 61`; `docs/BENCHMARK_RESULTS.md:5, 13, 77`; `src/bin/horizon_filter.rs:17`.

---

## Verdict Table

| Claim ID | Verdict | Summary of Inference Attack |
| :--- | :--- | :--- |
| **C1** | **NARROWED** | `horizon_filter` never calls `bocpd_gc_profile` (uses global contig GC or fixed `0.5`); `horizon_wga --graph` omits `-log2(out_degree)` DAG branching priors. |
| **C2** | **REFUTED** | `0.47s` `horizon_filter --net` is sequence-free greedy length-netting (no SDUST/BOCPD/sequence MDL); `1.6% -> 81.2%` is a record-count denominator artifact. |
| **C3** | **NARROWED** | `p_dist` and `jc69_dist` are computed only on surviving aligned bp (`~63%` of avian locus, survivorship bias), `jc69_distance` clamps at `0.74`, and `jaccard_cov` is Dice coverage. |
| **C4** | **NARROWED** | `--resolve-dups` only emits `# EG:Z:` tags without collapsing `282.2M` pairs; `124.16 Mbp` on `GRCh38.chr6` vs `110.7 Mbp` `CHM13.chr6` q-arm reflects many-to-many repeat alignments; speedup compares across pre-fix `8cd5db0`. |
| **C5** | **REFUTED** | `estimate_window_divergence` is a sequence-free topological proxy clamped to `[0.01, 0.25]`, making `MdlCostModel` divergence gating unreachable (`1.207 < 2.0` bits/bp); Tarjan `comp` is discarded. |
