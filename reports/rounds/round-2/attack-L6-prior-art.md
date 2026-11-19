# Round 2 Attack -- prior-art

## Literature, Method Cards, and Prior-Round Search
- **Held literature & method cards located:** The brief declares no appendix for literature or method cards. Searching the repository confirmed that `horizon-aligner` has **no** `refs/catalog_index.json` and **no** method-card directory (`refs/methods/`). Held bibliography and manuscript sources were located at `manuscript/references.bib` (106 lines) and `manuscript/horizon_paper.md` (112 lines), alongside `docs/STRATEGIC_PLAN.md`, `docs/BENCHMARK_RESULTS.md`, `docs/PANGENOME_SUITE.md`, and `src/`.
- **Prior rounds checked:** `./scripts/round history ""` returned Round 1 (`reports/rounds/round-1/adjudication-L6.md` and `reports/rounds/round-1/completeness-critic.md`). Prior Round 1 L6 adjudications relevant to Round 2 are cited explicitly below (`round 1`'s C2, C10, C11).

---

## F1 -- C1: `stream_gfa` implementation matches claim, though empty `L`-line orientation defaulting departs from GFA v1.0/v1.1 spec
- **Claim:** `C1` (`docs/STRATEGIC_PLAN.md:13`)
- **Verdict:** **SURVIVED**
- **Evidence:** Inspection of `src/gfa/mod.rs:98-184` and `219-318` confirms `stream_gfa` provides a buffered line-streaming GFA v1.0/v1.1 reader supporting `.gfa`, multi-member BGZF `.gfa.gz` via `MultiGzDecoder` (`src/gfa/mod.rs:122`), and `"-"` `stdin` (`src/gfa/mod.rs:113-116`), with unit tests for empty `L`-line orientations (`streams_all_record_types_and_handles_empty_link_orients`, lines 219-276) and multi-member `.gfa.gz` (`streams_gzipped_gfa_file`, lines 279-318). Furthermore, `docs/STRATEGIC_PLAN.md:13` cites `GBZ` / `GBZ-base` (Sirén et al.) and `gfaidx` and explicitly disclaims zero-copy. (Note: `src/gfa/mod.rs:167-175` silently defaults empty `L`-line orientations to `b'+'`, which departs from the GFA v1.0/v1.1 specification where `FromOrient`/`ToOrient` are mandatory `[+-]` fields.)
- **Correction:** None required for C1's factual existence claim.

## F2 -- C2: Adaptive Pair-HMM parameter estimation from anchors is prior art (`last-train`, Lunter 2007), and `identity_prior` is clamped to `[0.60, 0.72]`
- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:5`)
- **Verdict:** **NARROWED**
- **Flaw:**
  1. Estimating divergence-adaptive Pair-HMM match identity and affine gap open/extend probabilities from initial anchor alignments is established prior art in **`last-train`** (Hamada et al. 2017, *Bioinformatics* 33(6):926–928), **Lunter (2007)** (*Bioinformatics* 23(19):2485–2487), and **Durbin et al. (1998)** (*Biological Sequence Analysis*, Ch. 4), while the three-state MDL Pair-HMM and BOCPD are from **Allison, Wallace & Yee (1992)** and **Adams & MacKay (2007)** (already cited in `manuscript/references.bib:1-17` and adjudicated in `round 1`'s C2 at `adjudication-L6.md:16-19`, but omitted from `docs/BENCHMARK_RESULTS.md`). None of `last-train`, `Lunter (2007)`, or `Durbin et al. (1998)` appear in `manuscript/references.bib`.
  2. Additionally, `estimate_hmm_params` (`src/mdl/mod.rs:818-833`) clamps `identity_prior` to `[0.60, 0.72]` (narrower than `HmmParams`'s documented `0.50..=0.98` at `src/mdl/dp.rs:11`), `gap_open_bits` to `[-5.5, -4.0]`, and `gap_extend_bits` to `[-1.8, -0.5]`, and the graph/GFA benchmarks in `docs/BENCHMARK_RESULTS.md:54-64` (`horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_resolve`) do not invoke the Pair-HMM (`src/mdl/dp.rs`, `src/mdl/dag.rs`) at all.
- **Correction:** Narrow `docs/BENCHMARK_RESULTS.md:5` to state that pairwise/DAG alignment results use the clamped LCB-estimated Pair-HMM (`[0.60, 0.72]`) alongside Allison et al. (1992) MDL and Adams & MacKay (2007) BOCPD, and cite `last-train` (Hamada et al. 2017) as prior art for anchor-based Pair-HMM parameter estimation.

## F3 -- C3: Implementation and test-count claim makes no prior-art assertion
- **Claim:** `C3` (`docs/STRATEGIC_PLAN.md:3`)
- **Verdict:** **UNCHALLENGED**
- **Evidence:** `src/bin/` contains all nine named binary files (`bench_chr1.rs`, `horizon_clean.rs`, `horizon_eval.rs`, `horizon_gfa.rs`, `horizon_guide.rs`, `horizon_lift.rs`, `horizon_merge.rs`, `horizon_resolve.rs`, `horizon_wga.rs` with `--graph`). No prior-art or novelty claim is made in C3; runtime test-count verification belongs to the reproduction lens.

## F4 -- C4: `horizon_lift` benchmark omits established pangenome liftover prior art (`odgi position`, `halLiftover`, `vg surject`, `Liftoff`)
- **Claim:** `C4` (`docs/BENCHMARK_RESULTS.md:21,57`)
- **Verdict:** **NARROWED**
- **Flaw:** As already adjudicated in `round 1`'s C11 (`adjudication-L6.md:61-64`), coordinate liftover across pangenome paths is established prior art in **`odgi position`** (Guarracino et al. 2022, *Bioinformatics* 38(13):3319–3326), **`halLiftover`** (Hickey et al. 2013, *Bioinformatics* 29(10):1341–1342), **`vg surject`** (Garrison et al. 2018, *Nat. Biotechnol.* 36:875–879), and **`Liftoff`** (Shumate & Salzberg 2021, *Bioinformatics* 37(12):1639–1643). Round 2 adds the `horizon_lift` implementation (`src/bin/horizon_lift.rs`) and benchmark (`256.60s`, `619.2 MB` RSS for 500 `chr22` 1-kb windows), yet `docs/BENCHMARK_RESULTS.md` reports zero comparison against `odgi position`, `vg surject`, or `halLiftover`, and `odgi`, `vg`, and `Liftoff` remain absent from `manuscript/references.bib`.
- **Correction:** Narrow C4's presentation in `docs/BENCHMARK_RESULTS.md:21,57` by adding a baseline comparison (runtime, RSS, and lifted-interval concordance) against `odgi position` or `vg surject` on the same 500 `GRCh38.chr22` windows, and add `odgi`/`vg`/`Liftoff` to `manuscript/references.bib`.

## F5 -- C5: HPRC `chr22` pericentromere benchmark makes no prior-art assertion
- **Claim:** `C5` (`docs/BENCHMARK_RESULTS.md:23,59,63`)
- **Verdict:** **UNCHALLENGED**
- **Evidence:** C5 reports wall-clock, RSS, aligned-pair, and kept-link counts from a single run without making a priority or novelty claim. (Note: `docs/BENCHMARK_RESULTS.md` lacks precision/recall or edge-pruning comparisons against `minimap2`, `wfmash`, `vg clip`, or `smoothxg`, and `horizon_clean` severed `0` of `6,030` links because `src/bin/horizon_clean.rs:98` only scores GC-composition changepoints rather than repeat-depth support.)

## F6 -- C6: Synthetic E-Graph Resolver benchmark makes no prior-art assertion
- **Claim:** `C6` (`docs/BENCHMARK_RESULTS.md:64`)
- **Verdict:** **UNCHALLENGED**
- **Evidence:** C6 reports synthetic benchmark metrics (`200` segments, `20` walks x `225` steps, `24.0 MB`, `4.48s`) without making a priority or novelty claim. (Note: anchor-bounded snarl/bubble decomposition is prior art in `vg` snarls, Paten et al. 2018, and `smoothxg`, Garrison et al. 2024, neither cited in `manuscript/references.bib`.)

## F7 -- C7: Seeded DAG reachability + banded DP is prior art (`GraphAligner`, `vg map`/`GCSA2`, `vg giraffe`, `minigraph`), and node-only seeding misses junction k-mers
- **Claim:** `C7` (`docs/PANGENOME_SUITE.md:18`)
- **Verdict:** **NARROWED**
- **Flaw:** Extending `round 1`'s C10 (`adjudication-L6.md:56-59`, which established banded sequence-to-graph DP as prior art in `GraphAligner`, `minigraph`, and `vg giraffe`), `align_seq_to_dag_seeded` (`src/mdl/dag.rs:460-550`) adds spaced-seed (`MASK_W10`, `MASK_RY20`) hitting plus forward/backward graph reachability closure to isolate a subgraph before running `align_seq_to_dag_mdl`. Seed-and-reachability subgraph isolation prior to sequence-to-DAG DP is itself established prior art in **`GraphAligner`** (Rautiainen & Marschall 2020, *Genome Biology* 21:253), **`vg map` / `GCSA2`** (Garrison et al. 2018; Sirén 2017), **`vg giraffe`** (Sirén et al. 2021, *Science* 374:abg8871), and **`minigraph`** (Li et al. 2020, *Genome Biology* 21:265). Furthermore, unlike `GraphAligner` and `GCSA2` which index k-mers *across graph edges*, `src/mdl/dag.rs:475-479` only indexes within individual `GraphNode` (`S`-line) sequences, so `MASK_W10` (span 14) and `MASK_RY20` (span 25) cannot seed on GFA segments shorter than 14 bp (already flagged in `round 1` `completeness-critic.md:83-85`). None of `GraphAligner`, `vg`, `vg giraffe`, `minigraph`, or `GCSA2` are cited in `manuscript/references.bib`.
- **Correction:** Narrow `docs/PANGENOME_SUITE.md:18` to state that `align_seq_to_dag_seeded` adapts the standard seed-and-reachability subgraph-isolation pattern of `GraphAligner` / `vg` / `minigraph` using `MASK_W10` + `MASK_RY20` intra-node seeds, note the `<14 bp` segment limitation, and cite those papers in `manuscript/references.bib`.

## F8 -- C8: `ContigIndex` RAM bottleneck (`26.5 GB` on `chr6`) is confirmed and stems from omitting minimizer/syncmer sampling
- **Claim:** `C8` (`docs/BENCHMARK_RESULTS.md:26`)
- **Verdict:** **SURVIVED**
- **Evidence:** `src/seeding/mod.rs:37-120` confirms `canonical_spaced_seeds` allocates a dense per-position `HashMap<u64, Vec<usize>>` over the contig, directly explaining the reported `26.5 GB` RSS (~0.29 GB/Mbp) on `GRCh38.chr6` vs `CHM13.chr6`. From a prior-art standpoint, this scaling limit is a consequence of not using established minimizer sampling (Roberts et al. 2004, *Bioinformatics* 20(18):3363–3369; **`minimap2`**, Li 2018, *Bioinformatics* 34(18):3094–3100) or sorted flat seed arrays, neither of which is cited in `manuscript/references.bib` or `docs/BENCHMARK_RESULTS.md`.
- **Correction:** None required for C8's factual bottleneck report; consider citing `minimap2` (Li 2018) / Roberts et al. (2004) when explaining why per-position `HashMap` indexing hits 26.5 GB on a single chromosome.

## F9 -- Meta-finding: Round 2 claim list steers away from missing pangenome prior-art citations and baselines
- **Claim:** Claim-list selection (`BRIEF-prior-art.md:118-136`)
- **Verdict:** **NARROWED** *(meta-finding across C1–C8)*
- **Flaw:** Six of the eight claims in Round 2 (C1, C3, C4, C5, C6, C8) are narrow implementation-existence or single-run runtime/RSS statements, which steers attention away from the systemic prior-art gap in `docs/BENCHMARK_RESULTS.md`, `docs/PANGENOME_SUITE.md`, and `manuscript/references.bib`: none of the core pangenome and pairwise-alignment prior art (`GraphAligner`, `vg`, `vg giraffe`, `odgi`, `minigraph`, `PGGB`/`smoothxg`, `Minigraph-Cactus`, `minimap2`, `last-train`, `Liftoff`) is cited in `manuscript/references.bib` or benchmarked head-to-head in `docs/BENCHMARK_RESULTS.md`.

---

| C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SURVIVED | NARROWED | UNCHALLENGED | NARROWED | UNCHALLENGED | UNCHALLENGED | NARROWED | SURVIVED |
