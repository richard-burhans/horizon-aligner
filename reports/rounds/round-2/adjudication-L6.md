# Round 2 Adjudication -- L6 (prior-art)

- **Attack sealed digest verified:** `9cc5696df060c1a838e2c0f7f2e67d3aeab660843fea0ab4698a9c901fff30ac` (`reports/rounds/round-2/attack-L6-prior-art.md`) matches `9cc5696df060`.
- **Prior rounds checked:** `./scripts/round history "" --kind adjudication` returned Round 1 (`reports/rounds/round-1/adjudication-L6.md`, `reports/rounds/round-1/completeness-critic.md`). Round 1 L6 already settled prior art for MDL alignment (C2), pangenome liftover (C11), and banded sequence-to-graph DP (C10).
- **Method cards:** `horizon-aligner` declares no `refs/methods/` directory; no statistical method card governs these prior-art and implementation-scope findings.

---

## Finding-by-Finding Rulings

### F1 -- C1 (`stream_gfa` implementation and GFA v1.0/v1.1 reader)
- **Ruling:** **UPHELD**
- **What survives:** The attacker's `SURVIVED` assessment holds: `src/gfa/mod.rs:98-184` and `219-318` implement buffered streaming for `.gfa`, multi-member `.gfa.gz` via `MultiGzDecoder`, and `"-"` `stdin`, with the cited unit tests, and `docs/STRATEGIC_PLAN.md:13` cites `GBZ`/`GBZ-base` and `gfaidx` while disclaiming zero-copy. (The note that `src/gfa/mod.rs:167-175` defaults empty `L`-line orientations to `b'+'` is a spec-compliance detail for correctness/formats, not a prior-art defect in C1.)
- **Repair:** None required for C1.

### F2 -- C2 (Adaptive Pair-HMM parameter estimation prior art and `identity_prior` clamping)
- **Ruling:** **UPHELD**
- **What survives:** Both the prior-art omission (`last-train`, Hamada et al. 2017; Lunter 2007; Durbin et al. 1998 — all absent from `manuscript/references.bib`) and the scope overstatement in `docs/BENCHMARK_RESULTS.md:5` hold. `src/mdl/mod.rs:818-833` clamps `identity_prior` to `[0.60, 0.72]` (narrower than `HmmParams`'s `0.50..=0.98` doc comment at `src/mdl/dp.rs:11`), `gap_open_bits` to `[-5.5, -4.0]`, and `gap_extend_bits` to `[-1.8, -0.5]`, and the graph/GFA benchmarks (`horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_resolve` at `docs/BENCHMARK_RESULTS.md:54-64`) do not invoke the Pair-HMM.
- **Repair:** Endorsed (narrow `docs/BENCHMARK_RESULTS.md:5` to pairwise/DAG alignment results using the clamped `[0.60, 0.72]` LCB-estimated Pair-HMM alongside Allison et al. 1992 and Adams & MacKay 2007, and cite `last-train` in `manuscript/references.bib`). **The fix is free** (text and bibliography edit only; no re-measurement required).

### F3 -- C3 (Implemented binaries and `39` passing tests)
- **Ruling:** **UPHELD**
- **What survives:** The attacker's `UNCHALLENGED` assessment holds: `src/bin/` contains all nine binary files (`bench_chr1.rs`, `horizon_clean.rs`, `horizon_eval.rs`, `horizon_gfa.rs`, `horizon_guide.rs`, `horizon_lift.rs`, `horizon_merge.rs`, `horizon_resolve.rs`, `horizon_wga.rs` with `--graph`), and C3 makes no prior-art or novelty claim.
- **Repair:** None required under L6.

### F4 -- C4 (`horizon_lift` benchmark vs. `odgi position`, `halLiftover`, `vg surject`, `Liftoff`)
- **Ruling:** **NARROWED**
- **What survives:** As already ruled in Round 1 (`adjudication-L6.md:61-64`), pangenome coordinate liftover is established prior art in `odgi position`, `halLiftover`, `vg surject`, and `Liftoff`, and `odgi`, `vg`, and `Liftoff` remain absent from `manuscript/references.bib` alongside zero baseline comparison in `docs/BENCHMARK_RESULTS.md`.
- **Why narrowed:** C4 itself (`docs/BENCHMARK_RESULTS.md:21,57`) is a single-tool runtime/RSS/count measurement rather than a novelty or superiority claim, so the prior-art gap is an omission of citations and baselines in the surrounding document rather than a false assertion in C4's own sentence.
- **Repair:** Endorsed on the citation/scoping leg (add `odgi`/`vg`/`Liftoff` to `manuscript/references.bib` and label C4 as a single-tool measurement; **the text/citation fix is free**). Note that the attacker's additional request for a head-to-head runtime/RSS/concordance benchmark against `odgi position` or `vg surject` on the 500 `GRCh38.chr22` windows is **non-free** (requires running external tools on the 12 GB HPRC graph).

### F5 -- C5 (HPRC `chr22` pericentromere benchmark)
- **Ruling:** **UPHELD**
- **What survives:** The attacker's `UNCHALLENGED` assessment holds: C5 reports single-run runtime, RSS, aligned pairs, LCBs, Horizon Bubbles, and `horizon_clean` kept-link counts without a priority or novelty claim. (The parenthetical note that `src/bin/horizon_clean.rs:98-129` scores GC-composition changepoints rather than repeat-depth support and severed `0` of `6,030` links is verified.)
- **Repair:** None required under L6.

### F6 -- C6 (Synthetic E-Graph Resolver benchmark)
- **Ruling:** **UPHELD**
- **What survives:** The attacker's `UNCHALLENGED` assessment holds: C6 reports synthetic benchmark metrics without a priority or novelty claim. (The parenthetical note that anchor-bounded snarl/bubble decomposition is prior art in `vg` snarls and `smoothxg`, absent from `manuscript/references.bib`, is verified.)
- **Repair:** None required under L6.

### F7 -- C7 (Seeded DAG reachability + banded DP prior art and intra-node seeding)
- **Ruling:** **UPHELD**
- **What survives:** Extending Round 1 (`adjudication-L6.md:56-59`), `align_seq_to_dag_seeded` (`src/mdl/dag.rs:460-550`) uses `MASK_W10` (span 14, `src/seeding/mod.rs:10`) + `MASK_RY20` (span 25, `src/seeding/mod.rs:13`) intra-node seeding (`nd.seq` at lines 475–481) and forward/backward reachability closure (lines 490–518) before `align_seq_to_dag_mdl`. Seed-and-reachability subgraph isolation prior to sequence-to-DAG DP is established prior art in `GraphAligner`, `vg map`/`GCSA2`, `vg giraffe`, and `minigraph` (all absent from `manuscript/references.bib`), and intra-node seeding cannot seed on GFA segments shorter than 14 bp (falling back to full-graph DP when `seed_hits.is_empty()` at line 485).
- **Repair:** Endorsed (narrow `docs/PANGENOME_SUITE.md:18` to state that `align_seq_to_dag_seeded` adapts the seed-and-reachability subgraph-isolation pattern of `GraphAligner`/`vg`/`minigraph`, note the `<14 bp` intra-node seeding limitation, and cite those papers in `manuscript/references.bib`). **The fix is free**.

### F8 -- C8 (`ContigIndex` RAM bottleneck on `chr6`)
- **Ruling:** **UPHELD**
- **What survives:** The attacker's `SURVIVED` assessment holds: `src/seeding/mod.rs:37-130` allocates a per-position `HashMap<u64, Vec<usize>>` over the contig, directly explaining the `26.5 GB` RSS bottleneck reported in `docs/BENCHMARK_RESULTS.md:26` without minimizer/syncmer sampling (`minimap2`, Li 2018; Roberts et al. 2004).
- **Repair:** None required for C8's factual bottleneck report (citing `minimap2` / Roberts et al. 2004 is optional and free).

### F9 -- Meta-finding (Round 2 claim list and systemic prior-art/baseline gap)
- **Ruling:** **UPHELD**
- **What survives:** Six of the eight Round 2 claims (C1, C3, C4, C5, C6, C8) are narrow existence or single-run runtime/RSS statements, and `manuscript/references.bib` (11 entries) omits all of `GraphAligner`, `vg`, `vg giraffe`, `odgi`, `minigraph`, `PGGB`/`smoothxg`, `Minigraph-Cactus`, `minimap2`, `last-train`, and `Liftoff`, while `docs/BENCHMARK_RESULTS.md` includes no head-to-head baselines against them.
- **Repair:** Adding the missing prior-art citations to `manuscript/references.bib` is **free**; running head-to-head external tool benchmarks is **non-free**.

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `SURVIVED` assessment verified in `src/gfa/mod.rs:98-184,219-318` and `docs/STRATEGIC_PLAN.md:13`; no correction needed for C1 |
| F2 | UPHELD | `last-train` prior art omitted and `docs/BENCHMARK_RESULTS.md:5` overstates Pair-HMM scope across graph tools while clamping `identity_prior` to `[0.60, 0.72]`; text/citation fix is free |
| F3 | UPHELD | `UNCHALLENGED` assessment verified against `src/bin/`; C3 makes no prior-art assertion |
| F4 | NARROWED | `odgi position`/`halLiftover`/`vg surject`/`Liftoff` prior art and missing citations/baselines hold, but C4 itself is a single-run measurement; citation/scoping fix is free while head-to-head benchmarking is non-free |
| F5 | UPHELD | `UNCHALLENGED` assessment verified; C5 reports single-run metrics without a prior-art claim |
| F6 | UPHELD | `UNCHALLENGED` assessment verified; C6 reports synthetic metrics without a prior-art claim |
| F7 | UPHELD | Seed-and-reachability DAG isolation is prior art (`GraphAligner`/`vg`/`minigraph`, uncited) and intra-node seeding misses `<14 bp` segments (`src/mdl/dag.rs:475-481`); text/citation fix is free |
| F8 | UPHELD | `SURVIVED` assessment verified in `src/seeding/mod.rs:37-130`; `ContigIndex` per-position `HashMap` explains `26.5 GB` RSS |
| F9 | UPHELD | Verified `manuscript/references.bib` (11 entries) and `docs/BENCHMARK_RESULTS.md` omit core pangenome/aligner citations (free) and head-to-head baselines (non-free) |
