# Round 2 ADJUDICATION of L2 (Sampling Lens) -- Horizon Pangenome Suite

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner`
**Attack Reviewed:** `reports/rounds/round-2/attack-L2-sampling.md` (sha256 prefix `512b076b8229` verified)
**Prior Round & Method-Card Check:** `scripts/round history "horizon" --kind adjudication` reviewed (Round 1 left `R1-IMPL`, `R1-GFA`, and `R1-HMM` open; Round 1 completeness critic flagged `MASK_RY20` and `HB:i:1` `*`-sequence emission, but none of the Round 2 sampling defects in the newly implemented binaries/benchmarks were previously adjudicated). No `refs/methods/` directory exists in this repository.

---

## Finding-by-Finding Rulings

### F1 -- C2: `estimate_hmm_params` uses an unstratified prefix sample, ascertainment-biased anchors, and is bypassed by `horizon_wga --graph`
- **Verdict:** **NARROWED**
- **What survives:** `estimate_hmm_params` (`src/mdl/mod.rs:773-809`) iterates `lcbs_fwd` before `lcbs_rc` with `lcb.iter().take(64)`, `lcb.windows(2).take(64)`, and `break` at `total_bases >= 4000` (censoring `dx, dy > 500` and clamping `identity_prior` to `[0.60, 0.72]`), and `horizon_wga --graph` (`src/bin/horizon_wga.rs:195`) passes `&HmmParams::default()` rather than adaptive parameters.
- **Premise check:** The attacker's fourth leg (homologous-block-only MAF extraction in `scripts/benchmarks/extract_fasta_from_maf.py:28-45`) is **overturned** as an undisclosed defect because `docs/BENCHMARK_RESULTS.md:5` already states *"and MAF-extracted FASTAs sample homologous blocks."* C2 is therefore **NARROWED** (not `REFUTED` in full).
- **Repair:** Split. Qualifying `docs/BENCHMARK_RESULTS.md:5` to state that `horizon_wga --graph` uses `HmmParams::default()` and that `estimate_hmm_params` samples up to 4,000 anchor bases (`take(64)` per LCB, `dx, dy <= 500`) is **free** (doc edit only). Replacing prefix `take(64)`/`break` with stratified LCB/strand sampling is **non-free** (code change + full WGA/HPRC re-benchmarking).

### F2 -- C3: Four of nine binaries have zero unit tests and the five tested binaries sample one synthetic fixture each
- **Verdict:** **UPHELD**
- **What survives:** `cargo test` executes **40** passing tests (35 in `src/lib.rs`, 1 each in `horizon_clean`, `horizon_gfa`, `horizon_lift`, `horizon_resolve`, `horizon_wga`) plus 1 `#[ignore]` test (`src/hal/mod.rs:574`), not 39, and 4 of 9 binaries (`horizon_merge`, `horizon_guide`, `horizon_eval`, `bench_chr1`) have 0 unit tests.
- **Repair:** Endorsed. Updating `docs/STRATEGIC_PLAN.md:3` to report `40 passing tests (1 ignored)` and clarify that `horizon_merge`, `horizon_guide`, `horizon_eval`, and `bench_chr1` have 0 binary unit tests is **free** (doc edit only).

### F3 -- C4: `horizon_lift` 500-window sample excludes the entire `0–16.0 Mbp` `chr22` p-arm/centromere and all `-`-strand intervals
- **Verdict:** **UPHELD**
- **What survives:** Inspecting the entire `/tmp/grch38_chr22_windows.bed` (500 lines) confirms all 500 intervals are `+` strand spanning `16,000,000`–`49,933,000` with uniform `68,000` bp stride, omitting `0–16.0 Mbp` and all `-`-strand intervals.
- **Repair:** Endorsed. Narrowing C4 in `docs/BENCHMARK_RESULTS.md:21,57` to `+`-strand euchromatic q-arm windows (`GRCh38.chr22:16.0M–49.9M`) is **free** (doc edit only).

### F4 -- C5: `horizon_clean` benchmarked on sequence-less `*`-GFA (`0` links severed) and pre-paired 5 Mbp pericentromere slices
- **Verdict:** **NARROWED**
- **What survives:** All `6,031` `S` lines in `/tmp/chr22_peri_5m.gfa` (and all `2,777` `S` lines in `/tmp/chr6_peri_5m.gfa`) carry `*` placeholder sequences, so `src/bin/horizon_clean.rs:99,155` sets `len = 0 < 12` and returns `0.0` without ever executing `bocpd_gc_profile`.
- **Premise check:** The 5 Mbp slice scope (`GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M`) is already explicitly stated in C5, so only the undisclosed `*`-sequence / `0`-scorable-junctions defect survives.
- **Repair:** Split. Disclosing in `docs/BENCHMARK_RESULTS.md:23,24,63` that `/tmp/chr22_peri_5m.gfa` and `/tmp/chr6_peri_5m.gfa` are `*`-sequence GFAs with `0` junctions scorable by `bocpd_gc_profile` is **free** (doc edit only). Re-benchmarking `horizon_clean` on a sequence-bearing HPRC GFA subgraph is **non-free** (requires benchmark execution).

### F5 -- C6: `horizon_resolve` benchmark uses a synthetic 200-segment `ACGT`-repeat fixture, not HPRC
- **Verdict:** **UPHELD**
- **What survives:** `/tmp/bench_resolve_10k.gfa` consists of `200` identical-length (`200` bp) `ACGT`-repeat segments and `20` identical `225`-step walks, yet is listed under `## Pangenome Graph Benchmarks (12 GB hprc-v1.0-mc-grch38.gfa.gz & Subgraphs)` (`docs/BENCHMARK_RESULTS.md:50,64`).
- **Repair:** Endorsed. Labeling C6 at `docs/BENCHMARK_RESULTS.md:64` explicitly as a synthetic `ACGT`-repeat benchmark (`200` segs, `20` identical walks) rather than an HPRC subgraph benchmark is **free** (doc edit only).

### F6 -- C7: `align_seq_to_dag_seeded` drops seeds in short nodes (`<14`/`<25` bp) and high-copy/low-complexity windows, and over-expands reachability
- **Verdict:** **NARROWED**
- **What survives:** `align_seq_to_dag_seeded` (`src/mdl/dag.rs:466-518`) seeds per-node (missing `<14`/`<25` bp nodes and `L`-line junctions), drops `>8`-occurrence and low-complexity seeds (`src/seeding/mod.rs:20-30,102,136`), expands backward to all ancestors and forward to all descendants rather than a bounded neighborhood, and is bypassed completely by the 4-node benchmark (`nodes.len() <= 8` at `src/mdl/dag.rs:466`).
- **Premise check:** C7 (`docs/PANGENOME_SUITE.md:18`) accurately names the implemented primitives (`MASK_W10`, `MASK_RY20`, forward/backward reachability closure, `align_seq_to_dag_mdl`); what fails is the claim that this isolates relevant local components in general pangenome DAGs and that the benchmark exercises it.
- **Repair:** Split. Qualifying `docs/PANGENOME_SUITE.md:18` and `docs/BENCHMARK_RESULTS.md:65` to note per-node `<14`/`<25` bp/junction blind spots, unbounded ancestor+descendant closure, and the `<= 8`-node bypass in the 4-node benchmark is **free** (doc edit only). Implementing junction-spanning seeding with bounded closure and benchmarking on a `> 8`-node graph is **non-free** (code + benchmark work).

### F7 -- C1 & C8: Checked and survived under the sampling lens
- **Verdict:** **UPHELD**
- **What survives:** `stream_gfa` (`src/gfa/mod.rs:98-159,256-318`) streams all GFA records without subsampling and passes `.gfa.gz` and empty-`L` unit tests (C1), and `ContigIndex` (`src/seeding/mod.rs:100-137,329-354`) indexes every position in `HashMap<u64, Vec<usize>>` without minimizer/stride sampling (C8).
- **Repair:** None required.

### F8 -- Framing bias: Claim list omits homologous-block MAF filtering and pending B10K data
- **Verdict:** **NARROWED**
- **What survives:** The runner-authored claim list truncated C2's second sentence (*"and MAF-extracted FASTAs sample homologous blocks"*) and omitted `B10K 363-way` (`Pending Data` at `docs/BENCHMARK_RESULTS.md:48`) as well as the `*`-GFA / synthetic-fixture context of C5/C6.
- **Premise check:** `docs/BENCHMARK_RESULTS.md:5` itself *does* disclose homologous-block MAF sampling, so the framing omission is in the `round.json` / `BRIEF-sampling.md` claim excerpt rather than the underlying document.
- **Repair:** Covered by F1, F4, and F5; no additional document edit required.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | `estimate_hmm_params` prefix-4000-bp sampling, `dx,dy<=500` censoring, and `--graph` `HmmParams::default()` hold (doc fix is free; code fix is non-free), but MAF homologous-block sampling was already disclosed in `docs/BENCHMARK_RESULTS.md:5` |
| F2 | UPHELD | `cargo test` runs 40 passing tests (1 ignored), not 39, and 4/9 binaries have 0 unit tests; doc fix is free |
| F3 | UPHELD | All 500 windows in `/tmp/grch38_chr22_windows.bed` are `+` strand in `16.0M–49.933M`, omitting `0–16.0 Mbp` and `-` strand; doc fix is free |
| F4 | NARROWED | `/tmp/chr22_peri_5m.gfa` and `/tmp/chr6_peri_5m.gfa` have 0 non-`*` sequences so `horizon_clean` never runs BOCPD (doc fix is free; re-benchmarking is non-free), while the 5 Mbp slice scope was already stated in C5 |
| F5 | UPHELD | `/tmp/bench_resolve_10k.gfa` is a synthetic 200-segment `ACGT`-repeat fixture with 20 identical walks placed under the HPRC heading; doc fix is free |
| F6 | NARROWED | Per-node `<14/25` bp/junction blind spots, unbounded ancestor+descendant closure, and the `<=8`-node benchmark bypass hold (doc fix is free; code + benchmark fix is non-free), though C7's named primitives exist |
| F7 | UPHELD | C1 (`stream_gfa`) and C8 (`ContigIndex` dense per-position `HashMap`) survive verification under the sampling lens |
| F8 | NARROWED | Runner's claim excerpt truncated C2's MAF caveat and omitted pending B10K data, though `docs/BENCHMARK_RESULTS.md:5` itself disclosed MAF block filtering |
