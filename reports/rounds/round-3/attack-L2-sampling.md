# Round 3 Attack -- L2 (Sampling Lens)

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner`  
**Prior Round & Method-Card Check:** `scripts/round history "sampling" --kind adjudication`, `scripts/round history "cactus241way"`, and `scripts/round history "FlatSeedIndex"` reviewed. No `refs/methods/` directory exists in this repository. Prior Round 2 adjudications relevant to this lens (`round 2's F1`, `round 2's F2`, `round 2's F7`, `round 2's F8`) are explicitly disclosed below along with what changed in Round 3.

---

## F1 -- C2 REFUTED: `cactus241way chr22 1MB` samples only 15% homologous blocks, bifurcates seeding at `N_x * N_y < 10^9`, and leaves `estimate_hmm_params` prefix-biased

- **Claim:** C2 (`docs/BENCHMARK_RESULTS.md:16,47`) — *"UCSC cactus241way `chr22 1MB`: **62.3%** recall / **60.3%** precision (Dog, across 129 target scaffolds; **61.3% F1** in **17.16s**)."*
- **Verdict:** **REFUTED**
- **Prior-round disclosure:** `round 2's F1` (NARROWED) already established that `estimate_hmm_params` (`src/mdl/mod.rs:732-816`) samples only the first `4,000` anchor bases (`lcbs_fwd` before `lcbs_rc`, `lcb.iter().take(64)`), computes gap parameters from `lcb.windows(2).take(64)` with `dx <= 500 && dy <= 500` right-censoring, and clamps `identity_prior` to `[0.60, 0.72]`; `round 2's F8` noted that `docs/BENCHMARK_RESULTS.md:5` already mentions *"MAF-extracted FASTAs sample homologous blocks"*. **What changed in Round 3:** `cactus241way chr22 1MB` was re-run (`63.6% -> 62.3%` recall, `61.0% -> 60.3%` precision), `horizon_wga --graph` was wired to `estimate_dag_hmm_params`, **but** `estimate_hmm_params` in `src/mdl/mod.rs:732-816` was left unchanged and its prefix/censoring sampling remains undisclosed in `docs/BENCHMARK_RESULTS.md:5`.
- **New sampling defects:**
  1. **85% non-homologous sequence stripped & `+`-strand/case flattening:** Inspecting `/tmp/hg38_zoo_1m.fa` and `/tmp/dog_zoo_1m.fa` (produced by `scripts/benchmarks/extract_fasta_from_maf.py:28-45`) shows that the "`chr22 1MB`" reference FASTA is only **`154,483` bp** (1 contig, ~15% of 1 Mbp) and the target FASTA is **`132,089` bp** across **`31` contigs** (min `1` bp, max `30,551` bp; not `129` scaffolds in the workspace artifact), all concatenated onto `+` strand and uppercased (`.upper()`). Stripping ~85% of non-homologous/repetitive sequence inflates precision by removing decoy sequence where false positives occur. Moreover, `cactus241way` is absent from `scripts/benchmarks/run_all_evals.py:56-76`.
  2. **Scale-gated seeding discontinuity (`src/mdl/mod.rs:289,347`):** `align_contig_pair_with_config` gates 2-transition `LASTZ12`, `MASK_RY20`, and `w=10` 1-transition seeding behind `(seq_x.len() as u128) * (seq_y.len() as u128) < 1_000_000_000`. Against `154,483` bp of `hg38`, any dog scaffold `> 6,473` bp silently skips 3 of the 5 fallback seeding tiers, whereas scaffolds `<= 6,473` bp run them. On an actual 1 Mbp $\times$ 1 Mbp contig, all scaffolds $\ge 1,000$ bp skip those tiers.
- **Correction:** Re-label C2 in `docs/BENCHMARK_RESULTS.md:16,47` as a `154.5 kb` homologous-block MAF extract (not `1MB`), disclose the `4,000`-bp prefix / `take(64)` / `dx,dy<=500` sampling in `estimate_hmm_params` and the `10^9` cell-product seeding cutoff in `docs/BENCHMARK_RESULTS.md:5`, and benchmark on unconcatenated 1 Mbp genomic sequence with softmasking intact.

---

## F2 -- C3 NARROWED: All 9 binaries still sample only `n = 1` tiny synthetic fixture per binary (and `cargo test` runs `46` passing, not `45`)

- **Claim:** C3 (`docs/STRATEGIC_PLAN.md:3`) — *"All target pangenome binaries ... are now implemented in `src/bin/` with unit tests in every binary (`45` passing tests + `1` ignored; `cargo clippy --all-targets -- -D warnings` clean)."*
- **Verdict:** **NARROWED**
- **Prior-round disclosure:** `round 2's F2` (UPHELD) found that 4 of 9 binaries had 0 unit tests and the 5 tested binaries sampled 1 synthetic fixture each (`40` passing + `1` ignored). **What changed in Round 3:** Commit `488b020` added 1 unit test to each of the 4 previously untested binaries (`bench_chr1`, `horizon_eval`, `horizon_guide`, `horizon_merge`).
- **Evidence:** Running `cargo test` in `/home/rcb112/project/gemini/bench/horizon-aligner` executes **`46` passing tests + `1` ignored** (`37` passing + `1` ignored in `src/lib.rs`, plus `9` passing across `src/bin/*.rs`), not `45` passing. Under the sampling lens, every single one of the 9 binaries in `src/bin/` (`bench_chr1`, `horizon_clean`, `horizon_eval`, `horizon_gfa`, `horizon_guide`, `horizon_lift`, `horizon_merge`, `horizon_resolve`, `horizon_wga`) still has **`n = 1` unit test** sampling a single tiny synthetic fixture (e.g., 2–4 nodes, `<100` bp sequences), with zero coverage of real HPRC/UCSC inputs, multi-locus edge cases, or failure modes.
- **Correction:** Update `docs/STRATEGIC_PLAN.md:3` to state `46 passing tests + 1 ignored` and qualify that each binary's unit test covers a single minimal synthetic fixture (`n = 1` per binary).

---

## F3 -- C4 REFUTED: `FlatSeedIndex` censors all `>= 3`-copy seeds on `>= 100 kb` contigs and C4 samples sub-chromosome slices (`20 Mbp`, `5 Mbp`) instead of full `chr6` (`170.8 Mbp`)

- **Claim:** C4 (`docs/STRATEGIC_PLAN.md:11`) — *"`ContigIndex` (`src/seeding/mod.rs`) now uses a compact `FlatSeedIndex` (`12 bytes/entry` sorted `Vec<u64>` + `Vec<u32>` with on-the-fly reverse-complement coordinates, cutting `20 Mbp chr6 q-arm` peak RSS from `5,863.3 MB` to `2,409.7 MB` and `chr6` 5 Mbp pericentromere RSS from `3,084.5 MB` to `1,048.2 MB`)."*
- **Verdict:** **REFUTED**
- **Evidence:**
  1. **Aggressive occurrence censoring (`max_occ = 2`):** `prepare_contig_index_opts` (`src/seeding/mod.rs:567`) sets `let max_occ = if seq.len() >= 100_000 { 2 } else { 4 };`, and `FlatSeedIndex::build` (`src/seeding/mod.rs:433`) drops every seed key with `j - i > max_occ` (in addition to `is_low_complexity_window` at `src/seeding/mod.rs:390`). On any contig $\ge 100\text{ kb}$ (including `5 Mbp` pericentromere and `20 Mbp` q-arm), every spaced seed occurring $\ge 3$ times is excluded from the index—silently shrinking RSS by discarding repetitive/paralogous seeds without any precision/recall evaluation.
  2. **Sub-chromosome slice sampling instead of full `>100 Mbp` chromosome (`R2-INDEX-RAM`):** `R2-INDEX-RAM` (`BRIEF-sampling.md:134`) required bringing full `>100 Mbp` chromosome alignment well below `26.5 GB` RSS. C4 instead samples only a `20 Mbp` q-arm slice (`11.7%` of `170.8 Mbp` `GRCh38.chr6`) and a `5 Mbp` pericentromere slice (`2.9%` of `GRCh38.chr6`), both from a single chromosome (`chr6`) of a single pair (`GRCh38` vs `CHM13`), with no truth oracle. Meanwhile, the full `GRCh38.chr6` (`170.8 Mbp`) vs `CHM13.chr6` (`110.7 Mbp`) run did not complete (`/tmp/hprc_chr6.gfa` and `/tmp/hprc_chr6.maf` are `0` bytes; `docs/BENCHMARK_RESULTS.md:26` reports only an incomplete mid-run check at `13:59` of `10.09 GB` RSS).
  3. **Fallback `HashMap`s remain:** `src/mdl/mod.rs:248-261,312-325` still allocates up to four `HashMap<u64, Vec<usize>>`s (`MASK_LASTZ12` and `MASK_LASTZ11` fwd/rc) whenever `primary_cov < min_fallback_cov`.
- **Correction:** Qualify `docs/STRATEGIC_PLAN.md:11` to state that `FlatSeedIndex` drops all seeds with `occ > 2` on contigs `>= 100 kb`, that `LASTZ12`/`LASTZ11` fallbacks still use `HashMap`s, and that `20 Mbp` q-arm / `5 Mbp` pericentromere are sub-chromosome slices (full `170.8 Mbp chr6` remains incomplete).

---

## F4 -- C5 REFUTED: `chr22` Pericentromere 5 Mbp samples a pre-paired homologous window with `max_occ = 2` seed censoring, no truth oracle, and no baselines

- **Claim:** C5 (`docs/BENCHMARK_RESULTS.md:23,59`) — *"Real HPRC `chr22` Pericentromere 5 Mbp (`GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M`): **14.57s**, **1,861.9 MB** peak RSS, **5,106,541** many-to-many `(x, y)` aligned pairs (`unique_x_bp=4,387,655`, `unique_y_bp=4,549,858`), `16,083` Horizon Bubbles."*
- **Verdict:** **REFUTED**
- **Evidence:**
  1. **Pre-paired single-locus windowing:** C5 samples a single manually paired 5 Mbp window (`10%` of `50.8 Mbp` `GRCh38.chr22`, `14%` of `35.6 Mbp` `CHM13.chr22`), excluding the unalignable `0–15.5 Mbp` p-arm/satellite array and the remaining q-arm. Full `GRCh38.chr22` vs `CHM13.chr22` was not re-run on the current build (`docs/BENCHMARK_RESULTS.md:22,58` notes it is `pre-0e70693`).
  2. **Zero truth oracle and zero head-to-head baselines:** Unlike `R2-RESOLVE-LIFT` (`CARRIED.md:11`), which called for head-to-head pangenome baselines (`minimap2`/`odgi`/`vg`), C5 reports `5,106,541` pairs and `unique_x_bp`/`unique_y_bp` with no ground-truth pairs and no baseline comparison—leaving false-positive rate in a repetitive pericentromere completely unmeasured.
  3. **Silent algorithmic gating:** At `5 Mbp >= 100 kb`, `FlatSeedIndex` drops all seeds with `occ > 2` (`src/seeding/mod.rs:567`), and at $5\text{M} \times 5\text{M} = 2.5 \times 10^{13} \ge 10^9$, `src/mdl/mod.rs:289,347` skips 2-transition `LASTZ12`, `RY20`, and `w=10` 1-transition seeding.
- **Correction:** Label C5 in `docs/BENCHMARK_RESULTS.md:23,59` as an unvalidated single-locus `5 Mbp` pre-paired slice (no truth precision/recall or `minimap2`/`odgi`/`vg` baseline) and disclose `max_occ = 2` seed dropping.

---

## F5 -- C6 NARROWED: `horizon_gfa stats` streams all 81.4M segments without subsampling, but benchmarks `n = 1` graph/run with counter-only callbacks

- **Claim:** C6 (`docs/BENCHMARK_RESULTS.md:19,54`) — *"Full-graph streaming (`horizon_gfa stats`): **88.84s**, **97.1 MB** peak RSS across `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks."*
- **Verdict:** **NARROWED**
- **Evidence:** `cmd_stats` (`src/bin/horizon_gfa.rs:39-77`) calls `stream_gfa` and visits every `S`, `L`, `P`, and `W` record without subsampling. However, from the sampling lens, C6 is an **`n = 1` single-graph, single-run** measurement (no replicate variance, no other GFA graphs) whose callbacks (`|_| Ok(())` for `L`, `P`, `W`) only increment counters and never parse path/walk steps or retain sequence—whereas actual graph extraction (`horizon_gfa extract` on `chr6`, `docs/BENCHMARK_RESULTS.md:56`) requires `1,255.1 MB` RSS (`12.9x` higher) and `244.93s`.
- **Correction:** Qualify C6 in `docs/BENCHMARK_RESULTS.md:19,54` as an `n = 1` counter-only streaming pass over `hprc-v1.0-mc-grch38.gfa.gz` that does not reflect the memory/time profile of step-parsing or subgraph-extracting commands.

---

## F6 -- C1 SURVIVED: `stream_gfa` streams all GFA records without subsampling (confirming Round 2 L2 F7)

- **Claim:** C1 (`docs/STRATEGIC_PLAN.md:13`) — *"`stream_gfa` provides a buffered line-streaming GFA v1.0/v1.1 reader (supporting `.gfa`, multi-member BGZF `.gfa.gz` via `MultiGzDecoder`, and `"-"` `stdin`, with unit tests for empty `L`-line orientations and `.gfa.gz`)."*
- **Verdict:** **SURVIVED**
- **Prior-round disclosure:** `round 2's F7` (UPHELD) verified C1 under the sampling lens. No changes were made to `src/gfa/mod.rs` in Round 3.
- **Evidence:** Re-verified `src/gfa/mod.rs:97-219,229-329`: `stream_gfa` and `stream_gfa_reader` process all lines sequentially without sampling, and both unit tests (`streams_all_record_types_and_handles_empty_link_orients` and `streams_gzipped_gfa_file`) pass under `cargo test`.

---

## F7 -- Framing bias: Runner-authored claim list omits incomplete full-`chr6`, stale WGA rows, pending `B10K`, and no-op/synthetic pangenome benchmarks

- **Verdict:** **REFUTED** (framing)
- **Prior-round disclosure:** `round 2's F8` (NARROWED) flagged that the runner's Round 2 claim list truncated C2's MAF caveat and omitted `B10K 363-way` (`Pending Data`). **What changed in Round 3:** The runner's Round 3 claim list C1–C6 again omits unfavorable rows from `docs/BENCHMARK_RESULTS.md` and `CARRIED.md`:
  1. **Incomplete `R2-INDEX-RAM` full-chromosome target:** C4 highlights `20 Mbp` q-arm and `5 Mbp` pericentromere RSS while omitting that full `GRCh38.chr6` (`170.8 Mbp`) vs `CHM13.chr6` (`110.7 Mbp`) never finished (`/tmp/hprc_chr6.gfa` and `/tmp/hprc_chr6.maf` are `0` bytes; `docs/BENCHMARK_RESULTS.md:26` cites only a mid-run RSS check).
  2. **Stale pre-`FlatSeedIndex` WGA rows & `Pending Data`:** `docs/BENCHMARK_RESULTS.md:5,48` shows `multiz100way`, `cactus447way`, `multiz470way`, `cactus241way chr22_random`, and `GRCh38.chr22 vs CHM13.chr22` were **not re-run** on the current build (leaving `R2-HMM-XDROP`'s single-commit re-run incomplete) and `B10K 363-way` remains `Pending Data`. C2 selects only the single re-run WGA row (`cactus241way chr22 1MB`).
  3. **Dropped no-op/synthetic pangenome benchmarks & missing baselines:** `horizon_clean` (`0` links severed on `*`-seq GFA), `horizon_resolve` (synthetic `200`-seg `ACGT` repeat), `horizon_wga --graph` (`4`-node diamond `<=8` bypass), and the absence of `minimap2`/`odgi`/`vg` baselines (`R2-RESOLVE-LIFT`) are all omitted from C1–C6 while C3 presents all 9 binaries as implemented and unit-tested.

---

## Verdict Table

| Claim | Verdict | Primary Sampling Defect / Finding |
|---|---|---|
| C1 | SURVIVED | `stream_gfa` streams all GFA records without subsampling; Round 2 L2 F7 holds unchanged (`src/gfa/mod.rs:97-219,229-329`) |
| C2 | REFUTED | `chr22 1MB` extracts `154.5 kb` (+-strand, uppercased) homologous blocks, bifurcates fallback seeding at `N_x * N_y < 10^9` (`src/mdl/mod.rs:289,347`), and leaves `estimate_hmm_params` prefix-4000-bp / `dx,dy<=500` biased (`src/mdl/mod.rs:732-816`; `round 2's F1`) |
| C3 | NARROWED | `cargo test` runs `46` passing + `1` ignored (not `45`), and all 9 binaries still sample only `n = 1` tiny synthetic fixture each (`src/bin/*.rs`; `round 2's F2`) |
| C4 | REFUTED | `FlatSeedIndex` censors all `occ > 2` seeds on `>= 100 kb` contigs (`src/seeding/mod.rs:433,567`), `LASTZ12`/`11` still use `HashMap`s, and C4 samples `20 Mbp`/`5 Mbp` `chr6` slices while full `170.8 Mbp chr6` did not finish |
| C5 | REFUTED | Samples a single pre-paired `5 Mbp` pericentromere window with `max_occ = 2` seed censoring, zero truth oracle, and zero `minimap2`/`odgi`/`vg` baselines (`docs/BENCHMARK_RESULTS.md:23,59`) |
| C6 | NARROWED | `horizon_gfa stats` streams all 81.4M segments without subsampling, but is an `n = 1` counter-only pass (`src/bin/horizon_gfa.rs:39-77`) |
| Framing | REFUTED | Runner's claim list omits incomplete full-`chr6`, stale WGA rows, `B10K` pending data, no-op/synthetic pangenome rows, and missing external baselines |
