# Round 2 Attack -- Reproduction Lens (L7)

Target: `/home/rcb112/project/gemini/bench/horizon-aligner` at commit `b18f25b`.

- **Literature / Method Cards Check:** No appendix was declared in the brief. The repository contains `manuscript/` (no separate method cards directory exists). All reproduction checks below were executed directly against the workspace source, `cargo test`, `./target/release/*` binaries, git history, and cached `/tmp/` artifacts.
- **Prior Round Context:** Round 1 adjudicated `stream_gfa` and `horizon_clean` (`round 1's adjudication-L1.md`, `adjudication-L2.md`, `adjudication-L4.md`) when `stream_gfa` had zero callers/tests and `horizon_clean`/`horizon_resolve`/`horizon_lift` did not exist. Since Round 1, commits `b857377`, `7bf688c`, `0e70693`, and `b18f25b` implemented those binaries and wired `stream_gfa` with unit tests.

---

## F1 -- C2: "All results were generated out-of-the-box" is refuted by stale pre-change benchmarks and missing data/scripts

- **Claim:** C2 (`docs/BENCHMARK_RESULTS.md:5`) -- "All results were generated out-of-the-box using the divergence-adaptive Pair-HMM (`HmmParams`: `identity_prior`, `gap_open_bits`, and `gap_extend_bits` estimated from LCB anchors with Beta prior shrinkage in `src/mdl/mod.rs`), BOCPD compositional profiles (`src/bocpd.rs`), and normalized Pair-HMM transitions (`src/mdl/dp.rs`, `src/mdl/dag.rs`)."
- **Verdict:** **REFUTED**
- **Flaw & Evidence:**
  1. **Stale pre-change results:** Commit `0e70693` (`11:20:04`) activated the BOCPD jump gate in `src/mdl/dp.rs`, and commit `b18f25b` (`12:27:59`) made `gap_extend_bits` adaptive. `git diff 7bf688c..b18f25b -- docs/BENCHMARK_RESULTS.md` shows that `cactus241way (chr22 1MB)` was re-run and its output changed (`64.3%` recall / `3,362` bubbles -> `63.6%` recall / `3,411` bubbles), proving these code changes alter alignment results. However, at least 7 other alignment benchmarks in `docs/BENCHMARK_RESULTS.md` (`multiz100way` chrM Chimp/Mouse, `cactus447way` Horse/Macaque, `multiz470way` Marmoset, `cactus241way chr22_random`, and `GRCh38.chr22` vs `CHM13.chr22` `24,365,611` pairs) were **not** re-run and retain pre-`0e70693` / pre-`b18f25b` numbers.
  2. **Missing datasets and runner coverage:** `/home/rcb112/project/gemini/bench/data/benchmarks/` does not exist in the workspace, and `scripts/benchmarks/run_all_evals.py` omits `cactus241way`, all HPRC graph benchmarks, and all synthetic benchmarks (`hla_dup_5m.fa`, `bench_resolve_10k.gfa`), so none of those results can be reproduced out-of-the-box from the repository.
- **Correction:** Re-run all `docs/BENCHMARK_RESULTS.md` benchmarks on `HEAD` (`b18f25b`) or qualify the header at `docs/BENCHMARK_RESULTS.md:5` to identify which rows predate `0e70693`/`b18f25b` and document exact input-generation commands.

---

## F2 -- C3: `cargo test` runs 40 passing tests (not 39), and 4 of 9 binaries in `src/bin/` have 0 unit tests

- **Claim:** C3 (`docs/STRATEGIC_PLAN.md:3`) -- "All target pangenome binaries (`horizon_wga`, `horizon_wga --graph`, `horizon_merge`, `horizon_guide`, `horizon_eval`, `horizon_gfa`, `horizon_clean`, `horizon_resolve`, `horizon_lift`, `bench_chr1`) are now implemented in `src/bin/` and unit-tested (`39` passing tests)."
- **Verdict:** **REFUTED**
- **Flaw & Evidence:**
  Running `/home/agent/.cargo/bin/cargo test` on `HEAD` (`b18f25b`) produces:
  - `src/lib.rs`: `35 passed; 0 failed; 1 ignored` (`hal::tests::real_hal_decode_from_env`)
  - `src/bin/horizon_clean.rs`: `1 passed`
  - `src/bin/horizon_gfa.rs`: `1 passed`
  - `src/bin/horizon_lift.rs`: `1 passed`
  - `src/bin/horizon_resolve.rs`: `1 passed`
  - `src/bin/horizon_wga.rs`: `1 passed`
  - `src/bin/bench_chr1.rs`, `src/bin/horizon_eval.rs`, `src/bin/horizon_guide.rs`, `src/bin/horizon_merge.rs`: **`running 0 tests`** each.
  Total passing tests = **40** (with 1 ignored), not `39`, and **4 of the 9 binaries** in `src/bin/` have zero unit tests.
- **Correction:** Update `docs/STRATEGIC_PLAN.md:3` to state `40 passing tests (1 ignored)` and note that `bench_chr1`, `horizon_eval`, `horizon_guide`, and `horizon_merge` currently have `0` unit tests.

---

## F3 -- C4: `horizon_lift` counts match `/tmp/` cache, but run predates `0e70693` and lacks repo inputs

- **Claim:** C4 (`docs/BENCHMARK_RESULTS.md:21,57`) -- "Two-pass graph liftover (`horizon_lift --bridge-bp 10` on 500 `GRCh38.chr22` 1-kb windows -> `CHM13.chr22`): **256.60s**, **619.2 MB** peak RSS (`469` exact, `1` partial, `19` split, `11` unmapped)."
- **Verdict:** **NARROWED**
- **Flaw & Evidence:**
  - Cached outputs `/tmp/chm13_chr22_lifted_b10.bed` (`514` rows: `469 EXACT`, `1 PARTIAL`, `44 SPLIT` across `19` unique `win_*` IDs) and `/tmp/chm13_chr22_unmapped_b10.bed` (`11` rows) match the reported counts (`469 + 1 + 19 + 11 = 500` matching `/tmp/grch38_chr22_windows.bed`).
  - However, neither `hprc-v1.0-mc-grch38.gfa.gz` nor a script generating `/tmp/grch38_chr22_windows.bed` is in the repository. Moreover, the `/tmp/` files are timestamped `Nov 19 11:01`, whereas `src/bin/horizon_lift.rs` was modified in commit `0e70693` (`11:20:04`) to deduplicate GFA parsing via `stream_gfa`, so the `256.60s` / `619.2 MB` benchmark has not been reproduced on `HEAD`.
- **Correction:** Re-run `horizon_lift --bridge-bp 10` on `HEAD` and commit the window-generation command/script to `scripts/benchmarks/`.

---

## F4 -- C5: `chr22` pericentromere `horizon_wga` and `horizon_clean` counts reproduce on `/tmp/` files, but input generation is missing from repo

- **Claim:** C5 (`docs/BENCHMARK_RESULTS.md:23,59,63`) -- "Real HPRC `chr22` Pericentromere 5 Mbp (`GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M`): **23.45s**, **3,084.5 MB** peak RSS, **5,894,191** aligned pairs, `246` fwd / `184` rev LCBs, `17,789` Horizon Bubbles (`horizon_clean`: `6,030` links kept, **0.01s**, **7.8 MB** RSS)."
- **Verdict:** **NARROWED**
- **Flaw & Evidence:**
  - Re-running `./target/release/horizon_wga /tmp/grch38_chr22_peri_5m.fa /tmp/chm13_chr22_peri_5m.fa --resolve-dups` on `HEAD` reproduces `aligned_pairs=5894191`, `LF:i:246`, `LR:i:184`, `HB:i:17789` (`elapsed_s=49.0900` under shared load vs `23.45s` reported).
  - Re-running `./target/release/horizon_clean --gfa /tmp/chr22_peri_5m.gfa` reproduces `total_links=6030 kept=6030 severed=0 bubbles=0 split_walks=0 elapsed=0.05s peak_rss_mb=8.2` (though `/tmp/chr22_peri_5m.gfa` contains `*` for all `S` sequences, so `horizon_clean` cannot score sequence flanks and trivially keeps all links).
  - Neither the 5 Mbp pericentromere FASTAs nor `/tmp/chr22_peri_5m.gfa` have a generator script in the repository; on a clean checkout without `/tmp/`, C5 is not reproducible.
- **Correction:** Add a reproducible slice-extraction script to `scripts/benchmarks/` and note that `/tmp/chr22_peri_5m.gfa` has `*`-sequence segments.

---

## F5 -- C6: `horizon_resolve --collapse-dups` reproduces exact counts on `/tmp/bench_resolve_10k.gfa`, but the synthetic GFA has no generator script

- **Claim:** C6 (`docs/BENCHMARK_RESULTS.md:64`) -- "| GFA E-Graph Resolver | `horizon_resolve --collapse-dups` | 200 segs, 20 walks x 225 steps (`4,500` steps) | `280` windows, `840` dups collapsed, `308,760.9` bits saved | **24.0 MB** | 4.48s |"
- **Verdict:** **NARROWED**
- **Flaw & Evidence:**
  - Re-running `./target/release/horizon_resolve --collapse-dups /tmp/bench_resolve_10k.gfa` on `HEAD` reproduces `segments=200 anchors=180 paths=0 walks=20 windows=280 simplified=280 dups_detected=840 dups_collapsed=840 invs_cancelled=0 mdl_bits_saved=308760.9 elapsed_s=9.399 peak_rss_mb=23.3`.
  - However, `/tmp/bench_resolve_10k.gfa` exists only in `/tmp/` with no generator script or test fixture in the repository, so C6 cannot be reproduced from a clean checkout.
- **Correction:** Commit the generator script or fixture for `/tmp/bench_resolve_10k.gfa` into the repository.

---

## F6 -- C8: `ContigIndex` structure verified in code, but cancelled `chr6` OOM run has no saved log or reproducible harness

- **Claim:** C8 (`docs/BENCHMARK_RESULTS.md:26`) -- "(After replacing `HashSet<(usize, usize)>` with sorted/deduplicated `Vec<(usize, usize)>`, full 170.8 Mbp `GRCh38.chr6` vs 110.7 Mbp `CHM13.chr6` still reached **26.5 GB** RSS and `47:33` CPU before cancellation on the 30 GiB host due to `ContigIndex`'s per-position `HashMap<u64, Vec<usize>>` seed maps, ~0.29 GB/Mbp.)"
- **Verdict:** **NARROWED**
- **Flaw & Evidence:**
  - `ContigIndex` (`src/seeding/mod.rs:329-333`) stores `seeds_fwd: HashMap<u64, Vec<usize>>` and `seeds_rc: HashMap<u64, Vec<usize>>`, `src/mdl/mod.rs:456` uses `Vec<(usize, usize)>`, and `/tmp/hprc_chr6.gfa` / `/tmp/hprc_chr6.maf` are `0`-byte files from `11:05`.
  - However, no log file, script, or RSS measurement trace for the cancelled `26.5 GB` / `47:33` CPU run is stored in the repository, and the `~0.29 GB/Mbp` scaling rate has no documented denominator.
- **Correction:** Record the exact command, RSS measurement method, and Mbp denominator (e.g., `~91 Mbp` vs `281.5 Mbp` combined) if retaining this note.

---

## C1 and C7 (Survived Reproduction Checks)

- **C1 (`SURVIVED`):** `stream_gfa` (`src/gfa/mod.rs:98-128`) supports `.gfa`, `.gfa.gz` via `MultiGzDecoder`, and `"-"` for `stdin`, and both unit tests (`streams_all_record_types_and_handles_empty_link_orients` at `src/gfa/mod.rs:219-276` and `streams_gzipped_gfa_file` at `src/gfa/mod.rs:279-318`) pass in `cargo test`.
- **C7 (`SURVIVED`):** `align_seq_to_dag_seeded` (`src/mdl/dag.rs:460-540`) uses `MASK_W10` + `MASK_RY20` and forward/backward reachability closure (for `nodes.len() > 8` with non-empty `seed_hits`), its unit test `seeded_dag_aligns_ignoring_disconnected_decoy_component` (`src/mdl/dag.rs:621`) passes in `cargo test`, and `./target/release/horizon_wga --graph /tmp/test_dag.gfa /tmp/test_dag_q.fa` reproduces `>n0>n2>n3`, `kept_bp=96`, `score_bits=124.4`.

---

## Verdict Summary

| Claim | Verdict |
| :--- | :--- |
| C1 | SURVIVED |
| C2 | REFUTED |
| C3 | REFUTED |
| C4 | NARROWED |
| C5 | NARROWED |
| C6 | NARROWED |
| C7 | SURVIVED |
| C8 | NARROWED |
