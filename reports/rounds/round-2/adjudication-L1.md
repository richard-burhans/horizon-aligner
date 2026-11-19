# Round 2 Adjudication -- Correctness Lens (L1)

**Attack digest verified:** `9a7c332f2fe0` (`9a7c332f2fe00aec1d53679372a8e124c2f03d75ba3359e6fbc79de6e8cfff56` on `reports/rounds/round-2/attack-L1-correctness.md`).

**Method cards & prior-round check:** No method cards exist in `refs/methods/` (checked; none present in `horizon-aligner`). Prior adjudications (`scripts/round history "gfa" --kind adjudication`) settled `src/gfa/mod.rs` non-zero-copy allocation and `L`-line empty-orientation panic in Round 1 (`R1` `F4`–`F11`) and left open `R1-IMPL`, `R1-GFA`, and `R1-HMM`. None of Round 2's eight findings (`F1`–`F8`) merely re-file a settled Round 1 defect; they test the newly implemented binaries (`horizon_clean`, `horizon_resolve`, `horizon_lift`, `horizon_wga --graph`) and the newly reported HPRC benchmarks.

---

## F1 -- C3 (`docs/STRATEGIC_PLAN.md:3`): **UPHELD**

* **What survives:** `cargo test --no-fail-fast` reports **40 passing tests and 1 ignored** (`35 passed; 0 failed; 1 ignored` in `src/lib.rs`, plus `1 passed` each in `horizon_clean`, `horizon_gfa`, `horizon_lift`, `horizon_resolve`, and `horizon_wga`), while **4 of 9** binaries (`bench_chr1`, `horizon_eval`, `horizon_guide`, `horizon_merge`) have **0** binary unit tests (though `horizon_merge` logic is covered by library tests in `src/lib.rs`).
* **Repair:** Endorsed on `docs/STRATEGIC_PLAN.md:3` (report `40 passing unit tests (1 ignored)` and note which 4 binaries have `0` binary unit tests, optionally noting library-level coverage for `horizon_merge`). **The fix is free** (text-only edit, no re-measurement).

---

## F2 -- C2 (`docs/BENCHMARK_RESULTS.md:5`): **UPHELD**

* **What survives:** All four defects are verified in source:
  1. `mdl_xdrop_extend` (`src/mdl/dp.rs:369-532`) takes no `HmmParams`, applies no `mm_bits`/`im_bits` transition normalization, and hardcodes `IDENTITY_PRIOR` (`line 445`), `GAP_EXTEND_BITS` (`458, 473`), and `GAP_OPEN_BITS` (`460, 475`).
  2. `banded_mdl_gap_with_params` (`src/mdl/dp.rs:270-293`) has no backpointer matrix and selects predecessors at `(i, j)` by comparing `m_dp`, `ix_dp`, `iy_dp` without transition weights (`mm_bits`, `im_bits`, `go_bits`, `ge_bits`).
  3. `align_seq_to_dag_mdl` (`src/mdl/dag.rs:165-172, 258-315`) permits `Iy` transitions at both `i_u = nu` (end of predecessor `u`) and `i_s = 0` (start of successor `s`), double-counting boundary `(None, Some(j))` columns in Forward/Backward `Z`; additionally, `horizon_wga --graph` (`src/bin/horizon_wga.rs:195`) passes `&HmmParams::default()` rather than `estimate_hmm_params`.
  4. `estimate_hmm_params` (`src/mdl/mod.rs:828-836`) estimates `gap_extend_bits` via mean-length pseudo-count shrinkage `(total_indel_len + 2.5) / (indel_events + 1.0)` rather than Beta prior shrinkage.
* **Repair:** Endorsed on the documentation narrowing in `docs/BENCHMARK_RESULTS.md:5` (**the doc fix is free**, text-only edit). The code repairs (adding a Viterbi backpointer matrix in `src/mdl/dp.rs` and removing duplicate `Iy` boundary transitions in `src/mdl/dag.rs`) are **not free**—they require code changes and re-running the benchmark suite.

---

## F3 -- C4 (`docs/BENCHMARK_RESULTS.md:21`): **UPHELD**

* **What survives:** In `src/bin/horizon_lift.rs:256-270`, coalescing under `--bridge-bp 10` sets `prev.src_end = piece.src_end`, so `mapped_bp` at `line 490` includes `<= 10 bp` unmapped source gaps and `line 508` classifies bridged windows as `EXACT`; furthermore, `run_lift` (`lines 316, 419`) calls `stream_gfa` twice, so `"-"` `stdin` fails on Pass 2.
* **Repair:** Endorsed on the documentation clarification in `docs/BENCHMARK_RESULTS.md:21` noting that `469 exact` includes `<= 10 bp` bridged windows (**the doc fix is free**). Separating `mapped_bp` from bridged span in `LiftedPiece` and re-running the `256.60s` HPRC liftover benchmark is **not free**.

---

## F4 -- C5 (`docs/BENCHMARK_RESULTS.md:23`): **UPHELD**

* **What survives:** `5,894,191` aligned pairs on a `5,000,000 bp` locus (`117.9%`) is a many-to-many `(x, y)` pair count deduplicated only by exact `(x, y)` across `246 fwd / 184 rev` LCBs (`src/mdl/mod.rs:745-746`, with `--resolve-dups` only appending `# EG:Z:` comments at `src/bin/horizon_wga.rs:308-317`), and `horizon_clean` (`src/bin/horizon_clean.rs:152, 228-237, 246`) never severs walk-supported (`sup > 0`) links and fails on `"-"` `stdin` in Pass 2.
* **Repair:** Endorsed on `docs/BENCHMARK_RESULTS.md:23` to label `5,894,191` as many-to-many `(x, y)` aligned pairs rather than 1-to-1 aligned bases. **The fix is free** (text-only edit).

---

## F5 -- C6 (`docs/BENCHMARK_RESULTS.md:64`): **UPHELD**

* **What survives:** In `src/bin/horizon_resolve.rs:424-441`, `--collapse-dups` emits original `W`-line `start`/`end` coordinates (`cols[4]`, `cols[5]`) verbatim even when walk length shrinks (e.g., `0..80` on a collapsed `60 bp` walk at `line 507`) and copies `L` lines verbatim (`line 446`); without `--out`, `main` (`line 465`) writes to `io::sink()` instead of `stdout`, and `parse_args` (`line 103`) rejects `"-"` `stdin`.
* **Repair:** Endorsed on the documentation clarification in `docs/BENCHMARK_RESULTS.md:64` noting that `horizon_resolve` rewrites `P`/`W` walks without recomputing `W` coordinates or `L` edges (**the doc fix is free**). Updating `W`-line `end` in `src/bin/horizon_resolve.rs:424-441` is a code change (**not free**).

---

## F6 -- C7 (`docs/PANGENOME_SUITE.md:18`): **UPHELD**

* **What survives:** `align_seq_to_dag_seeded` (`src/mdl/dag.rs:466, 485, 503-518`) performs backward ancestor expansion followed by forward descendant expansion (activating every downstream branch from any shared ancestor, rather than an s-t / forward-backward intersection), and falls back to full `align_seq_to_dag_mdl` when `nodes.len() <= 8` or `seed_hits.is_empty()`.
* **Repair:** Endorsed on `docs/PANGENOME_SUITE.md:18` to describe the closure accurately as ancestor-then-descendant expansion with fallback to `align_seq_to_dag_mdl` for `nodes.len() <= 8` or empty `seed_hits`. **The fix is free** (text-only edit).

---

## F7 -- C8 (`docs/BENCHMARK_RESULTS.md:26`): **UPHELD**

* **What survives:** `26.5 GB` RSS on `170.8 + 110.7 Mbp` (`281.5 Mbp` combined) is `0.094 GB/Mbp` (`0.155 GB/Mbp` against `170.8 Mbp` alone), whereas `~0.29 GB/Mbp` matches the completed `20 Mbp` `chr6` q-arm run (`5,863.3 MB / 20 Mbp = 0.293 GB/Mbp`). Note that the same conflation also appears in `docs/STRATEGIC_PLAN.md:11`.
* **Repair:** Endorsed on `docs/BENCHMARK_RESULTS.md:26` (and `docs/STRATEGIC_PLAN.md:11`) to attribute `~0.29 GB/Mbp` to the `20 Mbp` `chr6` q-arm run rather than the cancelled full-`chr6` run. **The fix is free** (text-only edit).

---

## F8 -- C1 (`docs/STRATEGIC_PLAN.md:13`): **UPHELD**

* **What survives:** `stream_gfa_reader` (`src/gfa/mod.rs:198-203`) silently drops GFA v1.0 `C` and GFA v1.1 `J` (plus `U`/`O`) records and coerces unparseable `W` `start`/`end` fields (including `*`) to `0`, while `parse_path_steps` (`line 58`) can panic on multi-byte UTF-8 segment names (`split_at(tok.len() - 1)`), even though `S`/`L`/`P`/`W` streaming, `.gfa.gz` `MultiGzDecoder`, `"-"` `stdin`, and unit tests (`lines 98-128, 219-318`) are present.
* **Repair:** Endorsed on `docs/STRATEGIC_PLAN.md:13` to state that `stream_gfa` parses GFA v1.0/v1.1 `S`/`L`/`P`/`W` records (ignoring `C`/`J`/`U`/`O`). **The fix is free** (text-only edit).

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `cargo test` runs 40 passing + 1 ignored tests (not 39), and 4 of 9 binaries have 0 binary unit tests; doc fix is free |
| F2 | UPHELD | `mdl_xdrop_extend` uses hardcoded constants, `banded_mdl_gap_with_params` lacks backpointers, `align_seq_to_dag_mdl` double-counts boundary `Iy`, and `gap_extend_bits` uses mean-length shrinkage; doc narrowing is free, code fixes require re-benchmarking |
| F3 | UPHELD | `horizon_lift --bridge-bp 10` counts `<= 10 bp` bridged gaps as `EXACT` and `stdin` (`"-"`) fails on Pass 2; doc clarification is free |
| F4 | UPHELD | `5,894,191` on a 5 Mbp locus is a many-to-many `(x, y)` pair count, and `horizon_clean` never severs `sup > 0` links or supports `stdin` Pass 2; doc fix is free |
| F5 | UPHELD | `horizon_resolve --collapse-dups` leaves stale `W`-line `start`/`end` coordinates and `L` lines, drops output without `--out` to `io::sink()`, and rejects `"-"`; doc clarification is free |
| F6 | UPHELD | `align_seq_to_dag_seeded` uses ancestor-then-descendant expansion (not s-t intersection) and falls back to `align_seq_to_dag_mdl` for `<= 8` nodes or no hits; doc fix is free |
| F7 | UPHELD | `26.5 GB` over `170.8 + 110.7 Mbp` is `0.094 GB/Mbp`, whereas `~0.29 GB/Mbp` is from the `20 Mbp` `chr6` q-arm run (`5,863.3 MB / 20 Mbp`); doc fix in `BENCHMARK_RESULTS.md` and `STRATEGIC_PLAN.md` is free |
| F8 | UPHELD | `stream_gfa` parses `S`/`L`/`P`/`W` while dropping `C`/`J`/`U`/`O`, coerces invalid `W` coords to `0`, and `parse_path_steps` can panic on multi-byte UTF-8; doc fix is free |
