# Round 1 Attack — L4 (Provenance)

**Target:** `docs/STRATEGIC_PLAN.md` (`reports/rounds/round-1/ATTACKED-STRATEGIC_PLAN.md`) and `docs/PANGENOME_SUITE.md` (`reports/rounds/round-1/ATTACKED-PANGENOME_SUITE.md`), verified identical via `diff -u`.  
**Baseline checks (Rules 2–4):**
- **Prior rounds:** `./scripts/round history "STRATEGIC_PLAN"` returned 0 closed rounds and 0 prior mentions (Round 1 is the first open round).
- **Held literature:** `manuscript/references.bib` (lines 1–106) holds 11 BibTeX entries (`adams2007bayesian`, `allison1992finite`, `angelopoulos2024conformal`, `armstrong2020progressive`, `blanchette2004aligning`, `chiaromonte2002scoring`, `earl2014alignathon`, `foster2007combinators`, `harris2007improved`, `hickey2013hal`, `willsey2021egg`). None are `vg`, `odgi`, `liftoff`, or HPRC graph-crash papers.
- **Method cards:** No method-cards directory or card files exist in the repository (`round.toml`, `docs/`, `reports/`).

---

## F1 -- C1: Three of four named CLI tools do not exist and existing binaries do not compose via Unix pipes

- **Claim:** `C1` (`docs/STRATEGIC_PLAN.md:6`) — *"We are not building a monolithic 'run-all' executable. We are building a suite of specialized tools (`horizon_clean`, `horizon_wga`, `horizon_resolve`, `horizon_lift`) that compose perfectly via Unix pipes."*
- **Verdict:** **REFUTED**
- **Evidence:**
  - `src/bin/` contains five binaries: `bench_chr1.rs`, `horizon_eval.rs`, `horizon_guide.rs`, `horizon_merge.rs`, `horizon_wga.rs`. Neither `horizon_clean`, `horizon_resolve`, nor `horizon_lift` exists anywhere in `src/`, `Cargo.toml`, or git history (`rg -n "horizon_clean|horizon_resolve|horizon_lift" src/` returns 0 matches).
  - Existing binaries do not compose via Unix pipes (`stdin`/`stdout`): `rg -i "stdin" src/` returns 0 matches. `horizon_wga` (`src/bin/horizon_wga.rs:120-124`) requires positional file paths `<target.fa> <query.fa>` and writes outputs to `--maf <path>` / `--gfa <path>` files; `horizon_merge` (`src/bin/horizon_merge.rs:11-35`) reads `<sp=path.maf>` via `fs::read_to_string` and writes `<maf_out>` and `<gfa_out>` via `fs::write`.
- **Correction:** Replace `docs/STRATEGIC_PLAN.md:6` with an accurate inventory: *"Currently implemented binaries are `horizon_wga`, `horizon_merge`, `horizon_guide`, and `horizon_eval` (all file-path driven); `horizon_clean`, `horizon_resolve`, `horizon_lift`, and Unix-pipe (`stdin`/`stdout`) streaming are unimplemented proposals."*

---

## F2 -- C2: MDL replaces DP scoring matrices, not all heuristic thresholds, and has zero pangenome hairball provenance

- **Claim:** `C2` (`docs/STRATEGIC_PLAN.md:9`) — *"We are replacing arbitrary affine gap penalties and scoring matrices with Information-Theoretic bit-savings (MDL). This inherently prevents the 'hairball' over-alignment problem."*
- **Verdict:** **NARROWED**
- **Evidence:**
  - Pairwise DP (`src/mdl/dp.rs:81`, `banded_mdl_gap`) uses MDL bit-savings, but `src/mdl/mod.rs` still relies on heuristic/configurable thresholds: `let max_gap_bp = 2400;` (`src/mdl/mod.rs:123`), `hsp_threshold: Option<i32>` (`src/mdl/mod.rs:47`), `min_fallback_cov` (`src/mdl/mod.rs:184-191`), and `min_block_bits: f64` (`src/mdl/dp.rs:81`).
  - There is zero provenance in `src/` or `docs/` that MDL prevents pangenome "hairball" over-alignment: `src/mdl/` only aligns linear pairwise contigs, and `docs/BENCHMARK_RESULTS.md:39` leaves the sole pangenome entry unrun (`| HPRC (v1.1) | GRCh38 | CHM13 | Human Pangenome | *(Graph Eval)* | - | - |`).
- **Correction:** Narrow `docs/STRATEGIC_PLAN.md:9` to: *"We replace affine gap penalties and substitution matrices in pairwise DP (`src/mdl/dp.rs`) with MDL bit-savings (alongside heuristic `max_gap_bp` and HSP thresholds in `src/mdl/mod.rs`); whether this prevents pangenome 'hairball' over-alignment is untested (`docs/BENCHMARK_RESULTS.md:39`)."*

---

## F3 -- C3: E-graph resolution has no HLA or pangenome graph provenance

- **Claim:** `C3` (`docs/STRATEGIC_PLAN.md:10`) — *"Using congruence-closure to resolve complex, nested structural paralogs (like the HLA region). This natively models evolutionary trajectories without producing unnavigable cyclical graphs."*
- **Verdict:** **NARROWED**
- **Evidence:**
  - `src/egraph.rs:1-385` implements a congruence-closed `EGraph` (citing Willsey et al. 2021, `manuscript/references.bib:97-106`) used on pairwise/merged LCBs (`resolve_lcb_duplications` and `resolve_alignment_duplications`, called in `src/bin/horizon_wga.rs:254` and `src/merge.rs`).
  - No test, benchmark, or dataset in `src/`, `docs/BENCHMARK_DATASETS.md`, or `docs/BENCHMARK_RESULTS.md` evaluates the HLA region or demonstrates avoidance of cyclical pangenome graphs (`HPRC (v1.1)` at `docs/BENCHMARK_RESULTS.md:39` is `*(Graph Eval)*`).
- **Correction:** Narrow `docs/STRATEGIC_PLAN.md:10` to state that `src/egraph.rs` collapses tandem duplications and inversions in pairwise/merged LCBs, and that HLA and pangenome-graph evaluation remain unverified.

---

## F4 -- C4: Uncited `vg`/`odgi` crash claim and false "zero-copy" claim for unused `src/gfa/mod.rs`

- **Claim:** `C4` (`docs/STRATEGIC_PLAN.md:11`) — *"`vg` and `odgi` crash on 12GB pangenomes due to RAM exhaustion. Our custom zero-copy GFA streaming parser (`src/gfa/mod.rs`) solves the I/O bottleneck entirely."*
- **Verdict:** **REFUTED**
- **Evidence:**
  1. **Uncited external-tool claim:** Neither `docs/STRATEGIC_PLAN.md`, `docs/PANGENOME_SUITE.md`, `docs/BENCHMARK_RESULTS.md`, `docs/BENCHMARK_DATASETS.md`, nor `manuscript/references.bib` cites any benchmark, paper, or error log showing `vg` and `odgi` crashing on 12GB pangenomes.
  2. **Not zero-copy:** `src/gfa/mod.rs:75-120` (`stream_gfa`) reads lines into a `String` via `reader.read_line(&mut line)` (line 76) and heap-allocates owned `String` and `Vec<u8>` copies for every field of `GfaSegment` (`.to_string()`, `.to_vec()`, lines 87-89), `GfaLink` (lines 93-97), `GfaPath` (lines 102-104), and `GfaWalk` (lines 110-115).
  3. **Dead code & unbenchmarked:** `rg -n "stream_gfa" src/` shows `stream_gfa` has 0 callers anywhere in `src/` or `src/bin/`, and `src/gfa/mod.rs` contains 0 unit tests. No 12GB GFA I/O benchmark exists (`docs/BENCHMARK_RESULTS.md:39`).
- **Correction:** Remove the uncited `vg`/`odgi` 12GB crash assertion and revise `docs/STRATEGIC_PLAN.md:11` to describe `src/gfa/mod.rs` accurately as an unintegrated, allocating line-by-line GFA callback reader with no 12GB benchmark yet.

---

## F5 -- C5: `stdin`/`stdout` GFA/MAF/VCF communication does not exist, and `src/hal/` contradicts the anti-HAL framing

- **Claim:** `C5` (`docs/STRATEGIC_PLAN.md:16-17`) — *"We will **never invent a proprietary binary format**. Our tools will communicate purely by streaming standard GFA (v1.0/v1.1) and standard MAF/VCF over `stdin`/`stdout`."*
- **Verdict:** **REFUTED**
- **Evidence:**
  - **No VCF support:** `rg -i "vcf" src/` returns 0 matches.
  - **No `stdin` streaming:** `rg -i "stdin" src/` returns 0 matches; `horizon_wga` and `horizon_merge` read and write files via paths (`src/bin/horizon_wga.rs:120-124`, `src/bin/horizon_merge.rs:11-35`).
  - **Undisclosed `HAL` implementation in repo (claim-list omission):** `docs/STRATEGIC_PLAN.md:15-16` criticizes Cactus's `.hal` format, yet the repository itself implements a pure-Rust HDF5 `.hal` reader/writer in `src/hal/mod.rs`, `src/hal/hdf5_read.rs`, and `src/hal/hdf5_write.rs`, exposed via `horizon_wga --hal-info` (`src/bin/horizon_wga.rs:109-112, 150-169`).
- **Correction:** Remove the claim that tools communicate over `stdin`/`stdout` via GFA/MAF/VCF, and disclose that `horizon-aligner` currently uses file-based FASTA/MAF/GFA I/O plus an internal pure-Rust `.hal` (HDF5) reader/writer (`src/hal/`).

---

## F6 -- C6: Zero binaries read GFA or MAF from `stdin`

- **Claim:** `C6` (`docs/STRATEGIC_PLAN.md:23`) — *"Every tool must read GFA/MAF from `stdin` and write to `stdout`."*
- **Verdict:** **REFUTED**
- **Evidence:**
  - `rg -i "stdin" src/` returns 0 matches across all five binaries in `src/bin/` (`bench_chr1.rs`, `horizon_eval.rs`, `horizon_guide.rs`, `horizon_merge.rs`, `horizon_wga.rs`).
  - `horizon_wga` (`src/bin/horizon_wga.rs:120-124`) requires `<target.fa> <query.fa>` and `--maf`/`--gfa` file flags; `horizon_merge` (`src/bin/horizon_merge.rs:11-35`) requires `<maf_out> <gfa_out> <sp=path.maf>...`.
- **Correction:** Mark `docs/STRATEGIC_PLAN.md:23` explicitly as an unimplemented future requirement, noting that no current binary supports `stdin` GFA/MAF input.

---

## F7 -- C7: Linear GFA-stream-to-E-graph bubble stitching is unimplemented

- **Claim:** `C7` (`docs/STRATEGIC_PLAN.md:27`) — *"We will stream the pangenome linearly and only invoke the E-graph engine when we hit a complex, tangled sub-graph ('hairball'). The E-graph resolves the optimal path for that bubble, and we stitch the resolved path back into the linear stream."*
- **Verdict:** **REFUTED**
- **Evidence:**
  - `src/gfa/mod.rs` (`stream_gfa`) has 0 callers in `src/`.
  - `src/egraph.rs` (`resolve_lcb_duplications`, `resolve_alignment_duplications`) operates only on pairwise `ContigAlignment` or merged `MafBlock` structs (`src/formats.rs`, `src/merge.rs`), never on GFA sub-graphs or streamed bubbles. No GFA bubble-extraction or path-stitching code exists in `src/`.
- **Correction:** Label `docs/STRATEGIC_PLAN.md:27` as a future design proposal with no existing implementation in `src/`.

---

## F8 -- C8: Release strategy targets three non-existent tools and cites unrun HPRC benchmarks

- **Claim:** `C8` (`docs/STRATEGIC_PLAN.md:31`) — *"The HPRC produces massive graphs that researchers cannot use because standard tooling crashes or produces tangled paths. We will focus our benchmarking and release strategy entirely on `horizon_clean`, `horizon_resolve`, and `horizon_lift` first to capture the pangenome community immediately."*
- **Verdict:** **REFUTED**
- **Evidence:**
  - **Uncited premise:** No citation or benchmark in `docs/` or `manuscript/references.bib` supports the claim that HPRC graphs cannot be used because standard tooling crashes.
  - **Non-existent targets:** `horizon_clean`, `horizon_resolve`, and `horizon_lift` do not exist in `src/bin/` or `Cargo.toml`.
  - **Benchmark reality:** `docs/BENCHMARK_RESULTS.md:23-39` contains 7 populated pairwise/MSA benchmarks (`multiz100way`, `cactus447way`, `multiz470way`, `cactus241way`) and **zero** HPRC graph results (`*(Graph Eval)*` at line 39).
- **Correction:** Remove the uncited HPRC tooling-crash assertion and state that `horizon_clean`, `horizon_resolve`, and `horizon_lift` are unimplemented and have zero HPRC benchmark data.

---

## F9 -- C9: `horizon_clean` does not exist and `GfaLink` lacks sequences needed for BOCPD/MDL scoring

- **Claim:** `C9` (`docs/PANGENOME_SUITE.md:8`, under `horizon_clean`) — *"Stream an external `.gfa` file and evaluate every Link (`L` line) using the Horizon BOCPD/MDL scoring engine."*
- **Verdict:** **REFUTED**
- **Evidence:**
  - `horizon_clean` does not exist in `src/bin/` or `Cargo.toml`.
  - `src/gfa/mod.rs:50-124` (`stream_gfa`) has 0 callers in `src/` and never invokes `src/bocpd.rs` or `src/mdl/`.
  - `GfaLink` (`src/gfa/mod.rs:21-27`) stores only segment IDs, orientations, and `overlap` — not segment sequences — so evaluating an `L` line with `bocpd_gc_profile` / `banded_mdl_gap` is impossible from `L` lines alone without retaining `S` sequences in memory.
- **Correction:** Mark `docs/PANGENOME_SUITE.md:5-10` as a unimplemented proposal and note that sequence-level MDL/BOCPD scoring of `L` lines requires segment-sequence lookup, not pure stateless `L`-line streaming.

---

## F10 -- C10: `horizon_wga --graph` and graph-aware banded MDL DP do not exist

- **Claim:** `C10` (`docs/PANGENOME_SUITE.md:16`, under `horizon_wga --graph`) — *"Thread the query sequence through the graph, performing banded MDL dynamic programming across graph bifurcations."*
- **Verdict:** **REFUTED**
- **Evidence:**
  - `src/bin/horizon_wga.rs:32-145` (`parse_args`) has no `--graph` flag (supported flags are `--maf`, `--gfa`, `--min-shared-seeds`, `--min-top-lcb`, `--top-k`, `--resolve-dups`, `--rbest`, `--respect-softmask`, `--ry-seeds`, `--no-ry-seeds`, `--min-fallback-cov`, `--hsp-threshold`, `--hal-info`).
  - `src/mdl/dp.rs:81` (`banded_mdl_gap(seq_x: &[u8], seq_y: &[u8], ...)`) only accepts two linear byte slices; no graph-bifurcation DP exists anywhere in `src/mdl/`.
- **Correction:** Mark `docs/PANGENOME_SUITE.md:12-17` as unimplemented: `horizon_wga` currently performs pairwise FASTA-to-FASTA alignment only.

---

## F11 -- C11: `horizon_lift` does not exist and `coord::Interval` is not even compiled into `horizon_aligner`

- **Claim:** `C11` (`docs/PANGENOME_SUITE.md:29`, under `horizon_lift`) — *"Traverse the GFA Walk (`W`) lines and Path (`P`) lines using Horizon's strict `coord::Interval` zero-based half-open arithmetic."*
- **Verdict:** **REFUTED**
- **Evidence:**
  - `horizon_lift` does not exist in `src/bin/` or `Cargo.toml`.
  - `src/coord.rs:5-38` defines `Interval`, but `src/lib.rs:3-10` does **not** declare `pub mod coord;` (`rg -n "pub mod coord|coord::Interval" src/` returns 0 matches). Thus `coord::Interval` is dead, uncompiled code.
  - `src/gfa/mod.rs:100-117` parses `P` and `W` lines into `GfaPath` and `GfaWalk`, but `stream_gfa` has 0 callers in `src/` and never uses `coord::Interval`.
- **Correction:** Mark `docs/PANGENOME_SUITE.md:26-31` as unimplemented, and note that `src/coord.rs` is not currently exported in `src/lib.rs`.

---

## Verdict Summary

| C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 |
|---|---|---|---|---|---|---|---|---|---|---|
| REFUTED | NARROWED | NARROWED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED |
