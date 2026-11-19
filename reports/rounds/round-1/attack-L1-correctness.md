# Round 1 Adversarial Attack -- Correctness (L1)

All claims in `ATTACKED-STRATEGIC_PLAN.md` and `ATTACKED-PANGENOME_SUITE.md` were checked against the actual Rust codebase in `src/` and `Cargo.toml`.

Note on the claim list: both attacked documents present unimplemented tools (`horizon_clean`, `horizon_resolve`, `horizon_lift`), unimplemented CLI modes (`--graph`, `stdin`/`stdout` streaming), and an uncompiled module (`src/coord.rs`) as concrete architectural mechanisms or existing components ("Our custom zero-copy GFA streaming parser (`src/gfa/mod.rs`)", "Horizon's strict `coord::Interval`"). Furthermore, the codebase's core alignment engine contradicts the stated design principles (hardcoded affine gap penalties in `src/bocpd.rs`, HAL/HDF5 support in `src/hal/`).

## F1 -- C1: `horizon_clean`, `horizon_resolve`, `horizon_lift`, and Unix-pipe composition do not exist

- **Claim:** C1 (`docs/STRATEGIC_PLAN.md:6`) -- "We are not building a monolithic 'run-all' executable. We are building a suite of specialized tools (`horizon_clean`, `horizon_wga`, `horizon_resolve`, `horizon_lift`) that compose perfectly via Unix pipes."
- **Verdict:** REFUTED
- **Evidence:**
  - `src/bin/` contains only `bench_chr1.rs`, `horizon_eval.rs`, `horizon_guide.rs`, `horizon_merge.rs`, and `horizon_wga.rs`. None of `horizon_clean`, `horizon_resolve`, or `horizon_lift` exist anywhere in `src/` or `Cargo.toml`.
  - `grep -rn "stdin" src/` and `grep -rn "stdout" src/` both return zero matches: no binary in `src/bin/` reads from `stdin` or writes alignment streams to `stdout` for Unix-pipe composition.
- **Correction:** State that only `horizon_wga` and `horizon_merge` (plus evaluation/guide binaries `horizon_eval`, `horizon_guide`, `bench_chr1`) exist, and none of them compose via `stdin`/`stdout` Unix pipes.

## F2 -- C2: `src/bocpd.rs` and `src/mdl/dp.rs` use hardcoded arbitrary affine gap penalties

- **Claim:** C2 (`docs/STRATEGIC_PLAN.md:9`) -- "We are replacing arbitrary affine gap penalties and scoring matrices with Information-Theoretic bit-savings (MDL). This inherently prevents the 'hairball' over-alignment problem."
- **Verdict:** REFUTED
- **Evidence:**
  - `src/bocpd.rs:4-7` defines hardcoded affine gap constants and priors:
    `pub(crate) const GAP_OPEN_BITS: f64 = -4.5;`
    `pub(crate) const GAP_EXTEND_BITS: f64 = -0.8;`
    `pub(crate) const IDENTITY_PRIOR: f64 = 0.65;`
    `pub(crate) const MIN_POST_PROB: f64 = 0.15;`
  - `src/mdl/dp.rs:105-107` converts these fixed constants directly into DP transition penalties (`let go = GAP_OPEN_BITS * ln2; let ge = GAP_EXTEND_BITS * ln2;`). They are not derived from data or universal MDL code lengths.
- **Correction:** Replace "replacing arbitrary affine gap penalties" with "combining BOCPD local GC composition priors (`joint_pair_bits`) with fixed affine gap-open (`-4.5` bits) and gap-extend (`-0.8` bits) constants."

## F3 -- C3: `src/egraph.rs` cannot resolve nested structural paralogs or model evolutionary trajectories

- **Claim:** C3 (`docs/STRATEGIC_PLAN.md:10`) -- "Using congruence-closure to resolve complex, nested structural paralogs (like the HLA region). This natively models evolutionary trajectories without producing unnavigable cyclical graphs."
- **Verdict:** REFUTED
- **Evidence:**
  - `src/egraph.rs:165-273` (`saturate_synteny`) implements only four rewrite rules on linear 1D `seq`/`inv`/`dup` expressions: `inv(inv(x)) -> x` (lines 189-190), `inv(seq(a,b)) -> seq(inv(b),inv(a))` (lines 191-206), `seq(a,a) -> dup(a)` (lines 211-217), and `seq` associativity (lines 223-253).
  - It cannot collapse >2-copy tandem arrays (there is no `seq(a, dup(a)) -> dup(a)` rule) or non-tandem/interleaved paralogs (e.g., `A B A`), and it has no phylogenetic or evolutionary-trajectory model.
  - Its only caller, `resolve_alignment_duplications` (`src/formats.rs:73-103`), buckets pairwise alignment blocks into 1000-bp windows (`win in bins.chunks(12)`), collects `dup` labels, and emits `# EG:Z:...` comment lines; it never resolves or outputs a de-cycled graph.
- **Correction:** Narrow the claim to: `src/egraph.rs` provides a bounded 4-rule e-graph that detects adjacent tandem duplications and inverted tandem duplications in 12-block windows of 1000-bp binned pairwise alignments and emits `# EG:Z:` comment annotations.

## F4 -- C4: `src/gfa/mod.rs` is not zero-copy, is unused, and panics on empty `L`-line orientation fields

- **Claim:** C4 (`docs/STRATEGIC_PLAN.md:11`) -- "`vg` and `odgi` crash on 12GB pangenomes due to RAM exhaustion. Our custom zero-copy GFA streaming parser (`src/gfa/mod.rs`) solves the I/O bottleneck entirely."
- **Verdict:** REFUTED
- **Evidence:**
  - `src/gfa/mod.rs:85-117` allocates heap `String` and `Vec<u8>` on every record (`parts.next().unwrap_or("").to_string()`, `.as_bytes().to_vec()`), and reads into a growing `String` via `reader.read_line(&mut line)` (line 76). It is not zero-copy.
  - `src/gfa/mod.rs:94,96` indexes `.as_bytes()[0]` on `parts.next().unwrap_or("+")`: if an `L` line has an empty orientation field (`""`), `unwrap_or("+")` returns `""` and `.as_bytes()[0]` panics at runtime with an index-out-of-bounds error.
  - `src/gfa/mod.rs:113-114` silently swallows integer parse errors on `W` lines with `.parse().unwrap_or(0)`, conflating invalid fields or `*` with coordinate `0`.
  - `grep -rn "gfa::\|stream_gfa" src/` shows `stream_gfa` is never called anywhere in the codebase.
- **Correction:** Remove the "zero-copy" and "solves the I/O bottleneck" claims; `src/gfa/mod.rs` is an uncalled, heap-allocating line-by-line GFA reader with a panic bug on malformed `L` lines.

## F5 -- C5: No `stdin`/`stdout` or VCF support exists; `.hal` HDF5 code and custom GFA tags exist

- **Claim:** C5 (`docs/STRATEGIC_PLAN.md:17`) -- "We will **never invent a proprietary binary format**. Our tools will communicate purely by streaming standard GFA (v1.0/v1.1) and standard MAF/VCF over `stdin`/`stdout`."
- **Verdict:** REFUTED
- **Evidence:**
  - `grep -rn "stdin\|stdout" src/` returns zero matches, and `grep -rni "vcf" src/` returns zero matches.
  - `src/hal/mod.rs`, `src/hal/hdf5_read.rs`, `src/hal/hdf5_write.rs`, and `horizon_wga --hal-info` (`src/bin/horizon_wga.rs:150-169`) implement and inspect the HDF5 `.hal` binary format.
  - `src/formats.rs:15,43,50,102` emits custom non-standard GFA tags/comments (`HB:i:`, `XS:i:`, `YS:i:`, `ST:Z:`, `# EG:Z:`).
- **Correction:** State that existing tools read/write MAF, GFA (with custom tags), and inspect `.hal` files via file-path arguments, with no `stdin`/`stdout` streaming or VCF support.

## F6 -- C6: No tool in `src/bin/` reads GFA/MAF from `stdin` or writes to `stdout`

- **Claim:** C6 (`docs/STRATEGIC_PLAN.md:23`) -- "Every tool must read GFA/MAF from `stdin` and write to `stdout`."
- **Verdict:** REFUTED
- **Evidence:**
  - `grep -rn "stdin\|stdout" src/` returns zero matches.
  - `src/bin/horizon_wga.rs:16-20,122` requires file paths `<target.fa> <query.fa>` and `--maf out.maf` / `--gfa out.gfa`.
  - `src/bin/horizon_merge.rs:12-35` requires `<maf_out> <gfa_out> <sp=path.maf>...` as file paths on disk.
  - `src/gfa/mod.rs:50-64` (`stream_gfa`) only accepts `path: P` where `P: AsRef<Path>` and calls `File::open(path.as_ref())?`.
- **Correction:** State that all current CLI binaries take file paths as arguments and do not support `stdin`/`stdout` streaming.

## F7 -- C7: No linear GFA streaming, bubble extraction, or path-stitching pipeline exists

- **Claim:** C7 (`docs/STRATEGIC_PLAN.md:27`) -- "We will stream the pangenome linearly and only invoke the E-graph engine when we hit a complex, tangled sub-graph ('hairball'). The E-graph resolves the optimal path for that bubble, and we stitch the resolved path back into the linear stream."
- **Verdict:** REFUTED
- **Evidence:**
  - `gfa::stream_gfa` (`src/gfa/mod.rs:50`) is never called anywhere in `src/`.
  - `egraph::saturate_synteny` is only invoked via `resolve_alignment_duplications` (`src/formats.rs:73-94`) on fixed 12-element chunks of 1000-bp binned pairwise alignments, not on GFA subgraphs or bubbles, and outputs only `# EG:Z:` comments (`src/formats.rs:98-103`) rather than stitching paths back into a GFA stream.
- **Correction:** Remove or mark as unimplemented future work; no GFA bubble-extraction or path-stitching code exists.

## F8 -- C8: Benchmarking and implementation are entirely WGA/MAF, not `horizon_clean`/`horizon_resolve`/`horizon_lift`

- **Claim:** C8 (`docs/STRATEGIC_PLAN.md:31`) -- "The HPRC produces massive graphs that researchers cannot use because standard tooling crashes or produces tangled paths. We will focus our benchmarking and release strategy entirely on `horizon_clean`, `horizon_resolve`, and `horizon_lift` first to capture the pangenome community immediately."
- **Verdict:** REFUTED
- **Evidence:**
  - None of `horizon_clean`, `horizon_resolve`, or `horizon_lift` exist in `src/bin/` or `Cargo.toml`.
  - All existing benchmarks in `docs/BENCHMARK_RESULTS.md:9-39` and `src/bin/` (`horizon_wga`, `horizon_merge`, `horizon_guide`, `horizon_eval`, `bench_chr1`) are pairwise and multi-species linear WGA/MAF benchmarks; the sole HPRC row in `docs/BENCHMARK_RESULTS.md:39` has no results (`*(Graph Eval)*`).
- **Correction:** State that all existing binaries and benchmarks target pairwise/multi-species WGA (`horizon_wga`, `horizon_merge`, `horizon_guide`), while `horizon_clean`, `horizon_resolve`, and `horizon_lift` are unimplemented.

## F9 -- C9: No code streams external `.gfa` files or scores `L` lines with BOCPD/MDL

- **Claim:** C9 (`docs/PANGENOME_SUITE.md:8`) -- "Stream an external `.gfa` file and evaluate every Link (`L` line) using the Horizon BOCPD/MDL scoring engine."
- **Verdict:** REFUTED
- **Evidence:**
  - `gfa::stream_gfa` (`src/gfa/mod.rs:50`) has zero callers in `src/`.
  - Neither `src/bocpd.rs` nor `src/mdl/` contains any function accepting `GfaLink` or scoring GFA `L` lines.
  - `horizon_clean` does not exist in `src/bin/`.
- **Correction:** Mark `horizon_clean` and `L`-line BOCPD/MDL scoring as unimplemented.

## F10 -- C10: `horizon_wga` and `src/mdl/dp.rs` only align linear sequences, not graphs

- **Claim:** C10 (`docs/PANGENOME_SUITE.md:16`) -- "Thread the query sequence through the graph, performing banded MDL dynamic programming across graph bifurcations."
- **Verdict:** REFUTED
- **Evidence:**
  - `banded_mdl_gap` (`src/mdl/dp.rs:81`) takes two linear byte slices (`seq_x: &[u8], seq_y: &[u8]`) and has no graph or bifurcation representation.
  - `src/bin/horizon_wga.rs:48-123` has no `--graph` flag and only accepts two FASTA files (`<target.fa> <query.fa>`).
- **Correction:** Mark sequence-to-graph alignment (`horizon_wga --graph`) and graph-bifurcation DP as unimplemented.

## F11 -- C11: `src/coord.rs` is not even compiled in `src/lib.rs`, and `horizon_lift` does not exist

- **Claim:** C11 (`docs/PANGENOME_SUITE.md:29`) -- "Traverse the GFA Walk (`W`) lines and Path (`P`) lines using Horizon's strict `coord::Interval` zero-based half-open arithmetic."
- **Verdict:** REFUTED
- **Evidence:**
  - `src/lib.rs:3-9,394` declares `bocpd`, `egraph`, `formats`, `hal`, `mdl`, `merge`, `seeding`, and `gfa`, but omits `pub mod coord;` completely — `src/coord.rs` is not compiled into the library.
  - `src/gfa/mod.rs:29-46` (`GfaPath`, `GfaWalk`) stores unparsed strings (`segment_names: String`, `walk: String`) and never references `coord::Interval`.
  - `horizon_lift` does not exist in `src/bin/`.
- **Correction:** Mark `horizon_lift` and `W`/`P`-line traversal via `coord::Interval` as unimplemented, and note that `src/coord.rs` is currently unreferenced dead code omitted from `src/lib.rs`.

| C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED |
