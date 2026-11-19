# Round 1 Attack -- L7 (Reproduction)

Checked against `/home/rcb112/project/gemini/bench/horizon-aligner/` at commit `1ebc165`.
Prior round check: `./scripts/round history ""` returned 0 closed rounds (`no prior mention of '' in 1 round directories`).
Commands executed:
- `cargo check --all-targets`: exited 0 (with 2 unused-import warnings in `src/gfa/mod.rs:6-7`).
- `cargo test`: exited 0 (28 passed, 0 failed, 1 ignored in `src/lib.rs`; 0 tests in `src/gfa/mod.rs` and all 5 `src/bin/*.rs` binaries).
- `cargo run --bin horizon_clean`, `horizon_resolve`, `horizon_lift`: all failed with `error: no bin target named ...`.
- `echo -e "H\tVN:Z:1.1\nS\ts1\tACGT" | cargo run --bin horizon_wga -- --graph --query novel_genome.fa`: failed with `Error: No such file or directory (os error 2)`.

---

## F1 -- C1: `horizon_clean`, `horizon_resolve`, and `horizon_lift` do not exist, and existing binaries do not compose via Unix pipes

- **Claim:** C1 (`docs/STRATEGIC_PLAN.md`: suite of `horizon_clean`, `horizon_wga`, `horizon_resolve`, `horizon_lift` composing via Unix pipes).
- **Verdict:** REFUTED
- **Evidence:**
  - `src/bin/` contains only `bench_chr1.rs`, `horizon_eval.rs`, `horizon_guide.rs`, `horizon_merge.rs`, and `horizon_wga.rs`. Running `cargo run --bin horizon_clean`, `cargo run --bin horizon_resolve`, and `cargo run --bin horizon_lift` each fails with `error: no bin target named ...`.
  - None of the existing binaries compose via `stdin`/`stdout` pipes: `src/bin/horizon_wga.rs:120-126, 171-172, 206-217, 274-285` requires `<target.fa> <query.fa>` file paths, writes MAF/GFA to `--maf`/`--gfa` files, and prints summary metrics to `stdout`; `src/bin/horizon_merge.rs:11-36` requires `<ref_species> <bin_bp> <maf_out> <gfa_out> <sp=path.maf>...` file paths. Meanwhile, `src/bin/horizon_guide.rs:1` is a monolithic multi-species pipeline driver.
- **Correction:** Remove `horizon_clean`, `horizon_resolve`, `horizon_lift`, and Unix-pipe composability from present-tense claims; state that only `horizon_wga`, `horizon_merge`, `horizon_guide`, `horizon_eval`, and `bench_chr1` exist as file-based CLIs.

## F2 -- C2: Affine gap penalties are hard-coded in `src/bocpd.rs` and `src/mdl/dp.rs`, not replaced

- **Claim:** C2 (`docs/STRATEGIC_PLAN.md`: replacing affine gap penalties and scoring matrices with MDL bit-savings, inherently preventing the "hairball" over-alignment problem).
- **Verdict:** REFUTED
- **Evidence:**
  - `src/bocpd.rs:4-7` defines hard-coded affine gap penalties `GAP_OPEN_BITS: f64 = -4.5` and `GAP_EXTEND_BITS: f64 = -0.8` (alongside fixed `IDENTITY_PRIOR: f64 = 0.65` and `MIN_POST_PROB: f64 = 0.15`), which `src/mdl/dp.rs:105-107` converts directly into DP gap costs (`go = GAP_OPEN_BITS * ln2`, `ge = GAP_EXTEND_BITS * ln2`).
  - `cargo test` (`src/lib.rs:43-393`) contains no pangenome graph or "hairball" over-alignment test; only pairwise sequence alignment tests exist.
- **Correction:** State that the DP uses fixed affine gap penalties (`GAP_OPEN_BITS = -4.5`, `GAP_EXTEND_BITS = -0.8`) combined with BOCPD base-frequency priors, and no graph "hairball" benchmark has been reproduced.

## F3 -- C3: E-graph congruence closure is only tested on a 5-block synthetic tuple and pairwise LCB bins, not HLA or GFA graphs

- **Claim:** C3 (`docs/STRATEGIC_PLAN.md`: using congruence-closure to resolve complex, nested structural paralogs like the HLA region without cyclical graphs).
- **Verdict:** NARROWED
- **Evidence:**
  - `src/egraph.rs:1-385` implements `EGraph` with bounded `saturate_synteny` (`max_iters=4, max_nodes=800` in `src/formats.rs:66`) and passes `egraph_collapses_inverted_tandem_duplication` (`src/egraph.rs:375-384`) on a 5-element synthetic slice `[("A", 1), ("C", -1), ("B", -1), ("B", -1), ("D", 1)]`, plus `resolve_dups_and_merge_pairwise_mafs_work` (`src/lib.rs:281-298`).
  - `resolve_alignment_duplications` (`src/formats.rs:73-94`) only buckets pairwise `ContigAlignment` blocks into `bin_bp` windows (`chunks(12)`) and emits `# EG:Z:...` comment tags (`src/formats.rs:105`). No HLA test, real structural-paralog dataset, or GFA graph resolution exists.
- **Correction:** Narrow C3 to state that `src/egraph.rs` demonstrates bounded congruence-closure on synthetic block sequences and binned pairwise LCBs, not HLA or GFA graph resolution.

## F4 -- C4: `src/gfa/mod.rs` is neither zero-copy, nor called anywhere, nor tested, nor benchmarked on 12GB graphs

- **Claim:** C4 (`docs/STRATEGIC_PLAN.md`: `vg` and `odgi` crash on 12GB pangenomes; custom zero-copy GFA streaming parser `src/gfa/mod.rs` solves the I/O bottleneck).
- **Verdict:** REFUTED
- **Evidence:**
  - `src/gfa/mod.rs:50-124` (`stream_gfa`) is not zero-copy: it allocates a `String` line buffer via `reader.read_line(&mut line)` (`lines 75-76`) and allocates owned `String` and `Vec<u8>` on every `S`, `L`, `P`, and `W` line (`lines 87-88, 93-97, 102-104, 110-115`).
  - `grep -rn "stream_gfa" src/` shows `stream_gfa` is never called anywhere in `src/`.
  - `src/gfa/mod.rs` has 0 unit tests, `cargo check` emits 2 unused-import warnings (`lines 6-7`), and `lines 94, 96` (`parts.next().unwrap_or("+").as_bytes()[0]`) will panic on an `L` line with an empty orientation field.
  - No 12GB pangenome benchmark or `vg`/`odgi` comparison exists in the repository.
- **Correction:** Remove the "zero-copy" and "solves the I/O bottleneck" claims; `src/gfa/mod.rs` is an untested, unused, allocating line-by-line GFA reader.

## F5 -- C5: VCF and `stdin`/`stdout` GFA/MAF streaming do not exist, while HDF5 `.hal` I/O is implemented

- **Claim:** C5 (`docs/STRATEGIC_PLAN.md`: never invent a proprietary binary format; tools communicate purely by streaming standard GFA and MAF/VCF over `stdin`/`stdout`).
- **Verdict:** REFUTED
- **Evidence:**
  - `grep -rn "VCF\|vcf" src/` returns 0 matches; VCF is unsupported.
  - `grep -rn "stdin" src/` returns 0 matches; `src/gfa/mod.rs:64` only opens a `File::open(path.as_ref())`, and `src/bin/horizon_wga.rs` / `horizon_merge.rs` only read/write files on disk.
  - `src/hal/` (`src/hal/mod.rs`, `hdf5_read.rs`, `hdf5_write.rs`) actively implements HDF5 `.hal` reading and writing, contradicting the claimed pure GFA/MAF/VCF architecture.
- **Correction:** State that the existing tools communicate via on-disk FASTA/MAF/GFA/HAL files, with no VCF or `stdin`/`stdout` streaming support.

## F6 -- C6: No binary reads GFA or MAF from `stdin`

- **Claim:** C6 (`docs/STRATEGIC_PLAN.md`: every tool must read GFA/MAF from `stdin` and write to `stdout`).
- **Verdict:** REFUTED
- **Evidence:**
  - Running `echo -e "H\tVN:Z:1.1\nS\ts1\tACGT" | cargo run --bin horizon_wga -- --graph --query novel_genome.fa` failed with `Error: No such file or directory (os error 2)` because `src/bin/horizon_wga.rs:113-116, 125-126` treats unknown flags as positional FASTA file paths.
  - All 5 binaries in `src/bin/` (`horizon_wga.rs`, `horizon_merge.rs`, `horizon_guide.rs`, `horizon_eval.rs`, `bench_chr1.rs`) read from file paths via `read_fasta_contigs` or `fs::read_to_string`; none read GFA/MAF from `stdin`.
- **Correction:** Drop the `stdin`/`stdout` requirement as an implemented property; document that current binaries require file-path arguments.

## F7 -- C7: Linear pangenome streaming with bubble extraction into E-graphs and path stitching is unimplemented

- **Claim:** C7 (`docs/STRATEGIC_PLAN.md`: stream pangenome linearly, invoke E-graph only on "hairball" bubbles, stitch resolved path back into linear stream).
- **Verdict:** REFUTED
- **Evidence:**
  - `stream_gfa` (`src/gfa/mod.rs:50`) is never called, and `src/egraph.rs` has no GFA input/output, no bubble extraction, and no path-stitching logic.
  - The only callers of `egraph` (`src/formats.rs:63-94`) run on binned pairwise `ContigAlignment` blocks and emit `# EG:Z:...` comment lines (`src/formats.rs:105`).
- **Correction:** Reclassify C7 as an unimplemented future design, not a current capability.

## F8 -- C8: `horizon_clean`, `horizon_resolve`, and `horizon_lift` do not exist for HPRC benchmarking or release

- **Claim:** C8 (`docs/STRATEGIC_PLAN.md`: focus benchmarking and release strategy on `horizon_clean`, `horizon_resolve`, and `horizon_lift` first for HPRC graphs).
- **Verdict:** REFUTED
- **Evidence:**
  - `cargo run --bin horizon_clean`, `horizon_resolve`, and `horizon_lift` all fail (`no bin target named ...`).
  - No HPRC graph benchmarks exist in `src/bin/` or `tests/`; existing benchmarks (`src/bin/bench_chr1.rs`, `src/bin/horizon_eval.rs`) target Plasmodium Chr1 FASTA/MAF pairs.
- **Correction:** Note that none of the three named binaries exist yet in `Cargo.toml` or `src/bin/`.

## F9 -- C9: `horizon_clean` and BOCPD/MDL scoring of GFA `L` lines do not exist

- **Claim:** C9 (`docs/PANGENOME_SUITE.md`: stream an external `.gfa` file and evaluate every Link (`L` line) using the Horizon BOCPD/MDL scoring engine).
- **Verdict:** REFUTED
- **Evidence:**
  - `horizon_clean` does not exist in `src/bin/`.
  - `src/gfa/mod.rs:50` (`stream_gfa`) is unused and untested, and `src/bocpd.rs:58` / `src/mdl/dp.rs:81` only accept `&[u8]` nucleotide slices, not `GfaLink` records.
- **Correction:** Mark C9 as unimplemented.

## F10 -- C10: `horizon_wga --graph` does not exist and `banded_mdl_gap` is strictly linear 2D DP

- **Claim:** C10 (`docs/PANGENOME_SUITE.md`: `horizon_wga --graph` threads query sequence through graph, performing banded MDL DP across graph bifurcations).
- **Verdict:** REFUTED
- **Evidence:**
  - `cargo run --bin horizon_wga -- --graph --query novel_genome.fa` fails with `Error: No such file or directory (os error 2)` (`src/bin/horizon_wga.rs:49-117` has no `--graph` or `--query` flag).
  - `src/mdl/dp.rs:19-72, 81` (`BandGrid`, `banded_mdl_gap`) is strictly a 2D linear-sequence-to-linear-sequence DP grid (`n` x `m`) with no graph bifurcation support.
- **Correction:** Mark C10 as unimplemented.

## F11 -- C11: `horizon_lift` does not exist and `src/coord.rs` is dead, uncompiled code

- **Claim:** C11 (`docs/PANGENOME_SUITE.md`: `horizon_lift` traverses GFA `W` and `P` lines using Horizon's strict `coord::Interval` zero-based half-open arithmetic).
- **Verdict:** REFUTED
- **Evidence:**
  - `horizon_lift` does not exist (`cargo run --bin horizon_lift` fails with `no bin target named horizon_lift`).
  - `src/coord.rs` defines `pub struct Interval`, but `mod coord;` is not declared in `src/lib.rs:1-394`, so `coord::Interval` is not even compiled into `horizon_aligner`.
  - `src/gfa/mod.rs:100-117` stores `P` and `W` lines as raw `String`s without using `coord::Interval`.
- **Correction:** Mark C11 as unimplemented and note that `src/coord.rs` is currently uncompiled dead code.

---

## Claim Verdicts

| C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| REFUTED | REFUTED | NARROWED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED |
