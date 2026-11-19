# Round 1 ADJUDICATION of L7 (Reproduction) -- Strategic Plan & Pangenome Suite Roadmap

Adjudicating `reports/rounds/round-1/attack-L7-reproduction.md` (sha256 prefix `2a52c2bd0d7b`, verified) against `docs/STRATEGIC_PLAN.md`, `docs/PANGENOME_SUITE.md`, `Cargo.toml`, and `src/` at commit `1ebc165`.

Independent verification executed:
- `cargo check --all-targets`: exited 0 with 2 unused-import warnings in `src/gfa/mod.rs:6-7`.
- `cargo test`: exited 0 (28 passed, 0 failed, 1 ignored in `src/lib.rs`; 0 tests in `src/gfa/mod.rs` and all 5 `src/bin/*.rs` binaries).
- `cargo run --bin horizon_clean`, `horizon_resolve`, `horizon_lift`: all failed with `error: no bin target named ...`.
- `echo -e "H\tVN:Z:1.1\nS\ts1\tACGT" | cargo run --bin horizon_wga -- --graph --query novel_genome.fa`: failed with `Error: No such file or directory (os error 2)`.
- `grep -rn "stream_gfa" src/` returned only `src/gfa/mod.rs:50`; `grep -rn "VCF\|vcf" src/` and `grep -rn "stdin" src/` returned 0 matches; `grep -rn "mod gfa\|mod coord" src/` returned only `src/lib.rs:394:pub mod gfa;`.

Prior rounds & method cards:
- Prior round history: 0 closed rounds (attacker's check confirmed; `./scripts/round` is not present inside the `horizon-aligner/` subdirectory).
- Method cards (`refs/methods/`): F1–F11 are codebase/CLI reproduction checks rather than statistical disputes, so no method card governs them.

---

## Rulings

### F1 -- C1 (`horizon_clean`, `horizon_resolve`, `horizon_lift` & Unix pipes)
- **Verdict:** **UPHELD**
- **What survives:** `src/bin/` contains only `bench_chr1.rs`, `horizon_eval.rs`, `horizon_guide.rs`, `horizon_merge.rs`, and `horizon_wga.rs`; the three named pangenome binaries do not exist, and none of the five existing binaries compose via `stdin`/`stdout` pipes (`horizon_guide.rs` is a monolithic driver, while `horizon_wga.rs` and `horizon_merge.rs` take file paths).
- **Repair ruling:** Endorsed. Removing the missing binaries and Unix-pipe composability from present-tense claims and listing the five existing file-based CLIs is accurate; the fix is free (doc-only edit to `docs/STRATEGIC_PLAN.md`).

### F2 -- C2 (MDL bit-savings replacing affine gap penalties & scoring matrices)
- **Verdict:** **NARROWED**
- **What survives:** The defect is real on affine gap penalties and hairball testing—`src/bocpd.rs:4-5` hard-codes `GAP_OPEN_BITS = -4.5` and `GAP_EXTEND_BITS = -0.8` (used in `src/mdl/dp.rs:105-107, 130`), and `cargo test` has no graph "hairball" over-alignment test—but the attacker's blanket `REFUTED` verdict overstates the scope because `joint_pair_bits` (`src/bocpd.rs:28-39`) *does* replace traditional substitution scoring matrices with BOCPD base-frequency bit-savings.
- **Repair ruling:** Endorsed. The attacker's own proposed correction ("fixed affine gap penalties ... combined with BOCPD base-frequency priors, and no graph 'hairball' benchmark has been reproduced") already captures the exact narrowed scope; the fix is free (doc-only edit).

### F3 -- C3 (E-graph congruence closure on HLA / structural paralogs)
- **Verdict:** **UPHELD**
- **What survives:** The attacker's `NARROWED` finding is verified: `src/egraph.rs:375-384` tests `egraph_collapses_inverted_tandem_duplication` on a 5-element synthetic tuple, and `resolve_alignment_duplications` (`src/formats.rs:73-94`) buckets pairwise `ContigAlignment` blocks into `bin_bp` windows and emits `# EG:Z:...` tags (`src/formats.rs:102`, minor citation nit: attacker cited line `105`), with no HLA dataset or GFA graph resolution.
- **Repair ruling:** Endorsed. Narrowing C3 to synthetic block sequences and binned pairwise LCBs is accurate; the fix is free (doc-only edit).

### F4 -- C4 (Zero-copy GFA streaming parser `src/gfa/mod.rs`)
- **Verdict:** **UPHELD**
- **What survives:** `src/gfa/mod.rs:50-124` (`stream_gfa`) allocates a `String` line buffer via `read_line` and owned `String`/`Vec<u8>` fields per record, is never called anywhere in `src/`, has 0 unit tests and 2 unused-import warnings (`lines 6-7`), will panic on an `L` line with an empty orientation field (`lines 94, 96`), and has no 12GB pangenome benchmark against `vg`/`odgi`.
- **Repair ruling:** Endorsed. Removing "zero-copy" and "solves the I/O bottleneck" and describing `src/gfa/mod.rs` as an untested, unused, allocating line-by-line reader is accurate; the fix is free (doc-only edit).

### F5 -- C5 (No proprietary binary format; pure GFA/MAF/VCF over `stdin`/`stdout`)
- **Verdict:** **UPHELD**
- **What survives:** `src/` has 0 occurrences of `VCF`/`vcf` or `stdin`, while `src/hal/` actively implements HDF5 `.hal` reading and writing (with 5 passing tests in `cargo test` and `--hal-info` in `src/bin/horizon_wga.rs:109`).
- **Repair ruling:** Endorsed. Stating that existing tools communicate via on-disk FASTA/MAF/GFA/HAL files without VCF or `stdin`/`stdout` streaming is accurate; the fix is free (doc-only edit).

### F6 -- C6 (Every tool reads GFA/MAF from `stdin` and writes to `stdout`)
- **Verdict:** **UPHELD**
- **What survives:** None of the 5 binaries in `src/bin/` read GFA/MAF from `stdin`, and `horizon_wga -- --graph --query novel_genome.fa` fails with `No such file or directory (os error 2)` because `src/bin/horizon_wga.rs:113-116, 125-126` treats unknown flags as positional FASTA paths.
- **Repair ruling:** Endorsed. Dropping `stdin`/`stdout` as an implemented property and documenting file-path arguments is accurate; the fix is free (doc-only edit).

### F7 -- C7 (Linear pangenome streaming with bubble extraction into E-graphs)
- **Verdict:** **UPHELD**
- **What survives:** `stream_gfa` is uncalled and `src/egraph.rs` has no GFA I/O, bubble extraction, or path-stitching logic; `egraph` is only invoked on binned pairwise `ContigAlignment` blocks (`src/formats.rs:63-94`).
- **Repair ruling:** Endorsed. Reclassifying C7 as an unimplemented future design rather than a current capability is accurate; the fix is free (doc-only edit).

### F8 -- C8 (HPRC benchmarking & release strategy on `horizon_clean`, `horizon_resolve`, `horizon_lift`)
- **Verdict:** **UPHELD**
- **What survives:** None of `horizon_clean`, `horizon_resolve`, or `horizon_lift` exist in `Cargo.toml` or `src/bin/`, and no HPRC graph benchmarks exist in the repository (`bench_chr1.rs` and `horizon_eval.rs` target pairwise FASTA/MAF).
- **Repair ruling:** Endorsed. Noting that none of the three named binaries exist yet is accurate; the fix is free (doc-only edit).

### F9 -- C9 (`horizon_clean` BOCPD/MDL scoring of GFA `L` lines)
- **Verdict:** **UPHELD**
- **What survives:** `horizon_clean` does not exist, `stream_gfa` is unused/untested, and `bocpd_gc_profile` / `banded_mdl_gap` only accept `&[u8]` nucleotide slices rather than `GfaLink` records.
- **Repair ruling:** Endorsed. Marking C9 as unimplemented in the `docs/PANGENOME_SUITE.md` roadmap is accurate; the fix is free (doc-only edit).

### F10 -- C10 (`horizon_wga --graph` banded MDL DP across graph bifurcations)
- **Verdict:** **UPHELD**
- **What survives:** `horizon_wga` has no `--graph` or `--query` flag (`src/bin/horizon_wga.rs:49-117`), and `BandGrid` / `banded_mdl_gap` (`src/mdl/dp.rs:19-72, 81`) is strictly a 2D linear-sequence DP grid with no graph bifurcation support.
- **Repair ruling:** Endorsed. Marking C10 as unimplemented is accurate; the fix is free (doc-only edit).

### F11 -- C11 (`horizon_lift` GFA `W`/`P` traversal with `coord::Interval`)
- **Verdict:** **UPHELD**
- **What survives:** `horizon_lift` does not exist, `src/coord.rs` is not declared as a module in `src/lib.rs` (so `coord::Interval` is uncompiled dead code), and `src/gfa/mod.rs:100-117` stores `P` and `W` lines as raw `String`s.
- **Repair ruling:** Endorsed. Marking C11 as unimplemented and noting `src/coord.rs` is uncompiled dead code is accurate; the fix is free (doc-only edit).

---

## Summary Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `horizon_clean`, `horizon_resolve`, and `horizon_lift` do not exist and existing binaries use file paths, not Unix pipes; doc fix is free |
| F2 | NARROWED | `GAP_OPEN_BITS`/`GAP_EXTEND_BITS` are hard-coded and no hairball test exists, but `joint_pair_bits` does replace substitution matrices with BOCPD bit-savings; narrowed correction is free |
| F3 | UPHELD | Attacker's NARROWED verdict holds: `src/egraph.rs` is only tested on synthetic tuples and binned pairwise LCBs, not HLA or GFA; doc fix is free |
| F4 | UPHELD | `src/gfa/mod.rs` is allocating (not zero-copy), unused, untested, has 2 warnings, and lacks any 12GB benchmark; doc fix is free |
| F5 | UPHELD | VCF and `stdin`/`stdout` do not exist in `src/`, while `src/hal/` implements HDF5 `.hal` I/O; doc fix is free |
| F6 | UPHELD | None of the 5 binaries read GFA/MAF from `stdin`; `horizon_wga --graph --query` fails as a missing file path; doc fix is free |
| F7 | UPHELD | Linear GFA streaming with bubble extraction into E-graphs and path stitching is unimplemented; doc fix is free |
| F8 | UPHELD | None of the three named binaries or any HPRC graph benchmarks exist in `Cargo.toml` or `src/bin/`; doc fix is free |
| F9 | UPHELD | `horizon_clean` and BOCPD/MDL scoring of GFA `L` lines are unimplemented; doc fix is free |
| F10 | UPHELD | `horizon_wga --graph` does not exist and `banded_mdl_gap` is strictly linear 2D DP; doc fix is free |
| F11 | UPHELD | `horizon_lift` does not exist and `src/coord.rs` is uncompiled dead code not declared in `src/lib.rs`; doc fix is free |
