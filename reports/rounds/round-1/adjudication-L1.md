# Round 1 ADJUDICATION -- Correctness (L1)

Attack reviewed: `reports/rounds/round-1/attack-L1-correctness.md` (sealed digest `db1e0a921ab5`, verified).
Prior history: `scripts/round history "STRATEGIC_PLAN" --kind adjudication` returned no prior closed rounds.

## Finding-by-Finding Rulings

### F1 -- C1 (`horizon_clean`, `horizon_resolve`, `horizon_lift`, and Unix-pipe composition)
- **Verdict:** **NARROWED**
- **Ruling:** Diagnosis upheld on the factual gap (only `horizon_wga` exists in `src/bin/` alongside `horizon_merge`, `horizon_guide`, `horizon_eval`, and `bench_chr1`, and none of them read from `stdin` or write alignment streams to `stdout`), narrowed because C1 (`docs/STRATEGIC_PLAN.md:6`) uses progressive roadmap phrasing ("We are building...") rather than claiming the suite is already shipped. The proposed correction (explicitly listing which binaries exist today and that none compose via `stdin`/`stdout` pipes) is endorsed; the fix is free (text-only edit).

### F2 -- C2 (hardcoded affine gap penalties in `src/bocpd.rs` and `src/mdl/dp.rs`)
- **Verdict:** **UPHELD**
- **Ruling:** Both the diagnosis and the proposed correction hold: `docs/STRATEGIC_PLAN.md:9` makes a present-tense claim under `## 2. What We Are Doing Right` that Horizon replaces arbitrary affine gap penalties, whereas `src/bocpd.rs:4-7` hardcodes `GAP_OPEN_BITS = -4.5` and `GAP_EXTEND_BITS = -0.8`, and `src/mdl/dp.rs:105-107,129-146` uses them directly as affine gap-open/extend transition penalties. The proposed correction (describing BOCPD local GC composition priors combined with fixed `-4.5` / `-0.8` bit affine gap constants) is endorsed; the fix is free (text-only edit).

### F3 -- C3 (`src/egraph.rs` scope vs. nested structural paralogs and evolutionary trajectories)
- **Verdict:** **NARROWED**
- **Ruling:** Diagnosis upheld on the core defect (`src/egraph.rs:165-273` implements only four linear 1D rewrite rules on `seq`/`inv`/`dup`, cannot collapse >2-copy or non-tandem/interleaved paralogs, has no evolutionary trajectory model, and emits `# EG:Z:` comment annotations rather than de-cycled graphs), narrowed on the attacker's scope detail that `resolve_alignment_duplications` (`src/formats.rs:73-103`) with 1000-bp bins is the sole downstream user (`resolve_lcb_duplications` is also called by `src/merge.rs:206` with configurable `bin_bp`). The proposed correction is endorsed with that minor adjustment (`horizon_wga` and `horizon_merge`); the fix is free (text-only edit).

### F4 -- C4 (`src/gfa/mod.rs` zero-copy claim, unused status, and `L`-line panic bug)
- **Verdict:** **UPHELD**
- **Ruling:** Both the diagnosis and the proposed correction hold: `src/gfa/mod.rs:75-117` allocates `String` and `Vec<u8>` per record via `BufRead::read_line` (not zero-copy), panics with an out-of-bounds index on empty `L`-line orientation fields (`parts.next().unwrap_or("+").as_bytes()[0]` at lines 94 and 96), silently defaults `W`-line integer parse errors to `0` (lines 113-114), and has zero callers anywhere in `src/` despite being declared at `src/lib.rs:394`. Removing the "zero-copy" and "solves the I/O bottleneck" claims is endorsed; the fix is free (text-only edit).

### F5 -- C5 (no `stdin`/`stdout` or VCF support; `.hal` HDF5 code and GFA tags)
- **Verdict:** **NARROWED**
- **Ruling:** Diagnosis upheld on the missing `stdin`/`stdout` streaming and VCF support (and the active presence of `.hal` inspection/read/write code in `src/hal/` and `src/bin/horizon_wga.rs:150-169`), narrowed because C5 (`docs/STRATEGIC_PLAN.md:17`) is a future-tense commitment ("We will never invent a proprietary binary format") and neither reading/inspecting an existing external format (`.hal`) nor emitting standard GFA v1.0/v1.1 optional typed tags (`HB:i:`, `XS:i:`, `YS:i:`, `ST:Z:`) or `#` comments (`# EG:Z:`) constitutes inventing a proprietary binary format. The proposed correction (clarifying current file-path I/O, `.hal` support, and lack of `stdin`/`stdout` or VCF) is endorsed; the fix is free (text-only edit).

### F6 -- C6 (no binary reads GFA/MAF from `stdin` or writes to `stdout`)
- **Verdict:** **NARROWED**
- **Ruling:** Factual check upheld (`src/bin/horizon_wga.rs`, `src/bin/horizon_merge.rs`, and `src/gfa/mod.rs:50-64` all require file paths on disk, and `src/` has zero `stdin`/`stdout` streaming), narrowed because C6 (`docs/STRATEGIC_PLAN.md:23`) is a normative design requirement ("must") in a strategic mitigation section rather than a present-tense claim that current binaries already support pipes. Clarifying that current CLI binaries do not yet implement `stdin`/`stdout` streaming is endorsed; the fix is free (text-only edit).

### F7 -- C7 (no linear GFA streaming, bubble extraction, or path-stitching pipeline)
- **Verdict:** **NARROWED**
- **Ruling:** Factual check upheld (`gfa::stream_gfa` has zero callers in `src/`, and `egraph::saturate_synteny` is only called on 12-element chunks of binned alignment blocks in `src/formats.rs:73-94` and `src/merge.rs:206`, outputting `# EG:Z:` comments rather than stitched GFA paths), narrowed because C7 (`docs/STRATEGIC_PLAN.md:27`) is a future-tense mitigation proposal ("We will stream..."). Marking C7 explicitly as unimplemented future work is endorsed; the fix is free (text-only edit).

### F8 -- C8 (benchmarking and implementation are WGA/MAF, not `horizon_clean`/`horizon_resolve`/`horizon_lift`)
- **Verdict:** **NARROWED**
- **Ruling:** Factual check upheld (none of `horizon_clean`, `horizon_resolve`, or `horizon_lift` exist in `src/bin/` or `Cargo.toml`, and `docs/BENCHMARK_RESULTS.md:9-39` contains only pairwise/multi-species WGA/MAF benchmarks with an empty `*(Graph Eval)*` placeholder for HPRC), narrowed because C8 (`docs/STRATEGIC_PLAN.md:31`) states a future-tense release strategy ("We will focus... first") rather than claiming those benchmarks are already complete. Clarifying current WGA benchmark scope vs. unimplemented pangenome tools is endorsed; the fix is free (text-only edit).

### F9 -- C9 (`horizon_clean` and `L`-line BOCPD/MDL scoring)
- **Verdict:** **NARROWED**
- **Ruling:** Factual check upheld (`horizon_clean` does not exist in `src/bin/`, `gfa::stream_gfa` has zero callers, and neither `src/bocpd.rs` nor `src/mdl/` scores `GfaLink` `L` lines), narrowed because `docs/PANGENOME_SUITE.md:34` explicitly frames the document as a strategic roadmap rather than shipped functionality. Marking `horizon_clean` and `L`-line scoring as unimplemented roadmap items is endorsed; the fix is free (text-only edit).

### F10 -- C10 (`horizon_wga --graph` and graph-bifurcation DP)
- **Verdict:** **NARROWED**
- **Ruling:** Factual check upheld (`banded_mdl_gap` in `src/mdl/dp.rs:81` only accepts two linear `&[u8]` slices, and `src/bin/horizon_wga.rs:48-123` has no `--graph` flag), narrowed because `docs/PANGENOME_SUITE.md:34` frames the document as a strategic roadmap rather than shipped functionality. Marking `horizon_wga --graph` and graph DP as unimplemented is endorsed; the fix is free (text-only edit).

### F11 -- C11 (`src/coord.rs` uncompiled and `horizon_lift` unimplemented)
- **Verdict:** **NARROWED**
- **Ruling:** Factual check upheld (`src/coord.rs` is omitted from `src/lib.rs:3-9,394` and is uncompiled dead code, `GfaPath` and `GfaWalk` in `src/gfa/mod.rs:29-46` store unparsed strings, and `horizon_lift` does not exist in `src/bin/`), narrowed because `docs/PANGENOME_SUITE.md:34` frames the document as a strategic roadmap rather than shipped functionality (though citing `coord::Interval` as an active primitive is misleading without noting it is uncompiled). Marking `horizon_lift` as unimplemented and noting `src/coord.rs` is currently uncompiled is endorsed; the fix is free (text-only edit).

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | 3 of 4 named tools and pipe composition do not exist, but C1 uses progressive roadmap phrasing |
| F2 | UPHELD | present-tense claim contradicted by hardcoded `-4.5`/`-0.8` bit affine gap constants in `bocpd.rs` and `mdl/dp.rs` |
| F3 | NARROWED | 4-rule e-graph only emits `# EG:Z:` comments, but `merge.rs:206` also calls `resolve_lcb_duplications` |
| F4 | UPHELD | `src/gfa/mod.rs` allocates per record, panics on empty `L`-line orientations, and has zero callers |
| F5 | NARROWED | `stdin`/`stdout` and VCF are absent, but C5 is future-tense and `.hal` / GFA tags are not an invented proprietary format |
| F6 | NARROWED | current binaries require file paths on disk, but C6 states a normative design requirement |
| F7 | NARROWED | no GFA bubble extraction or path-stitching exists, but C7 is a future-tense mitigation proposal |
| F8 | NARROWED | all existing binaries and benchmarks are WGA/MAF, but C8 states a future-tense strategy |
| F9 | NARROWED | `horizon_clean` and `L`-line scoring are unimplemented, but `PANGENOME_SUITE.md` is explicitly a roadmap |
| F10 | NARROWED | `horizon_wga --graph` and graph DP are unimplemented, but `PANGENOME_SUITE.md` is explicitly a roadmap |
| F11 | NARROWED | `src/coord.rs` is uncompiled and `horizon_lift` is absent, but `PANGENOME_SUITE.md` is explicitly a roadmap |
