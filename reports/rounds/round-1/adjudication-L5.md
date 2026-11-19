# Round 1 Adjudication -- L5 (Inference)

Attack digest verified: `0535cbac162e` (`reports/rounds/round-1/attack-L5-inference.md`).
Prior adjudication history: `scripts/round history "" --kind adjudication` returned no prior closed rounds.
Method cards: `refs/methods/` does not exist in this repository; no method card governs these findings.

---

## F1 -- C1: Pipe composition is incompatible with graph indexing, out-of-order GFA records, and dual-input alignment

- **Verdict:** **UPHELD**
- **What survives:** Single-pass `stdin`/`stdout` Unix-pipe composition is structurally incompatible with GFA v1.0/v1.1 graph-global operations (`horizon_wga --query`, `horizon_resolve`) because GFA permits arbitrary record ordering, and three of the four named binaries do not exist.
- **Verification:** `docs/STRATEGIC_PLAN.md:6, 21-22` explicitly specifies `cat hprc.gfa | horizon_clean ... | horizon_wga --query novel_genome.fa | horizon_resolve > final.gfa`. GFA v1.0/v1.1 does not enforce `S`-before-`L` or topological ordering, so indexing all `S` lines and threading `L` bifurcations in a single-pass unidirectional pipe requires buffering the entire graph in RAM. In `src/bin/`, `horizon_clean`, `horizon_resolve`, and `horizon_lift` do not exist; `src/bin/horizon_wga.rs:15-17, 171-172` reads two FASTA file paths, and `src/bin/horizon_guide.rs:1` is a monolithic multi-species pipeline.
- **Diagnosis & Correction:** Diagnosis UPHELD. Correction UPHELD — replacing the pure-pipe claim with a two-phase/indexed graph-store architecture in `docs/STRATEGIC_PLAN.md` is a documentation edit (the fix is free).

## F2 -- C2: MDL bit-savings does not prevent paralogous or repetitive "hairball" over-alignment

- **Verdict:** **UPHELD**
- **What survives:** Horizon still uses fixed affine gap penalties and fixed thresholds, and MDL bit-savings against a random GC-matched null model abstains on non-homologous/divergent sequence rather than on high-identity repeats or paralogs.
- **Verification:** `src/bocpd.rs:4-7` hardcodes `GAP_OPEN_BITS = -4.5`, `GAP_EXTEND_BITS = -0.8`, `IDENTITY_PRIOR = 0.65`, and `MIN_POST_PROB = 0.15`; `src/mdl/mod.rs:62-64` hardcodes `band_w = 180`, `max_ext = 1500`, and `min_delta_bits = 10.0`. `joint_pair_bits` (`src/bocpd.rs:28-39`) computes log-odds relative to `p_null = pa * pb`, which abstains on random/divergent sequences (`src/lib.rs:105-131`), whereas repeats and paralogs share genuine sequence identity and score strongly positive. Hairball pruning in `src/bin/horizon_wga.rs:219-233` relies on `select_top_k_pairs` and `min_shared_seeds`.
- **Diagnosis & Correction:** Diagnosis UPHELD. Correction UPHELD — revising `docs/STRATEGIC_PLAN.md` to distinguish divergence abstention from repeat/paralog pruning is a documentation edit (the fix is free).

## F3 -- C3: Synteny e-graph rewrite rules do not model evolutionary trajectories or eliminate graph cycles

- **Verdict:** **UPHELD**
- **What survives:** `src/egraph.rs` is a bounded 1D term-rewriting pass over 12-block chunks of linear signed labels with fixed parsimony weights, not a phylogenetic or GFA cycle-resolution engine.
- **Verification:** `build_seq_expr` (`src/egraph.rs:134-161`) constructs a right-associated `seq` tree from `&[(&str, i8)]`. `saturate_synteny` (`src/egraph.rs:165-273`) applies four syntactic rules (`inv(inv(x)) = x`, `inv(seq(a,b)) = seq(inv(b),inv(a))`, `seq(a,a) = dup(a)`, and `seq` associativity), and `extract_min_cost` (`src/egraph.rs:291-306`) uses fixed weights (`1.0`, `+1.2`, `+1.3`, `+1.0`) to emit an `Expr` tree. `resolve_alignment_duplications` (`src/formats.rs:73-94`) chunks into 12-block windows with `max_iters = 4, max_nodes = 800` (`src/formats.rs:66`), which cannot collapse dispersed paralogs (`A B A`) or resolve GFA graph cycles.
- **Diagnosis & Correction:** Diagnosis UPHELD. Correction UPHELD — narrowing the claim in `docs/STRATEGIC_PLAN.md` to bounded 12-block tandem/inversion parsimony is a documentation edit (the fix is free).

## F4 -- C4: `src/gfa/mod.rs` is neither zero-copy nor a solution to graph RAM exhaustion

- **Verdict:** **UPHELD**
- **What survives:** `src/gfa/mod.rs` allocates owned `String` and `Vec<u8>` buffers for every parsed GFA record (and is unused in `src/lib.rs`), and line-by-line parsing does not eliminate the RAM required by random-access graph algorithms.
- **Verification:** `stream_gfa` (`src/gfa/mod.rs:75-117`) reads into a `String` via `read_line` and allocates owned `String`/`Vec<u8>` fields in `GfaSegment`, `GfaLink`, `GfaPath`, and `GfaWalk` (`src/gfa/mod.rs:14-46`). Moreover, `gfa` is not even declared as a module in `src/lib.rs:3-9` and `stream_gfa` has zero callers in `src/`.
- **Diagnosis & Correction:** Diagnosis UPHELD. Correction UPHELD — removing the "zero-copy" and RAM-exhaustion-solution claims in `docs/STRATEGIC_PLAN.md` is a documentation edit (the fix is free).

## F5 -- C5: Purely streaming standard GFA over `stdin`/`stdout` contradicts random-access and sidecar requirements

- **Verdict:** **NARROWED**
- **What survives:** Pure `stdin`/`stdout` streaming contradicts the same paragraph's optional sidecar-index requirement and is incompatible with random-access GFA operations, while `stream_gfa` (`src/gfa/mod.rs:50-73`) only opens filesystem `Path`s rather than `stdin`.
- **Verification & Overreach:** The `stdin`/`stdout` vs. sidecar-index contradiction (`docs/STRATEGIC_PLAN.md:17`), the impossibility of random-access graph queries on `stdin`, and `File::open` in `src/gfa/mod.rs:64` all hold. However, the attacker's secondary citation of `src/formats.rs:50, 102` (`# EG:Z:...` comments and `*`-sequence segments with `HB:i:1`) does *not* violate GFA v1.0/v1.1 syntax—comments (`#`), omitted sequences (`*` with `LN:i:`), and two-letter optional tags (`HB:i:1`) are valid standard GFA, even though `*`-sequence segments lack DNA bases for downstream sequence alignment.
- **Diagnosis & Correction:** Diagnosis NARROWED (drop the implication that `src/formats.rs:50, 102` emits non-standard GFA). Correction UPHELD — revising `docs/STRATEGIC_PLAN.md` to specify GFA/MAF/VCF as on-disk interchange formats with indexed/multi-pass access is a documentation edit (the fix is free).

## F6 -- C6: `stdin`/`stdout` GFA/MAF requirement is violated by every existing binary and impossible for dual-input alignment

- **Verdict:** **UPHELD**
- **What survives:** Every existing binary reads and writes filesystem paths rather than `stdin`/`stdout`, and multi-input alignment/merging (`horizon_wga`, `horizon_merge`) cannot receive multiple independent inputs on a single `stdin` stream without file arguments or multiplexing.
- **Verification:** `src/bin/horizon_wga.rs:121-127, 171-172, 206-217` takes `<target.fa> <query.fa>` and `--maf`/`--gfa` file paths; `src/bin/horizon_merge.rs:12, 20-35` takes `<maf_out> <gfa_out> <sp=path.maf>...`; `src/bin/horizon_guide.rs` reads FASTA files and writes output files.
- **Diagnosis & Correction:** Diagnosis UPHELD. Correction UPHELD — dropping the universal `stdin`/`stdout` mandate in `docs/STRATEGIC_PLAN.md` is a documentation edit (the fix is free).

## F7 -- C7: Linear streaming cannot detect hairballs or stitch resolved paths back into an already-emitted stream

- **Verdict:** **UPHELD**
- **What survives:** Single-pass linear streaming cannot detect cycles or tangles in arbitrary-order GFA files or retract/stitch already-emitted records in a Unix pipe.
- **Verification:** In GFA v1.0/v1.1, an `L` line near the end of a file can connect segments emitted at the start, forming a cycle or hairball after upstream lines have been written to `stdout`. Buffering until a connected component closes defeats linear streaming by requiring chromosome-scale component buffering in RAM. Additionally, `src/egraph.rs:134-161` takes a pre-built linear `&[(&str, i8)]`, not a GFA subgraph.
- **Diagnosis & Correction:** Diagnosis UPHELD. Correction UPHELD — replacing linear-stream bubble extraction with an explicit indexed/on-disk graph-partitioning pass in `docs/STRATEGIC_PLAN.md` is a documentation edit (the fix is free).

## F8 -- C8: Release strategy targets three non-existent tools while the codebase implements pairwise/progressive WGA

- **Verdict:** **UPHELD**
- **What survives:** Focusing immediate benchmarking and release "entirely on `horizon_clean`, `horizon_resolve`, and `horizon_lift` first" is unsupported by the repository, where none of those three binaries exist and all working binaries implement pairwise/progressive WGA.
- **Verification:** `src/bin/` contains `bench_chr1.rs`, `horizon_eval.rs`, `horizon_guide.rs`, `horizon_merge.rs`, and `horizon_wga.rs`. Neither `coord` nor `gfa` is declared in `src/lib.rs:3-9`.
- **Diagnosis & Correction:** Diagnosis UPHELD. Correction UPHELD — marking the pangenome binaries as unimplemented roadmap items and grounding immediate benchmarking in the existing WGA binaries in `docs/STRATEGIC_PLAN.md` is a documentation edit (the fix is free).

## F9 -- C9: `horizon_clean` cannot score GFA `L` lines via BOCPD/MDL in a stream

- **Verdict:** **UPHELD**
- **What survives:** GFA `L` lines carry only segment IDs, orientations, and overlap CIGARs (no DNA sequences), and Horizon's BOCPD/MDL engine scores pairwise sequence alignments rather than topological links, making single-pass `L`-line MDL scoring impossible.
- **Verification:** `GfaLink` (`src/gfa/mod.rs:21-27`) stores `from`, `from_orient`, `to`, `to_orient`, `overlap` with no sequence fields. `banded_mdl_gap` (`src/mdl/dp.rs:81`) takes `(seq_x: &[u8], seq_y: &[u8], ...)`. Even scoring alternative bubble paths requires an `S`-line sequence index because `L` lines can precede or appear gigabytes away from their `S` lines.
- **Diagnosis & Correction:** Diagnosis UPHELD. Correction UPHELD — specifying in `docs/PANGENOME_SUITE.md` that `horizon_clean` must index `S`-line sequences and `L`-line topology before pairwise MDL scoring of competing bubble paths is a documentation edit (the fix is free).

## F10 -- C10: Banded linear MDL DP (`BandGrid`) cannot run across graph bifurcations

- **Verdict:** **UPHELD**
- **What survives:** Horizon's existing `banded_mdl_gap` (`BandGrid`) and `bocpd_gc_profile` are strictly 1D/2D linear-sequence formulations and cannot traverse graph bifurcations or cycles without a new DAG/graph-HMM DP formulation, and `horizon_wga --graph` does not exist.
- **Verification:** `BandGrid` (`src/mdl/dp.rs:19-72`) indexes a 2D linear `(n, m)` grid around `center(i) = (i * m) / n`. `bocpd_gc_profile` (`src/bocpd.rs:58`) runs a linear left-to-right scan. `src/bin/horizon_wga.rs:32-120` has no `--graph` flag.
- **Diagnosis & Correction:** Diagnosis UPHELD. Correction UPHELD — replacing the claim in `docs/PANGENOME_SUITE.md` with an explicit graph-DP / POA formulation marked as unimplemented future work is a documentation edit (the fix is free).

## F11 -- C11: `W`/`P` line traversal with `coord::Interval` cannot work without segment lengths and mismatches `W`/`P` semantics

- **Verdict:** **UPHELD**
- **What survives:** `src/coord.rs` is not wired into `src/lib.rs`, `GfaWalk`/`GfaPath` lack segment lengths needed for coordinate projection, and `Interval::ucsc_start` implements UCSC reverse-strand arithmetic rather than GFA `W`/`P` forward-strand haplotype coordinates.
- **Verification:** `src/lib.rs:3-9` omits `coord`. `coord::Interval` (`src/coord.rs:4-38`) implements only `new`, `span`, and `ucsc_start`, with `ucsc_start` computing `src_size - end` on `'-'` strand. `GfaPath` and `GfaWalk` (`src/gfa/mod.rs:29-46`) store no segment lengths; GFA v1.1 `W` lines use 0-based half-open forward-strand haplotype coordinates regardless of `>`/`<` segment orientation, and GFA v1.0 `P` lines carry no `src_size` or haplotype bounds.
- **Diagnosis & Correction:** Diagnosis UPHELD. Correction NARROWED on cost — updating `docs/PANGENOME_SUITE.md` to distinguish GFA `W`/`P` coordinates from UCSC coordinates is free, whereas wiring `coord` into `src/lib.rs` and building a global `S`-line length index requires implementation work.

---

## Verdict Summary

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | Single-pass Unix pipes cannot support random-access GFA graph queries with arbitrary record order, and 3 of 4 named binaries do not exist |
| F2 | UPHELD | Fixed affine gap penalties/thresholds remain hardcoded, and MDL abstains on random/divergent pairs rather than high-identity repeats/paralogs |
| F3 | UPHELD | `src/egraph.rs` runs 4 syntactic rules on 12-block linear chunks with fixed weights and emits an `Expr` tree, not a phylogenetic or cycle-free GFA graph |
| F4 | UPHELD | `src/gfa/mod.rs` heap-allocates owned `String`/`Vec<u8>` per record, is unused in `src/lib.rs`, and line streaming does not solve graph RAM exhaustion |
| F5 | NARROWED | `stdin`/`stdout` streaming contradicts sidecar indexing and random-access graph ops (`File::open` used), though `# EG:Z:` and `*`/`HB:i:1` are valid GFA syntax |
| F6 | UPHELD | All existing binaries read/write filesystem paths, and dual-/multi-input tools cannot receive separate inputs from a single `stdin` without file args |
| F7 | UPHELD | Arbitrary-order GFA records prevent single-pass hairball detection or retracting/stitching already-emitted pipe records |
| F8 | UPHELD | `horizon_clean`, `horizon_resolve`, and `horizon_lift` do not exist in `src/bin/`; all working binaries implement pairwise/progressive WGA |
| F9 | UPHELD | GFA `L` lines contain no DNA sequences, and `banded_mdl_gap` scores pairwise sequence alignments rather than topological graph links |
| F10 | UPHELD | `BandGrid` and `bocpd_gc_profile` are strictly linear 1D/2D formulations with no graph-DP/POA implementation, and `--graph` does not exist |
| F11 | UPHELD | `coord` is absent from `src/lib.rs`, `GfaWalk`/`GfaPath` lack segment lengths, and `Interval::ucsc_start` assumes UCSC reverse-strand semantics |
