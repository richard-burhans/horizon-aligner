# Round 1 Attack -- L5 (Inference)

## F1 -- C1: Pipe composition is incompatible with graph indexing, out-of-order GFA records, and dual-input alignment

- **Claim:** C1 (`docs/STRATEGIC_PLAN.md` / `ATTACKED-STRATEGIC_PLAN.md:6, 21-22`)
- **Verdict:** REFUTED
- **Concrete flaw:** The inferred Unix-pipe architecture (`cat hprc.gfa | horizon_clean ... | horizon_wga --query novel_genome.fa | horizon_resolve > final.gfa`) cannot compose via single-pass `stdin`/`stdout` streams. First, `horizon_wga --query` must build a spaced-seed index over all GFA `S` lines and thread the query across `L`-line bifurcations (`ATTACKED-PANGENOME_SUITE.md:15-16`), which requires random access to the entire graph topology and segment sequences before alignment can proceed. Because GFA v1.0/v1.1 does not enforce topological or `S`-before-`L` ordering, a unidirectional pipe cannot satisfy random-access graph queries without buffering the entire graph in RAM (reproducing the `vg`/`odgi` memory bottleneck). Second, `horizon_clean`, `horizon_resolve`, and `horizon_lift` do not exist in `src/bin/`, while the existing `horizon_wga` (`src/bin/horizon_wga.rs:15-17, 171-172`) only reads two FASTA file paths and `horizon_guide` (`src/bin/horizon_guide.rs:1`) is itself a monolithic multi-species pipeline.
- **Evidence:** `ATTACKED-STRATEGIC_PLAN.md:21-23`; `ATTACKED-PANGENOME_SUITE.md:15-16`; `src/bin/horizon_wga.rs:15-17, 121-127, 171-172`; `src/bin/horizon_guide.rs:1`.
- **Correction:** Replace the claim of "compose perfectly via Unix pipes" with an explicit two-phase architecture: an indexed/disk-backed graph store (or multi-pass file reader) for graph-global operations (`horizon_wga --graph`, `horizon_resolve`, `horizon_lift`), restricting pipe streaming to linear MAF/VCF/summary outputs, and remove references to unimplemented binaries as current capabilities.

## F2 -- C2: MDL bit-savings does not prevent paralogous or repetitive "hairball" over-alignment

- **Claim:** C2 (`docs/STRATEGIC_PLAN.md` / `ATTACKED-STRATEGIC_PLAN.md:9`)
- **Verdict:** REFUTED
- **Concrete flaw:** The claim makes two false inferences. First, Horizon does not replace affine gap penalties and fixed parameters: `src/bocpd.rs:4-7` hardcodes `GAP_OPEN_BITS = -4.5`, `GAP_EXTEND_BITS = -0.8`, `IDENTITY_PRIOR = 0.65`, and `MIN_POST_PROB = 0.15`, and `src/mdl/mod.rs:62-67` hardcodes `band_w = 180`, `max_ext = 1500`, and `min_delta_bits = 10.0`. Second, MDL bit-savings (`joint_pair_bits` in `src/bocpd.rs:28-39`) tests whether two sequences are more similar than a random GC-matched null model (`p_null = pa * pb`). It abstains on unrelated/divergent sequences (`src/lib.rs:105-131`), **not** on repetitive elements or structural paralogs, which genuinely share high sequence identity and therefore yield strongly positive bit-savings across all paralogous copies. MDL alone would *produce* hairballs on repeats, which is why `horizon_wga` must explicitly prune pairs with `top_k`, `min_shared_seeds`, `fallback_max_occ`, and downstream `egraph.rs`.
- **Evidence:** `src/bocpd.rs:4-7, 28-39`; `src/mdl/mod.rs:62-67`; `src/lib.rs:105-131`; `src/bin/horizon_wga.rs:219-233`.
- **Correction:** State that MDL/BOCPD uses fixed affine gap bit-penalties (`-4.5` open, `-0.8` extend) and an identity prior (`0.65`) to abstain on *non-homologous/divergent* regions, while repetitive/paralogous over-alignment is controlled separately by seed-count/top-k filtering and downstream deduplication.

## F3 -- C3: Synteny e-graph rewrite rules do not model evolutionary trajectories or eliminate graph cycles

- **Claim:** C3 (`docs/STRATEGIC_PLAN.md` / `ATTACKED-STRATEGIC_PLAN.md:10`)
- **Verdict:** REFUTED
- **Concrete flaw:** The e-graph implementation in `src/egraph.rs` operates on a 1D linear sequence of signed block labels (`build_seq_expr`, `src/egraph.rs:134-161`) using four purely syntactic term-rewriting rules (`inv(inv(x)) = x`, `inv(seq(a,b)) = seq(inv(b),inv(a))`, `seq(a,a) = dup(a)`, and `seq` associativity in `src/egraph.rs:165-273`) with arbitrary fixed extraction weights (`blk: 1.0`, `seq: +1.0`, `inv: +1.2`, `dup: +1.3` in `src/egraph.rs:291-306`). It has no phylogenetic model, cannot collapse non-adjacent or dispersed paralogs (`A B A`), runs on fixed 12-block chunks (`src/formats.rs:82`) with `max_iters = 4, max_nodes = 800` (`src/formats.rs:66`), and produces a 1D `Expr` tree rather than a cycle-free GFA graph. It therefore neither models evolutionary trajectories nor resolves cyclical GFA topologies.
- **Evidence:** `src/egraph.rs:14-25, 134-161, 165-273, 277-327`; `src/formats.rs:63-94`.
- **Correction:** Narrow the claim to: "Bounded equality saturation (`src/egraph.rs`) collapses adjacent tandem duplications and nested inversions within 12-block windows of linear pairwise synteny chains using parsimony weights."

## F4 -- C4: `src/gfa/mod.rs` is neither zero-copy nor a solution to graph RAM exhaustion

- **Claim:** C4 (`docs/STRATEGIC_PLAN.md` / `ATTACKED-STRATEGIC_PLAN.md:11`)
- **Verdict:** REFUTED
- **Concrete flaw:** First, `src/gfa/mod.rs` is not zero-copy: `stream_gfa` allocates a `String` per line via `read_line` (`src/gfa/mod.rs:75-76`) and heap-allocates owned `String` and `Vec<u8>` copies for every field in `GfaSegment`, `GfaLink`, `GfaPath`, and `GfaWalk` (`src/gfa/mod.rs:14-46, 86-117`). Second, `vg` and `odgi` exhaust RAM on 12GB pangenomes not because they lack a line-by-line parser, but because graph alignment, bubble detection, and path traversal require random access to the graph's segments, links, and walks in memory (or an indexed store). Streaming GFA lines one at a time (and `stream_gfa` is never called anywhere in `src/`) does not eliminate the memory requirement of the downstream graph algorithms.
- **Evidence:** `src/gfa/mod.rs:14-46, 75-117`; `src/lib.rs:394` (no callers of `stream_gfa` in `src/`).
- **Correction:** Remove the "zero-copy" claim for `src/gfa/mod.rs` and state that streaming GFA parsing only bounds I/O buffer memory during sequential scans, whereas graph indexing, alignment, and topology resolution still require either in-memory graph structures or disk-backed indices.

## F5 -- C5: Purely streaming standard GFA over `stdin`/`stdout` contradicts random-access and sidecar requirements

- **Claim:** C5 (`docs/STRATEGIC_PLAN.md` / `ATTACKED-STRATEGIC_PLAN.md:17`)
- **Verdict:** REFUTED
- **Concrete flaw:** The claim that tools communicate *purely* by streaming standard GFA/MAF/VCF over `stdin`/`stdout` is internally contradicted in the same paragraph by the requirement for optional sidecar index files (`ATTACKED-STRATEGIC_PLAN.md:17`), and is operationally impossible for graph-global tools (`horizon_wga --graph`, `horizon_resolve`, `horizon_lift`) that require random access to `S`, `L`, `W`, and `P` records appearing in arbitrary order. Furthermore, `src/gfa/mod.rs:50-73` opens a `File` from a filesystem `Path`, not `stdin`, and `src/formats.rs:50, 102` emits custom `# EG:Z:...` comment lines and `*`-sequence bubble segments (`HB:i:1`) that downstream sequence-indexing tools cannot align against.
- **Evidence:** `ATTACKED-STRATEGIC_PLAN.md:17`; `src/gfa/mod.rs:50-73`; `src/formats.rs:50, 102`.
- **Correction:** State that GFA v1.0/v1.1 and MAF/VCF serve as the interchange formats on disk, while multi-pass or indexed file access (rather than pure `stdin`/`stdout` streaming) is required for graph-wide operations.

## F6 -- C6: `stdin`/`stdout` GFA/MAF requirement is violated by every existing binary and impossible for dual-input alignment

- **Claim:** C6 (`docs/STRATEGIC_PLAN.md` / `ATTACKED-STRATEGIC_PLAN.md:23`)
- **Verdict:** REFUTED
- **Concrete flaw:** None of the existing binaries obey this rule: `horizon_wga` (`src/bin/horizon_wga.rs:121-127, 171-172`) reads two FASTA file paths and writes MAF/GFA to `--maf`/`--gfa` file paths; `horizon_merge` (`src/bin/horizon_merge.rs:12, 20-35`) reads `<sp=path.maf>` files and writes `<maf_out>` and `<gfa_out>` files; `horizon_guide` (`src/bin/horizon_guide.rs`) reads FASTA files and writes output files. Logically, `horizon_wga` (pairwise or `--graph`) requires two distinct inputs (target graph/FASTA and query FASTA), which cannot both be streamed from a single `stdin` without out-of-band file arguments or buffering.
- **Evidence:** `src/bin/horizon_wga.rs:121-127, 171-172, 206-217`; `src/bin/horizon_merge.rs:12, 20-35`.
- **Correction:** Drop the universal `stdin`/`stdout` mandate and specify which tools are single-stream filters versus multi-input file/index tools.

## F7 -- C7: Linear streaming cannot detect hairballs or stitch resolved paths back into an already-emitted stream

- **Claim:** C7 (`docs/STRATEGIC_PLAN.md` / `ATTACKED-STRATEGIC_PLAN.md:27`)
- **Verdict:** REFUTED
- **Concrete flaw:** Detecting a "hairball" (high-degree tangle, cycle, or nested bubble) in a GFA stream is impossible in a single linear pass because GFA v1.0/v1.1 does not order `S`, `L`, `P`, or `W` lines by graph locality. An `L` line near the end of a 12GB GFA can link back to an `S` line at the beginning, forming a cycle or tangle after upstream `S`/`L` lines have already been streamed to `stdout`. Once emitted to a Unix pipe, upstream records cannot be retracted or stitched with a resolved e-graph path. Conversely, buffering until a connected component closes requires holding an entire chromosome-scale component in memory, negating linear streaming.
- **Evidence:** `ATTACKED-STRATEGIC_PLAN.md:27`; `src/gfa/mod.rs:76-119`; `src/egraph.rs:134-161` (takes a pre-built linear `&[(&str, i8)]`, not a GFA subgraph).
- **Correction:** Replace linear-stream bubble extraction with an explicit graph-partitioning pass (e.g., an on-disk or indexed connected-component / ultrabubble decomposition) prior to localized e-graph resolution.

## F8 -- C8: Release strategy targets three non-existent tools while the codebase implements pairwise/progressive WGA

- **Claim:** C8 (`docs/STRATEGIC_PLAN.md` / `ATTACKED-STRATEGIC_PLAN.md:31`)
- **Verdict:** REFUTED
- **Concrete flaw:** The strategic conclusion to focus benchmarking and release *entirely* on `horizon_clean`, `horizon_resolve`, and `horizon_lift` "first to capture the pangenome community immediately" does not follow from the repository's actual state. None of `horizon_clean`, `horizon_resolve`, or `horizon_lift` exist in `src/bin/`; `src/coord.rs` is not even included in `src/lib.rs`; and `src/gfa/mod.rs` is an unused 124-line parser with no graph index, adjacency table, or writer. All working binaries (`horizon_wga`, `horizon_guide`, `horizon_merge`, `horizon_eval`, `bench_chr1`) implement pairwise and progressive multi-species WGA.
- **Evidence:** `src/bin/` listing; `src/lib.rs:3-9, 394`; `src/gfa/mod.rs:1-124`.
- **Correction:** Acknowledge that `horizon_clean`, `horizon_resolve`, and `horizon_lift` are unimplemented roadmap items, and ground immediate benchmarking/release plans in the existing WGA pipeline (`horizon_wga`, `horizon_guide`, `horizon_merge`).

## F9 -- C9: `horizon_clean` cannot score GFA `L` lines via BOCPD/MDL in a stream

- **Claim:** C9 (`docs/PANGENOME_SUITE.md` / `ATTACKED-PANGENOME_SUITE.md:8`)
- **Verdict:** REFUTED
- **Concrete flaw:** Evaluating every GFA `L` line in a stream with the Horizon BOCPD/MDL scoring engine fails on two levels. First, a GFA `L` line (`src/gfa/mod.rs:21-27`) contains only segment IDs, orientations, and an overlap CIGAR (typically `0M` in pangenome graphs)—it contains no DNA sequences. Second, Horizon's BOCPD/MDL engine (`bocpd_gc_profile` in `src/bocpd.rs:58` and `banded_mdl_gap` in `src/mdl/dp.rs:81`) scores *pairwise homologous sequence alignments* `(seq_x, seq_y)`, not topological adjacency links between consecutive segments. Even if one attempted to score alternative paths across a bubble, doing so in a stream is impossible because `L` lines can appear before their referenced `S` lines or gigabytes away from them, so `stream_gfa` does not have the segment sequences in scope when an `L` line arrives.
- **Evidence:** `src/gfa/mod.rs:21-27, 91-99`; `src/bocpd.rs:58`; `src/mdl/dp.rs:81`.
- **Correction:** Specify that `horizon_clean` must first index `S`-line sequences and `L`-line topology to identify parallel/alternative bubble paths, and only then apply pairwise MDL scoring between competing path sequences.

## F10 -- C10: Banded linear MDL DP (`BandGrid`) cannot run across graph bifurcations

- **Claim:** C10 (`docs/PANGENOME_SUITE.md` / `ATTACKED-PANGENOME_SUITE.md:16`)
- **Verdict:** REFUTED
- **Concrete flaw:** Horizon's banded MDL DP (`banded_mdl_gap` using `BandGrid` in `src/mdl/dp.rs:19-72, 81-120`) is defined strictly over a 2D linear coordinate grid `(i, j)` with a single linear diagonal `center(i) = (i * m) / n`. At a graph bifurcation (or cycle), there is no single linear coordinate `i` or diagonal `center(i)`, and `bocpd_gc_profile` (`src/bocpd.rs:58`) depends on linear sequence history (`max_run = 80`). Extending forward-backward MDL DP across graph bifurcations requires a DAG/graph-HMM DP formulation (and cycle unrolling), which does not exist in `src/mdl/dp.rs`, nor does `horizon_wga --graph` exist in `src/bin/horizon_wga.rs`.
- **Evidence:** `src/mdl/dp.rs:19-72, 81-120`; `src/bocpd.rs:58`; `src/bin/horizon_wga.rs:32-120`.
- **Correction:** Replace "performing banded MDL dynamic programming across graph bifurcations" with an explicit graph-DP / Partial-Order-Alignment formulation (or linear path extraction followed by `banded_mdl_gap`) and mark it as unimplemented future work.

## F11 -- C11: `W`/`P` line traversal with `coord::Interval` cannot work without segment lengths and mismatches `W`/`P` semantics

- **Claim:** C11 (`docs/PANGENOME_SUITE.md` / `ATTACKED-PANGENOME_SUITE.md:29`)
- **Verdict:** REFUTED
- **Concrete flaw:** First, `src/coord.rs` is not even declared as a module in `src/lib.rs`, and `coord::Interval` (`src/coord.rs:4-38`) only implements `span()` and `ucsc_start()`—it has no GFA `W`/`P` traversal or projection logic. Second, neither `GfaWalk` (`W` lines) nor `GfaPath` (`P` lines) in `src/gfa/mod.rs:29-46` stores segment lengths; traversing a walk/path to project coordinates requires random access to a global `S`-line length index (and overlap CIGARs). Third, `Interval::ucsc_start` (`src/coord.rs:31-37`) assumes UCSC/MAF reverse-strand coordinates relative to `src_size`, whereas GFA `W` lines use forward-strand haplotype coordinates regardless of segment orientation (`>`/`<`) and GFA `P` lines carry no `src_size` or haplotype bounds at all.
- **Evidence:** `src/coord.rs:4-38`; `src/lib.rs:3-9`; `src/gfa/mod.rs:29-46`.
- **Correction:** Wire `coord` into `src/lib.rs`, build a global segment-length index from `S` lines before traversing `W`/`P` lines, and distinguish GFA forward-strand haplotype intervals from UCSC reverse-strand intervals.

## Verdict Summary

| C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 |
|---|---|---|---|---|---|---|---|---|---|---|
| REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED | REFUTED |
