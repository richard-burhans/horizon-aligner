# Round 1 Attack -- L6 (prior-art)

## Literature, Method Cards, and Prior-Round Search
- **Held literature & docs located:** `horizon-aligner` declares no `refs/catalog_index.json` or method cards. Held literature and internal manuscripts were located at `manuscript/references.bib` (106 lines), `manuscript/horizon_paper.md` (112 lines), `manuscript/guide_pure_rust_hal.md` (98 lines), `ARCHITECTURE.md` (28 lines), `src/gfa/mod.rs` (124 lines), `src/egraph.rs` (385 lines), and `src/bin/horizon_wga.rs`.
- **Prior rounds checked:** `scripts/round history "vg"` returned no closed rounds (Round 1 is the first round).

---

## F1 -- C1: Unix-pipe composition of specialized pangenome tools is established prior art in `odgi` and `vg`, and the named tools do not exist
- **Claim:** `C1` (`docs/STRATEGIC_PLAN.md:6` / `ATTACKED-STRATEGIC_PLAN.md:6`)
- **Verdict:** **NARROWED**
- **Flaw:** Presenting a suite of specialized Unix-pipe-composable pangenome tools (`horizon_clean`, `horizon_wga`, `horizon_resolve`, `horizon_lift`) as Horizon's differentiating vision overlooks that both **`odgi`** (Guarracino et al. 2022, *Bioinformatics* 38(13):3319–3326) and **`vg`** (Garrison et al. 2018, *Nature Biotechnology* 36:875–879), as well as **`gfatools`** (Li), already implement modular subcommands explicitly designed to compose via Unix `stdin`/`stdout` pipes (e.g., `odgi build -g - -o - | odgi sort ... | odgi view ...`). Moreover, checking `src/bin/` shows that `horizon_clean`, `horizon_resolve`, and `horizon_lift` do not exist; the implemented binaries are `horizon_wga`, `horizon_merge`, `horizon_guide`, `horizon_eval`, and `bench_chr1`.
- **Correction:** Narrow `C1` to state that Horizon adopts the established Unix-pipe subcommand pattern used by `odgi`/`vg` rather than introducing it, and restrict named tools to those actually implemented (`horizon_wga`, `horizon_merge`, `horizon_guide`).

## F2 -- C2: MDL alignment is from Allison et al. (1992) and does not prevent high-identity paralogous "hairballs"
- **Claim:** `C2` (`docs/STRATEGIC_PLAN.md:9` / `ATTACKED-STRATEGIC_PLAN.md:9`)
- **Verdict:** **REFUTED**
- **Flaw:** First, replacing affine gap penalties with three-state Minimum Description Length (MDL) alignment is prior art from **Allison, Wallace & Yee (1992)** (*J. Mol. Evol.* 35:77–89), which Horizon's own `manuscript/references.bib:8-17` and `manuscript/horizon_paper.md:36` cite. Second, the causal claim that MDL bit-savings "inherently prevents the 'hairball' over-alignment problem" conflates two distinct phenomena: MDL prevents *low-signal, compositionally biased random* over-alignment (e.g., 18% GC *Plasmodium*), whereas pangenome graph "hairballs" in HLA, centromeres, and segmental duplications—addressed by **PGGB / `smoothxg`** (Garrison et al. 2024, *Nature Methods*) and **Minigraph-Cactus** (Hickey et al. 2024, *Nature Biotechnology* 42:663–673)—are caused by *genuine high-identity multi-copy homology* (repeats and paralogs) that yields high alignment scores and positive MDL bit-savings.
- **Correction:** Replace the second sentence with: "MDL bit-savings suppresses low-signal, compositionally biased random over-alignment, though high-identity paralogous repeats still require graph-level decomposition."

## F3 -- C3: Congruence-closure on `seq`/`inv`/`dup` (`egg`) does not natively model HLA evolutionary trajectories
- **Claim:** `C3` (`docs/STRATEGIC_PLAN.md:10` / `ATTACKED-STRATEGIC_PLAN.md:10`)
- **Verdict:** **REFUTED**
- **Flaw:** E-graphs (**`egg`**, Willsey et al. 2021, *PACMPL* 5:1–29, cited at `manuscript/references.bib:97-106`) are algebraic term-rewriting structures, not evolutionary models. Inspection of `src/egraph.rs:16-25` shows Horizon's `Expr` only models four syntactic operators (`Block`, `Seq`, `Inv`, `Dup`) with tandem collapse `seq(X, X) -> dup(X)`. Real evolutionary trajectories of complex nested paralogs like the HLA region involve dispersed segmental duplications, gene conversion, interlocus recombination, and trans-species polymorphism—none of which are modeled by tandem `dup(X)` collapse. Existing pangenome tools already handle cyclical/tangled loci via graph smoothing (`smoothxg` in PGGB, Garrison et al. 2024) or hierarchical cactus graphs (Minigraph-Cactus, Hickey et al. 2024).
- **Correction:** State that `src/egraph.rs` collapses local tandem duplications and nested inversions syntactically via equality saturation, not that it natively models HLA evolutionary trajectories.

## F4 -- C4: `vg`/`odgi` do not crash on 12GB pangenomes, and `src/gfa/mod.rs` is neither zero-copy nor a graph-I/O solution
- **Claim:** `C4` (`docs/STRATEGIC_PLAN.md:11` / `ATTACKED-STRATEGIC_PLAN.md:11`)
- **Verdict:** **REFUTED**
- **Flaw:** Both halves of `C4` are factually contradicted by prior art and the codebase:
  1. **`vg` / `odgi` memory:** `libbdsg` (**Eizenga et al. 2020**, *Bioinformatics* 36(22–23):5139–5144) introduced `PackedGraph` (up to 20× less RAM than legacy `vg`) and `ODGI` specifically for gigabase-scale pangenomes. `odgi` (**Guarracino et al. 2022**) processes HPRC graphs at ~8× lower memory than legacy `vg`. `GBZ` (**Sirén & Paten 2022**, *Bioinformatics* 38(22):5012–5019) compresses the 90-haplotype HPRC v1.1 graph to ~3.4 GB, and **`GBZ-base`** (Sirén et al. 2026) and **`gfaidx`** provide disk-backed indexed access without loading whole graphs into RAM. **`gfatools`** (Li) already provides C-level GFA streaming.
  2. **`src/gfa/mod.rs`:** Inspection of `src/gfa/mod.rs:87-116` shows it is **not zero-copy**: it allocates a `String` line buffer and calls `.to_string()` and `.as_bytes().to_vec()` on every `S`, `L`, `P`, and `W` record. Nor does line-by-line streaming solve graph I/O for alignment or liftover, which require random access across `S` and `W`/`P` records.
- **Correction:** Remove the false claim that `vg`/`odgi` crash on 12GB pangenomes and remove "zero-copy" from the description of `src/gfa/mod.rs`.

## F5 -- C5: `.hal` characterization is inaccurate and Horizon's own codebase emits `.hal` and custom GFA tags
- **Claim:** `C5` (`docs/STRATEGIC_PLAN.md:17` / `ATTACKED-STRATEGIC_PLAN.md:17`, premise at line 16)
- **Verdict:** **NARROWED**
- **Flaw:** `STRATEGIC_PLAN.md:16` characterizes `.hal` (Hickey et al. 2013, *Bioinformatics* 29(10):1341–1342) as a walled-garden format that fractured the ecosystem, ignoring that `.hal` is an open schema on standard HDF5 with `hal2vg`, `hal2maf`, and `halLiftover`, and that `Minigraph-Cactus` (Hickey et al. 2024) natively exports GFA/GBZ/VCF. More directly, Horizon's own `manuscript/horizon_paper.md:46,92,98` and `manuscript/guide_pure_rust_hal.md:8-17,86-93` implement a pure-Rust `.hal` reader and HDF5 `v0` `.hal` writer (`horizon_guide --progressive --hal`), and `manuscript/horizon_paper.md:40,92` emits non-standard `# HB:Z:` and `# EG:Z:` GFA extensions.
- **Correction:** Narrow `C5` to acknowledge that Horizon already reads/writes `.hal` in `src/hal/` and emits custom GFA comment tags (`# HB:Z:`, `# EG:Z:`), rather than communicating purely via unmodified standard GFA/MAF/VCF.

## F6 -- C6: Universal `stdin`/`stdout` streaming is violated by `src/gfa/mod.rs` and `src/hal/`
- **Claim:** `C6` (`docs/STRATEGIC_PLAN.md:23` / `ATTACKED-STRATEGIC_PLAN.md:23`)
- **Verdict:** **REFUTED**
- **Flaw:** While `stdin`/`stdout` streaming is standard in `odgi`/`vg`/`gfatools`, Horizon's own `src/gfa/mod.rs:50-64` (`stream_gfa<P: AsRef<Path>>(path: P, ...)`) takes a filesystem `Path` and cannot read from `stdin`. Likewise, `src/hal/` requires `Read + Seek` (`manuscript/guide_pure_rust_hal.md:15,20`), which is impossible on `stdin`, and `src/bin/horizon_wga.rs:122` requires file paths for target/query FASTA inputs.
- **Correction:** Drop the universal "Every tool must read GFA/MAF from `stdin`" requirement or refactor `stream_gfa` to accept `BufRead` and document file-path requirements for `horizon_wga` and `src/hal/`.

## F7 -- C7: Single-pass linear GFA streaming cannot extract arbitrary bubbles, and local bubble decomposition is prior art
- **Claim:** `C7` (`docs/STRATEGIC_PLAN.md:27` / `ATTACKED-STRATEGIC_PLAN.md:27`)
- **Verdict:** **NARROWED**
- **Flaw:** Standard GFA v1.0/v1.1 does not guarantee topological or linear ordering of `S` and `L` lines; links (`L`) can reference segments (`S`) anywhere in the file. Single-pass linear streaming cannot detect or extract a tangled bubble without first building an adjacency index or requiring pre-sorted graphs (as `odgi sort` + `vg` snarls do, Paten et al. 2018, *J. Comput. Biol.*). Furthermore, decomposing a graph into local bubbles/blocks, running local resolution (`abPOA`), and stitching back is established prior art in **`PGGB` / `smoothxg`** (Garrison et al. 2024).
- **Correction:** Specify that bubble extraction requires a pre-sorted graph or an indexed graph pass rather than arbitrary single-pass GFA streaming, and cite `vg` snarls / `smoothxg` as prior art for local bubble decomposition.

## F8 -- C8: HPRC graphs are widely usable with existing tooling, and `horizon_clean`/`resolve`/`lift` are unimplemented
- **Claim:** `C8` (`docs/STRATEGIC_PLAN.md:31` / `ATTACKED-STRATEGIC_PLAN.md:31`)
- **Verdict:** **REFUTED**
- **Flaw:** HPRC pangenome graphs (Liao et al. 2023, *Nature*; Hickey et al. 2024, *Nature Biotechnology*) are actively used in production with `vg giraffe` (Sirén et al. 2021, *Science* 374:abg8871), `GBZ` (Sirén & Paten 2022), `odgi` (Guarracino et al. 2022), `GraphAligner` (Rautiainen & Marschall 2020), and `Panacus`. The premise that researchers cannot use HPRC graphs because standard tooling crashes is false. Additionally, the strategic focus on `horizon_clean`, `horizon_resolve`, and `horizon_lift` contradicts `src/bin/` and `manuscript/horizon_paper.md`, where none of those three binaries exist and all benchmarks were run on `horizon_wga`, `horizon_merge`, and `horizon_guide`.
- **Correction:** Remove the claim that standard HPRC tooling crashes, and align the release strategy with the binaries that actually exist (`horizon_wga`, `horizon_merge`, `horizon_guide`).

## F9 -- C9: Streaming GFA `L` lines cannot be scored with BOCPD/MDL without `S`-line sequences
- **Claim:** `C9` (`docs/PANGENOME_SUITE.md:8` / `ATTACKED-PANGENOME_SUITE.md:8`)
- **Verdict:** **REFUTED**
- **Flaw:** In GFA v1.0/v1.1, `L` (Link) lines contain only segment names, orientations, and an overlap CIGAR (typically `0M` in bluntified HPRC graphs from Minigraph-Cactus/PGGB)—no nucleotide sequences. Scoring an edge with BOCPD/MDL requires the sequences of the linked `S` segments, which cannot be obtained by streaming `L` lines in arbitrary GFA order without indexing/caching `S` lines in memory or on disk. Graph edge pruning based on depth/alignment support is already implemented in `Minigraph-Cactus` / `vg clip` (Hickey et al. 2024) and `PGGB`/`smoothxg` (Garrison et al. 2024).
- **Correction:** State that `L`-line MDL evaluation requires indexing or caching `S`-line sequences first, not pure `L`-line streaming.

## F10 -- C10: Banded DP across graph bifurcations is prior art in `GraphAligner` and `minigraph`
- **Claim:** `C10` (`docs/PANGENOME_SUITE.md:16` / `ATTACKED-PANGENOME_SUITE.md:16`)
- **Verdict:** **NARROWED**
- **Flaw:** Threading query sequences through a graph using banded dynamic programming across bifurcations is established prior art in **`GraphAligner`** (Rautiainen & Marschall 2020, *Genome Biology* 21:253, bit-parallel banded DP across graph bifurcations), **`minigraph`** (Li et al. 2020, *Genome Biology* 21:265), and **`vg giraffe`** (Sirén et al. 2021). Only the substitution of MDL forward-backward scoring inside graph DP is potentially novel, and `horizon_wga --graph` is not yet implemented in `src/bin/horizon_wga.rs` (`--gfa` is output-only).
- **Correction:** Narrow `C10` to credit `GraphAligner`/`minigraph`/`vg giraffe` for banded sequence-to-graph DP across bifurcations and mark `horizon_wga --graph` as unimplemented.

## F11 -- C11: 0-based half-open `W`-line coordinates are standard GFA v1.1, and graph liftover is prior art
- **Claim:** `C11` (`docs/PANGENOME_SUITE.md:29` / `ATTACKED-PANGENOME_SUITE.md:29`)
- **Verdict:** **NARROWED**
- **Flaw:** GFA v1.1 `W` (Walk) lines already specify `SeqStart` and `SeqEnd` as **0-based half-open `[SeqStart, SeqEnd)` intervals**—this is standard GFA v1.1 semantics, not a Horizon innovation. Moreover, coordinate/annotation liftover across pangenome paths and alignments is established prior art in **`odgi position`** (Guarracino et al. 2022), **`halLiftover`** (Hickey et al. 2013), **`vg surject`**, and **`Liftoff`** (Shumate & Salzberg 2021, *Bioinformatics* 37(12):1639–1643). `horizon_lift` does not yet exist in `src/bin/`, and `src/gfa/mod.rs:100-117` only parses raw `P`/`W` fields without `coord::Interval` traversal.
- **Correction:** Narrow `C11` to note that 0-based half-open `W`-line arithmetic follows the GFA v1.1 specification, cite `odgi position` and `halLiftover` as prior art, and mark `horizon_lift` as unimplemented.

---

| C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| NARROWED | REFUTED | REFUTED | REFUTED | NARROWED | REFUTED | NARROWED | REFUTED | REFUTED | NARROWED | NARROWED |
