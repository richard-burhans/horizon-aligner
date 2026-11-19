# Horizon Pangenome Suite: Strategic Plan and Risk Mitigation

## 1. The Vision: Composition over Monolith
The history of bioinformatics alignment is a history of composition: `lastz` (pairwise) feeding into `multiz` (MSA), or `cactus` (graph) feeding into `hal` (storage). The vision for `horizon-aligner` is to unify these fragmented tools under a single, theoretically grounded framework based on Minimum Description Length (MDL) and E-graphs. Horizon represents the unified composition of these concepts.

We are not building a monolithic "run-all" executable. We are building a suite of specialized tools (`horizon_clean`, `horizon_wga`, `horizon_resolve`, `horizon_lift`) that compose perfectly via Unix pipes.

## 2. What We Are Doing Right
* **MDL over Heuristics:** We are replacing arbitrary affine gap penalties and scoring matrices with Information-Theoretic bit-savings (MDL). This inherently prevents the "hairball" over-alignment problem.
* **E-graphs for Structural Variants:** Using congruence-closure to resolve complex, nested structural paralogs (like the HLA region). This natively models evolutionary trajectories without producing unnavigable cyclical graphs.
* **Native Rust Streaming (Bypass `vg`):** `vg` and `odgi` crash on 12GB pangenomes due to RAM exhaustion. Our custom zero-copy GFA streaming parser (`src/gfa/mod.rs`) solves the I/O bottleneck entirely.

## 3. Strategic Adjustments & Risk Mitigations

### A. Avoid the "Walled Garden" Format (Learn from `HAL`)
**Risk:** `Cactus` introduced the `.hal` format, which required heavy C-dependencies (`libhdf5`), fractured the ecosystem, and made downstream parsing difficult in other languages.
**Mitigation:** We will **never invent a proprietary binary format**. Our tools will communicate purely by streaming standard GFA (v1.0/v1.1) and standard MAF/VCF over `stdin`/`stdout`. Any necessary indexing must be a lightweight, optional sidecar file (like `.fai` or `.gzi`), never a locked central database.

### B. Embrace the Unix Philosophy
**Risk:** The temptation to build a monolithic CLI that does everything internally.
**Mitigation:** Ensure absolute composability via pipes. The ecosystem should look like:
`cat hprc.gfa | horizon_clean --mdl-threshold 5 | horizon_wga --query novel_genome.fa | horizon_resolve > final.gfa`
Every tool must read GFA/MAF from `stdin` and write to `stdout`.

### C. Localize the E-Graphs (Prevent Combinatorial Explosion)
**Risk:** E-graphs suffer from combinatorial explosion. Loading an entire chromosome into an e-graph to run saturation rules will never finish computing.
**Mitigation:** Strict **bubble extraction**. We will stream the pangenome linearly and only invoke the E-graph engine when we hit a complex, tangled sub-graph ("hairball"). The E-graph resolves the optimal path for that bubble, and we stitch the resolved path back into the linear stream.

### D. Target the Graph Market First
**Risk:** Competing against `lastz` or `minimap2` for simple pairwise alignment is an uphill battle because they are "good enough" for standard tasks.
**Mitigation:** The massive, urgent pain point right now is **Pangenome Graphs**. The HPRC produces massive graphs that researchers cannot use because standard tooling crashes or produces tangled paths. We will focus our benchmarking and release strategy entirely on `horizon_clean`, `horizon_resolve`, and `horizon_lift` first to capture the pangenome community immediately. Pairwise/MSA will follow naturally.
