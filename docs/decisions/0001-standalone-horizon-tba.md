# ADR-001: Standalone `horizon_tba` Reference-Free All-vs-All Aligner

## Status
Accepted (updated 2026-11-25 after Wave C scheduling, two-tier `--hinges`/`--auto-hinges`/`--extend-maf` [ADR-005], `--progressive-ancestors` [ADR-004], and `b131638` bounded-RAM multi-contig scheduling)

## Date
2026-11-22

## Context
Multi-species genome comparisons across divergent clades (e.g., 20 *Plasmodium* species, 17 *Candidozyma* species, 44–119 coronaviruses, or 7 *Saccharomyces* yeasts) lose clade-specific homologous blocks when anchored on a single hinge reference (`horizon_guide`). Users need a reference-free Threaded Blockset Aligner (TBA) workflow that:
- Computes all $\binom{N}{2}$ pairwise alignments, two-tier hinge + tip projections (`--hinges` / `--tip-map`, $O(H^2 + (N - H))$), post-hoc incremental tip extensions (`--extend-maf <backbone.maf>`), or merge-only runs (`--from-mafs <dir>`),
- Preserves clade-specific blocks absent from outgroups,
- Validates an explicit Newick guide tree (`--tree <newick|file>`) up front before compute or builds an automatic UPGMA tree using `[0, 1]`-normalized aligned-bp / MinHash Jaccard similarity,
- Supports atomic checkpointing/resuming of pairwise MAFs (`--pairwise-dir <dir>`) with bounded peak RAM on multi-contig genomes.

## Decision
Implement TBA as a standalone binary (`src/bin/horizon_tba/{main,cli,pipeline}.rs`, library core in `src/tba/{mod,model,merge,poa_refine,tree,two_tier,writer}.rs` and `src/tba/progressive/{mod,ancestral,hal_build}.rs`), dispatched through the unified CLI as `horizon tba`:
1. **Direct `RunAlignment -> TbaBlock` conversion, immediate `RunAlignment` drop, and Phase-1/Phase-2 `--pairwise-dir` disk spill (`b131638`):** Format each contig pair's MAF chunk inside the worker and drop `RunAlignment` immediately instead of holding all DP alignments across a species pair. When `--pairwise-dir` is set, Phase 1 streams finished pairwise MAFs atomically (`create_new(true)` + `0o600` + `rename`) without keeping `Vec<PairwiseTbaBlock>` in RAM, frees all `ContigIndex` seed tables completely, and only then loads `PairBlockMap` for Phase 2 tree merging—guaranteeing `ContigIndex` and `PairBlockMap` are never resident simultaneously.
2. **Non-nested adaptive Rayon scheduling (`b131638`):** Avoid nested `sp_pairs.par_iter()` × `contig_pairs.par_iter()` working-set amplification (which started 32 species pairs simultaneously and hit `29.1 GiB` OOM on 17-way *Candidozyma* before a single pair finished). Multi-contig panels (`avg >= 3` contigs/species) iterate species pairs sequentially while saturating all cores inside `contig_pairs.par_iter()` (`1×` active DP working set, `5.5–7.0 GiB` pairwise RSS at `26.0×` CPU efficiency on 32 cores); single-chromosome panels (e.g., single-chromosome sex-chromosome panels, 44-way coronaviruses) parallelize across `sp_pairs.par_iter()` with a serial inner loop.
3. **Borrowed tree merges, target-only tip indexing, & gapless GFA export:** `run_tba_tree_merge` threads cross-clade links using $O(\log M)$ `BTreeMap` `IntervalSet` tracking (`src/netting.rs`), two-tier tip projection uses `CladeIndex::build_for_species` to index only the single assigned hinge species per tip, and `write_tba_maf_and_gfa` reuses per-block row buffers (`row_offsets`, `mask`) while decomposing gapped blocks into gapless `S` sub-segments so every `W` walk satisfies $\sum \mathrm{LN} = W_{\mathrm{end}} - W_{\mathrm{start}}$.

## Alternatives Considered

### Mode flag (`horizon_guide --tba`)
- **Pros:** One fewer binary in `src/bin/`.
- **Cons:** `horizon_guide` is architected around a central hinge species, BFS transitive sister rescue onto hinge coordinates, and progressive HAL export; grafting an all-vs-all reference-free blockset merger and `--from-mafs` CLI mode onto `horizon_guide` would overload its CLI flags and data flow.
- **Rejected:** Violates the suite's single-responsibility binary pattern (`horizon_wga`, `horizon_guide`, `horizon_merge`, `horizon_gfa`, etc.).

### Chunked species-pair barriers (`sp_pairs.chunks(4)`)
- **Pros:** Initially tested to cap concurrent `DpScratch` allocations on small 8-species runs.
- **Cons:** Introduced global synchronization barriers between every 4 species pairs, stranding worker threads on long pairs. Removing the barrier (Wave C-opt-1) cut 44-way `wuhCor1` TBA (`-t 16`) wall time from `28.03 s` to `19.99 s` (`1.40×` faster), lifted CPU utilization from `8.98×` to `14.30×`, cut voluntary context switches `27.5×` (`43,097 → 1,567`), and *reduced* peak RSS (`646.3 → 629.2 MiB`) thanks to streaming pairwise emission and owned `ContigIndex` zero-copy reads.
- **Rejected:** Replaced by flat `par_iter()` + streaming pairwise MAF emission.

## Consequences
- Clean separation between reference-guided (`horizon_guide`, `horizon_merge`) and reference-free (`horizon_tba`) workflows.
- 44-way Sarbecovirus TBA (`-t 16`) runs in `18.78 s` (`14.35×` CPU, `621.8 MB` RSS) at `0.9805` exact F1 (`813` GFA segments, `96/96` exact walks); yeast `chrIII` TBA achieves `15.44×` CPU at `-t 16` (`15,368/15,368` exact walks).
- Users can reuse precomputed pairwise MAFs from external aligners or resumable HPC array jobs via `--from-mafs`.
