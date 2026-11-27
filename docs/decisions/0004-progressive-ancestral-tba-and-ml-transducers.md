# ADR-004: Bottom-Up Progressive Ancestral TBA, ML Indel Transducers, and POA Refinement

## Status
Accepted (implemented in `src/tba/progressive/{mod,ancestral,hal_build}.rs` and `src/tba/poa_refine.rs`, hardened in Passes A–C)

## Date
2026-11-25

## Context
All-vs-all Threaded Blockset Alignment (`horizon_tba`, ADR-001) scales as $O(N^2)$ pairwise alignments across $N$ genomes and retains only extant-species columns, whereas reference-guided progressive HAL export (`horizon_guide --progressive --hal`, ADR-002) anchors every ancestor on a single hinge coordinate frame. For deep multi-species panels (e.g., 20-species *Plasmodium* or vertebrate clades) users need:
1. $O(N)$ bottom-up Progressive Cactus-style scaling where each internal node aligns only its two children plus at most `K` nearest outgroup leaves (`--max-outgroups`, default `3`), skipping outgroup × outgroup pairs.
2. True reconstructed ancestral contigs (`Anc0..Anc_{n-2}`) emitted as FASTA (`--emit-ancestors <dir>`) and native hierarchical HAL (`--hal`, `build_progressive_hal`), with a choice between fast Fitch parsimony (`--ancestral-model fitch`) and probabilistic maximum-likelihood phylogenetic reconstruction (`--ancestral-model ml`).
3. Post-merge Partial Order Alignment (POA) block refinement (`src/tba/poa_refine.rs`) to close fragmented gap runs across threaded multi-species blocks.

## Decision
1. **Bottom-up progressive ancestral engine (`src/tba/progressive/`):**
   - `mod.rs`: Recursively traverses the guide tree (`run_progressive_ancestral_tba`), borrowing `leaf_genomes` directly via `Cow` slices to avoid redundant full-genome clones, selecting at most `max_outgroups` phylogenetically closest outgroup leaves (`select_node_outgroups`), indexing only the `2 + K` genomes active at the current internal node, and threading child clades onto the reconstructed ancestor.
   - `ancestral.rs`: Supports both Fitch parsimony (`reconstruct_ancestral_contigs`) and Ortheus/Holmes-style ML phylogenetic reconstruction (`reconstruct_ancestral_contigs_transducer`). The ML mode estimates empirical JC69 branch lengths from each subproblem (`estimate_subproblem_branches`), computes Felsenstein pruning likelihoods for substitutions, weights the outgroup clade by $1/K$ so $K$ correlated outgroups do not overwhelm the two children, counts all absent outgroups explicitly, and scores indels with an exact 2-state reversible continuous-time Markov chain ($\pi_1 = \lambda/(\lambda+\mu)$, $p_{\text{surv}} = \pi_1 + (1-\pi_1)e^{-rt}$, $p_{\text{ins}} = \pi_1(1-e^{-rt})$) that preserves exact detailed balance at root nodes.
   - `hal_build.rs`: Builds native hierarchical HAL graphs (`build_progressive_hal`) using `IntervalSet` availability (`src/netting.rs`), $O(\log R)$ `partition_point` lookup (`find_run`), and ascending child-coordinate sorting so split inversions and translocations tile cleanly without dropping half of a rearranged block.
2. **Rayon-parallel POA refinement (`src/tba/poa_refine.rs`):**
   - Refines qualifying multi-species TBA blocks in parallel (`blocks.par_iter_mut()`) by building a consensus POA DAG directly from the alignment matrix without spurious linear predecessor edges (`build_poa_dag_from_matrix`) and re-threading gapped rows.

## Alternatives Considered

### Keep only hinge-anchored progressive HAL (`horizon_guide`)
- **Cons:** Clade-specific sequences absent from the hinge species cannot be represented as native ancestral contigs in the HAL hierarchy.
- **Rejected:** Retained `horizon_guide --progressive --hal` (ADR-002) for reference-anchored $4\text{ MiB}$-windowed mammalian-scale streaming, while adding `horizon_tba --progressive-ancestors` for reference-free de novo ancestral reconstruction.

### Approximate indel survival/insertion with `1 - exp(-x)` vs `x / (1 + x)`
- **Cons:** Second-order $O(x^2/2)$ asymmetry caused root single-child insertions to be misclassified as ancestral even when $\lambda = \mu$.
- **Rejected:** Replaced in Pass A3 with the exact 2-state reversible CTMC transition probabilities.

## Consequences
- `horizon_tba --progressive-ancestors [--ancestral-model fitch|ml] [--max-outgroups K] [--emit-ancestors dir] [--hal]` delivers $O(N)$ scaling, de novo ancestral FASTAs, and multi-level hierarchical HAL files without holding more than `2 + K` `ContigIndex`es in memory.
- Cross-references ADR-001 (TBA core) and ADR-002 (hinge-anchored streaming HAL).
