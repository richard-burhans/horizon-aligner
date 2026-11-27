# ADR-007: Pangenome Interoperability — Incremental rGFA (`horizon_gfa rgfa`), Snarl VCF Export (`horizon_resolve --vcf`), and Lossless `.hnex` Archives

## Status
Accepted (implemented in `src/rgfa.rs`, `src/bin/horizon_resolve/vcf.rs`, `src/bin/horizon_clean/engine.rs`, `src/formats/hnex.rs`, and `src/tba/hnex.rs`, hardened in Pass A)

## Date
2026-11-25

## Context
Integrating `horizon-aligner` into modern pangenome workflows (Minigraph-Cactus, `vg`, `bcftools`, and comparative phylogenomics) requires three interoperability formats beyond standard MAF and GFA v1.1:
1. **Minigraph-style incremental rGFA (`SN`/`SO`/`SR` tags):** Constructing a stable reference-coordinate pangenome backbone by grafting structural variants ($\ge \text{min\_sv\_bp}$, $\le \text{DEFAULT\_MAX\_SV\_BP} = 1\text{ Mbp}$) from query assemblies onto rank-0 reference contigs one assembly rank at a time.
2. **Haplotype-phased VCFv4.2 variant projection from GFA snarls:** Calling SNPs, indels, and multiallelic structural variants directly from Tarjan bridge-delimited snarls in a GFA v1.0/v1.1 graph without external `vg deconstruct` dependencies, plus mapping-support and bubble-clipping filters in `horizon_clean`.
3. **Lossless probabilistic alignment exchange (`.hnex`):** Standard MAF discards per-column posterior probabilities and unaligned Horizon Bubble spans. Downstream phylogenomics and uncertainty-aware tools need a human-readable, bounded-memory, self-validating archive that round-trips losslessly with `TbaBlock`.

## Decision
1. **Incremental rGFA graph builder (`src/rgfa.rs`, `horizon_gfa rgfa`):**
   - `RgfaGraph::from_reference_contigs` seeds rank `0` segments with canonical `SN:Z:`, `SO:i:`, and `SR:i:0` tags.
   - `graft_assembly_contigs` assigns a single rank per grafted assembly, caches rank-0 `ContigIndex`es across queries, applies 1-to-1 netting (`net_one_to_one_runs_scored`), symmetrically trims overlapping flank coordinates, stores reverse-strand SV insertions in forward rGFA coordinates (`SN`/`SO`) linked with `-` orientations, deduplicates links via `HashSet<RgfaLink>`, and uses independent split guards in `split_rank0_at`.
   - `read_rgfa` enforces bounded line lengths, strict `+`/`-` orientations, `checked_add` on `SO:i:`, and max-ID tracking for `next_seg_num`.
2. **Tarjan-snarl VCFv4.2 export & `horizon_clean` pangenome filters (`src/bin/horizon_resolve/vcf.rs`, `src/bin/horizon_clean/engine.rs`):**
   - `horizon_resolve --vcf <out.vcf[.gz]> [--vcf-ref <sample>]` projects bridge-bounded snarls onto VCFv4.2 with phased diploid/haploid genotypes, skipping `*`-sequence anchor/reference snarls, emitting `.` for `*`-sequence alts, trimming closing link overlaps on non-`0M` graphs, and emitting `.` via `TraversalEntry::Conflict` when a haplotype traverses the same snarl multiple times with conflicting alleles.
   - `horizon_clean` adds `--min-haplotypes <K>` (parsing PanSN `(sample, hap)` pairs on both `W` and `P` lines) and `--clip-bubbles-bp <N>` for Minigraph-Cactus-style rare-allele/large-bubble clipping while preserving `S`-line optional tags.
3. **Lossless `.hnex` probabilistic archive (`src/formats/hnex.rs`, `horizon_locus --hnex`):**
   - Line-oriented text format (`#HNEX v1`) carrying `BLOCK`, `ROW`, `POST` (per-column posterior probabilities), and `BUBBLE` (`SPAN`) records with lossless `from_tba_blocks`/`to_tba_blocks` conversion (including bubble-only contigs).
   - Enforces a 256 MiB bounded line reader and strict fail-closed validation of strands, coordinate bounds, ASCII sequences, ungapped lengths, uniform column widths, finite `[0, 1]` posteriors, and control-character hygiene on write.

## Alternatives Considered

### Shell out to `minigraph` / `vg deconstruct`
- **Cons:** Introduces heavy external C++ toolchain dependencies and loses Horizon's MDL/Horizon Bubble semantics.
- **Rejected:** Pure-Rust implementations keep the suite self-contained and composable.

### Store posteriors in custom MAF comment lines
- **Cons:** Unstructured MAF comments break strict downstream MAF parsers and cannot cleanly represent bubble-only contigs lacking aligned rows.
- **Rejected:** Dedicated self-validating `.hnex` format (`src/formats/hnex.rs`) preserves full probabilistic fidelity alongside standard MAF/GFA output.

## Consequences
- `horizon_gfa rgfa`, `horizon_clean --min-haplotypes/--clip-bubbles-bp`, `horizon_resolve --vcf`, and `horizon_locus --hnex` form a complete, pure-Rust pangenome construction, cleaning, variant-calling, and probabilistic-archive pipeline.
- All parsers and writers enforce strict input validation and reject stdout collisions (e.g., `horizon_resolve --vcf -` without `--out`).
