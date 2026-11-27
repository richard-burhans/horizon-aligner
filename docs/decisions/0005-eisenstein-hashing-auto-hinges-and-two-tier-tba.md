# ADR-005: Eisenstein Sixth-Root Hashing, Auto-Hinges, and Two-Tier TBA

## Status
Accepted (implemented in `src/seeding/eisenstein/{mod,sketch,reconcile}.rs` and `src/tba/two_tier.rs`, streamed in Pass C)

## Date
2026-11-25

## Context
Scaling reference-free `horizon_tba` to large multi-genome panels (e.g., 20-species *Plasmodium* or 17-species *Candidozyma*) exposes two practical bottlenecks:
1. Full all-vs-all alignment requires $\binom{N}{2}$ pairwise jobs even when genomes cluster tightly into $H \ll N$ clades.
2. Users often want automatic, strand-invariant clade-hinge selection (`--auto-hinges H`) without manually curating a `--hinges` list or `--tip-map`, plus post-hoc incremental extension (`--extend-maf`) and low-disk gzipped pairwise checkpointing (`--pairwise-gz`).

At the same time, `kmer-polynomials` (see `docs/HORIZON_TIMELINE.md` Phase 1) proved two durable algebraic primitives over the Mersenne prime field $\mathbb{F}_q$ ($q = 2^{61} - 1$) using the sixth-root Eisenstein encoding ($A=1, C=\omega, G=-\omega, T=-1$):
- Strand-symmetric spaced-seed polynomial hashing (`KP-14`), and
- Exact algebraic power-sum symmetric-difference reconciliation (`KP-9`/`KP-10`, ~11.27 B/diff).

## Decision
1. **Eisenstein sixth-root hashing & streaming sketching (`src/seeding/eisenstein/`):**
   - `mod.rs`: Implements Mersenne-61 arithmetic (`fq_add`, `fq_sub`, `fq_mul`, `fq_inv`), canonical strand-symmetric `EisensteinHasher`, and zero-allocation streaming `for_each_kmer_hash` (Pass C) so genome-scale sketching never materializes a 16 B/bp `Vec<(usize, u64)>`.
   - `sketch.rs`: Implements `EisensteinSamplerSketch` (`from_contigs`, Jaccard similarity, and deterministic greedy k-center/DPP hinge selection `select_eisenstein_hinges`) used directly by `horizon_tba --auto-hinges <H>`.
   - `reconcile.rs`: Implements `PowerSumSketch` (Newton identities + Berlekamp–Massey / Chien search in $\mathbb{F}_q$) as a library primitive with fail-closed decoding when the error-locator degree is zero (Pass A2).
2. **Two-tier TBA & post-hoc MAF extension (`src/tba/two_tier.rs`, `src/bin/horizon_tba/`):**
   - `--hinges <sp1,sp2,...>` or `--auto-hinges <H>` restricts all-vs-all alignment to the $H$ hinge species ($\binom{H}{2}$ pairs) and projects the remaining $N - H$ tip species onto their assigned hinge via `CladeIndex::build_for_species` ($O(H^2 + (N - H))$ total work).
   - `--tip-map <tsv>` accepts explicit `<tip>\t<hinge>` mappings via bounded 64 KiB line reading (`read_bounded_line`) and rejects chained tip-to-tip projections (`tba_project_tips_onto_blockset`).
   - `--extend-maf <backbone.maf>` projects new tip genomes onto an existing multi-species TBA MAF without recomputing the backbone.
   - `--pairwise-gz` stores pairwise checkpoints in `--pairwise-dir` as `.maf.gz` and streams them lazily during Phase 2 merging (reducing 20-way *Plasmodium* peak merge RSS from 11.45 GiB to 3.79 GiB; see `reports/tba_panels_benchmark_report.md`).

## Alternatives Considered

### Use dense `FlatSeedIndex` for `--auto-hinges` sketching
- **Cons:** Requires building full 12 B/entry seed tables for all $N$ genomes up front before knowing which genomes are hinges vs tips, defeating the peak-RAM savings of two-tier scheduling.
- **Rejected:** Streaming `EisensteinSamplerSketch::from_contigs` sketches each genome in bounded memory and frees raw sequences immediately.

### Allow transitive tip -> tip -> hinge chains in `--tip-map`
- **Cons:** Introduces order-dependent coordinate drift and silently drops cycles or unprojected tips.
- **Rejected:** `tba_project_tips_onto_blockset` explicitly rejects chained tips and requires every tip to map directly to a backbone species.

## Consequences
- Two-tier `--auto-hinges` + `--pairwise-gz` cuts both pairwise compute ($O(H^2 + (N-H))$) and merge peak RSS (`3.0×` lower on 20-way *Plasmodium*), documented in `reports/tba_panels_benchmark_report.md`.
- `PowerSumSketch` remains a clean library-only algebraic reconciliation primitive under `horizon_aligner::seeding::eisenstein`.
