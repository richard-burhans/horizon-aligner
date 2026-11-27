# ADR-008: Regional Locus ProbCons Consistency, BioNJ Trees, and Library-Only Primitives

## Status
Accepted (implemented in `src/mdl/{consistency,squint,train}.rs`, `src/tba/tree.rs`, `src/seeding/{cdbg,tantan}.rs`, and `src/bin/horizon_locus.rs`, hardened in Passes A–C)

## Date
2026-11-25

## Context
Regional multi-sequence locus alignment (`horizon_locus`, targeting 25 kb VGP ILS windows and BUSCO gene intervals) benefits from ProbCons-style probabilistic consistency transformation and neighbor-joining phylogenetic trees, while our survey of whole-genome alignment architectures (`pathogen-genomes/manuscript/wga_methods.md`) identified several foundational algorithmic primitives—colored de Bruijn graphs, HMM tandem-repeat masking (`tantan`), directional Pair-HMM Baum–Welch training, and coding frameshift/fixed-share null scoring (`squint`)—that belong cleanly in the reusable `horizon_aligner` library crate even before exposing dedicated CLI flags.

## Decision
1. **ProbCons probabilistic consistency & BioNJ in `horizon_locus` (`src/mdl/consistency.rs`, `src/tba/tree.rs`, `src/bin/horizon_locus.rs`):**
   - `horizon_locus` emits `<out_prefix>.bionj.nwk` via `build_bionj_tree` (Gascuel 1997 variance-weighted neighbor joining using unclamped branch lengths in $d_{uk}$ distance reduction, zeroed diagonal, and strict non-empty species validation).
   - `horizon_locus --consistency` validates unique sequence names, orients reverse-complemented contigs onto the focal strand prior to MSA (marking `strand=-` on output), and executes `align_locus_probcons_msa`:
     - Rolling-row 8 B/cell Forward-Backward posterior calculation (`SparsePosteriorMatrix::from_pair`) with automatic banding (`half_band = |n - m| + 256` above 1.5M cells) and a hard `MAX_PAIR_HMM_CELLS = 16_000_000` cap (~102 MB peak on 25 kb loci vs ~35 GB unbanded).
     - Auto-symmetrized `triplet_consistency_transform` using a reusable dense accumulator and sparse row-scatter `merge_profiles_mea` guided by a precomputed $N \times N$ similarity matrix.
2. **Documented library-only primitives (exposed in `horizon_aligner`, intentionally not wired as standalone CLI binaries):**
   - **`src/seeding/cdbg.rs` (`ColoredDeBruijnGraph`):** Palindrome-aware (`KmerOrient`), deterministic (`BTreeMap`) colored compacted de Bruijn graph builder with $(k-1)$-overlap unitig chaining.
   - **`src/seeding/tantan.rs` (`tantan_posterior_probs`, `tantan_mask`):** Frith (2011) period-$P$ HMM tandem-repeat masker using compact 1 B/bp codes, single-loop backward pass, finite/NaN config guards, and 64 kb overlapping-window chunking (~12 MiB bounded RAM on arbitrarily long chromosomes).
   - **`src/mdl/train.rs` (`train_pair_hmm_em`):** Directional Pair-HMM Baum–Welch EM trainer with proper `pi_x`/`pi_y` renormalization and no self-reinforcing mismatch attenuation.
   - **`src/mdl/squint.rs` (`fixed_share_null_bits`, `repair_coding_frameshifts`):** Fixed-share compositional null bit-rate model and codon-aware 3-frame Pair-HMM frameshift repair (neutral `0.0` for `NNN`/`X`, exact `-1` nt deletion alignment, `16M` cell cap).
   - **`src/seeding/eisenstein/reconcile.rs` (`PowerSumSketch`):** Exact algebraic symmetric-difference set reconciliation over $\mathbb{F}_{2^{61}-1}$ (see ADR-005).

## Alternatives Considered

### Expose every library primitive as a separate CLI binary
- **Cons:** Fragments the 12-binary suite with overlapping utilities before real-world CLI contracts are validated (consistent with Round 8 adversarial triage in `docs/COMMUNITY_UTILITIES_PROPOSAL.md`).
- **Rejected:** Keep `cdbg`, `tantan`, `train`, `squint`, and `PowerSumSketch` as tested, documented library APIs under `horizon_aligner` until specific workflow integration justifies CLI surface.

## Consequences
- `horizon_locus` provides accurate, bounded-memory ProbCons consistency MSA (`--consistency`), BioNJ trees (`.bionj.nwk`), and lossless `.hnex` export (`--hnex`) on 25 kb regional loci.
- Library consumers have clear, explicit documentation distinguishing CLI-wired pipelines from reusable library primitives.
