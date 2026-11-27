# ADR-006: `k`-Copy E-Graph Repeat Saturation and Overlap-Clustered MDL Duplication Pruning

## Status
Accepted (implemented in `src/egraph.rs` and `src/formats/dup_resolve.rs`, hardened in Passes A–C)

## Date
2026-11-25

## Context
Pairwise whole-genome alignment and pangenome walk simplification frequently encounter multi-copy tandem arrays ($k \ge 3$ copies, such as satellite repeats or multigene clusters) as well as spurious alignments that overlap by chance across modulo bin boundaries. Two concrete failure modes were identified:
1. **Binary-only `Dup` rewrites:** Modeling only 2-copy `Dup(x)` missed $k$-copy tandem arrays (`Seq(x, Seq(x, x))` -> `Repeat(x, 3)`) or nested inverted arrays (`Inv(Repeat(x, k))`).
2. **Modulo-bin false collisions & `Dup(Seq(..))` MDL bypass:** Bucketing colinear run-groups by `x_start / bin_bp` could split a single duplication across a bin edge or falsely collide adjacent distinct loci; furthermore, evaluating `dup`/`rep:k` framing overhead as a flat constant allowed `Dup(Seq(B, B))` to beat `Seq(Seq(B, B), Seq(B, B))` even at JC69 saturation ($d = 0.75$) where duplication carries zero bit-savings.

## Decision
1. **`k`-copy `Repeat` operator & leaf-scaled MDL extraction (`src/egraph.rs`):**
   - Extend the congruence-closed e-graph (`ENode::Repeat(Id, usize)`, `Expr::Repeat(Box<Expr>, usize)`) with associative/distributive saturation rules collapsing $k$-copy tandem arrays (`k >= 3`) and inverted repeats (`Inv(Repeat(x, k))`), while keeping deterministic `BTreeMap`/`BTreeSet` storage.
   - Track leaf count (`leaves: usize`) during `extract_min_mdl_cost` and scale `dup`/`rep:k` structural framing bits by leaf count so composite `Dup(Seq(B, B))` expressions cannot bypass the $d = 0.75$ JC69 saturation MDL gate.
2. **Overlap-clustered empirical MDL paralog pruning (`src/formats/dup_resolve.rs`):**
   - Replace fixed modulo-bin duplication detection in `horizon_wga --resolve-dups` with `resolve_alignment_duplications_mdl` and `prune_resolved_dup_clusters`:
     - Cluster colinear run-groups by $\ge 50\%$ shorter-span reference overlap using an active-cluster sweep (`cluster_run_groups_by_ref_overlap`, $O(G \cdot D)$ when sorted by ascending `x_start`).
     - Estimate empirical divergence per cluster using the two lowest-divergence copies (`cluster_div_min1`, `cluster_div_min2`) so a noisy third ectopic copy does not poison a true high-identity 2-copy tandem duplication.
     - Feed cluster tokens and empirical `(mean_len_bp, est_div)` directly into `MdlCostModel`, then prune redundant paralog run-groups while recomputing accurate `lcbs_fwd`/`lcbs_rc` counts.
3. **Standalone `horizon_resolve` integration (`src/bin/horizon_resolve/resolve.rs`):**
   - Apply the same `k`-copy `Repeat` saturation to bounded GFA snarl windows when `--collapse-dups` or `--tag-walks` is enabled.

## Alternatives Considered

### Retain fixed `bin_bp` modulo bucketing
- **Cons:** Sensitive to arbitrary bin phase boundaries and unable to measure true per-cluster sequence divergence for `MdlCostModel`.
- **Rejected & removed:** Replaced by `resolve_alignment_duplications_mdl` + `prune_resolved_dup_clusters`, and the legacy modulo-bin helpers were removed from the public API.

## Consequences
- `horizon_wga --resolve-dups` and `horizon_resolve` handle arbitrary $k$-copy tandem arrays and inverted repeats under a unified, saturation-safe MDL cost model.
- Active-cluster sweeping eliminates quadratic $O(G^2)$ clustering overhead on dense whole-chromosome alignments.
