# Horizon Aligner Architecture Audit

**Verdict:** The architecture of `horizon-aligner` is structurally sound, exhibiting clean unidirectional dependency flow, well-encapsulated algorithm modules, and exceptional baseline hygiene (zero `clippy` warnings on strict settings). There are a few cross-tier boundary refinements and consistency proposals to address.

## 1. Boundary Integrity
- **Overall:** Excellent. `src/mdl/` correctly consumes `src/seeding/` without circular dependencies, and the `bin/` drivers consume everything cleanly.
- **Finding (Low):** `src/merge.rs` calls `resolve_lcb_duplications`, which is exported by `src/formats.rs`. However, structural duplication resolution is the domain of `src/egraph.rs`. Placing e-graph resolution wrappers inside the formatting tier forces `merge.rs` to import from `formats.rs` when it shouldn't need to.
  - *[proposal]* Move `resolve_lcb_duplications` and `resolve_alignment_duplications` out of `formats.rs` and into `egraph.rs` (or a dedicated `resolution.rs` module).

## 2. Abstraction Consistency
- **Overall:** `ContigAlignment` provides a solid universal representation for finalized alignments.
- **Finding (Medium):** The `hal.rs` module is asymmetric compared to `formats.rs`. While `formats.rs` provides generalized writers and block formatters (MAF/GFA) from `ContigAlignment`s, `hal.rs` is highly coupled to low-level HDF5 decompression and chunk reading, acting mostly as a decoder. It doesn't yet have a clean unifying adapter to convert HAL trees into `ContigAlignment` objects natively.
  - *[accepted]* This is intentional asymmetry for now, as HAL requires hierarchical tree resolution which maps poorly to a flat `ContigAlignment`.

## 3. Naming / Vocabulary Drift
- **Finding (Low):** "Horizon Bubbles" are referenced in `formats.rs` (`macro_bubbles` and `internal_bubbles`) but the term isn't strictly defined in the `mdl` DP layer; they are emitted as a consequence of gap properties. The vocabulary between the algorithmic "gap" and the formatted "bubble" is slightly disjointed.
  - *[proposal]* Add documentation in `mdl/mod.rs` connecting `GapResult` to the concept of "Horizon Bubbles" to bridge the terminology.

## 4. Contract-Enforcement Gaps
- **Finding (Low):** The 0-based half-open coordinate system (with reverse strands relative to forward length) is strictly adhered to, but coordinate projection is performed manually inline in `formats.rs` and `merge.rs` rather than through a dedicated coordinate-translation abstraction. This risks off-by-one errors if a new output format is added.
  - *[proposal]* Introduce a `Coordinate` or `Interval` struct that inherently manages strand inversion and length calculations to enforce this invariant systematically.

## 5. Duplication / Missed Reuse
- No significant code duplication detected. The logic for sequence reading (`read_fasta_contigs` vs HDF5 arrays) is distinctly specialized.

## 6. Dead / Reserved Surface
- **Finding (None):** The codebase passes strict `cargo clippy --all-targets -- -D warnings` with zero dead code or unused dependencies.

## 7. Doc / Code Agreement
- **Finding (Fixed):** `ARCHITECTURE.md` was missing, meaning there was no written baseline for these contracts.
  - *[fixed]* Authored `ARCHITECTURE.md` in the repository root to establish the tier map and cross-cutting contracts.