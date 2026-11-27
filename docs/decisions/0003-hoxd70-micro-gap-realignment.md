# ADR-003: HOXD70 Micro-Gap Realignment, `auN` Hinge Selection, and Global 1-to-1 Netting

## Status
Accepted (updated 2026-11-24 with Wave B–D `MAX_COLINEAR_GAP_BP = 40`, neutral `N`, `auN` hinge, $O(R \log R)$ global netting, and final `e095e6d` benchmark verification)

## Date
2026-11-24

## Context
Real unstripped multi-genome benchmarking against UCSC GoldenPath `wuhCor1` (`multiz44way`, `multiz119way`) and `sacCer3/multiz7way` (`chrIII` & `chrI`, 7 byte-exact draft yeast assemblies) exposed three concrete algorithmic issues:

1. **Balanced `X- / -Y` skip artifacts in MAF serialization (`src/formats/realign.rs`):** Unaligned equal-length mismatch stretches (`gx == gy > 0`) between consecutive colinear `AlignedRun`s were emitted as adjacent target-gap and query-gap blocks, inflating `StatSigMA-w` unsupported window rates (`28.84%` on 44-way TBA; `49.98%` on *S. kluyveri* TBA) and fragmenting GFA/HAL graphs (`46,570` segments on 44-way TBA; `284,456` on yeast `chrIII` TBA). An initial bounded Gotoh micro-realigner using ad-hoc weights (`+2 / -2 / -5 / -1`) eliminated skips but shattered diverged 4-bp coding windows and single 3-bp deletions into paired 1–2 bp frameshift indels, dropping 44-way CDS net in-frame indels from `41.0%` to `17.9%`.
2. **Unnormalized hinge selection on fragmented assemblies (`src/bin/horizon_guide/main.rs`):** Summing raw MinHash sketch overlaps across all contigs selected fragmented outgroups with larger flanking totals (e.g., `sacCas` on yeast `chrIII`, 45 contigs / `2.25 Mbp`), collapsing close-sister F1 from `0.89` to `0.30`.
3. **Per-contig-pair `--one-to-one` scoping (`src/bin/horizon_wga.rs`):** Applying 1-to-1 filtering independently to each `(contig_x, contig_y)` pair allowed duplicated loci on distinct draft contigs (e.g., `HMLα` on yeast `chrIII`) to retain `1.62×` cross-contig multiplicity.

## Decision
1. **Zero-alloc stack HOXD70 micro-gap realignment (`src/formats/realign.rs:align_micro_gap_affine`):** Use the standard evolutionary log-odds **HOXD70** DNA substitution matrix and affine gap penalties (Chiaromonte et al., 2002: `A/T = +91`, `C/G = +100`, transitions `-31`, transversions `-114 / -123 / -125`, `gap_open = -400`, `gap_ext = -30`, neutral `0` for ambiguous `N`/non-ACGT bases) with an inline `1×1` SNP fast path over all inter-run gaps up to `MAX_COLINEAR_GAP_BP = 40` (matching `src/netting.rs`).
2. **`auN`-based contiguity-aware hinge selection & explicit override (`src/bin/horizon_guide/main.rs`):** Add `--hinge <species>` and normalize sketch-overlap centrality by `1 + ln(n_eff)` where `n_eff = (sum L_i)^2 / sum(L_i^2) = total_bp / auN` (Li 2020) — identically `1.0` for single-contig chromosomes, resilient to minor unplaced scaffolds (e.g. `hg38`), and scaling smoothly with true fragmentation.
3. **Species-wide $O(R \log R)$ global 1-to-1 netting (`src/netting.rs:net_global_one_to_one_runs_scored`):** Enforce global non-overlapping target and query coordinates across all contigs in a species pair using `IntervalSet` (`O(R log R)` time, `O(R)` memory), ranking colinear run-groups by canonical ACGT identity score (`2 * matches - mismatches`) with tie-breaks by syntenic contig-pair total score and group bp.

## Alternatives Considered

### Keep ad-hoc `+2 / -2 / -5 / -1` micro-realignment weights
- **Cons:** Opening two compensating 1-bp gaps (`-12`) beat 4 mismatches (`-8`), introducing spurious frameshift indels inside diverged coding micro-bubbles (`17.86%` 44-way CDS net in-frame indels).
- **Rejected:** Violated first-principles evolutionary log-odds scaling and degraded CDS reading-frame preservation.

### Normalize hinge centrality by raw contig count `1 + 0.25 * ln(n_contigs)`
- **Cons:** Penalizes chromosome-level reference assemblies that carry hundreds of tiny unplaced contigs (e.g., `hg38` has 24 chromosomes + ~431 tiny unplaced scaffolds) as heavily as genuinely fragmented draft assemblies.
- **Rejected:** Replaced in Wave B with length-weighted `auN` effective contig count `total_bp / auN`.

### Keep raw `X- / -Y` skip emission without micro-realignment
- **Cons:** Left `64k–149k` balanced skips across 44-way Sarbecoviruses, inflated TBA `StatSigMA-w` unsupported rates to `28.8%–50.0%`, and bloated HAL/GFA artifacts `8×–36×`.
- **Rejected:** Unacceptable for downstream comparative genomics and pangenome graph consumers.

## Consequences
- **`wuhCor1` verified (`e095e6d`):** 44-way `guide_prog` CDS net in-frame indels rose to **`60.00%`** (**beating UCSC's `59.09%`**), balanced skips fell to **`4`** (`0.0` median, half of UCSC's `8`), exact F1 reached **`0.9679`** (`tba` **`0.9805`**, `±50 bp` **`0.9899`**), 44-way TBA GFA shrank to **`813` segments** (`96/96` exact walks), HAL shrank to **`5.13 MB`**, and 119-way `d ≥ 1.5` `StatSigMA-w` dropped to **`2.83%`** (`12.5×` below UCSC).
- **`sacCer3` verified (`e095e6d`):** `horizon_guide` automatically selects `sacCer3`, balanced skips across all 185 `chrIII` genes dropped to exact **`0 / 0 / 0 / 0 / 0 / 0`** across all 6 yeasts, `HMLα`/`MATα`/`HMRa` hit exact **`1.000×`** multiplicity across all 41 draft contigs, `wga --one-to-one` exact F1 improved on **6/6 yeasts** (`chrIII` `0.908`–`0.616`; `chrI` `+2.6` to `+4.6 pp`), and `StatSigMA-w` fell to **`0.12%–5.80%`** (`4.1×–15.1×` below UCSC).
