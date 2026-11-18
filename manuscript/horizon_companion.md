---
title: "Horizon Technical Companion: Data Representations, Mathematical Foundations, and Worked Examples"
author:
  - "Richard Burhans (ORCID 0009-0006-5949-0958)"
date: "2026-11-18"
bibliography: references.bib
---

# 1. Purpose and Scope

This companion document provides the data-format specifications, mathematical derivations, and implementation details underlying the Horizon whole-genome aligner (`https://github.com/richard-burhans/horizon-aligner`). While the main manuscript focuses on biological motivation, system architecture, and empirical benchmarks, this reference is designed for engineers, reviewers, and contributors inspecting the Rust kernels (`src/lib.rs`, `src/egraph.rs`, `src/hal.rs`, and `src/bin/`).

---

# 2. Data Representations and Wire Formats

## 2.1 Pairwise and Multi-Species `MAF`

Horizon emits standard `##maf version=1` alignment blocks:

```text
a score=184.2
s Pf3D7_01_v3   12040 42 + 638193 ATATATATATCCGGTTAA--CCGTAATATATATAT...
s Pr_01         11890 44 + 584285 ATATATATATCCGGTTAAAACCGTAATATATATAT...
```

- **Coordinates:** `start` is 0-based on the reported `strand` (`+` or `-`). For `-` strand rows, `start` is measured from the reverse-complemented sequence start (`seq_len - (fwd_start + size)`), matching UCSC/Galaxy `MAF` conventions.
- **Interval netting (`horizon_merge`):** When multiple pairwise blocks overlap on the reference coordinate space, `interval_match_score` scores each interval by exact matches minus mismatch/gap penalties and retains the highest-scoring non-overlapping tiling before threading query insertion columns (`thread_interval_columns`).

## 2.2 Reference-Anchored `rGFA` and Lens Metadata

Graph outputs (`horizon_wga --gfa` and `horizon_merge`) use GFA v1.1 / `rGFA` conventions:

- **Segment (`S`) lines:** Carry sequence and reference-coordinate tags:
  - `SN:Z:<seq_name>` — source contig or chromosome name
  - `SO:i:<offset>` — 0-based start offset on `SN`
  - `SR:i:0` — rank `0` for reference-backbone segments
- **Horizon Bubble (`# HB:Z:`) annotations:** Unaligned intervals where `ΔL < λ̂` or inter-LCB spans exceed the dynamic-programming ceiling are emitted as explicit bubble records carrying `(ref_start, ref_end, qry_start, qry_end, reason)` rather than forced base-by-base alignments.
- **E-Graph Paralog (`# EG:Z:`) annotations:** Collapsed tandem duplications and nested inversions resolved by `src/egraph.rs` are recorded as `# EG:Z:` tags so downstream tools can distinguish collapsed paralogs from single-copy synteny.
- **Bidirectional `GFA` ↔ `MAF` Lens:** Viewing `MAF` as an asymmetric view `A = get(C)` of the full homology graph `C` (which also stores bubbles, paralogs, and non-reference blocks), edits on the `MAF` side map back via `putback(A', C)` while obeying:
  1. **GETPUT:** `putback(get(C), C) == C`
  2. **PUTGET:** `get(putback(A', C)) == A'`
  3. **Constant Complement:** `# HB:Z:` bubbles, `# EG:Z:` paralogs, and non-reference blocks are preserved unchanged across round-trips.

## 2.3 Pure-Rust HDF5 and `HAL` Binary Layout (`src/hal.rs`)

`HAL` files (Hickey et al. 2013) are HDF5 containers storing a phylogenetic tree and per-genome sequence/segment tables. `src/hal.rs` parses HDF5 `v0`/`v1` containers natively in pure Rust without linking `libhdf5`:

1. **HDF5 Superblock (`read_hdf5_superblock`):** Checks the 8-byte magic signature (`89 48 44 46 0d 0a 1a 0a`), `offset_size`, `length_size`, and `root_group_addr`.
2. **Group Traversal (`read_hdf5_v0_group_entries`):** Reads the `SymbolTable` (`0x0011`) message in any group's `v1` object header to locate its `TREE` (B-tree `v1`, type `0`), `SNOD` symbol-table nodes (40-byte entries containing `(name_off, obj_hdr_addr)`), and `HEAP` local string heap.
3. **Dataset Reading (`read_hdf5_v0_dataset_bytes`):**
   - Inspects `Dataspace` (`0x0001`) for rank and dimensions, `Datatype` (`0x0003`) for fixed element width, and `FilterPipeline` (`0x000B`) for `zlib` deflate (`filter_id == 1`).
   - Reads either contiguous `DataLayout` (`0x0008`, class `1`) slices or walks chunk B-trees (`DataLayout` class `2`, `TREE` type `1`) and decompresses each chunk via `flate2::read::ZlibDecoder`.
4. **HAL Binary Record Decoders:**

| HAL Dataset | Record Width | Binary Field Layout (Little-Endian) | Decoder |
| :--- | :--- | :--- | :--- |
| `DNA_ARRAY` | 4 bits / base | High nibble first: `0..3` = `a,c,g,t`, `8..11` = `A,C,G,T`, `4` = `n`, else `N` | `decode_hal_dna` |
| `SEQIDX_ARRAY` | `24` bytes | `(start: u64, top_offset: u64, bottom_offset: u64)` | `decode_hal_seq_index` |
| `SEQNAME_ARRAY` | `stride` bytes | Fixed-width null-terminated ASCII/UTF-8 string | `decode_hal_seq_names` |
| `TOP_ARRAY` | `33` bytes | `(genome_idx: i64, bottom_idx: i64, paralogy_idx: i64, parent_idx: i64, rev: u8)` | `decode_hal_top_array` |
| `BOTTOM_ARRAY` | `24 + 9C` bytes | `(start: i64, len: u64, top_parse: i64)` + `C` × `(child_idx: i64, rev: u8)` | `decode_hal_bottom_array` |

---

# 3. Mathematical and Algorithmic Foundations

## 3.1 Strand-Symmetric `KP-14` Eisenstein Seeding

To avoid canonicalizing forward and reverse-complement *k*-mers separately, Horizon encodes nucleotides over the quadratic extension field **F**_{q²} (*q* = 2⁶¹ − 1) using sixth-root Eisenstein units:

```text
φ(A) = 1,   φ(T) = -1,   φ(C) = ω,   φ(G) = -ω
```

Because Watson–Crick complements map to exact algebraic negatives (`φ(b̄) = −φ(b)`), the forward polynomial `F` and reverse-complement polynomial `R` over any spaced-seed mask `M` satisfy `(F + λ)(R + λ)` invariance under strand reversal. Horizon pairs `KP-14` spaced seeds (`w = 12` and `w = 10`) with:
- **`LASTZ`-style transition seeds:** `MASK_LASTZ11` (1 transition) and `MASK_LASTZ12` (*T* ∈ {1, 2} transitions), gated by ungapped high-scoring pair (HSP) extension.
- **Reduced-alphabet nucleotide seeds:** `MASK_RY20` (purine/pyrimidine) and `MASK_W10` (weak/strong H-bond), activated on HSP-supported regions.
- **Six-frame biochemical amino-acid seeds (`5-of-6`):** Translates all six reading frames into 9 biochemical equivalence classes (`aa_group`) to recover coding homology where synonymous third-codon positions have saturated.

Seeds are chained in O(*M* log *M*) time via Longest Increasing Subsequence (LIS) dynamic programming and grouped into bidirectional locally collinear blocks (LCBs).

## 3.2 Online Compositional Tracking (`BOCPD`)

Instead of fixed sliding windows or vertebrate-trained matrices (`HOXD70`), Horizon models local AT/GC composition along the binary weak/strong channel `x_t ∈ {0, 1}` (`1` for `A/T`, `0` for `G/C`) using Beta-Bernoulli Bayesian Online Changepoint Detection (Adams & MacKay 2007). At each position *t*, the run-length distribution `P(r_t | x_{1:t})` updates recursively with constant hazard *h*:

```text
P(r_t = r_{t-1} + 1, x_{1:t}) ∝ (1 - h) · P(r_{t-1}, x_{1:t-1}) · π(x_t | r_{t-1})
P(r_t = 0,         x_{1:t}) ∝ h · Σ_{r_{t-1}} P(r_{t-1}, x_{1:t-1}) · π(x_t | r_{t-1})
```

where the Beta-Bernoulli predictive probability `π(x_t = 1 | r) = (α_0 + s_r) / (α_0 + β_0 + r)` adapts smoothly inside homogeneous isochores and resets immediately at sharp `18%` ↔ `45%` GC boundaries.

## 3.3 Banded Three-State MDL Forward-Backward Alignment

For each inter-anchor gap (`X`, `Y`), Horizon evaluates a three-state Minimum Description Length pair-HMM (`M` match/mismatch, `I` insertion, `D` deletion) in log-space (`logplus`), parameterized by the local BOCPD base frequencies. To eliminate Python's sparse-dictionary bottleneck, forward (α) and backward (β) tables are stored in a contiguous flat `BandGrid` (`Vec<f64>`) of stride `2W + 1` (`band_w = 180`), yielding:

- **Joint description length:** `L_joint = −log₂ P(X, Y | aligned)`
- **Local null description length:** `L_null = −log₂ P_BOCPD(X) − log₂ P_BOCPD(Y)`
- **Bit savings:** `ΔL = L_null − L_joint`
- **Per-column posterior:** `P(M_ij | X, Y) = exp(α_M(i, j) + β_M(i, j) − log P(X, Y))`

## 3.4 Conformal Risk Control and Horizon Bubbles

Horizon treats alignment emission as a calibrated statistical decision. Using Conformal Risk Control (Angelopoulos et al. 2024) at risk level `α = 0.07` (`n_cal = 15`), any block or gap whose bit savings falls below `λ̂ = 10.0` bits—or whose column posterior drops below `0.15`—is withheld from base-level alignment and emitted as a **Horizon Bubble**.

## 3.5 Congruence-Closed E-Graph Duplication Resolution

Horizon's `src/egraph.rs` models LCB sequences as e-classes in a congruence-closed e-graph (Willsey et al. 2021) and saturates four non-destructive algebraic rewrites:

1. `inv(inv(X)) <-> X` (double-inversion cancellation)
2. `inv(seq(A, B)) <-> seq(inv(B), inv(A))` (anti-homomorphic reversal)
3. `seq(seq(A, B), C) <-> seq(A, seq(B, C))` (associativity)
4. `seq(X, X) -> dup(X)` (tandem duplication collapse)

Bottom-up parsimony extraction then selects the minimum-cost ancestral block representation (`dup(X)`), recovering duplicated and inverted regions that greedy 1-to-1 synteny filters discard.

---

# 4. Worked Toy Example (20 bp Trace)

Consider a 20 bp *Plasmodium*-like pair with a conserved `15%`-GC left flank, a 6 bp unalignable divergent insert, and a conserved right flank:

- **Target `X` (20 bp):** `AAATAAATAA` `GGGGCC` `AAAT`
- **Query `Y` (20 bp):** `AAATAAATAA` `CCTTTT` `AAAT`

1. **Seeding & Chaining:** `KP-14` (`w = 10`/`12`) matches the identical `AAATAAATAA` prefix (`0..10`) and `AAAT` suffix (`16..20`), leaving an inter-anchor gap of `6 bp` (`X[10..16] = GGGGCC` vs. `Y[10..16] = CCTTTT`).
2. **BOCPD Tracking:** Over `X[0..10]` (`100%` A/T), `bocpd_gc_profile` estimates `P(W) ≈ 0.91`. At `X[10] = G`, the changepoint hazard shifts `P(W)` downward across the GC-rich insert.
3. **Banded MDL Scoring:**
   - On the conserved prefix (`10 bp`), matching `A`/`T` under `P(W) ≈ 0.91` still yields `ΔL ≈ +14.2` bits (`≥ λ̂ = 10.0` bits) with column posteriors `P(M_ij) > 0.95`.
   - Inside the `6 bp` gap (`GGGGCC` vs. `CCTTTT`), 0 of 6 bases match, giving `ΔL ≈ −6.4` bits (`< 10.0` bits) and `P(M_ij) < 0.15`.
4. **Conformal Gate & Output:** Horizon emits the aligned `MAF` blocks for `0..10` and `16..20`, and isolates `10..16` as an explicit `# HB:Z:` **Horizon Bubble** in the `GFA` graph rather than forcing 6 mismatched bases into a false alignment.

---

# 5. References {.unnumbered}

Full BibTeX entries for all cited works are maintained in [references.bib](references.bib).

- **Adams RP & MacKay DJC** (2007). Bayesian online changepoint detection. *arXiv:0710.3742*.
- **Allison L, Wallace CS, Yee CN** (1992). Finite-state models in the alignment of macromolecules. *Journal of Molecular Evolution*, 35(1):77–89.
- **Angelopoulos AN, Bates S, Fisch A, Lei L, Schuster T** (2024). Conformal risk control. *ICLR 2024* (`arXiv:2208.02814`).
- **Foster JN, Greenwald MB, Moore JT, Pierce BC, Schmitt A** (2007). Combinators for bidirectional tree transformations: A linguistic approach to the view-update problem. *ACM TOPLAS*, 29(3):17.
- **Hickey G, Paten B, Earl D, Zerbino D, Haussler D** (2013). HAL: a hierarchical format for storing and analyzing multiple genome alignments. *Bioinformatics*, 29(10):1341–1342.
- **Willsey M, Nandi C, Wang YR, Flatt O, Tatlock Z, Panchekha P** (2021). egg: Fast and extensible equality saturation. *PACMPL*, 5(POPL):1–29.
