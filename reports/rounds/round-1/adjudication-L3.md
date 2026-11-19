# Round 1 ADJUDICATION -- L3 (Statistics)

- **Attack verified:** `reports/rounds/round-1/attack-L3-statistics.md` (sha256 prefix `43eb59fcf98a` verified).
- **Prior rounds:** `scripts/round history "MDL" --kind adjudication` confirmed no closed prior rounds exist.
- **Method cards:** No `refs/methods/` directory exists in this repository; rulings are grounded directly in information theory, probability normalization, and the verified source code.

---

## F1 -- C2 (`docs/STRATEGIC_PLAN.md:9`)

- **Verdict:** **UPHELD**
- **What survives:** Claim C2's assertion that arbitrary affine gap penalties and scoring matrices are replaced by MDL bit-savings that inherently prevent hairballs is refuted by hardcoded gap/match constants, unnormalized HMM transition weights, a `0.55` search-space multiplier, and unbenchmarked HPRC graph evaluation.
- **Verification of diagnosis:**
  1. `src/bocpd.rs:4-7` hardcodes `GAP_OPEN_BITS = -4.5`, `GAP_EXTEND_BITS = -0.8`, `IDENTITY_PRIOR = 0.65`, and `MIN_POST_PROB = 0.15`. `src/mdl/mod.rs:566-567` scores `< 8`-anchor LCBs with a fixed linear rule `1.6 * matches - 2.4 * mismatches`, and `src/mdl/mod.rs:158-161, 180-182, 310, 344, 378, 398` uses ad-hoc linear HSP thresholds.
  2. In `src/mdl/dp.rs:105-107, 132-167`, `exp(go) = 2^(-4.5) ≈ 0.04419` and `exp(ge) = 2^(-0.8) ≈ 0.57435`. Transition `M -> M` (`dp.rs:162-167`) has weight `1.0`, while `M -> Ix` and `M -> Iy` each have weight `2^(-4.5)`, giving outgoing mass `1 + 2 * 2^(-4.5) ≈ 1.0884 > 1` from `M`; similarly `Ix -> M` (`1.0`) plus `Ix -> Ix` (`2^(-0.8)`) gives `1.5743 > 1` from `Ix` (and `Iy`), so `z_fwd` (`dp.rs:172`) is not a normalized partition function and `post_p` (`dp.rs:240`) is not a normalized posterior.
  3. `src/mdl/mod.rs:153-161` sets `effective_min_bits = cfg.min_delta_bits.max(0.55 * search_space_bits)`, charging `0.55 * log2(N * M)` rather than `log2(N * M)`.
  4. `docs/BENCHMARK_RESULTS.md:39` leaves HPRC as `*(Graph Eval)*` with no numbers, and `cactus241way` (chr22 1MB, line 37) reports `9.9%` recall with `494` Horizon Bubbles.
- **Ruling on proposed correction:** Endorsed with a minor addition: the proposed replacement text in `docs/STRATEGIC_PLAN.md:9` is accurate and installs no false claim, though it should also note the `1.6 / -2.4` short-LCB linear gate (`src/mdl/mod.rs:566`). **The fix is free** (documentation-only edit, no re-measurement needed).

---

## F2 -- C3 (`docs/STRATEGIC_PLAN.md:10`)

- **Verdict:** **UPHELD**
- **What survives:** `src/egraph.rs` uses fixed unitless parsimony weights (`1.0 / 1.2 / 1.3`) with no evolutionary trajectory probability model and zero HLA or paralog benchmark evaluation.
- **Verification of diagnosis:** `src/egraph.rs:291-308` (`extract_min_cost`) assigns `blk = 1.0`, `inv = c + 1.2`, `dup = c + 1.3`, `seq = c1 + c2 + 1.0`. `docs/BENCHMARK_RESULTS.md:20-40` contains no HLA or structural-paralog benchmark, and `src/egraph.rs:376-384` tests only a single 5-block synthetic example.
- **Ruling on proposed correction:** Endorsed. Narrowing `docs/STRATEGIC_PLAN.md:10` to describe `src/egraph.rs` as a bounded parsimony rewrite prototype not yet parameterized by an evolutionary model or benchmarked on HLA is accurate. **The fix is free** (documentation-only edit).

---

## F3 -- C4 (`docs/STRATEGIC_PLAN.md:11`)

- **Verdict:** **UPHELD**
- **What survives:** No empirical benchmark, RSS measurement, or crash reproduction exists anywhere in the repository for `vg`, `odgi`, or `src/gfa/mod.rs` on 12GB pangenomes, and `src/gfa/mod.rs` is not zero-copy.
- **Verification of diagnosis:** Repo-wide grep across `docs/` and `scripts/` confirms no 12GB benchmark or `vg`/`odgi` memory comparison exists. `src/gfa/mod.rs:86-116` allocates owned `String` and `Vec<u8>` heap buffers per record via `.to_string()` and `.as_bytes().to_vec()`.
- **Ruling on proposed correction:** Endorsed, with one addition: in addition to removing the unsubstantiated `12GB` / `vg` / `odgi` claim and `"solves the I/O bottleneck entirely"`, `docs/STRATEGIC_PLAN.md:11` should also replace `"zero-copy"` with `"line-streaming"` to match `src/gfa/mod.rs:86-116`. **The fix is free** (documentation-only edit).

---

## F4 -- C7 (`docs/STRATEGIC_PLAN.md:27`)

- **Verdict:** **NARROWED**
- **What survives:** `saturate_synteny(&mut eg, 4, 800)` (`src/egraph.rs:165-273, 381`; `src/formats.rs:66`) is hard-capped at `max_iters = 4` and `max_nodes = 800`, so global optimality is not guaranteed on arbitrary bubbles, and HPRC bubble sizes are unmeasured (`docs/BENCHMARK_RESULTS.md:39`).
- **Overturned leg (premise & arithmetic):** The attacker's claim that saturation truncates at `k >= 8` blocks because of Catalan growth (`C_7 = 429, C_8 = 1,430 > 800`) fails on two counts:
  1. By the attacker's own numbers, `C_7 = 429 < 800`, so `k = 8` would not exceed `800`.
  2. More fundamentally, an e-graph hash-conses shared sub-intervals in `eg.memo`. Under pure `seq` associativity on `k` distinct blocks, every contiguous sub-interval `[i..j]` of length `L >= 2` shares one e-class with `L - 1` `seq` e-nodes, yielding `sum_{L=2..k} (k - L + 1)(L - 1) + k = C(k+1, 3) + k` total `ENode`s (`92` at `k = 8`, `129` at `k = 9`, exceeding `800` only at `k = 17` with `833` nodes), rather than `C_{k-1}` e-nodes. (`inv`/`dup` rewrites and `max_iters = 4` can still limit or expand saturation depending on the bubble.)
  - *What the attacker should have checked instead:* The hash-consed `ENode` count `C(k+1, 3) + k` in `eg.memo` and the `max_iters = 4` iteration cap, plus `C_7 = 429 < 800`.
- **Ruling on proposed correction:** Rejected on the `typically k <= 7 blocks` wording, which installs a false threshold. Replace with a statement that `saturate_synteny` is bounded at `max_iters = 4, max_nodes = 800` and does not guarantee global optimality on large bubbles. **The fix is free** (documentation-only edit).

---

## F5 -- C9 (`docs/PANGENOME_SUITE.md:8`)

- **Verdict:** **UPHELD**
- **What survives:** Neither `src/bocpd.rs` nor `src/mdl/` defines a statistical likelihood or bit-savings function for GFA `L` lines (`GfaLink`), and `0M` links carry 0 aligned bases.
- **Verification of diagnosis:** `bocpd_gc_profile` (`src/bocpd.rs:58-119`) takes a 1D `&[u8]` sequence; `banded_mdl_gap` (`src/mdl/dp.rs:81-309`) aligns two `&[u8]` slices. `GfaLink` (`src/gfa/mod.rs:21-27`) stores segment IDs, orientations, and an `overlap` string (`0M` in `src/merge.rs:246` and `src/formats.rs:51-54`).
- **Ruling on proposed correction:** Endorsed. Clarifying in `docs/PANGENOME_SUITE.md:8` that evaluating `L` lines is unimplemented and requires specifying what flanking-sequence or read-support data would be scored is accurate. **The fix is free** (documentation-only edit).

---

## F6 -- C10 (`docs/PANGENOME_SUITE.md:16`)

- **Verdict:** **UPHELD**
- **What survives:** `BandGrid` (`src/mdl/dp.rs:19-72`) implements only a 2D linear-diagonal band `center(i) = (i * m) / n` and has no graph-topological DP across bifurcations.
- **Verification of diagnosis:** Verified in `src/mdl/dp.rs:19-72` (`center` at line 42, `slot` at lines 47-59) and `banded_mdl_gap` / `mdl_xdrop_extend`; no DAG or multi-predecessor aggregation exists.
- **Ruling on proposed correction:** Endorsed. Stating in `docs/PANGENOME_SUITE.md:16` that `src/mdl/dp.rs` currently supports only linear sequence-to-sequence `BandGrid` DP and requires a graph-topological formulation for bifurcations is accurate. **The fix is free** (documentation-only edit).

---

## F7 -- Scope note on `docs/BENCHMARK_RESULTS.md:5`

- **Verdict:** **UPHELD**
- **What survives:** `docs/BENCHMARK_RESULTS.md:5` claims `"parameter-free, Information-Theoretic Minimum Description Length (MDL) threshold defaults"`, which is contradicted by hardcoded priors and thresholds in `src/bocpd.rs:4-7, 20-25` and `src/mdl/mod.rs:52-70`.
- **Verification of diagnosis:** `src/bocpd.rs:6` fixes `IDENTITY_PRIOR = 0.65`, `src/bocpd.rs:20-25` assumes `G=C` and `A=T` strand parity, and `AlignConfig::default()` (`src/mdl/mod.rs:52-70`) fixes `band_w = 180`, `max_ext = 1500`, `min_delta_bits = 10.0`, `fallback_min_anchors = 4`, `fallback_min_span = 50`, `fallback_max_occ = 2`.
- **Ruling on proposed correction:** Removing `"parameter-free"` from `docs/BENCHMARK_RESULTS.md:5` is warranted. **The fix is free** (documentation-only edit).

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | Hardcoded gap/match constants, unnormalized HMM transitions (`1.0884`, `1.5743 > 1`), `0.55 * log2(N*M)`, and unbenchmarked HPRC refute C2; free text repair endorsed. |
| F2 | UPHELD | `src/egraph.rs:291-308` uses fixed unitless parsimony weights (`1.0/1.2/1.3`) with no evolutionary model or HLA evaluation; free text repair endorsed. |
| F3 | UPHELD | No 12GB/`vg`/`odgi` benchmark exists and `src/gfa/mod.rs:86-116` allocates owned `String`/`Vec<u8>` per record; free text repair endorsed (also drop `"zero-copy"`). |
| F4 | NARROWED | `max_iters = 4, max_nodes = 800` truncation and unmeasured HPRC bubble sizes hold, but `k >= 8` / `C_7 = 429 > 800` Catalan e-node bound is overturned (`C(k+1,3)+k` hash-consed e-nodes); `k <= 7` repair rejected. |
| F5 | UPHELD | `src/bocpd.rs` and `src/mdl/dp.rs` score 1D/2D nucleotide slices and have no statistical model for GFA `L` lines (`0M` overlaps); free text repair endorsed. |
| F6 | UPHELD | `BandGrid` (`src/mdl/dp.rs:19-72`) is a 2D linear-diagonal band undefined across graph bifurcations; free text repair endorsed. |
| F7 | UPHELD | `docs/BENCHMARK_RESULTS.md:5` `"parameter-free"` claim is contradicted by fixed priors and thresholds in `src/bocpd.rs` and `src/mdl/mod.rs`; free text repair endorsed. |
