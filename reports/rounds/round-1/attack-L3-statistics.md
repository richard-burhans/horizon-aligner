# Round 1 Attack -- L3 (Statistics)

**Scope & Provenance Checks:**
- Prior rounds: `scripts/round history "MDL"` confirmed Round 1 is the first round (no closed prior rounds).
- Method cards: No method-card directory exists in this repository (`docs/`, `reports/`, `manuscript/`).
- Lens: Statistics, information theory (MDL), probabilistic modeling (BOCPD/HMM), and quantitative/empirical validation. Claims outside this lens (architectural/CLI/format/coordinate-type claims **C1, C5, C6, C8, C11**) are marked `UNCHALLENGED` here.

---

## F1 -- C2 REFUTED: Hardcoded affine gap penalties, fixed priors, and an unnormalized HMM violate MDL

- **Claim:** **C2** (`docs/STRATEGIC_PLAN.md:9`) -- *"We are replacing arbitrary affine gap penalties and scoring matrices with Information-Theoretic bit-savings (MDL). This inherently prevents the 'hairball' over-alignment problem."*
- **Verdict:** `REFUTED`
- **Artifacts checked:**
  - `src/bocpd.rs:4-7, 28-39`
  - `src/mdl/dp.rs:105-173, 240`
  - `src/mdl/mod.rs:52-70, 153-162, 180-192, 566-567`
  - `docs/BENCHMARK_RESULTS.md:5, 37, 39`
- **Concrete flaws & evidence:**
  1. **Affine gap penalties and linear scoring matrices are hardcoded, not replaced:**
     - `src/bocpd.rs:4-7` defines fixed constants `GAP_OPEN_BITS = -4.5`, `GAP_EXTEND_BITS = -0.8`, `IDENTITY_PRIOR = 0.65`, and `MIN_POST_PROB = 0.15`.
     - `src/mdl/mod.rs:566-567` gates short LCBs with a fixed linear scoring matrix: `net_bits = 1.6 * (matches as f64) - 2.4 * (mismatches as f64)`.
     - `src/mdl/mod.rs:158-161, 180-182` sets integer HSP thresholds using ad-hoc linear scalings (`100.0 + effective_min_bits * 1.5`, `(search_space_bits - 20.0) * 3.8`, plus `+ 30` / `+ 40` at lines `310, 344, 378, 398`).
  2. **Forward-Backward HMM transition weights sum to > 1 (violating probability normalization / Kraft's inequality):**
     - In `src/mdl/dp.rs:105-107`, `go = GAP_OPEN_BITS * ln2` ($2^{-4.5} \approx 0.0442$) and `ge = GAP_EXTEND_BITS * ln2` ($2^{-0.8} \approx 0.5743$).
     - From match state `M` (`fm`), the transition `M -> M` (`dp.rs:162-167`) has weight $1.0$ (no $\log(1 - 2\delta)$ penalty), while `M -> Ix` (`dp.rs:132-136`) and `M -> Iy` (`dp.rs:142-146`) each have weight $2^{-4.5} \approx 0.0442$. Outgoing probability mass from `M` is $1 + 2 \cdot 2^{-4.5} \approx 1.0884 > 1$.
     - From gap state `Ix` (`fix`), `Ix -> M` (`dp.rs:164`) has weight $1.0$ and `Ix -> Ix` (`dp.rs:135`) has weight $2^{-0.8} \approx 0.5743$, giving outgoing mass $1.5743 > 1$.
     - Consequently, `z_fwd` (`dp.rs:172`) is not a normalized marginal likelihood under a valid generative HMM, invalidating both the MDL bit-savings interpretation and `post_p = (fm + bm - z_fwd).min(0.0).exp()` (`dp.rs:240`).
  3. **Under-penalized search-space description length:**
     - `src/mdl/mod.rs:153-161` computes `search_space_bits = log2(N * M)` and sets `effective_min_bits = cfg.min_delta_bits.max(0.55 * search_space_bits)`. Specifying a start coordinate pair in an $N \times M$ matrix requires $\log_2(NM)$ bits ($1.0 \times \text{search\_space\_bits}$); the `0.55` multiplier under-penalizes the search-space cost by 45%.
  4. **No empirical demonstration of preventing pangenome "hairball" over-alignment:**
     - `docs/BENCHMARK_RESULTS.md:39` lists HPRC graph evaluation as `*(Graph Eval)*` with no numbers; on `cactus241way` (chr22 1MB, Human vs Dog, line 37), recall drops to `9.9%` with `494` Horizon Bubbles.
- **Correction:** Replace the claim in `docs/STRATEGIC_PLAN.md:9` with: *"We replace traditional substitution matrices with a GC-aware BOCPD local-composition log-likelihood ratio (`IDENTITY_PRIOR = 0.65`) combined with fixed affine gap penalties (`-4.5` open, `-0.8` extend bits) and heuristic `0.55 * log2(N*M)` / HSP thresholds; graph hairball prevention remains unbenchmarked."*

---

## F2 -- C3 REFUTED: E-graph resolution uses arbitrary parsimony weights and lacks any HLA/paralog statistical model or evaluation

- **Claim:** **C3** (`docs/STRATEGIC_PLAN.md:10`) -- *"Using congruence-closure to resolve complex, nested structural paralogs (like the HLA region). This natively models evolutionary trajectories without producing unnavigable cyclical graphs."*
- **Verdict:** `REFUTED`
- **Artifacts checked:**
  - `src/egraph.rs:165-327, 381`
  - `docs/BENCHMARK_RESULTS.md:20-40`
- **Concrete flaws & evidence:**
  - `src/egraph.rs:291-308` (`extract_min_cost`) assigns arbitrary, unitless parsimony costs (`blk = 1.0`, `inv = +1.2`, `dup = +1.3`, `seq = c1 + c2 + 1.0`). There is no probabilistic evolutionary trajectory model, no branch-length or rearrangement-rate parameterization, and no MDL bit-cost.
  - `docs/BENCHMARK_RESULTS.md:20-40` contains zero quantitative evaluation on HLA or any nested structural paralog dataset.
- **Correction:** Narrow `docs/STRATEGIC_PLAN.md:10` to state that `src/egraph.rs` implements a bounded parsimony rewrite prototype (`inv`/`dup`/`seq` with fixed unitless weights `1.0/1.2/1.3`) not yet parameterized by an evolutionary model or benchmarked on HLA.

---

## F3 -- C4 REFUTED: Zero empirical evidence for the 12GB RAM exhaustion threshold or streaming throughput

- **Claim:** **C4** (`docs/STRATEGIC_PLAN.md:11`) -- *"`vg` and `odgi` crash on 12GB pangenomes due to RAM exhaustion. Our custom zero-copy GFA streaming parser (`src/gfa/mod.rs`) solves the I/O bottleneck entirely."*
- **Verdict:** `REFUTED`
- **Artifacts checked:**
  - `src/gfa/mod.rs:50-124`
  - `docs/BENCHMARK_RESULTS.md:1-40`
  - Repo-wide grep for `12GB`, `vg`, `odgi` across `docs/`, `manuscript/`, `scripts/`
- **Concrete flaws & evidence:**
  - No benchmark, memory profile, or crash reproduction for `vg`, `odgi`, or `src/gfa/mod.rs` on a 12GB pangenome exists anywhere in `docs/BENCHMARK_RESULTS.md` or `scripts/benchmarks/`.
  - `src/gfa/mod.rs:86-116` is not zero-copy: it allocates owned `String` and `Vec<u8>` heap buffers per record (`to_string()`, `as_bytes().to_vec()`), and no peak-RSS or GB/s I/O measurement is reported.
- **Correction:** Remove the unsubstantiated `12GB` / `vg` / `odgi` comparative claim and `"solves the I/O bottleneck entirely"` from `docs/STRATEGIC_PLAN.md:11` until a reproducible memory/throughput benchmark is recorded in `docs/BENCHMARK_RESULTS.md`.

---

## F4 -- C7 NARROWED: `seq` associativity induces Catalan combinatorial explosion truncated at `max_nodes = 800`

- **Claim:** **C7** (`docs/STRATEGIC_PLAN.md:27`) -- *"We will stream the pangenome linearly and only invoke the E-graph engine when we hit a complex, tangled sub-graph ('hairball'). The E-graph resolves the optimal path for that bubble, and we stitch the resolved path back into the linear stream."*
- **Verdict:** `NARROWED`
- **Artifacts checked:**
  - `src/egraph.rs:165-273, 381`
  - `src/formats.rs:66`
  - `docs/BENCHMARK_RESULTS.md:39`
- **Concrete flaws & evidence:**
  - `saturate_synteny` (`src/egraph.rs:223-253`) includes full `seq` associativity alongside `x seq x == dup(x)` (`lines 211-216`). For a bubble of $k$ blocks, the number of binary parenthesizations grows as the Catalan number $C_{k-1} = \Theta(4^k / k^{3/2})$.
  - Both production (`src/formats.rs:66`) and test (`src/egraph.rs:381`) hard-cap saturation at `saturate_synteny(&mut eg, 4, 800)`. For any bubble with $k \ge 8$ blocks ($C_7 = 429$, $C_8 = 1{,}430 > 800$, $C_9 = 4{,}862$), saturation truncates at `max_nodes = 800`, so `extract_min_cost` is not guaranteed to find the global optimum even under its own parsimony weights.
  - No bubble-size distribution on HPRC hairballs is measured (`docs/BENCHMARK_RESULTS.md:39` is unpopulated).
- **Correction:** Narrow `docs/STRATEGIC_PLAN.md:27` to specify that optimality holds only for bounded bubbles small enough not to hit `max_nodes = 800` under Catalan `seq`-associativity growth (typically $k \le 7$ blocks).

---

## F5 -- C9 REFUTED: BOCPD/MDL scores 1D/2D nucleotide sequences and has no statistical model for GFA `L` lines

- **Claim:** **C9** (`docs/PANGENOME_SUITE.md:8`) -- *"Stream an external `.gfa` file and evaluate every Link (`L` line) using the Horizon BOCPD/MDL scoring engine."*
- **Verdict:** `REFUTED`
- **Artifacts checked:**
  - `src/bocpd.rs:28-39, 58-119`
  - `src/mdl/dp.rs:81-309`
  - `src/gfa/mod.rs:21-27`
  - `src/merge.rs:246`
- **Concrete flaws & evidence:**
  - `bocpd_gc_profile` (`src/bocpd.rs:58-119`) is a 1D Beta-Bernoulli GC changepoint filter on a single `seq: &[u8]`, and `banded_mdl_gap` (`src/mdl/dp.rs:81-309`) scores pairwise alignments between two nucleotide slices `seq_x` and `seq_y`.
  - A GFA `L` line (`src/gfa/mod.rs:21-27`) represents a directed graph adjacency `(from, from_orient, to, to_orient, overlap)`, where `overlap` in pangenome graphs is typically `0M` (e.g., `src/merge.rs:246`). Neither `src/bocpd.rs` nor `src/mdl/` defines a statistical likelihood or bit-savings function for graph links, and a `0M` overlap carries 0 aligned bases (0 bits under `joint_pair_bits`).
- **Correction:** Replace `docs/PANGENOME_SUITE.md:8` with a specification of what data (e.g., flanking segment sequences or path-conditioned read support) would actually be scored, rather than claiming the existing BOCPD/MDL engine evaluates GFA `L` lines.

---

## F6 -- C10 REFUTED: `BandGrid` is a 2D linear-diagonal band undefined across graph bifurcations

- **Claim:** **C10** (`docs/PANGENOME_SUITE.md:16`) -- *"Thread the query sequence through the graph, performing banded MDL dynamic programming across graph bifurcations."*
- **Verdict:** `REFUTED`
- **Artifacts checked:**
  - `src/mdl/dp.rs:19-72`
- **Concrete flaws & evidence:**
  - `BandGrid` (`src/mdl/dp.rs:19-72`) indexes a strictly 2D linear band centered at `center(i) = (i * self.m) / self.n.max(1)` (`line 42`) with `|j - center(i)| <= band_w`.
  - Across graph bifurcations (where a node has multiple incoming/outgoing `L` edges and non-linear segment offsets), a single linear diagonal `(i * m) / n` is mathematically undefined, and `src/mdl/dp.rs` contains no DAG/graph predecessor aggregation.
- **Correction:** State in `docs/PANGENOME_SUITE.md:16` that `src/mdl/dp.rs` currently supports only linear sequence-to-sequence `BandGrid` DP and requires a graph-topological DP formulation to handle bifurcations.

---

## F7 -- Note on runner-authored claim list: `BENCHMARK_RESULTS.md` "parameter-free" claim and `IDENTITY_PRIOR`

- **Artifact checked:** `docs/BENCHMARK_RESULTS.md:5`, `src/bocpd.rs:6, 20-39`, `src/mdl/mod.rs:52-70`
- **Observation:** The claim list omits `docs/BENCHMARK_RESULTS.md:5` (*"parameter-free, Information-Theoretic Minimum Description Length (MDL) threshold defaults"*), which directly underpins **C2**. In addition to the constants in **F1**, `bocpd_gc_profile` models only 2-state GC/AT composition (`base_probs` assumes equal `G=C` and `A=T` strand parity, `src/bocpd.rs:20-25`) and `joint_pair_bits` fixes `IDENTITY_PRIOR = 0.65` (`src/bocpd.rs:6`) regardless of evolutionary divergence.

---

## Summary Verdicts

| C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| UNCHALLENGED | REFUTED | REFUTED | REFUTED | UNCHALLENGED | UNCHALLENGED | NARROWED | UNCHALLENGED | REFUTED | REFUTED | UNCHALLENGED |
