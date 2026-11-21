# Horizon Suite Updates & Additions for the VGP All-to-All Alignment Project

> **Target context:** [`pangenome/vgp-align-allocation/alignment-2026/main.md`](https://github.com/pangenome/vgp-align-allocation/blob/main/alignment-2026/main.md)
> **Horizon status (`42700e1` / `9dfec4a`):** **P0–P4 complete & benchmarked on real VGP/GenomeArk + HPRC data** (`120` passing active unit tests: `86` lib + `34` bin, `1` ignored; `11` suite CLI binaries + `bench_chr1`).

## 1. VGP Scale, Bottlenecks, and Where Horizon Fits

The VGP 2026 allocation scales an all-to-all vertebrate alignment from **581 assemblies (336,980 ordered pairs)** to **~1,381 assemblies (1,568,800 new pairs)** in the next tranche, on a path toward **4,467 Phase 1 assemblies (19.95M pairs)** spanning ~500 Myr of vertebrate divergence.

Key operational constraints from `main.md`:
* **Format currency is PAF (`.paf` / compressed PAF):** ~1.5 TB compressed at 581 assemblies, projecting to ~8.5 TB at 1,381 assemblies, hosted on GenomeArk (AWS Open Data).
* **High-divergence false positives vs. cost:** `FastGA` (0.263 node-h/pair) is **2.4× cheaper** than `wfmash` (0.627 node-h/pair), but produces high false-positive rates across ~500 Myr splits and currently requires `SweepGA` + CMA-ES calibration against expensive `wfmash` runs.
* **Implicit pangenome via `impg`:** Whole-genome graphs are not materialized globally; instead, `impg query` and `impg similarity` project 25 kb ILS windows and BUSCO gene intervals transitively through pairwise PAFs into local homologous sequence sets.
* **Batch execution model:** 48-hour wall-clock chunks under `pylauncher`/`ParaFly` with `TIMEOUT` checkpoint-and-resume, single-node memory bounds, and bounded scratch space.

Horizon's core strengths — **MDL/BOCPD statistical abstention (Horizon Bubbles)**, **6-frame translated `Aa5Index` + `MASK_RY20` sensitivity**, **E-graph structural paralog resolution**, and **streaming GFA/coord tooling** — map directly onto VGP's high-divergence filtering and regional `impg` locus analysis bottlenecks.

---

## 2. Immediate Utilization of Existing Horizon Tools

These workflows work today without new binary development:

1. **High-divergence validation & rescue on difficult pairs (`horizon_wga`, `horizon_eval`):**
   * Use [`src/bin/horizon_wga.rs`](../src/bin/horizon_wga.rs) (with `MASK_RY20` transition-tolerant seeds, `Aa5Index` translated seeds, and `estimate_hmm_params`) on bounded vertebrate validation subsets (e.g., mammal–teleost–lamprey or outgroup pairs) to audit where `FastGA` drops sensitivity or emits spurious hits.
   * Use [`src/bin/horizon_eval.rs`](../src/bin/horizon_eval.rs) to measure precision/recall/F1 and Horizon Bubble abstention rates against `wfmash` truth sets.
2. **Regional `impg` locus de-tangling & structural resolution (`horizon_clean`, `horizon_resolve`):**
   * Once `impg query` extracts homologous sequences for a 25 kb window or BUSCO locus and a local GFA is built (e.g., via `pggb`/`minigraph`), pipe it through [`horizon_clean`](../src/bin/horizon_clean.rs) to sever unsupported or compositionally discordant (`BOCPD + MDL`) spurious junctions into `HB:i:1` Horizon Bubbles.
   * Run [`horizon_resolve --collapse-dups --tag-walks`](../src/bin/horizon_resolve.rs) to simplify nested inversions and tandem duplications via MDL-calibrated E-graphs (`MdlCostModel`) before computing locus distance matrices or ILS trees.
3. **BUSCO & annotation liftover across assemblies (`horizon_lift`):**
   * Use [`src/bin/horizon_lift/`](../src/bin/horizon_lift/main.rs) on regional or chromosome-scale GFAs to project BED/GFF3 BUSCO coordinates with `context_jaccard` paralog disambiguation and `EXACT`/`PARTIAL`/`SPLIT_k/N`/`DUP` status reporting.
4. **Subtree MSA & optional HAL export (`horizon_merge`, `horizon_guide`):**
   * Use [`horizon_merge`](../src/bin/horizon_merge.rs) to thread pairwise MAFs onto a focal reference with exact insertion-column preservation, or [`horizon_guide`](../src/bin/horizon_guide/main.rs) for progressive star/guide-tree MSAs (and `.hal` export when Cactus-compatible downstream tools are used).

---

## 3. Recommended Updates & Additions to Horizon

### A. Native PAF Support & MDL Alignment Filtering (`horizon_filter` / PAF I/O) — *Highest Impact*
* **Gap:** Horizon currently reads/writes MAF, GFA, and HAL, whereas VGP's entire production pipeline and 1.5–8.5 TB GenomeArk archive are **PAF** (with `cg:Z:` CIGAR strings).
* **Addition:**
  1. Add a streaming `src/paf.rs` module (supporting `.paf`, `.paf.gz` via `MultiGzDecoder`, and optional `.zst` decompression) that parses standard 12-column PAF + `cg:Z:` CIGAR and `AS:i:`/`dv:f:` tags.
  2. Add a **`horizon_filter`** binary (or `horizon_wga --filter-paf`) that scores raw `FastGA` PAF records using Horizon's **BOCPD + Pair-HMM MDL bit-savings** (`joint_pair_bits`, `KARLIN_MDL_COEFF`, and `sdust_mask` repeat-aware gating) and prunes alignments with negative bit-savings.
  3. Emit filtered PAF + optional `HB:i:` / `md:f:` (MDL bits) tags on `stdout` so `horizon_filter` drops straight into VGP's `pylauncher`/`ParaFly` pipeline without needing CMA-ES/`wfmash` calibration per clade.

### B. Regional `impg` Locus Pipeline Mode (`horizon_wga --locus` / `impg` Bridge)
* **Gap:** VGP's `O3` objective runs `impg query` across thousands of 25 kb windows and BUSCO genes, then computes `impg similarity` / distance matrices. At ~500 Myr divergence, raw `FastGA` alignments inside a locus often mix paralogs, fragmented exons, and false-positive tails.
* **Addition:**
  * Add a lightweight **`horizon_locus`** (or `horizon_wga --paf-locus`) workflow taking a multi-FASTA from `impg query` (+ optional guide PAF), running SDUST-masked + `Aa5Index`/`MASK_RY20` anchored MDL alignment, resolving paralogs via `egraph::MdlCostModel`, and emitting:
    * Clean MAF / GFA for the locus,
    * Pairwise divergence/coverage matrix (`TSV`) for ILS and BUSCO fragmentation analysis, and
    * Explicit Horizon Bubble abstention coordinates where homology is statistically unsupported.

### C. Wire Translated `Aa5Index` Seeding into DAG Alignment & Liftover
* **Gap:** [`src/seeding/translated.rs`](../src/seeding/translated.rs) (`Aa5Index`) is implemented for pairwise contig alignment, while [`src/mdl/dag.rs`](../src/mdl/dag.rs) (`align_seq_to_dag_seeded`) and `horizon_lift` currently rely only on nucleotide `MASK_W10` + `MASK_RY20` seeds.
* **Addition:** Enable optional 6-frame `Aa5Index` seeding in `horizon_wga --graph` so divergent vertebrate coding exons (e.g., tetrapod–teleost–agnathan BUSCOs) anchor onto pangenome DAGs even when nucleotide identity drops below `MASK_RY20` thresholds.

### D. `horizon_clean` & `horizon_resolve` Hardening for Multi-Hundred-Species Vertebrate Graphs
* **Per-bubble divergence estimation in `horizon_resolve`:** Currently `horizon_resolve` uses a single global `--divergence` CLI parameter (`0.02`). A single vertebrate locus spanning primates, birds, and lamprey has order-of-magnitude variation in branch length; `MdlCostModel` should estimate `d` per snarl/window (or from path/segment divergence tags) rather than using a fixed global constant.
* **Snarl decomposition & `u32` `PackedGraph`:** Replace 1D walk windowing in `horizon_resolve` with linear-time biedged snarl/bubble decomposition and `u32`-interned CSR adjacency so regional graphs with 800–1,400 haplotypes scale without walk-order bias.
* **Link-overlap CIGAR accounting in `horizon_clean`:** When severing links and splitting `W` walk records in [`src/bin/horizon_clean.rs`](../src/bin/horizon_clean.rs), subtract link overlap bases (parsed from `L`-line CIGARs) when computing split `W` start/end coordinates on non-`0M` graphs.

### E. HPC / Stampede3 Execution Ergonomics
* **Bounded RSS & streaming:** Support `.zst` input/output streaming (matching VGP's compressed PAF storage), deterministic `--threads` control (Rayon pool sizing per `ParaFly` worker), and single-pass `stdin`/`stdout` emission for `horizon_filter` so scratch usage stays flat under 48-hour `TIMEOUT` checkpoint-and-resume cycles.

---

## 4. Implementation Status (P0–P4 Complete at `9dfec4a`)

| Priority | Item | Status | Real-Data Verification |
| :--- | :--- | :--- | :--- |
| **P0** | `src/paf.rs` + `horizon_filter` (SDUST+MDL PAF filter & `--net`) | **Implemented** | `chr6` FastGA (`56,857` recs): `0.45s–0.48s` seq-free (`565` recs / `110.00 Mbp`) / `4.86s–4.89s` SDUST+MDL (`553` recs; 2D HPRC re-scoring pending `numpy`); GenomeArk Koala × Amphioxus (`36,601` recs, pre-`9dfec4a` `/tmp/vgp_test/`): `0.10s` (`1,976` recs, `71.6%` overlap with `FastGA_11`) and rejects `FastGA_11`'s #1 longest hit (`3,548 bp` repeat, `md:f:-210.51 hb:i:1`) |
| **P1** | `horizon_clean` overlap-CIGAR `W` split fix + `horizon_resolve` topological window divergence | **Implemented** | Unit-tested & verified across synthetic and HPRC `10k` subgraphs |
| **P2** | `aa5_keys` 6-frame translated fallback in `horizon_wga --graph` (`src/mdl/dag_seed.rs`) | **Implemented** | Unit-tested & active alongside `MASK_W10` + `MASK_RY20` |
| **P3** | `horizon_locus` (`impg query` multi-FASTA -> `.similarity.tsv`, `.maf`, `.gfa`) | **Implemented** | Avian 2-way (`386.7 kb` × `349.4 kb`): `0.83s–0.86s` (`236.7–239.3 MB` RSS); Human 4-way 25 kb: `0.41s–0.42s` (`171.8 MB` RSS, `102` MAF / `92` multi-sp); Koala × Amphioxus 25 kb: `0.26s` (`40.1 MB` RSS, `9,697/25,000` bp) |
| **P4** | `u32` `PackedGraph` + $O(V+E)$ Tarjan biedged snarl decomposition in `horizon_resolve` + `.zst` streaming (`src/io_util.rs`) | **Implemented** | `chr6` `.paf.zst` → `.paf.zst` streaming filter in `0.50s`; HPRC `10k` `horizon_resolve` in `5.571s` (`23.6 MB` RSS) |
