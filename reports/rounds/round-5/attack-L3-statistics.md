# Round 5 Attack -- Statistics Lens

**Scope and rules disclosures (Rules 2–4):**
- **Literature (Rule 2):** `../library` is not present in this workspace, so all prior-work checks were performed against `refs/catalog_index.json` and `manuscript/references.bib` alone rather than from memory of held paper texts.
- **Method cards (Rule 3):** This project keeps **no method cards** yet (`BRIEF-statistics.md:142`). Statistical claims are checked against the codebase, benchmark data/scripts, and cited references.
- **Prior adjudications (Rule 4):** Checked `scripts/round history` for `n_cal`, `BOCPD`, `0.55`, and `precision`. Prior rulings on `r3 adjudication-L5.md:13 [F2]` (oracle `(ref_id, tgt_id)` pairing in `horizon_eval`), `r4 adjudication-L2.md:15 [F1]` (`n = 1` chromosome/run on `chr6`), and `r4 adjudication-L3.md:62 [F7]` (`n_cal = 15` CRC and `KARLIN_MDL_COEFF = 0.55`) are explicitly cited where they bear on Round 5.

---

## F1 -- C1: Single-run (`n = 1`) runtime comparison with unquantified variance and confounded task scope

- **Claim:** `C1` (`manuscript/horizon_paper.md:82`, `docs/BENCHMARK_RESULTS.md:25,76`)
- **Verdict:** `NARROWED`
- **Concrete flaw:** `C1` reports `horizon_gfa stats` at `0.08 s` (`22.6 MB` RSS) versus `vg stats -z -l -N` (`4.2984 s`, `77.0 MB` RSS) and `vg convert -g -p` (`1.8292 s`, `76.6 MB` RSS) from `n = 1` run per tool, with no replicate count, no mean/median/IQR/SD, no cold- vs. warm-cache control (critical for a `0.08 s` process), and asymmetric decimal precision (`0.08 s` to 2 decimals vs. `4.2984 s` / `1.8292 s` to 4 decimals). Moreover, `horizon_gfa stats` is a single-pass counter-only line scanner that discards records after counting (`ATTACKED-horizon_paper.md:47`), whereas `vg convert -g -p` constructs and serializes a `PackedGraph` and `vg stats -z -l -N` builds an in-memory graph representation; comparing wall times without equalizing the computational task confounds algorithmic efficiency with work avoided.
- **Evidence:** `ATTACKED-horizon_paper.md:47,82`, `docs/BENCHMARK_RESULTS.md:25,76` (single point estimates, no replicate statistics).
- **Correction:** Qualify `C1` (`manuscript/horizon_paper.md:82`) as a single-run (`n = 1`) measurement comparing streaming line-counting (`horizon_gfa stats`) against full graph construction/conversion (`vg stats` / `vg convert`), or report median and IQR across `n >= 5` cold/warm runs on matched subgraphs.

## F2 -- C2: Conditioned precision denominator (`48.63%`) omits `60.4M` (`21.4%`) emitted pairs, inconsistent with `~2.6x` multiplicity and `n = 1` locus

- **Claim:** `C2` (`manuscript/horizon_paper.md:10,76,90`, `docs/BENCHMARK_RESULTS.md:29,74`)
- **Verdict:** `REFUTED`
- **Concrete flaw:**
  1. **Undisclosed conditioned denominator:** `C2` reports `"48.63% many-to-many exact pair precision"`, which is computed in `manuscript/data/chr6_hprc_comparison.json:10,15,16,28` as `exact_coord_pairs / pairs_on_hprc_1to1_x = 107,836,669 / 221,734,653 = 0.486332`. That denominator **excludes `282,163,046 - 221,734,653 = 60,428,393` pairs (`21.42%` of Horizon's `282,163,046` emitted pairs)** whose `GRCh38.chr6` coordinate falls outside the `110,082,990` bp HPRC `1-to-1` shared segments (e.g., on the `GRCh38.chr6` p-arm `0-61 Mbp` or non-shared repeats). Across all emitted pairs (`horizon_total_pairs = 282,163,046`), Horizon's **unconditional** exact pair precision is `107,836,669 / 282,163,046 = 38.22%` — **10.41 percentage points lower** than `48.63%`. By contrast, `minimap2 -x asm5` (`manuscript/data/chr6_minimap2_summary.json:3,6,7,15`) emits `110,204,880` total pairs (`110,033,678` on HPRC `1-to-1` `x`), so its conditioned (`99.75%`) and unconditional (`99.60%`) exact pair precisions differ by only `0.15` percentage points.
  2. **Inconsistent multiplicity vs. precision denominator:** `C2` pairs `48.63%` with `(~2.6x repeat/paralog multiplicity)`. Using `horizon_total_pairs / horizon_and_hprc_x_bp = 282,163,046 / 109,575,507 = 2.575x` (`~2.6x`) includes all `282.2M` pairs in the numerator, whereas `48.63%` restricts the denominator to `221,734,653` pairs (`221,734,653 / 110,082,990 = 2.014x`).
  3. **Prior `n = 1` adjudication unresolved:** As already ruled `UPHELD` in `r4 adjudication-L2.md:15 [F1]`, `C2` still generalizes from `n = 1` chromosome (`chr6`), `n = 1` q-arm target (`CHM13.chr6:61401689-172126628`), and `n = 1` run (`282,163,046` pairs in `chr6_hprc_comparison.json:10` vs. `282,193,564` in `ATTACKED-horizon_paper.md:76`).
- **Evidence:** `manuscript/data/chr6_hprc_comparison.json:5-28`, `manuscript/data/chr6_minimap2_summary.json:3-15`, `ATTACKED-horizon_paper.md:76,90`.
- **Correction:** Replace `"48.63% many-to-many exact pair precision (~2.6x repeat/paralog multiplicity)"` in `manuscript/horizon_paper.md:76` (and Figure 7B) with `"38.22% unconditional exact pair precision across all 282.16M emitted pairs (48.63% conditioned on the 221.73M pairs within HPRC 1-to-1 GRCh38.chr6 segments, 2.01x local multiplicity; 2.58x total emitted multiplicity)"`, and qualify as `n = 1` chromosome/run.

## F3 -- C3: Single-run (`n = 1`) multi-species timings with zero accuracy or concordance statistics

- **Claim:** `C3` (`manuscript/horizon_paper.md:10,82`, `docs/BENCHMARK_RESULTS.md:10-11,57-58`)
- **Verdict:** `NARROWED`
- **Concrete flaw:** `C3` reports `n = 1` wall-clock runtime and block counts for `horizon_merge` (`21.78 s`, `26.6 MB` RSS, `460` multi-species blocks on `92` `chrM` species) and `horizon_guide` (`0.859 s`, `1,109` MAF blocks, `967` multi-species, `4` HAL bottom segments) with no replicate variance/CI and **zero alignment accuracy statistics** (no precision, recall, column-agreement, or tree-consistency comparison against the UCSC `100way` reference MAF or Cactus/MULTIZ). Block count (`460` or `1,109`) and `4` HAL bottom segments measure output volume, not statistical accuracy.
- **Evidence:** `ATTACKED-horizon_paper.md:82`, `docs/BENCHMARK_RESULTS.md:10-11,57-58`.
- **Correction:** Qualify `C3` (`manuscript/horizon_paper.md:82`) as `n = 1` throughput/memory measurements without alignment-accuracy validation against reference multi-species MAFs.

## F4 -- C4: Oracle `(ref_id, tgt_id)` pairing and MAF-only FASTA extraction inflate precision; `18%`-GC null claim rests on `n = 1` synthetic draw

- **Claim:** `C4` (`manuscript/horizon_paper.md:10,59-70`, `docs/BENCHMARK_RESULTS.md:13-22,36-51`)
- **Verdict:** `NARROWED`
- **Concrete flaw:**
  1. **Conditional precision on pre-extracted homologous blocks under oracle pairing:** As `r3 adjudication-L5.md:13 [F2]` already upheld, `src/bin/horizon_eval.rs:165-197` restricts evaluation strictly to `(ref_id, tgt_id)` keys present in `maf_by_pair`, never aligning non-homologous contig pairs (e.g., among the `31` `Canis_lupus_familiaris` contigs in `cactus241way chr22 1MB`). Moreover, `scripts/benchmarks/extract_fasta_from_maf.py:40-45` constructs `ref_fa` and `tgt_fa` by concatenating *only* the aligned MAF segments and stripping all unaligned intergenic/intronic background. While `ATTACKED-horizon_paper.md:59` now adds a note that `horizon_eval` evaluates `(ref_id, tgt_id)` pairs individually, `C4` (Abstract line 10) and Table 1 still present these conditional numbers (`97.4%`, `87.0%`, `95.1%`, `78.5%`) as unconditional `horizon_wga` precision.
  2. **`n = 1` random-DNA trial cannot establish `"eliminating all false-positive pairs"`:** `C4` claims `"eliminating all false-positive pairs on 18%-GC random DNA (0 pairs in 0.0069 s)"` from a single (`n = 1`) `10,000`-bp synthetic sequence. At `n = 1` trial, the 95% Clopper–Pearson upper bound on trial-level false-positive probability is `1 - 0.05^(1/1) = 95%`, and `30%`-GC and `45%`-GC single draws in the same sentence (`ATTACKED-horizon_paper.md:59`) emit `288` and `149` false-positive pairs respectively.
- **Evidence:** `src/bin/horizon_eval.rs:165-197`, `scripts/benchmarks/extract_fasta_from_maf.py:40-45`, `ATTACKED-horizon_paper.md:10,59-70`.
- **Correction:** Qualify `C4` (`manuscript/horizon_paper.md:10`, Table 1) as conditional precision on MAF-extracted homologous segments under oracle `(ref_id, tgt_id)` pairing, and replace `"eliminating all false-positive pairs on 18%-GC random DNA"` with `"emitting 0 pairs on a single 10,000-bp 18%-GC random DNA draw (vs. 288 and 149 pairs on single 30%-GC and 45%-GC draws)"`.

## F5 -- C5 and C6: Deterministic indexing and tag-emission claims

- **Claims:** `C5` (`manuscript/horizon_companion.md:139-143`), `C6` (`manuscript/horizon_companion.md:52`)
- **Verdict:** `UNCHALLENGED` (both)
- **Rationale:** Neither `C5` (flat `BandGrid` index formula `idx(i, j) = i * (2W + 1) + (j - i + W)`) nor `C6` (recording `# EG:Z:` tags from `src/egraph.rs` in GFA/MAF outputs, verified in `src/merge.rs:361`, `src/formats.rs:157`, `src/bin/horizon_resolve.rs:476`) makes a statistical claim. Under the statistics lens, both are `UNCHALLENGED`.

## F6 -- Claim-list steering: `C1`–`C6` omit uncalibrated `n_cal = 15` CRC, `0.55` Karlin coefficient, and companion toy-example discrepancy

- **Scope:** Meta-finding on `round.json` claim list vs. `ATTACKED-horizon_paper.md:40`, `ATTACKED-horizon_companion.md:134,153,181`
- **Concrete flaw:** The claim list `C1`–`C6` omits the manuscript's primary statistical calibration claims:
  1. **Conformal Risk Control (`n_cal = 15`, `α = 0.07`, `λ̂ = 10.0` bits) and `KARLIN_MDL_COEFF = 0.55`:** As adjudicated `NARROWED` in `r4 adjudication-L3.md:62 [F7]`, `n_cal = 15` only achieves `1 / (n_cal + 1) = 0.0625 <= 0.07` if empirical risk `R̂_15 = 0` (zero margin for calibration error), and in code `min_delta_bits = 10.0` (`src/mdl/mod.rs:71`) and `KARLIN_MDL_COEFF = 0.55` (`src/mdl/mod.rs:26`) remain hardcoded constants with no calibration script or dataset in the repository.
  2. **Worked toy-example bit-savings mismatch:** In `ATTACKED-horizon_companion.md:173-181`, for the `6 bp` zero-match gap (`GGGGCC` vs. `CCTTTT`) under `HmmParams::default()` (`identity_prior = 0.65`, `gap_open_bits = -4.5`, `mm_bits = log2(1 - 2 * 2^-4.5) = -0.1337`), `joint_pair_bits` (`src/bocpd.rs:41-52`) yields `6 * log2(1 - 0.65) + 5 * mm_bits = -9.09 - 0.67 = -9.75` bits (not `ΔL ≈ -6.4` bits as printed at line 181), and `bocpd_gc_profile` (`src/bocpd.rs:87-175`) gives `P(W) ≈ 0.684` at `t = 11` (mean `0.698` over `0..12`, vs. `P(W) ≈ 0.72` at line 177).
- **Evidence:** `src/mdl/mod.rs:26,71,140-148`, `src/bocpd.rs:41-52,87-175`, `ATTACKED-horizon_paper.md:40`, `ATTACKED-horizon_companion.md:134,153,177,181`.
- **Correction:** Include CRC calibration provenance (`n_cal = 15`, `0.55`) in the round's claim scope and recompute the worked toy example in `manuscript/horizon_companion.md:177-181` from the actual Rust `bocpd_gc_profile` / `banded_mdl_gap` output.

---

## Summary Table

| Claim ID | Verdict | Summary |
| :--- | :--- | :--- |
| `C1` | `NARROWED` | `n = 1` run, no variance/CI, asymmetric precision, confounded counter vs. graph-build task (`F1`) |
| `C2` | `REFUTED` | `48.63%` excludes `60.4M` (`21.4%`) emitted pairs; unconditional precision is `38.22%`; inconsistent `~2.6x` multiplicity (`F2`) |
| `C3` | `NARROWED` | `n = 1` run with zero alignment accuracy/concordance metrics (`F3`) |
| `C4` | `NARROWED` | Oracle `(ref_id, tgt_id)` pairing + MAF-only FASTA extraction inflate precision; `18%`-GC claim rests on `n = 1` draw (`F4`) |
| `C5` | `UNCHALLENGED` | Deterministic `BandGrid` indexing formula; no statistical claim (`F5`) |
| `C6` | `UNCHALLENGED` | Deterministic `# EG:Z:` tag emission; no statistical claim (`F5`) |
