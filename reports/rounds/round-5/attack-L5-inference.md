# Round 5 Attack — Inference Lens (`inference`)

**Target:** `horizon-aligner` (`manuscript/horizon_paper.md`, `manuscript/horizon_companion.md`)  
**Method-card & Literature Note (Brief Rules 2–3):** This project keeps **no method cards** yet, and `../library` is not present in this workspace, so all checks below are grounded directly in the repository source code, benchmark artifacts, and `manuscript/` texts.

---

## F1 -- C1: `horizon_gfa stats` vs. `vg` is an asymmetric counter-only vs. in-memory graph-construction comparison

- **Claim ID:** C1  
- **Verdict:** **NARROWED**  
- **Concrete flaw:** `horizon_paper.md:10, 82` and `Figure 6A` (`manuscript/scripts/generate_figures.py:768-769`, labeling `vg stats -z -l -N` as `168× slower`) infer that `horizon_gfa stats` outperforms `vg stats -z -l -N` and `vg convert -g -p` at pangenome graph analysis. In reality, `horizon_gfa stats` (`src/bin/horizon_gfa.rs:20-56`) is a counter-only line scanner that increments five `u64` counters (`s_count`, `s_bases`, `l_count`, `p_count`, `w_count`) and discards every record without building any graph topology, adjacency index, or validating link/path consistency. By contrast, `vg stats -z -l -N` and `vg convert -g -p` (`vg v1.60.0`) parse the GFA into a full in-memory bidirected graph (`PackedGraph`/`VG`) and, for `vg convert`, serialize the converted graph.  
- **Evidence:** `src/bin/horizon_gfa.rs:28-47` (`stream_gfa` callbacks only increment 5 integers).  
- **Correction:** Narrow C1 in `manuscript/horizon_paper.md:10, 82` and `Figure 6A` to state explicitly that `horizon_gfa stats` performs a counter-only streaming pass without constructing an in-memory graph topology, whereas `vg stats` and `vg convert` construct full in-memory graph representations; drop the unqualified `168× slower` / "outperforming `vg`" framing.

---

## F2 -- C2: `48.63%` "many-to-many exact pair precision" omits 21.4% of emitted pairs and misattributes 1-to-1 over-alignment to paralogs

- **Claim ID:** C2  
- **Verdict:** **NARROWED**  
- **Concrete flaw:** `manuscript/horizon_paper.md:76` (and `docs/BENCHMARK_RESULTS.md:29`) reports `48.63%` as Horizon's "many-to-many exact pair precision (`~2.6×` repeat/paralog multiplicity)" on `GRCh38.chr6` × `CHM13.chr6:61401689–172126628`. However, `manuscript/data/chr6_hprc_comparison.json:10, 15-16, 28` shows that `0.486332` is `exact_coord_precision_on_shared = exact_coord_pairs / pairs_on_hprc_1to1_x = 107,836,669 / 221,734,653`. This denominator **excludes** `282,163,046 - 221,734,653 = 60,428,393` aligned pairs (`21.4%` of Horizon's total output) that fall outside the `1-to-1` HPRC `x` mask. Against all emitted pairs (`282,163,046`), global exact-coordinate precision is **`38.22%`**. Moreover, Horizon was run with `--resolve-dups` and still emitted `221.7M` pairs on `110.1M` bp of **1-to-1** HPRC segments (`2.01×` multiplicity on `1-to-1` loci), so attributing the low precision to benign "`~2.6×` repeat/paralog multiplicity" is an unsupported post-hoc rationalization—especially alongside `minimap2 -x asm5`, which is `5.5×–6.5×` faster (`46.95 s` vs `258.88–304.53 s`), higher in exact recall (`99.71%` vs `97.96%`), and more than double in exact precision (`99.75%`).  
- **Evidence:** `manuscript/data/chr6_hprc_comparison.json:6, 10, 15-16, 28`.  
- **Correction:** Report both masked precision (`48.63%` on `221.7M` pairs overlapping HPRC `1-to-1` `x` coordinates) and unconditional exact-coordinate precision (`38.22%` across all `282.2M` emitted pairs) in `manuscript/horizon_paper.md:76` and `Figure 7B`, and remove the unverified claim that the `2.01×` multiplicity on `1-to-1` segments represents true paralogs.

---

## F3 -- C3: `horizon_guide --progressive --hal --rbest` writes a 4-segment contig-0 star-tree stub, not a verified 5-mammal progressive HAL archive

- **Claim ID:** C3  
- **Verdict:** **REFUTED**  
- **Concrete flaw:** C3 infers that `horizon_guide --progressive --hal --rbest` aligned 5 mammals (`1,109` MAF blocks, `967` multi-species) into a verified progressive `.hal` archive (`4` HAL bottom segments verified via `horizon_wga --hal-info`). In `src/bin/horizon_guide/main.rs:94-103, 251-296`, `align_species_pair` only records a HAL block when `xi == 0 && yi == 0` (the first contig of each species) and collapses that entire contig pair into a single bounding-box span `(x0, y0.min(y1), x1 - x0 + 1, y1 < y0)`. All remaining contigs and all blocks added by `run_progressive_sister_rescue` / `run_ancestral_consensus_rescue` (`main.rs:271-289`) are omitted from `hal_all_blocks`. `HalAlignment::from_star_blocks` (`src/hal/mod.rs:90-170`) then constructs a **star tree** rooted at `hinge_sp` (not a progressive guide tree) with hardcoded `.chr1` sequence names (`src/hal/mod.rs:144, 161`), yielding `4` bottom segments for `4` query species—not `1,109` blocks. Furthermore, `horizon_wga --hal-info` (`inspect_hal_genomes`, `src/hal/hdf5_read.rs:806`) only reads the HDF5 superblock and root-group member names, not alignment coordinates. Finally, neither `horizon_merge` (`92` `chrM` species) nor `horizon_guide` (`5` mammals) reports any alignment accuracy (precision/recall) against a reference.  
- **Evidence:** `src/bin/horizon_guide/main.rs:94-103, 251-296`; `src/hal/mod.rs:90-170`; `src/hal/hdf5_read.rs:806`.  
- **Correction:** Remove the claim in `manuscript/horizon_paper.md:10, 53, 82` that `horizon_guide --progressive --hal` builds and verifies a progressive 5-mammal `.hal` archive (it exports a contig-0 star-tree bounding-box stub of 4 segments), and qualify `horizon_merge` / `horizon_guide` metrics as runtime/block-count measurements without accuracy validation.

---

## F4 -- C4: "Held-out" MAF-slice accuracy and `18%`-GC null claims suffer from oracle pairing, iterative reuse, and cherry-picked GC nulls

- **Claim ID:** C4  
- **Verdict:** **NARROWED**  
- **Concrete flaw:**  
  1. **Prior-round context (Rule 4):** Round 1 (`F3`) and Round 3 (`F2`) already upheld that `extract_fasta_from_maf.py` filters FASTAs to homologous MAF blocks and `horizon_eval` evaluates oracle `(ref_id, tgt_id)` contig pairs. While `manuscript/horizon_paper.md:59` now discloses `(ref_id, tgt_id)` pair evaluation, C4 (`horizon_paper.md:10, 59`) still presents these numbers as unqualified `horizon_wga` accuracy and still labels the nuclear slices **"held-out"** (`horizon_paper.md:59`, `docs/BENCHMARK_RESULTS.md:42`) even though `Figure 6B` explicitly shows these exact datasets (horse, macaque, marmoset) were re-evaluated to measure precision gains from adding the `R1-HMM` BOCPD LCB gate.  
  2. **Cherry-picked random-DNA null:** C4 claims `horizon_wga` eliminates all false-positive pairs on `18%`-GC random DNA (`0` pairs in `0.0069 s`). However, `docs/BENCHMARK_RESULTS.md:13` and `horizon_paper.md:59` show that on the same `10,000`-bp random-DNA test, `horizon_wga` still emits **`288` false-positive pairs** at `30%`-GC (`0.4150 s`) and **`149` false-positive pairs** at `45%`-GC (`0.3829 s`). Testing a single `10 kb` pair at `18%`-GC (`n = 1`) does not license general claims of false-positive elimination across AT-rich or general genomes.  
- **Evidence:** `src/bin/horizon_eval.rs:187-205`; `scripts/benchmarks/run_all_evals.py:39, 65-75`; `docs/BENCHMARK_RESULTS.md:13, 42`; `manuscript/horizon_paper.md:59, 86` (`Figure 6B`).  
- **Correction:** Drop the "held-out" label on re-evaluated nuclear slices in `manuscript/horizon_paper.md:59` and `Table 1`, qualify the Abstract (`horizon_paper.md:10`) to note MAF-extracted `(ref_id, tgt_id)` slice evaluation, and report the `30%`-GC (`288` pairs) and `45%`-GC (`149` pairs) false positives alongside the `18%`-GC (`0` pairs) result.

---

## F5 -- C5: `BandGrid` index formula holds algebraically, but `W` clamping and diagonal centering limit asymmetric gaps

- **Claim ID:** C5  
- **Verdict:** **SURVIVED**  
- **Concrete flaw (non-refuting):** The indexing formula `idx(i, j) = i · (2W + 1) + (j − i + W)` matches `BandLayout::idx` in `src/mdl/dp.rs:164-169`. Two secondary companion-text caveats apply: (1) `BandLayout::new` (`src/mdl/dp.rs:142`) clamps `bw = band_w.min(n.max(m))`, so the worked 6-bp example in `manuscript/horizon_companion.md:179` inaccurately states `(6 + 1) × (2 · 180 + 1)` floats instead of `7 × 13`; and (2) centering the band on `j = i` rather than `j = i · (M / N)` makes the endpoint `(N, M)` unreachable whenever an inter-anchor gap has `|N − M| > W`.  
- **Evidence:** `src/mdl/dp.rs:141-169`; `manuscript/horizon_companion.md:139-145, 179`.  
- **Correction:** C5 survives as stated; optionally correct `manuscript/horizon_companion.md:179` (`(6 + 1) × (2 · 6 + 1)`).

---

## F6 -- C6: `# EG:Z:` tags omit nested inversions and cannot distinguish collapsed paralogs from single-copy synteny

- **Claim ID:** C6  
- **Verdict:** **REFUTED**  
- **Concrete flaw:** `manuscript/horizon_companion.md:52` claims that collapsed tandem duplications *and* nested inversions resolved by `src/egraph.rs` are recorded as `# EG:Z:` tags so downstream tools can distinguish collapsed paralogs from single-copy synteny. In fact:  
  1. `egraph::collect_dups` (`src/egraph.rs:518-538`) only collects leaves under `Expr::Dup` and **ignores `Expr::Inv`**, so nested inversions are never recorded in `# EG:Z:` tags.  
  2. `format_gfa_paralog_tags` (`src/formats.rs:153-158`) and `src/merge.rs:361` emit `# EG:Z:` as a **GFA comment line** (`# EG:Z:<ref>__<qry>\tB0,B3`) containing only coarse `bin_bp`-bucketed indices (`blk[0].0 / bin_bp`, `src/formats.rs:142`)—not segment IDs, walk coordinates, or per-record tags. Standard GFA parsers ignore `#` comment lines, and downstream tools cannot map `B0` back to specific `S` segments or `P`/`W` walks to separate collapsed paralogs from single-copy synteny. (Meanwhile, `horizon_resolve --tag-walks` in `src/bin/horizon_resolve.rs:476` writes `EG:Z:` as an optional tag on `W` lines, not `# EG:Z:` comments.)  
- **Evidence:** `src/egraph.rs:518-538`; `src/formats.rs:138-158`; `src/merge.rs:361`; `src/bin/horizon_resolve.rs:476`.  
- **Correction:** Replace C6 in `manuscript/horizon_companion.md:52` to state that `# EG:Z:` comment lines record only `bin_bp`-bucketed `Expr::Dup` block labels (omitting inversions and per-segment/walk mappings), and do not enable downstream per-segment separation of paralogs from single-copy synteny.

---

## F7 -- Claim-list steering note

- **Observation:** The runner-authored claim list steers attention toward a trivial algebraic index identity (C5) while omitting higher-risk inference claims in `manuscript/horizon_paper.md` and `manuscript/horizon_companion.md`: (1) Conformal Risk Control calibration of Horizon Bubbles (`α = 0.07`, `n_cal = 15`, hardcoded `λ̂ = 10.0` bits and `KARLIN_MDL_COEFF = 0.55`); (2) `egraph.rs` "minimum-description-length ancestral order" (`horizon_paper.md:44`) despite having no phylogenetic tree; and (3) `horizon_wga --graph` (`horizon_paper.md:82`) evaluated on only three short queries (`400`-bp all-`N`, `258`-bp, `356`-bp) with no precision/recall ground truth.

---

## Claim Verdict Summary

| Claim ID | Verdict |
| :--- | :--- |
| C1 | NARROWED |
| C2 | NARROWED |
| C3 | REFUTED |
| C4 | NARROWED |
| C5 | SURVIVED |
| C6 | REFUTED |
