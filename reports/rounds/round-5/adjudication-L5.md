# Round 5 ADJUDICATION — Inference Lens (`inference`)

**Attack reviewed:** `reports/rounds/round-5/attack-inference.md` (sealed digest `af4640c3c4a3`, verified via `sha256sum`).  
**Prior history checked:** `scripts/round history "horizon_eval" --kind adjudication` (Round 1 F3 / Round 3 F2 upheld MAF-slice filtering and `(ref_id, tgt_id)` pair evaluation).  
**Method cards:** None present in this repository (`refs/methods/` absent).

---

## Finding-by-Finding Rulings

### F1 — C1 (`horizon_gfa stats` vs. `vg` comparison)
- **Ruling:** **UPHELD** (attacker's **NARROWED** verdict on C1 holds).
- **What survives:** `src/bin/horizon_gfa.rs:20-56` (`cmd_stats`) only streams GFA lines to increment five `u64` counters (`s_count`, `s_bases`, `l_count`, `p_count`, `w_count`) without constructing any graph topology, adjacency index, or link/path validation, whereas `vg stats -z -l -N` and `vg convert -g -p` parse the GFA into a full in-memory graph representation. While `manuscript/horizon_paper.md:47, 74` notes "counter-only", `manuscript/horizon_paper.md:10, 82` and `manuscript/scripts/generate_figures.py:769` (`168× slower`) frame the comparison as an unqualified performance win on pangenome graph analysis.
- **Repair:** Endorsed (`manuscript/horizon_paper.md:10, 82` and `Figure 6A` label); **the fix is free** (text/label qualification only, no re-measurement).

### F2 — C2 (`48.63%` "many-to-many exact pair precision" on `GRCh38.chr6` × `CHM13.chr6`)
- **Ruling:** **UPHELD** (attacker's **NARROWED** verdict on C2 holds).
- **What survives:** `manuscript/data/chr6_hprc_comparison.json:10, 15-16, 28` confirms `0.486332` is `exact_coord_pairs / pairs_on_hprc_1to1_x = 107,836,669 / 221,734,653`, excluding `60,428,393` pairs (`21.4%` of all `282,163,046` emitted pairs) outside the `1-to-1` HPRC `x` mask; unconditional exact-coordinate precision across all emitted pairs is `38.22%` (`107,836,669 / 282,163,046`). Moreover, Horizon emits `221.7M` pairs on `110.1M` bp of `1-to-1` HPRC segments (`2.01×` multiplicity on `1-to-1` loci), so attributing the low precision to benign "`~2.6×` repeat/paralog multiplicity" (`282.2M / 110.1M`) is unsupported.
- **Repair:** Endorsed (report both masked `48.63%` and unconditional `38.22%` precision in `manuscript/horizon_paper.md:76` and `Figure 7B`, and drop the unverified paralog-multiplicity attribution); **the fix is free**.

### F3 — C3 (`horizon_guide --progressive --hal --rbest` and `horizon_merge`)
- **Ruling:** **NARROWED** (upheld on `horizon_guide --progressive --hal`; narrowed relative to the attacker's blanket **REFUTED** on C3 because `horizon_merge`'s runtime/block-count numbers stand once qualified).
- **What survives:** In `src/bin/horizon_guide/main.rs:94-103, 251-296`, `align_species_pair` only records a HAL block for `xi == 0 && yi == 0` as a single bounding box, omits all other contigs and all `run_progressive_sister_rescue` / `run_ancestral_consensus_rescue` blocks, and passes them to `HalAlignment::from_star_blocks` (`src/hal/mod.rs:90-171`), which builds a star tree with hardcoded `.chr1` sequence names (`4` bottom segments for `4` query species, not `1,109` progressive blocks). `inspect_hal_genomes` (`src/hal/hdf5_read.rs:806-817`, `--hal-info`) only reads the HDF5 superblock and v0 root member names. However, C3's `horizon_merge` `92`-species `chrM` runtime/block-count measurement (`21.78 s`, `26.6 MB` RSS, `460` blocks) is not refuted, only unvalidated for alignment accuracy.
- **Repair:** Endorsed for `horizon_guide --progressive --hal` and accuracy qualification (`manuscript/horizon_paper.md:10, 53, 82`); **the fix is free**.

### F4 — C4 ("Held-out" MAF-slice accuracy and `18%`-GC null)
- **Ruling:** **UPHELD** (attacker's **NARROWED** verdict on C4 holds).
- **What survives:** Prior rounds (Round 1 F3, Round 3 F2) already established that `extract_fasta_from_maf.py` and `horizon_eval` evaluate MAF-extracted `(ref_id, tgt_id)` pairs; `manuscript/horizon_paper.md:59` discloses this, but `horizon_paper.md:10` (Abstract) does not. Furthermore, `horizon_paper.md:59, 65-70` and `docs/BENCHMARK_RESULTS.md:42` label the nuclear slices "held-out" even though `Figure 6B` (`horizon_paper.md:86`) explicitly re-evaluates horse, macaque, and marmoset to measure `R1-HMM` BOCPD gate precision gains. Finally, `horizon_paper.md:10` highlights `0` false-positive pairs on `18%`-GC `10 kb` random DNA while omitting the `288` (`30%`-GC) and `149` (`45%`-GC) false-positive pairs reported at `horizon_paper.md:59` and `docs/BENCHMARK_RESULTS.md:13`.
- **Repair:** Endorsed (drop "held-out" on re-evaluated nuclear slices, qualify the Abstract on MAF-slice evaluation, and report the `30%`/`45%`-GC false positives alongside `18%`-GC); **the fix is free**.

### F5 — C5 (`BandGrid` index formula)
- **Ruling:** **OVERTURNED** on the attacker's code-inspection claims; **NARROWED** only on the companion-text typo.
- **What survives:** The attacker claimed `idx(i, j) = i · (2W + 1) + (j − i + W)` matches `BandLayout::idx` in `src/mdl/dp.rs:164-169` and that centering on `j = i` instead of `j = i · (M / N)` makes `(N, M)` unreachable when `|N − M| > W`. Both claims are false in the live code: `src/mdl/dp.rs:160-172` defines `BandLayout::slot` (not `BandLayout::idx`), and `self.center(i)` (`src/mdl/dp.rs:155-157, 244-246`) calls `band_center(i, m, n) = (i * m) / n.max(1)`, which **already** centers the band on `j = i · (M / N)`. The attacker should have checked `BandLayout::slot` and `band_center` in `src/mdl/dp.rs:155-172, 244-246`. Only the secondary companion-text discrepancy holds: `manuscript/horizon_companion.md:142` writes `j − i + W` (simplified `j = i` centering) rather than `j − center(i) + W`, and `manuscript/horizon_companion.md:179` writes `(6 + 1) × (2 · 180 + 1)` instead of `7 × 13` (because `BandLayout::new` at `src/mdl/dp.rs:142` clamps `bw = band_w.min(n.max(m))`).
- **Repair:** Update `manuscript/horizon_companion.md:142, 179` to match `BandLayout::slot` / `band_center` and `7 × 13`; **the fix is free**.

### F6 — C6 (`# EG:Z:` tags omit nested inversions and cannot distinguish collapsed paralogs from single-copy synteny)
- **Ruling:** **UPHELD** (attacker's **REFUTED** verdict on C6 holds).
- **What survives:** `egraph::collect_dups` (`src/egraph.rs:518-543`) only collects leaves under `Expr::Dup` and does not record `Expr::Inv`. Furthermore, `format_gfa_paralog_tags` (`src/formats.rs:153-158`) and `src/merge.rs:361` emit `# EG:Z:` as a `#` GFA comment line containing `bin_bp`-bucketed indices (`blk[0].0 / bin_bp`, `src/formats.rs:142`), which standard GFA parsers ignore and which lack segment/walk IDs needed to separate collapsed paralogs from single-copy synteny (distinct from `horizon_resolve --tag-walks` at `src/bin/horizon_resolve.rs:476`, which writes `EG:Z:` on `W` lines).
- **Repair:** Endorsed (`manuscript/horizon_companion.md:52`); **the fix is free**.

### F7 — Claim-list steering note
- **Ruling:** **UPHELD** (meta-observation).
- **What survives:** All three omitted inference claims cited by the attacker are verifiable in the repository: hardcoded `λ̂ = 10.0` bits and `KARLIN_MDL_COEFF = 0.55` (`src/mdl/mod.rs:26, 71, 148`, `manuscript/horizon_paper.md:40`), "minimum-description-length ancestral order" (`manuscript/horizon_paper.md:44`) without a phylogenetic tree in `src/egraph.rs`, and `horizon_wga --graph` (`manuscript/horizon_paper.md:82`, `docs/BENCHMARK_RESULTS.md:25`) evaluated on three short queries (`400`-bp all-`N`, `258`-bp, `356`-bp) without precision/recall ground truth.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `horizon_gfa stats` is counter-only (`horizon_gfa.rs:20-56`) while `vg` builds an in-memory graph; free text/figure-label fix |
| F2 | UPHELD | `48.63%` excludes `21.4%` (`60.4M`) of emitted pairs (`38.22%` unconditional) and `2.01×` multiplicity on `1-to-1` loci is unverified; free fix |
| F3 | NARROWED | `horizon_guide --progressive --hal` writes a 4-seg contig-0 star-tree stub (`horizon_guide/main.rs:94-103`, `hal/mod.rs:90-171`), but `horizon_merge` runtime/counts survive with qualification; free fix |
| F4 | UPHELD | "Held-out" nuclear slices were re-evaluated in `Figure 6B` and Abstract omits MAF-slice framing and `30%`/`45%`-GC false positives (`288`/`149`); free fix |
| F5 | OVERTURNED | `src/mdl/dp.rs:155-172, 244-246` defines `BandLayout::slot` (not `idx`) and already centers on `band_center(i, m, n) = i*m/n`; only companion doc `142, 179` has a free typo fix |
| F6 | UPHELD | `collect_dups` (`egraph.rs:518-543`) ignores `Expr::Inv` and `# EG:Z:` (`formats.rs:153-158`) emits coarse `bin_bp` comment lines, not per-segment/walk tags; free fix |
| F7 | UPHELD | Omitted CRC constants (`10.0`/`0.55`), tree-free `egraph` "ancestral order", and 3-query `--graph` evaluation are real and unexamined |
