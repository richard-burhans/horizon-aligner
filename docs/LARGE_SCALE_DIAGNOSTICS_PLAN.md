# Large-Scale Diagnostics, Bug-Finding & Performance Plan

**Status:** Executed (`B1`, `A1`, `A2`, `C1`, `B2`, `A3`, `A4`, `B3` complete; `C2` pending external tools / `.hal` fixture)  
**Scope:** Using existing large-scale alignment artifacts (`12 GB` HPRC v1.0 GFA, `chr6` MAF, `1.67 GB` `cactus241way_chr22_1MB.maf`, `97 MB` `gibbon_y_tba.maf`) to uncover correctness bugs, calibration gaps, and memory/runtime bottlenecks in `horizon-aligner` before whole-genome scaling.

---

## 1. Available Large-Scale Artifacts

| Artifact | Location | Size / Scale | Primary Diagnostic Use |
|---|---|---|---|
| HPRC v1.0 Minigraph-Cactus GFA | `../pg/data/benchmarks/hprc-v1.0-mc-grch38.gfa.gz` | `12 GB`, `81.4M` segments, `113M` links | Streaming GFA invariants, `10k→200k` subgraph scaling, `chr6` `1-to-1` & `multi` ground truth |
| Full `GRCh38.chr6 × CHM13.chr6` MAF | `/tmp/chr6_fixed.maf` | `624 MB`, `282.2M` pairs | Precision-gap stratification, `pairs_vec` RAM profiling, `-` strand coordinate checks |
| `cactus241way_chr22_1MB.maf` | `../pg/data/benchmarks/cactus241way_chr22_1MB.maf` | `1.67 GB`, `31` contigs | Multi-contig `horizon_guide --progressive --hal` & `horizon_merge` accuracy |
| `cactus447way` / `multiz470way` | `../pg/data/benchmarks/` | `69 MB` / `52 MB` | Multi-species coordinate & strand consistency |
| Gibbon Y 5-way TBA MAF | `../gibbon_y_tba.maf` | `97 MB`, `17,787` blocks | Repeat-rich / ampliconic Y-chromosome `horizon_resolve --collapse-dups` stress test |

---

## 2. Workstream A — Differential Bug-Finding & Invariant Checks

### A1. Stratify the `chr6` Precision Gap (`38.22%` unconditional vs. `48.63%` conditioned)
- **Target:** `scripts/benchmarks/compare_chr6_horizon_vs_hprc.py`
- **Action:** Extend the script (which already loads `hprc_multi_x` and computes `diff = |y_pred - y_truth|`) to bin the `61.78%` non-exact pairs into:
  1. Multi-copy / paralogous HPRC segments (`hprc_multi_x`),
  2. Small boundary / gap-placement shifts (`1–10 bp`, `11–100 bp`, `101–1000 bp`),
  3. True off-diagonal false positives (`> 1 kb` or outside HPRC shared segments, stratified by `'N'`-masked / low-complexity vs. unique sequence).
- **Goal:** Identify whether over-emission (`282.2M` pairs) is driven by `--resolve-dups` paralog retention, X-drop over-extension, or seeding in repetitive intervals.

### A2. Automated MAF / GFA / Liftover Invariant Suite (`scripts/benchmarks/validate_outputs.py`)
- **Motivation:** Negative-strand MAF coordinate bugs have occurred previously (`pg` `0c47500`, `chr6` coordinate helpers).
- **Invariants to assert on real outputs (`chr6`, `cactus241way`, `gibbon_y_tba.maf`):**
  1. **MAF–FASTA consistency:** For every `s` line, ungapped sequence length equals `size`, `start + size <= srcSize`, and the sequence matches the underlying FASTA slice (`+` strand) or its reverse complement (`-` strand).
  2. **GFA walk/path conservation:** `horizon_clean` and `horizon_resolve` preserve walk/path connectivity (every `L` link traversed by `W`/`P` exists with matching orientation; `S` segment lengths match walk spans).
  3. **Liftover roundtrip:** `horizon_lift` mappings `x → y → x` are idempotent on `1-to-1` aligned segments.

### A3. Multi-Contig `horizon_guide` / `horizon_merge` Validation on `cactus241way_chr22_1MB.maf`
- **Target:** `scripts/benchmarks/verify_merge_guide_maf.py`
- **Action:** Generalize `verify_merge_guide_maf.py` from `16.5 kb` `chrM.maf` to `cactus241way_chr22_1MB.maf` (`31` contigs), verifying both concatenated-coordinate pair recall/precision and multi-contig `.hal` `SEQNAME_ARRAY` / `BOTTOM_ARRAY` offsets produced by `HalAlignment::from_multi_contig_star_blocks`.

### A4. Ampliconic Repeat Stress Test on `gibbon_y_tba.maf`
- **Action:** Extract the 5 gibbon Y species sequences (`mHylPil1`, `mHylMol1`, `mHooLeu1`, `mNomLeu1`, `mSymSyn1`) from `../gibbon_y_tba.maf` and run `horizon_wga` + `horizon_resolve --collapse-dups`.
- **Goal:** Verify that `MdlCostModel` collapses tandem ampliconic repeats without collapsing distinct X-degenerate gametologs or introducing spurious `2-way` repeat matches.

---

## 3. Workstream B — Performance & Scaling Diagnostics

### B1. Replace `Vec<(usize, usize)>` with Run-Encoded Intervals (`src/mdl/mod.rs`, `src/mdl/dp.rs`)
- **Problem:** At `282.2M` pairs on `GRCh38.chr6 × CHM13.chr6`, `pairs_vec: Vec<(usize, usize)>` allocates **~4.5 GB** (`16 B/pair`) before `sort_unstable()` + `dedup()`.
- **Action:** Represent contiguous ungapped/colinear runs as `(x_start, y_start, len, is_rc)` in `ExtensionOutput` / `ContigAlignment`, expanding to per-base pairs only if a caller explicitly requests them.
- **Expected impact:** **10–50×** reduction in alignment-pair memory and sort time on full chromosomes, unblocking whole-genome (`chr1`–`chr22`) runs.

### B2. HPRC Subgraph Scaling Curve for `horizon_wga --graph` (`src/mdl/dag.rs`)
- **Observation:** On the `10k`-segment `chr22` subgraph, `horizon_gfa stats` used `22.7 MB` RSS while `horizon_wga --graph` used `154.5 MB` RSS (`DagBandGrid` + `DagSeedIndex`).
- **Action:** Extract `10k → 50k → 200k`-segment subgraphs from `hprc-v1.0-mc-grch38.gfa.gz` and measure wall time + peak RSS to locate any super-linear scaling in `DagBandGrid` allocation or `DagSeedIndex` lookups.

### B3. Automated `[profile]` Stage-Timing & RSS Regression Ledger
- **Action:** Record `horizon_wga` `[profile]` stage timings (`w12`, `w10`, `lz1`, `seed_end`, `dp+xd`) and peak RSS across `chr22_peri_5m`, `chr6_q20m`, and full `chr6` into a tracked JSONL ledger so quadratic chaining or DP band regressions are caught immediately.

---

## 4. Workstream C — Additional Improvements

### C1. LCB/Contig-Scale Conformal Risk Control (CRC) Calibration (`src/mdl/dp.rs`)
- **Observation:** `calibrate_null_crc_for_gc` currently calibrates `λ̂` (`n_cal = 15`, `alpha = 0.07`) on `64 bp` DP decoy windows, while `scripts/benchmarks/bench_random_gc_crc.py` still observes `0–502` false-positive pairs on `10 kb` random DNA (`18%–45%` GC) due to the larger LCB/X-drop search space.
- **Action:** Scale `λ_eff` with the search-space size (`log2(len_x * len_y)`) or calibrate `n_cal = 15` directly against decoy LCB/X-drop draws so multi-seed random-GC false positives drop to near zero across all seeds, not just `seed=999`.

### C2. Pending External Baselines & External `.hal` Audit
- **When installed in the environment:**
  - Run `FastGA`, `LASTZ`, `KegAlign`, `odgi`, and `gfatools` on identical hardware and inputs (`R2-RESOLVE-LIFT`).
  - Audit `src/hal/hdf5_read.rs` and `decode_hal_*` against a real external Cactus-generated `.hal` fixture and un-ignore `real_hal_decode_from_env` (`R2-HMM-XDROP`).

---

## 5. Recommended Execution Order

1. **B1 (Run-encoded `pairs_vec`)** — highest leverage for full-chromosome RAM and runtime.
2. **A1 + A2 (`chr6` error stratification + `validate_outputs.py`)** — surfaces exact causes of the `38.22%` unconditional precision gap and guards all MAF/GFA/liftover coordinates.
3. **C1 (LCB-scale CRC calibration)** — eliminates remaining multi-seed random-GC false positives.
4. **B2 + A3 + A4 (HPRC subgraph scaling, `cactus241way` multi-contig, `gibbon_y_tba` repeat stress)** — validates graph and multi-species scaling.
5. **B3 + C2 (`[profile]` ledger + external baselines/HAL when tools are available)**.
