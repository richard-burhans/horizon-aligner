# Round 5 Attack -- Sampling Lens (`sampling`)

**Target documents:** [`reports/rounds/round-5/ATTACKED-horizon_paper.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-5/ATTACKED-horizon_paper.md), [`reports/rounds/round-5/ATTACKED-horizon_companion.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-5/ATTACKED-horizon_companion.md), [`docs/BENCHMARK_RESULTS.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md), and [`manuscript/data/`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/data).

**Method cards & literature check (Rules 2–3):**
- This project keeps **no method cards** (`BRIEF-sampling.md` Rule 3).
- `../library` is not present on disk (`NO_LIBRARY`), so literature checks were performed against `refs/catalog_index.json`, `manuscript/references.bib`, and repository artifacts alone; no claims about prior papers are asserted from memory.

**Prior rounds check (Rule 4):**
- Checked `scripts/round history "sampling" --kind adjudication` and [`reports/rounds/round-4/adjudication-L2.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-4/adjudication-L2.md). Prior rounds already adjudicated: (a) `CHM13.chr6` q-arm-only sampling (`r4 F1`, now disclosed in `C2` as `61401689–172126628`), (b) `extract_fasta_from_maf.py` MAF-block concatenation (`r2 F8`, `r3 F1`, `r4 F4/F6`, now disclosed in `docs/BENCHMARK_RESULTS.md:5` but still omitted in `horizon_paper.md`), and (c) uncalibrated `n_cal = 15` (`r4 R1-HMM`). All findings below either address new Round 5 artifacts or explicitly note what changed since Round 4.

---

## F1 -- C1: The `10,000`-segment HPRC `chr22` subgraph is `96.09%` `N`-bases (`10.51 Mbp` in one segment) and mixes in `chr14` / `chrUn` contigs

- **Claim:** `C1` (`NARROWED`)
- **Concrete flaw:** The `10,000`-segment (`10,937,960` bp) HPRC `chr22` subgraph used for `horizon_gfa stats` vs. `vg v1.60.0` (and also `horizon_clean`, `horizon_resolve`, and `horizon_wga --graph` in Figure 6A) is an `n = 1` prefix sample whose sequence is `96.09%` `'N'`-gap bases and includes non-`chr22` unplaced contigs.
- **Evidence:**
  1. Direct inspection of [`/tmp/hprc_chr22_sub_real.gfa`](file:///tmp/hprc_chr22_sub_real.gfa) shows `10,000` `S`-lines totaling `10,937,960` bp, of which **`10,510,000` bp (`96.09%`)** are `'N'` bases inside a single segment (`S\t76067108\tNNNN...`), leaving only **`427,960` bp** of non-`N` ACGT sequence across the remaining `9,999` segments (`0.012%` of the `81,415,956`-segment HPRC graph).
  2. `SN:Z:` tags in `/tmp/hprc_chr22_sub_real.gfa` show only `3,991` segments belong to `GRCh38.chr22`; the prefix sample also includes `1,978` segments from `GRCh38.chrUn_GL000224v1`, `1,919` from `GRCh38.chr22_KI270736v1_random`, and `1,230` from `GRCh38.chr14_KI270722v1_random` (a different chromosome).
  3. The `P`-line (`GRCh38.chr22_sub`) has `2,000` steps and the `4` `W`-lines (`CHM13`) are each `18,001` characters, with `n = 1` subgraph and `n = 1` timed run.
- **Correction:** Narrow `C1` in [`ATTACKED-horizon_paper.md:10, 82`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-5/ATTACKED-horizon_paper.md#L82) and [`docs/BENCHMARK_RESULTS.md:25, 76`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L25) to state that this is an `n = 1` prefix-sampled `10,000`-segment subgraph (`427.96 kb` ACGT + `10.51 Mbp` `N`-gap, including `chrUn`/`chr14` unplaced contigs), not a representative `10.94 Mbp` euchromatic `chr22` graph.

---

## F2 -- C2: HPRC `chr6` precision (`48.63%`) conditions on `pairs_on_hprc_1to1_x`, excluding `60.43M` (`21.4%`) off-target pairs driven by `CHM13.chr6` q-arm-only sampling

- **Claim:** `C2` (`NARROWED`)
- **Prior round note:** Round 4 `F1` narrowed `C2` because `CHM13.chr6` was q-arm only (`61401689–172126628`). Round 5 now discloses those coordinates, **but** the new HPRC concordance metric introduces a denominator-conditioning bias from pairing full `GRCh38.chr6` (`170.8 Mbp`, including the `0–61 Mbp` p-arm and `HLA`) against q-arm-only `CHM13.chr6` (`110.7 Mbp`).
- **Concrete flaw:** `C2` reports `48.63%` many-to-many exact pair precision by restricting the precision denominator to `pairs_on_hprc_1to1_x` (`221,734,653` pairs), which drops **`60,428,393` (`21.4%`)** of Horizon's `282,163,046` emitted pairs that align `GRCh38.chr6` bases outside `hprc_1to1_x` (including `14,584,569` unique `GRCh38.chr6` bases, largely on the unmatched p-arm) onto `CHM13.chr6` q-arm.
- **Evidence:**
  1. [`/tmp/extract_hprc_chr6_truth.py:116-129`](file:///tmp/extract_hprc_chr6_truth.py#L116-L129) filters HPRC segments to `1,628,018` `1-to-1` segments (`110,082,990` bp; excluding `19` multi-copy segments / `38` pairs).
  2. [`/tmp/compare_chr6_horizon_vs_hprc.py:136-148, 213`](file:///tmp/compare_chr6_horizon_vs_hprc.py#L136-L148) and [`manuscript/data/chr6_hprc_comparison.json:10-28`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/data/chr6_hprc_comparison.json#L10-L28) compute `exact_coord_precision_on_shared = 107,836,669 / 221,734,653 = 0.486332` (`48.63%`), omitting `282,163,046 - 221,734,653 = 60,428,393` off-target pairs (`horizon_unique_x_bp = 124,160,076` vs. `horizon_and_hprc_x_bp = 109,575,507`, a `14,584,569` bp excess).
  3. Over all `282,163,046` emitted pairs, Horizon's unconditional exact pair precision is **`38.22%`** (`107,836,669 / 282,163,046`), whereas `minimap2 -x asm5` ([`manuscript/data/chr6_minimap2_summary.json:3-15`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/data/chr6_minimap2_summary.json#L3-L15)) has only `171,202` (`0.16%`) off-target pairs (`109,759,180 / 110,204,880 = 99.60%` unconditional exact precision).
- **Correction:** Narrow `C2` in [`ATTACKED-horizon_paper.md:10, 76`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-5/ATTACKED-horizon_paper.md#L76) and [`docs/BENCHMARK_RESULTS.md:29, 74`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L29) to report unconditional exact pair precision (`38.22%` across all `282.16M` pairs) alongside the `hprc_1to1_x`-conditioned `48.63%`, and state that `60.43M` (`21.4%`) aligned pairs fall outside `hprc_1to1_x` on `n = 1` chromosome pair.

---

## F3 -- C3: `horizon_merge` (`92/100` `chrM` species) samples only a `16.5 kb` organelle, and `horizon_guide` (`5` mammals) samples `11–33 KB` MAF-concatenated `chr22_random` slices

- **Claim:** `C3` (`NARROWED`)
- **Concrete flaw:** Neither multi-species benchmark samples full nuclear chromosomes: `horizon_merge` nets `92` of `100` `chrM` (`~16.5 kb` haploid, repeat-free mitochondrial) sequences, while `horizon_guide` aligns `11–33 KB` MAF-concatenated extracts from `chr22_KI270735v1_random` (`4` HAL bottom segments, `n = 1` run).
- **Evidence:**
  1. [`docs/BENCHMARK_RESULTS.md:10`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L10) confirms `horizon_merge` uses `chrM.maf` (`92` of `100` species; `8` excluded).
  2. [`scripts/benchmarks/run_all_evals.py:68-75`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/run_all_evals.py#L68-L75) and [`scripts/benchmarks/extract_fasta_from_maf.py:41-45, 78-85`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/extract_fasta_from_maf.py#L41-L45) extract the `5`-mammal inputs (`hg38`, `Macaca_mulatta`, `HLcalPym1`, `Equus_caballus`, `Canis_lupus`) from `chr22_KI270735v1_random` by concatenating only MAF-aligned blocks, yielding `11–33 KB` FASTAs (`cactus447way_Equus_caballus.fa` `11 KB`, `cactus447way_Macaca_mulatta.fa` `22 KB`, `multiz470way_HLcalPym1.fa` `33 KB`) and only `4` HAL bottom segments.
- **Correction:** Narrow `C3` in [`ATTACKED-horizon_paper.md:10, 82`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-5/ATTACKED-horizon_paper.md#L82) to specify `~16.5 kb` `chrM` (`92/100` species) for `horizon_merge` and `11–33 KB` MAF-concatenated `chr22_KI270735v1_random` slices (`4` HAL bottom segments) for `horizon_guide`.

---

## F4 -- C4: `C4` omits MAF-block concatenation (`11–387 KB` slices including new `B10K`), cherry-picks `18%`-GC (`0` pairs) over `30%`/`45%`-GC (`288`/`149` pairs) with no null script, and cites unimplemented `n_cal = 15`

- **Claim:** `C4` (`NARROWED`)
- **Prior round note:** Round 2 `F8`, Round 3 `F1`, and Round 4 `F4/F6` already established that `extract_fasta_from_maf.py` strips unaligned sequence (now disclosed in `docs/BENCHMARK_RESULTS.md:5`), and `r4 R1-HMM` carried `n_cal = 15`. In Round 5, `B10K 363-way` is newly added, and `horizon_paper.md` still omits the MAF-concatenation and random-GC sampling details.
- **Concrete flaw:**
  1. **MAF-concatenated slices:** `B10K 363-way` (`chr1:5.9M–6.4M`, `500 kb`) was extracted into [`/tmp/b10k_galGal.fa`](file:///tmp/b10k_galGal.fa) (`386,668` bp) and [`/tmp/b10k_taeGut.fa`](file:///tmp/b10k_taeGut.fa) (`349,383` bp), stripping `~113 kb` of unaligned sequence; `cactus447way` (`Equus_caballus` `11 KB`, `Macaca_mulatta` `22 KB`) and `multiz470way` (`HLcalPym1` `33 KB`) are tiny `chr22_KI270735v1_random` MAF-concatenated extracts.
  2. **Random GC null cherry-picking & missing script:** `ATTACKED-horizon_paper.md:10` (`C4`) cites only `18%`-GC random DNA (`0` pairs in `0.0069 s`), omitting `30%`-GC (`288` false-positive pairs in `0.4150 s`) and `45%`-GC (`149` false-positive pairs in `0.3829 s`) reported in [`docs/BENCHMARK_RESULTS.md:13`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L13). Moreover, `grep` across `src/`, `scripts/`, and `tests/` finds no script or test generating these `10,000`-bp random GC nulls (replicate count `n` and seed are untracked).
  3. **Unimplemented `n_cal = 15`:** `ATTACKED-horizon_paper.md:40` and `ATTACKED-horizon_companion.md:153` cite `α = 0.07` (`n_cal = 15`), but `n_cal` does not exist in `src/` or `tests/` (`min_delta_bits` is hardcoded to `10.0` in [`src/mdl/mod.rs:71, 142`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/mod.rs#L71)).
- **Correction:** Narrow `C4` in [`ATTACKED-horizon_paper.md:10, 40, 59`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-5/ATTACKED-horizon_paper.md#L59) to disclose MAF-block concatenation (`11–387 KB` slices), report all three random GC null bins (`18%`: `0`, `30%`: `288`, `45%`: `149` pairs) with a tracked script/replicates, and remove or implement `n_cal = 15`.

---

## F5 -- C5: `BandGrid` flat buffer indexing is a deterministic layout identity with no sampling component

- **Claim:** `C5` (`SURVIVED`)
- **Evidence:** Checked [`src/mdl/dp.rs:140-172`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/dp.rs#L140-L172) (`BandLayout::slot` computes `i * self.stride + col` with `stride = 2 * bw + 1`, `col = (j - center(i)) + bw`). This is a deterministic array-indexing formula, not a statistical sampling claim.

---

## F6 -- C6: `# EG:Z:` tag emission is a format claim outside the sampling lens

- **Claim:** `C6` (`UNCHALLENGED`)
- **Evidence:** Checked [`src/formats.rs:153-158`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/formats.rs#L153-L158) and [`src/bin/horizon_resolve.rs:476`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_resolve.rs#L476). `C6` is a structural format claim; the sampling flaw in `horizon_resolve`'s empirical `10k`-subgraph benchmark is already filed under `F1`.

---

## Verdict Table

| Claim | Verdict | Summary |
| :--- | :--- | :--- |
| C1 | NARROWED | `10,000`-segment `chr22` subgraph is `96.09%` `N`-bases (`10.51 Mbp` in 1 segment, `427.96 kb` ACGT) and mixes in `chrUn`/`chr14` contigs (`n = 1`). |
| C2 | NARROWED | `48.63%` exact precision conditions on `pairs_on_hprc_1to1_x`, dropping `60.43M` (`21.4%`) off-target pairs from full `GRCh38.chr6` vs. `CHM13.chr6` q-arm (`38.22%` unconditional). |
| C3 | NARROWED | `horizon_merge` samples `~16.5 kb` `chrM` (`92/100` species); `horizon_guide` samples `11–33 KB` MAF-concatenated `chr22_KI270735v1_random` slices (`4` HAL bottom segments). |
| C4 | NARROWED | `horizon_paper.md` omits MAF-block concatenation (`11–387 KB` slices incl. `B10K`), cherry-picks `18%`-GC (`0` pairs) over `30%`/`45%`-GC (`288`/`149` pairs) with no script, and cites unimplemented `n_cal = 15`. |
| C5 | SURVIVED | `BandGrid` indexing (`src/mdl/dp.rs:140-172`) is a deterministic layout identity with no sampling defect. |
| C6 | UNCHALLENGED | `# EG:Z:` tag emission is a structural format claim; `horizon_resolve` `10k`-subgraph sampling bias is covered in `F1`. |
