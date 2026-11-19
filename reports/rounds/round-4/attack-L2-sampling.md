# Round 4 Attack Report — L2 (`sampling`)

**Target:** `manuscript/horizon_paper.md` (claims `C1`–`C6`), verified against `HEAD` (`33bede8`), `docs/BENCHMARK_RESULTS.md`, `scripts/benchmarks/`, `manuscript/scripts/generate_figures.py`, `src/`, and `/tmp/` benchmark artifacts.
**Methodological disclosures (Rules 2–4):**
- **Literature (Rule 2):** `../library` is not present on disk (`ls -la ../library` failed), so I worked from `refs/catalog_index.json` alone and make no assertions about held paper contents from memory.
- **Method cards (Rule 3):** Per `BRIEF-sampling.md`, this project keeps **no method cards** yet.
- **Prior rounds (Rule 4):** Checked `scripts/round history` for `C1`–`C6`, `horizon_clean`, and `chr6`. Prior rounds already adjudicated: `horizon_clean` real-data no-op (`round 1's F9`, `round 2's R2-RESOLVE-LIFT`, `round 3's R2-RESOLVE-LIFT`), missing `minimap2`/`odgi`/`vg`/`LASTZ` baselines (`round 2's R2-RESOLVE-LIFT`, `round 3's F7`), and missing `/tmp/*.fa` extraction scripts / `CHM13.chr6` coordinate window (`round 3's F3`). Where those issues persist or take a new form on `HEAD` (`33bede8`), I cite the prior ruling and state what changed.

---

## F1 -- C2 samples only the CHM13.chr6 q-arm (`110.7 Mbp`), `n = 1` chromosome and `n = 1` run, with no accuracy metric and stale pre-`33bede8` timings

- **Claim:** `C2` (`manuscript/horizon_paper.md:101`) — *"full `GRCh38.chr6` (`170.8 Mbp`) × `CHM13.chr6` (`110.7 Mbp`) finished in `451.35 s` (`7m 31s`, `10,515.0 MB` RSS)—a `2.35×` total speedup over our `1,062.45 s` baseline (`12.3×` speedup in `LASTZ12` transition seeding from `588.67 s` to `47.87 s` via 16-bit prefix-bucketed `FlatSeedIndex` and Rayon parallelization; Figure 5)."*
- **Verdict:** `NARROWED`
- **Concrete flaws (sampling lens):**
  1. **Target chromosome is q-arm only, not "full `CHM13.chr6`":** Inspecting `/tmp/hprc_chm13_chr6.fa` shows the header `>CHM13.chr6:61401689-172126628` (`110,724,939` bp). The entire `0–61.4 Mbp` p-arm and centromeric/pericentromeric array (`35.7%` of `CHM13.chr6`) is absent from the HPRC `v1.0` walk. Unlike `chr22` (`manuscript/horizon_paper.md:101`, which qualifies `95.6%` of `CHM13.chr22` q-arm), `C2` labels this `110.7 Mbp` q-arm slice as "full `GRCh38.chr6` × `CHM13.chr6`".
  2. **`n = 1` chromosome, `n = 1` run, no variance, and no alignment accuracy/coverage:** Only `1` human chromosome pair (`chr6`) and `1` timed run (`451.35 s` vs. `1,062.45 s`) are sampled, with zero replicate runs, no confidence interval, and zero recall/precision/unique-base coverage metrics (so speedup cannot be separated from output quality).
  3. **Commit and workload confounding:** `docs/BENCHMARK_RESULTS.md:64-65` shows `1e47fe0` (`1,062.45 s`) and `53130b8` (`451.35 s`) did not even produce the same workload (`292,534` -> `301,156` fwd LCBs; `693,242` -> `778,471` DP jobs), and `HEAD` (`33bede8`) subsequently modified 20 files including `src/seeding/mod.rs`, `src/seeding/flat_index.rs`, `src/seeding/translated.rs`, `src/mdl/dp.rs`, and `src/mdl/mod.rs` without re-benchmarking `chr6`.
- **Evidence:**
  - `head -n 2 /tmp/hprc_chm13_chr6.fa /tmp/hprc_grch38_chr6.fa` -> `>CHM13.chr6:61401689-172126628` and `>GRCh38.chr6`.
  - `docs/BENCHMARK_RESULTS.md:27,64-65` (`53130b8` vs `1e47fe0`).
  - `git show --stat 33bede8` (modifies `src/seeding/mod.rs`, `src/seeding/flat_index.rs`, `src/seeding/translated.rs`, `src/mdl/dp.rs`, `src/mdl/mod.rs` after `53130b8`).
- **Correction:** Narrow `manuscript/horizon_paper.md:101` to state `GRCh38.chr6` (`170.8 Mbp`) × `CHM13.chr6:61.4M–172.1M` q-arm (`110.7 Mbp`), report `n = 1` single-run timing on commit `53130b8` (or re-run with replicates and unique-base coverage on `HEAD` `33bede8`).

---

## F2 -- C3 `horizon_gfa extract` samples only `2` chromosomes and `2` haplotypes out of `2,580` paths and `24,456` walks

- **Claim:** `C3` (`manuscript/horizon_paper.md:101`) — *"On the `12 GB` HPRC `v1.0` Minigraph-Cactus pangenome graph (`hprc-v1.0-mc-grch38.gfa.gz`), `horizon_gfa stats` streamed all `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks in `99.97 s` using `97.0 MB` peak RSS, while two-pass `horizon_gfa extract` isolated `GRCh38.chr22` + `CHM13.chr22` in `183.29 s` (`456.4 MB` RSS) and `GRCh38.chr6` + `CHM13.chr6` in `244.93 s` (`1,255.1 MB` RSS)."*
- **Verdict:** `NARROWED`
- **Concrete flaws (sampling lens):**
  1. **Subgraph extraction sampled on only `2` chromosomes (`chr22`, `chr6`) and `2` haplotypes (`GRCh38`, `CHM13`):** Out of `2,580` paths and `24,456` walks in `hprc-v1.0-mc-grch38.gfa.gz` (`11,947,915,887` bytes), `horizon_gfa extract` was only sampled on `2`-haplotype extractions for `chr22` (`921,252` cached segments, `456.4 MB` RSS) and `chr6` (`3,057,602` cached segments, `1,255.1 MB` RSS). Because `extract` caches all needed segments in RAM (`docs/BENCHMARK_RESULTS.md:56-57`), memory scales with extracted segment count rather than remaining `O(1)`; larger chromosomes (`chr1`) or multi-haplotype cohorts (e.g., all 90+ HPRC haplotypes) are unsampled.
  2. **Single-run timing (`n = 1`) and pre-`33bede8` measurement:** `chr6` extraction was pinned to commit `b857377` (`docs/BENCHMARK_RESULTS.md:57`), prior to `33bede8` editing `src/bin/horizon_gfa.rs` (19 lines) and `src/gfa/mod.rs` (16 lines).
- **Evidence:**
  - `docs/BENCHMARK_RESULTS.md:55-57`.
  - `git show --stat 33bede8` (`src/bin/horizon_gfa.rs`, `src/gfa/mod.rs`).
- **Correction:** Qualify `manuscript/horizon_paper.md:101` that `horizon_gfa extract` was measured on `n = 1` run per pair for `2`-haplotype (`GRCh38` + `CHM13`) `chr22` and `chr6` extractions (with RSS scaling with cached segment count), and re-measure on `33bede8`.

---

## F3 -- C4 `horizon_clean` has zero empirical sampling on sequence-bearing GFAs (only a `0`-severed `*-seq` no-op and 2 toy repeat tests)

- **Claim:** `C4` (`manuscript/horizon_paper.md:49`) — *"3. `horizon_clean` (`src/bin/horizon_clean.rs`) indexes `S`-segment flanks and scores candidate `L`-line junctions with `bocpd_gc_profile` and an MDL two-regime changepoint test (`--min-cp-bits`), severing unsupported (`walk_support == 0`) discordant edges into `HB:i:1` Horizon Bubbles."*
- **Verdict:** `NARROWED`
- **Prior round context:** `round 1's F9` and `round 2/3's R2-RESOLVE-LIFT` flagged `horizon_clean`'s real-data no-op. On `HEAD` (`33bede8`), the code implements `score_junction_cp_bits` (`src/bin/horizon_clean.rs:100-127`), **but** the benchmark gap remains unresolved.
- **Concrete flaws (sampling lens):**
  1. **Zero real-data or sequence-bearing GFA evaluation:** In `docs/BENCHMARK_RESULTS.md:66`, the **only** benchmark of `horizon_clean` in the repository is `horizon_clean --gfa -` on a `*-seq` pairwise `chr6` peri GFA (`2,776 kept, 0 severed (no-op on *-seq)`). When `S`-segment sequences are `*`, `score_junction_cp_bits` returns `0.0` for every link (`src/bin/horizon_clean.rs:154-176,252-260`).
  2. **Toy unit-test sample only:** The only tests exercising edge severing (`src/bin/horizon_clean.rs:559-626`) are two synthetic 60-bp 100%-AT vs. 100%-GC dinucleotide repeat (`ATATAT...` vs. `GCGCGC...`) fixtures. Neither real HPRC junctions nor realistic chimeric/discordant pangenome links have been sampled to measure precision, recall, or `--min-cp-bits` calibration.
- **Evidence:**
  - `docs/BENCHMARK_RESULTS.md:66`: `2,776 kept, 0 severed (no-op on *-seq)`.
  - `src/bin/horizon_clean.rs:100-127, 252-260, 559-626`.
- **Correction:** Clarify in `manuscript/horizon_paper.md:49` that `horizon_clean`'s junction-severing behavior is currently verified only on synthetic 60-bp AT/GC unit-test fixtures (`src/bin/horizon_clean.rs:559-626`) and has not yet been benchmarked on sequence-bearing HPRC or pangenome GFAs.

---

## F4 -- C5 rests on untracked *Plasmodium* / *Cannabis* samples with an internal denominator contradiction (`109,565` vs `119,582`) and MAF-block selection bias

- **Claim:** `C5` (`manuscript/horizon_paper.md:65`) — *"End-to-end `horizon_wga` abstained on `87.7%` of those MDL-rejected columns while recovering `54.4%` of MDL-supported pairs (`66.8%` with `w = 10` sparse-exon fallback and `band_w = 180`). Applied with zero parameter changes to intraspecific *Cannabis sativa* (`cs10` × `Cherry_Pie`, `45.5%` GC, `418,663` columns), Horizon retained `99.95%` (`418,465`) of reference `MAF` columns and `1.000` recall on scratch realignments."*
- **Verdict:** `REFUTED`
- **Concrete flaws (sampling lens):**
  1. **Untracked data and unreproducible numbers:** `grep -rn "Cherry_Pie\|418,663\|418465\|87\.7" src/ scripts/ docs/` returns **zero** matches in the codebase, `scripts/benchmarks/`, or `docs/BENCHMARK_RESULTS.md`. Neither the *P. falciparum* × *P. vivax* `chr1` dataset nor the *Cannabis sativa* `cs10` × `Cherry_Pie` dataset is present or scripted in the repository.
  2. **Contradictory sample denominators:** `manuscript/horizon_paper.md:65` reports `24,748 / 119,582` (`20.7%`) MDL-rejected columns for `Pf` × `Pv` `chr1`, whereas `manuscript/scripts/generate_figures.py:218-219` plots `n = 88,597` MDL-supported + `n = 20,968` MDL-rejected = `109,565` total pairs—a `10,017`-column/pair discrepancy with no explanation.
  3. **Condition-selected `n = 1` sampling bias:** Both evaluations sample `n = 1` locus/pair (`chr1` only for `Pf` × `Pv`; `1` intraspecific pair for *Cannabis*) and condition on columns **already inside** reference `MAF` (`LASTZ` `HoxD55`) blocks. As `scripts/benchmarks/extract_fasta_from_maf.py:41-45, 78-85` shows for the project's MAF evaluation pipeline, extracting only MAF-aligned blocks and concatenating them strips out unaligned/intergenic sequence, inflating retention/recall and preventing measurement of false-positive alignments outside reference blocks.
- **Evidence:**
  - `grep -rn "Cherry_Pie\|418,663\|418465\|87\.7" src/ scripts/ docs/` (0 hits).
  - `manuscript/horizon_paper.md:65` (`24,748 / 119,582`) vs. `manuscript/scripts/generate_figures.py:218-219` (`n=88,597` + `n=20,968` = `109,565`).
  - `scripts/benchmarks/extract_fasta_from_maf.py:41-45, 78-85`.
- **Correction:** Remove or re-run and track the *P. falciparum* × *P. vivax* and *Cannabis sativa* evaluations in `scripts/benchmarks/` and `docs/BENCHMARK_RESULTS.md` on `HEAD` (`33bede8`), reconcile the `119,582` vs. `109,565` denominator discrepancy, and evaluate on full unmasked chromosomes rather than MAF-extracted block concatenations.

---

## F5 -- C6 `30` randomized tandem-duplication / overlapping-inversion trials do not exist in the codebase

- **Claim:** `C6` (`manuscript/horizon_paper.md:99`) — *"Across `30` randomized tandem-duplication and overlapping-inversion trials, e-graph saturation recovered `1.00` ancestral block sets and `1.00` exact `dup(X)` duplications compared to `0.00` under greedy one-to-one synteny filtering."*
- **Verdict:** `REFUTED`
- **Concrete flaws (sampling lens):**
  1. **`0` of `30` randomized trials exist in the repository:** Searching `src/egraph.rs`, `src/bin/horizon_resolve.rs`, `tests/`, and `scripts/` reveals **no** 30-trial randomized harness or greedy-synteny comparison anywhere in the codebase.
  2. **Actual e-graph sample is `2` deterministic unit tests + `1` synthetic GFA run:** `src/egraph.rs:463-518` contains only `2` hardcoded unit tests (`egraph_collapses_inverted_tandem_duplication` on 5 blocks `[A, -C, -B, -B, D]` and `mdl_cost_model_gates_dup_by_sequence_divergence` on 2 blocks `[B, B]`). `docs/BENCHMARK_RESULTS.md:67` reports only `1` synthetic `horizon_resolve --collapse-dups` run (`200`-seg repeat, `20` identical walks, `--divergence 0.02`), not `30` randomized trials and not a `0.00` greedy baseline.
- **Evidence:**
  - `grep -rn "30" src/bin/horizon_resolve.rs tests/` (0 hits).
  - `src/egraph.rs:463-518` (2 unit tests only).
  - `docs/BENCHMARK_RESULTS.md:67`.
- **Correction:** Either add and run the 30-trial randomized simulation + greedy one-to-one baseline script in `scripts/benchmarks/` (and log it in `docs/BENCHMARK_RESULTS.md`), or remove the unsubstantiated `30`-trial sentence from `manuscript/horizon_paper.md:99`.

---

## F6 -- Claim-list steering: `C1`–`C6` omit unrun/missing benchmarks, `horizon_merge` provenance gap, and external baseline absences

- **Framing finding (`BRIEF-sampling.md`):** The runner-authored claim list `C1`–`C6` steers attention away from four major sampling gaps in `manuscript/horizon_paper.md` and `docs/BENCHMARK_RESULTS.md`:
  1. **Zero external pangenome/WGA head-to-head baselines:** `round 2's R2-RESOLVE-LIFT` and `round 3's F7` flagged the complete absence of `minimap2`, `odgi`, `vg`, and `LASTZ` baselines. `docs/BENCHMARK_RESULTS.md:5` still concedes *"head-to-head `minimap2`/`odgi`/`vg`/`LASTZ` baselines remain pending external tool availability on this host"* (even though `/home/rcb112/project/gemini/bench/vg` is present on disk).
  2. **`horizon_merge` scaling claims (`manuscript/horizon_paper.md:10, 63, 75, Table 1`) have no benchmark row in `docs/BENCHMARK_RESULTS.md`:** `docs/BENCHMARK_RESULTS.md:9` explicitly states *"no multi-species benchmark row yet"*, yet `horizon_paper.md` reports 7 multi-species merge timings (`406` mammals in `3.84 s`, `317` species in `2.53 s`, `60_amniotes` in `16.61 s`, etc.).
  3. **Pending dataset:** `docs/BENCHMARK_RESULTS.md:49` (`B10K 363-way` `galGal6` × `taeGut2`) remains `*(Pending Data)*`.
  4. **MAF-block-only extraction bias in `run_all_evals.py`:** `scripts/benchmarks/extract_fasta_from_maf.py:41-45, 78-85` strips out all unaligned intergenic/repetitive sequence before `horizon_eval` runs, so `Table 1` and `docs/BENCHMARK_RESULTS.md:36-48` sample only pre-aligned MAF blocks rather than whole chromosomes/contigs.

---

## Verdict Table

| Claim ID | Verdict | Summary |
| :--- | :--- | :--- |
| `C1` | `UNCHALLENGED` | Architectural/inventory claim; no empirical sampling assertion to refute under the sampling lens. |
| `C2` | `NARROWED` | `CHM13.chr6` is q-arm only (`61.4M–172.1M`), `n = 1` chromosome and `n = 1` run, no accuracy/coverage metric, and pre-`33bede8`. |
| `C3` | `NARROWED` | `horizon_gfa extract` sampled only `2` chromosomes (`chr22`, `chr6`) and `2` haplotypes out of `2,580` paths / `24,456` walks; `chr6` pre-`33bede8`. |
| `C4` | `NARROWED` | `horizon_clean` has `0` sequence-bearing real GFA benchmarks (only a `0`-severed `*-seq` no-op and `2` toy 60-bp repeat unit tests). |
| `C5` | `REFUTED` | Untracked *Plasmodium*/*Cannabis* datasets, contradictory denominators (`119,582` vs. `109,565`), and condition-selected `n = 1` MAF-block bias. |
| `C6` | `REFUTED` | `30` randomized duplication/inversion trials and `0.00` greedy baseline do not exist anywhere in the codebase or `BENCHMARK_RESULTS.md`. |
