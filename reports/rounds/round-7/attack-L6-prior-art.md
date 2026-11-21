# Round 7 Adversarial Attack -- L6 (prior-art)

## Context & Constraints (Rules 2, 3, 4)

- **Literature access (Rule 2):** `../library` is not present in this workspace. Accordingly, all prior-art checks were conducted strictly from [`refs/catalog_index.json`](../../../refs/catalog_index.json), [`manuscript/references.bib`](../../../manuscript/references.bib), [`reports/rounds/round-7/ATTACKED-horizon_paper.md`](ATTACKED-horizon_paper.md), [`reports/rounds/round-7/ATTACKED-BENCHMARK_RESULTS.md`](ATTACKED-BENCHMARK_RESULTS.md), and tracked source/benchmark files; no full-text claims are asserted from memory.
- **Method cards (Rule 3):** This project keeps no method cards (`refs/methods/` does not exist). Methodological/statistical prior-art claims (specifically `C3`) are checked against [`refs/catalog_index.json`](../../../refs/catalog_index.json) and the Rust implementation in [`src/mdl/`](../../../src/mdl/).
- **Prior-round adjudications (Rule 4):** Checked `scripts/round history` and [`reports/rounds/round-6/adjudication-L6.md`](../round-6/adjudication-L6.md) (`F1`–`F6`, all `NARROWED`). Commits `4fd2240` and `5cb83fb` discharged several Round 6 L6 items: the 9 held-catalog citations (`li_2026_finding`, `frith_2011_tantan`, `garrison_2024_building`, `andreace_2023_comparing`, `romain_2026_inversions-pangenome-graphs`, `kazemi_2022_nthash2`, `marcais_2024_knonical`, `wittler_2023_canonical-kmer-encoding`, `chen_2010_comparative`) plus `jukes_1969_evolution` were added to `manuscript/references.bib` and `manuscript/horizon_paper.md`, `foster2007combinators` was removed from `manuscript/references.bib`, and `minimap2`/`FastGA`/`odgi`/`gfatools`/`vg` baseline scripts were tracked in `scripts/benchmarks/`. Unresolved prior-art gaps and new defects are detailed below.

---

## F1 -- C3: Mischaracterization of CRC as offline/fixed (`λ̂ = 10.0`) and invalid exchangeability premise under Angelopoulos et al. (2024)

- **Claim:** `C3` ([`ATTACKED-horizon_paper.md:40`](ATTACKED-horizon_paper.md#L40); [`round.json:28-36`](round.json#L28-L36))
- **Verdict:** **REFUTED**
- **Concrete flaw:**
  1. **Draft contradicts code (`C3` summary vs. `ATTACKED-horizon_paper.md:40`):** `ATTACKED-horizon_paper.md:40` still states: *"Using fixed Conformal Risk Control defaults calibrated offline at risk level `α = 0.07` (`n_cal = 15`; Angelopoulos et al. 2024), Horizon enforces a constant bit-savings floor `λ̂ = 10.0` bits (`min_delta_bits = 10.0`)..."*, whereas `round.json:31` claims dynamic CRC is wired in `src/mdl/mod.rs` (confirmed at [`src/mdl/mod.rs:401-435`](../../../src/mdl/mod.rs#L401-L435)).
  2. **Invalid attribution to Angelopoulos et al. (2024) (`angelopoulos_2022_conformal`, `doi:10.48550/arXiv.2208.02814`):** Even in the dynamic implementation ([`src/mdl/dp.rs:82-165`](../../../src/mdl/dp.rs#L82-L165), `calibrate_null_crc_for_gc`), `λ̂` is calibrated on `n_cal = 15` **synthetic 64-bp i.i.d. PRNG decoy pairs** (`0x9E37_79B9_7F4A_7C15_u64`), not an exchangeable calibration sample from the target distribution of real genomic intervals. Conformal Risk Control requires exchangeability between calibration and test points; i.i.d. 64-bp Bernoulli draws are not exchangeable with real chromosomes containing repeats, composition shifts, and `>>64`-bp gaps (as shown by [`ATTACKED-horizon_paper.md:61`](ATTACKED-horizon_paper.md#L61), where `4/5` other `18%`-GC `10`-kb random seeds still emit `80–448` false-positive pairs and `30%`/`45%`-GC draws emit `288`/`149` pairs). Moreover, at `n_cal = 15` and `α = 0.07`, `(n_ge + 1) / (n_cal + 1) <= 0.07` holds only for `n_ge = 0` (`1/16 = 0.0625 <= 0.07`; `2/16 = 0.125 > 0.07`), reducing `calibrate_crc_threshold` to taking the maximum of 15 synthetic 64-bp decoy scores.
- **Evidence:** [`ATTACKED-horizon_paper.md:40, 61`](ATTACKED-horizon_paper.md#L40); [`src/mdl/dp.rs:82-165`](../../../src/mdl/dp.rs#L82-L165); [`src/mdl/mod.rs:401-435`](../../../src/mdl/mod.rs#L401-L435); [`refs/catalog_index.json:130-133`](../../../refs/catalog_index.json#L130-L133).
- **Correction:** Update [`ATTACKED-horizon_paper.md:40`](ATTACKED-horizon_paper.md#L40) to match the runtime behavior and replace the Conformal Risk Control (Angelopoulos et al. 2024) guarantee claim with an accurate description of the heuristic synthetic-decoy max-floor calibration (`n_cal = 15` 64-bp i.i.d. draws), unless calibrated on held-out exchangeable biological null intervals.

---

## F2 -- C2: Missing head-to-head `SweepGA` / PAF-netter baseline and `SweepGA` / `wfmash` / `impg` citations (carrying forward Round 6 L6 `F2`, `F3`)

- **Claim:** `C2` ([`ATTACKED-BENCHMARK_RESULTS.md:85-88`](ATTACKED-BENCHMARK_RESULTS.md#L85-L88); [`round.json:18-26`](round.json#L18-L26))
- **Verdict:** **NARROWED**
- **Prior-round context:** Round 6 L6 `F2` and `F3` ([`reports/rounds/round-6/adjudication-L6.md:17-25`](../round-6/adjudication-L6.md#L17-L25)) were adjudicated `NARROWED` because `horizon_filter --net` lacked a head-to-head `SweepGA` baseline and `horizon_locus` lacked an `impg` baseline and citations. In Round 7 (`4fd2240`), `horizon_paper.md:98` added a roadmap note that `SweepGA` and `KegAlign` comparisons remain future work, and `scripts/benchmarks/bench_external_graph_baselines.py` / `bench_external_wga_baselines.py` added `gfatools`/`odgi`/`vg` and `FastGA`/`minimap2`/`LASTZ` runs.
- **Concrete flaw:**
  1. `C2` (`ATTACKED-BENCHMARK_RESULTS.md:85-86`) still reports single-filter `horizon_filter --net` accuracy (`99.12%` exact recall, `99.65%` cond. / `99.32%` uncond. precision) on `FastGA` PAFs with **no head-to-head `SweepGA` (or `wfmash` / `minimap2` PAF-netting) baseline** on the identical `56,857`-record `chr6` input (`scripts/benchmarks/score_paf_vs_hprc.py` scores only `horizon_filter`).
  2. `SweepGA` (`ATTACKED-horizon_paper.md:98`), `wfmash` (`line 52`), and `impg` (`lines 10, 53`) remain **uncited** in `ATTACKED-horizon_paper.md` `# References` and `manuscript/references.bib` (none is in `refs/catalog_index.json`, as Round 6 L6 `F1`/`F3` already noted).
  3. In `bench_external_wga_baselines.py` (`5 Mbp` `chr6` pericentromere), `horizon_wga` is not compared against `KegAlign` (Gulhan et al. 2025, `doi:10.1186/s13059-025-03830-0`), and in `bench_external_graph_baselines.py` (`10k–200k` HPRC subgraphs), `horizon_resolve` and `horizon_wga --graph` have no external snarl-caller (`vg snarls`) or sequence-to-graph aligner (`vg giraffe` / `GraphAligner`) baseline.
- **Evidence:** [`ATTACKED-BENCHMARK_RESULTS.md:85-88`](ATTACKED-BENCHMARK_RESULTS.md#L85-L88); [`ATTACKED-horizon_paper.md:10, 52-53, 84, 98`](ATTACKED-horizon_paper.md#L84); [`scripts/benchmarks/bench_external_graph_baselines.py`](../../../scripts/benchmarks/bench_external_graph_baselines.py); [`scripts/benchmarks/bench_external_wga_baselines.py`](../../../scripts/benchmarks/bench_external_wga_baselines.py).
- **Correction:** Qualify `C2` in `docs/BENCHMARK_RESULTS.md:85-86` as a single-filter evaluation pending a head-to-head `SweepGA` baseline, and add bibliographic entries for `SweepGA`, `impg`, and `wfmash` to `manuscript/references.bib` and `manuscript/horizon_paper.md`.

---

## F3 -- C4: Progressive HAL/MAF benchmark lacks external `HAL` validation (`halValidate`/`hal2maf`), progressive WGA baselines, and `paten_2011_cactus` citation; stale `4`-segment stub persists in Abstract/Summary

- **Claim:** `C4` ([`ATTACKED-horizon_paper.md:84`](ATTACKED-horizon_paper.md#L84); [`round.json:38-46`](round.json#L38-L46))
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. `verify_merge_guide_maf.py:139-156` checks the emitted `.hal` only with Horizon's own `horizon_wga --hal-info` and scores only the `.maf` output, without external `HAL` toolkit verification (`halValidate` / `hal2maf`, Hickey et al. 2013, `doi:10.1093/bioinformatics/btt128`) or any head-to-head progressive WGA baseline (`Progressive Cactus`, `progressiveMauve`, `MULTIZ/TBA`, all cited in `references.bib`).
  2. `ATTACKED-horizon_paper.md:10` (Abstract) and `ATTACKED-BENCHMARK_RESULTS.md:11, 58` still report the superseded `"single-contig star-tree HDF5 v0 HAL stub (4 bottom segments)"` rather than the `2,103`-segment progressive ancestral-node HAL claimed in `C4` (`ATTACKED-horizon_paper.md:84`).
  3. `ATTACKED-horizon_paper.md` and `manuscript/references.bib` omit `paten_2011_cactus` (`doi:10.1101/gr.123356.111`) and `angiuoli_2011_mugsy` (`doi:10.1093/bioinformatics/btq665`) from `refs/catalog_index.json`, both foundational to progressive/Cactus multiple-genome alignment.
- **Evidence:** [`scripts/benchmarks/verify_merge_guide_maf.py:139-156`](../../../scripts/benchmarks/verify_merge_guide_maf.py#L139-L156); [`ATTACKED-horizon_paper.md:10, 84`](ATTACKED-horizon_paper.md#L10); [`ATTACKED-BENCHMARK_RESULTS.md:11, 58`](ATTACKED-BENCHMARK_RESULTS.md#L11); [`refs/catalog_index.json:135-138, 1110-1113`](../../../refs/catalog_index.json#L1110-L1113).
- **Correction:** Update `horizon_paper.md:10` and `docs/BENCHMARK_RESULTS.md:11, 58` to match the `2,103`-segment progressive HAL, qualify `C4` as lacking external `halValidate`/`hal2maf` and progressive-WGA baseline comparisons, and cite `paten_2011_cactus` (`doi:10.1101/gr.123356.111`).

---

## F4 -- C5 (`R5-BIB`): Primary 9-paper sync landed, but `impg`/`wfmash`/`SweepGA`, held-catalog prior art, and DOI hygiene remain incomplete

- **Claim:** `C5` ([`ATTACKED-horizon_paper.md:150`](ATTACKED-horizon_paper.md#L150); [`round.json:48-56`](round.json#L48-L56))
- **Verdict:** **NARROWED**
- **Prior-round context:** Round 6 L6 `F5` (`reports/rounds/round-6/adjudication-L6.md:32-35`, `NARROWED`) flagged orphan `foster2007combinators` in `manuscript/references.bib` and the unapplied Round 5 L6 citations. In `4fd2240`, `foster2007combinators` was removed and all 9 named `R5-BIB` catalog entries plus `jukes_1969_evolution` were synced (`49` entries in `references.bib` matching `49` in `ATTACKED-horizon_paper.md` `# References`). `wittler_2023_canonical-kmer-encoding` (`doi:10.24072/pcjournal.323`) matches `refs/catalog_index.json:1430-1433`.
- **Concrete flaw:**
  1. `impg`, `wfmash`, and `SweepGA` are explicitly named in `ATTACKED-horizon_paper.md:10, 52, 53, 98` (and `impg`/`wfmash` were endorsed in Round 6 L6 `F1`/`F3`) but remain absent from `manuscript/references.bib` and `# References`.
  2. Five `.bib` entries (`adams2007bayesian`, `angelopoulos2024conformal`, `chiaromonte2002scoring`, `harris2007improved`, `jukes_1969_evolution`) omit `doi` fields despite DOIs existing in `refs/catalog_index.json` for `adams_2007_bayesian` (`10.48550/arXiv.0710.3742`), `angelopoulos_2022_conformal` (`10.48550/arXiv.2208.02814`), and `chiaromonte_2001_scoring` (`10.1142/9789812799623_0012`); `frith_2011_tantan` lists `year = {2011}` in `references.bib:417` vs. `2010` in `refs/catalog_index.json:456`.
  3. Directly relevant held-catalog papers in `refs/catalog_index.json` remain uncited in `horizon_paper.md` and `references.bib`, including `paten_2011_cactus` (`10.1101/gr.123356.111`), `shumate_2021_liftoff` (`10.1093/bioinformatics/btaa1016`, liftover), `shivakumar_2025_mumemto` (`10.1186/s13059-025-03644-0`, multi-genome pangenome matching), `marcosola_2021_wfa` (`10.1093/bioinformatics/btaa777`), and `grootkoerkamp_2024_astarpa` (`10.1093/bioinformatics/btae032`).
- **Evidence:** [`manuscript/references.bib:1-489`](../../../manuscript/references.bib); [`ATTACKED-horizon_paper.md:10, 52-53, 98, 100-150`](ATTACKED-horizon_paper.md#L100); [`refs/catalog_index.json`](../../../refs/catalog_index.json).
- **Correction:** Add `impg`, `wfmash`, `SweepGA`, `paten_2011_cactus`, `shumate_2021_liftoff`, and `shivakumar_2025_mumemto` to `manuscript/references.bib` and `horizon_paper.md`, and populate missing DOIs from `refs/catalog_index.json`.

---

## F5 -- C1: Test/clippy inventory has no prior-art defect in the quoted claim

- **Claim:** `C1` ([`ATTACKED-horizon_paper.md:84`](ATTACKED-horizon_paper.md#L84); [`round.json:8-15`](round.json#L8-L15))
- **Verdict:** **UNCHALLENGED**
- **Concrete flaw:** None within the prior-art lens; `C1` is a Rust unit-test and `clippy` pass-count inventory (`122` unit tests + `1` `.hal` integration test) with no novelty or literature assertion. (The missing `impg`/`wfmash` citations for `horizon_locus`/`horizon_filter` are addressed under `F2` and `F4`.)

---

## Verdict Table

| claim_id | verdict | summary |
|---|---|---|
| C1 | UNCHALLENGED | Test/clippy inventory (`122 + 1`) makes no prior-art or novelty claim (`F5`). |
| C2 | NARROWED | `horizon_filter` lacks a head-to-head `SweepGA` (or PAF-netter) baseline, `SweepGA`/`wfmash`/`impg` lack citations, and `KegAlign`/`vg snarls` baselines are absent (`F2`). |
| C3 | REFUTED | `ATTACKED-horizon_paper.md:40` still claims offline/fixed `λ̂ = 10.0` (contradicting `src/mdl/mod.rs:401`), and `calibrate_null_crc_for_gc` (`src/mdl/dp.rs:120-165`) calibrates on `n_cal = 15` synthetic 64-bp i.i.d. PRNG decoys, violating Angelopoulos et al. (2024) exchangeability and degenerating to `max` (`F1`). |
| C4 | NARROWED | Progressive HAL (`2,103` bottom segments) lacks external `halValidate`/`hal2maf` and progressive-WGA baselines, `ATTACKED-horizon_paper.md:10` / `BENCHMARK_RESULTS.md:11,58` still cite the `4`-segment stub, and `paten_2011_cactus` is uncited (`F3`). |
| C5 | NARROWED | 9-paper `R5-BIB` sync and `foster2007combinators` removal landed, but `impg`/`wfmash`/`SweepGA`, held-catalog prior art (`paten_2011_cactus`, `shumate_2021_liftoff`, `shivakumar_2025_mumemto`), and 3 catalog DOIs remain missing (`F4`). |
