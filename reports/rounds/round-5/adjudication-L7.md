# Round 5 ADJUDICATION of reproduction -- Round 5: Full chr6 (8cd5db0/92e345c), HPRC & minimap2 Verification, B10K 363-way, and Technical Companion

Adjudicating `reports/rounds/round-5/attack-reproduction.md` (sha256 prefix `a97f95063960`, verified against `a97f9506396044bf90b049519be4ce7d6d8c3f3627499baee0392f451d98fce9`) against `manuscript/horizon_paper.md`, `manuscript/horizon_companion.md`, `docs/BENCHMARK_RESULTS.md`, and `horizon-aligner` `HEAD` (`92e345c`).

- **Prior round context (`scripts/round history "reproduction" --kind adjudication`):** Round 1 (`F3`), Round 3 (`F4`/`F5`), and Round 4 (`F3`) previously adjudicated `horizon_gfa stats` reproduction on the `12 GB` HPRC v1.0 GFA and flagged missing repository benchmark scripts in `scripts/benchmarks/`. Round 5 introduces new claims C1–C6 (`10k` HPRC subgraph vs `vg`, full `chr6` vs HPRC/`minimap2`, `horizon_merge`/`horizon_guide`, `B10K 363-way`, `BandGrid`, `# EG:Z:`).
- **Method cards (`refs/methods/`):** `refs/methods/` does not exist in this repository; F1–F6 are codebase, formula, and benchmark reproduction checks rather than statistical disputes governed by a method card.

---

## F1 -- C4 (`18%`-GC random DNA): `horizon_wga` emits false-positive pairs on 4 of 5 random seeds -- UPHELD

- **What survives:** The universal claim that `horizon_wga` eliminates all false-positive pairs on `10,000`-bp `18%`-GC random DNA (`0 pairs in 0.0069 s`) is refuted on `HEAD` (`92e345c`), while the four biological `horizon_eval` suites (`multiz100way chrM`, `cactus447way`, `B10K 363-way`, `cactus241way chr22 1MB`) reproduce as reported (with `B10K` and `cactus241way chr22 1MB` inputs residing only in `/tmp/` and absent from `scripts/benchmarks/run_all_evals.py`).
- **Verification:** Independent execution of `./target/release/horizon_wga <x.fa> <y.fa> --maf <out.maf>` on five `10,000`-bp i.i.d. `18%`-GC random pairs (`P(A)=P(T)=0.41`, `P(G)=P(C)=0.09`, `random.seed` in `{1, 7, 42, 123, 999}`) produced non-zero false-positive aligned pairs on `4/5` seeds (`seed=1`: `378` pairs; `seed=7`: `448` pairs; `seed=42`: `173` pairs; `seed=999`: `94` pairs; `seed=123`: `0` pairs), confirming the attacker's diagnosis that `0` pairs is a single-seed draw rather than a universal property.
- **Ruling on repair:** Endorsed on replacing the universal `0`-pair claim across `manuscript/horizon_paper.md:10,59,86` and `docs/BENCHMARK_RESULTS.md:13`, **narrowed/overturned** on hardcoding `0–189` pairs as the range (our independent 5-seed sample reached `448` pairs). The repair is **not free**: it requires running a larger multi-seed null sweep (e.g., `20+` seeds) and updating the manuscript text, Figure 6B, and `docs/BENCHMARK_RESULTS.md`.

## F2 -- C5 (`BandGrid` flat buffer index): wrong struct name, wrong row center, and wrong band width -- UPHELD

- **What survives:** `src/mdl/dp.rs:132-172,244-246` defines `BandLayout` (not `BandGrid`, which only appears as `DagBandGrid` in `src/mdl/dag.rs:11`), clamps `bw = band_w.min(n.max(m))`, centers row `i` at `floor(i * m / n)` (so `j = i` is not the center when `n != m`), and for the `6 bp` worked example (`n = m = 6`, `band_w = 180`) allocates `(6 + 1) * (2 * 6 + 1) = 91` floats with `idx(i, j) = i * 13 + (j - i + 6)` rather than `2,527` floats with `i * 361 + (j - i + 180)`.
- **Verification:** Direct inspection of `src/mdl/dp.rs:132-178,244-246` and `grep -rn "BandGrid\|BandLayout" src/` confirms every point of the attacker's diagnosis.
- **Ruling on repair:** Endorsed. Updating `BandGrid` to `BandLayout`, stating `bw = min(W, max(n, m))` and `c(i) = floor(i * m / n)`, and fixing the `6 bp` worked-example allocation to `7 × 13 = 91` floats in `manuscript/horizon_companion.md:19,122,139-145,147,149,178-179` and `manuscript/horizon_paper.md:22,36,80` requires no re-measurement: **the fix is free**.

## F3 -- C1 (`10k` HPRC `chr22` subgraph vs `vg`): metrics reproduce on `/tmp/` artifact, but no repo script generates or runs it -- UPHELD

- **What survives:** C1's counts and `vg v1.60.0` comparison reproduce on `/tmp/hprc_chr22_sub_real.gfa` (`horizon_gfa stats`: `segments=10000`, `segment_bp=10937960`, `links=10300`, `paths=1`, `walks=4`, `peak_rss_mb=22.7`, `elapsed=0.02s`; `/home/rcb112/project/gemini/bench/vg stats -z -l -N`: `nodes=10000`, `edges=10300`, `length=10937960` in `1.6867s`), but `/tmp/hprc_chr22_sub_real.gfa` is un-versioned and no script in `scripts/` extracts it or runs `vg`.
- **Ruling on repair:** Endorsed. Checking in an extraction and `vg` comparison script under `scripts/benchmarks/` is **not free**.

## F4 -- C2 (Full `chr6` vs HPRC graph and `minimap2`): summary artifacts match, but pipeline scripts and FASTAs live only in `/tmp/` -- UPHELD

- **What survives:** `manuscript/data/chr6_hprc_comparison.json` and `manuscript/data/chr6_minimap2_summary.json` are byte-for-byte identical (`diff -u` exit `0`) to `/tmp/chr6_comparison_summary.json` and `/tmp/chr6_minimap2_summary.json` and contain C2's exact reported metrics (`minimap2` `2.28-r1209` verified), while `/tmp/extract_hprc_chr6_truth.py`, `/tmp/compare_chr6_horizon_vs_hprc.py`, `/tmp/score_paf_vs_hprc.py`, `/tmp/hprc_grch38_chr6.fa`, `/tmp/hprc_chm13_chr6.fa`, and `/tmp/chr6_fixed.maf` exist only in `/tmp/` and are absent from `scripts/benchmarks/`.
- **Ruling on repair:** Endorsed. Committing the `chr6` truth-extraction, alignment, and scoring scripts under `scripts/benchmarks/` is **not free**.

## F5 -- C3 (`horizon_merge` and `horizon_guide`): `horizon_guide` reproduces exactly; `horizon_merge` 92-species split has no repo script -- UPHELD

- **What survives:** Running `./target/release/horizon_guide --progressive --hal --rbest <tmp>/out 1` on `hg38`, `Macaca_mulatta`, `HLcalPym1`, `Equus_caballus` (`../pg/data/benchmarks/fastas/`) and `Canis_lupus` (`/tmp/dog_zoo.fa`) reproduces `total_blocks=1109 multi_species_blocks=967 hal_bottom_segs=4` (`elapsed_s=0.9833`), and `horizon_wga --hal-info` confirms `hdf5_v=0 root_members=5`; however, `/tmp/dog_zoo.fa` lives only in `/tmp/` and no script in `scripts/` splits `../pg/data/benchmarks/chrM.maf` into `91` pairwise MAFs for `horizon_merge` or runs the `5`-mammal `horizon_guide` build.
- **Ruling on repair:** Endorsed. Adding reproducible scripts under `scripts/benchmarks/` for the `92`-species `chrM` `horizon_merge` split/run and the `5`-mammal `horizon_guide` run is **not free**.

## F6 -- C6 (`# EG:Z:` tags from `src/egraph.rs`): implemented and unit-tested -- UPHELD

- **What survives:** C6 survives as claimed: `src/formats.rs:153-158` (`format_gfa_paralog_tags`) and `src/merge.rs:358-363` emit `# EG:Z:<ref>__<qry>\t...` lines (and `src/bin/horizon_resolve.rs:476` emits `EG:Z:` tags on `W` lines), verified by passing unit tests under `cargo test --all-features` (`97` active passed + `1` ignored = `98` total; note `manuscript/horizon_paper.md:82` states `98 active` instead of `97 active`).
- **Ruling on repair:** No repair is needed for C6.

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `horizon_wga` emits false-positive pairs on `4/5` `18%`-GC random seeds (`94–448` in our run, `80–189` in attack), refuting the universal `0`-pair claim; multi-seed re-measurement is not free |
| F2 | UPHELD | `src/mdl/dp.rs` defines `BandLayout` (not `BandGrid`) with `center(i)=floor(i*m/n)` and `bw=min(W, max(n, m))` (`91` floats for `6 bp`, not `2,527`); the text/formula fix is free |
| F3 | UPHELD | `horizon_gfa stats` and `vg v1.60.0` reproduce on `/tmp/hprc_chr22_sub_real.gfa`, but no repo script generates it or runs `vg`; adding the script is not free |
| F4 | UPHELD | `manuscript/data/chr6_*.json` match `/tmp/` artifacts and C2 numbers (`minimap2 2.28-r1209`), but scripts and FASTAs exist only in `/tmp/`; committing scripts is not free |
| F5 | UPHELD | `horizon_guide --progressive --hal --rbest` (`1,109` blocks, `967` multi-sp, `4` bottom segs) reproduces, while `92`-species `horizon_merge` split and `5`-mammal run have no repo scripts; adding them is not free |
| F6 | UPHELD | `# EG:Z:` emission from `src/egraph.rs` is implemented in `src/formats.rs:157` / `src/merge.rs:361` and passes `cargo test --all-features` (`97` passed, `1` ignored); C6 survives |
