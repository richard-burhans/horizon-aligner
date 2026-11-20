# Round 5 Attack -- Reproduction Lens (L7)

- **Target artifacts:** `reports/rounds/round-5/ATTACKED-horizon_paper.md`, `reports/rounds/round-5/ATTACKED-horizon_companion.md`, `docs/BENCHMARK_RESULTS.md`, and `horizon-aligner` `HEAD` (`92e345c`).
- **Build & test baseline on `HEAD` (`92e345c`):**
  - `cargo test --all-features`: exited `0` (`97` active unit tests passed across `src/lib.rs` [`70` passed, `1` ignored] and `src/bin/*` [`27` passed]; note `ATTACKED-horizon_paper.md:82` states `98 active Rust unit tests (plus 1 environment-gated .hal integration test)`, whereas the actual count is `97` active + `1` ignored = `98` total).
  - `cargo clippy --all-targets --all-features -- -D warnings`: exited `0` with zero warnings.
- **Literature & method cards (`BRIEF-reproduction.md` Rules 2–3):** No `../library` checkout was needed for codebase/benchmark reproduction checks, and `refs/methods/` does not exist in this repository (`BRIEF-reproduction.md:142`).
- **Prior round context (`scripts/round history "reproduction" --kind adjudication`):** Round 1 (`F3`), Round 3 (`F4`/`F5`), and Round 4 (`F3`) previously adjudicated reproduction of `horizon_gfa stats` on the `12 GB` HPRC v1.0 GFA and flagged missing repository benchmark scripts (`run_all_evals.py`). Round 5 introduces new claims C1–C6 (`10k` HPRC subgraph vs `vg`, full `chr6` vs HPRC/`minimap2`, `horizon_merge`/`horizon_guide`, `B10K 363-way`, `BandGrid`, `# EG:Z:`).

---

## F1 -- C4 (`18%`-GC random DNA): `horizon_wga` produces 80–189 false-positive pairs on 4 of 5 random seeds (REFUTED)

- **Claim:** C4 (`ATTACKED-horizon_paper.md:10,59`): `"while eliminating all false-positive pairs on 18%-GC random DNA (0 pairs in 0.0069 s)"`.
- **Evidence:**
  1. Direct execution of `./target/release/horizon_eval` on `HEAD` (`92e345c`) reproduces the four biological benchmarks in C4:
     - `multiz100way chrM` (`hg38` × `panTro4`): `raw_maf_recall=95.6% raw_maf_precision=97.4% bubbles=103` (`elapsed=0.068s`).
     - `cactus447way` (`hg38` × `Macaca_mulatta`): `raw_maf_recall=91.4% raw_maf_precision=87.0% bubbles=267` (`elapsed=0.117s`).
     - `B10K 363-way` (`/tmp/b10k_galGal.fa`, `/tmp/b10k_taeGut.fa`, `/tmp/b10k_galGal_taeGut.maf`): `raw_maf_recall=71.1% raw_maf_precision=95.1% bubbles=3814` (`elapsed=0.564s`).
     - `cactus241way chr22 1MB` (`/tmp/repro_hg38_zoo_1m.fa`, `/tmp/repro_dog_zoo_1m.fa` [31 contigs], `/tmp/repro_pw_zoo_1m.maf`): `raw_maf_recall=59.8% raw_maf_precision=78.5% bubbles=2076` (`elapsed=0.914s` vs `0.852 s`).
  2. However, testing `./target/release/horizon_wga` on five independent `10,000`-bp `18%`-GC random DNA pairs (`random.seed` in `{1, 7, 42, 123, 999}`, i.i.d. `P(G)=P(C)=0.09`, `P(A)=P(T)=0.41`) yields:
     - `seed=1`: `aligned_pairs=145 unique_x_bp=145 unique_y_bp=145 total_bubbles=6 elapsed_s=0.0165`
     - `seed=7`: `aligned_pairs=80 unique_x_bp=80 unique_y_bp=80 total_bubbles=6 elapsed_s=0.0241`
     - `seed=42`: `aligned_pairs=171 unique_x_bp=163 unique_y_bp=161 total_bubbles=9 elapsed_s=0.0179`
     - `seed=123`: `aligned_pairs=189 unique_x_bp=180 unique_y_bp=189 total_bubbles=11 elapsed_s=0.0295`
     - `seed=999`: `pruned_pairs=1 aligned_pairs=0 elapsed_s=0.0045`
  3. Additionally, `B10K 363-way` and `cactus241way chr22 1MB` inputs reside only in `/tmp/` and are absent from `scripts/benchmarks/run_all_evals.py`.
- **Correction:** Replace `"eliminating all false-positive pairs on 18%-GC random DNA (0 pairs in 0.0069 s)"` in `manuscript/horizon_paper.md:10,59` with the observed multi-seed false-positive distribution (e.g., `0–189` pairs across `10,000`-bp `18%`-GC random pairs) rather than claiming `0` false-positive pairs universally.

## F2 -- C5 (`BandGrid` flat buffer index): wrong struct name, wrong row center, and wrong band width (REFUTED)

- **Claim:** C5 (`ATTACKED-horizon_companion.md:139-145`): `"For target index i and query index j within |j − i| ≤ W, the flat buffer index is: idx(i, j) = i · (2W + 1) + (j − i + W)"` (and `ATTACKED-horizon_companion.md:179`: `"For the 6 bp gap (band_w = 180), BandGrid allocates (6 + 1) × (2 · 180 + 1) floats and evaluates idx(i, j) = i · 361 + (j − i + 180)"`).
- **Evidence:**
  1. `grep -rn "BandGrid" src/` returns only `DagBandGrid` in `src/mdl/dag.rs:11`; `src/mdl/dp.rs:132` defines `BandLayout`, not `BandGrid`.
  2. In `src/mdl/dp.rs:141-172` and `244-246`, `BandLayout::new(n, m, band_w)` sets `bw = band_w.min(n.max(m))` and `stride = 2 * bw + 1`, and `BandLayout::slot(i, j)` centers row `i` at `band_center(i, m, n) = (i * m) / n.max(1)`:
     - The actual index formula is `i * (2 * bw + 1) + (j - floor(i * m / n) + bw)` for `|j - floor(i * m / n)| <= bw`.
     - When `n != m`, the main diagonal `j = i` is not the row center.
     - When `max(n, m) < W` (such as the `6 bp` worked example with `band_w = 180`), `bw = min(180, 6) = 6`, so `BandLayout` allocates `(6 + 1) * (2 * 6 + 1) = 91` floats, not `7 * 361 = 2,527` floats.
- **Correction:** Update `manuscript/horizon_companion.md:139-145,179` and `manuscript/horizon_paper.md:36` to name `BandLayout` (`src/mdl/dp.rs:132`), state the actual center `floor(i * m / n)` and clamped half-width `bw = min(W, max(n, m))`, and fix the worked-example allocation to `7 × 13 = 91` floats.

## F3 -- C1 (`10k` HPRC `chr22` subgraph vs `vg`): metrics reproduce on `/tmp/` artifact, but no repo script generates or runs it (NARROWED)

- **Claim:** C1 (`ATTACKED-horizon_paper.md:82`).
- **Evidence:**
  - Running `./target/release/horizon_gfa stats /tmp/hprc_chr22_sub_real.gfa` reproduces the exact counts (`segments=10000 segment_bp=10937960 links=10300 paths=1 walks=4`), `peak_rss_mb=22.6`, and `elapsed=0.03s` (`0.0173s` wall in Python).
  - Running `/home/rcb112/project/gemini/bench/vg` (`vg version v1.60.0 "Annicco"`) on `/tmp/hprc_chr22_sub_real.gfa` reproduces `nodes 10000`, `edges 10300`, `length 10937960` in `1.6888s` (`vg stats -z -l -N`) and `1.6836s` (`vg convert -g -p`).
  - However, `/tmp/hprc_chr22_sub_real.gfa` is an un-versioned `/tmp/` file, and no script in `scripts/` extracts this `10,000`-segment subgraph from `hprc-v1.0-mc-grch38.gfa.gz` or invokes `vg`.
- **Correction:** Check in the subgraph extraction and `vg` comparison script under `scripts/benchmarks/` so C1 is reproducible from the repository alone.

## F4 -- C2 (Full `chr6` vs HPRC graph and `minimap2`): summary artifacts match, but pipeline scripts and FASTAs live only in `/tmp/` (NARROWED)

- **Claim:** C2 (`ATTACKED-horizon_paper.md:76`).
- **Evidence:**
  - `manuscript/data/chr6_hprc_comparison.json` and `manuscript/data/chr6_minimap2_summary.json` match `/tmp/chr6_comparison_summary.json` and `/tmp/chr6_minimap2_summary.json` and contain the exact reported metrics (`hprc_1to1_segs=1,628,018`, `hprc_1to1_x_bp=110,082,990`, `horizon_blocks=370,937`, `horizon_total_pairs=282,163,046`, `x_base_coverage_recall_of_hprc=0.99539`, `exact_coord_recall_on_hprc_1to1=0.97959`, `within_1kb_recall_on_hprc_1to1=0.98144`, `exact_coord_precision_on_shared=0.48633`; `minimap2 v2.28-r1209`: `exact_coord_recall_on_hprc_1to1=0.99706`, `exact_coord_precision_on_shared=0.99751`). `/tmp/chr6_fixed.maf` contains `370,937` `a ` blocks and `/tmp/minimap2-2.28_x64-linux/minimap2 --version` reports `2.28-r1209`.
  - However, the truth-extraction and scoring scripts (`/tmp/extract_hprc_chr6_truth.py`, `/tmp/compare_chr6_horizon_vs_hprc.py`, `/tmp/score_paf_vs_hprc.py`) and input FASTAs (`/tmp/hprc_grch38_chr6.fa`, `/tmp/hprc_chm13_chr6.fa`) exist only in `/tmp/` and are not committed in `scripts/benchmarks/`.
- **Correction:** Commit the HPRC `chr6` truth-extraction, `horizon_wga`, `minimap2`, and scoring scripts under `scripts/benchmarks/`.

## F5 -- C3 (`horizon_merge` and `horizon_guide`): `horizon_guide` reproduces exactly; `horizon_merge` 92-species split has no repo script (NARROWED)

- **Claim:** C3 (`ATTACKED-horizon_paper.md:82`).
- **Evidence:**
  - Re-running `./target/release/horizon_guide --progressive --hal --rbest <tmp>/out 1` on `hg38`, `Macaca_mulatta`, `HLcalPym1`, `Equus_caballus` (from `pg/data/benchmarks/fastas/`) and `Canis_lupus` (`/tmp/dog_zoo.fa`) reproduces `total_blocks=1109 multi_species_blocks=967 hal_bottom_segs=4` in `1.011s`, and `horizon_wga --hal-info` confirms `hdf5_v=0 root_members=5`. Note that `/tmp/dog_zoo.fa` is only in `/tmp/` and the leftover `/tmp/guide_5sp.maf` on disk has `1,225` blocks from an earlier run.
  - For `horizon_merge`, `pg/data/benchmarks/chrM.maf` contains `92` unique species, but `horizon_merge` takes per-species pairwise MAFs (`<sp=path.maf>...`) and no script in `scripts/` splits `chrM.maf` into `91` pairwise MAFs to reproduce the `21.78 s` / `26.6 MB` / `460`-block run.
- **Correction:** Add a reproducible script under `scripts/benchmarks/` for both the `92`-species `chrM` `horizon_merge` split/run and the `5`-mammal `horizon_guide` run.

## F6 -- C6 (`# EG:Z:` tags from `src/egraph.rs`): implemented and unit-tested (SURVIVED)

- **Claim:** C6 (`ATTACKED-horizon_companion.md:52`).
- **Evidence:**
  - `src/egraph.rs` (`collect_windowed_duplications`), `src/formats.rs:153-158` (`format_gfa_paralog_tags`), and `src/merge.rs:358-363` emit `# EG:Z:<ref>__<qry>\t...` lines (while `src/bin/horizon_resolve.rs:476` emits `EG:Z:` optional tags on `W` lines).
  - Verified by passing unit tests `src/lib.rs:332` (`assert!(merged_gfa.contains("# EG:Z:Pf__Pr\tB0"))`) and `src/bin/horizon_resolve.rs:617-618` under `cargo test --all-features`, as well as `/tmp/chr6_8cd5db0.log`.
- **Correction:** None required for C6 itself (note that no downstream Horizon reader currently parses `# EG:Z:` comment lines).

---

## Claim Verdicts

| Claim | Verdict | Summary |
| :--- | :--- | :--- |
| C1 | NARROWED | `horizon_gfa stats` (`10k` segs, `22.6 MB` RSS) and `vg v1.60.0` speedup reproduce on `/tmp/hprc_chr22_sub_real.gfa`, but no repo script generates or runs it. |
| C2 | NARROWED | JSON/MAF/`minimap2 v2.28-r1209` artifacts match reported numbers, but extraction/scoring scripts and FASTAs exist only in `/tmp/`. |
| C3 | NARROWED | `horizon_guide --progressive --hal --rbest` (`1,109` blocks, `967` multi-sp, `4` bottom segs) reproduces; `92`-species `horizon_merge` has no split/run script in repo. |
| C4 | REFUTED | Biological `horizon_eval` suites reproduce, but `18%`-GC random DNA yields `80–189` false-positive pairs on `4/5` seeds (`0` only on `seed=999`). |
| C5 | REFUTED | `BandGrid` does not exist in `src/mdl/dp.rs` (`BandLayout` does); row center is `floor(i·m/n)` and half-width is `min(W, max(n, m))` (`91` floats for `6 bp`, not `2,527`). |
| C6 | SURVIVED | `# EG:Z:` emission from `src/egraph.rs` is implemented in `src/formats.rs:157` / `src/merge.rs:361` and passes `cargo test --all-features`. |
