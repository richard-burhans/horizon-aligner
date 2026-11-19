# Round 4 Attack — Reproduction (L7)

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `33bede8`, `2026-11-19`)
**Claims attacked:** `C1`–`C6` from `manuscript/horizon_paper.md`
**Literature & Method Cards (Rules 2–3):** `../library` does not exist on this host, so any literature context is from `refs/catalog_index.json` alone; per `docs/round_rubric.md`, this project keeps no method cards yet.
**Prior Round History (Rule 4):** Checked `scripts/round history` across closed Rounds 1–3. Round 3 (`adjudication-L7.md:38`) previously upheld `horizon_gfa stats` reproduction on `488b020` (Round 3 `C6`, now `C3`), and `adjudication-L4.md:19,35,43` upheld missing runner scripts/logs for several `docs/BENCHMARK_RESULTS.md` rows. `C2` (`451.35 s` chr6), `C5` (*Plasmodium* / *Cannabis*), and `C6` (`30` randomized e-graph trials) have no prior adjudications.

---

## F1 -- C1: 10 binaries in `src/bin/`, partial `stdin` composition, and 2 binaries outside the MDL/BOCPD core

- **Claim:** `C1` (`manuscript/horizon_paper.md:18`) — *"Horizon addresses these coupled bottlenecks in a unified, pure-Rust suite of nine composable executables sharing one information-theoretic core."*
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. **Executable count:** `src/bin/` contains **10** binary targets (`bench_chr1.rs`, `horizon.rs`, `horizon_clean.rs`, `horizon_eval.rs`, `horizon_gfa.rs`, `horizon_guide.rs`, `horizon_lift.rs`, `horizon_merge.rs`, `horizon_resolve.rs`, `horizon_wga.rs`), matching `manuscript/horizon_paper.md:101` ("all `10` CLI binaries"). `horizon` (`src/bin/horizon.rs:7-53`) is a 10th dispatcher binary invoking the 9 subcommands.
  2. **Partial pipe composability:** Only `horizon_gfa`, `horizon_lift`, `horizon_clean`, and `horizon_resolve` support `stdin` (`-`) / `stdout` streaming (`src/gfa/mod.rs:56,179`, `src/bin/horizon_resolve.rs:349-350`), whereas `horizon_wga`, `horizon_merge`, `horizon_guide`, `horizon_eval`, and `bench_chr1` require file paths.
  3. **Core sharing:** `horizon_gfa` (`src/bin/horizon_gfa.rs:4`) and `horizon_lift` (`src/bin/horizon_lift.rs:4-5`) import only `gfa` and `coord`, not the BOCPD/MDL/e-graph modules.
- **Evidence:**
  - `cargo test` at `33bede8`: `45` passed + `1` ignored (`real_hal_decode_from_env`) in `src/lib.rs`, plus `14` passed across all `10` binaries (`59` active + `1` ignored).
  - `cargo clippy --all-targets`: `0` warnings (`Finished dev profile` exit `0`).
  - `cargo build --release`: exit `0`.
  - `grep -rn "io::stdin\|\"-\"" src/bin/ src/gfa/`: `stdin` handling appears only in `src/gfa/mod.rs:56,179` and `src/bin/horizon_resolve.rs:350`.
- **Correction:** Change `manuscript/horizon_paper.md:18` from *"nine composable executables sharing one information-theoretic core"* to *"ten pure-Rust executables (a unified `horizon` dispatcher and nine specialized subcommands) built on a shared workspace crate, with streaming `stdin`/`stdout` composition across the GFA pangenome tools"*.

---

## F2 -- C2: Full `chr6` `451.35 s` speedup is pinned to stale commit `53130b8`, lacks a committed runner script, and uses a `110.7 Mbp` q-arm `CHM13.chr6` extract

- **Claim:** `C2` (`manuscript/horizon_paper.md:101`) — *"full `GRCh38.chr6` (`170.8 Mbp`) × `CHM13.chr6` (`110.7 Mbp`) finished in `451.35 s` (`7m 31s`, `10,515.0 MB` RSS)—a `2.35×` total speedup over our `1,062.45 s` baseline (`12.3×` speedup in `LASTZ12` transition seeding from `588.67 s` to `47.87 s` via 16-bit prefix-bucketed `FlatSeedIndex` and Rayon parallelization; Figure 5)."*
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. **Stale commit (`53130b8` vs `HEAD` `33bede8`):** `manuscript/horizon_paper.md:105` and `docs/BENCHMARK_RESULTS.md:27,65` attribute the `451.35 s` run to commit `53130b8`, which is 5 commits behind `HEAD` (`33bede8`). Commit `33bede8` modified `src/mdl/dp.rs`, `src/mdl/mod.rs`, `src/seeding/flat_index.rs`, `src/seeding/translated.rs`, and `src/seeding/mod.rs`. Re-running the `chr6` 20 Mbp q-arm slice at `HEAD` (`33bede8`) changed both runtime (`31.75 s` vs `38.23 s`) and alignment counts (`15,015,161` pairs / `16,417` bubbles vs `15,025,501` pairs / `16,877` bubbles in `docs/BENCHMARK_RESULTS.md:26,63`).
  2. **No committed runner script or log:** Figure 5's values are hardcoded in `manuscript/scripts/generate_figures.py:669-685`; no script in `scripts/` runs the `1e47fe0` baseline or `53130b8` `chr6` alignment.
  3. **`CHM13.chr6` scope:** `/tmp/hprc_chm13_chr6.fa` has FASTA header `>CHM13.chr6:61401689-172126628` (`110,724,939` bp), which is the q-arm walk extracted from `hprc-v1.0-mc-grch38.gfa.gz`, not the full ~171 Mbp T2T-CHM13 chromosome 6.
- **Evidence:**
  - `git diff 53130b8..HEAD -- src/seeding/ src/mdl/dp.rs src/mdl/mod.rs`: modifies `BandGrid::try_new`, `mdl_xdrop_extend_with_params`, `flat_index.rs`, `translated.rs`.
  - `HORIZON_PROFILE=1 target/release/horizon_wga /tmp/grch38_chr6_q20m.fa /tmp/chm13_chr6_q20m.fa --resolve-dups --maf /dev/null --gfa /dev/null`:
    `[profile] w12=0.24s (fwd=65,rc=16,cov=75.8%) w10=0.01s lz1=0.00s seed_end=0.00s (lcbs=30+14) dp=25.12s (jobs=126160) xdrop=0.00s (jobs=98)`
    `aligned_pairs=15015161 unique_x_bp=14989622 unique_y_bp=14987477 total_bubbles=16417 elapsed_s=31.7475`
  - `head -n 1 /tmp/hprc_chm13_chr6.fa` -> `>CHM13.chr6:61401689-172126628`; `grep -v "^>" /tmp/hprc_chm13_chr6.fa | tr -d '\n' | wc -c` -> `110724939`.
- **Correction:** Disclose in `manuscript/horizon_paper.md:101,105` that the `451.35 s` (`10,515.0 MB` RSS) measurement was recorded at commit `53130b8` prior to `33bede8` and aligns `GRCh38.chr6` (`170.8 Mbp`) against the HPRC-extracted `CHM13.chr6:61401689-172126628` (`110.7 Mbp`) walk, or re-run and log `HEAD` (`33bede8`).

---

## F3 -- C3: `horizon_gfa stats` on `12 GB` HPRC v1.0 GFA reproduces on `HEAD`

- **Claim:** `C3` (`manuscript/horizon_paper.md:101`) — *"On the `12 GB` HPRC `v1.0` Minigraph-Cactus pangenome graph (`hprc-v1.0-mc-grch38.gfa.gz`), `horizon_gfa stats` streamed all `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks in `99.97 s` using `97.0 MB` peak RSS, while two-pass `horizon_gfa extract` isolated `GRCh38.chr22` + `CHM13.chr22` in `183.29 s` (`456.4 MB` RSS) and `GRCh38.chr6` + `CHM13.chr6` in `244.93 s` (`1,255.1 MB` RSS)."*
- **Verdict:** **SURVIVED**
- **Concrete flaw:** None on the `horizon_gfa stats` graph counts or memory bound (which were previously upheld on `488b020` in Round 3 `adjudication-L7.md:38` and re-verified here on `HEAD` `33bede8`). Note that the two `horizon_gfa extract` timings were not re-run in this round (their extracted FASTAs `/tmp/hprc_grch38_chr22.fa`, `/tmp/hprc_chm13_chr22.fa`, `/tmp/hprc_grch38_chr6.fa`, `/tmp/hprc_chm13_chr6.fa` are present on disk).
- **Evidence:**
  - `target/release/horizon_gfa stats /home/rcb112/project/gemini/bench/pg/data/benchmarks/hprc-v1.0-mc-grch38.gfa.gz`:
    `elapsed=104.98s peak_rss_mb=97.1 segments=81415956 segment_bp=3287932785 links=112955105 paths=2580 walks=24456`
  - Every graph count (`81,415,956` segments, `3.29 Gbp`, `112,955,105` links, `2,580` paths, `24,456` walks) matches identically and peak RSS (`97.1 MB` vs `97.0 MB`) reproduces within `0.1 MB`.
- **Correction:** None required for `horizon_gfa stats`; optionally update `99.97 s` / `97.0 MB` to `104.98 s` / `97.1 MB` on `33bede8` if pinning all numbers to `HEAD`.

---

## F4 -- C4: `horizon_clean` severs `walk_support <= 1` (not only `== 0`) and is a `0`-severed no-op on real `*-seq` GFAs

- **Claim:** `C4` (`manuscript/horizon_paper.md:49`) — *"3. `horizon_clean` (`src/bin/horizon_clean.rs`) indexes `S`-segment flanks and scores candidate `L`-line junctions with `bocpd_gc_profile` and an MDL two-regime changepoint test (`--min-cp-bits`), severing unsupported (`walk_support == 0`) discordant edges into `HB:i:1` Horizon Bubbles."*
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. `src/bin/horizon_clean.rs:260` severs edges when `(cfg.prune_unsupported && sup == 0) || (sup <= 1 && cp >= cfg.min_cp_bits)`. Thus, discordant edges with `walk_support == 1` are also severed when `cp >= min_cp_bits`, not solely `walk_support == 0`.
  2. When `S` records carry placeholder `*` sequences (as emitted by `horizon_wga --gfa`), `score_junction_cp_bits` returns `0.0` (`src/bin/horizon_clean.rs:101,154-176`), so `horizon_clean` is a `0`-severed no-op on real `horizon_wga` GFAs (`docs/BENCHMARK_RESULTS.md:66`: `2,776 kept, 0 severed (no-op on *-seq)`). Only two synthetic 60-bp unit tests (`src/bin/horizon_clean.rs:559-626`) exercise severing.
- **Evidence:**
  - `src/bin/horizon_clean.rs:260`: `let sever = (cfg.prune_unsupported && sup == 0) || (sup <= 1 && cp >= cfg.min_cp_bits);`
  - `docs/BENCHMARK_RESULTS.md:66`: `GFA Junction Cleaner (stdin spool) | horizon_clean --gfa - | *-seq pairwise chr6 peri GFA (2,776 sup>0 links) | 2,776 kept, 0 severed (no-op on *-seq)`
- **Correction:** Update `manuscript/horizon_paper.md:49` to *"severing low-support (`walk_support <= 1`) discordant edges (and `walk_support == 0` edges under `--prune-unsupported`) into `HB:i:1` Horizon Bubbles on sequence-bearing GFAs"*.

---

## F5 -- C5: *Plasmodium* (`87.7%` / `54.4%` / `66.8%`) and *Cannabis* (`99.95%`) benchmarks ran Python Prototype 1, not `horizon_wga`, and their input files are absent on disk

- **Claim:** `C5` (`manuscript/horizon_paper.md:65`) — *"End-to-end `horizon_wga` abstained on `87.7%` of those MDL-rejected columns while recovering `54.4%` of MDL-supported pairs (`66.8%` with `w = 10` sparse-exon fallback and `band_w = 180`). Applied with zero parameter changes to intraspecific *Cannabis sativa* (`cs10` × `Cherry_Pie`, `45.5%` GC, `418,663` columns), Horizon retained `99.95%` (`418,465`) of reference `MAF` columns and `1.000` recall on scratch realignments."*
- **Verdict:** **REFUTED**
- **Concrete flaw:**
  1. **Wrong tool attributed:** Neither benchmark was produced by Rust `horizon_wga`. Both come from **Python Prototype 1** in the sibling `pg` directory: `../pg/prototypes/pf_pv_end_to_end.py:1` (*"End-to-end Prototype 1 (KP-14 + banded MDL/BOCPD) vs Galaxy LASTZ rbest_maf on Pf vs Pv Chr1"*) and `../pg/prototypes/cannabis_mdl_audit.py:1-129` (*"Out-of-domain generalization check of Prototype 1 on Cannabis roast MAF blocks"*).
  2. **Unreproducible (input data absent on disk):** Both prototype scripts hardcode paths under a deleted scratch directory that does not exist on disk, and neither the scripts nor the datasets exist in `horizon-aligner`.
- **Evidence:**
  - `ls -la /home/agent/.gemini/antigravity-cli/brain/c77b5ff5-88ad-43af-973a-c37dbadc186d/scratch/cannabis_test/cs10_small_sample.maf`: `No such file or directory` (referenced at `../pg/prototypes/cannabis_mdl_audit.py:19-21`).
  - `ls -la /home/agent/.gemini/antigravity-cli/brain/c77b5ff5-88ad-43af-973a-c37dbadc186d/scratch/plasmodium_fastas`: `No such file or directory` (referenced at `../pg/prototypes/pf_pv_crossclade_eval.py:22-24`).
- **Correction:** Remove or re-run the *Plasmodium* (`87.7%` / `54.4%` / `66.8%`) and *Cannabis* (`99.95%`) benchmarks using the Rust `horizon_wga` binary on committed/archived inputs, or explicitly label them as historical Python Prototype 1 results whose scratch inputs are not archived.

---

## F6 -- C6: The `30` randomized e-graph trials do not exist in `horizon-aligner` and test the superseded Python prototype's unitless cost model, not Rust `MdlCostModel`

- **Claim:** `C6` (`manuscript/horizon_paper.md:99`) — *"Across `30` randomized tandem-duplication and overlapping-inversion trials, e-graph saturation recovered `1.00` ancestral block sets and `1.00` exact `dup(X)` duplications compared to `0.00` under greedy one-to-one synteny filtering."*
- **Verdict:** **REFUTED**
- **Concrete flaw:**
  1. **Absent from `horizon-aligner`:** No script, benchmark, or Rust unit test in `horizon-aligner` executes `30` randomized e-graph trials.
  2. **Tests superseded Python cost model, not Rust `MdlCostModel`:** The `30`-trial benchmark resides in `../pg/prototypes/egraph_synteny_resolver.py:210-258` (`run_benchmark(n_trials=30, seed=42)`), which uses unitless parsimony weights `op_cost = {"seq": 1.0, "inv": 1.2, "dup": 1.3}` (`line 161`), not the Rust `MdlCostModel` in `src/egraph.rs` (which prices `dup(X)` via catalog pointer bits and sequence divergence residuals $H_2(d) + d \log_2 3$). Moreover, the Python greedy baseline (`line 235`: `if lbl not in seen_labels and st > 0`) discards all inverted blocks by construction when every trial inverts a 2–3-block span (`lines 225-228`).
- **Evidence:**
  - `python3 /home/rcb112/project/gemini/bench/pg/prototypes/egraph_synteny_resolver.py`:
    `{'trials': 30.0, 'greedy_uninverted_block_recall': 0.0, 'greedy_exact_dup_label_recall': 0.0, 'egraph_ancestral_block_recall': 1.0, 'egraph_exact_dup_label_recall': 1.0}`
  - `grep -rn "30.*trial\|randomized" src/egraph.rs scripts/`: `0` matches.
- **Correction:** Implement and run a `30`-trial randomized benchmark directly against Rust `src/egraph.rs` (`MdlCostModel` + `saturate_synteny` / `extract_min_mdl_cost`) with an orientation-aware greedy baseline, or drop the `30`-trial claim from `manuscript/horizon_paper.md:99`.

---

## Verdict Summary

| Claim ID | Verdict | Summary |
| :--- | :--- | :--- |
| `C1` | **NARROWED** | `10` binaries in `src/bin/` (`horizon` dispatcher + `9` subcommands); only GFA tools stream `stdin`; `horizon_gfa`/`horizon_lift` do not use MDL/BOCPD. |
| `C2` | **NARROWED** | `451.35 s` chr6 run was pinned at `53130b8` (pre-`33bede8`), lacks a committed runner script, and `CHM13.chr6` is `61401689-172126628` (`110.7 Mbp` q-arm). |
| `C3` | **SURVIVED** | `horizon_gfa stats` on `12 GB` HPRC `v1.0` reproduced on `HEAD` (`104.98 s`, `97.1 MB` RSS, exact `81,415,956` segs / `112,955,105` links / `2,580` paths / `24,456` walks). |
| `C4` | **NARROWED** | `horizon_clean.rs:260` severs `sup <= 1` (not only `walk_support == 0`) and is a `0`-severed no-op on `*-seq` GFAs. |
| `C5` | **REFUTED** | Attributed to `horizon_wga`, but actually Python Prototype 1 (`pf_pv_end_to_end.py`, `cannabis_mdl_audit.py`) with deleted/missing scratch input files. |
| `C6` | **REFUTED** | Absent from `horizon-aligner`; `../pg/prototypes/egraph_synteny_resolver.py` tests the superseded unitless Python cost model, not Rust `MdlCostModel`. |
