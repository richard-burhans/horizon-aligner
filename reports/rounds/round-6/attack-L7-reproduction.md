# Round 6 Attack -- L7 (Reproduction)

Checked on commit `42700e1` (`HEAD`). Per `docs/round_rubric.md`, this project keeps no method cards yet; prior round adjudications were checked via `scripts/round history` (`horizon_filter` and `horizon_locus` have no prior adjudications; `horizon_resolve` was adjudicated in earlier rounds prior to its current implementation). All findings below are grounded in direct execution of `cargo`, release binaries in `target/release/`, and inspection of `/tmp/` and `scripts/benchmarks/`.

## F1 -- C1 unit-test count is stale (`118` vs. `120` passing) and `cargo fmt --check` fails across 22 files
- **Claim:** `C1` (`docs/PANGENOME_SUITE.md:3`)
- **Verdict:** **NARROWED**
- **Evidence:**
  1. Running `PATH="$HOME/.cargo/bin:/tmp/ext_tools/env/bin:$PATH" cargo test --all-features` on `42700e1` passes **120** unit tests (`86` lib passed, `1` lib ignored `real_hal_decode_from_env`; `34` bin tests passed across `13` binary test runs: `1 + 2 + 6 + 2 + 1 + 2 + 4 + 11 + 1 + 1 + 1 + 2 + 0`), not `118` (`86` lib + `32` bin) as stated in `docs/PANGENOME_SUITE.md:3` (which is still tagged `2026-11-20, 9dfec4a`).
  2. `PATH="$HOME/.cargo/bin:/tmp/ext_tools/env/bin:$PATH" cargo clippy --all-targets --all-features -- -D warnings` passes with exit code `0`.
  3. `PATH="$HOME/.cargo/bin:/tmp/ext_tools/env/bin:$PATH" cargo fmt --check` exits with code `1` with formatting diffs across **22** Rust files (`src/bin/horizon_clean.rs`, `src/bin/horizon_eval.rs`, `src/bin/horizon_filter.rs`, `src/bin/horizon_gfa.rs`, `src/bin/horizon_guide/main.rs`, `src/bin/horizon_lift/features.rs`, `src/bin/horizon_lift/main.rs`, `src/bin/horizon_locus.rs`, `src/bin/horizon_resolve.rs`, `src/egraph.rs`, `src/gfa/mod.rs`, `src/hal/hdf5_write.rs`, `src/hal/tests.rs`, `src/io_util.rs`, `src/lib.rs`, `src/mdl/dp.rs`, `src/mdl/fallback.rs`, `src/mdl/mod.rs`, `src/mdl/runs.rs`, `src/merge.rs`, `src/netting.rs`, `src/paf.rs`).
- **Correction:** Update `docs/PANGENOME_SUITE.md:3` to `42700e1` and `120 passing unit tests: 86 lib + 34 bin, 1 ignored`, and run `cargo fmt` (or note that `cargo fmt --check` currently fails).

## F2 -- C2 `horizon_filter --net` keeps `565` records on HEAD, not `1,004`, and has no tracked runner or `numpy` environment
- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:85`)
- **Verdict:** **REFUTED**
- **Evidence:**
  1. Re-running `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --net --out /tmp/repro_c2_net.paf` on `42700e1` outputs:
     `[horizon_filter] total=56857 kept=565 filtered_mdl=1 filtered_len=0 elapsed=0.48s`
     (`565` kept records, not `1,004` records as claimed in `docs/BENCHMARK_RESULTS.md:85`).
  2. Running `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf.zst --net` without `/tmp/ext_tools/env/bin` in `PATH` fails with `spawning zstd -dc for '/tmp/chr6_fastga.paf.zst': No such file or directory (os error 2)`. With `/tmp/ext_tools/env/bin` in `PATH`, it completes in `0.50s` with `kept=565`.
  3. No benchmark script in `scripts/benchmarks/` or fixture in `data/` invokes `horizon_filter` or tracks `/tmp/chr6_fastga.paf`. Moreover, running `scripts/benchmarks/score_paf_vs_hprc.py` under both system `python3` and `/tmp/ext_tools/env/bin/python3` fails with `ModuleNotFoundError: No module named 'numpy'`, preventing reproduction of the claimed `81.2%` 2D HPRC truth hit rate and `99.7%` 1-to-1 coverage.
- **Correction:** Track the `horizon_filter` benchmark driver and document the `zstd`/`numpy` dependencies, and update `docs/BENCHMARK_RESULTS.md:85` (and `README.md:124`) to report the reproduced `565` kept records (or the exact CLI flags/commit that produced `1,004`).

## F3 -- C3 `horizon_locus` human 4-way 25 kb produces `102` MAF blocks (`92` multi-species), not `29` / `20`, and inputs/scripts are untracked
- **Claim:** `C3` (`README.md:126`)
- **Verdict:** **NARROWED**
- **Evidence:**
  1. Re-running `./target/release/horizon_locus --fasta /tmp/locus_avian.fa --out-prefix /tmp/repro_c3_avian --threads 4` on `42700e1` reproduces the avian 2-way numbers:
     `[horizon_locus] seqs=2 pairs=1 focal=galGal_Gallus_gallus.chr1 maf_blocks=1364 multi_blocks=0 elapsed=0.83s peak_rss_mb=239.3`
     with `aligned_bp=231444`, `p_dist=0.192552`, `jc69_dist=0.222528`.
  2. Re-running `./target/release/horizon_locus --fasta /tmp/locus_human_25k.fa --out-prefix /tmp/repro_c3_human --threads 4` on `42700e1` outputs:
     `[horizon_locus] seqs=4 pairs=6 focal=chm13_peri maf_blocks=102 multi_blocks=92 elapsed=0.42s peak_rss_mb=171.8`
     (`102` MAF blocks / `92` multi-species, contradicting the claimed `29` MAF blocks / `20` multi-species and `0.32s` / `165.8 MB` RSS in `README.md:126` and `docs/BENCHMARK_RESULTS.md:91`).
  3. Neither `/tmp/locus_avian.fa`, `/tmp/locus_human_25k.fa`, nor any `horizon_locus` benchmark script is tracked in `scripts/benchmarks/` or `data/`.
- **Correction:** Update `README.md:126` and `docs/BENCHMARK_RESULTS.md:91` to `102` MAF blocks (`92` multi-species) and `0.42s` (`171.8 MB` RSS) on HEAD, and track a reproducible `horizon_locus` benchmark script/fixture.

## F4 -- C4 whole-`chr6` numbers are pinned to historical commit `8cd5db0` (`14` commits behind HEAD) and rely on untracked `/tmp/` artifacts and missing `numpy`
- **Claim:** `C4` (`manuscript/horizon_paper.md:78`)
- **Verdict:** **NARROWED**
- **Evidence:**
  1. `manuscript/horizon_paper.md:78` cites `8cd5db0` (`258.88–267.86 s`), `bb43dd9` (`592.68 s`), and `1e47fe0` (`1,062.45 s`). `HEAD` (`42700e1`) is `14` commits ahead of `8cd5db0` (`git diff --stat 8cd5db0..HEAD` touches `139` files including `src/mdl/mod.rs`, `src/mdl/dp.rs`, `src/seeding/`, `src/egraph.rs`, `src/formats.rs`, `src/merge.rs`, `src/netting.rs`). Meanwhile `docs/BENCHMARK_RESULTS.md:5` asserts `"All benchmarks below were re-verified on HEAD"` even though line 73 still cites `8cd5db0`.
  2. The existing `/tmp/chr6_8cd5db0.log` matches the reported `8cd5db0` numbers (`aligned_pairs=282193564 unique_x_bp=124160423 unique_y_bp=109464447 total_bubbles=178737 elapsed_s=258.8827`), but the input FASTAs (`/tmp/hprc_grch38_chr6.fa`, etc.) and truth TSV (`/tmp/hprc_chr6_truth_segments.tsv`) are untracked in `/tmp/`, and `scripts/benchmarks/compare_chr6_horizon_vs_hprc.py` fails to run in this environment (`ModuleNotFoundError: No module named 'numpy'`).
- **Correction:** Re-run full `chr6` `horizon_wga --resolve-dups` on `42700e1` (or remove the `"re-verified on HEAD"` claim in `docs/BENCHMARK_RESULTS.md:5` for `8cd5db0` rows) and document/provide the `numpy` environment required by `scripts/benchmarks/compare_chr6_horizon_vs_hprc.py`.

## F5 -- C5 `horizon_resolve` functional claims reproduce on HEAD
- **Claim:** `C5` (`ARCHITECTURE.md:28`)
- **Verdict:** **SURVIVED**
- **Evidence:**
  1. `cargo test --all-features` passes `resolves_tandem_dup_and_cancels_double_inversion_in_gfa_walk` in `src/bin/horizon_resolve.rs:899-943`.
  2. Running `./target/release/horizon_resolve` on a synthetic GFA containing `S`, `L`, `P`, and `W` records with `--collapse-dups --tag-walks --min-anchor-bp 10 --divergence 0.05` verified:
     - `PackedGraph` + Tarjan bridge-finding (`compute_biedged_snarls` / `tarjan_bridges` in `src/bin/horizon_resolve.rs:210-375`),
     - `P` and `W` walk rewriting (`P p1 a+,b+,b+,c+ 0M,0M,0M` -> `P p1 a+,b+,c+ *`; `W HG002 1 chr6 10 90 >a>b>b>c` -> `W HG002 1 chr6 10 70 >a>b>c EG:Z:seq(a,seq(dup(b),c))`),
     - `W`-line `end` bound updating (`90` -> `70`),
     - Unused `L`-link pruning (`L b + b + 0M` and `L a + d + 0M` removed) and missing `L`-link synthesis (`L a + b - 0M` and `L b - c + 0M` emitted).

## Verdict Table

| Claim | Verdict | Summary |
| :--- | :--- | :--- |
| **C1** | **NARROWED** | `cargo test --all-features` passes `120` tests (`86` lib + `34` bin, `1` ignored), not `118`; `cargo clippy` passes; `cargo fmt --check` fails on `22` files. |
| **C2** | **REFUTED** | `horizon_filter --net` keeps `565` records on `42700e1`, not `1,004`; `.paf.zst` requires `zstd` in `PATH`; no tracked benchmark script; `score_paf_vs_hprc.py` fails (`No module named 'numpy'`). |
| **C3** | **NARROWED** | Avian `horizon_locus` metrics reproduce (`0.83s`, `239.3 MB` RSS, `231,444` bp, `p_dist=0.192552`), but human 4-way 25 kb emits `102` MAF blocks (`92` multi-sp) in `0.42s`, not `29` / `20` in `0.32s`; inputs/scripts untracked. |
| **C4** | **NARROWED** | `8cd5db0` log confirms historical numbers, but `42700e1` is `14` commits ahead (`139` files changed), contradicting `docs/BENCHMARK_RESULTS.md:5` (`re-verified on HEAD`); `/tmp/` inputs untracked and `compare_chr6_horizon_vs_hprc.py` fails (`No module named 'numpy'`). |
| **C5** | **SURVIVED** | `horizon_resolve` unit test and synthetic `P`/`W`/`L` GFA CLI run reproduce Tarjan snarl decomposition, `--divergence`, `P`/`W` rewriting, `W`-line `end` updating, and `L`-link pruning/synthesis. |
