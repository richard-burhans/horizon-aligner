# Round 6 Adjudication -- L7 (Reproduction)

Target commit: `42700e1` (`HEAD`).
Attack digest verified: `attack-L7-reproduction.md` SHA-256 prefix `2b54f1d71cb9` (`2b54f1d71cb9e5100ae7ba9906f11091c2543938f069dd7a92f27a0038c86295`).
Prior round history checked via `scripts/round history`: `horizon_filter` and `horizon_locus` have no prior adjudications; `horizon_resolve` has 112 prior mentions from earlier rounds predating its current `PackedGraph` + Tarjan snarl implementation.

---

## F1 -- C1 unit-test count is stale (`118` vs. `120` passing) and `cargo fmt --check` fails across 22 files

- **Ruling:** **UPHELD**
- **Verification:**
  - `PATH="$HOME/.cargo/bin:/tmp/ext_tools/env/bin:$PATH" cargo test --all-features --quiet` on `42700e1` passes **120** unit tests (`86` lib passed + `1` lib ignored `real_hal_decode_from_env`; `34` bin tests passed across `13` binary test targets: `1 + 2 + 6 + 2 + 1 + 2 + 4 + 11 + 1 + 1 + 1 + 2 + 0`), not `118` (`86` lib + `32` bin).
  - `cargo clippy --all-targets --all-features -- -D warnings` passes with exit code `0`.
  - `cargo fmt --check` exits with code `1` with `54` diff hunks across **22** unique Rust files. Note that `C1` (`docs/PANGENOME_SUITE.md:3`) is explicitly tagged `(2026-11-20, 9dfec4a)` and does not itself claim `cargo fmt --check` passes (that claim is in `docs/STRATEGIC_PLAN.md:3`), while the stale `118` count also appears in `docs/STRATEGIC_PLAN.md:3`, `docs/VGP_INTEGRATION_PLAN.md:4`, and `manuscript/horizon_paper.md:84`.
- **What survives:** The `118` test count in `docs/PANGENOME_SUITE.md:3` (and sibling status notes) is stale by `2` binary unit tests on `42700e1` (`120` passing: `86` lib + `34` bin, `1` ignored), and `cargo fmt --check` fails on `22` files.
- **Repair:** **ENDORSED** for updating the documented test count and commit tag to `42700e1` (`120 passing unit tests: 86 lib + 34 bin, 1 ignored`) — **the fix is free** (text-only update; `cargo fmt` is also free if formatting is normalized).

## F2 -- C2 `horizon_filter --net` keeps `565` records on HEAD, not `1,004`, and has no tracked runner or `numpy` environment

- **Ruling:** **NARROWED**
- **Verification:**
  - Re-running `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --net --out /tmp/adj_c2_net.paf` on `42700e1` outputs `[horizon_filter] total=56857 kept=565 filtered_mdl=1 filtered_len=0 elapsed=0.47s` (`565` kept records, not `1,004`). With `--max-occ 2`, it outputs `kept=7170`, confirming `1,004` is not reproduced on `42700e1`.
  - With `/tmp/ext_tools/env/bin` in `PATH`, `.paf.zst` input succeeds in `0.47s` with `kept=565`; without it, `zstd` is missing from default `PATH`.
  - No script in `scripts/benchmarks/` or fixture in `data/` invokes `horizon_filter` or tracks `/tmp/chr6_fastga.paf`.
  - `/tmp/ext_tools/env/bin/python3 -c "import numpy"` fails with `ModuleNotFoundError: No module named 'numpy'`, blocking execution of `scripts/benchmarks/score_paf_vs_hprc.py`.
- **What survives:** The `1,004` kept-record count in `C2` (`docs/BENCHMARK_RESULTS.md:85`) does not reproduce on `42700e1` (`565` records), the benchmark has no tracked driver/fixture, `.paf.zst` relies on `zstd` in `PATH`, and `score_paf_vs_hprc.py` cannot run without `numpy`. However, the `56,857`-record `0.47s` (`0.47s–0.50s` `.paf.zst`) runtime *does* reproduce, and the `81.2%` / `99.7%` (`110.22 Mbp`) accuracy metrics are **unverified** in this environment rather than proven false.
- **Repair:** **OVERTURNED** on simply replacing `1,004` with `565` in `docs/BENCHMARK_RESULTS.md:85` (`README.md:124`) while leaving `110.22 Mbp`, `81.2%`, and `99.7%` untouched, as that would pair a `42700e1` record count with un-re-scored metrics; re-scoring the `565`-record output requires `numpy`, or else the row must be pinned to the commit/flags that produced `1,004`.

## F3 -- C3 `horizon_locus` human 4-way 25 kb produces `102` MAF blocks (`92` multi-species), not `29` / `20`, and inputs/scripts are untracked

- **Ruling:** **UPHELD**
- **Verification:**
  - `./target/release/horizon_locus --fasta /tmp/locus_avian.fa --out-prefix /tmp/adj_c3_avian --threads 4` reproduces: `seqs=2 pairs=1 focal=galGal_Gallus_gallus.chr1 maf_blocks=1364 multi_blocks=0 elapsed=0.82s peak_rss_mb=240.9`, `aligned_bp=231444`, `p_dist=0.192552`, `jc69_dist=0.222528`.
  - `./target/release/horizon_locus --fasta /tmp/locus_human_25k.fa --out-prefix /tmp/adj_c3_human --threads 4` outputs `seqs=4 pairs=6 focal=chm13_peri maf_blocks=102 multi_blocks=92 elapsed=0.41s peak_rss_mb=172.2` (and `grch38_mhc` vs `hla_dup` `p_dist=0.041539`, vs `0.124507` in `docs/BENCHMARK_RESULTS.md:91`).
  - Neither `/tmp/locus_avian.fa`, `/tmp/locus_human_25k.fa`, nor any `horizon_locus` benchmark script is tracked in `scripts/benchmarks/` or `data/`.
- **What survives:** Avian 2-way `horizon_locus` metrics reproduce on `42700e1`, while human 4-way 25 kb emits `102` MAF blocks (`92` multi-species) in `0.41s–0.42s` (`~172 MB` RSS) instead of `29` / `20` in `0.32s` (`165.8 MB` RSS), and inputs/scripts are untracked.
- **Repair:** **ENDORSED** for updating `README.md:126` and `docs/BENCHMARK_RESULTS.md:91` to `102` MAF blocks (`92` multi-species) and the `42700e1` runtime/RSS (`0.41s–0.42s`, `~172 MB` RSS) — **the fix is free** for the doc text (tracking fixtures/scripts is separate).

## F4 -- C4 whole-`chr6` numbers are pinned to historical commit `8cd5db0` (`14` commits behind HEAD) and rely on untracked `/tmp/` artifacts and missing `numpy`

- **Ruling:** **NARROWED**
- **Verification:**
  - `manuscript/horizon_paper.md:78` explicitly labels the `258.88–267.86 s` / `282,193,564` pairs / `178,737` bubbles run as `on 8cd5db0`, and `/tmp/chr6_8cd5db0.log` confirms those exact figures (`aligned_pairs=282193564 unique_x_bp=124160423 unique_y_bp=109464447 total_bubbles=178737 elapsed_s=258.8827`).
  - `scripts/benchmarks/compare_chr6_horizon_vs_hprc.py` fails with `ModuleNotFoundError: No module named 'numpy'`, and `/tmp/hprc_grch38_chr6.fa` / `/tmp/hprc_chr6_truth_segments.tsv` are untracked in `/tmp/`.
  - **Attacker premise error:** `git rev-list --count 8cd5db0..HEAD` (and `--first-parent --count`) is **28** commits, not `14`.
- **What survives:** `C4` (`manuscript/horizon_paper.md:78`) accurately cites `8cd5db0` rather than claiming `HEAD`, whereas `docs/BENCHMARK_RESULTS.md:5` overclaims `"All benchmarks below were re-verified on HEAD"` for `8cd5db0` rows; `/tmp/` inputs are untracked and `compare_chr6_horizon_vs_hprc.py` is blocked by missing `numpy`. The attacker's `14`-commit gap is overturned (`28` commits).
- **Repair:** **ENDORSED** for scoping `docs/BENCHMARK_RESULTS.md:5` so it does not claim `8cd5db0` rows were re-verified on `HEAD` — **the fix is free**. Re-running full `chr6` on `42700e1` and installing `numpy` are non-free compute/environment tasks.

## F5 -- C5 `horizon_resolve` functional claims reproduce on HEAD

- **Ruling:** **UPHELD**
- **Verification:**
  - `cargo test --all-features` passes `resolves_tandem_dup_and_cancels_double_inversion_in_gfa_walk` in `src/bin/horizon_resolve.rs`.
  - Running `./target/release/horizon_resolve /tmp/adj_f5.gfa --out /tmp/adj_f5_out.gfa --collapse-dups --tag-walks --min-anchor-bp 10 --divergence 0.05` on a synthetic `S`/`L`/`P`/`W` GFA reproduces `P` rewriting (`a+,b+,b+,c+` -> `a+,b+,c+`), `W` rewriting (`>a>b>b>c` -> `>a>b>c EG:Z:seq(a,seq(dup(b),c))`), `W`-line `end` bound updating (`90` -> `70`), and unused `L`-link pruning (`b+ -> b+` and `a+ -> d+` removed).
- **What survives:** `C5` (`ARCHITECTURE.md:28`) survives all reproduction checks on `42700e1`.

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `cargo test --all-features` passes `120` tests (`86` lib + `34` bin, `1` ignored) rather than `118`, and `cargo fmt --check` fails on `22` files; doc update is free |
| F2 | NARROWED | `horizon_filter --net` keeps `565` records (not `1,004`) and `score_paf_vs_hprc.py` is blocked by missing `numpy`, but `0.47s` runtime reproduces and swapping `565` in without re-scoring `81.2%`/`99.7%` is rejected |
| F3 | UPHELD | Avian `horizon_locus` metrics reproduce, but human 4-way 25 kb emits `102` MAF blocks (`92` multi-species) in `0.41s–0.42s` rather than `29`/`20` in `0.32s`, and inputs/scripts are untracked; doc fix is free |
| F4 | NARROWED | `8cd5db0` log matches `C4` and `docs/BENCHMARK_RESULTS.md:5` overclaims `re-verified on HEAD` (free fix), `compare_chr6_horizon_vs_hprc.py` lacks `numpy`, and `42700e1` is `28` commits ahead of `8cd5db0` rather than `14` |
| F5 | UPHELD | `horizon_resolve` unit test and synthetic `P`/`W`/`L` GFA run reproduce Tarjan snarl decomposition, `--divergence`, `P`/`W` rewriting, `W`-line `end` updating, and `L`-link pruning |
