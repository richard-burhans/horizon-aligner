# Round 6 Attack — L4 (Provenance)

**Target revision:** `42700e142b70996c23cdf391ac29e1e5bcf4640d` (`42700e1`)
**Lens:** Provenance (frozen snapshots, commit pins, script paths, and benchmark artifacts)
**Method cards & prior history:** Per `docs/round_rubric.md`, this repository keeps no method cards. Prior rounds were checked via `scripts/round history "provenance" --kind adjudication` and `scripts/round history "scripts/benchmarks"`, plus `reports/adversarial_log.md` (Round 5's `R2-INDEX-RAM` noted `/tmp/` script provenance for `8cd5db0`/`2896229`, which was partially addressed by committing `compare_chr6_horizon_vs_hprc.py`, `extract_hprc_chr6_truth.py`, `score_paf_vs_hprc.py`, and `plot_chr6_hprc_fig.py` into `scripts/benchmarks/`). No remote or `../library` literature claims are made.

---

## F1 -- C1: Unit test count is 120 active (86 lib + 34 bin), not 118 (86 lib + 32 bin)
- **Claim:** `C1` (`reports/rounds/round-6/ATTACKED-PANGENOME_SUITE.md:3`, also in `docs/VGP_INTEGRATION_PLAN.md:4`, `docs/STRATEGIC_PLAN.md:3`, and `reports/rounds/round-6/ATTACKED-horizon_paper.md:84`)
- **Verdict:** **NARROWED**
- **Concrete flaw:** The claim states `118` passing unit tests (`86` lib + `32` bin, `1` ignored) on `9dfec4a` (and `42700e1`, whose commits `1ddc872` and `42700e1` touched only docs and `manuscript/horizon_paper.pdf`). Running `cargo test --all-features` at `42700e1` executes **120** active unit tests (`86` lib + `34` bin, plus `1` ignored in lib = `121` total).
- **Evidence:**
  - `cargo test --all-features` output at `42700e1`: `src/lib.rs` (`86 passed; 1 ignored`), `bench_chr1` (`1`), `horizon` (`2`), `horizon_clean` (`6`), `horizon_eval` (`2`), `horizon_filter` (`1`), `horizon_gfa` (`2`), `horizon_guide` (`4`), `horizon_lift` (`11`), `horizon_locus` (`1`), `horizon_merge` (`1`), `horizon_resolve` (`1`), `horizon_wga` (`2`). Sum of bin tests = `1+2+6+2+1+2+4+11+1+1+1+2 = 34` (`86 + 34 = 120` active).
  - Binary count (`11` suite binaries + `bench_chr1` = 12 files in `src/bin/`) and `0` clippy warnings (`cargo clippy --all-targets --all-features -- -D warnings`) are verified.
- **Correction:** Replace `**118** passing unit tests: \`86\` lib + \`32\` bin, \`1\` ignored` with `**120** passing unit tests: \`86\` lib + \`34\` bin, \`1\` ignored` across `docs/PANGENOME_SUITE.md:3`, `docs/VGP_INTEGRATION_PLAN.md:4`, `docs/STRATEGIC_PLAN.md:3`, and `manuscript/horizon_paper.md:84`.

## F2 -- C2: All `horizon_filter` VGP/chr6 artifacts in `/tmp/` predate `9dfec4a` and have zero tracked scripts or datasets
- **Claim:** `C2` (`reports/rounds/round-6/ATTACKED-BENCHMARK_RESULTS.md:5,81-88`, `reports/rounds/round-6/ATTACKED-README.md:123-125`, `reports/rounds/round-6/ATTACKED-VGP_INTEGRATION_PLAN.md:4,68,72`, `reports/rounds/round-6/ATTACKED-horizon_paper.md:84`)
- **Verdict:** **REFUTED**
- **Concrete flaw:** C2 attributes its `horizon_filter` benchmarks to commit `9dfec4a` (`2026-11-20 19:38:58 -0500`), and `docs/BENCHMARK_RESULTS.md:5` claims *"All benchmarks below were re-verified on `HEAD`"*. Every C2 output artifact in `/tmp/` predates `9dfec4a` (and `42700e1`), and no script or dataset fetcher in the repository reproduces any C2 benchmark.
- **Evidence:**
  - `git log -n 10 --format="%h %ad %s" --date=iso`: `39b5990` (first `horizon_filter.rs` commit) is `2026-11-20 19:12:27 -0500`; `4dba79e` (`19:25:28`) and `38e2d58` (`19:28:24`) modified `src/bin/horizon_filter.rs`, `src/netting.rs`, `src/paf.rs`, and `src/io_util.rs`; `9dfec4a` (`19:38:58`) modified `src/bin/horizon_filter.rs` and fixed `ZstdReadChild` deadlock in `src/io_util.rs`.
  - `stat -c "%y %n"` on `/tmp/` artifacts:
    - `/tmp/chr6_fastga_mdl_net.paf`: `2026-11-20 19:10:09 -0500` (28 min before `9dfec4a`, 2 min before `39b5990`; contains `553` lines vs. the claimed `1,004` / `848` records).
    - `/tmp/chr6_filtered.paf.zst`: `2026-11-20 19:19:58 -0500` (19 min before `9dfec4a`).
    - `/tmp/vgp_test/horizon_net.paf.gz`: `2026-11-20 19:30:26 -0500` (8.5 min before `9dfec4a`).
  - `grep -rnE "(horizon_filter|koala|amphioxus|GCA_003287225|GCA_015852565|FastGA_11)" scripts/ docs/BENCHMARK_DATASETS.md` returns zero matches in `scripts/` and `docs/BENCHMARK_DATASETS.md`. Neither `scripts/benchmarks/fetch_benchmarks.sh` nor any script in `scripts/benchmarks/` fetches or runs `chr6` `FastGA` PAFs or GenomeArk Koala × Amphioxus.
- **Correction:** Re-run all `horizon_filter` benchmarks on `42700e1` (or `HEAD`) via a tracked script in `scripts/benchmarks/` with documented dataset provenance in `docs/BENCHMARK_DATASETS.md`, or mark the `9dfec4a` attribution and `HEAD` re-verification claim in `docs/BENCHMARK_RESULTS.md:5,81` as unverified/pre-`9dfec4a` `/tmp/` runs.

## F3 -- C3: All `horizon_locus` outputs in `/tmp/` predate `9dfec4a` (avian/human predating `4dba79e` RC-complement fix) and have zero tracked scripts
- **Claim:** `C3` (`reports/rounds/round-6/ATTACKED-README.md:123,125-126`, `reports/rounds/round-6/ATTACKED-BENCHMARK_RESULTS.md:5,81,89-91`, `reports/rounds/round-6/ATTACKED-VGP_INTEGRATION_PLAN.md:4,68,75`, `reports/rounds/round-6/ATTACKED-horizon_paper.md:84`)
- **Verdict:** **REFUTED**
- **Concrete flaw:** C3 attributes all three `horizon_locus` benchmarks (Koala × Amphioxus 25 kb, avian 2-way `galGal` × `taeGut`, human 4-way 25 kb) to `9dfec4a` (`2026-11-20 19:38:58 -0500`), and `docs/BENCHMARK_RESULTS.md:5` claims re-verification on `HEAD`. In fact, all three `/tmp/` output files predate `9dfec4a`, and the avian and human 4-way outputs even predate `4dba79e` (`2026-11-20 19:25:28 -0500`), which fixed a reverse-complement match-counting bug in `src/bin/horizon_locus.rs`.
- **Evidence:**
  - `git show 4dba79e -- src/bin/horizon_locus.rs` (`2026-11-20 19:25:28 -0500`) changed `count_run_matches` to apply `complement_base(seq_y[y])` when `r.is_rc` is true, directly affecting `matches`, `p_dist`, and `jc69_dist`.
  - `stat -c "%y %n"` on `/tmp/` outputs:
    - `/tmp/locus_avian_out.similarity.tsv`: `2026-11-20 19:15:46 -0500` (10 min before `4dba79e`, 23 min before `9dfec4a`).
    - `/tmp/locus_human_out.similarity.tsv`: `2026-11-20 19:16:09 -0500` (9 min before `4dba79e`, 23 min before `9dfec4a`).
    - `/tmp/vgp_test/locus_out.similarity.tsv`: `2026-11-20 19:31:17 -0500` (7.7 min before `9dfec4a`).
  - `grep -rnE "(horizon_locus|galGal|taeGut|grch38_mhc|locus_avian|locus_human)" scripts/ docs/BENCHMARK_DATASETS.md` returns zero matches. All input FASTAs (`locus_avian.fa`, `locus_human_25k.fa`, `vgp_test/locus_koala_amphioxus.fa`) and outputs exist only in untracked `/tmp/`.
- **Correction:** Commit a reproducible benchmark script in `scripts/benchmarks/` that builds or fetches the three locus FASTAs and re-runs `horizon_locus` on `42700e1`, or remove the `9dfec4a` / `HEAD` provenance claims from `docs/BENCHMARK_RESULTS.md:5,81`, `README.md:123`, and `manuscript/horizon_paper.md:84`.

## F4 -- C4: `docs/BENCHMARK_RESULTS.md:5` falsely claims `HEAD` re-verification across 26 post-`8cd5db0` commits, and tracked scripts still hardcode `/tmp/`
- **Claim:** `C4` (`reports/rounds/round-6/ATTACKED-horizon_paper.md:78`, `reports/rounds/round-6/ATTACKED-BENCHMARK_RESULTS.md:5,28-29,73-75`, `reports/rounds/round-6/ATTACKED-README.md:110-122`)
- **Verdict:** **NARROWED**
- **Concrete flaw:** Commits `8cd5db0` (`2026-11-20 14:37:50 -0500`) and `2896229` (`2026-11-20 15:18:51 -0500`) exist in git history, and `scripts/benchmarks/` now tracks `extract_hprc_chr6_truth.py`, `compare_chr6_horizon_vs_hprc.py`, `score_paf_vs_hprc.py`, and `plot_chr6_hprc_fig.py` (partially addressing Round 5's `R2-INDEX-RAM`). However, `docs/BENCHMARK_RESULTS.md:5` asserts *"All benchmarks below were re-verified on `HEAD`"*, which contradicts both the table's own `8cd5db0` pin and git history: **26 commits** landed between `8cd5db0` and `42700e1`, modifying `src/mdl/mod.rs`, `src/mdl/dp.rs`, `src/mdl/fallback.rs`, `src/bocpd.rs`, `src/egraph.rs`, `src/formats.rs`, `src/seeding/`, `src/bin/horizon_wga.rs`, and `src/bin/horizon_eval.rs`.
- **Evidence:**
  - `git log --oneline 8cd5db0..42700e1` lists 26 commits (including `2d0b754`, `62402b8`, `4c91f71`, `87c11cf`, `d8339b7`, `4dba79e`, `38e2d58`, `9dfec4a`).
  - `scripts/benchmarks/compare_chr6_horizon_vs_hprc.py:9-14`, `extract_hprc_chr6_truth.py:8`, `score_paf_vs_hprc.py:10-16`, `plot_chr6_hprc_fig.py:43-44`, and `bench_profile_ledger.py:103-106` still default to `/tmp/...` paths, and `scripts/benchmarks/profile_ledger.jsonl` records only a single `chr6_peri_5m` run at `62402b8` (none at `8cd5db0`, `2896229`, or `42700e1`).
- **Correction:** Remove the blanket *"All benchmarks below were re-verified on `HEAD`"* claim in `docs/BENCHMARK_RESULTS.md:5` and explicitly state the pinned commit (`8cd5db0` / `2896229`) for WGA/HPRC rows until re-run on `42700e1`.

## F5 -- C5: `docs/VGP_INTEGRATION_PLAN.md` and `manuscript/horizon_paper.md` carry the wrong unit-test count and `9dfec4a` `HEAD` pin
- **Claim:** `C5` (`reports/rounds/round-6/ATTACKED-ARCHITECTURE.md:28`, `reports/rounds/round-6/ATTACKED-VGP_INTEGRATION_PLAN.md:4`, `reports/rounds/round-6/ATTACKED-horizon_paper.md:84`)
- **Verdict:** **NARROWED**
- **Concrete flaw:** `ARCHITECTURE.md:29` (`horizon` dispatcher subcommands matching `src/bin/horizon.rs:7-53`) and `ARCHITECTURE.md:28` (`horizon_resolve` `PackedGraph` + Tarjan biedged snarls + per-window divergence + `W`/`L` updates matching `src/bin/horizon_resolve.rs`) accurately reflect the source at `42700e1`, and `docs/manual.md:150-162` documents `horizon_filter` and `horizon_locus`. However, `docs/VGP_INTEGRATION_PLAN.md:4` and `manuscript/horizon_paper.md:84` both report `118` passing unit tests (`86` lib + `32` bin, `1` ignored) rather than the actual `120` active tests (`86` lib + `34` bin, `1` ignored), and `manuscript/horizon_paper.md:84` labels `9dfec4a` as `HEAD` instead of `42700e1`.
- **Evidence:**
  - `cargo test --all-features` at `42700e1`: `86` lib + `34` bin = `120` active (`1` ignored).
  - `reports/rounds/round-6/ATTACKED-VGP_INTEGRATION_PLAN.md:4` (`118 passing unit tests: 86 lib + 32 bin`); `reports/rounds/round-6/ATTACKED-horizon_paper.md:84` (`All 118 active Rust unit tests (86 library + 32 binary, plus 1 ignored .hal integration test; 119 total) pass on HEAD (9dfec4a)`).
- **Correction:** Update the test counts in `docs/VGP_INTEGRATION_PLAN.md:4` and `manuscript/horizon_paper.md:84` to `120` active (`86` lib + `34` bin, `1` ignored; `121` total) and update `HEAD (9dfec4a)` in `manuscript/horizon_paper.md:84` to `9dfec4a` (or re-verify and pin `42700e1`).

---

## Verdict Table

| Claim | Verdict |
| :--- | :--- |
| C1 | NARROWED |
| C2 | REFUTED |
| C3 | REFUTED |
| C4 | NARROWED |
| C5 | NARROWED |
