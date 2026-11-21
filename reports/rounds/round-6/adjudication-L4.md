# Round 6 ADJUDICATION of L4 -- Provenance Lens

**Target revision:** `42700e142b70996c23cdf391ac29e1e5bcf4640d` (`42700e1`)
**Target attack:** `reports/rounds/round-6/attack-L4-provenance.md` (SHA-256 `f678dd7df338dad7146607fddc71f0afa256c2e27fd9dbcf239d823348f3e768`, matches sealed `f678dd7df338`)
**Method cards & prior rounds:** `refs/methods/` does not exist in this repository. `scripts/round history "provenance" --kind adjudication` and `scripts/round history "scripts/benchmarks" --kind adjudication` confirm Round 2 (`F1`), Round 4 (`F6`), and Round 5 (`F2`, `F4` in L4; `F3`–`F5` in L7) previously adjudicated multi-commit conflation and `/tmp/`-only benchmark provenance.

---

## F1 -- C1: Unit test count is 120 active (86 lib + 34 bin), not 118 (86 lib + 32 bin)
- **Ruling:** **UPHELD** (attacker's `NARROWED` diagnosis on C1 stands).
- **What survives:** The C1 implementation inventory (`horizon_gfa` through `MdlCostModel`) and binary count (`11` suite binaries + `bench_chr1` = `12` entries in `src/bin/`) are verified at `42700e1`. The unit test count (`118`: `86` lib + `32` bin, `1` ignored) is refuted.
- **Independent verification:** `cargo test --all-features` at `42700e1` runs `86` lib tests (`1` ignored) and `34` bin tests (`bench_chr1` `1`, `horizon` `2`, `horizon_clean` `6`, `horizon_eval` `2`, `horizon_filter` `1`, `horizon_gfa` `2`, `horizon_guide` `4`, `horizon_lift` `11`, `horizon_locus` `1`, `horizon_merge` `1`, `horizon_resolve` `1`, `horizon_wga` `2`) = **120** active (`121` total). `cargo clippy --all-targets --all-features -- -D warnings` passes with `0` warnings. `git diff --stat 9dfec4a..42700e1` touches only docs and `manuscript/horizon_paper.pdf`, so `118` was already stale at `9dfec4a`.
- **Repair & cost:** Endorsed; **the fix is free** (text-only update of `118`/`32` to `120`/`34`, and `119 total` to `121 total` in `manuscript/horizon_paper.md:84`, across `docs/PANGENOME_SUITE.md:3`, `docs/VGP_INTEGRATION_PLAN.md:4`, `docs/STRATEGIC_PLAN.md:3`, and `manuscript/horizon_paper.md:84`, plus PDF rebuild).

## F2 -- C2: All `horizon_filter` VGP/chr6 artifacts in `/tmp/` predate `9dfec4a` and have zero tracked scripts or datasets
- **Ruling:** **UPHELD**.
- **What survives:** Nothing of C2's `9dfec4a` / `HEAD` provenance claim survives. Every `horizon_filter` `/tmp/` artifact predates `9dfec4a` (`2026-11-20 19:38:58 -0500`, which modified `src/bin/horizon_filter.rs` and `src/io_util.rs`), and no tracked script or dataset entry reproduces C2.
- **Independent verification:** Confirmed `stat -c "%y %n"`: `/tmp/chr6_fastga_mdl_net.paf` (`19:10:09 -0500`, `553` lines vs. `1,004`/`848` claimed), `/tmp/chr6_filtered.paf.zst` (`19:19:58 -0500`), `/tmp/vgp_test/horizon_net.paf.gz` (`19:30:26 -0500`). `grep -rnE "(horizon_filter|koala|amphioxus|GCA_003287225|GCA_015852565|FastGA_11)" scripts/ docs/BENCHMARK_DATASETS.md` returns zero matches.
- **Repair & cost:** Endorsed. Re-running `horizon_filter` benchmarks on `42700e1` and committing tracked scripts/dataset provenance in `scripts/benchmarks/` and `docs/BENCHMARK_DATASETS.md` is **not free**; alternatively, narrowing `docs/BENCHMARK_RESULTS.md:5,81`, `README.md:123`, `docs/VGP_INTEGRATION_PLAN.md:4,68`, and `manuscript/horizon_paper.md:84` to mark them as unverified pre-`9dfec4a` `/tmp/` runs is **free**.

## F3 -- C3: All `horizon_locus` outputs in `/tmp/` predate `9dfec4a` (avian/human predating `4dba79e` RC-complement fix) and have zero tracked scripts
- **Ruling:** **UPHELD** (diagnosis confirmed and strengthened by live re-run).
- **What survives:** Nothing of C3's `9dfec4a` / `HEAD` provenance survives, and the human 4-way metrics are actively stale at `42700e1`.
- **Independent verification:** Confirmed `/tmp/locus_avian_out.similarity.tsv` (`19:15:46 -0500`) and `/tmp/locus_human_out.similarity.tsv` (`19:16:09 -0500`) predate `4dba79e` (`19:25:28 -0500`, which fixed `complement_base(seq_y[y])` in `src/bin/horizon_locus.rs`), and `/tmp/vgp_test/locus_out.similarity.tsv` (`19:31:17 -0500`) predates `9dfec4a` (`19:38:58 -0500`). Re-running `./target/release/horizon_locus` at `42700e1` on `/tmp/locus_avian.fa` reproduces `231,444` aligned bp, `p_dist = 0.192552`, `jc69_dist = 0.222528`, and `1,364` MAF blocks (with `jaccard_cov` shifting from `0.458662` to `0.625374`), whereas `/tmp/locus_human_25k.fa` produces `102` MAF blocks (`92` multi-sp) instead of `29` (`20` multi-sp) and `grch38_mhc` vs `hla_dup` `p_dist = 0.041539` (`jc69_dist = 0.042734`) instead of `0.124507` (`jc69_dist = 0.136150`). Zero matches in `scripts/` and `docs/BENCHMARK_DATASETS.md`.
- **Repair & cost:** Endorsed; **not free** (requires committing reproducible `horizon_locus` benchmark scripts/dataset provenance and updating the stale human 4-way and avian `jaccard_cov` figures from `42700e1` runs).

## F4 -- C4: `docs/BENCHMARK_RESULTS.md:5` falsely claims `HEAD` re-verification across 26 post-`8cd5db0` commits, and tracked scripts still hardcode `/tmp/`
- **Ruling:** **UPHELD** (attacker's `NARROWED` diagnosis stands, with a minor commit-count correction).
- **What survives:** The `8cd5db0` / `2896229` git pins and the four tracked `chr6` HPRC evaluation scripts in `scripts/benchmarks/` survive; `docs/BENCHMARK_RESULTS.md:5`'s blanket *"All benchmarks below were re-verified on `HEAD`"* claim is refuted.
- **Independent verification:** `git log --oneline 8cd5db0..42700e1` lists **28** commits (the attacker's `26` is `8cd5db0..9dfec4a`), including code changes to `src/mdl/`, `src/bocpd.rs`, `src/egraph.rs`, `src/formats.rs`, `src/seeding/`, `horizon_wga.rs`, and `horizon_eval.rs`. Confirmed `/tmp/` default paths in `compare_chr6_horizon_vs_hprc.py:9-14`, `extract_hprc_chr6_truth.py:8`, `score_paf_vs_hprc.py:10-16`, `plot_chr6_hprc_fig.py:43-44`, and `bench_profile_ledger.py:103,106`, and only a single `62402b8` `chr6_peri_5m` entry in `scripts/benchmarks/profile_ledger.jsonl`.
- **Repair & cost:** Endorsed; **the fix is free** for removing the blanket *"All benchmarks below were re-verified on `HEAD`"* sentence in `docs/BENCHMARK_RESULTS.md:5` and pinning WGA/HPRC rows explicitly to `8cd5db0` / `2896229`. (Re-running full `chr6` on `42700e1` or refactoring `/tmp/` defaults is not free.)

## F5 -- C5: `docs/VGP_INTEGRATION_PLAN.md` and `manuscript/horizon_paper.md` carry the wrong unit-test count and `9dfec4a` `HEAD` pin
- **Ruling:** **UPHELD** (attacker's `NARROWED` diagnosis stands; overlaps F1 on test counts).
- **What survives:** `ARCHITECTURE.md:28` (`horizon_resolve` `PackedGraph` + Tarjan biedged snarls + per-window divergence + `P`/`W`/`L` updates) and `ARCHITECTURE.md:29` match `42700e1` source. `docs/VGP_INTEGRATION_PLAN.md:4` and `manuscript/horizon_paper.md:84` carry the stale `118` (`86` lib + `32` bin) test count, and `manuscript/horizon_paper.md:84` mislabels `9dfec4a` as `HEAD`.
- **Independent verification:** Verified `ARCHITECTURE.md:28-29`, `docs/VGP_INTEGRATION_PLAN.md:4`, `manuscript/horizon_paper.md:84`, and `cargo test --all-features` (`120` active: `86` lib + `34` bin, `1` ignored; `121` total).
- **Repair & cost:** Endorsed; **the fix is free** (update `118`/`32`/`119` to `120`/`34`/`121` and `HEAD (9dfec4a)` to `9dfec4a` or `42700e1` in `docs/VGP_INTEGRATION_PLAN.md:4` and `manuscript/horizon_paper.md:84`, plus PDF rebuild).

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `cargo test --all-features` at `42700e1` runs `120` active tests (`86` lib + `34` bin, `1` ignored), not `118` (`86` + `32`); free text fix endorsed |
| F2 | UPHELD | All `/tmp/` `horizon_filter` artifacts predate `9dfec4a` (`19:10`–`19:30` vs `19:38:58`) with zero tracked scripts or dataset entries in `scripts/` or `docs/BENCHMARK_DATASETS.md` |
| F3 | UPHELD | All `/tmp/` `horizon_locus` outputs predate `9dfec4a` (avian/human predate `4dba79e` RC fix); live re-run at `42700e1` shifts human 4-way `p_dist` from `0.124507` to `0.041539` and MAF blocks from `29` to `102` |
| F4 | UPHELD | `docs/BENCHMARK_RESULTS.md:5` claims `HEAD` re-verification across `28` post-`8cd5db0` commits (`26` to `9dfec4a`), and `scripts/benchmarks/*.py` still default to `/tmp/`; free wording fix endorsed |
| F5 | UPHELD | `ARCHITECTURE.md:28` matches `42700e1`, while `docs/VGP_INTEGRATION_PLAN.md:4` and `manuscript/horizon_paper.md:84` carry `118` tests and `manuscript/horizon_paper.md:84` mislabels `9dfec4a` as `HEAD`; free fix endorsed |
