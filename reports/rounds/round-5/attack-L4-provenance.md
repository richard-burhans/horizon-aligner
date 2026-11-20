# Round 5 Attack — Provenance Lens

**Target round:** Round 5 (`reports/rounds/round-5/round.json`, `HEAD` = `92e345c1fa938c30bff286ac2fd956adc9a2b4d6`)  
**Attacked snapshots:** `reports/rounds/round-5/ATTACKED-horizon_paper.md`, `reports/rounds/round-5/ATTACKED-horizon_companion.md`  
**Method cards & literature note:** This project keeps no method cards (`refs/methods/` does not exist; confirmed in `BRIEF-provenance.md:142`). All findings below were verified directly against repository commits, `docs/BENCHMARK_RESULTS.md`, `manuscript/data/`, `manuscript/scripts/generate_figures.py`, `scripts/benchmarks/run_all_evals.py`, and `cargo test` on `92e345c`.  
**Prior-round context:** Round 2 `F1` (`UPHELD`) flagged `docs/BENCHMARK_RESULTS.md:5` for claiming uniform provenance across three commits, and Round 4 `F7` (`UPHELD`) flagged `horizon_merge`, `horizon_guide`, and Table 1 for lacking repository provenance. In Round 5, benchmarks were refreshed across `8cd5db0` (`50dd576`) and `2896229`, but multi-commit conflation, stale Figure 6 data, and stale file paths remain.

---

## F1 — C1 (`NARROWED`): Figure 6A plots stale `fb8876f`/`bb43dd9` 10k-subgraph numbers and `horizon_paper.md:82` overcounts active unit tests (`97` active + `1` ignored = `98` total)

- **Claim:** `C1` (`round.json:13`, `ATTACKED-horizon_paper.md:10, 82`)
- **Concrete flaws:**
  1. **Figure 6A (`manuscript/scripts/generate_figures.py:765–772`, rendered to `manuscript/figures/fig6_suite_and_vg.pdf`) contradicts `ATTACKED-horizon_paper.md:82` and `docs/BENCHMARK_RESULTS.md:25, 76–79` across all four Horizon 10,000-segment HPRC `chr22` subgraph runs** because `generate_figures.py` was not updated from the pre-`8cd5db0` (`fb8876f`/`bb43dd9`) values:
     - `horizon_gfa stats`: `0.0256 s`, `22.7 MB` RSS in `generate_figures.py:765` (with `vg stats` annotated `"168× slower"` at `line 769`) vs. `0.08 s`, `22.6 MB` RSS in `C1`, `ATTACKED-horizon_paper.md:10, 82`, and `docs/BENCHMARK_RESULTS.md:25, 76` (which makes `vg stats` `4.2984 / 0.08 = 53.7×` slower, not `168×`).
     - `horizon_clean`: `0.0525 s`, `36.9 MB` RSS in `generate_figures.py:766` vs. `0.13 s`, `33.4 MB` RSS in `ATTACKED-horizon_paper.md:82` and `docs/BENCHMARK_RESULTS.md:25, 77`.
     - `horizon_wga --graph`: `0.821 s`, `404.4 MB` RSS in `generate_figures.py:767` vs. `2.98–6.09 s`, `153.9–154.6 MB` RSS in `ATTACKED-horizon_paper.md:82` and `docs/BENCHMARK_RESULTS.md:25, 79`.
     - `horizon_resolve`: `8.106 s`, `36.9 MB` RSS in `generate_figures.py:772` vs. `5.571 s`, `23.6 MB` RSS in `ATTACKED-horizon_paper.md:82` and `docs/BENCHMARK_RESULTS.md:25, 78`.
  2. **Unit-test count off-by-one:** `round.json:13` (`C1`) and `ATTACKED-horizon_paper.md:82` state `"All 98 active Rust unit tests (plus 1 environment-gated .hal integration test) pass on HEAD with zero clippy warnings."` Running `/home/agent/.cargo/bin/cargo test` on `HEAD` (`92e345c`) runs **97** passing tests (`70` lib + `27` across `src/bin/`) and **1** ignored test (`cargo test -- --list` lists `98` tests total), not `98` active + `1` ignored.
- **Required correction:** Update `manuscript/scripts/generate_figures.py:765–772` (and regenerate `fig6_suite_and_vg.*`) to match the `8cd5db0` values (`0.08 s` / `22.6 MB`, `0.13 s` / `33.4 MB`, `2.98–6.09 s` / `153.9–154.6 MB`, `5.571 s` / `23.6 MB`, and replace `"168× slower"` with `"54× slower"`), and change `98 active` to `97 active` in `ATTACKED-horizon_paper.md:82` and `C1`.

---

## F2 — C2 (`NARROWED`): Commit `8cd5db0` and `2896229` results are conflated, `chr6_*.json`/`.npz` lack a generating script, and `horizon_paper.md:59, 96` still claim `minimap2` is pending

- **Claim:** `C2` (`round.json:25`, `ATTACKED-horizon_paper.md:10, 59, 76, 96`)
- **Concrete flaws:**
  1. **Unacknowledged multi-commit mixing (`8cd5db0` vs. `2896229`):** `round.json:25` (`C2`) and `ATTACKED-horizon_paper.md:10` pin the `GRCh38.chr6 × CHM13.chr6` benchmark and HPRC concordance to `8cd5db0`. However, `8cd5db0` did not contain the `partition_colinear_blocks` fix in `src/formats.rs`—that fix was committed in `2896229` (`fix(formats): preserve overlapping colinear MAF chains`, two commits after `8cd5db0`):
     - Commit `50dd576` (the `8cd5db0` benchmark refresh) recorded `horizon_wga --resolve-dups` in `258.88–267.86 s` with `282,193,564` pairs (`124,160,423` X bp, `109,464,447` Y bp) and only `76.44%` exact recall (`84,301,290` pairs on shared bases).
     - Commit `2896229` modified `src/formats.rs` and re-ran full `chr6` in **`304.53 s`** (`docs/BENCHMARK_RESULTS.md:74`), emitting `370,937` blocks and `282,163,046` pairs (`124,160,076` X bp, `109,463,501` Y bp in `manuscript/data/chr6_hprc_comparison.json:9–12`) to achieve `99.54%` base-coverage recall, `97.96%` exact coordinate recall, and `48.63%` exact pair precision.
     - `ATTACKED-horizon_paper.md:76` mixes both commits in adjacent sentences without noting the commit change—reporting `8cd5db0`'s `258.88–267.86 s`, `282,193,564` pairs, `124,160,423` X bp, and `109,464,447` Y bp (omitting `304.53 s`), followed immediately by `2896229`'s `370,937` blocks and `282,163,046` pairs. Meanwhile, `ATTACKED-horizon_paper.md:10` and `C2` merge `258.88–304.53 s` and `99.54%` / `97.96%` under `(8cd5db0)`.
  2. **Missing provenance script for `manuscript/data/chr6_*`:** `manuscript/data/chr6_hprc_comparison.{json,npz}` and `chr6_minimap2_{summary.json,data.npz}` (added in `2896229`) are only read by `manuscript/scripts/generate_figures.py:831–834`; no script in `scripts/` or `manuscript/scripts/` documents or reproduces how they were generated from the MAF/PAF/GFA inputs.
  3. **Stale "minimap2 pending" text:** `ATTACKED-horizon_paper.md:59` still states `"head-to-head external pairwise aligner re-runs (minimap2, FastGA, LASTZ) remain pending external tool availability on the benchmark host"` and `ATTACKED-horizon_paper.md:96` lists `minimap2` as future work, contradicting `ATTACKED-horizon_paper.md:10, 76, 90` and Figure 7 where `minimap2 -x asm5` (`v2.28-r1209`) is reported.
- **Required correction:** Pin the `99.54%` / `97.96%` / `48.63%` concordance, `370,937` blocks, `282,163,046` pairs, and `304.53 s` runtime to `2896229` (not `8cd5db0`), reconcile the `8cd5db0` vs. `2896229` pair/bp counts in `ATTACKED-horizon_paper.md:76`, remove `minimap2` from the pending/future-work clauses in `ATTACKED-horizon_paper.md:59, 96`, and add the `chr6` comparison generation script to `scripts/benchmarks/`.

---

## F3 — C3 (`NARROWED`): `horizon_guide` `0.859 s` (`1,109` blocks) is from `2896229` rather than `8cd5db0`, and Figure 6A plots stale `fb8876f`/`bb43dd9` values

- **Claim:** `C3` (`round.json:37`, `ATTACKED-horizon_paper.md:10, 82`)
- **Concrete flaws:**
  1. **Commit-pinning mismatch (`8cd5db0` vs. `2896229`):** `round.json:37` (`C3`) and `ATTACKED-horizon_paper.md:10` attribute `horizon_guide --progressive --hal --rbest` (`0.859 s`, `1,109` MAF blocks, `967` multi-species) to `8cd5db0`. In `git show 50dd576` (the `8cd5db0` benchmark refresh), `horizon_guide` ran in **`1.460 s`** with **`1,177` MAF blocks (`966` multi-species)**. The `0.859 s` / `1,109`-block numbers only appeared in **`2896229`** after `src/formats.rs` was modified.
  2. **Stale Figure 6A values:** `manuscript/scripts/generate_figures.py:770–771` (`manuscript/figures/fig6_suite_and_vg.pdf`) still plots the older `fb8876f`/`bb43dd9` values for both multi-species tools: `horizon_guide` at **`2.527 s` / `610.1 MB` RSS** (vs. `0.859 s` in `C3` and `ATTACKED-horizon_paper.md:10, 82`) and `horizon_merge` at **`4.451 s` / `26.8 MB` RSS** (vs. `21.78 s` / `26.6 MB` in `C3` and `ATTACKED-horizon_paper.md:10, 82`).
- **Required correction:** Pin `horizon_guide` (`0.859 s`, `1,109` blocks) to `2896229` (or `92e345c`) rather than `8cd5db0`, and update `manuscript/scripts/generate_figures.py:770–771` to match `21.78 s` (`horizon_merge`) and `0.859 s` (`horizon_guide`).

---

## F4 — C4 (`NARROWED`): Figure 6B plots stale precision for Horse (`72.3%` vs. `67.6%`) and Macaque (`86.9%` vs. `87.0%`), and `run_all_evals.py` omits 3 of 8 Table 1 datasets

- **Claim:** `C4` (`round.json:50`, `ATTACKED-horizon_paper.md:10, 59–70`)
- **Concrete flaws:**
  1. **Figure 6B vs. Table 1 / `BENCHMARK_RESULTS.md` mismatch:** `manuscript/scripts/generate_figures.py:801` (`manuscript/figures/fig6_suite_and_vg.pdf`, Panel B) plots post-gate `HEAD` precision as **`72.3%`** for Horse and **`86.9%`** for Macaque (stale values from `fb8876f`/`bb43dd9`). Both `ATTACKED-horizon_paper.md:10, 59, 65, 67` and `docs/BENCHMARK_RESULTS.md:16, 19, 46, 47` report **`67.6%`** for Horse and **`87.0%`** for Macaque.
  2. **Incomplete evaluation harness:** `scripts/benchmarks/run_all_evals.py:56–76` only runs 5 of the 8 Table 1 datasets (`multiz100way` `mm10`/`panTro4`, `cactus447way` horse/macaque, `multiz470way` marmoset) and omits `B10K 363-way`, `cactus241way chr22 1MB`, and `cactus241way chr22_random`, leaving those three Table 1 rows without an in-repo reproducible runner.
- **Required correction:** Update `manuscript/scripts/generate_figures.py:801` to `[0.0, 67.6, 87.0, 87.7]` and regenerate `fig6_suite_and_vg.*`; add `B10K 363-way` and both `cactus241way` slices to `scripts/benchmarks/run_all_evals.py`.

---

## F5 — C5 (`NARROWED`): Stale file paths (`horizon_guide.rs`, `horizon_lift.rs`, `src/hal.rs`) and companion date (`2026-11-18`) after `8cd5db0` module decomposition

- **Claim:** `C5` (`round.json:60`, `ATTACKED-horizon_companion.md:5, 164`)
- **Concrete flaws:**
  1. `ATTACKED-horizon_companion.md:5` is dated `"2026-11-18"` (predating `8cd5db0` and `horizon_paper.md`'s `"2026-11-19"` / Round 5's `"2026-11-20"`).
  2. `ATTACKED-horizon_companion.md:164` references `src/bin/horizon_guide.rs`, which was deleted in `8cd5db0` and split into `src/bin/horizon_guide/main.rs` and `src/bin/horizon_guide/rescue.rs` (likewise `ATTACKED-horizon_paper.md:48, 53` still cites deleted `src/bin/horizon_lift.rs` and `src/bin/horizon_guide.rs`).
  3. `manuscript/scripts/generate_figures.py:586` (Figure G2 title) cites `src/hal.rs`, which does not exist (`src/hal/mod.rs` and submodules).
- **Required correction:** Update the companion date and replace `src/bin/horizon_guide.rs`, `src/bin/horizon_lift.rs`, and `src/hal.rs` with `src/bin/horizon_guide/`, `src/bin/horizon_lift/`, and `src/hal/` in `ATTACKED-horizon_companion.md:5, 164`, `ATTACKED-horizon_paper.md:48, 53`, and `manuscript/scripts/generate_figures.py:586`.

---

## F6 — C6 (`SURVIVED`): `# EG:Z:` paralog tagging, `partition_colinear_blocks`, and `# HB:Z:` / `HB:i:1` match the current Rust implementation

- **Claim:** `C6` (`round.json:70`, `ATTACKED-horizon_companion.md:51–52`, `ATTACKED-horizon_paper.md:40, 44, 76`)
- **Evidence checked:** Verified in `src/formats.rs:157, 161` (`format_gfa_paralog_tags` emitting `# EG:Z:` and `partition_colinear_blocks`), `src/merge.rs:361` (`# EG:Z:`), and `src/bin/horizon_resolve.rs:476` (`EG:Z:` walk tags). No provenance discrepancy found within the scope of `C6`.

---

## Verdict Summary

| Claim | Verdict | Primary Provenance Defect |
| :--- | :--- | :--- |
| **C1** | **NARROWED** | Figure 6A plots stale `fb8876f`/`bb43dd9` 10k-subgraph numbers; `98 active` unit tests overcounts `97 active + 1 ignored` (`98` total) |
| **C2** | **NARROWED** | Conflates `8cd5db0` (`258.88–267.86 s`, `282,193,564` pairs) and `2896229` (`304.53 s`, `282,163,046` pairs, `97.96%` recall) under `8cd5db0`; no script generates `manuscript/data/chr6_*`; `horizon_paper.md:59, 96` still say `minimap2` is pending |
| **C3** | **NARROWED** | Attributes `horizon_guide` (`0.859 s`, `1,109` blocks) from `2896229` to `8cd5db0` (`1.460 s`, `1,177` blocks); Figure 6A plots stale `2.527 s` / `4.451 s` |
| **C4** | **NARROWED** | Figure 6B plots stale precision for Horse (`72.3%` vs. `67.6%`) and Macaque (`86.9%` vs. `87.0%`); `run_all_evals.py` omits 3 of 8 Table 1 datasets |
| **C5** | **NARROWED** | Companion date (`2026-11-18`) and file paths (`src/bin/horizon_guide.rs`, `horizon_lift.rs`, `src/hal.rs`) are stale post-`8cd5db0` |
| **C6** | **SURVIVED** | Verified `# EG:Z:` / `EG:Z:` and `partition_colinear_blocks` in `src/formats.rs`, `src/merge.rs`, `src/bin/horizon_resolve.rs` |
