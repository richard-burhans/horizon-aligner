# Round 5 ADJUDICATION — Provenance Lens

**Attack reviewed:** `reports/rounds/round-5/attack-provenance.md` (sealed digest `a763d5f04745`, verified against `a763d5f047458a12feea3dccc9ba8f58ac5bbf1188f9e093b21e8bc444ed0d57`).  
**Prior-round context:** `scripts/round history "" --kind adjudication` shows Round 2 `F1` (`UPHELD`, multi-commit provenance conflation in `docs/BENCHMARK_RESULTS.md`) and Round 4 `F7` (`UPHELD`, missing repository provenance for `horizon_merge`, `horizon_guide`, and Table 1). Round 5 repeats the multi-commit conflation pattern across `8cd5db0` (`50dd576`) and `2896229`, alongside stale Figure 6 data and stale post-`8cd5db0` file paths.  
**Method cards:** `refs/methods/` does not exist in this repository; no method card governs these provenance checks.

---

## Finding-by-Finding Rulings

### F1 — C1 (`horizon_gfa stats` / 10k HPRC subgraph in Figure 6A and `98 active` unit-test count)

- **Verdict:** **UPHELD**
- **Diagnosis:** Verified against the live repository. First, `manuscript/scripts/generate_figures.py:765–772` still plots pre-`8cd5db0` values for all four 10,000-segment HPRC `chr22` subgraph runs (`horizon_gfa stats` `0.0256 s` / `22.7 MB` with `vg stats` annotated `"168× slower"`, `horizon_clean` `0.0525 s` / `36.9 MB`, `horizon_wga --graph` `0.821 s` / `404.4 MB`, `horizon_resolve` `8.106 s` / `36.9 MB`) instead of the `8cd5db0` values in `ATTACKED-horizon_paper.md:10, 82` and `docs/BENCHMARK_RESULTS.md:25, 76–79` (`0.08 s` / `22.6 MB` [`53.7×` slower], `0.13 s` / `33.4 MB`, `2.98–6.09 s` / `153.9–154.6 MB`, `5.571 s` / `23.6 MB`). Second, `/home/agent/.cargo/bin/cargo test` on `92e345c` runs **97** passing unit tests (`70` in `src/lib.rs` + `27` across `src/bin/`) and **1** ignored `.hal` test (`98` total), not `98 active` + `1` ignored.
- **Repair & cost:** Both proposed corrections are endorsed. Changing `98 active` to `97 active` in `horizon_paper.md:82` and `C1` is free, and updating `manuscript/scripts/generate_figures.py:765–772` (including `"54× slower"` at line 769) and regenerating `fig6_suite_and_vg.*` requires no re-measurement (the fix is free).
- **What survives:** The `8cd5db0` benchmark numbers and `97` active passing tests (+`1` ignored `.hal` test) in `C1` survive once Figure 6A and the test count are aligned.

---

### F2 — C2 (`8cd5db0` vs. `2896229` multi-commit conflation, `chr6_*` provenance, and stale `minimap2` pending text)

- **Verdict:** **UPHELD**
- **Diagnosis:** Verified against git history (`git show 50dd576`, `git show 2896229 --stat`) and `manuscript/data/chr6_hprc_comparison.json`. Commit `50dd576` (`8cd5db0`) recorded `258.88–267.86 s`, `282,193,564` pairs, `124,160,423` X bp, `109,464,447` Y bp, and `76.44%` exact recall (`110,750` blocks, `88,157,015` pairs). Commit `2896229` modified `partition_colinear_blocks` in `src/formats.rs` and re-ran full `chr6` in `304.53 s` (`docs/BENCHMARK_RESULTS.md:74`), emitting `370,937` blocks and `282,163,046` pairs (`124,160,076` X bp, `109,463,501` Y bp) to reach `99.54%` base-coverage recall, `97.96%` exact recall, and `48.63%` exact pair precision. `ATTACKED-horizon_paper.md:10, 76` and `C2` (`round.json:25`) conflate both commits under `8cd5db0`. Furthermore, no script in `scripts/` or `manuscript/scripts/` generates `manuscript/data/chr6_*`, and `ATTACKED-horizon_paper.md:59, 96` still claim `minimap2` is pending/future work despite `minimap2 -x asm5` appearing at lines `10, 76, 90` and Figure 7.
- **Repair & cost:** Pinning the `99.54%` / `97.96%` / `48.63%` concordance, `370,937` blocks, `282,163,046` pairs, and `304.53 s` runtime to `2896229`, reconciling the `8cd5db0` vs. `2896229` counts in `horizon_paper.md:76`, and removing `minimap2` from the pending/future-work clauses in `horizon_paper.md:59, 96` are all endorsed and free (no re-measurement needed). Adding the `chr6` comparison generation script to `scripts/benchmarks/` is endorsed but not free (requires reconstructing/writing the evaluation script).
- **What survives:** Both the `8cd5db0` runtime (`258.88–267.86 s`) and the `2896229` HPRC/minimap2 concordance (`304.53 s`, `97.96%` exact recall, `99.71%` minimap2 recall) survive once separated by commit and stripped of contradictory pending-tool text.

---

### F3 — C3 (`horizon_guide` `0.859 s` from `2896229` attributed to `8cd5db0`, and stale Figure 6A multi-species values)

- **Verdict:** **UPHELD**
- **Diagnosis:** Verified against `git show 50dd576:docs/BENCHMARK_RESULTS.md` and `git show 2896229 -- docs/BENCHMARK_RESULTS.md`. At `8cd5db0` (`50dd576`), `horizon_guide --progressive --hal --rbest` ran in `1.460 s` with `1,177` MAF blocks (`966` multi-species); the `0.859 s` / `1,109`-block (`967` multi-species) result appeared only in `2896229` after `src/formats.rs` changed. Additionally, `manuscript/scripts/generate_figures.py:770–771` still plots `horizon_guide` at `2.527 s` / `610.1 MB` RSS and `horizon_merge` at `4.451 s` / `26.8 MB` RSS instead of `0.859 s` and `21.78 s` / `26.6 MB`.
- **Repair & cost:** Pinning `horizon_guide` (`0.859 s`, `1,109` blocks) to `2896229` is free. Updating `generate_figures.py:771` (`horizon_merge` to `21.78 s`, `26.6 MB`) is free; updating `generate_figures.py:770` (`horizon_guide` to `0.859 s`) is free on the runtime axis, though note that `docs/BENCHMARK_RESULTS.md:58` records `horizon_guide` Peak RSS as `—` (unmeasured on `2896229`), so plotting `horizon_guide` on the RSS axis requires either measuring its `2896229` Peak RSS or omitting/annotating the RSS coordinate.
- **What survives:** The `horizon_merge` (`21.78 s`, `26.6 MB` at `8cd5db0`) and `horizon_guide` (`0.859 s`, `1,109` blocks at `2896229`) results survive once attributed to their respective commits and reflected in Figure 6A.

---

### F4 — C4 (Figure 6B stale precision for Horse/Macaque and `run_all_evals.py` omitting 3 of 8 Table 1 datasets)

- **Verdict:** **UPHELD**
- **Diagnosis:** Verified against `manuscript/scripts/generate_figures.py:801` and `scripts/benchmarks/run_all_evals.py:56–76`. Line 801 of `generate_figures.py` plots `after = [0.0, 72.3, 86.9, 87.7]` (stale `fb8876f`/`bb43dd9` precision for Horse and Macaque) rather than `[0.0, 67.6, 87.0, 87.7]` reported in `ATTACKED-horizon_paper.md:10, 59, 65, 67` and `docs/BENCHMARK_RESULTS.md:16, 19, 46, 47`. Meanwhile, `run_all_evals.py` only runs 5 of the 8 Table 1 datasets and omits `B10K 363-way`, `cactus241way chr22 1MB`, and `cactus241way chr22_random`.
- **Repair & cost:** Updating `generate_figures.py:801` to `[0.0, 67.6, 87.0, 87.7]` and regenerating `fig6_suite_and_vg.*` is free. Adding `B10K 363-way` and both `cactus241way` slices to `scripts/benchmarks/run_all_evals.py` is endorsed but not free (requires wiring dataset paths/extraction in the runner).
- **What survives:** The `8cd5db0` Table 1 precision/recall numbers in `C4` survive once Figure 6B is synced and the 3 omitted datasets are wired into `run_all_evals.py`.

---

### F5 — C5 (Stale companion date `2026-11-18` and deleted/non-existent file paths `horizon_guide.rs`, `horizon_lift.rs`, `src/hal.rs`)

- **Verdict:** **UPHELD**
- **Diagnosis:** Verified against the live filesystem (`ls src/bin/ src/hal/`). `ATTACKED-horizon_companion.md:5` is dated `"2026-11-18"`; `src/bin/horizon_guide.rs` and `src/bin/horizon_lift.rs` were decomposed into directories `src/bin/horizon_guide/` and `src/bin/horizon_lift/` in `8cd5db0`, yet are still cited in `ATTACKED-horizon_companion.md:164` and `ATTACKED-horizon_paper.md:48, 53`; and `manuscript/scripts/generate_figures.py:586` cites non-existent `src/hal.rs` (actually `src/hal/`).
- **Repair & cost:** Updating the companion date and replacing the three stale paths with `src/bin/horizon_guide/`, `src/bin/horizon_lift/`, and `src/hal/` (plus regenerating `fig_g2_hal_binary_layouts.*`) is endorsed and free.
- **What survives:** The technical specifications in `horizon_companion.md` survive once the date and decomposed module paths are updated.

---

### F6 — C6 (`# EG:Z:` paralog tagging, `partition_colinear_blocks`, and `# HB:Z:` / `HB:i:1` match implementation)

- **Verdict:** **UPHELD**
- **Diagnosis:** Spot-checked `src/formats.rs:157, 161` (`format_gfa_paralog_tags` emitting `# EG:Z:` and `partition_colinear_blocks`), `src/merge.rs:361` (`# EG:Z:`), and `src/bin/horizon_resolve.rs:476` (`EG:Z:` walk tags); all cited symbols and line numbers match the live codebase.
- **Repair & cost:** None required.
- **What survives:** `C6` survives provenance inspection intact.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `generate_figures.py:765-772` plots stale `fb8876f` 10k-subgraph values (`168×` vs `53.7×`) and `cargo test` runs `97` active + `1` ignored (`98` total); both fixes are free |
| F2 | UPHELD | `C2` and `horizon_paper.md:10,76` conflate `8cd5db0` (`50dd576`) and `2896229` (`304.53 s`), `manuscript/data/chr6_*` lacks a generating script, and `horizon_paper.md:59,96` still say `minimap2` is pending; text/pin fixes are free |
| F3 | UPHELD | `horizon_guide` (`0.859 s`, `1,109` blocks) is from `2896229` not `8cd5db0` (`1.460 s`, `1,177` blocks), and `generate_figures.py:770-771` plots stale `2.527 s`/`4.451 s` (note `horizon_guide` `2896229` RSS is unrecorded in `BENCHMARK_RESULTS.md:58`); commit-pin and `horizon_merge` plot fixes are free |
| F4 | UPHELD | `generate_figures.py:801` plots stale Horse (`72.3%` vs `67.6%`) and Macaque (`86.9%` vs `87.0%`) precision (fix is free), and `run_all_evals.py` omits 3 of 8 Table 1 datasets |
| F5 | UPHELD | Companion date (`2026-11-18`) and paths (`src/bin/horizon_guide.rs`, `horizon_lift.rs`, `src/hal.rs`) are stale post-`8cd5db0`; all text/figure-title fixes are free |
| F6 | UPHELD | Attacker's `SURVIVED` finding is verified against `src/formats.rs:157,161`, `src/merge.rs:361`, and `src/bin/horizon_resolve.rs:476` |
