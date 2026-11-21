# Round 6 Cold Post-Round Audit

**Auditor:** Cold Auditor (`aeebc041-88f7-4b19-bdbe-5031442991ff`)
**Date:** 2026-11-21
**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `42700e1`, working-tree repairs applied)
**Overall Verdict:** **FAIL** (code tests & clippy pass, but `scripts/round check 6 --full`, `scripts/round drift 6`, `cargo fmt --check`, and 3 doc-sync checks fail)

---

## 1. Gate & Command Verification

| Check | Command | Result | Evidence |
| :--- | :--- | :--- | :--- |
| Unit tests | `cargo test --all-features` | **PASS** | `120` active passed (`86` `src/lib.rs` + `34` `src/bin/*`), `0` failed, `1` ignored (`hal::tests::real_hal_decode_from_env`); `121` total. |
| Lints | `cargo clippy --all-targets --all-features -- -D warnings` | **PASS** | Exit `0`, zero warnings. |
| Formatting | `cargo fmt --check` | **FAIL** | Exit `1` (formatting diffs in Rust files including `src/bin/horizon_clean.rs:276,290,790`; tracked as `R6-FMT`). |
| Corrections | `scripts/round corrections 6` | **PASS** | `8 correction(s); 0 not applied or malformed` (`COR-1`..`COR-8` all `APPLIED`). |
| Drift | `scripts/round drift 6` | **FAIL** | Exit `1`. Flags `docs/manual.md` as `UNEXAMINED` (`CORRECTED BUT NEVER ATTACKED`) and `README.md`, `docs/PANGENOME_SUITE.md`, `docs/BENCHMARK_RESULTS.md`, `manuscript/horizon_paper.md` as `DRIFTED` (quoted claims changed post-snapshot). |
| Round ledger | `scripts/round check 6 --full` | **FAIL** | Exit `1`: `expected 28 outstanding, derived 29`. `scripts/round balance` reports `r6: expected 28, derived 29 (+1)`. |

---

## 2. Diff & Code Verification Findings

### Verified Accurate in Working-Tree Diff
- **Unit test counts (`120` active = `86` lib + `34` bin, `1` ignored; `121` total):** Verified against `cargo test --all-features` output (`docs/PANGENOME_SUITE.md:3`, `docs/STRATEGIC_PLAN.md:3`, `docs/VGP_INTEGRATION_PLAN.md:4`, `manuscript/horizon_paper.md:83`).
- **`horizon_resolve` (`src/bin/horizon_resolve.rs:458,476,843`):** Confirmed topological `[0.02, 0.15]` divergence heuristic (`(0.02 + 0.08 * non_anchor_frac + 0.05 * inv_frac).clamp(0.01, 0.25)`) and `0M` overlap assumption (`"L\t{}\t{}\t{}\t{}\t0M"`).
- **`horizon_filter` (`src/bin/horizon_filter.rs:56,123,167-168`):** Confirmed Pair-HMM MDL scoring uses `base_probs(0.5)` or contig-global `base_probs(qs.gc)` / `base_probs(ts.gc)`, with no BOCPD segmentation.
- **`horizon_locus` (`src/bin/horizon_locus.rs:67-71,130,141-142`):** Confirmed `jaccard_cov`, `p_dist`, and `jc69_dist` over aligned bases with `p_dist.clamp(0.0, 0.74)`.
- **`src/io_util.rs:92-105`:** Confirmed `.zst`/`.zstd` streaming spawns `zstd` CLI (`Command::new("zstd")`).
- **`docs/BENCHMARK_RESULTS.md:62`:** Confirmed `9,697 / 25,000` bp = `38.788%` (`38.8%`).

### Failures / Residual Defects Found
1. **Stale sibling number in `docs/VGP_INTEGRATION_PLAN.md:76`:**
   - The P4 row still reads `` `chr6` `.paf.zst` → `.paf.zst` streaming filter in `0.49s` ``, whereas `README.md:124`, `docs/BENCHMARK_RESULTS.md:58`, and `manuscript/horizon_paper.md:83` were updated to `0.50s` (`565` records / `110.00 Mbp`).
2. **Unmoved stale claims in `ARCHITECTURE.md:11,28`:**
   - `ARCHITECTURE.md` was in `attacked` and quoted by `C5`, yet remained unmoved when `COR-1`/`COR-7`/`COR-8` repaired the other docs:
     - `ARCHITECTURE.md:11` still describes `horizon_filter` as `(streaming MDL+BOCPD+SDUST PAF filter & netter)`.
     - `ARCHITECTURE.md:28` still states `horizon_resolve` `estimates per-window divergence (unless --divergence is passed)` without noting the topological `[0.02, 0.15]` heuristic or `0M` overlap assumption.
3. **Unmoved stale phrase in `docs/STRATEGIC_PLAN.md:49`:**
   - Still states `per-window divergence estimation` for `horizon_resolve.rs` (even though `docs/STRATEGIC_PLAN.md:3` was updated in the working tree).
4. **Untracked / unattacked working-tree edits:**
   - `docs/manual.md` was edited by `COR-7` without ever being in `attacked` (`scripts/round drift 6` flags `UNEXAMINED`).
   - `docs/STRATEGIC_PLAN.md`, `manuscript/references.bib`, and `manuscript/horizon_paper.pdf` are modified in `git status -s` outside `round.json`'s `attacked` and `drift_reviewed` lists.
5. **Round 6 ledger mismatch (`scripts/round check 6 --full`):**
   - `scripts/round balance` reports `r6: expected 28, derived 29 (+1)`, causing `scripts/round check 6 --full` to fail.

---

## 3. Explicit Coverage Statement (What Was Not Reached)

- **Did not read** `attack-*.md`, `adjudication-*.md`, or `completeness-critic.md` (to preserve cold-read independence).
- **Did not re-run** external/large-data benchmarks (`chr6` FastGA/HPRC/minimap2, `10k` HPRC subgraph/`vg`, Koala × Amphioxus, avian 2-way, human 4-way, or 92-species `chrM` / 5-mammal `horizon_guide`) or inspect binary `manuscript/horizon_paper.pdf`.
