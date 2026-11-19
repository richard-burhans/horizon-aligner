# Round 2 Cold Post-Round Audit

**Repository:** `horizon-aligner` (`/home/rcb112/project/gemini/bench/horizon-aligner`)  
**Pinned Base:** `b18f25b444985a5b7f1db38fd9dd312240230839` (`HEAD`)  
**Scope Audited:** Uncommitted working-tree diff (`git diff b18f25b44498`), `scripts/round` gates (`corrections 2`, `drift 2`, `check 2 --full`), `cargo test` / `cargo clippy`, and targeted code verification in `src/`.

---

## 1. Gate & Tree Status

- **Code & Tests — GREEN:**
  - `cargo test`: `40` passing (`35` in `src/lib.rs`, `1` each in `horizon_clean`, `horizon_gfa`, `horizon_lift`, `horizon_resolve`, `horizon_wga`), `0` failed, `1` ignored (`hal::tests::real_hal_decode_from_env`).
  - `cargo clippy --lib -- -D warnings`: exit `0`.
- **Round Gates — RED (3 failures):**
  1. `scripts/round corrections 2` (exit `1`): `COR-4` is **NOT APPLIED** (`docs/BENCHMARK_RESULTS.md` still contains `horizon_wga --rbest --resolve-dups` at line `58`).
  2. `scripts/round drift 2` (exit `1`): `docs/PANGENOME_SUITE.md` and `docs/STRATEGIC_PLAN.md` are **UNEXAMINED** (`CORRECTED BUT NEVER ATTACKED -- no snapshot exists`).
  3. `scripts/round check 2 --full` (exit `1`): `FAIL the runner wrote BOTH the text under attack and the claim list, and no externally chosen claim is recorded` (`round.json:157-159`: `runner_is_author: true`, `claim_list_author: "runner"`, `human_claim: ""`).

---

## 2. Findings

### F1. `COR-4` Incomplete: Table Row Missed (`docs/BENCHMARK_RESULTS.md:58`)
- **Artifact:** `docs/BENCHMARK_RESULTS.md:58` vs `docs/BENCHMARK_RESULTS.md:22`; `scripts/round corrections 2`.
- **Detail:** `COR-4` updated the Executive Summary bullet at `docs/BENCHMARK_RESULTS.md:22` (removing `--rbest`, adding `pre-0e70693`, and qualifying `24,365,611` as `many-to-many (x, y) aligned pairs` with `non-deduplicated ratio 68.5%`), but left the corresponding table row at `docs/BENCHMARK_RESULTS.md:58` unchanged:
  `| Chr22 Haplotype Alignment | horizon_wga --rbest --resolve-dups | GRCh38.chr22 vs CHM13.chr22:15729460-51324926 | 24,365,611 pairs (9,690 bubbles) | - | 154.92s |`
- **Note on `--rbest`:** `round.json:224` states `"removed unparsed --rbest flag"`, which is inaccurate — `src/bin/horizon_wga.rs:87` does parse `"--rbest"`. The actual issue is whether `--rbest` was used on that pre-`0e70693` run and that line `58` still omits the `many-to-many` and `pre-0e70693` caveats.

### F2. Unattacked Files Modified & Internal Contradictions in `docs/STRATEGIC_PLAN.md`
- **Artifact:** `scripts/round drift 2`; `docs/STRATEGIC_PLAN.md:3, 23, 31`; `src/bin/horizon_resolve.rs:346-347`.
- **Detail:** `round.json:111-119` only snapshotted `docs/BENCHMARK_RESULTS.md`, yet `COR-2` edited `docs/STRATEGIC_PLAN.md` and `COR-6` (plus additional edits) modified `docs/PANGENOME_SUITE.md`. In `docs/STRATEGIC_PLAN.md`:
  1. Line `3` (updated by `COR-2`) states: *"Only `horizon_gfa stats` supports single-pass `"-"` `stdin` streaming; two-pass tools (`horizon_clean`, `horizon_resolve`, `horizon_lift`) and `horizon_wga` require seekable file paths."*
  2. Line `23` (Section 3.B, untouched) still states: *"Support `stdin`/`stdout` GFA/MAF streaming on linear filter/transform tools (`horizon_clean`, `horizon_resolve`, `horizon_lift`)..."* — contradicting line `3`.
  3. Line `31` (Section 3.D, untouched) still describes HPRC `GRCh38` vs. `CHM13` evaluation as *"currently unrun"*, contradicting `docs/BENCHMARK_RESULTS.md:22-26, 58-62`.
  4. In code, `src/bin/horizon_resolve.rs:346-347` still accepts `"-"` (`io::stdin()`) in `open_reader()` despite `horizon_resolve` being a two-pass tool; passing `"-"` drains `stdin` on Pass 1 and hits EOF on Pass 2.

### F3. Derived-Quantity Denominator Mixing (`docs/BENCHMARK_RESULTS.md:26`)
- **Artifact:** `docs/BENCHMARK_RESULTS.md:26`.
- **Detail:** Line `26` reports `5,863.3 MB peak RSS (0.293 GB/Mbp)` for the `20 Mbp` `chr6` q-arm (`5.8633 / 20 = 0.293` **ref-only**), immediately followed by the cancelled full-`chr6` run reporting `0.094 GB/Mbp combined, 0.155 GB/Mbp ref-only`. Because `0.293 GB/Mbp` is not labeled `ref-only` (combined would be `0.147 GB/Mbp`), the two rates in the same bullet use different denominators without an explicit label on the first.

### F4. Code-Level Defects Remain Open (Doc-Only Round)
- **Artifact:** `git diff --stat b18f25b44498` (only `docs/BENCHMARK_RESULTS.md`, `docs/PANGENOME_SUITE.md`, `docs/STRATEGIC_PLAN.md` changed; `0` lines changed in `src/`).
- **Detail:** All Round 2 corrections (`COR-1`..`COR-7`) were documentation-only. The underlying code defects identified in the round remain in `src/`:
  - `mdl_xdrop_extend` (`src/mdl/dp.rs`) uses static defaults and `horizon_wga --graph` passes `&HmmParams::default()`.
  - `horizon_clean` (`src/bin/horizon_clean.rs:250`) strips all `S`-line optional tags (`S\t{}\t{}`).
  - `horizon_resolve` (`src/bin/horizon_resolve.rs:346-347`) accepts `"-"` `stdin` in a two-pass reader and leaves `W`-line `start`/`end` and `L` links unmodified when collapsing duplications.
  - `horizon_merge` (`src/bin/horizon_merge.rs:25`) loads all input MAFs into memory via `fs::read_to_string` and has `0` unit tests.

### F5. Minor Formatting Inconsistency (`docs/BENCHMARK_RESULTS.md:23, 24, 63`)
- **Artifact:** `docs/BENCHMARK_RESULTS.md:23, 24, 63` vs `docs/PANGENOME_SUITE.md:12`.
- **Detail:** `docs/BENCHMARK_RESULTS.md` writes `` `sup>0` `` (no spaces), whereas `docs/PANGENOME_SUITE.md:12` writes `` `sup > 0` ``.

---

## 3. What This Cold Audit Did Not Reach

- **No benchmark re-execution:** None of the HPRC 12 GB graph runs, `chr22`/`chr6` alignments, `horizon_lift`, `horizon_resolve`, or `cactus241way`/`multiz` WGA benchmarks were re-run; runtime, RSS, and accuracy figures in `docs/BENCHMARK_RESULTS.md` were checked only for internal arithmetic and code correspondence.
- **No full read of `attack-L1`..`L7` or `adjudication-L1`..`L7`:** Per `AUDIT-BRIEF.md`, the cold read was formed from the working-tree diff, `scripts/round` diagnostics, `cargo test`/`clippy`, `round.json`, and `adjudication-completeness.md` without reconstructing the seven individual lens reports.
