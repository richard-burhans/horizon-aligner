# Round 4 Cold Post-Round Audit

**Verdict:** PASS

**Scope & Baseline:**
- Repository: `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` `33bede8`, `M manuscript/horizon_paper.md`, `?? reports/rounds/round-4/`)
- Audited diff: `git diff manuscript/horizon_paper.md` against `docs/BENCHMARK_RESULTS.md` and `src/`

## 1. Build & Test Gate
All checks passed on `HEAD` (`33bede8` + working-tree `manuscript/horizon_paper.md` edits):
- `cargo fmt --check`: exited `0`.
- `cargo clippy --all-targets --all-features -- -D warnings`: exited `0` (0 warnings).
- `cargo test --all-features`: exited `0` — `59` active tests passed (`45` library + `14` across all `10` binaries in `src/bin/`: `bench_chr1` 1, `horizon` 1, `horizon_clean` 2, `horizon_eval` 1, `horizon_gfa` 1, `horizon_guide` 1, `horizon_lift` 3, `horizon_merge` 1, `horizon_resolve` 1, `horizon_wga` 2), `0` failed, `1` ignored (`hal::tests::real_hal_decode_from_env`, requires `HORIZON_TEST_HAL`).

## 2. Diff, Derived-Quantity, and Cross-File Checks (`manuscript/horizon_paper.md`)
- **Binary & test counts:** `10` binaries in `src/bin/` (`bench_chr1`, `horizon`, `horizon_clean`, `horizon_eval`, `horizon_gfa`, `horizon_guide`, `horizon_lift`, `horizon_merge`, `horizon_resolve`, `horizon_wga`) and `59` active tests (`45` lib + `14` bin + `1` env-gated) match `src/bin/` and `cargo test --all-features`.
- **Code constants & symbol names:** `max_gap_bp = 2400` (`src/mdl/mod.rs:128`), `KARLIN_MDL_COEFF = 0.55` (`src/mdl/mod.rs:26`), `band_w = 180` (`src/mdl/mod.rs:67`), `max_ext = 1500` (`src/mdl/mod.rs:68`), `score_junction_cp_bits` and `walk_support <= 1` / `prune_unsupported` (`src/bin/horizon_clean.rs:100,260`), `SpoolableGfaInput` (`src/gfa/mod.rs:20`), `MdlCostModel` (`src/egraph.rs:277`), `FlatSeedIndex` (`src/seeding/flat_index.rs:12`), `cluster_strand_lcbs` (`src/seeding/mod.rs:264`), `bocpd_gc_profile` (`src/bocpd.rs:58`), `mdl_xdrop_extend_with_params` (`src/mdl/dp.rs:448`), `interval_match_score` (`src/merge.rs:379`), and `rescue_transitive_maf` (`src/bin/horizon_guide.rs:184`) all match `src/`.
- **Benchmark numbers & derived quantities vs. `docs/BENCHMARK_RESULTS.md`:**
  - Table 1 (7 WGA rows) and Abstract highlights match `docs/BENCHMARK_RESULTS.md` lines 10–17, 36–37, 44–48; `cactus241way chr22 1MB` F1 `61.3%` matches `2 * 62.3 * 60.3 / (62.3 + 60.3) = 61.28%`.
  - `GRCh38.chr22` × `CHM13.chr22:15729460-51324926`: `34,016,799 / 35,595,466 = 95.56%` (`95.6%`).
  - `chr6` 20 Mbp q-arm RSS/Mbp: `2,409.7 MB / 20 Mbp` (`0.120 GB/Mbp` decimal) down from `5,863.3 MB` (`0.293 GB/Mbp`).
  - Full `chr6` (`53130b8` vs `1e47fe0`): `1,062.45 / 451.35 = 2.35×` total speedup; `588.67 / 47.87 = 12.30×` `lz1` speedup; `184.65 / 91.04 = 2.03×` `seed_end` speedup; `(778,471 - 693,242) / 693,242 = +12.29%` (`+12.3%`) DP jobs; `172,126,628 - 61,401,689 = 110,724,939` (`110.7 Mbp`) q-arm.
  - `horizon_gfa`, `horizon_lift`, `horizon_clean --gfa -`, `horizon_resolve --collapse-dups`, and `horizon_wga --graph` metrics match `docs/BENCHMARK_RESULTS.md` lines 55–68.
- **References:** All `38` entries in `# References` are cited in the manuscript body.

## 3. Round Harness Checks
- `scripts/round corrections 4`: `6 correction(s); 0 not applied or malformed` (`APPLIED` for `COR-1` through `COR-6` in `manuscript/horizon_paper.md`).
- `scripts/round drift 4`: exited `1` reporting `DRIFTED manuscript/horizon_paper.md: +56 / -65 line(s) since the snapshot this round attacked` (expected because `COR-1`..`COR-6` were applied after the round snapshot; `drift_reviewed` is already recorded in `round.json`).
- `scripts/round check 4 --full`: exited `0` (`manuscript/horizon_paper.md` drift is `REVIEWED` against 56 quoted lines).

## 4. What Was Not Reached
- External benchmark datasets (UCSC/Ensembl reference MAFs and the `12 GB` HPRC v1.0 GFA) were not re-executed during this cold audit; benchmark figures were verified for consistency against `docs/BENCHMARK_RESULTS.md`, commit/scope annotations, and arithmetic derivations.
- `attack-*.md`, `adjudication-*.md`, and `completeness-critic.md` were not read, per `AUDIT-BRIEF.md` cold-reader protocol.
