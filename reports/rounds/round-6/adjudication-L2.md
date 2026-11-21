# Round 6 ADJUDICATION of L2 (Sampling Lens)

**Attack evaluated:** `reports/rounds/round-6/attack-L2-sampling.md` (sealed sha256 prefix `93edb9286e96` verified; `HEAD` is `42700e1`).
**Prior-round & method-card check:** Checked `scripts/round history` for `"sampling"`, `"horizon_filter"`, `"horizon_locus"`, and `"horizon_resolve"` (`--kind adjudication`). No `refs/methods/` cards exist in this repository (`docs/round_rubric.md:36`), and no prior round adjudicated `horizon_filter`, `horizon_locus`, `horizon_resolve`'s `estimate_window_divergence`, or the `9dfec4a`/`42700e1` test-count and VGP/chr6 sampling claims.

---

## F1 -- C1: Test count is stale (`120` vs `118`) and P0–P4 binaries rest on 1–2 synthetic unit tests each

- **Verdict:** **UPHELD** (attacker's `NARROWED` ruling on C1 stands; diagnosis and repair both hold).
- **What survives:** `C1` (`docs/PANGENOME_SUITE.md:3`) is narrowed on its test-count figure (`120` passing unit tests on `42700e1`, not `118`) and on the fact that the new P0–P4 binaries each have only 1–2 synthetic unit tests.
- **Independent verification:**
  - Ran `cargo test --all-features` on `42700e1`: **120** passing unit tests (`86` lib + `34` bin: `1+2+6+2+1+2+4+11+1+1+1+2`) plus `1` ignored (`hal::tests::real_hal_decode_from_env`), vs `118` (`86` lib + `32` bin) pinned to `(2026-11-20, 9dfec4a)` at `docs/PANGENOME_SUITE.md:3` (and similarly at `manuscript/horizon_paper.md:84`).
  - Inspected `src/bin/horizon_filter.rs:484-517` (`1` synthetic 80-bp test, `0` tests for `--net`, `--max-occ`, `--keep-all`, or `.paf.zst`/`stdin` streaming), `src/bin/horizon_locus.rs:370-431` (`1` 160-bp LCG synthetic test), `src/bin/horizon_resolve.rs:900-944` (`1` synthetic test), and `src/bin/horizon_wga.rs:752-786` (`2` unit tests).
- **Repair:** **UPHELD — the fix is free.** Updating `docs/PANGENOME_SUITE.md:3` (and the matching `118` count at `manuscript/horizon_paper.md:84`) to `120` passing unit tests (`86` lib + `34` bin, `1` ignored) on `42700e1` and noting the single-test synthetic coverage of `horizon_filter`, `horizon_locus`, and `horizon_resolve` requires only text edits with no re-measurement.

---

## F2 -- C2: Unweighted record-count "truth hit rate" (`1.6%` -> `81.2%`) distorts repeat filtering, and `1,004` records does not reproduce on HEAD

- **Verdict:** **UPHELD** (attacker's `REFUTED` ruling on C2 stands; diagnosis and repair both hold, while noting `0.47s`–`0.48s` wall time and `~99.7%` 1-to-1 coverage reproduce).
- **What survives:** `C2` (`docs/BENCHMARK_RESULTS.md:85`) is refuted on its `1,004` records (`110.22 Mbp`) output count and its unweighted record-count `1.6%` -> `81.2%` 2D HPRC truth hit rate, which misrepresents `55.0%` p-arm/repeat fragments in `/tmp/chr6_fastga.paf`.
- **Independent verification:**
  - Ran `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --net --out /tmp/adj_l2_net.paf` on `42700e1`: `[horizon_filter] total=56857 kept=565 filtered_mdl=1 filtered_len=0 elapsed=0.48s` (`565` records, `109,995,785` bp), NOT `1,004` records (`110.22 Mbp`). Moreover, `git diff 9dfec4a 42700e1 -- src/bin/horizon_filter.rs src/netting.rs` is empty, confirming the `1,004` figure was already stale at `9dfec4a`.
  - Checked `/tmp/chr6_fastga.paf` and `/tmp/adj_l2_net.paf` against `/tmp/hprc_chr6_truth_segments.tsv` (`mult == "1"`): raw FastGA has `56,857` records (`292,357,161` bp), of which `31,298` (`55.0%`) lie in `GRCh38.chr6` p-arm (`<= 60 Mbp`) with no counterpart in the `CHM13.chr6` q-arm-only target (`61401689–172126628`). Unweighted record hit rate jumps from `~1.8%` to `~84.8%–85.0%`, whereas base-weighted aligned bp in HPRC-hitting records is already `~38.0%–38.1%` in raw FastGA (`62.59%` exact pair precision on shared HPRC bases per `/tmp/chr6_fastga_summary.json`) and rises to `~99.75%–99.85%` after `--net`.
- **Repair:** **UPHELD — the fix is free.** Replacing `1,004` records (`110.22 Mbp`, `81.2%` from `1.6%`) in `docs/BENCHMARK_RESULTS.md:85` (and the matching `1.6%` -> `81.2%` claim in `README.md:124` and `manuscript/horizon_paper.md:84`) with the reproduced `HEAD` output (`565` records, `110.00 Mbp`, both record-weighted and base-weighted hit rates, and `n = 1` q-arm-only target disclosure) requires no new runs.

---

## F3 -- C3: `horizon_locus` distance metrics suffer survivorship bias and benchmarks sample only `N=2` and `N=4` (25 kb)

- **Verdict:** **UPHELD** (attacker's `NARROWED` ruling on C3 stands; diagnosis and repair both hold).
- **What survives:** `C3` (`README.md:126`) is narrowed because `p_dist` and `jc69_dist` are conditioned exclusively on the `aligned_bp` subset (omitting `~33.8%–40.1%` of the avian locus), and the benchmarks sample only `N = 2` (`1` locus) and `N = 4` (`25 kb` windows).
- **Independent verification:**
  - Inspected `src/bin/horizon_locus.rs:135-145` (`p_dist = 1.0 - (matches / aligned_bp)`) and `:202-207` ($O(N^2)$ all-pairs loop).
  - Checked `/tmp/locus_avian_out.similarity.tsv`: `aligned_bp = 231,444` out of `386,668` bp (`galGal`, `40.1%` unaligned) and `349,383` bp (`taeGut`, `33.8%` unaligned).
  - Checked `/tmp/locus_human_out.similarity.tsv` and `/tmp/locus_human_25k.fa`: `N = 4` (`grch38_mhc`, `hla_dup`, `grch38_peri`, `chm13_peri`; `6` pairs of `25 kb` windows, with `hla_dup` `aligned_bp = 36,014 > 25,000` and `peri` `jaccard_cov = 0.0427–0.0801`).
- **Repair:** **UPHELD — the fix is free.** Qualifying `p_dist` and `jc69_dist` in `README.md:126` as conditioned on `aligned_bp` and noting the `N = 2` / `N = 4` (`25 kb`) benchmark scope is a text-only change.

---

## F4 -- C4: `CHM13.chr6` is q-arm only (`110.7 Mbp`), runs are `n = 1`, and `GRCh38.chr6` coverage exceeds target length

- **Verdict:** **UPHELD** (attacker's `NARROWED` ruling on C4 stands, noting that `manuscript/horizon_paper.md:78` already discloses `q-arm 61401689–172126628`, `n = 1 per run`, and `~2.58×` multiplicity in adjacent text).
- **What survives:** `C4` (`manuscript/horizon_paper.md:78`) is narrowed because `CHM13.chr6` omits `61.4 Mbp` of p-arm/centromeric sequence, `GRCh38.chr6`'s `124.2 Mbp` coverage exceeds the `110.7 Mbp` target via many-to-many paralog/repeat multiplicity, and whole-chromosome runs are `n = 1` human-only samples.
- **Independent verification:**
  - Verified `scripts/benchmarks/compare_chr6_horizon_vs_hprc.py:19-21` (`X_LEN = 170_805_979`, `Y_LEN = 110_724_939`, `Y_OFFSET = 61_401_689`).
  - Verified `manuscript/horizon_paper.md:78`: the excerpted sentence states `110.7 Mbp q-arm 61401689–172126628` and the surrounding paragraph notes `n = 1 per run` and `~2.58× overall multiplicity`, but does not explicitly flag the omitted `61.4 Mbp` p-arm/centromere or why `GRCh38.chr6` coverage (`124.2 Mbp`) exceeds target length (`110.7 Mbp`).
- **Repair:** **UPHELD — the fix is free.** Adding the brief clarification in `manuscript/horizon_paper.md:78` is a text-only edit.

---

## F5 -- C5: `horizon_resolve` "per-window divergence" samples topology step counts unweighted by base-pair length

- **Verdict:** **UPHELD** (attacker's `NARROWED` ruling on C5 stands; diagnosis and repair both hold).
- **What survives:** `C5` (`ARCHITECTURE.md:28`) is narrowed because `estimate_window_divergence` is a topology-based step-count heuristic bounded to `[0.02, 0.15]` unweighted by segment bp, not sequence-based divergence estimation.
- **Independent verification:**
  - Inspected `src/bin/horizon_resolve.rs:384, 402, 456-477`: `estimate_window_divergence` computes `(0.02 + 0.08 * non_anchor_frac + 0.05 * inv_frac).clamp(0.01, 0.25)` where each step has equal weight `1/n` regardless of `seg_lens[seg]` (so the output is strictly in `[0.02, 0.15]` and the clamp is a no-op). `resolve_window:384` skips `steps.len() <= 1`, and `:402` clamps segment lengths to `.max(100)` bp in `MdlCostModel`. (Minor precision note: a two-anchor window with 1 of 2 steps inverted receives `0.045`, not only `0.02` or `0.07`.)
- **Repair:** **UPHELD — the fix is free.** Updating `ARCHITECTURE.md:28` to describe `estimate_window_divergence` as a topology-based per-window divergence heuristic (`[0.02, 0.15]`, unweighted by segment bp) is a text-only edit.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `cargo test --all-features` runs `120` passing tests (`86` lib + `34` bin, `1` ignored) on `42700e1`, not `118`; P0–P4 binaries have 1–2 synthetic unit tests each; doc fix is free |
| F2 | UPHELD | `horizon_filter --net` outputs `565` records (`110.00 Mbp`) on `42700e1`, not `1,004` (`110.22 Mbp`), and unweighted `1.6%` -> `81.2%` distorts `55.0%` p-arm/repeat hits (`38.0%–38.1%` -> `99.8%` base-weighted); doc fix is free |
| F3 | UPHELD | `p_dist`/`jc69_dist` are conditioned on `aligned_bp` only (omitting `33.8%–40.1%` of the avian locus), and benchmarks sample only `N=2` and `N=4` (`25 kb`); doc fix is free |
| F4 | UPHELD | `CHM13.chr6` target omits `61.4 Mbp` p-arm/centromere, `GRCh38.chr6` `124.2 Mbp` coverage reflects `~2.58×` multiplicity on `110.7 Mbp` target, and runs are `n=1` human-only; doc fix is free |
| F5 | UPHELD | `estimate_window_divergence` computes an unweighted topology-step heuristic bounded to `[0.02, 0.15]`, not sequence divergence; doc fix is free |
