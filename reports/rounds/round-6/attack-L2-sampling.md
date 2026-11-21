# Round 6 Attack -- L2 (Sampling)

**Scope & Prior-Round / Method-Card Check:**
- As documented in `docs/round_rubric.md:36`, this project keeps **no method cards** yet; sampling claims are evaluated directly against the codebase, benchmark datasets/scripts, and live CLI outputs.
- Checked `scripts/round history "sampling" --kind adjudication`, `scripts/round history "horizon_filter" --kind adjudication` (0 prior mentions), `scripts/round history "horizon_locus" --kind adjudication` (0 prior mentions), and `scripts/round history "horizon_resolve" --kind adjudication`. Prior rounds (e.g., Round 2 L2 F1, Round 5 L2 F6) addressed `estimate_hmm_params` prefix sampling and `# EG:Z:` format claims; none adjudicated `horizon_filter`, `horizon_locus`, `horizon_resolve`'s `estimate_window_divergence`, or the `9dfec4a`/`42700e1` test-count and VGP/chr6 sampling claims below.

---

## F1 -- C1: Test count is stale (`120` vs `118`) and P0–P4 binaries rest on 1–2 synthetic unit tests each

- **Claim:** `C1` (`docs/PANGENOME_SUITE.md:3`)
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. `cargo test --all-features` on `HEAD` (`42700e1`) runs **120** passing unit tests (`86` lib + `34` bin) plus `1` ignored test, not `118` (`86` lib + `32` bin).
  2. Under the sampling lens, the newly claimed P0–P4 suite tools are barely sampled by the test suite: `horizon_filter` (`src/bin/horizon_filter.rs:484-517`) has **1** unit test on two 80-bp synthetic strings (`complex` vs `AT` repeat) and **0** tests for `--net`, `--max-occ`, `--keep-all`, `--min-nonrepeat-bp`, `--min-len`, `score_from_counts_only`, or `.paf.zst`/`stdin` streaming; `horizon_locus` (`src/bin/horizon_locus.rs:370-431`) has **1** unit test on a 160-bp LCG-generated synthetic sequence; `horizon_resolve` (`src/bin/horizon_resolve.rs`) has **1** unit test; `horizon_wga --graph` (`src/bin/horizon_wga.rs`) has **2** unit tests; and the **1 ignored** test (`hal::tests::real_hal_decode_from_env`) is the only test sampling a real external `.hal` file.
- **Evidence:** `cargo test --all-features` output on `42700e1` (`86` lib + `1+2+6+2+1+2+4+11+1+1+1+2 = 34` bin = `120` passing, `1` ignored); `docs/PANGENOME_SUITE.md:3` still pins `(2026-11-20, 9dfec4a)`.
- **Correction:** Update `docs/PANGENOME_SUITE.md:3` to `120` passing unit tests (`86` lib + `34` bin, `1` ignored) on `42700e1`, and note that `horizon_filter`, `horizon_locus`, and `horizon_resolve` each have only a single synthetic unit test (`--net` and multi-record streaming are untested in `cargo test`).

---

## F2 -- C2: Unweighted record-count "truth hit rate" (`1.6%` -> `81.2%`) distorts repeat filtering, and `1,004` records does not reproduce on HEAD

- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:85`)
- **Verdict:** **REFUTED**
- **Concrete flaw:**
  1. **Unweighted record-count sampling bias:** Measuring "2D HPRC truth hit rate" per PAF record (`hit_records / total_records`) treats a 50-bp repeat fragment and a 1-Mbp syntenic block as equal `n = 1` units. In `/tmp/chr6_fastga.paf` (`56,857` records, `292,357,161` aln bp), `31,298` records (`55.0%`) lie in `GRCh38.chr6` p-arm (`<= 60 Mbp`), which has **no** counterpart in the `CHM13.chr6:61401689-172126628` q-arm-only target, and tens of thousands of short repeat hits inflate the denominator. Base-weighted, `38.12%` of raw FastGA aligned bases already belong to HPRC-hitting records (`62.59%` exact pair precision on shared HPRC bases per `/tmp/chr6_fastga_summary.json`).
  2. **Stale / unreproduced output count on HEAD:** Running `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --net` on `HEAD` (`42700e1`) produces **565** records (`109,995,785` bp; `84.96%` 2D record hit rate; `99.85%` of emitted bases in HPRC-hitting records), NOT `1,004` records (`110.22 Mbp`, `81.2%`), with `filtered_mdl=1` (i.e. greedy interval netting removed 56,291 records, while MDL filtering removed 1).
  3. **Single-chromosome q-arm sample:** The benchmark evaluates `n = 1` human chromosome pair (`GRCh38.chr6` × `CHM13.chr6` q-arm only).
- **Evidence:**
  - `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --net --out /tmp/test_l2_net.paf`: `[horizon_filter] total=56857 kept=565 filtered_mdl=1 filtered_len=0 elapsed=0.47s`.
  - Overlap check against `/tmp/hprc_chr6_truth_segments.tsv` (`mult == "1"`): `raw_fastga: n=56857, bp=292,357,161, p_arm(<=60M)=31298 (55.0%), hit_2d=1008 (1.77%), bp_in_hit_recs=38.12%`; `horizon_filter_net: n=565, bp=109,995,785, p_arm(<=60M)=76 (13.5%), hit_2d=480 (84.96%), bp_in_hit_recs=99.85%`.
  - `/tmp/chr6_fastga_summary.json`: `exact_coord_precision_on_shared = 0.625885312833943`, `exact_coord_pairs = 109245821` out of `total_pairs = 285646096` (`38.25%`).
- **Correction:** Replace `1,004` records (`110.22 Mbp`), `81.2%` 2D HPRC truth hit rate (from `1.6%`) in `docs/BENCHMARK_RESULTS.md:85` with the HEAD (`42700e1`) output (`565` records, `110.00 Mbp`, `85.0%` record hit rate from `1.8%` unweighted / `99.85%` base-weighted from `38.1%`), and disclose that the benchmark is `n = 1` chromosome pair with `CHM13.chr6` restricted to the q-arm.

---

## F3 -- C3: `horizon_locus` distance metrics suffer survivorship bias and benchmarks sample only `N=2` and `N=4` (25 kb)

- **Claim:** `C3` (`README.md:126`)
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. **Survivorship / ascertainment bias in `p_dist` and `jc69_dist`:** In `src/bin/horizon_locus.rs:135-145`, `p_dist = 1.0 - (matches / aligned_bp)` and `jc69_dist` are computed exclusively over `aligned_bp` (bases surviving MDL + Horizon Bubble gating). In `/tmp/locus_avian_out.similarity.tsv`, `aligned_bp = 231,444` out of `349,383` bp (`taeGut`) and `386,668` bp (`galGal`), omitting `~34%–40%` of the locus (unaligned/divergent/bubble sequence). Reporting `p_dist = 0.192552` and `jc69_dist = 0.222528` without stating they are conditioned on the aligned subset understates whole-locus divergence.
  2. **Non-representative small-$N$ sampling:** `horizon_locus` executes $O(N^2)$ all-pairs alignment (`src/bin/horizon_locus.rs:202-207`), yet C3 samples only `N = 2` (`1` pair, `n = 1` avian locus) and `N = 4` (`6` pairs of `25 kb` windows in `/tmp/locus_human_25k.fa`: `grch38_mhc`, synthetic `hla_dup` with `aligned_bp = 36,014 > 25,000`, and unrelated `grch38_peri`/`chm13_peri` with `jaccard_cov = 0.0427–0.0801`), leaving `impg`-style multi-haplotype pangenome scaling (e.g., 50–100 haplotypes) unsampled.
- **Evidence:** `src/bin/horizon_locus.rs:135-145, 202-207`; `/tmp/locus_avian_out.similarity.tsv`; `/tmp/locus_human_out.similarity.tsv`; `/tmp/locus_human_25k.fa`.
- **Correction:** Qualify `p_dist = 0.192552` and `jc69_dist = 0.222528` in `README.md:126` as conditioned on the `231,444` aligned bp (`~60–66%` of locus length), and note that `horizon_locus` is benchmarked only on `N = 2` (`1` locus) and `N = 4` (`25 kb` composite windows), not large-$N$ pangenome cohorts.

---

## F4 -- C4: `CHM13.chr6` is q-arm only (`110.7 Mbp`), runs are `n = 1`, and `GRCh38.chr6` coverage exceeds target length

- **Claim:** `C4` (`manuscript/horizon_paper.md:78`)
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. **Partial-chromosome target sampling:** C4 frames the benchmark as "full `GRCh38.chr6` (`170.8 Mbp`) × `CHM13.chr6` (`110.7 Mbp` q-arm `61401689–172126628`)". While `GRCh38.chr6` is full-length, `CHM13.chr6` is truncated to its q-arm, omitting `61.4 Mbp` of p-arm and centromeric sequence.
  2. **Many-to-many inflation beyond target length:** `GRCh38.chr6` covers `124,160,423` bp (`72.7%`) against a `110,724,939` bp target (`13.4 Mbp` more query bases than the entire target), with `282,193,564` aligned pairs (`~2.58×` multiplicity on `109.46 Mbp` of target).
  3. **`n = 1` chromosome and timing sample:** As `manuscript/horizon_paper.md:78` notes, timings are `n = 1` per run (a two-run range `258.88–267.86 s` on `8cd5db0` vs single runs on `bb43dd9` and `1e47fe0`), and whole-chromosome evaluation is restricted to `chr6` (plus `chr22` q-arm) with zero non-human whole-chromosome samples.
- **Evidence:** `manuscript/horizon_paper.md:78`; `scripts/benchmarks/compare_chr6_horizon_vs_hprc.py:19-21` (`X_LEN = 170_805_979`, `Y_LEN = 110_724_939`, `Y_OFFSET = 61_401_689`).
- **Correction:** Clarify in `manuscript/horizon_paper.md:78` that the `CHM13.chr6` target excludes the `61.4 Mbp` p-arm/centromere (so `GRCh38.chr6`'s `124.2 Mbp` coverage reflects many-to-many paralog/repeat mappings onto the `110.7 Mbp` q-arm) and that speedups rest on single-chromosome `n = 1` runs.

---

## F5 -- C5: `horizon_resolve` "per-window divergence" samples topology step counts unweighted by base-pair length

- **Claim:** `C5` (`ARCHITECTURE.md:28`)
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. `estimate_window_divergence` (`src/bin/horizon_resolve.rs:456-477`) does **not** sample sequence bases or alignments; it computes `0.02 + 0.08 * non_anchor_frac + 0.05 * inv_frac` (which is bounded to `[0.02, 0.15]` despite `.clamp(0.01, 0.25)`).
  2. Under the sampling lens, every GFA segment step in the window receives equal weight (`1/n`) regardless of segment length in bp: a 1-bp SNP node has the same weight as a 100-kb anchor. Windows with `steps.len() <= 1` are skipped (`resolve_window:384`), two-anchor windows always receive `0.02` (or `0.07` if inverted), and `resolve_window:402` clamps all segment lengths to `.max(100)` bp in `MdlCostModel`.
- **Evidence:** `src/bin/horizon_resolve.rs:384, 402, 456-477`.
- **Correction:** Update `ARCHITECTURE.md:28` to state that `horizon_resolve` uses a topology-based per-window divergence heuristic (fraction of non-anchor and inverted steps, bounded to `[0.02, 0.15]`, unweighted by segment bp) rather than sequence-based divergence estimation.

---

## Verdict Summary

| Claim ID | Verdict | Summary |
| :--- | :--- | :--- |
| C1 | NARROWED | `120` passing tests (`86` lib + `34` bin, `1` ignored), not `118`; P0–P4 binaries have only 1–2 synthetic unit tests each |
| C2 | REFUTED | Unweighted record-count hit rate (`1.6%` -> `81.2%`) distorts `55%` p-arm/repeat fragments; HEAD outputs `565` records, not `1,004` |
| C3 | NARROWED | `p_dist`/`jc69_dist` conditioned on aligned subset only (survivorship bias); only `N=2` and `N=4` (`25 kb`) sampled |
| C4 | NARROWED | `CHM13.chr6` target is q-arm only (`110.7 Mbp`), `GRCh38.chr6` `124.2 Mbp` coverage reflects `~2.58×` multiplicity, `n=1` runs |
| C5 | NARROWED | `estimate_window_divergence` samples unweighted GFA step topology (`[0.02, 0.15]`), not sequence divergence |
