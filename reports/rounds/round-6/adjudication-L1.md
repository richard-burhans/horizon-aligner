# Round 6 Adjudication — L1 (Correctness)

- **Target:** `horizon-aligner` (`HEAD` = `42700e1`)
- **Attack file verified:** `reports/rounds/round-6/attack-L1-correctness.md` (`sha256`: `f3193ca7ef4bb98086722bce208ddd4279e6b83e7cfa3c7ee9a9a4c094490c4a`, matching sealed prefix `f3193ca7ef4b`)
- **Method cards:** `refs/methods/` does not exist in this repository; rulings rest on live code inspection, `cargo test --no-fail-fast`, and re-execution of compiled release binaries.

---

## F1 — C1 (`docs/PANGENOME_SUITE.md:3`): Unit test count is stale (`120` active, not `118`) and `horizon_filter` does not run BOCPD

- **Ruling on diagnosis:** **UPHELD.**
  1. Re-running `cargo test --no-fail-fast` on `42700e1` passes **120** active unit tests (`86` in `src/lib.rs` + `34` across the 12 binaries in `src/bin/`: `bench_chr1` 1, `horizon` 2, `horizon_clean` 6, `horizon_eval` 2, `horizon_filter` 1, `horizon_gfa` 2, `horizon_guide` 4, `horizon_lift` 11, `horizon_locus` 1, `horizon_merge` 1, `horizon_resolve` 1, `horizon_wga` 2; plus `1` ignored `.hal` test in `src/lib.rs`). Both `docs/PANGENOME_SUITE.md:3` and `manuscript/horizon_paper.md:84` still report `118` (`86` lib + `32` bin).
  2. `src/bin/horizon_filter.rs:4-6` imports only `GAP_EXTEND_BITS`, `GAP_OPEN_BITS`, `IDENTITY_PRIOR`, `base_probs`, and `joint_pair_bits` from `horizon_aligner::bocpd`. It never calls `bocpd_gc_profile`, `bocpd_gc_profile_into`, or `BocpdScratch`; it computes a contig-global GC fraction (`IndexedSeq::from_dna`, lines `38-58`) when FASTAs are supplied or fixed `0.5` GC in `score_from_counts_only` (lines `123-130`).
  3. `src/io_util.rs:92-114, 136-156` delegates `.zst`/`.zstd` to an external `zstd` subprocess, and `which zstd` confirms `zstd` is absent on this host, so `zstd_roundtrip_and_early_drop_do_not_deadlock` (`src/io_util.rs:198-200`) returns early without exercising `.zst` I/O.
- **What survives:** `C1` overstates `horizon_filter` as running BOCPD and cites a stale active unit-test count (`118` instead of `120`).
- **Ruling on correction:** **ENDORSED (text update); the fix is free.** Updating `docs/PANGENOME_SUITE.md:3` and `manuscript/horizon_paper.md:84` to `120` passing unit tests (`86` lib + `34` bin, `1` ignored) and describing `horizon_filter` as a streaming PAF SDUST + global-GC Pair-HMM MDL filter & netter requires no re-measurement. (Wiring `bocpd_gc_profile_into` into `horizon_filter.rs` instead would be a non-free code change requiring re-testing and re-benchmarking.)

---

## F2 — C2 (`docs/BENCHMARK_RESULTS.md:85`): `horizon_filter --net` produces `565` records (`110.00 Mbp`), not `1,004` records (`110.22 Mbp`)

- **Ruling on diagnosis:** **UPHELD.** Re-running `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --net --out /tmp/adj_c2_net.paf` on `42700e1` outputs `[horizon_filter] total=56857 kept=565 filtered_mdl=1 filtered_len=0 elapsed=0.45s`, with `565` records and `109,995,785` `aln_len` bp (`110.00 Mbp`), not `1,004` records (`110.22 Mbp`). Because `zstd` is not installed on this host, the `0.49s` `.paf.zst` timing in `C2` could not be independently re-run here.
- **What survives:** `C2`'s `1,004` records (`110.22 Mbp`) figure is stale against `HEAD` (`565` records, `110.00 Mbp` `aln_len`).
- **Ruling on correction:** **ENDORSED WITH COST (not free).** Replacing `1,004` records (`110.22 Mbp`) with `565` records (`110.00 Mbp` `aln_len`) in `docs/BENCHMARK_RESULTS.md:85` is required, and the `81.2%` 2D HPRC truth hit rate and `99.7%` 1-to-1 coverage must be re-measured on the `565`-record output rather than carried over from the earlier `1,004`-record run.

---

## F3 — C3 (`README.md:126`): Human 4-way `horizon_locus` block counts (`29 MAF blocks / 20 multi-species`) do not match `HEAD` (`102` / `92` or `320` / `110`), and `jaccard_cov` is not Jaccard similarity

- **Ruling on diagnosis:** **UPHELD.**
  1. Re-running `./target/release/horizon_locus --fasta /tmp/locus_avian.fa --out-prefix /tmp/adj_c3_avian --threads 4` reproduces `231,444` aligned bp, `p_dist = 0.192552`, `jc69_dist = 0.222528`, and `1,364` MAF blocks (`0.89s`, `223.4 MB` RSS).
  2. Re-running `./target/release/horizon_locus --fasta /tmp/locus_human_25k.fa --out-prefix /tmp/adj_c3_human --threads 4` selects `focal=chm13_peri` (`max_by_key` at `src/bin/horizon_locus.rs:190-195`) and emits `102` MAF blocks / `92` multi-blocks (`0.42s`, `170.8 MB` RSS); passing `--ref grch38_mhc` emits `320` MAF blocks / `110` multi-blocks (`0.30s`, `182.4 MB` RSS). Neither matches `29 MAF blocks / 20 multi-species`, and in `/tmp/locus_human_out.maf` all `29` blocks have `>= 2` `s`-lines (`min = 2`), not `20`.
  3. `src/bin/horizon_locus.rs:128-131` computes `jaccard_cov` as `(cov_x + cov_y) as f64 / (seq_x.len() + seq_y.len()) as f64` (average covered fraction), not `|A ∩ B| / |A ∪ B|`.
- **What survives:** `C3`'s avian 2-way metrics reproduce, but its human 4-way block counts (`29` / `20`) do not match `HEAD` (`102`/`92` default or `320`/`110` with `--ref grch38_mhc`) and `jaccard_cov` is misnamed.
- **Ruling on correction:** **ENDORSED; the text fix is free** (using the already-reproduced `HEAD` block counts across `README.md:126`, `docs/BENCHMARK_RESULTS.md:91`, and `manuscript/horizon_paper.md:84`, and specifying `--ref`). Note that `docs/BENCHMARK_RESULTS.md:91` also reports `grch38_mhc` vs `hla_dup` `jaccard_cov = 0.999400` and `p_dist = 0.124507`, which should be checked when updating that row.

---

## F4 — C4 (`manuscript/horizon_paper.md:78`): `8cd5db0` full `chr6` metrics and speedup ratios are internally consistent

- **Ruling on diagnosis:** **UPHELD.** Independent arithmetic check confirms `172126628 - 61401689 = 110,724,939` bp (`110.7 Mbp`), `124,160,423 / 170,805,979 = 72.69%` (`72.7%`), `109,464,447 / 110,724,939 = 98.86%` (`98.9%`), `592.68 / 267.86 = 2.21×`, and `1,062.45 / 267.86 = 3.97×`. `C4` is explicitly pinned to `8cd5db0`, and the following sentence in `manuscript/horizon_paper.md:78` reports the `2896229` multi-chain run matching `manuscript/data/chr6_hprc_comparison.json`.
- **What survives:** `C4` stands as written; no repair is needed.

---

## F5 — C5 (`ARCHITECTURE.md:28`): `horizon_resolve` divergence estimate is topological-only and `W`-line `end` bounds ignore `L`-line overlaps

- **Ruling on diagnosis:** **UPHELD.**
  1. `src/bin/horizon_resolve.rs:456-477` (`estimate_window_divergence`) computes `(0.02 + 0.08 * non_anchor_frac + 0.05 * inv_frac).clamp(0.01, 0.25)` solely from non-anchor and reverse-step fractions without inspecting segment DNA sequences.
  2. `src/bin/horizon_resolve.rs:722-727` (`emit_resolved_walk`) sets `end_pos = w.start.saturating_add(resolved_bp)` by summing `seg_lens` without subtracting `L`-line CIGAR overlaps (unlike `src/bin/horizon_clean.rs:444`), and `pass2_rewrite` (`src/bin/horizon_resolve.rs:841-848`) hardcodes `0M` on synthesized `L` links.
- **What survives:** `C5` omits that `horizon_resolve`'s per-window divergence is a topological heuristic and that its `W`-line `end` bounds and synthesized `L` links assume `0M` overlaps.
- **Ruling on correction:** **ENDORSED; the fix is free.** Qualifying `ARCHITECTURE.md:28` to state both assumptions is a text-only clarification requiring no re-measurement.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `cargo test` passes 120 active tests (86 lib + 34 bin), not 118, and `horizon_filter.rs` uses global/fixed GC without calling BOCPD; text fix is free |
| F2 | UPHELD | `horizon_filter --net` on `/tmp/chr6_fastga.paf` outputs 565 records (110.00 Mbp `aln_len`), not 1,004 (110.22 Mbp); truth-hit/cov metrics require re-measurement |
| F3 | UPHELD | Human 4-way `horizon_locus` emits 102/92 (default) or 320/110 (`--ref grch38_mhc`), not 29/20, and `jaccard_cov` is average covered fraction; text fix is free |
| F4 | UPHELD | Attacker's `SURVIVED` assessment of `C4` is verified: all spans, percentages, speedups, and the `8cd5db0` pin hold |
| F5 | UPHELD | `estimate_window_divergence` is topological-only and `emit_resolved_walk` / `pass2_rewrite` assume `0M` overlaps; text fix is free |
