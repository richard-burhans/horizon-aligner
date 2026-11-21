# Round 6 Attack — L1 (Correctness)

- **Target:** `horizon-aligner` (`42700e1` / `9dfec4a` in `src/`)
- **Method cards:** This project keeps **no method cards** yet (`docs/round_rubric.md:36`).
- **Literature:** No literature claims were asserted from memory; this correctness attack checked the live Rust codebase (`src/`), `cargo test --no-fail-fast`, compiled release binaries, and local benchmark artifacts.

---

## F1 — C1 (`docs/PANGENOME_SUITE.md:3`): Unit test count is stale (`120` active, not `118`) and `horizon_filter` does not run BOCPD

- **Claim:** `C1` (`docs/PANGENOME_SUITE.md:3`)
- **Verdict:** **NARROWED**
- **Prior-round context:** Round 5 `F1` (`UPHELD`, `adjudication-L4.md:14-16`) previously caught a stale active unit-test count (`98` vs `97` + `1` ignored). That count was updated, but `C1` is now stale again after subsequent binary tests were added.
- **Concrete flaws:**
  1. **Wrong unit-test count:** Running `cargo test --no-fail-fast` on `HEAD` (`42700e1`) executes **120** passing unit tests (`86` lib in `src/lib.rs` + `34` bin across the 12 binaries in `src/bin/`, plus `1` ignored `.hal` test in `src/lib.rs`), not `118` (`86` lib + `32` bin). Specifically, the binary test counts are: `bench_chr1` (`1`), `horizon` (`2`), `horizon_clean` (`6`), `horizon_eval` (`2`), `horizon_filter` (`1`), `horizon_gfa` (`2`), `horizon_guide` (`4`), `horizon_lift` (`11`), `horizon_locus` (`1`), `horizon_merge` (`1`), `horizon_resolve` (`1`), `horizon_wga` (`2`), summing to `34`. (The same stale `118` / `32` count appears in `manuscript/horizon_paper.md:84`.)
  2. **`horizon_filter` does not use BOCPD:** `C1` describes `horizon_filter` as `(streaming PAF SDUST+BOCPD+MDL filter & netter)`. In `src/bin/horizon_filter.rs:4-6`, it imports only `GAP_EXTEND_BITS`, `GAP_OPEN_BITS`, `IDENTITY_PRIOR`, `base_probs`, and `joint_pair_bits` from `horizon_aligner::bocpd`. It never calls `bocpd_gc_profile`, `bocpd_gc_profile_into`, or `BocpdScratch`; instead it computes a single contig-global GC fraction (`IndexedSeq::from_dna`, lines `38-58`) when FASTAs are provided, or fixed `0.5` GC in sequence-free mode (`score_from_counts_only`, lines `123-130`).
  3. **`.zst` support depends on external `zstd`:** `src/io_util.rs:92-114, 136-156` delegates `.zst`/`.zstd` compression and decompression to an external `zstd` subprocess rather than pure-Rust decoding. On this host (`zstd: command not found`), `zstd_roundtrip_and_early_drop_do_not_deadlock` (`src/io_util.rs:198`) silently returns early without testing `.zst` I/O.
- **Correction:** Update `docs/PANGENOME_SUITE.md:3` (and `manuscript/horizon_paper.md:84`) to report **120** passing unit tests (`86` lib + `34` bin, `1` ignored), and describe `horizon_filter` accurately as a streaming PAF SDUST + global-GC Pair-HMM MDL filter & netter (or wire `bocpd_gc_profile_into` into `horizon_filter.rs`).

---

## F2 — C2 (`docs/BENCHMARK_RESULTS.md:85`): `horizon_filter --net` produces `565` records (`110.00 Mbp`), not `1,004` records (`110.22 Mbp`)

- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:85`)
- **Verdict:** **REFUTED**
- **Prior-round context:** `scripts/round history "horizon_filter" --kind adjudication` returned zero prior mentions.
- **Concrete flaw:** Re-running the release binary built from `HEAD` (`42700e1` / `9dfec4a`) on the exact `93 MB` `56,857`-record `GRCh38.chr6` × `CHM13.chr6` `FastGA` PAF (`/tmp/chr6_fastga.paf`):
  ```text
  target/release/horizon_filter --paf /tmp/chr6_fastga.paf --net --out /tmp/test_c2_net2.paf
  [horizon_filter] total=56857 kept=565 filtered_mdl=1 filtered_len=0 elapsed=0.44s
  ```
  The output contains **`565` records** (with `109,995,785` `aln_len` bp = **`110.00 Mbp`**, `109,932,124` query bp, `109,934,990` target bp), **not** `1,004` records (`110.22 Mbp`) as claimed in `docs/BENCHMARK_RESULTS.md:85`. Even the earlier `/tmp/chr6_fastga_mdl_net.paf` on disk contains `553` records (`109,990,091` `aln_len` bp), not `1,004`. (Note also that `zstd` is absent on this host, so the `.paf.zst` timing could not be re-run here.)
- **Correction:** Replace `1,004` records (`110.22 Mbp`) in `docs/BENCHMARK_RESULTS.md:85` with **`565` records (`110.00 Mbp` `aln_len` / `109.93 Mbp` span)** and re-verify the `81.2%` 2D HPRC truth hit rate and `99.7%` 1-to-1 coverage on the `565`-record output.

---

## F3 — C3 (`README.md:126`): Human 4-way `horizon_locus` block counts (`29 MAF blocks / 20 multi-species`) do not match `HEAD` (`102` / `92` or `320` / `110`), and `jaccard_cov` is not Jaccard similarity

- **Claim:** `C3` (`README.md:126`)
- **Verdict:** **NARROWED**
- **Prior-round context:** `scripts/round history "horizon_locus" --kind adjudication` returned zero prior mentions.
- **Concrete flaws:**
  1. **Avian 2-way (`galGal` × `taeGut`) reproduces, with minor metric drift:** Running `target/release/horizon_locus --fasta /tmp/locus_avian.fa --out-prefix /tmp/test_c3_avian --threads 4` on `HEAD` reproduces `231,444` aligned bp, `p_dist = 0.192552`, `jc69_dist = 0.222528`, and `1,364` MAF blocks in `0.80s` (`226.9 MB` RSS).
  2. **Human 4-way 25 kb MAF/multi-species counts do not match `HEAD`:** Running `target/release/horizon_locus --fasta /tmp/locus_human_25k.fa --out-prefix /tmp/test_c3_human --threads 4` on `HEAD` selects `focal=chm13_peri` (because `max_by_key` in `src/bin/horizon_locus.rs:190-195` picks the last equal-length 25 kb contig) and emits **`102` MAF blocks / `92` multi-blocks** (`0.42s`, `174.7 MB` RSS); passing `--ref grch38_mhc` emits **`320` MAF blocks / `110` multi-blocks** (`0.32s`, `170.9 MB` RSS). Neither matches `29 MAF blocks / 20 multi-species`. Moreover, even in the older `/tmp/locus_human_out.maf` (`29` blocks), all `29` blocks have `>= 2` `s`-lines (`[2, 3, 3, 4, ..., 2]`), not `20`.
  3. **`jaccard_cov` formula is not Jaccard similarity:** In `src/bin/horizon_locus.rs:128-131`, `jaccard_cov` is computed as `(cov_x + cov_y) as f64 / (seq_x.len() + seq_y.len()) as f64` (average covered fraction), not `|A ∩ B| / |A ∪ B|`.
- **Correction:** Update `README.md:126` (and `docs/BENCHMARK_RESULTS.md:91`, `manuscript/horizon_paper.md:84`) with the `HEAD` human 4-way `horizon_locus` block counts (`102` MAF / `92` multi-species for default focal `chm13_peri`, or `320` / `110` for `--ref grch38_mhc`) and specify which `--ref` was used.

---

## F4 — C4 (`manuscript/horizon_paper.md:78`): `8cd5db0` full `chr6` metrics and speedup ratios are internally consistent

- **Claim:** `C4` (`manuscript/horizon_paper.md:78`)
- **Verdict:** **SURVIVED**
- **Evidence checked:**
  - Coordinate span: `CHM13.chr6` q-arm `61401689–172126628` is `110,724,939` bp (`110.7 Mbp`), and `GRCh38.chr6` is `170,805,979` bp (`170.8 Mbp`).
  - Coverage percentages: `124,160,423 / 170,805,979 = 72.69%` (`72.7%`) and `109,464,447 / 110,724,939 = 98.86%` (`98.9%`).
  - Speedup ratios: `592.68 / 267.86 = 2.21×` and `1,062.45 / 267.86 = 3.97×`.
  - Provenance pin: `C4` is explicitly pinned to `8cd5db0`, and the immediately following sentence in `manuscript/horizon_paper.md:78` reports the updated `2896229` multi-chain numbers (`304.53 s`, `282,163,046` pairs), matching `manuscript/data/chr6_hprc_comparison.json`.

---

## F5 — C5 (`ARCHITECTURE.md:28`): `horizon_resolve` divergence estimate is topological-only and `W`-line `end` bounds ignore `L`-line overlaps

- **Claim:** `C5` (`ARCHITECTURE.md:28`)
- **Verdict:** **NARROWED**
- **Prior-round context:** Round 2 `R2-RESOLVE-LIFT` and Round 3 `R2-RESOLVE-LIFT` required updating `W`-line coordinates and `L` links in `horizon_resolve`; `src/bin/horizon_resolve.rs` now implements both, but with two remaining correctness defects.
- **Concrete flaws:**
  1. **`estimate_window_divergence` does not inspect DNA sequence:** `src/bin/horizon_resolve.rs:456-477` computes `(0.02 + 0.08 * non_anchor_frac + 0.05 * inv_frac).clamp(0.01, 0.25)` purely from the fraction of non-anchor and reverse-oriented steps in the window, without reading segment sequences.
  2. **`W`-line `end` bounds and synthesized `L` links ignore `L`-line overlaps:** In `emit_resolved_walk` (`src/bin/horizon_resolve.rs:710-727`), when `resolved != steps`, `end_pos` is set to `w.start.saturating_add(resolved_bp)` where `resolved_bp` is the raw sum of `seg_lens` — unlike `horizon_clean.rs:444`, it does not subtract `L`-line CIGAR overlaps, overestimating `end_pos` on any graph with non-`0M` overlaps. Similarly, `pass2_rewrite` (`src/bin/horizon_resolve.rs:841-848`) hardcodes `0M` overlap on synthesized `L` links, and `stream_gfa_records` silently drops any `C`/`J` GFA records.
- **Correction:** Qualify `ARCHITECTURE.md:28` to state that `horizon_resolve` uses a topological step-orientation/non-anchor heuristic for per-window divergence and assumes `0M` overlaps when updating `W`-line `end` bounds and synthesizing `L` links.

---

## Verdict Table

| Claim | Verdict | Summary |
| :--- | :--- | :--- |
| **C1** | **NARROWED** | `cargo test` passes `120` active tests (`86` lib + `34` bin, `1` ignored), not `118` (`32` bin); `horizon_filter` does not call BOCPD (`bocpd_gc_profile*`), using global/fixed GC. |
| **C2** | **REFUTED** | `horizon_filter --net` on `/tmp/chr6_fastga.paf` outputs `565` records (`110.00 Mbp` `aln_len`) in `0.44–0.47s`, not `1,004` records (`110.22 Mbp`). |
| **C3** | **NARROWED** | Avian 2-way reproduces (`231,444` bp, `p_dist = 0.192552`, `jc69_dist = 0.222528`), but human 4-way on `HEAD` produces `102`/`92` (default) or `320`/`110` (`--ref grch38_mhc`), not `29`/`20`; `jaccard_cov` is not Jaccard similarity. |
| **C4** | **SURVIVED** | Arithmetic (`2.21×`, `3.97×`, `72.7%`, `98.9%`), `110.7 Mbp` q-arm span, and `8cd5db0` pin all hold. |
| **C5** | **NARROWED** | `estimate_window_divergence` is a topological heuristic (no sequence comparison), and `W`-line `end` bounds / synthesized `L` links assume `0M` overlaps. |
