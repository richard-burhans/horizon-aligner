# Adversarial Attack — Round 3 (Lens: Inference)

**Target:** `horizon-aligner` (`488b020`) — Claims `C1`–`C6` from `docs/STRATEGIC_PLAN.md` and `docs/BENCHMARK_RESULTS.md`.
**Prior Round Disclosure:** `scripts/round history` was queried before filing (`FlatSeedIndex`: 0 prior hits; `stream_gfa`: 27 prior hits across Rounds 1–2; `cactus241way`: 10 prior hits across Rounds 1–2). Round 1 and Round 2 previously adjudicated `stream_gfa` (`NARROWED`/`UPHELD` on dead-code/non-zero-copy/stdin gaps) and `cactus241way` (`UPHELD` on multi-commit staleness and MAF-extracted sampling). Below, only *new* or *post-`488b020`* inference defects are filed, with prior adjudications noted where relevant.

---

## F1 -- C1: `stream_gfa` tests do not license GFA v1.0/v1.1 interoperability or `stdin` `.gfa.gz` streaming

* **Claim:** `C1` (`docs/STRATEGIC_PLAN.md:13`)
* **Verdict:** **NARROWED**
* **Concrete flaw:** While `stream_gfa` exists and passes the two named unit tests (`src/gfa/mod.rs:229-329`), `docs/STRATEGIC_PLAN.md:13` infers from this that *"native streaming gives Horizon direct GFA v1.0/v1.1 interoperability without external `vg`/`odgi` dependencies."* That inference is unsupported on three counts:
  1. `stream_gfa_reader` (`src/gfa/mod.rs:213`) silently discards standard GFA v1.0/v1.1 record types (`_ => {} // Ignore custom tags/lines like C, J, U, O`), where `C` (Containment) and `J` (Jump) are part of the GFA v1.0/v1.1 specifications.
  2. `stdin` (`"-"`, `src/gfa/mod.rs:112-116`) uses plain `BufReader::new(stdin.lock())` without `MultiGzDecoder`, so piping compressed `.gfa.gz` streams via `"-"` fails.
  3. Neither `"-"` `stdin` nor `C`/`J` record handling is covered by any unit test in `src/gfa/mod.rs:229-329` (the tests only exercise `Cursor` and a 2-member `GzEncoder` file). (Note: Round 1 and Round 2 already adjudicated `stream_gfa`'s non-zero-copy allocation and two-pass `stdin` incompatibility in downstream tools.)
* **Evidence:** `src/gfa/mod.rs:112-127, 213, 229-329`; `docs/STRATEGIC_PLAN.md:13`.
* **Correction:** Narrow `docs/STRATEGIC_PLAN.md:13` to state that `stream_gfa` supports a GFA v1.0/v1.1 subset (`H`/`S`/`L`/`P`/`W`, silently ignoring `C`/`J`), and that `"-"` `stdin` is uncompressed-only and untested.

---

## F2 -- C2: Oracle `(ref_id, tgt_id)` pairing in `horizon_eval` and repeated re-runs invalidate precision and held-out inferences

* **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:16, 47`)
* **Verdict:** **NARROWED**
* **Concrete flaw:**
  1. **Oracle contig-pairing confounder:** `horizon_eval` (`src/bin/horizon_eval.rs:151-162`) iterates *only* over `((ref_id, tgt_id), maf_set)` entries already present in `maf_by_pair` (the truth MAF), and `docs/BENCHMARK_RESULTS.md:5` notes that *"MAF-extracted FASTAs sample homologous blocks."* Because non-homologous scaffold pairs across the 129 Dog target scaffolds are never aligned against `hg38.chr22`, the reported `60.3%` precision (`61.3%` F1) is conditioned on an oracle truth-MAF contig pairing and cannot be inferred as realistic multi-scaffold alignment precision.
  2. **Held-out contamination & regression:** `docs/BENCHMARK_RESULTS.md:38-47` places `cactus241way (chr22 1MB)` under **Held-Out Test Sets** (*"completely locked away during algorithm development to prevent overfitting"*), yet it has been repeatedly re-evaluated across development commits (`b857377`, `b18f25b`, `488b020`). On `488b020`, `cactus241way (chr22 1MB)` *regressed* from `63.6%` recall / `61.0%` precision (`62.3%` F1, `15.11s` on `b18f25b`) to `62.3%` recall / `60.3%` precision (`61.3%` F1, `17.16s`), while the other five WGA rows were left at pre-`FlatSeedIndex` commits.
* **Evidence:** `src/bin/horizon_eval.rs:148-162`; `docs/BENCHMARK_RESULTS.md:5, 16, 38-47`; `git show 488b020 -- docs/BENCHMARK_RESULTS.md`.
* **Correction:** Qualify `C2` in `docs/BENCHMARK_RESULTS.md:16, 47` as oracle-paired, MAF-extracted block precision/recall (not all-vs-all 129-scaffold precision), remove the "completely locked away" held-out claim for `cactus241way (chr22 1MB)`, and note the `-1.0 pp` F1 / `+2.05s` regression on `488b020`.

---

## F3 -- C3: Helper-only binary unit tests do not license pipeline correctness or `R2-RESOLVE-LIFT` discharge

* **Claim:** `C3` (`docs/STRATEGIC_PLAN.md:3`)
* **Verdict:** **NARROWED**
* **Concrete flaw:** While all 9 binaries in `src/bin/` now have at least one `#[test]` (`45` passing + `1` ignored; `clippy` clean), using this to infer that the pangenome suite is verified or that `R2-RESOLVE-LIFT` is discharged (commit `488b020`) is unsupported:
  1. The newly added binary tests in `src/bin/bench_chr1.rs:129-137` (`revcomp` + `read_first_contig` on 8 bp), `src/bin/horizon_eval.rs:188-212` (`read_named_contigs` + `parse_maf_pairs_by_contig` on 3 bp), and `src/bin/horizon_merge.rs:52-66` (`merge_pairwise_mafs` on 10 bp) only test small helper functions on toy strings; none of them invoke `main()`, CLI argument parsing, or end-to-end binary pipelines.
  2. `docs/STRATEGIC_PLAN.md:3` itself concedes that `horizon_clean`, `horizon_resolve`, and `horizon_lift` still require seekable file paths rather than `stdin` streaming (contradicting Section 3B, `docs/STRATEGIC_PLAN.md:23`), `horizon_merge` has no multi-species benchmark (`docs/BENCHMARK_RESULTS.md:9`), and `R2-RESOLVE-LIFT`'s required head-to-head baselines (`minimap2`/`odgi`/`vg`) remain completely absent.
* **Evidence:** `src/bin/bench_chr1.rs:129-137`; `src/bin/horizon_eval.rs:188-212`; `src/bin/horizon_merge.rs:52-66`; `docs/STRATEGIC_PLAN.md:3, 23`; `docs/BENCHMARK_RESULTS.md:9`.
* **Correction:** Qualify `docs/STRATEGIC_PLAN.md:3` to state that the binary unit tests are helper-level smoke tests (not end-to-end CLI/pipeline tests) and that `R2-RESOLVE-LIFT` baselines (`minimap2`/`odgi`/`vg`) and `stdin` streaming for `horizon_clean`/`horizon_resolve`/`horizon_lift` remain open.

---

## F4 -- C4: RSS reduction is confounded by bundled `488b020` algorithmic changes, and `10.09 GB` is an unfinished mid-run sample

* **Claim:** `C4` (`docs/STRATEGIC_PLAN.md:11`; `docs/BENCHMARK_RESULTS.md:26, 60, 62`)
* **Verdict:** **NARROWED**
* **Concrete flaw:**
  1. **Confounded causal attribution (no ablation):** `C4` attributes the `20 Mbp chr6 q-arm` RSS drop (`5,863.3 MB -> 2,409.7 MB`) and `chr6` 5 Mbp pericentromere RSS drop (`3,084.5 MB -> 1,048.2 MB`) solely to `FlatSeedIndex`. However, commit `488b020` bundled `FlatSeedIndex` with `rev_offsets` fixes (`src/seeding/mod.rs:96`), `bocpd.rs` N-base handling, `mdl_xdrop_extend_with_params` (`src/mdl/dp.rs`), and `horizon_wga.rs` changes. Across those two commits (`b18f25b` -> `488b020`), `chr6` pericentromere aligned pairs dropped by **11.3%** (`4,785,880 -> 4,247,511`), `chr22` pericentromere pairs dropped by **13.4%** (`5,894,191 -> 5,106,541`), and `chr6 q-arm 20 Mbp` bubbles dropped by **21.4%** (`21,475 -> 16,877`). Without an ablation holding the alignment algorithm fixed, the RSS reduction cannot be attributed purely to `FlatSeedIndex`.
  2. **Transient build footprint & incomplete whole-chromosome evidence:** `FlatSeedIndex::build` (`src/seeding/mod.rs:385-450`) allocates `pairs: Vec<(u64, u32)>` (`16 bytes/entry` with alignment padding) simultaneously with `keys: Vec<u64>` + `pos: Vec<u32>` (`12 bytes/entry`) before dropping `pairs`. Furthermore, `docs/BENCHMARK_RESULTS.md:26` cites *"a mid-run check at 13:59 used 10.09 GB RSS"* on full `GRCh38.chr6` vs `CHM13.chr6`—an unfinished mid-run sample cannot establish peak RSS or discharge `R2-INDEX-RAM`.
* **Evidence:** `src/seeding/mod.rs:96, 385-450`; `git show 488b020 --stat`; `git show 488b020 -- docs/BENCHMARK_RESULTS.md`; `docs/BENCHMARK_RESULTS.md:26, 59-62`.
* **Correction:** Attribute the RSS reduction to the `488b020` build (not `FlatSeedIndex` in isolation) alongside the `-11.3%` to `-13.4%` pair-count shift, and remove the unfinished `10.09 GB` mid-run check as evidence of full-chromosome peak RSS until a completed run is measured.

---

## F5 -- C5: Pericentromere coverage/pair counts license no inference of alignment accuracy or structural resolution

* **Claim:** `C5` (`docs/BENCHMARK_RESULTS.md:23, 59`)
* **Verdict:** **NARROWED**
* **Concrete flaw:** Reporting `5,106,541` many-to-many `(x, y)` aligned pairs (`unique_x_bp=4,387,655`, `unique_y_bp=4,549,858`) and `16,083` Horizon Bubbles on a pre-sliced 5 Mbp-vs-5 Mbp window (`GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M`) licenses no inference of alignment accuracy, orthology, or pericentromere resolution quality:
  1. There is no truth alignment, no precision/recall/F1, and no baseline comparison (`minimap2`, `wfmash`, `lastz`, `odgi`, `vg`). High base coverage (`87.8%` / `91.0%`) in a satellite-rich pericentromere cannot distinguish true orthologous/paralogous alignment from spurious repeat-induced alignment.
  2. Between `b18f25b` and `488b020`, `chr22` pericentromere aligned pairs dropped by `787,650` (`-13.4%`), yet without a truth oracle it is impossible to infer whether this was a false-positive reduction or a sensitivity loss (mirroring the `-1.3 pp` recall drop on `cactus241way` in the same commit).
  3. Pre-slicing homologous 5 Mbp windows bypasses the genome-wide localization problem in pericentromeres.
* **Evidence:** `docs/BENCHMARK_RESULTS.md:23, 59`; `git show 488b020 -- docs/BENCHMARK_RESULTS.md`.
* **Correction:** Explicitly label `C5` in `docs/BENCHMARK_RESULTS.md:23, 59` as unvalidated throughput/coverage on a pre-localized 5 Mbp window with no accuracy oracle or baseline.

---

## F6 -- C6: `horizon_gfa stats` (`97.1 MB` RSS) does not license pangenome graph-processing scalability

* **Claim:** `C6` (`docs/BENCHMARK_RESULTS.md:19, 54`; `docs/STRATEGIC_PLAN.md:13`)
* **Verdict:** **NARROWED**
* **Concrete flaw:** `horizon_gfa stats` (`src/bin/horizon_gfa.rs:39-77`) only increments five `u64` counters (`s_count`, `s_bases`, `l_count`, `p_count`, `w_count`) and discards every parsed record immediately. Citing its `97.1 MB` peak RSS in `docs/STRATEGIC_PLAN.md:13` to infer HPRC-scale GFA interoperability without `vg`/`odgi` is a category error: as soon as records are actually used, RSS grows to `456.1 MB` (`horizon_gfa extract` `chr22`), `1,255.1 MB` (`extract` `chr6`), `619.2 MB` (`horizon_lift` on 500 `chr22` windows), and `2,409.7 MB` (`horizon_wga` `chr6` 20 Mbp), while `horizon_wga --graph` has only been tested on a 4-node diamond GFA (`docs/BENCHMARK_RESULTS.md:65`). Moreover, `88.84s` has no baseline (`gfatools stat`, `odgi stats`, or `zcat | wc -l`).
* **Evidence:** `src/bin/horizon_gfa.rs:39-77`; `docs/BENCHMARK_RESULTS.md:19, 54-65`; `docs/STRATEGIC_PLAN.md:13`.
* **Correction:** Narrow `docs/STRATEGIC_PLAN.md:13` and `docs/BENCHMARK_RESULTS.md:19` to clarify that `97.1 MB` measures record-discarding line counting (`stats`), not graph indexing, extraction, or alignment memory.

---

## F7 -- Claim-list steering: runner omitted no-op/synthetic pangenome rows, missing baselines, and stale pre-`488b020` WGA rows

* **Claim:** Framing / claim selection (`BRIEF-inference.md`)
* **Verdict:** **NARROWED** *(framing finding)*
* **Concrete flaw:** The runner-authored claim list `C1`–`C6` steers attention toward passing smoke tests (`C3`), record-counting RSS (`C6`), and unvalidated pericentromere pair counts (`C5`), while omitting the weakest rows in `docs/BENCHMARK_RESULTS.md`:
  1. `horizon_clean` on real `chr22`/`chr6` pericentromere GFAs (`docs/BENCHMARK_RESULTS.md:63`) is a **complete no-op** (`6,030 / 2,776 kept, 0 severed (no-op)`).
  2. `horizon_resolve --collapse-dups` (`docs/BENCHMARK_RESULTS.md:64`) is still benchmarked only on a synthetic 200-segment `ACGT` repeat with 20 identical walks.
  3. `horizon_wga --graph` (`docs/BENCHMARK_RESULTS.md:65`) is still benchmarked only on a 4-node diamond GFA + 96-bp query.
  4. `minimap2`/`odgi`/`vg` baselines (`R2-RESOLVE-LIFT`) remain completely absent, and 5 of 6 WGA rows (`multiz100way`, `cactus447way`, `multiz470way`, `cactus241way chr22_random`, `GRCh38.chr22 vs CHM13.chr22`) remain stale on pre-`488b020` commits (`docs/BENCHMARK_RESULTS.md:5`).
* **Evidence:** `docs/BENCHMARK_RESULTS.md:5, 9, 63-65`; `reports/rounds/round-3/CARRIED.md`.
* **Correction:** Require future claim lists and `docs/BENCHMARK_RESULTS.md` summaries to surface the no-op/synthetic status of `horizon_clean`, `horizon_resolve`, and `horizon_wga --graph`, and to re-run stale WGA rows on a single pinned commit.

---

## Verdict Summary

| Claim | Verdict |
| :--- | :--- |
| C1 | NARROWED |
| C2 | NARROWED |
| C3 | NARROWED |
| C4 | NARROWED |
| C5 | NARROWED |
| C6 | NARROWED |
