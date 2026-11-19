# Round 4 Adjudication — L2 (`sampling`)

**Target attack:** [`reports/rounds/round-4/attack-L2-sampling.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-4/attack-L2-sampling.md) (sealed SHA-256 prefix `c99917ac727c`, verified).
**Target manuscript & codebase:** [`manuscript/horizon_paper.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md) (claims `C1`–`C6`, `HEAD` `33bede8`), [`docs/BENCHMARK_RESULTS.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md), [`manuscript/scripts/generate_figures.py`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/scripts/generate_figures.py), [`scripts/benchmarks/`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks), and `src/`.
**Method cards & prior rounds (Rules 5–6):** This repository has no method cards in `refs/methods/`. Checked `scripts/round history` for `C1`, `horizon_clean`, `chr6`, and `30`; prior rounds adjudicated `horizon_clean`'s real-data no-op (`r1` `F9`, `r2`/`r3` `R2-RESOLVE-LIFT`) and missing `minimap2`/`odgi`/`vg`/`LASTZ` baselines (`r2` `R2-RESOLVE-LIFT`, `r3` `F7`), which the attack properly discloses.

---

## F1 — C2 samples only the CHM13.chr6 q-arm (`110.7 Mbp`), `n = 1` chromosome and `n = 1` run, with no accuracy metric and stale pre-`33bede8` timings

- **Ruling:** **UPHELD** (attacker's `NARROWED` verdict on `C2` stands).
- **What survives:** `C2` overstates a single-run (`53130b8`) q-arm-only `CHM13.chr6` benchmark as "full `GRCh38.chr6` × `CHM13.chr6`" without replicates, accuracy/coverage metrics, or re-measurement on `HEAD` (`33bede8`).
- **Verification:**
  1. `/tmp/hprc_chm13_chr6.fa` begins with `>CHM13.chr6:61401689-172126628` (`110,724,939` bp; `35.7%` of `172.1 Mbp` `CHM13.chr6` is absent), whereas `manuscript/horizon_paper.md:10, 101, 105` calls it "full `GRCh38.chr6` × `CHM13.chr6`" (unlike `chr22` at line 101, which qualifies `CHM13.chr22` q-arm).
  2. [`docs/BENCHMARK_RESULTS.md:27, 64-65`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L64-L65) confirms `n = 1` chromosome pair and `n = 1` run with no recall/precision/unique-base coverage, and shows `1e47fe0` (`292,534` fwd / `291,150` rev LCBs, `693,242` DP jobs) vs. `53130b8` (`301,156` fwd / `298,178` rev LCBs, `778,471` DP jobs) had different workloads.
  3. `git show --stat 33bede8` confirms `HEAD` subsequently modified `src/seeding/mod.rs`, `src/seeding/flat_index.rs`, `src/seeding/translated.rs`, `src/mdl/dp.rs`, and `src/mdl/mod.rs` without re-benchmarking `chr6`.
- **Repair:** Endorse the text-only qualification in `manuscript/horizon_paper.md:10, 101, 105` to `GRCh38.chr6` (`170.8 Mbp`) × `CHM13.chr6:61.4M–172.1M` q-arm (`110.7 Mbp`) and `n = 1` single-run timing on `53130b8` (**the fix is free**); re-benchmarking on `33bede8` with replicates and unique-base coverage is non-free.

---

## F2 — C3 `horizon_gfa extract` samples only `2` chromosomes and `2` haplotypes out of `2,580` paths and `24,456` walks

- **Ruling:** **UPHELD** (attacker's `NARROWED` verdict on `C3` stands).
- **What survives:** `horizon_gfa extract` was measured on `n = 1` run per pair for only `2`-haplotype (`GRCh38` + `CHM13`) extractions (`chr22`, `chr6`), with RSS scaling with cached segments and `chr6` pinned prior to `33bede8` edits to `cmd_extract`.
- **Verification:**
  1. [`docs/BENCHMARK_RESULTS.md:55-57`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L55-L57) confirms only `GRCh38` + `CHM13` `chr22` (`921,252` cached segs, `456.4 MB` RSS) and `chr6` (`3,057,602` cached segs, `1,255.1 MB` RSS) were sampled, and `chr6` is pinned to `b857377`.
  2. `git show 33bede8 -- src/bin/horizon_gfa.rs src/gfa/mod.rs` confirms `33bede8` modified `cmd_extract` (adding `SpoolableGfaInput::new` and `matches_sample_prefix`).
- **Repair:** Endorse qualifying `manuscript/horizon_paper.md:101` to state `n = 1` run per pair for `2`-haplotype (`GRCh38` + `CHM13`) `chr22` and `chr6` extractions with RSS scaling with cached segment count (**the fix is free**); re-measuring on `33bede8` is non-free.

---

## F3 — C4 `horizon_clean` has zero empirical sampling on sequence-bearing GFAs (only a `0`-severed `*-seq` no-op and 2 toy repeat tests)

- **Ruling:** **UPHELD** (attacker's `NARROWED` verdict on `C4` stands).
- **What survives:** `horizon_clean` has zero empirical evaluation on sequence-bearing real or pangenome GFAs, resting only on a `0`-severed `*-seq` no-op and two synthetic 60-bp AT/GC repeat unit tests.
- **Verification:**
  1. [`docs/BENCHMARK_RESULTS.md:66`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L66) records only `2,776 kept, 0 severed (no-op on *-seq)`.
  2. [`src/bin/horizon_clean.rs:100-127, 154-176, 252-260`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_clean.rs#L100-L127) returns `0.0` for every link when `S`-segment sequences are `*`, and [`src/bin/horizon_clean.rs:559-626`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_clean.rs#L559-L626) exercises edge severing only on two 60-bp `ATATAT...` vs. `GCGCGC...` unit-test fixtures.
- **Repair:** Endorse clarifying in `manuscript/horizon_paper.md:49` that `horizon_clean`'s junction-severing behavior is currently verified only on synthetic 60-bp AT/GC unit-test fixtures and has not yet been benchmarked on sequence-bearing HPRC or pangenome GFAs (**the fix is free**).

---

## F4 — C5 rests on untracked *Plasmodium* / *Cannabis* samples with an internal denominator contradiction (`109,565` vs `119,582`) and MAF-block selection bias

- **Ruling:** **UPHELD** (attacker's `REFUTED` verdict on `C5` stands).
- **What survives:** `C5`'s *Plasmodium* (`87.7%` abstention, `54.4%`/`66.8%` recovery) and *Cannabis* (`418,663` columns, `99.95%` retention) figures are absent from `src/`, `scripts/`, and `docs/`, contradict `generate_figures.py` by `10,017` units (`119,582` vs. `109,565`), and condition on MAF-aligned blocks.
- **Verification:**
  1. `grep -rn "Cherry_Pie\|418,663\|418465\|87\.7" src/ scripts/ docs/` returns `0` hits.
  2. [`manuscript/horizon_paper.md:65`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L65) reports `24,748 / 119,582` (`20.7%`) MDL-rejected columns, whereas [`manuscript/scripts/generate_figures.py:218-219`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/scripts/generate_figures.py#L218-L219) plots `n = 88,597` + `n = 20,968` = `109,565` pairs (`10,017` discrepancy).
  3. [`scripts/benchmarks/extract_fasta_from_maf.py:41-45, 78-85`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/extract_fasta_from_maf.py#L41-L45) confirms the project's MAF pipeline concatenates only MAF-aligned blocks, excluding unaligned intergenic/repetitive sequence.
- **Repair:** Re-running, tracking, and reconciling the *Plasmodium* and *Cannabis* benchmarks on `33bede8` is non-free (removing the unverified claims from `manuscript/horizon_paper.md:10, 65` and `fig3` is a free text removal if the authors choose removal over re-execution).

---

## F5 — C6 `30` randomized tandem-duplication / overlapping-inversion trials do not exist in the codebase

- **Ruling:** **UPHELD** (attacker's `REFUTED` verdict on `C6` stands).
- **What survives:** No 30-trial randomized simulation or `0.00` greedy one-to-one baseline exists anywhere in the repository; the codebase has only two deterministic unit tests and one synthetic GFA run.
- **Verification:**
  1. `grep -rn "greedy\|30" src/egraph.rs src/bin/horizon_resolve.rs` and `grep -rn "tandem\|egraph\|greedy" scripts/` return `0` hits.
  2. [`src/egraph.rs:463-518`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/egraph.rs#L463-L518) contains only `2` deterministic unit tests (`[A, -C, -B, -B, D]` and `[B, B]`), and [`docs/BENCHMARK_RESULTS.md:67`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L67) reports only `1` synthetic `horizon_resolve --collapse-dups` run.
- **Repair:** Endorse removing the unsubstantiated 30-trial / `0.00` greedy-baseline sentence from `manuscript/horizon_paper.md:99` (**the fix is free**); implementing and running a 30-trial randomized harness is non-free.

---

## F6 — Claim-list steering: `C1`–`C6` omit unrun/missing benchmarks, `horizon_merge` provenance gap, and external baseline absences

- **Ruling:** **UPHELD** (framing finding).
- **What survives:** `C1`–`C6` omit `horizon_merge`'s missing benchmark row (`docs/BENCHMARK_RESULTS.md:9` vs. 7 merge timings in `horizon_paper.md`), the pending `B10K 363-way` dataset, and recurring external-baseline (`vg` present at `/home/rcb112/project/gemini/bench/vg`) and MAF-block-extraction gaps.
- **Verification:**
  1. `/home/rcb112/project/gemini/bench/vg` exists (`46,259,960` bytes) while [`docs/BENCHMARK_RESULTS.md:5`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L5) states head-to-head baselines remain pending external tool availability.
  2. [`docs/BENCHMARK_RESULTS.md:9`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L9) states *"no multi-species benchmark row yet"*, while `manuscript/horizon_paper.md:10, 63, 75, Table 1` reports 7 `horizon_merge` timings.
  3. [`docs/BENCHMARK_RESULTS.md:49`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L49) is `*(Pending Data)*`, and [`scripts/benchmarks/run_all_evals.py:39`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/run_all_evals.py#L39) invokes `extract_fasta_from_maf.py`.
- **Note:** Items 1 and 4 are recurring from rounds 2–3 (`r2` `R2-RESOLVE-LIFT`, `r3` `F7`) and disclosed in `docs/BENCHMARK_RESULTS.md:5`; item 2 (`horizon_merge` provenance gap) is a new Round 4 manuscript-vs-benchmark discrepancy.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `/tmp/hprc_chm13_chr6.fa` is q-arm only (`61.4M–172.1M`), `n=1` run/chromosome with no accuracy metric, and `33bede8` edited seeding/MDL without re-benchmarking; free text repair endorsed. |
| F2 | UPHELD | `horizon_gfa extract` sampled only `2` chromosomes and `2` haplotypes with RSS scaling with cached segments, and `chr6` predates `33bede8` `cmd_extract` edits; free text repair endorsed. |
| F3 | UPHELD | `horizon_clean` has `0` sequence-bearing real GFA benchmarks (`docs/BENCHMARK_RESULTS.md:66` is a `*-seq` no-op; unit tests are 2 synthetic 60-bp AT/GC fixtures); free text repair endorsed. |
| F4 | UPHELD | *Plasmodium*/*Cannabis* evaluation numbers are absent from `src/`/`scripts/`/`docs/`, `horizon_paper.md:65` (`119,582`) contradicts `generate_figures.py:218-219` (`109,565`), and MAF-block extraction strips unaligned sequence. |
| F5 | UPHELD | `30` randomized duplication/inversion trials and `0.00` greedy baseline do not exist in `src/egraph.rs`, `src/bin/horizon_resolve.rs`, or `scripts/`; free text removal endorsed. |
| F6 | UPHELD | `horizon_merge` timings in `horizon_paper.md` have no benchmark row in `docs/BENCHMARK_RESULTS.md:9`, `B10K` is pending, `vg` is on disk unused, and `run_all_evals.py` uses MAF-block extraction. |
