# Round 1 Adversarial Attack — L2 (Sampling)

**Attacker lens:** `sampling`
**Targets:** `docs/STRATEGIC_PLAN.md`, `docs/PANGENOME_SUITE.md` (claims `C1`–`C11`)
**Prior rounds / method cards checked:**
- `scripts/round history "" --kind adjudication`: no prior closed rounds exist (`round 1` is the first round).
- No statistical method cards or held-literature appendix were declared in `BRIEF-sampling.md` or found in the repository (`docs/`, `scripts/`).

---

## F1 -- C8 REFUTED: Target-market and benchmark-strategy claims contradict actual benchmarked samples and missing pangenome binaries

* **Claim:** `C8` (`docs/STRATEGIC_PLAN.md:31`) — *"The HPRC produces massive graphs that researchers cannot use because standard tooling crashes or produces tangled paths. We will focus our benchmarking and release strategy entirely on `horizon_clean`, `horizon_resolve`, and `horizon_lift` first to capture the pangenome community immediately."*
* **Verdict:** `REFUTED`
* **Concrete flaw & evidence:**
  1. **Zero HPRC or pangenome-tool benchmarks exist:** In `docs/BENCHMARK_RESULTS.md:20-39`, 100% of reported benchmark numbers (7 rows across lines `25-26` and `33-37`) are pairwise whole-genome alignment (`horizon_wga` / `horizon_eval`) runs on UCSC MAF slices (`multiz100way`, `cactus447way`, `multiz470way`, `cactus241way`). The only pangenome row (`docs/BENCHMARK_RESULTS.md:39`, `HPRC (v1.1) | GRCh38 | CHM13 | Human Pangenome | *(Graph Eval)* | - | -`) has no measured recall, precision, bubble count, or runtime, and `B10K 363-way` (line `38`) is `*(Pending Data)*`.
  2. **The evaluation harness excludes HPRC and all three named pangenome tools:** `scripts/benchmarks/run_all_evals.py:56-76` runs only `horizon_eval` on pairwise FASTA/MAF pairs (`multiz100way` chrM, `cactus447way` chr22_random, `multiz470way` chr22_random). None of `horizon_clean`, `horizon_resolve`, or `horizon_lift` exist in `src/bin/` (which contains only `bench_chr1.rs`, `horizon_eval.rs`, `horizon_guide.rs`, `horizon_merge.rs`, and `horizon_wga.rs`).
  3. **HPRC dataset version mismatch & missing local data:** `docs/BENCHMARK_DATASETS.md:46` specifies `hprc-v1.1-mc-grch38.gfa.gz`, whereas `scripts/benchmarks/fetch_benchmarks.sh:27` downloads `hprc-v1.0-mc-grch38.gfa.gz`, and `data/benchmarks/` is absent from the workspace.
  4. **Zero sample of standard-tooling failures:** No benchmark or dataset in `docs/BENCHMARK_DATASETS.md` or `docs/BENCHMARK_RESULTS.md` measures `vg` or `odgi` crashes or tangled-path counts on HPRC.
* **Correction:** Replace `docs/STRATEGIC_PLAN.md:31` with an accurate scope statement acknowledging that all current benchmarks sample pairwise UCSC MAF slices via `horizon_wga`/`horizon_eval`, and that `horizon_clean`, `horizon_resolve`, `horizon_lift`, and HPRC graph evaluations are unimplemented/unbenchmarked.

---

## F2 -- C4 REFUTED: Zero 12GB pangenome or vg/odgi crash samples are benchmarked, and stream_gfa has zero callers

* **Claim:** `C4` (`docs/STRATEGIC_PLAN.md:11`) — *"`vg` and `odgi` crash on 12GB pangenomes due to RAM exhaustion. Our custom zero-copy GFA streaming parser (`src/gfa/mod.rs`) solves the I/O bottleneck entirely."*
* **Verdict:** `REFUTED`
* **Concrete flaw & evidence:**
  1. **No 12GB pangenome sample tested:** Although `scripts/benchmarks/fetch_benchmarks.sh:26-27` includes a curl command for `hprc-v1.0-mc-grch38.gfa.gz` ("12 GB"), `docs/BENCHMARK_RESULTS.md:39` reports no runtime, memory usage, or completion result for HPRC (`*(Graph Eval)* | - | -`), and `scripts/benchmarks/run_all_evals.py` never invokes any GFA parser.
  2. **Parser is uncalled across the entire codebase:** `src/gfa/mod.rs:50-124` defines `stream_gfa`, but searching `src/` shows zero callers in any binary, library function, or test (it is only declared via `pub mod gfa;` in `src/lib.rs:394`).
  3. **No comparative failure sample:** There is no benchmark artifact in `docs/BENCHMARK_RESULTS.md` or `scripts/benchmarks/` sampling `vg` or `odgi` RAM exhaustion on 12GB graphs.
* **Correction:** Remove the unverified claim in `docs/STRATEGIC_PLAN.md:11` that `src/gfa/mod.rs` solves the 12GB pangenome I/O bottleneck until `stream_gfa` is wired into a CLI binary and benchmarked on an actual multi-gigabyte GFA sample alongside `vg`/`odgi`.

---

## F3 -- C2 NARROWED: Benchmark sampling pre-filters to homologous MAF blocks on tiny slices and omits graph/hairball evaluations

* **Claim:** `C2` (`docs/STRATEGIC_PLAN.md:9`) — *"We are replacing arbitrary affine gap penalties and scoring matrices with Information-Theoretic bit-savings (MDL). This inherently prevents the "hairball" over-alignment problem."*
* **Verdict:** `NARROWED`
* **Concrete flaw & evidence:**
  1. **Selection bias in pairwise input extraction:** `scripts/benchmarks/extract_fasta_from_maf.py:28-45` builds the test FASTAs (`ref_fa`, `tgt_fa`) by iterating over the truth MAF and retaining **only** blocks where both `ref_prefix` and `target_prefix` are present (`if ref_row and tgt_row:`), concatenating those homologous segments into contiguous sequences starting at coordinate `0`. Non-homologous, repetitive, and unaligned inter-block regions—the primary substrate where over-alignment ("hairball") occurs—are excluded from the benchmark sample before `horizon_eval` (`src/bin/horizon_eval.rs:86-115`) runs.
  2. **Tiny, unrepresentative genomic slices:** `scripts/benchmarks/fetch_benchmarks.sh:7-24` and `docs/BENCHMARK_DATASETS.md:8-26` sample only mitochondrial DNA (`multiz100way` `chrM.maf.gz`, ~16.5 kb) and a ~38 kb unplaced contig (`chr22_KI270735v1_random` for `cactus447way`, `multiz470way`, and `cactus241way`).
  3. **Zero graph/hairball metrics:** `docs/BENCHMARK_RESULTS.md:23-39` reports only pairwise MAF recall/precision and bubble counts, with no graph topology or over-alignment/hairball benchmark on any pangenome.
* **Correction:** Narrow `docs/STRATEGIC_PLAN.md:9` to state that MDL bit-savings has only been evaluated on pre-filtered homologous MAF blocks from `chrM` and `chr22_KI270735v1_random` slices, not on full unmasked chromosomes or pangenome graphs.

---

## F4 -- C3 NARROWED: Named HLA and centromere paralog loci are never sampled in benchmark datasets or tests

* **Claim:** `C3` (`docs/STRATEGIC_PLAN.md:10`) — *"Using congruence-closure to resolve complex, nested structural paralogs (like the HLA region). This natively models evolutionary trajectories without producing unnavigable cyclical graphs."*
* **Verdict:** `NARROWED`
* **Concrete flaw & evidence:**
  1. **Zero HLA (`chr6`) or centromere sampling:** `docs/BENCHMARK_DATASETS.md:8-48` and `scripts/benchmarks/run_all_evals.py:59-75` sample `chrM`, `chr22_KI270735v1_random`, `chr22 1MB` (`docs/BENCHMARK_RESULTS.md:37`), Ensembl `92_mammals` (`1_14`), `sars-cov-2.hal`, `b10k_363way`, and an unrun HPRC GFA. Neither `chr6` (the HLA locus) nor any centromeric region is sampled or evaluated in `docs/BENCHMARK_RESULTS.md:23-39`.
  2. **No HLA/centromere test in `src/`:** Searching `src/` for `hla` or `centromer` returns zero matches (`docs/STRATEGIC_PLAN.md:10` and `docs/PANGENOME_SUITE.md:24` are the only occurrences in the repository).
* **Correction:** Narrow `docs/STRATEGIC_PLAN.md:10` (and `docs/PANGENOME_SUITE.md:24`) to remove the unverified assertion that the pipeline resolves HLA or centromeric paralogs until an HLA/`chr6` or centromere dataset is actually sampled and benchmarked.

---

## F5 -- Claim-list steering: Runner's claim list omits distant-relative and chr22 1MB benchmark-collapse sampling flaws

* **Evidence:** `BRIEF-sampling.md:135-138` instructs attackers to flag when the runner-authored claim list steers away from the weakest point.
* **Concrete flaw & evidence:**
  1. `docs/PANGENOME_SUITE.md:12-17` claims that `horizon_wga --graph` *"Enables the direct mapping of distant relatives (e.g., archaic hominids, distant primates) into established pangenomes"*, and `docs/BENCHMARK_RESULTS.md:14-16` highlights *"Distant divergence"* results (`cactus447way` Horse 47.3% recall; `multiz100way` chrM Mouse 79.1% recall) in its Executive Summary.
  2. However, the held-out table in `docs/BENCHMARK_RESULTS.md:36-37` shows that on `cactus241way` (Human vs Dog, `Distant Mammal`), moving from the 38 kb `chr22_random` slice (`45.3%` recall) to `chr22 1MB` causes recall to collapse to **`9.9%`** (with `494` Horizon Bubbles)—a result omitted from the Executive Summary (`docs/BENCHMARK_RESULTS.md:14-16`) and absent from the `C1`–`C11` claim list.
  3. Similarly, `scripts/benchmarks/run_all_evals.py:59-63` calibrates (`run_tuning`) exclusively on `multiz100way` `chrM` (mitochondrial DNA, ~16.5 kb), while `run_testing` (lines `68-75`) evaluates on nuclear `chr22_KI270735v1_random` (~38 kb), creating a compartment/scale mismatch between tuning and test samples.

---

## Claim Verdicts (Sampling Lens)

| C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| UNCHANGED | NARROWED | NARROWED | REFUTED | UNCHANGED | UNCHANGED | UNCHANGED | REFUTED | UNCHANGED | UNCHANGED | UNCHANGED |
