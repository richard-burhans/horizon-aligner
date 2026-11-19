# Round 3 Attack -- Statistics Lens (L3)

- **Method cards (`refs/methods/`):** Checked `/home/rcb112/project/gemini/bench/horizon-aligner/refs/methods/` and `/home/rcb112/project/gemini/bench/refs/methods/`; neither exists, so no project method card governs these statistical claims.
- **Prior round history:** Checked `reports/adversarial_log.md`, `reports/rounds/round-2/adjudication-L3.md`, and `scripts/round history` for `cactus241way`, `FlatSeedIndex`, `horizon_gfa`, `horizon_eval`, and `extract_fasta_from_maf`. Prior adjudications relevant to each finding are disclosed below alongside what changed in commit `488b020`.

---

## F1 -- C2: `horizon_eval` oracle-pairs target scaffolds with zero cross-scaffold false positives, and `estimate_hmm_params` remains truncated/clamped while `cactus241way` regressed

- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:16, 47`) -- *"UCSC cactus241way `chr22 1MB`: **62.3%** recall / **60.3%** precision (Dog, across 129 target scaffolds; **61.3% F1** in **17.16s**)."*
- **Verdict:** **NARROWED**
- **Prior round context:** Round 2 L3-F1 (`UPHELD`) documented that `estimate_hmm_params` (`src/mdl/mod.rs:732-816`) truncates at `4000` bases, clamps `identity_prior` to `[0.60, 0.72]`, and collapses multi-indel intervals (`dx != dy`), and Round 1/2 L2 noted `extract_fasta_from_maf.py` homologous-block extraction. In `488b020`, `mdl_xdrop_extend_with_params` was wired in, but `estimate_hmm_params` itself was not changed.
- **Concrete flaws:**
  1. **Oracle contig pairing masks cross-scaffold false positives (`src/bin/horizon_eval.rs:148-162`):** Although C2 cites `"across 129 target scaffolds"`, `horizon_eval` does not align the query against all 129 scaffolds competitively. Instead, `for ((ref_id, tgt_id), maf_set) in &pairs` iterates exclusively over `(ref_id, tgt_id)` pairs present in the truth MAF (`maf_by_pair`) and aligns each known-homologous pair in isolation via `align_contig_pair_with_config(idx_r, idx_t, &cfg)`. Non-paired target scaffolds are never tested, making cross-scaffold false positives `0` by construction.
  2. **Unweighted micro-averaging without per-scaffold dispersion (`src/bin/horizon_eval.rs:165-176`):** `recall` (`total_overlap / total_maf`) and `precision` (`total_overlap / total_proto`) are micro-averaged totals across the 129 scaffolds with no per-scaffold median, IQR, or confidence interval; `61.3% F1` is the harmonic mean of those two point estimates (`2 * 62.3 * 60.3 / (62.3 + 60.3) = 61.28%`).
  3. **Unchanged `estimate_hmm_params` biases & post-`488b020` regression (`src/mdl/mod.rs:747-809`, `docs/BENCHMARK_RESULTS.md:5, 16`):** `estimate_hmm_params` still samples `lcb.iter().take(64)` up to `total_bases >= 4000`, clamps `identity_prior` to `[0.60, 0.72]` around `Beta(195, 105)`, and counts `dx != dy` over `lcb.windows(2).take(64)` (`dx, dy <= 500`). After `488b020` wired `HmmParams` into `mdl_xdrop_extend_with_params`, `cactus241way chr22 1MB` regressed from `63.6%` recall / `61.0%` precision (`62.3%` F1, `15.1s`) in Round 2 to `62.3%` recall / `60.3%` precision (`61.3%` F1, `17.16s`).
- **Correction:** Narrow `docs/BENCHMARK_RESULTS.md:16, 47` to state that `horizon_eval` evaluates oracle-paired `(ref_id, tgt_id)` contigs from the truth MAF in isolation (zero cross-scaffold competition) using micro-averaged base-pair totals without per-scaffold variance, and note that `estimate_hmm_params` retains its `[0.60, 0.72]` clamp and `4000`-bp cap.

---

## F2 -- C3: Unit test count is `46` (not `45`) and binary coverage is `n = 1` smoke test per binary

- **Claim:** `C3` (`docs/STRATEGIC_PLAN.md:3`) -- *"All target pangenome binaries (`horizon_wga`, `horizon_wga --graph`, `horizon_merge`, `horizon_guide`, `horizon_eval`, `horizon_gfa`, `horizon_clean`, `horizon_resolve`, `horizon_lift`, `bench_chr1`) are now implemented in `src/bin/` with unit tests in every binary (`45` passing tests + `1` ignored; `cargo clippy --all-targets -- -D warnings` clean)."*
- **Verdict:** **NARROWED**
- **Prior round context:** Round 2 (L1/L2/L4/L5/L7) established that `b18f25b` had `40` passing tests (`35` in `src/lib.rs` + `1` each in `5` binaries, `0` in `4` binaries). Commit `488b020` added `1` test each to the remaining `4` binaries and `2` tests to `src/lib.rs`.
- **Concrete flaw:** Running `/home/agent/.cargo/bin/cargo test --no-fail-fast` on `HEAD` (`488b020`) executes **`46` passing tests + `1` ignored** (`37 passed; 0 failed; 1 ignored` in `src/lib.rs`, plus `1 passed` in each of the `9` binaries in `src/bin/`), not `45`. Statistically, `37/46` (`80.4%`) of the suite resides in `src/lib.rs`, while every binary in `src/bin/` has sample size `n = 1` synthetic unit test.
- **Correction:** Update `docs/STRATEGIC_PLAN.md:3` from `45 passing tests + 1 ignored` to `46 passing tests (37 library + 1 per binary across all 9 binaries in src/bin/) + 1 ignored`.

---

## F3 -- C4: `FlatSeedIndex` peaks at `28 bytes/entry` during `build`, and `5,863.3 -> 2,409.7 MB` RSS attribution is confounded by simultaneous `488b020` changes while `4` `HashMap`s remain active

- **Claim:** `C4` (`docs/STRATEGIC_PLAN.md:11`) -- *"`ContigIndex` (`src/seeding/mod.rs`) now uses a compact `FlatSeedIndex` (`12 bytes/entry` sorted `Vec<u64>` + `Vec<u32>` with on-the-fly reverse-complement coordinates, cutting `20 Mbp chr6 q-arm` peak RSS from `5,863.3 MB` to `2,409.7 MB` and `chr6` 5 Mbp pericentromere RSS from `3,084.5 MB` to `1,048.2 MB`)."*
- **Verdict:** **NARROWED**
- **Concrete flaws:**
  1. **Peak allocation during `FlatSeedIndex::build` is `28 bytes/entry`, not `12 bytes/entry` (`src/seeding/mod.rs:386, 424-441, 542-546`):** While the steady-state `FlatSeedIndex` stores `keys: Vec<u64>` (`8 B`) + `pos: Vec<u32>` (`4 B`) = `12 bytes/entry`, `FlatSeedIndex::build` first allocates `pairs: Vec<(u64, u32)>` (`16 bytes/entry` with 4-byte struct alignment padding) and then allocates both `keys` (`8 B`) and `pos` (`4 B`) at lines `424-425` **before** `drop(pairs)` at line `441`, yielding a peak transient footprint of **`28 bytes/entry`** (`2.33x` the claimed `12 bytes/entry`). In addition, `ContigIndex` (`lines 542-546`) still stores both `seq` and `seq_rc` (`2 bytes/bp`).
  2. **Confounded before/after RSS attribution (`src/mdl/mod.rs:247-343`, `git show 488b020`):** The `5,863.3 MB -> 2,409.7 MB` (`20 Mbp chr6 q-arm`) and `3,084.5 MB -> 1,048.2 MB` (`5 Mbp chr6 peri`) comparisons are uncontrolled pre/post measurements across `488b020`, which simultaneously: (a) replaced `MASK_W10` `HashMap`s with `FlatSeedIndex`, (b) eliminated redundant `MASK_LASTZ12` `HashMap` rebuilds (`src/mdl/mod.rs:247-293`), (c) changed `mdl_xdrop_extend_with_params` and `bocpd_gc_profile` (altering alignment pair counts), and (d) fixed `rev_offsets` in `canonical_spaced_seeds_impl` (`src/seeding/mod.rs:96`) as well as `FlatSeedIndex::build` (`line 385`).
  3. **`4` `HashMap<u64, Vec<usize>>`s remain active whenever `primary_cov < min_fallback_cov` (`src/mdl/mod.rs:247-343`):** On `20 Mbp chr6 q-arm` (`75.0%` deduplicated coverage `< 85%`), `MASK_LASTZ12` (`iy_lz_fwd`, `iy_lz_rc`) and `MASK_LASTZ11` (`iy_11_fwd`, `iy_11_rc`) still construct `HashMap<u64, Vec<usize>>` indexes, explaining why peak RSS remains `2,409.7 MB` (`5.0x` larger than `40M * 12 B = 480 MB`).
- **Correction:** Qualify `12 bytes/entry` in `docs/STRATEGIC_PLAN.md:11` as steady-state index size (`28 bytes/entry` peak during `FlatSeedIndex::build`, plus `seq` + `seq_rc` in `ContigIndex`), and attribute the `20 Mbp chr6 q-arm` and `chr6` pericentromere RSS reductions to the combined `488b020` changes (noting that `MASK_LASTZ12`/`MASK_LASTZ11` `HashMap` fallbacks remain active when `primary_cov < 0.85`).

---

## F4 -- C5: `5,106,541` pairs still includes `14.1%` many-to-many multiplicity (`718,886` duplicate pairs) with zero accuracy/truth metric or pangenome baseline

- **Claim:** `C5` (`docs/BENCHMARK_RESULTS.md:23, 59`) -- *"Real HPRC `chr22` Pericentromere 5 Mbp (`GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M`): **14.57s**, **1,861.9 MB** peak RSS, **5,106,541** many-to-many `(x, y)` aligned pairs (`unique_x_bp=4,387,655`, `unique_y_bp=4,549,858`), `16,083` Horizon Bubbles."*
- **Verdict:** **NARROWED**
- **Prior round context:** Round 2 L3-F2 (`NARROWED`) flagged that `5,894,191` pairs on `5 Mbp` was a many-to-many count without unique-base coverage. Commit `488b020` added `unique_x_bp` and `unique_y_bp` in `src/bin/horizon_wga.rs:350-360` and re-ran the locus.
- **Concrete flaws:**
  1. **Residual many-to-many multiplicity:** Even with `unique_x_bp = 4,387,655` (`87.8%` of 5 Mbp) and `unique_y_bp = 4,549,858` (`91.0%` of 5 Mbp), `5,106,541` aligned pairs exceeds the `5,000,000`-bp locus length (`102.1%`) and represents `1.164x` `unique_x_bp` (`718,886` duplicate `(x, y)` pairings, `14.1%` of total pairs).
  2. **Zero precision/false-positive rate or head-to-head baseline:** On a repetitive pericentromeric locus, high base coverage (`87.8%` / `91.0%`) without a truth oracle or comparison to `minimap2`/`odgi`/`vg` (still missing despite carried `R2-RESOLVE-LIFT`) cannot distinguish true orthologous alignment from false-positive tandem-repeat or satellite multi-mapping.
  3. **Unexplained `3.1x` RSS/Mbp variance (`n = 1`):** `1,861.9 MB` on `5 Mbp` (`372.4 MB/Mbp`) is `3.1x` higher per Mbp than `20 Mbp chr6 q-arm` (`120.5 MB/Mbp`) and `1.78x` higher than `5 Mbp chr6 peri` (`209.6 MB/Mbp`), reported from a single `n = 1` run with no repeat-run variance.
- **Correction:** Label `C5` in `docs/BENCHMARK_RESULTS.md:23` as an unvalidated `n = 1` coverage/resource measurement (`14.1%` many-to-many pair multiplicity, no orthology precision/truth oracle, and no `minimap2`/`odgi`/`vg` baseline).

---

## F5 -- C6: `horizon_gfa stats` reports an `n = 1` no-op counting pass with no variance or baseline

- **Claim:** `C6` (`docs/BENCHMARK_RESULTS.md:19, 54`) -- *"Full-graph streaming (`horizon_gfa stats`): **88.84s**, **97.1 MB** peak RSS across `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks."*
- **Verdict:** **NARROWED**
- **Concrete flaw:** In `src/bin/horizon_gfa.rs:39-77`, `cmd_stats` is a record-discarding line counter that increments five `u64` counters and reads `/proc/self/status` `VmHWM:`. The `88.84s` wall time and `97.1 MB` peak RSS are `n = 1` single-run point estimates with no repeat-run variance and no head-to-head baseline (`gfatools stat`, `odgi stats`, or `vg stats`), and reflect record-discarding counting rather than graph extraction/analysis (which uses `456.1–1,255.1 MB` in `horizon_gfa extract` at `docs/BENCHMARK_RESULTS.md:55-56`).
- **Correction:** Annotate `docs/BENCHMARK_RESULTS.md:19, 54` that `88.84s` / `97.1 MB` is a single-run (`n = 1`) measurement for record-discarding counting without variance or external tool comparison.

---

## F6 -- C1: No statistical content under the statistics lens

- **Claim:** `C1` (`docs/STRATEGIC_PLAN.md:13`) -- *"`stream_gfa` provides a buffered line-streaming GFA v1.0/v1.1 reader (supporting `.gfa`, multi-member BGZF `.gfa.gz` via `MultiGzDecoder`, and `"-"` `stdin`, with unit tests for empty `L`-line orientations and `.gfa.gz`)."*
- **Verdict:** **UNCHALLENGED**
- **Note:** Checked `src/gfa/mod.rs`; `C1` states a software/API capability and two unit-test cases with no statistical estimator, parameter, or benchmark distribution.

---

## F7 -- Claim-list framing steers away from stale WGA rows, synthetic/no-op graph rows, and missing baselines

- **Observation (`docs/BENCHMARK_RESULTS.md:5, 33-47, 58, 63-65`):** The runner-authored claim list (`C1-C6`) omits:
  1. Five of six pairwise WGA benchmark rows (`multiz100way` Chimp/Mouse, `cactus447way` Horse/Macaque, `multiz470way` Marmoset, `cactus241way chr22_random`, and `GRCh38.chr22 vs CHM13.chr22`), which `docs/BENCHMARK_RESULTS.md:5` admits were recorded prior to `FlatSeedIndex` / `rev_offsets` fixes and remain stale on `488b020`.
  2. The `horizon_clean` (`0` severed, no-op), `horizon_resolve` (`308,760.9` analytical `2*len` bits on a synthetic `200`-seg repeat at `--divergence 0.02`), and `horizon_wga --graph` (`4`-node toy diamond GFA) rows (`docs/BENCHMARK_RESULTS.md:63-65`).
  3. The continued absence of any head-to-head `minimap2`/`odgi`/`vg`/`lastz` baseline across the entire benchmark suite despite carried item `R2-RESOLVE-LIFT`.

---

## Verdict Table

| claim | verdict | why |
|---|---|---|
| C1 | UNCHALLENGED | API/reader capability claim in `src/gfa/mod.rs` with no statistical content |
| C2 | NARROWED | `horizon_eval` (`src/bin/horizon_eval.rs:148-176`) oracle-pairs truth MAF contigs with zero cross-scaffold false positives and no per-scaffold variance; `estimate_hmm_params` (`src/mdl/mod.rs:747-809`) retains its `4000`-bp cap and `[0.60, 0.72]` clamp while `cactus241way` regressed on `488b020` |
| C3 | NARROWED | `cargo test --no-fail-fast` runs `46` passing tests (`37` lib + `9` binary, `n = 1` per binary) + `1` ignored, not `45` |
| C4 | NARROWED | `FlatSeedIndex::build` (`src/seeding/mod.rs:386, 424-441`) peaks at `28 bytes/entry` (not `12`), RSS drop across `488b020` is confounded by simultaneous changes, and `4` `HashMap`s (`MASK_LASTZ12`/`MASK_LASTZ11`) remain active at `primary_cov < 0.85` (`src/mdl/mod.rs:247-343`) |
| C5 | NARROWED | `5,106,541` pairs still carries `14.1%` (`718,886`) many-to-many multiplicity (`1.164x` `unique_x_bp`), with zero precision/truth oracle, zero `minimap2`/`odgi`/`vg` baseline, and `n = 1` `372.4 MB/Mbp` RSS |
| C6 | NARROWED | `horizon_gfa stats` (`src/bin/horizon_gfa.rs:39-77`) is an `n = 1` record-discarding counting pass with no repeat-run variance or external tool baseline |
