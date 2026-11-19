# Round 3 Adjudication -- Statistics Lens (L3)

- **Attack digest verified:** `f112216c6692` (`reports/rounds/round-3/attack-L3-statistics.md`).
- **Method cards (`refs/methods/`):** Confirmed absent under both `/home/rcb112/project/gemini/bench/horizon-aligner/refs/methods/` and `/home/rcb112/project/gemini/bench/refs/methods/`; no project method card governs these statistical claims.
- **Prior round history (`scripts/round history --kind adjudication`):** Round 2 L3-F1 (`UPHELD`) already documented `estimate_hmm_params` truncation (`4000` bp), `identity_prior` clamping (`[0.60, 0.72]`), and `dx != dy` indel-event counting; Round 2 L1/L2/L4/L5/L7 (`UPHELD`) documented the `40`-test count on `b18f25b`; Round 2 L3-F2 (`NARROWED`) documented many-to-many pair counts without unique-base coverage on `chr22` pericentromere. The attacker appropriately disclosed this prior context and scoped Round 3 findings to what changed or persisted in `488b020`.

---

## F1 -- C2: `horizon_eval` oracle-pairs target scaffolds with zero cross-scaffold false positives, and `estimate_hmm_params` remains truncated/clamped while `cactus241way` regressed

- **Diagnosis: UPHELD.**
  - **What survives:** Verified in `src/bin/horizon_eval.rs:148-176` that `horizon_eval` iterates exclusively over `(ref_id, tgt_id)` pairs present in the truth MAF (`maf_by_pair`), aligning each pair in isolation without cross-scaffold competition, and computes micro-averaged base-pair totals (`total_overlap / total_maf`, `total_overlap / total_proto`) with no per-scaffold dispersion (median, IQR, or CI). Verified in `src/mdl/mod.rs:747-809` that `estimate_hmm_params` retains its `take(64)` / `4000`-bp cap, `[0.60, 0.72]` `identity_prior` clamp, and `lcb.windows(2).take(64)` indel-event heuristic. Verified in `git show 488b020 -- docs/BENCHMARK_RESULTS.md` that `cactus241way chr22 1MB` shifted from `63.6%` recall / `61.0%` precision (`62.3%` F1, `15.11s`) to `62.3%` recall / `60.3%` precision (`61.3%` F1, `17.16s`).
- **Repair: UPHELD.**
  - Narrowing `docs/BENCHMARK_RESULTS.md:16, 47` to state that `horizon_eval` evaluates oracle-paired `(ref_id, tgt_id)` contigs from the truth MAF in isolation using micro-averaged base-pair totals without per-scaffold variance, and noting the retained `estimate_hmm_params` `[0.60, 0.72]` clamp and `4000`-bp cap, accurately qualifies C2; **the fix is free** (documentation-only edit requiring no re-measurement).

---

## F2 -- C3: Unit test count is `46` (not `45`) and binary coverage is `n = 1` smoke test per binary

- **Diagnosis: UPHELD.**
  - **What survives:** Running `/home/agent/.cargo/bin/cargo test --no-fail-fast` on `HEAD` (`488b020`) executes **`46` passing tests + `1` ignored** (`37 passed; 0 failed; 1 ignored` in `src/lib.rs`, plus `1 passed` in each of the `9` binaries in `src/bin/`), so `45` undercounts by `1` (`37/46 = 80.4%` in `src/lib.rs`, and `n = 1` unit test per binary in `src/bin/`). Note that `horizon_wga --graph` is a CLI mode of `horizon_wga` rather than a 10th binary file.
- **Repair: UPHELD.**
  - Updating `docs/STRATEGIC_PLAN.md:3` from `45 passing tests + 1 ignored` to `46 passing tests (37 library + 1 per binary across all 9 binaries in src/bin/) + 1 ignored` is exact; **the fix is free**.

---

## F3 -- C4: `FlatSeedIndex` peaks at `28 bytes/entry` during `build`, and `5,863.3 -> 2,409.7 MB` RSS attribution is confounded by simultaneous `488b020` changes while `4` `HashMap`s remain active

- **Diagnosis: UPHELD.**
  - **What survives:** Verified in `src/seeding/mod.rs:386, 424-441, 542-546` that `FlatSeedIndex::build` allocates `pairs: Vec<(u64, u32)>` (`16 bytes/entry` with alignment padding) alongside `keys: Vec<u64>` (`8 bytes/entry`) and `pos: Vec<u32>` (`4 bytes/entry`) before `drop(pairs)` at line `441`, giving a `28 bytes/entry` transient build peak (`2.33x` the `12 bytes/entry` steady-state index), while `ContigIndex` also retains both `seq` and `seq_rc` (`2 bytes/bp`). Verified in `git show 488b020` and `src/mdl/mod.rs:247-343` that the `5,863.3 -> 2,409.7 MB` (`20 Mbp chr6 q-arm`) and `3,084.5 -> 1,048.2 MB` (`5 Mbp chr6 peri`) RSS drops coincide with multiple simultaneous changes (eliminating redundant `MASK_LASTZ12` rebuilds, `mdl_xdrop_extend_with_params`, `bocpd_gc_profile`, and `rev_offsets` fixes), and that `4` `HashMap<u64, Vec<usize>>` fallback indexes (`iy_lz_fwd`, `iy_lz_rc`, `iy_11_fwd`, `iy_11_rc`) remain active whenever `primary_cov < min_fallback_cov` (`0.85`), as on `20 Mbp chr6 q-arm` (`75.0%` deduplicated coverage).
- **Repair: UPHELD.**
  - Qualifying `12 bytes/entry` in `docs/STRATEGIC_PLAN.md:11` as steady-state index size (`28 bytes/entry` transient peak during `FlatSeedIndex::build`, plus `seq`/`seq_rc` in `ContigIndex`) and attributing the RSS reductions to the combined `488b020` changes (noting active `MASK_LASTZ12`/`MASK_LASTZ11` `HashMap` fallbacks when `primary_cov < 0.85`) is accurate; **the fix is free**.

---

## F4 -- C5: `5,106,541` pairs still includes `14.1%` many-to-many multiplicity (`718,886` duplicate pairs) with zero accuracy/truth metric or pangenome baseline

- **Diagnosis: NARROWED.**
  - **What survives:** Verified in `docs/BENCHMARK_RESULTS.md:23, 59` that `5,106,541` pairs on the `5 Mbp` `chr22` pericentromere exceeds `unique_x_bp = 4,387,655` (`87.8%`) by `718,886` pairs (`14.1%` of total pairs, `1.164x` `unique_x_bp`) and `unique_y_bp = 4,549,858` (`91.0%`) by `556,683` pairs, with no orthology precision/truth oracle, no `minimap2`/`odgi`/`vg` baseline, and an unexplained `n = 1` peak RSS of `1,861.9 MB` (`372.4 MB/Mbp`, `3.1x` `chr6 q-arm` per-Mbp RSS). **Narrowing:** `718,886` is excess many-to-many multiplicity relative to `unique_x_bp` (one-to-many/many-to-many mappings), not duplicate identical `(x, y)` tuples (as `aligned_pairs` is deduplicated in `horizon_wga`).
- **Repair: UPHELD.**
  - Labeling `C5` in `docs/BENCHMARK_RESULTS.md:23` as an unvalidated `n = 1` coverage/resource measurement (`14.1%` many-to-many pair multiplicity over `unique_x_bp`, no orthology precision/truth oracle, and no `minimap2`/`odgi`/`vg` baseline) is accurate; **the fix is free**.

---

## F5 -- C6: `horizon_gfa stats` reports an `n = 1` no-op counting pass with no variance or baseline

- **Diagnosis: UPHELD.**
  - **What survives:** Verified in `src/bin/horizon_gfa.rs:39-77` that `cmd_stats` discards every parsed GFA record after incrementing five `u64` counters and reading `/proc/self/status` `VmHWM:`, so `88.84s` and `97.1 MB` (`docs/BENCHMARK_RESULTS.md:19, 54`) are `n = 1` single-run estimates for record-discarding counting with no repeat-run variance and no `gfatools`/`odgi`/`vg` baseline (whereas `horizon_gfa extract` uses `456.1–1,255.1 MB` at `docs/BENCHMARK_RESULTS.md:55-56`).
- **Repair: UPHELD.**
  - Annotating `docs/BENCHMARK_RESULTS.md:19, 54` that `88.84s` / `97.1 MB` is a single-run (`n = 1`) measurement for record-discarding counting without variance or external tool comparison is accurate; **the fix is free**.

---

## F6 -- C1: No statistical content under the statistics lens

- **Diagnosis: UPHELD.**
  - **What survives:** Verified in `src/gfa/mod.rs` and `docs/STRATEGIC_PLAN.md:13` that `C1` describes `stream_gfa`'s streaming API and two unit-test cases with no statistical estimator, parameter, or benchmark distribution (`UNCHALLENGED` under L3).

---

## F7 -- Claim-list framing steers away from stale WGA rows, synthetic/no-op graph rows, and missing baselines

- **Diagnosis: UPHELD.**
  - **What survives:** Verified in `docs/BENCHMARK_RESULTS.md:5, 33-47, 58, 63-65` that five of six pairwise WGA rows plus `GRCh38.chr22 vs CHM13.chr22` predate `FlatSeedIndex`/`rev_offsets` fixes and remain stale on `488b020`, that `horizon_clean` (`0` severed, no-op), `horizon_resolve` (`308,760.9` analytical `2*len` bits on synthetic `200`-seg repeat), and `horizon_wga --graph` (`4`-node toy diamond GFA) were omitted from `C1-C6`, and that no `minimap2`/`odgi`/`vg`/`lastz` baseline exists anywhere in `docs/BENCHMARK_RESULTS.md`.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `horizon_eval` (`src/bin/horizon_eval.rs:148-176`) oracle-pairs truth MAF contigs with zero cross-scaffold competition and unweighted micro-averaging; `estimate_hmm_params` (`src/mdl/mod.rs:747-809`) retains its `4000`-bp cap and `[0.60, 0.72]` clamp while `cactus241way` regressed on `488b020`; doc fix is free |
| F2 | UPHELD | `cargo test --no-fail-fast` runs `46` passing tests (`37` lib + `9` binary, `n = 1` per binary) + `1` ignored, not `45`; doc fix is free |
| F3 | UPHELD | `FlatSeedIndex::build` (`src/seeding/mod.rs:386, 424-441`) peaks at `28 bytes/entry`, `488b020` RSS attribution is confounded by simultaneous changes, and `4` `HashMap`s (`MASK_LASTZ12`/`MASK_LASTZ11`) remain active at `primary_cov < 0.85` (`src/mdl/mod.rs:247-343`); doc fix is free |
| F4 | NARROWED | `5,106,541` pairs has `14.1%` (`718,886`) excess many-to-many multiplicity over `unique_x_bp` (not duplicate identical `(x, y)` tuples), with zero truth/precision oracle, zero baseline, and `n = 1` `372.4 MB/Mbp` RSS; doc fix is free |
| F5 | UPHELD | `horizon_gfa stats` (`src/bin/horizon_gfa.rs:39-77`) is an `n = 1` record-discarding counting pass with no variance or external tool baseline; doc fix is free |
| F6 | UPHELD | `C1` (`src/gfa/mod.rs`, `docs/STRATEGIC_PLAN.md:13`) states an API/reader capability with no statistical content (`UNCHALLENGED` under L3) |
| F7 | UPHELD | `C1-C6` omit stale pre-`488b020` WGA rows, synthetic/no-op graph rows (`horizon_clean`, `horizon_resolve`, `horizon_wga --graph`), and the absence of any `minimap2`/`odgi`/`vg`/`lastz` baseline |
