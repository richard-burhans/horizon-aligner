# Round 3 Adjudication — Lens L5 (Inference)

**Target:** `horizon-aligner` (`488b020`) — `reports/rounds/round-3/attack-L5-inference.md` (seal digest `c3f8bb931ad0` verified: `c3f8bb931ad00042a46e45089245254b459014ecfdf773c455984afd49bb157e`).
**Prior Round Context:** `scripts/round history` confirms Round 1–2 previously adjudicated `stream_gfa` (27 hits) and `cactus241way` (10 hits, including Round 2 `adjudication-L4.md` `UPHELD` on multi-commit staleness), with 0 prior hits for `FlatSeedIndex`. No `refs/methods/` cards exist in this repository.

---

## F1 — C1: `stream_gfa` tests do not license GFA v1.0/v1.1 interoperability or `stdin` `.gfa.gz` streaming
* **Ruling:** **UPHELD**
* **What survives:** `src/gfa/mod.rs:112-127, 213, 229-329` confirms that `stream_gfa_reader` silently discards `C` (GFA v1.0/v1.1) and `J`/`U`/`O` records, `stdin` (`"-"`) uses plain `BufReader` without `MultiGzDecoder`, and neither `"-"` `stdin` nor `C`/`J` records are covered by unit tests (though `C1`'s own parenthetical only claims `"-"` `stdin` and `.gfa.gz` separately, it is the `"direct GFA v1.0/v1.1 interoperability"` inference in `docs/STRATEGIC_PLAN.md:13` that overreaches; minor note: `J` is GFA v1.2, while `C` is v1.0/v1.1).
* **Repair:** **ENDORSED** — narrowing `docs/STRATEGIC_PLAN.md:13` to a GFA subset (`H`/`S`/`L`/`P`/`W`, ignoring `C`/`J`) with uncompressed-only, untested `stdin` requires no re-measurement; **the fix is free**.

## F2 — C2: Oracle `(ref_id, tgt_id)` pairing in `horizon_eval` and repeated re-runs invalidate precision and held-out inferences
* **Ruling:** **UPHELD**
* **What survives:** `src/bin/horizon_eval.rs:151-162` iterates strictly over `((ref_id, tgt_id), maf_set)` keys present in the truth MAF, so non-homologous scaffold pairs among the 129 Dog scaffolds never enter the precision denominator; furthermore, `git show 488b020 -- docs/BENCHMARK_RESULTS.md` confirms `cactus241way (chr22 1MB)` was re-evaluated across development commits and regressed by `-1.0 pp` F1 (`62.3%` -> `61.3%`) and `+2.05s` (`15.11s` -> `17.16s`) while the other 5 WGA rows remained on pre-`488b020` commits.
* **Repair:** **ENDORSED** — qualifying `docs/BENCHMARK_RESULTS.md:16, 47` as oracle-paired MAF-block precision/recall, removing the `"completely locked away"` held-out claim for `cactus241way (chr22 1MB)`, and noting the `-1.0 pp` F1 / `+2.05s` shift is doc-only; **the fix is free**.

## F3 — C3: Helper-only binary unit tests do not license pipeline correctness or `R2-RESOLVE-LIFT` discharge
* **Ruling:** **NARROWED**
* **What survives:** Checking the whole set of 9 binaries refutes the attacker's blanket `"helper-only"` label—`horizon_clean` (`clean_gfa`, `src/bin/horizon_clean.rs:504-540`), `horizon_lift` (`lift_intervals`, `src/bin/horizon_lift.rs:665-710`), `horizon_resolve` (`pass1_index`/`pass2_rewrite`, `src/bin/horizon_resolve.rs:525-564`), and `horizon_wga` (`load_gfa_dag`, `src/bin/horizon_wga.rs:582`) test core domain pipelines, and `cargo test` passes `46` tests (`37` lib + `9` bin) + `1` ignored (not `45`), while `C3` (`docs/STRATEGIC_PLAN.md:3`) explicitly concedes two-pass tools require seekable paths; what survives is that none of the 9 binary tests invoke `main()`/CLI parsing, `docs/STRATEGIC_PLAN.md:23` (Section 3B) still conflicts with line 3 on `stdin` streaming, and commit `488b020`'s `R2-RESOLVE-LIFT` discharge lacks `minimap2`/`odgi`/`vg` baselines.
* **Repair:** **SPLIT** — qualifying `docs/STRATEGIC_PLAN.md:3, 23` on untested `main()` CLI paths and open `stdin` streaming is **free**; adding `minimap2`/`odgi`/`vg` baselines is **non-free** (requires external benchmark execution).

## F4 — C4: RSS reduction is confounded by bundled `488b020` algorithmic changes, and `10.09 GB` is an unfinished mid-run sample
* **Ruling:** **UPHELD**
* **What survives:** `git show 488b020 --stat` and `docs/BENCHMARK_RESULTS.md:26, 59-62` confirm `FlatSeedIndex` was bundled with `rev_offsets`, `bocpd.rs`, `mdl_xdrop_extend_with_params`, and `horizon_wga.rs` changes without an ablation (shifting `chr22`/`chr6` pericentromere pairs by `-13.4%`/`-11.3%` and `chr6 q-arm` bubbles by `-21.4%`, even though `chr6 q-arm` pairs changed only `+0.02%`), `FlatSeedIndex::build` (`src/seeding/mod.rs:386-443`) holds `pairs: Vec<(u64, u32)>` (`16 B/entry`) alongside `keys` + `pos` (`12 B/entry`) before dropping `pairs`, and the `10.09 GB` check (`docs/BENCHMARK_RESULTS.md:26`) is an unfinished mid-run sample.
* **Repair:** **ENDORSED** — attributing the RSS reduction to the `488b020` build alongside the pair/bubble shifts and removing the unfinished `10.09 GB` mid-run sample as peak-RSS evidence requires no re-measurement; **the fix is free**.

## F5 — C5: Pericentromere coverage/pair counts license no inference of alignment accuracy or structural resolution
* **Ruling:** **UPHELD**
* **What survives:** `docs/BENCHMARK_RESULTS.md:23, 59` reports `5,106,541` many-to-many pairs, `unique_x_bp`/`unique_y_bp`, and `16,083` bubbles on a pre-localized 5 Mbp window with no truth oracle or baseline (`-13.4%` pairs vs `b18f25b` cannot distinguish false-positive reduction from sensitivity loss), and omits the explicit `"no external coordinate truth oracle"` caveat that line 21 includes for `horizon_lift`.
* **Repair:** **ENDORSED** — adding the `"no truth oracle / no baseline"` qualifier to `docs/BENCHMARK_RESULTS.md:23, 59` is doc-only; **the fix is free**.

## F6 — C6: `horizon_gfa stats` (`97.1 MB` RSS) does not license pangenome graph-processing scalability
* **Ruling:** **UPHELD**
* **What survives:** `src/bin/horizon_gfa.rs:39-77` confirms `cmd_stats` only increments five `u64` counters and immediately discards each parsed record, whereas actual graph operations scale to `456.1 MB`–`2,409.7 MB` RSS (`docs/BENCHMARK_RESULTS.md:55-62`) and `horizon_wga --graph` is only tested on a 4-node diamond GFA (`docs/BENCHMARK_RESULTS.md:65`).
* **Repair:** **ENDORSED** — clarifying in `docs/STRATEGIC_PLAN.md:13` and `docs/BENCHMARK_RESULTS.md:19` that `97.1 MB` measures record-discarding line counting rather than graph indexing or alignment memory is doc-only; **the fix is free**.

## F7 — Claim-list steering: runner omitted no-op/synthetic pangenome rows, missing baselines, and stale pre-`488b020` WGA rows
* **Ruling:** **UPHELD**
* **What survives:** `docs/BENCHMARK_RESULTS.md:5, 9, 63-65` and `reports/rounds/round-3/CARRIED.md` confirm that `C1`–`C6` omitted the `horizon_clean` no-op (`0 severed`), the synthetic-only `horizon_resolve` and 4-node `horizon_wga --graph` benchmarks, absent `minimap2`/`odgi`/`vg` baselines (`R2-RESOLVE-LIFT`), and 5 stale pre-`488b020` WGA rows (`R2-HMM-XDROP`).
* **Repair:** **SPLIT** — surfacing no-op/synthetic status and staleness in future claim lists and `docs/BENCHMARK_RESULTS.md` summaries is **free**; re-running the 5 stale WGA rows on a single pinned commit and adding external baselines is **non-free**.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `stream_gfa` ignores `C` (v1.0/v1.1) and `J`/`U`/`O`, `stdin` lacks `MultiGzDecoder` and tests, and the doc narrowing fix is free |
| F2 | UPHELD | `horizon_eval` pairs only oracle truth-MAF contigs and `cactus241way chr22 1MB` was re-run across commits with a `-1.0 pp` F1 regression; the doc fix is free |
| F3 | NARROWED | `horizon_clean`/`horizon_lift`/`horizon_resolve`/`horizon_wga` tests exercise core pipelines (not just helpers), though `main()` CLI parsing, Section 3B `stdin`, and `R2-RESOLVE-LIFT` baselines remain open |
| F4 | UPHELD | `488b020` bundled `FlatSeedIndex` with alignment-changing edits without ablation, `build` transiently allocates `16+12 B/entry`, and `10.09 GB` is an unfinished mid-run sample; the doc fix is free |
| F5 | UPHELD | Pericentromere pair/coverage counts on a pre-localized 5 Mbp window have no truth oracle or baseline; adding the caveat is free |
| F6 | UPHELD | `horizon_gfa stats` (`97.1 MB`) only increments five `u64` counters and discards records, while real graph tasks use `456 MB`–`2.41 GB`; the doc clarification is free |
| F7 | UPHELD | Runner's `C1`–`C6` omitted the `horizon_clean` no-op, synthetic-only `horizon_resolve`/`--graph` rows, missing baselines, and 5 stale WGA rows |
