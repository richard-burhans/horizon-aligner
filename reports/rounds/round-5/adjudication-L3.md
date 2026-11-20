# Round 5 Adjudication -- Statistics Lens

**Target attack:** [`reports/rounds/round-5/attack-statistics.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-5/attack-statistics.md) (sealed SHA-256 `aba1a127eb078890912cb9d038f685d95206058af53f9760dbafdad2e6f24bf8`, verified).
**Target surfaces:** [`manuscript/horizon_paper.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md), [`manuscript/horizon_companion.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_companion.md), [`docs/BENCHMARK_RESULTS.md`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md), [`manuscript/data/chr6_hprc_comparison.json`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/data/chr6_hprc_comparison.json), [`manuscript/data/chr6_minimap2_summary.json`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/data/chr6_minimap2_summary.json), `src/`.
**Method cards & prior rounds (Rules 5–6):** This repository has no `refs/methods/` directory, so no method card governs these disputes. Checked `scripts/round history` for `n_cal`, `chr6`, and `horizon_eval`; prior rulings `r3 adjudication-L5.md:13 [F2]`, `r4 adjudication-L2.md:15 [F1]`, and `r4 adjudication-L3.md:62 [F7]` are properly disclosed by the attack.

---

## F1 -- C1: Single-run (`n = 1`) runtime comparison with unquantified variance and confounded task scope

- **Ruling:** **UPHELD** (attacker's `NARROWED` verdict on `C1` stands).
- **What survives:** `C1` compares a single (`n = 1`) run of a counter-only line scanner (`horizon_gfa stats`, `0.08 s`, `22.6 MB` RSS) against graph-constructing `vg` commands (`vg stats -z -l -N` at `4.2984 s`, `vg convert -g -p` at `1.8292 s`) with no replicate variance/CI, no cold/warm cache control, and asymmetric decimal precision.
- **Verification:**
  1. [`manuscript/horizon_paper.md:82`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L82) and [`docs/BENCHMARK_RESULTS.md:25,76`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L25) report single point estimates with no replicate count or dispersion metric and mixed precision (`0.08 s` vs. `4.2984 s` / `1.8292 s`).
  2. [`src/bin/horizon_gfa.rs:20-45`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_gfa.rs#L20-L45) confirms `cmd_stats` only increments five `u64` counters (`s_count`, `s_bases`, `l_count`, `p_count`, `w_count`) in a single streaming pass and discards every record, whereas `vg convert -g -p` and `vg stats -z -l -N` build/serialize an in-memory/packed graph representation.
- **Repair:** Qualifying [`manuscript/horizon_paper.md:82`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L82) as a single-run (`n = 1`) measurement comparing counter-only streaming line-counting (`horizon_gfa stats`) against full graph construction/conversion (`vg stats` / `vg convert`), and matching decimal precision, is endorsed (**the fix is free**). Running `n >= 5` cold/warm replicates is non-free.

---

## F2 -- C2: Conditioned precision denominator (`48.63%`) omits `60.4M` (`21.4%`) emitted pairs, inconsistent with `~2.6x` multiplicity and `n = 1` locus

- **Ruling:** **NARROWED** (**UPHELD** on the undisclosed conditioned denominator, multiplicity mismatch, and `n = 1` qualification; **OVERTURNED** on the attacker's `REFUTED` verdict for `C2` as a whole, which is `NARROWED` rather than `REFUTED`).
- **What survives:** `48.63%` is `exact_coord_pairs / pairs_on_hprc_1to1_x` (`107,836,669 / 221,734,653`), silently excluding `60,428,393` (`21.42%`) emitted pairs outside HPRC `1-to-1` segments (across all `282,163,046` pairs, exact pair precision against HPRC `1-to-1` segments is `38.22%`), while `~2.6x` multiplicity (`282,163,046 / 109,575,507 = 2.575x`) uses all emitted pairs rather than the conditioned `221.73M` subset (`2.014x`).
- **Verification:**
  1. [`manuscript/data/chr6_hprc_comparison.json:10,13,15,16,28`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/data/chr6_hprc_comparison.json#L10-L28) confirms `107,836,669 / 221,734,653 = 0.486332` (`exact_coord_precision_on_shared`), `282,163,046 - 221,734,653 = 60,428,393` (`21.416%`), and `107,836,669 / 282,163,046 = 0.382179` (`38.22%`). [`manuscript/data/chr6_minimap2_summary.json:3,6,7,15`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/data/chr6_minimap2_summary.json#L3-L15) confirms `minimap2` differs by only `0.15` pp (`99.75%` conditioned vs. `99.60%` unconditional).
  2. Multiplicity arithmetic is verified (`282,163,046 / 109,575,507 = 2.575x` vs. `221,734,653 / 110,082,990 = 2.014x`), and `n = 1` chromosome/q-arm/run was already ruled `UPHELD` in `r4 adjudication-L2.md:15 [F1]`.
  3. **Why `C2` is `NARROWED`, not `REFUTED`:** `C2`'s recall claims (`99.54%` base-coverage recall, `97.96%` exact coordinate recall, `98.14%` within `1 kb`) and the `minimap2` metrics (`99.71%` exact recall, `99.75%` exact precision) are verified in the JSON files. Moreover, HPRC `1-to-1` segments (`110.08 Mbp`) are a `1-to-1` orthology reference, not a catalog of true paralogous/many-to-many pairs, so `38.22%` measures exact precision *against `1-to-1` HPRC segments* rather than total false-positive rate.
- **Repair:** Endorse replacing `"48.63% many-to-many exact pair precision (~2.6× repeat/paralog multiplicity)"` in [`manuscript/horizon_paper.md:76`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L76) and Figure 7B with `"38.22% unconditional exact pair precision against HPRC 1-to-1 segments across all 282.16M emitted pairs (48.63% conditioned on the 221.73M pairs within HPRC 1-to-1 GRCh38.chr6 segments, 2.01× local multiplicity; 2.58× total emitted multiplicity)"`, qualified as `n = 1` chromosome/q-arm/run (**the fix is free**).

---

## F3 -- C3: Single-run (`n = 1`) multi-species timings with zero accuracy or concordance statistics

- **Ruling:** **UPHELD** (attacker's `NARROWED` verdict on `C3` stands).
- **What survives:** `C3` reports `n = 1` wall-clock timings and block/segment counts (`horizon_merge`: `21.78 s`, `26.6 MB`, `460` blocks; `horizon_guide`: `0.859 s`, `1,109` blocks, `4` HAL bottom segments) with no replicate variance and zero alignment-accuracy or column-concordance metrics against reference multi-species MAFs.
- **Verification:** [`manuscript/horizon_paper.md:10,82`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L82) and [`docs/BENCHMARK_RESULTS.md:10-11,57-58`](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/BENCHMARK_RESULTS.md#L10-L11) contain only single-run runtime, RSS, and output-volume counts.
- **Repair:** Qualifying [`manuscript/horizon_paper.md:82`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L82) as `n = 1` throughput/memory measurements without alignment-accuracy validation against reference multi-species MAFs is endorsed (**the fix is free**).

---

## F4 -- C4: Oracle `(ref_id, tgt_id)` pairing and MAF-only FASTA extraction inflate precision; `18%`-GC null claim rests on `n = 1` synthetic draw

- **Ruling:** **UPHELD** (attacker's `NARROWED` verdict on `C4` stands).
- **What survives:** `horizon_eval` evaluates only oracle truth-MAF `(ref_id, tgt_id)` contig pairs on MAF-extracted homologous segments stripped of unaligned background, and `"eliminating all false-positive pairs on 18%-GC random DNA"` rests on `n = 1` `10,000`-bp synthetic draw (95% Clopper–Pearson upper bound `95%`), alongside `288` and `149` pairs on single `30%`-GC and `45%`-GC draws.
- **Verification:**
  1. [`src/bin/horizon_eval.rs:165-197`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon_eval.rs#L165-L197) iterates strictly over `(ref_id, tgt_id)` keys in `maf_by_pair` (already upheld in `r3 adjudication-L5.md:13 [F2]`), and [`scripts/benchmarks/extract_fasta_from_maf.py:40-45`](file:///home/rcb112/project/gemini/bench/horizon-aligner/scripts/benchmarks/extract_fasta_from_maf.py#L40-L45) concatenates only MAF-aligned blocks.
  2. [`manuscript/horizon_paper.md:10,61-70`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L10) omits both qualifications in the Abstract and Table 1, and [`manuscript/horizon_paper.md:59`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L59) generalizes a single `10,000`-bp `18%`-GC draw as eliminating all false-positive pairs.
- **Repair:** Qualifying [`manuscript/horizon_paper.md:10`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_paper.md#L10) and Table 1 as conditional precision on MAF-extracted homologous segments under oracle `(ref_id, tgt_id)` pairing, and replacing `"eliminating all false-positive pairs on 18%-GC random DNA"` with `"emitting 0 pairs on a single 10,000-bp 18%-GC random DNA draw (vs. 288 and 149 pairs on single 30%-GC and 45%-GC draws)"` is endorsed (**the fix is free**).

---

## F5 -- C5 and C6: Deterministic indexing and tag-emission claims

- **Ruling:** **UPHELD** (`UNCHALLENGED` assessment for both `C5` and `C6` stands).
- **What survives:** Neither [`manuscript/horizon_companion.md:139-143`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_companion.md#L139-L143) (`BandGrid` index formula) nor [`manuscript/horizon_companion.md:52`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_companion.md#L52) (`# EG:Z:` tag emission) makes a statistical claim.
- **Repair:** None needed.

---

## F6 -- Claim-list steering: `C1`–`C6` omit uncalibrated `n_cal = 15` CRC, `0.55` Karlin coefficient, and companion toy-example discrepancy

- **Ruling:** **UPHELD**.
- **What survives:**
  1. `C1`–`C6` omit `n_cal = 15` (`α = 0.07`, `λ̂ = 10.0` bits) CRC and `KARLIN_MDL_COEFF = 0.55` ([`src/mdl/mod.rs:26,71`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/mdl/mod.rs#L26)), already adjudicated `NARROWED` in `r4 adjudication-L3.md:62 [F7]` and still lacking a calibration dataset/script in the repository.
  2. In [`manuscript/horizon_companion.md:173-181`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_companion.md#L173-L181), the worked `22 bp` toy example prints `P(W) ≈ 0.72` and `ΔL ≈ −6.4` bits for the `6 bp` zero-match gap (`GGGGCC` vs. `CCTTTT`), whereas direct evaluation of `bocpd_gc_profile_into` ([`src/bocpd.rs:87-175`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bocpd.rs#L87-L175)) and `joint_pair_bits` ([`src/bocpd.rs:41-52`](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bocpd.rs#L41-L52)) under `identity = 0.65`, `gap_open_bits = -4.5` yields `P(W) = 0.687` at `t = 11` (mean `0.699` over `0..11`) and `6 * log2(1 - 0.65) + 5 * mm_bits = -9.755` bits (not `-6.4` bits).
- **Repair:** Recomputing the worked toy example in [`manuscript/horizon_companion.md:177-181`](file:///home/rcb112/project/gemini/bench/horizon-aligner/manuscript/horizon_companion.md#L177-L181) from the actual Rust `bocpd_gc_profile` / `banded_mdl_gap` output and including CRC calibration provenance in the claim scope is endorsed (**the fix is free**).

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `C1` is `n = 1` with no variance, asymmetric precision, and counter-only `horizon_gfa stats` vs. graph-building `vg`; free text qualification endorsed |
| F2 | NARROWED | `48.63%` omits `60.43M` (`21.4%`) pairs (`38.22%` unconditional vs. `1-to-1`) and mismatches `~2.6x`, but `C2`'s recall and `minimap2` metrics hold (`C2` is `NARROWED`, not `REFUTED`); free text fix endorsed |
| F3 | UPHELD | `C3` reports `n = 1` timings and block counts with zero alignment-accuracy metrics; free text qualification endorsed |
| F4 | UPHELD | `horizon_eval` uses oracle `(ref_id, tgt_id)` pairing on MAF-only FASTAs, and `18%`-GC null rests on `n = 1` draw; free text qualification endorsed |
| F5 | UPHELD | `C5` (`BandGrid` indexing) and `C6` (`# EG:Z:` tags) are deterministic, non-statistical claims (`UNCHALLENGED`) |
| F6 | UPHELD | `C1`–`C6` omit `n_cal = 15` CRC / `0.55` (`r4 [F7]`), and `horizon_companion.md:177-181` prints `-6.4` bits / `0.72` vs. `-9.755` bits / `0.687` from `src/bocpd.rs`; free text fix endorsed |
