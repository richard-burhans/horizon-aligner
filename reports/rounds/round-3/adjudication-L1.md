# Round 3 Adjudication -- Correctness Lens (L1)

**Attack digest verified:** `71d2ec66776b` (`71d2ec66776b67ca40c0f189a75736c3ec01ac2edeeba6f78e4501695d368dd0` on `reports/rounds/round-3/attack-L1-correctness.md`).

**Method cards & prior-round check:** No method cards exist in `refs/methods/`. Prior adjudications (`./scripts/round history` for `stream_gfa`, `FlatSeedIndex`, `cactus241way`) show `stream_gfa` (`R2` `F8` in `reports/rounds/round-2/adjudication-L1.md`) was previously adjudicated and endorsed for a free doc narrowing that was not applied; `FlatSeedIndex` has no prior adjudication.

---

## F1 -- C1 (`docs/STRATEGIC_PLAN.md:13`, `stream_gfa`): **UPHELD**

* **What survives:** `stream_gfa_reader` (`src/gfa/mod.rs:148-216`) silently drops GFA v1.0 `C` and GFA v1.1 `J`/`U`/`O` records (`_ => {}` at line 213) and coerces missing or non-numeric `W`-line `start`/`end` fields (including valid GFA v1.1 `*` wildcards) to `0` (`parse().unwrap_or(0)` at lines 208-209), `"-"` `stdin` (`lines 112-116`) does not decompress gzipped input, and Round 2's endorsed doc repair (`adjudication-L1.md` `F8`) to `docs/STRATEGIC_PLAN.md:13` was never applied. (Minor note: `C1` claimed unit tests for empty `L`-line orientations and `.gfa.gz`, both of which exist in `src/gfa/mod.rs:229-330`, not for `stdin`.)
* **Repair:** Endorsed on `docs/STRATEGIC_PLAN.md:13` to state that `stream_gfa` parses `S`/`L`/`P`/`W` records (ignoring `C`/`J`/`U`/`O`), coerces non-numeric `W` coordinates to `0`, and reads uncompressed `"-"` `stdin`. **The fix is free** (text-only edit, no re-measurement).

---

## F2 -- C2 (`docs/BENCHMARK_RESULTS.md:16`, `cactus241way chr22 1MB`): **NARROWED**

* **What survives:** The attacker's diagnosis is verified on `/tmp/dog_zoo_1m.fa` (`31` FASTA headers, not `129`), `/tmp/pw_zoo_1m.maf` (`31` unique `Canis_lupus_familiaris` scaffolds), and `src/bin/horizon_eval.rs:151-162` (which iterates only over `(ref_id, tgt_id)` pairs present in the truth MAF); however, re-running `horizon_eval` reproduces **62.3%** recall, **60.3%** precision, and **61.3% F1** (`elapsed=14.310s proto_pairs=136348 maf_pairs=132089 overlap=82228 bubbles=3490`), so the `REFUTED` label overstates the scope.
* **Repair:** Endorsed on `docs/BENCHMARK_RESULTS.md:16` to replace `across 129 target scaffolds` with `across 31 truth-paired target scaffolds` (noting wall-clock variance `14.3s`–`22.0s`). **The fix is free** (text-only edit, no re-measurement).

---

## F3 -- C3 (`docs/STRATEGIC_PLAN.md:3`, binaries and test count): **UPHELD**

* **What survives:** `cargo test` runs **46** passing tests + `1` ignored (`37` in `src/lib.rs` + `1` in each of the `9` binaries in `src/bin/`), not `45` passing tests (`cargo clippy --all-targets -- -D warnings` exits `0`), and `horizon_clean` (`src/bin/horizon_clean.rs:156, 224-234, 265-274`) is a no-op on `*`-sequence GFAs and emits unlinked `S` bubble records (`HB:i:1`) when it does sever a link.
* **Repair:** Endorsed on `docs/STRATEGIC_PLAN.md:3` to update `45 passing tests` to `46 passing tests` and qualify `horizon_clean` on `*`-sequence GFAs. **The fix is free** (text-only edit, no re-measurement).

---

## F4 -- C4 (`docs/STRATEGIC_PLAN.md:11`, `FlatSeedIndex`): **UPHELD**

* **What survives:** `FlatSeedIndex` (`src/seeding/mod.rs:332-337, 525-532`) achieves `12 bytes/entry` post-build storage and on-the-fly reverse-complement coordinates (`other.seq_len - other.span - y_pos`), but `FlatSeedIndex::build` (`lines 386-441`) peaks at up to `28 bytes/entry` before `drop(pairs)`, `ContigIndex` (`lines 542-546, 570-574`) still materializes `seq_rc`, and `align_contig_pair_with_config` (`src/mdl/mod.rs:210-261`) unconditionally builds two `MASK_W10` `FlatSeedIndex`es per call and falls back to `HashMap<u64, Vec<usize>>` (`MASK_LASTZ12` on `seq_y` and `seq_y_rc`) when `primary_cov < min_fallback_cov`.
* **Repair:** Endorsed on `docs/STRATEGIC_PLAN.md:11` to qualify `12 bytes/entry` as post-build storage (`28 bytes/entry` peak during `FlatSeedIndex::build`) and note `ContigIndex.seq_rc` and the `MASK_W10`/`MASK_LASTZ12` paths. **The fix is free** (text-only edit, no re-measurement).

---

## F5 -- C5 (`docs/BENCHMARK_RESULTS.md:23`, `chr22` Pericentromere 5 Mbp): **UPHELD**

* **What survives:** Re-running `./target/release/horizon_wga /tmp/grch38_chr22_peri_5m.fa /tmp/chm13_chr22_peri_5m.fa --resolve-dups` yields `contigs_x=1 contigs_y=1 pruned_pairs=0 aligned_pairs=5106541 unique_x_bp=4387655 unique_y_bp=4549858 total_bubbles=16083 elapsed_s=13.6094`, confirming the attacker's `SURVIVED` ruling.

---

## F6 -- C6 (`docs/BENCHMARK_RESULTS.md:19`, `horizon_gfa stats`): **UPHELD**

* **What survives:** `hprc-v1.0-mc-grch38.gfa.gz` is absent on disk in the workspace and `/tmp/`, so the counts remain `UNCHALLENGED`, while `src/bin/horizon_gfa.rs:28-76` confirms streaming without retaining records in RAM and `VmHWM:` KiB-to-MiB reporting.

---

## F7 -- Framing bias (`C1`–`C6` selection): **UPHELD**

* **What survives:** The claim list selects only the single re-run WGA row (`C2`) while omitting the 5 stale WGA rows acknowledged in `docs/BENCHMARK_RESULTS.md:5` (`R2-HMM-XDROP`), the missing `minimap2`/`odgi`/`vg` pangenome baselines (`R2-RESOLVE-LIFT`), the `horizon_clean` no-op (`docs/BENCHMARK_RESULTS.md:63`), and the unapplied Round 2 repair to `docs/STRATEGIC_PLAN.md:13` (`C1`).

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `stream_gfa` drops `C`/`J`/`U`/`O`, coerces invalid `W` coords to `0`, does not gunzip `stdin`, and Round 2's endorsed doc fix was never applied; doc fix is free |
| F2 | NARROWED | `129 target scaffolds` is false (`31` in `/tmp/dog_zoo_1m.fa` and `/tmp/pw_zoo_1m.maf`) and `horizon_eval` filters by truth MAF, but `62.3%` recall / `60.3%` precision / `61.3%` F1 reproduce; doc fix is free |
| F3 | UPHELD | `cargo test` runs `46` passing tests + `1` ignored (not `45`), and `horizon_clean` is a no-op on `*`-seq GFAs; doc fix is free |
| F4 | UPHELD | `12 bytes/entry` is post-build (`28 bytes/entry` peak in `build`), `ContigIndex` retains `seq_rc`, and `align_contig_pair_with_config` builds `MASK_W10` and `MASK_LASTZ12` `HashMap` fallbacks; doc fix is free |
| F5 | UPHELD | `horizon_wga --resolve-dups` reproduced `5,106,541` pairs, `unique_x_bp=4,387,655`, `unique_y_bp=4,549,858`, and `16,083` bubbles |
| F6 | UPHELD | `hprc-v1.0-mc-grch38.gfa.gz` is absent on disk (`UNCHALLENGED`), and `cmd_stats` streaming/`VmHWM` code has no correctness flaw |
| F7 | UPHELD | Claim list omits 5 stale WGA rows (`R2-HMM-XDROP`), missing `minimap2`/`odgi`/`vg` baselines (`R2-RESOLVE-LIFT`), `horizon_clean` no-op, and Round 2's unapplied `C1` doc fix |
