# Round 3 Attack -- Correctness (L1)

Prior round history checked via `./scripts/round history` for `FlatSeedIndex`, `stream_gfa`, `cactus241way`, and `horizon_clean` (`--kind adjudication`).

## F1 -- C1 (`stream_gfa`): Round 2's endorsed narrowing to `S`/`L`/`P`/`W` was never applied
- **Claim:** `C1` (`docs/STRATEGIC_PLAN.md:13`)
- **Verdict:** **NARROWED**
- **Prior adjudication:** Round 2 `adjudication-L1.md` (`F8`) **UPHELD** narrowing `C1` because `stream_gfa_reader` (`src/gfa/mod.rs:148-216`) silently drops GFA v1.0 `C` and GFA v1.1 `J` (plus `U`/`O`) records (`_ => {}` at `src/gfa/mod.rs:213`) and coerces missing or unparseable `W`-line `start`/`end` coordinates (including valid GFA v1.1 `*` wildcards) to `0` (`parse().unwrap_or(0)` at `src/gfa/mod.rs:208-209`). Round 2's adjudicator endorsed repairing `docs/STRATEGIC_PLAN.md:13` to state that `stream_gfa` parses `S`/`L`/`P`/`W` records (ignoring `C`/`J`/`U`/`O`).
- **What changed:** Nothing changed in `docs/STRATEGIC_PLAN.md:13` or `src/gfa/mod.rs:202-214`. Additionally, `"-"` `stdin` has no unit test in `src/gfa/mod.rs:223-330` and does not decompress gzipped `stdin`.
- **Correction:** Update `docs/STRATEGIC_PLAN.md:13` to state that `stream_gfa` parses GFA v1.0/v1.1 `S`/`L`/`P`/`W` records (ignoring `C`/`J`/`U`/`O`) and coerces non-numeric `W` coordinates to `0`.

## F2 -- C2 (`cactus241way chr22 1MB`): target scaffold count is `31`, not `129`, and `horizon_eval` uses truth-MAF pair filtering
- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:16`)
- **Verdict:** **REFUTED**
- **Evidence:**
  - Re-running `/home/rcb112/project/gemini/bench/horizon-aligner/target/release/horizon_eval hg38 Canis_lupus_familiaris /tmp/hg38_zoo_1m.fa /tmp/dog_zoo_1m.fa /tmp/pw_zoo_1m.maf` outputs `elapsed=21.987s proto_pairs=136348 maf_pairs=132089 overlap=82228 raw_maf_recall=62.3% raw_maf_precision=60.3% bubbles=3490` (`61.3% F1`).
  - However, the claim `"across 129 target scaffolds"` is factually false: `/tmp/dog_zoo_1m.fa` contains **31** target scaffolds (`grep -c "^>" /tmp/dog_zoo_1m.fa` = `31`), and `/tmp/pw_zoo_1m.maf` contains **31** distinct `Canis_lupus_familiaris` scaffolds (`grep "^s Canis_lupus_familiaris" /tmp/pw_zoo_1m.maf | awk '{print $2}' | sort -u | wc -l` = `31`).
  - Moreover, `src/bin/horizon_eval.rs:151-162` iterates only over `(ref_id, tgt_id)` pairs already present in `maf_by_pair` (the truth MAF), rather than discovering homologous scaffolds across an unannotated assembly.
- **Correction:** Replace `across 129 target scaffolds` in `docs/BENCHMARK_RESULTS.md:16` with `across 31 truth-paired target scaffolds` (and note the `21.99s` / `17.16s` run-to-run variance).

## F3 -- C3 (pangenome binaries and test count): `cargo test` runs `46` passing tests, not `45`, and `horizon_clean` is a no-op on `*`-seq GFAs
- **Claim:** `C3` (`docs/STRATEGIC_PLAN.md:3`)
- **Verdict:** **NARROWED**
- **Evidence:**
  - Running `cargo test` in `/home/rcb112/project/gemini/bench/horizon-aligner` executes **46** passing tests + `1` ignored (`37` in `src/lib.rs` + `1` in each of the `9` binaries in `src/bin/`), not `45` passing tests. (`cargo clippy --all-targets -- -D warnings` exits `0`.)
  - While all 9 binaries in `src/bin/` and `horizon_wga --graph` compile and pass their unit tests, `horizon_clean` (`docs/BENCHMARK_RESULTS.md:63`) severed `0` of `6,030` / `2,776` links on real HPRC pairwise GFAs (a no-op because `horizon_wga` emits `*` placeholder sequences), and when `horizon_clean` does sever a link (`src/bin/horizon_clean.rs:265-274`), it emits an orphan `S` bubble record with no connecting `L` links.
- **Correction:** Update `45 passing tests` to `46 passing tests` in `docs/STRATEGIC_PLAN.md:3`, and qualify that `horizon_clean` is a no-op on `*`-sequence GFAs.

## F4 -- C4 (`FlatSeedIndex`): `12 bytes/entry` is post-build steady state, not build-time peak (`28 bytes/entry`), and `ContigIndex` retains `seq_rc` and `HashMap` fallbacks
- **Claim:** `C4` (`docs/STRATEGIC_PLAN.md:11`)
- **Verdict:** **NARROWED**
- **Evidence:**
  - `FlatSeedIndex` (`src/seeding/mod.rs:332-337`) stores `Vec<u64>` (8 B) + `Vec<u32>` (4 B) = `12 bytes/entry` once built, and `FlatSeedIndex::intersect` (`src/seeding/mod.rs:525-532`) computes reverse-complement coordinates as `other.seq_len - other.span - y_pos`.
  - However, `FlatSeedIndex::build` (`src/seeding/mod.rs:386-441`) allocates `pairs: Vec<(u64, u32)>` (`16 bytes/entry` with alignment padding) concurrently with `keys` (`8 B`) + `pos` (`4 B`) before `drop(pairs)` at `line 441` (`28 bytes/entry` peak during build, excluding `Vec` growth).
  - Additionally, `ContigIndex` (`src/seeding/mod.rs:542-546, 570-574`) still materializes both `seq` and `seq_rc`, and `align_contig_pair_with_config` (`src/mdl/mod.rs:210-261`) unconditionally builds two extra `FlatSeedIndex`es (`MASK_W10`) on every call and falls back to `HashMap<u64, Vec<usize>>` (`canonical_spaced_seeds_with_opts` on `seq_y` and `seq_y_rc` for `MASK_LASTZ12`) when `primary_cov < min_fallback_cov`.
- **Correction:** Qualify `12 bytes/entry` in `docs/STRATEGIC_PLAN.md:11` as post-build storage (`28 bytes/entry` during `FlatSeedIndex::build`), and note that `ContigIndex` still retains `seq_rc` and `HashMap` fallback paths.

## F5 -- C5 (`chr22` Pericentromere 5 Mbp): reproduced exact pair, unique-bp, and bubble counts
- **Claim:** `C5` (`docs/BENCHMARK_RESULTS.md:23`)
- **Verdict:** **SURVIVED**
- **Evidence:** Running `/home/rcb112/project/gemini/bench/horizon-aligner/target/release/horizon_wga /tmp/grch38_chr22_peri_5m.fa /tmp/chm13_chr22_peri_5m.fa --resolve-dups` output `contigs_x=1 contigs_y=1 pruned_pairs=0 aligned_pairs=5106541 unique_x_bp=4387655 unique_y_bp=4549858 total_bubbles=16083 elapsed_s=14.1242`, matching every count in `C5` (`horizon_wga.rs:329-338` emits `# EG:Z:` paralog tags while keeping many-to-many `aligned_pairs`, which `C5` explicitly states).

## F6 -- C6 (`horizon_gfa stats` on 12 GB HPRC graph): input `.gfa.gz` absent on disk; no code defect found
- **Claim:** `C6` (`docs/BENCHMARK_RESULTS.md:19`)
- **Verdict:** **UNCHALLENGED**
- **Evidence:** `hprc-v1.0-mc-grch38.gfa.gz` is not present on disk in the workspace or `/tmp/`, so the counts could not be re-executed. Inspection of `src/bin/horizon_gfa.rs:28-76` confirms `cmd_stats` streams via `stream_gfa` without retaining records in RAM and reads `VmHWM:` from `/proc/self/status` divided by `1024.0` (MiB, labeled `MB`).

## F7 -- Framing bias: claim list omits carried failures (`R2-HMM-XDROP`, `R2-RESOLVE-LIFT`) and stale WGA rows
- **Claim:** Framing of `C1`–`C6`
- **Verdict:** **REFUTED**
- **Evidence:**
  1. `docs/BENCHMARK_RESULTS.md:5` admits that 5 earlier WGA rows (`multiz100way`, `cactus447way`, `multiz470way`, `cactus241way chr22_random`, `GRCh38.chr22 vs CHM13.chr22`) were **not** re-run on the current build (`R2-HMM-XDROP`), yet the runner's claim list selects only the single re-run WGA row (`C2`).
  2. `docs/BENCHMARK_RESULTS.md` still has **zero** head-to-head pangenome baselines (`minimap2`/`odgi`/`vg`) required by `R2-RESOLVE-LIFT`, `horizon_clean` is a `0`-severed no-op on real HPRC pairwise GFAs (`docs/BENCHMARK_RESULTS.md:63`), and `horizon_resolve` has only a synthetic 200-segment repeat benchmark (`docs/BENCHMARK_RESULTS.md:64`).
  3. Round 2's endorsed repair to `docs/STRATEGIC_PLAN.md:13` (`C1`) was left unapplied.

| Claim | Verdict | Summary |
| :--- | :--- | :--- |
| C1 | NARROWED | `stream_gfa` ignores `C`/`J`/`U`/`O` and coerces invalid `W` coords to `0`; Round 2's endorsed doc fix was never applied |
| C2 | REFUTED | `129 target scaffolds` is false (`31` in `/tmp/dog_zoo_1m.fa` and `/tmp/pw_zoo_1m.maf`); `horizon_eval` filters pairs via truth MAF |
| C3 | NARROWED | `cargo test` runs `46` passing tests (not `45`); `horizon_clean` is a no-op on `*`-seq GFAs |
| C4 | NARROWED | `12 bytes/entry` is post-build (`28 bytes/entry` peak in `build`); `ContigIndex` retains `seq_rc` and `HashMap` fallbacks |
| C5 | SURVIVED | `horizon_wga --resolve-dups` reproduced `5,106,541` pairs, `unique_x_bp=4,387,655`, `unique_y_bp=4,549,858`, `16,083` bubbles |
| C6 | UNCHALLENGED | `hprc-v1.0-mc-grch38.gfa.gz` absent on disk; `cmd_stats` streaming/VmHWM code has no correctness flaw |
| Framing | REFUTED | Omits 5 stale WGA rows (`R2-HMM-XDROP`), missing `minimap2`/`odgi`/`vg` baselines (`R2-RESOLVE-LIFT`), and `horizon_clean` no-op |
