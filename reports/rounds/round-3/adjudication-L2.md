# Round 3 ADJUDICATION of L2 (Sampling Lens) -- Horizon Pangenome Suite

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner`  
**Attack Reviewed:** `reports/rounds/round-3/attack-L2-sampling.md` (sha256 prefix `5b713c187962` verified)  
**Prior Round & Method-Card Check:** `scripts/round history "sampling" --kind adjudication`, `scripts/round history "cactus241way" --kind adjudication`, and `scripts/round history "FlatSeedIndex" --kind adjudication` reviewed (Round 2's `F1`, `F2`, `F4`, `F7`, and `F8` bear directly on this lens). No `refs/methods/` directory exists in this repository.

---

## Finding-by-Finding Rulings

### F1 -- C2: `cactus241way chr22 1MB` samples `154.5 kb` homologous blocks (`31` contigs, not `129`), bifurcates seeding at `N_x * N_y < 10^9`, and leaves `estimate_hmm_params` prefix-biased
- **Verdict:** **NARROWED**
- **What survives:** Inspecting `/tmp/hg38_zoo_1m.fa` and `/tmp/dog_zoo_1m.fa` confirms the "`chr22 1MB`" reference is `154,483` bp (`1` contig) and the target is `132,089` bp across `31` contigs (min `1` bp, max `30,551` bp; not `129` scaffolds), `+`-strand concatenated and `.upper()`-cased (`scripts/benchmarks/extract_fasta_from_maf.py:28-45`), while `src/mdl/mod.rs:289,347` gates 2-transition `LASTZ12`, `MASK_RY20`, and `w=10` 1-transition seeding behind `(seq_x.len() as u128) * (seq_y.len() as u128) < 1_000_000_000` (splitting dog scaffolds at `6,473` bp), and `estimate_hmm_params` (`src/mdl/mod.rs:732-816`) retains the `4,000`-bp prefix / `take(64)` / `dx,dy<=500` censoring carried from Round 2 `F1`.
- **Premise check:** As Round 2 `F8` already ruled, `docs/BENCHMARK_RESULTS.md:5` already discloses *"and MAF-extracted FASTAs sample homologous blocks"*, and the reported `62.3%` recall / `60.3%` precision numbers are not shown to be false on the extracted `154.5 kb` subset itself—so C2 is narrowed rather than wholly invalidated.
- **Repair:** Split. Re-labeling C2 in `docs/BENCHMARK_RESULTS.md:16,47` as a `154.5 kb` homologous-block extract (`31` target contigs, `+`-strand, uppercased) and disclosing the `4,000`-bp prefix/`dx,dy<=500` `estimate_hmm_params` sampling and `10^9` cell-product seeding cutoff in `docs/BENCHMARK_RESULTS.md:5` is **free** (doc edit only). Benchmarking on unconcatenated 1 Mbp genomic sequence with softmasking intact is **non-free** (requires external data and benchmark execution).

### F2 -- C3: All 9 binaries sample `n = 1` tiny synthetic fixture per binary, and `cargo test` runs `46` passing + `1` ignored (not `45`)
- **Verdict:** **UPHELD**
- **What survives:** Running `cargo test` executes **46** passing unit tests (`37` in `src/lib.rs` plus `1` in each of the `9` binaries in `src/bin/*.rs`) and `1` ignored test, not `45`, and every binary in `src/bin/` still has `n = 1` unit test over a minimal synthetic fixture.
- **Repair:** Endorsed. Updating `docs/STRATEGIC_PLAN.md:3` to report `46 passing tests + 1 ignored` and qualify that each binary's unit test covers a single minimal synthetic fixture (`n = 1` per binary) is **free** (doc edit only).

### F3 -- C4: `FlatSeedIndex` censors `>= 3`-copy seeds on `>= 100 kb` contigs and C4 samples sub-chromosome `chr6` slices while full `chr6` did not finish
- **Verdict:** **NARROWED**
- **What survives:** `prepare_contig_index_opts` (`src/seeding/mod.rs:567`) sets `max_occ = 2` for `seq.len() >= 100_000` (vs `4` below `100 kb`), and `FlatSeedIndex::build` (`src/seeding/mod.rs:433`) drops all seeds with `occ > max_occ`, silently contributing to RSS reduction without precision/recall evaluation; meanwhile, `LASTZ11` also uses `HashMap`s (`src/mdl/mod.rs:312-325`) alongside `LASTZ12`, and `/tmp/hprc_chr6.gfa` and `/tmp/hprc_chr6.maf` are `0` bytes (full `170.8 Mbp` `chr6` did not finish).
- **Premise check:** `docs/STRATEGIC_PLAN.md:11` already explicitly labels the RSS numbers as `20 Mbp chr6 q-arm` and `chr6 5 Mbp pericentromere` slices (not full `chr6`), and `docs/STRATEGIC_PLAN.md:3` already mentions `LASTZ12 fallback HashMaps` (though not `LASTZ11`, and line 3 says "documented below" when line 11 omits it).
- **Repair:** Endorsed on the doc qualification. Qualifying `docs/STRATEGIC_PLAN.md:11` (that `FlatSeedIndex` drops `occ > 2` seeds on `>= 100 kb` contigs and `LASTZ12`/`LASTZ11` fallbacks still use `HashMap`s) and `docs/BENCHMARK_RESULTS.md:26` (that full `170.8 Mbp chr6` did not finish) is **free** (doc edit only).

### F4 -- C5: `chr22` Pericentromere 5 Mbp lacks a truth oracle and baselines, and silently triggers `max_occ = 2` and `10^9` seeding gating
- **Verdict:** **NARROWED**
- **What survives:** C5 (`docs/BENCHMARK_RESULTS.md:23,59`) reports `5,106,541` aligned pairs (`4.39M`/`4.55M` unique bp) in a repetitive pericentromere with zero truth oracle and zero `minimap2`/`odgi`/`vg` baseline comparisons, while silently hitting `max_occ = 2` seed censoring (`src/seeding/mod.rs:567`, `5 Mbp >= 100 kb`) and the `10^9` cell-product seeding cutoff (`src/mdl/mod.rs:289,347`, $5\text{M} \times 5\text{M} = 2.5 \times 10^{13}$).
- **Premise check:** Just as Round 2 `F4` ruled for the `5 Mbp` pericentromere slice scope, C5 already explicitly states the `5 Mbp` coordinates (`GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M`), so windowing is not undisclosed.
- **Repair:** Split. Qualifying `docs/BENCHMARK_RESULTS.md:23,59` to note no truth oracle or external baseline (`minimap2`/`odgi`/`vg`) and disclose `max_occ = 2` / `10^9` seeding gating is **free** (doc edit only). Running external baselines is **non-free**.

### F5 -- C6: `horizon_gfa stats` streams all 81.4M segments without subsampling, but is an `n = 1` counter-only pass
- **Verdict:** **UPHELD**
- **What survives:** `cmd_stats` (`src/bin/horizon_gfa.rs:39-77`) visits every `S`, `L`, `P`, and `W` record without subsampling, but is an `n = 1` single-graph, single-run measurement whose `L`/`P`/`W` callbacks (`|_| Ok(())`) only increment counters without parsing steps or retaining sequence (as reflected by `horizon_gfa extract` taking `1,255.1 MB` RSS on `chr6`).
- **Repair:** Endorsed. Qualifying `docs/BENCHMARK_RESULTS.md:19,54` as an `n = 1` counter-only streaming pass over `hprc-v1.0-mc-grch38.gfa.gz` is **free** (doc edit only).

### F6 -- C1: `stream_gfa` streams all GFA records without subsampling (confirming Round 2 L2 F7)
- **Verdict:** **UPHELD**
- **What survives:** `stream_gfa` (`src/gfa/mod.rs:97-219,229-329`) streams all GFA records sequentially without subsampling, and both unit tests (`streams_all_record_types_and_handles_empty_link_orients` and `streams_gzipped_gfa_file`) pass under `cargo test`.
- **Repair:** None required.

### F7 -- Framing bias: Runner-authored claim list omits incomplete full-`chr6`, stale WGA rows, pending `B10K`, and no-op/synthetic pangenome benchmarks
- **Verdict:** **NARROWED**
- **What survives:** The runner-authored `C1`–`C6` claim list omits the incomplete full-`chr6` run (`/tmp/hprc_chr6.gfa` and `/tmp/hprc_chr6.maf` are `0` bytes), the stale pre-`FlatSeedIndex` WGA rows, `B10K 363-way` (`Pending Data`), and the no-op/synthetic/bypassed pangenome rows (`horizon_clean`, `horizon_resolve`, `horizon_wga --graph`).
- **Premise check:** Consistent with Round 2 `F8`, `docs/BENCHMARK_RESULTS.md` itself (`lines 5, 22, 48, 58, 63-65`) already discloses the stale WGA rows, `Pending Data`, `pre-0e70693`, `*`-seq no-op, synthetic `200`-seg repeat, and `<= 8`-node bypass; only `docs/BENCHMARK_RESULTS.md:26` omits that the full-`chr6` run did not finish.
- **Repair:** Covered by `F3` (noting incomplete full-`chr6` in `docs/BENCHMARK_RESULTS.md:26` is **free**); no additional document edit required.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | `154.5 kb` / `31`-contig (`+`-strand, uppercased) extract, `10^9` cell-product seeding cutoff (`src/mdl/mod.rs:289,347`), and `estimate_hmm_params` prefix/`dx,dy<=500` censoring hold (doc fix is free; 1 Mbp genomic re-benchmark is non-free), while MAF homologous-block sampling was already disclosed in `docs/BENCHMARK_RESULTS.md:5` |
| F2 | UPHELD | `cargo test` runs `46` passing tests + `1` ignored (not `45`), and all `9` binaries in `src/bin/` have `n = 1` minimal synthetic fixture test each; doc fix is free |
| F3 | NARROWED | `FlatSeedIndex` censors `occ > 2` seeds on `>= 100 kb` contigs (`src/seeding/mod.rs:433,567`), `LASTZ11` also uses `HashMap`s, and full `chr6` did not finish (doc fix is free), though C4 already labeled the `20 Mbp`/`5 Mbp` slice scope |
| F4 | NARROWED | C5 lacks a truth oracle or `minimap2`/`odgi`/`vg` baseline and hits `max_occ = 2` / `10^9` seeding gating (doc fix is free; external baseline is non-free), while the `5 Mbp` window coordinates were already stated in C5 |
| F5 | UPHELD | `horizon_gfa stats` streams all `81.4M` segments without subsampling, but is an `n = 1` counter-only pass (`src/bin/horizon_gfa.rs:39-77`); doc qualification is free |
| F6 | UPHELD | C1 (`stream_gfa`) streams all GFA records without subsampling and passes both unit tests under `cargo test`, confirming Round 2 L2 F7 |
| F7 | NARROWED | Runner's `C1`–`C6` claim list omits incomplete full-`chr6`, stale WGA rows, pending `B10K`, and no-op/synthetic pangenome rows, though `docs/BENCHMARK_RESULTS.md` itself discloses most of them except full-`chr6` incompleteness |
