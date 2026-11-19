# Round 3 Adjudication -- L6 (prior-art)

- **Attack sealed digest verified:** `2817e0412c42c6e842e83dfa30ea70a11de4f2a34701c9af81839758f44507a4` (`reports/rounds/round-3/attack-L6-prior-art.md`) matches `2817e0412c42`.
- **Prior rounds checked:** `scripts/round history "" --kind adjudication` returned Round 1 (`reports/rounds/round-1/adjudication-L6.md`) and Round 2 (`reports/rounds/round-2/adjudication-L6.md`). Round 2 L6 already settled `C1` (`SURVIVED`), `C3` (`UNCHALLENGED`), `C5` (`UNCHALLENGED`), and the systemic citation/baseline gap (`F9`, `UPHELD`).
- **Method cards:** `horizon-aligner` declares no `refs/methods/` directory; no statistical method card governs these prior-art and implementation-scope findings.

---

## Finding-by-Finding Rulings

### F1 -- C1 (`stream_gfa` implementation and GFA v1.0/v1.1 reader)
- **Ruling:** **UPHELD**
- **What survives:** The attacker's `SURVIVED` assessment holds: `C1` is verbatim Round 2's `C1` (already adjudicated `SURVIVED` at `reports/rounds/round-2/adjudication-L6.md:11-14`), `src/gfa/mod.rs:97-219` implements `stream_gfa` / `stream_gfa_reader` with `MultiGzDecoder` (`line 121`) and `"-"` `stdin` (`lines 112-116`) plus the two cited unit tests (`lines 229-287`, `289-329`), and `docs/STRATEGIC_PLAN.md:13` cites `GBZ`/`GBZ-base` and `gfaidx` while disclaiming zero-copy.
- **Repair:** None required for `C1`.

### F2 -- C2 (`cactus241way` "Truth MAF" framing and omitted `LASTZ`/`minimap2` baselines)
- **Ruling:** **UPHELD**
- **What survives:** Both the diagnosis and the attacker's `NARROWED` verdict on `C2` hold: `docs/BENCHMARK_RESULTS.md:41,47` evaluates `cactus241way (chr22 1MB)` under the column header `Recall vs Truth MAF` without noting that Progressive Cactus (Armstrong et al. 2020, `manuscript/references.bib:27-36`) is a 241-way progressive estimator rather than simulated ground truth (Earl et al. 2014, `manuscript/references.bib:57-66`), and `docs/BENCHMARK_RESULTS.md` omits `LASTZ` / `minimap2` baselines on the same `chr22 1MB` pair.
- **Repair:** Endorsed on the text/citation leg — renaming the table column header (`docs/BENCHMARK_RESULTS.md:33,41`) to `Recall vs Reference MAF` (since the tables also span `cactus447way`, `multiz470way`, and `multiz100way`) and citing Armstrong et al. (2020) / Earl et al. (2014) is **free**; running pairwise `LASTZ` / `minimap2` baselines on `chr22 1MB` is **non-free**.

### F3 -- C3 (Binary inventory and test/clippy counts)
- **Ruling:** **UPHELD**
- **What survives:** The attacker's `UNCHALLENGED` assessment holds: `src/bin/` contains all nine binary files (`bench_chr1.rs`, `horizon_clean.rs`, `horizon_eval.rs`, `horizon_gfa.rs`, `horizon_guide.rs`, `horizon_lift.rs`, `horizon_merge.rs`, `horizon_resolve.rs`, `horizon_wga.rs` with `--graph`), and consistent with Round 2's `C3` (`reports/rounds/round-2/adjudication-L6.md:21-24`), `C3` makes no prior-art or novelty assertion.
- **Repair:** None required under L6.

### F4 -- C4 (`FlatSeedIndex` prior art and omitted minimizer/syncmer sampling)
- **Ruling:** **UPHELD**
- **What survives:** Both the diagnosis and the attacker's `NARROWED` verdict on `C4` hold: `FlatSeedIndex` (`src/seeding/mod.rs:332-538`, alongside `ContigIndex` at `541-546`) implements standard sorted parallel `Vec<u64>` + `Vec<u32>` canonical `min(fwd, rc)` indexing (`line 419`), on-the-fly reverse-complement coordinates (`line 530`), and two-pointer sort-merge intersection (`lines 497-537`)—established prior art in `minimap`/`minimap2` (Li 2016, 2018), `lastdb` (Kiełbasa et al. 2011), `Kraken` (Wood & Salzberg 2014), and `Jellyfish` (Marçais & Kingsford 2011), none cited in `manuscript/references.bib`—while `FlatSeedIndex::build` (`lines 388-420`) still indexes every valid position without minimizer (Roberts et al. 2004) or syncmer (Edgar 2021) subsampling, leaving `20 Mbp chr6 q-arm` at `2,409.7 MB` RSS (`docs/BENCHMARK_RESULTS.md:26`).
- **Repair:** Endorsed on the text/citation leg — clarifying in `docs/STRATEGIC_PLAN.md:11` and `docs/BENCHMARK_RESULTS.md:26` that `FlatSeedIndex` adapts standard sorted flat `(u64, u32)` canonical k-mer indexing without minimizer/syncmer sampling and adding `minimap2` / `lastdb` / Roberts et al. (2004) to `manuscript/references.bib` is **free**; adding a `minimap2` RSS/runtime baseline on `chr6` is **non-free**.

### F5 -- C5 (HPRC `chr22` pericentromere benchmark)
- **Ruling:** **UPHELD**
- **What survives:** The attacker's `UNCHALLENGED` assessment holds: consistent with Round 2's `C5` (`reports/rounds/round-2/adjudication-L6.md:32-35`), `C5` (`docs/BENCHMARK_RESULTS.md:23,59`) reports single-run runtime, RSS, aligned pairs, unique-base coverage, and Horizon Bubbles without a priority or novelty assertion.
- **Repair:** None required for `C5`'s factual report; citing Liao et al. (2023) in `manuscript/references.bib` is **free**, while adding `minimap2`/`wfmash` pericentromere baselines is **non-free**.

### F6 -- C6 (`horizon_gfa stats` benchmark)
- **Ruling:** **UPHELD**
- **What survives:** The attacker's `UNCHALLENGED` assessment holds: `C6` (`docs/BENCHMARK_RESULTS.md:19,54`) reports single-run wall time, peak RSS, and GFA record counts without a priority or novelty claim, and `docs/STRATEGIC_PLAN.md:8,13` already cites `gfatools`, `odgi`, `GBZ`/`GBZ-base`, and `gfaidx`.
- **Repair:** None required for `C6`'s factual report; adding `odgi`/`vg`/`Minigraph-Cactus` citations to `manuscript/references.bib` is **free**, while benchmarking `gfatools stat` on `hprc-v1.0-mc-grch38.gfa.gz` is **non-free**.

### F7 -- Meta-finding (Round 3 claim-list framing and persistent `round 2` `F9` citation/baseline gaps)
- **Ruling:** **UPHELD**
- **What survives:** Verified that four of the six Round 3 claims (`C1`, `C3`, `C5`, `C6`) make no prior-art assertions (`C1` and `C3` repeating settled Round 2 claims), while the systemic prior-art and baseline gaps already upheld in Round 2's `F9` (`reports/rounds/round-2/adjudication-L6.md:52-55`) remain unresolved: `manuscript/references.bib` (106 lines, 11 entries) still omits all core pangenome/aligner references (`minimap2`, `last-train`, `lastdb`, Roberts et al. 2004, `GraphAligner`, `vg`, `vg giraffe`, `odgi`, `minigraph`, `Minigraph-Cactus`, `smoothxg`, `Liftoff`, HPRC Liao et al. 2023), and `docs/BENCHMARK_RESULTS.md` still has zero head-to-head external baselines despite carried item `R2-RESOLVE-LIFT`.
- **Repair:** Adding the missing prior-art citations to `manuscript/references.bib` is **free**; running head-to-head external tool benchmarks is **non-free**.

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `SURVIVED` assessment verified in `src/gfa/mod.rs:97-219,229-329` and `docs/STRATEGIC_PLAN.md:13`; verbatim Round 2 `C1` |
| F2 | UPHELD | `Recall vs Truth MAF` header (`docs/BENCHMARK_RESULTS.md:33,41,47`) treats Progressive Cactus as truth without citing Alignathon qualification and omits `LASTZ`/`minimap2` baselines; text/citation fix is free, external baselines are non-free |
| F3 | UPHELD | `UNCHALLENGED` assessment verified against `src/bin/`; `C3` makes no prior-art assertion |
| F4 | UPHELD | `FlatSeedIndex` (`src/seeding/mod.rs:332-538`) adapts standard sorted flat `(u64, u32)` canonical indexing (`minimap2`/`lastdb`, uncited) without minimizer/syncmer sampling; text/citation fix is free, `minimap2` baseline is non-free |
| F5 | UPHELD | `UNCHALLENGED` assessment verified in `docs/BENCHMARK_RESULTS.md:23,59`; `C5` makes no prior-art assertion |
| F6 | UPHELD | `UNCHALLENGED` assessment verified in `docs/BENCHMARK_RESULTS.md:19,54`; `C6` makes no prior-art assertion |
| F7 | UPHELD | Verified `manuscript/references.bib` (106 lines, 11 entries) and `docs/BENCHMARK_RESULTS.md` still omit Round 2 `F9` core citations (free) and `R2-RESOLVE-LIFT` head-to-head baselines (non-free) |
