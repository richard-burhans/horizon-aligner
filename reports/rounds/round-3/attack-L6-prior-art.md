# Round 3 Attack -- prior-art

## Literature, Method Cards, and Prior-Round Search
- **Held literature & method cards located:** `BRIEF-prior-art.md` declares no appendix for literature or method cards. Searching `/home/rcb112/project/gemini/bench/horizon-aligner` confirmed that `horizon-aligner` has **no** `refs/catalog_index.json` and **no** method-card directory (`refs/methods/`). Held bibliography and manuscript sources were inspected at `manuscript/references.bib` (106 lines, 11 entries) and `manuscript/horizon_paper.md` (112 lines), alongside `docs/STRATEGIC_PLAN.md`, `docs/BENCHMARK_RESULTS.md`, `src/gfa/mod.rs`, `src/seeding/mod.rs`, and `src/bin/`.
- **Prior rounds checked:** `scripts/round history "" --kind adjudication` and targeted queries (`cactus241way`, `Alignathon`, `81,415,956`) returned Round 1 (`reports/rounds/round-1/adjudication-L6.md`) and Round 2 (`reports/rounds/round-2/adjudication-L6.md`). Prior Round 2 L6 adjudications relevant to Round 3 (`round 2`'s C1, C3, C5, C6, C8, and F9) are cited explicitly below.

---

## F1 -- C1: `stream_gfa` matches `round 2`'s C1 with no new prior-art defect
- **Claim:** `C1` (`docs/STRATEGIC_PLAN.md:13`)
- **Verdict:** **SURVIVED**
- **Evidence:** This claim is verbatim `round 2`'s C1, which was adjudicated **SURVIVED** (`reports/rounds/round-2/adjudication-L6.md:11-14`). Re-checking `src/gfa/mod.rs:97-219` and unit tests `streams_all_record_types_and_handles_empty_link_orients` (`lines 229-287`) and `streams_gzipped_gfa_file` (`lines 290-329`) confirms nothing changed: `stream_gfa` implements buffered line-streaming for `.gfa`, multi-member `.gfa.gz` via `MultiGzDecoder` (`line 121`), and `"-"` `stdin` (`lines 112-116`), while `docs/STRATEGIC_PLAN.md:13` cites `GBZ` / `GBZ-base` (Sirén et al.) and `gfaidx` and explicitly disclaims zero-copy.
- **Correction:** None required for C1 under `prior-art`.

## F2 -- C2: `cactus241way` "Truth MAF" framing treats Progressive Cactus as ground truth and omits `LASTZ`/`minimap2` baselines
- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:16,47`)
- **Verdict:** **NARROWED**
- **Flaw:** `docs/BENCHMARK_RESULTS.md:16,41,47` evaluates pairwise `horizon_wga` on `UCSC cactus241way chr22 1MB` (`62.3%` recall / `60.3%` precision, `61.3% F1` in `17.16s`) under the column header `Recall vs Truth MAF`. Treating the UCSC 241-way **Progressive Cactus** alignment (Armstrong et al. 2020, *Nature* 587:246–251, cited in `manuscript/references.bib:27-36`) as a ground-truth oracle ("Truth MAF") without qualification contradicts **Alignathon** (**Earl et al. 2014**, *Genome Research* 24(12):2077–2089, cited in `manuscript/references.bib:57-66`), which established that real-genome WGA outputs are *estimator* consensus/hypotheses rather than ground truth, and conflates a 241-species progressive ancestral-graph alignment with pairwise Human-vs-Dog alignment. Moreover, `docs/BENCHMARK_RESULTS.md:16,47` omits any head-to-head comparison on the same `chr22 1MB` Human-vs-Dog input against established pairwise aligners (**`LASTZ`**, Harris 2007; **`minimap2`**, Li 2018, *Bioinformatics* 34(18):3094–3100, absent from `manuscript/references.bib`).
- **Correction:** Narrow `docs/BENCHMARK_RESULTS.md:16,41,47` to label the reference as `vs. UCSC cactus241way MAF (Armstrong et al. 2020)` rather than `Truth MAF`, note that Progressive Cactus is a 241-way progressive estimator rather than a simulated ground truth (Earl et al. 2014), and add a pairwise `LASTZ` / `minimap2` baseline on the same `chr22 1MB` pair.

## F3 -- C3: Binary inventory and test/clippy counts make no prior-art assertion
- **Claim:** `C3` (`docs/STRATEGIC_PLAN.md:3`)
- **Verdict:** **UNCHALLENGED**
- **Evidence:** Extending `round 2`'s C3 (adjudicated **UNCHALLENGED** at `reports/rounds/round-2/adjudication-L6.md:21-24`), `src/bin/` contains all nine binary files (`bench_chr1.rs`, `horizon_clean.rs`, `horizon_eval.rs`, `horizon_gfa.rs`, `horizon_guide.rs`, `horizon_lift.rs`, `horizon_merge.rs`, `horizon_resolve.rs`, `horizon_wga.rs` with `--graph`), and C3 makes no prior-art or novelty assertion (test-count and clippy verification belong to reproduction/correctness).
- **Correction:** None required under `prior-art`.

## F4 -- C4: `FlatSeedIndex` adapts standard sorted flat `(key, pos)` canonical indexing (`minimap2`, `lastdb`) while still omitting minimizer/syncmer sampling
- **Claim:** `C4` (`docs/STRATEGIC_PLAN.md:11`, `docs/BENCHMARK_RESULTS.md:26`)
- **Verdict:** **NARROWED**
- **Flaw:** Addressing `round 2`'s C8 (`reports/rounds/round-2/adjudication-L6.md:47-50`), `FlatSeedIndex` (`src/seeding/mod.rs:332-546`) replaces the per-position `HashMap<u64, Vec<usize>>` with sorted parallel `Vec<u64>` + `Vec<u32>` arrays, canonical `min(fwd, rc)` keys (`line 419`), on-the-fly reverse-complement coordinates (`seq_len - span - y_pos`, `line 530`), and two-pointer sort-merge intersection (`intersect`, `lines 497-537`). From a prior-art standpoint:
  1. Flat sorted `(hash, pos)` seed tables with two-pointer/radix-sort lookup and canonical `min(fwd, rc)` encoding are established prior art in **`minimap` / `minimap2`** (Li 2016, *Bioinformatics* 32(14):2103–2110; Li 2018, *Bioinformatics* 34(18):3094–3100), **`lastdb` / `LAST`** (Kiełbasa et al. 2011, *Genome Research* 21(3):487–493), **`Kraken`** (Wood & Salzberg 2014, *Genome Biology* 15:R46), and **`Jellyfish`** (Marçais & Kingsford 2011, *Bioinformatics* 27(6):764–770) — none of which are cited in `manuscript/references.bib` or `docs/STRATEGIC_PLAN.md`.
  2. `FlatSeedIndex::build` (`src/seeding/mod.rs:388-420`) still hashes *every* valid position `0..=(seq.len() - span)` rather than using **minimizer** sampling (Roberts et al. 2004, *Bioinformatics* 20(18):3363–3369; `minimap2`, Li 2018) or **syncmers** (Edgar 2021, *PeerJ* 9:e10805). Consequently, `20 Mbp chr6 q-arm` still requires `2,409.7 MB` RSS (`0.120 GB/Mbp`, or `10.09 GB` mid-run on full `chr6`, `docs/BENCHMARK_RESULTS.md:26`), with no RSS/runtime comparison against `minimap2` or `LAST`.
- **Correction:** Narrow `docs/STRATEGIC_PLAN.md:11` and `docs/BENCHMARK_RESULTS.md:26` to state that `FlatSeedIndex` adopts standard sorted flat `(u64, u32)` canonical k-mer indexing (`minimap2`, `lastdb`) without minimizer/syncmer subsampling, cite `minimap2` (Li 2018), `lastdb` (Kiełbasa et al. 2011), and Roberts et al. (2004) in `manuscript/references.bib`, and add a `minimap2` RSS/runtime baseline on `chr6`.

## F5 -- C5: HPRC `chr22` pericentromere benchmark makes no prior-art assertion
- **Claim:** `C5` (`docs/BENCHMARK_RESULTS.md:23,59`)
- **Verdict:** **UNCHALLENGED**
- **Evidence:** Consistent with `round 2`'s C5 (`reports/rounds/round-2/adjudication-L6.md:32-35`), C5 reports single-run wall time (`14.57s`), peak RSS (`1,861.9 MB`), aligned pairs (`5,106,541`), unique-base coverage (`unique_x_bp=4,387,655`, `unique_y_bp=4,549,858`), and Horizon Bubbles (`16,083`) without a priority or novelty assertion. (Note: `docs/BENCHMARK_RESULTS.md` still omits HPRC v1.0, Liao et al. 2023 *Nature* 617:312–324, and head-to-head pericentromere baselines against `minimap2`, `wfmash`, or `Minigraph-Cactus`.)
- **Correction:** None required for C5's own factual report; cite Liao et al. (2023) and add `minimap2`/`wfmash` baselines in `docs/BENCHMARK_RESULTS.md`.

## F6 -- C6: `horizon_gfa stats` benchmark makes no prior-art assertion (though `gfatools stat` baseline is omitted)
- **Claim:** `C6` (`docs/BENCHMARK_RESULTS.md:19,54`)
- **Verdict:** **UNCHALLENGED**
- **Evidence:** C6 reports single-run wall time (`88.84s`), peak RSS (`97.1 MB`), and GFA record counts (`81,415,956` segments / `3.29 Gbp`, `112,955,105` links, `2,580` paths, `24,456` walks) on `hprc-v1.0-mc-grch38.gfa.gz` without making a priority or novelty claim, and `docs/STRATEGIC_PLAN.md:8,13` already cites `gfatools`, `odgi`, `GBZ`/`GBZ-base`, and `gfaidx`. (Note: streaming GFA statistics in low constant memory is prior art in `gfatools stat`, Li; `odgi stats`, Guarracino et al. 2022; and `vg stats`, Garrison et al. 2018 — none benchmarked head-to-head in `docs/BENCHMARK_RESULTS.md` or cited in `manuscript/references.bib`.)
- **Correction:** None required for C6's own factual report; consider adding a `gfatools stat` comparison on `hprc-v1.0-mc-grch38.gfa.gz` and citing `odgi`/`vg`/`Minigraph-Cactus` in `manuscript/references.bib`.

## F7 -- Meta-finding: Round 3 claim-list framing steers away from `round 2`'s unresolved `F9` prior-art and baseline gaps
- **Claim:** Claim-list framing (`BRIEF-prior-art.md:120-131`)
- **Verdict:** **REFUTED** *(meta-finding on claim-list framing)*
- **Flaw:** Four of the six Round 3 claims (`C1`, `C3`, `C5`, `C6`) are narrow implementation-existence, test-count, or single-run runtime/RSS statements that make zero prior-art claims—and `C1` and `C3` repeat claims already settled in Round 2 (`round 2`'s C1 `SURVIVED`, `round 2`'s C3 `UNCHALLENGED`). Meanwhile, the systemic prior-art and baseline-comparison gap already upheld in `round 2`'s `F9` (`reports/rounds/round-2/adjudication-L6.md:52-55`) remains untouched in Round 3:
  1. `manuscript/references.bib` (106 lines, 11 entries) is unchanged and still omits **every** core pangenome and pairwise-alignment prior-art paper (`minimap2` Li 2018, `last-train` Hamada et al. 2017, `lastdb` Kiełbasa et al. 2011, Roberts et al. 2004 minimizers, `GraphAligner` Rautiainen & Marschall 2020, `vg` Garrison et al. 2018, `vg giraffe` Sirén et al. 2021, `odgi` Guarracino et al. 2022, `minigraph` Li et al. 2020, `Minigraph-Cactus` Hickey et al. 2024, `smoothxg` Garrison et al. 2024, `Liftoff` Shumate & Salzberg 2021, HPRC Liao et al. 2023).
  2. `docs/BENCHMARK_RESULTS.md` still contains **zero** head-to-head comparisons against external tools (`LASTZ`, `minimap2`, `wfmash`, `gfatools`, `odgi`, `vg`), even though carried item `R2-RESOLVE-LIFT` explicitly requires `head-to-head pangenome baselines (minimap2/odgi/vg)`.

---

| C1 | C2 | C3 | C4 | C5 | C6 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SURVIVED | NARROWED | UNCHALLENGED | NARROWED | UNCHALLENGED | UNCHALLENGED |
