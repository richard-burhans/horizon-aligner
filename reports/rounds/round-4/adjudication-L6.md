# Round 4 Adjudication -- L6 (prior-art)

- **Attack sealed digest verified:** `a5bf242ef800` (`reports/rounds/round-4/attack-L6-prior-art.md`) matches `sha256sum` prefix `a5bf242ef800`.
- **Prior rounds checked:** `scripts/round history "prior-art" --kind adjudication` returned Round 1 (`reports/rounds/round-1/adjudication-L6.md`), Round 2 (`reports/rounds/round-2/adjudication-L6.md`), and Round 3 (`reports/rounds/round-3/adjudication-L6.md`). Round 1 `F1` (CLI pipe composition), Round 2 `F2`/`F4`/`F9` (`last-train`, `Lunter`, external baselines, systemic `references.bib` omissions), and Round 3 `F2`/`F4`/`F6`/`F7` (`FlatSeedIndex` prior art, HPRC/Minigraph-Cactus/`odgi`/`vg` citations, persistent `F9` gaps) already established the distinction between **free** text/citation repairs and **non-free** external benchmark runs.
- **Method cards:** `horizon-aligner` declares no `refs/methods/` directory; no statistical method card governs these prior-art and citation-hygiene findings.
- **Manuscript vs. attacked snapshot:** `diff -u reports/rounds/round-4/ATTACKED-horizon_paper.md manuscript/horizon_paper.md` is empty; `../library` is absent on disk, and all 24 cited bioinformatics keys and DOIs in `attack-L6-prior-art.md` were independently verified against `refs/catalog_index.json`.

---

## Finding-by-Finding Rulings

### F1 -- C1: Binary count inconsistency (10 vs. 9) and uncited composable genomic/pangenome toolkit prior art
- **Ruling:** **UPHELD**
- **What survives:** Both legs of the attacker's `NARROWED` diagnosis hold: `src/bin/` contains 10 `.rs` executables (`horizon.rs:7-53` dispatches 9 subcommands while itself being a 10th binary), `manuscript/horizon_paper.md:101` states "across all `10` CLI binaries" in contradiction to "nine composable executables" at `line 18`, and `guarracino_2022_odgi` appears in `# References` (`line 122`) without any body-text citation while `vg` (`garrison_2018_variation`) and `BEDTools` (`quinlan-hall_2010_bedtools`) are omitted.
- **Repair:** Reconciling `line 18` with `line 101` ("ten CLI binaries, nine of which are dispatched via `horizon <subcommand>`"), citing `guarracino_2022_odgi` in the body text, and adding `garrison_2018_variation` and `quinlan-hall_2010_bedtools` to `manuscript/references.bib` and the Introduction/Discussion is **free** (endorsed; text/citation edit only, no re-measurement).

### F2 -- C2: `FlatSeedIndex` prefix-bucketed Rayon speedup omits `minimap2`/`FastGA`/`KegAlign`/`KMC 3` prior art and external baselines
- **Ruling:** **UPHELD** (diagnosis holds; repair split into free and non-free legs)
- **What survives:** Consistent with Round 3 `F4`, `FlatSeedIndex::build_with_sdust` (`src/seeding/flat_index.rs:12-120`) indexes every unmasked k-mer (`line 81`) into a 16-bit prefix-bucketed sorted `(u64, u32)` array (`buckets: Vec<u32>`, `line 17`), `C2` (`manuscript/horizon_paper.md:101`) reports only an internal commit-over-commit speedup (`53130b8` vs. `1e47fe0`, explicitly labeled in Figure 5's caption at `line 105` but phrased as "our `1,062.45 s` baseline" at `line 101`), and `li_2018_minimap2`, `kokot_2017_counting`, `myers_2024_fastga`, and `gulhan_2025_kegalign` are all present in `refs/catalog_index.json` yet absent from `manuscript/references.bib` and `manuscript/horizon_paper.md`.
- **Repair:** Qualifying `line 101` as an internal commit-over-commit speedup (`53130b8` vs. `1e47fe0`) and citing `li_2018_minimap2`, `myers_2024_fastga`, `gulhan_2025_kegalign`, and `kokot_2017_counting` is **free** (endorsed); running head-to-head `minimap2`/`LASTZ`/`FastGA` baselines on `GRCh38.chr6` × `CHM13.chr6` is **non-free**.

### F3 -- C3: HPRC `v1.0` Minigraph-Cactus benchmark omits HPRC/Minigraph-Cactus/`vg` citations and `odgi`/`vg`/`gfatools` baselines
- **Ruling:** **UPHELD** (diagnosis holds; repair split into free and non-free legs)
- **What survives:** Consistent with Round 3 `F6` and `F7`, `C3` (`manuscript/horizon_paper.md:101`) reports single-tool runtime, RSS, and record counts without making a priority claim, but the manuscript benchmarks directly on `hprc-v1.0-mc-grch38.gfa.gz` while omitting `liao_2023_draft` (HPRC), `hickey_2023_pangenome` (Minigraph-Cactus), `li_2020_design` (`minigraph`), and `garrison_2018_variation` (`vg`), leaving `guarracino_2022_odgi` uncited in the body text (`line 122`), and providing no external `odgi`/`vg`/`gfatools` comparison.
- **Repair:** Citing `liao_2023_draft`, `hickey_2023_pangenome`, `guarracino_2022_odgi`, and `garrison_2018_variation` in the body text and `manuscript/references.bib` is **free** (endorsed); running head-to-head `odgi`/`gfatools`/`vg` stats and extraction baselines on `hprc-v1.0-mc-grch38.gfa.gz` is **non-free**.

### F4 -- C4: Fabricated `Sena L, Paten B, Eizenga JM (2026)` citation and ungrounded "Horizon Bubble" vs. topological bubble decomposition
- **Ruling:** **UPHELD**
- **What survives:** Verified that `manuscript/horizon_paper.md:130` fabricates `Sena L, Paten B, Eizenga JM (2026). Linear-time bubble decomposition via SPQR trees. Journal of Computational Biology` whereas both `refs/catalog_index.json` and `manuscript/references.bib:241-247` hold `Sena, F., Politov, A., Moumard, C., Cairo, M., Rizzi, R., Cáceres, M., Schmidt, S., Harviainen, J., & Tomescu, A. I. (2026). Identifying bubble-like subgraphs in linear-time via a unified SPQR-tree framework. arXiv (DOI: 10.48550/arXiv.2604.08071)`, and neither `Paten et al. (2018)` (`line 127`) nor `Sena et al. (2026)` (`line 130`) is cited in the body text to distinguish topological superbubbles/SPQR decompositions from `horizon_clean`'s local GC/MDL `walk_support == 0` edge severing (`HB:i:1`).
- **Repair:** Replacing `manuscript/horizon_paper.md:130` with the exact `Sena et al. (2026)` citation from `manuscript/references.bib:241-247`, citing `paten_2018_superbubbles-ultrabubbles-cacti` and `sena_2026_bubble-subgraphs-spqr` in the body text, and clarifying that `HB:i:1` marks severed discordant `L`-edges rather than topological superbubble/SPQR decompositions is **free** (endorsed).

### F5 -- C5: Posterior-support abstention and parameter-free transfer omit `Lunter et al. (2008)` and `last-train` (`Hamada et al. 2017`)
- **Ruling:** **UPHELD**
- **What survives:** Repeating the settled omissions from Round 2 `F2` and Round 3 `F7`, `manuscript/horizon_paper.md:16` frames conventional aligners as emitting deterministic blocks without per-interval posterior support while omitting `lunter_2008_uncertainty` (Pair-HMM posterior decoding) and `hamada_2017_training` (`last-train` automatic parameter estimation) from both `manuscript/horizon_paper.md` and `manuscript/references.bib` (and omitting `lynch_2025_cannabis-pangenome-synthases` from the bibliography).
- **Repair:** Citing `lunter_2008_uncertainty` and `hamada_2017_training` in `manuscript/horizon_paper.md:16,36,65` and `manuscript/references.bib` and qualifying `line 16` so it does not imply posterior abstention or automatic parameter adaptation are absent from prior aligners is **free** (endorsed).

### F6 -- C6: E-graph duplication recovery (`1.00` vs. `0.00`) benchmarks only against a strawman greedy 1-to-1 filter and omits duplication-aware WGA/synteny prior art
- **Ruling:** **UPHELD**
- **What survives:** `C6` (`manuscript/horizon_paper.md:99`) compares e-graph duplication recovery (`1.00` `dup(X)`) solely against "greedy one-to-one synteny filtering" (`0.00`, which discards paralogs by definition) while omitting all duplication-aware WGA and synteny prior art held in `refs/catalog_index.json` (`hou_2006_controlling-size`, `hou_2009_aligning`, `paten_2008_enredo`, `darling_2010_progressivemauve`, `minkin_2020_scalable`, `wang_2012_mcscanx`, `kent_2003_evolutions-cauldron`).
- **Repair:** Citing `hou_2009_aligning`, `paten_2008_enredo`, `darling_2010_progressivemauve`, `minkin_2020_scalable`, and `wang_2012_mcscanx` in `manuscript/horizon_paper.md:44,99` and `manuscript/references.bib` and qualifying `C6` as a comparison against 1-to-1 greedy filtering rather than duplication-aware synteny methods is **free** (endorsed).

### F7 -- Meta-finding: Uncited `# References` entries and `references.bib` / `horizon_paper.md` drift
- **Ruling:** **UPHELD**
- **What survives:** Cross-checking `manuscript/horizon_paper.md` against `manuscript/references.bib` confirms all three counts: seven `# References` entries are never cited in the body text (`Guarracino 2022` `L122`, `Heumos 2024` `L123`, `Morgulis 2006` `L126`, `Paten 2018` `L127`, `Sena 2026` `L130`, `Foster 2007` `L134`, `Hickey 2013` `L136`), two `.bib` entries (`jiang_2008_universal-dna-representation`, `miller_2020_population-structure-clustering`) are absent from `manuscript/horizon_paper.md`, and `L130` is fabricated.
- **Repair:** Citing every `# References` entry in the body text (or removing unused ones), removing the two orphan `.bib` entries, and fixing `manuscript/horizon_paper.md:130` to match `sena_2026_bubble-subgraphs-spqr` is **free** (endorsed).

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | Verified 10 binaries in `src/bin/` vs. 9 in `C1` (`line 18`) and 10 at `line 101`, plus uncited `odgi`/`vg`/`BEDTools` prior art; text/citation fix is free |
| F2 | UPHELD | Verified `FlatSeedIndex::build_with_sdust` (`src/seeding/flat_index.rs:12-120`), internal commit baseline (`1e47fe0` vs. `53130b8`), and uncited `minimap2`/`KMC 3`/`FastGA`/`KegAlign`; text/citation fix is free, external baselines are non-free |
| F3 | UPHELD | Verified `liao_2023_draft`, `hickey_2023_pangenome`, `li_2020_design`, and `garrison_2018_variation` omitted and `guarracino_2022_odgi` uncited in body; citations are free, external GFA baselines are non-free |
| F4 | UPHELD | Verified fabricated `Sena L, Paten B, Eizenga JM (2026)` at `horizon_paper.md:130` vs. `references.bib:241-247`, uncited `Paten 2018`/`Sena 2026` in body, and edge-severing vs. topological bubble conflation; fix is free |
| F5 | UPHELD | Verified `lunter_2008_uncertainty` and `hamada_2017_training` still omitted from `horizon_paper.md` and `references.bib` (repeating Round 2 `F2` / Round 3 `F7`); text/citation fix is free |
| F6 | UPHELD | Verified `C6` benchmarks only against 1-to-1 greedy filtering (`0.00` by definition) and omits `Enredo`/`Hou et al.`/`progressiveMauve`/`SibeliaZ`/`MCScanX`; text/citation qualification is free |
| F7 | UPHELD | Verified 7 uncited `# References` entries, 2 orphan `.bib` entries (`jiang_2008`, `miller_2020`), and fabricated `L130`; bibliography hygiene fix is free |
