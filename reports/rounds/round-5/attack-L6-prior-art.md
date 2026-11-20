# Round 5 Attack -- L6 (prior-art)

- **Scope & constraints:** `../library` is not present on disk, so this attack worked from `refs/catalog_index.json` (pinned at `refs/LIBRARY_VERSION` = `d4df9122eddd84694b2c3215e5478c46f4314044`), `manuscript/references.bib`, `reports/rounds/round-5/ATTACKED-horizon_paper.md`, `reports/rounds/round-5/ATTACKED-horizon_companion.md`, and `reports/rounds/round-5/round.json` alone; no paper's internal claims are asserted from memory.
- **Method cards:** This project keeps **no method cards** (`refs/methods/` does not exist).
- **Prior rounds checked:** Checked `scripts/round history` for `prior-art`, `C1`, `minimap2`, and `FastGA`, and read `reports/rounds/round-4/adjudication-L6.md`. Round 4 `F1`–`F7` were all `UPHELD`; between Round 4 and Round 5, `ATTACKED-horizon_paper.md` added body-text and `# References` entries for `vg`, `BEDTools`, `minimap2`, `KMC 3`, `FastGA`, `KegAlign`, `Liao et al. (2023)`, `Hickey et al. (2023)`, `minigraph`, `Lunter et al. (2008)`, `last-train`, `Enredo`, `Hou et al. (2009)`, `progressiveMauve`, `SibeliaZ`, and `MCScanX`, and corrected `Sena et al. (2026)`. However, **none** of the 15 newly added paper references were added to `manuscript/references.bib`, the two orphan `.bib` entries remain, and a wrong DOI was introduced for `Gulhan et al. (2025)`.

---

## F1 -- C1: `vg` subgraph comparison is apples-to-oranges and omits `odgi` / `gfatools` baselines and pangenome graph construction prior art

- **Claim:** `C1` (`ATTACKED-horizon_paper.md:10, 18, 46-51, 82`; `round.json:13-17`)
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. `C1` compares `horizon_gfa stats` (`0.08 s`, `22.6 MB` RSS)—which `"performs single-pass counter-only graph stats (discarding records after counting)"` (`ATTACKED-horizon_paper.md:47`)—against `vg stats -z -l -N` (`4.2984 s`, `77.0 MB` RSS) and `vg convert -g -p` (`1.8292 s`, `76.6 MB` RSS), which construct in-memory bidirected handle-graph structures rather than performing streaming line counting. Meanwhile, the actual streaming/compact pangenome baselines `odgi` (`guarracino_2022_odgi`) and `gfatools` remain unbenchmarked (acknowledged as pending in `ATTACKED-horizon_paper.md:96`).
  2. `C1`'s summary (`round.json:13`) and `ATTACKED-horizon_paper.md:82` also report `horizon_clean`, `horizon_resolve`, and `horizon_wga --graph` on the `10,000`-segment `chr22` subgraph with zero external baselines (`vg map`/`giraffe`, `minigraph`, `abPOA`, `odgi`).
  3. In `ATTACKED-horizon_paper.md:46`, pangenome graph construction/evaluation prior art held in `refs/catalog_index.json` is omitted: `garrison_2024_building` (*Building pangenome graphs*, Nature Methods 2024, DOI `10.1038/s41592-024-02430-3`), `andreace_2023_comparing` (*Comparing methods for constructing and representing human pangenome graphs*, Genome Biology 2023, DOI `10.1186/s13059-023-03098-2`), `siren_2024_personalized-pangenome` (Nature Methods 2024, DOI `10.1038/s41592-024-02407-2`), and `kopalli_2025_toward` (GigaScience 2025, DOI `10.1093/gigascience/giaf121`).
- **Prior-round status:** Round 4 `F1` and `F3` (`adjudication-L6.md:12-15, 22-25`) upheld missing `vg`/`odgi`/`BEDTools` citations and missing external graph baselines. Round 5 added `vg stats`/`vg convert` on the `10k` subgraph and cited `vg`/`odgi`/`BEDTools` in `horizon_paper.md`, but `vg` vs `horizon_gfa stats` is apples-to-oranges and `odgi`/`gfatools` remain unbenchmarked.
- **Correction:** Qualify `ATTACKED-horizon_paper.md:10, 82` to note that `vg stats`/`vg convert` build in-memory graph representations whereas `horizon_gfa stats` performs counter-only line scanning, add `odgi stats` and `gfatools stat` baselines, and cite `garrison_2024_building` and `andreace_2023_comparing`.

---

## F2 -- C2: `minimap2 -x asm5` outperforms Horizon on speed, exact recall, and precision while `FastGA`/`LASTZ`/`KegAlign` and `ntHash2`/`k`-nonical prior art are missing

- **Claim:** `C2` (`ATTACKED-horizon_paper.md:10, 28, 76, 90`; `round.json:25-29`)
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. `C2` compares **many-to-many** `horizon_wga --resolve-dups` (`258.88–304.53 s`, `97.96%` exact recall, `48.63%` exact pair precision, `~2.6×` multiplicity) against **1-to-1 primary** `minimap2 -x asm5` (`46.95 s`, `99.71%` exact recall, `99.75%` exact precision). `minimap2` is `5.5×–6.5×` faster with higher exact recall and precision, yet `C2` omits secondary/many-to-many `minimap2` modes and still omits head-to-head `FastGA` (`myers_2024_fastga`), `LASTZ` (`harris_2007_improved`), and `KegAlign` (`gulhan_2025_kegalign`) runs (`ATTACKED-horizon_paper.md:59, 96`).
  2. `ATTACKED-horizon_paper.md:28` discusses `FlatSeedIndex` strand-symmetric spaced-seed indexing (`KP-14`) without citing directly relevant held prior art in `refs/catalog_index.json`: `mohamadi_2016_nthash` (`10.1093/bioinformatics/btw397`), `kazemi_2022_nthash2` (`10.1093/bioinformatics/btac564`, recursive spaced-seed hashing), `marcais_2024_knonical` (`10.1093/bioinformatics/btae629`, sketching with reverse complements), `wittler_2023_canonical-kmer-encoding` (`10.24072/pcjournal.323`), and `shivakumar_2025_mumemto` (`10.1186/s13059-025-03644-0`).
- **Prior-round status:** Round 4 `F2` (`adjudication-L6.md:17-20`) upheld missing `minimap2`/`FastGA`/`KegAlign`/`KMC 3` citations and baselines. Round 5 added those four citations to `horizon_paper.md:28` and ran `minimap2 -x asm5`, but did not add them to `manuscript/references.bib` (see `F6`) and still lacks `FastGA`/`LASTZ`/`KegAlign` baselines.
- **Correction:** Qualify `C2` (`ATTACKED-horizon_paper.md:76`) as slower (`258.88–304.53 s` vs `46.95 s`) and lower-precision (`48.63%` vs `99.75%`) than `1-to-1` `minimap2 -x asm5`, add `FastGA`/`LASTZ`/`KegAlign` or secondary-alignment `minimap2` baselines, and cite `kazemi_2022_nthash2` and `marcais_2024_knonical`.

---

## F3 -- C3: Multi-species `horizon_merge` / `horizon_guide` and `horizon_lift` lack external baselines and omit UCSC netting (`Kent et al. 2003`), `Mugsy`, and `Liftoff`

- **Claim:** `C3` (`ATTACKED-horizon_paper.md:10, 53, 74, 82`; `round.json:37-41`)
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. `C3` reports single-tool runtime/RSS/block counts for `horizon_merge` (`92` `chrM` species, `21.78 s`, `26.6 MB`) and `horizon_guide --progressive --hal --rbest` (`5` mammals, `0.859 s`, `1,109` MAF blocks, `4` HAL bottom segments) with **zero** head-to-head comparison against existing multi-species aligners/mergers (`MULTIZ`/`TBA` `blanchette_2004_threaded-blockset-aligner`, `Progressive Cactus` `armstrong_2020_progressive`, `progressiveMauve` `darling_2010_progressivemauve`, `SibeliaZ` `minkin_2020_scalable`, or `hal2maf` `hickey_2013_hierarchical`). Similarly, `horizon_lift` (`round.json:37`, `ATTACKED-horizon_paper.md:74`) has no baseline against `halLiftover`, `odgi position`, `vg surject`, or `Liftoff`.
  2. `ATTACKED-horizon_paper.md:53, 82` uses the UCSC term `"netting"` (`horizon_merge` netted `92` vertebrate `chrM` species; `"best-score interval netting"`) while omitting `kent_2003_evolutions-cauldron` (*Evolution's cauldron: Duplication, deletion, and rearrangement in the mouse and human genomes*, PNAS 2003, DOI `10.1073/pnas.1932072100`)—held in `refs/catalog_index.json:815-818`. It also omits `angiuoli_2011_mugsy` (`10.1093/bioinformatics/btq665`), `dubchak_2009_multiple` (`10.1101/gr.081778.108`), and `shumate_2021_liftoff` (`10.1093/bioinformatics/btaa1016`).
- **Prior-round status:** Round 4 `F3` and Round 2 `F4` (`adjudication-L6.md:22-25, 29`) established that single-tool runtime/RSS measurements without external baselines are narrowed rather than refuted, and `kent_2003_evolutions-cauldron` was noted in Round 4 `F6`.
- **Correction:** Qualify `C3` as single-tool resource/count measurements pending external `MULTIZ`/`Cactus`/`halLiftover` baselines, and cite `kent_2003_evolutions-cauldron`, `angiuoli_2011_mugsy`, and `shumate_2021_liftoff`.

---

## F4 -- C4: Pairwise WGA benchmarks omit external aligner comparisons (`LASTZ`, `minimap2`, `FastGA`, `last-train`) and `tantan` repeat-masking prior art

- **Claim:** `C4` (`ATTACKED-horizon_paper.md:10, 59-70`; `round.json:50-54`)
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. `ATTACKED-horizon_paper.md:59` explicitly acknowledges that `"head-to-head external pairwise aligner re-runs (minimap2, FastGA, LASTZ) remain pending external tool availability on the benchmark host."` Consequently, none of the eight Table 1 benchmarks (`100-way chrM`, `cactus447way`, `multiz470way`, `cactus241way`, `B10K 363-way`) or the `18%`-GC random DNA null include a baseline comparison against `LASTZ`, `minimap2`, `FastGA`, or `last-train`.
  2. While `ATTACKED-horizon_paper.md:16, 61` now cites `earl_2014_alignathon` and labels Table 1 `Recall vs Reference MAF (%)` (addressing Round 3 `F2`), it omits `chen_2010_comparative` (`10.1038/nbt.1637`), `pollard_2004_benchmarking-tools-alignment` (`10.1186/1471-2105-5-6`), `wu_2022_multiple` (`10.1002/tpg2.20204`), and `frith_2011_tantan` (`10.1093/nar/gkq1212`, directly relevant to suppressing AT-rich/simple-repeat false positives alongside `morgulis_2006_dustmasker`).
- **Prior-round status:** Round 3 `F2` (`scripts/round history` `adjudication-L6.md:51`) and Round 4 `F5` (`adjudication-L6.md:32-35`) upheld missing `Alignathon` qualification and `last-train`/`Lunter` citations. `horizon_paper.md` now cites `earl_2014_alignathon`, `hamada_2017_training`, and `lunter_2008_uncertainty` in the body text, though the latter two remain absent from `manuscript/references.bib`.
- **Correction:** Qualify `C4` as reference-MAF concordance without external aligner baselines (`LASTZ`/`minimap2`/`FastGA`/`last-train`), and cite `frith_2011_tantan` and `chen_2010_comparative`.

---

## F5 -- C5: `horizon_companion.md` leaves two `# 5. References` entries uncited in its body and omits prior art for Eisenstein DNA codes, banded DP, and Karlin–Altschul scaling

- **Claim:** `C5` (`ATTACKED-horizon_companion.md:53-56, 87-96, 124, 131-145, 186-193`; `round.json:60-64`)
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. `ATTACKED-horizon_companion.md` `# 5. References` (`lines 189, 191`) lists `Allison L, Wallace CS, Yee CN (1992)` and `Foster JN, Greenwald MB, Moore JT, Pierce BC, Schmitt A (2007)`, neither of which is cited anywhere in the companion body text (even though Section 2.2 `lines 53-56` states the `GETPUT`/`PUTGET`/`Constant Complement` lens laws and Section 3.3 `line 124` describes the three-state MDL pair-HMM).
  2. Section 3.1 (`lines 87-96`) presents `KP-14` Eisenstein nucleotide encoding (`{1, -1, ω, -ω}`) without citing `sajjad_2023_eisenstein-bch` (`10.3934/math.20231508`), `sajjad_2026_eisenstein-dna-bch` (`10.3934/math.20261378`), `abualrub_2006_cyclic-codes-gf4-dna`, or `rykov_2001_quaternary-cyclic-dna-codes` from `refs/catalog_index.json`.
  3. Section 3.3 (`lines 131-145`) invokes `"Continuous Karlin–Altschul search-space scaling"` (`line 131`) without citing Karlin & Altschul or `rissanen_1978_shortest-data-description`/`grunwald_2007_mdl-principle` in `horizon_companion.md`, and presents `BandGrid` diagonal flat-buffer indexing (`lines 139-145`) without citing banded/diagonal DP prior art (`ukkonen_1985_approximate-string-matching`, `chao_1994_linear-space`).
- **Correction:** Cite `Allison et al. (1992)` in Section 3.3 and `Foster et al. (2007)` in Section 2.2 of `horizon_companion.md`, and add prior-art citations for `sajjad_2026_eisenstein-dna-bch` and `ukkonen_1985_approximate-string-matching` / `chao_1994_linear-space`.

---

## F6 -- C6: `C6`'s claim of "accurate prior-art attribution" is refuted by a wrong DOI for `Gulhan et al. (2025)`, 15 `horizon_paper.md` references missing from `manuscript/references.bib`, and uncited companion prior art

- **Claim:** `C6` (`ATTACKED-horizon_paper.md:113, 100-137`; `ATTACKED-horizon_companion.md:157-162, 189, 191`; `manuscript/references.bib:1-247`; `round.json:70-74`)
- **Verdict:** **REFUTED**
- **Concrete flaw:** `round.json:70` explicitly asserts that `"The companion and main manuscript match the current Rust implementation on ... accurate prior-art attribution."` That assertion is false on four concrete counts:
  1. **Wrong DOI:** `ATTACKED-horizon_paper.md:113` cites `Gulhan A, et al. (2025). KegAlign: optimizing pairwise alignments with diagonal partitioning. Genome Biology (10.1186/s13059-025-03468-0)`. In `refs/catalog_index.json:560-563`, `gulhan_2025_kegalign` has DOI `10.1186/s13059-025-03830-0` (not `10.1186/s13059-025-03468-0`).
  2. **15 paper references absent from `manuscript/references.bib`:** `ATTACKED-horizon_paper.md` `# References` contains `38` entries, whereas `manuscript/references.bib` contains only `25` entries. All `15` entries added to `horizon_paper.md` after Round 4 (`Darling 2010`, `Garrison 2018`, `Gulhan 2025`, `Hamada 2017`, `Hickey 2023`, `Hou 2009`, `Kokot 2017`, `Li 2018`, `Li 2020`, `Liao 2023`, `Lunter 2008`, `Myers 2025`, `Paten 2008`, `Quinlan 2010`, `Wang 2012`) are missing from `manuscript/references.bib`.
  3. **2 orphan `.bib` entries retained:** `jiang_2008_universal-dna-representation` (`manuscript/references.bib:153-162`) and `miller_2020_population-structure-clustering` (`manuscript/references.bib:186-195`) remain in `manuscript/references.bib` despite being cited in neither `horizon_paper.md` nor `horizon_companion.md` (and despite Round 4 `F7` upholding their removal).
  4. **Inversion/duplication prior-art gaps:** `ATTACKED-horizon_paper.md` omits `romain_2026_inversions-pangenome-graphs` (*Representation and detection of inversions in pangenome graphs*, NAR Genomics and Bioinformatics 2026, DOI `10.1093/nargab/lqag125`), and `ATTACKED-horizon_companion.md:157-162` describes e-graph saturation for tandem duplications and nested inversions without citing `Paten et al. (2008)`, `Hou et al. (2009)`, `Darling et al. (2010)`, `Minkin & Medvedev (2020)`, or `Wang et al. (2012)`.
- **Prior-round status:** Round 4 `F7` (`adjudication-L6.md:42-45`) upheld `references.bib` / `horizon_paper.md` drift and orphan `.bib` entries. In Round 5, the author updated `horizon_paper.md` `# References` without syncing `manuscript/references.bib` and introduced the wrong `Gulhan et al. (2025)` DOI.
- **Correction:** Fix `ATTACKED-horizon_paper.md:113` to DOI `10.1186/s13059-025-03830-0`, add all `15` missing entries to `manuscript/references.bib`, remove `jiang_2008_universal-dna-representation` and `miller_2020_population-structure-clustering` from `manuscript/references.bib`, and cite `romain_2026_inversions-pangenome-graphs` (and the duplication-aware WGA references in `horizon_companion.md:157-162`).

---

## F7 -- Meta-finding: Runner-authored claim list steers away from missing external baselines and bibliography drift

- **Claim:** Meta-finding on `round.json` (`C1`–`C6`, `round.json:111-112`)
- **Verdict:** **UPHELD** (framing observation)
- **Concrete flaw:** `round.json:111-112` notes `runner_is_author: true` and `claim_list_author: "runner"`. `C1` highlights `horizon_gfa stats` vs `vg stats`/`vg convert` (an apples-to-oranges counter-only vs in-memory graph comparison) while omitting `odgi` and `gfatools`; `C2` highlights an internal `3.97×` commit-over-commit speedup (`1e47fe0` → `8cd5db0`) while burying that `minimap2 -x asm5` is `5.5×–6.5×` faster (`46.95 s` vs `258.88–304.53 s`) with higher exact recall (`99.71%` vs `97.96%`) and precision (`99.75%` vs `48.63%`); `C3` and `C4` report single-tool numbers with zero external baselines; and `C6` asserts `"accurate prior-art attribution"` despite 15 missing `references.bib` entries and a wrong DOI.

---

| Claim | Verdict | Summary |
|---|---|---|
| C1 | NARROWED | `horizon_gfa stats` vs `vg stats`/`vg convert` is counter-only vs in-memory graph construction; `odgi`/`gfatools` baselines and `garrison_2024_building`/`andreace_2023_comparing` omitted |
| C2 | NARROWED | `minimap2 -x asm5` is 5.5×–6.5× faster with higher exact recall/precision; `FastGA`/`LASTZ`/`KegAlign` baselines and `kazemi_2022_nthash2`/`marcais_2024_knonical` omitted |
| C3 | NARROWED | Single-tool `horizon_merge`/`horizon_guide`/`horizon_lift` timings lack `MULTIZ`/`Cactus`/`halLiftover` baselines and omit `kent_2003_evolutions-cauldron`, `angiuoli_2011_mugsy`, `shumate_2021_liftoff` |
| C4 | NARROWED | All 8 Table 1 WGA slices and the 18%-GC null lack external aligner baselines (`LASTZ`/`minimap2`/`FastGA`/`last-train`); `frith_2011_tantan` and `chen_2010_comparative` omitted |
| C5 | NARROWED | `horizon_companion.md` leaves `Allison 1992` and `Foster 2007` uncited in body, invokes Karlin–Altschul uncited, and omits `sajjad_2026_eisenstein-dna-bch` and `ukkonen_1985_approximate-string-matching`/`chao_1994_linear-space` |
| C6 | REFUTED | `"accurate prior-art attribution"` is false: wrong `Gulhan 2025` DOI (`10.1186/s13059-025-03468-0` vs `10.1186/s13059-025-03830-0`), 15 `horizon_paper.md` refs missing from `references.bib`, 2 orphan `.bib` entries, and uncited `romain_2026_inversions-pangenome-graphs` |
