# Round 5 Adjudication -- prior (prior-art)

- **Attack sealed digest verified:** `da4bd66bb3a1` (`reports/rounds/round-5/attack-prior-art.md`) matches `sha256sum` prefix `da4bd66bb3a1`.
- **Prior rounds checked:** `scripts/round history "prior" --kind adjudication` returned Round 1–4 `adjudication-L6.md`. Round 4 `F1`–`F7` were all `UPHELD`, establishing the split between **free** text/citation/bibliography repairs and **non-free** external benchmark runs (`odgi`, `gfatools`, `FastGA`, `LASTZ`, `KegAlign`, `MULTIZ`, `Cactus`, `halLiftover`). Round 5 resolved several Round 4 citation omissions in `manuscript/horizon_paper.md` (including fixing `Sena et al. 2026` and adding `vg`, `odgi`, `BEDTools`, `minimap2`, `KMC 3`, `FastGA`, `KegAlign`, `Liao 2023`, `Hickey 2023`, `minigraph`, `Lunter 2008`, `last-train`, `Enredo`, `Hou 2009`, `progressiveMauve`, `SibeliaZ`, `MCScanX`) and added a `vg v1.60.0` subgraph comparison (`C1`) and `minimap2 -x asm5` `chr6` baseline (`C2`), but did not sync `manuscript/references.bib`.
- **Method cards:** `refs/methods/` does not exist in this repository; no statistical method card governs these prior-art and bibliography-hygiene findings.
- **Manuscript vs. attacked snapshot:** `diff -u reports/rounds/round-5/ATTACKED-horizon_paper.md manuscript/horizon_paper.md` and `diff -u reports/rounds/round-5/ATTACKED-horizon_companion.md manuscript/horizon_companion.md` are both empty. All cited slugs and DOIs were verified against `refs/catalog_index.json` and `manuscript/references.bib` (`harris_2007_improved` is absent from `refs/catalog_index.json` but already present in `manuscript/references.bib:79-84` and `manuscript/horizon_paper.md:16, 115`).

---

## Finding-by-Finding Rulings

### F1 -- C1: `vg` subgraph comparison is apples-to-oranges and omits `odgi` / `gfatools` baselines and pangenome graph construction prior art
- **Ruling:** **NARROWED**
- **What survives:** `manuscript/horizon_paper.md:10, 82` states that `horizon_gfa stats` (`0.08 s`, `22.6 MB` RSS) outperforms `vg stats -z -l -N` (`4.2984 s`, `77.0 MB` RSS) and `vg convert -g -p` (`1.8292 s`, `76.6 MB` RSS) without noting that `horizon_gfa stats` is a single-pass counter-only scanner (`line 47`) whereas `vg` builds in-memory handle-graph structures, and `odgi`/`gfatools` (plus external baselines for `horizon_clean`/`horizon_resolve`/`horizon_wga --graph`) remain unbenchmarked (acknowledged as pending at `line 96`).
- **Scope narrowed:** `manuscript/horizon_paper.md:46` already cites primary pangenome graph construction frameworks (`minigraph`, `Minigraph-Cactus`, `PGGB`/`odgi`, `vg`); `siren_2024_personalized-pangenome` and `kopalli_2025_toward` are peripheral to streaming GFA `stats`/`extract`, whereas the evaluation/comparison surveys `garrison_2024_building` and `andreace_2023_comparing` are relevant omissions.
- **Repair:** Qualifying `manuscript/horizon_paper.md:10, 82` to distinguish counter-only scanning from `vg` in-memory graph construction and citing `garrison_2024_building` and `andreace_2023_comparing` is **free** (endorsed); running `odgi stats` and `gfatools stat` baselines on the `10k` subgraph is **non-free**.

### F2 -- C2: `minimap2 -x asm5` outperforms Horizon on speed, exact recall, and precision while `FastGA`/`LASTZ`/`KegAlign` and `ntHash2`/`k`-nonical prior art are missing
- **Ruling:** **NARROWED**
- **What survives:** `C2` lacks a many-to-many `minimap2` mode or `FastGA`/`LASTZ`/`KegAlign` baselines (acknowledged as pending at `manuscript/horizon_paper.md:59, 96`), and `manuscript/horizon_paper.md:28` omits `kazemi_2022_nthash2`, `marcais_2024_knonical`, and `wittler_2023_canonical-kmer-encoding` on strand-symmetric spaced-seed / canonical *k*-mer hashing.
- **Scope narrowed:** `manuscript/horizon_paper.md:10, 76, 90` and `round.json:25, 28` already report `minimap2 -x asm5`'s faster runtime (`46.95 s` vs. `258.88–304.53 s`) and higher `1-to-1` exact recall/precision (`99.71%`/`99.75%` vs. `97.96%`/`48.63%`) side-by-side with Horizon's many-to-many numbers rather than concealing them; `shivakumar_2025_mumemto` is peripheral to `KP-14` spaced-seed hashing.
- **Repair:** Citing `kazemi_2022_nthash2`, `marcais_2024_knonical`, and `wittler_2023_canonical-kmer-encoding` in `manuscript/horizon_paper.md:28` and explicitly noting in `line 76` prose that `1-to-1` `minimap2 -x asm5` is faster and higher-precision is **free** (endorsed); running secondary-mode `minimap2`, `FastGA`, `LASTZ`, or `KegAlign` baselines is **non-free**.

### F3 -- C3: Multi-species `horizon_merge` / `horizon_guide` and `horizon_lift` lack external baselines and omit UCSC netting (`Kent et al. 2003`), `Mugsy`, and `Liftoff`
- **Ruling:** **UPHELD**
- **What survives:** Consistent with Round 4 `F3`, `C3` (`manuscript/horizon_paper.md:74, 82`) reports single-tool resource/count measurements for `horizon_merge`, `horizon_guide`, and `horizon_lift` without external baselines (`MULTIZ`/`Cactus`/`halLiftover`/`Liftoff`), and `manuscript/horizon_paper.md:53, 82` uses `"netting"` while omitting `kent_2003_evolutions-cauldron` (plus `shumate_2021_liftoff` and `angiuoli_2011_mugsy`).
- **Repair:** Qualifying `C3` as single-tool resource/count measurements and citing `kent_2003_evolutions-cauldron`, `shumate_2021_liftoff`, and `angiuoli_2011_mugsy` is **free** (endorsed); running external `MULTIZ`/`Cactus`/`halLiftover`/`Liftoff` baselines is **non-free**.

### F4 -- C4: Pairwise WGA benchmarks omit external aligner comparisons (`LASTZ`, `minimap2`, `FastGA`, `last-train`) and `tantan` repeat-masking prior art
- **Ruling:** **NARROWED**
- **What survives:** None of the eight Table 1 benchmarks or the `18%`-GC random DNA null include head-to-head external aligner baselines (`LASTZ`, `minimap2`, `FastGA`, `last-train`), and `manuscript/horizon_paper.md` omits `frith_2011_tantan` (directly relevant to suppressing AT-rich/simple-repeat false positives alongside `morgulis_2006_dustmasker`) and `chen_2010_comparative`.
- **Scope narrowed:** `manuscript/horizon_paper.md:59, 61` already qualifies Table 1 as reference-MAF concordance (`Recall vs Reference MAF (%)`) and explicitly states that head-to-head external pairwise aligner re-runs remain pending, so the body-text qualification leg is already largely met.
- **Repair:** Citing `frith_2011_tantan` and `chen_2010_comparative` (and echoing reference-MAF concordance in `line 10`) is **free** (endorsed); running external `LASTZ`/`minimap2`/`FastGA`/`last-train` baselines on Table 1 slices and the `18%`-GC null is **non-free**.

### F5 -- C5: `horizon_companion.md` leaves two `# 5. References` entries uncited in its body and omits prior art for Eisenstein DNA codes, banded DP, and Karlin–Altschul scaling
- **Ruling:** **UPHELD**
- **What survives:** Verified in `manuscript/horizon_companion.md` that `Allison et al. (1992)` (`line 189`) and `Foster et al. (2007)` (`line 191`) are listed in `# 5. References` with zero body citations, Section 3.3 (`line 131`) invokes `"Continuous Karlin–Altschul search-space scaling"` uncited and omits `Rissanen (1978)`/`Grünwald (2007)` and banded DP prior art (`ukkonen_1985_approximate-string-matching`, `chao_1994_linear-space`), and Section 3.1 (`lines 87-96`) omits `sajjad_2026_eisenstein-dna-bch` (related algebraic DNA encoding over the Eisenstein field).
- **Repair:** Citing `Allison et al. (1992)` in Section 3.3, `Foster et al. (2007)` in Section 2.2, `Rissanen (1978)`/`Grünwald (2007)` and `ukkonen_1985_approximate-string-matching`/`chao_1994_linear-space` in Section 3.3, and `sajjad_2026_eisenstein-dna-bch` in Section 3.1 (plus adding Karlin & Altschul 1990 if retained in prose) is **free** (endorsed).

### F6 -- C6: `C6`'s claim of "accurate prior-art attribution" is refuted by a wrong DOI for `Gulhan et al. (2025)`, 15 `horizon_paper.md` references missing from `manuscript/references.bib`, and uncited companion prior art
- **Ruling:** **UPHELD**
- **What survives:** All four counts refuting `C6`'s `"accurate prior-art attribution"` assertion (`round.json:70`) are verified: (1) `manuscript/horizon_paper.md:113` gives `Gulhan et al. (2025)` DOI `10.1186/s13059-025-03468-0` instead of `10.1186/s13059-025-03830-0` (`refs/catalog_index.json`); (2) `manuscript/horizon_paper.md` `# References` has `38` entries while `manuscript/references.bib` has `25` entries, omitting all `15` entries added after Round 4; (3) orphan `.bib` entries `jiang_2008_universal-dna-representation` and `miller_2020_population-structure-clustering` remain in `manuscript/references.bib:153-162, 186-195` despite Round 4 `F7`; and (4) `romain_2026_inversions-pangenome-graphs` and duplication-aware WGA references (`Paten 2008`, `Hou 2009`, `Darling 2010`, `Minkin & Medvedev 2020`, `Wang 2012`) are omitted from `horizon_paper.md` and `horizon_companion.md:157-162`, respectively.
- **Repair:** Fixing `manuscript/horizon_paper.md:113` to DOI `10.1186/s13059-025-03830-0`, adding all `15` missing entries to `manuscript/references.bib`, removing the `2` orphan `.bib` entries, and citing `romain_2026_inversions-pangenome-graphs` and the duplication-aware WGA references in `manuscript/horizon_companion.md:157-162` is **free** (endorsed).

### F7 -- Meta-finding: Runner-authored claim list steers away from missing external baselines and bibliography drift
- **Ruling:** **NARROWED**
- **What survives:** `round.json:119-120` confirms `runner_is_author: true` (`claim_list_author: "runner"`), and the framing critique holds for `C1` (highlighting `vg` without noting counter-only vs. in-memory construction or missing `odgi`/`gfatools`), `C3`/`C4` (single-tool numbers without external baselines), and `C6` (asserting `"accurate prior-art attribution"` despite a wrong DOI and `15` missing `.bib` entries).
- **Scope narrowed:** Overturned on the `C2` burying leg—`round.json:25, 28` explicitly quotes `minimap2 -x asm5`'s `46.95 s` runtime and `99.71%`/`99.75%` exact recall/precision alongside Horizon's numbers.

---

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | Counter-only vs. `vg` in-memory qualification and `odgi`/`gfatools` baselines missing (text + `garrison_2024_building`/`andreace_2023_comparing` citations free; runs non-free); `line 46` already cites primary graph tools |
| F2 | NARROWED | `minimap2 -x asm5` numbers are already reported side-by-side; `kazemi_2022_nthash2`/`marcais_2024_knonical`/`wittler_2023_canonical-kmer-encoding` citations free, `FastGA`/`LASTZ`/`KegAlign`/secondary `minimap2` runs non-free |
| F3 | UPHELD | Single-tool `C3` timings lack external baselines (non-free) and omit `kent_2003_evolutions-cauldron`, `shumate_2021_liftoff`, and `angiuoli_2011_mugsy` (text/citations free) |
| F4 | NARROWED | `horizon_paper.md:59, 61` already qualifies Table 1 as reference-MAF concordance; `frith_2011_tantan`/`chen_2010_comparative` citations free, external aligner runs non-free |
| F5 | UPHELD | Verified `Allison 1992` and `Foster 2007` uncited in `horizon_companion.md` body plus omitted `sajjad_2026_eisenstein-dna-bch` and `ukkonen_1985_approximate-string-matching`/`chao_1994_linear-space`; fix is free |
| F6 | UPHELD | Refutation of `C6`'s `"accurate prior-art attribution"` holds: wrong `Gulhan 2025` DOI, 15 refs missing from `references.bib`, 2 orphan `.bib` entries, and uncited `romain_2026_inversions-pangenome-graphs`; fix is free |
| F7 | NARROWED | `runner_is_author: true` framing critique holds for `C1`, `C3`, `C4`, and `C6`, but `C2` (`round.json:25, 28`) explicitly quotes `minimap2 -x asm5`'s faster runtime and higher 1-to-1 precision |
