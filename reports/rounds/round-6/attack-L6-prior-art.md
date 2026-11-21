# Round 6 Attack -- L6 (prior-art)

- **Pinned revision & snapshots:** `42700e142b70996c23cdf391ac29e1e5bcf4640d`; attacked frozen snapshots in `reports/rounds/round-6/ATTACKED-*.md`, `refs/catalog_index.json` (`301` entries, pinned at `refs/LIBRARY_VERSION` `d4df9122eddd84694b2c3215e5478c46f4314044`), and `manuscript/references.bib` (`40` entries).
- **Library availability (Rule 2):** `../library` is not present on disk (`ls: cannot access '/home/rcb112/project/gemini/bench/library': No such file or directory`). Per `docs/round_rubric.md`, this attack works from `refs/catalog_index.json` and `manuscript/references.bib` alone and does not assert paper contents from memory.
- **Method cards (Rule 3):** This project keeps **no method cards** (`refs/methods/` does not exist).
- **Prior rounds checked (Rule 4):** `scripts/round history "prior-art" --kind adjudication` and `reports/rounds/round-5/adjudication-L6.md` were reviewed. Round 5 `F1`–`F7` established the split between **free** text/citation/bibliography repairs and **non-free** external benchmark runs (`odgi`, `gfatools`, `FastGA`, `LASTZ`, `KegAlign`, `MULTIZ`, `Cactus`, `halLiftover`). Round 6 (`1ddc872` / `42700e1`) applied several Round 5 L6 repairs (`Gulhan et al. 2025` DOI fixed to `10.1186/s13059-025-03830-0`, `15` missing `horizon_paper.md` entries added to `manuscript/references.bib`, `2` Round 4/5 orphan `.bib` entries removed, `Kent et al. 2003` added to `horizon_paper.md:55, 122`, and `Allison et al. 1992` cited in `horizon_companion.md:117`), while leaving other Round 5 L6-endorsed citations unapplied and introducing new P0–P4 prior-art gaps in `C1`–`C3`.

---

## F1 -- C1: P0–P4 suite (`horizon_filter`, `horizon_locus`, `horizon_resolve`) omits held low-complexity/inversion prior art (`longdust`, `tantan`, `Romain et al. 2026`) and `impg`/`wfmash` citations

- **Claim:** `C1` (`docs/PANGENOME_SUITE.md`, `reports/rounds/round-6/ATTACKED-PANGENOME_SUITE.md:3`).
- **Verdict:** **NARROWED**
- **Prior-round context:** Discloses Round 5 L6 `F1` (`NARROWED`), which found that `manuscript/horizon_paper.md` lacked `garrison_2024_building` and `andreace_2023_comparing` as well as `odgi`/`gfatools` baselines on the `10k` HPRC subgraph. Both citations remain absent in `reports/rounds/round-6/ATTACKED-horizon_paper.md`, and `odgi`/`gfatools` remain pending (`ATTACKED-horizon_paper.md:84, 98`).
- **Concrete flaw (new in Round 6):**
  1. `C1` introduces `horizon_filter` (streaming PAF SDUST+BOCPD+MDL filter) and `horizon_locus` (regional `impg`-style multi-FASTA locus aligner), and `ATTACKED-horizon_paper.md:28, 52-53, 84` cites `Morgulis et al. (2006)` (`morgulis_2006_dustmasker`) for SDUST repeat masking while omitting two held low-complexity and tandem-repeat masking papers in `refs/catalog_index.json`: `li_2026_finding` (*Finding low-complexity DNA sequences with longdust*, `doi:10.1093/bioinformatics/btag112`) and `frith_2011_tantan` (*A new repeat-masking method enables specific detection of homologous sequences*, `doi:10.1093/nar/gkq1212`).
  2. `ATTACKED-horizon_paper.md:10, 50, 52-53, 84` describes `horizon_locus` as `"impg-style"`, `horizon_filter` as filtering `wfmash` PAFs, and `horizon_resolve` as using `PackedGraph` + $O(|V| + |E|)$ Tarjan biedged snarl decomposition, without citing `impg` or `wfmash` (both absent from `manuscript/references.bib` and `refs/catalog_index.json`), without tying `Paten et al. (2018)` / `Sena et al. (2026)` (`line 40`) to `horizon_resolve` (`line 50`), and without citing `romain_2026_inversions-pangenome-graphs` (`doi:10.1093/nargab/lqag125`, held in `refs/catalog_index.json` and already endorsed in Round 5 L6 `F6`).
- **Scope:** The Rust implementation inventory and `118` passing unit tests in `C1` are outside the prior-art lens, and `ATTACKED-horizon_paper.md:40, 46` already cites primary pangenome graph frameworks (`minigraph`, `Minigraph-Cactus`, `PGGB`/`odgi`, `vg`) and `Paten et al. (2018)` / `Sena et al. (2026)`.
- **Correction:** Cite `li_2026_finding`, `frith_2011_tantan`, `garrison_2024_building`, `andreace_2023_comparing`, and `romain_2026_inversions-pangenome-graphs` in `manuscript/horizon_paper.md` and `manuscript/references.bib`, and add citations for `impg` and `wfmash`.

---

## F2 -- C2: `horizon_filter` benchmarks omit head-to-head `SweepGA` / PAF-netter baselines and tandem-repeat masking prior art

- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md`, `reports/rounds/round-6/ATTACKED-BENCHMARK_RESULTS.md:85-88`).
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. `C2` evaluates `horizon_filter --net` only against raw `FastGA` (`56,857` `chr6` records; `36,601` Koala × Amphioxus records) and VGP's pre-filtered `FastGA_11`, with no head-to-head comparison against `SweepGA` (explicitly named in `ATTACKED-VGP_INTEGRATION_PLAN.md:12` as VGP's current PAF filter) or other PAF netting pipelines.
  2. `C2` highlights flagging `FastGA_11`'s #1 longest hit (`3,548 bp`) as a `GA`/`AT`/`CA` simple tandem repeat (`md:f:-210.51 hb:i:1`) using `SDUST + BOCPD + MDL`, yet `ATTACKED-horizon_paper.md` omits `frith_2011_tantan` (tandem-repeat masking) and `li_2026_finding` (`longdust`), both held in `refs/catalog_index.json`, and neither `SweepGA` nor `wfmash` is in `manuscript/references.bib` or `refs/catalog_index.json`.
- **Scope:** `C2`'s `0.47s`/`4.73s` (`chr6`) and `0.10s` (Koala × Amphioxus) numbers are single-tool / raw-`FastGA` measurements rather than a claimed superiority over `SweepGA`.
- **Correction:** Qualify `C2` in `docs/BENCHMARK_RESULTS.md` and `manuscript/horizon_paper.md:84` as lacking a head-to-head `SweepGA` baseline (or run `SweepGA` on identical inputs), and cite `frith_2011_tantan` and `li_2026_finding`.

---

## F3 -- C3: `"impg-style"` `horizon_locus` benchmarks omit `impg` / MSA baselines and `impg` / `JC69` citations

- **Claim:** `C3` (`README.md`, `reports/rounds/round-6/ATTACKED-README.md:126`).
- **Verdict:** **NARROWED**
- **Concrete flaw:**
  1. `C3` (and `ATTACKED-horizon_paper.md:10, 53, 84`) frames `horizon_locus` as an `"impg-style"` multi-FASTA locus aligner, and `ATTACKED-VGP_INTEGRATION_PLAN.md:13, 47` notes that VGP runs `impg query` and `impg similarity` on `25 kb` windows. Yet `C3` reports single-tool runtimes (`0.25s`, `0.86s`, `0.32s`), RSS, and similarity metrics (`jaccard_cov`, `p_dist`, `jc69_dist`) with zero head-to-head comparison against `impg` (or `MAFFT`/`progressiveMauve`/`wfmash`) on the identical loci.
  2. `impg` is never cited in `manuscript/horizon_paper.md` or `manuscript/references.bib` (and is absent from `refs/catalog_index.json`), and `jc69_dist` (Jukes–Cantor 1969 distance) is reported in `C3` and `ATTACKED-horizon_paper.md:84` without any citation in `manuscript/horizon_paper.md` or `manuscript/references.bib`.
- **Scope:** `C3` reports single-tool resource and distance measurements without asserting a speedup or accuracy claim over `impg`.
- **Correction:** Cite `impg` and Jukes & Cantor (1969) in `manuscript/horizon_paper.md` and `manuscript/references.bib`, and qualify `C3` as single-tool measurements pending a head-to-head `impg` / MSA comparison.

---

## F4 -- C4: Core WGA benchmarks (`8cd5db0`/`2896229`) still omit Round 5 L6-endorsed spaced-seed, repeat-masking, and WGA-assessment citations

- **Claim:** `C4` (`manuscript/horizon_paper.md`, `reports/rounds/round-6/ATTACKED-horizon_paper.md:61-78`).
- **Verdict:** **NARROWED**
- **Prior-round context:** Discloses Round 5 L6 `F2` (`NARROWED`) and `F4` (`NARROWED`). Round 5 already established that `manuscript/horizon_paper.md:10, 78, 92` reports `minimap2 -x asm5` side-by-side and `lines 59, 61` qualify Table 1 as reference-MAF concordance with external aligner runs (`FastGA`, `LASTZ`, `KegAlign`) pending.
- **Concrete flaw (remaining from Round 5):** None of the held-catalog citations endorsed in Round 5 L6 `F2` and `F4` were added in Round 6 (`42700e1`):
  - `kazemi_2022_nthash2` (`doi:10.1093/bioinformatics/btac564`), `marcais_2024_knonical` (`doi:10.1093/bioinformatics/btae629`), and `wittler_2023_canonical-kmer-encoding` (`doi:10.24072/pcjournal.323`) remain omitted from `ATTACKED-horizon_paper.md:28` on canonical spaced-seed hashing;
  - `frith_2011_tantan` remains omitted on AT-rich/simple-repeat masking (`ATTACKED-horizon_paper.md:28, 61`);
  - `chen_2010_comparative` (`doi:10.1038/nbt.1637`) remains omitted on comparative WGA assessment (`ATTACKED-horizon_paper.md:61`).
- **Scope:** The `8cd5db0`/`2896229` empirical numbers and `minimap2 -x asm5` comparison are already accurately qualified in `ATTACKED-horizon_paper.md:59-61, 78, 98`.
- **Correction:** Add `kazemi_2022_nthash2`, `marcais_2024_knonical`, `wittler_2023_canonical-kmer-encoding`, `frith_2011_tantan`, and `chen_2010_comparative` to `manuscript/horizon_paper.md` and `manuscript/references.bib`.

---

## F5 -- C5 (`R5-BIB`): `manuscript/references.bib` retains orphan `foster2007combinators`, while `horizon_companion.md` and `horizon_paper.md` leave Round 5 L6-endorsed citations unapplied

- **Claim:** `C5` (`ARCHITECTURE.md` et al., carrying `R5-BIB`, `reports/rounds/round-6/round.json:60-72`).
- **Verdict:** **NARROWED**
- **Prior-round context:** Discloses Round 5 L6 `F5` (`UPHELD`) and `F6` (`UPHELD`). Round 6 (`1ddc872`) partially resolved `R5-BIB` and Round 5 L6 `F5`/`F6`: it fixed the `Gulhan et al. (2025)` DOI (`ATTACKED-horizon_paper.md:115`), added the `15` missing `horizon_paper.md` entries to `manuscript/references.bib`, removed `jiang_2008_universal-dna-representation` and `miller_2020_population-structure-clustering`, cited `Allison et al. (1992)` in `ATTACKED-horizon_companion.md:117`, and added `Kent et al. (2003)` to `ATTACKED-horizon_paper.md:55, 122`.
- **Concrete flaw (remaining from Round 5):**
  1. `manuscript/references.bib` (`40` entries) and `ATTACKED-horizon_paper.md` `# References` (`39` entries) remain out of sync: `foster2007combinators` (`manuscript/references.bib:68-77`, `doi:10.1145/1232420.1232424`) is an orphan `.bib` entry absent from `horizon_paper.md`, while `ATTACKED-horizon_companion.md:55` explicitly discusses bidirectional lens combinators (`get`/`putback`) **without citing `Foster et al. (2007)`** or listing it in `# 5. References` (`lines 181-185`) — leaving Round 5 L6 `F5`'s endorsed repair unapplied.
  2. `ATTACKED-horizon_companion.md` also leaves the other Round 5 L6 `F5`/`F6` endorsed citations unapplied: `Rissanen (1978)` / `Grünwald (2007)`, `ukkonen_1985_approximate-string-matching`, `chao_1994_linear-space`, and a Karlin–Altschul citation in Section 3.3 (`lines 115-138`, where `line 124` still invokes `"Continuous Karlin–Altschul search-space scaling"` uncited); `sajjad_2026_eisenstein-dna-bch` in Section 3.1 (`lines 86-98`); and duplication-aware WGA references (`Paten 2008`, `Hou 2009`, `Darling 2010`, `Minkin & Medvedev 2020`, `Wang 2012`) in Section 3.5 (`lines 150-156`).
  3. `ATTACKED-horizon_paper.md` still omits `10` of the `11` Round 5 L6-endorsed held-catalog citations (`garrison_2024_building`, `andreace_2023_comparing`, `kazemi_2022_nthash2`, `marcais_2024_knonical`, `wittler_2023_canonical-kmer-encoding`, `shumate_2021_liftoff`, `angiuoli_2011_mugsy`, `frith_2011_tantan`, `chen_2010_comparative`, `romain_2026_inversions-pangenome-graphs`).
- **Scope:** Code-to-doc consistency in `ARCHITECTURE.md`, `README.md`, and `docs/manual.md` is outside the prior-art lens, and the primary `R5-BIB` diff (`15` `.bib` additions, `2` `.bib` removals, `Gulhan 2025` DOI fix) landed.
- **Correction:** Cite `Foster et al. (2007)` in `manuscript/horizon_companion.md:55` (and `# 5. References`) or remove `foster2007combinators` from `manuscript/references.bib`, and apply the remaining Round 5 L6-endorsed citations in `horizon_companion.md` and `horizon_paper.md`.

---

## F6 -- Meta-finding: Runner-authored claim list steers away from missing `SweepGA` / `impg` baselines and unapplied Round 5 L6 citations

- **Verdict:** **NARROWED** (claim-list framing critique)
- **Prior-round context:** Discloses Round 5 L6 `F7` (`NARROWED`), which flagged `"runner_is_author": true` / `"claim_list_author": "runner"` for steering away from missing external baselines and bibliography drift.
- **Evidence:** `reports/rounds/round-6/round.json:144-145` again records `"runner_is_author": true` and `"claim_list_author": "runner"`. `C2` and `C3` highlight single-tool timings and raw `FastGA`/`FastGA_11` overlap while omitting `SweepGA` and `impg` (the exact tools named in `ATTACKED-VGP_INTEGRATION_PLAN.md:12-13`), and `C5` carries `R5-BIB` without noting that `foster2007combinators` remains orphaned and `10` Round 5 L6-endorsed `horizon_paper.md` citations plus `horizon_companion.md` `F5`/`F6` repairs remain unapplied.

---

## Verdict Table

| Claim | Verdict | Summary |
|---|---|---|
| C1 | NARROWED | `horizon_filter`/`horizon_locus`/`horizon_resolve` omit held `li_2026_finding`, `frith_2011_tantan`, `garrison_2024_building`, `andreace_2023_comparing`, `romain_2026_inversions-pangenome-graphs`, and `impg`/`wfmash` citations |
| C2 | NARROWED | `horizon_filter` benchmarks omit head-to-head `SweepGA` baseline and held `frith_2011_tantan`/`li_2026_finding` repeat-masking prior art |
| C3 | NARROWED | `"impg-style"` `horizon_locus` benchmarks report single-tool numbers without `impg`/MSA baselines or `impg`/`JC69` citations |
| C4 | NARROWED | `8cd5db0`/`2896229` numbers and `minimap2` are qualified, but Round 5 L6-endorsed citations (`ntHash2`, `k`-nonical, `Wittler`, `tantan`, `Chen & Tompa`) remain omitted |
| C5 | NARROWED | `15` `.bib` entries and `Gulhan 2025` DOI fixed, but `foster2007combinators` remains an orphan `.bib` entry uncited in `horizon_companion.md:55` and Round 5 L6 `F5`/`F6` repairs remain incomplete |
