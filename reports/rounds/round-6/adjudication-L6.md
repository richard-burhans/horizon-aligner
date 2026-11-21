# Round 6 Adjudication -- L6 (prior-art)

- **Attack sealed digest verified:** `9e2d00896b2f` (`reports/rounds/round-6/attack-L6-prior-art.md` `sha256sum` = `9e2d00896b2f3b835084d65f304b5478093299d091094e4482e8900494679fc9`).
- **Prior rounds checked:** `scripts/round history "prior-art" --kind adjudication` and `reports/rounds/round-5/adjudication-L6.md` (`F1`–`F7`) were reviewed. Round 5 established the split between **free** text/citation/bibliography repairs and **non-free** external benchmark runs (`odgi`, `gfatools`, `FastGA`, `LASTZ`, `KegAlign`, `MULTIZ`, `Cactus`, `halLiftover`, plus Round 6's `SweepGA`, `impg`, `MAFFT`/`progressiveMauve`/`wfmash`).
- **Method cards:** `refs/methods/` does not exist in this repository; no statistical method card governs these prior-art and bibliography-hygiene findings.
- **Snapshots vs. live files & catalog verification:** `diff -u` is empty across all attacked snapshots (`ATTACKED-horizon_paper.md`, `ATTACKED-horizon_companion.md`, `ATTACKED-PANGENOME_SUITE.md`, `ATTACKED-BENCHMARK_RESULTS.md`, `ATTACKED-README.md`, `ATTACKED-ARCHITECTURE.md` vs. their live files). `../library` is absent; all cited slugs and DOIs were independently verified against `refs/catalog_index.json` (`301` entries) and `manuscript/references.bib` (`40` entries).

---

## Finding-by-Finding Rulings

### F1 -- C1: P0–P4 suite (`horizon_filter`, `horizon_locus`, `horizon_resolve`) omits held low-complexity/inversion prior art (`longdust`, `tantan`, `Romain et al. 2026`) and `impg`/`wfmash` citations
- **Ruling:** **NARROWED**
- **What survives:** `docs/PANGENOME_SUITE.md:3` (`C1`) is a Rust implementation and unit-test inventory (`118` passing tests, `11` binaries + `bench_chr1`), and `manuscript/horizon_paper.md:40, 46` already cites primary pangenome graph frameworks (`minigraph`, `Minigraph-Cactus`, `PGGB`/`odgi`, `vg`) and `Paten et al. (2018)` / `Sena et al. (2026)`. However, `horizon_paper.md:28, 52-53, 84` cites `Morgulis et al. (2006)` for SDUST repeat masking while omitting `li_2026_finding` (`doi:10.1093/bioinformatics/btag112`) and `frith_2011_tantan` (`doi:10.1093/nar/gkq1212`), both verified in `refs/catalog_index.json`; `horizon_paper.md:10, 50, 52-53, 84` also omits `impg` and `wfmash` citations (neither is in `catalog_index.json` or `references.bib`), does not tie `Paten et al. (2018)` / `Sena et al. (2026)` (`line 40`) to `horizon_resolve` (`line 50`), and leaves `romain_2026_inversions-pangenome-graphs` (`doi:10.1093/nargab/lqag125`), `garrison_2024_building`, and `andreace_2023_comparing` (endorsed in Round 5 L6 `F1`, `F6`) omitted.
- **Repair:** Adding `li_2026_finding`, `frith_2011_tantan`, `garrison_2024_building`, `andreace_2023_comparing`, `romain_2026_inversions-pangenome-graphs`, `impg`, and `wfmash` to `manuscript/horizon_paper.md` and `manuscript/references.bib` (and linking `Paten 2018`/`Sena 2026` at `line 50`) is **free** (endorsed; bibliographic details for `impg`/`wfmash` must be verified from primary sources since neither is in `refs/catalog_index.json`).

### F2 -- C2: `horizon_filter` benchmarks omit head-to-head `SweepGA` / PAF-netter baselines and tandem-repeat masking prior art
- **Ruling:** **NARROWED**
- **What survives:** Verified in `docs/BENCHMARK_RESULTS.md:85-88` and `manuscript/horizon_paper.md:84` that `horizon_filter --net` is evaluated only against raw `FastGA` and VGP's pre-filtered `FastGA_11`, with no head-to-head `SweepGA` baseline (`docs/VGP_INTEGRATION_PLAN.md:12` names `SweepGA` as VGP's current PAF filter), and `horizon_paper.md:84` flags `FastGA_11`'s `3,548 bp` simple tandem repeat while omitting `frith_2011_tantan` and `li_2026_finding` (both verified in `refs/catalog_index.json`). As narrowed, `C2` reports single-tool / raw-`FastGA` measurements without claiming superiority over `SweepGA`.
- **Repair:** Qualifying `C2` in `docs/BENCHMARK_RESULTS.md` and `manuscript/horizon_paper.md:84` as lacking a head-to-head `SweepGA` baseline and citing `frith_2011_tantan` and `li_2026_finding` is **free** (endorsed); running `SweepGA` on identical inputs is **non-free**.

### F3 -- C3: `"impg-style"` `horizon_locus` benchmarks omit `impg` / MSA baselines and `impg` / `JC69` citations
- **Ruling:** **NARROWED**
- **What survives:** Verified that `README.md:126` and `manuscript/horizon_paper.md:10, 53, 84` frame `horizon_locus` as `"impg-style"` and report `jc69_dist` without citing `impg` or Jukes & Cantor (1969) in `horizon_paper.md` or `references.bib` (both absent from `catalog_index.json`), and without head-to-head `impg` or MSA (`MAFFT`/`progressiveMauve`/`wfmash`) baselines on the identical loci (`docs/VGP_INTEGRATION_PLAN.md:13, 47`). As narrowed, `C3` reports single-tool resource and distance measurements without asserting superiority over `impg`.
- **Repair:** Citing `impg` and Jukes & Cantor (1969) in `manuscript/horizon_paper.md` and `manuscript/references.bib` and qualifying `C3` as single-tool measurements pending a head-to-head `impg`/MSA comparison is **free** (endorsed); running external `impg`/MSA baselines is **non-free**.

### F4 -- C4: Core WGA benchmarks (`8cd5db0`/`2896229`) still omit Round 5 L6-endorsed spaced-seed, repeat-masking, and WGA-assessment citations
- **Ruling:** **NARROWED**
- **What survives:** Verified against `manuscript/horizon_paper.md`, `manuscript/references.bib`, and `refs/catalog_index.json` that all five Round 5 L6-endorsed held-catalog citations (`kazemi_2022_nthash2`, `marcais_2024_knonical`, `wittler_2023_canonical-kmer-encoding` for `line 28`; `frith_2011_tantan` for `lines 28, 61`; `chen_2010_comparative` for `line 61`) remain omitted. As narrowed, the `8cd5db0`/`2896229` empirical numbers and `minimap2 -x asm5` comparison are already accurately qualified in `horizon_paper.md:59-61, 78, 98`.
- **Repair:** Adding `kazemi_2022_nthash2`, `marcais_2024_knonical`, `wittler_2023_canonical-kmer-encoding`, `frith_2011_tantan`, and `chen_2010_comparative` to `manuscript/horizon_paper.md` and `manuscript/references.bib` is **free** (endorsed).

### F5 -- C5 (`R5-BIB`): `manuscript/references.bib` retains orphan `foster2007combinators`, while `horizon_companion.md` and `horizon_paper.md` leave Round 5 L6-endorsed citations unapplied
- **Ruling:** **NARROWED**
- **What survives:** While the primary `R5-BIB` sync (`15` `.bib` additions, `2` orphan `.bib` removals, `Gulhan 2025` DOI fix, `Allison et al. 1992` body citation in `horizon_companion.md:117`, `Kent et al. 2003` in `horizon_paper.md:55, 122`) landed, verified that: (1) `manuscript/references.bib` has `40` entries vs. `39` in `horizon_paper.md` `# References`, leaving `foster2007combinators` (`references.bib:68-77`) orphaned while `horizon_companion.md:55` explicitly discusses `"formal bidirectional lens combinators such as get/putback"` without citing `Foster et al. (2007)` or listing it in `# 5. References` (`lines 181-185`); (2) `horizon_companion.md` still omits Round 5 L6 `F5`/`F6`-endorsed citations (`Rissanen 1978`, `Grünwald 2007`, `ukkonen_1985_approximate-string-matching`, `chao_1994_linear-space`, and Karlin–Altschul at `line 124`; `sajjad_2026_eisenstein-dna-bch` in Section 3.1; `Paten 2008`, `Hou 2009`, `Darling 2010`, `Minkin & Medvedev 2020`, `Wang 2012` in Section 3.5); and (3) `horizon_paper.md` still omits `10` of the `11` Round 5 L6-endorsed held-catalog citations.
- **Repair:** Citing `Foster et al. (2007)` in `manuscript/horizon_companion.md:55` (and `# 5. References`) or removing `foster2007combinators` from `manuscript/references.bib`, and applying the remaining Round 5 L6-endorsed citations in `horizon_companion.md` and `horizon_paper.md` is **free** (endorsed).

### F6 -- Meta-finding: Runner-authored claim list steers away from missing `SweepGA` / `impg` baselines and unapplied Round 5 L6 citations
- **Ruling:** **NARROWED**
- **What survives:** `reports/rounds/round-6/round.json:152-153` (cited as `144-145` in the attack, which are `"sealed"` entries) confirms `"runner_is_author": true` and `"claim_list_author": "runner"`. The framing critique holds: `C2` and `C3` foreground single-tool timings and raw `FastGA`/`FastGA_11` overlap while omitting `SweepGA` and `impg` (`docs/VGP_INTEGRATION_PLAN.md:12-13`), and `C5` carries `R5-BIB` without noting the remaining `foster2007combinators` orphan or the unapplied Round 5 L6-endorsed citations.

---

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | `C1` is an implementation inventory and `line 46` cites primary graph tools, but `horizon_paper.md` omits `li_2026_finding`, `frith_2011_tantan`, `garrison_2024_building`, `andreace_2023_comparing`, `romain_2026_inversions-pangenome-graphs`, and `impg`/`wfmash` (citations free) |
| F2 | NARROWED | `C2` makes no superiority claim over `SweepGA`, but lacks a `SweepGA` baseline (run non-free; qualification free) and omits held `frith_2011_tantan`/`li_2026_finding` (citations free) |
| F3 | NARROWED | `C3` makes no superiority claim over `impg`, but reports single-tool `"impg-style"`/`jc69_dist` numbers without `impg`/MSA baselines (runs non-free; qualification and `impg`/JC69 citations free) |
| F4 | NARROWED | `horizon_paper.md:59-61, 78, 98` already qualifies `8cd5db0`/`2896229` and `minimap2`, but 5 Round 5 L6-endorsed held-catalog citations (`ntHash2`, `k`-nonical, `Wittler`, `tantan`, `Chen & Tompa`) remain omitted (citations free) |
| F5 | NARROWED | Primary `R5-BIB` sync landed, but `foster2007combinators` remains an orphan `.bib` entry uncited in `horizon_companion.md:55` and Round 5 L6 `F5`/`F6` citations remain incomplete (fix is free) |
| F6 | NARROWED | `round.json:152-153` confirms `runner_is_author: true`; claim framing omits `SweepGA`/`impg` baselines in `C2`/`C3` and unapplied Round 5 L6 repairs in `C5` |
