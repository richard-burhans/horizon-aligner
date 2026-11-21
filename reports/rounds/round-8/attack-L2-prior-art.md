# Round 8 Attack — `prior-art` Lens

**Target:** `docs/COMMUNITY_UTILITIES_PROPOSAL.md` (frozen snapshot `reports/rounds/round-8/ATTACKED-COMMUNITY_UTILITIES_PROPOSAL.md`, commit `0f31c26cf24e52a865cdcd9963abe9e6ec7ad499`).  
**Literature & standards provenance:**
- Checked `refs/catalog_index.json` and `manuscript/references.bib`, and read held paper summaries at the pinned `refs/LIBRARY_VERSION` revision `d4df9122eddd84694b2c3215e5478c46f4314044` in `/home/rcb112/project/claude/bench/repos/library`.
- Checked `scripts/round history` for `COMMUNITY_UTILITIES_PROPOSAL`, `horizon hal`, and `horizon gfa`: no prior mentions in closed rounds 1–7.
- Per the project map, this project keeps **no method cards** yet.

---

## F1 -- `horizon hal` (C1) 1:1 duplicates `halTools` and `hal2vg`; only zero-`libhdf5` static packaging is non-duplicative

- **Claim ID:** `C1`
- **Verdict:** `NARROWED`
- **Concrete flaw:** All four proposed subcommands (`stats`, `to-maf`, `liftover`, `to-gfa`) replicate existing `halTools` and `hal2vg` functionality without adding a new algorithmic capability, and the claim that deploying `halTools` in containers or HPC environments is "brittle" ignores existing prebuilt Bioconda and Cactus container distributions.
- **Evidence:**
  - `hickey_2013_hierarchical` (Table 1, summary citing lines 126–138 and 188–204) already ships `halStats`, `hal2maf`, and `halLiftover` in `ComparativeGenomicsToolkit/hal`.
  - `romain_2026_inversions-pangenome-graphs` (lines 465–468) uses `hal2vg --noAncestors` + `vg convert -f -W` for HAL-to-GFA conversion.
  - `hal`, `hal2vg`, and `comparative-annotation-toolkit` are distributed as prebuilt Bioconda packages and in official Cactus Docker/Singularity images (`quay.io/comparative-genomics-toolkit/cactus`).
  - For reference-free coordinate liftover outside `.hal` files, `pangenome/impg` already performs fast PAF/1aln/TPA coordinate liftover without `libhdf5`.
- **Correction:** Narrow C1 in `docs/COMMUNITY_UTILITIES_PROPOSAL.md` (`### 1. horizon hal`) from a missing-capability or container-deployment gap to a **zero-`libhdf5` static-binary / embedded-Rust packaging convenience** for environments where Bioconda or Cactus containers are unavailable.

---

## F2 -- `horizon gfa-qc` (C2) duplicates `odgi`, `panacus`, and `gfatools`, and misrepresents `kopalli_2025_toward`'s `minimap2`+`Liftoff` glue

- **Claim ID:** `C2`
- **Verdict:** `REFUTED`
- **Concrete flaw:** The proposal claims `horizon gfa-qc` fills a retention, clipping, component-bridge, and node-depth/AF diagnostic gap not met without `odgi paths` + `minimap2` + `samtools` + `Liftoff` glue. In reality: (1) `odgi`, `panacus`, and `gfatools` already compute node depth, AF/coverage histograms, connected components, and path lengths without `minimap2`/`samtools`/`Liftoff`; and (2) `kopalli_2025_toward` used `minimap2`+`samtools`+`Liftoff` for *base-level assembly-to-path re-alignment* and *gene-loss annotation*, neither of which `horizon gfa-qc`'s `--fasta-manifest` length comparison replaces.
- **Evidence:**
  - `guarracino_2022_odgi` (summary lines 254–303, 451, 229–230) documents `odgi stats` (nodes, edges, paths, weakly connected components `-W`, MultiQC YAML), `odgi depth` (per-node path depth up to `10^5–10^6`-fold), and `odgi paths` (per-path retained bp).
  - `panacus` (explicitly run in `kopalli_2025_toward`, summary citing line 227) computes GFA coverage-depth histograms, quorum/AF tables, and pangenome growth directly from blunt GFA (`panacus histgrowth`, `panacus report`). `gfatools stat` reports GFA segment/link/component stats.
  - `kopalli_2025_toward` (lines 298–301) extracted paths with `odgi paths`, aligned assemblies back to their paths with `minimap2`+`samtools` (a metric `kopalli_2025_toward` itself notes is self-referential), and lifted gene annotations with `Liftoff` v1.6.3 to count missing genes. Comparing path lengths against `.fai` lengths in `horizon gfa-qc` cannot detect missing genes or sequence divergence.
- **Correction:** Refute C2 as an unmet analytical gap. At most, `horizon gfa-qc` could be reframed as a **single-pass, ~97 MB RSS convenience wrapper** that avoids `odgi build`'s `.og` conversion memory (`10 GB` peak on chr6 in `guarracino_2022_odgi`, line 660).

---

## F3 -- `horizon inversions` (C3) reimplements `INVPG_annot` and cites SPQR trees (`sena_2026`) while proposing Tarjan snarls instead of `BubbleFinder`

- **Claim ID:** `C3`
- **Verdict:** `NARROWED`
- **Concrete flaw:** `INVPG_annot` (`romain_2026_inversions-pangenome-graphs`) already implements the exact two-class (`path_explicit` + `alignment_rescued`) algorithm with the same defaults (`--min-len 50`, `--min-cov 0.5`), and `BubbleFinder` (`algbio/BubbleFinder`, companion to `sena_2026_bubble-subgraphs-spqr`) already implements linear-time SPQR-tree snarl/ultrabubble enumeration. Moreover, `horizon inversions` cites `sena_2026_bubble-subgraphs-spqr` for SPQR-tree complexity while proposing Tarjan biedged snarl decomposition (`horizon_resolve.rs`).
- **Evidence:**
  - `romain_2026_inversions-pangenome-graphs` (lines 305–377, 479–481; Suppl. Table 3 `suppl` L174–183): `INVPG_annot` v1.0.0 (`SandraLouise/INVPG_annot`) already classifies bubbles >50 bp with `-m 0.5` into path-explicit and alignment-rescued inversions in `2m12s–35m11s` and `174–669 MB` RAM on human chr7+X graphs. `bedtools intersect` in `romain_2026` (lines 599–606) was used only to score against the truth set, not as pipeline glue.
  - `sena_2026_bubble-subgraphs-spqr` (summary, Theorems 1, 2, 5, 11) and `algbio/BubbleFinder` already provide $O(|V|+|E|)$ SPQR-tree snarl/ultrabubble decomposition.
  - `romain_2026_inversions-pangenome-graphs` Table 2 (lines 903–906, 927–931) shows real human inversion recall is **9.8%–53.7%** due to graph builders fragmenting or collapsing inversions, which a different snarl caller cannot fix.
- **Correction:** Narrow C3 to **eliminating the two-step `vg deconstruct` → `INVPG_annot` handoff and external `minimap2` binary dependency** via an integrated GFA-native binary, rather than claiming an unsolved algorithmic or `bedtools`-glue gap.

---

## F4 -- `horizon gfa-prune` (C4) duplicates `vg clip` (which does not require `XG`) and misattributes `kopalli_2025_toward`'s Giraffe failure

- **Claim ID:** `C4`
- **Verdict:** `NARROWED`
- **Concrete flaw:** `vg clip` already performs allele-frequency filtering (`-d`) and deletion-bubble capping (`-D`) without requiring an `XG` index, and `kopalli_2025_toward`'s Giraffe indexing failure on PGGB graphs was a Giraffe/GBWT indexing bottleneck rather than an AF-filtering or `XG` conversion gap.
- **Evidence:**
  - `liao_2023_draft` (Methods lines 3959–3960, 4045, 4357) and `hickey_2023_pangenome` (lines 345–355, 1827–1828) already use `vg clip -d 9 -m 10000` (AF/haplotype-support filtering) and deletion-edge filtering (`-D 10000000`), alongside `vcfbub -l 0 -r 100000` for VCF bubble filtering. `vg clip` operates on `HandleGraph` (`PackedGraph`/`HashGraph`/GFA) without `XG`.
  - `kopalli_2025_toward` (lines 586–588, 772–773) states: *"We were unable to generate mapping indexes for the PGGB pangenome graph due to heavy computational and memory requirements"*; pruning GFA walks does not build a Giraffe index.
- **Correction:** Narrow C4 to a **standalone, streaming two-pass Rust CLI alternative to `vg clip`** that avoids installing `vg` and reduces peak RAM on large GFAs, removing the `XG` and Giraffe-indexing motivation claims.

---

## F5 -- `horizon eval-null` (C5) duplicates `tantan`/`LASTEX`, `Alignathon`/`mafComparator`, and `bedtools`, and relies on uncalibrated CRC

- **Claim ID:** `C5`
- **Verdict:** `REFUTED`
- **Concrete flaw:** Reversed/shuffled-sequence empirical FPR benchmarking and cross-aligner concordance are already established and implemented in existing tools (`LAST`/`LASTEX`, `uShuffle`, `mafComparator`/`Alignathon`, `bedtools`), leaving no general community gap beyond `horizon-aligner`'s own internal CRC calibration.
- **Evidence:**
  - `frith_2011_tantan` (lines 362, 379–386) already defines and runs the reversed-sequence and shuffled-sequence empirical FPR protocol using `LAST` and `LASTEX` (expected random alignments for sequence length and composition), and explicitly notes (summary Limits) that reversed-sequence tests measure false positives only and require masking both sequences.
  - `earl_2014_alignathon` (lines 402–427, 826–873, 953–956) already implements alignment-relation precision/recall/F-score, Jaccard cross-aligner concordance, and PSAR statistical assessment (and `mafTools`/`mafComparator` compares MAF alignments).
  - `quinlan-hall_2010_bedtools` (`bedtools shuffle`, `bedtools intersect`, `bedtools jaccard`) already provides interval permutation and concordance.
  - `reports/rounds/round-8/round.json` (`R7-CRC-BIO`) records that `calibrate_null_crc_from_seqs` itself is not yet calibrated on biological empirical negatives beyond 64-bp synthetic decoys.
- **Correction:** Refute C5 as a community utility; keep reversed-sequence / CRC calibration as an internal `horizon-aligner` evaluation script (`horizon_eval.rs`).

---

## F6 -- `horizon synteny` (C6) misreads `shivakumar_2025_mumemto`'s RAM, ignores held RagTag/Panagram literature, and conflates pairwise chaining with N-way multi-MUMs

- **Claim ID:** `C6`
- **Verdict:** `REFUTED`
- **Concrete flaw:** Three independent prior-art errors invalidate C6: (1) `shivakumar_2025_mumemto`'s `707–800 GB` figure is the *8-thread total memory footprint*, whereas serial Mumemto runs 320 human assemblies in `139 GB` (and chr19 in `16.2 GB`); (2) `horizon synteny`'s pairwise `FlatSeedIndex` + LIS chaining does not compute N-way multi-MUMs/multi-MEMs and cannot catch consensus orientation flips that `shivakumar_2025_mumemto` showed pairwise alignment (`RagTag`) misses; and (3) two held papers in `refs/catalog_index.json`—`alonge_2022_automated` (`RagTag`) and `jenike_2026_panagram` (`Panagram`)—are uncited in `manuscript/references.bib` and `COMMUNITY_UTILITIES_PROPOSAL.md`, yet directly address assembly chimera/orientation correction (`RagTag correct`) and alignment-free k-mer pangenome profiling (`Panagram`).
- **Evidence:**
  - `shivakumar_2025_mumemto` (lines 128–130, 361–362; Table 2; `moesm1` L31–37): `707.5 GB` (89 human haplotypes) and `800 GB` (320 human assemblies, 8 threads) are *"the total memory footprint across all parallel threads"*; **serially**, 320 human assemblies completed *"under a week within 139 GB of memory"*.
  - `alonge_2022_automated` (held in `refs/catalog_index.json`, uncited in `manuscript/references.bib` and proposal; summary citing lines 309–343): `RagTag correct` already uses `minimap2`/`nucmer` PAF/delta alignments to detect and break `--inter`/`--intra` chimeric contigs at merged-alignment boundaries in low RAM.
  - `shivakumar_2025_mumemto` (lines 331–336) specifically showed that `RagTag` (pairwise `minimap2`-based) produced chr8 orientation-flip errors with confidence `1.0` that only *multi-genome* consensus MUMs caught—demonstrating that pairwise reference-vs-assembly chaining (`horizon synteny`) does not solve the orientation-flip problem.
  - `jenike_2026_panagram` (held in `refs/catalog_index.json`, uncited in `manuscript/references.bib` and proposal) already provides k-mer-based pangenome containment/multiplicity across 160–605 genomes without WGA or suffix arrays.
  - Existing pairwise tools (`minimap2` + `SyRI`, `MUMmer4`/`nucmer`, `QUAST`) already perform low-RAM pairwise chimeric-join and orientation-flip scanning.
- **Correction:** Refute C6. Do not build `horizon synteny` as a low-RAM Mumemto replacement; cite `alonge_2022_automated` and `jenike_2026_panagram` if discussing assembly-QC or k-mer pangenome tools.

---

## Claim Verdict Summary

| Claim ID | Verdict | Summary of Prior-Art Ruling |
|---|---|---|
| `C1` | `NARROWED` | Duplicates `halTools` (`halStats`, `hal2maf`, `halLiftover`) + `hal2vg` (available via Bioconda/Cactus containers); narrows to zero-`libhdf5` static-binary portability. |
| `C2` | `REFUTED` | `odgi` (`stats`, `depth`, `paths`), `panacus`, and `gfatools` already perform depth/AF/component/path-bp QC; `horizon gfa-qc` does not replace `kopalli_2025_toward`'s `minimap2`+`Liftoff` gene-loss/alignment checks. |
| `C3` | `NARROWED` | `INVPG_annot` already implements the exact path-explicit + alignment-rescued algorithm (`--min-len 50 --min-cov 0.5`), and `BubbleFinder` already implements SPQR-tree snarls; narrows to single-binary UX without external `vg deconstruct` + `minimap2`. |
| `C4` | `NARROWED` | `vg clip` (`-d`, `-D`, `-m`) already performs AF and deletion-bubble pruning without `XG`, and `kopalli_2025_toward`'s Giraffe failure is a PGGB/Giraffe index bottleneck; narrows to standalone streaming Rust alternative to `vg clip`. |
| `C5` | `REFUTED` | `frith_2011_tantan` (`LAST`/`LASTEX`), `earl_2014_alignathon` (`Alignathon`/`mafComparator`), and `quinlan-hall_2010_bedtools` already cover reversed/shuffled FPR and concordance; only internal CRC calibration remains. |
| `C6` | `REFUTED` | Misreads Mumemto's 8-thread `707–800 GB` total vs `139 GB` serial RAM; pairwise `FlatSeedIndex`+LIS cannot match N-way multi-MUMs and duplicates `minimap2`+`SyRI`/`QUAST`/`RagTag` (`alonge_2022_automated`, held & uncited, alongside `jenike_2026_panagram`). |
