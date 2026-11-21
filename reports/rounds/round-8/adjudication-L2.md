# Round 8 Adjudication — `prior` (`prior-art` Lens)

**Target:** [attack-prior-art.md](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-8/attack-prior-art.md) (verified sha256 prefix `01886b808825`, matching the sealed digest in [ADJ-BRIEF-prior.md](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-8/ADJ-BRIEF-prior.md)).  
**Prior rounds:** `scripts/round history "COMMUNITY_UTILITIES_PROPOSAL" --kind adjudication` confirms no prior mentions in closed rounds 1–7. No method cards apply to this prior-art proposal triage.

---

## C1 (`F1` — `horizon hal`) — **UPHELD**

What survives: `horizon hal` (`stats`, `to-maf`, `liftover`, `to-gfa`) survives only as a **zero-`libhdf5` static-binary / embedded-Rust packaging convenience**, not as an unmet algorithmic or container-deployment gap.

Verified against `hickey_2013_hierarchical` (Table 1, lines 126–138, 188–204) and `romain_2026_inversions-pangenome-graphs` (lines 465–468): `halStats`, `hal2maf`, `halLiftover`, and `hal2vg --noAncestors` + `vg convert -f -W` already implement all four subcommands, and Bioconda (`gruning_2018_bioconda` in `refs/catalog_index.json`) and official Cactus containers already cover standard container/HPC deployment. (One minor caveat in the attack: `pangenome/impg` operates on PAF/1aln/TPA rather than `.hal`, so it is not `.hal` prior art, but `halLiftover` already covers `.hal` liftover.) The narrowing fix to [COMMUNITY_UTILITIES_PROPOSAL.md](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/COMMUNITY_UTILITIES_PROPOSAL.md) is free.

---

## C2 (`F2` — `horizon gfa-qc`) — **UPHELD**

What survives: C2 is refuted as an unmet analytical or `minimap2`+`Liftoff` glue gap; at most, `horizon gfa-qc` could be reframed as a low-RSS (`~97 MB`) single-pass GFA streaming wrapper avoiding `odgi build`'s `10 GB` `.og` conversion peak on chr6.

Verified against `guarracino_2022_odgi` (lines 229–230, 254–303, 451, 660) and `kopalli_2025_toward` (lines 227, 298–301): `odgi` (`stats`, `depth`, `paths`), `panacus`, and `gfatools` already compute node depth, AF/coverage histograms, connected components, and path lengths. Furthermore, `kopalli_2025_toward` used `minimap2`+`samtools`+`Liftoff` for base-level assembly-to-path re-alignment and gene-loss annotation, neither of which `horizon gfa-qc`'s `--fasta-manifest` length comparison replaces. Dropping or reframing C2 in [COMMUNITY_UTILITIES_PROPOSAL.md](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/COMMUNITY_UTILITIES_PROPOSAL.md) is free.

---

## C3 (`F3` — `horizon inversions`) — **UPHELD**

What survives: `horizon inversions` survives only as a **single-binary UX convenience** eliminating the two-step `vg deconstruct` → `INVPG_annot` handoff and external `minimap2` binary dependency.

Verified against `romain_2026_inversions-pangenome-graphs` (lines 305–377, 479–481, 599–606, 903–906, 927–931; `suppl` L174–183) and `sena_2026_bubble-subgraphs-spqr` (Theorems 1, 2, 5, 11): `INVPG_annot` v1.0.0 already implements the exact `path_explicit` + `alignment_rescued` classification with the same defaults (`--min-len 50`, `--min-cov 0.5`), `bedtools intersect` was used only for truth-set scoring rather than pipeline glue, and real human recall (`9.8%–53.7%`) is limited by graph builders fragmenting or collapsing inversions rather than by snarl enumeration. Moreover, the proposal cites `sena_2026_bubble-subgraphs-spqr` for SPQR-tree complexity while proposing Tarjan biedged snarls (`horizon_resolve.rs`). (Note: `sena_2026_bubble-subgraphs-spqr` is theory-only in the held `arXiv:2604.08071v1` version, with `BubbleFinder` living in the earlier `arXiv:2511.21919` preprint/repo, which does not affect the contradiction in the proposal.) The narrowing fix is free.

---

## C4 (`F4` — `horizon gfa-prune`) — **UPHELD**

What survives: `horizon gfa-prune` survives only as a **standalone, streaming two-pass Rust CLI alternative to `vg clip`** that avoids installing `vg` and lowers peak RAM, without the `XG` or Giraffe-indexing motivation claims.

Verified against `liao_2023_draft` (lines 3953–3960, 4045, 4357), `hickey_2023_pangenome` (lines 345–355, 1827–1828), and `kopalli_2025_toward` (lines 586–588, 772–773): `vg clip` (`-d 9 -m 10000`, `-D 10000000`) already performs AF and deletion-bubble pruning on `HandleGraph`/GFA without requiring an `XG` index, and `kopalli_2025_toward`'s Giraffe indexing failure on PGGB was a Giraffe/GBWT indexing memory/compute bottleneck rather than an AF-filtering or `XG` gap. The narrowing fix is free.

---

## C5 (`F5` — `horizon eval-null`) — **UPHELD**

What survives: C5 is refuted as a community utility; reversed-sequence / CRC calibration remains strictly an internal `horizon-aligner` evaluation script (`horizon_eval.rs`).

Verified against `frith_2011_tantan` (lines 362, 379–386, 971–972), `earl_2014_alignathon` (lines 402–427, 826–873, 953–956), `quinlan-hall_2010_bedtools` (Table 1: `shuffleBed`, `intersectBed`), and `reports/rounds/round-8/round.json` (`R7-CRC-BIO`): reversed/shuffled empirical FPR (`LAST`/`LASTEX`) and cross-aligner concordance (`Alignathon`/`mafComparator`, `bedtools`) are already established, and `calibrate_null_crc_from_seqs` is not yet calibrated on biological empirical negatives beyond 64-bp synthetic decoys. (Minor note: `bedtools jaccard` was added after the 2010 paper, though `shuffleBed` and `intersectBed` are in Table 1.) Dropping C5 from the community-utility roadmap is free.

---

## C6 (`F6` — `horizon synteny`) — **UPHELD**

What survives: C6 is refuted; `horizon synteny` should not be built as a low-RAM Mumemto replacement, and any assembly-QC or k-mer pangenome discussion should cite `alonge_2022_automated` and `jenike_2026_panagram`.

Verified against `shivakumar_2025_mumemto` (lines 128–130, 331–336, 361–362; `moesm1` L31–37), `alonge_2022_automated` (lines 309–343), `jenike_2026_panagram`, `refs/catalog_index.json`, and `manuscript/references.bib`: (1) Mumemto's `707–800 GB` is the 8-thread total footprint while serial Mumemto runs 320 human assemblies in `139 GB` (and chr19 in `16.2 GB` using PFP, not full SA); (2) pairwise `FlatSeedIndex` + LIS chaining cannot compute N-way multi-MUMs/multi-MEMs or catch consensus orientation flips that `shivakumar_2025_mumemto` showed pairwise alignment (`RagTag`) misses; and (3) `alonge_2022_automated` (`RagTag correct`) and `jenike_2026_panagram` (`Panagram`) are held in `refs/catalog_index.json` yet uncited in `manuscript/references.bib` and `COMMUNITY_UTILITIES_PROPOSAL.md`. Dropping C6 and adding the missing citations if needed is free.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| C1 | UPHELD | `halTools` (`halStats`, `hal2maf`, `halLiftover`) + `hal2vg` and Bioconda/Cactus containers already cover functionality and deployment; only zero-`libhdf5` static-binary packaging survives. |
| C2 | UPHELD | `odgi`, `panacus`, and `gfatools` already compute depth/AF/components/path-bp, and `--fasta-manifest` length checks do not replace `kopalli_2025_toward`'s `minimap2`+`Liftoff` gene/alignment checks. |
| C3 | UPHELD | `INVPG_annot` already implements the exact two-class algorithm (`--min-len 50 --min-cov 0.5`), `bedtools` was truth-set scoring only, and the proposal cites `sena_2026` SPQR while using Tarjan snarls; only single-binary UX survives. |
| C4 | UPHELD | `vg clip` (`-d`, `-D`, `-m`) already runs without `XG`, and `kopalli_2025_toward`'s Giraffe failure was PGGB indexing cost; only a standalone streaming Rust alternative to `vg clip` survives. |
| C5 | UPHELD | `frith_2011_tantan` (`LAST`/`LASTEX`), `earl_2014_alignathon`, and `quinlan-hall_2010_bedtools` already cover reversed/shuffled FPR and concordance, and `R7-CRC-BIO` shows CRC is uncalibrated on biological negatives. |
| C6 | UPHELD | `707–800 GB` is 8-thread Mumemto vs `139 GB` serial, pairwise chaining cannot catch N-way consensus orientation flips, and held `alonge_2022_automated` (`RagTag`) and `jenike_2026_panagram` (`Panagram`) are uncited. |
