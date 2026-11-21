# Round 8 Attack — `provenance`

- **Target:** `reports/rounds/round-8/ATTACKED-COMMUNITY_UTILITIES_PROPOSAL.md` (`docs/COMMUNITY_UTILITIES_PROPOSAL.md`)
- **Literature provenance:** Checked against held summaries in `/home/rcb112/project/claude/bench/repos/library/papers_md/summaries/` at pinned revision `d4df9122eddd84694b2c3215e5478c46f4314044` (`refs/LIBRARY_VERSION`; `garrison_2024_building.md` read via `git show` at `d4df9122`, remaining 15 summaries verified identical between `d4df9122` and `HEAD`).
- **Method cards (Rule 3):** This project keeps **no method cards** yet (`docs/round_rubric.md`).
- **Prior-round adjudications (Rule 4):** `scripts/round history` returned zero prior mentions of `COMMUNITY_UTILITIES_PROPOSAL`, `horizon hal`, or `gfa-qc`; prior rounds 1, 4, 5, and 6 adjudicated `src/egraph.rs` inversion/duplication resolution, not this proposal.

---

## F1 -- C1: Unverified multi-GiB HAL benchmarks and overstated `to_maf_for_genome` / HDF5 streaming capabilities

- **Claim ID:** C1
- **Verdict:** REFUTED
- **Concrete flaw:**
  1. `ATTACKED-COMMUNITY_UTILITIES_PROPOSAL.md:41, 42, 117` claims `src/hal/hdf5_read.rs` is *"already verified against 5.5 GiB `4-t2t-orangs-mc-2023v2.full.hal` and 56.4 GiB `jarvis_lab_birds.hal`, including stacked multi-layer zlib deflate filters"*. Neither filename appears anywhere in the repository (`src/`, `docs/`, `manuscript/`, `scripts/`, `tests/`, or `git log` prior to commit `0f31c26`). The only `.hal` fixture in the repo is `tests/data/halTest.hal` (`src/hal/tests.rs:242`).
  2. `src/hal/hdf5_read.rs` only provides low-level HDF5 v0 group/dataset readers (`inspect_hal_genomes`, `read_hdf5_v0_dataset_bytes_seek`) and `tests/data/halTest.hal` only tests decoding `SEQNAME_ARRAY` (`src/hal/tests.rs:240-280`). There is **no** HDF5-to-`HalAlignment` loader in `src/hal/`. Moreover, `read_hdf5_v0_dataset_bytes_seek` loads entire datasets into a `Vec<u8>` up to `MAX_CHUNK_DECOMPRESSED_BYTES` (`256 MiB`, `src/hal/hdf5_read.rs:9`), which is not streaming.
  3. `HalAlignment::to_maf_for_genome` (`src/hal/mod.rs:416-470`) only projects a single parent-child level in memory (`root.bottom_segments` -> `root.children`, or `leaf -> parent`), requires `root.dna` in memory, and emits `root.name` with global concatenated coordinates rather than per-contig `HalSequence` names/offsets.
- **Evidence:**
  - `grep -rn "orangs\|jarvis" src/ docs/ manuscript/ scripts/ tests/` returns only `docs/COMMUNITY_UTILITIES_PROPOSAL.md:41, 42, 117`.
  - `src/hal/tests.rs:240-280` (`real_hal_decode_from_env` on `tests/data/halTest.hal`).
  - `src/hal/hdf5_read.rs:9, 683-843`; `src/hal/mod.rs:416-470`.
  - Literature citations (`hickey_2013_hierarchical`, `liao_2023_draft`, `armstrong_2020_progressive`) accurately motivate HAL usage.
- **Correction:** Remove the unverified claims about `4-t2t-orangs-mc-2023v2.full.hal` and `jarvis_lab_birds.hal` (or record reproducible verification logs in `docs/BENCHMARK_RESULTS.md`), and state accurately that `src/hal/hdf5_read.rs` currently provides low-level HDF5 v0 metadata/dataset readers (verified on `tests/data/halTest.hal`), while `HalAlignment::to_maf_for_genome` is an in-memory single-level projector requiring an HDF5-to-`HalAlignment` loader and streaming multi-level tree traversal to match `hal2maf`.

---

## F2 -- C2: Citations and `~97 MB` counter-only RSS verified, with minor provenance scope notes

- **Claim ID:** C2
- **Verdict:** SURVIVED
- **Concrete flaw:** None that invalidates C2; two minor provenance qualifications apply:
  1. The `15–22%` clipping range for *Arabidopsis* / *E. coli* / gibbon (`ATTACKED-COMMUNITY_UTILITIES_PROPOSAL.md:48`) combines text figures from `garrison_2024_building` Supplementary E (`22.1%` *E. coli* at Sup L237-238; `>20%` gibbon `mSymSyn` at Sup L413-414; non-centromeric gene clipping in *A. thaliana* at Sup L362-363) with `15–17%` *A. thaliana* read off the rendered Supplementary Fig. 2 rather than printed text.
  2. The `~97 MB` peak RSS (`ATTACKED-COMMUNITY_UTILITIES_PROPOSAL.md:59, 118`) is verified in `docs/BENCHMARK_RESULTS.md:24, 64` (`96.9 MB` RSS in `196.93s` for `horizon_gfa stats` on `hprc-v1.0-mc-grch38.gfa.gz`), but `horizon_gfa stats` (`src/bin/horizon_gfa.rs:20-56`) is a stateless counter-only pass and `stream_gfa` (`docs/PANGENOME_SUITE.md:3`) skips `C`/`J` records.
- **Evidence:**
  - `liao_2023_draft.md:49-50` (`188.6 Mb` per non-reference assembly).
  - `garrison_2024_building.md` (at `d4df9122`, Sup L237-238, L362-363, L413-414).
  - `kopalli_2025_toward.md:44-48` (`odgi paths`, `minimap2`, `samtools`, `Liftoff`).
  - `andreace_2023_comparing.md:65-68` (giant components on H10).
  - `guarracino_2022_odgi.md:56-58` (`10^5–10^6`-fold node depth).
  - `docs/BENCHMARK_RESULTS.md:24, 64`; `src/bin/horizon_gfa.rs:20-56`.
- **Correction:** None required for C2's provenance; optionally note that `96.9 MB` RSS is for `horizon_gfa stats`'s stateless counter pass, whereas `gfa-qc` component/depth tracking across `81.4M` segments will require additional state.

---

## F3 -- C3: `horizon_resolve.rs` does not perform Tarjan snarl decomposition, and `sena_2026` is mischaracterized

- **Claim ID:** C3
- **Verdict:** REFUTED
- **Concrete flaw:**
  1. **Overstated codebase capability:** `ATTACKED-COMMUNITY_UTILITIES_PROPOSAL.md:25, 73, 119` claims `horizon_resolve.rs` provides `PackedGraph + Tarjan snarls` / `Tarjan biedged snarl decomposition`. In `src/bin/horizon_resolve.rs:210-250, 483-527`, `tarjan_bridges` only counts `num_bridges` and `num_components`, while `decompose_and_resolve` (`lines 483-527`) slices `P`/`W` walks into 1D sliding windows between high-visit `is_anchor` segments (`j - i < max_bubble_steps`) without using Tarjan bridges or 2-edge cuts to extract snarls or their allele subgraphs.
  2. **Overstated `FlatSeedIndex` capability:** `ATTACKED-COMMUNITY_UTILITIES_PROPOSAL.md:74` attributes `FlatSeedIndex strand-symmetric 11-of-15 + banded extension` as a reusable unit. `FlatSeedIndex` (`src/seeding/flat_index.rs:11`) is a generic spaced-seed index without built-in banded extension (which lives in `src/mdl/xdrop.rs` / `src/mdl/dp.rs`).
  3. **Misattributed literature claim (`sena_2026_bubble-subgraphs-spqr`):** `ATTACKED-COMMUNITY_UTILITIES_PROPOSAL.md:67` states that `sena_2026_bubble-subgraphs-spqr` *"avoids Paten (2018)'s quadratic worst case when enumerating all snarls/ultrabubbles"*. In `sena_2026_bubble-subgraphs-spqr.md:12-14, 34-37`, the total number of snarls can still be quadratic in `|V|` (so enumerating all snarls is not linear-time); Sena (2026) gives a linear-size *representation* of all snarls (Theorem 2) and linear-time enumeration of *ultrabubbles* (Theorem 5/11), and is an unreviewed arXiv preprint with no implementation or benchmark (`sena_2026_bubble-subgraphs-spqr.md:7, 51-54`).
- **Evidence:**
  - `src/bin/horizon_resolve.rs:210-250, 483-527`.
  - `src/seeding/flat_index.rs:11-29`.
  - `sena_2026_bubble-subgraphs-spqr.md:7, 12-14, 34-37, 51-54`.
  - `romain_2026_inversions-pangenome-graphs.md:65-67, 86-94` (`path-explicit` and `alignment-rescued` topologies are accurately cited).
- **Correction:** State that `horizon_resolve.rs` currently provides `u32` `PackedGraph` + `tarjan_bridges` bridge counting and 1D anchor-window walk slicing (actual 2-edge-cut / SPQR snarl extraction must be built), and narrow the `sena_2026_bubble-subgraphs-spqr` citation to linear-time *ultrabubble* enumeration and a linear-size *representation* of snarls.

---

## F4 -- C4: `vcfbub` VCF filtering conflated with GFA graph pruning, and `horizon_clean.rs` splits rather than bridges walks

- **Claim ID:** C4
- **Verdict:** NARROWED
- **Concrete flaw:**
  1. **Conflation of VCF post-processing with GFA pruning:** `ATTACKED-COMMUNITY_UTILITIES_PROPOSAL.md:80` cites `liao_2023_draft` as using `vcfbub -l 0 -r 100000` *"to drop >100 kb / >10 Mb spurious deletion edges"*. In `liao_2023_draft.md:70-88`, `vcfbub` (`-l 0 -r 100000` and `-l 0 -r 10000000`) was run on **`vg deconstruct` VCF output** (`lines 4042-4045, 4088-4089, 4506`), whereas graph-level >10 Mb deletion removal used `vg clip -D 10000000` (`lines 3953-3957`). Additionally, `liao_2023_draft.md:124-125` notes a printed discrepancy between Methods (`vg clip -d 9 -m 10000`, line 3959) and the mapping graph name (`clip.d9.m1000.D10M.m1000`, line 4357).
  2. **Overstated `horizon_clean.rs` capability:** `ATTACKED-COMMUNITY_UTILITIES_PROPOSAL.md:87` describes *"splitting or bridging cleanly using the walk-splitting and CIGAR-overlap logic already implemented in `horizon_clean.rs`"*. `src/bin/horizon_clean.rs:410-466` implements walk-splitting (`emit_clean_walk`, `emit_clean_path`) and `parse_cigar_overlap_bp` (`lines 151-168`), but has **no** walk-bridging logic.
- **Evidence:**
  - `liao_2023_draft.md:70-88, 124-125`; `kopalli_2025_toward.md:22, 93-95` (PGGB Giraffe indexing failure verified).
  - `src/bin/horizon_clean.rs:151-168, 410-466`.
- **Correction:** Separate `vg clip -D 10000000` (graph deletion pruning) from `vcfbub` (VCF post-processing) in `liao_2023_draft`, and state that `horizon_clean.rs` provides walk-splitting and CIGAR-overlap tracking while walk-bridging across pruned bubbles must be newly implemented.

---

## F5 -- C5: `quinlan-hall_2010_bedtools` misattributed as WGA ground-truth critique; `calibrate_null_crc_from_seqs` scope overstated

- **Claim ID:** C5
- **Verdict:** NARROWED
- **Concrete flaw:**
  1. **Misattributed citation:** `ATTACKED-COMMUNITY_UTILITIES_PROPOSAL.md:27, 92` includes `quinlan-hall_2010_bedtools` in the list of papers emphasizing that *"real-genome WGA lacks ground truth"*. `quinlan-hall_2010_bedtools` is a two-page BEDTools applications note on genomic interval overlaps; it never discusses WGA ground truth (`quinlan-hall_2010_bedtools.md:49-80` only mentions `shuffleBed` for BED feature permutation enrichment).
  2. **Overstated `calibrate_null_crc_from_seqs` capability:** `ATTACKED-COMMUNITY_UTILITIES_PROPOSAL.md:98-99` cites `calibrate_null_crc_from_seqs` (`src/mdl/dp.rs`) for *"GC-stratified bit-score distributions"* and dinucleotide-preserving shuffles. In `src/mdl/dp.rs:202-260`, `calibrate_null_crc_from_seqs` samples `64`-bp (`WIN = 64`) reversed `y_rev` windows at a single estimated GC fraction and returns one scalar `f64` threshold, without dinucleotide shuffling or multi-bin GC distribution reporting.
- **Evidence:**
  - `quinlan-hall_2010_bedtools.md:49-80`; `earl_2014_alignathon.md:25-30`; `frith_2011_tantan.md:47-55`; `garrison_2024_building.md` (at `d4df9122`, lines 1027-1028); `kopalli_2025_toward.md:113-115`.
  - `src/mdl/dp.rs:202-260`; `src/bin/horizon_eval.rs:1-120`.
- **Correction:** Remove `quinlan-hall_2010_bedtools` from the WGA ground-truth sentence (or cite `shuffleBed` specifically for interval-permutation nulls), and clarify that `calibrate_null_crc_from_seqs` currently returns a scalar `64`-bp reversed-window CRC floor rather than GC-stratified distributions or dinucleotide shuffles.

---

## F6 -- C6: `FlatSeedIndex` has a hard 4 Gbp `u32` single-sequence limit

- **Claim ID:** C6
- **Verdict:** NARROWED
- **Concrete flaw:** `shivakumar_2025_mumemto` (`707.5–800 GB` RAM, `HG02080` `chr17–chr19` chimeric join, RagTag `chr8` orientation flips, pericentromeric LTR spikes), `liu_2018_inferring` (up to `40%` underestimation on fragmented assemblies), `FlatSeedIndex` (`12 bytes/entry`), `lis_chain` (`src/seeding/chain.rs`), and `sdust_mask` (`src/seeding/mod.rs`) are all verified. However, `FlatSeedIndex::build_with_sdust` (`src/seeding/flat_index.rs:41-44`) indexes a single sequence (`seq: &[u8]`) and explicitly rejects sequences exceeding `u32::MAX` (`4 Gbp`) with `"sequence length exceeds 4 Gbp u32 index limit"`, so multi-assembly pangenome scanning requires per-chromosome/per-contig sharding or a multi-sequence index wrapper.
- **Evidence:**
  - `shivakumar_2025_mumemto.md:55-57, 67-69, 79-88`; `liu_2018_inferring.md:69-74`.
  - `src/seeding/flat_index.rs:6, 41-44`; `src/seeding/chain.rs:9-42`; `src/seeding/mod.rs:54`.
- **Correction:** Note that `FlatSeedIndex` is a single-sequence index capped at `4 Gbp` (`u32` positions), so `horizon synteny` must index per-reference/per-contig rather than a concatenated multi-genome pangenome.

---

## Claim Verdict Summary

| Claim | Verdict | Summary of Provenance Finding |
|---|---|---|
| C1 | REFUTED | `4-t2t-orangs-mc-2023v2.full.hal` & `jarvis_lab_birds.hal` verification absent from repo; no HDF5->`HalAlignment` loader; `to_maf_for_genome` is in-memory single-level. |
| C2 | SURVIVED | All 5 citations & `96.9 MB` `horizon_gfa stats` RSS verified (`docs/BENCHMARK_RESULTS.md:24, 64`); minor caveat that `stats` is stateless counter-only. |
| C3 | REFUTED | `horizon_resolve.rs` only counts Tarjan bridges and slices 1D walk windows (no snarl extraction); `sena_2026` does not enumerate all snarls in linear time. |
| C4 | NARROWED | `liao_2023_draft` used `vcfbub` on VCFs (graph deletion pruning was `vg clip -D 10000000`); `horizon_clean.rs` splits walks but does not bridge them. |
| C5 | NARROWED | `quinlan-hall_2010_bedtools` does not discuss WGA ground truth; `calibrate_null_crc_from_seqs` returns a scalar `64`-bp reversed-window threshold. |
| C6 | NARROWED | Literature and `12 bytes/entry` / `lis_chain` / `sdust_mask` verified, but `FlatSeedIndex` is single-sequence with a hard `4 Gbp` `u32` limit. |
