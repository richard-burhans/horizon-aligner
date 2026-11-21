# Round 8 Attack — Lens: `correctness`

*Note on rules 2 & 3: This project keeps no method cards, and no `../library` checkout was consulted; all findings below test C1–C6 directly against the `horizon-aligner` source tree and algorithmic/memory invariants.*

---

## F1 -- C1 (`horizon hal`): NARROWED — no HDF5-to-`HalAlignment` loader exists, `256 MiB` dataset cap breaks on real multi-GiB HALs, and `to_maf_for_genome` is single-level only

- **Claim:** C1 (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:23,34-42`) proposes `horizon hal` (`stats`, `to-maf`, `liftover`, `to-gfa`) reusing `src/hal/hdf5_read.rs` (claimed verified on 5.5 GiB `4-t2t-orangs-mc-2023v2.full.hal` and 56.4 GiB `jarvis_lab_birds.hal`), `src/hal/mod.rs` (`HalAlignment::to_maf_for_genome`), and `src/coord.rs`.
- **Verdict:** **NARROWED**
- **Evidence:**
  1. **No HDF5-to-`HalAlignment` loader exists:** `HalAlignment::to_maf_for_genome` (`src/hal/mod.rs:416-495`) operates solely on an in-memory `HalAlignment` struct constructed by `horizon_guide` (`src/bin/horizon_guide/main.rs:315-322`), never on an HDF5 `.hal` file on disk. `inspect_hal_genomes` (`src/hal/hdf5_read.rs:832-843`) only reads the HDF5 superblock and v0 root-group member names (`read_hdf5_v0_root_names_seek`).
  2. **256 MiB dataset cap blocks real multi-GiB HAL datasets:** `read_hdf5_v0_dataset_bytes_seek` (`src/hal/hdf5_read.rs:9,758,790-804`) enforces `MAX_CHUNK_DECOMPRESSED_BYTES = 256 * 1024 * 1024` and materializes the entire dataset into a single `Vec<u8>`, returning `Ok(None)` if decompressed output exceeds 256 MiB. In a 56.4 GiB avian HAL (`jarvis_lab_birds.hal`) or 5.5 GiB primate HAL (`4-t2t-orangs-mc-2023v2.full.hal`), a single genome's `DNA_ARRAY` (~1.5 GB for 3 Gbp) or `TOP_ARRAY`/`BOTTOM_ARRAY` exceeds 256 MiB, so `read_hdf5_v0_dataset_bytes_seek` fails silently (`Ok(None)`) without streaming chunk iteration.
  3. **Single-level MAF projection & missing subcommands:** `HalAlignment::to_maf_for_genome` (`src/hal/mod.rs:432-493`) only projects one tree level (`root.bottom_segments` onto `root.children`, or `write_leaf_maf` to immediate `parent`), lacking multi-level clade traversal (`--target-genomes`) and `--region` slicing. Neither `to-gfa` nor `liftover` exists (`src/coord.rs:5-44` is only a 60-line `Interval` struct).
- **Correction:** `horizon hal` is worth building as P1, but must be scoped as requiring: (a) a streaming chunked HDF5 dataset reader without the 256 MiB `Vec<u8>` cap, (b) an HDF5-to-segment/tree loader, and (c) multi-level phylogenetic MAF/liftover traversal.

---

## F2 -- C2 (`horizon gfa-qc`): NARROWED — `~97 MB` counter architecture cannot compute connected components or per-node depth/AF in a single pass

- **Claim:** C2 (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:24,46-59`) proposes single-pass `horizon gfa-qc` (`<out>.paths.tsv`, `<out>.components.tsv`, `<out>.depth_hist.tsv`) reusing `src/gfa/mod.rs`, `src/io_util.rs`, and "the low-memory counter architecture of `src/bin/horizon_gfa.rs` (`~97 MB` peak RSS on the `12 GB` HPRC graph)".
- **Verdict:** **NARROWED**
- **Evidence:**
  1. **Stateless counters vs. stateful graph topology/depth:** `horizon_gfa stats` (`src/bin/horizon_gfa.rs:20-56`) achieves `96.9 MB` peak RSS only because it maintains five scalar `u64` counters (`s_count`, `s_bases`, `l_count`, `p_count`, `w_count`) and zero per-segment or per-link state.
  2. **Single-pass component & depth/AF tracking requires multi-GB RAM:** In GFA, `S`, `L`, `P`, and `W` lines reference segments by arbitrary string names. Computing connected components across `81.4M` segments and `113M` links plus per-node depth and haplotype-support (`AF`) histograms across `24,456` walks in a single pass requires a segment-name-to-`u32` interner (`HashMap<String, u32>` for 81.4M strings: several GB), Union-Find parent/rank arrays (`~325 MB`), and per-node depth/AF counters (`~325–650 MB`).
  3. **Per-record heap allocation:** `stream_gfa` (`src/gfa/mod.rs:195-211`) allocates owned `String`/`Vec<u8>` fields per record, and `parse_walk_steps` (`src/gfa/mod.rs:167-189`) allocates a `Vec` per walk.
- **Correction:** `horizon gfa-qc` is worth building as P1 (multi-GB RAM is still far below `vg`/`odgi`'s 70–154 GB), but the proposal must drop the `~97 MB` single-pass RSS claim and budget for segment-ID interning + Union-Find + depth arrays.

---

## F3 -- C3 (`horizon inversions`): REFUTED — `PackedGraph` stores no snarl tree/sequences and `FlatSeedIndex::intersect` conflates forward and RC hits

- **Claim:** C3 (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:25,63-75`) proposes `horizon inversions` reusing `PackedGraph` + Tarjan snarls (`src/bin/horizon_resolve.rs`) and `FlatSeedIndex` (`src/seeding/flat_index.rs`) to enumerate snarls/bubbles and classify `path_explicit` vs `alignment_rescued` inversions without `minimap2`.
- **Verdict:** **REFUTED**
- **Evidence:**
  1. **`PackedGraph` has no snarl enumeration, no adjacency list, and no sequences:** In `src/bin/horizon_resolve.rs:170-177`, `PackedGraph` stores only `seg_to_id`, `id_to_seg`, `seg_lens`, `is_anchor: Vec<bool>`, `num_bridges`, and `num_components`. `compute_biedged_snarls` (`src/bin/horizon_resolve.rs:272-374`) runs `tarjan_bridges` + 2-ECC BFS solely to set `is_anchor: Vec<bool>`, not to build a snarl tree, enumerate `(start_side, end_side)` ultrabubbles, or implement SPQR decomposition (`sena_2026_bubble-subgraphs-spqr`). `horizon_resolve.rs` only inspects linear `P`/`W` windows capped at `max_bubble_steps = 16` and discards `S` sequences.
  2. **`FlatSeedIndex::intersect` cannot distinguish forward vs. reverse-complement matches:** `FlatSeedIndex::build_with_sdust` (`src/seeding/flat_index.rs:73`) stores canonical keys `fwd.min(rc)` without recording a strand bit. Consequently, `FlatSeedIndex::intersect` (`src/seeding/flat_index.rs:246-247`) pushes **every** shared seed to **both** `hits_fwd` and `hits_rc`. Using `FlatSeedIndex` alone to rescue inversions on allele sequences produces spurious RC hits for every forward-strand match unless base-level orientation verification is added. Moreover, `FlatSeedIndex` is `pub(crate)` (`src/seeding/flat_index.rs:11`), not accessible to binaries.
- **Correction:** Do not build `horizon inversions` on these primitives as proposed; snarl/ultrabubble enumeration and strand-aware allele alignment would have to be built from scratch.

---

## F4 -- C4 (`horizon gfa-prune`): REFUTED — oversized deletion shortcuts cannot be detected in two streaming passes, and `horizon_clean.rs` lacks AF counting and walk bridging

- **Claim:** C4 (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:26,79-87`) proposes a two-pass streaming `horizon gfa-prune` (`--min-af 0.10`, `--max-del-bp 100000`) reusing two-pass `W`/`P` walk-splitting and CIGAR-overlap logic in `src/bin/horizon_clean.rs` to prune unsupported segments/links and oversized deletion shortcuts without `vg`/`XG`.
- **Verdict:** **REFUTED**
- **Evidence:**
  1. **Deletion-bubble size (`--max-del-bp 100000`) is undefined without a graph/reference topology index:** A GFA link `L u + v + 0M` carries no genomic span. Determining that `u -> v` is a >100 kb deletion shortcut (like `vcfbub -l 0 -r 100000`) requires knowing the length of the alternative path(s) between `u` and `v` in a snarl/reference coordinate index — which two streaming passes over `L` and `W` lines without a graph index cannot compute.
  2. **`horizon_clean.rs` counts raw edge traversals, not haplotype AF:** `record_walk_support` (`src/bin/horizon_clean.rs:178-194`) increments `walk_support[e]` for every adjacent step pair in `P`/`W` lines without deduplicating by `(sample, haplotype)`. A repetitive segment traversed multiple times by a single haplotype inflates support.
  3. **`horizon_clean.rs` only splits walks, never bridges pruned segments:** `emit_clean_path` and `emit_clean_walk` (`src/bin/horizon_clean.rs:326-466`) split walks into `_part1`, `_part2` at severed `L` links; they contain no logic to bridge/splice across pruned segments or synthesize valid replacement `L` edges without fragmenting haplotypes.
- **Correction:** `horizon gfa-prune` cannot work as a two-pass streaming extension of `horizon_clean.rs`; deletion-bubble capping requires at least a reference-coordinate or snarl-boundary index.

---

## F5 -- C5 (`horizon eval-null`): REFUTED — `calibrate_null_crc_from_seqs` only calibrates a single 64-bp scalar threshold and `horizon_eval.rs` cannot score external aligner MAFs

- **Claim:** C5 (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:27,91-99`) proposes `horizon eval-null` (`--aln out.maf`, `--decoy-aln reversed.maf`, `--n-shuffles 20`) to compute empirical FPR on reversed/dinucleotide-shuffled decoys, GC-stratified bit-score distributions, and cross-aligner concordance by reusing `calibrate_null_crc_from_seqs` (`src/mdl/dp.rs`) and `src/bin/horizon_eval.rs`.
- **Verdict:** **REFUTED**
- **Evidence:**
  1. **`calibrate_null_crc_from_seqs` does not do dinucleotide shuffles or GC-stratified distributions:** In `src/mdl/dp.rs:202-287`, `calibrate_null_crc_from_seqs` estimates a single global `obs_gc` from the first 4,096 bp (`estimate_sampled_gc`, lines 180-199), draws i.i.d. mononucleotide 64-bp windows (`calibrate_null_crc_for_gc`, lines 120-174) plus `n_cal` reversed 64-bp genomic windows (`const WIN: usize = 64`, lines 233-286), runs Horizon's own `banded_mdl_gap_with_scratch`, and returns a single `f64` scalar threshold `λ̂`. It implements neither dinucleotide-preserving shuffling nor GC-stratified bit-score distributions nor empirical FPR (aligned bp/kb) reporting.
  2. **`horizon_eval.rs` hardcodes running `horizon-aligner` itself:** `src/bin/horizon_eval.rs:186-208` calls `align_contig_pair_with_config` directly on `(ref_fa, tgt_fa)` and compares against `truth_maf`. It has no parser or scorer for an external aligner's `--aln out.maf` or `--decoy-aln reversed.maf`.
- **Correction:** Neither `calibrate_null_crc_from_seqs` nor `horizon_eval.rs` provides the claimed decoy-generation, dinucleotide-shuffle, or cross-aligner MAF/PAF scoring primitives.

---

## F6 -- C6 (`horizon synteny`): REFUTED — dense 1-bp `FlatSeedIndex` takes ~37–50 GB RAM for one human genome, `lis_chain` is single-contig monotonic, and SDUST misses interspersed repeats

- **Claim:** C6 (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:28,103-111`) proposes `horizon synteny` using `FlatSeedIndex` (`12 bytes/entry`) + SDUST (`src/seeding/mod.rs`) + LIS chaining (`src/seeding/chain.rs`) to detect inter-chromosomal chimeric joins (`chr17–chr19`), orientation flips, and pericentromeric LTR repeat spikes in "a few GB of RAM rather than hundreds of GB".
- **Verdict:** **REFUTED**
- **Evidence:**
  1. **Memory is ~37–50 GB for a single human genome, not "a few GB":** `FlatSeedIndex::build_with_sdust` (`src/seeding/flat_index.rs:57-74`) indexes every 1-bp offset `0..=(seq.len() - span)` with no minimizer/syncmer thinning. At `12 bytes/entry` (`keys: Vec<u64>` + `pos: Vec<u32>`) plus `16 bytes/entry` (`pairs: Vec<(u64, u32)>`) during sorting, a single 3.1 Gbp human reference requires **37.2 GB** steady-state and **~50 GB** peak RAM, and `seq.len() > u32::MAX` (`src/seeding/flat_index.rs:41-44`) rejects >4 Gbp concatenated inputs. Detecting inter-chromosomal chimeric joins (`chr17–chr19`) requires all reference chromosomes indexed simultaneously.
  2. **Algorithmic mismatch (`lis_chain` and `FlatSeedIndex::intersect`):** `lis_chain` (`src/seeding/chain.rs:9-42`) computes a single monotonic increasing `(x, y)` chain per contig pair, which cannot represent multi-chromosome chimeric splits or inverted segments without multi-contig clustering (`cluster_strand_lcbs`), and `FlatSeedIndex::intersect` (`src/seeding/flat_index.rs:246-247`) emits every hit on both `hits_fwd` and `hits_rc`.
  3. **SDUST does not detect LTR repeat spikes:** `sdust_mask` (`src/seeding/mod.rs:54-120`) masks only simple low-complexity triplet repeats (`W in {16, 24, 32}`), not interspersed LTR/Alu/LINE repeats.
- **Correction:** Building `horizon synteny` in "a few GB of RAM" requires minimizer/syncmer thinning (not dense 1-bp `FlatSeedIndex`), multi-contig/strand-aware chaining, and copy-number/multiplicity tracking rather than SDUST.

---

## Summary Table

| Claim ID | Proposed Utility | Verdict |
|---|---|---|
| **C1** | `horizon hal` | **NARROWED** |
| **C2** | `horizon gfa-qc` | **NARROWED** |
| **C3** | `horizon inversions` | **REFUTED** |
| **C4** | `horizon gfa-prune` | **REFUTED** |
| **C5** | `horizon eval-null` | **REFUTED** |
| **C6** | `horizon synteny` | **REFUTED** |
