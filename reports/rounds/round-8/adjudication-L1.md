# Round 8 ADJUDICATION — Lens: `correctness`

- **Attack reviewed:** `reports/rounds/round-8/attack-correctness.md` (sealed digest `d0178606a809`, verified via `sha256sum`).
- **Prior history:** `scripts/round history "" --kind adjudication` shows only Round 1 (`STRATEGIC_PLAN`, `adjudication-L1.md`); no prior closed rounds adjudicate `docs/COMMUNITY_UTILITIES_PROPOSAL.md`.
- **Method cards:** None exist in this repository; all rulings below are grounded directly in `src/` and `docs/COMMUNITY_UTILITIES_PROPOSAL.md`.

---

## Finding-by-Finding Rulings

### F1 -- C1 (`horizon hal`): **UPHELD**
- **Ruling:** **UPHELD** (attacker verdict `NARROWED` stands).
- **What survives:** `horizon hal` remains worth building as P1, but the proposal must drop the implication that `src/hal/` can already stream or project multi-GiB HAL files and instead scope in a chunk-streaming HDF5 reader, an HDF5-to-`HalAlignment` loader, and multi-level tree traversal.
- **Verification:**
  1. `inspect_hal_genomes` (`src/hal/hdf5_read.rs:832-843`) only reads the superblock and v0 root-group names; no HDF5-to-`HalAlignment` loader exists.
  2. `read_hdf5_v0_dataset_bytes_seek` (`src/hal/hdf5_read.rs:9,758,790-804`) enforces `MAX_CHUNK_DECOMPRESSED_BYTES = 256 * 1024 * 1024`, materializes datasets into a single `Vec<u8>`, and returns `Ok(None)` above 256 MiB, which breaks on multi-GiB HAL `DNA_ARRAY`/`TOP_ARRAY`/`BOTTOM_ARRAY` datasets.
  3. `HalAlignment::to_maf_for_genome` (`src/hal/mod.rs:416-495`) projects only one tree level (immediate children or immediate parent), and `src/coord.rs:5-44` is only a 60-line `Interval` struct with no `liftover` or `to-gfa` implementation.

### F2 -- C2 (`horizon gfa-qc`): **UPHELD**
- **Ruling:** **UPHELD** (attacker verdict `NARROWED` stands).
- **What survives:** `horizon gfa-qc` remains worth building as P1 (multi-GB RAM is still much smaller than `vg`/`odgi` at 70–154 GB), but the proposal cannot claim the `~97 MB` peak RSS of `horizon_gfa stats`.
- **Verification:**
  1. `cmd_stats` (`src/bin/horizon_gfa.rs:20-56`) achieves `~97 MB` RSS only because it keeps five `u64` scalar counters and zero per-segment or per-link state.
  2. Computing connected components and per-node depth/AF histograms across `81.4M` string-named segments in a single pass requires segment-name interning, Union-Find arrays, and per-node counters totaling multiple GB of RAM.
  3. `stream_gfa` (`src/gfa/mod.rs:195-211`) and `parse_walk_steps` (`src/gfa/mod.rs:167-189`) allocate owned per-record strings and per-walk `Vec`s.

### F3 -- C3 (`horizon inversions`): **UPHELD**
- **Ruling:** **UPHELD** (attacker verdict `REFUTED` stands).
- **What survives:** `horizon inversions` cannot be built by reusing `PackedGraph` + `compute_biedged_snarls` (`src/bin/horizon_resolve.rs`) and `FlatSeedIndex` (`src/seeding/flat_index.rs`) as proposed; snarl/ultrabubble enumeration and sequence storage must be built from scratch.
- **Verification:**
  1. `PackedGraph` (`src/bin/horizon_resolve.rs:170-177`) stores no adjacency list, no snarl tree, and no sequences; `compute_biedged_snarls` (`src/bin/horizon_resolve.rs:272-374`) computes `tarjan_bridges` + 2-ECC BFS solely to populate `is_anchor: Vec<bool>`.
  2. `FlatSeedIndex` (`src/seeding/flat_index.rs:11,73,246-247`) is `pub(crate)`, stores canonical `fwd.min(rc)` keys without a strand bit, and pushes every hit to both `hits_fwd` and `hits_rc`. (Note: C3's design does pair `FlatSeedIndex` with banded extension, which would filter spurious RC hits, but the lack of snarl enumeration, sequences, and public visibility still refutes C3's primitive-reuse claim.)

### F4 -- C4 (`horizon gfa-prune`): **UPHELD**
- **Ruling:** **UPHELD** (attacker verdict `REFUTED` stands).
- **What survives:** `horizon gfa-prune` cannot work as a two-pass streaming extension of `src/bin/horizon_clean.rs`; `--max-del-bp` requires a snarl or reference-coordinate topology index, and AF filtering requires `(sample, haplotype)` deduplication plus walk bridging.
- **Verification:**
  1. GFA `L` links carry no genomic span of the bypassed allele(s), so `--max-del-bp 100000` is undefined in two streaming passes without a graph/reference index.
  2. `record_walk_support` (`src/bin/horizon_clean.rs:178-194`) increments `walk_support[e]` per step pair without deduplicating by `(sample, haplotype)`.
  3. `emit_clean_path` and `emit_clean_walk` (`src/bin/horizon_clean.rs:326-466`) only split walks into `_part1`, `_part2` at severed links and have no bridging/splicing logic.

### F5 -- C5 (`horizon eval-null`): **UPHELD**
- **Ruling:** **UPHELD** (attacker verdict `REFUTED` stands).
- **What survives:** Neither `calibrate_null_crc_from_seqs` (`src/mdl/dp.rs`) nor `src/bin/horizon_eval.rs` provides dinucleotide-shuffled decoys, GC-stratified bit-score distributions, or external MAF/PAF scoring.
- **Verification:**
  1. `calibrate_null_crc_from_seqs` (`src/mdl/dp.rs:202-287`) estimates a single global `obs_gc` from the first 4,096 bp (`estimate_sampled_gc`, lines 180-199), evaluates 64-bp i.i.d. mononucleotide and reversed windows with `banded_mdl_gap_with_scratch`, and returns one scalar `f64` threshold.
  2. `align_and_score_pairs` (`src/bin/horizon_eval.rs:186-208`) hardcodes calling `align_contig_pair_with_config` directly and has no parser or scorer for external `--aln` or `--decoy-aln` MAFs.

### F6 -- C6 (`horizon synteny`): **UPHELD**
- **Ruling:** **UPHELD** (attacker verdict `REFUTED` stands).
- **What survives:** `horizon synteny` cannot run in "a few GB of RAM" using dense 1-bp `FlatSeedIndex`, single-contig `lis_chain`, and `sdust_mask`; it requires minimizer/syncmer thinning, multi-contig/strand-aware clustering, and repeat-multiplicity tracking.
- **Verification:**
  1. `FlatSeedIndex::build_with_sdust` (`src/seeding/flat_index.rs:41-44,55-80`) indexes every 1-bp offset (`16 bytes/entry` `pairs` + `12 bytes/entry` `keys`/`pos` simultaneously during build), requiring tens of GB for a single 3.1 Gbp human genome and rejecting `seq.len() > u32::MAX` (>4 Gbp).
  2. `lis_chain` (`src/seeding/chain.rs:9-42`) produces a single monotonic `(x, y)` chain per contig pair rather than multi-chromosome chimeric splits or inversions, and `FlatSeedIndex::intersect` (`src/seeding/flat_index.rs:246-247`) duplicates every hit into `hits_fwd` and `hits_rc`.
  3. `sdust_mask` (`src/seeding/mod.rs:54-120`) masks only low-complexity triplet repeats (`W in {16, 24, 32}`), not interspersed LTR/Alu/LINE repeats.

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `hdf5_read.rs` has a 256 MiB `Vec<u8>` cap and no HDF5-to-`HalAlignment` loader, `to_maf_for_genome` is single-level, and `coord.rs` is only `Interval`; narrowing C1 holds |
| F2 | UPHELD | `horizon_gfa stats` uses `~97 MB` only because it has 5 scalar counters; single-pass components and depth/AF require multi-GB interning and arrays; narrowing C2 holds |
| F3 | UPHELD | `PackedGraph` has no adjacency list, snarl tree, or sequences, and `FlatSeedIndex` is `pub(crate)` with canonical keys lacking strand bits; refuting C3 holds |
| F4 | UPHELD | Two-pass streaming cannot compute `--max-del-bp` without a topology index, and `horizon_clean.rs` lacks haplotype-deduplicated AF and walk bridging; refuting C4 holds |
| F5 | UPHELD | `calibrate_null_crc_from_seqs` returns one scalar `f64` from 64-bp windows, and `horizon_eval.rs` hardcodes `align_contig_pair_with_config` with no external MAF scorer; refuting C5 holds |
| F6 | UPHELD | Dense 1-bp `FlatSeedIndex` takes tens of GB on 3.1 Gbp and rejects >4 Gbp, `lis_chain` is single-contig monotonic, and `sdust_mask` misses interspersed LTRs; refuting C6 holds |
