# Round 8 Adjudication — Reproduction Lens (`horizon-aligner`)

**Target:** `reports/rounds/round-8/attack-reproduction.md` (verified SHA-256 prefix `30c0d85ad544`).
**Method-card & prior-round check:** `refs/methods/` contains no method cards, so no statistical method card governs these findings. Prior adjudications via `./scripts/round history` confirm that Round 3 (`adjudication-L1.md`, `adjudication-L2.md`, `adjudication-L3.md`) already settled `FlatSeedIndex` peak memory (`28 bytes/entry` during `build`, `12 bytes/entry` steady-state) and `max_occ = 2` censoring on `>= 100 kb` contigs, and Round 6 (`adjudication-L5.md`) already settled that `horizon_resolve.rs` computes 2-edge-connected components and linear `P`/`W` windows rather than full snarl/ultrabubble trees.

---

## F1 -- C1 (`horizon hal`): **NARROWED**

- **What survives:** `read_hdf5_v0_dataset_bytes_seek` (`src/hal/hdf5_read.rs:9, 683-816`) hard-caps every contiguous or chunked dataset at `MAX_CHUNK_DECOMPRESSED_BYTES = 256 * 1024 * 1024` (`256 MiB`) with no chunk-level random-access slice reader, and `HalAlignment::to_maf_for_genome` (`src/hal/mod.rs:416-495`) only projects onto direct parent or children while requiring all genomes' arrays in RAM simultaneously.
- **What is overturned:** The attacker's first premise—that `4-t2t-orangs-mc-2023v2.full.hal` and `jarvis_lab_birds.hal` are absent on the host—is **FALSE**. Both files are present at `/home/rcb112/project/claude/bench/data/horizon_hal_bench/4-t2t-orangs-mc-2023v2.full.hal` (`5.94 GB`) and `/home/rcb112/project/claude/bench/data/horizon_hal_bench/jarvis_lab_birds.hal` (`60.60 GB`), and `./target/release/horizon_wga --hal-info` succeeds on both (`root_members=10` and `root_members=69`, respectively). The attacker's `find /home/rcb112 -maxdepth 4` stopped two directories too shallow; they should have searched without `-maxdepth 4` or checked `/home/rcb112/project/claude/bench/data/horizon_hal_bench/`.
- **Repair:** **ENDORSED** on the technical scoping fix (explicitly scoping `horizon hal` as requiring a streaming/chunk-indexed HDF5 slice reader that removes the `256 MiB` cap and a tree-traversal MAF projector rather than reusing `read_hdf5_v0_dataset_bytes_seek` and `HalAlignment::to_maf_for_genome` as-is), **OVERTURNED** on removing the benchmark files as "absent". Qualifying `docs/COMMUNITY_UTILITIES_PROPOSAL.md` is free (doc-only); implementing the slice reader and tree projector is non-free.

---

## F2 -- C2 (`horizon gfa-qc`): **UPHELD**

- **What survives:** Verified in `src/bin/horizon_gfa.rs:20-56` and `docs/BENCHMARK_RESULTS.md:64` that `horizon_gfa stats` achieves `96.9 MB` peak RSS on the `81.4M`-segment HPRC graph solely because it increments five `u64` scalar counters (`s_count`, `s_bases`, `l_count`, `p_count`, `w_count`) with zero segment-ID interning or graph state; single-pass connected components and node-depth/AF histograms require multi-GB RAM (`~6–8 GB` RSS) for `81.4M` segment-name interning and per-segment arrays.
- **Repair:** **ENDORSED — the fix is free.** Narrowing the `~97 MB` RSS baseline in `docs/COMMUNITY_UTILITIES_PROPOSAL.md:59` to counter-only `horizon_gfa stats` and specifying an explicit multi-GB RAM budget (or integer-ID compaction / two-pass external counting) for `horizon gfa-qc` is a documentation-only edit.

---

## F3 -- C3 (`horizon inversions`): **UPHELD**

- **What survives:** Verified in `src/bin/horizon_resolve.rs:170-198, 272-374, 558-636` and `docs/BENCHMARK_RESULTS.md:87` (`95.5 MB` at `200k` segments, projecting to `~38.9 GB` at `81.4M` segments, exceeding the `30 GiB` host RAM) that `PackedGraph` does not scale to whole-genome HPRC GFAs in low RAM, stores no segment DNA sequences for `alignment_rescued`, and `compute_biedged_snarls` computes 2-edge-connected components and `is_anchor` rather than a full snarl/ultrabubble tree (as previously adjudicated in Round 6 `adjudication-L5.md`).
- **Repair:** **ENDORSED — the fix is free** for the documentation qualification in `docs/COMMUNITY_UTILITIES_PROPOSAL.md:72-75` (scoping `horizon inversions` to pre-extracted locus/chromosome subgraphs or noting that a new bit-packed/external-memory graph representation and snarl-tree builder are required).

---

## F4 -- C4 (`horizon gfa-prune`): **UPHELD**

- **What survives:** Verified in `src/bin/horizon_clean.rs:170-259, 389-466` and `docs/BENCHMARK_RESULTS.md:87` (`130.6 MB` at `200k` segments, `86.9 MB` measured by the attacker, projecting to `~35–53 GB` at `81.4M` segments) that `collect_clean_pass1` loads the full segment/link/flank/support tables into `HashMap`s rather than streaming, and `emit_clean_walk` splits walks into disjoint `W` chunks (`chunks.push(...)`) on severed edges rather than bridging across pruned bubbles.
- **Repair:** **ENDORSED — the fix is free** for the documentation qualification in `docs/COMMUNITY_UTILITIES_PROPOSAL.md:85-87` (replacing the streaming claim with an explicit multi-GB / disk-backed requirement and noting that bridging requires new graph-traversal walk rewriting).

---

## F5 -- C5 (`horizon eval-null`): **UPHELD**

- **What survives:** Verified in `src/mdl/dp.rs:202-287` and `src/bin/horizon_eval.rs:42, 46-107, 193-199` that `calibrate_null_crc_from_seqs` samples `n_cal` 64-bp windows (`WIN = 64`) in bounded memory, whereas `horizon_eval.rs` materializes every aligned base pair in `MafPairMap = HashMap<(String, String), HashSet<(usize, usize)>>`, which does not scale to whole-genome `--n-shuffles 20` runs.
- **Repair:** **ENDORSED — the fix is free** for the documentation qualification in `docs/COMMUNITY_UTILITIES_PROPOSAL.md:97-99` (scoping `horizon eval-null` to interval/block-level MAF/PAF overlap and `calibrate_null_crc_from_seqs` window sampling).

---

## F6 -- C6 (`horizon synteny`): **UPHELD**

- **What survives:** Verified in `src/seeding/flat_index.rs:41-45, 57, 242-249` and `src/seeding/mod.rs:260-288` (and previously adjudicated in Round 3 `adjudication-L1.md`, `adjudication-L2.md`, `adjudication-L3.md`) that `FlatSeedIndex::build_with_sdust` hard-errors above `u32::MAX` (`4 Gbp`) per sequence, indexes every 1-bp position at `28 bytes/entry` peak during `build` (plus `seq`, `seq_rc`, and `sdust`, totaling `~93+ GB` peak on a `3.1 Gbp` human genome), and uses `PRIMARY_SEED_WEIGHT = 11` with `max_occ = 2` on `>= 100 kb` contigs, which saturates at chromosome/genome scale.
- **Repair:** **ENDORSED — the fix is free** for the documentation qualification in `docs/COMMUNITY_UTILITIES_PROPOSAL.md:109-111` (withdrawing the "few GB of RAM" claim for dense `FlatSeedIndex` and specifying minimizer/syncmer sampling with a larger k-mer weight).

---

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | `256 MiB` `hdf5_read.rs` cap and direct-parent/child-only `to_maf_for_genome` hold, but both benchmark HAL files ARE present on host and pass `--hal-info` |
| F2 | UPHELD | `~97 MB` RSS only holds for counter-only `horizon_gfa stats`; single-pass components and depth/AF histograms require multi-GB interning on `81.4M` segments |
| F3 | UPHELD | `PackedGraph` projects to `~38.9 GB` at `81.4M` segments, lacks sequence storage, and `compute_biedged_snarls` does not build a snarl tree |
| F4 | UPHELD | `horizon_clean.rs` loads full `HashMap` state in pass 1 (OOMs at `81.4M` segments) and `emit_clean_walk` splits walks instead of bridging |
| F5 | UPHELD | `calibrate_null_crc_from_seqs` is bounded-memory, but `horizon_eval.rs` materializes per-base `HashSet<(usize, usize)>` |
| F6 | UPHELD | Dense 1-bp `FlatSeedIndex` hard-caps at `4 Gbp`, peaks at `~93+ GB` on `3.1 Gbp`, and weight-11 `max_occ = 2` saturates at genome scale |
