# Round 8 Attack — Reproduction Lens (`horizon-aligner`)

**Method-card & prior-round check:** Per the project rubric, `refs/methods/` keeps no method cards. Prior rounds were searched via `./scripts/round history`: Round 3 (`adjudication-L1.md`, `adjudication-L2.md`, `adjudication-L3.md`) already adjudicated `FlatSeedIndex` (`12 bytes/entry` post-build vs. `28 bytes/entry` peak during `build`, and `max_occ = 2` censoring on `>= 100 kb` contigs); Round 5 (`adjudication-completeness.md` Gap 3) flagged that `src/hal/hdf5_read.rs` and `src/hal/mod.rs` had not been audited at scale; Round 6 (`adjudication-L5.md`) established that `horizon_resolve.rs` computes 2-edge-connected components and linear `P`/`W` windows rather than full snarl/ultrabubble trees. All host experiments below were run on this machine (`30 GiB` RAM, `0B` swap).

---

## F1 -- C1 (`horizon hal`): benchmark HAL files are absent on host, `hdf5_read.rs` aborts above 256 MiB, and `to_maf_for_genome` is not a tree-wide projector

- **Claim:** C1 (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:41-42`) — `src/hal/hdf5_read.rs` is already verified against `4-t2t-orangs-mc-2023v2.full.hal` (5.5 GiB) and `jarvis_lab_birds.hal` (56.4 GiB), and `horizon hal` (`stats`, `to-maf`, `liftover`) can reuse `hdf5_read.rs` and `HalAlignment::to_maf_for_genome` to match `halTools` on those datasets.
- **Verdict:** **REFUTED**
- **Concrete flaw:**
  1. Neither `4-t2t-orangs-mc-2023v2.full.hal` nor `jarvis_lab_birds.hal` is present on the host (`find /home/rcb112 -maxdepth 4 -name "*.hal"` returned empty; `/home/rcb112/project/gemini/bench/pg/data/benchmarks/sars-cov-2.hal` is `196` bytes and fails `./target/release/horizon_wga --hal-info` with `not a valid HDF5 HAL file`). The only real-file test in `src/hal/tests.rs:241-243` (`real_hal_decode_from_env`) defaults to the `55 KiB` fixture `tests/data/halTest.hal` (all 17 `hal::tests` pass in `0.01s`).
  2. `read_hdf5_v0_dataset_bytes_seek` (`src/hal/hdf5_read.rs:9, 683-816`) hard-caps every dataset at `MAX_CHUNK_DECOMPRESSED_BYTES = 256 * 1024 * 1024` (`256 MiB`), returning `Ok(None)` if a contiguous or chunked dataset exceeds `256 MiB` (lines `729, 767, 785, 795, 800`), and has **no** chunk-level random-access slice reader (it decompresses every chunk in the B-tree into a single `Vec<u8>`). In a 5.5–56.4 GiB vertebrate/ape HAL, a single genome's 4-bit nibble-packed `DNA_ARRAY` (~600 MiB–1.6 GiB) or `TOP_ARRAY` (`33 bytes/seg`) exceeds `256 MiB` and fails with `Ok(None)`.
  3. `HalAlignment::to_maf_for_genome` (`src/hal/mod.rs:416-495`) only projects onto direct children (`bseg.child_tops`) or immediate parent (`write_leaf_maf`), does not traverse multi-branch phylogenetic trees, and requires loading all genomes' `dna`, `top_segments`, and `bottom_segments` into RAM simultaneously (infeasible for `56.4 GiB` `jarvis_lab_birds.hal` on `30 GiB` RAM).
- **Evidence:** `src/hal/hdf5_read.rs:6-9, 683-816`; `src/hal/mod.rs:416-495`; `src/hal/tests.rs:241-276`; `cargo test --lib hal::tests` (`17 passed` in `0.01s`); `./target/release/horizon_wga --hal-info .../sars-cov-2.hal`.
- **Correction:** Remove the claim that `4-t2t-orangs-mc-2023v2.full.hal` and `jarvis_lab_birds.hal` are reproducibly verified on this codebase, and explicitly scope `horizon hal` as requiring a new streaming/chunk-indexed HDF5 slice reader (removing the `256 MiB` `MAX_CHUNK_DECOMPRESSED_BYTES` cap) and tree-traversal MAF projector rather than reusing `read_hdf5_v0_dataset_bytes_seek` and `HalAlignment::to_maf_for_genome` as-is.

---

## F2 -- C2 (`horizon gfa-qc`): `~97 MB` RSS only holds for `horizon_gfa stats` scalar counters, not connected components or node-depth/AF histograms

- **Claim:** C2 (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:56-59`) — `horizon gfa-qc` can perform a single-pass connected-component, chromosome-bridge, and node-depth/AF audit by reusing `src/bin/horizon_gfa.rs` (`~97 MB` peak RSS on the `12 GB` HPRC graph).
- **Verdict:** **NARROWED**
- **Concrete flaw:** `horizon_gfa stats` (`src/bin/horizon_gfa.rs:20-56`) achieves `96.9 MB` peak RSS (`196.93s` on `81.4M` segments / `113.0M` links in `docs/BENCHMARK_RESULTS.md:64`) **only because** it increments five `u64` scalar counters (`s_count`, `s_bases`, `l_count`, `p_count`, `w_count`) with zero segment-ID interning or graph state. Computing connected components (Union-Find over `81.4M` segments) and per-segment depth/AF histograms in a single pass requires interning `81.4M` string segment names (`HashMap<String, u32>` ≈ `4.5–6 GB` alone) plus Union-Find/length/depth/haplotype arrays (`~1.3+ GB`), totaling **`~6–8 GB` RSS** (60–80× higher than `96.9 MB`).
- **Evidence:** `src/bin/horizon_gfa.rs:20-56`; `docs/BENCHMARK_RESULTS.md:64`.
- **Correction:** Narrow the `~97 MB` RSS baseline to counter-only `horizon_gfa stats` and specify an explicit `~6–10 GB` RAM budget (or integer-ID compaction / two-pass external counting) for `horizon gfa-qc` on `81.4M`-segment HPRC graphs.

---

## F3 -- C3 (`horizon inversions`): `PackedGraph` OOMs on whole-genome HPRC GFA and `horizon_resolve.rs` has neither snarl enumeration nor sequence storage

- **Claim:** C3 (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:72-75`) — `horizon inversions` can build `PackedGraph` + Tarjan snarls (reusing `horizon_resolve.rs`) and `FlatSeedIndex` to call path-explicit and alignment-rescued inversions on pangenome graphs.
- **Verdict:** **REFUTED**
- **Concrete flaw:**
  1. Measuring `horizon_resolve` directly on `/tmp/hprc_scaling_cache/hprc_sub_200000.gfa` (`200,000` segments, `278,464` links, `0` paths, `0` walks) gave **`95.7 MB` peak RSS** (`~478 bytes/segment`; `docs/BENCHMARK_RESULTS.md:87` reports `95.5 MB`). Scaling `407×` to the full `81.4M`-segment HPRC graph requires **`~38.9 GB` RAM** even with **zero** `P`/`W` lines—exceeding the host's `30 GiB` RAM—and `pass1_index` (`src/bin/horizon_resolve.rs:552, 570, 623-624`) additionally buffers all un-deduplicated walk/path step pairs in `gray_edges: Vec<(u32, u32)>` before sorting.
  2. Prior Round 6 (`adjudication-L5.md:61-62`) already established that `compute_biedged_snarls` (`src/bin/horizon_resolve.rs:272-374`) only computes 2-edge-connected components and `is_anchor: Vec<bool>`, then `decompose_and_resolve` walks linear `P`/`W` windows (`max_bubble_steps = 16`). `PackedGraph` (`src/bin/horizon_resolve.rs:170-177`) stores no snarl/ultrabubble tree and no segment DNA sequences for `alignment_rescued`.
  3. On `hprc_chr22_sub_real.gfa` (`10,000` segments, `1` path, `4` walks, `749` windows), `horizon_resolve` took `2.881s` (`3.8 ms/window`).
- **Evidence:** `src/bin/horizon_resolve.rs:170-198, 272-374, 558-636`; `docs/BENCHMARK_RESULTS.md:87`; host runs of `./target/release/horizon_resolve` on `hprc_sub_200000.gfa` (`95.7 MB` RSS, `0.237s`) and `hprc_chr22_sub_real.gfa` (`25.7 MB` RSS, `2.881s`).
- **Correction:** Do not claim `PackedGraph` in `horizon_resolve.rs` scales to whole-genome HPRC GFAs in low RAM or provides snarl/ultrabubble enumeration; scope `horizon inversions` either to pre-extracted locus/chromosome subgraphs or to a new bit-packed/external-memory graph representation with actual snarl-tree construction.

---

## F4 -- C4 (`horizon gfa-prune`): `horizon_clean.rs` OOMs on whole-genome HPRC GFA and splits walks instead of bridging

- **Claim:** C4 (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:85-87`) — `horizon gfa-prune` provides streaming AF-filtering and deletion-bubble capping, rewriting affected `W`/`P` walks cleanly using the walk-splitting and CIGAR-overlap logic in `horizon_clean.rs`.
- **Verdict:** **REFUTED**
- **Concrete flaw:**
  1. Measuring `horizon_clean` directly on `/tmp/hprc_scaling_cache/hprc_sub_200000.gfa` (`200,000` segments, `278,464` links, `0` paths, `0` walks) gave **`86.9 MB` peak RSS** (`docs/BENCHMARK_RESULTS.md:87` reports `130.6 MB`). Scaling `407×` to `81.4M` segments requires **`~35.4–53.2 GB` RAM** before `walk_support: HashMap<DirectedEdge, usize>` (`src/bin/horizon_clean.rs:174`) is even populated from `2,580` `P` and `24,456` `W` records—exceeding the host's `30 GiB` RAM.
  2. `emit_clean_walk` (`src/bin/horizon_clean.rs:389-466`) **splits** walks into disjoint `W` chunks whenever an edge is severed (`chunks.push(...)`). It cannot bridge across pruned unsupported segments or deleted >100 kb shortcut edges, so reusing it would fragment haplotype walks.
- **Evidence:** `src/bin/horizon_clean.rs:170-259, 389-466`; `docs/BENCHMARK_RESULTS.md:87`; host run of `./target/release/horizon_clean --gfa /tmp/hprc_scaling_cache/hprc_sub_200000.gfa --out /dev/null` (`86.9 MB` RSS, `0.24s`).
- **Correction:** Replace the `horizon_clean.rs` streaming claim with an explicit multi-GB memory requirement (or external/disk-backed edge table) and note that bridging pruned bubbles requires new graph-traversal walk rewriting rather than `emit_clean_walk`'s walk-splitting logic.

---

## F5 -- C5 (`horizon eval-null`): `calibrate_null_crc_from_seqs` is lightweight, but `horizon_eval.rs` materializes every aligned base pair in a `HashSet`

- **Claim:** C5 (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:97-99`) — `horizon eval-null` can run `--n-shuffles 20` and score decoy/concordance MAF/PAF outputs alongside `horizon_eval.rs` recall/precision.
- **Verdict:** **NARROWED**
- **Concrete flaw:** `calibrate_null_crc_from_seqs` (`src/mdl/dp.rs:202-287`) samples `n_cal` 64-bp windows (`WIN = 64`) in bounded memory and reproduces cleanly. However, `horizon_eval.rs` (`src/bin/horizon_eval.rs:42, 46-107, 193-199`) materializes every aligned base pair into `MafPairMap = HashMap<(String, String), HashSet<(usize, usize)>>` (~24–32+ bytes/bp), which requires ~4–5 GB for a single human chromosome (`chr6`, `170.8 Mbp`) and >70 GB for whole-genome mammalian alignments, making `--n-shuffles 20` on whole-genome outputs infeasible via `horizon_eval.rs`.
- **Evidence:** `src/mdl/dp.rs:202-287`; `src/bin/horizon_eval.rs:42, 67-74, 193-199`.
- **Correction:** Narrow `horizon eval-null` to interval/block-level MAF/PAF overlap and `calibrate_null_crc_from_seqs` window sampling, rather than reusing `horizon_eval.rs`'s per-base `HashSet<(usize, usize)>` representation at whole-genome scale.

---

## F6 -- C6 (`horizon synteny`): `FlatSeedIndex` hard-errors above 4 Gbp, peaks at ~93 GB on one human genome, and weight-11 `max_occ = 2` saturates

- **Claim:** C6 (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:109-111`) — `horizon synteny` can index `<ref.fa> <assemblies...>` with `FlatSeedIndex` (`12 bytes/entry`) + SDUST + LIS chaining (`src/seeding/chain.rs`) using "a few GB of RAM rather than hundreds of GB".
- **Verdict:** **REFUTED**
- **Concrete flaw:**
  1. `FlatSeedIndex::build_with_sdust` (`src/seeding/flat_index.rs:41-45, 57`) hard-errors above `u32::MAX` (`4 Gbp`) per sequence and indexes **every** 1-bp position (`for i in 0..=(seq.len() - span)` — no minimizers). As Round 3 (`adjudication-L1.md:30-33`, `adjudication-L3.md:27-32`) already established, `FlatSeedIndex::build` peaks at **`28 bytes/entry`** (`16 B/entry` `pairs` + `12 B/entry` `keys`/`pos`), plus `seq` (`1 B/bp`) and `sdust` (`1 B/bp`) in `prepare_contig_index_opts` (`src/seeding/mod.rs:260-288`, which also allocates `seq_rc` = `1 B/bp`). Indexing a single `3.1 Gbp` human genome requires **`~93 GB` peak RAM** (`~43 GB` post-build), not "a few GB of RAM".
  2. `prepare_contig_index_opts` (`src/seeding/mod.rs:266-272`) uses `PRIMARY_SEED_WEIGHT = 11` (`2,097,152` canonical keys) with `max_occ = 2` on contigs `>= 100 kb` (adjudicated in Round 3 `adjudication-L2.md:22-26`). On a `250 Mbp` chromosome or `3.1 Gbp` genome, random Poisson multiplicity is ~119–1,478× per key, so `max_occ = 2` drops nearly all seeds, while raising `max_occ` causes $O(\text{occ}_x \times \text{occ}_y)$ blowup in `FlatSeedIndex::intersect` (`src/seeding/flat_index.rs:242-249`).
- **Evidence:** `src/seeding/flat_index.rs:41-95, 242-249`; `src/seeding/mod.rs:260-288`; Round 3 `adjudication-L1.md:30-33`, `adjudication-L2.md:22-26`, `adjudication-L3.md:27-32`.
- **Correction:** Withdraw the "few GB of RAM" claim for dense `FlatSeedIndex` on `<ref.fa> <assemblies...>`; `horizon synteny` requires minimizer/syncmer sampling and a larger k-mer weight (e.g., `k >= 19–21`) to fit whole-genome multi-assembly indexing in a few GB of RAM.

---

## Claim Verdicts

| Claim | Verdict | Summary |
|---|---|---|
| C1 | REFUTED | Benchmark HALs absent on host; `hdf5_read.rs` caps datasets at 256 MiB with no slice reader; `to_maf_for_genome` only projects direct parent/children. |
| C2 | NARROWED | `~97 MB` RSS only holds for counter-only `horizon_gfa stats`; single-pass components + depth/AF require ~6–8 GB on `81.4M` HPRC segments. |
| C3 | REFUTED | `horizon_resolve` (`95.7 MB` at `200k` segs) projects to `~38.9 GB` at `81.4M` segs (OOM on 30 GiB host) and lacks snarl-tree/sequence support. |
| C4 | REFUTED | `horizon_clean` (`86.9 MB` at `200k` segs) projects to `~35.4+ GB` at `81.4M` segs (OOM on 30 GiB host) and splits walks instead of bridging. |
| C5 | NARROWED | `calibrate_null_crc_from_seqs` is lightweight, but `horizon_eval.rs` materializes per-base `HashSet<(usize, usize)>` (unscalable to whole genomes). |
| C6 | REFUTED | Dense 1-bp `FlatSeedIndex` hard-caps at 4 Gbp, peaks at `~93 GB` on 3.1 Gbp, and weight-11 `max_occ = 2` saturates at chromosome/genome scale. |
