# Round 8 ADJUDICATION — `provenance`

- **Attack reviewed:** `reports/rounds/round-8/attack-provenance.md` (sealed digest `77d810004c3e`, verified via `sha256sum`).
- **Target surface:** `docs/COMMUNITY_UTILITIES_PROPOSAL.md` (`C1`–`C6`).
- **Prior-round adjudications (Rule 5):** `scripts/round history "" --kind adjudication` confirms no prior closed round adjudicated `COMMUNITY_UTILITIES_PROPOSAL.md`, `horizon hal`, or `gfa-qc`.
- **Method cards (Rule 6):** This project keeps no method cards (`refs/methods/` is absent).

---

### F1 — C1 (`horizon hal`)
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; proposed correction **UPHELD** (text repair in `docs/COMMUNITY_UTILITIES_PROPOSAL.md` is **free**). While `4-t2t-orangs-mc-2023v2.full.hal` (`5.6G`) and `jarvis_lab_birds.hal` (`57G`) exist on host at `/home/rcb112/project/claude/bench/data/horizon_hal_bench/`, neither filename appears anywhere in `src/`, `docs/`, `manuscript/`, `scripts/`, `tests/`, or `git log` prior to `0f31c26`, and `src/hal/tests.rs:240-280` only tests `SEQNAME_ARRAY` decoding on `tests/data/halTest.hal`. Moreover, `src/hal/hdf5_read.rs:9, 683-816` buffers entire datasets up to `MAX_CHUNK_DECOMPRESSED_BYTES` (`256 MiB`) in `Vec<u8>` (not streaming), no HDF5-to-`HalAlignment` loader exists in `src/hal/`, and `HalAlignment::to_maf_for_genome` (`src/hal/mod.rs:416-470`) is an in-memory single-level projector using concatenated `root.name` coordinates.
- **What survives:** `hickey_2013_hierarchical`, `liao_2023_draft`, and `armstrong_2020_progressive` accurately motivate `horizon hal`, and `src/hal/hdf5_read.rs` + `src/coord.rs` exist as low-level building blocks rather than a verified multi-GiB streaming HAL-to-MAF pipeline.

---

### F2 — C2 (`horizon gfa-qc`)
- **Verdict:** **UPHELD**
- **Ruling:** Survival ruling **UPHELD**; optional caveat note is **free**. `docs/BENCHMARK_RESULTS.md:24, 64` verifies `96.9 MB` peak RSS in `196.93s` for `horizon_gfa stats` (`src/bin/horizon_gfa.rs:20-56`) on `hprc-v1.0-mc-grch38.gfa.gz`, and `docs/PANGENOME_SUITE.md:3` confirms `stream_gfa` skips `C`/`J` records. The attacker's provenance notes—that `horizon_gfa stats` is a stateless counter-only pass (so `gfa-qc` component/depth tracking will require additional state) and that `15–22%` combines `garrison_2024_building` Sup E text with rendered Sup Fig 2—are accurate.
- **What survives:** C2's literature motivation and `~97 MB` counter-only `horizon_gfa stats` RSS baseline survive intact.

---

### F3 — C3 (`horizon inversions`)
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; proposed correction **UPHELD** (text repair in `docs/COMMUNITY_UTILITIES_PROPOSAL.md` is **free**). In `src/bin/horizon_resolve.rs:210-268, 326-360, 483-527`, `tarjan_bridges` and `compute_biedged_snarls` only count bridges/components and mark `is_anchor` booleans, while `decompose_and_resolve` slices `P`/`W` walks into 1D windows between anchors (`j - i < max_bubble_steps`) without extracting 2-edge-cut snarls or allele subgraphs. `FlatSeedIndex` (`src/seeding/flat_index.rs:11-29`) has no built-in banded extension, and `sena_2026_bubble-subgraphs-spqr.md:7, 12-14, 34-37, 51-54` is an unreviewed theory preprint giving a linear-size *representation* of all snarls (whose count can still be quadratic) and linear-time enumeration of *ultrabubbles*, not linear-time enumeration of all snarls.
- **What survives:** `romain_2026_inversions-pangenome-graphs`'s `path-explicit` vs `alignment-rescued` inversion topology motivation survives; `horizon_resolve.rs` contributes `PackedGraph` + bridge-based anchor flagging rather than full snarl extraction.

---

### F4 — C4 (`horizon gfa-prune`)
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; proposed correction **UPHELD** (text repair in `docs/COMMUNITY_UTILITIES_PROPOSAL.md` is **free**). `liao_2023_draft.md:70-88, 124-125` confirms `vcfbub` (`-l 0 -r 100000` / `-r 10000000`) filtered `vg deconstruct` VCFs whereas graph-level >10 Mb deletion removal used `vg clip -D 10000000`, and `src/bin/horizon_clean.rs:151-168, 410-466` implements `parse_cigar_overlap_bp` and walk-splitting without walk-bridging.
- **What survives:** C4's literature motivation and `horizon_clean.rs`'s two-pass walk-splitting and CIGAR-overlap primitives survive; walk-bridging remains to be built.

---

### F5 — C5 (`horizon eval-null`)
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; proposed correction **UPHELD** (text repair in `docs/COMMUNITY_UTILITIES_PROPOSAL.md` is **free**). `quinlan-hall_2010_bedtools.md:49-80` is a two-page interval-overlap applications note that never discusses WGA ground truth (only `shuffleBed` for BED permutation enrichment), and `calibrate_null_crc_from_seqs` (`src/mdl/dp.rs:202-260`) samples `WIN = 64`-bp reversed `y_rev` windows at a single estimated GC fraction to return one scalar `f64` threshold, without dinucleotide shuffles or GC-stratified distributions.
- **What survives:** C5's WGA ground-truth motivation (`earl_2014_alignathon`, `frith_2011_tantan`, `garrison_2024_building`, `kopalli_2025_toward`) and the scalar reversed-window CRC calibration primitive survive.

---

### F6 — C6 (`horizon synteny`)
- **Verdict:** **UPHELD**
- **Ruling:** Diagnosis **UPHELD**; proposed correction **UPHELD** (text repair in `docs/COMMUNITY_UTILITIES_PROPOSAL.md` is **free**). `shivakumar_2025_mumemto.md`, `liu_2018_inferring.md`, `FlatSeedIndex` (`12 bytes/entry`), `lis_chain` (`src/seeding/chain.rs:9-42`), and `sdust_mask` (`src/seeding/mod.rs:54`) are all verified, and `FlatSeedIndex::build_with_sdust` (`src/seeding/flat_index.rs:41-44`) explicitly enforces a single-sequence `u32::MAX` (`4 Gbp`) limit.
- **What survives:** C6's literature motivation and low-RAM chaining/masking primitives survive, scoped to per-reference/per-contig indexing under `4 Gbp`.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | Repo lacks any `4-t2t-orangs`/`jarvis_lab_birds` verification log, HDF5->`HalAlignment` loader, or streaming dataset reader (`256 MiB` cap), and `to_maf_for_genome` is in-memory single-level; text fix is free |
| F2 | UPHELD | All 5 citations and `96.9 MB` `horizon_gfa stats` RSS are verified, with accurate caveats that `stats` is stateless counter-only and `stream_gfa` skips `C`/`J` |
| F3 | UPHELD | `horizon_resolve.rs` only counts Tarjan bridges and slices 1D anchor windows (no snarl extraction), `FlatSeedIndex` lacks extension, and `sena_2026` is a theory preprint without linear-time all-snarl enumeration; text fix is free |
| F4 | UPHELD | `liao_2023_draft` ran `vcfbub` on VCFs and `vg clip -D 10000000` on graphs, and `horizon_clean.rs` splits walks without bridging them; text fix is free |
| F5 | UPHELD | `quinlan-hall_2010_bedtools` does not discuss WGA ground truth and `calibrate_null_crc_from_seqs` returns a scalar `64`-bp reversed-window threshold; text fix is free |
| F6 | UPHELD | Literature and `12 bytes/entry` / `lis_chain` / `sdust_mask` are verified, and `FlatSeedIndex::build_with_sdust` enforces a hard `4 Gbp` `u32` single-sequence limit; text fix is free |
