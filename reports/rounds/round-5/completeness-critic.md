# Round 5 Completeness Critique

**Reviewed artifacts:**
- `reports/rounds/round-5/CRITIC-BRIEF.md`, `reports/rounds/round-5/round.json`, `reports/rounds/round-5/CARRIED.md`
- `reports/rounds/round-5/ATTACKED-horizon_paper.md`, `reports/rounds/round-5/ATTACKED-horizon_companion.md`
- All 7 attack files (`attack-correctness.md`, `attack-inference.md`, `attack-prior-art.md`, `attack-provenance.md`, `attack-reproduction.md`, `attack-sampling.md`, `attack-statistics.md`) and all 7 adjudication files (`adjudication-correctness.md`, `adjudication-inference.md`, `adjudication-prior.md`, `adjudication-provenance.md`, `adjudication-reproduction.md`, `adjudication-sampling.md`, `adjudication-statistics.md`).

---

## 1. Lens Execution & Coverage Audit (Rules 1–2)

- **UNRUN lenses:** **0**. All 7 declared lenses (`correctness`, `inference`, `prior-art`, `provenance`, `reproduction`, `sampling`, `statistics`) filed complete attacks and were adjudicated.
- **Per-lens counts and thinness:**
  - `correctness`: **4 findings** (`F1`–`F4`), with **`C2` and `C4` `UNCHALLENGED`** (the `12 GB` HPRC `chr6` comparison script and Table 1 datasets were not re-run under this lens). `correctness` is the thinnest lens on empirical claims.
  - `inference`: **7 findings** (`F1`–`F7`, including meta-finding `F7`; `F5` on `C5` was `OVERTURNED` on code inspection and `NARROWED` on companion text).
  - `prior-art`: **7 findings** (`F1`–`F7`, including meta-finding `F7`, `NARROWED`).
  - `provenance`: **6 findings** (`F1`–`F6`; `F6` on `C6` `SURVIVED`).
  - `reproduction`: **6 findings** (`F1`–`F6`; `F6` on `C6` `SURVIVED`).
  - `sampling`: **6 findings** (`F1`–`F6`; `F5` on `C5` `SURVIVED`, `F6` on `C6` `UNCHALLENGED`).
  - `statistics`: **6 findings** (`F1`–`F6`, including meta-finding `F6`; `F5` marked `C5` and `C6` `UNCHALLENGED`).

Note on `UNCHALLENGED` vs. `UNRUN`: `C2` (`correctness`), `C4` (`correctness`), `C5` (`statistics`), and `C6` (`sampling`, `statistics`) are `UNCHALLENGED` (the lens inspected the claim and found no defect within its scope), not `UNRUN`.

---

## 2. Framing & Claim-List Critique (Rules 3–4)

`round.json:119-120` records `runner_is_author: true` and `claim_list_author: "runner"`. Three framing defects shaped what the 7 lenses examined:

1. **Overloaded composite `C5` with a `BandGrid`-only quote (`round.json:57-65`):** `C5`'s summary bundles **seven** distinct technical specifications (`KP-14 Eisenstein seeding`, `BOCPD beta-bernoulli`, `search_space_bits`, `BandGrid idx(i,j)`, `GFA<->MAF lens laws`, `HDF5/HAL binary record decoders`, and `22 bp worked toy example`), but its `source.quote` (`round.json:63`) points solely to `BandGrid` (`ATTACKED-horizon_companion.md:139`). As a result, all 7 lenses treated `C5` primarily as a `BandGrid` check, leaving `KP-14` Eisenstein seeding, `GFA ↔ MAF` lens laws, and `HDF5/HAL` decoders completely unexamined by `correctness`, `reproduction`, and `inference`.
2. **Overloaded composite `C1` (`round.json:8-18`):** `C1` bundles unit-test counts, `horizon_gfa stats` vs. `vg`, `horizon_clean`, `horizon_resolve`, and `horizon_wga --graph` into one claim. Every lens focused on `horizon_gfa stats` vs. `vg` (and `provenance`/`sampling` checked the `10k` subgraph), leaving `horizon_wga --graph` (`src/mdl/dag.rs`) and `horizon_clean` (`src/bin/horizon_clean.rs`) largely unexamined.
3. **Omitted high-risk manuscript assertions (`inference F7`, `statistics F6`):** `C1`–`C6` omitted Conformal Risk Control (`α = 0.07`, `n_cal = 15`, `λ̂ = 10.0`, `0.55`), `egraph` "ancestral order" without a phylogenetic tree (`ATTACKED-horizon_paper.md:44`), and the `34.2%` prototype Python runtime figure (`ATTACKED-horizon_paper.md:36`, `ATTACKED-horizon_companion.md:139`).

---

## 3. Concrete Unattacked Gaps

### Gap 1 — `KP-14` Eisenstein field `F_{q^2}` (`q = 2^61 - 1`), `(F + λ)(R + λ)` invariance, and `MASK_W12`/`MASK_W10` specs (`C5`)
- **Artifacts:** `reports/rounds/round-5/ATTACKED-horizon_companion.md:87-101`, `reports/rounds/round-5/round.json:60`, `src/seeding/mod.rs:9, 13, 197-202`, `src/seeding/flat_index.rs:67-70`.
- **What went unattacked:** Every lens skipped `C5`'s `"KP-14 Eisenstein seeding"` claim. `ATTACKED-horizon_companion.md:89-95` claims that *"To avoid canonicalizing forward and reverse-complement k-mers separately, Horizon encodes nucleotides over the quadratic extension field F_{q²} (q = 2⁶¹ − 1) using sixth-root Eisenstein units: φ(A) = 1, φ(T) = -1, φ(C) = ω, φ(G) = -ω"* with `(F + λ)(R + λ)` invariance, and `line 101` lists `MASK_W12` (weight 12, span 15) and `MASK_W10` (weight 10, span 13).
- **Why it matters:** Direct inspection of `src/seeding/mod.rs:197-202` and `src/seeding/flat_index.rs:67-70` shows **zero** Eisenstein field arithmetic, **zero** `2^61 - 1`, and **zero** `(F + λ)(R + λ)` polynomial—`encode_spaced_pair` computes two 2-bit packed `u64` integers and canonicalizes via `fwd.min(rc)` (the exact separate forward/reverse canonicalization the companion claims to avoid). Moreover, `src/seeding/mod.rs:9, 13` defines `MASK_PRIMARY_11_15` (weight **11**, span 15—not `MASK_W12` weight 12) and `MASK_W10` (`"11011011011011"`, span **14**, not span 13).
- **What would close it:** Replace `manuscript/horizon_companion.md:87-101` (and `round.json:60`) with the actual `fwd.min(rc)` 2-bit packed `u64` canonicalization in `src/seeding/mod.rs` / `flat_index.rs` and the true mask names/spans (`MASK_PRIMARY_11_15` weight 11 / span 15; `MASK_W10` weight 10 / span 14).

### Gap 2 — Bidirectional `GFA ↔ MAF` Lens Laws (`GETPUT`, `PUTGET`, `Constant Complement`) (`C5`)
- **Artifacts:** `reports/rounds/round-5/ATTACKED-horizon_companion.md:21, 53-56`, `reports/rounds/round-5/ATTACKED-horizon_paper.md:18, 22`, `reports/rounds/round-5/round.json:60`, `src/formats.rs`.
- **What went unattacked:** `prior-art` (`F5`) noted that `Foster et al. (2007)` was uncited in the companion body, but **no** lens (`correctness`, `inference`, `reproduction`, `provenance`) checked whether `src/formats.rs` actually implements `get(C)`, `putback(A', C)`, or the three lens laws (`GETPUT`, `PUTGET`, `Constant Complement`).
- **Why it matters:** `grep -rn "putback\|GETPUT\|PUTGET" src/` returns **zero matches**; `src/formats.rs` has no `get` or `putback` functions and no GFA-to-MAF round-trip logic. The companion's Section 2.2 lens specification (`ATTACKED-horizon_companion.md:53-56`) describes unimplemented functionality, and `ATTACKED-horizon_paper.md:18, 22` claims `"lossless MAF, rGFA, and HAL-compatible interoperability"` without a working bidirectional lens.
- **What would close it:** Either implement and unit-test `get` / `putback` satisfying `GETPUT`, `PUTGET`, and `Constant Complement` in `src/formats.rs`, or remove/mark as unimplemented future work the bidirectional lens laws in `manuscript/horizon_companion.md:21, 53-56` and qualify `manuscript/horizon_paper.md:18, 22`.

### Gap 3 — Pure-Rust HDF5 `Read + Seek` & `HAL` binary record decoders (`C5`)
- **Artifacts:** `reports/rounds/round-5/ATTACKED-horizon_companion.md:58-81`, `reports/rounds/round-5/ATTACKED-horizon_paper.md:53, 96`, `reports/rounds/round-5/round.json:60`, `src/hal/mod.rs:332-425`, `src/hal/hdf5_read.rs`, `src/hal/hdf5_write.rs`.
- **What went unattacked:** Only `inference` (`F3`, checking `inspect_hal_genomes` `--hal-info`) and `provenance` (`F5`, noting `src/hal.rs` vs. `src/hal/` in `generate_figures.py:586`) touched `src/hal/`. No lens audited `src/hal/hdf5_read.rs`, `src/hal/hdf5_write.rs`, or `decode_hal_*` (`src/hal/mod.rs:332-425`) against the binary layout table in `ATTACKED-horizon_companion.md:66-81`, and the only `.hal` integration test on `HEAD` (`92e345c`) is `#[ignore]`d (`provenance F1`).
- **Why it matters:** The companion's binary decoder widths (`DNA_ARRAY` 4-bit, `SEQIDX_ARRAY` 24B, `TOP_ARRAY` 33B, `BOTTOM_ARRAY` `24 + 9C` B) and the paper's readiness claim for multi-hundred-gigabyte `.hal` archives (`ATTACKED-horizon_paper.md:96`) remain unverified against real external `.hal` files.
- **What would close it:** Audit `src/hal/hdf5_read.rs`, `src/hal/hdf5_write.rs`, and `src/hal/mod.rs:332-425` against `manuscript/horizon_companion.md:66-81`, and run the environment-gated `.hal` test on a real `.hal` fixture.

### Gap 4 — Unverified `34.2%` prototype Python runtime & unexamined `horizon_wga --graph` (`src/mdl/dag.rs`)
- **Artifacts:** `reports/rounds/round-5/ATTACKED-horizon_paper.md:36, 51, 82`, `reports/rounds/round-5/ATTACKED-horizon_companion.md:139`, `src/mdl/dag.rs`.
- **What went unattacked:**
  1. Both manuscripts state that pointer-heavy sparse dictionaries *"accounted for `34.2%` of prototype Python runtime"* (`ATTACKED-horizon_paper.md:36`, `ATTACKED-horizon_companion.md:139`) with no source script, commit, or benchmark artifact, and zero lens coverage.
  2. While `inference` (`F7`) noted `horizon_wga --graph` was omitted from `C1`–`C6` evaluation and `provenance` (`F1`) noted stale Figure 6A numbers, no lens (`correctness`, `reproduction`) inspected `src/mdl/dag.rs` (`DagBandGrid`, `align_seq_to_dag_seeded`) or reproduced its `2.98–6.09 s` / `153.9–154.6 MB` run on the `10k` HPRC subgraph.
- **Why it matters:** The `34.2%` statistic is an ungrounded historical number, and `horizon_wga --graph` is a flagship abstract claim (`ATTACKED-horizon_paper.md:10`) whose implementation (`DagBandGrid` in `src/mdl/dag.rs`) was never verified by `correctness` or `reproduction`.
- **What would close it:** Cite or remove the `34.2%` Python prototype statistic, and have `correctness`/`reproduction` verify `src/mdl/dag.rs` and reproduce the `horizon_wga --graph` timing on the `10k` subgraph.

---

## Summary Action List

| Gap | Artifact | Why It Matters | What Would Close It |
| :--- | :--- | :--- | :--- |
| **1. `KP-14` Eisenstein field & mask specs (`C5`)** | `ATTACKED-horizon_companion.md:87-101`, `src/seeding/mod.rs`, `flat_index.rs` | Companion specifies `F_{q^2}` Eisenstein hashing & `MASK_W12`/`MASK_W10` (span 13); code actually uses `fwd.min(rc)` 2-bit `u64` canonicalization, `MASK_PRIMARY_11_15` (weight 11), and `MASK_W10` (span 14) | Rewrite `horizon_companion.md:87-101` to match `src/seeding/mod.rs` and `flat_index.rs` |
| **2. `GFA ↔ MAF` Lens Laws (`C5`)** | `ATTACKED-horizon_companion.md:21, 53-56`, `ATTACKED-horizon_paper.md:18, 22`, `src/formats.rs` | `get`/`putback`/`GETPUT`/`PUTGET`/`Constant Complement` do not exist anywhere in `src/formats.rs` | Implement and test `get`/`putback` in `src/formats.rs`, or remove/qualify the lens-law claims |
| **3. HDF5 `Seek` & `HAL` binary decoders (`C5`)** | `ATTACKED-horizon_companion.md:58-81`, `src/hal/mod.rs`, `hdf5_read.rs` | Binary record decoders and multi-GB `.hal` readiness (`horizon_paper.md:96`) were never audited or run on an external `.hal` file | Audit `src/hal/` against `horizon_companion.md:66-81` and run the `#[ignore]`d `.hal` integration test |
| **4. `34.2%` Python claim & `horizon_wga --graph` (`src/mdl/dag.rs`)** | `ATTACKED-horizon_paper.md:36, 82`, `horizon_companion.md:139`, `src/mdl/dag.rs` | `34.2%` has no provenance; `DagBandGrid` / `horizon_wga --graph` was not inspected or reproduced | Cite or remove `34.2%`; verify `src/mdl/dag.rs` and reproduce `horizon_wga --graph` |
