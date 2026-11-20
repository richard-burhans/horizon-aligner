# Round 5 Adjudication of Completeness (`completeness-critic.md`)

**Target:** `reports/rounds/round-5/completeness-critic.md` (unsealed; quoted per finding below).
**Verification:** Checked directly against `src/seeding/mod.rs`, `src/seeding/flat_index.rs`, `src/formats.rs`, `src/hal/mod.rs`, `src/mdl/dag.rs`, `manuscript/horizon_paper.md`, `manuscript/horizon_companion.md`, `reports/rounds/round-5/round.json`, and all 7 Round 5 `attack-*.md` / `adjudication-*.md` files.

---

## Context Audit (Sections 1–2 of `completeness-critic.md`)

- **Lens execution (Section 1):** Verified. All 7 declared lenses (`correctness`, `inference`, `prior-art`, `provenance`, `reproduction`, `sampling`, `statistics`) filed attacks and were adjudicated (`0` `UNRUN` lenses).
- **Framing & claim-list critique (Section 2):** Verified in `reports/rounds/round-5/round.json:8-18, 57-65, 119-120`: `runner_is_author: true`, `claim_list_author: "runner"`, `C5` bundles 7 distinct specifications under a `BandGrid`-only `source.quote`, and `C1` bundles unit tests, `horizon_gfa stats` vs. `vg`, `horizon_clean`, `horizon_resolve`, and `horizon_wga --graph`.

---

## Finding-by-Finding Rulings

### Gap 1 — `KP-14` Eisenstein field `F_{q^2}` (`q = 2^61 - 1`), `(F + λ)(R + λ)` invariance, and `MASK_W12`/`MASK_W10` specs (`C5`)

> **Quoted critic claim:** *"Every lens skipped `C5`'s `"KP-14 Eisenstein seeding"` claim. `ATTACKED-horizon_companion.md:89-95` claims that ... Horizon encodes nucleotides over the quadratic extension field F_{q²} (q = 2⁶¹ − 1) using sixth-root Eisenstein units ... with `(F + λ)(R + λ)` invariance, and `line 101` lists `MASK_W12` (weight 12, span 15) and `MASK_W10` (weight 10, span 13). Direct inspection of `src/seeding/mod.rs:197-202` and `src/seeding/flat_index.rs:67-70` shows **zero** Eisenstein field arithmetic ... `encode_spaced_pair` computes two 2-bit packed `u64` integers and canonicalizes via `fwd.min(rc)` ... `src/seeding/mod.rs:9, 13` defines `MASK_PRIMARY_11_15` (weight **11**, span 15—not `MASK_W12` weight 12) and `MASK_W10` (`"11011011011011"`, span **14**, not span 13)."*

- **Verdict:** **UPHELD**
- **Diagnosis:** Verified against the codebase and all 7 Round 5 lenses. Only `prior-art` (`F5`) mentioned `KP-14` Eisenstein encoding (for an uncited reference `sajjad_2026_eisenstein-dna-bch`); no lens checked whether Eisenstein field arithmetic exists in `src/seeding/`. `grep -rn -i "eisenstein\|mersenne\|2\^61\|extension field" src/` returns zero matches. `src/seeding/mod.rs:197-202` (`canonical_spaced_seeds_impl`) and `src/seeding/flat_index.rs:67-73` (`FlatSeedIndex`) compute 2-bit packed `u64` pairs via `encode_spaced_pair` and take `fwd.min(rc)`. Furthermore, `src/seeding/mod.rs:9, 13` defines `MASK_PRIMARY_11_15 = "111011011011011"` (weight 11, span 15) and `MASK_W10 = "11011011011011"` (weight 10, span 14)—neither `MASK_W12` nor a span-13 `MASK_W10` exists. (Additionally, `horizon_companion.md:99` describes `di_count_01 >= len * 3 / 8` in `is_low_complexity_window`, whereas `src/seeding/mod.rs:31-32` checks `distinct <= 2` for `sub.len() >= 14`.)
- **Prescribed repair:** Endorsed. Replacing `manuscript/horizon_companion.md:87-101` (and `round.json:60`) with the actual `fwd.min(rc)` 2-bit packed `u64` canonicalization and true mask constants (`MASK_PRIMARY_11_15` weight 11 / span 15; `MASK_W10` weight 10 / span 14) requires no re-measurement; **the fix is free**.

---

### Gap 2 — Bidirectional `GFA ↔ MAF` Lens Laws (`GETPUT`, `PUTGET`, `Constant Complement`) (`C5`)

> **Quoted critic claim:** *"`prior-art` (`F5`) noted that `Foster et al. (2007)` was uncited in the companion body, but **no** lens (`correctness`, `inference`, `reproduction`, `provenance`) checked whether `src/formats.rs` actually implements `get(C)`, `putback(A', C)`, or the three lens laws (`GETPUT`, `PUTGET`, `Constant Complement`). `grep -rn "putback\|GETPUT\|PUTGET" src/` returns **zero matches**; `src/formats.rs` has no `get` or `putback` functions and no GFA-to-MAF round-trip logic."*

- **Verdict:** **UPHELD**
- **Diagnosis:** Verified. No lens in Round 5 checked `src/formats.rs` for bidirectional lens functions or laws. Direct inspection of `src/formats.rs` confirms zero `get`, `putback`, `GETPUT`, `PUTGET`, or `Constant Complement` implementations and no GFA-to-MAF round-trip logic. Both `manuscript/horizon_companion.md:21, 53-56` and `manuscript/horizon_paper.md:18, 22` (`"lossless MAF, rGFA, and HAL-compatible interoperability"`) overclaim unimplemented bidirectional round-tripping.
- **Prescribed repair:** Endorsed on the text-removal/qualification option: removing or marking as unimplemented future work the bidirectional lens claims in `manuscript/horizon_companion.md:21, 53-56` and qualifying `"lossless"` in `manuscript/horizon_paper.md:18, 22` requires no code or benchmark run; **the fix is free**. (The alternative—implementing and unit-testing `get`/`putback` in `src/formats.rs`—is a non-free engineering task.)

---

### Gap 3 — Pure-Rust HDF5 `Read + Seek` & `HAL` binary record decoders (`C5`)

> **Quoted critic claim:** *"Only `inference` (`F3`, checking `inspect_hal_genomes` `--hal-info`) and `provenance` (`F5`, noting `src/hal.rs` vs. `src/hal/` in `generate_figures.py:586`) touched `src/hal/`. No lens audited `src/hal/hdf5_read.rs`, `src/hal/hdf5_write.rs`, or `decode_hal_*` (`src/hal/mod.rs:332-425`) against the binary layout table in `ATTACKED-horizon_companion.md:66-81`, and the only `.hal` integration test on `HEAD` (`92e345c`) is `#[ignore]`d (`provenance F1`)."*

- **Verdict:** **UPHELD**
- **Diagnosis:** Verified as a genuine Round 5 coverage gap. Grep across all 7 `attack-*.md` and `adjudication-*.md` files confirms zero references to `hdf5_read`, `hdf5_write`, or `decode_hal_*`. (Note that `src/hal/mod.rs:386` documents `decode_hal_top_array`'s first `i64` field as `start`, whereas `manuscript/horizon_companion.md:80` labels it `genome_idx: i64`, illustrating why the table was unverified.) Per the adjudicator rules, this is an absence of Round 5 verification coverage rather than a refutation of the decoders themselves.
- **Prescribed repair:** Endorsed. Auditing `src/hal/hdf5_read.rs`, `src/hal/hdf5_write.rs`, and `src/hal/mod.rs:332-425` against `manuscript/horizon_companion.md:66-81` and running the `#[ignore]`d `.hal` integration test on a real external `.hal` fixture is **not free** (requires an external `.hal` fixture and code/test execution).

---

### Gap 4 — Unverified `34.2%` prototype Python runtime & unexamined `horizon_wga --graph` (`src/mdl/dag.rs`)

> **Quoted critic claim:** *"1. Both manuscripts state that pointer-heavy sparse dictionaries *"accounted for `34.2%` of prototype Python runtime"* (`ATTACKED-horizon_paper.md:36`, `ATTACKED-horizon_companion.md:139`) with no source script, commit, or benchmark artifact, and zero lens coverage. 2. While `inference` (`F7`) noted `horizon_wga --graph` was omitted from `C1`–`C6` evaluation and `provenance` (`F1`) noted stale Figure 6A numbers, no lens (`correctness`, `reproduction`) inspected `src/mdl/dag.rs` (`DagBandGrid`, `align_seq_to_dag_seeded`) or reproduced its `2.98–6.09 s` / `153.9–154.6 MB` run on the `10k` HPRC subgraph."*

- **Verdict:** **UPHELD**
- **Diagnosis:** Both parts are verified. First, `grep -rn "34\.2%" docs/ scripts/ src/ manuscript/` returns hits only in `manuscript/horizon_paper.md:36` and `manuscript/horizon_companion.md:139`, with zero provenance in `docs/`, `scripts/`, or `src/` and zero coverage by any Round 5 lens. Second, while `reproduction` (`F2`) noted that `DagBandGrid` exists at `src/mdl/dag.rs:11`, `inference` (`F7`) flagged the 3-query `--graph` evaluation, and `provenance` (`F1`) flagged stale Figure 6A numbers, neither `correctness` nor `reproduction` inspected `src/mdl/dag.rs` or re-ran `horizon_wga --graph` on the `10k` HPRC subgraph.
- **Prescribed repair:** Endorsed with split cost: (a) removing the ungrounded `34.2%` claim from `manuscript/horizon_paper.md:36` and `manuscript/horizon_companion.md:139` is a **free** text fix; (b) inspecting `src/mdl/dag.rs` and reproducing `horizon_wga --graph` on the `10k` subgraph is **not free** (requires code audit and benchmark execution).

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| Gap 1 | UPHELD | `src/seeding/mod.rs` and `flat_index.rs` use `fwd.min(rc)` 2-bit `u64` keys with `MASK_PRIMARY_11_15` (w=11, span 15) and `MASK_W10` (span 14), zero Eisenstein `F_{q^2}` arithmetic; text fix in `horizon_companion.md:87-101` is free |
| Gap 2 | UPHELD | `src/formats.rs` has no `get`, `putback`, `GETPUT`, `PUTGET`, or `Constant Complement` logic and no lens checked it; text removal/qualification in `horizon_companion.md:21,53-56` and `horizon_paper.md:18,22` is free |
| Gap 3 | UPHELD | No Round 5 lens audited `src/hal/hdf5_read.rs`, `hdf5_write.rs`, or `decode_hal_*` (`src/hal/mod.rs:332-425`) or ran the `#[ignore]`d `.hal` test; repair requires external `.hal` fixture and audit (not free) |
| Gap 4 | UPHELD | `34.2%` Python runtime has no provenance in `docs/`/`scripts/`/`src/` (removing it is free), and `src/mdl/dag.rs` / `horizon_wga --graph` (`2.98-6.09 s`) was neither audited by `correctness` nor reproduced by `reproduction` (not free) |
