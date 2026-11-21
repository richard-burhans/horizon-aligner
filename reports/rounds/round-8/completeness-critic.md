# Round 8 Completeness Critic — `docs/COMMUNITY_UTILITIES_PROPOSAL.md`

**Files read:** `reports/rounds/round-8/CRITIC-BRIEF.md`, `reports/rounds/round-8/round.json`, `docs/COMMUNITY_UTILITIES_PROPOSAL.md`, `src/bin/horizon.rs`, and all eight round-8 artifacts (`attack-correctness.md`, `adjudication-correctness.md`, `attack-prior-art.md`, `adjudication-prior.md`, `attack-provenance.md`, `adjudication-provenance.md`, `attack-reproduction.md`, `adjudication-reproduction.md`).

---

## 1. Lens Coverage: UNRUN vs. UNCHALLENGED, and Filing Counts

Per Rule 1, **UNRUN** and **UNCHALLENGED** are kept strictly separate:

- **UNRUN lenses (3):** `sampling`, `statistics`, `inference` — explicitly dropped in `round.json:89-96` (*"Proposal triage round: evaluating engineering feasibility, prior art, source-grounded motivation, and benchmark criteria for 6 proposed CLI utilities rather than a statistical survey"*). Nobody ran them.
- **UNCHALLENGED claims (0 across the 4 run lenses):** Every claim `C1`–`C6` was attacked by all 4 declared lenses (`correctness`, `prior-art`, `provenance`, `reproduction`). Only `C2` (`horizon gfa-qc`) received a `SURVIVED` verdict in one lens (`provenance`, `adjudication-provenance.md:17-20`), while being `NARROWED` by `correctness` and `reproduction` and `REFUTED` by `prior-art`.
- **Carried items (`round.json:101-134`):** `R7-SWEEPGA`, `R7-CRC-BIO`, and `R7-HAL-ANC` were explicitly dropped from Round 8 as out of scope for proposal triage and remain open for their respective future rounds.

### Per-Lens Filing & Adjudication Counts

| Lens | Findings Filed | Claims Covered | Attack Verdicts | Adjudication Rulings | Thin? |
|---|---|---|---|---|---|
| `correctness` | 6 (`F1`–`F6`) | `C1`–`C6` (1 each) | 2 `NARROWED`, 4 `REFUTED` | 6 `UPHELD` | No |
| `prior-art` | 6 (`F1`–`F6`) | `C1`–`C6` (1 each) | 3 `NARROWED`, 3 `REFUTED` | 6 `UPHELD` | No |
| `provenance` | 6 (`F1`–`F6`) | `C1`–`C6` (1 each) | 1 `SURVIVED`, 3 `NARROWED`, 2 `REFUTED` | 6 `UPHELD` | No |
| `reproduction` | 6 (`F1`–`F6`) | `C1`–`C6` (1 each) | 2 `NARROWED`, 4 `REFUTED` | 5 `UPHELD`, 1 `NARROWED` (`F1`) | No |

All four run lenses filed 6 substantive, source-grounded findings (24 total). None of the run lenses is thin on `C1`–`C6` individually; the gaps below are **framing, claim-structure, and cross-cutting constraint** blind spots that *all four lenses* missed.

---

## 2. Claims Nobody Could Attack (Claim-List Defects)

### Gap 1 — Compound ("triple-barrel") claims + undefined `P1`/`P2`/`P3` priority criteria (`round.json:6-61`, `docs/COMMUNITY_UTILITIES_PROPOSAL.md:19-29, 115-120`)
- **What went unattacked:** Every claim `C1`–`C6` in `round.json` bundles three distinct assertions into one sentence: (1) *fills a real community gap*, (2) *can be built on specific existing `horizon-aligner` primitives*, and (3) *at `P1`/`P2`/`P3` priority*. No lens attacked the **Recommended Implementation Order** or the **`P1`/`P2`/`P3` priority ranking** (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:19-29, 115-120`), because "at P1 priority" has no falsifiable criterion (no effort estimate, user-demand metric, or build-vs-drop threshold).
- **Why it matters:** Because the claims were compound, `correctness` and `reproduction` refuted `C3`–`C6` on *primitive-reuse* grounds (e.g., `PackedGraph` lacks a snarl tree, `FlatSeedIndex` is `pub(crate)` and 1-bp dense), leaving unsettled whether the *underlying utility* (e.g., an AF pruner or inversion caller built on *new* primitives) is still worth building. Conversely, `prior-art` refuted `C2` (`horizon gfa-qc`, P1) as duplicating `odgi`+`panacus`+`gfatools`, while `correctness` and `reproduction` kept `C2` as "worth building as P1" and `provenance` marked it `SURVIVED`. Without separating community value, primitive feasibility, and priority criteria, the round's headline question (*"which proposed community utilities ... are worth building"*) remains ambiguous for `C2`, `C3`, and `C4`.
- **What would close it:** Split each utility in `docs/COMMUNITY_UTILITIES_PROPOSAL.md` into separate claims for (a) *unmet community need*, (b) *reusable vs. new engineering scope*, and (c) *explicit triage decision (Build / Re-scope / Drop)* with effort-vs-impact justification.

---

## 3. Framing, Lens-Drop, and Missing-Constraint Gaps

### Gap 2 — Dropping `statistics` / `inference` left statistical validity in `C5` (`horizon eval-null`) unexamined (`round.json:89-96`, `docs/COMMUNITY_UTILITIES_PROPOSAL.md:15, 98-99`)
- **What went unexamined:** `statistics` and `inference` were dropped for the whole round as *"not a statistical survey"*. While `C5` was refuted by `correctness`, `prior-art`, and `reproduction` on code-reuse and prior-art grounds, no lens evaluated the statistical assumptions of `C5`'s proposed methodology—specifically, whether reversed-sequence decoys inflate false positives on palindromic/tandem-repeat DNA (where a sequence and its reverse share local similarity even after `sdust` or `tantan`), or whether `calibrate_null_crc_from_seqs`'s Conformal Risk Control (CRC) exchangeability assumption holds across GC-stratified genomic bins.
- **Why it matters:** Even if `horizon eval-null` is dropped as a standalone community CLI, `calibrate_null_crc_from_seqs` remains an active internal primitive (tracked in `R7-CRC-BIO`, `round.json:109-113`). Treating statistical methodology as out-of-scope for `C5` meant neither `C5` nor its underlying CRC null-calibration assumptions were audited by a statistical lens.
- **What would close it:** Confirm `C5` is dropped from the community roadmap (so no user-facing statistical claims ship), and ensure `R7-CRC-BIO` explicitly includes a `statistics`/`inference` check on decoy exchangeability and GC-bin sample sizes.

### Gap 3 — Unchecked overlap and fragmentation with the existing `horizon` CLI (`src/bin/horizon.rs:7-63`, `src/bin/horizon_lift/`, `docs/COMMUNITY_UTILITIES_PROPOSAL.md:21-29, 34-60, 79-87`)
- **What every lens missed:** `src/bin/horizon.rs` already dispatches 11 subcommands (`wga`, `merge`, `guide`, `gfa`, `lift`, `clean`, `resolve`, `filter`, `locus`, `eval`, `bench-chr1`). None of the four lenses checked `src/bin/horizon.rs` or `src/bin/horizon_lift/`:
  1. `C2` (`horizon gfa-qc`) and `C4` (`horizon gfa-prune`) propose new top-level subcommands (`gfa-qc`, `gfa-prune`) alongside existing `horizon gfa` (`horizon_gfa`: `stats` and `extract`) and `horizon clean` (`horizon_clean`), rather than sub-actions under `horizon gfa` / `horizon clean`.
  2. `C1` (`horizon hal liftover`) proposes a new liftover command without specifying its relationship to the existing `horizon lift` (`horizon_lift`: *"Graph-aware BED/GFF3 coordinate liftover with micro-bubble coalescing"*).
- **Why it matters:** Shipping redundant top-level subcommands (`gfa-qc` vs. `gfa`, `gfa-prune` vs. `clean`, `hal liftover` vs. `lift`) fragments the `horizon` CLI UX and risks duplicating CLI parsing, I/O, and maintenance across parallel binaries.
- **What would close it:** Map any surviving utilities (`C1`, and re-scoped `C2`/`C4` if kept) onto the existing `src/bin/horizon.rs` subcommand hierarchy (e.g., `horizon gfa qc`, `horizon hal ...`) and clarify how `horizon hal liftover` relates to `horizon lift`.

### Gap 4 — Missing GFA dialect and HAL-to-GFA interoperability constraints (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:39-40, 51-59, 70-75, 81-87`, `attack-provenance.md:33`)
- **What every lens missed:** `provenance` F2 noted in passing that `stream_gfa` skips `C`/`J` records (`attack-provenance.md:33`), but no lens systematically checked GFA dialect and output-schema compatibility constraints:
  1. **Blunt vs. non-zero CIGAR overlaps:** Minigraph-Cactus emits blunt GFAs (`0M` links), whereas PGGB/minigraph GFAs frequently contain non-zero link overlaps. `C2` (`gfa-qc` retained-bp and depth) and `C4` (`gfa-prune` deletion-bubble span) give no specification for overlap-adjusted lengths, risking silent miscounts on non-blunt GFAs.
  2. **GFA v1.0 `P`-only vs. v1.1 `W`-only vs. rGFA (`SN`/`SO`/`SR`):** `C2`–`C4` assume `W`/`P` paths without specifying behavior on rGFA graphs lacking `W`/`P` lines or mixed `P`+`W` graphs (where counting both double-counts haplotypes).
  3. **`horizon hal to-gfa` (`C1`):** No lens checked what GFA v1.1 conventions (`W`-line sample/haplotype/seqid formatting, `--no-ancestors` semantics, blunt segment splitting) downstream `vg`/`odgi`/`panacus` consumers require.
- **Why it matters:** A community utility that silently miscounts overlap bases, double-counts `P`+`W` haplotypes, or emits non-standard `W` lines will fail on real PGGB/minigraph/vg workflows even if its Rust implementation is memory-efficient.
- **What would close it:** Add explicit GFA-dialect support matrices (blunt vs. overlap CIGARs, `P`/`W`/`J`/`C`/rGFA handling) and `vg`/`odgi` schema-compatibility acceptance criteria for any surviving GFA or `hal to-gfa` utility.

### Gap 5 — Missing opportunity-cost / unproposed community utilities check (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:1-120`, `src/bin/horizon.rs`)
- **What every lens missed:** All four lenses evaluated only the 6 *proposed* tools (`C1`–`C6`) and never asked whether higher-leverage community utilities exist in the cited literature or on `horizon-aligner`'s *actual* verified primitives. After Round 8's attacks, `C3`–`C6` are refuted, `C2` is refuted by `prior-art` and narrowed by `correctness`/`reproduction`, and `C1` is narrowed to static-binary packaging convenience requiring a new HDF5 slice reader and tree-wide MAF projector. Yet `horizon-aligner` already has verified primitives (e.g., `horizon_merge` `O(N log N)` MAF+rGFA merger, `horizon_lift` micro-bubble-coalescing liftover, `horizon_gfa` 97 MB streaming stats/extract) that were not compared against `C1`–`C6` for community ROI.
- **Why it matters:** Without checking alternatives, the team risks spending weeks building a new HDF5 slice reader for `C1` (which duplicates Bioconda `halTools`) or re-architecting `C2`–`C4` (which duplicate `odgi`/`panacus`/`vg clip`/`INVPG_annot`) instead of hardening or exposing higher-leverage existing capabilities.
- **What would close it:** Add a brief "Build / Re-scope / Drop" synthesis and alternative-utility comparison to `docs/COMMUNITY_UTILITIES_PROPOSAL.md` before scheduling implementation work.

---

## 4. Summary of Actionable Gaps

1. **Compound claims & unattacked `P1`/`P2`/`P3` priority (`round.json:6-61`, `docs/COMMUNITY_UTILITIES_PROPOSAL.md:19-29, 115-120`):**
   - *Why it matters:* Conflates "wrong internal primitives" with "no community need" and leaves the `P1`/`P2`/`P3` roadmap arbitrary.
   - *Close:* Separate community-need, primitive-feasibility, and priority/triage claims per utility.
2. **Dropped `statistics`/`inference` lens on `C5` (`round.json:89-96`, `docs/COMMUNITY_UTILITIES_PROPOSAL.md:98-99`):**
   - *Why it matters:* Leaves reversed-sequence repeat-inflation and CRC GC-bin exchangeability unexamined while `R7-CRC-BIO` remains active internally.
   - *Close:* Keep `C5` dropped from the community roadmap and audit CRC null assumptions under `R7-CRC-BIO`.
3. **Unchecked `horizon` CLI overlap (`src/bin/horizon.rs:7-63`, `src/bin/horizon_lift/`):**
   - *Why it matters:* `gfa-qc`/`gfa-prune`/`hal liftover` fragment and duplicate existing `horizon gfa`, `horizon clean`, and `horizon lift` subcommands.
   - *Close:* Map any surviving utilities into the existing `horizon.rs` subcommand hierarchy.
4. **Missing GFA-dialect & `hal to-gfa` interoperability specs (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:39-40, 51-59, 70-75, 81-87`):**
   - *Why it matters:* Non-zero CIGAR overlaps, mixed `P`/`W` lines, rGFA, and `vg`/`odgi` `W`-line conventions can silently break community workflows.
   - *Close:* Specify dialect support matrices and `vg`/`odgi` compatibility tests for any surviving GFA/HAL tools.
5. **No comparison against unproposed / existing-primitive alternatives (`docs/COMMUNITY_UTILITIES_PROPOSAL.md:1-120`):**
   - *Why it matters:* With `C3`–`C6` refuted and `C1`/`C2` heavily narrowed, engineering effort could go to low-ROI wrappers duplicating `halTools`/`odgi`/`vg clip`.
   - *Close:* Add a Build / Re-scope / Drop decision table weighing surviving proposals against existing `horizon` tools.
