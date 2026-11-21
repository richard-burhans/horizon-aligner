# Round 8 Adjudication of Completeness — `completeness-critic.md`

- **Critic filing reviewed:** [completeness-critic.md](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-8/completeness-critic.md)
- **Round artifacts checked:** [round.json](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-8/round.json), [COMMUNITY_UTILITIES_PROPOSAL.md](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/COMMUNITY_UTILITIES_PROPOSAL.md), [horizon.rs](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon.rs), `attack-*.md`, and `adjudication-*.md`.
- **Prior rounds (Rule 5):** `./scripts/round history "R7-CRC-BIO" --kind adjudication` confirms no prior closed-round adjudication of `R7-CRC-BIO` or `COMMUNITY_UTILITIES_PROPOSAL.md`.
- **Method cards (Rule 6):** No method cards exist in `refs/methods/`.

The critic's Section 1 lens-coverage accounting (3 `UNRUN` lenses: `sampling`, `statistics`, `inference`; 0 `UNCHALLENGED` claims across the 4 run lenses; 24 findings total) is verified accurate against `round.json:89-96` and the four `attack-*.md` / `adjudication-*.md` files. Rulings on the five filed gaps (`Gap 1`–`Gap 5`) follow.

---

## Finding-by-Finding Rulings

### `Gap 1` — Compound ("triple-barrel") claims + undefined `P1`/`P2`/`P3` priority criteria — **UPHELD**
> Critic's claim: *"Every claim `C1`–`C6` in `round.json` bundles three distinct assertions into one sentence: (1) fills a real community gap, (2) can be built on specific existing `horizon-aligner` primitives, and (3) at `P1`/`P2`/`P3` priority. No lens attacked the Recommended Implementation Order or the `P1`/`P2`/`P3` priority ranking..."*

- **Ruling:** **UPHELD** (both diagnosis and documentation-restructuring repair).
- **Verification:** Verified in [round.json](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-8/round.json#L6-L61) and [COMMUNITY_UTILITIES_PROPOSAL.md](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/COMMUNITY_UTILITIES_PROPOSAL.md#L19-L29). Across all four `attack-*.md` and `adjudication-*.md` files, no lens evaluated the `P1`/`P2`/`P3` priority criteria or the Recommended Implementation Order (`correctness` F1/F2 merely echoed `"as P1"`). Consequently, `C3`–`C6` were refuted on primitive-reuse grounds without settling whether any re-scoped version has unmet community value, and `C2` (`prior-art` REFUTED vs. `correctness`/`reproduction` NARROWED) remains ambiguous.
- **Repair cost:** Separating community need, reusable vs. new engineering scope, and explicit `Build / Re-scope / Drop` triage criteria in `docs/COMMUNITY_UTILITIES_PROPOSAL.md` is a documentation-only edit; **the fix is free**.

### `Gap 2` — Dropping `statistics` / `inference` left statistical validity in `C5` (`horizon eval-null`) unexamined — **NARROWED**
> Critic's claim: *"no lens evaluated the statistical assumptions of `C5`'s proposed methodology—specifically, whether reversed-sequence decoys inflate false positives on palindromic/tandem-repeat DNA ... or whether `calibrate_null_crc_from_seqs`'s Conformal Risk Control (CRC) exchangeability assumption holds across GC-stratified genomic bins. ... Even if `horizon eval-null` is dropped as a standalone community CLI, `calibrate_null_crc_from_seqs` remains an active internal primitive (tracked in `R7-CRC-BIO`, `round.json:109-113`)."*

- **Ruling:** **NARROWED** — **OVERTURNED** as a Round 8 defect on `C5`; **UPHELD** only as a forward-looking scope note for `R7-CRC-BIO`.
- **Verification:** `C5` was already **REFUTED** by `correctness`, `prior-art`, and `reproduction` (and narrowed/refuted by `provenance`), with all four adjudicators upholding the refutations. Dropping `statistics`/`inference` in Round 8 therefore did **not** allow any false statistical claim in `C5` to survive. Moreover, `round.json:109-113, 126-129` already tracks internal CRC biological calibration under `R7-CRC-BIO` and explicitly deferred it to dedicated empirical CRC calibration rounds.
- **What survives:** Noting that `R7-CRC-BIO` should include `statistics`/`inference` checks on reversed-decoy repeat inflation and GC-bin exchangeability is a useful forward-looking scope note (documentation-only; **the fix is free**).

### `Gap 3` — Unchecked overlap and fragmentation with the existing `horizon` CLI — **UPHELD**
> Critic's claim: *"`src/bin/horizon.rs` already dispatches 11 subcommands (`wga`, `merge`, `guide`, `gfa`, `lift`, `clean`, `resolve`, `filter`, `locus`, `eval`, `bench-chr1`). None of the four lenses checked `src/bin/horizon.rs` or `src/bin/horizon_lift/`..."*

- **Ruling:** **UPHELD** (diagnosis and documentation repair).
- **Verification:** Verified in [horizon.rs](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon.rs#L7-L63) and `src/bin/horizon_lift/` (`main.rs`, `lift.rs`, `features.rs`). Zero `attack-*.md` or `adjudication-*.md` files mention `horizon.rs` or `horizon_lift`. Proposing top-level `gfa-qc`, `gfa-prune`, and `hal liftover` alongside existing `horizon gfa`, `horizon clean`, and `horizon lift` creates genuine CLI fragmentation and maintenance duplication if any utility survives.
- **Repair cost:** Mapping any surviving utilities (`C1`, and re-scoped `C2` if kept) into the existing `horizon.rs` subcommand hierarchy in `docs/COMMUNITY_UTILITIES_PROPOSAL.md` is a documentation-only edit; **the fix is free**.

### `Gap 4` — Missing GFA dialect and HAL-to-GFA interoperability constraints — **UPHELD**
> Critic's claim: *"`provenance` F2 noted in passing that `stream_gfa` skips `C`/`J` records (`attack-provenance.md:33`), but no lens systematically checked GFA dialect and output-schema compatibility constraints: 1. Blunt vs. non-zero CIGAR overlaps ... 2. GFA v1.0 `P`-only vs. v1.1 `W`-only vs. rGFA ... 3. `horizon hal to-gfa` (`C1`)..."*

- **Ruling:** **UPHELD** (diagnosis and specification repair).
- **Verification:** Verified in [attack-provenance.md](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-8/attack-provenance.md#L33) and [COMMUNITY_UTILITIES_PROPOSAL.md](file:///home/rcb112/project/gemini/bench/horizon-aligner/docs/COMMUNITY_UTILITIES_PROPOSAL.md#L39-L87); no lens audited non-zero CIGAR overlaps, mixed `P`+`W` haplotype double-counting, rGFA-without-paths behavior, or `vg`/`odgi` `W`-line conventions for `hal to-gfa`. For any surviving GFA or `hal to-gfa` tool, these omissions cause silent miscounts or downstream tool failures on real PGGB/minigraph/vg graphs.
- **Repair cost:** Adding GFA-dialect support matrices and `vg`/`odgi` schema-compatibility acceptance criteria to `docs/COMMUNITY_UTILITIES_PROPOSAL.md` is **free** (doc-only; implementing the corresponding code and tests later is non-free).

### `Gap 5` — Missing opportunity-cost / unproposed community utilities check — **OVERTURNED**
> Critic's claim: *"All four lenses evaluated only the 6 proposed tools (`C1`–`C6`) and never asked whether higher-leverage community utilities exist in the cited literature or on `horizon-aligner`'s actual verified primitives. ... `horizon-aligner` already has verified primitives (e.g., `horizon_merge` `O(N log N)` MAF+rGFA merger, `horizon_lift` micro-bubble-coalescing liftover, `horizon_gfa` 97 MB streaming stats/extract) that were not compared against `C1`–`C6` for community ROI."*

- **Ruling:** **OVERTURNED** (a gap that is not a defect, and subsumed by `Gap 1` + `Gap 3`).
- **Verification:** Round 8's declared scope in [round.json](file:///home/rcb112/project/gemini/bench/horizon-aligner/reports/rounds/round-8/round.json#L4) is *"Evaluate which proposed community utilities in docs/COMMUNITY_UTILITIES_PROPOSAL.md are worth building"*. Faulting the four lenses for evaluating the six proposed utilities rather than inventing unproposed alternatives is outside the round's verification target. Furthermore, `horizon_merge`, `horizon_lift`, and `horizon_gfa` are already shipped subcommands in [horizon.rs](file:///home/rcb112/project/gemini/bench/horizon-aligner/src/bin/horizon.rs#L13-L32), not unproposed utilities, and the need for a post-round `Build / Re-scope / Drop` synthesis weighing surviving proposals against existing `horizon` subcommands is already covered by `Gap 1` and `Gap 3`.
- **What the critic should have checked instead (Rule 7):** Whether `Gap 1` (explicit `Build / Re-scope / Drop` triage) and `Gap 3` (reconciliation with existing `horizon.rs` subcommands) already capture the actionable decision gap without filing a separate open-ended "unproposed utilities" defect.

---

| finding | verdict | why |
|---|---|---|
| Gap 1 | UPHELD | Compound claims (`C1`–`C6`) and unattacked `P1`/`P2`/`P3` priority criteria leave `Build / Re-scope / Drop` triage ambiguous; doc restructuring is free |
| Gap 2 | NARROWED | Overturned as a Round 8 defect on `C5` (`C5` is already refuted and `R7-CRC-BIO` tracks empirical CRC separately); upheld only as a free scope note for `R7-CRC-BIO` |
| Gap 3 | UPHELD | No lens checked `src/bin/horizon.rs` (11 existing subcommands) or `src/bin/horizon_lift/` for `gfa-qc`/`gfa-prune`/`hal liftover` overlap; doc mapping is free |
| Gap 4 | UPHELD | Only a passing `C`/`J` note existed; non-zero CIGAR overlaps, mixed `P`+`W`/rGFA handling, and `hal to-gfa` `W`-line specs were unexamined; doc spec is free |
| Gap 5 | OVERTURNED | Evaluating the 6 proposed utilities matches `round.json` scope, cited tools (`horizon_merge`/`lift`/`gfa`) already exist, and triage synthesis is subsumed by `Gap 1` + `Gap 3` |
