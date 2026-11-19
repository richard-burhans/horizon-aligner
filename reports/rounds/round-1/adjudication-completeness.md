# Round 1 ADJUDICATION -- Completeness Critic

**Target filing:** `reports/rounds/round-1/completeness-critic.md` (unsealed; key passages quoted below per `ADJ-BRIEF-completeness.md`).
**Prior rounds & method cards:** `scripts/round history` confirms 0 prior closed rounds. No `refs/methods/` directory exists in this repository; rulings are grounded in direct verification of `src/`, `docs/`, `reports/rounds/round-1/round.json`, `CRITIC-BRIEF.md`, and `adjudication-L1.md` through `adjudication-L7.md`.

---

## Finding-by-Finding Rulings

### S1-S2 -- Execution Status, `CRITIC-BRIEF.md:51` Discrepancy, and Per-Lens Coverage
> *"Zero lenses were UNRUN... `CRITIC-BRIEF.md:51` states `attack-L5-inference.md (inference) — NO ADJUDICATION FILED`, but `reports/rounds/round-1/adjudication-L5.md` (102 lines) exists on disk... Zero claims were UNCHALLENGED across the round... L2 (`sampling`) is the thinnest lens (5 findings; 7 of 11 claims unchallenged)... L3 (`statistics`) is the second thinnest (7 findings; 5 claims unchallenged: C1, C5, C6, C8, C11)."*

- **Verdict:** **UPHELD**
- **What survives:** All 7 declared lenses filed both an attack and an adjudication, `CRITIC-BRIEF.md:51` indeed reflects stale state (`NO ADJUDICATION FILED` despite `adjudication-L5.md` being 102 lines on disk), and every per-lens finding count, challenged/unchallenged claim list, and verdict tally in the Section 2 table matches `adjudication-L1.md` through `adjudication-L7.md` exactly.
- **Repair:** Refreshing `round.json` / `CRITIC-BRIEF.md` registration for `adjudication-L5.md` is endorsed; **the fix is free** (metadata/text-only update).

---

### S3 -- Unfalsifiable Framing of Future-Tense Commitments and Normative Rules
> *"No codebase or benchmark evidence can falsify a future intention. As a result, `adjudication-L1.md`, `adjudication-L2.md`, and `adjudication-L4.md` repeatedly had to NARROW `REFUTED` findings on these claims—not because the technical criticisms were wrong, but because the claims were unfalsifiable promises rather than testable present-state or architectural-feasibility claims."*

- **Verdict:** **NARROWED**
- **What survives:** `C5`, `C6`, `C7`, `C8`, and `PANGENOME_SUITE.md:34` (`C9`–`C11`) are indeed phrased as future-tense commitments, normative rules, or a strategic roadmap, which forced `adjudication-L1.md` (F1, F5–F11), `adjudication-L2.md` (F1), and `adjudication-L4.md` (F1, F5–F7, F9–F11) to narrow existence-based refutations.
- **Overturned leg:** The critic's blanket assertion that *"No codebase or benchmark evidence can falsify a future intention"* overreaches. Technical-feasibility and internal-consistency attacks *did* falsify parts of these commitments (e.g., `L5` on `C6` dual-input `stdin`, `C7` single-pass arbitrary-order GFA stitching, and `C9` stateless `L`-line scoring; `L4` F8 on `C8`'s uncited present-tense HPRC crash premise).

---

### Gap 1 -- Omitted Claim: `.fai`/`.gzi`-Style Sidecar Indexing (`ATTACKED-STRATEGIC_PLAN.md:17`)
> *"`round.json:46` truncated the second sentence of `C5`: 'Any necessary indexing must be a lightweight, optional sidecar file (like `.fai` or `.gzi`), never a locked central database.' ... no lens attacked the technical claim that `.fai`/`.gzi`-style sidecars suffice for graph operations."*

- **Verdict:** **UPHELD**
- **What survives:** Verified in `round.json:44-51` and `ATTACKED-STRATEGIC_PLAN.md:17`, and verified across `attack-L1-correctness.md` through `attack-L7-reproduction.md`: while `L5` (F5) noted the contradiction between `stdin`/`stdout` and optional sidecars and `L6` (F4) cited `gfaidx`/`GBZ-base`, no lens attacked the technical sufficiency of linear-offset `.fai`/`.gzi`-style sidecars for 2D `L`-line topology, bifurcation seeding, or `W`/`P` walk indexing.
- **Repair:** Endorsed (specify an explicit graph-index format or drop the `.fai`/`.gzi` analogy in `docs/STRATEGIC_PLAN.md:17`). **The fix is free** (text-only edit).

---

### Gap 2 -- Omitted Claim: `MASK_RY20` and `AA5` Indexing on GFA `S` Lines (`ATTACKED-PANGENOME_SUITE.md:15`)
> *"`C10` (`round.json:91-94`) quoted only line 16... omitting line 15... No lens examined `MASK_RY20` or AA masks in `horizon_wga --graph`. In `src/seeding/mod.rs:13`, `MASK_RY20` is a 25-bp spaced seed (`\"1110110110111011011011011\"`), and `find_translated_aa5_hits` (`src/seeding/translated.rs:72`) requires 5-codon (15-bp) windows."*

- **Verdict:** **NARROWED**
- **What survives (Gap UPHELD):** Verified in `round.json:89-96` and `ATTACKED-PANGENOME_SUITE.md:15`: `C10.source.quote` omitted line 15, and zero `attack-L*.md` files examined `MASK_RY20` or `AA5` seeding on fragmented GFA `S` lines. `src/seeding/mod.rs:13` indeed defines `MASK_RY20` as a 25-bp mask, so indexing `S` lines individually without traversing `L` junctions drops seeds spanning short variant nodes.
- **Overturned detail:** `find_translated_aa5_hits` (`src/seeding/translated.rs:72`) calls `index_aa_spaced_4of6_frames` (`src/seeding/translated.rs:27-41`), which uses a spaced 5-of-6 mask (`"111011"`, **18-bp** span, `while pos + 18 <= seq.len()`), not 15-bp windows (15 bp appears only in the doc comment on `extend_aa_ungapped_hsp` at line 93). This strengthens rather than weakens the biological gap.
- **Repair:** Endorsed (specify path/k-mer traversal across `L` edges in `docs/PANGENOME_SUITE.md:15` and mark unimplemented). **The fix is free** (text-only edit).

---

### Gap 3 -- Omitted Claim: "Horizon Bubble" Emission in `horizon_clean` (`ATTACKED-PANGENOME_SUITE.md:9`)
> *"`C9` (`round.json:82-85`) quoted only line 8, omitting line 9: 'If a path through the graph does not yield sufficient bit-savings, explicitly sever the edge and emit a Horizon Bubble.' No lens examined what a 'Horizon Bubble' is in `src/formats.rs:15, 50`... Severing an `L` line in an external pangenome GFA and emitting a `*` `HB:i:1` node breaks existing `W`/`P` haplotype walks that traverse the severed link."*

- **Verdict:** **UPHELD**
- **What survives:** Verified in `round.json:80-87`, `ATTACKED-PANGENOME_SUITE.md:9`, and `src/formats.rs:15, 50`: `HB:i:1` is a `*`-sequence placeholder `S` segment for pairwise unaligned gaps. While `L5` (F5) briefly noted under `C5` that `src/formats.rs:50` emits `*`-sequence segments that cannot be aligned against, no lens examined `horizon_clean` (`C9`) severing `L` lines in external GFAs or the resulting breakage of `W`/`P` haplotype walks.
- **Repair:** Endorsed (define `W`/`P` walk/path handling when `L` edges are severed and clarify `HB:i:1` in `src/formats.rs` as a pairwise unaligned-gap placeholder). **The fix is free** (text-only edit).

---

### Gap 4 -- Omitted Claim: Graceful Abstention vs. Scale-Dependent Recall Collapse (`ATTACKED-PANGENOME_SUITE.md:17`)
> *"no lens audited whether `horizon_wga`'s abstention mechanism distinguishes true divergence abstention from scale-dependent recall collapse (recall drops 4.5x from `45.3%` on 38 kb `chr22_random` to `9.9%` on `chr22 1MB`)."*

- **Verdict:** **UPHELD**
- **What survives:** Verified in `docs/BENCHMARK_RESULTS.md:36-37` (`cactus241way` `chr22_random`: `45.3%` recall, `233` Horizon Bubbles; `chr22 1MB`: `9.9%` recall, `494` Horizon Bubbles—a `4.58x` drop). `L2` (F5) and `L3` (F1) cited the `9.9%` figure and `L5` (F2) noted MDL abstains on random sequence, but no lens audited why recall collapses `4.5x` from 38 kb to 1 MB on the same clade (`Human` vs `Dog`) or whether "graceful abstention" masks scale-dependent failure.
- **Repair:** Endorsed on diagnosis; split on cost. Documenting the `9.9%` collapse and qualifying "graceful abstention" in `docs/BENCHMARK_RESULTS.md` and `docs/PANGENOME_SUITE.md` **is free**, whereas auditing `horizon_wga` seeding/chaining on `cactus241way` `chr22 1MB` requires code/benchmark investigation.

---

### Gap 5 -- Omitted Claim: Strict Structural-Path Annotation Liftover (`ATTACKED-PANGENOME_SUITE.md:30`)
> *"`C11` (`round.json:100-103`) quoted line 29... omitting line 30: 'Project functional annotations from a reference path (e.g., `GRCh38`) to any alternate haplotype strictly via high-confidence structural paths.' No lens attacked the biological failure mode of strict structural-path projection..."*

- **Verdict:** **UPHELD**
- **What survives:** Verified in `round.json:98-105` and `ATTACKED-PANGENOME_SUITE.md:30`: although `L6` (F11) cited `Liftoff` as prior art, no lens attacked the biological failure mode of projecting annotations strictly along shared structural paths (dropping or fragmenting exons/genes in haplotype-specific insertions, gene conversions, or structural rearrangements).
- **Repair:** Endorsed (document fallback alignment/splice handling or scope limitation for non-syntenic haplotype regions in `docs/PANGENOME_SUITE.md:30`). **The fix is free** (text-only edit).

---

### Gap 6a -- Contradiction: `pub mod gfa;` Declaration in `src/lib.rs` (`adjudication-L5.md` vs. `L1`/`L7`)
> *"`adjudication-L5.md:34` asserts: 'Moreover, `gfa` is not even declared as a module in `src/lib.rs:3-9`' (repeated at line 62: 'Neither `coord` nor `gfa` is declared in `src/lib.rs:3-9`'). In fact, `src/lib.rs:394` does declare `pub mod gfa;`..."*

- **Verdict:** **UPHELD**
- **What survives:** Verified directly in `adjudication-L5.md:34, 62`, `src/lib.rs:394` (`pub mod gfa;` appears at the end of the file after `mod tests`), `adjudication-L1.md:22`, and `adjudication-L7.md:10`. (`gfa` still has zero callers in `src/`, so `L5`'s F4/F8 verdicts stand, but the factual misstatement in `adjudication-L5.md` must be corrected.)
- **Repair:** Endorsed (correct `adjudication-L5.md:34, 62` to cite `src/lib.rs:394` and note `gfa` is declared but unused). **The fix is free** (text-only edit).

---

### Gap 6b -- Contradiction: Verdict Standard on Roadmap/Future-Tense Claims (`L1`/`L2`/`L4` vs. `L5`/`L6`/`L7`)
> *"`adjudication-L1.md` (F1, F5–F11), `adjudication-L2.md` (F1 on C8), and `adjudication-L4.md` (F1, F5–F7, F9–F11) NARROWED `REFUTED` findings because `STRATEGIC_PLAN.md` uses future-tense/normative phrasing and `PANGENOME_SUITE.md:34` labels itself a strategic roadmap. `adjudication-L5.md` (F1, F6–F11), `adjudication-L6.md` (F5–F11), and `adjudication-L7.md` (F1, F5–F11) UPHELD `REFUTED` / `NARROWED` verdicts on those same claims without applying the roadmap/future-tense narrowing rule."*

- **Verdict:** **UPHELD**
- **What survives:** Verified across all six adjudication files: `L1`, `L2`, and `L4` systematically narrowed existence-based refutations on future-tense/roadmap claims, whereas `L7` (F1, F5–F11) and `L5`/`L6` (F8) upheld `REFUTED` verdicts on the same existence checks without applying that rule (note that `L5`'s `UPHELD` rulings on F1, F6, F7, F9–F11 target technical infeasibility rather than mere unimplemented status).
- **Repair:** Endorsed (standardize how roadmap-tense vs. present-tense claims are scored when closing Round 1). **The fix is free** (protocol/scoring edit).

---

### Gap 6c -- Contradiction: Standard vs. Non-Standard GFA Tags/Comments (`L1`/`L5` vs. `L6` on `C5`)
> *"`adjudication-L1.md:26` and `adjudication-L5.md:41` rule that `# EG:Z:` comments, `*` sequences with `LN:i:`, and `HB:i:1` / `XS:i:` / `YS:i:` / `ST:Z:` tags emitted by `src/formats.rs:15, 43, 50, 102` are valid standard GFA v1.0/v1.1 syntax. `adjudication-L6.md:33-34` upholds the characterization that Horizon 'emits non-standard GFA tags/comments'."*

- **Verdict:** **UPHELD**
- **What survives:** Verified in `adjudication-L1.md:26`, `adjudication-L5.md:41`, and `adjudication-L6.md:33-34`. Under GFA v1.0/v1.1, `#` comment lines, `*` sequences with `LN:i:`, and two-letter `TAG:TYPE:VALUE` optional fields are valid standard GFA syntax (`L1`/`L5` are correct; `L6` conflated Horizon-specific tag names with non-standard syntax).
- **Repair:** Endorsed (clarify in Round 1 closeout that `src/formats.rs` tags/comments are valid GFA v1.0/v1.1 syntax with Horizon-specific semantics). **The fix is free** (text-only edit).

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| S1-S2 | UPHELD | All 7 lenses ran and filed, `CRITIC-BRIEF.md:51` is stale on `adjudication-L5.md`, and all per-lens counts/verdicts match disk; metadata fix is free |
| S3 | NARROWED | Future-tense/roadmap phrasing forced L1/L2/L4 narrowing, but L4 (F8), L5, and L6 still falsified embedded premises and technical feasibility |
| Gap1 | UPHELD | `round.json:46` omitted `.fai`/`.gzi` sidecars and no lens attacked their technical insufficiency for 2D graph/walk indexing; doc fix is free |
| Gap2 | NARROWED | No lens examined `MASK_RY20`/`AA5` on short GFA `S` lines, though `find_translated_aa5_hits` requires 18-bp (6-codon) spans, not 15-bp; doc fix is free |
| Gap3 | UPHELD | No lens examined `horizon_clean` (`C9`) severing `L` lines and emitting `*`/`HB:i:1` segments that break `W`/`P` walks; doc fix is free |
| Gap4 | UPHELD | No lens audited why `cactus241way` recall collapses 4.5x (`45.3%` at 38 kb to `9.9%` at 1 MB) vs. "graceful abstention"; doc fix is free, code audit is not |
| Gap5 | UPHELD | `round.json:100-103` omitted strict structural-path projection (`PANGENOME_SUITE.md:30`) and no lens attacked its exon/rearrangement failure mode; doc fix is free |
| Gap6a | UPHELD | `adjudication-L5.md:34, 62` falsely states `gfa` is not declared in `src/lib.rs`, whereas `pub mod gfa;` is at `src/lib.rs:394`; text fix is free |
| Gap6b | UPHELD | L1/L2/L4 narrowed existence attacks on future-tense/roadmap claims while L7 (F1, F5–F11) and L5/L6 (F8) upheld them; scoring standardization is free |
| Gap6c | UPHELD | L1/L5 correctly rule `# EG:Z:` and `HB:i:`/`XS:i:`/`YS:i:`/`ST:Z:` valid GFA v1.0/v1.1 syntax while L6 calls them non-standard; closeout clarification is free |
