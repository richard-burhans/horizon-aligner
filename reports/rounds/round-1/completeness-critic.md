# Round 1 Completeness Critic — Strategic Plan & Pangenome Suite Roadmap

**Target round:** `reports/rounds/round-1/` (`round.json`, `ATTACKED-STRATEGIC_PLAN.md`, `ATTACKED-PANGENOME_SUITE.md`)
**Artifacts read:** `CRITIC-BRIEF.md`, `round.json`, `ATTACKED-STRATEGIC_PLAN.md`, `ATTACKED-PANGENOME_SUITE.md`, all 7 attack files (`attack-L1-correctness.md` through `attack-L7-reproduction.md`), all 7 adjudication files (`adjudication-L1.md` through `adjudication-L7.md`), `src/lib.rs`, `src/formats.rs`, `src/seeding/mod.rs`, `src/seeding/translated.rs`, `docs/BENCHMARK_RESULTS.md`, and `Cargo.toml`.

---

## 1. Execution Status: UNRUN vs. UNCHALLENGED (and Brief Discrepancy)

* **Zero lenses were UNRUN:** All 7 declared lenses (`correctness`, `sampling`, `statistics`, `provenance`, `inference`, `prior-art`, `reproduction`) filed both an attack (`attack-L*.md`) and an adjudication (`adjudication-L*.md`).
* **Metadata discrepancy in `CRITIC-BRIEF.md:51`:** `CRITIC-BRIEF.md:51` states `attack-L5-inference.md (inference) — NO ADJUDICATION FILED`, but `reports/rounds/round-1/adjudication-L5.md` (102 lines) exists on disk in the round directory. `L5` was **not** unrun or unadjudicated; `CRITIC-BRIEF.md` reflects stale state prior to `adjudication-L5.md` being registered.
* **Zero claims were UNCHALLENGED across the round:** Every claim `C1`–`C11` was attacked by at least 5 of the 7 lenses.

---

## 2. Per-Lens Coverage and Thin Lenses

| Lens | Attack File | Findings Filed | Claims Challenged | Claims UNCHALLENGED by This Lens | Adjudication Rulings |
| :--- | :--- | :---: | :--- | :--- | :--- |
| L1 (`correctness`) | `attack-L1-correctness.md` | 11 | C1–C11 | None | 2 UPHELD, 9 NARROWED |
| L2 (`sampling`) | `attack-L2-sampling.md` | 5 | C2, C3, C4, C8 (+ framing note F5) | C1, C5, C6, C7, C9, C10, C11 | 4 UPHELD, 1 NARROWED |
| L3 (`statistics`) | `attack-L3-statistics.md` | 7 | C2, C3, C4, C7, C9, C10 (+ scope note F7) | C1, C5, C6, C8, C11 | 6 UPHELD, 1 NARROWED |
| L4 (`provenance`) | `attack-L4-provenance.md` | 11 | C1–C11 | None | 4 UPHELD, 7 NARROWED |
| L5 (`inference`) | `attack-L5-inference.md` | 11 | C1–C11 | None | 10 UPHELD, 1 NARROWED |
| L6 (`prior-art`) | `attack-L6-prior-art.md` | 11 | C1–C11 | None | 10 UPHELD, 1 NARROWED |
| L7 (`reproduction`) | `attack-L7-reproduction.md` | 11 | C1–C11 | None | 10 UPHELD, 1 NARROWED |

* **Thin lenses:**
  * **L2 (`sampling`)** is the thinnest lens (5 findings; 7 of 11 claims unchallenged). While C1, C5, C6, C7, C9, C10, and C11 largely describe unimplemented architecture/CLI contracts, L2 did not examine whether existing pairwise/MAF benchmarks sample `W`/`P` coordinate-liftover edge cases (`C11`) or whether `horizon_wga` abstention behavior (`PANGENOME_SUITE.md:17`) holds across chromosome-scale inputs (`cactus241way` `chr22 1MB` was only flagged as a framing omission in F5, not audited for sampling/abstention).
  * **L3 (`statistics`)** is the second thinnest (7 findings; 5 claims unchallenged: C1, C5, C6, C8, C11). C1, C5, and C6 are CLI/format claims outside statistics, but C8 (HPRC graph usability/tangling) and C11 (lossy vs. exact liftover) have quantitative failure modes that L3 left to other lenses.

---

## 3. Claims Nobody Could Falsify (Unfalsifiable Framing)

* **Future-tense commitments and normative rules (`C5`, `C6`, `C7`, `C8`, and roadmap `C9`–`C11`):**
  * `round.json` extracts claims from two roadmap/strategy documents (`docs/STRATEGIC_PLAN.md` and `docs/PANGENOME_SUITE.md`) phrased in future tense or normative imperatives (*"We will never invent..."* in `C5`, *"Every tool must..."* in `C6`, *"We will stream..."* in `C7`, *"We will focus..."* in `C8`, *"strategic roadmap"* in `PANGENOME_SUITE.md:34`).
  * No codebase or benchmark evidence can falsify a future intention. As a result, `adjudication-L1.md`, `adjudication-L2.md`, and `adjudication-L4.md` repeatedly had to **NARROW** `REFUTED` findings on these claims—not because the technical criticisms were wrong, but because the claims were unfalsifiable promises rather than testable present-state or architectural-feasibility claims.

---

## 4. Single-Party Framing Bias (`runner_is_author=true`, `claim_list_author="runner"`)

`round.json:149-151` records `runner_is_author: true`, `claim_list_author: "runner"`, and `human_claim: ""`. Beyond the `PANGENOME_SUITE.md:12-17` and `BENCHMARK_RESULTS.md:37` omissions already flagged by `L2` (`F5`) and `L3` (`F7`), the runner's `C1`–`C11` claim list omitted five concrete technical claims present in the attacked documents themselves:

1. **`.fai`/`.gzi`-style sidecar indexing for pangenome graphs (`ATTACKED-STRATEGIC_PLAN.md:17`):**
   * `round.json:46` truncated the second sentence of `C5`: *"Any necessary indexing must be a lightweight, optional sidecar file (like `.fai` or `.gzi`), never a locked central database."*
   * L5 (`F5`) noted the contradiction between pure `stdin`/`stdout` and optional sidecars, but **no lens attacked the technical claim that `.fai`/`.gzi`-style sidecars suffice for graph operations**. `.fai` indexes linear sequence byte offsets and `.gzi` indexes BGZF blocks; neither indexes 2D `L`-line graph topology, spaced seeds across bifurcations, or `W`/`P` path intervals without either an actual graph index or loading the graph into RAM.
2. **Spaced-seed (`MASK_RY20`) and Amino Acid (`AA5`) indexing on GFA `S` lines (`ATTACKED-PANGENOME_SUITE.md:15`):**
   * `C10` (`round.json:91-94`) quoted only line 16 (banded DP across bifurcations), omitting line 15: *"Ingest the sequence (`S`) lines of a GFA and build a graph-aware spaced-seed index, utilizing transition-tolerant (`MASK_RY20`) and Amino Acid masks."*
   * **No lens examined `MASK_RY20` or AA masks in `horizon_wga --graph`.** In `src/seeding/mod.rs:13`, `MASK_RY20` is a 25-bp spaced seed (`"1110110110111011011011011"`), and `find_translated_aa5_hits` (`src/seeding/translated.rs:72`) requires 5-codon (15-bp) windows. In pangenome GFAs (e.g., HPRC Minigraph-Cactus/PGGB), `S` segments at SNPs and small indels are frequently 1–20 bp; indexing `S` lines individually without traversing `L` junctions drops all seeds spanning variant nodes.
3. **"Horizon Bubble" emission in `horizon_clean` (`ATTACKED-PANGENOME_SUITE.md:9`):**
   * `C9` (`round.json:82-85`) quoted only line 8, omitting line 9: *"If a path through the graph does not yield sufficient bit-savings, explicitly sever the edge and emit a Horizon Bubble."*
   * **No lens examined what a "Horizon Bubble" is in `src/formats.rs:15, 50`:** it is a placeholder `S` segment with `*` sequence and `HB:i:1` tag representing an unaligned pairwise gap, not a GFA ultrabubble/snarl. Severing an `L` line in an external pangenome GFA and emitting a `*` `HB:i:1` node breaks existing `W`/`P` haplotype walks that traverse the severed link.
4. **Graceful abstention at low-confidence divergence thresholds (`ATTACKED-PANGENOME_SUITE.md:17`):**
   * `L2` (`F5`) flagged the omission of `PANGENOME_SUITE.md:12-17` from `C1`–`C11`, and `L3` (`F1`) cited `cactus241way` `chr22 1MB` (`9.9%` recall, `494` Horizon Bubbles in `docs/BENCHMARK_RESULTS.md:37`). However, **no lens audited whether `horizon_wga`'s abstention mechanism distinguishes true divergence abstention from scale-dependent recall collapse** (recall drops 4.5x from `45.3%` on 38 kb `chr22_random` to `9.9%` on `chr22 1MB`).
5. **Strict structural-path annotation liftover (`ATTACKED-PANGENOME_SUITE.md:30`):**
   * `C11` (`round.json:100-103`) quoted line 29 (`W`/`P` traversal with `coord::Interval`), omitting line 30: *"Project functional annotations from a reference path (e.g., `GRCh38`) to any alternate haplotype strictly via high-confidence structural paths."*
   * **No lens attacked the biological failure mode of strict structural-path projection:** when an alternate haplotype contains novel insertions, gene conversions, or structural rearrangements not shared as a contiguous structural path with `GRCh38`, strict path projection drops or fragments exons/genes (unlike alignment-based liftover such as `Liftoff`).

---

## 5. Contradictions Between Adjudications

1. **`pub mod gfa;` declaration in `src/lib.rs` (`adjudication-L5.md` vs. `adjudication-L1.md` / `adjudication-L7.md`):**
   * `adjudication-L5.md:34` asserts: *"Moreover, `gfa` is not even declared as a module in `src/lib.rs:3-9`"* (repeated at line 62: *"Neither `coord` nor `gfa` is declared in `src/lib.rs:3-9`"*).
   * In fact, `src/lib.rs:394` **does** declare `pub mod gfa;` (at the end of the file, after `mod tests`), as correctly verified by `adjudication-L1.md:22` (`src/lib.rs:394`) and `adjudication-L7.md:10` (`src/lib.rs:394:pub mod gfa;`).
2. **Verdict standard on roadmap/future-tense claims (`L1`/`L2`/`L4` vs. `L5`/`L6`/`L7`):**
   * `adjudication-L1.md` (F1, F5–F11), `adjudication-L2.md` (F1 on C8), and `adjudication-L4.md` (F1, F5–F7, F9–F11) **NARROWED** `REFUTED` findings because `STRATEGIC_PLAN.md` uses future-tense/normative phrasing and `PANGENOME_SUITE.md:34` labels itself a strategic roadmap.
   * `adjudication-L5.md` (F1, F6–F11), `adjudication-L6.md` (F5–F11), and `adjudication-L7.md` (F1, F5–F11) **UPHELD** `REFUTED` / `NARROWED` verdicts on those same claims without applying the roadmap/future-tense narrowing rule.
3. **Standard vs. non-standard GFA tags/comments (`L1`/`L5` vs. `L6` on `C5`):**
   * `adjudication-L1.md:26` and `adjudication-L5.md:41` rule that `# EG:Z:` comments, `*` sequences with `LN:i:`, and `HB:i:1` / `XS:i:` / `YS:i:` / `ST:Z:` tags emitted by `src/formats.rs:15, 43, 50, 102` are **valid standard GFA v1.0/v1.1 syntax**.
   * `adjudication-L6.md:33-34` upholds the characterization that Horizon *"emits non-standard GFA tags/comments"*.

---

## 6. Open Gaps Summary

1. **Gap: `.fai`/`.gzi`-style sidecar indexing for graphs (`ATTACKED-STRATEGIC_PLAN.md:17`).**
   * *Why it matters:* Treating linear-offset sidecars (`.fai`, `.gzi`) as sufficient for graph indexing hides the need for an actual graph topology/k-mer/walk index, which reintroduces the RAM or disk-index problem `C4`/`C5` claim to avoid.
   * *What would close it:* Specify an explicit graph-index format (e.g., GBZ/gfaidx-style segment + adjacency + walk index) or drop the `.fai`/`.gzi` analogy in `docs/STRATEGIC_PLAN.md:17`.

2. **Gap: `MASK_RY20` (25-bp) and `AA5` (15-bp) seeding on fragmented GFA `S` lines (`ATTACKED-PANGENOME_SUITE.md:15`, `src/seeding/mod.rs:13`, `src/seeding/translated.rs:72`).**
   * *Why it matters:* In pangenome GFAs, variant `S` nodes are often shorter than 25 bp (or 15 bp), so `S`-line-only spaced-seed indexing silently drops seeds that cross `L`-line bifurcations.
   * *What would close it:* Specify path/k-mer traversal across `L` edges during graph seed indexing in `docs/PANGENOME_SUITE.md:15` (and mark as unimplemented).

3. **Gap: `horizon_clean` severing `L` lines and emitting `*`-sequence Horizon Bubbles (`ATTACKED-PANGENOME_SUITE.md:9`, `src/formats.rs:15, 50`).**
   * *Why it matters:* Severing an `L` edge and inserting a `*`-sequence `HB:i:1` segment breaks existing `W`/`P` haplotype walks in external GFAs unless walk/path repair is defined.
   * *What would close it:* Define how `W`/`P` lines are updated or split when `L` edges are severed, and clarify that `HB:i:1` in `src/formats.rs` is currently a pairwise unaligned-gap placeholder.

4. **Gap: Scale-dependent recall collapse vs. "graceful abstention" (`ATTACKED-PANGENOME_SUITE.md:17`, `docs/BENCHMARK_RESULTS.md:37`).**
   * *Why it matters:* `cactus241way` recall drops from `45.3%` (`chr22_random`, 38 kb) to `9.9%` (`chr22 1MB`) with `494` Horizon Bubbles; without auditing this collapse, false-negative failure at 1 MB scale is misframed as "graceful abstention."
   * *What would close it:* Audit `horizon_wga` seeding/chaining behavior on `cactus241way` `chr22 1MB` and document the recall collapse in `docs/BENCHMARK_RESULTS.md` and `docs/PANGENOME_SUITE.md`.

5. **Gap: Strict structural-path liftover failure on divergent/rearranged haplotypes (`ATTACKED-PANGENOME_SUITE.md:30`).**
   * *Why it matters:* Projecting annotations strictly along shared graph structural paths drops exons/genes in haplotype-specific insertions, gene conversions, or rearrangements.
   * *What would close it:* Document fallback alignment/splice handling (or scope limitation) for non-syntenic haplotype regions in `docs/PANGENOME_SUITE.md:30`.

6. **Gap: Adjudication factual error and cross-lens verdict inconsistency (`adjudication-L5.md:34, 62`, `CRITIC-BRIEF.md:51`).**
   * *Why it matters:* `adjudication-L5.md` falsely states `gfa` is not declared in `src/lib.rs` (it is at `src/lib.rs:394`), `CRITIC-BRIEF.md:51` misreports `L5` as having no adjudication, and `L1`/`L2`/`L4` vs. `L5`/`L6`/`L7` apply inconsistent standards to future-tense/roadmap claims.
   * *What would close it:* Correct `adjudication-L5.md:34, 62` (`pub mod gfa;` is at `src/lib.rs:394`), refresh `round.json`/brief status for `adjudication-L5.md`, and standardize how roadmap-tense vs. present-tense claims are scored when closing Round 1.
