# Round 5 Attack — Correctness Lens

**Scope & Disclosures (Rules 2–4):**
- **Literature (Rule 2):** No external papers from `../library` were needed or asserted; all findings below are verified directly against `src/`, `docs/BENCHMARK_RESULTS.md`, and `reports/rounds/round-5/ATTACKED-horizon_*.md`.
- **Method cards (Rule 3):** Per the project map, this project keeps **no method cards**.
- **Prior-round history (Rule 4):** Checked via `scripts/round history "horizon_gfa" --kind adjudication`, `scripts/round history "chr6" --kind adjudication`, and `scripts/round history "EG:Z" --kind adjudication`. Rounds 3 and 4 (`adjudication-L5.md:24 [F3]`, `adjudication-L7.md:32 [F3]`) already upheld that `horizon_gfa stats` is a counter-only record-discarding pass; what changed in Round 5 (`C1`) is the new head-to-head claim that `horizon_gfa stats` "outperforms" `vg stats -z -l -N` and `vg convert -g -p` on a `10,000`-segment `chr22` subgraph without disclosing that `horizon_gfa stats` performs no graph construction, validation, or format conversion.

---

## C1 — `horizon_gfa stats` vs `vg stats -z -l -N` and `vg convert -g -p` (`manuscript/horizon_paper.md:82`)

### Verdict: NARROWED

### F1 -- `horizon_gfa stats` is a counter-only line stream and is not functionally equivalent to `vg stats -z -l -N` or `vg convert -g -p`

- **Prior-round context:** Rounds 3 and 4 (`adjudication-L5.md:24 [F3]`) upheld that `horizon_gfa stats` (`12 GB` HPRC graph) is a counter-only streaming pass that discards records after counting. Round 5's `C1` introduces a new comparison against `vg v1.60.0` (`vg stats -z -l -N` and `vg convert -g -p`) on a `10,000`-segment `chr22` subgraph (`10,937,960` bp, `10,300` links, `1` path, `4` walks).
- **Evidence:** In `src/bin/horizon_gfa.rs:20-56` (`cmd_stats`), `horizon_gfa stats` calls `stream_gfa(...)` and only increments five `u64` counters (`s_count`, `s_bases`, `l_count`, `p_count`, `w_count`), immediately discarding every `S`, `L`, `P`, and `W` record. It never constructs an in-memory graph, never checks whether `L`-line endpoints exist in `S`, never validates `P`/`W` walk steps, and never serializes a packed/indexed graph (unlike `vg convert -g -p`, which builds and writes a PackedGraph, and `vg stats -z -l -N`, which loads and validates the variation graph).
- **Correction:** Narrow `manuscript/horizon_paper.md:82` (and the Abstract) to state that `horizon_gfa stats` (`0.08 s`, `22.6 MB` RSS) is a counter-only streaming pass that does not construct, validate, or convert the graph, rather than claiming it outperforms `vg stats -z -l -N` and `vg convert -g -p` on equivalent graph-loading/conversion work.

---

## C2 — Full `GRCh38.chr6` × `CHM13.chr6` concordance against HPRC Minigraph-Cactus and `minimap2` (`manuscript/horizon_paper.md:76`)

### Verdict: UNCHALLENGED

- **Checks performed:**
  - Verified exact coordinate recall arithmetic: `107,836,669 / 110,082,990 = 0.97959429...` (`97.96%`).
  - Verified why MAF pair count (`282,163,046`) is `30,518` pairs lower than `horizon_wga`'s raw `282,193,564` aligned pairs (`manuscript/horizon_paper.md:76`): `partition_colinear_blocks` (`src/formats.rs:161-243`) drops colinear runs shorter than `12` pairs (`if chain.len() >= 12`, lines `181, 228, 237`), which is accurately described as "after multi-chain colinear partitioning."
  - Note: The `12 GB` HPRC graph truth-comparison script and `chr6` MAF output were not re-executed here, so `C2` is marked `UNCHALLENGED`.

---

## C3 — `horizon_merge` (`92` `chrM` species) and `horizon_guide --progressive --hal --rbest` (`5` mammals) (`manuscript/horizon_paper.md:82`)

### Verdict: NARROWED

### F2 -- `horizon_guide --progressive --hal --rbest` emits a flat star-tree HAL with at most one bounding-box segment per species pair, discarding `1,105+` MAF blocks, all indels, and all progressive/ancestral rescue blocks

- **Evidence:**
  1. In `src/bin/horizon_guide/main.rs:94-103`, `align_species_pair` constructs `HalBlock` ONLY when `xi == 0 && yi == 0`, collapsing the entire contig alignment into a single ungapped bounding-box tuple `(x0, y0.min(y1), x1 - x0 + 1, y1 < y0)` from `aln.aligned_pairs.first()` to `aln.aligned_pairs.last()`. All other MAF blocks, internal gaps/indels, inversions, secondary contigs (`xi > 0 || yi > 0`), and all blocks added later by `run_progressive_sister_rescue` / `run_ancestral_consensus_rescue` (`src/bin/horizon_guide/main.rs:271-289`, which only mutate `pairwise_mafs`) are completely omitted from `hal_all_blocks`. That is why a `5`-mammal alignment with `1,109` MAF blocks (`967` multi-species) yields only **`4` HAL bottom segments** (`1` per query species).
  2. In `src/bin/horizon_guide/main.rs:291-296` and `src/hal/mod.rs:90-151`, `HalAlignment::from_star_blocks` builds a flat **star tree** rooted at `hinge_sp` (`hg38`) rather than a progressive guide-tree HAL, and hardcodes `format!("{ref_name}.chr1")` using only `contigs.first()`.
- **Correction:** Narrow `C3` in `manuscript/horizon_paper.md:82` (and `docs/BENCHMARK_RESULTS.md:11,58`) to disclose that `--hal` currently writes a `1`-contig star-tree summary (`4` bounding-box bottom segments, one per query species) rather than a progressive HAL encoding of the `1,109` MAF blocks, or fix `horizon_guide` to populate HAL segments from the progressive/rescued blocks before claiming HAL equivalence. (`horizon_merge` on `92` `chrM` species was not refuted.)

---

## C4 — `horizon_wga` pairwise accuracy across UCSC/B10K slices and `18%`-GC random DNA (`manuscript/horizon_paper.md:59-70`)

### Verdict: UNCHALLENGED

- **Checks performed:**
  - Inspected `src/bin/horizon_eval.rs:17-91, 154-239`: confirmed `horizon_eval` evaluates each `(ref_id, tgt_id)` pair present in `truth.maf` individually (`align_contig_pair_with_config` with `AlignConfig::default()` and `prepare_contig_index_opts(seq, false)`), rather than running all-vs-all scaffold alignment; this is already disclosed in `manuscript/horizon_paper.md:59`.
  - Inspected `src/bocpd.rs:87-126` and `src/mdl/mod.rs:654-667`: confirmed `bocpd_gc_profile_into` skips non-ACGT bases and `joint_pair_bits` / `GAP_OPEN_BITS` gating is wired in `src/mdl/mod.rs`.
  - Benchmark datasets were not re-run under this lens, so `C4` is marked `UNCHALLENGED`.

---

## C5 — `BandGrid` flat buffer index formula (`manuscript/horizon_companion.md:139-145`)

### Verdict: NARROWED

### F3 -- `BandLayout` in `src/mdl/dp.rs` centers the band on `c(i) = floor(i · m / max(n, 1))`, not `i`, and clamps `W` to `min(band_w, max(n, m))`

- **Evidence:** `manuscript/horizon_companion.md:139-145` states that for the forward (`α`) and backward (`β`) tables, for target index `i` and query index `j` within `|j − i| ≤ W`, the flat buffer index is `idx(i, j) = i · (2W + 1) + (j − i + W)`, mapping `j = i` to column `W`.
  - In `src/mdl/dp.rs:141-172` (`BandLayout`), `self.band_w` is clamped to `band_w.min(n.max(m))` (line `142`), and `slot(i, j)` (lines `160-172`) computes `diff = (j as isize) - (self.center(i) as isize)` where `center(i) = band_center(i, self.m, self.n) = floor(i * m / max(n, 1))` (`src/mdl/dp.rs:244-246`). Thus the actual condition is `|j − c(i)| ≤ W` and the actual index is `i · (2W + 1) + (j − c(i) + W)`.
  - The unscaled `j − i` formula in `C5` only holds in `BandLayout` when `n == m` (square gaps), or in `src/mdl/xdrop.rs:166-177` (X-drop flank extension, which uses `diff = j - i`).
- **Correction:** Update `manuscript/horizon_companion.md:139-145` (and Figure C3A's caption if needed) to state `c(i) = ⌊i · M / max(N, 1)⌋`, `|j − c(i)| ≤ W` with `W = min(band_w, max(N, M))`, and `idx(i, j) = i · (2W + 1) + (j − c(i) + W)` for `BandLayout` (`src/mdl/dp.rs`), noting that `c(i) = i` is used for X-drop extension (`src/mdl/xdrop.rs`) and square gaps (`N = M`).

---

## C6 — E-graph collapsed duplications and nested inversions recorded as `# EG:Z:` tags (`manuscript/horizon_companion.md:52`)

### Verdict: REFUTED

### F4 -- `# EG:Z:` lines are skipped as comments by `src/gfa/mod.rs`, `collect_dups` omits non-`Dup` nested inversions, and `horizon_resolve` emits `EG:Z:` on `W` lines without `#`

- **Evidence:**
  1. **Ignored by GFA parsers (including Horizon's own):** `src/formats.rs:153-158` (`format_gfa_paralog_tags`) and `src/merge.rs:361` emit `# EG:Z:<ref>__<qry>\t<bins>` comment lines. In `src/gfa/mod.rs:317-323`, Horizon's own GFA parser explicitly skips every line starting with `#` (`if trimmed.is_empty() || trimmed.starts_with('#') ... { continue; }`). Downstream GFA tools cannot read `# EG:Z:` lines as GFA tags.
  2. **Nested inversions without `Dup` are not recorded:** In `src/bin/horizon_wga.rs:401-410` and `src/formats.rs:101-149`, `resolve_alignment_duplications` runs `collect_dups(&expr)` (`src/egraph.rs:518-542`), which only collects block labels inside `Expr::Dup(...)`. Nested inversions resolved solely by `inv(inv(X)) -> X` or `inv(seq(A, B)) -> seq(inv(B), inv(A))` without a `Dup` wrapper produce no `# EG:Z:` output. Moreover, `horizon_wga --resolve-dups` still writes the uncollapsed `aln` to both MAF and GFA `S`/`L` records (`src/bin/horizon_wga.rs:387-400`).
  3. **`horizon_resolve` uses `W`-line `EG:Z:` tags, not `# EG:Z:`:** `src/bin/horizon_resolve.rs:473-485` records `EG:Z:<tags>` as an optional field on `W` lines (only when `--tag-walks` is passed), with no `# ` prefix.
- **Correction:** Replace `manuscript/horizon_companion.md:52` to state accurately that: (a) `horizon_wga --resolve-dups` and `horizon_merge` append `# EG:Z:` comment lines listing `Expr::Dup`-collapsed 1-kb bin labels (which GFA parsers skip as comments, while leaving `S`/`L`/MAF records uncollapsed, and omitting non-`Dup` inversions), whereas (b) `horizon_resolve --tag-walks` writes `EG:Z:` as a standard optional tag on `W` lines.

---

## Summary Table

| Claim | Verdict | Findings |
| :--- | :--- | :--- |
| C1 | NARROWED | F1 (`horizon_gfa stats` is counter-only record-discarding; not equivalent to `vg stats`/`vg convert`) |
| C2 | UNCHALLENGED | — |
| C3 | NARROWED | F2 (`horizon_guide --hal` emits 1-contig star-tree with 1 bounding-box block per pair, discarding 1,105+ MAF blocks) |
| C4 | UNCHALLENGED | — |
| C5 | NARROWED | F3 (`BandLayout` in `src/mdl/dp.rs` centers at `floor(i*m/max(n,1))` with `W = min(band_w, max(n,m))`, not `j - i`) |
| C6 | REFUTED | F4 (`# EG:Z:` skipped as `#` comment by `src/gfa/mod.rs`; `collect_dups` drops non-`Dup` inversions; `horizon_resolve` uses `W`-line `EG:Z:`) |
