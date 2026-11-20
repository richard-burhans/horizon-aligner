# Round 5 Adjudication — Correctness Lens

**Attack digest verified:** `3d478f9488c1c5b7cf95076fce23d58f888a8551681a04e3119a737297f29264` (`reports/rounds/round-5/attack-correctness.md` matches `ADJ-BRIEF-correctness.md`).
**Method cards:** None exist in this repository (`refs/methods/` is absent).
**Prior-round history:** Checked via `scripts/round history` for `horizon_gfa`, `EG:Z`, `BandLayout`, and `horizon_guide`. Rounds 3 and 4 already upheld that `horizon_gfa stats` is a counter-only record-discarding pass; Round 5's `C1` introduces a new head-to-head comparison against `vg v1.60.0` (`vg stats -z -l -N` and `vg convert -g -p`) without carrying that qualification into `manuscript/horizon_paper.md:10,82`.

---

## F1 — `horizon_gfa stats` is a counter-only line stream and is not functionally equivalent to `vg stats -z -l -N` or `vg convert -g -p` (`C1`)

- **Diagnosis:** **UPHELD.** Verified in `src/bin/horizon_gfa.rs:20-56`: `cmd_stats` calls `stream_gfa(...)` and only increments five `u64` counters (`s_count`, `s_bases`, `l_count`, `p_count`, `w_count`), discarding every `S`, `L`, `P`, and `W` record without constructing an in-memory graph, validating `L`/`P`/`W` references, or serializing a packed graph. While `manuscript/horizon_paper.md:47,74` already label `horizon_gfa stats` as "counter-only", the new Round 5 head-to-head claim against `vg stats -z -l -N` and `vg convert -g -p` in `manuscript/horizon_paper.md:10` (Abstract) and `:82` (Results) omits that disclosure.
- **Repair:** **UPHELD — the fix is free.** Narrowing `manuscript/horizon_paper.md:10,82` to state that `horizon_gfa stats` is a counter-only streaming pass that does not construct, validate, or convert the graph requires only a small text edit and no re-measurement.
- **What survives:** The `0.08 s` / `22.6 MB` RSS timing on the `10,000`-segment `chr22` subgraph stands as a counter-only streaming measurement, not an equivalent-work comparison to `vg stats` or `vg convert`.

---

## F2 — `horizon_guide --progressive --hal --rbest` emits a flat star-tree HAL with at most one bounding-box segment per species pair, discarding `1,105+` MAF blocks, all indels, and all progressive/ancestral rescue blocks (`C3`)

- **Diagnosis:** **UPHELD.** Verified in `src/bin/horizon_guide/main.rs:94-103`, `:271-296`, and `src/hal/mod.rs:90-151`:
  1. `align_species_pair` constructs a `HalBlock` only when `xi == 0 && yi == 0`, collapsing the first contig pair into a single ungapped bounding-box tuple `(x0, y0.min(y1), x1 - x0 + 1, y1 < y0)` from `aln.aligned_pairs.first()` to `.last()`, omitting all other MAF blocks, internal indels, secondary contigs (`xi > 0 || yi > 0`), and all blocks added by `run_progressive_sister_rescue` / `run_ancestral_consensus_rescue` (which mutate only `pairwise_mafs`).
  2. `HalAlignment::from_star_blocks` builds a flat star tree rooted at `hinge_sp` using only `contigs.first()` (`format!("{ref_name}.chr1")`), explaining why `1,109` MAF blocks (`967` multi-species) yield only `4` HAL bottom segments (`1` per query species).
- **Repair:** **UPHELD — the fix is free.** Narrowing `manuscript/horizon_paper.md:82` (and `docs/BENCHMARK_RESULTS.md:11,58`) to disclose that `--hal` currently writes a `1`-contig star-tree summary (`4` bounding-box bottom segments, one per query species) rather than a progressive HAL encoding of the `1,109` MAF blocks requires only text edits and no re-measurement.
- **What survives:** `horizon_merge` (`92` `chrM` species, `21.78 s`, `460` multi-species blocks) and the MAF output (`1,109` MAF blocks) of `horizon_guide` survive; the `.hal` output is a `4`-segment star-tree bounding-box summary rather than a progressive HAL encoding of those blocks.

---

## F3 — `BandLayout` in `src/mdl/dp.rs` centers the band on `c(i) = floor(i · m / max(n, 1))`, not `i`, and clamps `W` to `min(band_w, max(n, m))` (`C5`)

- **Diagnosis:** **UPHELD.** Verified in `src/mdl/dp.rs:141-172, 244-246`: `BandLayout` (the forward/backward `α`/`β` table layout) clamps `band_w` to `band_w.min(n.max(m))` and computes `slot(i, j)` using `diff = (j as isize) - (self.center(i) as isize)` where `center(i) = band_center(i, m, n) = floor(i * m / max(n, 1))`. Thus `manuscript/horizon_companion.md:139-145`'s `|j − i| ≤ W` and `idx(i, j) = i · (2W + 1) + (j − i + W)` hold only for square gaps (`N = M`) or X-drop flank extension (`src/mdl/xdrop.rs:166-177`).
- **Repair:** **UPHELD — the fix is free.** Updating `manuscript/horizon_companion.md:139-145` (and Figure C3A's caption if needed) to state `c(i) = ⌊i · M / max(N, 1)⌋`, `|j − c(i)| ≤ W` with `W = min(band_w, max(N, M))`, and `idx(i, j) = i · (2W + 1) + (j − c(i) + W)` for `BandLayout`, while noting `c(i) = i` for X-drop extension and square gaps, is a text-only edit requiring no re-measurement.
- **What survives:** The contiguous flat `Vec<f64>` band layout of stride `2W + 1` survives with the row-scaled center `c(i)` and clamped `W` in `BandLayout`.

---

## F4 — `# EG:Z:` lines are skipped as comments by `src/gfa/mod.rs`, `collect_dups` omits non-`Dup` nested inversions, and `horizon_resolve` emits `EG:Z:` on `W` lines without `#` (`C6`)

- **Diagnosis:** **UPHELD.** Verified across all cited files:
  1. `src/formats.rs:153-158` (`format_gfa_paralog_tags`) and `src/merge.rs:361` emit `# EG:Z:<ref>__<qry>\t<bins>` comment lines, and `src/gfa/mod.rs:317-323` explicitly skips every line starting with `#` (`trimmed.starts_with('#')`), so downstream GFA parsers (including Horizon's own) cannot read `# EG:Z:` lines as GFA tags.
  2. `src/egraph.rs:518-542` (`collect_dups`) only collects block labels wrapped in `Expr::Dup(...)`, omitting nested inversions resolved without a `Dup` wrapper, and `src/bin/horizon_wga.rs:387-410` still writes uncollapsed `aln` to MAF and GFA `S`/`L` records.
  3. `src/bin/horizon_resolve.rs:473-485` writes `EG:Z:<tags>` as a standard optional field on `W` lines (when `--tag-walks` is passed), without a `# ` prefix.
- **Repair:** **UPHELD — the fix is free.** Replacing `manuscript/horizon_companion.md:52` to accurately distinguish `# EG:Z:` comment lines (`horizon_wga --resolve-dups`, `horizon_merge`, `Expr::Dup`-collapsed 1-kb bin labels only, uncollapsed `S`/`L`/MAF) from `W`-line `EG:Z:` optional tags (`horizon_resolve --tag-walks`) is a text-only edit requiring no re-measurement.
- **What survives:** E-graph paralog/walk annotation exists in two distinct forms (`# EG:Z:` comment lines for `Expr::Dup` bins in `horizon_wga`/`horizon_merge`, and `W`-line `EG:Z:` tags in `horizon_resolve --tag-walks`), rather than `# EG:Z:` GFA tags recording all collapsed duplications and nested inversions for downstream tools.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `src/bin/horizon_gfa.rs:20-56` is counter-only and record-discarding, not equivalent to `vg stats`/`vg convert`; narrowing `horizon_paper.md:10,82` is free |
| F2 | UPHELD | `src/bin/horizon_guide/main.rs:94-103,291-296` and `src/hal/mod.rs:90-151` emit a 1-contig star-tree HAL with at most 1 bounding-box block per query pair (`4` segments total); narrowing `horizon_paper.md:82` and `BENCHMARK_RESULTS.md:11,58` is free |
| F3 | UPHELD | `src/mdl/dp.rs:141-172,244-246` centers `BandLayout` on `floor(i*m/max(n,1))` and clamps `W = min(band_w, max(n,m))`, unlike `horizon_companion.md:139-145`'s `j - i`; updating the formula is free |
| F4 | UPHELD | `src/gfa/mod.rs:317-323` skips `# EG:Z:` as comments, `collect_dups` omits non-`Dup` inversions, and `horizon_resolve` emits `W`-line `EG:Z:` tags; updating `horizon_companion.md:52` is free |
