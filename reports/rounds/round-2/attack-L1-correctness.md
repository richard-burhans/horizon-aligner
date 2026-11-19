# Round 2 Attack Report -- Correctness Lens

**Method cards & literature check:** No method cards exist in `horizon-aligner` (checked via workspace search for `*card*`, `*lit*`, `*bib*`). Prior adjudications were checked via `scripts/round history "gfa" --kind adjudication` and `reports/adversarial_log.md` (Round 1 settled `src/gfa/mod.rs` zero-copy and `L`-line panic, and left open `R1-IMPL`, `R1-GFA`, `R1-HMM`).

---

## F1 -- C3: `cargo test` has 40 passing + 1 ignored tests (not 39), and 4 of 9 binaries have 0 unit tests

- **Claim:** C3 (`docs/STRATEGIC_PLAN.md:3`)
- **Verdict:** **REFUTED**
- **Concrete flaw & evidence:**
  Running `cargo test --no-fail-fast` and `cargo test -- --list` in `/home/rcb112/project/gemini/bench/horizon-aligner` produces:
  - `src/lib.rs`: `35 passed; 0 failed; 1 ignored` (36 tests total)
  - `src/bin/horizon_clean.rs`: `1 passed`
  - `src/bin/horizon_gfa.rs`: `1 passed`
  - `src/bin/horizon_lift.rs`: `1 passed`
  - `src/bin/horizon_resolve.rs`: `1 passed`
  - `src/bin/horizon_wga.rs`: `1 passed`
  - `src/bin/bench_chr1.rs`: `0 tests`
  - `src/bin/horizon_eval.rs`: `0 tests`
  - `src/bin/horizon_guide.rs`: `0 tests`
  - `src/bin/horizon_merge.rs`: `0 tests`
  Total passing tests is **40** (with **1 ignored**), not `39`, and four of the claimed binaries (`bench_chr1`, `horizon_eval`, `horizon_guide`, `horizon_merge`) have zero unit tests.
- **Correction:** Update `docs/STRATEGIC_PLAN.md:3` to report `40 passing unit tests (1 ignored)` and explicitly note that `bench_chr1`, `horizon_eval`, `horizon_guide`, and `horizon_merge` currently have `0` unit tests.

---

## F2 -- C2: `mdl_xdrop_extend` uses unnormalized hardcoded constants; `banded_mdl_gap_with_params` has a stateless Viterbi traceback bug; `align_seq_to_dag_mdl` double-counts boundary `Iy` transitions

- **Claim:** C2 (`docs/BENCHMARK_RESULTS.md:5`)
- **Verdict:** **REFUTED**
- **Concrete flaw & evidence:**
  1. **Unnormalized & non-adaptive X-drop extension:** `mdl_xdrop_extend` (`src/mdl/dp.rs:369-532`) does not accept `HmmParams` and does not apply `mm_bits`/`im_bits` transition normalization. It hardcodes `IDENTITY_PRIOR` (`src/mdl/dp.rs:445`), `GAP_EXTEND_BITS` (`458, 473`), `GAP_OPEN_BITS` (`460, 475`), and adds `sb` directly (`v = dp[ps] + sb` at `line 447`).
  2. **Stateless Viterbi traceback in `banded_mdl_gap_with_params`:** `src/mdl/dp.rs:270-293` stores no backpointer matrix and instead compares `m_dp.get(i, j)`, `ix_dp.get(i, j)`, `iy_dp.get(i, j)` without transition weights (`mm_bits`, `im_bits`, `go_bits`, `ge_bits`). After stepping from `M(i, j)` to `(i-1, j-1)`, it re-picks `max(m_dp, ix_dp, iy_dp)` at `(i-1, j-1)` independently of which predecessor state actually maximized `m_dp(i, j)`.
  3. **Boundary `Iy` double-counting & default params in DAG DP:** In `src/mdl/dag.rs:165-172` and `258-315`, `align_seq_to_dag_mdl` permits `Iy` transitions at both `i_u = nu` (end of predecessor `u`) and `i_s = 0` (start of successor `s`), double-counting identical `(None, Some(j))` gap columns across node boundaries in Forward/Backward `Z`. Furthermore, `horizon_wga --graph` (`src/bin/horizon_wga.rs:195`) passes `HmmParams::default()`, not `estimate_hmm_params`.
  4. **Shrinkage formula:** `estimate_hmm_params` (`src/mdl/mod.rs:828-836`) estimates `gap_extend_bits` via geometric mean-length shrinkage `(total_indel_len + 2.5) / (indel_events + 1.0)`, not Beta prior shrinkage.
- **Correction:** Narrow `docs/BENCHMARK_RESULTS.md:5` to state that `HmmParams` (with Beta shrinkage for `identity_prior`, empirical/Laplace shrinkage for `gap_open_bits`/`gap_extend_bits`) and normalized transitions apply to `banded_mdl_gap_with_params`, while `mdl_xdrop_extend` and `horizon_wga --graph` use fixed `HmmParams::default()` / `IDENTITY_PRIOR` constants, and fix the Viterbi backpointer and DAG `Iy` boundary bugs.

---

## F3 -- C4: `horizon_lift --bridge-bp 10` absorbs unmapped gaps into `EXACT` counts, and `stdin` (`"-"`) fails on Pass 2

- **Claim:** C4 (`docs/BENCHMARK_RESULTS.md:21`)
- **Verdict:** **NARROWED**
- **Concrete flaw & evidence:**
  1. **`EXACT` inflation under `--bridge-bp 10`:** In `src/bin/horizon_lift.rs:256-270`, when `piece.src_start - prev.src_end <= bridge_bp`, coalescing sets `prev.src_end = piece.src_end` (`line 267`). At `line 490`, `mapped_bp` sums `p.src_end - p.src_start`, which now includes any `<= 10 bp` unmapped source gap between the bridged segments. Consequently, `mapped_bp == total_bp` (`line 508`) classifies bridged/gapped windows as `EXACT` (`469 exact` includes bridged intervals).
  2. **Broken `stdin` (`"-"`) support:** `run_lift` (`src/bin/horizon_lift.rs:316, 419`) calls `stream_gfa(cfg.gfa_path, ...)` twice. When `gfa_path == "-"` (advertised at `line 608` and `docs/STRATEGIC_PLAN.md:3`), Pass 1 exhausts `stdin` and Pass 2 reads `0` segments, leaving `seg_lens` empty.
- **Correction:** Track true `mapped_bp` separately from bridged span in `LiftedPiece` (so bridged intervals with missing bases are reported as `PARTIAL` or `BRIDGED`, not `EXACT`), and clarify in `docs/BENCHMARK_RESULTS.md:21` that `469 exact` includes `<= 10 bp` bridged windows.

---

## F4 -- C5: `5,894,191` aligned pairs on a 5 Mbp locus reflects many-to-many LCB multiplicity, and `horizon_clean` never severs supported links or works on `stdin`

- **Claim:** C5 (`docs/BENCHMARK_RESULTS.md:23`)
- **Verdict:** **NARROWED**
- **Concrete flaw & evidence:**
  1. **Pair count > locus length (`5,894,191 > 5,000,000`):** Both `GRCh38.chr22:15.5M-20.5M` and `CHM13.chr22:15.73M-20.73M` are `5,000,000 bp`. Because `horizon_wga --resolve-dups` was run without `--rbest`, `src/mdl/mod.rs:745-746` only deduplicates exact `(x, y)` pairs across `246 fwd / 184 rev` LCBs, and `--resolve-dups` (`src/bin/horizon_wga.rs:308-315`) only appends `# EG:Z:` comments without filtering paralogous pairs. Thus `5,894,191` (`117.9%` of 5 Mbp) is a many-to-many pair count, not 1-to-1 aligned bases.
  2. **`horizon_clean` behaviour:** `src/bin/horizon_clean.rs:228-237` only severs links with `sup == 0`, so walk-induced links in `/tmp/chr22_peri_5m.gfa` (`6,030 kept, 0 severed`) can never be severed. Like `horizon_lift`, `clean_gfa` (`src/bin/horizon_clean.rs:152, 246`) calls `stream_gfa` twice and fails on `"-"` `stdin`.
- **Correction:** Label `5,894,191` as many-to-many `(x, y)` aligned pairs (not 1-to-1 aligned bases) in `docs/BENCHMARK_RESULTS.md:23`.

---

## F5 -- C6: `horizon_resolve --collapse-dups` corrupts GFA v1.1 `W`-line `start`/`end` coordinates and discards stdout when `--out` is omitted

- **Claim:** C6 (`docs/BENCHMARK_RESULTS.md:64`)
- **Verdict:** **NARROWED**
- **Concrete flaw & evidence:**
  1. **Stale `W`-line coordinates:** In `src/bin/horizon_resolve.rs:424-441`, when `--collapse-dups` shortens a walk (e.g., `>a>b>b>c` -> `>a>b>c`), `pass2_rewrite` writes the original `cols[4]` (`start`) and `cols[5]` (`end`) verbatim (e.g., `0` and `80` in the unit test at `line 507` even though the collapsed walk is `60 bp`), breaking the GFA v1.1 invariant `end - start == sum(segment lengths)`. It also copies `L` lines verbatim (`line 446`) without updating topology.
  2. **No `stdout` or `stdin` streaming:** Without `--out`, `main` (`src/bin/horizon_resolve.rs:465`) writes to `io::sink()` instead of `stdout`, and `pass1_index` + `pass2_rewrite` (`lines 301, 367`) read the file twice (and `parse_args` at `line 103` rejects `"-"`).
- **Correction:** Update `W`-line `end` (or emit `end = start + sum(collapsed_seg_lens)`) in `horizon_resolve.rs:424-441` and note in `docs/BENCHMARK_RESULTS.md:64` that `horizon_resolve` rewrites `P`/`W` walks without recomputing `W` coordinates or `L` edges.

---

## F6 -- C7: `align_seq_to_dag_seeded` reachability closure expands to full ancestor-descendant components rather than s-t bubble subgraphs

- **Claim:** C7 (`docs/PANGENOME_SUITE.md:18`)
- **Verdict:** **NARROWED**
- **Concrete flaw & evidence:**
  In `src/mdl/dag.rs:503-518`, `align_seq_to_dag_seeded` runs a backward pass marking all ancestors of `seed_hits` active, followed by a forward pass marking all descendants of `active`. Any seed hit whose ancestors reach a shared root therefore activates **every** downstream branch from that root (not just nodes on paths through `seed_hits`), and for `nodes.len() <= 8` or `seed_hits.is_empty()` (`lines 466, 485`) it falls back to unseeded `align_seq_to_dag_mdl` (which also has the `Iy` boundary double-counting bug from F2).
- **Correction:** Clarify in `docs/PANGENOME_SUITE.md:18` that the current closure is an ancestor-then-descendant expansion (not an s-t / forward-backward intersection) with fallback to full `align_seq_to_dag_mdl` when `nodes.len() <= 8` or no seeds hit.

---

## F7 -- C8: `~0.29 GB/Mbp` does not match `26.5 GB` on `170.8 + 110.7 Mbp` (`0.094 GB/Mbp`)

- **Claim:** C8 (`docs/BENCHMARK_RESULTS.md:26`)
- **Verdict:** **NARROWED**
- **Concrete flaw & evidence:**
  `src/mdl/mod.rs:745-747` confirms `pairs_vec.sort_unstable(); pairs_vec.dedup();`, and `src/seeding/mod.rs:329-333` confirms `ContigIndex` stores `HashMap<u64, Vec<usize>>` seed maps. However, `26.5 GB / (170.8 + 110.7 Mbp) = 0.094 GB/Mbp` (or `0.155 GB/Mbp` against `170.8 Mbp` alone). The `~0.29 GB/Mbp` rate instead comes from the 20 Mbp `chr6` q-arm run (`5,863.3 MB / 20 Mbp = 0.293 GB/Mbp`) on the same line, not the cancelled full-`chr6` run.
- **Correction:** Attribute `~0.29 GB/Mbp` to the completed 20 Mbp `chr6` q-arm run (`5,863.3 MB / 20 Mbp`) rather than the cancelled `170.8 vs 110.7 Mbp` run.

---

## F8 -- C1: `stream_gfa` ignores GFA v1.0 `C` and GFA v1.1 `J` records and silently coerces invalid `W` coordinates to `0`

- **Claim:** C1 (`docs/STRATEGIC_PLAN.md:13`)
- **Verdict:** **NARROWED**
- **Concrete flaw & evidence:**
  `src/gfa/mod.rs:98-128` implements `stream_gfa` supporting `.gfa`, multi-member `.gfa.gz` via `MultiGzDecoder`, and `"-"` `stdin`, with unit tests for empty `L`-line orientations (`lines 219-276`) and multi-member `.gfa.gz` (`lines 279-318`). However, `stream_gfa_reader` (`src/gfa/mod.rs:203`) silently drops GFA v1.0 `C` (Containment) and GFA v1.1 `J` (Jump) lines, `lines 198-199` silently coerce `*` or unparseable `W` `start`/`end` fields to `0`, and `parse_path_steps` (`line 58`) can panic on multi-byte UTF-8 segment names (`split_at(tok.len() - 1)`).
- **Correction:** State in `docs/STRATEGIC_PLAN.md:13` that `stream_gfa` parses GFA v1.0/v1.1 `S`/`L`/`P`/`W` records (ignoring `C`/`J`/`U`/`O`).

---

## Verdict Table

| Claim | Verdict |
| :--- | :--- |
| **C1** | NARROWED |
| **C2** | REFUTED |
| **C3** | REFUTED |
| **C4** | NARROWED |
| **C5** | NARROWED |
| **C6** | NARROWED |
| **C7** | NARROWED |
| **C8** | NARROWED |
