# Round 4 Adjudication — Reproduction (L7)

- **Attack answered:** `reports/rounds/round-4/attack-L7-reproduction.md` (digest `8411a3f3e86a`, verified against `8411a3f3e86a45eb6cfe02a999368c1d6aaccc8f295745c9a2d3173aafddde88`)
- **Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `33bede8`)
- **Prior round context:** `scripts/round history "horizon_gfa" --kind adjudication` shows Round 3 (`adjudication-L7.md:38`) previously upheld `horizon_gfa stats` reproduction on `488b020`, while Round 3's completeness critic (`adjudication-completeness.md:49`) flagged that `horizon_gfa extract` was left unreproduced. Per `docs/round_rubric.md`, this project keeps no method cards in `refs/methods/`.

---

## F1 -- C1: 10 binaries in `src/bin/`, partial `stdin` composition, and 2 binaries outside the MDL/BOCPD core

- **Ruling:** **NARROWED** (**UPHELD** on the defect; **NARROWED** on the prescribed repair).
- **What survives:** All three factual premises hold on `HEAD` (`33bede8`):
  1. `src/bin/` contains **10** binary targets (`bench_chr1.rs`, `horizon.rs`, `horizon_clean.rs`, `horizon_eval.rs`, `horizon_gfa.rs`, `horizon_guide.rs`, `horizon_lift.rs`, `horizon_merge.rs`, `horizon_resolve.rs`, `horizon_wga.rs`), where `horizon` (`src/bin/horizon.rs:7-53`) dispatches the 9 subcommands.
  2. `grep -rn 'io::stdin\|"-"' src/bin/ src/gfa/` confirms `stdin` (`-`) streaming is implemented only in `src/gfa/mod.rs:56,179` (`stream_gfa`, `SpoolableGfaInput`) and `src/bin/horizon_resolve.rs:349-350`, covering the 4 GFA tools (`horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_resolve`), while `horizon_wga`, `horizon_merge`, `horizon_guide`, `horizon_eval`, and `bench_chr1` take file paths.
  3. `horizon_gfa` (`src/bin/horizon_gfa.rs:4`) and `horizon_lift` (`src/bin/horizon_lift.rs:4-5`) import only `gfa` and `coord`, not BOCPD/MDL/e-graph modules.
  4. `/home/agent/.cargo/bin/cargo test --quiet` confirms `45` passed + `1` ignored in `src/lib.rs` plus `14` passed across the `10` binaries (`59` active + `1` ignored).
- **Repair ruling & cost:** The attacker's prescribed wording (`"built on a shared workspace crate"`) is slightly inaccurate because `Cargo.toml:1-5` defines a single `[package]` (`horizon-aligner`), not a Cargo `[workspace]`. Replacing `"shared workspace crate"` with `"shared library crate"` fixes the repair; **the fix is free** (text-only edit in `manuscript/horizon_paper.md:18`).

---

## F2 -- C2: Full `chr6` `451.35 s` speedup is pinned to stale commit `53130b8`, lacks a committed runner script, and uses a `110.7 Mbp` q-arm `CHM13.chr6` extract

- **Ruling:** **UPHELD** (with a minor commit-count correction: `53130b8` is **6** commits behind `HEAD` `33bede8`, not 5).
- **What survives:**
  1. `git log --oneline 53130b8..HEAD` lists 6 commits (`44d281d`, `2695375`, `8c600a0`, `a4e1546`, `f777935`, `33bede8`), and `git diff --stat 53130b8..HEAD -- src/seeding/ src/mdl/dp.rs src/mdl/mod.rs` confirms `509` insertions / `356` deletions across all 5 files (`src/mdl/dp.rs`, `src/mdl/mod.rs`, `src/seeding/flat_index.rs`, `src/seeding/mod.rs`, `src/seeding/translated.rs`).
  2. `scripts/` contains only `benchmarks/`, `biblio/`, and `round`; no committed runner script reproduces the `1e47fe0` baseline or `53130b8` full-`chr6` run, and `manuscript/scripts/generate_figures.py:669-685` hardcodes Figure 5's values.
  3. `head -n 1 /tmp/hprc_chm13_chr6.fa` is `>CHM13.chr6:61401689-172126628` with `110,724,939` bp—the HPRC-extracted q-arm walk, not full ~171 Mbp T2T-CHM13 chromosome 6.
- **Repair ruling & cost:** Disclosing in `manuscript/horizon_paper.md:101,105` that the `451.35 s` (`10,515.0 MB` RSS) run was recorded at `53130b8` and aligns `GRCh38.chr6` (`170.8 Mbp`) against the HPRC q-arm walk `CHM13.chr6:61401689-172126628` (`110.7 Mbp`) **is free**; re-running full `chr6` on `HEAD` (`33bede8`) is not free (~7.5 min).

---

## F3 -- C3: `horizon_gfa stats` on `12 GB` HPRC v1.0 GFA reproduces on `HEAD`

- **Ruling:** **NARROWED** (the `SURVIVED` verdict holds only for `horizon_gfa stats`; the two `horizon_gfa extract` sub-claims in `C3` were not re-executed and remain unreproduced on `HEAD`).
- **What survives:** L7's `horizon_gfa stats` reproduction (`104.98 s`, `97.1 MB` RSS, `81,415,956` segments / `3.29 Gbp` / `112,955,105` links / `2,580` paths / `24,456` walks) confirms the stats portion of `C3` (consistent with Round 3 `adjudication-L7.md:38`). However, `C3` also asserts `horizon_gfa extract` runtimes and peak RSS for `GRCh38.chr22` + `CHM13.chr22` (`183.29 s`, `456.4 MB` RSS) and `GRCh38.chr6` + `CHM13.chr6` (`244.93 s`, `1,255.1 MB` RSS, pinned to `b857377` in `docs/BENCHMARK_RESULTS.md:57`), which neither Round 3 (flagged in Round 3 `adjudication-completeness.md:49`) nor Round 4 L7 re-ran.
- **Repair ruling & cost:** Optionally updating `horizon_gfa stats` to `104.98 s` / `97.1 MB` on `33bede8` and noting `horizon_gfa extract` was pinned at `b857377` **is free**; re-running both `horizon_gfa extract` jobs on `HEAD` is not free (~7 min).

---

## F4 -- C4: `horizon_clean` severs `walk_support <= 1` (not only `== 0`) and is a `0`-severed no-op on real `*-seq` GFAs

- **Ruling:** **UPHELD**.
- **What survives:**
  1. `src/bin/horizon_clean.rs:260` literally reads `let sever = (cfg.prune_unsupported && sup == 0) || (sup <= 1 && cp >= cfg.min_cp_bits);`, so `walk_support == 1` discordant edges are also severed when `cp >= cfg.min_cp_bits`.
  2. `src/bin/horizon_clean.rs:101-102,154-176` shows that when `S` records have `*` placeholder sequences (as emitted by `horizon_wga --gfa`), `prefix`/`suffix` are empty and `score_junction_cp_bits` returns `0.0`, yielding a `0`-severed no-op (`docs/BENCHMARK_RESULTS.md:66`). Only two synthetic 60-bp unit tests (`src/bin/horizon_clean.rs:559-626`) exercise severing.
- **Repair ruling & cost:** Updating `manuscript/horizon_paper.md:49` to *"severing low-support (`walk_support <= 1`) discordant edges (and `walk_support == 0` edges under `--prune-unsupported`) into `HB:i:1` Horizon Bubbles on sequence-bearing GFAs"* is accurate and **the fix is free**.

---

## F5 -- C5: *Plasmodium* (`87.7%` / `54.4%` / `66.8%`) and *Cannabis* (`99.95%`) benchmarks ran Python Prototype 1, not `horizon_wga`, and their input files are absent on disk

- **Ruling:** **UPHELD**.
- **What survives:** Both `../pg/prototypes/pf_pv_end_to_end.py:1` and `../pg/prototypes/cannabis_mdl_audit.py:1-21` confirm these benchmarks were run with Python Prototype 1, not Rust `horizon_wga`, and both reference `/home/agent/.gemini/antigravity-cli/brain/c77b5ff5-88ad-43af-973a-c37dbadc186d/scratch/...` paths that do not exist on disk (`No such file or directory`).
- **Repair ruling & cost:** Removing the misattributed claims or explicitly relabeling them in `manuscript/horizon_paper.md:65` (and Abstract/Figure 3) as historical Python Prototype 1 results with unarchived scratch inputs **is free**; re-running them with Rust `horizon_wga` on archived inputs is not free.

---

## F6 -- C6: The `30` randomized e-graph trials do not exist in `horizon-aligner` and test the superseded Python prototype's unitless cost model, not Rust `MdlCostModel`

- **Ruling:** **UPHELD**.
- **What survives:** `grep -rn "30.*trial\|randomized" src/egraph.rs scripts/` returns `0` matches in `horizon-aligner`. Re-running `python3 /home/rcb112/project/gemini/bench/pg/prototypes/egraph_synteny_resolver.py` reproduces the exact dictionary (`{'trials': 30.0, 'greedy_uninverted_block_recall': 0.0, 'greedy_exact_dup_label_recall': 0.0, 'egraph_ancestral_block_recall': 1.0, 'egraph_exact_dup_label_recall': 1.0}`), and `egraph_synteny_resolver.py:161,210-258` confirms it uses unitless parsimony weights `op_cost = {"seq": 1.0, "inv": 1.2, "dup": 1.3}` (not Rust `MdlCostModel`) and a greedy baseline (`line 235`: `if lbl not in seen_labels and st > 0`) that drops inverted blocks by construction when every trial inverts a 2–3-block span (`lines 225-228`).
- **Repair ruling & cost:** Dropping the `30`-trial claim from `manuscript/horizon_paper.md:99` **is free**; implementing and running a `30`-trial benchmark against Rust `src/egraph.rs` (`MdlCostModel`) with an orientation-aware greedy baseline is not free.

---

## Verdict Table

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | Defect upheld (`10` binaries, `stdin` only on 4 GFA tools, `gfa`/`lift` outside MDL/BOCPD); prescribed repair narrowed from `"workspace crate"` to `"library crate"` (free). |
| F2 | UPHELD | `451.35 s` chr6 run is pinned at `53130b8` (`6` commits behind `HEAD`), lacks a committed runner script, and aligns `110.7 Mbp` q-arm `CHM13.chr6:61401689-172126628` (disclosure fix is free). |
| F3 | NARROWED | `SURVIVED` holds only for `horizon_gfa stats` (`104.98 s`, `97.1 MB` RSS, exact counts); both `horizon_gfa extract` timings/RSS in `C3` were not re-run on `HEAD`. |
| F4 | UPHELD | `horizon_clean.rs:260` severs `sup <= 1` (not only `== 0`) and is a `0`-severed no-op on `*-seq` GFAs; prescribed repair is accurate and free. |
| F5 | UPHELD | `87.7%`/`54.4%`/`66.8%` *Plasmodium* and `99.95%` *Cannabis* claims ran Python Prototype 1, not `horizon_wga`, and scratch inputs are absent on disk (relabeling/removal is free). |
| F6 | UPHELD | `30`-trial benchmark is absent from `horizon-aligner`; `../pg/prototypes/egraph_synteny_resolver.py` tests unitless parsimony weights and a strand-filtered greedy baseline (dropping claim is free). |
