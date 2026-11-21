# Round 7 Completeness Critique

**Target:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `5cb83fb6a1346a8193f9e4048055638d0a422b0a`)
**Artifacts read:** `reports/rounds/round-7/CRITIC-BRIEF.md`, `reports/rounds/round-7/round.json`, `reports/rounds/round-7/CARRIED.md`, all 5 `ATTACKED-*.md` snapshots (`ATTACKED-README.md`, `ATTACKED-ARCHITECTURE.md`, `ATTACKED-BENCHMARK_RESULTS.md`, `ATTACKED-VGP_INTEGRATION_PLAN.md`, `ATTACKED-horizon_paper.md`), all 7 attack reports (`attack-L1-correctness.md` through `attack-L7-reproduction.md`), and `scripts/benchmarks/bench_random_gc_crc.py`.
**Adjudication status:** Per `CRITIC-BRIEF.md`, no `adjudication-*.md` files have been filed yet in `reports/rounds/round-7/`.

---

## 1. Lens Execution & Filing Audit (`UNRUN` vs `UNCHALLENGED`)

All **7 declared lenses** (`correctness`, `sampling`, `statistics`, `provenance`, `inference`, `prior-art`, `reproduction`) **RAN** and filed substantive attack files (`10,004`–`14,909` bytes each). **Zero lenses are `UNRUN`.**

| Lens | File | Bytes | Findings | C1 | C2 | C3 | C4 | C5 | Meta / Steering |
| :--- | :--- | :---: | :---: | :--- | :--- | :--- | :--- | :--- | :--- |
| L1 Correctness | `attack-L1-correctness.md` | 10,004 | 6 (`F1`–`F6`) | `NARROWED` | `NARROWED` | `REFUTED` | `NARROWED` | `SURVIVED` | `F6` (stale Abstract/`BENCHMARK_RESULTS` HAL stub) |
| L2 Sampling | `attack-L2-sampling.md` | 10,465 | 5 (`F1`–`F5`) | `NARROWED` | `NARROWED` | `REFUTED` | `NARROWED` | `UNCHALLENGED` | — |
| L3 Statistics | `attack-L3-statistics.md` | 12,522 | 4 (`F1`–`F4`) | `UNCHALLENGED` | `NARROWED` | `REFUTED` | `NARROWED` | `UNCHALLENGED` | — |
| L4 Provenance | `attack-L4-provenance.md` | 11,758 | 5 (`F1`–`F5`) | `NARROWED` | `NARROWED` | `REFUTED` | `NARROWED` | `SURVIVED` | — |
| L5 Inference | `attack-L5-inference.md` | 14,909 | 6 (`F1`–`F6`) | `NARROWED` | `NARROWED` | `REFUTED` | `NARROWED` | `SURVIVED` | `F6` (`n=1`, MAF-concat, `minimap2` baseline) |
| L6 Prior-Art | `attack-L6-prior-art.md` | 13,589 | 5 (`F1`–`F5`) | `UNCHALLENGED` | `NARROWED` | `REFUTED` | `NARROWED` | `NARROWED` | — |
| L7 Reproduction | `attack-L7-reproduction.md` | 13,240 | 5 (`F1`–`F5`) | `NARROWED` | `NARROWED` | `REFUTED` | `NARROWED` | `UNCHALLENGED` | — |

- **`UNCHALLENGED` (critic looked, found no concrete defect in-lens):**
  - `C1` under **L3** (no statistical claim in unit-test/clippy count) and **L6** (no prior-art claim in unit-test/clippy count).
  - `C5` under **L2** (no sampling claim in `Wittler 2023` citation), **L3** (no statistical claim), and **L7** (author/year/title/journal match `refs/catalog_index.json`; `3:e87` uncheckable without `../library`).
- **Thinness assessment:** No lens failed to engage, and `C1`–`C4` were uniformly challenged (`NARROWED` or `REFUTED`) across all applicable lenses. However, **L3 (`4` findings)** and **L2/L6/L7 (`5` findings)** stayed strictly inside the quoted fragments of `C1`–`C5`; only L1 (`F6`) and L5 (`F6`) flagged claim-list steering.

---

## 2. Framing Critique & Claim-Bundling Defects

`round.json` has `runner_is_author: true` and `claim_list_author: "runner"`. Four of the five claims (`C1`, `C2`, `C3`, `C5`) suffer from **summary-vs-quote bundling**: the `summary` asserts multiple deliverables or experiments, while `source.quote` points to a single narrow sentence. Attackers anchored on `source.quote`, leaving the remaining bundled assertions completely unexamined:

1. **`C1` (`round.json:14-18`):** Bundles `cargo fmt`, `cargo clippy`, the unified 11-subcommand `horizon` dispatcher, `aa5_keys` DAG fallback test, `horizon_clean` `10M` L-overlap W-split test, and `122+1` unit tests—but quotes only the `122+1` test sentence (`ATTACKED-horizon_paper.md:84`). Attackers verified `cargo test`, `cargo clippy`, and `cargo fmt --check`, while `horizon` CLI dispatch, `aa5_keys`, and `horizon_clean` `10M` L-overlap W-splits were never individually checked.
2. **`C2` (`round.json:27-31`):** Bundles `horizon_filter --net` `chr6` scoring (`553`/`565`) **and** `10k`–`200k` HPRC subgraph baselines (`horizon` vs `gfatools`/`odgi`/`vg`) **and** `5 Mbp` `chr6` WGA baselines (`horizon_wga` vs `FastGA`/`minimap2`/`LASTZ`)—but quotes only the `horizon_filter` `553`-record line (`ATTACKED-BENCHMARK_RESULTS.md:86`). Zero attackers ran `bench_external_graph_baselines.py` or `bench_external_wga_baselines.py`.
3. **`C3` (`round.json:39-43`):** Bundles online CRC (`n_cal = 15`), `estimate_hmm_params` adaptive gap-bit penalties, **and** a multi-seed random-GC null sweep (`18%/30%/45%` GC, `10 kb`, `median=0` across 6 seeds on `5cb83fb`)—but quotes only the obsolete offline `λ̂ = 10.0` paragraph (`ATTACKED-horizon_paper.md:40`). All 7 attackers refuted the stale quote, and zero attackers ran `scripts/benchmarks/bench_random_gc_crc.py` or audited `estimate_hmm_params`.
4. **`C5` (`round.json:63-67`):** Bundles syncing all 9 catalog references + `JC69` and rebuilding `manuscript/horizon_paper.pdf` cleanly—but quotes only `- **Wittler R** (2023)...`. Six of seven attackers checked only `Wittler R (2023)` (L6 checked the full bibliography), and zero attackers verified `manuscript/horizon_paper.pdf`.

---

## 3. Concrete Completeness Gaps

| # | Gap & Artifact | Why It Matters | What Would Close It |
| :--- | :--- | :--- | :--- |
| **G1** | **Refuted `λ̂ = 10.0` / fixed CRC floor persists unchallenged in Figure 1 and `README.md`** (`ATTACKED-horizon_paper.md:22`, `ATTACKED-README.md:8`). | All 7 attackers refuted `C3` at `ATTACKED-horizon_paper.md:40`, but none checked `horizon_paper.md:22` (*"gated by Conformal Risk Control (`λ̂ = 10.0` bits)"*) or `README.md:8` (*"`0.55 log2(N*M)` + fixed `10.0`-bit conformal floor"*, which also writes `+` instead of `max`). Patching only line 40 leaves contradictory, false claims in Figure 1 and `README.md`. (Similarly, `ATTACKED-ARCHITECTURE.md` and `ATTACKED-VGP_INTEGRATION_PLAN.md` were not inspected by any attacker.) | Update `horizon_paper.md:22` and `README.md:8` alongside `horizon_paper.md:40` to match `src/mdl/mod.rs:151,428`, and confirm `ARCHITECTURE.md` / `VGP_INTEGRATION_PLAN.md` consistency. |
| **G2** | **External baseline scripts in `C2` (`scripts/benchmarks/bench_external_graph_baselines.py`, `bench_external_wga_baselines.py`) were never executed or audited on `5cb83fb`** (`ATTACKED-horizon_paper.md:84`, `ATTACKED-BENCHMARK_RESULTS.md:87-88`). | `horizon_paper.md:84` and `BENCHMARK_RESULTS.md:87-88` report head-to-head speed/RSS wins (`horizon_gfa stats` `0.058s` vs `vg`/`odgi`/`gfatools`; `horizon_wga` `1.520s` vs `FastGA` `2.894s` and `minimap2` `2.062s`; `LASTZ` 15s timeout). L6 noted missing `KegAlign`/`vg snarls` baselines, but **no attacker** (not even L7 Reproduction) ran either script on `HEAD` to verify that the numbers reproduce or that `horizon_wga` / `horizon_resolve` / `horizon_clean` produce valid outputs on those inputs. | Run `bench_external_graph_baselines.py` and `bench_external_wga_baselines.py` on `5cb83fb` under L7 and record exact commands, outputs, and `n` replicates. |
| **G3** | **Random-GC `median=0` sweep (`scripts/benchmarks/bench_random_gc_crc.py`) and `estimate_hmm_params` in `C3` were never run or audited** (`round.json:39`, `ATTACKED-horizon_paper.md:61`, `ATTACKED-BENCHMARK_RESULTS.md:13`). | `C3.summary` claims `bench_random_gc_crc.py` (`18%/30%/45%` GC, `10 kb`) yields `median=0` false-positive pairs across 6 seeds on `5cb83fb`, whereas `ATTACKED-horizon_paper.md:61` and `ATTACKED-BENCHMARK_RESULTS.md:13` report `288` pairs (`30%` GC) and `149` pairs (`45%` GC) on single draws. Neither L3 nor L7 ran `bench_random_gc_crc.py` to test whether `median=0` actually holds at `30%` and `45%` GC on `5cb83fb`. | Run `python3 scripts/benchmarks/bench_random_gc_crc.py` on `5cb83fb`, audit `estimate_hmm_params` in `src/mdl/mod.rs`, and reconcile `horizon_paper.md:61` / `BENCHMARK_RESULTS.md:13`. |
| **G4** | **Unverified `C1`/`C5` sub-claims: `horizon` 11-subcommand CLI dispatch (`src/bin/horizon.rs`), `aa5_keys` test, `horizon_clean` `10M` L-overlap W-split test, and `manuscript/horizon_paper.pdf`**. | While `cargo test` passed overall, no attacker verified that `horizon <subcommand>` dispatches all 11 subcommands (including `filter` and `locus`) or inspected `aa5_keys` / `10M` L-overlap test coverage, and no attacker checked `manuscript/horizon_paper.pdf` against `horizon_paper.md` (which still contains the refuted `C3` text and stale `4`-segment HAL stub). | Smoke-test `horizon --help` / all 11 subcommands, cite the exact `aa5_keys` and `10M` unit test names/paths, and rebuild `manuscript/horizon_paper.pdf` after fixing `horizon_paper.md`. |
| **G5** | **Entire Table 1 (`horizon_paper.md:63-72`), `horizon_gfa` / `horizon_lift` / `horizon_merge` / full `chr6` WGA (`horizon_paper.md:76-78`), VGP Koala × Amphioxus / `horizon_locus` (`BENCHMARK_RESULTS.md:89-93`), and carried `R2-INDEX-RAM` / `estimate_dag_hmm_params` items omitted from `C1`–`C5`**. | None of the 8 Table 1 WGA accuracy benchmarks, whole-`chr6` `horizon_wga` (`258.88–304.53 s`), `12 GB` HPRC `horizon_gfa`/`horizon_lift`, `92`-species `horizon_merge`, or `horizon_locus` (`0.26–0.86 s`) were re-verified on `5cb83fb` in Round 7, nor was `R2-INDEX-RAM` (moving full `chr6`/`minimap2` scripts out of `/tmp/`) discharged. | Either re-run/pin these benchmarks on `5cb83fb` in a follow-up round or explicitly scope `horizon_paper.md` / `BENCHMARK_RESULTS.md` to their pinned commits (`8cd5db0`, `2896229`, `42700e1`) and carry `R2-INDEX-RAM`. |
