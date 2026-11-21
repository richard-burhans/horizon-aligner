# Round 7 Adversarial Attack — L1 Correctness

**Repository:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `5cb83fb`)
**Lens:** L1 Correctness
**Note on Rules 2 & 3:** `../library` is not present on disk (`NO_LIBRARY_DIR`), so prior-work checks worked from `refs/catalog_index.json` and `manuscript/references.bib` alone without asserting held paper contents from memory. Per the project map, this repository keeps **no method cards**.

---

## F1 — C3 (`manuscript/horizon_paper.md:40`): Horizon does NOT calibrate CRC offline or enforce a constant `λ̂ = 10.0` bits floor

- **Claim:** C3 (`manuscript/horizon_paper.md:40`) — *"Using fixed Conformal Risk Control defaults calibrated offline at risk level `α = 0.07` (`n_cal = 15`; Angelopoulos et al. 2024), Horizon enforces a constant bit-savings floor `λ̂ = 10.0` bits (`min_delta_bits = 10.0`) and a minimum column posterior floor (`0.15`)."*
- **Verdict:** **REFUTED**
- **Prior-round context:** Round 4 (`adjudication-L3.md [F7]`), Round 5 (`adjudication-L3.md [F6]`), and Round 6 (`adjudication-L3.md [F6]`) previously noted uncalibrated `n_cal = 15` / `min_delta_bits = 10.0` as claim-list steering omissions. On `HEAD` (`5cb83fb`), `src/mdl/mod.rs` and `src/mdl/dp.rs` now implement online GC-conditioned calibration, making C3's "calibrated offline" and "constant bit-savings floor `λ̂ = 10.0` bits" statement factually false.
- **Evidence:**
  1. `src/mdl/mod.rs:401-429`: `extend_anchors` computes `observed_gc` online from the first `4,096` bp of `seq_x` and `seq_y` and calls `dp::calibrate_null_crc_for_gc(observed_gc, 15, 0.07, ctx.effective_min_bits, &hmm_params)` at runtime (not offline).
  2. `src/mdl/mod.rs:151`: `let effective_min_bits = cfg.min_delta_bits.max(KARLIN_MDL_COEFF * search_space_bits);` (with `KARLIN_MDL_COEFF = 0.55`). For any sequence pair with `N * M > 2^(10 / 0.55) ≈ 297,298` bp², `effective_min_bits` exceeds `10.0` bits (e.g., `~14.6` bits at `10 kb × 10 kb`, `~29.8` bits on `chr6`), and `calibrate_null_crc_for_gc` (`src/mdl/dp.rs:86-165`) returns `λ̂ >= base_floor_bits`, which is neither offline nor constant at `10.0` bits.
- **Correction:** Replace `manuscript/horizon_paper.md:40` with: *"At runtime, Horizon estimates `observed_gc` from the first `4,096` bp of each sequence pair and runs online Conformal Risk Control calibration (`n_cal = 15`, `α = 0.07`, `64`-bp synthetic decoy windows; Angelopoulos et al. 2024) floored by `max(min_delta_bits, 0.55 * log2(N * M))` (`min_delta_bits = 10.0`), alongside a minimum column posterior floor (`0.15`)."*

## F2 — C4 (`manuscript/horizon_paper.md:84`): Commit conflation, Abstract/BENCHMARK_RESULTS contradiction, and progressive HAL pectinate-tree / reference-clone defects

- **Claim:** C4 (`manuscript/horizon_paper.md:84`) — *"`horizon_guide --progressive --hal --rbest` (`2896229` / `HEAD`) aligned `5` mammals ... in `0.859 s` (`1,109` MAF blocks, `967` multi-species) and on 5-mammal UCSC `chrM.maf` (`verify_merge_guide_maf.py`, `0.258 s`) emitted a progressive ancestral-node HDF5 `v0` HAL (`((((hg38,panTro4)Anc0,rheMac3)Anc1,mm10)Anc2,canFam3)Anc3;`, `2,103` bottom segments) with `95.5%/99.8%` (`panTro4`), `89.2%/97.1%` (`rheMac3`), `71.8%/97.5%` (`mm10`), and `74.1%/98.2%` (`canFam3`) MAF recall/precision."*
- **Verdict:** **NARROWED**
- **Prior-round context:** Round 5 (`adjudication-L7.md [F5]`) and Round 6 (`adjudication-completeness.md [2.2]`) flagged the missing `chrM` MAF accuracy verification and star-tree HAL stub; `5cb83fb` (`HEAD`) added `verify_merge_guide_maf.py` and `from_multi_contig_progressive_blocks`.
- **Evidence:**
  1. **Reproduced `chrM.maf` metrics:** Running `python3 scripts/benchmarks/verify_merge_guide_maf.py` on `HEAD` (`5cb83fb`) reproduces `hal_tree=((((hg38,panTro4)Anc0,rheMac3)Anc1,mm10)Anc2,canFam3)Anc3;`, `hal_bottom_segs=2103`, `elapsed_s=0.3404`, and exact MAF recall/precision: `hg38 x panTro4` (`0.955`/`0.998`), `hg38 x rheMac3` (`0.892`/`0.971`), `hg38 x mm10` (`0.718`/`0.975`), `hg38 x canFam3` (`0.741`/`0.982`).
  2. **Commit conflation & internal contradiction:** `from_multi_contig_progressive_blocks` was introduced in `5cb83fb` (`HEAD`), not `2896229`. Meanwhile, `manuscript/horizon_paper.md:10` (Abstract) and `docs/BENCHMARK_RESULTS.md:11,58` still report the obsolete `2896229` *"single-contig star-tree HDF5 `v0` HAL stub (`4` bottom segments)"*, directly contradicting C4.
  3. **Correctness defects in progressive HAL (`src/hal/mod.rs:250-365`):** `from_multi_contig_progressive_blocks` builds a hardcoded pectinate (caterpillar) tree in **CLI argument order** (`((...((ref_name,q0)Anc0,q1)Anc1...),q_{k-1})Anc_{k-1};`) rather than an inferred guide tree, and every ancestral node `Anc0..Anc_{n_q-1}` clones `ref_genome.dna` and `ref_genome.bottom_segments` (`src/hal/mod.rs:302-324`) with identity `1-to-1` `top_segments`. Moreover, `hal_bottom_segs=2103` is the pre-export unpartitioned `bottom_segments.len()` from overlapping pairwise MAF runs (`src/bin/horizon_guide/main.rs:339-410`), and `verify_merge_guide_maf.py:139-156` only scores the `.maf` output (never decoding `.hal` beyond `horizon_wga --hal-info`).
- **Correction:** Update `manuscript/horizon_paper.md:10` and `docs/BENCHMARK_RESULTS.md:11,58` to match `HEAD` (`5cb83fb`), separate the `2896229` `5`-mammal run from the `5cb83fb` `chrM.maf` run in C4, and qualify the HAL output as a reference-framed pectinate ancestral-node HDF5 `v0` HAL (`Anc0..Anc3` cloning hinge coordinates/DNA) whose MAF accuracy numbers (`95.5%/99.8%`, etc.) are measured on the emitted `.maf`, not decoded from `.hal`.

## F3 — C1 (`manuscript/horizon_paper.md:84`): `halTest.hal` is absent from the repository so `real_hal_decode_from_env` is ignored on clean `HEAD`

- **Claim:** C1 (`manuscript/horizon_paper.md:84`) — *"All `122` active Rust unit tests (`88` library + `34` binary, plus `1` `.hal` integration test verified against `halTest.hal` via `HORIZON_TEST_HAL`; `123` total) pass on `HEAD` with zero `clippy` warnings."*
- **Verdict:** **NARROWED**
- **Prior-round context:** Round 3 (`F3`) and Round 6 (`adjudication-L4.md [F1]`, `adjudication-L7.md [F1]`) flagged stale unit-test counts (`46`, `120`) and `cargo fmt` diffs. On `HEAD` (`5cb83fb`), `/home/agent/.cargo/bin/cargo test --all-features --no-fail-fast` passes `88` lib + `34` bin = `122` active tests (`1` ignored), `/home/agent/.cargo/bin/cargo clippy --all-targets --all-features -- -D warnings` exits `0`, and `/home/agent/.cargo/bin/cargo fmt --check` exits `0`.
- **Evidence:** `find . -name "halTest.hal"` returns nothing in the workspace, and `src/hal/tests.rs:240-243` marks `real_hal_decode_from_env` as `#[ignore = "requires HORIZON_TEST_HAL=/path/to/file.hal"]`. On a clean checkout of `HEAD`, `122` active tests pass and `1` `.hal` test is skipped (`1 ignored`).
- **Correction:** Qualify C1 to note that `real_hal_decode_from_env` is `#[ignore]`d by default and requires an external `HORIZON_TEST_HAL` fixture not tracked in the repository.

## F4 — C2 (`docs/BENCHMARK_RESULTS.md:86`): Pre-netting `436` MDL count conflated with post-netting `12`-record delta, and `score_paf_vs_hprc.py` fails without `numpy`

- **Claim:** C2 (`docs/BENCHMARK_RESULTS.md:86`) — *"`553` records (`436` filtered by MDL; removes `12` noisy/non-1-to-1 records while preserving `100%` of exact `109,119,624` bp coordinates: **99.12%** exact recall, **99.65%** cond. / **99.32%** uncond. exact precision)"*
- **Verdict:** **NARROWED**
- **Evidence:**
  1. Running `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --query /tmp/hprc_grch38_chr6.fa --target /tmp/hprc_chm13_chr6.fa --net` produces `[horizon_filter] total=56857 kept=553 filtered_mdl=436 filtered_len=0 elapsed=4.82s`.
  2. Pure-Python scoring of `/tmp/chr6_fastga_mdl_net_head.paf` (`553` records) against `/tmp/hprc_chr6_truth_segments.tsv` (`110,082,990` bp 1-to-1) confirms `exact_pairs=109,119,624`, `exact_recall=0.991249` (`99.12%`), `cond_prec=0.996497` (`99.65%`), and `uncond_prec=0.993205` (`99.32%`).
  3. However, `553 records (436 filtered by MDL; removes 12 noisy/non-1-to-1 records...)` conflates pre-netting `filtered_mdl=436` (out of `56,857` raw records) with the post-netting reduction (`565 - 553 = 12`). Additionally, `python3 scripts/benchmarks/score_paf_vs_hprc.py` fails on this host with `ModuleNotFoundError: No module named 'numpy'`.
- **Correction:** Clarify that `436` records are filtered by MDL *prior to netting* (out of `56,857` raw records), yielding `553` netted records (`12` fewer than sequence-free `--net`'s `565`).

## F5 — C5 (`manuscript/horizon_paper.md:150`): Citation matches `refs/catalog_index.json` and `manuscript/references.bib`

- **Claim:** C5 (`manuscript/horizon_paper.md:150`) — *"- **Wittler R** (2023). General encoding of canonical *k*-mers. *Peer Community Journal*, 3:e87."*
- **Verdict:** **SURVIVED**
- **Evidence:** `refs/catalog_index.json:1430-1434` (`slug: "wittler_2023_canonical-kmer-encoding"`, `cite: "Wittler, R. (2023). *General encoding of canonical k-mers*. Peer Community Journal."`, `doi: "10.24072/pcjournal.323"`), `manuscript/references.bib:481-489`, and `manuscript/horizon_paper.md:28,150` all match. (Worked from catalog index alone because `../library` is absent.)

## F6 — Meta-finding: Claim-list steering omits stale `4`-segment HAL stub claims in Abstract and `BENCHMARK_RESULTS.md`

- **Evidence:** The runner-authored claim list `C1`–`C5` targets `manuscript/horizon_paper.md:84` (`C4`) while omitting `manuscript/horizon_paper.md:10` (Abstract) and `docs/BENCHMARK_RESULTS.md:11,58`, which still claim `horizon_guide --progressive --hal` emits a *"single-contig star-tree HDF5 `v0` HAL stub (`4` bottom segments)"* that no longer exists on `HEAD` (`5cb83fb`).

---

| claim_id | verdict | finding |
| :--- | :--- | :--- |
| C1 | NARROWED | F3 |
| C2 | NARROWED | F4 |
| C3 | REFUTED | F1 |
| C4 | NARROWED | F2, F6 |
| C5 | SURVIVED | F5 |
