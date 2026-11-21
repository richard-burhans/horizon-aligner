# Round 7 Adjudication — L1 Correctness

**Repository:** `/home/rcb112/project/gemini/bench/horizon-aligner` (`HEAD` = `5cb83fb`)  
**Attack Verified:** `reports/rounds/round-7/attack-L1-correctness.md` (`sha256` prefix `a2ee89b3a596`, matching sealed digest)  
**Method Cards:** None present in this repository (`refs/methods/` does not exist); rulings rest on direct code/test execution and `refs/catalog_index.json` / `manuscript/references.bib`.

---

## F1 — C3 (`manuscript/horizon_paper.md:40`): Horizon does NOT calibrate CRC offline or enforce a constant `λ̂ = 10.0` bits floor

- **Ruling:** **UPHELD** (diagnosis and proposed correction both upheld).
- **What survives:** C3's claims of offline CRC calibration and a constant `λ̂ = 10.0` bits floor are factually false on `HEAD` (`5cb83fb`), because `src/mdl/mod.rs:401-429` computes `observed_gc` at runtime and calls `dp::calibrate_null_crc_for_gc` (`src/mdl/dp.rs:120-165`), floored by `effective_min_bits = cfg.min_delta_bits.max(0.55 * search_space_bits)` (`src/mdl/mod.rs:151`).
- **Prior-round context:** Rounds 4 (`adjudication-L3.md [F7]`), 5 (`adjudication-L3.md [F6]`), and 6 (`adjudication-L3.md [F6]`) flagged uncalibrated `n_cal = 15` / `min_delta_bits = 10.0` as claim-list steering omissions; `HEAD` now implements online GC-conditioned calibration, invalidating C3's static description.
- **Whole-set check & repair cost:** `manuscript/horizon_paper.md:22` (Figure 1 caption) also repeats `λ̂ = 10.0` bits and should be updated alongside line 40. **The fix is free** (text-only update from verified source code, no re-measurement needed).

## F2 — C4 (`manuscript/horizon_paper.md:84`): Commit conflation, Abstract/BENCHMARK_RESULTS contradiction, and progressive HAL pectinate-tree / reference-clone defects

- **Ruling:** **UPHELD** (endorsing the attacker's **NARROWED** assessment of C4 and its proposed repair).
- **What survives:** Re-running `python3 scripts/benchmarks/verify_merge_guide_maf.py` on `HEAD` (`5cb83fb`) reproduces `hal_tree=((((hg38,panTro4)Anc0,rheMac3)Anc1,mm10)Anc2,canFam3)Anc3;`, `hal_bottom_segs=2103`, `total_blocks=266`, `multi_species_blocks=255` (`elapsed_s=0.3256`), and exact MAF recall/precision (`panTro4` `0.955`/`0.998`, `rheMac3` `0.892`/`0.971`, `mm10` `0.718`/`0.975`, `canFam3` `0.741`/`0.982`); however, `from_multi_contig_progressive_blocks` was added in `5cb83fb` (not `2896229`), `manuscript/horizon_paper.md:10` and `docs/BENCHMARK_RESULTS.md:11,58` still state the obsolete `4`-segment star-tree stub, `src/hal/mod.rs:250-365` builds a CLI-order pectinate tree with `Anc0..Anc_{n_q-1}` cloning `ref_genome`'s coordinates/DNA from unpartitioned pairwise MAF runs (`src/bin/horizon_guide/main.rs:339-420`), and `verify_merge_guide_maf.py:139-156` scores `.maf` rather than decoding `.hal`.
- **Prior-round context:** Round 5 (`adjudication-L5.md [F3]`, `adjudication-L7.md [F5]`) and Round 6 (`adjudication-completeness.md`) flagged the `4`-segment star-tree stub and missing `chrM` MAF accuracy check; `5cb83fb` implemented progressive HAL emission and `verify_merge_guide_maf.py`, leaving the stale Abstract/`BENCHMARK_RESULTS.md` text and structural/scoring caveats.
- **Repair cost:** **The fix is free** (text-only updates to `manuscript/horizon_paper.md:10,84` and `docs/BENCHMARK_RESULTS.md:11,58`).

## F3 — C1 (`manuscript/horizon_paper.md:84`): `halTest.hal` is absent from the repository so `real_hal_decode_from_env` is ignored on clean `HEAD`

- **Ruling:** **UPHELD** (endorsing the attacker's **NARROWED** assessment of C1 and its proposed qualification).
- **What survives:** `/home/agent/.cargo/bin/cargo test --all-features --no-fail-fast` passes `88` lib + `34` bin = `122` active tests with `1` ignored (`real_hal_decode_from_env` at `src/hal/tests.rs:240-243`), `cargo clippy --all-targets --all-features -- -D warnings` exits `0`, and `find . -name "halTest.hal"` confirms `halTest.hal` is not tracked in the workspace.
- **Prior-round context:** Round 3 (`F3`) and Round 6 (`adjudication-L4.md [F1]`, `adjudication-L7.md [F1]`) flagged stale test counts (`46`, `120`) and `cargo fmt`; those are resolved on `HEAD`.
- **Repair cost:** **The fix is free** (minor parenthetical clarification in `manuscript/horizon_paper.md:84`).

## F4 — C2 (`docs/BENCHMARK_RESULTS.md:86`): Pre-netting `436` MDL count conflated with post-netting `12`-record delta, and `score_paf_vs_hprc.py` fails without `numpy`

- **Ruling:** **UPHELD** (endorsing the attacker's **NARROWED** assessment of C2 and its proposed clarification).
- **What survives:** Executing `./target/release/horizon_filter --paf /tmp/chr6_fastga.paf --query /tmp/hprc_grch38_chr6.fa --target /tmp/hprc_chm13_chr6.fa --net` yields `total=56857 kept=553 filtered_mdl=436 filtered_len=0 elapsed=4.88s`, while sequence-free `--net` yields `total=56857 kept=565 filtered_mdl=1`, confirming `436` is pre-netting (out of `56,857`) whereas `12` (`565 - 553`) is post-netting; `python3 scripts/benchmarks/score_paf_vs_hprc.py` also fails on this host with `ModuleNotFoundError: No module named 'numpy'`.
- **Whole-set check & repair cost:** `manuscript/horizon_paper.md:84` already avoids the `436` conflation (`553 records; MDL removed 12 noisy/non-1-to-1 records`), so only `docs/BENCHMARK_RESULTS.md:86` needs clarification. **The fix is free**.

## F5 — C5 (`manuscript/horizon_paper.md:150`): Citation matches `refs/catalog_index.json` and `manuscript/references.bib`

- **Ruling:** **UPHELD** (C5 **SURVIVED**).
- **What survives:** `refs/catalog_index.json:1430-1434`, `manuscript/references.bib:481-489`, and `manuscript/horizon_paper.md:28,150` all match for Wittler (2023), *Peer Community Journal* 3:e87 (`doi:10.24072/pcjournal.323`), with `../library` absent on disk.

## F6 — Meta-finding: Claim-list steering omits stale `4`-segment HAL stub claims in Abstract and `BENCHMARK_RESULTS.md`

- **Ruling:** **UPHELD**.
- **What survives:** `manuscript/horizon_paper.md:10` and `docs/BENCHMARK_RESULTS.md:11,58` still report the obsolete `4`-segment star-tree HAL stub while the claim list targeted only `manuscript/horizon_paper.md:84` (`C4`).
- **Repair cost:** **The fix is free** (text-only update aligned with `HEAD` `5cb83fb`).

---

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | `src/mdl/mod.rs:151,401-429` and `src/mdl/dp.rs:120-165` run online GC-conditioned CRC floored by `0.55 * log2(N * M)`, refuting C3's offline/constant-`10.0` claim; fix (including Fig. 1 caption at line 22) is free |
| F2 | UPHELD | `chrM.maf` metrics reproduce, but `5cb83fb` vs `2896229` provenance, Abstract/`BENCHMARK_RESULTS.md` `4`-seg contradiction, CLI-order pectinate reference-clone HAL, and `.maf`-only scoring are verified; fix is free |
| F3 | UPHELD | `122` active (`88` lib + `34` bin) pass and `clippy` exits `0`, while `halTest.hal` is absent from the repo and `real_hal_decode_from_env` (`src/hal/tests.rs:240-243`) is `#[ignore]`d; fix is free |
| F4 | UPHELD | `horizon_filter` confirms `436` is pre-netting (vs `565 - 553 = 12` post-netting) and `score_paf_vs_hprc.py` fails with `No module named 'numpy'`; `docs/BENCHMARK_RESULTS.md:86` clarification is free |
| F5 | UPHELD | Wittler (2023) citation in `manuscript/horizon_paper.md:150` matches `refs/catalog_index.json:1430-1434` and `manuscript/references.bib:481-489` |
| F6 | UPHELD | `manuscript/horizon_paper.md:10` and `docs/BENCHMARK_RESULTS.md:11,58` still carry the superseded `4`-segment HAL stub claim omitted from `C1`–`C5`; fix is free |
