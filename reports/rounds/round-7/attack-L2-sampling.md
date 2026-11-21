# Round 7 Adversarial Attack -- L2 (Sampling)

**Scope & Environment Notes:**
- **Literature (`../library`):** No `../library` checkout is present on this host, so all literature checks were performed from `refs/catalog_index.json` and `manuscript/references.bib` alone.
- **Method cards:** As documented in `BRIEF-sampling.md`, this project keeps no method cards yet.
- **Prior rounds checked:** Queried `scripts/round history` for `n_cal`, `553`, `chrM`, and `Wittler` prior to filing; prior adjudications are cited explicitly below with what changed on `4fd2240` / `5cb83fb`.

---

## F1 -- `C1`: `123` total tests rely on a single `55 KB` external `.hal` fixture (`HORIZON_TEST_HAL`) and synthetic/toy unit fixtures
- **Claim:** `C1` (`ATTACKED-horizon_paper.md:84`)
- **Verdict:** `NARROWED`
- **Concrete flaw:** Passing `122` active unit tests plus `1` `#[ignore]`d `.hal` test (`src/hal/tests.rs:240-260`) samples only `n = 1` external `55 KB` HDF5 `.hal` file (`halTest.hal`, `4` synthetic genomes `Genome_0..Genome_3`, `root_members=8`), which is skipped under default `cargo test` unless `HORIZON_TEST_HAL` is set to an external path outside the repository. None of the `123` tests validate `horizon_guide`'s own progressive `.hal` output (`Anc0..Anc3`) or any multi-contig/chromosome-scale real `.hal` file.
- **Evidence:**
  - `/home/agent/.cargo/bin/cargo test --quiet` ran `88` passed, `0` failed, `1` ignored (`src/lib.rs`) + `34` binary tests (`122` active).
  - `HORIZON_TEST_HAL=/home/rcb112/project/claude/bench/data/iuc-port/tools-iuc/tools/haltools/test-data/halTest.hal /home/agent/.cargo/bin/cargo test --lib real_hal_decode_from_env -- --ignored --nocapture` passed `1` test (`src/hal/tests.rs:242`).
  - `./target/release/horizon_wga --hal-info .../halTest.hal` reports `hdf5_v=0 offset_bytes=8 len_bytes=8 root_members=8`, `sample_members=Genome_0,Genome_1,Genome_2,Genome_3,Genomes,Meta,Phylogeny,Verison` (`55K` file).
- **Correction:** Qualify `C1` in `manuscript/horizon_paper.md:84` to state that the `1` `.hal` integration test (`#[ignore]`d without `HORIZON_TEST_HAL`) exercises a single `55 KB` 4-genome fixture (`halTest.hal`) rather than multi-contig or `horizon_guide`-emitted `.hal` files.

## F2 -- `C2`: `chr6` `horizon_filter` scoring (`100%` preservation of `109,119,624` bp) evaluates `n = 1` q-arm-only `1-to-1` truth (`60.23M–170.74M`), where all `15` removed records lie in the unmeasured p-arm (`0–45.7 Mbp`)
- **Claim:** `C2` (`ATTACKED-BENCHMARK_RESULTS.md:86`, `ATTACKED-horizon_paper.md:84`)
- **Verdict:** `NARROWED`
- **Concrete flaw:** The claim that SDUST+MDL removes `12` noisy/non-1-to-1 records while "preserving `100%` of exact `109,119,624` bp coordinates" is a sampling artifact of the `n = 1` truth set (`/tmp/hprc_chr6_truth_segments.tsv`). First, `/tmp/hprc_chr6_truth_segments.tsv` only covers `GRCh38.chr6:60,229,947–170,743,762` (the q-arm, because `CHM13.chr6` in HPRC v1.0 is `61,401,689–172,126,628`), and `scripts/benchmarks/score_paf_vs_hprc.py:34` filters to `mult == "1"` (`1,628,018` of `1,628,056` segments, `110,082,990` bp). Second, comparing `/tmp/chr6_fastga_net_head.paf` (`565` records) and `/tmp/chr6_fastga_mdl_net_head.paf` (`553` records) by `(qname, qstart, qend, strand, tname, tstart, tend)` shows `15` removed coordinate intervals (not `12`, as `3` records changed during netting), and **every single removed record** has `GRCh38.chr6` `qstart` between `147,635` and `45,643,002` (p-arm). Because the truth sample has zero coverage below `60,229,947`, removing p-arm alignments cannot reduce exact recall on the `109,119,624` bp q-arm truth. Additionally, `score_paf_vs_hprc.py` cannot be re-executed on this host (`ModuleNotFoundError: No module named 'numpy'`).
- **Evidence:**
  - `/tmp/hprc_chr6_truth_segments.tsv`: `GRCh38.chr6` range `60,229,947 .. 170,743,762`; `mult==1`: `1,628,018` (`110,082,990` bp), `mult!=1`: `38` (`310` bp).
  - Coordinate diff between `/tmp/chr6_fastga_net_head.paf` (`565`) and `/tmp/chr6_fastga_mdl_net_head.paf` (`553`): `15` removed records, all with `qstart` in `147,635 .. 45,643,002` (`< 60,229,947`).
  - `python3 scripts/benchmarks/score_paf_vs_hprc.py --paf /tmp/chr6_fastga_mdl_net_head.paf ...` failed with `ModuleNotFoundError: No module named 'numpy'`.
- **Correction:** Narrow `C2` in `docs/BENCHMARK_RESULTS.md:86` and `manuscript/horizon_paper.md:84` to state that `100%` coordinate preservation applies only to the `1-to-1` `GRCh38.chr6` q-arm (`60.23M–170.74M`) truth subset (`n = 1` chromosome pair, `n = 1` aligner `FastGA`), and that all MDL-removed netted records fall in the unmeasured `GRCh38.chr6` p-arm (`< 45.7 Mbp`).

## F3 -- `C3`: `horizon_paper.md:40` claims an offline constant `λ̂ = 10.0` bits, whereas `src/mdl/mod.rs:401–428` runs online `calibrate_null_crc_for_gc` with `n_cal = 15` 64-bp i.i.d. draws from a `4,096`-bp prefix sample that is `100%` `'N'` on `GRCh38.chr6`
- **Claim:** `C3` (`ATTACKED-horizon_paper.md:40`)
- **Verdict:** `REFUTED`
- **Prior round context:** Round 4 (`adjudication-L3.md:62-67`, `F7`) narrowed `n_cal = 15` CRC as underpowered (`1/(15+1) = 0.0625`). **What changed:** `4fd2240` / `5cb83fb` added `calibrate_null_crc_for_gc` (`src/mdl/dp.rs:120-165`) and wired it into `src/mdl/mod.rs:401-428`, while `ATTACKED-horizon_paper.md:40` still describes an offline constant `λ̂ = 10.0` bits.
- **Concrete flaw:**
  1. `ATTACKED-horizon_paper.md:40` is factually false against `HEAD`: `src/mdl/mod.rs:428` no longer enforces a constant offline `λ̂ = 10.0` bits; it calls `dp::calibrate_null_crc_for_gc(observed_gc, 15, 0.07, ctx.effective_min_bits, &hmm_params)` dynamically.
  2. The implemented GC sampler (`src/mdl/mod.rs:407-412`) reads only `ctx.seq_x.iter().take(4096).chain(ctx.seq_y.iter().take(4096))`. On `/tmp/hprc_grch38_chr6.fa` (`170.8 Mbp`), the first `4,096` bp are **`4,096` `'N'`s and `0` ACGT bases** (`0.000%` of `GRCh38.chr6`), so `observed_gc` is sampled solely from the first `4,096` bp of `CHM13.chr6` (`0.0037%` of `CHM13.chr6`, at the pericentromeric q-arm boundary).
  3. `calibrate_null_crc_for_gc` (`src/mdl/dp.rs:120-165`) draws `n_cal = 15` synthetic i.i.d. windows of fixed `win_len = 64` bp from a deterministic xorshift64 PRNG. Because `(n_ge + 1) / (15 + 1) <= 0.07` forces `n_ge = 0` (`1/16 = 0.0625`), `λ̂` is always the **sample maximum** of `15` short i.i.d. draws, giving maximum sensitivity to a single draw and ignoring variable gap lengths and genomic autocorrelation/repeats.
- **Evidence:**
  - `src/mdl/mod.rs:401-428` and `src/mdl/dp.rs:86-165`.
  - Prefix check on `/tmp/hprc_grch38_chr6.fa`: `GRCh38.chr6 first 4096 bp: N count = 4096 ACGT = 0`; `/tmp/hprc_chm13_chr6.fa`: `N count = 0 ACGT = 4096`.
- **Correction:** Replace `manuscript/horizon_paper.md:40` to document the actual online `calibrate_null_crc_for_gc` procedure (`n_cal = 15`, `win_len = 64` i.i.d. draws, `4,096`-bp prefix GC sample), and replace the `4,096`-bp prefix sampler with chromosome-wide ACGT sampling (and increase `n_cal` beyond the `n_ge = 0` boundary).

## F4 -- `C4`: Progressive HAL accuracy is evaluated only on `n = 1` `~16.5 kb` mitochondrial `chrM.maf` (`15–16 kb` per pair) with `hg38` recall dropping to `71.8%` (`mm10`) and `74.1%` (`canFam3`), while the `5`-mammal nuclear run (`0.859 s`) has no accuracy verification
- **Claim:** `C4` (`ATTACKED-horizon_paper.md:84`)
- **Verdict:** `NARROWED`
- **Prior round context:** Round 5 (`adjudication-L2.md:35-40`, `F4`) upheld that `horizon_paper.md` omitted MAF-block concatenation (`11–387 KB` slices). **What changed:** `5cb83fb` added `scripts/benchmarks/verify_merge_guide_maf.py` and the `5`-mammal `chrM.maf` progressive HAL run (`2,103` bottom segments).
- **Concrete flaw:** The only accuracy-verified progressive HAL run samples a single `~16.5 kb` intronless mitochondrial genome (`../pg/data/benchmarks/chrM.maf`, `n = 1` tree, `n = 1` contig, `15,025–15,980` truth pairs per species pair) lacking nuclear repeats, rearrangements, and multi-contig scaffolding, and recall already drops to `71.8%` (`mm10`) and `74.1%` (`canFam3`). Meanwhile, the `0.859 s` `5`-mammal nuclear run (`hg38`, `Macaca_mulatta`, `HLcalPym1`, `Equus_caballus`, `Canis_lupus`, `1,109` MAF blocks) uses concatenated MAF-extracted slices with no MAF accuracy verification.
- **Evidence:**
  - `python3 scripts/benchmarks/verify_merge_guide_maf.py` output: `hg38 x panTro4: truth=15980 recall=0.955 prec=0.998`; `hg38 x rheMac3: truth=15025 recall=0.892 prec=0.971`; `hg38 x mm10: truth=15700 recall=0.718 prec=0.975`; `hg38 x canFam3: truth=15575 recall=0.741 prec=0.982` (`elapsed_s=0.3237`).
- **Correction:** Narrow `C4` in `manuscript/horizon_paper.md:84` to state explicitly that progressive HAL accuracy is verified only on `n = 1` single-contig `~16.5 kb` `chrM` (`15–16 kb` truth pairs) and that the `5`-mammal `0.859 s` nuclear slice has no ground-truth recall/precision measurement.

## F5 -- `C5`: Pure bibliographic citation entry (`Wittler R (2023)`) verified in `references.bib` and `catalog_index.json`
- **Claim:** `C5` (`ATTACKED-horizon_paper.md:150`)
- **Verdict:** `UNCHALLENGED`
- **Concrete flaw:** None under the sampling lens; `C5` is a bibliographic reference entry (`manuscript/references.bib:481-489`, `refs/catalog_index.json:1430-1433`, `doi:10.24072/pcjournal.323`) making no sampling claim. (Full text not checked because `../library` is absent.)

---

## Claim Verdict Table

| claim_id | verdict | summary |
| :--- | :--- | :--- |
| C1 | NARROWED | `122+1` tests rely on `n = 1` external `55 KB` `halTest.hal` (`#[ignore]`d by default) and do not test `horizon_guide` HAL output. |
| C2 | NARROWED | `100%` coordinate preservation holds only on `1-to-1` `GRCh38.chr6` q-arm truth (`60.23M–170.74M`); all `15` removed records lie in the unmeasured p-arm (`< 45.7 Mbp`). |
| C3 | REFUTED | `horizon_paper.md:40` claims offline constant `λ̂ = 10.0` bits, whereas `src/mdl/mod.rs:407-428` runs online `n_cal = 15` 64-bp i.i.d. CRC using a `4,096`-bp prefix GC sample (`100%` `'N'` on `GRCh38.chr6`). |
| C4 | NARROWED | Progressive HAL accuracy is verified only on `n = 1` `~16.5 kb` `chrM.maf` (recall drops to `71.8%`/`74.1%` on mouse/dog); `5`-mammal nuclear run (`0.859 s`) has no accuracy check. |
| C5 | UNCHALLENGED | Pure bibliographic citation (`Wittler R 2023`, verified in `references.bib` and `catalog_index.json`); no sampling claim. |
