# Round 3 Attack -- Reproduction (L7)

Target claims: `C1`–`C6` from `reports/rounds/round-3/BRIEF-reproduction.md`, tested at commit `488b020` (`HEAD`).
Prior round history checked via `./scripts/round history "" --kind adjudication` and `./scripts/round history "129"`: none of the specific reproduction discrepancies below (`129 target scaffolds` vs `31`, `45` vs `46` passing tests, missing `/tmp/*.fa` slice-extraction scripts, and broken `run_all_evals.py` data path) were adjudicated in Rounds 1 or 2.

---

## F1 -- C2 claims `129 target scaffolds` when the extracted Dog FASTA contains `31` contigs, and `cactus241way_chr22_1MB.maf` is absent from benchmark scripts

- **Claim:** `C2` (`docs/BENCHMARK_RESULTS.md:16`) -- "UCSC cactus241way `chr22 1MB`: **62.3%** recall / **60.3%** precision (Dog, across 129 target scaffolds; **61.3% F1** in **17.16s**)."
- **Verdict:** **NARROWED**
- **Evidence:**
  1. Re-running extraction and evaluation on `/home/rcb112/project/gemini/bench/pg/data/benchmarks/cactus241way_chr22_1MB.maf`:
     - `python3 scripts/benchmarks/extract_fasta_from_maf.py /home/rcb112/project/gemini/bench/pg/data/benchmarks/cactus241way_chr22_1MB.maf hg38 Canis_lupus_familiaris /tmp/repro_hg38_zoo_1m.fa /tmp/repro_dog_zoo_1m.fa /tmp/repro_pw_zoo_1m.maf`
       Output: `Extracted 1 ref contigs and 31 target contigs.`
     - `grep -c "^>" /tmp/dog_zoo_1m.fa /tmp/repro_dog_zoo_1m.fa` -> both files have `31` contigs, not `129`.
     - `./target/release/horizon_eval hg38 Canis_lupus_familiaris /tmp/repro_hg38_zoo_1m.fa /tmp/repro_dog_zoo_1m.fa /tmp/repro_pw_zoo_1m.maf`
       Output: `elapsed=21.520s proto_pairs=136348 maf_pairs=132089 overlap=82228 raw_maf_recall=62.3% raw_maf_precision=60.3% bubbles=3490`.
  2. The recall (`62.3%`), precision (`60.3%`), F1 (`61.3%`), and bubble count (`3,490`) reproduce exactly, but the `129 target scaffolds` figure does not exist in either the pre-existing `/tmp/dog_zoo_1m.fa` or a fresh extraction from `cactus241way_chr22_1MB.maf`.
  3. Neither `scripts/benchmarks/fetch_benchmarks.sh` nor `scripts/benchmarks/run_all_evals.py` fetches or runs `cactus241way_chr22_1MB.maf` (the only reference is an ad-hoc parent-directory script `/home/rcb112/project/gemini/bench/run_zoo_1mb.py` pointing to a stale path).
- **Correction:** Change `across 129 target scaffolds` in `docs/BENCHMARK_RESULTS.md:16` to `across 31 target contigs`, and add the `cactus241way_chr22_1MB.maf` slice command to `scripts/benchmarks/fetch_benchmarks.sh` / `run_all_evals.py`.

---

## F2 -- C3 miscounts the passing unit test suite as `45` instead of `46`

- **Claim:** `C3` (`docs/STRATEGIC_PLAN.md:3`) -- "All target pangenome binaries (`horizon_wga`, `horizon_wga --graph`, `horizon_merge`, `horizon_guide`, `horizon_eval`, `horizon_gfa`, `horizon_clean`, `horizon_resolve`, `horizon_lift`, `bench_chr1`) are now implemented in `src/bin/` with unit tests in every binary (`45` passing tests + `1` ignored; `cargo clippy --all-targets -- -D warnings` clean)."
- **Verdict:** **NARROWED**
- **Evidence:**
  - `cargo clippy --all-targets -- -D warnings` exited `0` with no warnings.
  - `cargo test` exited `0` with:
    - `src/lib.rs`: `37 passed; 0 failed; 1 ignored` (`hal::tests::real_hal_decode_from_env`)
    - 9 binaries in `src/bin/` (`bench_chr1`, `horizon_clean`, `horizon_eval`, `horizon_gfa`, `horizon_guide`, `horizon_lift`, `horizon_merge`, `horizon_resolve`, `horizon_wga`): `1 passed` each (`9` total).
  - Total passing tests = `37 + 9 = 46` passing + `1` ignored (`cargo test -- --list | grep ": test$" | wc -l` = `47`), not `45` passing + `1` ignored.
- **Correction:** Update `docs/STRATEGIC_PLAN.md:3` from `` `45` passing tests + `1` ignored `` to `` `46` passing tests + `1` ignored ``.

---

## F3 -- C4 and C5 rely on untracked `/tmp/*.fa` slices with no extraction script, and C4 misstates the `CHM13.chr6` coordinate window

- **Claims:**
  - `C4` (`docs/STRATEGIC_PLAN.md:11`) -- "`ContigIndex` (`src/seeding/mod.rs`) now uses a compact `FlatSeedIndex` (`12 bytes/entry` sorted `Vec<u64>` + `Vec<u32>` with on-the-fly reverse-complement coordinates, cutting `20 Mbp chr6 q-arm` peak RSS from `5,863.3 MB` to `2,409.7 MB` and `chr6` 5 Mbp pericentromere RSS from `3,084.5 MB` to `1,048.2 MB`)."
  - `C5` (`docs/BENCHMARK_RESULTS.md:23`) -- "Real HPRC `chr22` Pericentromere 5 Mbp (`GRCh38.chr22:15.5M-20.5M` vs `CHM13.chr22:15.73M-20.73M`): **14.57s**, **1,861.9 MB** peak RSS, **5,106,541** many-to-many `(x, y)` aligned pairs (`unique_x_bp=4,387,655`, `unique_y_bp=4,549,858`), `16,083` Horizon Bubbles."
- **Verdicts:** **NARROWED** (`C4`), **NARROWED** (`C5`)
- **Evidence:**
  1. Re-running `./target/release/horizon_wga <x.fa> <y.fa> --resolve-dups` on the pre-existing `/tmp/` FASTAs reproduces the exact deterministic counts and closely matches peak RSS:
     - `chr6_q20m` (`/tmp/grch38_chr6_q20m.fa`, `/tmp/chm13_chr6_q20m.fa`): `dt=50.17s`, `rss=2425.3 MiB` (`2543.1 MB`) vs claimed `2,409.7 MB`; `aligned_pairs=15025501`, `unique_x_bp=14999823`, `unique_y_bp=14998488`, `total_bubbles=16877`.
     - `chr6_peri_5m` (`/tmp/grch38_chr6_peri_5m.fa`, `/tmp/chm13_chr6_peri_5m.fa`): `dt=12.98s`, `rss=1077.1 MiB` (`1129.4 MB`) vs claimed `1,048.2 MB`; `aligned_pairs=4247511`, `unique_x_bp=4012961`, `unique_y_bp=3998723`, `total_bubbles=23299`.
     - `chr22_peri_5m` (`/tmp/grch38_chr22_peri_5m.fa`, `/tmp/chm13_chr22_peri_5m.fa`): `dt=17.90s` vs `14.57s`, `max_rss=1788.6 MiB` (`1831528 KB`) vs `1,861.9 MB`; `aligned_pairs=5106541`, `unique_x_bp=4387655`, `unique_y_bp=4549858`, `total_bubbles=16083`.
  2. However, `grep -rn "chr6_q20m\|chr22_peri_5m\|chr6_peri_5m" src/ scripts/` returns `No matches`: the 6 FASTA inputs exist only as untracked files in `/tmp/` with no reproducible extraction script from `hprc-v1.0-mc-grch38.gfa.gz`.
  3. Moreover, `head -n 1 /tmp/chm13_chr6_q20m.fa` is `>CHM13.chr6_61M_81M`, whereas `docs/BENCHMARK_RESULTS.md:26,62` reports `CHM13.chr6:61.4M-81.4M`.
  4. Note also that `FlatSeedIndex::build` (`src/seeding/mod.rs:386-441`) temporarily holds `pairs: Vec<(u64, u32)>` (16 bytes/entry) concurrently with `keys: Vec<u64>` and `pos: Vec<u32>` (12 bytes/entry) before `drop(pairs)` at line 441, so peak build-time index memory is ~28 bytes/entry, not 12 bytes/entry.
- **Correction:** Add a slice-extraction script under `scripts/benchmarks/` that reproduces the `/tmp/*_{peri_5m,q20m}.fa` inputs from `hprc-v1.0-mc-grch38.gfa.gz`, and fix `CHM13.chr6:61.4M-81.4M` to `CHM13.chr6:61M-81M` in `docs/BENCHMARK_RESULTS.md:26,62`.

---

## F4 -- `C1` and `C6` survive direct reproduction

- **Claims:**
  - `C1` (`docs/STRATEGIC_PLAN.md:13`) -- `stream_gfa` buffered line-streaming GFA v1.0/v1.1 reader (`.gfa`, `.gfa.gz` via `MultiGzDecoder`, `"-"` `stdin`, unit tests for empty `L`-line orientations and `.gfa.gz`).
  - `C6` (`docs/BENCHMARK_RESULTS.md:19,54`) -- `horizon_gfa stats`: `88.84s`, `97.1 MB` peak RSS across `81,415,956` segments (`3.29 Gbp`), `112,955,105` links, `2,580` paths, and `24,456` walks.
- **Verdicts:** **SURVIVED** (`C1`), **SURVIVED** (`C6`)
- **Evidence:**
  - `C1`: Verified in `src/gfa/mod.rs:97-127,232-329`; `cargo test` passes both `gfa::tests::streams_all_record_types_and_handles_empty_link_orients` and `gfa::tests::streams_gzipped_gfa_file`; `printf "H\tVN:Z:1.0\nS\ts1\tACGT\nS\ts2\tTGCA\nL\ts1\t+\ts2\t+\t0M\n" | ./target/release/horizon_gfa stats -` succeeded (`segments=2 segment_bp=8 links=1`).
  - `C6`: `./target/release/horizon_gfa stats /home/rcb112/project/gemini/bench/pg/data/benchmarks/hprc-v1.0-mc-grch38.gfa.gz` produced `elapsed=244.01s peak_rss_mb=96.7 segments=81415956 segment_bp=3287932785 links=112955105 paths=2580 walks=24456` (exact counts and `96.7 MB` vs `97.1 MB` RSS; slower wall time `244.01s` vs `88.84s` under concurrent benchmark load).

---

## F5 -- Framing & steering: `run_all_evals.py` points to a nonexistent data directory, and `C2` omits the 5 stale pre-`488b020` WGA benchmarks

- **Artifact:** `scripts/benchmarks/run_all_evals.py:10`, `docs/BENCHMARK_RESULTS.md:5,10-17,35-47`
- **Evidence:**
  1. `scripts/benchmarks/run_all_evals.py:10` sets `BENCH_DIR = ROOT.parent / "data" / "benchmarks"`, which does not exist (`/home/rcb112/project/gemini/bench/data/benchmarks: No such file or directory`); the actual downloaded benchmarks reside in `/home/rcb112/project/gemini/bench/pg/data/benchmarks`. Running `python3 scripts/benchmarks/run_all_evals.py` out of the box skips every dataset.
  2. `docs/BENCHMARK_RESULTS.md:5` notes that `multiz100way`, `cactus447way`, `multiz470way`, `cactus241way chr22_random`, and `GRCh38.chr22 vs CHM13.chr22` were recorded prior to `FlatSeedIndex` / `rev_offsets` fixes. The runner-authored claim list selected only `cactus241way chr22 1MB` (`C2`)—the single WGA benchmark re-run on `488b020`—steering away from the 5 un-reproduced WGA rows.

---

## Verdict Table

| Claim | Verdict |
| :--- | :--- |
| C1 | SURVIVED |
| C2 | NARROWED |
| C3 | NARROWED |
| C4 | NARROWED |
| C5 | NARROWED |
| C6 | SURVIVED |
