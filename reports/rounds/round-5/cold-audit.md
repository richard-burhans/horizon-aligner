# Round 5 Cold Post-Round Audit

**Verdict:** PASS

**Evidence:** `cargo fmt --check` exit 0; `cargo clippy --all-targets --all-features -- -D warnings` exit 0; `cargo test --all-features` exit 0 (97 passed, 0 failed, 1 ignored); `rg -n "KP-14|BandGrid|1110110110111" manuscript/ README.md` 0 matches; `scripts/round corrections 5` exit 0 (8 applied, 0 not applied); `scripts/round check 5 --full` exit 0 (`drift_reviewed` covers both `manuscript/horizon_paper.md` and `manuscript/horizon_companion.md`); `scripts/round drift 5` exit 1 (`UNEXAMINED manuscript/scripts/generate_figures.py`, corrected without a pre-round snapshot).

## 1. Tree & Gate Checks

All Rust toolchain checks were run in `/home/rcb112/project/gemini/bench/horizon-aligner` with `PATH="$HOME/.cargo/bin:$PATH"`:

- `cargo fmt --check`: **PASS** (`FMT_EXIT:0`).
- `cargo clippy --all-targets --all-features -- -D warnings`: **PASS** (`CLIPPY_EXIT:0`).
- `cargo test --all-features`: **PASS** (`TEST_EXIT:0`; 97 passed, 0 failed, 1 ignored [`hal::tests::real_hal_decode_from_env`] across `src/lib.rs` and 10 binaries).
- `python3 -m py_compile manuscript/scripts/generate_figures.py`: **PASS** (`PY_COMPILE_OK`).

## 2. Banned Identifier Removal

- `rg -n "KP-14|BandGrid|1110110110111" manuscript/ README.md`: **PASS** (`NO_MATCHES`, 0 hits across `manuscript/horizon_paper.md`, `manuscript/horizon_companion.md`, `manuscript/scripts/generate_figures.py`, and `README.md`).

## 3. Round Script Checks

- `scripts/round corrections 5`: **PASS** (`EXIT:0`; `COR-1`..`COR-8` all `APPLIED`, `8 correction(s); 0 not applied or malformed`).
- `scripts/round check 5 --full`: **PASS** (`EXIT:0`; notes that `manuscript/horizon_paper.md` [+17/-16] and `manuscript/horizon_companion.md` [+22/-32] are `REVIEWED` against 17 and 22 quoted lines respectively, and notes open Round 4 drift binding on `manuscript/horizon_paper.md`).
- `scripts/round drift 5`: **EXIT:1** — `drift_reviewed` covers both `manuscript/horizon_paper.md` and `manuscript/horizon_companion.md`, so neither manuscript blocks `check 5 --full`; `drift 5` exits 1 solely because `manuscript/scripts/generate_figures.py` is flagged `UNEXAMINED` (corrected without a pre-round snapshot).

## 4. Diff & Derived-Quantity Inspection

Inspected `git diff` across `README.md` (`+7/-7`), `manuscript/horizon_paper.md` (`+17/-16`), `manuscript/horizon_companion.md` (`+22/-32`), and `manuscript/scripts/generate_figures.py` (`+18/-18`):

- `KP-14` → canonical `fwd.min(rc)` / `w=11/12` / `MASK_W10` across all four files; `1110110110111` → `MASK_W10` (`11011011011011`, span 14).
- `BandGrid` → `BandLayout` across all four files, with `BandLayout` clamping `W = min(band_w, max(N, M))` and proportional center `c(i) = floor(i · M / max(N, 1))` reflected in `manuscript/horizon_companion.md` (lines 87–101, 126–129) and `manuscript/scripts/generate_figures.py` (`N = M` noted in Fig C2A title).
- Derived quantities verified against live text:
  - Unit test count in `manuscript/horizon_paper.md` (`97` active + `1` ignored = `98` total) matches `cargo test --all-features`.
  - `GRCh38.chr6` × `CHM13.chr6` precision/multiplicity figures (`107,836,669 / 221,734,653 = 48.63%`, `107,836,669 / 282,163,046 = 38.22%`, `60,428,393 / 282,163,046 = 21.42%`) are arithmetically consistent.
  - HPRC `10k`-segment subgraph `'N'` fraction (`10,510,000 / 10,937,960 = 96.09%`, `427,960` bp ACGT) and `vg stats` speedup ratio (`4.2984 / 0.08 = 53.73×` → `53.7×`) are consistent between prose and `generate_figures.py`.
  - `18%`-GC null suppression (`seed=999`: `18,077` → `0`, `4/5` other seeds `80–448`) and `horse` (`65.4%` → `67.6%`), `macaque` (`84.6%` → `87.0%`), `marmoset` (`86.1%` → `87.7%`) match between `manuscript/horizon_paper.md` and `manuscript/scripts/generate_figures.py`.

## 5. Scope & What Was Not Run

- Did not re-run `generate_figures.py` to regenerate PNG/PDF/SVG assets (only syntax-checked with `py_compile`).
- Did not re-run whole-chromosome `chr6` or HPRC `12 GB` benchmarks; relied on `cargo test --all-features` and diff/arithmetic verification.
- Did not read `attack-*.md`, `adjudication-*.md`, or `completeness-critic.md`.
