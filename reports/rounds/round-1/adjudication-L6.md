# Round 1 Adjudication -- L6 (prior-art)

- **Attack sealed digest verified:** `6ca3f716b537781f67fb917de02f62b2e85d5865c5a86903ce57506400fc815b` matches `6ca3f716b537`.
- **Prior rounds checked:** `./scripts/round history "vg" --kind adjudication` returned no closed rounds (Round 1 is the first round).
- **Method cards:** `horizon-aligner` declares no `refs/methods/` directory; this lens adjudicates prior-art and codebase-implementation claims.

---

## F1 -- C1: Unix-pipe composition prior art and unimplemented binaries
- **Ruling:** **NARROWED** (UPHELD on the defect; OVERTURNED on the prescribed repair).
- **What survives:** `odgi` (Guarracino et al. 2022), `vg` (Garrison et al. 2018), and `gfatools` (Li) already implement modular Unix-pipe composition, and `horizon_clean`, `horizon_resolve`, and `horizon_lift` do not exist in `src/bin/`.
- **Why:**
  - *Diagnosis:* `docs/STRATEGIC_PLAN.md:6` frames pipe composition as an architectural goal ("We are building...") rather than claiming to have invented Unix pipes, so the prior-art omission is real but narrower than a priority claim. Verification of `src/bin/` confirms `horizon_clean`, `horizon_resolve`, and `horizon_lift` do not exist (only `horizon_wga`, `horizon_merge`, `horizon_guide`, `horizon_eval`, and `bench_chr1` exist).
  - *Prescribed repair (OVERTURNED):* The attacker's fix—restricting C1's named tools to `horizon_wga`, `horizon_merge`, and `horizon_guide`—would install a **new false claim**, because none of those three binaries read from `stdin` or compose via Unix pipes (`src/bin/horizon_wga.rs:122` requires file paths; `src/bin/horizon_merge.rs:12,25,32-35` requires `<maf_out>`, `<gfa_out>`, and `<sp=path.maf>` file paths). Any repair (free, text-only edit) must mark `horizon_clean`/`horizon_resolve`/`horizon_lift` as planned/unimplemented and not claim existing binaries compose via pipes.

## F2 -- C2: MDL alignment prior art and paralogous "hairballs"
- **Ruling:** **UPHELD**.
- **What survives:** Three-state MDL alignment is prior art from Allison, Wallace & Yee (1992) (already cited at `manuscript/references.bib:8-17` and `manuscript/horizon_paper.md:36`), and MDL bit-savings does not prevent pangenome "hairballs" caused by genuine high-identity multi-copy homology.
- **Why:** Verified against `manuscript/references.bib:8-17` and `docs/STRATEGIC_PLAN.md:9`. Positive MDL bit-savings (`ΔL = L_null - L_joint`) suppresses low-signal random AT-rich alignment, but genuine multi-copy repeats and paralogs (e.g., HLA, centromeres) yield positive bit-savings and still require graph-level decomposition (`PGGB`/`smoothxg`, Garrison et al. 2024; `Minigraph-Cactus`, Hickey et al. 2024). The proposed text repair is free (and should also cite Allison et al. 1992).

## F3 -- C3: `egg` congruence-closure vs. HLA evolutionary trajectories
- **Ruling:** **UPHELD**.
- **What survives:** `src/egraph.rs` implements syntactic equality saturation (`egg`, Willsey et al. 2021) over `Block`, `Seq`, `Inv`, and `Dup` with tandem collapse `seq(X, X) -> dup(X)`, which does not natively model HLA evolutionary trajectories.
- **Why:** Verified `src/egraph.rs:15-25,165-273`: only double-inversion cancellation, anti-homomorphic reversal, associativity, and tandem `dup(X)` collapse are implemented—no dispersed segmental duplications, gene conversion, or interlocus recombination. The proposed text repair is free.

## F4 -- C4: `vg`/`odgi` RAM usage and `src/gfa/mod.rs` zero-copy claim
- **Ruling:** **UPHELD**.
- **What survives:** `libbdsg` (`PackedGraph`, Eizenga et al. 2020), `odgi` (Guarracino et al. 2022), `GBZ` (Sirén & Paten 2022), and `gfatools` handle gigabase pangenomes without crashing, and `src/gfa/mod.rs` is neither zero-copy nor a complete graph-I/O solution.
- **Why:** Verified `src/gfa/mod.rs:75-116`: `stream_gfa` allocates a `String` line buffer (`read_line`) and calls `.to_string()` and `.as_bytes().to_vec()` on every `S`, `L`, `P`, and `W` record. The proposed text repair is free (and should also drop "solves the I/O bottleneck entirely" and update the `Zero-copy` doc comment in `src/gfa/mod.rs:1`).

## F5 -- C5: `.hal` characterization and Horizon's own `.hal` / custom GFA tags
- **Ruling:** **UPHELD**.
- **What survives:** `.hal` (Hickey et al. 2013) is an open HDF5 schema with GFA/MAF converters and `Minigraph-Cactus` exports GFA/GBZ/VCF; Horizon itself implements a pure-Rust `.hal` reader/writer (`src/hal/`) and emits non-standard GFA tags/comments.
- **Why:** Verified `src/hal/hdf5_read.rs`, `src/hal/hdf5_write.rs`, `manuscript/horizon_paper.md:40,46,92,98`, and `src/formats.rs:15,43,50,102` (which emit custom `HB:i:`, `XS:i:`, `YS:i:`, `ST:Z:` tags and `# EG:Z:` comments). The proposed text repair is free.

## F6 -- C6: Universal `stdin`/`stdout` requirement vs. `src/gfa/mod.rs`, `src/hal/`, and `horizon_wga`
- **Ruling:** **UPHELD**.
- **What survives:** Universal `stdin`/`stdout` streaming is contradicted by `src/gfa/mod.rs:50-64` (`AsRef<Path>`), `src/hal/` (`Read + Seek`), `src/bin/horizon_wga.rs:122` (file paths), and `src/bin/horizon_merge.rs:12,25,32-35` (file paths).
- **Why:** Verified directly in `src/gfa/mod.rs`, `src/bin/horizon_wga.rs`, `src/bin/horizon_merge.rs`, and `manuscript/guide_pure_rust_hal.md:15,20`. Dropping the universal `stdin` requirement in `docs/STRATEGIC_PLAN.md:23` is free (refactoring `stream_gfa` to `BufRead` would be a small code change).

## F7 -- C7: Single-pass linear GFA streaming for bubble extraction
- **Ruling:** **UPHELD**.
- **What survives:** Arbitrary GFA v1.0/v1.1 has no topological ordering guarantee for single-pass bubble extraction, and local bubble decomposition/resolution/stitching is prior art in `vg` snarls (Paten et al. 2018) and `PGGB`/`smoothxg` (Garrison et al. 2024).
- **Why:** Verified against the GFA v1.0/v1.1 specification: `L` lines can reference `S` lines anywhere in the file, so bubble extraction requires pre-sorting (`odgi sort`) or adjacency indexing. The proposed text repair is free.

## F8 -- C8: HPRC graph usability and strategic focus on unimplemented binaries
- **Ruling:** **UPHELD**.
- **What survives:** HPRC graphs are actively used in production with `vg giraffe`, `GBZ`, `odgi`, `GraphAligner`, and `Panacus`, and `horizon_clean`, `horizon_resolve`, and `horizon_lift` do not exist in `src/bin/`.
- **Why:** Verified against `src/bin/` and `manuscript/horizon_paper.md` (all benchmarks use `horizon_wga`, `horizon_merge`, and `horizon_guide`). The proposed text repair is free.

## F9 -- C9: Streaming GFA `L` lines for BOCPD/MDL scoring without `S`-line sequences
- **Ruling:** **UPHELD**.
- **What survives:** GFA `L` lines carry no nucleotide sequences (`src/gfa/mod.rs:20-27`), so BOCPD/MDL edge scoring requires indexing or caching `S`-line sequences first.
- **Why:** Verified `src/gfa/mod.rs:20-27`: `GfaLink` stores only segment names, orientations, and overlap CIGAR (`0M` in bluntified graphs). Depth/support-based edge pruning is also prior art in `Minigraph-Cactus`/`vg clip` and `PGGB`/`smoothxg`. The proposed text repair is free.

## F10 -- C10: Banded DP across graph bifurcations prior art and unimplemented `--graph`
- **Ruling:** **UPHELD**.
- **What survives:** Banded sequence-to-graph DP across bifurcations is prior art in `GraphAligner` (Rautiainen & Marschall 2020), `minigraph` (Li et al. 2020), and `vg giraffe` (Sirén et al. 2021), and `horizon_wga --graph` is not implemented.
- **Why:** Verified `src/bin/horizon_wga.rs:48-118`: `--graph` does not exist and `--gfa` is output-only. Only MDL forward-backward scoring inside graph DP is potentially novel. The proposed text repair is free.

## F11 -- C11: 0-based half-open `W`-line coordinates and graph liftover prior art
- **Ruling:** **UPHELD**.
- **What survives:** 0-based half-open `[SeqStart, SeqEnd)` intervals are standard GFA v1.1 `W`-line semantics, coordinate liftover is prior art in `odgi position`, `halLiftover`, `vg surject`, and `Liftoff`, and `horizon_lift` does not exist.
- **Why:** Verified `src/bin/` (no `horizon_lift`), `src/gfa/mod.rs:100-117` (raw field parsing only), and `src/coord.rs:1-38` (`Interval` is not wired into GFA traversal). The proposed text repair is free.

---

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | `odgi`/`vg`/`gfatools` prior art and missing binaries hold, but the proposed fix would falsely claim `horizon_wga`/`horizon_merge`/`horizon_guide` compose via Unix pipes |
| F2 | UPHELD | Allison et al. (1992) MDL prior art is omitted and MDL does not prevent high-identity paralogous hairballs; text repair is free |
| F3 | UPHELD | `src/egraph.rs` only performs syntactic `seq`/`inv`/`dup` rewrites, not HLA evolutionary modeling; text repair is free |
| F4 | UPHELD | `libbdsg`/`odgi`/`GBZ` handle gigabase graphs and `src/gfa/mod.rs:75-116` allocates per record rather than being zero-copy; text repair is free |
| F5 | UPHELD | Horizon implements `.hal` I/O in `src/hal/` and emits custom GFA tags (`HB:i:`, `# EG:Z:`); text repair is free |
| F6 | UPHELD | `src/gfa/mod.rs`, `src/hal/`, `horizon_wga`, and `horizon_merge` all require file paths rather than `stdin`; dropping the universal `stdin` rule is free |
| F7 | UPHELD | Arbitrary GFA cannot be bubble-extracted in a single linear pass, and local bubble resolution is prior art in `vg` snarls and `smoothxg`; text repair is free |
| F8 | UPHELD | HPRC graphs are used in production with existing tools, and `horizon_clean`/`resolve`/`lift` do not exist in `src/bin/`; text repair is free |
| F9 | UPHELD | GFA `L` lines contain no sequences (`src/gfa/mod.rs:20-27`), so BOCPD/MDL scoring requires indexing/caching `S` lines first; text repair is free |
| F10 | UPHELD | Banded DP across graph bifurcations is prior art in `GraphAligner`/`minigraph`/`vg giraffe`, and `horizon_wga --graph` is unimplemented; text repair is free |
| F11 | UPHELD | 0-based half-open `W`-line coordinates are standard GFA v1.1, liftover is prior art (`odgi position`/`halLiftover`), and `horizon_lift` is unimplemented; text repair is free |
