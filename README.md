# horizon-aligner

A pure-Rust **information-theoretic genome alignment and pangenome graph suite** combining **strand-symmetric spaced-seed indexing (`FlatSeedIndex`)**, **Bayesian Online Change-Point Detection (BOCPD)** compositional segmentation, **Minimum Description Length (MDL)** banded & DAG dynamic programming with explicit epistemic abstention (**Horizon Bubbles**), **MDL-calibrated e-graph paralog resolution**, **inversion-aware VCF (`SVTYPE=INV`) & bidirected `rGFA` construction**, and **streaming `GFA v1.0/v1.1` + pure-Rust zero-`libhdf5` `HAL`** interoperability.

## Key Capabilities

- **Pairwise & DAG/POA MDL alignment (`horizon_wga`):** Canonical `fwd.min(rc)` spaced-seed indexing (`FlatSeedIndex`, 12 bytes/entry + 16-bit prefix buckets) with SDUST + optional `--tantan` tandem-repeat masking and `--respect-softmask`, Rayon-parallel `LASTZ12`/`LASTZ11` transition seeding, `w=10`/`w=9`, `MASK_RY18`, and 6-frame `AA5` fallback seeding, optional Baum–Welch `--train-hmm`, flat `BandLayout` Pair-HMM Forward–Backward MEA alignment with continuous `md:f:`/`mb:f:` bit-score PAF tags, HOXD70 affine micro-gap realignment, $O(R \log R)$ species-wide `--one-to-one` netting, native UCSC `--chain` and hierarchical `--net` output, and seeded DAG/POA alignment (`horizon_wga --graph`).
- **Calibrated epistemic abstention (Horizon Bubbles):** Intervals falling below the online biological + GC-matched Conformal Risk Control floor or `min_post_prob` abstain explicitly as `# HB:Z:` / `HB:i:1` Horizon Bubbles rather than forcing low-confidence base pairings.
- **Streaming `12 GB+` GFA, rGFA & PAF pangenome toolkit:**
  - Transparent plain, `.gz`/`.bgz`, `.zst`/`.zstd`, and `stdin`/`stdout` streaming across GFA, PAF, FASTA, MAF, VCF, and `.hnex`.
  - `horizon_gfa`: Single-pass counter-only BGZF/gzip/zstd `stats`, lightweight `qc`, `list-walks`, streaming two-pass chromosome `extract`, streaming MAF-to-variation-GFA `from-maf`, and parallel incremental bidirected `rgfa` construction (with balanced-inversion strand-flip edges).
  - `horizon_lift`: Graph-aware (`--gfa`) and UCSC `.chain` (`--chain`) BED/GFF3 coordinate liftover with `--bridge-bp` micro-bubble coalescing and reverse-strand arithmetic.
  - `horizon_clean`: Two-pass `S`-flank-indexed BOCPD+MDL junction cleaner (`--min-cp-bits`, `--min-haplotypes`, `--clip-bubbles-bp`).
  - `horizon_resolve`: Standalone MDL-calibrated `k`-copy e-graph resolution (`MdlCostModel`) over `u32` `PackedGraph` + $O(V + E)$ iterative Tarjan biedged bridge decomposition with `--collapse-dups`, `--tag-walks`, and phased haplotype `--vcf` / `--vcf-ref` / `--nested-snarls` (`LV=0`/`LV=1`, `SVTYPE=INV`) VCFv4.2 export.
  - `horizon_filter`: Streaming SDUST + Pair-HMM MDL PAF filter with biological CRC (`md:f:`/`hb:i:` tags) and symmetric reciprocal-best netting (`--net --net-score mdl|log-length-ani`).
  - `horizon_locus`: Regional multi-FASTA all-pairs locus aligner emitting `.similarity.tsv`, `.bionj.nwk`, `.maf`, `.gfa`, optional lossless `.hnex`, ProbCons consistency MSA (`--consistency`), and codon frameshift repair (`--repair-frameshifts`).
- **Multi-species merging, progressive guide-tree WGA & reference-free TBA (`horizon_merge`, `horizon_guide`, `horizon_tba`):** $O(N \log N)$ sweep-line MAF/rGFA merging, `auN`-weighted MinHash hinge selection, BFS UPGMA sister rescue, and reference-free Threaded Blockset Alignment with all-vs-all, two-tier `--hinges`/`--auto-hinges`/`--extend-maf`, bottom-up `--progressive-ancestors` (`--ancestral-model fitch|ml`, `--hal`), packed-RAM `--variation-gfa`, direct gzip output (`--out-gz`), Rayon POA block refinement, and `--pairwise-dir` atomic checkpointing.
- **Pure-Rust zero-`libhdf5` HDF5 / HAL streaming toolkit & writer (`horizon_hal`):** Fail-closed streaming `Read + Seek` HDF5 `v0`/`v1` parser with $O(\log C)$ chunk B-tree range descent (`stats [--chroms]`, arbitrary-depth `liftover`, `to-maf`, `to-gfa`) and `libhdf5`-compatible `v0` `4 MiB`-windowed streaming Fitch progressive HAL writer.

## Installation & Verification

Requires Rust (`edition = "2024"`, Rust 1.85+):

```bash
cargo build --release --bins --locked
cargo test --all-features --locked
cargo clippy --all-targets --all-features --locked -- -D warnings
RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features --locked
cargo fmt --check
```

## CLI Quick Reference

All tools can be invoked directly (`target/release/horizon_<tool>`) or through the unified `horizon <subcommand>` dispatcher:

```bash
# Pairwise WGA (MAF + GFA + PAF + native UCSC .chain/.net, 1-to-1 netting, e-graph paralog resolution)
horizon wga target.fa query.fa 0 --one-to-one --respect-softmask \
  --maf out.maf --gfa out.gfa --paf out.paf --chain out.chain --net out.net --resolve-dups

# UCSC .chain or GFA-walk BED/GFF3 coordinate liftover
horizon lift --chain out.chain --bed regions.bed --out lifted.bed --unmapped unmapped.bed
horizon lift --gfa pangenome.gfa.gz --bed regions.bed --from GRCh38#0#chr22 --to CHM13#0#chr22 --bridge-bp 10

# Sequence-to-graph DAG/POA alignment
horizon wga --graph target.gfa.gz query.fa

# Streaming GFA stats, QC, chromosome extraction, MAF-to-variation-GFA, & incremental bidirected rGFA
horizon gfa stats pangenome.gfa.gz
horizon gfa qc pangenome.var.gfa.gz
horizon gfa extract pangenome.gfa.gz GRCh38#0:chr6:chr6.fa
horizon gfa from-maf multi.maf.gz --out pangenome.var.gfa.gz
horizon gfa rgfa --ref ref.fa --query q1.fa --query q2.fa --min-sv-bp 50 --respect-softmask -t 8 --out out.rgfa

# Pure-Rust zero-libhdf5 streaming multi-level HAL toolkit
horizon hal stats alignment.hal --chroms RefSpecies
horizon hal liftover alignment.hal --src SpeciesA --dst SpeciesB --bed region.bed --out lifted.bed
horizon hal to-maf alignment.hal --ref SpeciesA --seq chr1 --start 1000000 --length 10000 --out slice.maf
horizon hal to-gfa alignment.hal --ref SpeciesA --seq chr1 --start 1000000 --length 10000 --out slice.gfa

# Junction cleaning, e-graph resolution, and phased nested-snarl VCFv4.2 (+ SVTYPE=INV) calling
horizon clean --gfa in.gfa.gz --min-cp-bits 12.0 --min-haplotypes 2 --out clean.gfa
horizon resolve in.gfa.gz --collapse-dups --tag-walks --vcf variants.vcf.gz --nested-snarls --out resolved.gfa

# Streaming PAF MDL filtering & regional locus alignment
horizon filter --paf raw.paf.zst --query q.fa --target t.fa --net --net-score log-length-ani --out filtered.paf.zst
horizon locus --fasta locus.fa --out-prefix locus_out --threads 4 --hnex --consistency --repair-frameshifts

# Multi-species progressive guide-tree WGA & reference-free TBA
horizon merge ref_species 1000 out.maf out.gfa sp1=pair1.maf sp2=pair2.maf
horizon guide --progressive --hinge sp1 --rbest --hal -t 16 out_prefix 0 sp1=sp1.fa sp2=sp2.fa sp3=sp3.fa
horizon tba --pairwise-dir pw_mafs/ --pairwise-gz --variation-gfa --out-gz -t 16 out_prefix 0 sp1=sp1.fa sp2=sp2.fa sp3=sp3.fa
horizon tba --progressive-ancestors --ancestral-model ml --variation-gfa --hal -t 16 out_prefix 0 sp1=sp1.fa sp2=sp2.fa sp3=sp3.fa
```

## Documentation

- [Technical Manual](docs/manual.md) — Processing stages, MDL/BOCPD theory, and CLI usage for all 13 suite binaries.
- [Architecture Map](ARCHITECTURE.md) & [ADRs](docs/decisions/) — Tier map, module contracts, dependency flow, and architecture decision records (`ADR-001`–`ADR-008`).
- [Pangenome Suite](docs/PANGENOME_SUITE.md) — Design and contracts of `horizon_gfa`, `horizon_lift`, `horizon_clean`, `horizon_resolve`, `horizon_filter`, `horizon_locus`, and `horizon_wga --graph`.

## License & Citation

Released under the [MIT License](LICENSE). See [CITATION.cff](CITATION.cff) for citation metadata.
