# Build genome alignments with Horizon

Horizon builds pairwise and multi-species genome alignments from FASTA files. You can build the release binaries from source, align two genomes, merge pairwise alignments, or run a multi-species alignment in one command.

## Build the tools

1. Install Rust (edition 2024).
2. Clone the repository and enter it.
3. Run the release build:

```bash
cargo build --release --bins
```

This creates `horizon_wga`, `horizon_merge`, and `horizon_guide` in `target/release/`.

## Align two genomes

Use `horizon_wga` to align a query FASTA against a target FASTA:

```bash
target/release/horizon_wga target.fa query.fa --maf out.maf --gfa out.gfa
```

By default, `horizon_wga` emits MAF blocks and prints a GFA summary. Add `--maf` or `--gfa` to write those files to disk.

### Common flags for `horizon_wga`

| Flag | Effect |
|---|---|
| `max_contigs` | Optional third positional argument. Limits how many contigs Horizon reads from each FASTA (`0` reads all). |
| `--rbest` | Keeps only the best-scoring query contig per reference contig. |
| `--top-k K` | Aligns each reference contig only against its top `K` query contigs by shared seed count. |
| `--resolve-dups` | Resolves overlapping query blocks using net-style filtering. |
| `--respect-softmask` | Skips lowercase bases during fallback `w=10` and `R/Y` seeding. |
| `--no-ry-seeds` | Disables 1-bit `R/Y` fallback seeding. |

## Merge pairwise alignments

If you already have pairwise MAF files against one reference species, merge them with `horizon_merge`:

```bash
target/release/horizon_merge homo_sapiens 50000 merged.maf merged.gfa \
  pan_troglodytes=chimp.maf \
  mus_musculus=mouse.maf
```

`homo_sapiens` is the reference species name. `50000` is the window size in base pairs for chaining and bubble detection.

## Run a multi-species alignment

Use `horizon_guide` to align three or more genomes directly from FASTA files:

```bash
target/release/horizon_guide --progressive --hal out_prefix 0 \
  homo_sapiens=human.fa \
  pan_troglodytes=chimp.fa \
  papio_anubis=baboon.fa \
  mus_musculus=mouse.fa
```

`out_prefix` sets the output file prefix (`out_prefix.maf`, `out_prefix.gfa`, and `out_prefix.hal`). `0` tells Horizon to read all contigs.

### What each `horizon_guide` flag does

- `--progressive`: Aligns non-hinge sister pairs and projects transitive blocks onto the central hinge species. It also runs two iterative ancestral consensus passes for distant species with under 50% hinge coverage.
- `--hal`: Writes a pure-Rust HDF5 `v0` HAL archive (`out_prefix.hal`) alongside the MAF and GFA files.
- `--rbest`: Default for multi-contig inputs (`max_contigs != 1`). Uses top-1 reciprocal best contig pairing.
- `--top-k K`: Sets a custom top-`K` contig filter.
- `--all-pairs`: Disables top-`K` filtering and tests all contig pairs that meet the seed threshold.
- `--respect-softmask`: Skips lowercase softmasked bases during fallback seeding.

## Inspect a HAL file

You can inspect any HAL file (from Horizon or Cactus) with `horizon_wga`:

```bash
target/release/horizon_wga --hal-info out_prefix.hal
```

Because `horizon_guide --hal` writes standard HDF5 `v0` files, you can also inspect them with `h5ls -r out_prefix.hal`.
