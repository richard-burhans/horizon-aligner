#!/usr/bin/env python3
"""Validate MAF-FASTA coordinate/sequence invariants and GFA walk/link invariants."""
import argparse
import gzip
from pathlib import Path

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")


def load_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def validate_maf(maf_path: Path, fastas: dict[str, str], max_errors: int = 20) -> tuple[int, int]:
    errors = 0
    rows_checked = 0
    with open_text(maf_path) as f:
        for lineno, line in enumerate(f, 1):
            if not line.startswith("s "):
                continue
            parts = line.split()
            if len(parts) < 7:
                errors += 1
                if errors <= max_errors:
                    print(f"[MAF ERROR] line {lineno}: malformed s-line")
                continue
            src, start_s, size_s, strand, src_size_s, seq = parts[1:7]
            start, size, src_size = int(start_s), int(size_s), int(src_size_s)
            ungapped = seq.replace("-", "").upper()
            rows_checked += 1
            if len(ungapped) != size:
                errors += 1
                if errors <= max_errors:
                    print(f"[MAF ERROR] line {lineno} ({src}): ungapped len {len(ungapped)} != size {size}")
                continue
            if start < 0 or start + size > src_size:
                errors += 1
                if errors <= max_errors:
                    print(f"[MAF ERROR] line {lineno} ({src}): {start}+{size} > srcSize {src_size}")
                continue
            # Match FASTA by exact name or suffix after first '.'
            contig_seq = fastas.get(src)
            if contig_seq is None and "." in src:
                contig_seq = fastas.get(src.split(".", 1)[1])
            if contig_seq is None and len(fastas) == 1:
                contig_seq = next(iter(fastas.values()))
            if contig_seq is not None:
                if len(contig_seq) != src_size:
                    errors += 1
                    if errors <= max_errors:
                        print(f"[MAF ERROR] line {lineno} ({src}): fasta len {len(contig_seq)} != srcSize {src_size}")
                    continue
                if strand == "+":
                    expected = contig_seq[start:start + size]
                else:
                    fwd_end = src_size - start
                    fwd_start = fwd_end - size
                    expected = revcomp(contig_seq[fwd_start:fwd_end])
                if ungapped != expected:
                    errors += 1
                    if errors <= max_errors:
                        print(f"[MAF ERROR] line {lineno} ({src}): sequence mismatch at {start}..{start+size} ({strand})")
    return rows_checked, errors


def validate_gfa(gfa_path: Path, max_errors: int = 20) -> tuple[int, int, int]:
    segs: set[str] = set()
    links: set[tuple[str, str, str, str]] = set()
    walks: list[tuple[int, str, list[tuple[str, str]]]] = []
    errors = 0
    with open_text(gfa_path) as f:
        for lineno, line in enumerate(f, 1):
            if not line or line.startswith("#") or line.startswith("H\t"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            tag = parts[0]
            if tag == "S":
                segs.add(parts[1])
            elif tag == "L":
                u, uo, v, vo = parts[1], parts[2], parts[3], parts[4]
                links.add((u, uo, v, vo))
                rev_uo = "-" if uo == "+" else "+"
                rev_vo = "-" if vo == "+" else "+"
                links.add((v, rev_vo, u, rev_uo))
            elif tag == "W" and len(parts) >= 7:
                walk_str = parts[6]
                steps = []
                i = 0
                while i < len(walk_str):
                    if walk_str[i] in "><":
                        orient = "+" if walk_str[i] == ">" else "-"
                        j = i + 1
                        while j < len(walk_str) and walk_str[j] not in "><":
                            j += 1
                        steps.append((walk_str[i + 1:j], orient))
                        i = j
                    else:
                        i += 1
                walks.append((lineno, parts[1], steps))
    for lineno, sample, steps in walks:
        for seg, _ in steps:
            if seg not in segs:
                errors += 1
                if errors <= max_errors:
                    print(f"[GFA ERROR] line {lineno} ({sample}): missing segment {seg}")
        for (u, uo), (v, vo) in zip(steps, steps[1:]):
            if (u, uo, v, vo) not in links:
                errors += 1
                if errors <= max_errors:
                    print(f"[GFA ERROR] line {lineno} ({sample}): missing link {u}{uo} -> {v}{vo}")
    return len(segs), len(walks), errors


def main():
    parser = argparse.ArgumentParser(description="Validate MAF-FASTA and GFA invariants")
    parser.add_argument("--maf", type=Path, help="MAF file to validate")
    parser.add_argument("--fastas", type=Path, nargs="*", default=[], help="FASTA files for MAF slice verification")
    parser.add_argument("--gfa", type=Path, help="GFA file to validate")
    args = parser.parse_args()

    total_errors = 0
    if args.maf:
        fastas = {}
        for fp in args.fastas:
            fastas.update(load_fasta(fp))
        rows, errs = validate_maf(args.maf, fastas)
        print(f"MAF validated ({args.maf}): rows={rows:,} errors={errs}")
        total_errors += errs
    if args.gfa:
        n_segs, n_walks, errs = validate_gfa(args.gfa)
        print(f"GFA validated ({args.gfa}): segs={n_segs:,} walks={n_walks} errors={errs}")
        total_errors += errs
    if total_errors > 0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
