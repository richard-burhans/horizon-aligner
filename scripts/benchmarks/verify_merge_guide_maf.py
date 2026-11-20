#!/usr/bin/env python3
"""Verify horizon_guide (--progressive --hal) and horizon_merge MAF accuracy against UCSC chrM.maf."""
import argparse
import os
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def parse_maf_blocks(path: Path):
    cur = []
    with open(path) as f:
        for line in f:
            if line.startswith("s "):
                cur.append(line.split())
            elif line.strip() == "" and cur:
                yield cur
                cur = []
    if cur:
        yield cur


def extract_species_fastas_and_truth(
    maf_path: Path, species: list[str], out_dir: Path
) -> tuple[dict[str, Path], dict[tuple[str, str], set[tuple[int, int]]]]:
    seqs = {sp: [] for sp in species}
    offsets = {sp: 0 for sp in species}
    ref_sp = species[0]
    truth_pairs: dict[tuple[str, str], set[tuple[int, int]]] = {
        (ref_sp, sp): set() for sp in species[1:]
    }
    for blk in parse_maf_blocks(maf_path):
        by_sp = {}
        for row in blk:
            sp = row[1].split(".")[0]
            if sp in seqs and sp not in by_sp:
                by_sp[sp] = row[6].upper()
        # Record truth pairs in concatenated coordinates when ref_sp and tgt_sp are both in block
        if ref_sp in by_sp:
            r_seq = by_sp[ref_sp]
            for tgt_sp in species[1:]:
                if tgt_sp not in by_sp:
                    continue
                t_seq = by_sp[tgt_sp]
                rp = offsets[ref_sp]
                tp = offsets[tgt_sp]
                for rb, tb in zip(r_seq, t_seq):
                    if rb != "-" and tb != "-":
                        truth_pairs[(ref_sp, tgt_sp)].add((rp, tp))
                    if rb != "-":
                        rp += 1
                    if tb != "-":
                        tp += 1
        for sp, aln_s in by_sp.items():
            ug = aln_s.replace("-", "")
            seqs[sp].append(ug)
            offsets[sp] += len(ug)
    paths = {}
    for sp, chunks in seqs.items():
        if not chunks:
            continue
        p = out_dir / f"{sp}.fa"
        p.write_text(f">{sp}.chrM\n{''.join(chunks)}\n")
        paths[sp] = p
    return paths, truth_pairs


def maf_pairs(maf_path: Path, ref_sp: str, tgt_sp: str) -> set[tuple[int, int]]:
    pairs = set()
    for blk in parse_maf_blocks(maf_path):
        r_row = next((r for r in blk if r[1].split(".")[0] == ref_sp), None)
        t_row = next((r for r in blk if r[1].split(".")[0] == tgt_sp), None)
        if not r_row or not t_row:
            continue
        r_pos = int(r_row[2])
        r_strand = r_row[4]
        r_size = int(r_row[5])
        t_pos = int(t_row[2])
        t_strand = t_row[4]
        t_size = int(t_row[5])
        for rb, tb in zip(r_row[6], t_row[6]):
            rc = r_pos if r_strand == "+" else (r_size - 1 - r_pos)
            tc = t_pos if t_strand == "+" else (t_size - 1 - t_pos)
            if rb != "-" and tb != "-":
                pairs.add((rc, tc))
            if rb != "-":
                r_pos += 1
            if tb != "-":
                t_pos += 1
    return pairs


def main():
    parser = argparse.ArgumentParser(description="Verify horizon_guide/horizon_merge MAF accuracy on chrM.maf")
    parser.add_argument("--maf", default=str(ROOT.parent / "pg" / "data" / "benchmarks" / "chrM.maf"))
    parser.add_argument("--species", nargs="+", default=["hg38", "panTro4", "rheMac3", "mm10", "canFam3"])
    args = parser.parse_args()

    maf_path = Path(args.maf)
    if not maf_path.exists():
        raise FileNotFoundError(f"Missing {maf_path}")

    cargo = os.environ.get("CARGO", "/home/agent/.cargo/bin/cargo")
    subprocess.run(
        [cargo, "build", "--release", "--bin", "horizon_guide", "--bin", "horizon_wga"],
        cwd=ROOT, check=True,
    )
    rel = ROOT / "target" / "release"

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        fastas, truth_pairs = extract_species_fastas_and_truth(maf_path, args.species, tmp_dir)
        specs = [f"{sp}={p}" for sp, p in fastas.items()]
        out_prefix = tmp_dir / "guide_out"
        res = subprocess.run(
            [str(rel / "horizon_guide"), "--progressive", "--hal", "--rbest", str(out_prefix), "0", *specs],
            capture_output=True, text=True, check=True,
        )
        print("horizon_guide:", res.stdout.strip())

        hal_res = subprocess.run(
            [str(rel / "horizon_wga"), "--hal-info", f"{out_prefix}.hal"],
            capture_output=True, text=True, check=True,
        )
        print("horizon_wga --hal-info:", " | ".join(hal_res.stdout.strip().splitlines()))

        pred_maf = Path(f"{out_prefix}.maf")
        ref_sp = args.species[0]
        for tgt_sp in args.species[1:]:
            if tgt_sp not in fastas:
                continue
            truth = truth_pairs[(ref_sp, tgt_sp)]
            pred = maf_pairs(pred_maf, ref_sp, tgt_sp)
            tp = len(truth & pred)
            rec = tp / len(truth) if truth else 0.0
            prec = tp / len(pred) if pred else 0.0
            print(f"  {ref_sp} x {tgt_sp}: truth={len(truth)} pred={len(pred)} tp={tp} recall={rec:.3f} prec={prec:.3f}")


if __name__ == "__main__":
    main()
