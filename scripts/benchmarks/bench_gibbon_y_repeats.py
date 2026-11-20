#!/usr/bin/env python3
"""A4: Ampliconic repeat stress test on gibbon_y_tba.maf (5 gibbon Y species)."""
import argparse
import os
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_GIBBON_MAF = ROOT.parent / "gibbon_y_tba.maf"
SPECIES = [
    "mNomLeu1_NLE_Y_pri",
    "mSymSyn1_SSY_Y_pri",
    "mHooLeu1_HLE_Y_pri",
    "mHylPil1_HPI_Y_pri",
    "mHylMol1_HMO_Y_pri",
]


def parse_maf_blocks(path: Path, max_blocks: int = 0):
    cur = []
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith("s "):
                cur.append(line.split())
            elif line.strip() == "" and cur:
                yield cur
                cur = []
                n += 1
                if max_blocks > 0 and n >= max_blocks:
                    return
    if cur and (max_blocks == 0 or n < max_blocks):
        yield cur


def extract_gibbon_fastas(maf_path: Path, out_dir: Path, max_blocks: int) -> tuple[dict[str, Path], set[tuple[int, int]]]:
    seqs = {sp: [] for sp in SPECIES}
    offsets = {sp: 0 for sp in SPECIES}
    truth_pairs: set[tuple[int, int]] = set()
    ref_sp, tgt_sp = SPECIES[0], SPECIES[1]
    for blk in parse_maf_blocks(maf_path, max_blocks=max_blocks):
        by_sp = {}
        for row in blk:
            sp = row[1].split(".")[0]
            if sp in seqs and sp not in by_sp:
                by_sp[sp] = row[6].upper()
        if ref_sp in by_sp and tgt_sp in by_sp:
            rp = offsets[ref_sp]
            tp = offsets[tgt_sp]
            for rb, tb in zip(by_sp[ref_sp], by_sp[tgt_sp]):
                if rb != "-" and tb != "-":
                    truth_pairs.add((rp, tp))
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
        if chunks:
            p = out_dir / f"{sp}.fa"
            p.write_text(f">{sp}\n{''.join(chunks)}\n")
            paths[sp] = p
    return paths, truth_pairs


def maf_pairs(maf_path: Path, ref_sp: str, tgt_sp: str) -> set[tuple[int, int]]:
    pairs: set[tuple[int, int]] = set()
    for blk in parse_maf_blocks(maf_path):
        r_row = next((r for r in blk if r[1].split(".")[0] == ref_sp), None)
        t_row = next((r for r in blk if r[1].split(".")[0] == tgt_sp), None)
        if not r_row or not t_row:
            continue
        r_pos, r_strand, r_size = int(r_row[2]), r_row[4], int(r_row[5])
        t_pos, t_strand, t_size = int(t_row[2]), t_row[4], int(t_row[5])
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


def main() -> None:
    parser = argparse.ArgumentParser(description="A4: Gibbon Y ampliconic repeat stress test")
    parser.add_argument("--maf", type=Path, default=DEFAULT_GIBBON_MAF)
    parser.add_argument("--max-blocks", type=int, default=500, help="Number of TBA blocks to extract (0 = all)")
    args = parser.parse_args()

    cargo = os.environ.get("CARGO", "/home/agent/.cargo/bin/cargo")
    subprocess.run(
        [cargo, "build", "--release", "--bin", "horizon_wga", "--bin", "horizon_resolve"],
        cwd=ROOT,
        check=True,
    )
    rel = ROOT / "target" / "release"

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        fastas, truth = extract_gibbon_fastas(args.maf, tmp_dir, args.max_blocks)
        ref_fa, tgt_fa = fastas[SPECIES[0]], fastas[SPECIES[1]]
        for label, extra in [("default", []), ("one-to-one", ["--one-to-one"]), ("rbest+resolve", ["--rbest", "--resolve-dups"])]:
            maf_out = tmp_dir / f"gibbon_{label}.maf"
            gfa_out = tmp_dir / f"gibbon_{label}.gfa"
            res = subprocess.run(
                [
                    str(rel / "horizon_wga"),
                    str(ref_fa),
                    str(tgt_fa),
                    "1",
                    *extra,
                    "--maf",
                    str(maf_out),
                    "--gfa",
                    str(gfa_out),
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            pred = maf_pairs(maf_out, SPECIES[0], SPECIES[1])
            tp = len(truth & pred)
            rec = tp / len(truth) if truth else 0.0
            prec = tp / len(pred) if pred else 0.0
            print(f"[{label}] truth={len(truth):,} pred={len(pred):,} tp={tp:,} recall={rec:.4f} prec={prec:.4f}")
            print("  wga:", res.stdout.strip().splitlines()[0])

            subprocess.run(
                [
                    "python3",
                    str(ROOT / "scripts" / "benchmarks" / "validate_outputs.py"),
                    "--maf",
                    str(maf_out),
                    "--fastas",
                    str(ref_fa),
                    str(tgt_fa),
                    "--gfa",
                    str(gfa_out),
                ],
                check=True,
            )

            res_gfa = tmp_dir / f"gibbon_{label}_collapsed.gfa"
            r_out = subprocess.run(
                [str(rel / "horizon_resolve"), str(gfa_out), "--out", str(res_gfa), "--collapse-dups", "--tag-walks"],
                capture_output=True,
                text=True,
                check=True,
            )
            print("  resolve:", " | ".join(r_out.stdout.strip().splitlines()))


if __name__ == "__main__":
    main()
