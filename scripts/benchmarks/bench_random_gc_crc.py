#!/usr/bin/env python3
"""Multi-seed random-GC false-positive sweep for horizon_wga (18%, 30%, 45% GC)."""
import argparse
import os
import random
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def gen_random_dna(length: int, gc: float, rng: random.Random) -> str:
    p_gc = gc / 2.0
    p_at = (1.0 - gc) / 2.0
    bases = ["A", "T", "G", "C"]
    weights = [p_at, p_at, p_gc, p_gc]
    return "".join(rng.choices(bases, weights=weights, k=length))


def main():
    parser = argparse.ArgumentParser(description="Multi-seed random-GC FP sweep for horizon_wga")
    parser.add_argument("--length", type=int, default=10_000, help="Sequence length (bp)")
    parser.add_argument("--seeds", type=int, nargs="+", default=[999, 101, 202, 303, 404, 505])
    parser.add_argument("--gcs", type=float, nargs="+", default=[0.18, 0.30, 0.45])
    args = parser.parse_args()

    cargo = os.environ.get("CARGO", "/home/agent/.cargo/bin/cargo")
    subprocess.run([cargo, "build", "--release", "--bin", "horizon_wga"], cwd=ROOT, check=True)
    bin_path = ROOT / "target" / "release" / "horizon_wga"

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for gc in args.gcs:
            counts = []
            t0 = time.time()
            for seed in args.seeds:
                rng = random.Random(seed)
                sx = gen_random_dna(args.length, gc, rng)
                sy = gen_random_dna(args.length, gc, rng)
                fa_x = tmp_dir / f"x_{gc}_{seed}.fa"
                fa_y = tmp_dir / f"y_{gc}_{seed}.fa"
                maf_out = tmp_dir / f"out_{gc}_{seed}.maf"
                fa_x.write_text(f">rand_x\n{sx}\n")
                fa_y.write_text(f">rand_y\n{sy}\n")
                res = subprocess.run(
                    [str(bin_path), str(fa_x), str(fa_y), "1", "--maf", str(maf_out)],
                    capture_output=True, text=True, check=True,
                )
                pairs = 0
                for line in res.stdout.splitlines():
                    if "aligned_pairs=" in line:
                        for tok in line.split():
                            if tok.startswith("aligned_pairs="):
                                pairs = int(tok.split("=")[1])
                counts.append((seed, pairs))
            dt = time.time() - t0
            vals = [p for _, p in counts]
            print(
                f"GC={gc*100:.0f}% len={args.length}bp ({dt:.2f}s): "
                f"min={min(vals)} median={sorted(vals)[len(vals)//2]} max={max(vals)} "
                f"per_seed={dict(counts)}"
            )


if __name__ == "__main__":
    main()
