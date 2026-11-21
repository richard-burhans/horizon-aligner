#!/usr/bin/env python3
"""Reproducible VGP/HPRC benchmark driver for `horizon_filter` and `horizon_locus`."""
import argparse
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark horizon_filter on chr6 FastGA PAF against 2D HPRC truth")
    parser.add_argument("--paf", type=Path, default=Path("/tmp/chr6_fastga.paf"))
    parser.add_argument("--query-fa", type=Path, default=Path("/tmp/hprc_grch38_chr6.fa"))
    parser.add_argument("--target-fa", type=Path, default=Path("/tmp/hprc_chm13_chr6.fa"))
    parser.add_argument("--truth", type=Path, default=Path("/tmp/hprc_chr6_truth_segments.tsv"))
    parser.add_argument("--locus-fa", type=Path, default=None, help="Optional multi-FASTA for horizon_locus")
    args = parser.parse_args()

    cargo = os.environ.get("CARGO", "/home/agent/.cargo/bin/cargo")
    subprocess.run(
        [cargo, "build", "--release", "--bin", "horizon_filter", "--bin", "horizon_locus"],
        cwd=ROOT,
        check=True,
    )
    rel = ROOT / "target" / "release"
    scorer = ROOT / "scripts" / "benchmarks" / "score_paf_vs_hprc.py"

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        if args.paf.exists():
            net_paf = tmp_dir / "chr6_net.paf"
            t0 = time.perf_counter()
            subprocess.run(
                [str(rel / "horizon_filter"), "--paf", str(args.paf), "--net", "--out", str(net_paf)],
                check=True,
            )
            print(f"horizon_filter --net: {time.perf_counter() - t0:.2f}s")
            if args.truth.exists():
                subprocess.run(
                    [sys.executable, str(scorer), "--paf", str(net_paf), "--truth", str(args.truth), "--all-records"],
                    check=True,
                )

            if args.query_fa.exists() and args.target_fa.exists():
                mdl_paf = tmp_dir / "chr6_mdl_net.paf"
                t0 = time.perf_counter()
                subprocess.run(
                    [
                        str(rel / "horizon_filter"),
                        "--paf", str(args.paf),
                        "--net",
                        "--query", str(args.query_fa),
                        "--target", str(args.target_fa),
                        "--out", str(mdl_paf),
                    ],
                    check=True,
                )
                print(f"horizon_filter --net (SDUST+MDL): {time.perf_counter() - t0:.2f}s")
                if args.truth.exists():
                    subprocess.run(
                        [sys.executable, str(scorer), "--paf", str(mdl_paf), "--truth", str(args.truth), "--all-records"],
                        check=True,
                    )

        if args.locus_fa and args.locus_fa.exists():
            t0 = time.perf_counter()
            subprocess.run(
                [str(rel / "horizon_locus"), "--fasta", str(args.locus_fa), "--out-prefix", str(tmp_dir / "locus")],
                check=True,
            )
            print(f"horizon_locus: {time.perf_counter() - t0:.2f}s")


if __name__ == "__main__":
    main()
