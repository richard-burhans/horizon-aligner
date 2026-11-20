#!/usr/bin/env python3
"""Head-to-head benchmark on a 10,000-segment HPRC chr22 subgraph: Horizon vs vg."""
import argparse
import os
import resource
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run_timed(cmd: list[str], cwd: Path | None = None, discard_stdout: bool = False) -> tuple[float, float, str]:
    t0 = time.perf_counter()
    res = subprocess.run(
        ["/usr/bin/env", "python3", "-c",
         "import resource, subprocess, sys; "
         "discard = sys.argv[1] == '1'; "
         "p = subprocess.run(sys.argv[2:], stdout=(subprocess.DEVNULL if discard else subprocess.PIPE), stderr=subprocess.PIPE, check=True); "
         "r = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss; "
         "print(f'MAXRSS_KB={r}'); "
         "sys.stdout.buffer.write(b'' if discard else p.stdout)",
         "1" if discard_stdout else "0",
         *cmd],
        cwd=cwd, capture_output=True, text=True, check=True,
    )
    dt = time.perf_counter() - t0
    rss_mb = 0.0
    lines = []
    for line in res.stdout.splitlines():
        if line.startswith("MAXRSS_KB="):
            rss_mb = int(line.split("=", 1)[1]) / 1024.0
        else:
            lines.append(line)
    return dt, rss_mb, "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Benchmark 10k-segment HPRC chr22 subgraph")
    parser.add_argument("--gfa", default="/tmp/hprc_chr22_sub_clean.gfa")
    parser.add_argument("--query", default="/tmp/hprc_chr22_sub_query.fa")
    parser.add_argument("--vg", default=str(ROOT.parent / "vg"))
    args = parser.parse_args()

    gfa = Path(args.gfa)
    query = Path(args.query)
    if not gfa.exists():
        raise FileNotFoundError(f"Missing GFA: {gfa}")

    cargo = os.environ.get("CARGO", "/home/agent/.cargo/bin/cargo")
    subprocess.run(
        [cargo, "build", "--release", "--bin", "horizon_gfa", "--bin", "horizon_clean",
         "--bin", "horizon_resolve", "--bin", "horizon_wga"],
        cwd=ROOT, check=True,
    )
    rel = ROOT / "target" / "release"

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        dt, rss, out = run_timed([str(rel / "horizon_gfa"), "stats", str(gfa)])
        print(f"horizon_gfa stats: {dt:.3f}s | {rss:.1f} MB RSS\n  {out.strip()}")

        clean_out = tmp_dir / "clean.gfa"
        dt, rss, out = run_timed([str(rel / "horizon_clean"), "--gfa", str(gfa), "--out", str(clean_out), "--min-cp-bits", "10.0"])
        print(f"horizon_clean:     {dt:.3f}s | {rss:.1f} MB RSS")

        res_out = tmp_dir / "resolved.gfa"
        dt, rss, out = run_timed([str(rel / "horizon_resolve"), str(gfa), "--out", str(res_out), "--collapse-dups", "--tag-walks"])
        print(f"horizon_resolve:   {dt:.3f}s | {rss:.1f} MB RSS")

        if query.exists():
            dt, rss, out = run_timed([str(rel / "horizon_wga"), "--graph", str(gfa), str(query), "0"])
            print(f"horizon_wga --graph: {dt:.3f}s | {rss:.1f} MB RSS\n  {out.strip().splitlines()[-1]}")

        vg_bin = Path(args.vg)
        if vg_bin.exists():
            dt, rss, out = run_timed([str(vg_bin), "stats", "-z", "-l", "-N", str(gfa)])
            print(f"vg stats:          {dt:.3f}s | {rss:.1f} MB RSS")
            dt, rss, _ = run_timed([str(vg_bin), "convert", "-g", "-p", str(gfa)], discard_stdout=True)
            print(f"vg convert -g -p:  {dt:.3f}s | {rss:.1f} MB RSS")
        else:
            print(f"vg not found at {vg_bin}; skipping vg baseline")


if __name__ == "__main__":
    main()
