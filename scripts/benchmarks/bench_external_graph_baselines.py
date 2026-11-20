#!/usr/bin/env python3
"""C2: External graph tool baselines (gfatools, odgi, vg) vs Horizon on 10k->200k HPRC subgraphs."""
import argparse
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run_timed(cmd: list[str], env: dict[str, str] | None = None, discard_stdout: bool = True) -> tuple[float, float]:
    t0 = time.perf_counter()
    res = subprocess.run(
        [
            "/usr/bin/env",
            "python3",
            "-c",
            "import resource, subprocess, sys; "
            "discard = sys.argv[1] == '1'; "
            "p = subprocess.run(sys.argv[2:], stdout=(subprocess.DEVNULL if discard else subprocess.PIPE), stderr=subprocess.DEVNULL, check=True); "
            "r = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss; "
            "print(f'MAXRSS_KB={r}')",
            "1" if discard_stdout else "0",
            *cmd,
        ],
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    dt = time.perf_counter() - t0
    rss_mb = 0.0
    for line in res.stdout.splitlines():
        if line.startswith("MAXRSS_KB="):
            rss_mb = int(line.split("=", 1)[1]) / 1024.0
    return dt, rss_mb


def main() -> None:
    parser = argparse.ArgumentParser(description="C2: External graph baselines on HPRC subgraphs")
    parser.add_argument("--cache-dir", type=Path, default=Path("/tmp/hprc_scaling_cache"))
    parser.add_argument("--ext-bin", type=Path, default=Path("/tmp/ext_tools/env/bin"))
    parser.add_argument("--vg", type=Path, default=ROOT.parent / "vg")
    args = parser.parse_args()

    cargo = os.environ.get("CARGO", "/home/agent/.cargo/bin/cargo")
    subprocess.run(
        [cargo, "build", "--release", "--bin", "horizon_gfa", "--bin", "horizon_clean", "--bin", "horizon_resolve", "--bin", "horizon_wga"],
        cwd=ROOT,
        check=True,
    )
    rel = ROOT / "target" / "release"
    env = dict(os.environ)
    env["PATH"] = f"{args.ext_bin}:{env.get('PATH', '')}"

    gfatools = shutil.which("gfatools", path=env["PATH"])
    odgi = shutil.which("odgi", path=env["PATH"])
    vg = str(args.vg) if args.vg.exists() else shutil.which("vg", path=env["PATH"])
    query_fa = args.cache_dir / "hprc_sub_query_500bp.fa"

    print(f"\n{'Scale':>8} | {'Tool':<24} | {'Time (s)':>9} | {'Peak RSS (MB)':>13}")
    print("-" * 64)
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for n in (10_000, 50_000, 200_000):
            gfa = args.cache_dir / f"hprc_sub_{n}.gfa"
            if not gfa.exists():
                print(f"Skipping {n}: missing {gfa} (run bench_hprc_subgraph_scaling.py first)")
                continue

            for label, cmd in [
                ("horizon_gfa stats", [str(rel / "horizon_gfa"), "stats", str(gfa)]),
                ("horizon_clean", [str(rel / "horizon_clean"), "--gfa", str(gfa), "--out", str(tmp_dir / f"c_{n}.gfa"), "--min-cp-bits", "10.0"]),
                ("horizon_resolve", [str(rel / "horizon_resolve"), str(gfa), "--out", str(tmp_dir / f"r_{n}.gfa"), "--collapse-dups", "--tag-walks"]),
                ("horizon_wga --graph", [str(rel / "horizon_wga"), "--graph", str(gfa), str(query_fa), "0"]),
            ]:
                dt, rss = run_timed(cmd, env=env)
                print(f"{n:>8,} | {label:<24} | {dt:>9.3f} | {rss:>13.1f}")

            if gfatools:
                dt, rss = run_timed([gfatools, "stat", str(gfa)], env=env)
                print(f"{n:>8,} | {'gfatools stat':<24} | {dt:>9.3f} | {rss:>13.1f}")
            if odgi:
                og_out = tmp_dir / f"sub_{n}.og"
                dt_b, rss_b = run_timed([odgi, "build", "-g", str(gfa), "-o", str(og_out), "-t", "8"], env=env)
                dt_s, rss_s = run_timed([odgi, "stats", "-i", str(og_out), "-S"], env=env)
                print(f"{n:>8,} | {'odgi build':<24} | {dt_b:>9.3f} | {rss_b:>13.1f}")
                print(f"{n:>8,} | {'odgi stats -S':<24} | {dt_s:>9.3f} | {rss_s:>13.1f}")
            if vg:
                dt, rss = run_timed([vg, "stats", "-z", "-l", "-N", str(gfa)], env=env)
                print(f"{n:>8,} | {'vg stats':<24} | {dt:>9.3f} | {rss:>13.1f}")
                dt, rss = run_timed([vg, "convert", "-g", "-p", str(gfa)], env=env)
                print(f"{n:>8,} | {'vg convert -g -p':<24} | {dt:>9.3f} | {rss:>13.1f}")


if __name__ == "__main__":
    main()
