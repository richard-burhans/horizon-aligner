#!/usr/bin/env python3
"""C2: External WGA baselines (FastGA, minimap2, LASTZ) vs horizon_wga."""
import argparse
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run_timed(cmd: list[str], env: dict[str, str], stdout_path: Path | None = None) -> tuple[float, float]:
    t0 = time.perf_counter()
    res = subprocess.run(
        [
            "/usr/bin/env",
            "python3",
            "-c",
            "import resource, subprocess, sys; "
            "out_p = sys.argv[1]; "
            "f = open(out_p, 'wb') if out_p != '-' else subprocess.DEVNULL; "
            "p = subprocess.run(sys.argv[2:], stdout=f, stderr=subprocess.DEVNULL, check=True); "
            "r = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss; "
            "print(f'MAXRSS_KB={r}')",
            str(stdout_path) if stdout_path else "-",
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
    parser = argparse.ArgumentParser(description="C2: External WGA baselines vs horizon_wga")
    parser.add_argument("--fa-x", type=Path, default=Path("/tmp/grch38_chr6_peri_5m.fa"))
    parser.add_argument("--fa-y", type=Path, default=Path("/tmp/chm13_chr6_peri_5m.fa"))
    parser.add_argument("--ext-bin", type=Path, default=Path("/tmp/ext_tools/env/bin"))
    parser.add_argument("--skip-lastz", action="store_true", help="Skip LASTZ (recommended for >10 Mbp inputs)")
    args = parser.parse_args()

    cargo = os.environ.get("CARGO", "/home/agent/.cargo/bin/cargo")
    subprocess.run([cargo, "build", "--release", "--bin", "horizon_wga"], cwd=ROOT, check=True)
    rel = ROOT / "target" / "release"

    env = dict(os.environ)
    env["PATH"] = f"{args.ext_bin}:{env.get('PATH', '')}"
    fastga = shutil.which("FastGA", path=env["PATH"])
    minimap2 = shutil.which("minimap2", path=env["PATH"])
    lastz = shutil.which("lastz", path=env["PATH"])

    print(f"\nDataset: {args.fa_x.name} vs {args.fa_y.name}")
    print(f"{'Aligner':<24} | {'Time (s)':>9} | {'Peak RSS (MB)':>13}")
    print("-" * 54)

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        dt, rss = run_timed([str(rel / "horizon_wga"), str(args.fa_x), str(args.fa_y), "1", "--maf", str(tmp_dir / "h_def.maf")], env)
        print(f"{'horizon_wga (default)':<24} | {dt:>9.3f} | {rss:>13.1f}")

        dt, rss = run_timed([str(rel / "horizon_wga"), str(args.fa_x), str(args.fa_y), "1", "--one-to-one", "--maf", str(tmp_dir / "h_1to1.maf")], env)
        print(f"{'horizon_wga --one-to-one':<24} | {dt:>9.3f} | {rss:>13.1f}")

        if fastga:
            paf = tmp_dir / "fastga.paf"
            dt, rss = run_timed([fastga, "-pafx", "-T8", f"-P{tmp_dir}", str(args.fa_x), str(args.fa_y)], env, stdout_path=paf)
            print(f"{'FastGA (-pafx -T8)':<24} | {dt:>9.3f} | {rss:>13.1f}")

        if minimap2:
            paf = tmp_dir / "mm2.paf"
            dt, rss = run_timed([minimap2, "-x", "asm5", "-c", "-t", "8", str(args.fa_x), str(args.fa_y)], env, stdout_path=paf)
            print(f"{'minimap2 (-x asm5 -t8)':<24} | {dt:>9.3f} | {rss:>13.1f}")

        if lastz and not args.skip_lastz:
            maf = tmp_dir / "lastz.maf"
            dt, rss = run_timed([lastz, str(args.fa_x), str(args.fa_y), "--format=maf", "--ambiguous=iupac"], env, stdout_path=maf)
            print(f"{'LASTZ (--format=maf)':<24} | {dt:>9.3f} | {rss:>13.1f}")


if __name__ == "__main__":
    main()
