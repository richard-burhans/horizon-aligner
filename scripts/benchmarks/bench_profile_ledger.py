#!/usr/bin/env python3
"""B3: Automated [profile] stage-timing and peak RSS regression ledger for horizon_wga."""
import argparse
import datetime
import json
import os
import re
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LEDGER = ROOT / "scripts" / "benchmarks" / "profile_ledger.jsonl"

PROFILE_RE = re.compile(
    r"\[profile\]\s+w12=([0-9.]+)s\s+\(fwd=(\d+),rc=(\d+),cov=([0-9.]+)%\)\s+"
    r"w10=([0-9.]+)s\s+lz1=([0-9.]+)s\s+seed_end=([0-9.]+)s\s+\(lcbs=(\d+)\+(\d+)\)\s+"
    r"dp\+xd=([0-9.]+)s\s+\(dp=(\d+),xd=(\d+)\)"
)
SUMMARY_RE = re.compile(
    r"aligned_pairs=(\d+)\s+unique_x_bp=(\d+)\s+unique_y_bp=(\d+)\s+total_bubbles=(\d+)\s+elapsed_s=([0-9.]+)"
)


def run_profiled(bin_path: Path, fa_x: Path, fa_y: Path, max_contigs: int = 1) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        maf_out = Path(tmp) / "out.maf"
        t0 = time.perf_counter()
        res = subprocess.run(
            [
                "/usr/bin/env",
                "python3",
                "-c",
                "import os, resource, subprocess, sys; "
                "env = dict(os.environ, HORIZON_PROFILE='1'); "
                "p = subprocess.run(sys.argv[1:], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True); "
                "r = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss; "
                "print(f'MAXRSS_KB={r}'); "
                "sys.stdout.buffer.write(p.stdout); sys.stderr.buffer.write(p.stderr)",
                str(bin_path),
                str(fa_x),
                str(fa_y),
                str(max_contigs),
                "--maf",
                str(maf_out),
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        wall_s = time.perf_counter() - t0
    combined = res.stdout + "\n" + res.stderr
    rss_kb = 0
    for line in combined.splitlines():
        if line.startswith("MAXRSS_KB="):
            rss_kb = int(line.split("=", 1)[1])
    prof_m = PROFILE_RE.search(combined)
    sum_m = SUMMARY_RE.search(combined)
    rec: dict[str, object] = {
        "wall_s": round(wall_s, 4),
        "max_rss_mb": round(rss_kb / 1024.0, 2),
    }
    if prof_m:
        rec["stages_s"] = {
            "w12": float(prof_m.group(1)),
            "w10": float(prof_m.group(5)),
            "lz1": float(prof_m.group(6)),
            "seed_end": float(prof_m.group(7)),
            "dp_xd": float(prof_m.group(10)),
        }
        rec["counts"] = {
            "fwd_lcbs": int(prof_m.group(2)),
            "rc_lcbs": int(prof_m.group(3)),
            "w12_cov_pct": float(prof_m.group(4)),
            "dp_jobs": int(prof_m.group(11)),
            "xdrop_jobs": int(prof_m.group(12)),
        }
    if sum_m:
        rec["aligned_pairs"] = int(sum_m.group(1))
        rec["unique_x_bp"] = int(sum_m.group(2))
        rec["unique_y_bp"] = int(sum_m.group(3))
        rec["total_bubbles"] = int(sum_m.group(4))
        rec["elapsed_s"] = float(sum_m.group(5))
    return rec


def main() -> None:
    parser = argparse.ArgumentParser(description="B3: [profile] stage-timing & RSS regression ledger")
    parser.add_argument("--record", action="store_true", help="Append results to profile_ledger.jsonl")
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    args = parser.parse_args()

    cargo = os.environ.get("CARGO", "/home/agent/.cargo/bin/cargo")
    subprocess.run([cargo, "build", "--release", "--bin", "horizon_wga"], cwd=ROOT, check=True)
    bin_path = ROOT / "target" / "release" / "horizon_wga"

    sha = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True).strip()
    ts = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")

    cases = [
        ("chr6_peri_5m", Path("/tmp/grch38_chr6_peri_5m.fa"), Path("/tmp/chm13_chr6_peri_5m.fa"), 1),
    ]
    rows = []
    for name, fx, fy, mc in cases:
        if not (fx.exists() and fy.exists()):
            print(f"Skipping {name}: missing {fx} or {fy}")
            continue
        metrics = run_profiled(bin_path, fx, fy, mc)
        row = {"timestamp": ts, "git_sha": sha, "benchmark": name, **metrics}
        rows.append(row)
        print(json.dumps(row))

    if args.record and rows:
        args.ledger.parent.mkdir(parents=True, exist_ok=True)
        with args.ledger.open("a") as f:
            for row in rows:
                f.write(json.dumps(row) + "\n")
        print(f"Appended {len(rows)} entries to {args.ledger}")


if __name__ == "__main__":
    main()
