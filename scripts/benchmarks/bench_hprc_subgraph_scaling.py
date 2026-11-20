#!/usr/bin/env python3
"""B2: HPRC subgraph scaling curve (10k -> 50k -> 200k segments) for Horizon graph tools."""
import argparse
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_HPRC_GFA_GZ = ROOT.parent / "pg" / "data" / "benchmarks" / "hprc-v1.0-mc-grch38.gfa.gz"


def run_timed(cmd: list[str], cwd: Path | None = None) -> tuple[float, float, str]:
    t0 = time.perf_counter()
    res = subprocess.run(
        [
            "/usr/bin/env",
            "python3",
            "-c",
            "import resource, subprocess, sys; "
            "p = subprocess.run(sys.argv[1:], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True); "
            "r = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss; "
            "print(f'MAXRSS_KB={r}'); "
            "sys.stdout.buffer.write(p.stdout)",
            *cmd,
        ],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=True,
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


def ensure_subgraphs(gfa_gz: Path, cache_dir: Path, max_segs: int = 200_000) -> tuple[dict[int, Path], Path]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    gfa_200k = cache_dir / f"hprc_sub_{max_segs}.gfa"
    query_fa = cache_dir / "hprc_sub_query_500bp.fa"

    if not gfa_200k.exists():
        print(f"Extracting {max_segs:,}-segment subgraph from {gfa_gz} (one-time)...", flush=True)
        decomp = shutil.which("pigz") or shutil.which("gzip")
        if decomp is None:
            raise RuntimeError("Neither pigz nor gzip found")
        proc = subprocess.Popen([decomp, "-dc", str(gfa_gz)], stdout=subprocess.PIPE, text=True, bufsize=1024 * 1024)
        assert proc.stdout is not None
        kept_ids: set[str] = set()
        seen_links = False
        with gfa_200k.open("w") as out:
            out.write("H\tVN:Z:1.1\n")
            for line in proc.stdout:
                tag = line[0]
                if tag == "S":
                    if len(kept_ids) < max_segs:
                        parts = line.split("\t", 3)
                        kept_ids.add(parts[1])
                        out.write(line)
                elif tag == "L":
                    seen_links = True
                    parts = line.split("\t", 5)
                    if parts[1] in kept_ids and parts[3] in kept_ids:
                        out.write(line)
                elif seen_links and tag in ("P", "W"):
                    break
        proc.terminate()
        proc.wait()

    out_paths: dict[int, Path] = {}
    for n in (10_000, 50_000, max_segs):
        p = cache_dir / f"hprc_sub_{n}.gfa"
        if not p.exists():
            kept: set[str] = set()
            with gfa_200k.open("r") as inp, p.open("w") as out:
                out.write("H\tVN:Z:1.1\n")
                for line in inp:
                    tag = line[0]
                    if tag == "S" and len(kept) < n:
                        seg_id = line.split("\t", 3)[1]
                        kept.add(seg_id)
                        out.write(line)
                    elif tag == "L":
                        parts = line.split("\t", 5)
                        if parts[1] in kept and parts[3] in kept:
                            out.write(line)
        out_paths[n] = p

    if not query_fa.exists():
        seq_parts: list[str] = []
        total_len = 0
        with out_paths[10_000].open("r") as inp:
            for line in inp:
                if line.startswith("S\t"):
                    seq = line.rstrip("\n").split("\t", 3)[2]
                    if seq != "*" and "N" not in seq.upper():
                        seq_parts.append(seq)
                        total_len += len(seq)
                        if total_len >= 500:
                            break
        query_seq = "".join(seq_parts)[:500]
        query_fa.write_text(f">hprc_sub_q{len(query_seq)}\n{query_seq}\n")

    return out_paths, query_fa


def main() -> None:
    parser = argparse.ArgumentParser(description="B2: HPRC subgraph scaling benchmark (10k -> 50k -> 200k)")
    parser.add_argument("--gfa-gz", type=Path, default=DEFAULT_HPRC_GFA_GZ)
    parser.add_argument("--cache-dir", type=Path, default=Path("/tmp/hprc_scaling_cache"))
    args = parser.parse_args()

    cargo = os.environ.get("CARGO", "/home/agent/.cargo/bin/cargo")
    subprocess.run(
        [
            cargo,
            "build",
            "--release",
            "--bin",
            "horizon_gfa",
            "--bin",
            "horizon_clean",
            "--bin",
            "horizon_resolve",
            "--bin",
            "horizon_wga",
        ],
        cwd=ROOT,
        check=True,
    )
    rel = ROOT / "target" / "release"
    subgraphs, query_fa = ensure_subgraphs(args.gfa_gz, args.cache_dir)

    print(f"\n{'Scale':>8} | {'Tool':<20} | {'Time (s)':>9} | {'Peak RSS (MB)':>13}")
    print("-" * 60)
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for n in sorted(subgraphs):
            gfa = subgraphs[n]
            dt, rss, _ = run_timed([str(rel / "horizon_gfa"), "stats", str(gfa)])
            print(f"{n:>8,} | {'horizon_gfa stats':<20} | {dt:>9.3f} | {rss:>13.1f}")

            clean_out = tmp_dir / f"clean_{n}.gfa"
            dt, rss, _ = run_timed(
                [str(rel / "horizon_clean"), "--gfa", str(gfa), "--out", str(clean_out), "--min-cp-bits", "10.0"]
            )
            print(f"{n:>8,} | {'horizon_clean':<20} | {dt:>9.3f} | {rss:>13.1f}")

            res_out = tmp_dir / f"res_{n}.gfa"
            dt, rss, _ = run_timed(
                [str(rel / "horizon_resolve"), str(gfa), "--out", str(res_out), "--collapse-dups", "--tag-walks"]
            )
            print(f"{n:>8,} | {'horizon_resolve':<20} | {dt:>9.3f} | {rss:>13.1f}")

            dt, rss, out = run_timed([str(rel / "horizon_wga"), "--graph", str(gfa), str(query_fa), "0"])
            last = out.strip().splitlines()[-1] if out.strip() else ""
            print(f"{n:>8,} | {'horizon_wga --graph':<20} | {dt:>9.3f} | {rss:>13.1f}  ({last})")


if __name__ == "__main__":
    main()
