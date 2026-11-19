#!/usr/bin/env python3
import subprocess
import time
import sys
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ALIGNER_DIR = ROOT
BENCH_DIR = ROOT.parent / "pg" / "data" / "benchmarks"
FASTA_DIR = BENCH_DIR / "fastas"
EXTRACT_SCRIPT = ROOT / "scripts" / "benchmarks" / "extract_fasta_from_maf.py"

def run_cmd(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess:
    print(f"\n> {' '.join(cmd)}")
    return subprocess.run(cmd, cwd=cwd, check=True, text=True, capture_output=True)

def evaluate_pair(dataset_name: str, maf_gz: Path, ref_prefix: str, tgt_prefix: str, extra_args: list[str] = None):
    extra_args = extra_args or []
    print(f"\n=== Evaluating {dataset_name} ({ref_prefix} vs {tgt_prefix}) ===")
    
    # 1. Unzip the MAF if necessary
    unzipped_maf = maf_gz.with_suffix('')
    if maf_gz.suffix == '.gz' and not unzipped_maf.exists():
        print(f"Unzipping {maf_gz.name}...")
        run_cmd(["gunzip", "-k", str(maf_gz)])
    
    target_maf = unzipped_maf if maf_gz.suffix == '.gz' else maf_gz
    if not target_maf.exists():
        print(f"Missing {target_maf}, skipping.")
        return

    # 2. Extract specific pairs (include tgt_prefix in ref_fa so multiple targets don't overwrite)
    ref_fa = FASTA_DIR / f"{dataset_name}_{ref_prefix.replace('.', '_')}_vs_{tgt_prefix.replace('.', '_')}.fa"
    tgt_fa = FASTA_DIR / f"{dataset_name}_{tgt_prefix.replace('.', '_')}.fa"
    pw_maf = BENCH_DIR / f"{dataset_name}_{ref_prefix.replace('.', '_')}_{tgt_prefix.replace('.', '_')}.maf"
    
    print(f"Extracting {ref_prefix} and {tgt_prefix} from {target_maf.name}...")
    run_cmd([str(EXTRACT_SCRIPT), str(target_maf), ref_prefix, tgt_prefix, str(ref_fa), str(tgt_fa), str(pw_maf)])
    
    # 3. Run horizon_eval
    print(f"Running horizon_eval on {dataset_name}...")
    t0 = time.time()
    try:
        cmd = [
            "/home/agent/.cargo/bin/cargo", "run", "--release", "--bin", "horizon_eval", "--", 
            ref_prefix, tgt_prefix, str(ref_fa), str(tgt_fa), str(pw_maf)
        ] + extra_args
        res = run_cmd(cmd, cwd=ALIGNER_DIR)
        print(res.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running evaluator:\n{e.stderr}")
    t1 = time.time()
    print(f"Eval completed in {t1-t0:.2f}s")

def run_tuning(extra_args):
    print("\n--- RUNNING TUNING DATASETS ---")
    print("These datasets are used for calibrating algorithm thresholds.")
    # Dataset 3: multiz100way chrM (Human vs Mouse - Distant Vertebrate)
    evaluate_pair("multiz100way", BENCH_DIR / "chrM.maf.gz", "hg38", "mm10", extra_args)

    # Dataset 4: multiz100way chrM (Human vs Chimp - Very Close)
    evaluate_pair("multiz100way", BENCH_DIR / "chrM.maf.gz", "hg38", "panTro4", extra_args)

def run_testing(extra_args):
    print("\n--- RUNNING HELD-OUT TESTING DATASETS ---")
    print("These datasets MUST NOT be used for tuning parameters to prevent overfitting.")
    # Dataset 1: cactus447way (Human vs Horse - Distant Mammal)
    evaluate_pair("cactus447way", BENCH_DIR / "cactus447way_chr22_random.maf.gz", "hg38", "Equus_caballus", extra_args)

    # Dataset 2: cactus447way (Human vs Macaque - Close Primate)
    evaluate_pair("cactus447way", BENCH_DIR / "cactus447way_chr22_random.maf.gz", "hg38", "Macaca_mulatta", extra_args)

    # Dataset 5: multiz470way (Human vs Marmoset - Medium Primate)
    evaluate_pair("multiz470way", BENCH_DIR / "multiz470way_chr22_random.maf", "hg38", "HLcalPym1", extra_args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate Horizon Aligner")
    parser.add_argument("--mode", choices=["tune", "test", "all"], default="all", help="Which split to evaluate")
    parser.add_argument("--min-fallback-cov", type=float, help="Explicitly pass --min-fallback-cov to aligner")
    parser.add_argument("--hsp-threshold", type=int, help="Explicitly pass --hsp-threshold to aligner")
    args = parser.parse_args()
    
    print("Building horizon-aligner (release mode)...")
    run_cmd(["/home/agent/.cargo/bin/cargo", "build", "--release"], cwd=ALIGNER_DIR)
    
    extra_aligner_args = []
    if args.min_fallback_cov is not None:
        extra_aligner_args.extend(["--min-fallback-cov", str(args.min_fallback_cov)])
    if args.hsp_threshold is not None:
        extra_aligner_args.extend(["--hsp-threshold", str(args.hsp_threshold)])

    if args.mode in ["tune", "all"]:
        run_tuning(extra_aligner_args)
    if args.mode in ["test", "all"]:
        run_testing(extra_aligner_args)
    
    print("\nAll requested evaluations complete!")
