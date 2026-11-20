#!/usr/bin/env python3
"""Score minimap2 PAF (with CIGAR) against HPRC GRCh38.chr6 <-> CHM13.chr6 1-to-1 segments."""
import json
import re
import time
import numpy as np

TRUTH_TSV = "/tmp/hprc_chr6_truth_segments.tsv"
PAF_PATH = "/tmp/chr6_minimap2_asm5.paf"
OUT_NPZ = "/tmp/chr6_minimap2_data.npz"
OUT_JSON = "/tmp/chr6_minimap2_summary.json"

X_LEN = 170_805_979
Y_LEN = 110_724_939
BIN_BP = 1_000_000

t0 = time.time()
hprc_y_for_x = np.full(X_LEN, -1, dtype=np.int32)
hprc_cov_x = np.zeros(X_LEN, dtype=np.bool_)
hprc_cov_y = np.zeros(Y_LEN, dtype=np.bool_)

with open(TRUTH_TSV) as f:
    f.readline()
    for line in f:
        sid, xs, xe, ys, ye, strand, mult = line.rstrip("\n").split("\t")
        if mult != "1":
            continue
        xs, xe, ys, ye = int(xs), int(xe), int(ys), int(ye)
        hprc_cov_x[xs:xe] = True
        hprc_cov_y[ys:ye] = True
        if strand == "+":
            hprc_y_for_x[xs:xe] = np.arange(ys, ye, dtype=np.int32)
        else:
            hprc_y_for_x[xs:xe] = np.arange(ye - 1, ys - 1, -1, dtype=np.int32)

print(f"HPRC loaded in {time.time()-t0:.1f}s", flush=True)

mm_cov_x = np.zeros(X_LEN, dtype=np.bool_)
mm_cov_y = np.zeros(Y_LEN, dtype=np.bool_)
mm_exact_x = np.zeros(X_LEN, dtype=np.bool_)
mm_100bp_x = np.zeros(X_LEN, dtype=np.bool_)
mm_1kb_x = np.zeros(X_LEN, dtype=np.bool_)

total_pairs = 0
pairs_on_shared = 0
exact_pairs = 0
within_10bp = 0
within_100bp = 0
within_1kb = 0
n_records = 0

cigar_re = re.compile(r"(\d+)([MIDNSHP=X])")

t1 = time.time()
with open(PAF_PATH) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 12:
            continue
        # Only score primary alignments? Let's check tp:A:P vs all
        is_primary = any(tag == "tp:A:P" for tag in parts[12:])
        if not is_primary:
            continue
        n_records += 1
        qstart = int(parts[2])
        qend = int(parts[3])
        strand = parts[4]
        tstart = int(parts[7])
        cg = next((tag[5:] for tag in parts[12:] if tag.startswith("cg:Z:")), None)
        if cg is None:
            continue

        x_cur = tstart
        y_cur = qstart if strand == "+" else qend - 1
        y_step = 1 if strand == "+" else -1

        xs_list = []
        ys_list = []
        for num_s, op in cigar_re.findall(cg):
            n = int(num_s)
            if op in ("M", "=", "X"):
                xs_list.append(np.arange(x_cur, x_cur + n, dtype=np.int32))
                if y_step == 1:
                    ys_list.append(np.arange(y_cur, y_cur + n, dtype=np.int32))
                else:
                    ys_list.append(np.arange(y_cur, y_cur - n, -1, dtype=np.int32))
                x_cur += n
                y_cur += y_step * n
            elif op in ("I", "S"):
                y_cur += y_step * n
            elif op in ("D", "N"):
                x_cur += n

        if not xs_list:
            continue
        xs = np.concatenate(xs_list)
        ys = np.concatenate(ys_list)
        total_pairs += len(xs)
        mm_cov_x[xs] = True
        mm_cov_y[ys] = True

        ty = hprc_y_for_x[xs]
        mask = ty >= 0
        if np.any(mask):
            pairs_on_shared += int(np.count_nonzero(mask))
            diff = np.abs(ys[mask] - ty[mask])
            m0 = diff == 0
            m10 = diff <= 10
            m100 = diff <= 100
            m1k = diff <= 1000
            exact_pairs += int(np.count_nonzero(m0))
            within_10bp += int(np.count_nonzero(m10))
            within_100bp += int(np.count_nonzero(m100))
            within_1kb += int(np.count_nonzero(m1k))
            xm = xs[mask]
            mm_exact_x[xm[m0]] = True
            mm_100bp_x[xm[m100]] = True
            mm_1kb_x[xm[m1k]] = True

print(f"PAF (primary) scored in {time.time()-t1:.1f}s: records={n_records}, pairs={total_pairs:,}", flush=True)

hprc_1to1_x_bp = int(hprc_cov_x.sum())
hprc_1to1_y_bp = int(hprc_cov_y.sum())
n_xbins = (X_LEN + BIN_BP - 1) // BIN_BP
pad = n_xbins * BIN_BP - X_LEN
x_mm_bp = np.pad(mm_cov_x, (0, pad)).reshape(n_xbins, BIN_BP).sum(axis=1)
x_mm_exact_bp = np.pad(mm_exact_x, (0, pad)).reshape(n_xbins, BIN_BP).sum(axis=1)

summary = {
    "primary_records": n_records,
    "total_pairs": total_pairs,
    "unique_x_bp": int(mm_cov_x.sum()),
    "unique_y_bp": int(mm_cov_y.sum()),
    "pairs_on_hprc_1to1_x": pairs_on_shared,
    "exact_coord_pairs": exact_pairs,
    "within_10bp_pairs": within_10bp,
    "within_100bp_pairs": within_100bp,
    "within_1kb_pairs": within_1kb,
    "exact_coord_x_bp": int(mm_exact_x.sum()),
    "exact_coord_recall_on_hprc_1to1": int(mm_exact_x.sum()) / hprc_1to1_x_bp,
    "within_100bp_recall_on_hprc_1to1": int(mm_100bp_x.sum()) / hprc_1to1_x_bp,
    "within_1kb_recall_on_hprc_1to1": int(mm_1kb_x.sum()) / hprc_1to1_x_bp,
    "exact_coord_precision_on_shared": exact_pairs / pairs_on_shared if pairs_on_shared else 0.0,
    "within_10bp_precision_on_shared": within_10bp / pairs_on_shared if pairs_on_shared else 0.0,
    "within_100bp_precision_on_shared": within_100bp / pairs_on_shared if pairs_on_shared else 0.0,
    "within_1kb_precision_on_shared": within_1kb / pairs_on_shared if pairs_on_shared else 0.0,
}
with open(OUT_JSON, "w") as f:
    json.dump(summary, f, indent=2)
np.savez_compressed(OUT_NPZ, x_mm_bp=x_mm_bp, x_mm_exact_bp=x_mm_exact_bp)
print(json.dumps(summary, indent=2))
