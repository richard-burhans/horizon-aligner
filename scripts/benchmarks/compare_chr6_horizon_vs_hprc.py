#!/usr/bin/env python3
"""Compare Horizon GRCh38.chr6 vs CHM13.chr6 (q-arm) alignment against HPRC Minigraph-Cactus graph."""
import json
import numpy as np
import time

TRUTH_TSV = "/tmp/hprc_chr6_truth_segments.tsv"
MAF_PATH = "/tmp/chr6_fixed.maf"
OUT_NPZ = "/tmp/chr6_comparison_data.npz"
OUT_JSON = "/tmp/chr6_comparison_summary.json"

X_LEN = 170_805_979
Y_LEN = 110_724_939
Y_OFFSET = 61_401_689
BIN_BP = 1_000_000

t0 = time.time()
print("Loading HPRC truth segments...", flush=True)
hprc_y_for_x = np.full(X_LEN, -1, dtype=np.int32)
hprc_cov_x = np.zeros(X_LEN, dtype=np.bool_)
hprc_cov_y = np.zeros(Y_LEN, dtype=np.bool_)
hprc_multi_x = np.zeros(X_LEN, dtype=np.bool_)
hprc_multi_y = np.zeros(Y_LEN, dtype=np.bool_)

n_1to1 = 0
n_multi = 0
inv_1to1_bp = 0
truth_segs_sample = []  # (xs, xe, ys, ye, strand) for segments >= 500 bp for dotplot

with open(TRUTH_TSV) as f:
    header = f.readline()
    for line in f:
        sid, xs, xe, ys, ye, strand, mult = line.rstrip("\n").split("\t")
        xs, xe, ys, ye = int(xs), int(xe), int(ys), int(ye)
        seg_len = xe - xs
        if mult == "1":
            n_1to1 += 1
            hprc_cov_x[xs:xe] = True
            hprc_cov_y[ys:ye] = True
            if strand == "+":
                hprc_y_for_x[xs:xe] = np.arange(ys, ye, dtype=np.int32)
            else:
                hprc_y_for_x[xs:xe] = np.arange(ye - 1, ys - 1, -1, dtype=np.int32)
                inv_1to1_bp += seg_len
            if seg_len >= 500:
                truth_segs_sample.append((xs, xe, ys, ye, 1 if strand == "+" else -1))
        else:
            n_multi += 1
            hprc_multi_x[xs:xe] = True
            hprc_multi_y[ys:ye] = True

print(f"HPRC loaded in {time.time()-t0:.1f}s: 1to1_segs={n_1to1}, x_bp={hprc_cov_x.sum():,}, y_bp={hprc_cov_y.sum():,}, inv_bp={inv_1to1_bp:,}", flush=True)

t1 = time.time()
print("Streaming Horizon MAF...", flush=True)
hor_cov_x = np.zeros(X_LEN, dtype=np.bool_)
hor_cov_y = np.zeros(Y_LEN, dtype=np.bool_)
hor_match_exact_x = np.zeros(X_LEN, dtype=np.bool_)
hor_match_100bp_x = np.zeros(X_LEN, dtype=np.bool_)
hor_match_1kb_x = np.zeros(X_LEN, dtype=np.bool_)
hor_blocks_sample = []  # (xs, xe, ys, ye, strand) for blocks >= 500 bp

total_aligned_pairs = 0
pairs_on_hprc_x = 0
exact_pairs = 0
within_10bp = 0
within_100bp = 0
within_1kb = 0
n_blocks = 0

def parse_maf_blocks(path):
    cur = []
    with open(path) as f:
        for line in f:
            if line.startswith("s "):
                cur.append(line.split())
                if len(cur) == 2:
                    yield cur[0], cur[1]
                    cur = []
            elif line == "\n":
                cur = []

for s1, s2 in parse_maf_blocks(MAF_PATH):
    n_blocks += 1
    # s1 = GRCh38, s2 = CHM13
    x_start = int(s1[2])
    x_size = int(s1[3])
    x_strand = s1[4]
    x_src_size = int(s1[5])
    x_seq = s1[6]

    y_start = int(s2[2])
    y_size = int(s2[3])
    y_strand = s2[4]
    y_src_size = int(s2[5])
    y_seq = s2[6]

    # Convert MAF coordinates to forward-strand 0-based coordinates
    if x_strand == "+":
        xs_fwd = x_start
        xe_fwd = x_start + x_size
    else:
        xe_fwd = x_src_size - x_start
        xs_fwd = xe_fwd - x_size

    if y_strand == "+":
        ys_fwd = y_start
        ye_fwd = y_start + y_size
    else:
        ye_fwd = y_src_size - y_start
        ys_fwd = ye_fwd - y_size

    if x_size >= 500:
        hor_blocks_sample.append((xs_fwd, xe_fwd, ys_fwd, ye_fwd, 1 if x_strand == y_strand else -1))

    # Fast vectorized/C-speed coordinate extraction when ungapped or general
    if "-" not in x_seq and "-" not in y_seq:
        xs = np.arange(x_start, x_start + x_size, dtype=np.int32) if x_strand == "+" else np.arange(x_src_size - 1 - x_start, x_src_size - 1 - x_start - x_size, -1, dtype=np.int32)
        ys = np.arange(y_start, y_start + y_size, dtype=np.int32) if y_strand == "+" else np.arange(y_src_size - 1 - y_start, y_src_size - 1 - y_start - y_size, -1, dtype=np.int32)
    else:
        # Walk alignment columns
        xa = np.frombuffer(x_seq.encode("ascii"), dtype=np.uint8)
        ya = np.frombuffer(y_seq.encode("ascii"), dtype=np.uint8)
        both = (xa != 45) & (ya != 45)
        x_step = np.cumsum(xa != 45) - 1
        y_step = np.cumsum(ya != 45) - 1
        xi = x_step[both]
        yi = y_step[both]
        xs = (x_start + xi) if x_strand == "+" else (x_src_size - 1 - (x_start + xi))
        ys = (y_start + yi) if y_strand == "+" else (y_src_size - 1 - (y_start + yi))

    total_aligned_pairs += len(xs)
    hor_cov_x[xs] = True
    hor_cov_y[ys] = True

    ty = hprc_y_for_x[xs]
    mask = ty >= 0
    if np.any(mask):
        pairs_on_hprc_x += int(np.count_nonzero(mask))
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
        hor_match_exact_x[xm[m0]] = True
        hor_match_100bp_x[xm[m100]] = True
        hor_match_1kb_x[xm[m1k]] = True

print(f"MAF processed in {time.time()-t1:.1f}s: blocks={n_blocks:,}, pairs={total_aligned_pairs:,}", flush=True)

# 1-Mbp binned profiles along GRCh38.chr6 and CHM13.chr6 q-arm
n_xbins = (X_LEN + BIN_BP - 1) // BIN_BP
n_ybins = (Y_LEN + BIN_BP - 1) // BIN_BP

def bin_sums(arr, nbins):
    pad = nbins * BIN_BP - len(arr)
    if pad:
        arr = np.pad(arr, (0, pad))
    return arr.reshape(nbins, BIN_BP).sum(axis=1)

x_hprc_bp = bin_sums(hprc_cov_x, n_xbins)
x_hor_bp = bin_sums(hor_cov_x, n_xbins)
x_exact_bp = bin_sums(hor_match_exact_x, n_xbins)
x_100bp_bp = bin_sums(hor_match_100bp_x, n_xbins)
x_1kb_bp = bin_sums(hor_match_1kb_x, n_xbins)

y_hprc_bp = bin_sums(hprc_cov_y, n_ybins)
y_hor_bp = bin_sums(hor_cov_y, n_ybins)

# Base-level recall & precision of Horizon against HPRC 1-to-1 shared bases on GRCh38.chr6
hprc_1to1_x_bp = int(hprc_cov_x.sum())
hprc_1to1_y_bp = int(hprc_cov_y.sum())
hor_x_bp = int(hor_cov_x.sum())
hor_y_bp = int(hor_cov_y.sum())
hor_exact_x_bp = int(hor_match_exact_x.sum())
hor_100bp_x_bp = int(hor_match_100bp_x.sum())
hor_1kb_x_bp = int(hor_match_1kb_x.sum())
hor_and_hprc_x_bp = int((hor_cov_x & hprc_cov_x).sum())
hor_and_hprc_y_bp = int((hor_cov_y & hprc_cov_y).sum())

summary = {
    "x_len_bp": X_LEN,
    "y_len_bp": Y_LEN,
    "y_offset_bp": Y_OFFSET,
    "hprc_1to1_segs": n_1to1,
    "hprc_1to1_x_bp": hprc_1to1_x_bp,
    "hprc_1to1_y_bp": hprc_1to1_y_bp,
    "hprc_inv_1to1_bp": inv_1to1_bp,
    "horizon_blocks": n_blocks,
    "horizon_total_pairs": total_aligned_pairs,
    "horizon_unique_x_bp": hor_x_bp,
    "horizon_unique_y_bp": hor_y_bp,
    "horizon_and_hprc_x_bp": hor_and_hprc_x_bp,
    "horizon_and_hprc_y_bp": hor_and_hprc_y_bp,
    "pairs_on_hprc_1to1_x": pairs_on_hprc_x,
    "exact_coord_pairs": exact_pairs,
    "within_10bp_pairs": within_10bp,
    "within_100bp_pairs": within_100bp,
    "within_1kb_pairs": within_1kb,
    "exact_coord_x_bp": hor_exact_x_bp,
    "within_100bp_x_bp": hor_100bp_x_bp,
    "within_1kb_x_bp": hor_1kb_x_bp,
    "x_base_coverage_recall_of_hprc": hor_and_hprc_x_bp / hprc_1to1_x_bp if hprc_1to1_x_bp else 0.0,
    "y_base_coverage_recall_of_hprc": hor_and_hprc_y_bp / hprc_1to1_y_bp if hprc_1to1_y_bp else 0.0,
    "exact_coord_recall_on_hprc_1to1": hor_exact_x_bp / hprc_1to1_x_bp if hprc_1to1_x_bp else 0.0,
    "within_100bp_recall_on_hprc_1to1": hor_100bp_x_bp / hprc_1to1_x_bp if hprc_1to1_x_bp else 0.0,
    "within_1kb_recall_on_hprc_1to1": hor_1kb_x_bp / hprc_1to1_x_bp if hprc_1to1_x_bp else 0.0,
    "exact_coord_precision_on_shared": exact_pairs / pairs_on_hprc_x if pairs_on_hprc_x else 0.0,
    "within_100bp_precision_on_shared": within_100bp / pairs_on_hprc_x if pairs_on_hprc_x else 0.0,
    "within_1kb_precision_on_shared": within_1kb / pairs_on_hprc_x if pairs_on_hprc_x else 0.0,
}

with open(OUT_JSON, "w") as f:
    json.dump(summary, f, indent=2)

np.savez_compressed(
    OUT_NPZ,
    x_hprc_bp=x_hprc_bp,
    x_hor_bp=x_hor_bp,
    x_exact_bp=x_exact_bp,
    x_100bp_bp=x_100bp_bp,
    x_1kb_bp=x_1kb_bp,
    y_hprc_bp=y_hprc_bp,
    y_hor_bp=y_hor_bp,
    truth_segs=np.array(truth_segs_sample, dtype=np.int64),
    hor_blocks=np.array(hor_blocks_sample, dtype=np.int64),
)

print(json.dumps(summary, indent=2))
