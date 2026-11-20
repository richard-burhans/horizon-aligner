#!/usr/bin/env python3
"""Extract exact GRCh38.chr6 <-> CHM13.chr6 shared-segment coordinates from HPRC GFA."""
import re
import subprocess
import time

GFA_GZ = "/home/rcb112/project/gemini/bench/pg/data/benchmarks/hprc-v1.0-mc-grch38.gfa.gz"
OUT_TSV = "/tmp/hprc_chr6_truth_segments.tsv"

t0 = time.time()
grch_steps = None  # list of (seg_id, orient)
chm13_steps = None
grch_info = None
chm13_info = None

def parse_p_steps(s: str):
    out = []
    for tok in s.split(","):
        if not tok:
            continue
        out.append((tok[:-1], tok[-1]))
    return out

def parse_w_steps(s: str):
    # >seg1<seg2...
    return [(m.group(2), "+" if m.group(1) == ">" else "-") for m in re.finditer(r"([><])([^><]+)", s)]

print("Pass 1: scanning P/W lines for GRCh38 chr6 and CHM13 chr6...", flush=True)
p = subprocess.Popen(["zcat", GFA_GZ], stdout=subprocess.PIPE, bufsize=16 * 1024 * 1024)
for raw in p.stdout:
    if raw.startswith(b"S\t") or raw.startswith(b"L\t"):
        continue
    if raw.startswith(b"P\t"):
        parts = raw.rstrip(b"\r\n").decode("ascii").split("\t")
        name = parts[1]
        if grch_steps is None and (name == "GRCh38.chr6" or name == "GRCh38#0#chr6"):
            grch_steps = parse_p_steps(parts[2])
            grch_info = f"P:{name}"
            print(f"  Found {grch_info}: {len(grch_steps)} steps", flush=True)
        elif chm13_steps is None and ("CHM13" in name and "chr6" in name):
            chm13_steps = parse_p_steps(parts[2])
            chm13_info = f"P:{name}"
            print(f"  Found {chm13_info}: {len(chm13_steps)} steps", flush=True)
    elif raw.startswith(b"W\t"):
        parts = raw.rstrip(b"\r\n").decode("ascii").split("\t")
        sample, hap, seq_id, start, end, walk = parts[1], parts[2], parts[3], parts[4], parts[5], parts[6]
        if grch_steps is None and sample == "GRCh38" and seq_id == "chr6":
            grch_steps = parse_w_steps(walk)
            grch_info = f"W:{sample}.{seq_id}:{start}-{end}"
            print(f"  Found {grch_info}: {len(grch_steps)} steps", flush=True)
        elif chm13_steps is None and sample == "CHM13" and seq_id == "chr6" and start == "61401689":
            chm13_steps = parse_w_steps(walk)
            chm13_info = f"W:{sample}.{seq_id}:{start}-{end}"
            print(f"  Found {chm13_info}: {len(chm13_steps)} steps", flush=True)
    if grch_steps is not None and chm13_steps is not None:
        break
p.stdout.close()
p.terminate()
p.wait()

assert grch_steps is not None and chm13_steps is not None, f"Missing paths: {grch_info}, {chm13_info}"
needed_segs = {sid for sid, _ in grch_steps} | {sid for sid, _ in chm13_steps}
print(f"Pass 1 done in {time.time()-t0:.1f}s: {len(needed_segs)} unique segments needed", flush=True)

t1 = time.time()
seg_lens = {}
print("Pass 2: scanning S lines for segment lengths...", flush=True)
p = subprocess.Popen(["zcat", GFA_GZ], stdout=subprocess.PIPE, bufsize=16 * 1024 * 1024)
for raw in p.stdout:
    if not raw.startswith(b"S\t"):
        continue
    # S\t<id>\t<seq>...
    tab1 = raw.find(b"\t", 2)
    sid = raw[2:tab1].decode("ascii")
    if sid in needed_segs:
        tab2 = raw.find(b"\t", tab1 + 1)
        seq_len = (tab2 if tab2 != -1 else len(raw.rstrip(b"\r\n"))) - (tab1 + 1)
        seg_lens[sid] = seq_len
        if len(seg_lens) == len(needed_segs):
            break
p.stdout.close()
p.terminate()
p.wait()
print(f"Pass 2 done in {time.time()-t1:.1f}s: got {len(seg_lens)}/{len(needed_segs)} lengths", flush=True)

# Build coordinate maps
grch_map = {}  # sid -> list of (start, end, orient)
pos = 0
for sid, orient in grch_steps:
    l = seg_lens[sid]
    grch_map.setdefault(sid, []).append((pos, pos + l, orient))
    pos += l
grch_total = pos

chm_map = {}
pos = 0
for sid, orient in chm13_steps:
    l = seg_lens[sid]
    chm_map.setdefault(sid, []).append((pos, pos + l, orient))
    pos += l
chm_total = pos

print(f"Total bp: GRCh38={grch_total}, CHM13={chm_total}", flush=True)

# Emit uniquely shared segments (1-to-1) and multi-copy shared segments
shared_1to1_bp = 0
shared_multi_bp = 0
n_1to1 = 0
n_multi = 0
with open(OUT_TSV, "w") as out:
    out.write("seg_id\tx_start\tx_end\ty_start\ty_end\tstrand\tmultiplicity\n")
    for sid, g_occs in grch_map.items():
        c_occs = chm_map.get(sid)
        if not c_occs:
            continue
        if len(g_occs) == 1 and len(c_occs) == 1:
            xs, xe, go = g_occs[0]
            ys, ye, co = c_occs[0]
            strand = "+" if go == co else "-"
            out.write(f"{sid}\t{xs}\t{xe}\t{ys}\t{ye}\t{strand}\t1\n")
            shared_1to1_bp += (xe - xs)
            n_1to1 += 1
        else:
            n_multi += 1
            shared_multi_bp += seg_lens[sid]
            for xs, xe, go in g_occs:
                for ys, ye, co in c_occs:
                    strand = "+" if go == co else "-"
                    out.write(f"{sid}\t{xs}\t{xe}\t{ys}\t{ye}\t{strand}\t{len(g_occs)}x{len(c_occs)}\n")

print(f"Wrote {OUT_TSV}: 1-to-1 segs={n_1to1} ({shared_1to1_bp} bp), multi-copy segs={n_multi} ({shared_multi_bp} bp)", flush=True)
