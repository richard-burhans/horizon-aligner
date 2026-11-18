#!/usr/bin/env python3
import sys
import gzip
from pathlib import Path
from collections import defaultdict

def extract(maf_path: Path, ref_prefix: str, target_prefix: str, out_ref_fa: Path, out_tgt_fa: Path, out_pw_maf: Path):
    open_func = gzip.open if maf_path.suffix == '.gz' else open
    
    ref_seqs = defaultdict(list)
    tgt_seqs = defaultdict(list)
    
    with open_func(maf_path, 'rt') as f, open(out_pw_maf, 'w') as out_maf:
        out_maf.write("##maf version=1\n")
        
        current_block = []
        for line in f:
            if line.startswith('a'):
                if current_block:
                    # process block
                    ref_row = None
                    tgt_row = None
                    for row in current_block:
                        src = row[1]
                        if src.startswith(ref_prefix) and not ref_row:
                            ref_row = row
                        elif src.startswith(target_prefix) and not tgt_row:
                            tgt_row = row
                    
                    if ref_row and tgt_row:
                        
                        r_len = sum(len(s) for s in ref_seqs.get(ref_row[1], []))
                        t_len = sum(len(s) for s in tgt_seqs.get(tgt_row[1], []))
                        
                        # rewrite coordinates for the glued fasta
                        ref_row[2] = str(r_len)
                        ref_row[5] = "0" # dummy src len, we don't know total
                        
                        tgt_row[2] = str(t_len)
                        tgt_row[5] = "0"
                        
                        out_maf.write("a\n")
                        out_maf.write(" ".join(ref_row) + "\n")
                        out_maf.write(" ".join(tgt_row) + "\n\n")
                        
                        r_src = ref_row[1]
                        r_seq = ref_row[6].replace('-', '').upper()
                        ref_seqs[r_src].append(r_seq)
                        
                        t_src = tgt_row[1]
                        t_seq = tgt_row[6].replace('-', '').upper()
                        tgt_seqs[t_src].append(t_seq)
                        
                current_block = [line.strip().split()] if line.strip() != 'a' else []
            elif line.startswith('s '):
                current_block.append(line.strip().split())
                
        # don't forget the last block
        if current_block:
            ref_row = None
            tgt_row = None
            for row in current_block:
                src = row[1]
                if src.startswith(ref_prefix) and not ref_row:
                    ref_row = row
                elif src.startswith(target_prefix) and not tgt_row:
                    tgt_row = row
            
            if ref_row and tgt_row:
                out_maf.write("a\n")
                out_maf.write(" ".join(ref_row) + "\n")
                out_maf.write(" ".join(tgt_row) + "\n\n")
                
                r_src = ref_row[1]
                r_seq = ref_row[6].replace('-', '').upper()
                ref_seqs[r_src].append(r_seq)
                
                t_src = tgt_row[1]
                t_seq = tgt_row[6].replace('-', '').upper()
                tgt_seqs[t_src].append(t_seq)

    with open(out_ref_fa, 'w') as f:
        for src, seqs in ref_seqs.items():
            f.write(f">{src}\n{''.join(seqs)}\n")
            
    with open(out_tgt_fa, 'w') as f:
        for src, seqs in tgt_seqs.items():
            f.write(f">{src}\n{''.join(seqs)}\n")
            
    print(f"Extracted {len(ref_seqs)} ref contigs and {len(tgt_seqs)} target contigs.")

if __name__ == "__main__":
    maf = Path(sys.argv[1])
    ref_pre = sys.argv[2]
    tgt_pre = sys.argv[3]
    out_r = Path(sys.argv[4])
    out_t = Path(sys.argv[5])
    out_pw = Path(sys.argv[6])
    extract(maf, ref_pre, tgt_pre, out_r, out_t, out_pw)
