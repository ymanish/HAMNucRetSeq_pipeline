
from pyfaidx import Fasta                      
from pathlib import Path
import pandas as pd
from src.config.path import DATA_DIR
import pybedtools as pbt
import numpy as np
import sys

ANNOT_DIR   = DATA_DIR/"annotations"

MACS3_PEAKS_DIR = DATA_DIR/ "macs3_mapq30_keepdup"
OUT_DIR     = DATA_DIR / "bound_regions"

if not OUT_DIR.exists():
    print(f"[+] Creating output directory: {OUT_DIR}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

hgfa        = ANNOT_DIR / "hg38.fa"
peaks_path       = MACS3_PEAKS_DIR / "pooled_peaks.narrowPeak"
out_fa = OUT_DIR / "pooled_peaks_bound.fa"



if not OUT_DIR.exists():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

if not hgfa.with_suffix(".fa.fai").exists():
    sys.exit("FASTA index .fai not found; run 'samtools faidx hg38.fa' first.")

valid_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

fasta = Fasta(hgfa, one_based_attributes=False)

peaks = pd.read_csv(peaks_path, 
                    sep="\t", 
                    header=None,
                    names=['chr','start','end','id','score','strand', 'fe','logP','logQ','rel_summit'])
valid_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

peaks = peaks[peaks.chr.isin(valid_chroms)]

with open(out_fa, "w") as out_fh:
    for row in peaks.itertuples(index=False):
        seq = fasta[row.chr][row.start:row.end].seq.upper()
        header = f"{row.chr}:{row.start}-{row.end}|{row.id}|{row.score}" 
        out_fh.write(f">{header}\n{seq}\n")

print(f"[+] FASTA written: {out_fa}")