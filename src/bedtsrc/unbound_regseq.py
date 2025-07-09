"""
Create length-matched un-bound control intervals that exclude 
- assembly gaps, 
- RepeatMasker repeats,
- and the original peak set;
then extract their FASTA sequences.

Author: Manish 
"""

import sys, pathlib, tempfile
from src.config.path import DATA_DIR
from pathlib import Path
import pandas as pd
import pybedtools as pbt
from pyfaidx import Fasta



MACS3_PEAKS_DIR = DATA_DIR/ "macs3_mapq30_keepdup"
ANNOT_DIR   = DATA_DIR/"annotations"
chrom_sizes  = ANNOT_DIR / "hg38.chrom.sizes"

peaks_path       = MACS3_PEAKS_DIR / "pooled_peaks.narrowPeak"
hgfa        = ANNOT_DIR / "hg38.fa"
OUT_DIR     = DATA_DIR / "unbound_regions"
OUT_DIR.mkdir(parents=True, exist_ok=True)


gap_bed      = ANNOT_DIR / "gap_hg38.bed"
rmsk_bed     = ANNOT_DIR / "rmsk_hg38.bed"
excl_bed     = OUT_DIR / "combined_exclusion.bed"
safe_bed     = OUT_DIR / "genome_nonpeak_regions.bed"
unbound_bed  = OUT_DIR / "unbound_seq.bed"
unbound_fa   = OUT_DIR / "unbound_seq.fa"

tmp = tempfile.TemporaryDirectory()
tmpdir = pathlib.Path(tmp.name)

# -------- sanity checks ------------------------------------------------------

if not hgfa.with_suffix(".fa.fai").exists():
    sys.exit("FASTA index .fai not found; run 'samtools faidx hg38.fa' first.")


peaks = pd.read_csv(peaks_path, 
                    sep="\t", 
                    header=None,
                    names=['chr','start','end','id','score','strand', 'fe','logP','logQ','rel_summit'])
valid_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

peaks = peaks[peaks.chr.isin(valid_chroms)]
peaks_bed = pbt.BedTool.from_dataframe(peaks[['chr','start','end']])

allowed = {line.split()[0] for line in open(chrom_sizes)}
keep_primary = lambda f: f.chrom in allowed  
# --------combine exclusions ---------------------------------------
gap_bt   = pbt.BedTool(gap_bed ).filter(keep_primary).saveas()
rmsk_bt  = pbt.BedTool(rmsk_bed).filter(keep_primary).saveas()
peaks_bt = peaks_bed.filter(keep_primary).saveas()
excl_bt = (
        gap_bt
        # .cat(rmsk_bt, postmerge=False)
        .cat(peaks_bt, postmerge=False)
        .sort(g=str(chrom_sizes)) 
        .merge()
)
excl_bt.saveas(str(excl_bed))

# -------- complement (safe regions) --------------------------------
safe_bt = excl_bt.complement(g=str(chrom_sizes))       
safe_bt = safe_bt.filter(lambda f: f.chrom in valid_chroms).saveas(safe_bed)
print(f"[+] Safe regions BED written ({safe_bed})")

# -------- shuffle ---------------------------------------------------
shuffled_bt = peaks_bed.shuffle(
                g=str(chrom_sizes),
                incl=str(safe_bed),
                chrom=True,
                noOverlapping=True,
                seed=42)                               
shuffled_bt.saveas(unbound_bed)
print(f"[+] Un-bound BED written  ({unbound_bed})")

# -------- extract sequences --------------------------------------
fasta = Fasta(hgfa, one_based_attributes=False)  
with open(unbound_fa, "w") as out_fh:
    for feat in shuffled_bt:
        header  = f"{feat.chrom}:{feat.start}-{feat.end}"
        seq_str = fasta[feat.chrom][feat.start:feat.end].seq.upper()
        out_fh.write(f">{header}\n{seq_str}\n")

print(f"[+] Un-bound FASTA written ({unbound_fa})")
print(f"Un-bound BED : {unbound_bed}")
print(f"Un-bound FA  : {unbound_fa}")