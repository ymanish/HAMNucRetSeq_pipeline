

from __future__ import annotations
from pathlib import Path
import pandas as pd
import pybedtools as pbt
from pyfaidx import Fasta
from src.config.path import DATA_DIR
import matplotlib.pyplot as plt
import numpy as np
from Bio.Seq import Seq


import logging
import sys
from src.utils.logger_util import get_logger

from src.bedtsrc.nub.GetUniquePromo import parse_gtf_promoters, PromoterBedGenerator



# ------------------------------------------------------------------
def promoter_peak_stats(bound_bed: Path,
                        out_dir: Path, 
                        total_peaks: int) -> None:


    cols = ["p_chr","p_start","p_end","transcript_id", "TSS", "strand",
            "peak_chr","peak_start","peak_end", "peak_name","peak_score",
            "overlap_bp"]
    df = pbt.BedTool(str(bound_bed)).to_dataframe(names=cols)

    # overlap length ------------------------------------------------
    plt.figure(figsize=(6,4))
    plt.hist(df.overlap_bp, bins=40)
    plt.xlabel("peak-promoter overlap length (bp)")
    plt.ylabel("Count")
    plt.title("Distribution of peak-promoter overlap size")
    plt.tight_layout()
    plt.savefig(out_dir / "overlap_length_hist.png", dpi=300)

    # peaks per promoter -------------------------------------------
    peaks_per_prom = df.transcript_id.value_counts()
    plt.figure(figsize=(6,4))
    plt.hist(peaks_per_prom, bins=range(1, peaks_per_prom.max()+2))
    plt.xlabel("# peaks per promoter")
    plt.ylabel("Promoters")
    plt.title("How many peaks overlap each promoter")
    plt.tight_layout()
    plt.savefig(out_dir / "peaks_per_promoter_hist.png", dpi=300)

    # peak distribution around TSS --------------------------------------
    df['peak_mid'] = ((df.peak_start + df.peak_end) // 2)
    def signed_dist(r):
        return r.peak_mid - r.TSS if r.strand=='+' else r.TSS - r.peak_mid
    df['dist_TSS'] = df.apply(signed_dist, axis=1)

    plt.figure(figsize=(6,4))
    plt.hist(df.dist_TSS, bins=100, range=(-2100,2100))
    plt.xlabel("Peak midpoint - TSS (bp)  (upstream < 0, downstream > 0)")
    plt.ylabel("Count")
    plt.title("Peak distribution around TSS")
    plt.tight_layout()
    plt.savefig(out_dir / "peak_distance_hist.png", dpi=300)


    # numeric summary ----------------------------------------------
    summary = pd.DataFrame({"metric": ["n_pairs",
                                       "percent_of_totalpeaks_on_promoters",
                                        "median_overlap_bp","mean_overlap_bp",
                                        "median_peaks_per_prom","mean_peaks_per_prom"],
                            "value": [len(df),
                                    round((len(df) / total_peaks * 100), 2),
                                    df.overlap_bp.median(), df.overlap_bp.mean().round(2),
                                    peaks_per_prom.median().round(2), peaks_per_prom.mean().round(2)]
    })
    
    logger.info("[+] Numeric summary:")
    for _, row in summary.iterrows():
        logger.info(f"{row['metric']}: {row['value']}")

    logger.info(f"[+] statistics and PNG plots written to {out_dir}")
# ------------------------------------------------------------------


def drop_duplicates_bedtool(bedtool: pbt.BedTool) -> pbt.BedTool:

    colnames = ["chrom", "start", "end", "name", "score", "strand"]
    df = bedtool.to_dataframe(names=colnames + list(range(6, len(bedtool[0].fields))))
    df_clean = df[colnames].drop_duplicates()
    return pbt.BedTool.from_dataframe(df_clean)

def bed_to_fa(hg_fasta: Path, bedtool: pbt.BedTool, out_path: Path):
    fasta = Fasta(hg_fasta, one_based_attributes=False)
    with out_path.open("w") as out:
        for f in bedtool:
            header = f"{f.chrom}:{f.start}-{f.end}|{f.name}|{f.strand}|{f.score}"
            seq = fasta[f.chrom][f.start:f.end].seq.upper()
            if f.strand == '-':
                seq = str(Seq(seq).reverse_complement())
            out.write(f">{header}\n{seq}\n")



if __name__ == "__main__":

    ANNOT_DIR      = DATA_DIR / "annotations"
    PEAK_DIR       = DATA_DIR / "macs3_mapq30_keepdup"
    PROM_DIR       = DATA_DIR / "promoter_regions"

    if not PROM_DIR.exists():
        PROM_DIR.mkdir(parents=True, exist_ok=True)

    log_file = PROM_DIR / "promoter_regseq.log"
    logger = get_logger(__name__, log_file=log_file, level=logging.INFO)
    logger.info(f"[+] Log file configured at {log_file}")
    logger.info("Starting promoter_regseq analysis ...")

    HGFA   = ANNOT_DIR / "hg38.fa"
    CHROMS = ANNOT_DIR / "hg38.chrom.sizes"

    PEAKS  = PEAK_DIR  / "pooled_peaks.narrowPeak"

    PROM_BED        = PROM_DIR / "UNIQUE_promoter_regions.bed"
    PROM_BOUND_BED  = PROM_DIR / "bound_promoters.bed"
    PROM_UNBOUND_BED= PROM_DIR / "unbound_promoters.bed"

    PROM_BOUND_FA   = PROM_DIR / "promoter_bound.fa"
    PROM_UNBOUND_FA = PROM_DIR / "promoter_unbound.fa"


    # ---------------------------------------------------------------------------
    # read GTF -> UNIQUE promoter BED
    UP = 2000
    DOWN = 2000
    logger.info(f"[+] Parsing GTF to extract promoter regions (up={UP}, down={DOWN}) ...")
    pc_tx = parse_gtf_promoters(up=UP, down=DOWN)

    logger.info(f"[+] {len(pc_tx)} protein-coding transcripts extracted from GTF")

    PromoterBedGenerator(annot_dir=ANNOT_DIR, logger=logger).generate_unique_promoter_bed(
                                                                        pc_tx,
                                                                        method="canonical",
                                                                        out_bed=PROM_BED
                                                                    )


    # ---------------------------------------------------------------------------
    # intersect with peaks
    logger.info("[+] Intersecting with peaks ....")
    peak_bt   = pbt.BedTool(PEAKS).cut([0,1,2,3,4])
    TOTAL_PEAKS = len(peak_bt)
    logger.info(f"[+] Total peaks: {TOTAL_PEAKS}")

    prom_bt   = pbt.BedTool(PROM_BED)
    bound_bt  = prom_bt.intersect(peak_bt, wo=True)               # promoters with peaks
    unbound_bt= prom_bt.intersect(peak_bt, v=True)               # promoters without peaks

    bound_bt.saveas(PROM_BOUND_BED)
    unbound_bt.saveas(PROM_UNBOUND_BED)


    logger.info(f"[+] promoter peaks BED -> {PROM_BOUND_BED}")
    logger.info(f"[+] unbound promoters BED -> {PROM_UNBOUND_BED}")

    # ---------------------------------------------------------------------------
    # promoter peak stats
    logger.info("[+] Calculating promoter-peak stats ...")
    promoter_peak_stats(bound_bed=PROM_BOUND_BED,
                        out_dir=PROM_DIR, 
                        total_peaks=TOTAL_PEAKS)

    # ---------------------------------------------------------------------------
    # write FASTA files for bound and unbound promoters
    logger.info("[+] Writing FASTA files for bound and unbound promoters ...")
    unique_boundprom_bt = drop_duplicates_bedtool(bound_bt)

    bed_to_fa(hg_fasta=HGFA, bedtool=unique_boundprom_bt, out_path=PROM_BOUND_FA)
    bed_to_fa(hg_fasta=HGFA, bedtool=unbound_bt, out_path=PROM_UNBOUND_FA)
    logger.info(f"[+] FASTA written {PROM_BOUND_FA}")
    logger.info(f"[+] FASTA written {PROM_UNBOUND_FA}")
