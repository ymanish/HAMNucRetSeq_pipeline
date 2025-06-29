

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

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)

def bed_to_fa(hg_fasta: Path, bedtool: pbt.BedTool, out_path: Path):
    fasta = Fasta(hg_fasta, one_based_attributes=False)
    with out_path.open("w") as out:
        for f in bedtool:
            header = f"{f.chrom}:{f.start}-{f.end}|{f.name}|{f.strand}|{f.score}"
            seq = fasta[f.chrom][f.start:f.end].seq.upper()
            if f.strand == '-':
                seq = str(Seq(seq).reverse_complement())
            out.write(f">{header}\n{seq}\n")

def nonredundant_promoters(df: pd.DataFrame) -> pd.DataFrame:
    """Keep ≤1 transcript per nearly identical promoter ( ≥50 % overlap )."""
    keep = []
    intervals = df[['prom_start','prom_end','transcript_id']].sort_values('prom_start')
    for _, row in intervals.iterrows():
        s, e, tid = row.prom_start, row.prom_end, row.transcript_id
        if all(max(0, min(e, pe)-max(s, ps)) < 0.5*min(e-s, pe-ps)
               for ps, pe, _ in keep):
            keep.append((s, e, tid))
    keep_ids = {tid for _,_,tid in keep}
    return df[df.transcript_id.isin(keep_ids)]


def canonical_promoters(prom: pd.DataFrame,
                        mane_file: Path,
                        appris_file: Path) -> pd.DataFrame:
    """
    Return exactly ONE promoter per gene, ranked by:
      1) MANE Select / MANE Plus Clinical       (priority 1)
      2) APPRIS PRINCIPAL:1 to PRINCIPAL:5      (priority 2-6)
      3) longest CDS length                     (priority 7)

    prom dataframe must already contain: seqname, prom_start, prom_end, gene_id, transcript_id, cds_len
    """
    # -----------------------------------------------------------------
    mane = (pd.read_csv(mane_file, sep='\t', usecols=['Ensembl_Gene',
                                                                    'Ensembl_nuc',
                                                                    'MANE_status']).rename(columns={'Ensembl_Gene':'gene_id',
                                                                                                    'Ensembl_nuc':'transcript_id'}))
    mane['mane_pri'] = 1

    # -----------------------------------------------------------------
    # read APPRIS annotations (principal tags)
    appris = (pd.read_csv(appris_file, sep='\t',
                            usecols=['Gene ID',
                                    'Transcript ID',
                                    'APPRIS Annotation'])
                .rename(columns={'Gene ID':'gene_id',
                                    'Transcript ID':'transcript_id', 
                                    'APPRIS Annotation':'APPRIS_annotation'}))

    # keep only PRINCIPAL 
    appris = appris[appris.APPRIS_annotation.str.startswith('PRINCIPAL')]
    appris['appris_rank'] = (appris.APPRIS_annotation.str.extract(r'PRINCIPAL:(\d)').astype(int))
    appris['priority'] = appris['appris_rank'] + 1 
    # -----------------------------------------------------------------
    # merge priorities onto promoter table
    prom = prom.merge(
        mane[['gene_id','transcript_id','mane_pri']],
        how='left'
    ).merge(
        appris[['gene_id','transcript_id','priority']],
        how='left'
    )

    # choose MANE priority if present, else APPRIS, else 7
    prom['priority'] = prom['mane_pri'].fillna(prom['priority']).fillna(7)

    # -----------------------------------------------------------------
    # sort by gene, then priority, then CDS length  
    prom = (prom.sort_values(['gene_id','priority','cds_len'], ascending=[True,True,False]).groupby('gene_id', as_index=False).first())

    return prom


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

def calc_promoter_coords(r, up, down):
    # For positive strand: upstream = up, downstream = down
    if r.strand == '+':
        start = max(r.TSS - up - 1, 0)
        end = r.TSS + down
    else:  # negative strand: upstream is toward higher coordinates.
        start = max(r.TSS - down - 1, 0)
        end = r.TSS + up
    return pd.Series({'prom_start': start, 'prom_end': end})



def parse_gtf_promoters(gtf_file: Path, up: int = 2000, down: int = 2000) -> pd.DataFrame:

    gtf_cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    gtf = pd.read_csv(gtf_file, sep='\t', comment='#', names=gtf_cols)
    
    # Select protein-coding transcripts.
    pc_tx = gtf[(gtf.feature == 'transcript') & 
                gtf.attributes.str.contains('transcript_type "protein_coding"')].copy()
    # Extract gene and transcript IDs.
    pc_tx.loc[:, 'gene_id'] = pc_tx.attributes.str.extract(r'gene_id "([^"]+)"')
    pc_tx.loc[:, 'transcript_id'] = pc_tx.attributes.str.extract(r'transcript_id "([^"]+)"')
    pc_tx.loc[:, 'cds_len'] = (pc_tx.end - pc_tx.start)
    
    # Define TSS: for plus strand, TSS is start; for minus strand, TSS is end.
    pc_tx.loc[:, 'TSS'] = pc_tx.apply(lambda r: r.start if r.strand == '+' else r.end, axis=1)
    # Calculate promoter coordinates using a helper function.
    pc_tx[['prom_start', 'prom_end']] = pc_tx.apply(lambda r: calc_promoter_coords(r, up, down), axis=1)
    
    return pc_tx

def generate_unique_promoter_bed(pc_tx: pd.DataFrame,
                          trans_choice_method: str,
                          OUT_PROM_BED_FILE: Path) -> pd.DataFrame:

    # Filter promoters based on the provided method.
    if trans_choice_method == "canonical":
        logger.info("[+] Using canonical promoter selection ...")
        prom_tx_filtered = canonical_promoters(
                                pc_tx,
                                mane_file   = ANNOT_DIR / "MANE.GRCh38.v1.4.summary.txt",
                                appris_file = ANNOT_DIR / "appris_data.principal.txt"
                            )
    else:
        logger.info("[+] Using non-redundant promoter selection ...")
        prom_tx_filtered = (pc_tx.groupby('gene_id', group_keys=False)
                              .apply(nonredundant_promoters)
                              .reset_index(drop=True)
                              .drop_duplicates())

    logger.info(f"[+] After Filtering : {len(prom_tx_filtered)} protein-coding transcripts found in GTF")

    # Select key columns and filter for canonical chromosomes.
    prom = prom_tx_filtered[['seqname','prom_start','prom_end','gene_id', 'transcript_id','cds_len','TSS','strand']].copy()
    canonical_chr = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    prom = prom[prom.seqname.isin(canonical_chr)]
    
    # Write out the BED file
    prom[['seqname','prom_start','prom_end','transcript_id','TSS','strand']].to_csv(OUT_PROM_BED_FILE,
                                                                                        sep='\t',
                                                                                        header=False,
                                                                                        index=False)
    logger.info(f"[+] Promoter BED written to {OUT_PROM_BED_FILE}")

    return prom

def drop_duplicates_bedtool(bedtool: pbt.BedTool) -> pbt.BedTool:

    colnames = ["chrom", "start", "end", "name", "score", "strand"]
    df = bedtool.to_dataframe(names=colnames + list(range(6, len(bedtool[0].fields))))
    df_clean = df[colnames].drop_duplicates()
    return pbt.BedTool.from_dataframe(df_clean)
    
if __name__ == "__main__":
    # ---------------------------------------------------------------------------
    logger.info("Starting promoter_regseq analysis ...")

    ANNOT_DIR      = DATA_DIR / "annotations"
    PEAK_DIR       = DATA_DIR / "macs3_mapq30_keepdup"
    PROM_DIR       = DATA_DIR / "promoter_regions"

    log_file = PROM_DIR / "promoter_regseq.log"
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.info(f"[+] Log file configured at {log_file}")


    if not PROM_DIR.exists():
        logger.info(f"[+] Creating output directory: {PROM_DIR}")
        PROM_DIR.mkdir(parents=True, exist_ok=True)

    HGFA   = ANNOT_DIR / "hg38.fa"
    GTF    = ANNOT_DIR / "gencode.v47.annotation.gtf"
    CHROMS = ANNOT_DIR / "hg38.chrom.sizes"


    PEAKS  = PEAK_DIR  / "pooled_peaks.narrowPeak"

    PROM_BED        = PROM_DIR / "UNIQUE_promoter_regions.bed"
    PROM_BOUND_BED  = PROM_DIR / "bound_promoters.bed"
    PROM_UNBOUND_BED= PROM_DIR / "unbound_promoters.bed"

    PROM_BOUND_FA   = PROM_DIR / "promoter_bound.fa"
    PROM_UNBOUND_FA = PROM_DIR / "promoter_unbound.fa"



    # ---------------------------------------------------------------------------
    # read GTF -> promoter BED
    UP = 2000
    DOWN = 2000
    logger.info("[+] Parsing GTF to extract promoter regions (up={}, down={}) ...".format(UP, DOWN))

    pc_tx = parse_gtf_promoters(GTF, up=UP, down=DOWN)

    logger.info(f"[+] {len(pc_tx)} protein-coding transcripts extracted from GTF")

    unique_pc_tx = generate_unique_promoter_bed(pc_tx, 
                                                    trans_choice_method="canonical", 
                                                    OUT_PROM_BED_FILE=PROM_BED)


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
