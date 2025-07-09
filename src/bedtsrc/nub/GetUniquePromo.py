
from pathlib import Path
import pandas as pd
import pybedtools as pbt
from typing import Literal
from src.utils.logger_util import get_logger
from src.config.path import DATA_DIR


class PromoterBedGenerator:
    """
    Encapsulate logic for filtering and writing unique promoter BEDs.
    """

    def __init__(self,
                 annot_dir: Path,
                 logger=None):
        self.ANNOT_DIR = annot_dir
        self.logger = logger or get_logger(__name__, log_file=None, level="INFO")

    def nonredundant_promoters(self, df: pd.DataFrame) -> pd.DataFrame:
        """Keep <=1 transcript per nearly identical promoter ( >=50 % overlap )."""
        keep = []
        intervals = df[['prom_start','prom_end','transcript_id']].sort_values('prom_start')
        for _, row in intervals.iterrows():
            s, e, tid = row.prom_start, row.prom_end, row.transcript_id
            if all(max(0, min(e, pe)-max(s, ps)) < 0.5*min(e-s, pe-ps)
                   for ps, pe, _ in keep):
                keep.append((s, e, tid))
        keep_ids = {tid for _,_,tid in keep}
        return df[df.transcript_id.isin(keep_ids)]

    def canonical_promoters(self,
                            prom: pd.DataFrame,
                            mane_file: Path,
                            appris_file: Path) -> pd.DataFrame:
        """
        Return exactly ONE promoter per gene, ranked by:
        1) MANE Select / MANE Plus Clinical       (priority 1)
        2) APPRIS PRINCIPAL:1 to PRINCIPAL:5      (priority 2-6)
        3) longest CDS length                     (priority 7)

        prom dataframe must already contain: seqname, prom_start, prom_end, gene_id, transcript_id, cds_len
        """
        mane = (pd.read_csv(mane_file, sep='\t',
                             usecols=['Ensembl_Gene','Ensembl_nuc','MANE_status'])
                    .rename(columns={'Ensembl_Gene':'gene_id',
                                     'Ensembl_nuc':'transcript_id'}))
        mane['mane_pri'] = 1

        appris = (pd.read_csv(appris_file, sep='\t',
                              usecols=['Gene ID','Transcript ID','APPRIS Annotation'])
                      .rename(columns={'Gene ID':'gene_id',
                                       'Transcript ID':'transcript_id',
                                       'APPRIS Annotation':'APPRIS_annotation'}))
        appris = appris[appris.APPRIS_annotation.str.startswith('PRINCIPAL')]
        appris['appris_rank'] = (appris.APPRIS_annotation
                                    .str.extract(r'PRINCIPAL:(\d)')
                                    .astype(int))
        appris['priority'] = appris['appris_rank'] + 1

        merged = (prom.merge(mane[['gene_id','transcript_id','mane_pri']], how='left')
                      .merge(appris[['gene_id','transcript_id','priority']], how='left'))
        merged['priority'] = merged['mane_pri'].fillna(merged['priority']).fillna(7)

        out = (merged.sort_values(['gene_id','priority','cds_len'],
                                  ascending=[True,True,False])
                      .groupby('gene_id', as_index=False)
                      .first())
        return out

    def generate_unique_promoter_bed(self,
                                    pc_tx: pd.DataFrame,
                                    method: Literal['canonical','nonredundant'],
                                    out_bed: Path) -> pd.DataFrame:
        """
        Filter pc_tx by either 'canonical' or 'nonredundant' rules,
        write BED to out_bed, and return the filtered DataFrame.
        """
        if method == "canonical":
            self.logger.info("[+] Using canonical promoter selection …")
            prom_tx = self.canonical_promoters(
                pc_tx,
                mane_file   = self.ANNOT_DIR/"MANE.GRCh38.v1.4.summary.txt",
                appris_file = self.ANNOT_DIR/"appris_data.principal.txt"
            )
        else:
            self.logger.info("[+] Using non-redundant promoter selection …")
            prom_tx = (pc_tx.groupby('gene_id', group_keys=False)
                          .apply(self.nonredundant_promoters)
                          .reset_index(drop=True)
                          .drop_duplicates())

        self.logger.info(f"[+] After filtering: {len(prom_tx)} transcripts")

        # restrict to canonical chromosomes
        canonical = {f"chr{i}" for i in range(1,23)} | {"chrX","chrY"}
        prom_tx = prom_tx[prom_tx.seqname.isin(canonical)]

        # write BED
        cols = ['seqname','prom_start','prom_end','transcript_id','TSS','strand']
        prom_tx[cols].to_csv(out_bed, sep='\t', header=False, index=False)
        self.logger.info(f"[+] Promoter BED written to {out_bed}")
        return prom_tx




def calc_promoter_coords(r, up, down):
    # For positive strand: upstream = up, downstream = down
    if r.strand == '+':
        start = max(r.TSS - up - 1, 0)
        end = r.TSS + down
    else:  # negative strand: upstream is toward higher coordinates.
        start = max(r.TSS - down - 1, 0)
        end = r.TSS + up
    return pd.Series({'prom_start': start, 'prom_end': end})

def parse_gtf_promoters(up: int = 2000, down: int = 2000, 
                        gtf_file: Path=DATA_DIR / "annotations/gencode.v47.annotation.gtf" 
                        ) -> pd.DataFrame:

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



if __name__ == "__main__":

    from src.config.path import DATA_DIR

    parent_logger = get_logger(__name__, log_file=None, level="INFO")


    GTF    = DATA_DIR / "annotations/gencode.v47.annotation.gtf"
    UP     = 2000
    DOWN   = 2000
    pc_tx = parse_gtf_promoters(up=UP, down=DOWN, gtf_file=GTF)

    prombedgen = PromoterBedGenerator(annot_dir=Path(DATA_DIR/"annotations"), logger=parent_logger)
    print ("Generating unique promoter BEDs...")
    prombedgen.generate_unique_promoter_bed(
                                                            pc_tx = pc_tx,
                                                            method = "canonical", 
                                                            out_bed = "example_promoter.bed"
                                                        )

    # parent_logger.info(f"Generated unique promoter BED with {len(prom_bed)} entries.")