

from pathlib import Path
import pandas as pd
import pybedtools as pbt
from pyfaidx import Fasta
from src.config.path import DATA_DIR
import matplotlib.pyplot as plt
import logging
from src.utils.logger_util import get_logger



def extract_minpoint_sequences(minima_dict: dict,
                               fasta_path: Path,
                               length: int = 147) -> dict:

    fa = Fasta(str(fasta_path))
    result = {}
    for seq_id, pos_list in minima_dict.items():
        if seq_id not in fa:
            print(f"Sequence ID {seq_id} not found in FASTA file {fasta_path}. Skipping.")
            continue
        rec = fa[seq_id]
        seqs = []
        for p in pos_list:
            start = p
            end = start + length
            subseq = rec[start:end].seq
            # print(f"Extracted subsequence for {seq_id} at position {p} and of length {len(subseq)}: {subseq}")
            seqs.append(subseq)
        result[seq_id] = seqs
    return result


def write_minpoints_fasta(minima_dict: dict,
                          seq_dict: dict,
                          output_path: Path):
    """
    Write out sequences in seq_dict to FASTA, header=">seq_id|pos:{p}"
    """
    with open(output_path, "w") as out_fa:
        for seq_id, seqs in seq_dict.items():
            positions = minima_dict.get(seq_id, [])
            for pos, seq in zip(positions, seqs):
                header = f">{seq_id}|pos:{pos}\n"
                out_fa.write(header)
                out_fa.write(seq + "\n")





if __name__ == "__main__":

    PROM_DIR       = DATA_DIR / "promoter_regions"

    log_file = None
    logger = get_logger(__name__, log_file=log_file, level=logging.INFO)
    logger.info(f"[+] Log file configured at {log_file}")
    logger.info("Starting promoter_regseq analysis ...")

    PROM_BOUND_FA   = PROM_DIR / "promoter_bound.fa"
    PROM_UNBOUND_FA = PROM_DIR / "promoter_unbound.fa"

    MINPOINT_PROM_DIR = DATA_DIR / "nbdata/promoters"



    # ---------------------------------------------------------------------------
    # read MINPOINT -> UNIQUE promoter BED


    import pickle

    # load unbound minima dict
    unbound_pkl = MINPOINT_PROM_DIR / "minpoints_unbound_dict.pkl"
    with open(unbound_pkl, "rb") as f:
        minima_unbound_dict = pickle.load(f)
    logger.info(f"Loaded unbound minima dict ({len(minima_unbound_dict)} ids) from {unbound_pkl}")

    # load bound minima dict
    bound_pkl = MINPOINT_PROM_DIR / "minpoints_bound_dict.pkl"
    with open(bound_pkl, "rb") as f:
        minima_bound_dict = pickle.load(f)
    logger.info(f"Loaded bound minima dict ({len(minima_bound_dict)} ids) from {bound_pkl}")

    # print (f"Unbound minima dict: {minima_unbound_dict}")
    # print (f"Bound minima dict: {minima_bound_dict}")


    # extract sequences
    bound_seqs   = extract_minpoint_sequences(minima_bound_dict, PROM_BOUND_FA)
    unbound_seqs = extract_minpoint_sequences(minima_unbound_dict, PROM_UNBOUND_FA)

    # write to FASTA
    out_bound = MINPOINT_PROM_DIR/"boundprom_minpoints.fa"
    out_unbound = MINPOINT_PROM_DIR/"unboundprom_minpoints.fa"
    write_minpoints_fasta(minima_bound_dict, bound_seqs, out_bound)
    write_minpoints_fasta(minima_unbound_dict, unbound_seqs, out_unbound)
    logger.info(f"Wrote bound minpoints FASTA to {out_bound}")
    logger.info(f"Wrote unbound minpoints FASTA to {out_unbound}")