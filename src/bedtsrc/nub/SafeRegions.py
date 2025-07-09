import pandas as pd
import pybedtools as pbt
from pyfaidx import Fasta
from pathlib import Path
from typing import Union, List, Optional
import warnings


def get_safe_regions(
    gap_bed: Union[str, Path],
    rmsk_bed: Union[str, Path],
    peaks_bt: pbt.BedTool,
    chrom_sizes: Union[str, Path],
    valid_chroms: List[str], 
    out_bed: Optional[Union[str, Path]] = None
) -> pbt.BedTool:
    """
    Get safe regions of the genome (complement of gaps and peaks).
    
    :param gap_bed: Path to the gaps BED file.
    :param rmsk_bed: Path to the repeats/masked BED file.
    :param peaks_bt: BedTool object of peaks.
    :param chrom_sizes: Path to the chromosome sizes file.
    :param valid_chroms: List of valid chromosomes.
    :return: Safe regions as a BedTool object.
    """

    warnings.warn("Not using the rmsk_bed in the exclusion set, but it can be included if needed.")

    allowed = {line.split()[0] for line in open(chrom_sizes)}
    keep_primary = lambda f: f.chrom in allowed
    gap_bt_filtered = pbt.BedTool(gap_bed).filter(keep_primary).saveas()
    # rmsk_bt_filtered = pbt.BedTool(rmsk_bed).filter(keep_primary).saveas()
    peaks_bt_filtered = peaks_bt.filter(keep_primary).saveas()
    excl_bt = (
        gap_bt_filtered
        # Optionally include rmsk_bt_filtered by cat()
        .cat(peaks_bt_filtered, postmerge=False)
        .sort(g=str(chrom_sizes))
        .merge()
    )
    safe_bt = excl_bt.complement(g=str(chrom_sizes))
    safe_bt = safe_bt.filter(lambda f: f.chrom in valid_chroms).saveas(out_bed)
    return safe_bt