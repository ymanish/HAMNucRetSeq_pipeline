import pandas as pd
from pathlib import Path
from typing import Union, List


def load_narrowpeaks(narrowpeaks_path: Union[str, Path], valid_chroms: List[str]) -> pd.DataFrame:

    if not str(narrowpeaks_path).endswith('.narrowPeak'):
        raise ValueError(f"Expected narrowPeak file, got {narrowpeaks_path}")

    peaks = pd.read_csv(
        narrowpeaks_path, sep="\t", header=None,
        names=['chr','start','end','id','score','strand','fe','logP','logQ','rel_summit']
    )
    return peaks[peaks.chr.isin(valid_chroms)]