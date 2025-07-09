from pathlib import Path
from typing import Union, List, Tuple, Optional, Dict
import numpy as np
import pandas as pd
import pybedtools as pbt
from pyfaidx import Fasta
from collections import defaultdict
import warnings
import sys


class GcLengthMatcher:
    """
    A class to build GC and length distribution-based sequences.
    
    The required inputs:
      - A human genome FASTA file.
      - A safe region file (as BED or BedTool) from which to sample windows.
      - Customizable bin sizes for sequence length and GC content.
      - Output file locations for the resulting BED and FASTA files.
    
    Example:
        >>> matcher = GcLengthMatcher(
                fasta_file="hg38.fa",
                safe_bed_file="safe_regions.bed",
                length_bin_size=50,
                gc_bin_factor=50,
                output_bed="matched.bed",
                output_fa="matched.fa"
            )
        >>> # peaks_bt must be provided as a BedTool object from peaks, e.g., from pooled_peaks.narrowPeak.
        >>> matcher.run(peaks_bt)
    """
    def __init__(self,
                 fasta_file: Union[str, Path],
                 safe_bed_file: Union[str, Path],
                 length_bin_size: int,
                 gc_bin_factor: int,
                 output_bed: Union[str, Path],
                 output_fa: Union[str, Path],
                 max_tries: int = 10_000_000,
                 seed:int = 42) -> None:
        
        self.rng = np.random.default_rng(seed)
        self.fasta: Fasta = self.load_fasta(fasta_file)

        self.safe_bt: pbt.BedTool = pbt.BedTool(str(safe_bed_file))
        self.safe_df: pd.DataFrame = self.prepare_safe_dataframe(self.safe_bt)
        self.length_bin_size: int = length_bin_size
        self.gc_bin_factor: int = gc_bin_factor
        self.output_bed: Union[str, Path] = output_bed
        self.output_fa: Union[str, Path] = output_fa
        self.max_tries: int = max_tries

    def load_fasta(self, hgfa_path: Union[str, Path]) -> Fasta:
        return Fasta(hgfa_path, one_based_attributes=False)

    def prepare_safe_dataframe(self, safe_bt: pbt.BedTool) -> pd.DataFrame:
        df = safe_bt.to_dataframe(names=["chr", "start", "end"])
        df['len'] = df['end'] - df['start']
        return df

    def gc_frac(self, seq: str) -> float:
        seq = seq.upper()
        return (seq.count('G') + seq.count('C')) / len(seq)

    def bin_key(self, L: int, gc: float) -> Tuple[int, int]:
        """
        Compute the bin key using custom bin sizes.
        
        :param L: Length.
        :param gc: GC fraction (should be between 0 and 1).
        :return: Tuple (length_bin, gc_bin) representing the computed bin.
        """
        if gc < 0 or gc > 1:
            raise ValueError(f"gc value must be between 0 and 1, got {gc}")
        return (L // self.length_bin_size, int(gc * self.gc_bin_factor))


    def build_target_histogram(self, peaks_bt: pbt.BedTool) -> Dict[Tuple[int, int], int]:
        """
        Build the target histogram based on length and GC content of bound peaks.
        
        :param peaks_bt: BedTool object of peaks.
        :return: Dictionary with keys as (length_bin, gc_bin) and counts as values.
        
        Example:
            Suppose peaks_bt has two intervals:
              - Interval 1 of length 100 and corresponding GC fraction 0.4.
              - Interval 2 of length 120 and GC fraction 0.4.
            With default parameters:
              length_bin = 100//50 = 2 and 120//50 = 2,
              gc_bin = int(0.4*50) = 20.
            The histogram will be: {(2, 20): 2}
        """
        df = peaks_bt.to_dataframe()
        lengths = df['end'] - df['start']
        gc_values = [self.gc_frac(self.fasta[r.chrom][r.start:r.end].seq) for r in peaks_bt]
        need = defaultdict(int)
        for L, G in zip(lengths, gc_values):
            key = self.bin_key(L, G)
            need[key] += 1
        return need

    def match_windows(self, need: Dict[Tuple[int, int], int]) -> List[pbt.Interval]:
        """
        Attempt to match windows from safe regions that satisfy the target histogram bins.
        
        :param need: Dictionary of target bins with counts.
        :return: List of BedTool Interval objects that match the bins.
        """
        matched_feats: List[pbt.Interval] = []
        tries = 0
        while any(v > 0 for v in need.values()) and tries < self.max_tries:
            remaining_bins = [k for k, v in need.items() if v > 0]
            key = self.rng.choice(remaining_bins)
            length_bin, gc_bin = key

            target_len = (length_bin * self.length_bin_size) + self.rng.integers(0, self.length_bin_size)
            tpl = self.random_window(target_len)
            if tpl is None:
                tries += 1
                continue
            chrom, start, end = tpl
            seq = self.fasta[chrom][start:end].seq
            G = self.gc_frac(seq)
            if self.bin_key(int(target_len), G) == tuple(key):
                feat = pbt.create_interval_from_list(
                    [chrom, str(start), str(end), f"{chrom}:{start}-{end}"]
                )
                matched_feats.append(feat)
                need[tuple(key)] -= 1
            tries += 1
        if any(v > 0 for v in need.values()):
            missing = sum(need.values())
            warnings.warn(f"[!] Warning: could not fill {missing} bins after {tries} draws")
        else:
            print(f"[+] All bins satisfied after {tries} draws")
        return matched_feats

    def random_window(self, target_len: int) -> Optional[Tuple[str, int, int]]:
        """
        Randomly pick a window from safe regions that can accommodate target_len.
        
        :param target_len: Target window length.
        :return: Tuple (chrom, start, end) if a candidate is found, else None.
        """
        candidates = self.safe_df[self.safe_df['len'] >= target_len]
        if candidates.empty:
            return None
        row = candidates.sample(n=1, random_state=self.rng).iloc[0]
        remaining = row['len'] - target_len
        offset = self.rng.integers(0, remaining) if remaining > 0 else 0
        return row['chr'], row['start'] + offset, row['start'] + offset + target_len

    def write_outputs(self, matched_feats: List[pbt.Interval]) -> None:
        """
        Write the matched intervals to BED and FASTA files.
        
        :param matched_feats: List of matched BedTool Interval objects.
        """
        pbt.BedTool(matched_feats).saveas(self.output_bed)
        with open(self.output_fa, "w") as fh:
            for f in matched_feats:
                fh.write(f">{f.chrom}:{f.start}-{f.end}|{f.name}\n{self.fasta[f.chrom][f.start:f.end].seq.upper()}\n")
        # print(f"[+] GC+length-matched sequence sampled from Genome: BED  ({self.output_bed})")
        # print(f"[+] GC+length-matched sequence sampled from Genome: FASTA ({self.output_fa})")

    def run(self, peaks_bt: pbt.BedTool) -> None:
        """
        Run the full pipeline: build target histogram from peaks,
        match windows from safe regions, and write output files.
        
        :param peaks_bt: BedTool object of peaks.
        """
        need = self.build_target_histogram(peaks_bt)
        matched_feats = self.match_windows(need)
        if not matched_feats:
            sys.exit("No matched features found. Check parameters or input files.")
        elif len(matched_feats) < len(peaks_bt):
            warnings.warn(f"[!] Only {len(matched_feats)} matched features found out of {len(peaks_bt)} peaks."
                           f" To get more matches, try to increase max_tries or (increase the bin size for length_bin_size and decrease the gc_bin_factor parameters).")
        else :
            print(f"[+] Successfully matched {len(matched_feats)} features.")
        
        self.write_outputs(matched_feats)