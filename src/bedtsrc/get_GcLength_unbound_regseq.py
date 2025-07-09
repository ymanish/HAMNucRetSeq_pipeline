import pandas as pd
import pybedtools as pbt
from pathlib import Path
import tempfile
import concurrent.futures
import numpy as np

from src.utils.bdtload_files import load_narrowpeaks
from src.bedtsrc.nub.GcLengthMatcher import GcLengthMatcher
from src.bedtsrc.nub.SafeRegions import get_safe_regions


def combine_batch_files(output_dir: Path, combined_bed: Path, combined_fa: Path) -> None:
   
    # Combine BED files.
    bed_files = sorted(output_dir.glob("matched_batch_*.bed"))
    with combined_bed.open("w") as out_bed:
        for bed in bed_files:
            with bed.open("r") as infile:
                out_bed.write(infile.read())
            bed.unlink()  # delete the batch file

    # Combine FASTA files.
    fa_files = sorted(output_dir.glob("matched_batch_*.fa"))
    with combined_fa.open("w") as out_fa:
        for fa in fa_files:
            with fa.open("r") as infile:
                out_fa.write(infile.read())
            fa.unlink() 


def sample_batch(batch_peaks_df: pd.DataFrame,
                 hgfa_file: Path,
                 safe_bed_file: Path,
                 length_bin_size: int,
                 gc_bin_factor: int,
                 max_tries: int,
                 output_bed: Path,
                 output_fa: Path,
                 seed: int) -> None:
    """
    Function to run sampling for a single batch.
    """
    batch_bed = pbt.BedTool.from_dataframe(batch_peaks_df[['chr', 'start', 'end']])
    matcher = GcLengthMatcher(
        fasta_file=hgfa_file,
        safe_bed_file=safe_bed_file,
        length_bin_size=length_bin_size,
        gc_bin_factor=gc_bin_factor,
        output_bed=output_bed,
        output_fa=output_fa,
        max_tries=max_tries,
        seed=seed
    )
    matcher.run(peaks_bt=batch_bed)


def parallel_sample(peaks_file: Path,
                    hgfa_file: Path,
                    gap_bed_file: Path,
                    rmsk_bed_file: Path,
                    chrom_sizes_file: Path,
                    output_dir: Path,
                    valid_chroms: list,
                    length_bin_size: int,
                    gc_bin_factor: int,
                    max_tries: int,
                    n_batches: int = 100,
                    seed: int = 42) -> None:
    """
    Split peaks into batches and run sampling in parallel.
    """
    # Load peaks.
    peaks_df = load_narrowpeaks(narrowpeaks_path=peaks_file, valid_chroms=valid_chroms)
    print(f"[+] Loaded {len(peaks_df)} peaks from {peaks_file}")

    
    # Create a temporary safe regions file.
    outbed = Path(tempfile.NamedTemporaryFile(suffix=".bed", delete=False).name)
    safe_bt = get_safe_regions(gap_bed=gap_bed_file,
                               rmsk_bed=rmsk_bed_file,
                               peaks_bt=pbt.BedTool.from_dataframe(peaks_df[['chr','start','end']]),
                               chrom_sizes=chrom_sizes_file,
                               valid_chroms=valid_chroms,
                               out_bed=outbed)
    total_bp = sum(f.end - f.start for f in safe_bt)
    print(f"[+] Safe regions span a total of {total_bp/1_000_000:.2f} million base pairs.")

    # Split peaks_df into n_batches roughly equally.
    batches = np.array_split(peaks_df, n_batches)
    
    # Prepare unique output files for each batch.
    tasks = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=21) as executor:
        for i, batch_df in enumerate(batches):
            batch_bed_out = output_dir / f"matched_batch_{i}.bed"
            batch_fa_out = output_dir / f"matched_batch_{i}.fa"
            tasks.append(executor.submit(sample_batch,
                                         batch_df,
                                         hgfa_file,
                                         outbed, # using the safe regions file
                                         length_bin_size,
                                         gc_bin_factor,
                                         max_tries,
                                         batch_bed_out,
                                         batch_fa_out,
                                         seed))
        # Wait for all tasks to complete.
        for future in concurrent.futures.as_completed(tasks):
            future.result() 

    print(f"[+] Sampling completed for all {n_batches} batches.")




if __name__ == "__main__":
    from src.config.path import DATA_DIR
    import sys
    
    ### Define input files.
    ANNOT_DIR = DATA_DIR / "annotations"
    chrom_sizes_file = ANNOT_DIR / "hg38.chrom.sizes"
    gap_bed_file = ANNOT_DIR / "gap_hg38.bed"
    rmsk_bed_file = ANNOT_DIR / "rmsk_hg38.bed"
    hgfa_file = ANNOT_DIR / "hg38.fa"
    
    # Check FASTA index.
    if not hgfa_file.with_suffix(".fa.fai").exists():
        sys.exit("FASTA index .fai not found; run 'samtools faidx hg38.fa' first.")
    
    # MACS3 peaks input.
    MACS3_PEAKS_DIR = DATA_DIR / "macs3_mapq30_keepdup"
    peaks_file = MACS3_PEAKS_DIR / "pooled_peaks.narrowPeak"
    valid_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    
    # Output directory.
    out_dir = DATA_DIR / "gclength_matched_unbound_regions"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    length_bin_size = 20
    gc_bin_factor = 50
    max_tries = 10_000_000
    
    print(f"[+] Using length bin size: {length_bin_size} bp")
    print(f"[+] Using GC bin factor: {gc_bin_factor} (bins of width 1/{gc_bin_factor} ≈ {1/gc_bin_factor:.2f})")
    print(f"[+] Using max tries: {max_tries}")
    print(f"[+] Valid chromosomes: {valid_chroms}")
    
    parallel_sample(peaks_file=peaks_file,
                    hgfa_file=hgfa_file,
                    gap_bed_file=gap_bed_file,
                    rmsk_bed_file=rmsk_bed_file,
                    chrom_sizes_file=chrom_sizes_file,
                    output_dir=out_dir,
                    valid_chroms=valid_chroms,
                    length_bin_size=length_bin_size,
                    gc_bin_factor=gc_bin_factor,
                    max_tries=max_tries,
                    n_batches=100)

    combined_bed_path = out_dir / "combined_matched.bed"
    combined_fa_path = out_dir / "combined_matched.fa"
    combine_batch_files(out_dir, combined_bed_path, combined_fa_path)
