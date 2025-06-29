from src.modules.NucFreeEnergy import NucleosomeBreath

import pandas as pd
from src.config.path import DATA_DIR, RESULTS_DIR
from src.config.custom_types import FreeEnergyResult, ProcessedSequence
from typing import List, Tuple, Iterator, Optional
import itertools
import tqdm 
from concurrent.futures import ProcessPoolExecutor, as_completed

from src.utils.fasta_checks import contains_non_canonical
from pathlib import Path
from Bio import SeqIO 
import csv

def calc_batch_energy(batch:List[ProcessedSequence],
                    method_nuc:str, *,  
                    free_dna_method:Optional[str]=None,
                      hard:bool=False, 
                      style:str="b_index", 
                      style_sites:Tuple[int, int]=(0,13))-> List[FreeEnergyResult]:

    nuc_breath = NucleosomeBreath(nuc_method=method_nuc, free_dna_method=free_dna_method)
    results:List[FreeEnergyResult]=[]

    if hard:
        print("Using hard method")
        compute_energy = lambda subseq, seq_id, sub_id: nuc_breath.calculate_free_energy_hard(seq147=subseq,
                                                                                            left=style_sites[0], 
                                                                                            right=style_sites[1], 
                                                                                            id=seq_id, 
                                                                                            subid=sub_id)
    else:
        print("Using soft method")
        compute_energy = lambda subseq, seq_id, sub_id: nuc_breath.calculate_free_energy_soft(seq601=subseq, 
                                                                                              left=style_sites[0], 
                                                                                            right=style_sites[1], 
                                                                                              id=seq_id, 
                                                                                              subid=sub_id,
                                                                                                style=style)
        
    for rec in batch:

        res = compute_energy(subseq = rec.sequence, seq_id = rec.id, sub_id = rec.subid)
        results.append(res)
    return results

def sliding_windows(seq: str, win: int, step: int) -> Iterator[tuple]:
    for off in range(0, len(seq) - win + 1, step):
        yield off, off + win, seq[off:off + win]

def process_fasta(path: Path, win: int = 147, step: int = 1) -> Iterator[ProcessedSequence]:
        for rec in SeqIO.parse(path, "fasta"):          
            header = rec.id                          
            s = str(rec.seq).upper()

            if contains_non_canonical(s):
                print(f"Skipping sequence {header} due to non-canonical characters.")
                continue


            sub_seq= 0
            for start, end, slice_ in sliding_windows(s, win, step):
                subid = str(sub_seq)
                sub_seq += 1
                yield ProcessedSequence(header, subid, slice_, start, end)

def main(fasta_path: Path, nuc_method:str, *, batch_size: int, 
         n_workers: int, outfile:Path = None, 
         freedna_method:Optional[str] = None)-> None:


    # generator of all windows
    windows = process_fasta(fasta_path, win=147, step=1)
    # windows = itertools.islice(process_fasta(fasta_path, win=147, step=1), 21)

    # chunk generator into batches of BATCH_SIZE
    def batcher(it, size):
        it = iter(it)
        for first in it:
            batch = list(itertools.chain([first], itertools.islice(it, size-1)))
            yield batch
    
    with outfile.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t") 
        
        with ProcessPoolExecutor(max_workers=n_workers) as pool:

            w.writerow(["id", "subid", "F", "F_entropy", "F_enthalpy", "F_freedna"])

            futures = [pool.submit(calc_batch_energy, 
                                                    batch,
                                                    method_nuc=nuc_method,
                                                    free_dna_method=freedna_method,
                                                    hard=False,
                                                    style="b_index", 
                                                    style_sites=(0,13)) for batch in batcher(windows, batch_size)
                        ]

            for fut in tqdm.tqdm(as_completed(futures),
                                total=len(futures),
                                desc="Processing batches"):
                for res in fut.result():
                    w.writerow([res.id, res.subid,
                                res.F, res.F_entropy,
                                res.F_enthalpy, res.F_freedna])
                    # print(f"id:{res.id} sub:{res.subid}  F:{res.F:.2f}")

    print(f"All done - results written to {outfile}")

    return None

def arg_parser():
    import argparse
    parser = argparse.ArgumentParser(description="Calculate nucleosome free energy from FASTA sequences.")
    parser.add_argument("--infile", type=Path, help="Path to the input FASTA file.")
    parser.add_argument("--batch_size", type=int, default=100, help="Number of sequences per batch.")
    parser.add_argument("--n_workers", type=int, default=21, help="Number of parallel workers.")
    parser.add_argument("--outfile", type=Path, help="Output file path for results.")
    parser.add_argument("--nuc_method", type=str, default="hybrid", help="Method for nucleosome free energy calculation.")
    parser.add_argument("--freedna_method", type=str, default=None, help="Method for free DNA energy calculation.")
    return parser.parse_args()

if __name__ == "__main__":
    import time 
    start = time.perf_counter()

    # FASTA = Path(DATA_DIR/"bound_regions/pooled_peaks_bound.fa")
    # OUTFILE = Path(RESULTS_DIR/"nucfe/bound_free_energy_results.txt")

    args = arg_parser()
    if args.infile:
        FASTA = args.infile
    if args.batch_size:
        batch_size = args.batch_size
    if args.n_workers:
        n_workers = args.n_workers
    if args.outfile:
        OUTFILE = args.outfile
    if args.nuc_method:
        NUC_METHOD = args.nuc_method
    if args.freedna_method:
        FREEDNA_METHOD = args.freedna_method
    if not FASTA.exists():
        raise FileNotFoundError(f"FASTA file {FASTA} does not exist. Please check the path.")
    
    if not OUTFILE.parent.exists():
        OUTFILE.parent.mkdir(parents=True, exist_ok=True)  


    print(f"Processing FASTA file: {FASTA}")
    print(f"Output will be written to: {OUTFILE}")

    # main(FASTA, batch_size=100, n_workers=21, outfile=Path(RESULTS_DIR/"nucfe/bound_free_energy_results.txt"))
    main(FASTA, batch_size=batch_size, nuc_method=NUC_METHOD, 
         n_workers=n_workers, outfile=Path(OUTFILE), 
         freedna_method=FREEDNA_METHOD)

    end = time.perf_counter()
    print(f"Total time taken: {end - start:.2f} seconds")
