# src/core/compute_weighted_fe.py
# Created on 2025-04-06

import numpy as np
from Bio import SeqIO
from src.config import path
from src.config.custom_types import FreeEnergyResultSimple 
from src.config.var import MAX_WORKERS, BATCH_SIZE, PARAM_TYPE
from src.modules.MarkovNucFe import calc_markov_fe
from src.utils.fasta_checks import contains_non_canonical
from typing import List
from concurrent.futures import ProcessPoolExecutor
import sys


def F_weighted_markov_model(id:str, seq:str)->FreeEnergyResultSimple:
    """
    This function calculates the weighted free energy for a given sequence.
    It uses the Markov model to compute the free energy based on the sequence.
    """
    # Placeholder for actual implementation
    P,ID = calc_markov_fe(id, 
                   seq,
                    order=3, 
                    file_l=path.PARAM_DIR / f"markov_prob/{PARAM_TYPE}_RoomTemp.trinucdist", 
                    file_s=path.PARAM_DIR / f"markov_prob/{PARAM_TYPE}_RoomTemp.dinucdist")

    # Calculate the weighted free energy
    P_min = np.min(P)
    P_rel = P/P_min
    Z_s = np.sum(P_rel)
    P_log = np.log(P)
    P_log_row = P_log.reshape(-1, 1)  # Convert to row vector
    P_rel_row = P_rel.reshape(-1, 1)  
    bF_weighted = -(P_log_row.T @ P_rel_row)[0, 0]/Z_s ### F_weighted*beta and beta is 1


    # print(P_log_col.shape)
    # print(P_rel_col.shape)
    return FreeEnergyResultSimple(bF_weighted, ID)


def batch_calculate_free_energy(seq_records:List[SeqIO.SeqRecord]) -> List[FreeEnergyResultSimple]:
    
    results: List[FreeEnergyResultSimple] = []

    for record in seq_records:

        if contains_non_canonical(str(record.seq)):
            print(f"Skipping sequence {record.id} due to non-canonical characters.")
            results.append(FreeEnergyResultSimple(-99999.0, record.id))

        elif len(record.seq) > 10000:
            print(f"Skipping sequence {record.id} because the length > 10k.")
            results.append(FreeEnergyResultSimple(-99999.0, record.id))

        else:
            # Calculate free energy using the Markov model
            F_weight_seq = F_weighted_markov_model(record.id, str(record.seq).upper())
            results.append(F_weight_seq)
            # print(f"Free energy for sequence {F_weight_seq.id}: {F_weight_seq.energy}")


    return results


def run_batches(seq_records:List[SeqIO.SeqRecord], batch_size:int = 100, max_physical_core=11) -> List[FreeEnergyResultSimple]:
    """
    This function divides the sequences into batches and processes each batch in parallel.
    """
    # Divide the sequences into batches
    batches = [seq_records[i:i + batch_size] for i in range(0, len(seq_records), batch_size)]
    print(f"Total batches: {len(batches)}")

    # Use ProcessPoolExecutor for parallel processing
    all_results = []
    with ProcessPoolExecutor(max_workers=max_physical_core) as executor:
        for batch_result in executor.map(batch_calculate_free_energy, batches):
            all_results.extend(batch_result)

    return all_results


if __name__ == '__main__': 
    import time 
    start = time.perf_counter()

    max_workers = int(sys.argv[1]) if len(sys.argv) > 1 else MAX_WORKERS
    batch_size = int(sys.argv[2]) if len(sys.argv) > 2 else BATCH_SIZE

    # Load the sequences from the FASTA file
    fasta_file = path.DATA_DIR / "bound_regions/pooled_peaks_bound.fa"


    seq_records = list(SeqIO.parse(fasta_file, "fasta"))
    all_results = run_batches(seq_records, batch_size=batch_size, max_physical_core=max_workers)


    # output_file = path.RESULTS_DIR / f"bound_promoter_free_energy_results_{PARAM_TYPE}.txt"
    output_file = path.RESULTS_DIR / f"bound_free_energy_results_{PARAM_TYPE}.txt"

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        for result in all_results:
            f.write(f"{result.id}\t{result.F}\n")
    print(f"Results written to {output_file}")

    del all_results, seq_records



    fasta_file = path.DATA_DIR / "unbound_regions/unbound_seq.fa"

    seq_records = list(SeqIO.parse(fasta_file, "fasta"))

    all_results = run_batches(seq_records, batch_size=batch_size, max_physical_core=max_workers)
    output_file = path.RESULTS_DIR / f"unbound_free_energy_results_{PARAM_TYPE}.txt"


    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        for result in all_results:
            f.write(f"{result.id}\t{result.F}\n")
    print(f"Results written to {output_file}")



    end = time.perf_counter()
    print(f"Finished in {round(end - start, 2)} seconds.")