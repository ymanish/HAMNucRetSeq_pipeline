import os
if os.environ.get("IMPORT_ENV_SETTINGS", "1") == "1":
    from src.config.env_settings import *  # Triggers env_settings import

from src.config.path import DATA_DIR, RESULTS_DIR
from typing import List, Tuple, Optional
import itertools
import tqdm 
from concurrent.futures import ProcessPoolExecutor, as_completed

from pathlib import Path
from src.core.helper.breathstates_calculator import total_states_index, total_open_states
import time, logging, datetime as dt
from src.utils.logger_util import get_logger
import shutil

from src.core.helper.bkeep import init_worker

import src.core.nucbreath_core as nc

def _init_BREATHER(nuc_method:str="hybrid", free_dna_method:str=None):
    from src.modules.NucFreeEnergy import NucleosomeBreath
    nuc_breath = NucleosomeBreath(nuc_method=nuc_method,
                                    free_dna_method=free_dna_method)
    return nuc_breath

def batcher(it, size):
    it = iter(it)
    for first in it:
        batch = list(itertools.chain([first], itertools.islice(it, size-1)))
        yield batch

def main(fasta_path: Path, nuc_method:str, *, batch_size: int, 
         n_workers: int, outfile:Path = None, 
         freedna_method:Optional[str] = None, 
         style:str="b_index",
         states: Optional[List[Tuple[int, int]]] = None, flush_every: int = 10000) -> None:


    # generator of all windows
    # windows = nc.process_fasta(fasta_path, win=147, step=1)
    windows = itertools.islice(nc.process_fasta(fasta_path, win=147, step=1), 40)
   
    # chunk generator into batches of BATCH_SIZE

    nuc_breath = _init_BREATHER(nuc_method=nuc_method, free_dna_method=freedna_method)
    nc.BREATHER = nuc_breath

    temp_file_paths = []
    with ProcessPoolExecutor(max_workers=n_workers,
                                initializer=init_worker, 
                                initargs=(flush_every,)) as pool:

        futures = [pool.submit(nc.calc_batch_breathing_energy,
                                                batch,
                                                hard=False,
                                                style=style, 
                                                style_states=states) for batch in batcher(windows, batch_size)
                    ]

        for fut in tqdm.tqdm(as_completed(futures),
                            total=len(futures),
                            desc="Processing batches"):
        
            temp_file_paths.append(fut.result())
            print(f"Temporary file created: {fut.result()}")


    HEADER = ["id", "subid", "sequence", "left_index", "right_index", "F", "F_entropy", "F_enthalpy", "F_freedna"]

    with open(outfile, "wb") as final:
        final.write(("\t".join(HEADER) + "\n").encode())            
        for path in temp_file_paths:
            with open(path, "rb") as src:
                shutil.copyfileobj(src, final)                      
            os.remove(path)                   

    parent_logger.info(f"All done - results written to {outfile}")

    return None



def arg_parser():
    import argparse
    parser = argparse.ArgumentParser(description="Calculate nucleosome free energy from FASTA sequences.")
    parser.add_argument("--infile", type=Path, help="Path to the input FASTA file.")
    parser.add_argument("--batch_size", type=int, default=1, help="Number of sequences per batch.")
    parser.add_argument("--n_workers", type=int, default=2, help="Number of parallel workers.")
    parser.add_argument("--outfile", type=Path, help="Output file path for results.")
    parser.add_argument("--nuc_method", type=str, default="hybrid", help="Method for nucleosome free energy calculation.")
    parser.add_argument("--freedna_method", type=str, default=None, help="Method for free DNA energy calculation.")
    parser.add_argument("--flush_every", type=int, default=10000, help="Number of rows to flush to disk per batch.")
    return parser.parse_args()

if __name__ == "__main__":

    import time 
    start = time.perf_counter()

    parent_logger = get_logger(__name__, level=logging.INFO)


    FASTA = Path(DATA_DIR/"nucbreath_NBdata/promoters/boundprom_minpoints.fa")
    OUTFILE = Path(RESULTS_DIR/"nucbreathfe_minpromoter/bound_free_energy_results.txt")
    tmp_dir = os.path.join(Path(__file__).parent.parent.parent, "temps")
    os.makedirs(tmp_dir, exist_ok=True)
    os.environ["TMPDIR"] = tmp_dir
    parent_logger.info(f"Using temporary directory: {tmp_dir}")



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
    else:
        FREEDNA_METHOD = None

    if args.flush_every:
        flush_every = args.flush_every

    if not FASTA.exists():
        raise FileNotFoundError(f"FASTA file {FASTA} does not exist. Please check the path.")

    if not OUTFILE.parent.exists():
        OUTFILE.parent.mkdir(parents=True, exist_ok=True)

    STYLE = "b_index"  # or "open_sites" or "ph_index"

    if STYLE == "open_sites":
        states = total_open_states()
    elif STYLE == "b_index":
        states = total_states_index(length=14)
    else:
        states = total_states_index(length=28)


    parent_logger.info(f"Processing FASTA file: {FASTA}")
    parent_logger.info(f"Output will be written to: {OUTFILE}")
    parent_logger.info(f"Using nucleosome method: {NUC_METHOD}")
    parent_logger.info(f"Using free DNA method: {FREEDNA_METHOD}")
    parent_logger.info(f"Using style: {STYLE} with states: {len(states)}")
    parent_logger.info(f"Batch size: {batch_size}, Number of workers: {n_workers}")
    parent_logger.info(f"Flush every: {flush_every}")

    # main(FASTA, nuc_method="crystal", freedna_method=FREEDNA_METHOD, batch_size=100, n_workers=21, outfile=OUTFILE, style=STYLE,
    #      states=states)
    main(FASTA, batch_size=batch_size,
          nuc_method=NUC_METHOD, 
         n_workers=n_workers, 
         outfile=Path(OUTFILE), 
         freedna_method=FREEDNA_METHOD, 
        style=STYLE,
         states=states, flush_every=flush_every)

    end = time.perf_counter()
    print(f"Total time taken: {end - start:.2f} seconds")
