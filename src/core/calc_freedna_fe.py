import os
if os.environ.get("IMPORT_ENV_SETTINGS", "1") == "1":
    from src.config.env_settings import *  # Triggers env_settings import

import csv, itertools
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from src.modules.AddFreeDNAfe import AddFreeDNAEnergy  
from src.utils.logger_util import get_logger

from src.config.path import DATA_DIR, RESULTS_DIR, EXAMPLES_DATA, EXAMPLES_OUT
from typing import Iterator, List, Tuple
import src.core.helper.bkeep as bk
from src.config.custom_types import FREEDNAResult, SequenceTask
import shutil
import tqdm
import datetime as dt
import psutil

FREEDNAER = None  # Global variable to hold the FreeDNA energy calculator instance

def load_sequence_tasks(tsv_path: Path) -> Iterator[SequenceTask]:
    """
    Read a TSV with columns (id, subid, sequence) and yield one tuple per row.
    """
    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield SequenceTask(
                id=row["id"],
                subid=row["subid"],
                sequence=row["sequence"]
            )

def run_batch_freedna_energy(batch: List[SequenceTask], style: str = "b_index") -> Path:
    """
    Calculate free energy for a batch of sequences.
    """

    start_g = dt.datetime.now()
    proc = psutil.Process(os.getpid())
    
    bk.WORKER_LOGGER.info("Worker PID %d using FREEDNAER: %s", proc.pid, FREEDNAER)
    bk.WORKER_LOGGER.info("Processing batch of %d sequences", len(batch))

    tmpfile, writer = bk.new_batch_writer()
    bk.WORKER_LOGGER.info("Temporary file created: %s", tmpfile)


    results: List[FREEDNAResult] = []
    for task in batch:
        # Calculate free energy for each sequence in the batch
        fb, fu, starts, ends = FREEDNAER.energies_for_all_states(task.sequence, style=style, length=14)
        for i,j,k,l in zip(starts, ends, fb, fu):
            result = FREEDNAResult(id=task.id, subid=task.subid, leftbind_indx=i, rightbind_indx=j, Ffree_bound=k, Ffree_unbound=l)   
            results.append(result)

        if len(results) >= bk.FLUSH_EVERY:
            writer.writerows(results)
            results.clear()

    if results:
        writer.writerows(results)
        results.clear()

    tmpfile.close()

    rss = proc.memory_info().rss / 2**20 ## CONVERT FROM BYTES TO MB
    bk.WORKER_LOGGER.info("Batch of %d done by %s; RSS %.1f MB; t %.1fs",
                len(batch), os.getpid(),
                rss, (dt.datetime.now() - start_g).total_seconds())
    bk.WORKER_LOGGER.info("Temporary file %s written with %d rows", tmpfile.name, len(results))
    return Path(tmpfile.name)



def _batcher(it, size):
    it = iter(it)
    for first in it:
        batch = list(itertools.chain([first], itertools.islice(it, size-1)))
        yield batch

def main(tsv_path: Path, dna_method:str, *, batch_size: int, 
         n_workers: int, outfile:Path = None, 
         flush_every: int = 10000,  # Number of rows to flush to disk
         style:str="b_index") -> None:


    # generator of all windows
    # _sequences = load_sequence_tasks(tsv_path)
    _sequences = itertools.islice(load_sequence_tasks(tsv_path), 10)
   
    # chunk generator into batches of BATCH_SIZE
    global FREEDNAER
    nuc_freedna = AddFreeDNAEnergy(dna_method=dna_method)
    FREEDNAER = nuc_freedna

    temp_file_paths = []
    with ProcessPoolExecutor(max_workers=n_workers,
                                initializer=bk.init_worker,
                                initargs=(flush_every,)) as pool:

        futures = [pool.submit(run_batch_freedna_energy,
                                                batch,
                                                style=style) for batch in _batcher(_sequences, batch_size)
                    ]

        for fut in tqdm.tqdm(as_completed(futures),
                            total=len(futures),
                            desc="Processing batches"):
        
            temp_file_paths.append(fut.result())
            print(f"Temporary file created: {fut.result()}")


    HEADER = ["id", "subid", "left_index", "right_index", "Ffree_bound", "Ffree_unbound"]

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
    parser = argparse.ArgumentParser(description="Calculate free energy of 147bp" \
    "sequences treating them as Free DNA for different combinations.")
    parser.add_argument("--infile", type=Path, help="Input TSV file with format <id> <subid> <sequence>.")
    parser.add_argument("--batch_size", type=int, default=1, help="Number of sequences per batch.")
    parser.add_argument("--n_workers", type=int, default=2, help="Number of parallel workers.")
    parser.add_argument("--outfile", type=Path, help="Output file path for results.")
    parser.add_argument("--dna_method", type=str, default="md", help="Method for DNA free energy calculation.")
    parser.add_argument("--flush_every", type=int, default=10000, help="Number of rows to flush to disk per batch.")
    return parser.parse_args()



if __name__ == "__main__":

    parent_logger = get_logger(__name__)

    TSV_FILE = EXAMPLES_DATA / "unique_001.tsv"
    OUTFILE = EXAMPLES_OUT / "unique_001_freednafe.tsv"

    tmp_dir = os.path.join(Path(__file__).parent.parent.parent, "temps")
    os.makedirs(tmp_dir, exist_ok=True)
    os.environ["TMPDIR"] = tmp_dir
    parent_logger.info(f"Using temporary directory: {tmp_dir}")


    args = arg_parser()
    if args.infile:
        TSV_FILE = args.infile
    if args.batch_size:
        batch_size = args.batch_size
    if args.n_workers:
        n_workers = args.n_workers
    if args.outfile:
        OUTFILE = args.outfile
    if args.dna_method:
        DNA_METHOD = args.dna_method
    if args.flush_every:
        flush_every = args.flush_every

    if not TSV_FILE.exists():
        raise FileNotFoundError(f"TSV file {TSV_FILE} does not exist. Please check the path.")

    if not OUTFILE.parent.exists():
        OUTFILE.parent.mkdir(parents=True, exist_ok=True)

    STYLE = "b_index"  # or "open_sites" or "ph_index"


    parent_logger.info(f"Processing TSV file: {TSV_FILE}")
    parent_logger.info(f"Output will be written to: {OUTFILE}")
    parent_logger.info(f"Using DNA method: {DNA_METHOD}")
    parent_logger.info(f"Using style: {STYLE}")
    parent_logger.info(f"Batch size: {batch_size}, Number of workers: {n_workers}")
    parent_logger.info(f"Flushing every: {flush_every} rows")

    main(tsv_path=TSV_FILE, dna_method=DNA_METHOD, batch_size=batch_size,
         n_workers=n_workers, outfile=OUTFILE, flush_every=flush_every, style=STYLE)