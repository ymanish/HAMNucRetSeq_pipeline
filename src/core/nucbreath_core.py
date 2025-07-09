import os
import csv
import psutil
import datetime as dt
from pathlib import Path
from typing import List, Tuple, Iterator
from Bio import SeqIO
from src.config.custom_types import ProcessedSequence
from src.utils.fasta_checks import contains_non_canonical
from src.utils.logger_util import get_logger
# from line_profiler import profile   # When running kernprof this decorator is active
import logging
from src.core.helper.fn_wraps import timeit
import tempfile

BREATHER = None
@timeit
def _init_worker():
    """Run once in every child before it executes any task."""

    global WORKER_LOGGER, SCRATCH_PATH, TSV_WRITER, FLUSH_EVERY

    FLUSH_EVERY = 1000
    proc = psutil.Process(os.getpid())
    WORKER_LOGGER = get_logger(f"Worker-{os.getpid()}", log_file=None, level=logging.INFO)


    scratch_tmpdir = os.environ.get("SLURM_TMPDIR") or os.environ.get("TMPDIR")
    if not scratch_tmpdir or not Path(scratch_tmpdir).is_dir():
        raise RuntimeError("Scratch directory not found")

    tmpfile = tempfile.NamedTemporaryFile(mode="w", dir=scratch_tmpdir,
                                          prefix=f"wbatch-{os.getpid()}-", suffix=".tsv",
                                          delete=False, newline='')
    
    SCRATCH_PATH = Path(tmpfile.name)
    TSV_WRITER = csv.writer(tmpfile, delimiter="\t",
                            lineterminator="\n")
               
    WORKER_LOGGER.info("Worker %s online - RSS %.1f MB",
                 os.getpid(), proc.memory_info().rss / 2**20)
    WORKER_LOGGER.info("Temporary file created at %s", SCRATCH_PATH)    

@timeit    
def calc_batch_breathing_energy(batch:List[ProcessedSequence],
                                        *,
                                    hard:bool=False, 
                                    style:str="b_index", 
                                    style_states:List[Tuple[int, int]])-> Path:

    global BREATHER
    start_g = dt.datetime.now()
    proc = psutil.Process(os.getpid())

    # results:List[NuclBreathingResult]=[]
    row_cache = []
    if hard:
        WORKER_LOGGER.info("Using hard method")
        compute_energy = lambda subseq, seq_id, sub_id, left_s, right_s: BREATHER.calculate_free_energy_hard(seq147=subseq,
                                                                                            left=left_s, 
                                                                                            right=right_s, 
                                                                                            id=seq_id, 
                                                                                            subid=sub_id)
    else:
        WORKER_LOGGER.info("Using soft method")
        compute_energy = lambda subseq, seq_id, sub_id, left_s, right_s: BREATHER.calculate_free_energy_soft(seq601=subseq, 
                                                                                              left=left_s, 
                                                                                            right=right_s, 
                                                                                              id=seq_id, 
                                                                                              subid=sub_id,
                                                                                                style=style)

    for rec in batch:
        # WORKER_LOGGER.info(f"Processing sequence: {rec.id} subid: {rec.subid} with length {len(rec.sequence)}")
        for left_s, right_s in style_states:
            # start = dt.datetime.now()

            res = compute_energy(subseq = rec.sequence, seq_id = rec.id, sub_id = rec.subid, left_s=left_s, right_s=right_s)

            # rss = proc.memory_info().rss / 2**20 ## CONVERT FROM BYTES TO MB
            # WORKER_LOGGER.info("Single Iteration of %d done by %s; RSS %.1f MB; t %.3fs",
            #     len(batch), os.getpid(),
            #     rss, (dt.datetime.now() - start).total_seconds())
        
            #     res_nucbreath = NuclBreathingResult(
            #     id=rec.id,
            #     subid=rec.subid,
            #     sequence=rec.sequence,
            #     leftbind_indx=left_s,
            #     rightbind_indx=right_s,
            #     F_vals=res
            # )
            #     # print(f"id:{rec.id} sub:{rec.subid}  F:{res.F:.2f}")
            #     results.append(res_nucbreath)


            row_cache.append([rec.id, rec.subid, rec.sequence,
                              left_s, right_s,
                            res.F, res.F_entropy, res.F_enthalpy, res.F_freedna])
            
        if len(row_cache) >= FLUSH_EVERY:
            TSV_WRITER.writerows(row_cache)
            row_cache.clear()

    if row_cache:
        TSV_WRITER.writerows(row_cache)
        row_cache.clear()

    rss = proc.memory_info().rss / 2**20 ## CONVERT FROM BYTES TO MB
    WORKER_LOGGER.info("Batch of %d done by %s; RSS %.1f MB; t %.1fs",
                len(batch), os.getpid(),
                rss, (dt.datetime.now() - start_g).total_seconds())
    return SCRATCH_PATH

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

