from src.modules.NucFreeEnergy import NucleosomeBreath
from src.core.helper.breathstates_calculator import total_states_index, total_open_states
from src.core.nucbreath_core import calc_batch_breathing_energy
from src.config.custom_types import ProcessedSequence
from pathlib import Path
from typing import List, Tuple, Optional
import random
from src.utils.logger_util import get_logger
import logging
import src.core.nucbreath_core as nc
import os

def dummy_sequence_generator(n:int) -> List[ProcessedSequence]:
    """Generate dummy sequences for testing."""
    for i in range(n):
        yield ProcessedSequence(
            id=f"seq{i+1}",
            subid=f"sub{i+1}",
            sequence=''.join(random.choice(['A', 'T', 'C', 'G']) for _ in range(147)),
            start_site=0,
            end_site=147
        )

if __name__ == "__main__":
    parent_logger = get_logger(__name__, level=logging.INFO)

    tmp_dir = os.path.join(Path(__file__).parent.parent.parent, "temps")
    os.makedirs(tmp_dir, exist_ok=True)
    os.environ["TMPDIR"] = tmp_dir
    parent_logger.info(f"Using temporary directory: {tmp_dir}")




    states = total_states_index(length=14)
    batch = list(dummy_sequence_generator(1))

    parent_logger.info(f"Total states: {len(states)}")
    parent_logger.info(f"Batch: {batch}")

    nc._init_worker(method_nuc="crystal", free_dna_method=None)

    print(nc.BREATHER)
    print(nc.SCRATCH_PATH)
    print(nc.TSV_WRITER)

    nc.calc_batch_breathing_energy(batch,
                                    hard=False,
                                    style="b_index",
                                    style_states=states)