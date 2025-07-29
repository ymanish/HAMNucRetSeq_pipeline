import os
if os.environ.get("IMPORT_ENV_SETTINGS", "1") == "1":
    from src.config.env_settings import * 

from src.modules import NUC_STATE_PATH ###import just to set the path from init.py otherwise need to explicitly set it
from methods import GenStiffness
import numpy as np
import numba as nb
from functools import lru_cache
from src.core.helper.breathstates_calculator import total_states_index, total_open_states
from line_profiler import profile  # When running kernprof this decorator is active

@nb.jit(nopython=True, cache=True)
def _batch_logdet(stiff: np.ndarray,
                starts: np.ndarray,
                ends:   np.ndarray,
                log_2pi: float)-> tuple[np.ndarray, np.ndarray]:

    n = starts.size
    Fb = np.empty(n)
    Fu = np.empty(n)
    N = stiff.shape[0]

    for k in range(n):
        s = starts[k]
        e = ends[k]

        ##### ----- bound block ------------------------------------------------ #
        B = stiff[s:e, s:e]
        Lb = np.linalg.cholesky(B)
        logdet_b = 2.0 * np.sum(np.log(np.diag(Lb)))
        Fb[k] = -0.5 * B.shape[0] * log_2pi + 0.5 * logdet_b

        #### ----- unbound block ------------------------------------------ #
        idx = np.concatenate((np.arange(s), np.arange(e, N))).astype(np.int32)

        tmp = stiff[idx, :]
        UB  = tmp[:, idx]
        Lu = np.linalg.cholesky(UB)
        logdet_ub = 2.0 * np.sum(np.log(np.diag(Lu)))

        Fu[k] = -0.5 * UB.shape[0] * log_2pi + 0.5 * logdet_ub
    return Fb, Fu







class AddFreeDNAEnergy:
    def __init__(self, dna_method:str='hybrid'):
        self.dna_method = dna_method
        print(f"Using free DNA method: {self.dna_method}")

        self.genstiff_nuc = GenStiffness(method=self.dna_method)
        self.phosphate_bind_sites = self._select_phosphate_bind_sites()

    @staticmethod
    def _select_phosphate_bind_sites(left=0, right=13):
        
        phosphate_bind_sites = [2, 6, 14, 17, 24, 29, 34, 38, 
                                    45, 49, 55, 59, 65, 69, 76, 
                                    80, 86, 90, 96, 100, 107, 111, 
                                    116, 121, 128, 131, 139, 143]
        
        return phosphate_bind_sites[left*2:(right*2)+2]
    

    @lru_cache(maxsize=1)
    def _cached_gen_params(self, seq_str: str):
        """Pure-python wrapper around the expensive call; result is memoised."""
        return self.genstiff_nuc.gen_params(seq_str, use_group=True, sparse=False)

    @staticmethod
    def _get_left_right_open(left: int, right: int, *, style: str = "b_index") -> tuple[int, int]:
        if style == "b_index":     # bound-site indices
            return 2 * left, 28 - (2 * right) - 2
        if style == "ph_index":    # phosphate-site indices
            return left, 28 - right - 1
        if style == "open_sites":  # already counts of open sites
            return left, right
        raise ValueError("style must be 'b_index', 'ph_index' or 'open_sites'")

    def get_freedna_energy_single_state(self, fullseq: str, left:int, right:int, style:str="b_index" ) -> tuple[float, float]:
        """ Calculate energy of the free DNA-region not bound to the nucleosome.
        This is used in the case where I run the binding_model_free_energy on the
        regions that is bound to the nucleosome and the free DNA region is not included in the calculation.
        This is done to increase the speed of the calculation."""
        
        stiff, _ = self._cached_gen_params(fullseq)
        
        l_open, r_open = self._get_left_right_open(left=left, right=right, style=style)
        print(f"Left open: {l_open}, Right open: {r_open}, Style: {style}")



        # bound_locs = self._select_phosphate_bind_sites(l_open, r_open)
        bound_locs = self.phosphate_bind_sites[l_open:len(self.phosphate_bind_sites)-r_open]
        start = 6 * bound_locs[0]
        end = 6 * (bound_locs[-1] + 1)
        print(f"Start: {start}, End: {end}, Length: {stiff.shape[0]}")
        bound_block_stiff = stiff[start:end, start:end]

        unbound_indices = np.r_[0:start, end:stiff.shape[0]]
        unbound_stiffness = stiff[unbound_indices][:, unbound_indices]


        logdet_sign_b, logdet_b = np.linalg.slogdet(bound_block_stiff)
        Fe_bound_dna = -0.5*len(bound_block_stiff)*np.log(2*np.pi) + 0.5*logdet_b
        
        
        logdet_sign_ub, logdet_ub = np.linalg.slogdet(unbound_stiffness)
        Fe_unbound_dna = -0.5*len(unbound_stiffness)*np.log(2*np.pi) + 0.5*logdet_ub

        print(f"Free DNA energy: {Fe_unbound_dna}, Bound DNA energy: {Fe_bound_dna}")

        return Fe_bound_dna, Fe_unbound_dna


    def energies_for_all_states(self,
                                fullseq: str,
                                style: str = "b_index",
                                length: int = 14) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

        stiff, _ = self._cached_gen_params(fullseq)

        if style == "open_sites":
            states = total_open_states()
        elif style == "b_index":
            states = total_states_index(length=14)
        else:
            states = total_states_index(length=28)

        starts, ends = [], []
        lindx, rindx = [], []
        for left, right in states:
            l_open, r_open = self._get_left_right_open(left, right, style=style)
            bound_locs = self.phosphate_bind_sites[l_open : len(self.phosphate_bind_sites) - r_open]
            
            starts.append(6 * bound_locs[0])
            ends.append(6 * (bound_locs[-1] + 1))


            lindx.append(left)
            rindx.append(right)

        starts = np.asarray(starts, dtype=np.int32)
        ends   = np.asarray(ends,   dtype=np.int32)

        log_2pi = np.log(2.0 * np.pi)
        Fb, Fu = _batch_logdet(stiff, starts, ends, log_2pi)
        return Fb, Fu, lindx, rindx




if __name__ == "__main__":
    import time
    start = time.perf_counter()
    Seq601 = "CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT"

    free_dna_energy = AddFreeDNAEnergy(dna_method='hybrid')
    # for i in range(0, 14):
    #     for j in range(0, 14):
    #         if j >= i:
    #             result = free_dna_energy.get_freedna_energy_single_state(fullseq=Seq601, left=i, right=j, style="b_index")
    #             print(f"L:{i}, R:{j} ------ Free DNA energy: {result}")

    
    energies_for_all_states = free_dna_energy.energies_for_all_states(fullseq=Seq601, style="b_index")
    print(f"Total states: {len(energies_for_all_states[0])}")
    # for fb, fu, starts, ends in zip(*energies_for_all_states):
    #     print(f"Fb: {fb}, Fu: {fu}, Starts: {starts}, Ends: {ends}")

    end = time.perf_counter()
    print(f"Time taken: {end - start:.2f} seconds")