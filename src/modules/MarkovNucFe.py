# src/utils/markov_fe.py
# Created on 2025-04-06

import numpy as np
from src.modules.nucleosome_positioning import NPBackend
def calc_markov_fe(ID, Sequence, order, file_l, file_s=None):


  # Load in the probability tensor for the 'long' oligonucleotides
  rshptuplong = (148-order,) + (4,)*order

  Pl = np.genfromtxt(file_l).reshape(rshptuplong)

  if ( order > 1 ):

    fileshrt = file_s
    rshptupshrt = (149-order,) + (4,)*(order-1)
    Ps = np.genfromtxt(fileshrt).reshape(rshptupshrt)
    NP = NPBackend(order, 147, Pl, Ps)
  else:
    NP = NPBackend(order, 147, Pl)

  # Calculate the landscape
  p = NP.ProbLandscape(Sequence)
  
  return p, ID

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################




if __name__=='__main__':

  import numpy as np
  import math
  import time
  import concurrent.futures
  import pandas as pd
  from src.config.path import PARAM_DIR, DATA_DIR
  from Bio import SeqIO
  from tqdm import tqdm
  start = time.perf_counter()

  fasta_file = DATA_DIR / "bound_regions/pooled_peaks_bound.fa"
  seq_records = []
  for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
      if i >= 10:
          break
      seq_records.append(record)





  file_l_path =PARAM_DIR / "markov_prob/MD_RoomTemp.trinucdist"
  file_s_path = PARAM_DIR / "markov_prob/MD_RoomTemp.dinucdist"
  markov_order = 3

  results = []
  with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
      futures = []
      for rec in seq_records:
          # rec.id will be used as ID, and rec.seq (converted to str) as Sequence.
          futures.append(
              executor.submit(
                  calc_markov_fe,
                  ID=rec.id,
                  Sequence=str(rec.seq).upper(),
                  order=markov_order,
                  file_l=file_l_path,
                  file_s=file_s_path
              )
          )
      total = len(futures)
      for count, future in enumerate(tqdm(concurrent.futures.as_completed(futures), total=total, desc="Processing sequences")):
          res = future.result()
          print(f"Simulation instance {count} completed:", res)
          results.append(res)

  # Optionally, combine the results in a DataFrame.
  Prob_land_df = pd.DataFrame(results, columns=['id', 'prob_land'])
  print(Prob_land_df)

  end = time.perf_counter()
  print(f'Finished in {round(end - start, 2)} second(s)')




