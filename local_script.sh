NUMBA_CACHE_DIR="$PWD/numba_cache"
mkdir -p "$NUMBA_CACHE_DIR"
export NUMBA_CACHE_DIR
export NUMBA_CPU_NAME=generic 


echo "-> Warming up Numba kernels (one-time per node)...."
singularity exec \
    --env NUMBA_CACHE_DIR="$NUMBA_CACHE_DIR" \
    --bind $PWD:/project \
    hamnucret.sif python3 /project/src/modules/NucFreeEnergy.py 
echo "[+] Warm-up done - native .so files in \$NUMBA_CACHE_DIR"


singularity exec \
  --env NUMBA_CACHE_DIR="$NUMBA_CACHE_DIR" \
  --bind $PWD:/project \
  hamnucret.sif python3 /project/src/core/compute_nucbreath.py 







# echo "-> Warming up Numba kernels (one-time per node)...."
# singularity exec \
#     --env NUMBA_CACHE_DIR="$NUMBA_CACHE_DIR" \
#     --bind $PWD:/project \
#     hamnucret.sif python3 /project/src/modules/AddFreeDNAfe.py 
# echo "[+] Warm-up done - native .so files in \$NUMBA_CACHE_DIR"


# singularity exec \
#   --env NUMBA_CACHE_DIR="$NUMBA_CACHE_DIR" \
#   --bind $PWD:/project \
#   hamnucret.sif python3 /project/src/core/calc_freedna_fe.py 