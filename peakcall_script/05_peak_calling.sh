set -euo pipefail
ulimit -n 8192  

threads=20                       # same as upstream
MACS_DIR=../data/macs3_mapq30_keepdup
BAM_DIR=../data/bam_groups
# BAM_DIR=../data/bam_groups_notrim

LOG_DIR=../logs
mkdir -p "${MACS_DIR}" "${LOG_DIR}"


###############################################################################
# Peak calling with MACS3                                        
###############################################################################

declare -A groups2=(
  [donor1]="donor1_hist pooled_input"
  [pooled]="pooled_hist pooled_input"
  # [donor2]="donor2_hist d2d3d4_input"
  # [donor3]="donor3_hist d2d3d4_input"
  # [donor4]="donor4_hist d2d3d4_input"
  # [d2d3d4]="d2d3d4_hist d2d3d4_input"
)


for tag in "${!groups2[@]}"; do
  IFS=' ' read -r CHIP CTRL <<< "${groups2[$tag]}"
  echo "[macs3] Calling peaks for ${tag} (ChIP: ${CHIP}, Control: ${CTRL})"
  macs3 callpeak \
      -t "${BAM_DIR}/${CHIP}.mkdup.mapq30.keepdup.bam" \
      -c "${BAM_DIR}/${CTRL}.mkdup.mapq30.keepdup.bam" \
      --nomodel --extsize 150 \
      --min-length 150 \
      --llocal 10000 --slocal 1000 \
      -f BAM --gsize hs -q 0.01 \
      --keep-dup all -B --SPMR \
      --outdir "${MACS_DIR}" -n "${tag}" 
done
