############################################################### Optional Script #####################
################################ Just Jump to 05_peak_calling.sh ######################################

set -o pipefail ## Exit on any error in a pipeline

threads=20                        
ulimit -n 8192    ## INCREASE MAX OPEN FILES

BAM_DIR_GRPS=../data/bam_groups
LOG_DIR=../logs
BAM_SAMPLED_DIR=../data/bam_sampled
REFSIZE=../data/annotations/hg38.chrom.sizes                 

mkdir -p "${BAM_SAMPLED_DIR}" "${LOG_DIR}"

###############################################################################
#  Sampling of the merged BAM groups for Peack calling with MACS3                                      
###############################################################################


######### D1 Sampling for treatment and input #########


# macs3 randsample -i "${BAM_DIR_GRPS}/donor1_hist.mkdup.mapq30.bam" \
#                  -n 9200000 \
#                  -s 42 -f BAM | awk '$2>=0' > "${BAM_SAMPLED_DIR}/donor1_hist9M.bed"

# bedtools sort -i "${BAM_SAMPLED_DIR}/donor1_hist9M.bed" > "${BAM_SAMPLED_DIR}/donor1_hist9M.sorted.bed"
# rm "${BAM_SAMPLED_DIR}/donor1_hist9M.bed"




# bedtools bamtobed -i "${BAM_DIR_GRPS}/donor1_input.mkdup.mapq30.bam" | awk '$2>=0' \
#                             | bedtools sort -i - > "${BAM_SAMPLED_DIR}/donor1_input.sorted.bed"

# rm "${BAM_SAMPLED_DIR}/donor1_input.bed"


######### Pool Sampling for treatment #########
macs3 randsample -i "${BAM_DIR_GRPS}/donor1_hist.mkdup.mapq30.keepdup.bam" \
                 -n 5500000 \
                 -s 42 -f BAM | awk '$2>=0' > "${BAM_SAMPLED_DIR}/donor1_hist5M.bed"

bedtools bamtobed -i "${BAM_DIR_GRPS}/donor2_hist.mkdup.mapq30.keepdup.bam" | awk '$2>=0' > "${BAM_SAMPLED_DIR}/donor2_hist.bed"
bedtools bamtobed -i "${BAM_DIR_GRPS}/donor3_hist.mkdup.mapq30.keepdup.bam" | awk '$2>=0' > "${BAM_SAMPLED_DIR}/donor3_hist.bed"
bedtools bamtobed -i "${BAM_DIR_GRPS}/donor4_hist.mkdup.mapq30.keepdup.bam" | awk '$2>=0' > "${BAM_SAMPLED_DIR}/donor4_hist.bed"

cat "${BAM_SAMPLED_DIR}/donor1_hist5M.bed" \
    "${BAM_SAMPLED_DIR}/donor2_hist.bed" \
    "${BAM_SAMPLED_DIR}/donor3_hist.bed" \
    "${BAM_SAMPLED_DIR}/donor4_hist.bed" | bedtools sort -i - > "${BAM_SAMPLED_DIR}/pooled_hist20M.sorted.bed"

rm "${BAM_SAMPLED_DIR}/donor1_hist5M.bed" "${BAM_SAMPLED_DIR}/donor2_hist.bed" "${BAM_SAMPLED_DIR}/donor3_hist.bed" "${BAM_SAMPLED_DIR}/donor4_hist.bed"


bedtools bamtobed -i "${BAM_DIR_GRPS}/pooled_input.mkdup.mapq30.keepdup.bam" | awk '$2>=0' \
                            | bedtools sort -i - > "${BAM_SAMPLED_DIR}/pooled_input.sorted.bed"




# macs3 callpeak \
#     -t "${BAM_SAMPLED_DIR}/pooled_hist20M.sorted.bed" \
#     -c "${BAM_SAMPLED_DIR}/pooled_input.sorted.bed" \
#     --nomodel --extsize 150 \
#     --min-length 150 \
#     --llocal 10000 --slocal 1000 \
#     -f BED --gsize hs -q 0.01 \
#     --keep-dup all -B --SPMR \
#     --outdir "${BAM_SAMPLED_DIR}" -n "pooled" 




# macs3 callpeak \
#     -t "${BAM_SAMPLED_DIR}/donor1_hist9M.sorted.bed" \
#     -c "${BAM_SAMPLED_DIR}/pooled_input.sorted.bed" \
#     --nomodel --extsize 200 \
#     --min-length 150 \
#     --llocal 10000 --slocal 1000 \
#     -f BED --gsize hs -q 0.01 \
#     --keep-dup all -B --SPMR \
#     --outdir "${BAM_SAMPLED_DIR}" -n "donor1" \
