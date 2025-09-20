set -o pipefail ## Exit on any error in a pipeline

threads=20                        
ulimit -n 8192    ## INCREASE MAX OPEN FILES

REF=../data/annotations/hg38.fa 
REFSIZE=../data/annotations/hg38.chrom.sizes                 
IDX=../data/annotations/hg38
# MET=../data/metrics_notrim
MET=../data/metrics

BAM_DIR=../data/bam
BAM_DIR_GRPS=../data/bam_groups

# BAM_DIR=../data/bam_notrim
# BAM_DIR_GRPS=../data/bam_groups_notrim
LOG_DIR=../logs

mkdir -p "${BAM_DIR}" "${LOG_DIR}" "${MET}" "${BAM_DIR_GRPS}"
# STATS_FILE="${MET}/group_depths_mapq13_nodup.tsv"
STATS_FILE="${MET}/group_depths_mapq30_keepdup.tsv"

echo -e "group\treads_M" > "$STATS_FILE"    # header line
##############################################################################
#  Merging the Lanes                                             
###############################################################################
SRR_D1_HIST=(SRR019958 SRR019959 SRR019960 SRR019961 SRR019962 SRR019963 SRR019964)
SRR_D2_HIST=(SRR019965 SRR019966)
SRR_D3_HIST=(SRR019967 SRR019968)
SRR_D4_HIST=(SRR019969 SRR019970)
SRR_POOL_HIST="${SRR_D1_HIST[*]} ${SRR_D2_HIST[*]} ${SRR_D3_HIST[*]} ${SRR_D4_HIST[*]}"

SRR_D1_INPUT=(SRR019971 SRR019972 SRR019973)
SRR_POOL_INPUT=(SRR019974 SRR019975 SRR019976 SRR019977)

declare -A groups1=(
  [donor1_hist]="${SRR_D1_HIST[*]}"
  [donor1_input]="${SRR_D1_INPUT[*]}"
  [pooled_hist]="${SRR_D1_HIST[*]} ${SRR_POOL_HIST[*]}"
  [pooled_input]="${SRR_D1_INPUT[*]} ${SRR_POOL_INPUT[*]}"
  [donor2_hist]="${SRR_D2_HIST[*]}"
  [donor3_hist]="${SRR_D3_HIST[*]}"
  [donor4_hist]="${SRR_D4_HIST[*]}"
  [d2d3d4_hist]="${SRR_D2_HIST[*]} ${SRR_D3_HIST[*]} ${SRR_D4_HIST[*]}"
  [d2d3d4_input]="${SRR_POOL_INPUT[*]}"
)
echo "Samples to process: ${#groups1[@]}"

for grp in "${!groups1[@]}"; do
    echo "Processing group: $grp"
    echo "Members: ${groups1[$grp]}"

    fout="${BAM_DIR_GRPS}/${grp}.mkdup.mapq30.keepdup.bam"
    
    [[ -f $fout ]] || {
        echo "[merge] $grp"
        
        bam_list=()
        for id in ${groups1[$grp]}; do
            bam_list+=("${BAM_DIR}/${id}.mkdup.mapq30.keepdup.bam")
        done
        
        samtools merge -@${threads} -o "$fout" "${bam_list[@]}"
        samtools index -@${threads} "$fout"
    }
    
    
    # ----------  stats  ----------
    total=$(samtools idxstats "$fout" | awk '{s+=$3} END{print s}') 
    reads_M=$(awk -v n="$total" 'BEGIN{printf "%.3f", n/1e6}')
    # tag_bp=$(samtools stats "$fout" | awk '$1=="RL" {print $3}')
    # echo -e "${grp}\t${reads_M}\t${tag_bp}" >> "$STATS_FILE"
    echo -e "${grp}\t${reads_M}" >> "$STATS_FILE"

done


###############################################################################
################## mapq30 no duplicates ##########################
######## group	reads_M

######## donor1_input	9.298
######## pooled_input	20.978
######## d2d3d4_input	11.680

######## donor1_hist	21.027
######## donor2_hist	5.531
######## donor3_hist	5.040
######## donor4_hist	4.813
######## pooled_hist	57.438
######## d2d3d4_hist	15.384


################## mapq13 all duplicates ##########################


######## group	reads_M

######## donor1_input	9.391
######## pooled_input	21.173
######## d2d3d4_input	11.783

######## donor1_hist	21.559
######## donor2_hist	5.657
######## donor3_hist	5.139
######## donor4_hist	4.918
######## pooled_hist	58.833
######## d2d3d4_hist	15.714

################## mapq30 all duplicates ##########################
######### group	reads_M

######### donor1_input	9.298
######### pooled_input	20.978
######### d2d3d4_input	11.680

######### donor1_hist	21.027
######### donor2_hist	5.531
######### donor3_hist	5.040
######### donor4_hist	4.813
######### pooled_hist	57.438
######### d2d3d4_hist	15.384

################## mapq13 no duplicates ##########################

######### group	reads_M

######### donor1_input	9.391
######### pooled_input	21.173
######### d2d3d4_input	11.783

######### donor1_hist	21.559
######### donor2_hist	5.657
######### donor3_hist	5.139
######### donor4_hist	4.918
######### pooled_hist	58.833
######### d2d3d4_hist	15.714



