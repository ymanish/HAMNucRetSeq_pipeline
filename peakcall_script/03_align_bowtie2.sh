set -o pipefail ## Exit on any error in a pipeline

threads=21                        
# ulimit -n 8192    ## INCREASE MAX OPEN FILES

REF=../data/annotations/hg38.fa                  
IDX=../data/annotations/hg38

FASTQ_DIR=../data/raw/SRA_trimmed_1/     # *.trim.fastq.gz from fastp step
# FASTQ_DIR=../data/raw/SRA/     # *.trim.fastq.gz from fastp step

BAM_DIR=../data/bam
LOG_DIR=../logs
MET=../data/metrics


# BAM_DIR=../data/bam_notrim
# MET=../data/metrics_notrim


mkdir -p "${BAM_DIR}" "${LOG_DIR}" "${MET}"

###############################################################################
#Build Bowtie2 index (do once)                                            
###############################################################################
# if [[ ! -f "${IDX}.1.bt2" ]]; then
#   echo "[index] Building Bowtie2 index..."
#   bowtie2-build --threads "${threads}" "${REF}" "${IDX}"
# fi

###############################################################################
# Define sample groups                                                     
###############################################################################
SRR_D1_HIST=(SRR019958 SRR019959 SRR019960 SRR019961 SRR019962 SRR019963 SRR019964)
SRR_POOL_HIST=(SRR019965 SRR019966 SRR019967 SRR019968 SRR019969 SRR019970)
SRR_D1_INPUT=(SRR019971 SRR019972 SRR019973)
SRR_POOL_INPUT=(SRR019974 SRR019975 SRR019976 SRR019977)
SRR_ACCESSIONS=("${SRR_D1_HIST[@]}" "${SRR_POOL_HIST[@]}" "${SRR_D1_INPUT[@]}" "${SRR_POOL_INPUT[@]}")

echo "Samples to process: ${#SRR_ACCESSIONS[@]}"

###############################################################################
# Loop over every sample                                                   
###############################################################################


for sample in "${SRR_ACCESSIONS[@]}"; do

  # fq="${FASTQ_DIR}/${sample}.trim.fastq.gz" 
  fq="${FASTQ_DIR}/${sample}.fastq.gz" 

  bam_sort="${BAM_DIR}/${sample}.sort.bam"       # temp sorted BAM
  log="${LOG_DIR}/${sample}.bowtie2.log"
  bam_dup="${BAM_DIR}/${sample}.mkdup.bam"       # duplicate-marked

  ## Skip finished samples ---------------------------------------------------
  if [[ -f "${bam_dup}" ]]; then
    echo "[skip] ${sample} alignment & dup marking already done."
    continue
  fi

  ## 4a. Align & sort --------------------------------------------------------
  echo "[map ] ${sample}"
  bowtie2                                \
      --very-sensitive                   \
      -p "${threads}"                    \
      -x "${IDX}"                        \
      -U "${fq}"                         \
      2> "${log}"                        | samtools sort -@ "${threads}" -o "${bam_sort}" -


  ### Mark duplicates --------------------------------------------------------
  samtools markdup -@"${threads}" -s "${bam_sort}" "${bam_dup}"
  samtools index "${bam_dup}"

  total=$(samtools view -c "${bam_dup}")
  dups=$(samtools view -c -f 1024 "${bam_dup}")
  printf "%s\t%d\t%d\t%.2f%%\n" "$sample" "$total" "$dups" "$(bc -l <<< "$dups*100/$total")" >> "${MET}/dup_summary.tsv"

  # optional: remove intermediates to save space
  rm "${bam_sort}"
  # rm "${bam_sort}" "${bam_dup}" "${bam_dup}.bai"
done


###### Filtering low-quality reads loop
for sample in "${SRR_ACCESSIONS[@]}"; do
  bam_dup="${BAM_DIR}/${sample}.mkdup.bam"                    # input for filtering
  bam_filtered="${BAM_DIR}/${sample}.mkdup.mapq30.keepdup.bam"    # output filtered BAM

  # Only process if the duplicate-marked file exists but the filtered file doesn't
  if [[ ! -f "${bam_dup}" ]]; then
    echo "[skip] ${sample} duplicate-marked BAM not found."
    continue
  fi

  if [[ -f "${bam_filtered}" ]]; then
    echo "[skip] ${sample} already filtered."
    continue
  fi

  echo "[filter] MAPQ≥30 primary alignments"
  # samtools view -@ "${threads}" -b -F 4 -F 256 -q 30 "${bam_dup}" > "${bam_filtered}"
  samtools view -@ "${threads}" -b -F 4 -q 30 "${bam_dup}" > "${bam_filtered}"
  # -F 4: exclude unmapped reads
  # -F 256: exclude secondary alignments
  # -q 30: exclude reads with MAPQ < 30
  samtools index "${bam_filtered}"

done

echo "All samples done.  BAMs in ${BAM_DIR}"




