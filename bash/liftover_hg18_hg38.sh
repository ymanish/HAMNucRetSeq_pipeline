###############################################################################
#  liftover_and_compare.sh
#  1) Convert Hammoud 2009 peak table (hg18) to BED6
#  2) Lift to hg38 with UCSC liftOver
#  3) Sort & merge to ensure ≥300-bp windows
#  4) Compare against your 1-bp dyad peaks with bedtools intersect
###############################################################################
set -euo pipefail
threads=4

# paths -----------------------------------------------------------------------
HAM_DIR=../data/hammoud_peaks
HAM_RAW_D1=$HAM_DIR/GSE15690_Seq_HistoneMnase_D1.txt       
HAM_HG18_D1=$HAM_DIR/Hammoud_hg18_d1.bed
HAM_LFT_D1=$HAM_DIR/Hammoud_hg38_d1.bed
HAM_UNM_D1=$HAM_DIR/Hammoud_unmapped_d1.txt
CHAIN=../data/annotations/hg18ToHg38.over.chain.gz        

HAM_RAW_POOLED=$HAM_DIR/GSE15690_Seq_HistoneMnase_PooledDonor.txt
HAM_HG18_POOLED=$HAM_DIR/Hammoud_hg18_pooled.bed
HAM_LFT_POOLED=$HAM_DIR/Hammoud_hg38_pooled.bed
HAM_UNM_POOLED=$HAM_DIR/Hammoud_unmapped_pooled.txt

if [[ -f "${HAM_LFT_D1%.bed}_sorted.bed" && -f "${HAM_LFT_POOLED%.bed}_sorted.bed" ]]; then
    echo "Both sorted liftover files exist. Skipping liftover steps."
else
    echo "Performing liftover for donor1 peaks..."
    awk 'NR>1 {printf "%s\t%s\t%s\tpeak%s\t0\t+\n", $2, $3, $4, NR-1}' "$HAM_RAW_D1" > "$HAM_HG18_D1"
    liftOver -bedPlus=6 "$HAM_HG18_D1" "$CHAIN" "$HAM_LFT_D1" "$HAM_UNM_D1"
    bedtools sort -i "$HAM_LFT_D1" > "${HAM_LFT_D1%.bed}_sorted.bed"

    echo "Performing liftover for pooled peaks..."
    awk 'NR>1 {printf "%s\t%s\t%s\tpeak%s\t0\t+\n", $2, $3, $4, NR-1}' "$HAM_RAW_POOLED" > "$HAM_HG18_POOLED"
    liftOver -bedPlus=6 "$HAM_HG18_POOLED" "$CHAIN" "$HAM_LFT_POOLED" "$HAM_UNM_POOLED"
    bedtools sort -i "$HAM_LFT_POOLED" > "${HAM_LFT_POOLED%.bed}_sorted.bed"

    rm "$HAM_UNM_D1" "$HAM_UNM_POOLED"  # remove unmapped files
fi



#### Compare your dyad peaks with Hammoud's lifted peaks #######################

MY_PEAKS_D1=../data/macs3_mapq30_nodup/donor1_peaks.narrowPeak  #  1-bp peaks
STATS_D1=compare_stats_d1.txt
MY_PEAKS_POOLED=../data/macs3_mapq30_nodup/pooled_peaks.narrowPeak  #  1-bp peaks
STATS_POOLED=compare_stats_pooled.txt


tot_my=$(wc -l < $MY_PEAKS_D1)
ovl=$(bedtools intersect -u -a $MY_PEAKS_D1 -b ${HAM_LFT_D1%.bed}_sorted.bed | wc -l)
echo -e "My_peaks\t$tot_my\tOverlap_with_Hammoud\t$ovl\t$(awk 'BEGIN{printf "%.2f",('$ovl'/'$tot_my')*100}')%" > $STATS_D1
cat $STATS_D1


tot_my=$(wc -l < $MY_PEAKS_POOLED)
ovl=$(bedtools intersect -u -a $MY_PEAKS_POOLED -b ${HAM_LFT_POOLED%.bed}_sorted.bed | wc -l)
echo -e "My_peaks\t$tot_my\tOverlap_with_Hammoud\t$ovl\t$(awk 'BEGIN{printf "%.2f",('$ovl'/'$tot_my')*100}')%" > $STATS_POOLED
cat $STATS_POOLED

