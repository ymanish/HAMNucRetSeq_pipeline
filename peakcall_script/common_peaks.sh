set -euo pipefail
threads=4

# ---------------------------------- paths -----------------------------------
CHAIN=../data/annotations/hg18ToHg38.over.chain.gz                     # UCSC chain
HAM_DIR=../data/hammoud_peaks
# MY_DIR=../data/macs3_mapq13_keepdup
MY_DIR=../data/macs3_mapq30_keepdup_notrim
# MY_DIR=../data/bam_sampled


declare -A HAM_RAW=(
    [d1]=$HAM_DIR/GSE15690_Seq_HistoneMnase_D1.txt
    [pooled]=$HAM_DIR/GSE15690_Seq_HistoneMnase_PooledDonor.txt
)


declare -A MY_SUMMITS=(
    [d1]=$MY_DIR/donor1_summits.bed
    [pooled]=$MY_DIR/pooled_summits.bed
)


# declare -A MY_SUMMITS=(
#     [pooled]=../data/bam_sampled/pooled_summits.bed
# )

txt2bed6 () {               # $1 txt  $2 bed6_out
    awk 'NR>1 {
            OFS="\t";
            id="ham_peak" NR-1;
            score=int($12);  # -10*log10q value (column 5)
            print $2,$3,$4,id,score,"."
         }' "$1" > "$2"
}

# narrow_to_summit () {       # $1 narrowPeak  $2 bed6_out
#     awk '{
#             summit=$2+$10;            # abs summit (column 10)
#             q=$5;               # -10*log10q value (column 5)
#             OFS="\t";
#             print $1,summit,summit+1,$4,q,"."
#          }' "$1" > "$2"
# }


summit_to_bed6 () {       # $1 _summit.bed  $2 bed6_out
    awk '{ q=10*$5;               # -10*log10q value (column 5)
            OFS="\t";
            print $1,$2,$3,$4,q,"."
         }' "$1" > "$2"
}

echo -e "Building overlap tables ....\n"

for tag in "${!HAM_RAW[@]}"; do
    echo "Dataset: $tag"
    ham_hg18=$HAM_DIR/${tag}.bed
    ham_hg38=$HAM_DIR/${tag}_hg38.bed
    ham_unm=$HAM_DIR/${tag}_unmapped.bed

    my_summit=$MY_DIR/${tag}_summits.bed6

    out_tsv=$MY_DIR/${tag}_overlap.tsv

    txt2bed6 "${HAM_RAW[$tag]}" "$ham_hg18"
    liftOver -bedPlus=6 "$ham_hg18" "$CHAIN" "$ham_hg38" "$ham_unm"  ### 6+ coulumns will be copied directly
    rm "$ham_hg18" "$ham_unm"

    summit_to_bed6 "${MY_SUMMITS[$tag]}" "$my_summit"

    bedtools intersect -wa -wb -a "$my_summit" -b "$ham_hg38" | awk 'BEGIN{OFS="\t"}{
             print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
           }' > "$out_tsv"

    tot_my=$(wc -l < "${MY_SUMMITS[$tag]}")
    hmd_peaks=$(wc -l < "$ham_hg38")
    ovl=$(wc -l < "$out_tsv")
    echo -e "My_peaks\t$tot_my\tHammoud Peaks\t$hmd_peaks\tOverlap in them $ovl\t$(awk 'BEGIN{printf "%.2f",('$ovl'/'$tot_my')*100}')%"
done


