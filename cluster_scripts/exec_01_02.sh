#!/bin/bash
module load apps/singularity
# Run the 01_split_fasta.sh script and wait for it to finish
FASTA_FILE="/group/pol_schiessel/Manish/HAMNucRetSeq_pipeline/data/bound_regions/pooled_peaks_bound.fa"
OUTPUT_DIR="/group/pol_schiessel/Manish/HAMNucRetSeq_pipeline/data/bound_regions/chunks"
PER_FILE_SEQ=1000  # Number of sequences per chunk

if [ ! -f "$FASTA_FILE" ]; then
    echo "ERROR: FASTA file $FASTA_FILE not found. Please ensure the file exists before running this script."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Starting 01_split_fasta.sh..."

# singularity exec --bind $PWD:/project hamnucret.sif bash /project/01_split_fasta.sh -i "$FASTA_FILE" -o "$OUTPUT_DIR" -c "$PER_FILE_SEQ"

if [ $? -eq 0 ]; then
    echo "01_split_fasta.sh completed successfully."
    CHUNK_FILES=(${OUTPUT_DIR}/*.fa)
    TOTAL_CHUNKS=${#CHUNK_FILES[@]}

    echo "Total chunks created: $TOTAL_CHUNKS"
else
    echo "ERROR: 01_split_fasta.sh encountered an error. Aborting submission."
    exit 1
fi

if [ $? -eq 0 ]; then
    echo "Submitting 01_compute_nucfe.job..."

    sbatch --array=1-${TOTAL_CHUNKS}%20 --time=30:00:00 02_compute_nucfe.job 
    # sbatch --time=30:00:00 02_compute_nucfe.job 

else
    echo "ERROR: 01_split_fasta.sh encountered an error. Aborting submission."
    exit 1
fi


# awk 'BEGIN {FS=OFS="\t"} FNR==1 && NR!=1{next} {print}' *.tsv > combined_file.tsv