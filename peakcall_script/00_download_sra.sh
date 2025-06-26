#!/usr/bin/bash

SRA_DIR=../data/raw/SRA
mkdir -p "${SRR_DIR}"

SRR_D1_ACCESSIONS=("SRR019958" "SRR019959" "SRR019960" "SRR019961" "SRR019962" "SRR019963" "SRR019964") ## LANE 1 to 7 all are D1
SRR_POOLED_ACCESSIONS=("SRR019965" "SRR019966" "SRR019967" "SRR019968" "SRR019969" "SRR019970") ##LANE 8 and 9 are D2, LANE 10 and 11 are D3, LANE 12 and 13 are D4.
SRR_D1_INPUT_ACCESSIONS=("SRR019971" "SRR019972" "SRR019973")
SRR_POOLED_INPUT_ACCESSIONS=("SRR019974" "SRR019975" "SRR019976" "SRR019977")

# SRR_H2az_ACCESSION=("SRR019952") ## H2az for Pooled Donors
# SRR_H3K4Me3_ACCESSION=("SRR019953" "SRR019954" "SRR019955") ## 	H3K4Me3 for Pooled Donors 
# SRR_H3K27Me3_ACCESSION=("SRR019956" "SRR019957") ## H3K27Me3 for Pooled Donors

SRR_ACCESSIONS=("${SRR_D1_ACCESSIONS[@]}" "${SRR_POOLED_ACCESSIONS[@]}" "${SRR_D1_INPUT_ACCESSIONS[@]}" "${SRR_POOLED_INPUT_ACCESSIONS[@]}")


THREADS=20    

for SRR in "${SRR_ACCESSIONS[@]}"; do
  echo -e "\n>> Processing ${SRR}"


  if [[ ! -f "${SRA_DIR}/${SRR}/${SRR}.sra" ]]; then
    prefetch ${SRR} -O "${SRA_DIR}"
  else
    echo "SRA file already present."
  fi

  echo "Converting to FASTQ (gzipped)..."
  fasterq-dump \
      --threads "${THREADS}" \
      --outdir  "${SRA_DIR}" \
      --progress \
      "${SRA_DIR}/${SRR}/${SRR}.sra"

   echo "Compressing FASTQ..."
   gzip "${SRA_DIR}/${SRR}.fastq"  # Compress after fasterq-dump

  rm -rf "${SRA_DIR}/${SRR}"
done