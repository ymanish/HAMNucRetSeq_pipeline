#!/usr/bin/bash

###############################################################################
# Trim raw sperm MNase-seq reads (GAII, 26 nt input & 36 nt test) and rerun QC
###############################################################################
THREADS=16

FASTQ_IN_DIR=../data/raw/SRA           
FASTQ_OUT_DIR=../data/raw/SRA_trimmed_1
QC_DIR=../data/raw/qc_reports_trimmed_1
mkdir -p "${FASTQ_OUT_DIR}" "${QC_DIR}"


for fq in "${FASTQ_IN_DIR}"/*.fastq.gz; do
    base=$(basename "${fq}" .fastq.gz)
    out="${FASTQ_OUT_DIR}/${base}.trim.fastq.gz"
    html="${FASTQ_OUT_DIR}/${base}.fastp.html"
    json="${FASTQ_OUT_DIR}/${base}.fastp.json"

    if [[ -f "${out}" ]]; then
        echo "[skip] ${base} already trimmed"; continue
    fi

    echo "[trim] ${base}"
    fastp -i "${fq}"                \
          -o "${out}"               \
          --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20 \
          --trim_front1 1 \
          --thread "${THREADS}" \
          --qualified_quality_phred 20 \
          --verbose \
          --html "${html}" --json "${json}"
done

echo -e "\n[qc] FastQC -> MultiQC on trimmed reads"
fastqc -t "${THREADS}" "${FASTQ_OUT_DIR}"/*.fastq.gz -o "${QC_DIR}"
multiqc "${QC_DIR}" -o "${QC_DIR}"

echo "Done. Inspect ${QC_DIR}/multiqc_report.html – all modules should be green."
