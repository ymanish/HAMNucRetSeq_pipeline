QC_DIR="../data/raw/qc_reports"
FASTQ_DIR="../data/raw/SRA"

THREADS=20

echo -e "\n>> Running FastQC on all gz FASTQ ..."
fastqc -t "${THREADS}" "${FASTQ_DIR}"/*.fastq.gz -o "${QC_DIR}"

echo -e "\n>> Aggregating QC reports with MultiQC ..."
multiqc "${QC_DIR}" -o "${QC_DIR}"

echo -e "\nAll downloads, conversions and QC complete!"