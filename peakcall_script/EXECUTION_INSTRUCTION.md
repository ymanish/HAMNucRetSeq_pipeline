# MNase-seq Peak-Calling Pipeline (Hammoud et al., 2009)

This repository contains a reproducible shell-script pipeline for processing the **human sperm MNase-seq dataset from Hammoud et al. (2009)**.  
It downloads raw data from the Sequence Read Archive (SRA), performs quality control, trimming, alignment, duplicate marking, group merging, and finally peak calling with **MACS3**.  


---
## Directory Layout

data/
    annotations/      # hg38.fa, hg38.chrom.sizes, cpg_islands.bed, ...
    raw/
        SRA/            # *.sra and *.fastq.gz after step 00
        SRA_trimmed_1/  # optional trimmed FASTQs (step 02)
        qc_reports/     # FastQC + MultiQC (raw)
        qc_reports_trimmed_1/  # FastQC + MultiQC (trimmed)
    bam/             # per-SRR aligned BAMs
    bam_groups/ # merged group BAMs
    metrics/     # duplication/read-depth summaries
    macs3_mapq30_keepdup/    # MACS3 outputs
logs/                           # mapping / MACS logs

---

**Recommended conda environment:**

conda create -n mnase \
  python=3.10 \
  sra-tools=3.2.1 fastp=0.23.4 fastqc=0.12.1 multiqc=1.15 \
  bowtie2=2.5.4 samtools=1.21 macs3=3.0.3 \
  bedtools=2.31.1 deeptools=3.5.5 \
  ucsc-liftover ucsc-bedgraphtobigwig trimmomatic=0.39 \
  -c conda-forge -c bioconda


---

## Workflow Overview


0. **Reference annotations**  
   - `get_annotations.sh`  
     Downloads the hg38 reference genome, chromosome sizes, CpG islands, assembly gaps, and repeat masks from UCSC.


1. **Data acquisition**  
   - `00_download_sra.sh`  
     Downloads raw `.sra` files, converts them into gzipped FASTQ files, and cleans up temporary directories.

2. **Quality control (raw reads)**  
   - `01_QC.sh`  
     Runs **FastQC** on raw FASTQ files and aggregates results with **MultiQC**.

3. **Optional trimming**  
   - `02_trim_QC.sh`  
     Uses **fastp** to trim low-quality bases, followed by FastQC/MultiQC to confirm improvements.

4. **Alignment**  
   - `03_align_bowtie2.sh`  
     Aligns reads to the **hg38** reference genome with **Bowtie2**, sorts, marks duplicates with **samtools markdup**, and filters alignments at **MAPQ ≥ 30** (retaining duplicates for biological reasons).

5. **Group merging**  
   - `04_merge_groups.sh`  
     Merges donor-specific BAM files into group BAMs (e.g. donor1, pooled), indexes them, and records total read depth in a summary table.

6. **Peak calling**  
   - `05_peak_calling.sh`  
     Runs **MACS3** with parameters tuned for mononucleosome-sized fragments (~147 bp):  
     - `--nomodel --extsize 150 --min-length 150`  
     - FDR threshold `-q 0.01`  
     - Keeps duplicates (`--keep-dup all`) to preserve true biological enrichments.


---
## Usage

bash get_annotations.sh
bash 00_download_sra.sh
.
.
.
bash 05_peak_calling.sh
 



