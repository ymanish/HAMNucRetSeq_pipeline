#!/usr/bin/bash

cd ../data/annotations/
# Download annotation files (hg38 example)
# wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz"

if [ ! -f "cpg_islands.bed" ]; then
    url_cpg="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz"
    echo "Downloading and converting cpgIslandExt.txt.gz..."
    wget -qO- "$url_cpg" | gunzip -c | awk '{print $2"\t"$3"\t"$4"\t"$5" "$6}' > cpg_islands.bed
fi


if [ ! -f "gap_hg38.bed" ]; then
    url_gap="https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz"
    echo "Downloading and converting gap.txt.gz..."
    wget -qO- "$url_gap" | gunzip -c | awk '{print $2"\t"$3"\t"$4}' > gap_hg38.bed
fi

if [ ! -f "rmsk_hg38.bed" ]; then
    url_rmsk="https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
    echo "Downloading and converting rmsk.txt.gz..."
    wget -qO- "$url_rmsk" | gunzip -c | awk '{print $6"\t"$7"\t"$8}' > rmsk_hg38.bed
fi


# wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz 
# # wget https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/e106v47/appris_data.appris.txt

# wget https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/e106v47/appris_data.principal.txt
# # # Uncompress files
# gunzip *.txt.gz

# # Process annotation files
# awk 'BEGIN {OFS="\t"} {print $2,$3,$4,$5}' cpgIslandExt.txt > cpg_islands.bed
# awk 'BEGIN {OFS="\t"} {print $6,$7,$8,$11}' rmsk.txt > repeats.bed
# awk 'BEGIN {OFS="\t"} {print $3,$5,$6,$13}' refGene.txt > genes.bed

# echo "Annotation files downloaded and processed."