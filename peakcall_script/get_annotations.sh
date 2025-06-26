#!/usr/bin/bash

cd ../data/annotations/
# Download annotation files (hg38 example)
# wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz"
# wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
# wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"

# # Uncompress files
# gunzip *.txt.gz

# # Process annotation files
awk 'BEGIN {OFS="\t"} {print $2,$3,$4,$5}' cpgIslandExt.txt > cpg_islands.bed
awk 'BEGIN {OFS="\t"} {print $6,$7,$8,$11}' rmsk.txt > repeats.bed
awk 'BEGIN {OFS="\t"} {print $3,$5,$6,$13}' refGene.txt > genes.bed

echo "Annotation files downloaded and processed."