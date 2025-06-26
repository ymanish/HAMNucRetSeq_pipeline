

#   Split a FASTA file into equal-sized chunks.
#   You can specify:
#     –c <N>   number of sequences per chunk           (mutually exclusive)#
#   Usage examples
#   -------------
#   # 1) 100 sequences per file
#   ./split_fasta.sh -i big.fa -o chunks/ -c 100
#
# ---------------------------------------------------------------------------

set -euo pipefail
THREADS=2
# ----------------------------- CLI -----------------------------------------
usage() {
  echo "Usage: $0 -i input.fa -o output_dir -c seqs_per_chunk"
  exit 1
}

# getopts processes command-line arguments
# :i:o:c: means it expects options: -i, -o, -c, (colons mean they require values)
# For each option found:
# -i stores value in IN (input file)
# -o stores value in OUT (output directory)
# -c stores value in PER_SEQ (sequences per chunk)
# *) handles any invalid options by showing help

while getopts ":i:o:c:" opt; do
  case $opt in
    i) IN=$OPTARG ;;
    o) OUT=$OPTARG ;;
    c) PER_SEQ=$OPTARG ;;
    *) usage ;;
  esac
done

[[ -z ${IN:-} || -z ${OUT:-} ]] && usage
[[ ! -e $IN ]] && { echo "Error: FASTA $IN not found"; exit 1; }

mkdir -p "$OUT"

# ----------------------------- split logic ---------------------------------
if [[ -n ${PER_SEQ:-} ]]; then
  echo "[+] Splitting $IN into chunks of $PER_SEQ sequences ..."
  num_seqs=$(seqkit stats -T "$IN" | awk 'NR==2{print $4}')
  echo "The number of sequences is $num_seqs"
  parts=$(((num_seqs / PER_SEQ) + 1 ))
  echo "Splitting into $parts parts ..."

  seqkit split "$IN" --threads $THREADS -p $parts -O "$OUT" 
else
  usage
fi

echo "Done — chunks are in $OUT"
