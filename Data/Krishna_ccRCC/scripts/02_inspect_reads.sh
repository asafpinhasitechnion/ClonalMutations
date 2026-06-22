#!/usr/bin/env bash
# =============================================================================
# Krishna ccRCC — inspect FASTQ read structure (lengths, barcode/UMI layout)
# =============================================================================
# Records how read lengths and R1/R2 roles were determined for this dataset.
# Works on partial downloads; can also peek at remote files via HTTP range
# requests without downloading the whole file.
# =============================================================================
set -euo pipefail

DATASET_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LABEL="${LABEL:-IpiNivo_Complete_Center}"
DIR="$DATASET_DIR/fastq/$LABEL"

len_dist() {  # $1 = fastq.gz path ; prints length distribution over first 40k reads
  zcat "$1" 2>/dev/null | awk 'NR%4==2{print length($0)} NR>160000{exit}' | sort -n | uniq -c
}

for f in "$DIR"/RNA/*_1.fastq.gz "$DIR"/RNA/*_2.fastq.gz \
         "$DIR"/TCR/*_1.fastq.gz "$DIR"/TCR/*_2.fastq.gz; do
  [ -e "$f" ] || continue
  echo "=== $f ==="; len_dist "$f"
done

# Peek at a remote file's read length WITHOUT full download (HTTP range request):
#   URL=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/063/SRR13806063/SRR13806063_2.fastq.gz
#   curl -s -r 0-500000 "$URL" | zcat 2>/dev/null | awk 'NR%4==2{print length($0)}' | sort | uniq -c
#
# Sanity-check barcode sharing between RNA-R1 and TCR-R1 (first 16bp = 10x CB):
#   zcat RNA/*_1.fastq.gz | awk 'NR%4==2{print substr($0,1,16)} NR>200000{exit}' | sort -u > /tmp/rna_bc
#   zcat TCR/*_1.fastq.gz | awk 'NR%4==2{print substr($0,1,16)} NR>200000{exit}' | sort -u > /tmp/tcr_bc
#   comm -12 /tmp/rna_bc /tmp/tcr_bc | wc -l   # >0 confirms shared cells / barcode in TCR-R1[0:16]
