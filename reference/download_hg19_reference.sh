#!/usr/bin/env bash
# Download the UCSC hg19 FASTA needed for BAM-first runs on submitted cSCC BAMs.
set -euo pipefail

OUT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/hg19"
mkdir -p "$OUT_DIR"

FA_GZ="$OUT_DIR/hg19.fa.gz"
FA="$OUT_DIR/hg19.fa"
MD5S="$OUT_DIR/md5sum.txt"

curl -L --fail --retry 5 -o "$MD5S" \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/md5sum.txt"

if [ -f "$FA_GZ" ]; then
  echo "[$(date)] resuming/checking $FA_GZ"
else
  echo "[$(date)] downloading $FA_GZ"
fi

curl -L --fail --retry 5 --continue-at - -o "$FA_GZ" \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"

(cd "$OUT_DIR" && grep ' hg19.fa.gz$' md5sum.txt | md5sum --check -)

if [ ! -f "$FA" ]; then
  gunzip -k "$FA_GZ"
fi

samtools faidx "$FA"
