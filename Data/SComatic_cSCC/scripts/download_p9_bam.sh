#!/usr/bin/env bash
# Download the submitted Cell Ranger BAM for the P9 cSCC tumor scRNA sample.
# Source: ENA submitted files for SRR11832857 / GSM4284244.
set -euo pipefail

DATASET_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$DATASET_DIR/bam/P9_scRNA"
mkdir -p "$OUT"

download_one() {
  local url="$1"
  local outfile="$2"
  local expected_bytes="$3"
  local expected_md5="$4"

  if [ -f "$outfile" ] && [ "$(stat -c '%s' "$outfile")" -eq "$expected_bytes" ]; then
    echo "[$(date)] already complete $outfile"
  else
    echo "[$(date)] downloading $(basename "$outfile")"
    curl -L --fail --retry 5 --continue-at - -o "$outfile" "$url"
  fi

  local observed_bytes
  observed_bytes="$(stat -c '%s' "$outfile")"
  if [ "$observed_bytes" -ne "$expected_bytes" ]; then
    echo "Unexpected byte count for $outfile" >&2
    echo "expected=$expected_bytes observed=$observed_bytes" >&2
    exit 1
  fi

  echo "$expected_md5  $outfile" | md5sum --check -
}

download_one \
  "https://ftp.sra.ebi.ac.uk/vol1/run/SRR118/SRR11832857/P9_cSCC_scRNA.bam" \
  "$OUT/P9_cSCC_scRNA.bam" \
  "30913461365" \
  "62d20a5397b892363a3ad447ae1a94df"

download_one \
  "https://ftp.sra.ebi.ac.uk/vol1/run/SRR118/SRR11832857/P9_cSCC_scRNA.bam.bai" \
  "$OUT/P9_cSCC_scRNA.bam.bai" \
  "11721680" \
  "630718745a700236853e70518aea4b49"
