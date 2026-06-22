#!/usr/bin/env bash
# Download the submitted Cell Ranger BAM for the P3 cSCC tumor sample.
# Source: ENA submitted files for SRR11832842 / GSM4284229.
set -euo pipefail

DATASET_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$DATASET_DIR/bam/P3_scRNA"
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
  "https://ftp.sra.ebi.ac.uk/vol1/run/SRR118/SRR11832842/P3_cSCC_scRNA_1.bam" \
  "$OUT/P3_cSCC_scRNA_1.bam" \
  "3632920739" \
  "1ec67c1a261d5b5cf80bff96ec28df4e"

download_one \
  "https://ftp.sra.ebi.ac.uk/vol1/run/SRR118/SRR11832842/P3_cSCC_scRNA_1.bam.bai" \
  "$OUT/P3_cSCC_scRNA_1.bam.bai" \
  "5704208" \
  "5cfec86031055fffff654112c4129fb7"
