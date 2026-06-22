#!/usr/bin/env bash
# Download the small P3 cSCC scRNA example from the SComatic validation dataset.
# Source series: GSE144236 / SRP244706.
set -euo pipefail

INCLUDE_NORMAL="${INCLUDE_NORMAL:-0}"

DATASET_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$DATASET_DIR/sra/P3_scRNA"
mkdir -p "$OUT"

download_one() {
  local run="$1"
  local label="$2"
  local url="$3"
  local expected_bytes="$4"
  local outfile="$OUT/${run}.${label}.sralite"

  if [ -f "$outfile" ] && [ "$(stat -c '%s' "$outfile")" -eq "$expected_bytes" ]; then
    echo "[$(date)] already complete $outfile"
    return
  fi

  echo "[$(date)] downloading $run ($label)"
  curl -L --fail --retry 5 --continue-at - -o "$outfile" "$url"
  if [ "$(stat -c '%s' "$outfile")" -ne "$expected_bytes" ]; then
    echo "Unexpected byte count for $outfile" >&2
    echo "expected=$expected_bytes observed=$(stat -c '%s' "$outfile")" >&2
    exit 1
  fi
  echo "[$(date)] done $outfile"
}

download_one \
  "SRR11832842" \
  "P3_cSCC_1_scRNA" \
  "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR011/832/SRR11832842.sralite.1" \
  "591287409"

if [ "$INCLUDE_NORMAL" = "1" ]; then
  download_one \
    "SRR11832843" \
    "P3_normal_scRNA" \
    "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-16/SRR011/832/SRR11832843.sralite.1" \
    "861675640"
fi
