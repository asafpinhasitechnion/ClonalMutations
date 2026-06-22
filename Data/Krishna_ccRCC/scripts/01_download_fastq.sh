#!/usr/bin/env bash
# =============================================================================
# Krishna ccRCC (PRJNA705464) — FASTQ download from ENA over HTTPS
# =============================================================================
# Records the exact commands used to fetch raw reads for this dataset, so the
# download is reproducible. Each dataset keeps its own scripts/ dir because
# download/processing details differ between datasets.
#
# Usage:
#   ./01_download_fastq.sh                 # default pair (IpiNivo Complete Center)
#   ./01_download_fastq.sh RNA=SRRxxxxxxx TCR=SRRyyyyyyy LABEL=Some_Region
#
# Notes:
#   - ENA fastq_ftp paths are ftp.sra.ebi.ac.uk/... ; prefix https:// for HTTPS.
#   - wget -q -c resumes partial files (safe to re-run).
#   - 10x 5' chemistry. Read structure is documented per-run in the README that
#     this dataset's fastq dir carries (RNA: R1=27bp barcode+UMI / R2=92bp cDNA;
#     TCR: R1=R2=151bp, barcode+UMI in first 26bp of R1).
# =============================================================================
set -euo pipefail

# ---- parameters (override via KEY=VALUE args) -------------------------------
RNA="${RNA:-SRR13806063}"   # 5' gene-expression run
TCR="${TCR:-SRR13806031}"   # 5' V(D)J / TCR run ("-" to skip when using published clonotypes)
LABEL="${LABEL:-IpiNivo_Complete_Center}"
for kv in "$@"; do export "$kv"; done   # allow RNA=.. TCR=.. LABEL=.. on CLI

DATASET_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"   # Data/Krishna_ccRCC
OUT="$DATASET_DIR/fastq/$LABEL"
mkdir -p "$OUT/RNA" "$OUT/TCR"
LOG="$OUT/download.log"

# ENA HTTPS path helper: vol1/fastq/<first6>/<3-digit-suffix>/<acc>/<acc>_N.fastq.gz
ena_url() {                       # $1 = run accession, $2 = read (1|2)
  local acc="$1" read="$2"
  printf 'https://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%03d/%s/%s_%s.fastq.gz' \
    "${acc:0:6}" "$((10#${acc: -3}))" "$acc" "$acc" "$read"
}

echo "[$(date)] download start: RNA=$RNA TCR=$TCR -> $OUT" | tee -a "$LOG"
wget -q -c "$(ena_url "$RNA" 1)" -O "$OUT/RNA/${RNA}_1.fastq.gz" 2>>"$LOG"
wget -q -c "$(ena_url "$RNA" 2)" -O "$OUT/RNA/${RNA}_2.fastq.gz" 2>>"$LOG"
if [ "$TCR" != "-" ]; then   # skip raw TCR when published clonotypes are used
  wget -q -c "$(ena_url "$TCR" 1)" -O "$OUT/TCR/${TCR}_1.fastq.gz" 2>>"$LOG"
  wget -q -c "$(ena_url "$TCR" 2)" -O "$OUT/TCR/${TCR}_2.fastq.gz" 2>>"$LOG"
else
  rmdir "$OUT/TCR" 2>/dev/null || true
  echo "[$(date)] TCR=- : skipped raw TCR (using published clonotype table)" | tee -a "$LOG"
fi
echo "[$(date)] download complete" | tee -a "$LOG"

# ---- optional integrity check against ENA md5 -------------------------------
# Fetch expected md5s and verify (uncomment to run):
# curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RNA},${TCR}&result=read_run&fields=run_accession,fastq_md5,fastq_ftp&format=tsv"
# md5sum "$OUT"/RNA/*.fastq.gz "$OUT"/TCR/*.fastq.gz
