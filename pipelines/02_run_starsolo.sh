#!/usr/bin/env bash
# Align a 10x 5' scRNA-seq sample with STARsolo -> tagged BAM (CB/UB/nM/NH) + count
# matrices. Configured for 10x 5' v1/v2 chemistry with CellRanger>=4 emulation so
# output matches the published Cell Ranger processing closely. Run in `align` env.
#
# Usage:
#   THREADS=32 ./02_run_starsolo.sh <R1.fastq.gz> <R2.fastq.gz> <out_dir> [sample_prefix]
# 10x layout: R1 = barcode+UMI (16+10), R2 = cDNA. STAR wants cDNA first.
set -euo pipefail

R1="$1"; R2="$2"; OUT="$3"; PREFIX="${4:-sample_}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IDX="$REPO/reference/star_index_GRCh38_gencode_v32"
WL="$REPO/reference/737K-august-2016.txt"      # 10x 5' v1/v2 whitelist
THREADS="${THREADS:-32}"
mkdir -p "$OUT"

# Chemistry: CB = R1[1..16], UMI = R1[17..26]. (Our R1 is 27bp; the extra base is
# ignored. --soloBarcodeReadLength 0 disables the strict R1-length==26 check.)
STAR --runMode alignReads \
  --runThreadN "$THREADS" \
  --genomeDir "$IDX" \
  --readFilesIn "$R2" "$R1" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$OUT/$PREFIX" \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist "$WL" \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 \
  --soloBarcodeReadLength 0 \
  --soloStrand Forward \
  --soloFeatures Gene GeneFull \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM CB UB GX GN \
  --outSAMprimaryFlag AllBestScore \
  --clipAdapterType CellRanger4 \
  --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloCellFilter EmptyDrops_CR \
  --limitBAMsortRAM 40000000000

samtools index -@ "$THREADS" "$OUT/${PREFIX}Aligned.sortedByCoord.out.bam"
echo "DONE -> $OUT/${PREFIX}Aligned.sortedByCoord.out.bam"
echo "matrices -> $OUT/${PREFIX}Solo.out/"
