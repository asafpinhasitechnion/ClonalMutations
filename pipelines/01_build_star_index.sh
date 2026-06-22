#!/usr/bin/env bash
# Build a STAR genome index from the 10x GRCh38-2020-A reference (genome + GTF).
# One-time, reused for every sample. ~25-40 min, needs ~30 GB RAM. Run in `align` env.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="$REPO/reference/gencode_v32"               # GENCODE v32 (= basis of 10x 2020-A)
IDX="$REPO/reference/star_index_GRCh38_gencode_v32"
THREADS="${THREADS:-32}"
SJDB_OVERHANG="${SJDB_OVERHANG:-91}"            # cDNA(R2) read length 92 - 1

FA="$REF/GRCh38.primary_assembly.genome.fa"
GTF="$REF/gencode.v32.primary_assembly.annotation.gtf"
# STAR needs uncompressed FASTA/GTF
[ -f "$FA" ]  || gunzip -kf "$FA.gz"
[ -f "$GTF" ] || gunzip -kf "$GTF.gz"
mkdir -p "$IDX"

STAR --runMode genomeGenerate \
  --runThreadN "$THREADS" \
  --genomeDir "$IDX" \
  --genomeFastaFiles "$FA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang "$SJDB_OVERHANG" \
  --limitGenomeGenerateRAM 60000000000

echo "STAR index built at: $IDX"
