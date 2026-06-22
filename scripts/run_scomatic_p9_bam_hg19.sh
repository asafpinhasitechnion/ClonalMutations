#!/usr/bin/env bash
# BAM-first SComatic run for P9 cSCC submitted BAM.
#
# This is a pragmatic pre-realignment pass. The submitted BAM is hg19/GRCh37
# style, whereas the paper's exact benchmark used Cell Ranger v6 on GRCh38.
# Because local SComatic RNA-editing/PoN resources are hg38, step2 is run without
# those files here. Use this for plumbing and approximate call/support checks;
# use the later GRCh38 realignment for paper-faithful results.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCOMATIC="$ROOT/external/SComatic"

sample="P9_cSCC"
bam="$ROOT/Data/SComatic_cSCC/bam/P9_scRNA/P9_cSCC_scRNA.bam"
expected_bytes="30913461365"
meta="$ROOT/results/SComatic_cSCC/P9_cSCC/P9_cSCC.scomatic_meta.tsv"
ref="$ROOT/reference/hg19/hg19.fa"
out="$ROOT/results/SComatic_cSCC/P9_cSCC/scomatic_hg19"

if [ ! -f "$bam" ] || [ "$(stat -c '%s' "$bam")" -ne "$expected_bytes" ]; then
  echo "P9 BAM is not complete yet: $bam" >&2
  exit 1
fi
if [ ! -f "$bam.bai" ]; then
  echo "P9 BAM index is missing: $bam.bai" >&2
  exit 1
fi
if [ ! -f "$ref.fai" ]; then
  echo "hg19 reference index is missing: $ref.fai" >&2
  exit 1
fi
if [ ! -f "$meta" ]; then
  echo "SComatic metadata is missing: $meta" >&2
  exit 1
fi

mkdir -p "$out"

o1="$out/Step1_BamCellTypes"
mkdir -p "$o1"
python "$SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py" \
  --bam "$bam" \
  --meta "$meta" \
  --id "$sample" \
  --n_trim 5 \
  --max_nM 5 \
  --max_NH 1 \
  --min_MQ 255 \
  --outdir "$o1"

for split_bam in "$o1"/*.bam; do
  samtools sort -o "${split_bam}.fixed" "$split_bam"
  mv "${split_bam}.fixed" "$split_bam"
  samtools index "$split_bam"
done

o2="$out/Step2_BaseCellCounts"
mkdir -p "$o2"
for split_bam in "$o1"/*.bam; do
  ct="$(basename "$split_bam" | awk -F'.' '{print $(NF-1)}')"
  tmp="$o2/temp_${ct}"
  rm -rf "$tmp"
  mkdir -p "$tmp"
  python "$SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py" \
    --bam "$split_bam" \
    --ref "$ref" \
    --chrom all \
    --out_folder "$o2" \
    --min_bq 30 \
    --min_mq 255 \
    --tmp_dir "$tmp" \
    --nprocs "${NPROCS:-8}"
  rm -rf "$tmp"
done

o3="$out/Step3_BaseCellCountsMerged"
mkdir -p "$o3"
python "$SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py" \
  --tsv_folder "$o2" \
  --outfile "$o3/${sample}.BaseCellCounts.AllCellTypes.tsv"

o4="$out/Step4_VariantCalling"
mkdir -p "$o4"
python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py" \
  --infile "$o3/${sample}.BaseCellCounts.AllCellTypes.tsv" \
  --outfile "$o4/${sample}" \
  --ref "$ref"

python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py" \
  --infile "$o4/${sample}.calling.step1.tsv" \
  --outfile "$o4/${sample}"

awk '$1 ~ /^#/ || $6 == "PASS"' "$o4/${sample}.calling.step2.tsv" \
  > "$o4/${sample}.calling.step2.PASS.tsv"

echo "Done: $o4/${sample}.calling.step2.tsv"
