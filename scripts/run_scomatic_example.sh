#!/usr/bin/env bash
# Validate the SComatic base pipeline end-to-end on its bundled example
# (small chr10 region, 5 known SNVs). Establishes a known-good base before we
# build the clone-aware front-end. Run inside the `scomatic` conda env.
set -euo pipefail

SCOMATIC="$(cd "$(dirname "${BASH_SOURCE[0]}")/../external/SComatic" && pwd)"
out="$SCOMATIC/example_data/results"
mkdir -p "$out"

# Reference: unpack + index if needed
[ -f "$SCOMATIC/example_data/chr10.fa" ] || gunzip -k "$SCOMATIC/example_data/chr10.fa.gz"
[ -f "$SCOMATIC/example_data/chr10.fa.fai" ] || samtools faidx "$SCOMATIC/example_data/chr10.fa"
REF="$SCOMATIC/example_data/chr10.fa"
sample=Example

# Step 1 — split BAM by cell type
o1="$out/Step1_BamCellTypes"; mkdir -p "$o1"
python "$SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py" \
  --bam "$SCOMATIC/example_data/Example.scrnaseq.bam" \
  --meta "$SCOMATIC/example_data/Example.cell_barcode_annotations.tsv" \
  --id "$sample" --n_trim 5 --max_nM 5 --max_NH 1 --outdir "$o1"
# SComatic's pysam-written split BAMs can be subtly malformed (bad BGZF EOF /
# index that pysam.pileup cannot read). Rewrite each via samtools sort, then index.
for bam in "$o1"/*.bam; do
  samtools sort -o "${bam}.fixed" "$bam" && mv "${bam}.fixed" "$bam" && samtools index "$bam"
done

# Step 2 — base/cell counts per cell type
o2="$out/Step2_BaseCellCounts"; mkdir -p "$o2"
for bam in "$o1"/*.bam; do
  ct=$(basename "$bam" | awk -F'.' '{print $(NF-1)}')
  tmp="$o2/temp_${ct}"; rm -rf "$tmp"; mkdir -p "$tmp"
  python "$SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py" \
    --bam "$bam" --ref "$REF" --chrom all --out_folder "$o2" \
    --min_bq 30 --tmp_dir "$tmp" --nprocs 8
  rm -rf "$tmp"
done

# Step 3 — merge
o3="$out/Step3_BaseCellCountsMerged"; mkdir -p "$o3"
python "$SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py" \
  --tsv_folder "$o2" --outfile "$o3/${sample}.BaseCellCounts.AllCellTypes.tsv"

# Step 4.1 — beta-binomial calling ; 4.2 — RNA-editing/PoN/clustered filters
o4="$out/Step4_VariantCalling"; mkdir -p "$o4"
python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py" \
  --infile "$o3/${sample}.BaseCellCounts.AllCellTypes.tsv" --outfile "$o4/${sample}" --ref "$REF"
python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py" \
  --infile "$o4/${sample}.calling.step1.tsv" --outfile "$o4/${sample}" \
  --editing "$SCOMATIC/RNAediting/AllEditingSites.hg38.txt" \
  --pon "$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv"

echo "=== PASS calls ==="
awk '$1 ~ /^#/ || $6 == "PASS"' "$o4/${sample}.calling.step2.tsv" | grep -v '^##' | head -20
