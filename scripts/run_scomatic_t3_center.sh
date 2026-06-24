#!/usr/bin/env bash
# Clone-aware SComatic run for the t3_Center ccRCC sample (GRCh38 / GENCODE v32).
#
# Inputs (all already staged on the server):
#   bam   STARsolo CB/UB-tagged alignment (results/t3_Center/star/...)
#   meta  tiered split meta from build_split_meta.py: one group per expanded
#         clone, pooled small clones, singletons -> internal "normal" background
#   ref   GENCODE v32 primary assembly (the exact reference the BAM was aligned to)
#   PoN + RNA-editing: hg38 resources bundled with the SComatic submodule
#
# Unlike the hg19 P9 run, this IS the paper-faithful hg38 path: step2 applies the
# hg38 RNA-editing and PoN filters. Run inside the `scomatic` conda env.
#
# Usage:  NPROCS=16 bash scripts/run_scomatic_t3_center.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCOMATIC="$ROOT/external/SComatic"

sample="t3_Center"
bam="$ROOT/results/t3_Center/star/t3_Center_Aligned.sortedByCoord.out.bam"
meta="$ROOT/results/t3_Center/t3_Center.split_meta.tsv"
ref="$ROOT/reference/gencode_v32/GRCh38.primary_assembly.genome.fa"
out="$ROOT/results/t3_Center/scomatic"
NPROCS="${NPROCS:-8}"

# hg38 filter resources ship gzipped in the submodule; SComatic wants them plain.
editing_gz="$SCOMATIC/RNAediting/AllEditingSites.hg38.txt.gz"
pon_gz="$SCOMATIC/PoNs/PoN.scRNAseq.hg38.tsv.gz"
editing="${editing_gz%.gz}"
pon="${pon_gz%.gz}"
[ -f "$editing" ] || gunzip -k "$editing_gz"
[ -f "$pon" ]     || gunzip -k "$pon_gz"

# ---- input guards: fail early with a clear message ----
for f in "$bam" "$bam.bai" "$ref" "$ref.fai" "$meta"; do
  [ -f "$f" ] || { echo "MISSING input: $f" >&2; exit 1; }
done
mkdir -p "$out"

# Step 1 — split the BAM into one BAM per group (clone tier)
o1="$out/Step1_BamCellTypes"; mkdir -p "$o1"
python "$SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py" \
  --bam "$bam" --meta "$meta" --id "$sample" \
  --n_trim 5 --max_nM 5 --max_NH 1 --min_MQ 255 --outdir "$o1"
# SComatic's pysam-written split BAMs can carry a bad BGZF EOF / index that
# pysam.pileup cannot read. Rewrite each via samtools sort, then index.
for split_bam in "$o1"/*.bam; do
  samtools sort -o "${split_bam}.fixed" "$split_bam"
  mv "${split_bam}.fixed" "$split_bam"
  samtools index "$split_bam"
done

# Step 2 — per-group base/cell counts (the heavy step: one whole-genome pass
# per group; bump NPROCS on a cluster).
o2="$out/Step2_BaseCellCounts"; mkdir -p "$o2"
for split_bam in "$o1"/*.bam; do
  ct="$(basename "$split_bam" | awk -F'.' '{print $(NF-1)}')"
  tmp="$o2/temp_${ct}"; rm -rf "$tmp"; mkdir -p "$tmp"
  python "$SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py" \
    --bam "$split_bam" --ref "$ref" --chrom all --out_folder "$o2" \
    --min_bq 30 --min_mq 255 --tmp_dir "$tmp" --nprocs "$NPROCS"
  rm -rf "$tmp"
done

# Step 3 — merge per-group counts into one matrix
o3="$out/Step3_BaseCellCountsMerged"; mkdir -p "$o3"
python "$SCOMATIC/scripts/MergeCounts/MergeBaseCellCounts.py" \
  --tsv_folder "$o2" --outfile "$o3/${sample}.BaseCellCounts.AllCellTypes.tsv"

# Step 4 — calling: 4.1 beta-binomial; 4.2 RNA-editing + PoN filtering
o4="$out/Step4_VariantCalling"; mkdir -p "$o4"
python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step1.py" \
  --infile "$o3/${sample}.BaseCellCounts.AllCellTypes.tsv" \
  --outfile "$o4/${sample}" --ref "$ref"
python "$SCOMATIC/scripts/BaseCellCalling/BaseCellCalling.step2.py" \
  --infile "$o4/${sample}.calling.step1.tsv" --outfile "$o4/${sample}" \
  --editing "$editing" --pon "$pon"

awk '$1 ~ /^#/ || $6 == "PASS"' "$o4/${sample}.calling.step2.tsv" \
  > "$o4/${sample}.calling.step2.PASS.tsv"

echo "Done: $o4/${sample}.calling.step2.tsv"
echo "PASS calls -> $o4/${sample}.calling.step2.PASS.tsv"
