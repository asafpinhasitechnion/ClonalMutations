# clonemut.extract — molecule-level extraction (Tier A)

The extraction half of the First Deliverable: turn an aligned scRNA-seq BAM into a
**molecule-level** coverage summary plus a permissive **candidate universe**, written
as Parquet for the explorer / modeling arms.

## Why this exists (vs SComatic)
SComatic counts *reads* per cell type. We collapse each `(CB, UB)` UMI to one
**consensus molecule** first, so PCR over-amplification of one erroneous cDNA counts
once, not N times. The molecule is the atomic unit everywhere downstream. Grouping is
generic — the meta maps barcode → group, where group is a `clone_id` (real target) or
a `cell_type` (SComatic stand-in). Identical machinery either way.

## Files
- `molecules.py` — core: pysam pileup → consensus molecules. SComatic-compatible read
  filtering (NH/nM/MQ/BQ/trim, drop secondary/dup/supp, strip CB `-N` suffix), plus a
  required `UB`. Importable + unit-testable. Holds the (v1, swappable) `eps_m` model.
- `extract_tier_a.py` — CLI: aggregates molecules into the Tier-A Parquet data layer
  and the candidate list.

## Run (in the `clonemut` conda env)
```
conda run -n clonemut python src/clonemut/extract/extract_tier_a.py \
  --bam   external/SComatic/example_data/Example.scrnaseq.bam \
  --fasta external/SComatic/example_data/chr10.fa \
  --meta  external/SComatic/example_data/Example.cell_barcode_annotations.tsv \
  --out-prefix results/example_extract/Example --region chr10
```
Outputs `<prefix>.group_summary.parquet`, `<prefix>.candidates.parquet`,
`<prefix>.coverage_by_site.parquet` (skip the last with `--no-site-coverage`).

## Validation
`scripts/validate_extract_example.py` asserts our candidate universe recovers all of
SComatic's `PASS` calls on the bundled example. Current result: **5/5 recovered, 0
missing, 0 extra** — and SComatic's one `Cell_type_noise`-filtered site is correctly
absent (read-level noise collapses below the molecule threshold).

## Environment note (important)
The env is pinned to **python 3.7 + pysam 0.21.0 + htslib 1.17** on purpose: every
py3.10 pysam conda build (0.21–0.24) segfaults / corrupts the heap during pileup
iteration on this WSL2 host (they bundle their own htslib). The py3.7 build links the
external htslib 1.17 and is stable (verified over 2.38M pileup columns). The explorer
reads only Parquet, so its modern stack lives in a separate env.

## Next (Tier B)
Per-candidate × cell × molecule extraction (richer than Tier A) for the per-site
evidence view and the hierarchical arm, including germline calibration sites.
