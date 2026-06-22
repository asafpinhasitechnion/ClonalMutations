# t3_Center — barcode → clone map (T2)

Clone grouping for the `t3_Center` pilot, built from the **published** TCR table
(no raw TCR processing). Clone identity = exact paired α+β clonotype string; cells
sharing it are one clone (the published table gives one clonotype per cell).

## Command
```
python src/clonemut/frontend/build_clone_meta.py \
  --tcr Data/Krishna_ccRCC/ccRCC_TCRs.txt \
  --region _t3_Center \
  --out-prefix results/t3_Center/t3_Center
```

## Outputs
- `t3_Center.clone_meta.tsv` — `Index<TAB>Cell_type` (Cell_type = clone_id). **This is
  the SComatic grouping file** — feed as `--meta` to SplitBam to group reads by clone.
  `Index` = bare 16 bp barcode (matches a STARsolo `CB` tag; add `-1` if the aligner does).
- `t3_Center.clone_dict.tsv` — clone_id ↔ full clonotype string, n_cells, TRB, TRA.
- `t3_Center.clone_sizes.tsv` — size distribution.

## Numbers
- 3,986 cells, **1,436 distinct clones**.
- Coverage by clone size: ≥2 → 355 clones / 2,905 cells; ≥5 → 108 / 2,298; ≥10 → 53 / 1,952; ≥50 → 8 / 1,105.
- Largest clones: 442, 171, 134, 128, 61, 59, 56, 54 cells.

## Validation
- **100% (3,986/3,986)** of clone barcodes appear in the downloaded RNA run
  (`SRR13806069_1`, 5 M-read sample) → grouping keys will match the BAM.

## Scope note for downstream (T5)
The meta currently contains **all** clones (`--min-size 1`). Singletons (1,081 clones
of size 1) add no clonal-pooling power and would spawn ~1,400 tiny per-clone BAMs.
When we split by clone, restrict to expanded clones (`--min-size 2`, or 5 for real
power) — this is the plan's "min_cells restricts to expanded clones; embrace as a
scope statement." Candidate *generation* still uses pooled/clone-agnostic reads.
