# SComatic cSCC Validation Data

This folder tracks public cSCC data used by SComatic for matched scRNA/WES validation.

Current BAM-first comparison target:

| Sample | GEO | SRA run | Label | Submitted BAM size |
|---|---|---|---|---:|
| P9 tumor | GSM4284244 | SRR11832857 | P9_cSCC_scRNA | 30.91 GB |

Metadata downloaded:

- `GSE144236_patient_metadata_new.txt.gz`: corrected cell metadata.
- `GSE144236_family.soft.gz`: GEO sample labels and SRA experiment mapping.
- `SRP244706_RunInfo.csv`: SRA scRNA run table.
- `P9_scRNA_ENA_submitted.tsv`: submitted BAM/BAI locations and checksums.
- `GSE144237_WES_RunInfo.csv`: matched WES run table, for later validation.
- `paper/`: SComatic paper source tables and benchmark-count extraction.

Run the P9 submitted BAM download:

```bash
bash Data/SComatic_cSCC/scripts/download_p9_bam.sh
```

Scope note: P3 sample payload data was removed because P3 has no non-true-negative
SComatic entries in the paper's cSCC benchmark table. P9 is a better one-sample
BAM-first target: the paper table reports 43 `scRNA-seq & WES` SComatic calls and
4 `scRNA-seq-specific` calls for `P9_cSCC`.
