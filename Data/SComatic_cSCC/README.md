# SComatic cSCC Validation Data

This folder tracks public cSCC data used by SComatic for matched scRNA/WES validation.

Initial fast example:

| Sample | GEO | SRA run | Label | Submitted BAM size |
|---|---|---|---|---:|
| P3 tumor | GSM4284229 | SRR11832842 | P3_cSCC_1_scRNA | 3.63 GB |
| P3 normal | GSM4284230 | SRR11832843 | P3_normal_scRNA | 4.20 GB |

Metadata downloaded:

- `GSE144236_patient_metadata_new.txt.gz`: corrected cell metadata.
- `GSE144236_family.soft.gz`: GEO sample labels and SRA experiment mapping.
- `SRP244706_RunInfo.csv`: SRA scRNA run table.
- `P3_scRNA_ENA_submitted.tsv`: submitted BAM/BAI locations and checksums.
- `GSE144237_WES_RunInfo.csv`: matched WES run table, for later validation.

Run the tumor-only submitted BAM fast example:

```bash
bash Data/SComatic_cSCC/scripts/download_p3_bam.sh
```

The SRA-lite download path is retained as a fallback:

```bash
bash Data/SComatic_cSCC/scripts/download_p3_sra.sh
```

Resume/download the matched normal scRNA later:

```bash
INCLUDE_NORMAL=1 bash Data/SComatic_cSCC/scripts/download_p3_sra.sh
```
