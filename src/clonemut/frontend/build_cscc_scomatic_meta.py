#!/usr/bin/env python3
"""Build SComatic cell-type metadata from the cSCC patient metadata table."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create Index/Cell_type metadata for SComatic from GSE144236 cSCC annotations."
    )
    parser.add_argument("--metadata", required=True, help="GSE144236_patient_metadata_new.txt.gz")
    parser.add_argument("--patient", required=True, help="patient id, for example P9")
    parser.add_argument("--tum-norm", default="Tumor", choices=["Tumor", "Normal"])
    parser.add_argument(
        "--celltype-level",
        default="level1_celltype",
        choices=["level1_celltype", "level2_celltype", "level3_celltype"],
    )
    parser.add_argument(
        "--cb-suffix",
        default="-1",
        help="suffix to append to bare 10x barcodes so they match BAM CB tags",
    )
    parser.add_argument("--out-prefix", required=True)
    return parser.parse_args()


def clean(value: str) -> str:
    return value.strip().strip('"')


def main() -> None:
    args = parse_args()
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(args.metadata, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = [clean(x) for x in next(reader)]
        columns = ["cell_id"] + header
        level_idx = columns.index(args.celltype_level)
        patient_idx = columns.index("patient")
        tum_norm_idx = columns.index("tum.norm")

        rows: list[tuple[str, str, str]] = []
        for raw in reader:
            if not raw:
                continue
            row = [clean(x) for x in raw]
            if len(row) != len(columns):
                raise ValueError(f"Unexpected field count: expected {len(columns)}, got {len(row)}")
            if row[patient_idx] != args.patient or row[tum_norm_idx] != args.tum_norm:
                continue

            cell_id = row[0]
            parts = cell_id.split("_")
            if len(parts) < 3:
                raise ValueError(f"Cannot parse cell id: {cell_id}")
            barcode = parts[-1] + args.cb_suffix
            cell_type = row[level_idx].replace(" ", "_")
            rows.append((barcode, cell_type, cell_id))

    if not rows:
        raise SystemExit(f"No rows matched patient={args.patient} tum.norm={args.tum_norm}")

    meta_path = out_prefix.with_suffix(".scomatic_meta.tsv")
    full_path = out_prefix.with_suffix(".metadata.tsv")
    counts_path = out_prefix.with_suffix(".celltype_counts.tsv")

    with meta_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Index", "Cell_type"])
        for barcode, cell_type, _ in rows:
            writer.writerow([barcode, cell_type])

    with full_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Index", "Cell_type", "source_cell_id"])
        for row in rows:
            writer.writerow(row)

    counts = Counter(cell_type for _, cell_type, _ in rows)
    with counts_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Cell_type", "n_cells"])
        for cell_type, n_cells in counts.most_common():
            writer.writerow([cell_type, n_cells])

    print(f"wrote {meta_path} ({len(rows)} cells)")
    print(f"wrote {counts_path} ({len(counts)} cell types)")


if __name__ == "__main__":
    main()
