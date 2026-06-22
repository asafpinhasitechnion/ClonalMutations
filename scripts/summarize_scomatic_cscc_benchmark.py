#!/usr/bin/env python3
"""Summarize the SComatic cSCC WES benchmark tables downloaded from the paper."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "Data" / "SComatic_cSCC" / "paper" / "csv"


def read_csv(name: str) -> list[list[str]]:
    with (PAPER / name).open(newline="") as handle:
        return list(csv.reader(handle))


def print_fig2d() -> None:
    rows = read_csv("Figure2d.csv")
    print("Fig. 2d source-data counts")
    print(",".join(rows[0]))
    for row in rows[1:]:
        print(",".join(row))


def print_table3_scomatic_status() -> None:
    rows = read_csv("Supplementary_table_3.csv")
    header_row = next(i for i, row in enumerate(rows) if row and row[0] == "sample_ID")
    header = rows[header_row]
    index = {name: i for i, name in enumerate(header)}

    counts: Counter[tuple[str, str]] = Counter()
    for row in rows[header_row + 1 :]:
        if len(row) < len(header) or not row or not row[0]:
            continue
        status = row[index["Stat.SComatic"]]
        if status and status not in {"True negative", "NA", "."}:
            counts[(row[index["sample_ID"]], status)] += 1

    statuses = [
        "WES-specific",
        "WES-specific (no alt. reads in WES)",
        "WES-specific (no alt. reads in scRNA-seq)",
        "scRNA-seq & WES",
        "scRNA-seq (other alleles found in WES)",
        "scRNA-seq also found in WES matched normal",
        "scRNA-seq with read support in WES",
        "scRNA-seq-specific",
    ]
    samples = [f"P{i}_cSCC" for i in range(2, 11)]

    print()
    print("Supplementary Table 3 SComatic benchmark-status counts")
    print("sample," + ",".join(statuses))
    for sample in samples:
        print(sample + "," + ",".join(str(counts[(sample, status)]) for status in statuses))


def main() -> None:
    print_fig2d()
    print_table3_scomatic_status()


if __name__ == "__main__":
    main()
