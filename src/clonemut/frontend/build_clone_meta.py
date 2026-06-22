#!/usr/bin/env python3
"""Build a barcode -> clone map from a published TCR clonotype table.

Produces SComatic-shaped grouping metadata where the "Cell_type" column is a
clean clone_id, so SComatic's SplitBam/BaseCellCounter group reads BY CLONE
instead of by cell type.

Inputs
------
A TCR table (TSV) with a cell-barcode column and a clonotype-string column.
For Krishna ccRCC: `cell` = e.g. AAACCTGAGGGCTTGA-1_t3_Center, `TCR_clone` =
e.g. TRBV24-1:TRBJ2-7:CATSE..._TRAV26-1:TRAJ49:CIATG... (paired beta_alpha).

The clone identity is the *exact clonotype string* (paired alpha+beta); cells
sharing it are one clone. The published table gives one clonotype per cell.

Outputs (written next to --out_prefix)
--------------------------------------
  <prefix>.clone_meta.tsv   Index<TAB>Cell_type  (Cell_type = clone_id) -> SComatic
  <prefix>.clone_dict.tsv   clone_id, tcr_clone, n_cells, trb, tra
  <prefix>.clone_sizes.tsv  clone_id, n_cells (descending)

Notes
-----
* clone_id is `clone_<rank>` (1 = largest), deterministic by (size desc, string).
* Index is the bare cell barcode (16 bp by default), to match the BAM CB tag.
  The exact CB format (bare 16bp vs with a `-1` suffix) depends on the aligner;
  use --keep-suffix / --bc-regex if the BAM tags differ.
* --min-size restricts the META (per-clone split) to expanded clones; singletons
  add no clonal-pooling power and would create thousands of tiny BAMs. The dict
  and sizes tables always list ALL clones for transparency.
"""
import argparse
import csv
import re
import sys
from collections import Counter


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tcr", required=True, help="TCR clonotype TSV")
    p.add_argument("--out-prefix", required=True, help="output path prefix")
    p.add_argument("--region", default=None,
                   help="keep only cells whose barcode contains this token (e.g. _t3_Center)")
    p.add_argument("--cell-col", default="cell")
    p.add_argument("--clone-col", default="TCR_clone")
    p.add_argument("--bc-regex", default=r"^([ACGTN]+)",
                   help="regex whose group(1) extracts the BAM-matching barcode from the cell id")
    p.add_argument("--min-size", type=int, default=1,
                   help="min cells per clone to include in clone_meta (default 1 = all)")
    return p.parse_args()


def main():
    a = parse_args()
    bc_re = re.compile(a.bc_regex)

    # 1. read table, optional region filter, extract (barcode, clonotype)
    cell_to_clone = {}
    n_read = n_region = n_badbc = 0
    with open(a.tcr) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            n_read += 1
            cell = row[a.cell_col]
            if a.region and a.region not in cell:
                continue
            n_region += 1
            m = bc_re.match(cell)
            if not m:
                n_badbc += 1
                continue
            bc = m.group(1)
            clonotype = row[a.clone_col].strip()
            # published table = one clonotype per cell; guard anyway
            if bc in cell_to_clone and cell_to_clone[bc] != clonotype:
                sys.stderr.write(f"WARN: barcode {bc} has two clonotypes; keeping first\n")
                continue
            cell_to_clone[bc] = clonotype

    if not cell_to_clone:
        sys.exit("ERROR: no cells matched (check --region / --cell-col)")

    # 2. clone sizes + deterministic clone_id (size desc, then string)
    sizes = Counter(cell_to_clone.values())
    ranked = sorted(sizes.items(), key=lambda kv: (-kv[1], kv[0]))
    width = max(4, len(str(len(ranked))))
    clone_id = {ct: f"clone_{i+1:0{width}d}" for i, (ct, _) in enumerate(ranked)}

    def split_ab(ct):
        # "TRB...:..._TRA...:..." -> (beta, alpha); robust to missing '_'
        parts = ct.split("_", 1)
        return (parts[0], parts[1] if len(parts) > 1 else "")

    # 3. write dict + sizes (ALL clones)
    with open(f"{a.out_prefix}.clone_dict.tsv", "w") as fh:
        fh.write("clone_id\ttcr_clone\tn_cells\ttrb\ttra\n")
        for ct, n in ranked:
            trb, tra = split_ab(ct)
            fh.write(f"{clone_id[ct]}\t{ct}\t{n}\t{trb}\t{tra}\n")
    with open(f"{a.out_prefix}.clone_sizes.tsv", "w") as fh:
        fh.write("clone_id\tn_cells\n")
        for ct, n in ranked:
            fh.write(f"{clone_id[ct]}\t{n}\n")

    # 4. write SComatic meta (Index<TAB>Cell_type), restricted to expanded clones
    n_meta = 0
    with open(f"{a.out_prefix}.clone_meta.tsv", "w") as fh:
        fh.write("Index\tCell_type\n")
        for bc, ct in sorted(cell_to_clone.items()):
            if sizes[ct] < a.min_size:
                continue
            fh.write(f"{bc}\t{clone_id[ct]}\n")
            n_meta += 1

    # 5. report
    def coverage(thr):
        cl = [n for n in sizes.values() if n >= thr]
        return len(cl), sum(cl)
    sys.stderr.write(
        f"rows={n_read} region_cells={n_region} bad_barcodes={n_badbc}\n"
        f"cells_with_clone={len(cell_to_clone)} distinct_clones={len(sizes)}\n"
        f"meta written (min_size={a.min_size}): {n_meta} cells\n"
        "clone/cell coverage by size threshold:\n"
    )
    for thr in (1, 2, 5, 10, 50):
        nc, ncell = coverage(thr)
        sys.stderr.write(f"  size>= {thr:<3d}: {nc:5d} clones, {ncell:5d} cells\n")
    sys.stderr.write(f"largest clones: {[n for _, n in ranked[:8]]}\n")


if __name__ == "__main__":
    main()
