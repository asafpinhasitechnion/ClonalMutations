#!/usr/bin/env python3
"""Build a tiered SComatic split/group meta from a clone map.

`build_clone_meta.py --min-size N` *drops* clones below N, which would discard
most cells. For somatic calling we instead want to KEEP every cell but group it
by clone size, so that:

  * expanded clones (>= --expanded-min cells) each get their OWN group -> these
    are where we generate/call candidates with real per-clone power;
  * small clones (--small-min .. expanded-min-1) are pooled into one group -> real
    but underpowered clonal signal, kept distinct for visualization;
  * singletons (< --small-min, i.e. size 1 by default) are pooled into a "normal"
    group -> an internal background/pseudo-PoN: a true somatic mutation in an
    expanded clone should be absent across these unrelated cells.

This loses no cells (every barcode in the input meta is emitted) while keeping
the per-group BAM count manageable (~110 instead of ~1,400).

Inputs
------
  --clone-meta   <prefix>.clone_meta.tsv   Index<TAB>Cell_type (Cell_type=clone_id)
  --clone-sizes  <prefix>.clone_sizes.tsv  clone_id<TAB>n_cells

Output
------
  --out  split meta, Index<TAB>Cell_type, where Cell_type is the clone_id for
         expanded clones, --small-label for small clones, --normal-label for
         singletons. This is the file you pass to SComatic SplitBamCellTypes
         as --meta.

Notes
-----
* STARsolo writes bare 16 bp CB tags (no `-1`); the clone meta Index is already
  bare, so they match directly. Use --cb-suffix to append one if a different
  aligner's CB tags carry it.
"""
import argparse
import sys
from collections import Counter


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clone-meta", required=True, help="<prefix>.clone_meta.tsv (Index, Cell_type=clone_id)")
    p.add_argument("--clone-sizes", required=True, help="<prefix>.clone_sizes.tsv (clone_id, n_cells)")
    p.add_argument("--out", required=True, help="output split meta path")
    p.add_argument("--expanded-min", type=int, default=5,
                   help="min cells for a clone to get its own group (default 5)")
    p.add_argument("--small-min", type=int, default=2,
                   help="min cells to count as a 'small' clone vs a singleton (default 2)")
    p.add_argument("--small-label", default="small_clones", help="group name for small clones")
    p.add_argument("--normal-label", default="normal", help="group name for singletons (internal background)")
    p.add_argument("--cb-suffix", default="", help="append to each Index to match BAM CB tags (e.g. -1)")
    return p.parse_args()


def main():
    a = parse_args()
    if a.small_min > a.expanded_min:
        sys.exit("ERROR: --small-min must be <= --expanded-min")

    # 1. clone_id -> size
    size = {}
    with open(a.clone_sizes) as fh:
        header = fh.readline()
        for line in fh:
            cid, n = line.rstrip("\n").split("\t")
            size[cid] = int(n)

    # 2. relabel each cell by its clone's tier
    n_in = 0
    tier_groups = Counter()      # group label -> n cells
    tier_clones = {"expanded": set(), "small": set(), "normal": set()}
    rows = []
    with open(a.clone_meta) as fh:
        header = fh.readline()
        for line in fh:
            bc, cid = line.rstrip("\n").split("\t")
            n_in += 1
            n = size.get(cid)
            if n is None:
                sys.exit(f"ERROR: clone {cid} in meta but absent from sizes table")
            if n >= a.expanded_min:
                group, tier = cid, "expanded"
            elif n >= a.small_min:
                group, tier = a.small_label, "small"
            else:
                group, tier = a.normal_label, "normal"
            rows.append((bc + a.cb_suffix, group))
            tier_groups[group] += 1
            tier_clones[tier].add(cid)

    # 3. write split meta
    with open(a.out, "w") as fh:
        fh.write("Index\tCell_type\n")
        for bc, group in sorted(rows):
            fh.write(f"{bc}\t{group}\n")

    # 4. report
    n_groups = len(tier_groups)
    sys.stderr.write(
        f"cells in: {n_in}  cells out: {len(rows)}  (no cell dropped)\n"
        f"groups written: {n_groups}\n"
        f"  expanded (>= {a.expanded_min}): {len(tier_clones['expanded'])} clones, "
        f"{sum(v for g, v in tier_groups.items() if g not in (a.small_label, a.normal_label))} cells\n"
        f"  small ({a.small_min}..{a.expanded_min - 1}) -> '{a.small_label}': "
        f"{len(tier_clones['small'])} clones, {tier_groups.get(a.small_label, 0)} cells\n"
        f"  singletons (< {a.small_min}) -> '{a.normal_label}': "
        f"{len(tier_clones['normal'])} clones, {tier_groups.get(a.normal_label, 0)} cells\n"
        f"wrote {a.out}\n"
    )


if __name__ == "__main__":
    main()
