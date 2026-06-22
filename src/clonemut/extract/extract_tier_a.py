#!/usr/bin/env python3
"""Tier-A extraction: BAM + grouping meta -> coverage summary + candidate universe.

This is the *cheap, genome-wide* pass of the First Deliverable. It produces the
data layer the explorer's Gate-1 triage reads, plus the permissive candidate list
that both modeling arms must later adjudicate over the *same* universe.

Outputs (Parquet, written next to --out-prefix):
  <prefix>.group_summary.parquet     per group (clone/cell-type):
      group_id, n_cells, n_callable_sites, mean_cov,
      median_mol_per_cell_per_site, mol_per_cell_hist (list)
  <prefix>.coverage_by_site.parquet  per site x group (omit with --no-site-coverage):
      site_id, chrom, pos, group_id, depth, n_cells, n_molecules
  <prefix>.candidates.parquet        permissive candidate universe (pooled):
      site_id, chrom, pos, ref, alt, depth, n_molecules, n_alt_molecules,
      n_cells, n_alt_cells, vaf

Candidate rule (permissive, by design — NOT a significance call): a site enters the
universe if pooled molecule depth >= min_dp AND alt molecules >= min_ac_mol. Keeping
this loose is deliberate: the pseudobulk and hierarchical callers must both be able
to recover the full universe, otherwise the arm comparison is rigged.

Run in the `clonemut` conda env.
  python src/clonemut/extract/extract_tier_a.py \
      --bam BAM --fasta REF.fa --meta GROUPS.tsv --out-prefix OUT [--region chr10]
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import Counter, defaultdict

import pyarrow as pa
import pyarrow.parquet as pq

# allow running as a plain script (no install): import the sibling module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from molecules import Filters, iter_sites, load_meta  # noqa: E402

_ACGT = frozenset("ACGT")


def write_parquet(rows, path, schema):
    """Write a list of row dicts to Parquet under an explicit schema.

    The explicit schema keeps column types stable even when ``rows`` is empty
    (e.g. a region with no candidates), so downstream readers never see a
    schema-less file.
    """
    if rows:
        table = pa.Table.from_pylist(rows, schema=schema)
    else:
        table = schema.empty_table()
    pq.write_table(table, path)


class GroupAcc:
    """Per-group Tier-A accumulator."""
    __slots__ = ("cells", "n_callable_sites", "sum_depth", "mpc_hist")

    def __init__(self):
        self.cells: set[str] = set()
        self.n_callable_sites = 0
        self.sum_depth = 0
        self.mpc_hist: Counter[int] = Counter()  # molecules-per-(cell, callable site)


def median_from_hist(hist: Counter[int]) -> float:
    total = sum(hist.values())
    if total == 0:
        return 0.0
    items = sorted(hist.items())
    # lower & upper median positions (1-based)
    lo_i, hi_i = (total + 1) // 2, (total + 2) // 2
    lo = hi = None
    cum = 0
    for val, cnt in items:
        cum += cnt
        if lo is None and cum >= lo_i:
            lo = val
        if cum >= hi_i:
            hi = val
            break
    return (lo + hi) / 2.0


def parse_region(region: str | None, bam_path: str):
    """Yield (chrom, start0, end) intervals. None -> every reference in the BAM."""
    import pysam
    bam = pysam.AlignmentFile(bam_path)
    if not region:
        for chrom in bam.references:
            yield chrom, None, None
        return
    if ":" in region:
        chrom, span = region.split(":", 1)
        a, b = span.replace(",", "").split("-")
        yield chrom, int(a) - 1, int(b)   # 1-based inclusive -> 0-based half-open
    else:
        yield region, None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True, help="reference FASTA (indexed .fai)")
    ap.add_argument("--meta", required=True, help="Index<TAB>group grouping file")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--region", default=None,
                    help="chrom or chrom:start-end (1-based); default = all references")
    ap.add_argument("--min-dp", type=int, default=5,
                    help="min molecule depth for a callable/candidate site (default 5)")
    ap.add_argument("--min-ac-mol", type=int, default=2,
                    help="min alt molecules for a candidate (default 2)")
    ap.add_argument("--no-site-coverage", action="store_true",
                    help="skip the per-site x group coverage table (large genome-wide)")
    # baked filters (SComatic-compatible defaults)
    ap.add_argument("--max-nh", type=int, default=1)
    ap.add_argument("--max-nm", type=int, default=5)
    ap.add_argument("--min-mq", type=int, default=255)
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--n-trim", type=int, default=5)
    a = ap.parse_args()

    meta = load_meta(a.meta)
    if not meta:
        sys.exit("ERROR: empty meta")
    filt = Filters(max_nh=a.max_nh, max_nm=a.max_nm, min_mq=a.min_mq,
                   min_bq=a.min_bq, n_trim=a.n_trim)

    groups: dict[str, GroupAcc] = defaultdict(GroupAcc)
    site_rows: list[dict] = []
    cand_rows: list[dict] = []
    n_sites = 0

    for chrom, start0, end in parse_region(a.region, a.bam):
        for pos1, ref, molecules in iter_sites(a.bam, a.fasta, chrom, start0, end,
                                               meta, filt):
            n_sites += 1
            site_id = f"{chrom}:{pos1}"

            # bucket molecules by group, and by cell within group
            by_group: dict[str, list] = defaultdict(list)
            for m in molecules:
                by_group[m.group].append(m)

            for g, mols in by_group.items():
                acc = groups[g]
                cells_here: Counter[str] = Counter(m.cb for m in mols)
                acc.cells.update(cells_here)
                depth = len(mols)
                if depth >= a.min_dp:
                    acc.n_callable_sites += 1
                    acc.sum_depth += depth
                    for c, n_mol in cells_here.items():
                        acc.mpc_hist[n_mol] += 1
                    if not a.no_site_coverage:
                        site_rows.append(dict(
                            site_id=site_id, chrom=chrom, pos=pos1, group_id=g,
                            depth=depth, n_cells=len(cells_here), n_molecules=depth))

            # pooled candidate test
            depth = len(molecules)
            alt_counts = Counter(m.base for m in molecules
                                 if m.base in _ACGT and m.base != ref)
            if depth >= a.min_dp and alt_counts:
                alt, n_alt_mol = alt_counts.most_common(1)[0]
                if n_alt_mol >= a.min_ac_mol:
                    alt_cells = {m.cb for m in molecules if m.base == alt}
                    all_cells = {m.cb for m in molecules}
                    cand_rows.append(dict(
                        site_id=site_id, chrom=chrom, pos=pos1, ref=ref, alt=alt,
                        depth=depth, n_molecules=depth, n_alt_molecules=n_alt_mol,
                        n_cells=len(all_cells), n_alt_cells=len(alt_cells),
                        vaf=round(n_alt_mol / depth, 4)))

    # ---- write group_summary ----
    gs_rows = []
    for g, acc in sorted(groups.items()):
        mean_cov = acc.sum_depth / acc.n_callable_sites if acc.n_callable_sites else 0.0
        hist = [int(acc.mpc_hist.get(k, 0)) for k in range(1, (max(acc.mpc_hist) + 1) if acc.mpc_hist else 1)]
        gs_rows.append(dict(
            group_id=g, n_cells=len(acc.cells), n_callable_sites=acc.n_callable_sites,
            mean_cov=round(mean_cov, 3),
            median_mol_per_cell_per_site=median_from_hist(acc.mpc_hist),
            mol_per_cell_hist=hist))

    gs_schema = pa.schema([
        ("group_id", pa.string()), ("n_cells", pa.int64()),
        ("n_callable_sites", pa.int64()), ("mean_cov", pa.float64()),
        ("median_mol_per_cell_per_site", pa.float64()),
        ("mol_per_cell_hist", pa.list_(pa.int64()))])
    cand_schema = pa.schema([
        ("site_id", pa.string()), ("chrom", pa.string()), ("pos", pa.int64()),
        ("ref", pa.string()), ("alt", pa.string()), ("depth", pa.int64()),
        ("n_molecules", pa.int64()), ("n_alt_molecules", pa.int64()),
        ("n_cells", pa.int64()), ("n_alt_cells", pa.int64()), ("vaf", pa.float64())])
    site_schema = pa.schema([
        ("site_id", pa.string()), ("chrom", pa.string()), ("pos", pa.int64()),
        ("group_id", pa.string()), ("depth", pa.int64()),
        ("n_cells", pa.int64()), ("n_molecules", pa.int64())])

    write_parquet(gs_rows, f"{a.out_prefix}.group_summary.parquet", gs_schema)
    write_parquet(cand_rows, f"{a.out_prefix}.candidates.parquet", cand_schema)
    if not a.no_site_coverage:
        write_parquet(site_rows, f"{a.out_prefix}.coverage_by_site.parquet", site_schema)

    # ---- legible report ----
    n_cand = len(cand_rows)
    sys.stderr.write(
        f"sites_covered={n_sites}  groups={len(groups)}  candidates={n_cand}\n"
        f"group_summary -> {a.out_prefix}.group_summary.parquet\n"
        f"candidates    -> {a.out_prefix}.candidates.parquet\n")
    if not a.no_site_coverage:
        sys.stderr.write(f"coverage_by_site -> {a.out_prefix}.coverage_by_site.parquet\n")
    # Gate-1 headline: do cells carry >=2 molecules at callable sites?
    if gs_rows:
        meds = [r["median_mol_per_cell_per_site"] for r in gs_rows]
        sys.stderr.write(
            "GATE-1 (molecules/cell/callable-site) median across groups: "
            f"min={min(meds)} median={sorted(meds)[len(meds)//2]} max={max(meds)}\n")


if __name__ == "__main__":
    main()
