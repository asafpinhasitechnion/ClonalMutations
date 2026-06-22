#!/usr/bin/env python3
"""Molecule-level pileup core for the clone-aware extractor.

This is the project's reusable substrate: it turns an aligned scRNA-seq BAM into
**consensus molecules** at each genomic position, grouping reads by ``(CB, UB)``.

Why molecules and not reads (the central adaptation vs SComatic):
  SComatic's BaseCellCounter counts *reads* per cell-type. Three PCR copies of one
  over-amplified erroneous cDNA look like three votes. The independent unit is the
  *molecule* (one UMI), so we collapse a UMI's reads to a single consensus base and
  summarise read support into a per-molecule error ``eps_m``. Everything downstream
  counts molecules.

Read/base filtering mirrors SComatic (BaseCellCounter): drop secondary / duplicate /
supplementary reads; require unique mapping (``NH<=max_nh``, default 1 — critical in
T cells); ``nM<=max_nm``; mapping quality ``>=min_mq``; base quality ``>=min_bq``;
trim ``n_trim`` bases from each read end; strip the ``-N`` suffix from the CB tag.
We additionally require a ``UB`` tag (a read without a UMI cannot form a molecule).

The grouping is generic: the meta maps barcode -> group, where group is a clone_id
(real target) or a cell_type (SComatic stand-in). Identical machinery either way.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import reduce

import pysam

_ACGT = frozenset("ACGT")


@dataclass(frozen=True)
class Filters:
    """Baked read/base filters (SComatic-compatible defaults)."""
    max_nh: int = 1        # unique mappers only (NH tag)
    max_nm: int = 5        # max mismatches per read (nM tag)
    min_mq: int = 255      # STAR reports 255 for unique alignments
    min_bq: int = 20       # min base quality at the queried position
    n_trim: int = 5        # ignore this many bases at each read end


@dataclass(frozen=True)
class Molecule:
    """One UMI collapsed to a consensus base at a single site."""
    cb: str            # cell barcode (suffix stripped)
    group: str         # clone_id / cell_type
    base: str          # consensus base, one of A/C/G/T (or 'N' if unresolved)
    n_reads: int       # reads supporting this molecule at this site
    agree_frac: float  # fraction of reads agreeing with the consensus base
    eps_m: float       # P(consensus base is wrong); read-count + agreement summary


def consensus(calls: list[tuple[str, int]]) -> tuple[str, int, float, float]:
    """Collapse a UMI's (base, base_quality) read calls into a consensus.

    Returns (base, n_reads, agree_frac, eps_m).

    eps_m model (v1, deliberately simple and swappable — it is *stored* now but only
    *consumed* by the hierarchical arm / per-molecule evidence views later, so the
    exact form is not yet load-bearing):
      - supporting reads contribute their independent error prob 10^(-bq/10);
        more agreeing reads -> lower eps (product),
      - any dissenting reads inflate eps by their fraction.
    Clamped to [1e-6, 0.5].  TODO(modeling): replace with a proper consensus error.
    """
    n = len(calls)
    by_base: Counter[str] = Counter()
    qsum: dict[str, int] = defaultdict(int)
    for base, bq in calls:
        by_base[base] += 1
        qsum[base] += bq
    # consensus = most reads; tie broken by summed base quality
    top = max(by_base.items(), key=lambda kv: (kv[1], qsum[kv[0]]))[0]
    n_agree = by_base[top]
    agree_frac = n_agree / n

    e_support = [10.0 ** (-bq / 10.0) for base, bq in calls if base == top]
    e_combined = reduce(lambda x, y: x * y, e_support, 1.0) if e_support else 0.5
    n_dissent = n - n_agree
    eps_m = min(0.5, max(1e-6, e_combined + n_dissent / n))
    return top, n, agree_frac, eps_m


def _tag(aln: pysam.AlignedSegment, name: str):
    """Read an optional tag, tolerating pysam's occasional get_tag failures."""
    try:
        return aln.get_tag(name)
    except (KeyError, TypeError):
        return None


def _passes_read(aln: pysam.AlignedSegment, f: Filters) -> bool:
    if aln.is_secondary or aln.is_duplicate or aln.is_supplementary:
        return False
    nh = _tag(aln, "NH")
    if nh is not None and nh > f.max_nh:
        return False
    nm = _tag(aln, "nM")
    if nm is not None and nm > f.max_nm:
        return False
    return True


def iter_sites(bam_path: str, fasta_path: str, chrom: str,
               start: int | None, end: int | None,
               meta: dict[str, str], f: Filters, max_depth: int = 1_000_000):
    """Yield ``(pos1, ref_base, [Molecule, ...])`` for each covered site in a region.

    ``pos1`` is 1-based. Sites with reference base N are skipped. Molecules whose CB
    is not in ``meta`` are dropped (we only assess annotated cells). ``max_depth`` is
    pysam's pileup depth cap (keep it bounded — multi-billion values segfault pysam).
    """
    bam = pysam.AlignmentFile(bam_path)
    fasta = pysam.FastaFile(fasta_path)
    if start is None:
        start = 0
    if end is None:
        end = bam.get_reference_length(chrom)

    pileup = bam.pileup(
        chrom, start, end,
        truncate=True,                 # only columns within [start, end)
        min_base_quality=0,            # we filter bq ourselves for full control
        min_mapping_quality=f.min_mq,
        ignore_overlaps=False,
        max_depth=max_depth,
    )
    for col in pileup:
        pos0 = col.reference_pos
        ref = fasta.fetch(chrom, pos0, pos0 + 1).upper()
        if ref == "N" or ref == "":
            continue

        # gather read calls per (cb, ub) -> molecule
        mol_calls: dict[tuple[str, str], list[tuple[str, int]]] = defaultdict(list)
        mol_group: dict[tuple[str, str], str] = {}
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            aln = pr.alignment
            if not _passes_read(aln, f):
                continue
            qpos = pr.query_position
            qlen = aln.query_length
            if qpos < f.n_trim or qpos >= qlen - f.n_trim:
                continue
            bq = aln.query_qualities[qpos]
            if bq < f.min_bq:
                continue
            base = aln.query_sequence[qpos].upper()
            if base not in _ACGT:
                continue
            cb = _tag(aln, "CB")
            ub = _tag(aln, "UB")
            if not cb or not ub:
                continue
            cb = cb.split("-")[0]
            group = meta.get(cb)
            if group is None:
                continue
            key = (cb, ub)
            mol_calls[key].append((base, bq))
            mol_group[key] = group

        if not mol_calls:
            continue
        molecules = []
        for (cb, ub), calls in mol_calls.items():
            base, n_reads, agree_frac, eps_m = consensus(calls)
            molecules.append(Molecule(cb, mol_group[(cb, ub)], base,
                                      n_reads, agree_frac, eps_m))
        yield pos0 + 1, ref, molecules


def load_meta(path: str) -> dict[str, str]:
    """Read an ``Index<TAB>Cell_type`` grouping file -> {barcode: group}.

    Barcode is stripped of any ``-N`` suffix to match the BAM CB tag handling.
    """
    meta: dict[str, str] = {}
    with open(path) as fh:
        header = fh.readline()  # skip header (Index<TAB>Cell_type)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            meta[parts[0].split("-")[0]] = parts[1]
    return meta
