# clonemut.explorer — Gate-1 triage UI (Explorer v1)

The exploration half of the First Deliverable. A lightweight, reactive
[marimo](https://marimo.io) app that reads **only the Tier-A Parquet** (via DuckDB) —
never raw reads — and answers the Phase-0 / **Gate-1** question:

> Do groups (clones — or cell-types as a stand-in) carry enough molecules per cell at
> callable sites for the hierarchical arm to have purchase over pseudobulk?

## Views
- **Overview + Gate-1 headline** — % of `(cell, callable-site)` pairs with ≥2 molecules.
- **Group sizes** and **callable sites per group** (top 40).
- **Gate-1 histogram** — molecules per cell per callable site.
- **Candidate universe** — filterable table (min VAF / min alt molecules) + VAF distribution.

## Run (in the `explorer` conda env)
```
conda run -n explorer marimo edit src/clonemut/explorer/explorer_v1.py
```
The prefix dropdown auto-discovers every `*.group_summary.parquet` under `results/`, so
the same app serves the SComatic example today and the real Krishna clone extract once
it lands. Headless validation (runs every cell):
```
conda run -n explorer marimo export html src/clonemut/explorer/explorer_v1.py -o /tmp/x.html
```

## Status / caveats
- v1 reads the **Tier-A** layer (`group_summary`, `candidates`, `coverage_by_site`).
  The per-site molecule→cell→clone evidence panel and germline-calibration view are
  **Tier-B** (v2) and not built yet.
- On the **SComatic example** the Gate-1 number is *not meaningful* (cell-type grouping,
  shallow toy BAM). It becomes meaningful on the real clone extract (Krishna `t3_Center`).
- "Coverage by gene" from the plan is deferred — Tier-A has no gene annotation yet
  (add in Tier-B).

## Environment
Separate `explorer` env (python 3.11 + marimo/duckdb/polars/altair) with **no pysam** —
deliberately decoupled from the py3.7 extraction env (see `src/clonemut/extract/README.md`).
