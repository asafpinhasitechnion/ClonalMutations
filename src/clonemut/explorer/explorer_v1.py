#!/usr/bin/env python3
"""Explorer v1 — Gate-1 triage over the Tier-A molecule extraction.

A lightweight, reactive instrument (the plan's "simpler app"): it reads ONLY the
Tier-A Parquet via DuckDB — never raw reads — and answers the Phase-0 / Gate-1
question: do groups (clones, or cell-types as a stand-in) carry enough molecules per
cell at callable sites for the hierarchical arm to have purchase over pseudobulk?

Views:
  - overview + the headline Gate-1 number (% of (cell, callable-site) pairs with >=2
    molecules),
  - group-size distribution, callable-sites-per-group,
  - molecules-per-cell-per-callable-site histogram (the Gate-1 view),
  - filterable candidate table + candidate VAF distribution.

Run in the `explorer` conda env:
    conda run -n explorer marimo edit src/clonemut/explorer/explorer_v1.py
    # or headless export for a quick check:
    conda run -n explorer marimo export html src/clonemut/explorer/explorer_v1.py -o /tmp/x.html
"""
import marimo

app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import duckdb
    import polars as pl
    import altair as alt
    from pathlib import Path
    return Path, alt, duckdb, mo, pl


@app.cell
def _(Path, mo):
    # Discover extraction prefixes (anything with a *.group_summary.parquet).
    found = sorted(str(p)[: -len(".group_summary.parquet")]
                   for p in Path("results").rglob("*.group_summary.parquet"))
    options = found or ["results/example_extract/Example"]
    prefix = mo.ui.dropdown(options=options, value=options[0], label="Extraction prefix")
    mo.md("# Clone-mutation explorer — v1 (Tier A)\n"
          "Molecule-level Gate-1 triage. Pick an extraction below.")
    return options, prefix


@app.cell
def _(prefix):
    prefix
    return


@app.cell
def _(duckdb, prefix):
    con = duckdb.connect()
    pfx = prefix.value
    gs = con.execute(
        f"SELECT * FROM read_parquet('{pfx}.group_summary.parquet')"
    ).pl()
    cand = con.execute(
        f"SELECT * FROM read_parquet('{pfx}.candidates.parquet')"
    ).pl()
    # molecules-per-cell histogram, unnested from the per-group list column
    mpc = con.execute(f"""
        SELECT group_id, idx AS n_molecules, freq
        FROM (
          SELECT group_id, UNNEST(mol_per_cell_hist) AS freq,
                 generate_subscripts(mol_per_cell_hist, 1) AS idx
          FROM read_parquet('{pfx}.group_summary.parquet')
        ) WHERE freq > 0
    """).pl()
    return cand, con, gs, mpc, pfx


@app.cell
def _(gs, mo, mpc):
    total_pairs = int(mpc["freq"].sum()) if mpc.height else 0
    ge2 = int(mpc.filter(mpc["n_molecules"] >= 2)["freq"].sum()) if mpc.height else 0
    pct_ge2 = (100.0 * ge2 / total_pairs) if total_pairs else 0.0
    gate = "✅ hierarchy has purchase" if pct_ge2 >= 10 else "⚠️ thin — pseudobulk/targeted likely"
    mo.md(f"""
    ## Overview
    - **groups:** {gs.height} &nbsp;|&nbsp; **cells (max group):** {int(gs['n_cells'].max()) if gs.height else 0}
      &nbsp;|&nbsp; **callable sites (total):** {int(gs['n_callable_sites'].sum()) if gs.height else 0}
    - **Gate-1:** {pct_ge2:.1f}% of (cell, callable-site) pairs carry **≥2 molecules**
      &nbsp; → &nbsp; {gate}
    """)
    return (pct_ge2,)


@app.cell
def _(alt, gs, mo):
    mo.stop(gs.height == 0, mo.md("_no groups_"))
    top = gs.sort("n_cells", descending=True).head(40).to_pandas()
    sizes = alt.Chart(top).mark_bar().encode(
        x=alt.X("group_id:N", sort="-y", title=None,
                axis=alt.Axis(labels=top.shape[0] <= 25)),
        y=alt.Y("n_cells:Q", title="cells"),
        tooltip=["group_id", "n_cells", "n_callable_sites", "mean_cov"],
    ).properties(title="Group sizes (top 40 by cells)", height=220)
    callable_ch = alt.Chart(top).mark_bar(color="#5a8").encode(
        x=alt.X("group_id:N", sort="-y", title=None,
                axis=alt.Axis(labels=top.shape[0] <= 25)),
        y=alt.Y("n_callable_sites:Q", title="callable sites"),
        tooltip=["group_id", "n_callable_sites", "mean_cov"],
    ).properties(title="Callable sites per group", height=220)
    mo.hstack([sizes, callable_ch], widths=[1, 1])
    return


@app.cell
def _(alt, mo, mpc, pl):
    mo.stop(mpc.height == 0, mo.md("_no coverage_"))
    agg = (mpc.group_by("n_molecules").agg(pl.col("freq").sum())
           .sort("n_molecules").to_pandas())
    hist = alt.Chart(agg).mark_bar().encode(
        x=alt.X("n_molecules:O", title="molecules per cell at a callable site"),
        y=alt.Y("freq:Q", title="(cell, site) pairs"),
        tooltip=["n_molecules", "freq"],
    ).properties(
        title="Gate-1: molecules per cell per callable site", height=260)
    mo.md("## Gate-1 view")
    return (hist,)


@app.cell
def _(hist):
    hist
    return


@app.cell
def _(cand, mo):
    mo.md(f"## Candidate universe ({cand.height} sites)")
    min_vaf = mo.ui.slider(0.0, 1.0, value=0.0, step=0.01, label="min VAF")
    min_alt = mo.ui.slider(1, 20, value=1, step=1, label="min alt molecules")
    mo.hstack([min_vaf, min_alt], justify="start", gap=2)
    return min_alt, min_vaf


@app.cell
def _(cand, min_alt, min_vaf, mo):
    filt = cand.filter(
        (cand["vaf"] >= min_vaf.value) & (cand["n_alt_molecules"] >= min_alt.value)
    ).sort("n_alt_molecules", descending=True)
    mo.ui.table(filt, selection=None, page_size=15)
    return (filt,)


@app.cell
def _(alt, filt, mo):
    mo.stop(filt.height == 0, mo.md("_no candidates pass the filter_"))
    vaf = alt.Chart(filt.to_pandas()).mark_bar().encode(
        x=alt.X("vaf:Q", bin=alt.Bin(maxbins=25), title="VAF (alt molecules / depth)"),
        y=alt.Y("count():Q", title="candidates"),
    ).properties(title="Candidate VAF distribution", height=200)
    vaf
    return


if __name__ == "__main__":
    app.run()
