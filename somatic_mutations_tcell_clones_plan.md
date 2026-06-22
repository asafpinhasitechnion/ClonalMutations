# Project Plan — Detecting Somatic Mutations in T-cell Clones from Paired scRNA + TCR Sequencing

> **One-line summary.** Detect somatic mutations in (primarily) T-cell clones from paired single-cell RNA and TCR sequencing, using the TCR as a *clonal barcode* to aggregate sparse, error-prone reads across clonally related cells. Compare clone-level pooling ("pseudobulk") against a cell-resolved hierarchical model. Build on **SComatic** for the clone-agnostic candidate-generation and filtering layer; innovate on clone-resolved calling and cross-clone classification.

This document is self-contained: a reader who has not seen the originating discussion should be able to act on it.

---

## SUPERVISION & WAYS OF WORKING — READ FIRST

> This is a **supervised** research project. The project owner stays on top of all consequential decisions. Any agent or session working from this plan must operate within the protocol below. **Default posture: do routine, reversible work autonomously, but stop and check in at every decision point and before any consequential or hard-to-reverse choice.** When in doubt, ask rather than proceed.

1. **Stop at phase/stage gates.** Do not advance from one Stage to the next without (a) presenting the current stage's outputs in a reviewable form and (b) getting explicit sign-off. **Phase 0 → direction is the critical gate:** do not commit to de novo vs targeted/pooled, or to a leading setting (healthy vs disease), until the owner has seen the triage numbers and approved.

2. **Surface assumptions before acting on them.** Whenever you set a parameter, threshold, dataset inclusion/exclusion, reference/resource, or modeling assumption, state it explicitly and label it as a *documented default* vs. a *new choice*. Proceed on defaults for routine steps; **pause for sign-off** on anything that affects scientific direction or results.

3. **Never silently resolve an open question.** The Part III open questions (coverage/power, the Stage E cross-clone model, expected yield, T-vs-B scope, FDR at scale, NMD/editing handling) are brought to the owner with options and a recommendation — never decided unilaterally.

4. **Proceed-vs-check, explicitly.** *Proceed autonomously* on routine reversible engineering (running the pipeline, intermediate files, standard QC, refactors, reproducible re-runs). *Stop and check first* on: committing to a methodological choice that affects results; deviating from this plan; deleting/overwriting prior results; large or costly compute runs; selecting or excluding datasets; anything touching an open question.

5. **Calibrate before trusting.** Do not report or build on somatic calls until the germline-recovery calibration (Stage 0 / Stage F) has been run and shown to the owner. Calibration is a gate, not an afterthought.

6. **Report legibly and honestly.** Present intermediate outputs (triage numbers, candidate counts, calibration metrics, the Arm 1 vs Arm 2 gap) in an auditable form alongside the assumptions that produced them. **Report negative results plainly** — e.g., if the hierarchical model does not beat pseudobulk, that is a finding to surface, not a shortfall to hide. No tuning of thresholds to manufacture a desired outcome; no over-claiming; state limitations.

7. **Flag deviations immediately.** If reality forces a departure from this plan (a resource is missing, an assumption fails, results contradict expectations), stop, surface the situation and the proposed adaptation, and wait for direction rather than quietly adapting.

8. **Keep the decision log current.** Record each consequential decision — what, why, and that the owner approved it — in Part III's "Settled decisions," so context accumulates across sessions. Treat this document as the living source of truth and propose updates to it as the project evolves.

**Default check-in cadence:** at every stage gate; before any consequential or irreversible choice; whenever results contradict expectations; and whenever an open question is reached.

---

## PART I — CONCEPTUAL PLAN

### 1. Problem and gap
Detecting somatic mutations at single-cell resolution is limited. Single-cell DNA and targeted assays can do it but are neither standardized nor abundant. scRNA-seq, by contrast, is ubiquitous and exists at enormous scale — a largely untapped substrate for mutation detection. No scalable approach currently measures somatic mutation *and* gene expression in the same cells, which is what is needed to link genotype to transcriptional state in real human data.

### 2. Why immune cells
Somatic mutation in immune cells has been understudied relative to cancer. T (and B) lymphocytes undergo clonal expansion — a receptor-defined ancestor gives rise to large clones of receptor-identical descendants — and mutations can arise and be inherited by the whole clone. This is an opportunity hiding in existing data.

### 3. Core idea
The TCR is a **somatically generated clonal barcode**: clone members descend from one ancestor and therefore share its *truncal* mutations. Pooling reads across a clone is thus justified **by descent**, not merely by statistical convenience. This raises effective depth at expressed positions and overcomes scRNA-seq sparsity. (Contrast: SComatic pools by *cell type*, which groups genetically unrelated cells; the clone is a genetically principled group.)

### 4. The two modeling approaches — and why both
- **Pseudobulk (Arm 1).** Collapse a clone's reads at a site; test alt-allele fraction against a background error model (beta-binomial). Simple, fast, battle-tested (this is essentially SComatic's core).
- **Hierarchical / cell-resolved (Arm 2).** A generative model over `read → molecule → cell → clone`. Evidence accumulates across *molecules within a cell* (rejecting molecule-level artifacts) and *cells within a clone* (raising power and confirming clonality).
- **Key relationship.** Pseudobulk is the **degenerate limit** of the hierarchical model when every cell contributes ≤1 molecule per site. As coverage rises they diverge. **The sensitivity/specificity gap between them, as a function of coverage and clone size, is the central methodological result** ("how much does resolving cells buy you, and when").

### 5. Strategy
- **Reuse SComatic** for the *clone-agnostic* layer: read/base QC, RNA-editing and panel-of-normals (PoN) filtering, clustered-variant flagging, candidate generation.
- **Innovate** on three things SComatic does not handle for this setting: (a) **clone** (not cell-type) grouping; (b) **cell-resolved hierarchical calling**; (c) **cross-clone classification** (germline vs recurrent vs clone-specific vs pre-thymic).

### 6. Validation philosophy
**No dependence on matched DNA.** Calibrate sensitivity with germline-heterozygous-site recovery; check the null with germline-homozygous-reference sites; corroborate orthogonally with **mtDNA lineage variants** (do TCR clones share them; are somatic calls consistent); use known **driver recurrence** in disease data as positive controls; use matched WES/WGS as validation-only where it happens to exist.

### 7. Phases (conceptual)
- **Phase 0 — Feasibility triage** (do first; gates everything).
- **Phase 1 — Shared pipeline** (align, assign clones, generate + filter candidates).
- **Phase 2 — Pseudobulk arm.**
- **Phase 3 — Hierarchical arm.**
- **Phase 4 — Cross-clone classification layer.**
- **Phase 5 — Benchmarking + validation.**
- **Phase 6 — Cohort application + biology.**

### 8. Scope and standing caveats
- **Coverage decides the headline.** If most cells contribute ≤1 molecule/site, the hierarchy collapses to pseudobulk and adds nothing; Phase 0 measures this directly.
- **T cells are the clean target** (no somatic hypermutation; stable barcode). **B cells are a hard extension**: AID-driven SHM injects mutations to be separated from the signal, and the BCR barcode itself mutates (clone definition needs lineage clustering, not exact match). Decide scope consciously.
- **Detectable space is limited** to expressed, exonic, and (with the TCR-compatible 5′ chemistry) 5′-biased positions.
- **Low somatic burden in healthy T cells** (most mutations intronic/intergenic, invisible to RNA) means de novo discovery may require a disease / clonal-expansion context; a **targeted** formulation (known driver loci) is the fallback and is much easier.

---

## PART II — DETAILED PLAN (implementation base)

## FIRST DELIVERABLE — EXTRACTION + EXPLORER (BUILD THIS FIRST)

**Build this before any modeling.** It is the shared substrate for both arms, it answers the feasibility and trust questions, and the explorer's findings decide *which* modeling arm to build. It consolidates Stage 0, the extraction half of Stage B, and the visualization data layer into one build target; the per-stage detail below is the component reference. Two decision gates are built in.

### Build sequence
1. **Front-end pipeline.** Align (STAR) → assign clones → read/base QC → **UMI-collapse to consensus molecules** (see Stage A) → bake-in blacklist (TCR/Ig) → permissive candidate detection.
2. **Tier-A extraction** (cheap, genome-wide summary) → Parquet.
3. **Explorer v1** = distribution views over Tier A. This *is* the Phase 0 triage.
   → **GATE 1:** do enough cells carry ≥2 molecules at candidate sites? Yes → plan the hierarchical arm. No → drop it; pseudobulk/targeted only.
4. **Tier-B extraction** (rich, candidate sites only; include germline calibration sites) → Parquet + raw drill-down.
5. **Explorer v2** = per-site evidence view + germline calibration.
   → **GATE 2:** germline recovery vs expectation (trust)? dominant artifact modes (which filters to strengthen)? where signal sits — healthy vs disease, which clones/genes (scope/focus)?
6. **Then build the modeling arm** the explorer points to — pseudobulk regardless, hierarchical only if Gate 1 = yes — reusing the same extracted data.

### Filters: baked vs flagged
- **Baked into the front-end (irreversible — only "certain non-usable"):** multimapping reads (`max_NH 1`); a conservative mapping-quality floor; the TCR/Ig rearranged loci (reference wrong by construction).
- **Computed as flags on candidates, adjustable in the explorer (never used to exclude from extraction):** RNA editing (DB membership + A>G/T>C signature); PoN/germline; clustered-variant flag (adjustable distance); HLA region; base-quality and mismatch as retained tags; the candidate-detection threshold (kept permissive).
- **Mechanic:** candidate generation uses **only** the baked filters + the permissive threshold to decide *where* to extract; the judgment filters are annotations on every extracted candidate. The explorer holds all candidates + all flags and filters interactively, so nothing borderline is lost before it can be inspected.

### Data layer — Parquet schemas (indicative)
*Tier A:*
- `clone_summary.parquet` (per clone): `clone_id, n_cells, n_callable_sites, mean_cov, median_mol_per_cell_per_site, mol_per_cell_hist[]`
- `coverage_by_site.parquet` (per site×clone, optional): `site_id, clone_id, depth, n_cells_covered, n_molecules`

*Tier B (candidate sites only; `class ∈ {somatic_candidate, germline_het, germline_homref}`):*
- `candidate_clone_site.parquet` (per candidate × clone): `site_id, chrom, pos, ref, alt, gene, class, clone_id, depth, n_molecules, n_alt_molecules, n_cells, n_alt_cells, vaf, strand_alt_frac, mean_consensus_q, flag_editing, flag_pon, flag_clustered, flag_hla`
- `candidate_cell_site.parquet` (per candidate × cell): `site_id, cell_id, clone_id, n_molecules, n_alt_molecules, consensus_q, carrier_post` *(carrier_post filled after modeling)*
- raw drill-down (per candidate, on-disk, single-site inspection only): `site_id, cell_id, clone_id, umi, consensus_base, n_reads, agree_frac`

### Explorer app (the "simpler app")
- **Principle:** never touch raw reads at interaction time — interactive views read Parquet via DuckDB; raw molecule pileup is fetched only when drilling into one site.
- **Stack:** DuckDB (+ Polars) over Parquet; **marimo** or **Panel** for the reactive UI; **Datashader** for large heatmaps/embeddings; optional **igv.js** (or link to IGV) for raw single-site pileup. Local-only, no server.
- **v1 views (Tier A) = triage / Gate 1:** molecules-per-cell-per-site histogram (per clone & overall); clone-size distribution; callable-sites-per-clone; coverage by gene.
- **v2 views (Tier B) = Gate 2:** per-site evidence panel (molecule→cell→clone breakdown for a selected candidate + all flags); germline calibration panel (het recovery vs expectation; hom-ref null); signal-by-clone/gene; filterable candidate table with live flag toggles.
- **Later (as modeling produces outputs):** cell×mutation carrier heatmap; cross-clone pattern (Stage E classes); Arm 1 vs Arm 2 diff; mutation↔expression on the GEX embedding.
- **Keep it lightweight** — an exploratory instrument, not a product. Until interaction becomes repetitive, DuckDB queries + notebook plots may suffice.

### Inputs / data model
- **Aligned reads (BAM)** with: cell barcode tag `CB`, UMI tag `UB`, and `nM` (mismatches) + `NH` (hits) tags. 5′ single-cell chemistry (the TCR-V(D)J–compatible kit), stranded.
- **TCR clonotype assignments** per cell (Cell Ranger `vdj`, scRepertoire, or TRUST4).
- **Barcode → clone map** and **clone-size table**.
- **Reference genome** (+ `.fai`), trinucleotide context.
- **Resources:** RNA-editing site list (REDIportal-derived), PoN, germline AF database (gnomAD), blacklist BED (see Stage B).
- Exclude doublets. Prefer clones defined by paired α+β.

### Stage 0 — Feasibility triage (RUN FIRST)
Outputs three numbers, each a decision gate:
1. **Molecules per cell per callable site within expanded clones** → does the hierarchy have purchase over pseudobulk? (If modal cell ≤1 molecule everywhere → reframe to pooled/targeted.)
2. **Callable footprint × plausible per-cell mutation density** → expected detectable clonal mutations per clone; compare healthy vs disease → decides the paper's leading setting.
3. **Germline-het recovery rate** vs model expectation, and **germline-hom-ref false-call rate** → first calibration of the read/molecule pipeline and error model.
> SComatic's "callable sites per cell / per cell type" utilities can be repurposed with clone groupings for (1)–(2).

### Stage A — Preprocessing & clone assignment
- Align with a splice-aware aligner (STAR) if needed; ensure `CB`/`UB`/`nM`/`NH` tags present.
- **UMI collapse (do this early).** Collapse the reads of each UMI into a **consensus molecule**, retaining read support as a per-molecule confidence: store the consensus base + a consensus error `ε_m` (a function of read count and agreement). This removes PCR-duplicate redundancy *without* discarding sequencing-error information — read multiplicity is summarized into `ε_m`. **The molecule, not the read, is the atomic unit everywhere downstream**, including the explorer, which never displays multiple reads per molecule. (Avoid a naive bare consensus that drops read support.)
- **Clone definition:** nucleotide-level CDR3 + V/J genes; paired α+β preferred; handle dual-TCR (allelic inclusion), single-chain dropout, and doublets explicitly.
- Output: `barcode → clone` table; clone-size distribution (drives power expectations).

### Stage B — Candidate generation & filtering (CLONE-AGNOSTIC; SComatic-based)
Generate candidates from **pooled/merged** data, not per-clone (the user's point: candidate detection is separate from clone assessment).
- **Read/base QC (inherit SComatic defaults, tune to aligner):** `--max_NH 1` (unique mappers — *critical in T cells*), `--max_nM ~5`, `--min_MQ` (aligner-dependent; 255 for STAR-unique), `--n_trim` (zero read-end qualities), `--min_bq ~20`.
- **Blacklist (ADDITION — not in SComatic):** exclude **TCR, Ig, and HLA loci** via `--bed_out`. TCR/Ig are somatically rearranged (reference wrong by construction); HLA is hyperpolymorphic with chronic mapping artifacts.
- **Permissive candidate generation:** keep any clean site with ≥ *k* alt **molecules** across the pooled data. **Do NOT use the pseudobulk significance call as the gate** — otherwise Arm 2 can only ever recover a subset of Arm 1 and can never demonstrate added sensitivity. The two callers must adjudicate the *same* candidate universe.
- **External filters:** RNA editing (REDIportal **+ augment with an A>G/T>C strand-aware signature flag**, because activation-induced ADAR1 in expanded T cells creates editing sites absent from static DBs); custom **PoN** (build from the same 5′ chemistry and, ideally, non-T cells of the same donors); **clustered-variant flag** (`--min_distance 5`).
- **UMI awareness (KEY ADAPTATION):** collapse reads → molecules; all support thresholds count **molecules, not reads** (three reads can be one over-amplified erroneous molecule).
- **Output:** (i) clean candidate site list; (ii) **per-clone, per-cell, per-UMI, per-read pileups** at candidate sites. *(SComatic emits collapsed count tables, which are sufficient for Arm 1 but NOT for Arm 2 — Arm 2 needs this richer extraction.)*

### Stage C — Pseudobulk caller (Arm 1) ≈ SComatic core with clone groups
- For each `clone × candidate site`: pool molecules; run the **read-count beta-binomial** (alt vs background error) and the **cell-count beta-binomial** (number of alt-bearing cells vs background).
- **Re-parameterize** α/β for this chemistry (SComatic ships estimation scripts). Feed **UMI-collapsed** counts.
- **FDR** across the `clone × site` space — at hundreds of clones the test count is 10–100× SComatic's, so use explicit FDR / empirical-Bayes calibration (see Stage D's outer loop); SComatic defaults will not be calibrated for this scale.
- Reuse SComatic code where possible; the main changes are grouping label, UMI-collapse, re-parameterization, and FDR.

### Stage D — Hierarchical caller (Arm 2)
Generative model over `read → molecule → cell → clone`. Run on the **same candidate set + UMI pileups** from Stage B.

**Latent variables** (per candidate site *j*): `Z_cj ∈ {0,1}` clone carries mutation; `π_cj ∈ [0,1]` within-clone carrier cell-fraction (truncal→1, subclonal<1); `S_ij ∈ {0,1}` cell carrier (`S_ij | Z_cj=1 ~ Bern(π_cj)`); `φ_ij ∈ [0,1]` fraction of a carrier cell's molecules from the mutant allele (ASE + transcriptional bursting); `G_j ∈ {0,1}` germline competitor.

**Background parameters:** `η_j` site-specific rate a ref-source molecule is *encoded* alt (RNA editing + RT/PCR), context/strand-specific; `δ` mutant→ref molecule flip; `ε` per-read sequencing error (from base quality); `ρ_j` small somatic prior (rarity).

**Likelihood (bottom-up):**
- Per-read: `alt-reads | content=alt ~ Bin(R, 1-ε)`, `| content=ref ~ Bin(R, ε)`.
- Per-molecule content probability: `q = η_j` (non-carrier) or `q(φ) = φ(1-δ) + (1-φ)η_j` (carrier / germline-het).
- Per-molecule marginal: `L_m(q) = q·Bin(a;R,1-ε) + (1-q)·Bin(a;R,ε)`.
- Per cell: carrier → `∫ p(φ) Π_m L_m(q(φ)) dφ`; non-carrier → `Π_m L_m(η_j)`.
- Per clone: `Π_i [ (1-π_cj)·P(cell|S=0) + π_cj·P(cell|S=1) ]` (Z=1); `Π_i P(cell|S=0)` (Z=0).
- Germline `G_j=1` overrides every cell of every clone to the het emission regardless of Z; the site model weighs **germline (signal across all clones) vs somatic (clone-restricted) vs nothing**.

> **With UMI consensus (Stage A)** the per-read binomial collapses to one consensus base per molecule with error `ε_m` (`R→1`, `a∈{0,1}`), so `L_m(q)` becomes a simple 2-term mixture over the consensus base. The molecule/cell/clone structure above is unchanged; read support lives inside `ε_m`.

**Priors:** `Z ~ Bern(ρ_j)` (rarity = the multiple-testing control); `φ ~ Beta`/3-component mixture with heavy mass near 0/1 (bursting); `π` skewed toward 1 (truncal) with a subclonal tail; `η_j` estimated from background + editing DBs + context; `ε` from (recalibrated) qualities; `G_j` from population AF.

**Inference (no MCMC required for a first version):**
- Sum out molecule content and `S_ij` analytically (finite mixtures).
- `φ` integral is 1-D → quadrature, or closed-form **Beta-Binomial** under a Beta prior (folding the small `η` cross-term into background). `π` likewise 1-D.
- Enumerate `Z_cj`, `G_j`.
- **Empirical-Bayes outer loop (EM):** E-step computes posteriors over `Z, S, G` (with the φ, π integrals); M-step learns the **shared hyperparameters** `η` (context model), `ε` recalibration, `ρ`, and the φ/π prior shapes across all sites. *This is the "null is calibrated by the genome's emptiness" mechanism.* Swap EM→MCMC only if full posterior uncertainty is needed.

**Outputs:** `P(Z_cj=1 | data)` (the call); posterior on `π_cj` (clonal prevalence / timing); per-cell `P(S_ij=1 | data)` (the **cell × mutation carrier matrix**, input to within-clone phylogeny later).

**Decision rule:** rank by `P(Z=1|·)`; admit until running mean of `1 - P(Z=1|·)` (local FDR) hits target. No separate frequentist correction.

### Stage E — Cross-clone classification (REDESIGN; no off-the-shelf equivalent)
SComatic's `--max_cell_types 1` rule **breaks** with hundreds of clones (germline now appears in *hundreds*, and the rule would discard all recurrent drivers). Replace with a **graded model over the fraction of covered clones carrying a candidate**, resolving four classes:
- **Germline** — present in a large fraction of covered clones.
- **Recurrent somatic** — present in a few clones (candidate driver hotspots).
- **Clone-specific somatic** — present in one clone.
- **Pre-thymic / clonal hematopoiesis** — shared across several clones with *different* TCRs (cannot share by TCR descent ⇒ a mutation predating thymic rearrangement, in a common progenitor). *This is a novel readout the clone resolution uniquely enables.*
> This module needs design work; it is the main statistical novelty beyond SComatic. Parameterizing "fraction of covered clones" into these classes (accounting for coverage per clone and recurrence vs shared-descent) is an open task.

### Stage F — Benchmarking & validation
- **Headline:** Arm 1 vs Arm 2 sensitivity/specificity **as a function of coverage and clone size**.
- **Calibration:** germline-het recovery vs expectation; germline-hom-ref false-call rate (null).
- **Orthogonal:** mtDNA lineage variants (scMitoMut-style beta-binomial) — concordance of TCR clones with shared mtDNA variants; consistency with somatic calls.
- **Positive controls:** known driver recurrence in disease datasets.
- **Validation-only:** matched WES/WGS where available (not a dependency).

### Stage G — Application & biology
- Assemble the cohort (target *X* datasets / conditions).
- Link mutation ↔ transcriptional state; reconstruct within-clone structure from the carrier matrix; mutational signatures (trinucleotide context); characterize pre-thymic/CH events.

### SComatic filter audit — fit in the clone setting (summary table)
| SComatic component | Default | Fit for clones | Adaptation |
|---|---|---|---|
| `--max_NH`, `--max_nM`, `--min_MQ`, `--n_trim`, `--min_bq` (read/base QC) | varies | **Transfers wholesale** (grouping-agnostic) | Tune `min_MQ` to aligner; `max_NH 1` essential |
| TCR/Ig/HLA blacklist | — (absent) | **Missing — must add** | Add via `--bed_out` |
| RNA-editing filter | REDIportal | **Transfers; more important here** | Augment with A>G signature (activation-induced ADAR1) |
| PoN | provided | **Transfers** | Build custom from 5′ chemistry + non-T normals |
| Clustered-variant flag `--min_distance` | 5 | **Transfers untouched** | Keep |
| `min_dp`, `min_cells` (coverage floors) | 5, 5 | Transfers but **restricts to expanded clones** | Embrace as scope statement |
| `min_ac_cells` (≥ cells supporting alt) | 2 | **Good fit** (clonal consistency) | Possibly raise for stringency |
| `min_ac_reads` (≥ reads supporting alt) | 3 | **Trap with UMIs** | Recast as ≥ *k* **molecules** |
| Read-count beta-binomial (α1/β1) | fitted | Transfers as Arm 1 core | **Re-parameterize**; feed UMI-collapsed counts |
| Cell-count beta-binomial (α2/β2) | fitted | Coarse precursor to Arm 2 | Re-estimate; expect weak for small clones |
| Strand bias `--fisher_cutoff` | off | Low priority (RNA strand semantics) | Leave off; molecule/cell consistency is better |
| `--max_cell_types 1` (cross-group germline) | 1 | **BREAKS** | Replace with graded cross-clone model (Stage E) |
| `--min_cell_types 2` (callability) | 2 | Reinterpret | "min clones covered to assess cross-clone pattern" |

### Dependencies (indicative)
STAR, samtools, bedtools, pysam; SComatic (filtering/PoN/beta-binomial scaffold); Cell Ranger `vdj` / scRepertoire / TRUST4 (clonotypes); a beta-binomial/optimization stack (e.g. `scipy`, or `VGAM` in R as SComatic uses); REDIportal, gnomAD; mtDNA tooling (e.g. scMitoMut) for validation.

---

## PART III — DECISIONS, RATIONALE, OPEN QUESTIONS

### Settled decisions (with rationale)
- **Group by TCR clone, not cell type.** Clone members share ancestry → pooling justified by descent; finer and more principled than SComatic's cell-type grouping.
- **Build both arms; the gap is the result.** Pseudobulk is the degenerate limit of the hierarchical model; comparing them measures the value of cell resolution and de-risks the project (clone-pseudobulk alone is already novel and publishable).
- **Build on SComatic for the clone-agnostic layer.** Its filters/resources encode years of artifact handling and dominate the genome-scale false-positive rate; reinventing them would be worse and slower.
- **Separate candidate generation (pooled, clone-agnostic) from clone assessment.** Avoids redundant per-clone scans and keeps the Arm 1 vs Arm 2 comparison fair.
- **Count molecules, not reads.** Molecule is the independent unit; read multiplicity within a molecule only rejects sequencing error.
- **No matched-DNA dependency; calibrate via germline recovery + mtDNA orthogonal check.**
- **Primary scope = T cells.** B cells deferred due to SHM + mutating barcode.
- **Extract + explore before modeling.** The shared extraction + a thin explorer is the first deliverable; the explorer's findings (coverage, calibration, artifacts) decide which modeling arm to build, so modeling sits behind two decision gates.
- **UMI consensus collapse, with retained read support.** Collapse early to one consensus molecule per UMI, keeping read count/agreement as the molecule's confidence (`ε_m`); the molecule is the atomic unit. Naive bare-consensus (discarding read support) is explicitly avoided.
- **Filters: baked vs flagged.** Bake in only "certain non-usable" filters (multimapping, MQ floor, TCR/Ig loci); carry editing, PoN, clustered, HLA, and quality/mismatch cutoffs as adjustable flags in the explorer (extract permissively, filter interactively).

### Open questions / risks
1. **Coverage (Phase 0).** Does the hierarchy beat pseudobulk on real data? Determines whether the de novo ambition is powered or we lead with targeted/pooled.
2. **Cross-clone model (Stage E).** How exactly to parameterize "fraction of covered clones" into germline / recurrent / clone-specific / pre-thymic, accounting for per-clone coverage and distinguishing recurrence from shared descent. **No off-the-shelf solution.**
3. **Yield in healthy tissue.** May be near zero in expressed space → may require disease/clonal-expansion datasets as the leading context.
4. **NMD bias.** Loss-of-function (PTC) mutations under-detected (mutant transcript degraded) → state as a known sensitivity bias.
5. **Editing beyond databases.** Activation-induced editing sites; how aggressive the A>G-signature augmentation should be without discarding real C>T/G>A somatic signal.
6. **FDR at hundreds-of-clones scale + group-size heterogeneity.** Confidence should scale with clone size, not pass a uniform threshold.
7. **B-cell extension** (if pursued): lineage-aware clone definition + SHM separation.

### Suggested build order for a coding agent
**Build the extraction + explorer first (see "FIRST DELIVERABLE" in Part II); defer modeling until the explorer answers which arm to build.**
1. Front-end: align → clone-assign → read/base QC → **UMI-collapse to consensus molecules** → bake-in blacklist (TCR/Ig) → permissive candidate detection.
2. **Tier-A extraction** (cheap, genome-wide summary) → Parquet.
3. **Explorer v1** (distribution views = Phase 0 triage). **Gate 1:** enough cells with ≥2 molecules at candidates? → decides whether to build the hierarchical arm.
4. **Tier-B extraction** (rich, candidate sites only; include germline calibration sites) → Parquet + raw drill-down.
5. **Explorer v2** (per-site evidence + germline calibration). **Gate 2:** trust, dominant artifacts, scope/focus.
6. **Modeling arm** the explorer points to — pseudobulk (Stage C) always; hierarchical (Stage D) only if Gate 1 = yes — reusing the same extracted data.
7. Stage E (cross-clone classification), then Stages F–G (benchmark, validate, apply).
