# Runbook: UCSF Perturb-seq Analysis (2026-04-17)

## Goal

Define a practical analysis sequence for the UCSF perturb-seq production data so
we do not jump from release packaging straight into ad hoc plots.

This runbook is organized into three tiers:

1. **Tier 1: QC / release gate**
2. **Tier 2: exploratory analysis**
3. **Tier 3: biological validation**

The intent is:

- Tier 1 answers: "is this a usable perturb-seq dataset?"
- Tier 2 answers: "what biology changes under perturbation?"
- Tier 3 answers: "do those changes make biological sense?"

## Scope / dataset surface

Authoritative UCSF analysis surface:

- guarded rerun release root:
  `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009`
- release/readme summary:
  `scripts/ucsf_release_guarded_README.md`

Default analysis inputs:

- per-sample `filtered_counts.h5ad`
- optionally per-sample `default_singlet_filtered_counts.h5ad` if we want the
  stricter singlet-only surface

Do **not** use `final_counts.h5ad` as the default exploratory input. It is the
full barcode surface and is much larger than the filtered analysis surface.

Current guarded rerun scale:

- `16` samples total
- filtered cells across all samples: about `156,648`
- default singlet-filtered cells across all samples: about `144,761`
- per-sample filtered cell range: about `6.8k` to `11.9k`

This is comfortably small enough for Scanpy/Seurat exploratory analysis.

## Principle

Treat perturb-seq analysis as two separate questions:

1. **Did the perturbation work?**
2. **What cellular programs changed after it worked?**

Those are related but not interchangeable.

Leiden/UMAP is useful for question 2. It is **not** the first-line answer for
question 1.

## Tier 1: QC / release gate

Tier 1 is split into three parts, run in order:

1. **technical QC** (per-sample) — Tier 1A
2. **dataset-level tenor** (pooled, perturbation-agnostic) — Tier 1B
3. **perturbation-efficacy QC** (per-target) — Tier 1C

1A confirms each sample is individually usable. 1B then confirms the collection
of 16 samples hangs together as a single coherent dataset before we commit
compute to 1C's per-target efficacy pass. If any of the three fails, the
dataset is not analysis-ready in the perturb-seq sense — and the specific
failure mode tells us whether to re-run a sample, re-run integration, or
investigate the library.

### Tier 1A: technical QC

Minimum checks:

- cell counts per sample
- genes / UMIs / mitochondrial fraction distributions
- doublet rate per sample
- guide-assignment rate per sample
- fraction of `single` / `multi` / `ambiguous` / `none` guide calls
- top-guide abundance and long-tail guide coverage

Minimum plots:

- histogram / violin of `n_genes`
- histogram / violin of `total_counts`
- histogram / violin of `mt_pct`
- barplot of guide call categories
- barplot of top guides / targets by cell count
- cell-count summary table across all 16 samples

Existing surfaces already available:

- `gene_quantile_histogram.html`
- `adaptive_qc_threshold.json`
- `summary.txt`
- CRISPR call summaries under `run/outs/crispr_analysis/`

Tier 1A exit criteria:

- no sample has obviously broken cell-level QC
- no sample is dominated by `none` or pathological `ambiguous` guide calls
- guide coverage is broad enough that downstream per-target summaries are worth
  doing

### Tier 1B: dataset-level tenor (perturbation-agnostic)

Tier 1A tells us that each sample is individually OK. Tier 1B asks a different
and equally important question before we commit compute to per-target work:
**does the 16-sample collection hang together as a single coherent perturb-seq
dataset?** This is the Replogle-style Figure 1 orientation pass — it is
perturbation-agnostic on purpose, because we want to catch gross batch failure,
integration collapse, or sample-dominated clusters before 48 target-level QCs
are running against them.

Tier 1B runs a light first integration pass and produces the dataset-orientation
card *without* perturbation overlays. The same `sample_overview_umap` card is
re-emitted in Tier 2A with Mixscape and feature-call overlays populated once
Tier 1C has run — this is an intentional two-pass pattern, not duplicated work.

Required compute:

- run a Harmony-integrated cross-lineage embedding of all 16 samples on the
  filtered counts (PCA → Harmony → neighbors → UMAP → Leiden)
- this is explicitly a first-pass integration, not the analytical product;
  the lineage-specific integrations in Tier 2A are the product

Minimum metrics (perturbation-agnostic):

- total cell counts — overall, per lineage, per sample
- per-sample cell-count range and coefficient of variation
- **integration mixing** — sample-entropy per cluster or kBET / LISI per
  cluster, to detect sample-dominated clusters
- **lineage separation** — silhouette / ASW of EBs vs iPSC on the integrated
  embedding
- **cluster × sample composition** — matrix, with entropy or Shannon
  diversity per cluster
- **NTC pool adequacy** — number of NTC cells per lineage, to confirm that
  the per-target NTC comparisons in Tier 1C are well-powered
- overall fraction of cells with a `single` guide call vs `multi` /
  `ambiguous` / `none`, aggregated
- (optional) doublet-rate summary across samples, if doublet calls are
  available

Required plot:

- `sample_overview_umap` in light mode (no Mixscape panel), as the Tier 1B
  summary slide:
  - main integrated UMAP colored by Leiden with cluster labels
  - 4 × 4 per-sample small-multiples band on shared coordinates
  - cluster × sample composition strip underneath
  - per-sample QC table (n_cells, median n_genes, median mt_pct, n_NTC_cells)

Biological questions this figure must answer at a glance:

- do EBs and iPSC separate cleanly on the shared manifold?
- are the samples within a lineage reproducible, or does one sample dominate?
- is any Leiden cluster an obvious batch artifact (present in one sample only
  or dominated by a single sample)?
- is the NTC pool adequate for per-target comparison in each lineage?

Tier 1B exit criteria (quantitative gate):

- ≥`6000` cells per sample after per-sample QC (guarded rerun range is
  ~`6.8k`–`11.9k`, so this is a sanity floor)
- ≥`1000` NTC cells per lineage
- no Leiden cluster is >`70%` dominated by a single sample (sample-entropy
  threshold), **or** such clusters are identified and annotated as
  batch-suspect before Tier 1C runs against them
- EBs and iPSC are clearly separable on the integrated UMAP (silhouette /
  ASW above a reviewed threshold; exact threshold is an open decision — see
  Open questions)
- guide-call categorical distribution is roughly consistent across samples
  (no single sample catastrophically low on `single` calls)

### Tier 1C: perturbation-efficacy QC

This is the missing piece if we stop at generic scRNA QC.

A simple "target-guide mean vs NTC mean" comparison is necessary but not
sufficient. It cannot distinguish "50% of cells strongly knocked down" from
"100% of cells weakly knocked down," and it misses perturbations that act
strongly in representation space without self-regulating the target mRNA. Tier
1B therefore runs three complementary checks, not one.

Primary question:

- for each targeted gene, do cells carrying guides for that target actually show
  lower expression of that target than NTC / non-targeting controls?

#### Tier 1C.1: target-expression knockdown (mean-level)

- compare target-gene expression in target-guide cells vs NTC cells
- summarize per target, per sample, **and** per lineage
- track:
  - number of perturbed cells
  - median / mean target expression
  - log2 fold change vs NTC
  - p-value / FDR
  - whether the target gene is expressed at a usable baseline in that lineage
- report both **guide-level** and **target-level** summaries:
  - target-level is the primary efficacy call
  - guide-level is the consistency overlay (guide-to-guide concordance per
    target, e.g., Spearman of guide-level DE signatures)
  - this closes the "guide vs target" open question rather than leaving it
    implicit

Recommended plots:

- target-knockdown dotplot across samples
  - dot size = number of perturbed cells
  - color = target-gene log2FC vs NTC
  - significance ring or label
- guide-level knockdown strip plot per target (all guides side by side) as the
  guide-concordance QC visual
- power-ranked target summary table
- plot restricted to genes with adequate baseline expression

Relevant existing script:

- `scripts/perturbseq_qc_dotplot.py`

That script already encodes the right idea for Tier 1C.1:

- per-sample target vs NTC comparisons
- cell-count-aware summaries
- expression-aware filtering

#### Tier 1C.2: per-cell perturbation labeling (Mixscape)

Mean knockdown alone cannot see escapers / incompletely perturbed cells. Run a
per-cell labeling pass using the Mixscape approach (Papalexi et al. 2021):

1. Compute a perturbation signature per target by removing cell-cycle and other
   confounders from the expression of target-guide cells relative to NTC.
2. Fit a two-component Gaussian mixture per target and assign each cell a
   `mixscape_class` of `perturbed` (`KO`) or `non_perturbed` (`NP`).
3. Report per target:
   - fraction of guide-carrying cells labeled `perturbed`
   - per-cell `mixscape_score`
   - number of `NP` escapers to exclude from downstream effect estimation

Downstream Tier 2 DE should default to the Mixscape-`perturbed` subset when
Mixscape runs successfully; Tier 2 contrasts should fall back to the raw
target-guide set only when Mixscape cannot fit (too few cells, degenerate
signature).

#### Tier 1C.3: perturbation magnitude in representation space (E-distance)

Target mRNA readout misses regulators that do not self-regulate and can
underweight regulators with strong downstream effect. Complement Tier 1C.1 with
a representation-space summary:

- compute an **E-distance** (scPerturb, Peidli et al. 2023) between target-guide
  cells and NTC cells in PCA space (top ~50 PCs), per target, per lineage
- permute NTC vs NTC to establish a null distribution
- report per target:
  - E-distance
  - E-distance permutation p-value
  - rank among all targets

A target with high E-distance and weak target-mRNA knockdown is a flag, not a
failure — often a true regulator whose transcript is not itself the primary
readout.

#### Tier 1C.4: CRISPRoff temporal-control / selection-bias QC

Because this experiment is **CRISPRoff**, one additional QC question is
mandatory:

- was perturbation effectively delayed until near assay time, or do we already
  see CRISPRi-like proliferative selection during expansion?

Important distinction:

- ideal CRISPRoff does **not** imply uniform target-mRNA knockdown magnitude
  across genes; promoter context, baseline expression, and target biology will
  make log2FC heterogeneous
- ideal CRISPRoff **does** imply approximately uniform **guide / target
  representation at readout**, after conditioning on library composition and
  sampling depth, because fitness-sensitive knockdowns should not have had much
  time to expand or deplete before profiling

Strong enrichment of tumor suppressor knockdowns together with depletion of
essential / pro-growth knockdowns is therefore a red flag for:

- constitutive CRISPRi rather than true CRISPRoff
- leaky induction
- or CRISPRoff induction occurring too early relative to the expansion window

Required summaries:

- perturbed-cell count distribution per **guide** and per **target**
  - overall
  - per lineage
  - per sample where useful
- rank-ordered target-representation plot colored by expected fitness class
  - tumor suppressor
  - neutral / unknown
  - essential / pro-growth
- class-specific median perturbed-cell counts
- **tumor-suppressor / essential ratio** as the headline temporal-control QC
  metric
- guide / target abundance **Gini** or coefficient of variation across the
  library
- dropout rate by fitness class
  - fraction of guides / targets with `<5`, `<10`, and `<30` cells
- correlation between guide / target abundance and external fitness annotation
  when available
  - e.g. essentiality class, DepMap-style dependency score, literature
    annotation
- lineage-specific representation ratios, because some targets may only be
  strongly fitness-sensitive in iPSC or in EBs

Recommended plots:

- ranked barplot of target representation across the whole library, colored by
  fitness class
- class-median comparison plot
  - tumor suppressor vs neutral vs essential / pro-growth
- per-lineage representation violin / box plot by fitness class
- cumulative dropout curve by fitness class

Interpretation rule:

- **uniformity should be judged on representation, not on knockdown depth**
- if representation is strongly fitness-skewed, the dataset may still be usable
  for neutral-fitness targets, but it is not behaving like ideal CRISPRoff
  temporal control

Relevant existing script / precedent:

- `scripts/generate_crispoff_report.py`

That script already encodes the key logic:

- uniform representation is the CRISPRoff expectation
- tumor-suppressor enrichment plus essential-gene depletion is evidence against
  effective temporal control

#### Tier 1C exit criteria (quantitative gate)

"A substantial fraction" is not a gate. The prespecified pass criterion for the
guarded rerun is:

- **≥50%** of evaluable, baseline-expressed targets show target-level log2FC
  ≤ −0.5 with FDR < 0.05 vs NTC in their primary lineage **and**
- **≥50%** of evaluable targets show E-distance permutation p < 0.05 vs NTC
  **and**
- **≥60%** median guide-to-guide Spearman concordance within target (for
  targets with ≥2 guides and ≥30 cells/guide)

CRISPRoff-specific temporal-control readout:

- report the **tumor-suppressor / essential median representation ratio** and
  class-specific medians as a required QC output
- values near `1x` are the CRISPRoff ideal
- large deviations above `1x` are a red flag for proliferative selection during
  expansion
- exact pass / fail threshold is still an open decision and should be reviewed
  against the full library annotation rather than guessed in code

Evaluability rules (record these explicitly rather than silently dropping):

- `not_expressed_in_lineage` — baseline expression below threshold in the
  relevant lineage; lack of knockdown signal is not evidence of guide failure
- `low_power` — too few perturbed cells (default threshold `n_target_cells < 30`)
- `no_mixscape_fit` — Mixscape could not separate KO/NP
- `suspicious` — adequate power and baseline but directionally wrong or null

Tier 1C should be condition-aware: EBs and iPSC should be summarized separately
because some targets will only be evaluable in one lineage.

### Tier 1 summary

Tier 1 deliverables should be enough for a short "analysis readiness" report.

Suggested outputs:

- one per-sample QC table
- one cross-sample QC summary table
- one target knockdown summary table
- one CRISPRoff selection-bias summary figure / table
- one or two compact perturbation-efficacy figures

## Tier 2: exploratory analysis

Tier 2 starts only after Tier 1 shows the dataset is working as perturb-seq.

This is where Leiden clustering belongs.

### Tier 2A: embeddings and clustering

Framing: **per-sample embeddings are a QC diagnostic; the analytical product is
a lineage-integrated embedding.** 16 samples at ~10k cells each is well inside
the regime where integration is straightforward and necessary — the biology we
care about (perturbation × lineage × state) requires a shared coordinate
system across samples.

Recommended default sequence:

1. **Per-sample** `PCA -> neighbors -> UMAP -> Leiden` on `filtered_counts.h5ad`
   - diagnostic only: confirm that major states are real within a sample before
     trusting integrated views
   - used to catch sample-specific pathology (failed library, dominant debris
     cluster, etc.) before it pollutes a pooled run
2. **Lineage-pooled integrated** analysis as the analytical product:
   - pooled EBs, pooled iPSC
   - integrate by sample with **Harmony** (fast, good default) or **scVI**
     (model-based, preferred when we want latent-space DE later)
   - rerun neighbors / UMAP / Leiden on the integrated latent space
3. **Cross-lineage integrated** (all 16 samples) only after (2), and only to
   answer cross-lineage questions; the primary Tier 2 views remain
   lineage-specific

Do **not** start with a single 16-sample pooled Leiden without explicit
integration. EBs vs iPSC and sample-specific effects will dominate the
embedding otherwise.

Minimum Tier 2A plots:

- UMAP colored by sample (integration diagnostic)
- UMAP colored by Leiden cluster
- UMAP colored by `n_genes`, `total_counts`, `mt_pct`
- UMAP colored by `feature_call_category`
- UMAP colored by Mixscape `mixscape_class` (`KO` vs `NP` vs `NTC`)
- UMAP colored by selected marker genes
- UMAP colored by selected perturbation targets

#### Dataset-orientation summary figure (16-sample overview, full pass)

This is the **second pass** of the `sample_overview_umap` card. A light,
perturbation-agnostic first pass was already produced in Tier 1B (see above);
Tier 2A re-emits the same card with Mixscape per-cell labels, feature-call
categories, and per-target MELD overlays populated — now that Tier 1C has run.
Two passes of the same card is intentional: the Tier 1B pass is the "is this
dataset worth drilling into" gate, and the Tier 2A pass is the analytical
product that collaborators will actually interact with.

Recommended composition (single composite card, full pass):

- a main cross-lineage integrated UMAP colored by Leiden cluster, with
  cluster labels in the margin
- a secondary coloring of the same UMAP by `lineage_label` (EBs vs iPSC)
- a **4 × 4 small-multiples band** with one mini-UMAP per sample on shared
  coordinates, each showing that sample's cells highlighted against a gray
  backdrop of all cells — this is the sample-reproducibility view
- a cluster × sample composition heatmap strip underneath for context
- a compact per-sample cell-count / QC table (n cells, % Mixscape-perturbed,
  median n_genes, n NTC cells)
- Mixscape-class overlay toggle (`KO` vs `NP` vs `NTC` vs `unassigned`)
- feature-call-category overlay toggle (`single` vs `multi` vs `ambiguous` vs
  `none`)

Biological questions this figure should answer at a glance:

- do EBs and iPSC separate cleanly on the shared manifold?
- are the EBs samples consistent with each other (differentiation reproducible)
  or does one sample dominate / drop out?
- are there Leiden clusters represented in only a subset of samples?
- is any Leiden cluster an obvious batch artifact (present in one sample only)?
- is perturbed-cell coverage uniform across samples, or biased toward a few?
- are Mixscape-`NP` escapers concentrated in particular clusters or spread
  evenly across the manifold?

This figure corresponds to the `sample_overview_umap` dynamic plot class (see
graphics compat runbook), which is specified to support both the Tier 1B light
pass and the Tier 2A full pass from a single manifest with different data
populations.

Useful follow-up summaries:

- cluster sizes
- cluster-by-sample composition (also serves as an integration sanity check)
- cluster-by-guide / cluster-by-target composition
- enrichment of particular perturbations in particular clusters

### Tier 2B: differential expression

Exploratory DE questions:

1. guide/target vs NTC across all filtered cells
2. guide/target vs NTC within a lineage
3. guide/target vs NTC within a Leiden cluster
4. cluster markers independent of perturbation identity

Recommended DE stack (pinned, not optional):

- **Exploratory layer — single-cell Wilcoxon.**
  Use `scanpy.tl.rank_genes_groups` (Wilcoxon) on target vs NTC, per lineage.
  Fast, useful for iteration, but known to overstate significance by ignoring
  pseudoreplication (Squair et al. 2021) — never the final reportable number.
- **Serious-claim layer — pseudobulk.**
  Aggregate counts by (sample, target) (or (sample, target, lineage)) and run
  **DESeq2** (or limma-voom) on the resulting matrix. This is the quotable DE
  number and the default for any claim leaving this group.
- **Perturb-seq-specific calibration layer — Sceptre.**
  For target-vs-NTC hypothesis testing with correct perturb-seq null
  calibration, run **Sceptre** (Barry et al. 2024) on the best-powered targets.
  Sceptre resamples NTC cells to build a correctly calibrated null and avoids
  the parametric-FDR inflation that trips up naive single-cell DE on CRISPR
  screens.

When Mixscape runs cleanly, DE contrasts at every layer default to
Mixscape-`perturbed` cells vs NTC. When Mixscape fit fails for a target, fall
back to all target-guide cells vs NTC and record `mixscape_status = failed`
on the DE record so the interpretation is explicit.

Recommended first-pass contrasts:

- per target gene:
  - all target-guide cells vs NTC (exploratory only)
  - **EBs Mixscape-perturbed vs EBs NTC (primary)**
  - **iPSC Mixscape-perturbed vs iPSC NTC (primary)**
  - optional: per-cluster target vs NTC inside interesting clusters

Recommended first-pass outputs:

- DE result tables per target (Wilcoxon + pseudobulk side-by-side)
- Sceptre-calibrated call set for the top-powered targets
- top-up / top-down gene heatmaps
- single-target volcano plots as drill-downs
- **volcano-matrix small-multiples** across selected targets as the DE summary
  figure
- target-by-gene effect matrix for the strongest, best-powered targets

Important caution:

- guide-level effects should be collapsed to target-level summaries before
  interpretation (see Tier 1C.1 guide/target posture)
- multiple guides against the same target are a consistency check, not noise to
  ignore

### Tier 2C: perturbation-state relationships

This is the analysis layer that Leiden makes possible. Cluster composition bars
are a starting point, not the end — they bake in a discretization (the Leiden
cut) that can swamp the signal. Complement them with two continuous methods.

Useful questions:

- do some perturbations shift cells into or out of specific states?
- do perturbations create new subpopulations or just move existing proportions?
- are effects global across the sample or concentrated in one state?
- for differentiating samples (EBs in particular), do perturbations shift
  developmental progression rather than steady-state identity?

Required methods:

- **Cluster-composition summaries.** Target-by-cluster composition tables,
  observed-vs-expected enrichment per cluster (chi-square / hypergeometric), and
  stacked barplots of cluster proportions by perturbation group. Useful as a
  discrete summary; known to be sensitive to the Leiden resolution.
- **Augur** (Skinnider et al. 2021). Trains a classifier to separate perturbed
  from control cells within each cell type / cluster and reports an AUC per
  cluster. This is strictly more informative than chi-square on cluster×target
  for ranking which states are responsive to which perturbations, because it
  uses the full expression vector, not cluster membership alone.
- **MELD** (Burkhardt et al. 2021). Computes a smooth "relative likelihood of
  perturbation" scalar field on the manifold, giving a *continuous* answer to
  "where in state space does perturbation X act" that does not depend on the
  Leiden cut. This is the correct display primitive for perturbation density on
  a UMAP.

Recommended for EBs specifically:

- **Pseudotime / trajectory analysis.** EBs are differentiating; the most
  biologically interesting perturbation effects in this dataset are likely
  shifts in developmental progression rather than steady-state cluster
  identity. Run a pseudotime inference (`scanpy.tl.dpt` for a quick pass,
  **CellRank 2** with a suitable kernel for the serious pass) and summarize
  pseudotime distribution per target vs NTC (KS test + violin plots).

Useful summaries:

- target-by-cluster composition tables (discrete)
- Augur AUC per (cluster, target) (continuous, classifier-based)
- MELD scalar field per target, rendered as a UMAP density overlay
- pseudotime distribution shift per target (EBs), with KS statistic and p-value
- UMAP overlays for selected high-confidence perturbations

## Tier 3: biological validation

Tier 3 asks whether the observed changes are coherent with prior knowledge.

This is not just "does the target go down?" It is:

- do downstream transcriptional effects align with known target biology?
- are pathway-level shifts directionally sensible?
- are related targets producing related signatures?
- do our effects agree with published perturb-seq atlases on the targets we
  share with them?

**Headline deliverable for Tier 3 is the target × target signature-correlation
heatmap** (the Replogle 2022 figure-4 analogue). Everything else in Tier 3
feeds into that figure or qualifies it. Pathway enrichments and per-target
validation notes are important but secondary.

### Tier 3A: target-level validation

For each high-confidence target:

- confirm target knockdown direction and magnitude (from Tier 1C.1)
- confirm per-cell perturbation fraction from Mixscape (Tier 1C.2)
- inspect canonical downstream markers
- inspect known pathway readouts
- compare multiple guides targeting the same gene
- compare the target-level signature to an external reference dataset where one
  exists (e.g., Replogle 2022 K562 essential-gene CRISPRi, Norman 2019, Jost
  2020) — see Tier 3D below

Examples of useful checks:

- transcription factor perturbation should affect its known regulon or lineage
  program
- pathway perturbation should affect expected pathway members, not random genes

### Tier 3B: pathway / program enrichment

Use **decoupler-py** (Badia-i-Mompel et al. 2022) as the default single-cell
pathway / TF activity layer so pathway, TF, and hallmark activities come out of
a single consistent framework rather than a patchwork of GSEA variants.

Recommended collections, in order of priority:

- **PROGENy** — 14 curated signaling pathway signatures, designed for
  single-cell activity inference
- **CollecTRI / DoRothEA** — TF–target regulons for TF activity scoring
  (essential for CRISPRi of transcription factors)
- **MSigDB Hallmark** — 50 Hallmark gene sets as a cross-study reference
- Reactome / KEGG / GO as supplementary, not primary

Run decoupler in two modes:

1. per-cell activity (applied to Mixscape-`perturbed` vs NTC distributions)
2. per-target activity (applied to ranked pseudobulk DE lists from Tier 2B)

Expected outputs:

- pathway / TF activity tables per target, per lineage
- compact pathway × target heatmap across selected perturbations
- TF × target heatmap for transcription-factor perturbations specifically
- target-to-pathway summary sheet

### Tier 3C: coherence / reproducibility checks

Questions:

- do replicate guides for the same target agree?
- do related targets cluster together by transcriptional effect?
- are EBs and iPSC effects shared or lineage-specific in a biologically
  sensible way?

Required analyses:

- guide-to-guide correlation within target (Spearman on pseudobulk log2FC
  vectors; already reported from Tier 1C.1)
- **target × target signature correlation matrix** on pseudobulk log2FC
  vectors across the top N highly variable / highly responsive genes; cluster
  hierarchically and annotate functional modules — this is the Tier 3 headline
  figure
- pathway concordance between lineages (EBs vs iPSC activity vectors per
  target)

### Tier 3D: external-atlas consistency

Perturb-seq datasets are rarely interpretable in isolation. For every target
that overlaps a published CRISPRi perturb-seq atlas in a comparable cell
context, compute and report a cross-dataset signature correlation.

Recommended reference datasets (union of what is practically available):

- Replogle et al. 2022 (genome-wide CRISPRi in K562 / RPE1)
- Norman et al. 2019 (dual-guide combinatorial CRISPRa in K562)
- Jost et al. 2020

Expected outputs:

- per-target Spearman correlation of our pseudobulk log2FC signature vs the
  external dataset's signature on shared genes
- a scatter plot of external log2FC vs UCSF log2FC for each high-confidence
  target
- a flag on any target where correlation is low despite adequate power in both
  datasets — that is either a genuine lineage-specific effect or an artifact,
  and downstream interpretation needs to say which

## Recommended immediate order

If starting from the current UCSF release surface, do the work in this order:

1. Finalize Tier 1A technical QC summary from the guarded rerun outputs.
2. Tier 1B — run a cross-lineage Harmony integration and emit the light
   `sample_overview_umap` card; evaluate the Tier 1B quantitative gate and
   produce a `tier1b_readiness.md` report.
3. Tier 1C.1 — run target-level and guide-level knockdown QC.
4. Tier 1C.2 — run Mixscape per target; record per-cell perturbation labels.
5. Tier 1C.3 — compute per-target E-distance vs NTC with permutation nulls.
6. Evaluate the Tier 1C quantitative gate; generate the Tier 1C readiness
   report.
7. Generate per-sample UMAP/Leiden diagnostics.
8. Generate **lineage-integrated** (Harmony or scVI) EBs and iPSC UMAP/Leiden
   as the analytical product.
9. Re-emit the `sample_overview_umap` card in full mode (Tier 2A) with
   Mixscape and feature-call overlays populated.
10. Tier 2C — run Augur and MELD on the integrated embeddings; run pseudotime
    for EBs.
11. Tier 2B — run Wilcoxon + pseudobulk DESeq2 + Sceptre contrasts for the
    best-powered targets.
12. Tier 3A / 3B — per-target validation notes and decoupler-py pathway/TF
    activity.
13. Tier 3C — target × target signature-correlation heatmap (headline figure).
14. Tier 3D — external-atlas consistency check for overlapping targets.

## Concrete first iteration

The lowest-effort, highest-value first iteration is:

1. Run a cross-lineage Harmony integration on all 16 samples' filtered counts
   and emit the light Tier 1B `sample_overview_umap` card, producing:
   - integrated `.h5ad` (first-pass, cross-lineage)
   - Tier 1B composite PNG (pooled UMAP + small-multiples + composition strip)
   - `sample_summary.tsv` (per-sample cell counts, n_NTC, median_n_genes,
     median_mt_pct, dominant Leiden cluster)
   - `cluster_by_sample.tsv`
   - integration-mixing metrics table (sample-entropy per cluster, silhouette
     of EBs vs iPSC)
   - `tier1b_readiness.md` stating pass/fail against the Tier 1B criteria
2. Run `scripts/perturbseq_qc_dotplot.py` on the guarded rerun surface to
   produce Tier 1C.1 outputs:
   - `cell_count_distribution.png`
   - `power_ranked_summary.tsv`
   - `knockdown_dotplot.png`
   - `knockdown_expressed_only.png`
3. Add Mixscape (Tier 1C.2) and E-distance (Tier 1C.3) passes, producing:
   - per-cell `mixscape_class` / `mixscape_score` columns appended to the
     filtered AnnData
   - `target_mixscape_summary.tsv` (fraction perturbed per target)
   - `target_edistance_summary.tsv` (E-distance, permutation p, rank)
   - `guide_concordance.tsv` (guide-to-guide Spearman per target)
4. Evaluate the Tier 1C quantitative gate and emit a one-page
   `tier1c_readiness.md` stating pass/fail with the three criteria.
5. Run a per-sample Scanpy workflow for `PCA -> neighbors -> UMAP -> Leiden` on
   one sample (recommended: `EBs2_2`) as an integration-independent diagnostic.
6. Run a **lineage-pooled integrated** Scanpy + Harmony workflow for EBs,
   producing:
   - integrated `.h5ad`
   - UMAP colored by sample, Leiden, `feature_call_category`, `mixscape_class`
   - cluster-by-target composition table
   - MELD density overlay per top-5 targets (PNG + payload)
7. Re-emit the `sample_overview_umap` card in full mode (Tier 2A) with
   Mixscape and feature-call overlays populated on the integrated embedding.
8. Pick one clean, high-confidence target (e.g., `PAX6` in pooled EBs) and run
   one full DE cascade (Wilcoxon + pseudobulk DESeq2 + Sceptre) to produce
   one real volcano and one real DE table.

That gives a usable Tier 1B + Tier 1C + Tier 2A (both passes of the orientation
card) + one real DE example without committing to full pathway / cross-dataset
validation yet.

## Graphics packaging

If we want the UCSF outputs to drop into the MorPhiC gene-page surface, the
figure contract should be frozen now rather than after plots start to
accumulate.

See:

- `docs/RUNBOOK_UCSF_MORPHIC_GRAPHICS_COMPAT_20260417.md`

If the dynamic target is the Shiny / Plotly surface, the first contract-forcing
plot set is five plots:

1. `sample_overview_umap` — dataset-orientation summary; pooled integrated
   UMAP with a per-sample small-multiples band (the 16-sample summary slide)
2. `umap_scatter` — cell-level embedding data (with MELD density overlay mode)
3. `target_qc_dotplot` — target-by-sample aggregate QC (with Mixscape and
   E-distance fields)
4. `de_volcano` — gene-level DE (single-target and volcano-matrix modes)
5. `target_correlation_heatmap` — target-by-target signature correlation (the
   Tier 3 headline)

The first plot is the "open this first" orientation figure. The fifth plot is
what forces the schema to handle target × target grain, hierarchical ordering,
and on-demand gene-list payloads at v1 instead of retrofitting them later.

## What Leiden is and is not

Leiden is recommended here, but only in the right place.

Use Leiden to answer:

- what cell states are present?
- where do perturbations accumulate?
- are there state-specific perturbation effects?

Do not use Leiden as the primary answer to:

- did the perturbation work?

That answer comes first from target-expression and guide-level efficacy checks.

## Decisions closed in this revision

These were previously listed as open. They are now closed and should be
followed unless explicitly reopened:

- **Guide vs target-level Tier 1C summaries** — report both. Target-level is
  the primary efficacy call; guide-level is the consistency overlay.
- **Pooled exploratory scope** — lineage-integrated (EBs, iPSC) with Harmony
  or scVI is the analytical product; per-sample is a diagnostic only.
- **DE framework** — exploratory Wilcoxon + serious-claim pseudobulk DESeq2 +
  perturb-seq-calibrated Sceptre for top-powered targets. This is a stack, not
  a choice.
- **Default pathway collections** — PROGENy + CollecTRI/DoRothEA + MSigDB
  Hallmark via decoupler-py, with Reactome/KEGG/GO as supplementary.

## Open questions for next iteration

- what exact thresholds should Tier 1B use for integration-mixing and
  lineage-separation metrics (sample-entropy per cluster, EBs-vs-iPSC
  silhouette)? The structural gate is prespecified; the numeric thresholds
  still need a first empirical pass.
- which external atlases (Replogle 2022 vs Norman 2019 vs Jost 2020) should be
  the required Tier 3D reference set, and which are optional?
- should Mixscape failure automatically downgrade a target's DE record, or
  should we retain both Mixscape-filtered and unfiltered DE pairs by default?
- should the integrated EBs embedding use Harmony (fast, simple) or scVI
  (model-based, enables latent-space DE) as the primary; or both?
- should the first dynamic Shiny view read parquet directly or go through a
  thinner derived JSON layer?
- how much of the MorPhiC card metadata should be generated by the analysis
  pipeline vs a separate packaging step?

## Current recommendation

Current default recommendation:

- Tier 1 is a **three-part** staged gate:
  - Tier 1A — per-sample technical QC
  - Tier 1B — dataset-level tenor (perturbation-agnostic integration pass;
    emits the light `sample_overview_umap` summary card)
  - Tier 1C — perturbation-efficacy QC: target-level knockdown + Mixscape
    per-cell labels + E-distance in representation space, with a quantitative
    pass threshold
- Tier 2 embeddings use **lineage-integrated** Harmony/scVI as the analytical
  product; per-sample embeddings are diagnostic only
- Tier 2C goes beyond cluster composition: add Augur, MELD, and (for EBs)
  pseudotime
- Tier 2B DE is a three-layer stack (Wilcoxon / pseudobulk DESeq2 / Sceptre),
  not a single framework choice
- Tier 3 headline is the target × target signature-correlation heatmap; pathway
  activity comes through decoupler-py; external-atlas consistency is required
  for overlapping targets

That gives a disciplined sequence:

- first prove the perturbations are real at both cell-level and signature-level
- then map their effects onto an integrated cell-state manifold
- then test whether the resulting biology is coherent internally and
  externally
