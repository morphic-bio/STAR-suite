# Handoff: UCSF Perturb-seq Analysis Direction (2026-04-17)

## Goal

Summarize the current state of the UCSF perturb-seq analysis planning work so a
coding agent can continue from here without re-deriving the decisions.

This handoff covers:

- what has been reviewed
- what has already been written
- the recommended analysis sequence
- the graphics / product contract direction
- the next concrete implementation steps

## Current objective

We are moving from:

- release packaging
- simple cell-level QC

to:

- perturbation-aware QC
- exploratory perturb-seq analysis
- biological validation
- graphics that are compatible with the current MorPhiC site while also
  pushing toward a Plotly / Shiny future

The immediate product goal is **not** to build a full polished portal yet. The
goal is to define the real contracts by instantiating a small number of useful
plots and result bundles.

## What has been reviewed

### UCSF production / guarded rerun surface

Authoritative analysis surface:

- `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009`

Important repo references already reviewed:

- `docs/HANDOFF_UCSF_DATA_ORGANIZATION_AND_RELEASE_NOTES_20260331.md`
- `docs/PAPER_UCSF_CORRECTED_PRODUCTION_WORKFLOW_20260329.md`
- `docs/HANDOFF_UCSF_PRODUCTION_RUN_20260329.md`
- `docs/HANDOFF_UCSF_ADAPTIVE_DOWNSTREAM_REPAIR_20260331.md`
- `scripts/ucsf_release_guarded_README.md`
- `scripts/ucsf_release_README.md`

Important findings:

- default exploratory input should be `filtered_counts.h5ad`
- `default_singlet_filtered_counts.h5ad` is a stricter option
- `final_counts.h5ad` should **not** be the default analysis input
- no saved PCA / UMAP / Leiden outputs were found in the UCSF production
  surface
- `filtered_counts.h5ad` already carries the key cell metadata needed for
  downstream perturb-seq analysis

### Existing repo capabilities

Relevant existing scripts / precedents:

- `scripts/perturbseq_qc_dotplot.py`
  - target-guide vs NTC perturbation-efficacy QC
- `scripts/generate_qc_histogram.py`
  - Plotly precedent in Python
- `scripts/paper/project_flex_leiden_to_cr9.py`
  - Scanpy precedent for `PCA -> neighbors -> UMAP -> Leiden`
- `docs/HEATMAP_REFACTOR_SUMMARY.md`
  - existing repo precedent for Plotly HTML + JSON outputs

### MorPhiC public surface

Reviewed on `2026-04-17`:

- `https://morphic.bio/genes/HGNC:8620/`
- `https://morphic.bio/methods/`
- MorPhiC client JS bundle for the gene page
- MorPhiC gene API for `HGNC:8620`

Important findings:

- the current public MorPhiC gene page is **static-first**
- the page renders from API records like `Analysis_Results` and
  `Enrichment_Analysis`
- current live cards use:
  - `title`
  - optional `s3_tsv_key`
  - either `svg` or `s3_png_key`
  - `analysis_context`
- current public bundles inspected contained **no Plotly references**
- there is no public Plotly contract on the current gene page
- however, the intended direction is toward dynamic Plotly / Shiny views

## Documents and artifacts created

### 1. Main analysis runbook

- `docs/RUNBOOK_UCSF_PERTURB_ANALYSIS_20260417.md`

Purpose:

- defines the overall analysis program
- separates technical QC from perturbation-efficacy QC
- places Leiden/UMAP in the correct stage

Key structure:

1. Tier 1: QC / release gate
2. Tier 2: exploratory analysis
3. Tier 3: biological validation

### 2. Graphics compatibility / future contract runbook

- `docs/RUNBOOK_UCSF_MORPHIC_GRAPHICS_COMPAT_20260417.md`

Purpose:

- defines how UCSF outputs should map onto MorPhiC-style cards
- documents current MorPhiC behavior
- proposes the dynamic Plotly / Shiny manifest contract

### 3. Example dynamic manifest files

- `docs/examples/ucsf_plot_contracts/sample_overview_umap.manifest.example.json`
- `docs/examples/ucsf_plot_contracts/umap_scatter.manifest.example.json`
- `docs/examples/ucsf_plot_contracts/target_qc_dotplot.manifest.example.json`
- `docs/examples/ucsf_plot_contracts/de_volcano.manifest.example.json`
- `docs/examples/ucsf_plot_contracts/target_correlation_heatmap.manifest.example.json`

Purpose:

- make the proposed dynamic contract concrete
- define real schema shape for the first five dynamic plot classes:
  composite-layout dataset orientation, cell-level exploration,
  target-by-sample aggregate, gene-level DE, and target-by-target correlation

### 4. Card mockup page

- `docs/examples/ucsf_plot_contracts/mockup_cards.html`

Purpose:

- local standalone visual mockup of the proposed MorPhiC-style cards
- shows the first three card types
- includes a contract inspector pane showing manifest fields and data grain

## Recommended overall analysis direction

The current recommended direction is:

1. prove the perturbations worked
2. then do clustering / embedding / exploratory analysis
3. then run DE and biological validation
4. package outputs in a card-oriented way that works for both static MorPhiC
   pages and future Plotly / Shiny views

This is intentionally disciplined. A perturb-seq dataset is not analysis-ready
just because generic scRNA QC looks fine.

## Analysis runbook summary

### Tier 1: QC / release gate

Tier 1 is split into:

1. technical QC
2. perturbation-efficacy QC

#### Tier 1A technical QC

Core checks:

- cell counts per sample
- genes / UMIs / mt distributions
- doublet rate
- guide-call category distributions
- guide coverage and target abundance

#### Tier 1B dataset-level tenor (perturbation-agnostic)

Core question:

- does the 16-sample collection hang together as a single coherent perturb-seq
  dataset, before we commit compute to per-target efficacy analysis?

Required compute:

- first-pass cross-lineage Harmony integration on filtered counts
- `sample_overview_umap` card in light mode (no Mixscape overlays) as the
  Tier 1B summary slide

Quantitative pass gate (prespecified; exact thresholds open):

- \u22656000 cells per sample after per-sample QC
- \u22651000 NTC cells per lineage
- no Leiden cluster >70% dominated by a single sample, or such clusters are
  annotated batch-suspect before Tier 1C runs against them
- EBs and iPSC are clearly separable on the integrated UMAP

Expected deliverables:

- cross-lineage integrated `.h5ad`
- Tier 1B composite PNG
- `sample_summary.tsv`, `cluster_by_sample.tsv`
- integration-mixing metrics table
- `tier1b_readiness.md`

#### Tier 1C perturbation-efficacy QC

Core question:

- for each target, do target-guide cells actually show reduced target
  expression relative to NTC **and** contain a population of truly perturbed
  cells **and** move meaningfully in representation space?

This is the first real perturb-seq gate, and it is a *three-part* gate:

1. **Tier 1C.1 — target-level mean knockdown.** Guide-level and target-level
   log2FC vs NTC, per sample and per lineage, with guide-to-guide Spearman
   concordance as the consistency overlay. Target-level is the primary
   efficacy call.
2. **Tier 1C.2 — Mixscape per-cell labeling** (Papalexi 2021). Fit a two-component
   mixture per target and report fraction of `perturbed` (`KO`) vs `non_perturbed`
   (`NP`) cells. Escapers are excluded from downstream effect estimation when
   Mixscape fits; Tier 2 DE defaults to Mixscape-`perturbed` vs NTC.
3. **Tier 1C.3 — E-distance** (scPerturb, Peidli 2023) between target-guide
   cells and NTC cells in PCA space, with NTC-vs-NTC permutation null. Catches
   regulators with strong downstream effect but no target self-knockdown.

Quantitative pass gate (prespecified):

- \u226550% of evaluable, baseline-expressed targets show target-level log2FC
  \u2264 \u22120.5 at FDR < 0.05 in their primary lineage
- \u226550% of evaluable targets show E-distance permutation p < 0.05 vs NTC
- \u226560% median guide-to-guide Spearman concordance for multi-guide targets

Recommended first output:

- run `scripts/perturbseq_qc_dotplot.py` (Tier 1C.1)
- add Mixscape and E-distance passes (Tier 1C.2 / 1C.3)

Expected deliverables:

- `cell_count_distribution.png`
- `power_ranked_summary.tsv`
- `knockdown_dotplot.png`
- `knockdown_expressed_only.png`
- `target_mixscape_summary.tsv`
- `target_edistance_summary.tsv`
- `guide_concordance.tsv`
- `tier1c_readiness.md` (one-page pass/fail report against the quantitative gate)

### Tier 2: exploratory analysis

Per-sample embeddings are a QC diagnostic; the analytical product is a
**lineage-integrated** embedding.

Recommended order:

1. per-sample `PCA -> neighbors -> UMAP -> Leiden` (diagnostic)
2. lineage-pooled integrated EBs (Harmony or scVI; analytical product)
3. lineage-pooled integrated iPSC
4. cross-lineage integrated (all 16 samples) only for cross-lineage questions

Tier 2C adds continuous perturbation-state methods on top of cluster
composition:

- **Augur** (Skinnider 2021) — classifier AUC per (cluster, target) for
  ranking responsive states
- **MELD** (Burkhardt 2021) — smooth perturbation-likelihood scalar field on
  the manifold; correct primitive for UMAP density overlays
- **Pseudotime / CellRank 2** — required for EBs, where perturbation-driven
  shifts in developmental progression are often the biology

Expected outputs:

- per-sample diagnostic UMAP PNGs
- integrated `.h5ad` per lineage
- cluster tables, cluster-by-target composition
- Augur AUC tables, MELD scalar overlays, pseudotime violin/KS summaries
- dynamic UMAP card bundles with MELD overlay mode

### Tier 2B differential expression

Pinned three-layer DE stack (not a framework choice):

1. **Exploratory — single-cell Wilcoxon** (`scanpy.tl.rank_genes_groups`) per
   lineage. Fast, but overstates significance due to pseudoreplication
   (Squair 2021); never the final reportable number.
2. **Serious-claim — pseudobulk DESeq2** (or limma-voom) aggregated by
   (sample, target, lineage). Default for any quotable DE result.
3. **Perturb-seq calibration — Sceptre** (Barry 2024) on top-powered targets
   for correctly calibrated null distributions on CRISPR-screen contrasts.

When Mixscape fits cleanly, DE defaults to Mixscape-`perturbed` vs NTC at
every layer; failed fits fall back to raw target-guide vs NTC with
`mixscape_status = failed` recorded on the DE record.

First-pass contrasts:

- all target-guide cells vs NTC (exploratory only)
- EBs Mixscape-perturbed vs EBs NTC (primary)
- iPSC Mixscape-perturbed vs iPSC NTC (primary)
- optional per-cluster target vs NTC inside interesting clusters

### Tier 3: biological validation

Headline deliverable is the **target \u00d7 target signature-correlation heatmap**
(Replogle 2022 figure-4 analogue). Pathway enrichments and per-target
validation notes are secondary.

Questions:

- do downstream signatures make sense biologically?
- do related targets produce coherent effects and cluster into recognizable
  functional modules?
- are pathway shifts plausible?
- do multiple guides for the same target agree?
- do target signatures agree with published perturb-seq atlases (Replogle
  2022, Norman 2019, Jost 2020) for shared targets? (Tier 3D)

Default pathway / TF activity stack:

- **decoupler-py** (Badia-i-Mompel 2022) as the single framework
- PROGENy (signaling), CollecTRI/DoRothEA (TFs), MSigDB Hallmark as the
  required collections; Reactome / KEGG / GO as supplementary

Expected outputs:

- target \u00d7 target signature-correlation heatmap (headline)
- decoupler-py pathway \u00d7 target and TF \u00d7 target activity heatmaps
- target-level validation notes for the best-supported hits
- external-atlas consistency scatter and flags for overlapping targets

## Graphics / product direction

### Current MorPhiC compatibility target

For immediate compatibility, each result should still be expressible as a
MorPhiC-style static card:

- `title`
- `study_label`
- `dataset_id`
- `result_type`
- `role`
- `analysis_context`
- `s3_png_key` or `svg`
- optional `s3_tsv_key`

For DE cards, preserve current MorPhiC-style fields where relevant:

- `de_summary`
- `de_thresholds`
- `de_conditions`
- `default_condition_id`

### Dynamic target architecture

The dynamic contract should be:

1. card record
2. interactive manifest
3. data payloads
4. static fallback

The dynamic unit should **not** be raw Plotly traces alone.

Why:

- Shiny needs controls, defaults, and data-grain awareness
- versioning is needed
- provenance is needed
- not every consumer will directly render Plotly

### Proposed schema version

- `morphic-plotly-shiny-v1`

### Proposed first five dynamic plot classes

These were chosen specifically to force the major schema decisions. The set is
five plots. The first (`sample_overview_umap`) is the "open this first"
dataset-orientation card and forces composite layouts with facet + highlight
semantics. The last (`target_correlation_heatmap`) carries the target \u00d7 target
grain that is the biological headline of the dataset and forces hierarchical
ordering, interactive reordering, and on-demand gene-list payloads.

1. `sample_overview_umap`
   - dataset-orientation summary card (featured role)
   - composite layout: pooled integrated UMAP + per-sample small-multiples
     band (4 \u00d7 4) + cluster \u00d7 sample composition strip + per-sample QC table
   - layout modes: `composite`, `pooled_only`, `small_multiples_only`,
     `lineage_split`
   - sample-highlight control (foreground vs gray backdrop)
   - answers at a glance: do EBs/iPSC separate, are samples reproducible, is
     any cluster a batch artifact, is perturbed-cell coverage uniform
2. `umap_scatter`
   - cell-level grain (exploration, not orientation)
   - embeddings, overlays, hover fields, filters
   - includes MELD density overlay mode and continuous gene-expression lookup
3. `target_qc_dotplot`
   - aggregate target-by-sample grain
   - size + color encoding, thresholds, evaluability semantics
   - carries Mixscape fraction-perturbed and E-distance fields alongside
     mean-level log2FC
4. `de_volcano`
   - gene-level grain
   - thresholds, labels, table linkage
   - supports single-target and volcano-matrix small-multiples modes; records
     `de_method` (Wilcoxon / pseudobulk_deseq2 / sceptre) and `mixscape_status`
5. `target_correlation_heatmap`
   - target \u00d7 target grain
   - hierarchical ordering, interactive reordering, on-demand gene-list
     payloads per pair
   - the Tier 3 headline figure

These five are sufficient to define the main data grains and layout semantics
needed by the UCSF analysis product:

- composite-layout dataset orientation (with sample faceting and highlight)
- cell-level exploration (with continuous-scalar overlay semantics)
- aggregate perturbation QC (with controlled evaluability vocabulary)
- gene-level DE (with multi-method and small-multiples semantics)
- target \u00d7 target correlation (with on-demand gene-list payloads)

## Recommended implementation order from here

### Phase 1: make Tier 1B real (dataset tenor, light pass)

1. run a cross-lineage Harmony integration on all 16 samples' filtered counts
2. produce:
   - cross-lineage integrated `.h5ad` (first-pass, not analytical product)
   - per-sample summary TSV (n_cells, n_NTC, median_n_genes, median_mt_pct,
     dominant Leiden cluster)
   - cluster \u00d7 sample composition TSV
   - integration-mixing metrics (sample-entropy per cluster, EBs-vs-iPSC
     silhouette)
3. instantiate the first `sample_overview_umap` manifest + payload bundle in
   `dataset_tenor_light` mode (no Mixscape overlays populated)
4. emit `tier1b_readiness.md` with pass/fail against the Tier 1B quantitative
   gate
5. package all Tier 1B outputs as static card records (composite card +
   readiness report)

### Phase 2: make Tier 1C real (three-part efficacy gate)

1. run `scripts/perturbseq_qc_dotplot.py` (Tier 1C.1)
2. add a Mixscape pass (Tier 1C.2), producing per-cell `mixscape_class` /
   `mixscape_score` and a `target_mixscape_summary.tsv`
3. add an E-distance pass (Tier 1C.3), producing `target_edistance_summary.tsv`
4. add `guide_concordance.tsv` from pseudobulk guide-level log2FC vectors
5. emit `tier1c_readiness.md` with pass/fail against the three-part quantitative
   gate
6. package all Tier 1C outputs as static card records (target dotplot card +
   guide strip plot card + readiness report)

### Phase 3: make the first embedding real (including Tier 2A full pass)

1. build a per-sample Scanpy workflow as a diagnostic (start with `EBs2_2`)
2. build a lineage-pooled **integrated** (Harmony) Scanpy workflow for EBs as
   the analytical product
3. generate:
   - integrated per-lineage UMAP coordinates
   - Leiden labels
   - Mixscape class overlay on the integrated embedding
   - MELD scalar overlay for top-5 targets
   - Augur AUC per (cluster, target)
   - per-sample summary table (n_cells, % Mixscape-perturbed, median n_genes)
   - cluster \u00d7 sample composition table
   - a static PNG set
   - cluster-by-target composition table
4. **re-emit** the `sample_overview_umap` card in `dataset_orientation_full`
   mode on the cross-lineage embedding, now with Mixscape and feature-call
   overlays populated (second pass of the Phase 1 card)
5. instantiate the `umap_scatter` manifest bundle for the EBs exploration card
   with MELD overlay mode

### Phase 4: make one DE card real (full stack)

1. choose one clean, high-confidence target-vs-NTC comparison (e.g., `PAX6` in
   pooled EBs)
2. run the full DE cascade: Wilcoxon + pseudobulk DESeq2 + Sceptre
3. generate a real DE table (all three methods, side-by-side)
4. generate a static single-target volcano and a volcano-matrix PNG across
   the top-powered targets
5. instantiate the first `de_volcano` manifest bundle with `de_method` and
   `mixscape_status` fields populated

### Phase 5: make the Tier 3 headline real

1. build the per-target pseudobulk log2FC vectors across the top 2000 HVGs
2. compute the target \u00d7 target Spearman / Pearson correlation matrix and
   hierarchically cluster to identify functional modules
3. generate a static clustered-heatmap PNG
4. instantiate the first `target_correlation_heatmap` manifest bundle with
   per-pair gene-list payloads resolved on-demand
5. run decoupler-py (PROGENy + CollecTRI + Hallmark) and produce a pathway \u00d7
   target heatmap card
6. run Tier 3D external-atlas consistency for overlapping targets

### Phase 6: connect the five together

Once one real example exists for each of the five dynamic classes:

- the Shiny implementation can be built against a concrete contract
- the MorPhiC-style static cards can be populated with real assets
- the `sample_overview_umap` card becomes the featured landing card on the
  UCSF study page; drill-downs flow from it into `umap_scatter`,
  `target_qc_dotplot`, `de_volcano`, and `target_correlation_heatmap`
- coordinated linking (UMAP \u2194 DE table \u2194 guide strip; dotplot drill-down;
  correlation-heatmap drill-down) can be implemented against working data
- future plot classes can be added as extensions rather than new schema
  categories

## Immediate instructions for a coding agent

If assigning this to another coding agent, the recommended immediate tasks are:

1. Read:
   - `docs/RUNBOOK_UCSF_PERTURB_ANALYSIS_20260417.md`
   - `docs/RUNBOOK_UCSF_MORPHIC_GRAPHICS_COMPAT_20260417.md`
2. Inspect:
   - `docs/examples/ucsf_plot_contracts/mockup_cards.html`
   - the five example manifest JSON files (including
     `sample_overview_umap.manifest.example.json` and
     `target_correlation_heatmap.manifest.example.json`)
3. Treat `filtered_counts.h5ad` as the default analysis input.
4. Do **not** use `final_counts.h5ad` as the default exploratory surface.
5. Implement Tier 1B first (dataset-level tenor, perturbation-agnostic):
   - run a cross-lineage Harmony integration on all 16 samples
   - emit the light `sample_overview_umap` card (`dataset_tenor_light` role)
   - compute integration-mixing metrics and lineage separation
   - emit `tier1b_readiness.md` against the Tier 1B quantitative gate
     (per-sample cell counts, NTC pool, cluster sample-entropy, lineage
     separability)
6. Then implement Tier 1C as a three-part efficacy gate:
   - Tier 1C.1: run / adapt `scripts/perturbseq_qc_dotplot.py`
   - Tier 1C.2: add Mixscape per-cell labeling (Papalexi 2021)
   - Tier 1C.3: add scPerturb E-distance with NTC-vs-NTC permutation null
   - emit a `tier1c_readiness.md` against the prespecified quantitative gate
     (\u226550% targets log2FC \u2264 \u22120.5 at FDR < 0.05; \u226550% E-distance p < 0.05;
     \u226560% guide concordance)
7. Then implement one **lineage-integrated** (Harmony) Scanpy workflow for
   pooled EBs:
   - PCA
   - neighbors
   - UMAP
   - Leiden
   - MELD scalar overlay per top-5 targets
   - Augur AUC per (cluster, target)
   - re-emit the `sample_overview_umap` card in `dataset_orientation_full` mode
     with Mixscape and feature-call overlays populated
8. Then pick one clean target (e.g., `PAX6` in pooled EBs) and run the full DE
   cascade: Wilcoxon + pseudobulk DESeq2 + Sceptre.
9. Use the `sample_overview_umap`, `umap_scatter`, `de_volcano`, and
   `target_correlation_heatmap` manifest examples as the target contract for
   the first real dynamic bundles.

## Open decisions to keep explicit

These are still open and should remain explicit rather than getting decided
accidentally in code:

- which external atlases (Replogle 2022 vs Norman 2019 vs Jost 2020) should be
  the required Tier 3D reference set, and which are optional
- whether Mixscape failure should automatically downgrade a target's DE record
  or whether we retain both Mixscape-filtered and unfiltered DE pairs by default
- whether the integrated EBs embedding should use Harmony (fast, simple) or
  scVI (model-based, latent-space DE) as primary, or both
- whether the first dynamic Shiny view should read parquet directly or use a
  thinner derived JSON layer
- how much of the MorPhiC card metadata should be generated by the analysis
  pipeline vs a separate packaging step

### Decisions now closed (do not revisit without explicit reason)

- guide-level vs target-level Tier 1C: **both** — target-level primary,
  guide-level consistency overlay
- pooled exploratory scope: **lineage-integrated (Harmony/scVI) as the
  analytical product**; per-sample as diagnostic only
- default DE framework: **three-layer stack** (Wilcoxon / pseudobulk DESeq2 /
  Sceptre), not a single choice
- default pathway collections: **decoupler-py with PROGENy + CollecTRI +
  Hallmark**; Reactome/KEGG/GO supplementary

## Bottom line

The project direction is now:

- **analysis first, not just graphics**
- **Tier 1 is a three-part staged gate**: 1A per-sample technical, 1B
  dataset-level tenor (perturbation-agnostic integration pass that emits the
  light orientation card), 1C per-target efficacy
- **perturbation efficacy as a three-part gate** (target-mean + Mixscape +
  E-distance) with a prespecified quantitative threshold, **plus a
  CRISPRoff-specific temporal-control / selection-bias readout** based on
  guide / target representation
- **lineage-integrated embeddings (Harmony / scVI) as the analytical product**;
  per-sample as diagnostic
- **continuous perturbation-state methods** (Augur, MELD, pseudotime) layered
  on top of cluster composition
- **three-layer DE stack** (Wilcoxon / pseudobulk DESeq2 / Sceptre) on
  Mixscape-perturbed cells
- **target \u00d7 target signature-correlation heatmap as the Tier 3 headline**;
  decoupler-py for pathway / TF activity; external-atlas consistency required
- **Plotly / Shiny introduced through a manifest contract**
- **five concrete dynamic plot classes** used to force the contract
  (`sample_overview_umap`, `umap_scatter`, `target_qc_dotplot`, `de_volcano`,
  `target_correlation_heatmap`); the first is the featured dataset-orientation
  card that summarizes all 16 samples on one composite view
- **static MorPhiC-compatible cards retained throughout**

The next best step is to stop writing abstract schema prose and instantiate one
real UCSF sample bundle, starting with the Tier 1B light `sample_overview_umap`
card, followed by the three-part Tier 1C efficacy output, the Tier 2A full-pass
orientation card, one real DE cascade, and one real target-correlation heatmap.
