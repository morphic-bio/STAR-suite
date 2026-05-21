# Runbook: UCSF MorPhiC Graphics Compatibility (2026-04-17)

## Goal

Define a graphics/result contract for UCSF perturb-seq outputs that:

- drops cleanly into the current MorPhiC gene-page surface
- preserves static fallbacks
- leaves a clean upgrade path for interactive Plotly views

This is intentionally conservative. The current MorPhiC page is static-first.
We should match that now, then extend it without breaking the existing card
model.

## What MorPhiC currently does

Observed on `2026-04-17` from:

- `https://morphic.bio/genes/HGNC:8620/`
- `https://morphic.bio/methods/`
- client bundle:
  `https://morphic.bio/component---src-pages-genes-gene-id-js-ec0898c96d7614318fe8.js`
- gene API:
  `https://46ucfedadd.execute-api.us-east-1.amazonaws.com/api/gene/HGNC%3A8620`

Important findings:

- the gene page HTML itself is not the authoritative payload; the page is
  populated client-side from the API
- the public site bundles we inspected contain **no Plotly references**
- the gene page renders analysis cards from `Analysis_Results`
- the gene page renders enrichment cards from `Enrichment_Analysis`
- for each analysis card, the frontend behavior is:
  - render `title`
  - if `s3_tsv_key` exists, render a TSV download link
  - if `svg` exists, render the inline SVG
  - else if `s3_png_key` exists, render the PNG
  - else render "No image available"
- `analysis_context` is already used as metadata and is the right place for
  lineage / assay / target annotations
- DE-style cards can additionally carry:
  - `de_summary`
  - `de_thresholds`
  - `de_conditions`
  - `default_condition_id`
- a current UMAP-style card exists as a PNG-only result item:
  - title: `MSK Pooled scRNA-seq: PAX6_umap`
  - metadata includes `analysis_context.plot_kind = "umap"`
- the current CSS uses a white image card style and a nominal card width of
  about `400px`
- one live MorPhiC UMAP PNG is `1920 x 1280` pixels, which is a good static
  export target

### Immediate consequence

The compatibility target today is:

- static PNG or inline SVG
- optional TSV download
- structured metadata

Plotly is **not** currently the public MorPhiC gene-page standard.

## Design principle

Treat every figure as two layers:

1. a **card record** that the gene page can list and render
2. a **dynamic manifest** that Shiny / Plotly can open

Do **not** make the dynamic manifest the only representation. The static card
must remain sufficient on its own.

## Target architecture for Plotly / Shiny

The dynamic target should not be "raw Plotly JSON somewhere in S3".

That is too renderer-specific and does not give Shiny enough structure to build
controls, filters, downloads, or linked tables.

The better contract is a four-part bundle:

1. **card record**
   - what the MorPhiC gene page indexes and renders
2. **interactive manifest**
   - what the Shiny app reads to understand plot type, datasets, encodings,
     controls, defaults, and downloads
3. **data payloads**
   - the actual cell-level or aggregate tables, ideally in parquet
4. **static fallback**
   - the PNG or SVG that works even when no dynamic viewer is present

This keeps the dynamic view declarative without pretending the Plotly trace list
is itself the long-term product contract.

## Recommended UCSF result families

Use a small controlled vocabulary so the frontend can group or filter results
later without parsing free text.

| Analysis tier | `result_type` | Typical `role` values |
|---|---|---|
| Tier 1A technical QC | `QC` | `technical_qc`, `guide_call_qc`, `cell_count_qc` |
| Tier 1B perturbation QC | `PERTURB_QC` | `knockdown_qc`, `expressed_target_qc` |
| Tier 2 embeddings | `EMBEDDING` | `canonical_umap`, `leiden_clusters`, `target_overlay`, `qc_overlay` |
| Tier 2 compositions | `COMPOSITION` | `cluster_by_sample`, `cluster_by_target`, `cluster_enrichment` |
| Tier 2 differential expression | `DE` | `canonical_de`, `cluster_specific_de` |
| Tier 3 validation / pathways | `ENRICHMENT` | `pathway_enrichment`, `marker_validation`, `program_validation` |

The exact strings are less important than keeping them stable once adopted.

## Static card contract

This is the minimum contract I would target for UCSF outputs.

```json
{
  "title": "UCSF CRISPRi EBs2_2: UMAP colored by Leiden",
  "study_label": "UCSF_CRISPRi_PerturbSeq",
  "study_id": "ucsf_crispri_perturbseq_v1",
  "dataset_id": "EBs2_2",
  "result_type": "EMBEDDING",
  "role": "leiden_clusters",
  "s3_png_key": "png/ucsf/ebs2_2_leiden_umap.png",
  "s3_tsv_key": "tsv/ucsf/ebs2_2_leiden_cluster_counts.tsv.gz",
  "analysis_context": {
    "plot_kind": "umap",
    "lineage_label": "Embryoid body",
    "sample_label": "EBs2_2",
    "study_assay": "Perturb-seq",
    "target_gene_symbols": [],
    "cluster_key": "leiden",
    "embedding_key": "X_umap",
    "n_cells": 10897
  },
  "summary": {
    "n_cells": 10897,
    "n_clusters": 14
  }
}
```

### Required fields

- `title`
- `study_label`
- `result_type`
- `role`
- `analysis_context`
- one of:
  - `svg`
  - `s3_png_key`

### Strongly recommended fields

- `dataset_id`
- `s3_tsv_key` when a natural table exists
- `summary` for compact card metadata

### DE-specific fields

For DE cards, follow the current MorPhiC shape rather than inventing a new one.

Recommended DE extras:

- `de_summary`
- `de_thresholds`
- `de_conditions`
- `default_condition_id`

This lets UCSF DE outputs fit the same pattern as existing MorPhiC bulk DE
cards.

## Proposed Plotly / Shiny contract

Because MorPhiC does not yet expose a public Plotly contract, the right move is
to define a viewer-oriented manifest now instead of pushing raw Plotly state
directly into the card record.

### Card-level extension

Add this only when an interactive view exists:

```json
{
  "dynamic_plot": {
    "schema_version": "morphic-plotly-shiny-v1",
    "viewer": "plotly_shiny",
    "manifest_key": "plotly/ucsf/ebs2_2_leiden_umap.manifest.json",
    "fallback_png_key": "png/ucsf/ebs2_2_leiden_umap.png",
    "default_state": {
      "color_by": "leiden"
    }
  }
}
```

### Interactive manifest contents

The JSON file referenced by `manifest_key` should look like:

```json
{
  "schema_version": "morphic-plotly-shiny-v1",
  "viewer": "plotly_shiny",
  "plot_id": "ucsf_ebs2_2_umap_leiden",
  "plot_kind": "umap_scatter",
  "title": "UCSF CRISPRi EBs2_2: UMAP colored by Leiden",
  "datasets": [
    {
      "dataset_id": "cells",
      "format": "parquet",
      "s3_key": "plotly/ucsf/ebs2_2_cells.parquet",
      "grain": "cell"
    }
  ],
  "encodings": {
    "x": {
      "field": "umap_1",
      "type": "quantitative"
    },
    "y": {
      "field": "umap_2",
      "type": "quantitative"
    },
    "color": {
      "field": "leiden",
      "type": "categorical"
    },
    "hover": [
      "cell_id",
      "best_feature",
      "feature_call_category",
      "n_genes",
      "mt_pct"
    ]
  },
  "controls": {
    "color_by": [
      "leiden",
      "sample_label",
      "feature_call_category"
    ]
  },
  "plotly_defaults": {
    "trace_type": "scattergl",
    "layout": {},
    "config": {}
  },
  "summary": {
    "n_cells": 10897,
    "n_clusters": 14
  },
  "analysis_context": {
    "lineage_label": "Embryoid body",
    "sample_label": "EBs2_2",
    "cluster_key": "leiden",
    "embedding_key": "X_umap"
  },
  "data_dictionary": {
    "cell_id": "cell barcode",
    "umap_1": "UMAP dimension 1",
    "umap_2": "UMAP dimension 2",
    "leiden": "Leiden cluster label"
  },
  "provenance": {
    "pipeline": "scanpy",
    "pipeline_version": "v1",
    "generated_at": "2026-04-17T00:00:00Z"
  }
}
```

### Why use a manifest instead of storing raw Plotly only

- we need a stable place for `analysis_context`
- we need schema versioning
- we need provenance
- we need Shiny-readable controls and defaults
- we need named datasets rather than opaque traces
- we need a clean way to evolve without breaking old cards
- not every consumer of the JSON will be a Plotly renderer

## Contract-forcing starter set

If the goal is to force real definition work now, the first dynamic set should
be intentionally small but structurally demanding.

The starter set is **five plot kinds**. Four kinds covered cell, aggregate,
gene, and target × target grains; adding a fifth, `sample_overview_umap`, gives
the dataset a dedicated "open this first" orientation figure that combines a
pooled integrated UMAP with a per-sample small-multiples band. The summary
layout is distinct enough from a single colored scatter that treating it as a
mode of `umap_scatter` would obscure it; keeping it separate keeps the
biological question ("does this dataset hang together across 16 samples?")
addressable as a top-level card.

### 1. `sample_overview_umap`

Purpose:

- serve as the dataset-orientation summary slide — the one figure a reviewer
  or collaborator opens first to understand the shape of the UCSF data before
  any perturbation-aware analysis
- show pooled integrated UMAP coordinates with Leiden structure and
  simultaneously expose per-sample coverage, dropout, and reproducibility
  through a small-multiples band

**Two-pass usage pattern.** The same `sample_overview_umap` manifest is
instantiated twice on a UCSF study page, driven by different data populations
rather than a different plot kind:

- **Tier 1B light pass** (perturbation-agnostic). Emitted right after per-sample
  technical QC, before any per-target efficacy analysis. The cell-level
  dataset omits `mixscape_class` population (column may exist but values are
  all `unassigned`); controls that reference Mixscape are hidden. This is the
  "is this dataset worth drilling into" gate.
- **Tier 2A full pass.** Emitted after Tier 1C (Mixscape + E-distance + DE)
  has run. Same manifest, now with `mixscape_class`, `feature_call_category`,
  and optional per-target MELD overlays populated. This is the analytical
  landing card on the study page.

Cards produced by the two passes should carry distinct `dataset_id`s
(e.g., `..._tenor` and `..._orientation`) so the frontend can distinguish
them even though they share the plot kind. A `card_role` of
`dataset_tenor_light` or `dataset_orientation_full` is recommended on each
card record.

Why this plot forces the contract (beyond what `umap_scatter` covers):

- requires a **composite layout** (main panel + small-multiples band + context
  strip) rather than a single panel
- requires **faceting semantics** on a categorical field (`sample_label`)
  with shared coordinates and gray-backdrop highlighting
- requires **sample-highlight** as a first-class control (pick a sample →
  that sample's cells foreground, all others backdrop)
- requires an auxiliary aggregate dataset for the per-sample QC strip /
  cluster × sample composition panel, alongside the cell-level dataset
- distinguishes "orientation view" vs "exploration view" card roles so the
  frontend can feature the orientation card on the gene page without it being
  confused with the drill-down UMAP

Minimum cell-level columns (same grain as `umap_scatter`):

- `cell_id`
- `sample_label`
- `lineage_label`
- `umap_1`
- `umap_2`
- `leiden`
- `mixscape_class`
- `feature_call_category`

Minimum per-sample aggregate columns (for the QC strip and cluster composition):

- `sample_label`
- `lineage_label`
- `n_cells`
- `n_mixscape_perturbed`
- `frac_mixscape_perturbed`
- `median_n_genes`
- `median_mt_pct`
- `dominant_leiden_cluster`

Layout modes to support:

- `composite` — main integrated UMAP + small-multiples band + QC strip (the
  summary-slide default)
- `pooled_only` — single integrated UMAP colored by a chosen categorical
  (Leiden, sample, or lineage)
- `small_multiples_only` — 4 × 4 grid of per-sample mini-UMAPs with shared
  coordinates, each on a gray backdrop of all other samples
- `lineage_split` — two pooled UMAPs side by side, EBs and iPSC, on separate
  integrated coordinate systems

### 2. `umap_scatter`

Purpose:

- explore state structure
- color by cluster, sample, guide-call category, Mixscape class, or selected
  target metadata
- render MELD density overlays per target (continuous scalar field on the
  manifold) as a primary display mode, not an afterthought

Why this plot forces the contract:

- requires cell-level grain
- requires categorical and numeric color semantics
- requires hover fields
- requires filter state
- requires a clear distinction between coordinates, overlays, and
  density/scalar fields
- requires a density/hex mode for cases where categorical color has too many
  levels to read

Minimum cell-level columns:

- `cell_id`
- `sample_label`
- `lineage_label`
- `umap_1`
- `umap_2`
- `leiden`
- `feature_call_category`
- `best_feature`
- `mixscape_class` (`KO` / `NP` / `NTC` / `unassigned`)
- `mixscape_score`
- `n_genes`
- `total_counts`
- `mt_pct`

Supplementary cell-level columns for overlays (one column per target on the
dynamic view):

- `meld_score__<target_symbol>` — MELD relative likelihood for cells given
  the target's perturbation vs NTC; present for the top-N dynamic targets

Display modes to support:

- point scatter colored by categorical or continuous field
- density / hex mode (for crowded categorical colorings)
- MELD scalar overlay mode per target
- gene-expression lookup mode (lazy-loaded per-gene vectors)

### 3. `target_qc_dotplot`

Purpose:

- show whether each targeted gene is actually reduced vs NTC
- summarize power, effect size, and evaluability across sample or lineage
- expose per-cell perturbation fraction (Mixscape) and E-distance alongside
  the mean-level knockdown so a single glance captures all three Tier 1B
  signals

Why this plot forces the contract:

- requires aggregate rather than cell-level grain
- requires simultaneous size and color encoding
- requires explicit threshold metadata
- requires "not evaluable" semantics rather than silently dropping rows

Minimum aggregate columns:

- `target_gene`
- `sample_label`
- `lineage_label`
- `n_target_cells`
- `n_ntc_cells`
- `baseline_expression`
- `log2fc_vs_ntc`
- `fdr`
- `mixscape_fraction_perturbed`
- `edistance`
- `edistance_pvalue`
- `guide_concordance_spearman`
- `is_evaluable`
- `non_evaluable_reason`
  (controlled vocabulary: `not_expressed_in_lineage`, `low_power`,
  `no_mixscape_fit`, `suspicious`, `ok`)

### 4. `de_volcano`

Purpose:

- inspect per-comparison differential expression results
- label top genes and support direct table download
- support both single-target volcano mode and a **volcano-matrix
  small-multiples mode** across a configurable set of targets as the DE summary
  figure

Why this plot forces the contract:

- requires gene-level grain
- requires threshold semantics
- requires multiple comparisons or condition metadata
- requires linked tabular output
- requires a small-multiples layout specification (rows × columns of facets
  keyed by comparison)

Minimum gene-level columns:

- `gene_id`
- `symbol`
- `log2fc`
- `padj`
- `pvalue`
- `base_mean`
- `is_target_gene`
- `condition_id`
- `de_method` (`wilcoxon`, `pseudobulk_deseq2`, `sceptre`)
- `mixscape_status` (`used`, `failed`, `na`)

### 5. `target_correlation_heatmap`

Purpose:

- render the target × target signature-correlation matrix (Tier 3C headline)
- cluster hierarchically, allow interactive row/column reordering, and expose
  the shared gene list on hover/click between any two targets

Why this plot forces the contract:

- introduces a target × target grain that is different from cell / aggregate /
  gene
- requires a primary matrix payload plus per-pair gene-list payloads that are
  resolved on demand
- requires hierarchical clustering / row-and-column order state to be part of
  the manifest defaults
- requires interactive reordering, row pinning, and module annotation
- exercises the hover / drill-down contract harder than volcano or dotplot do

Minimum primary-matrix columns (long format):

- `target_a`
- `target_b`
- `spearman`
- `pearson`
- `n_shared_genes`
- `module_id_a`
- `module_id_b`

Minimum per-pair payload columns (resolved on hover/click, typically
parquet-keyed by `target_a,target_b`):

- `gene_id`
- `symbol`
- `log2fc_a`
- `log2fc_b`
- `shared_rank`

These five plots together define the hard parts of the contract:

- composite-layout dataset orientation (main + small-multiples + QC strip,
  with sample faceting and highlight semantics)
- cell-level embedding data (with continuous-scalar overlay semantics)
- aggregate matrix-style QC data (with controlled evaluability vocabulary)
- gene-level DE data (with multi-method and small-multiples semantics)
- target × target correlation data (with on-demand gene-list payloads)

Once these five exist, cluster composition, pathway heatmaps, and guide-level
inspectors are extensions within existing grains rather than new contract
classes.

## Example manifest files

Concrete examples live here:

- `docs/examples/ucsf_plot_contracts/sample_overview_umap.manifest.example.json`
- `docs/examples/ucsf_plot_contracts/umap_scatter.manifest.example.json`
- `docs/examples/ucsf_plot_contracts/target_qc_dotplot.manifest.example.json`
- `docs/examples/ucsf_plot_contracts/de_volcano.manifest.example.json`
- `docs/examples/ucsf_plot_contracts/target_correlation_heatmap.manifest.example.json`
- `docs/examples/ucsf_plot_contracts/mockup_cards.html`

Those files are not just illustrations. They are the proposed v1 shape for the
first dynamic UCSF outputs.

## Which UCSF plots should be dynamic

The first dynamic candidates are the ones where interactivity is actually worth
the complexity — where a second question follows naturally from the first and
a static image cannot answer it without regenerating.

Recommended dynamic (v1):

- dataset-orientation summary composite (pooled integrated UMAP + 16-sample
  small-multiples band + cluster × sample composition strip), with
  sample-highlight and layout-mode controls
- UMAP colored by Leiden / sample / guide-call / Mixscape class
- UMAP with MELD density overlay per target
- UMAP with continuous gene-expression lookup (lazy-loaded per-gene vectors)
- target knockdown QC dotplot (with Mixscape and E-distance fields)
- DE volcano for selected comparisons (single-target)
- DE volcano-matrix small-multiples across the top-powered targets
- target × target signature-correlation heatmap (with hover-revealed shared
  gene lists)

Cross-view linking worth building in v1:

- linked UMAP ↔ DE table ↔ guide-strip — lasso on UMAP recomputes DE on the
  selection; clicking a gene in the DE table recolors the UMAP; clicking a
  target reveals the guide-level strip plot for concordance QC
- target dotplot drill-down — clicking a target cell opens the volcano + UMAP
  overlay + guide-concordance strip for that target in one coordinated panel
- correlation-heatmap drill-down — clicking a pair reveals the per-pair gene
  list and the scatter of log2FC_a vs log2FC_b

Keep static-only, at least initially:

- cell-count summaries
- guide-call barplots
- pathway / Hallmark / TF barplots and heatmaps (static rendering from
  decoupler-py outputs is fine for v1; promote to dynamic only if users need
  to reorder or drill into gene lists)
- guide-level knockdown strip plots (emit as a static PNG inside the Tier 1B
  card; promote to an interactive mode only when wired into the UMAP linking
  above)

Those are already understandable as static figures and do not gain enough from
Plotly to justify the payload and manifest overhead at v1.

## UCSF-specific packaging recommendation

For each produced figure, write a small bundle:

1. static image:
   - `png/...`
2. optional table:
   - `tsv/...`
3. optional interactive manifest:
   - `plotly/...manifest.json`
4. optional interactive data payloads:
   - `plotly/...parquet`
5. one metadata record that points to the files above

That is a much better contract than forcing the frontend to infer what a file
means from its filename.

## Naming and metadata recommendations

### Study label

Use a stable study label. I would start with:

- `UCSF_CRISPRi_PerturbSeq`

If lineage-specific grouping becomes important, keep that in
`analysis_context.lineage_label` rather than splitting the study label
immediately.

### Title format

Keep titles short and front-load the meaningful distinction.

Good examples:

- `UCSF CRISPRi EBs2_2: UMAP colored by Leiden`
- `UCSF CRISPRi iPSC2_1: target knockdown summary`
- `UCSF CRISPRi pooled EBs: PAX6 vs NTC volcano`

Avoid titles that require the frontend to parse long free text.

### Analysis context

For UCSF, `analysis_context` should usually carry:

- `plot_kind`
- `lineage_label`
- `sample_label` or `pool_label`
- `study_assay`
- `target_gene_symbols`
- `cluster_key` when relevant
- `embedding_key` when relevant
- `comparison_label` for DE / validation outputs
- `n_cells`
- `parse_status`
- `parse_confidence`

## Styling defaults

To fit the current MorPhiC page:

- use a white-background theme
- use `plotly_white` for dynamic figures
- export static PNGs at `1920 x 1280` when practical
- assume the image may be viewed at about `400px` width on the card
- avoid tiny legends and dense text blocks
- prefer clean subtitles or annotations over paragraph text inside the figure

## Recommended first implementation

For the first UCSF dynamic pass, I would do this:

1. Emit MorPhiC-compatible static metadata records for:
   - Tier 1 technical QC figures
   - Tier 1 knockdown dotplots
   - Tier 2 UMAP/Leiden figures
   - Tier 2 DE volcanoes
2. Keep `s3_png_key` and `s3_tsv_key` as the canonical live fields.
3. Add `dynamic_plot` for exactly three plot kinds:
   - `umap_scatter`
   - `target_qc_dotplot`
   - `de_volcano`
4. Back those views with explicit manifest + dataset files, not raw Plotly only.
5. Reuse existing MorPhiC-style DE fields for DE cards.

That gives immediate compatibility with the current site while preserving a
clean path to interactive embeddings later.

## Bottom line

For MorPhiC compatibility today:

- **static PNG + optional TSV + metadata is the standard**
- **Plotly/Shiny should be introduced through a versioned manifest contract**
- the right move is to make three dynamic plot types real now, not to wait for
  a general-purpose frontend abstraction to appear later
