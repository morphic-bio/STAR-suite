# Deterministic quality ranking for bootstrap OrdMag

Bootstrap OrdMag estimates how many barcodes to retain. The final selection
keeps the rounded estimate, clamped to the number of eligible barcodes above
the configured inclusive UMI floor (500 by default).
It no longer expands UMI ties and then reverts when expansion exceeds 20%.

## Ranking protocol

Within the OrdMag input window, sort by:

1. Total UMIs, descending.
2. Non-mitochondrial UMIs, descending, when an annotation mask is supplied.
3. Detected genes, descending; exclude mitochondrial genes when a mask is supplied.
4. Full barcode identity, lexicographically ascending.

At equal total UMIs, step 2 is equivalent to preferring a lower mitochondrial
fraction, without floating-point division. Detected genes are distinct feature
IDs with positive counts; zero entries and duplicate coordinates do not inflate
the score. The final barcode key includes the entire Flex `CB16+TAG8` identity.
The original index is only a fallback for missing or duplicate identities.

The selected indices, rather than `UMI >= retainThreshold`, define OrdMag calls.
`retainThreshold` records the last selected UMI and can split a UMI tie. Simple
cells remain the prefix of the EmptyDrops candidate list. The shared OrdMag
implementation enforces the inclusive UMI floor before constructing primary
calls, their median, and tail candidates. This applies to both bootstrap and
fallback paths, preventing low-count unused tags from inflating occupancy.

The bootstrap still estimates cell number. There is no bootstrap for tie-breaking.
Given fixed scores, moving the rounded target by one adds or removes one cell.
This does not bound how much the bootstrap target itself can move after a count
perturbation. Bootstrap streams also still depend on thread count.

## Shared implementation and assay annotations

`core/features/libscrna/include/OrdMagRank.h` provides sparse scoring and ranking.
The matrix-aware `runCRSimpleFilterBootstrap` overload accepts detected genes,
full barcode IDs, and optional non-MT UMI scores. Existing count-only callers keep
their three-argument signature, but need metadata for quality and identity ties.

The C API `scrna_emptydrops_run` computes detected-gene scores when sparse data is
available. `scrna_emptydrops_run_with_rank_mask` additionally accepts a byte mask
of length `input.n_features`, with 1 for mitochondrial features. Existing C ABI
struct layouts are preserved. The additive `scrna_emptydrops_run_with_rank_options`
also accepts an explicit bootstrap worker count. STAR Solo's bootstrap path uses
this shared C API. STAR loads the optional annotation with
`--soloCellFilterMitochondrialGenes PATH`; `--soloCellFilterBootstrapThreads N`
sets the bootstrap workers (0 inherits `runThreadN`). These options also apply to
the native tag-aware Flex caller.

The rule can be reused across RNA modules. Flex maps its targeted-gene feature
IDs to the mitochondrial annotation; other RNA assays use their own reference
feature annotations. No species-specific gene list is hardcoded in libscrna.
A matrix with probe IDs as rows must first map those IDs to the intended gene
features; counting probes is not the same as counting distinct genes.

## Standalone caller

Add these arguments to an existing full, bootstrap caller command:

```text
--features /path/to/raw/features.tsv
--mitochondrial-genes /path/to/mitochondrial_gene_ids.tsv
```

Use `--mode full --use-bootstrap`, or `--flex-tag-aware`, which enables both.
Both files must be plain text. The first whitespace-delimited column contains
exact feature IDs; gene-list blank lines and comments beginning with `#` are
ignored. Features must match matrix row order and row count. Empty lists and
lists with no matching IDs are rejected; partial matches are reported. The caller
does not guess gene IDs from display names or strip gene-version suffixes.

The MT mask affects tie scores only. Original matrix counts remain in the
bootstrap, ambient profile, EmptyDrops likelihoods, and outputs. Ambient and
retain-window UMI ties use barcode identity only, so quality scores do not select
low-quality cells preferentially into the background profile.

The September 8 Flex validation derives its list from `gene` records on the
mitochondrial chromosome in the actual GRCh38-2024-A reference GTF. That reference
has 13 such genes in the CR feature set and 12 in the STAR Flex feature set;
STAR's target set lacks MT-ATP8. The list, GTF hash, and matched counts are saved
with the validation artifacts. These are dataset-specific facts, not universal
mitochondrial gene counts.

## Verification and boundaries

Unit tests cover realistic UMI ties, quality priority, full barcode identity,
input permutation, explicit zeros, duplicate gene entries, the whole-tie midpoint
counterexample, and isolation of ambient selection from MT scores. Standalone
fixtures check both the C API and custom Flex ambient paths against independently
ranked expected barcode sets.

The non-bootstrap legacy knee method retains its threshold-based protocol.
Quality ranking does not introduce an MT filter or replace bootstrap cell-number
estimation. The separate primary-floor fix excludes one-UMI barcodes under the
default 500-UMI floor. Cell-call accuracy is measured separately; the complete
320K validation is documented in
`docs/HANDOFF_FULL320K_STAR_DEPRECATED_20260909.md`.

Validation artifacts and comparison results:
`/mnt/pikachu/star_suite_paper/analysis/ordmag_quality_ties_20260908/`.

## Native tag-aware Flex integration

`--soloFlexCellCaller tag-aware` selects the shared `runSimpleEDWithAmbient`
implementation used by `scrna_simpleed --flex-tag-aware`. This is the new STAR
default when the Flex filter runs; `--soloFlexCellCaller legacy` retains the
previous per-tag pipeline and its expected-cell settings. Tag-aware mode rejects
legacy count, ambient, occupancy, and fallback overrides rather than ignoring
them. Remove `--soloFlexExpectedCellsPerTag 3000` from the new-mode command.

Repeated biological sample labels in `--soloSampleWhitelist` and
`--soloFlexAllowedTags` group their tag sequences into one call and one output
MEX. Full `CB16+TAG8` identities are retained. Per tag, the protocol uses an
OrdMag window of 90,000, a maximum expected-cell search bound of 22,500, an
ambient rank start of 45,000 and end of 90,000. These four values scale by the
number of tags in the sample. The ambient UMI target is zero: the fixed rank
window does not extend to reach a target mass. The EmptyDrops candidate cap
remains 100,000 per sample. OrdMag is primary, with a BH-adjusted EmptyDrops
tail at default FDR 0.01 and 100,000 simulations. Observed-tag occupancy is
applied across the complete library after all grouped samples have been called.
Modeling-only features remain available to both calling stages and are omitted
from filtered exports using `--soloFlexFilteredGeneList`; see
`docs/FLEX_MODEL_FEATURES.md`.

Enable reproducible caller diagnostics with, for example:

```text
--soloRunFlexFilter yes
--soloFlexCellCaller tag-aware
--soloCellFilterMitochondrialGenes /path/to/mitochondrial_gene_ids.tsv
--soloCellFilterBootstrapThreads 48
--soloFlexDebugTagLog yes
--soloFlexDebugOutputDir /path/to/fresh/diagnostics
--soloFlexInvariantChecks yes
--soloFlexFatalOnError yes
```

Each sample diagnostic directory contains `caller_diagnostics.json`,
`ordmag_rank.tsv`, `ambient_profile.tsv`, `sample_tags.tsv`, and the standard
EmptyDrops results ledger. The JSON records the bootstrap mean, standard
deviation and workers, selected count, candidate floor, ambient mass, FDR,
and final calls. `feature_rank_mask.tsv` records the annotation mapping.
The ledger's `is_simple_cell` now reflects actual OrdMag membership, and
the tail count excludes primary cells; these correct diagnostic reporting
without changing the external caller's numerical decisions.

Integration fixtures compare internal and external per-cell ledgers for a
single tag and a two-tag sample, including shared CB16 values across tags,
MT ties, the 500-UMI boundary, and duplicate-tag rejection. Synthetic ambient
profiles need at least five distinct nonzero gene-count frequencies for the
existing Simple Good-Turing implementation; smaller degenerate profiles can
yield an empty simulation profile. This pre-existing limitation is not fixed
by the integration.

L004 integration artifacts:
`/mnt/pikachu/star_suite_paper/analysis/flex_internal_l004_20260908/`.
