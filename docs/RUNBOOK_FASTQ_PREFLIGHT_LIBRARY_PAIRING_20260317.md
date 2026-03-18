# Runbook: FASTQ Preflight for Chemistry Detection and GEX/Guide Pairing

**Date**: 2026-03-17  
**Status**: proposed iteration plan  
**Primary motivation**: avoid UCSF-style confusion where files are valid FASTQs but
are mislabeled, mixed-chemistry, or routed into the wrong assignment path.

## Goal

Add a lightweight preflight stage before expensive STAR/pf-multi runs that:

1. samples `R1/R2` directly from each FASTQ set
2. detects likely library type from `R2`
3. detects barcode chemistry from `R1` (`TRU` vs `NXT`)
4. normalizes barcodes to a canonical namespace (`TRU`)
5. pairs GEX and guide libraries by barcode overlap
6. emits a report and suggested config before the main run starts

This is intended to prevent:

- guide FASTQs being treated like GEX
- GEX FASTQs being treated like guide capture
- `NXT`/`TRU` namespace mismatches
- incorrect feature offset assumptions hiding the real issue

Related context:

- [HANDOFF_UCSF_FASTQ_MISLABEL_INVESTIGATION_20260316.md](/mnt/pikachu/STAR-suite/docs/HANDOFF_UCSF_FASTQ_MISLABEL_INVESTIGATION_20260316.md)
- [HANDOFF_UCSF_AALG1_TRU_AUTODETECT_20260316.md](/mnt/pikachu/STAR-suite/docs/HANDOFF_UCSF_AALG1_TRU_AUTODETECT_20260316.md)
- [RUNBOOK_STRICT_BARCODE_NAMESPACE_NORMALIZATION_20260303.md](/mnt/pikachu/STAR-suite/docs/RUNBOOK_STRICT_BARCODE_NAMESPACE_NORMALIZATION_20260303.md)
- [HANDOFF_BOOTSTRAP_TIERED_HASH_ACTIVATION_20260317.md](/mnt/pikachu/STAR-suite/docs/HANDOFF_BOOTSTRAP_TIERED_HASH_ACTIVATION_20260317.md)

## Inputs

For each candidate library or FASTQ set:

- `R1` FASTQ path(s)
- `R2` FASTQ path(s)
- optional expected sample count `k`
- whitelist resources:
  - TRU whitelist
  - NXT translation whitelist
- optional feature reference for guide libraries

Recommended sampling budget:

- `200000` reads from `R1` (configurable via `--sample-reads`)
- `200000` reads from `R2`

This is large enough to stabilize top-barcode overlap while remaining cheap
(~1 s/file for gzipped FASTQs). If pairing S/N is thin, increase to `500000`.
See [Tuning](#tuning-parameters) below.

## Canonical outputs

The preflight should write a machine-readable report per run, for example:

- `preflight_summary.tsv`
- `preflight_pairing_matrix.tsv`
- `preflight_pairing_report.json`

Per FASTQ set, record at minimum:

- inferred library type: `GEX`, `GUIDE`, `AMBIGUOUS`, `UNKNOWN`
- inferred chemistry: `TRU`, `NXT`, `AMBIGUOUS`, `UNKNOWN`
- chemistry confidence metrics
- top corrected barcodes in canonical `TRU`
- detected guide-anchor / scaffold evidence in `R2`
- observed guide offset histogram if guide-like
- suggested paired partner
- pairing confidence

## Proposed preflight algorithm

### 1. Sample reads

For each FASTQ set:

- stream the first `100000` `R1` reads
- stream the first `100000` `R2` reads
- do not require full decompression of the file

If multiple lanes belong to one logical library, aggregate the sample across
lanes before computing chemistry or pairing.

### 2. Detect chemistry from `R1`

For sampled `R1` barcodes:

1. extract the candidate barcode bases using the known read structure
2. score against:
   - TRU whitelist
   - NXT whitelist / translation surface
3. compute:
   - exact whitelist hit rate
   - corrected hit rate
   - translated hit rate to canonical `TRU`

Canonical rule:

- internal comparison namespace = `TRU`
- if file looks `NXT`, translate sampled barcodes to `TRU` immediately

Suggested chemistry call:

- `TRU` if TRU hit rate strongly dominates NXT
- `NXT` if NXT hit rate strongly dominates TRU
- `AMBIGUOUS` if both are plausible
- `UNKNOWN` if neither is plausible

This stage should reuse existing namespace normalization logic where possible,
not invent a new translation path.

### 3. Detect library type from `R2`

For sampled `R2` reads, score evidence for guide capture vs non-guide content.

Guide indicators can include:

- guide scaffold / anchor motifs from the feature reference
- common guide-capture adapter context
- recurrent fixed preambles such as the UCSF TSO case
- strong exact / `H1` support for guide sequences at a consistent offset band

Record:

- fraction of reads with guide-like anchor evidence
- dominant guide offset histogram
- top guide-like sequences or anchors

Suggested library-type call:

- `GUIDE` if guide-anchor/scaffold evidence is strong and coherent
- `GEX` if guide evidence is absent and barcode chemistry looks like ordinary
  GEX
- `AMBIGUOUS` if mixed
- `UNKNOWN` if too weak to classify

Important:

- preflight should capture the *observed* guide offset band
- it should not assume offset `0`

### 4. Build canonical top-barcode signatures

For each file:

1. correct sampled `R1` barcodes using the inferred chemistry
2. normalize them to canonical `TRU`
3. count barcode frequencies
4. keep the top `N` (configurable via `--top-barcodes`, default `500`)

Validated signature widths (UCSF 16-sample Perturb-seq, 10x Chromium V3):

| `top_n` | Pairing result | Notes |
|---------|---------------|-------|
| 50 | classification ✓, pairing ✗ | Too narrow; GEX/GUIDE rank distributions differ enough that top-50 sets barely overlap |
| 200 | adequate for clean data | |
| **500** | **robust (default)** | **Worst-case pairing S/N = 50x on UCSF dataset** |
| 1000 | diminishing returns | Only needed if S/N is still thin with 500 |

The signature also records diagnostic slices at top-50, top-100, and top-500
for the output report.

### 4b. Merge duplicate / multi-lane files

Before cross-type pairing, the preflight automatically detects and merges files
that represent the same logical library (e.g., duplicate copies from different
delivery directories, or multi-lane FASTQs from the same library).

Merge criteria (all must hold):

- same inferred chemistry (TRU/NXT)
- same inferred library type (GEX/GUIDE)
- pairwise set-Jaccard on barcode signatures ≥ `MERGE_JACCARD_MIN` (0.20)

Merged groups pool their barcode counters (union of counts, trimmed to `top_n`)
to create a single stronger signature for pairing. This is critical when the
same library appears in multiple directories — without merging, these duplicates
absorb the "second-best" slot and destroy pairing S/N.

Validated on UCSF dataset: 64 individual files → 32 merged logical libraries
(each original directory pair merged). Worst-case merge S/N = 288x.

### 5. Pair files by barcode overlap

Construct pairwise similarity between all files after chemistry normalization.

Primary metric:

- top-barcode Jaccard on canonical `TRU`

Secondary metrics:

- weighted Jaccard using sampled counts
- overlap of top-10 and top-50 ranks
- chemistry compatibility
- library-type compatibility (`GEX` vs `GUIDE`)

### 6. Infer groups and assign files to samples

Do not require `k` up front.

The default behavior should infer sample groups from the barcode-overlap graph.

#### Iteration 1: overlap graph first

After chemistry normalization and library-type detection:

1. build a pairwise similarity matrix across files
2. keep only strong, chemistry-compatible edges
3. require reciprocal-best-match or a clear best-vs-second-best margin
4. form connected components on the retained graph

Each connected component is the inferred file group for one biological sample
or sample-family.

This should be the default because it naturally handles:

- unknown sample count
- missing mates
- extra lanes
- orphan libraries
- mislabeled files
- partial runs where not every sample has every library type

If the graph is clean, use it directly and do not invoke clustering.

#### Iteration 2: deterministic matching within components

Inside each clean graph component:

1. split files into likely `GEX` and likely `GUIDE`
2. compute all cross-group Jaccards
3. greedily or bipartite-match the highest-confidence pairs
4. reject weak or conflicting matches

This remains the simplest first version and should cover most cases.

#### Iteration 3: clustering only for ambiguous cases

If the graph does not separate cleanly, escalate to one of:

- manual review
- optional clustering summary

If clustering is needed:

- if `k` is known, use `k = number_of_samples`
- if `k` is unknown, estimate it from the graph or from a dendrogram cut, then
  cluster
- cluster on barcode-set similarity after canonical normalization

`k`-medoids is still a reasonable second-phase option because:

- the distance can be `1 - Jaccard`
- it is robust to outliers
- it works directly on pairwise distances

But if `k` is unknown, hierarchical clustering is usually a better first
fallback than fixed-`k` medoids because it does not force an up-front sample
count.

## Confidence and failure rules

The preflight should be willing to say “unclear”.

Hard-stop cases:

- chemistry `UNKNOWN`
- guide/GEX type `UNKNOWN`
- two candidate pairings with near-identical Jaccard
- graph structure that does not admit a confident component split
- no candidate pair above a minimum overlap threshold

Soft-warning cases:

- chemistry `AMBIGUOUS`
- low barcode diversity in the sample
- guide offset histogram is broad instead of concentrated
- one file looks like a mixed library

The report should explicitly separate:

- `safe auto-pair`
- `auto-pair with warning`
- `manual review required`

## Suggested implementation phases

### Phase 1: offline script

Build a standalone script first.

Inputs:

- FASTQ paths
- whitelist paths
- optional feature reference
- optional expected sample count `k`

Outputs:

- summary TSV/JSON
- pairing matrix
- suggested library grouping

This lets us validate on UCSF/MSK style problem cases without touching the main
STAR run path.

### Phase 2: STAR preflight mode

Add a dry-run or preflight-only mode in STAR/pf-multi that:

- runs the same sampling logic
- prints the report
- exits before the expensive assignment stage

Possible interface:

- `--pfPreflight yes`
- `--pfPreflightSampleReads 100000`
- `--pfPreflightExpectedSamples K`
- `--pfPreflightOut <dir>`

### Phase 3: optional automatic pairing/config synthesis

Once validated, allow preflight to emit a suggested config file:

- corrected library labels
- inferred chemistry per library
- inferred sample groups
- proposed GEX/guide pairing

The main run should still require either:

- explicit user acceptance
- or a confidence threshold high enough to trust automation

## Validation plan

Use real confusion cases first.

Priority fixtures:

1. UCSF mixed `AALG1/AALG2` TRU/NXT sets
2. small UCSF guide fixture where guide offset is around `31`
3. MSK datasets where chemistry or feature routing was previously confusing

Acceptance criteria for Phase 1:

- chemistry classification is correct on known TRU/NXT fixtures
- guide vs GEX classification is correct on known mixed fixtures
- canonical barcode normalization reproduces known same-sample pairings
- top-barcode Jaccard separates true pairs from false pairs with a clear gap
- inferred overlap-graph components match known sample structure when `k` is
  not provided

Acceptance criteria for Phase 2:

- no expensive STAR alignment is launched before preflight report is available
- preflight adds only a small fixed startup cost
- failures are explicit and actionable

## Tuning parameters

The two key knobs are `--sample-reads` and `--top-barcodes`. Both are printed
in the `PREFLIGHT SUMMARY` header and recorded in the JSON report under
`parameters`.

| Parameter | Default | When to increase | Cost |
|-----------|---------|-----------------|------|
| `--sample-reads` | 200,000 | Pairing S/N < 10x; very low barcode diversity | ~1 s/file per 200K reads (gzip) |
| `--top-barcodes` | 500 | Pairing S/N < 10x after increasing reads | Negligible (in-memory counter trim) |

**Decision rule**: if the `PAIRING CONFIDENCE` table shows S/N below 10x on any
pair, double `--sample-reads` first (cheap), then widen `--top-barcodes` to
1000. If S/N is still below 5x, the libraries likely lack sufficient shared
barcode population and manual review is warranted.

The `MERGE CONFIDENCE` table typically shows S/N > 100x because same-library
duplicates share nearly identical barcode profiles. If merge S/N is low (< 20x),
the files being merged may not actually be from the same library — inspect the
`NAME_MISMATCH` warnings.

## Validated results

### UCSF 16-sample Perturb-seq (10x Chromium V3, AALG1=GEX/TRU, AALG2=GUIDE/NXT)

- **Input**: 64 individual FASTQ sets (32 from `GEX/` dir, 32 from `guides/` dir)
- **Parameters**: `--sample-reads 200000 --top-barcodes 500`
- **Runtime**: ~65 s total (~1 s/file for gzip decompression + scoring)

| Metric | Value |
|--------|-------|
| Chemistry detection | 32/32 TRU (AALG1), 32/32 NXT (AALG2) — 100% correct |
| Library type detection | 32/32 GEX, 32/32 GUIDE — 100% correct |
| Merge phase | 64 → 32 logical libraries (all duplicates detected) |
| Merge S/N (worst) | 288x |
| Pairing | 16/16 correct pairs (matches corrected symlink set) |
| Pairing S/N (worst) | 50x |
| Pairing S/N (median) | 107x |
| Name warnings | 16× PAIR_NAME_MISMATCH (AALG1≠AALG2, expected for this dataset) |

## Implementation status

Phase 1 is implemented as `scripts/preflight_library_pairing.py`.

Pipeline (all barcode-driven, no filename assumptions for pairing):

1. sample `200K` reads from each `R1/R2` (`--sample-reads`)
2. detect chemistry from `R1` (TRU/NXT/AMBIGUOUS/UNKNOWN)
3. normalize sampled barcodes to canonical `TRU`
4. detect guide evidence and offset band from `R2`
5. build top-N barcode signatures (`--top-barcodes`, default 500)
6. merge same-type/chemistry duplicates (Jaccard ≥ 0.20) into logical libraries
7. pairwise cross-type Jaccard scoring on merged libraries
8. connected-component detection + greedy bipartite GEX↔GUIDE matching
9. confidence tables (S/N for both merge and pairing decisions)
10. post-hoc name-based sanity checks (anomalous filenames, missing R1/R2,
    empty GEX/GUIDE sides)

The overlap graph is the primary algorithm; fixed-`k` clustering is not used.
Hierarchical clustering / `k`-medoids remain available as a refinement for
noisy multi-file runs if the graph does not separate cleanly.

This gives a cheap (~1 min for 64 files), explainable gate before the main
pipeline and would have prevented the specific confusion that triggered this
runbook.
