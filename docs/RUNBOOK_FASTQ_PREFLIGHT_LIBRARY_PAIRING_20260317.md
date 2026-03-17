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
- expected sample count `k`
- whitelist resources:
  - TRU whitelist
  - NXT translation whitelist
- optional feature reference for guide libraries

Recommended sampling budget:

- `100000` reads from `R1`
- `100000` reads from `R2`

This is large enough to stabilize top-barcode overlap while remaining cheap.

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
4. keep the top `N`

Recommended initial `N`:

- start with `50`
- also record `100` and `500` for diagnostics

The main signature for pairing should be the top-50 canonical barcode set plus
their counts.

### 5. Pair files by barcode overlap

Construct pairwise similarity between all files after chemistry normalization.

Primary metric:

- top-barcode Jaccard on canonical `TRU`

Secondary metrics:

- weighted Jaccard using sampled counts
- overlap of top-10 and top-50 ranks
- chemistry compatibility
- library-type compatibility (`GEX` vs `GUIDE`)

### 6. Assign files to samples

Two viable approaches:

#### Iteration 1: deterministic matching

Use simple, auditable rules first:

1. split files into likely `GEX` and likely `GUIDE`
2. compute all cross-group Jaccards
3. greedily or bipartite-match the highest-confidence pairs
4. reject weak or conflicting matches

This is the simplest first version and should cover most cases.

#### Iteration 2: clustering

If the run has more files than simple pairs, or if labels are noisy:

- use `k = number_of_samples`
- cluster on barcode-set similarity after canonical normalization
- prefer a medoid-style clustering because the sample count is small and the
  metric is not Euclidean

`k`-medoids is a reasonable second-phase improvement because:

- the distance can be `1 - Jaccard`
- it is robust to outliers
- it works directly on pairwise distances

But it should not block the first implementation. The deterministic matcher is
easier to debug and explain.

## Confidence and failure rules

The preflight should be willing to say “unclear”.

Hard-stop cases:

- chemistry `UNKNOWN`
- guide/GEX type `UNKNOWN`
- two candidate pairings with near-identical Jaccard
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
- expected sample count `k`

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

Acceptance criteria for Phase 2:

- no expensive STAR alignment is launched before preflight report is available
- preflight adds only a small fixed startup cost
- failures are explicit and actionable

## Initial recommendation

Implement the deterministic preflight first.

Specifically:

1. sample `100K` reads from each `R1/R2`
2. detect chemistry from `R1`
3. normalize sampled barcodes to `TRU`
4. detect guide evidence and offset band from `R2`
5. compute top-50 barcode Jaccard between files
6. emit a pairing report

Do **not** start with clustering as the only algorithm. Add `k`-medoids after
the deterministic surface is validated, as a refinement for noisy multi-file
runs.

That gives us a cheap, explainable gate before the main pipeline and would have
prevented the specific confusion that triggered this runbook.
