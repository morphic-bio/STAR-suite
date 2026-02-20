# Handoff: CR-Compatible Multimap Policy (UCSF/A375)

Date: 2026-02-19  
Scope owner: next agent  
Goal: mirror Cell Ranger multimap resolution behavior in a STAR compatibility mode.

## Objective

Improve STAR-vs-Cell Ranger GEX parity by harmonizing multimapper resolution in a
compatibility path, while keeping current STAR defaults unchanged.

## What Was Completed

### 1) Script infrastructure updated for multimap-policy comparisons

Updated:
- `tests/compare_feature_mex.py`
- `scripts/report_additional_parity_metrics.py`
- `scripts/run_gex_feature_parity_checks.sh`
- `scripts/ucsf_parity/README.md`
- `docs/UCSF_PARITY_RECOVERY_RUNBOOK_20260218.md`

New capabilities:
- Side-specific matrix selection (CR vs STAR):
  - `--matrix-basename-a`
  - `--matrix-basename-b`
- Run-level matrix controls:
  - `--cr-raw-matrix-basename`
  - `--cr-filtered-matrix-basename`
  - `--star-raw-matrix-basename`
  - `--star-filtered-matrix-basename`
- Real-valued MatrixMarket support for STAR `UniqueAndMult-*` matrices.

This enables direct policy probes (for example STAR `Unique` vs
`UniqueAndMult-Rescue`) without changing analysis code each run.

### 2) Methodology documentation aligned with current workflow

Runbook now includes:
- explicit multimap compatibility workstream
- commands for Unique baseline and Rescue probe
- interpretation guardrails for raw vs filtered matrix availability.

## Findings To Carry Forward

1. Namespace harmonization is mandatory before parity analysis.
- Use 2-column whitelist mapping NXT->TRU.
- Avoid CR-side remapping when CR is already TRU.

2. Feature-call nondeterminism from PF matching was fixed separately.
- Remaining STAR-vs-CR differences are not due to that nondeterminism.

3. Multimapper handling explains a meaningful part of GEX discordance.
- Differences are not only total-count under/over counting; there is
  classification/policy divergence on ambiguous loci.

4. `UniqueAndMult-*` inspection is required for policy debugging.
- `matrix.mtx` alone cannot answer all multimap-policy questions.

## Implementation Target (Compatibility Mode)

Add a compatibility policy that prefers gene-annotated candidates over
non-gene candidates among equal-best alignments in Solo counting logic.

Proposed behavior (mode ON only):
1. For a read with multiple best-score alignments:
   - if at least one candidate is gene-annotated and others are non-gene,
     keep gene-annotated candidate set, drop non-gene contenders.
2. Continue existing multimapper handling (`Unique`, `Rescue`, etc.) on the
   post-filtered candidate set.
3. Preserve default behavior when compatibility flag is OFF.

## Candidate Code Touchpoints

- `core/legacy/source/SoloReadFeature_record_base.cpp`
  - read-level unique/multigene decision path (`readGe->size()` handling).
- `core/legacy/source/SoloFeature_collapseUMIall.cpp`
  - multimapper UMI redistribution (`Uniform/Rescue/PropUnique/EM`).
- `core/legacy/source/SoloFeature_outputResults.cpp`
  - `UniqueAndMult-*` output generation and auditability.
- `core/legacy/source/Parameters.h`
- `core/legacy/source/Parameters.cpp`
- `core/legacy/source/parametersDefault`
  - add compatibility flag(s), default OFF.

## Validation Checklist

1. UCSF parity regression:
- run Unique baseline and Rescue probe using the updated scripts.
- compare:
  - cell-level GEX correlations
  - gene-level GEX correlations
  - feature-call parity (guide/target/gene collapse).

2. A375 control:
- rerun with same analysis scripts and confirm no regression in expected high
  parity behavior.

3. Determinism:
- two repeated STAR runs with identical settings must match exactly on MEX
  files and CRISPR call CSVs.

4. Backward compatibility:
- with compatibility mode OFF, outputs must remain unchanged vs current baseline.

## Immediate Next Step

Implement the compatibility gate in Solo read classification (gene-over-non-gene
tie policy), then rerun UCSF Unique/Rescue parity with the updated script
switches to quantify impact before broader refactors.
