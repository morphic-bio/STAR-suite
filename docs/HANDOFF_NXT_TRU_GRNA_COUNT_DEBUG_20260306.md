# Handoff: gRNA Parity in MSK 30-KO (Canonical-TRU)

**Date**: 2026-03-06
**Branch**: `CR-Larry-perturb`
**Commit**: `f16f30995d36d077f8c999a92a959fa4221d8b4b`

## Summary

Canonical-TRU barcode normalization resolves the MSK 30-KO mixed-chemistry
namespace mismatch. The MSK full-sample STAR run achieves **1.044x
CellRanger's guide UMIs** in the merged filtered matrix, with per-guide
Pearson **0.9999** and per-cell total guide UMI Pearson **0.9997**.

With the corrected 30-guide CellRanger CRISPR run
`/storage/MSK-perturb-comparison/cr_full_grna_30crispr_20260306_173247/`,
call-level concordance is also strong in the full 30-guide space:

- exact called guide-set match: **0.9944**
- exact singleton match: **0.99995**
- exact multi-guide match: **0.9506**
- mean call-set Jaccard: **0.9973**
- CR-in-STAR containment: **0.9963**
- STAR-in-CR containment: **0.9918**

The previously reported "63% undercount" was an artifact of comparing
the wrong STAR output stage (gRNA-only EmptyDrops, 5,643 cells) against
CellRanger's merged filtered matrix (32,256 cells). The three output
stages and their correct comparators are documented below.

## How To Compare This Correctly

Use the following rules:

1. **Guide count parity**: compare merged matrices, not gRNA-only
   EmptyDrops output.
   - STAR comparator: `outs/filtered_feature_bc_matrix/`
   - CellRanger comparator: `outs/per_sample_outs/.../count/sample_filtered_feature_bc_matrix/`
   - Match features by **feature ID**, not by `feature_type` string.

2. **Raw assignment sensitivity**: compare STAR `cr_assign/.../PolyIII/`
   against CellRanger `raw_feature_bc_matrix/` restricted to the same
   30 guide IDs.

3. **CRISPR calling parity**: compare STAR `outs/crispr_analysis/`
   against CellRanger `outs/.../count/crispr_analysis/`, but only after
   putting both runs in the same CRISPR feature space.
   - **Primary comparator**:
     `/storage/MSK-perturb-comparison/cr_full_grna_30crispr_20260306_173247/`
     where all 30 guides are typed `CRISPR Guide Capture`.
   - **Historical comparator**:
     `/storage/MSK-perturb-comparison/cr_full_grna_withcalls_20260306_095454/`
     where `SIX2` is typed `Custom`, so call parity must be projected to
     the shared 29-guide space by dropping `SIX2`.

4. **Do not compare** STAR `cr_assign/.../PolyIII/filtered/` to any
   CellRanger merged matrix. That directory is a feature-library-only
   EmptyDrops artifact, not the final parity output.

## Methodology

### Set Metrics

Use symmetric set metrics for overlap:

- **Filtered-barcode Jaccard**:
  `|STAR_filtered ∩ CR_filtered| / |STAR_filtered ∪ CR_filtered|`
- **Per-cell guide-set Jaccard**:
  for cells where both methods call at least one guide,
  `|STAR_call_set ∩ CR_call_set| / |STAR_call_set ∪ CR_call_set|`

Jaccard is appropriate for overlap, but not for count parity.

### Directional Set Metrics

Because CellRanger often behaves like a subset of STAR, also report
directional containment:

- **CR in STAR**:
  `|STAR_called ∩ CR_called| / |CR_called|`
- **STAR in CR**:
  `|STAR_called ∩ CR_called| / |STAR_called|`

Interpretation:

- `CR in STAR` near 1.0 means STAR recovers nearly all CellRanger calls.
- `STAR in CR` lower than `CR in STAR` means STAR is more permissive.
- This asymmetry is not necessarily bad; it is a concise way to show skew.

For call sets on cells where both methods call at least one guide, also
report mean cellwise containment:

- `mean_cellwise_CR_in_STAR = mean(|STAR_call_set ∩ CR_call_set| / |CR_call_set|)`
- `mean_cellwise_STAR_in_CR = mean(|STAR_call_set ∩ CR_call_set| / |STAR_call_set|)`

### Count Metrics

For abundance parity, use:

- merged raw / filtered total UMI ratios
- per-guide Pearson on merged totals
- per-cell total guide UMI Pearson on shared called cells

Do not use Jaccard as a replacement for these count-based statistics.

## Stage-by-Stage Guide UMI Accounting

### Stage 1: Raw gRNA Assignment Output

Raw `assignBarcodes` output before any cell calling.

| Metric | STAR | CellRanger | STAR/CR |
|--------|------|------------|---------|
| Total deduped guide UMIs | 4,751,051 | 4,515,296 | 1.052 |
| Barcodes with >=1 guide UMI | 238,726 | 231,128 | 1.033 |

STAR is slightly higher at raw assignment, likely due to Hamming-distance
correction differences. No undercount at this stage.

### Stage 2: gRNA-Only EmptyDrops Filtered Output

Feature-library standalone EmptyDrops cell calling on the gRNA library alone.
This is **not** the final output and should **not** be compared to CR's
merged matrix.

| Metric | STAR gRNA EmptyDrops |
|--------|---------------------|
| Total deduped guide UMIs | 1,958,729 |
| Called barcodes | 5,643 |
| Excluded barcodes | 233,083 |

The low barcode count is expected: EmptyDrops on a feature-only library
calls only the small number of barcodes with enough guide UMIs to pass
the ambient filter. This is an intermediate artifact, not the parity
comparator.

### Stage 3: Merged GEX-Filtered Guide Output (Correct Comparator)

Guide UMIs from the raw gRNA matrix intersected with the GEX-called
filtered cell set, then merged into `filtered_feature_bc_matrix/`.
Guide loss between raw gRNA output and final merged cells is **zero**.

| Metric | STAR merged filtered | CR merged filtered | STAR/CR |
|--------|---------------------|-------------------|---------|
| Total guide UMIs (30 guides) | 3,247,179 | 3,111,361 | 1.044 |
| Guide-positive cells | 30,263 | 32,017 | 0.945 |
| Filtered cells (all) | 30,497 | 32,256 | 0.945 |

Per-guide breakdown (merged filtered, sorted by name):

| Guide | CR UMIs | STAR UMIs | STAR/CR |
|-------|---------|-----------|---------|
| CARM1 | 38,540 | 40,231 | 1.044 |
| CBX4 | 199,910 | 209,246 | 1.047 |
| DPY30 | 72,489 | 74,967 | 1.034 |
| EOMES | 251 | 271 | 1.080 |
| FEV | 183,155 | 190,992 | 1.043 |
| GCM1 | 219,624 | 229,884 | 1.047 |
| HES1 | 6,666 | 6,950 | 1.043 |
| HNF1A | 147,197 | 152,705 | 1.037 |
| HNF4A | 140,795 | 147,545 | 1.048 |
| INSM1 | 81,069 | 84,132 | 1.038 |
| ISL1 | 179,755 | 187,267 | 1.042 |
| JMJD6 | 10,706 | 11,263 | 1.052 |
| MAFB | 107,181 | 111,695 | 1.042 |
| MEIS1 | 66,973 | 69,660 | 1.040 |
| MEN1 | 116,078 | 121,125 | 1.043 |
| NCOA5 | 110,363 | 115,507 | 1.047 |
| NEUROG3 | 133,577 | 138,035 | 1.033 |
| NKX6-1 | 430,248 | 455,247 | 1.058 |
| NOTCH1 | 16,337 | 17,053 | 1.044 |
| NOTCH2 | 86,376 | 90,477 | 1.047 |
| ONECUT1 | 86,111 | 88,953 | 1.033 |
| PAX6 | 151,081 | 157,267 | 1.041 |
| PBX2 | 89,711 | 93,260 | 1.040 |
| PDX1 | 200,938 | 207,148 | 1.031 |
| RFX3 | 27,247 | 28,825 | 1.058 |
| RFX7 | 6,143 | 6,492 | 1.057 |
| RREB1 | 6,067 | 6,293 | 1.037 |
| SIM1 | 98,822 | 102,993 | 1.042 |
| SIX2 | 74,503 | 76,593 | 1.028 |
| SIX3 | 23,673 | 25,103 | 1.060 |

Per-guide Pearson (merged filtered totals): **0.9999**. Ratio range:
1.028 (SIX2) to 1.080 (EOMES).

## Feature Call Concordance

### Primary 30-Guide Call Comparator

### Cell Set Overlap

| Metric | Value |
|--------|-------|
| STAR filtered cells | 30,497 |
| CR filtered cells | 32,256 |
| Shared filtered cells | 30,417 |
| Filtered-barcode Jaccard | 0.9407 |

### Guide Presence on Shared Cells

| Category | Count |
|----------|-------|
| Both uncalled | 6,798 |
| STAR-only called | 193 |
| CR-only called | 86 |
| Both called | 23,340 |

### Per-Cell Guide Assignment (23,340 shared called cells)

| Metric | Value |
|--------|-------|
| Same `0 / 1 / >1 guides` bucket | 0.9870 |
| Same `num_features` | 0.9866 |
| Exact called guide-set match | 23,210 / 23,340 (0.9944) |
| Exact singleton match | 0.99995 |
| Exact multi-guide set match | 0.9506 |
| Mean per-cell guide-set Jaccard | 0.9973 |
| Mean cellwise CR-in-STAR containment | 0.9991 |
| Mean cellwise STAR-in-CR containment | 0.9981 |
| Per-cell total guide UMI Pearson | 0.9997 |

Directional called-cell containment:

| Metric | Value |
|--------|-------|
| CR-in-STAR called-cell containment | 23,340 / 23,426 = 0.9963 |
| STAR-in-CR called-cell containment | 23,340 / 23,533 = 0.9918 |
| Overlap coefficient | 0.9963 |

The largest raw mismatches are:

- `108` cells: `STAR=RREB1`, `CR=None`
- `50` cells: `STAR=None`, `CR=NOTCH1`
- `13` cells: `STAR=RFX3`, `CR=None`
- a small tail of `STAR=<guide>|RREB1`, `CR=<guide>` and similar
  secondary-guide disagreements

This is the paper-quality 30-guide CRISPR call comparison.

### Historical 29-Guide Call Comparison

If STAR calls are projected into CellRanger's current 29-guide CRISPR
space by dropping `SIX2` before comparison, concordance improves to:

| Metric | Value |
|--------|-------|
| Both uncalled | 7,606 |
| STAR-only called | 194 |
| CR-only called | 87 |
| Both called | 22,530 |
| Exact called guide-set match (both called) | 22,405 / 22,530 (0.9945) |
| Exact singleton match | 20,101 / 20,102 (0.99995) |
| Exact multi-guide set match | 2,304 / 2,428 (0.9489) |

This historical 29-guide projection remains useful for interpreting the
older misconfigured CellRanger call run.

Directional containment in this 29-guide space:

| Metric | Value |
|--------|-------|
| CR-in-STAR called-cell containment | 22,530 / 22,617 = 0.9962 |
| STAR-in-CR called-cell containment | 22,530 / 22,724 = 0.9915 |
| Overlap coefficient | 0.9962 |
| Mean call-set Jaccard (both-called cells) | 0.9973 |
| Mean cellwise CR-in-STAR containment | 0.9991 |
| Mean cellwise STAR-in-CR containment | 0.9982 |

This captures the skew directly: CellRanger calls are almost entirely
contained within STAR, while STAR is slightly more permissive.

### Historical 29-Guide Count Parity (Old Full-Call CellRanger Run)

Dropping `SIX2` from STAR to match the current CellRanger CRISPR feature
space gives:

| Metric | STAR | CellRanger | STAR/CR |
|--------|------|------------|---------|
| Merged raw UMIs (29 guides) | 4,612,544 | 4,403,015 | 1.0476 |
| Merged raw guide-positive cells | 219,323 | 230,103 | 0.9532 |
| Merged filtered UMIs (29 guides) | 3,170,586 | 3,037,083 | 1.0440 |
| Merged filtered guide-positive cells | 30,243 | 31,997 | 0.9452 |

Per-guide Pearson in this 29-guide shared space: **0.999926**.

## CRISPR Calling Summary

### STAR Protospacer Calls (full run)

| Category | Cells | Pct |
|----------|-------|-----|
| 1 called feature | 20,968 | 68.7% |
| >1 called features | 2,573 | 8.4% |
| 0 called features | 6,956 | 22.8% |

### CellRanger 30-Guide Call Comparator

From `/storage/MSK-perturb-comparison/cr_full_grna_30crispr_20260306_173247/`:

- all 30 guides are typed `CRISPR Guide Capture`
- merged filtered guide UMIs: **3,111,361**
- merged filtered guide-positive cells: **32,017**
- filtered cell set: **32,256**
- called-cell containment versus STAR:
  - `CR in STAR = 0.9963`
  - `STAR in CR = 0.9918`

The old `/storage/MSK-perturb-comparison/cr_full_grna_withcalls_20260306_095454/`
run remains useful only as a historical 29-guide call-space comparator.

## Residual GEX Watchlist and Decision

The remaining GEX divergence is small and concentrated in a short watchlist of
ambiguity-prone genes, not in the broad GEX or gRNA signal:

- **CR-high**: `RPL39`, `RPL9`
- **STAR-high**: `STEAP1B`, `ENSG00000283674`, `ENSG00000280441`

We audited CellRanger's molecule-level provenance from:

- `outs/per_sample_outs/cr9_gex_grna_30crispr/count/sample_molecule_info.h5`
- `outs/multi/count/raw_molecule_info.h5`

Saved audit artifact:

```
/storage/MSK-perturb-comparison/paper_artifacts/msk_cr_molecule_provenance_20260307/
```

Key result: the CR-vs-STAR difference is mostly **not** driven by extra cells
or robust unique molecules that STAR completely drops. It is mainly in the
same shared filtered cells where both methods already have the gene.

Examples on the `30,417` shared filtered cells:

- `RPL39`: CR `302,865` vs STAR `44,218`
  - delta from cells where both methods already have `RPL39`: `209,607`
  - delta from cells where CR has `RPL39` and STAR has zero: `49,042`
- `RPL9`: CR `475,441` vs STAR `381,695`
  - delta from cells where both methods already have `RPL9`: `93,416`
  - delta from cells where CR has `RPL9` and STAR has zero: `331`
- `STEAP1B`: the opposite pattern; STAR is much higher in the same cells,
  not because CR is missing unique signal

The molecule support is weak in exactly the way expected for ambiguous short
genes:

- median reads per molecule: `1`
- p90 reads per molecule: `2`

Interpretation:

- this does **not** look like a cell-calling problem
- this does **not** look like a large UMI-collapse policy effect
- it is most consistent with subtle multimapper candidate-set / assignment
  differences on short ambiguity-prone genes
- differences in STAR version / alignment behavior remain a plausible residual
  source here

We tested a protein-coding-preference rescue heuristic in the CR rescue path.
It was reasonable in principle, but it made STAR slightly more permissive
without improving CellRanger compatibility on the MSK replay. We therefore
leave the multimapper rescue policy **unchanged**.

Conclusion: the residual watchlist difference is real but small, localized,
and not worth extra rescue complexity. Preserve the audit path, but keep the
current policy.

## Context: What Was Fixed

The canonical-TRU patch (this branch) ensures:
- NXT read barcodes are canonicalized to TRU at ingress in `assignBarcodes.c`
- NXT whitelists are canonicalized to TRU in `PfMultiAssign.cpp`
- Internal feature barcode state is always TRU
- Filtered-barcode normalization targets TRU
- Feature MEX entering merge is TRU (verified: namespace-assert shows 100% overlap)

This fixed the previous **zero** gRNA assignments from the namespace mismatch.

## UCSF Regression Gate

The UCSF iPSC2 full dataset is **stable** under canonical-TRU:
- GEX Jaccard: 0.9898
- Gene Pearson: 0.998
- Feature call set-equivalence: 99.44%

UCSF uses TRU-native chemistry, so the canonical-TRU change is a no-op there.

## Output Files

### STAR full MSK 3-lib run
```
/storage/MSK-perturb-comparison/canonical_tru_seq_20260306_052040/star_3lib/
```
- **GEX filtered**: `Solo.out/GeneFull/filtered/`
- **gRNA raw MEX**: `cr_assign/CRISPR_Guide_Capture/grna_de/PolyIII/` (stats.txt: 4,751,051 UMIs)
- **gRNA EmptyDrops filtered**: `cr_assign/CRISPR_Guide_Capture/grna_de/PolyIII/filtered/` (1,958,729 UMIs, 5,643 cells)
- **LARRY filtered MEX**: `cr_assign/Custom/larry_de/LARRY/filtered/`
- **Merged filtered MEX (GEX+gRNA+LARRY)**: `outs/filtered_feature_bc_matrix/` (3,247,179 guide UMIs, 30,497 cells)
- **Merged raw MEX**: `outs/raw_feature_bc_matrix/`
- **CRISPR calling**: `outs/crispr_analysis/protospacer_calls_per_cell.csv`
- **Timing**: 41:40 wall, 73.4 GB RSS

### CellRanger full MSK GEX+gRNA (with CRISPR calling, primary 30-guide comparator)
```
/storage/MSK-perturb-comparison/cr_full_grna_30crispr_20260306_173247/cr9_gex_grna_30crispr/
```
- **Filtered MEX**: `outs/per_sample_outs/cr9_gex_grna_30crispr/count/sample_filtered_feature_bc_matrix/`
- **Raw MEX**: `outs/multi/count/raw_feature_bc_matrix/`
- **CRISPR analysis**: `outs/per_sample_outs/cr9_gex_grna_30crispr/count/crispr_analysis/`
- **Feature ref**: `30 x CRISPR Guide Capture`
- **Guide counts**: 4,515,296 raw UMIs, 3,111,361 filtered UMIs

### CellRanger full MSK GEX+gRNA (historical 29-guide call-space run)
```
/storage/MSK-perturb-comparison/cr_full_grna_withcalls_20260306_095454/cr9_gex_grna_calls/
```
- **Filtered MEX**: `outs/per_sample_outs/cr9_gex_grna_calls/count/sample_filtered_feature_bc_matrix/`
- **Raw MEX**: `outs/multi/count/raw_feature_bc_matrix/`
- **CRISPR analysis**: `outs/per_sample_outs/cr9_gex_grna_calls/count/crispr_analysis/`
- **Feature ref**: `ref_feature_geneBC_crispr.csv` (29 `CRISPR Guide Capture` + 1 `Custom` feature, `SIX2`)
- **Timing**: 1:01:02 wall, 14.6 GB RSS
- **Config**: `create-bam,false`; no `no-secondary-analysis` (full pipeline including CRISPR calling)

## Reproducible Analysis Scripts

Saved under:

```
scripts/msk_grna_parity/
```

Files:

- `compare_guide_mex.py`
  - stage-by-stage guide-count parity on a common feature-ID set
- `compare_protospacer_calls.py`
  - per-cell call-set concordance, Jaccard, and directional containment
- `README.md`
  - methodology, formulas, primary 30-guide commands, and historical 29-guide commands

These scripts are intended to be the reproducible source for any
publication tables or parity figures derived from this handoff.

Saved 30-guide paper artifacts:

```
/storage/MSK-perturb-comparison/paper_artifacts/msk_grna_parity_30guide_20260306/
```

- `guide_counts/guide_stage_summary.tsv`
- `guide_counts/guide_per_feature_filtered.tsv`
- `guide_counts/guide_parity_summary.json`
- `call_concordance/call_concordance_summary.json`
- `call_concordance/call_pair_counts.tsv`
- `call_concordance/call_mismatch_rows.tsv`

Saved molecule-provenance artifact for the residual GEX watchlist:

```
/storage/MSK-perturb-comparison/paper_artifacts/msk_cr_molecule_provenance_20260307/
```

- `cr_molecule_provenance_watchlist.json`
- `README.txt`

### CellRanger full MSK GEX+gRNA (no secondary analysis)
```
/storage/MSK-perturb-comparison/canonical_tru_seq_cr_grna_20260306_052040/cr9_gex_grna/
```
- **Filtered MEX**: `outs/per_sample_outs/cr9_gex_grna/count/sample_filtered_feature_bc_matrix/`
- **No crispr_analysis** (disabled by `no-secondary-analysis,true`)
- **Timing**: 58:16 wall, 86.3 GB true RSS

### CellRanger full MSK GEX+LARRY
```
/storage/MSK-perturb-comparison/canonical_tru_seq_cr_larry_20260306_052040/cr9_gex_larry/
```
- **Filtered MEX**: `outs/per_sample_outs/cr9_gex_larry/count/sample_filtered_feature_bc_matrix/`
- **Timing**: 1:49:45 wall, 87.8 GB true RSS

### UCSF full STAR (regression gate)
```
/storage/ucsf-full/paper_benchmarks/canonical_tru_seq_20260306_052040/
```
- **Parity report**: `parity/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
- **Timing**: 11:03 wall, 51.2 GB RSS

## Timing & Memory Summary

| Run | Wall Clock | Peak RSS | CPU % |
|-----|-----------|----------|-------|
| **STAR 3-lib** (GEX+gRNA+LARRY) | 41:40 | 73.4 GB | 1919% |
| **UCSF STAR** (GEX+guides) | 11:03 | 51.2 GB | 2437% |
| **CR GEX+gRNA** (historical 29-guide calls) | 1:01:02 | 14.6 GB | 1884% |
| **CR GEX+gRNA** (no secondary) | 58:16 | 86.3 GB | 1941% |
| **CR GEX+LARRY** | 1:49:45 | 87.8 GB | 1897% |

Machine: pikachu, i9-13900KF, 126 GB, 32 threads.
All runs sequential (no contention), `/usr/bin/time -v` for STAR,
Martian `_perf` `highmem.rss` for CR true peak.

## Errata: Previous "63% Undercount"

The original version of this document reported STAR at 0.63x CR guide UMIs.
That number compared STAR's gRNA-only EmptyDrops output (1,958,729 UMIs on
5,643 cells) against CellRanger's merged filtered matrix (3,111,361 UMIs on
32,256 cells). These are different output stages and not comparable.

The correct merged-filtered comparison is **1.044x** (3,247,179 / 3,111,361).
The abundance-dependent ratio pattern in the original table was entirely an
artifact of the gRNA-only EmptyDrops cell set being biased toward
high-expression guides.

## Actionable Guidance For Future Agents

When debugging MSK 30-KO gRNA parity:

- First decide whether you are checking **raw assignment**, **merged guide
  counts**, or **CRISPR calling**. These are different outputs.
- For merged guide counts, ignore `cr_assign/.../filtered/` and use the
  merged matrices.
- Use `/storage/MSK-perturb-comparison/cr_full_grna_30crispr_20260306_173247/`
  as the primary 30-guide CellRanger comparator.
- Keep `/storage/MSK-perturb-comparison/cr_full_grna_withcalls_20260306_095454/`
  only as a historical 29-guide call-space reference.
