# MSK 30polyKO CRISPR Parity Regression — Root Cause Analysis

**Date**: 2026-03-17
**Status**: ROOT CAUSE IDENTIFIED — bad squash merge; `master` missing 15 commits from `CR-Larry-perturb`

## Problem Statement

Running the MSK 30polyKO paper benchmark (`run_msk_30polyko_benchmark.sh`) produced
catastrophically low CRISPR feature calling: 28,550 of 30,520 filtered cells had
**zero** CRISPR molecules, and only 394 cells received a single-feature call (vs.
20,968 in the reference run).

## Root Cause

The paper benchmark script ran the **`master` branch binary** (`2f0bd9cce7`), while
the original good MSK benchmark ran the **`CR-Larry-perturb` branch binary**
(`f16f30995d`). The `master` branch is missing 15 commits that contain the entire
feature assignment pipeline overhaul:

```
92efead docs: record MSK residual watchlist conclusion
463004e docs: add local agent reference override
f16f309 CPU-pressure permit bypass, benchmark scripts, NXT/TRU namespace bug analysis
5898302 Phase 0 async overlap, fix synthetic readId CB conflict, review hardening
9cc1670 MEX merge optimization: vector remap, streaming gzip, in-memory GEX cache
aa99976 Skip heatmaps in CR-compat mode, add ratio-based lineage calling for Custom libraries
b0db247 Bootstrap replay buffer and tiered search integration for feature assignment
a917996 FeatureN=1 default, ratio-based producer/consumer thread apportioning
ec04452 Add per-library star_max_hamming column for independent Hamming distance control
faa9f3b Integer-key prehash refactor: fix OOM on large feature sets, add tiered search order
384de1d Fix translateNxt namespace for filtered barcodes, enable MC rescue in bootstrap mode, and harden Flex loader
09444ab Strict namespace normalization: fix union bypass, partial-spec hard error, feature_barcodes opt-in, and E2E test
19f24fe Fix prehash ambiguity bookkeeping, add union whitelist support, harden regression harness
c5b9169 process_features: uncap PF consumers, fix d2 prehash budget, and clarify limit_search hard-bound CLI
2e0b93c docs: strengthen clean-rebuild policy in AGENTS.md
```

Key code-level changes missing from `master` (total ~3,500 lines in feature assignment):

| File | Lines changed |
|---|---|
| `core/features/process_features/src/assignBarcodes.c` | +1784 |
| `core/legacy/source/PfMultiProcess.cpp` | +1000 |
| `core/legacy/source/PfMultiMerge.cpp` | +531 |
| `core/legacy/source/PfMultiProcess.h` | +80 |
| `core/legacy/source/PfMultiMerge.h` | +55 |
| `core/legacy/source/PfMultiAssign.cpp` | +25 |

## Evidence: assignBarcodes stats comparison

| Metric | GOOD (`CR-Larry-perturb`) | BAD (`master`) | Delta |
|---|---|---|---|
| Total feature counts | 100,505,763 | 33,734,677 | **-66%** |
| Total deduped feature counts | 4,751,051 | 2,088,149 | -56% |
| Reads assigned to barcode % | 95.14% | 86.76% | -8.4pp |
| Total barcodes | 238,726 | 198,633 | -17% |
| Resolve no_hit | 0 | 190,963 | **new failures** |
| Raw matrix entries | 1,312,315 | 631,093 | -52% |

## Evidence: CRISPR calling comparison

| Metric | GOOD (`CR-Larry-perturb`) | BAD (`master`) |
|---|---|---|
| Cells with 0 molecules | 234 | 28,550 |
| Cells with no confident call | 6,722 | 1,279 |
| Cells with 1 feature | 20,968 | 394 |
| Cells with >1 features | 2,573 | 297 |
| Total filtered cells | 30,497 | 30,520 |

## Additional Parameter Differences (secondary)

The paper script also introduced these parameters not in the original run:

| Parameter | Original | Paper script |
|---|---|---|
| `--genomeDir` | `refdata-gex-GRCh38-autoindex11044-crstar/star` | `bulk_index` |
| `--clip3pPolyG` | (not set) | `yes` |
| `--dynamicThreadConstMapPermits` | (not set) | `32` |
| `--crAssignSkipQcOutputs` | (not set) | `1` |

These are **not** the cause of the CRISPR regression — the genome index and poly-G
trimming affect GEX mapping, not feature barcode assignment. The `crAssignSkipQcOutputs`
only skips histogram generation. The core issue is the binary itself.

## Resolution Plan

1. **Switch to `CR-Larry-perturb` branch**, clean rebuild, re-run MSK benchmark
2. OR: Squash-merge `CR-Larry-perturb` → `master` properly, clean rebuild, re-run
3. Verify CRISPR stats match the reference run before proceeding to parity

## Verification Command

After rebuilding on the correct branch:
```bash
# Check the binary reports the correct branch
head -2 Log.out | grep "STAR git"
# Should show: On branch CR-Larry-perturb (or master with all commits merged)

# Quick sanity: stats.txt should show ~100M total feature counts
cat cr_assign/CRISPR_Guide_Capture/grna_de/PolyIII/stats.txt | head -3
```

## Run Paths

- **GOOD run**: `/storage/MSK-perturb-comparison/canonical_tru_seq_20260306_052040/star_3lib/`
- **BAD run**: `/storage/MSK-perturb-comparison/paper_bench_20260317_152518/`
- **Paper script**: `/mnt/pikachu/STAR-suite/scripts/paper/run_msk_30polyko_benchmark.sh`
