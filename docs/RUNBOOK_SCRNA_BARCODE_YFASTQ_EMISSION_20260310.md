# Runbook: scRNA Barcode Y/noY FASTQ Emission

Date: March 10, 2026

## Goal

Extend Y/noY FASTQ emission for scRNA/Flex runs so that the corresponding
barcode read is emitted alongside the cDNA read:

- cDNA/GEX read goes to `Y` or `noY` based on the current Y-routing decision
- barcode read goes to the same side as its corresponding cDNA read

Target behavior:

- bulk PE remains unchanged
- scRNA/Flex with separate barcode read emits both `R2` and barcode `R1`
- existing Y/noY BAM behavior remains unchanged

## Why This Is A Core Change

This is not a `libscrna` task.

Reason:

- `libscrna` starts after FASTQ ingestion and alignment output decisions
- Y/noY FASTQ routing is decided while STAR still owns:
  - original input ends
  - per-read routing state
  - output path derivation
  - chunk-level FASTQ writers

The current omission is caused by core Solo logic excluding the barcode read
from `readNmates`, and the Y/noY emission path iterating over `readNmates`.

## Current Root Cause

Solo setup intentionally excludes the barcode read from aligned mates:

- `ParametersSolo.cpp`
  - `readNmates = readNends - 1`
  - barcode read is the last read for separate-barcode Solo modes

Y/noY FASTQ emission then uses `readNmates` in all major lifecycle stages:

- path derivation in `Parameters.cpp`
- batch reset path derivation in `Parameters_batchReset.cpp`
- chunk file open/reopen in `ReadAlignChunk.cpp`
- final gzip concatenation in `STAR.cpp`
- actual write loop in `ReadAlign_outputAlignments.cpp`

Result:

- bulk PE emits both mates
- scRNA/Flex emits only the aligned cDNA read

## Scope

In scope:

- Solo/Flex runs with separate barcode read
- emitting barcode-read `R1` plus cDNA `R2` under `y_separated/`
- batch-mode compatibility
- gzip and uncompressed emission paths

Out of scope for this change:

- `libscrna`
- BAM routing logic
- `YReadNames` logic
- pfMulti feature-library FASTQ emission policy beyond the main aligned
  cDNA/barcode pair

## Proposed Design

### 1. Introduce a Core Helper For Y/noY FASTQ Emit End Count

Add one core helper or derived count, conceptually:

- bulk PE: emit `readNmates`
- Solo/Flex separate-barcode mode: emit `readNends`

This should be the single source of truth for Y/noY FASTQ emission loops.

Suggested placement:

- `Parameters` or `ParametersSolo`

The value must be available to:

- path derivation
- batch reset
- chunk open/reopen
- final concat
- actual write loop

### 2. Use That Helper Everywhere Y/noY FASTQ Loops Currently Use `readNmates`

Touch points:

- `core/legacy/source/Parameters.cpp`
- `core/legacy/source/Parameters_batchReset.cpp`
- `core/legacy/source/ReadAlignChunk.cpp`
- `core/legacy/source/STAR.cpp`
- `core/legacy/source/ReadAlign_outputAlignments.cpp`

Important detail:

- arrays are already sized to `MAX_N_MATES=3`
- this change is about loop bounds and path derivation, not storage expansion

### 3. Route Barcode Read With Its Corresponding cDNA Read

For Solo/Flex route-per-read mode:

- keep the Y/noY decision based on the aligned cDNA read
- write all emitted ends for that record to the same side

That means:

- if current cDNA read is Y, emit both cDNA `R2` and barcode `R1` to `Y`
- otherwise emit both to `noY`

This preserves the current routing model and avoids introducing any new
barcode-specific Y decision.

### 4. Preserve Existing Header Behavior

Current scRNA/Flex emitted FASTQs strip header comments.

Do not change that in this patch.

The new barcode outputs should follow the existing emitted FASTQ formatting,
not the standalone `remove_y_reads` formatting.

## Test Plan

### Extend Existing New Smokes

1. `tests/run_flex_yremove_smoke.sh`

Add expectations for:

- `*_Y_R1_001.fastq.gz`
- `*_noY_R1_001.fastq.gz`

Reference:

- run `remove_y_reads` on both cDNA `R2` and barcode `R1`
- aggregate outputs in original lane order
- compare emitted `R1` and `R2` against reference with:
  - normalized header comparison
  - exact record order
  - exact sequence
  - exact quality

2. `tests/run_perturb_yremove_smoke.sh`

Add the same check for the main GEX barcode read paired with the emitted GEX
cDNA stream.

Keep the existing perturb guard:

- `outs/crispr_analysis/*` must still exist

### Keep Existing Invariants

The following must remain true:

- `noY.bam` contains no `chrY` alignments
- Y/noY BAM split still passes
- existing bulk PE Y/noY tests still pass
- batch multisample PE equivalence still passes

## Validation Sequence

1. Clean rebuild

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

2. Run targeted smokes

```bash
bash tests/run_flex_yremove_smoke.sh
bash tests/run_perturb_yremove_smoke.sh
```

3. Re-run existing PE regressions

```bash
bash tests/run_public_bulk_pe_smoke.sh
bash tests/run_pe_bulk_regression_replay.sh
bash tests/run_bulk_pe_multisample_equivalence.sh
```

## Acceptance Criteria

The patch is complete when all of the following are true:

1. Flex emits barcode-read Y/noY FASTQs.
2. Perturb emits barcode-read Y/noY FASTQs for the main GEX pair.
3. Emitted scRNA/Flex Y/noY FASTQs match `remove_y_reads` on normalized
   header, order, sequence, and quality.
4. No bulk PE behavior changes.
5. No new batch-mode regressions.

## Main Risk

The only real risk is mixing:

- aligned mate count (`readNmates`)
- physical input end count (`readNends`)

If the change is only applied in one layer, output paths and chunk concatenation
will diverge.

This is why the helper/count must be adopted consistently across all five core
touch points listed above.

