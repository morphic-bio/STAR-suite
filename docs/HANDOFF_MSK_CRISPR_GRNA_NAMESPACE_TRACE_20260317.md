# MSK gRNA Namespace Trace Handoff

**Date**: 2026-03-17  
**Status**: TARGETED ROOT CAUSE IDENTIFIED FOR BAD POLYIII PATH

## Summary

The remaining catastrophic MSK PolyIII regression is not just a late merge/filter
namespace bug.

A targeted barcode-level trace shows an earlier failure in
`checkAndCorrectBarcode()`:

- the raw guide barcode enters in `NXT` space
- the code can see that the translated `TRU` barcode is in the whitelist
- but it does **not** translate before exact whitelist lookup
- instead it falls into ordinary 1-mismatch rescue
- that rescue rewrites the barcode to the wrong `TRU` neighbor

This is a silent **misassignment**, not just a dropout.

## Concrete Traced Example

Target barcode pair:

- raw FASTQ barcode (`NXT`): `AACGAAAAGAGCGTCG`
- true translated whitelist barcode (`TRU`): `AACGAAATCAGCGTCG`

Observed wrong rescue target:

- wrong internal rescue barcode (`TRU`): `AACGAAAAGAGCATCG`
- wrong final output barcode after `translate_NXT` output rewrite: `AACGAAATCAGCATCG`

What the trace showed:

- `raw_hit=0`
- `translated_hit=1`
- `n_matches=1`
- chosen rescue candidate: `AACGAAAAGAGCATCG`

So the code knows the translated barcode would hit the whitelist, but still
chooses a generic 1MM rescue candidate instead of doing deterministic
namespace normalization first.

## Evidence

### Trace Run

Focused bad-path reproduction:

- output: `/tmp/msk_bad_trace_ns_20260317`
- config: `/tmp/msk_bad_trace_ns_20260317.multi_config.csv`
- shape: `GEX + PolyIII`, no LARRY, same bad whitelist/control surface

Key trace file:

- `/tmp/msk_bad_trace_ns_20260317/namespace_trace.log`

Representative lines:

```text
NS_CHECK raw=AACGAAAAGAGCGTCG translated=AACGAAATCAGCGTCG translate_NXT=1 feat=18 match_pos=20 raw_hit=0 translated_hit=1
NS_RESCUE raw=AACGAAAAGAGCGTCG translated=AACGAAATCAGCGTCG n_matches=1
  cand[0]=AACGAAAAGAGCATCG mismatch_pos=12
NS_DECISION raw=AACGAAAAGAGCGTCG decision=unique_rescue chosen=AACGAAAAGAGCATCG
```

Counts from the trace run:

- `NS_CHECK` seen for this barcode: `24`
- `unique_rescue chosen=AACGAAAAGAGCATCG`: `24`
- corrected reads written with `bc=AACGAAAAGAGCATCG`: `24`

Supporting files:

- namespace trace: `/tmp/msk_bad_trace_ns_20260317/namespace_trace.log`
- corrected read trace: `/tmp/msk_bad_trace_ns_20260317/reads_trace.log`
- note: `/tmp/msk_bad_trace_ns_20260317/rescue_trace.log` is empty because this
  never reaches posterior/pending rescue; the error is earlier in unique 1MM rescue

### Output-Level Confirmation

Good run:

- `/tmp/msk_100k_ucsf_guideflags_nolarry_20260317_1912/cr_assign/CRISPR_Guide_Capture/grna_de/PolyIII/feature_per_cell.csv`

Bad traced run:

- `/tmp/msk_bad_trace_ns_20260317/cr_assign/CRISPR_Guide_Capture/grna_de/PolyIII/feature_per_cell.csv`

Observed result:

- good run contains `AACGAAAAGAGCGTCG` with `(num_features=2, top_feature=18, total_deduped_umi=25)`
- bad run does **not** contain `AACGAAAAGAGCGTCG`
- bad run instead contains `AACGAAATCAGCATCG` with `(num_features=2, top_feature=18, total_deduped_umi=24)`

This is the same logical cell getting rewritten to the wrong barcode.

## Why This Matters

This trace explains why the bad MSK PolyIII path can still recover some guide
signal while remaining catastrophically wrong:

- it is not performing deterministic `NXT -> TRU` normalization for assignment
- it is letting the large `TRU` whitelist act like a giant near-neighbor search space
- many `NXT` guide reads can then be uniquely or ambiguously rescued to the wrong
  `TRU` barcodes
- that matches the bad run’s huge `Resolve_calls_total` and poor barcode parity

## Relevant Code

Instrumentation added here:

- `core/features/process_features/src/assignBarcodes.c`
  - trace helpers: around lines `439-479`
  - barcode trace inside `checkAndCorrectBarcode()`: around lines `4002-4099`

Key behavioral point:

- `checkAndCorrectBarcode()` currently does exact lookup on the raw barcode first
- if it misses, it immediately falls into `find_closest_barcodes()`
- it does **not** first use the known `translate_NXT` relation to transform the
  barcode into whitelist space for exact matching

## Reproduction Command

This is the focused command used for the trace:

```bash
OUT=/tmp/msk_bad_trace_ns_20260317
CFG=/tmp/msk_bad_trace_ns_20260317.multi_config.csv
awk 'NR<=4' /tmp/msk30ko_100k_repair_20260317_182315/multi_config.csv > "$CFG"

export PF_TRACE_NAMESPACE=1
export PF_TRACE_NAMESPACE_BARCODES='AACGAAAAGAGCGTCG,AACGAAATCAGCGTCG,AACGAAAAGAGCATCG,AACGAAATCAGCATCG'
export PF_TRACE_NAMESPACE_LOG="$OUT/namespace_trace.log"
export PF_TRACE_RESCUE=1
export PF_TRACE_RESCUE_BARCODES='AACGAAAAGAGCGTCG,AACGAAATCAGCGTCG,AACGAAAAGAGCATCG,AACGAAATCAGCATCG'
export PF_TRACE_RESCUE_LOG="$OUT/rescue_trace.log"
export PF_TRACE_READS=1
export PF_TRACE_READS_BARCODES='AACGAAAAGAGCGTCG,AACGAAATCAGCGTCG,AACGAAAAGAGCATCG,AACGAAATCAGCATCG'
export PF_TRACE_READS_LOG="$OUT/reads_trace.log"

./core/legacy/source/STAR \
  --runMode alignReads \
  --runThreadN 32 \
  --dynamicThreadInterface 1 \
  --genomeDir /storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-cellranger/star \
  --readFilesIn /storage/MSK-perturb-comparison/fixtures/msk_30ko_100k_l001/mRNA/mRNA_L001_R2_001.fastq.gz /storage/MSK-perturb-comparison/fixtures/msk_30ko_100k_l001/mRNA/mRNA_L001_R1_001.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix "$OUT/" \
  --outSAMtype None \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 \
  --soloUMIstart 17 \
  --soloCBlen 16 \
  --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist /storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt \
  --soloFeatures Gene GeneFull \
  --soloUMIdedup 1MM_CR \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloCellFilter EmptyDrops_CR \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloMultiMappers Rescue \
  --soloCbUbRequireTogether no \
  --soloCrGexFeature GeneFull \
  --soloCrMultimapRescue yes \
  --pfMultiConfig "$CFG" \
  --crChemistry auto \
  --crMinUmi 2 \
  --crAssignMaxHamming 1 \
  --crAssignConsumerThreads -1 \
  --crAssignSearchThreads 1 \
  --defaultCrCompat yes
```

## Recommended Next Step

Do **not** spend more time on downstream namespace plumbing for this specific
failure until this is fixed.

The next code change should be one of:

1. Ensure PolyIII assignment uses an `NXT` whitelist in the bad/default path.
2. Or, if `translateNxtForAssign` remains part of the design, add a deterministic
   `NXT -> TRU` normalization step **before** exact barcode lookup and before
   1MM rescue in `checkAndCorrectBarcode()`.

The key constraint:

- generic mismatch rescue must not stand in for namespace translation

## Validation After the Fix

Re-run the focused trace first and verify:

- `raw_hit=0 translated_hit=1` no longer falls into wrong `unique_rescue`
- the correct barcode `AACGAAAAGAGCGTCG` appears in
  `feature_per_cell.csv`
- the wrong barcode `AACGAAATCAGCATCG` disappears

Then re-run the broader 100K MSK guide fixture and watch:

- `Total barcodes`
- `Total feature counts`
- `Resolve_calls_total`
- final `protospacer_calls_summary.csv`

## Related Docs

- runbook: `docs/RUNBOOK_MSK_CRISPR_MASTER_REPAIR_20260317.md`
- original regression handoff:
  `docs/HANDOFF_MSK_CRISPR_REGRESSION_DEBUG_20260317.md`
