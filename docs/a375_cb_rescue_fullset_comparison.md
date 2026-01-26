# A375 CB-Rescue Parity (Full Set) Comparison

This documents the full set `process_features` run with the same settings used in the 1k L001 pair parity check, plus a CR comparison derived from the Cell Ranger molecule file.

## Inputs
- CRISPR FASTQs (full set, all lanes): `/storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/crispr`
- Feature reference: `/storage/A375/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv`
- Whitelist: `/storage/A375/3M-5pgex-jan-2023.txt`

## Settings (same as 1k parity)
- `maxHammingDistance=1`
- `feature_n=1`
- `barcode_n=2`
- `max_barcode_mismatches=5`
- `feature_constant_offset=0`
- `limit_search=-1`
- `min_counts=0`

## CR comparison source (full set)
Cell Ranger molecules file used to derive CRISPR-only stats and nonzero barcodes:
- `/storage/A375/outputs/1k_CRISPR_5p_gemx_count_sample_molecule_info.h5`

Derived CR CRISPR-only stats:
- CRISPR features: 11
- CRISPR UMIs (total counts): 3,170,423
- Nonzero CRISPR barcodes: 108,867
- Nonzero barcode list: `/tmp/a375_crispr_full_nonzero_barcodes.tsv`

## `process_features` full set (unfiltered)
Command:
```
MIN_COUNTS=0 \
FASTQ_DIR=/storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/crispr \
MAX_HAMMING_LIST="1" FEATURE_N_LIST="1" BARCODE_N_LIST="2" \
FEATURE_OFFSET_LIST="0" LIMIT_SEARCH_LIST="-1" \
MAX_BARCODE_MISMATCHES_LIST="5" \
OUT_ROOT=/tmp/a375_feature_full_m1_mbm5 \
bash scripts/a375_feature_sweep.sh
```

Summary:
- Output: `/tmp/a375_feature_full_m1_mbm5/summary.tsv`
- barcodes: 110,904
- features: 11
- total_counts: 3,225,308
- whitelisted_barcodes: 110,916
- unmatched_reads: 345,083
- percent_assigned: 96.4142

## `process_features` full set (CR nonzero barcodes)
Command:
```
MIN_COUNTS=0 \
FASTQ_DIR=/storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/crispr \
FILTERED_BARCODES=/tmp/a375_crispr_full_nonzero_barcodes.tsv \
MAX_HAMMING_LIST="1" FEATURE_N_LIST="1" BARCODE_N_LIST="2" \
FEATURE_OFFSET_LIST="0" LIMIT_SEARCH_LIST="-1" \
MAX_BARCODE_MISMATCHES_LIST="5" \
OUT_ROOT=/tmp/a375_feature_full_m1_mbm5_filtered \
bash scripts/a375_feature_sweep.sh
```

Summary:
- Output: `/tmp/a375_feature_full_m1_mbm5_filtered/summary.tsv`
- barcodes: 108,817
- features: 11
- total_counts: 3,222,884
- whitelisted_barcodes: 110,916
- unmatched_reads: 345,083
- percent_assigned: 96.4142

## Delta vs CR (CRISPR-only)
- CR barcodes: 108,867
- STAR barcodes: 108,817 (delta -50)
- CR UMIs: 3,170,423
- STAR counts: 3,222,884 (delta +52,461; ratio 1.0165)
