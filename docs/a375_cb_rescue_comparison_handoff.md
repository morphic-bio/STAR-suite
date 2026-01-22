# A375 CB-Rescue Parity (1k Read Pair) Handoff

This captures the **correct comparison method** for Cell Ranger vs `process_features` on the **1k-read L001 pair** and documents the key details that were previously missed.

## Correct comparison setup (apples-to-apples)

### Use a single L001 pair (not multi-lane)
Using all lanes can exceed 1,000 counts even when “downsampled_1000” is intended. Use the single-pair downsample:

- CRISPR FASTQs (single pair, 1k reads):
  - `/storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/crispr/downsampled_1000_L001_pair`
- GEX FASTQs (single pair, 1k reads):
  - `/storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/gex/downsampled_1000_L001_pair`

### CR run and CRISPR-only MEX
- Cell Ranger config: `a375_crispr_1000_pair_multi_config.csv`
- Cell Ranger run: `/tmp/a375_crispr_1000_pair`
- CRISPR-only raw MEX (separated):
  - `/tmp/a375_crispr_1000_pair/outs/multi/count/raw_feature_bc_matrix_crispr_only`

### Eliminate zero-count barcodes before comparing
CR raw matrices include many zero-count barcodes; `process_features` only emits barcodes with nonzero counts. Create a “nonzero-only” CR matrix before comparing:

- `/tmp/a375_crispr_1000_pair/outs/multi/count/raw_feature_bc_matrix_crispr_only_nonzero`

Raw CRISPR-only MEX stats (raw):
- barcodes = 914 total
- barcodes with counts = 475
- zero-count barcodes = 439

Resulting CR stats (nonzero-only):
- features = 11
- barcodes = 475
- total_counts = 772

## `process_features` run (parity settings)
Use `limit_search=-1` so offsets are not limiting and set `MIN_COUNTS=0` to allow single-read UMIs.

Expected command shape:

```
MIN_COUNTS=0 \
FASTQ_DIR=/storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/crispr/downsampled_1000_L001_pair \
FILTERED_BARCODES=/tmp/a375_crispr_1000_pair/outs/multi/count/raw_feature_bc_matrix_crispr_only_nonzero/barcodes.tsv \
MAX_HAMMING_LIST="1" FEATURE_N_LIST="1" BARCODE_N_LIST="1" \
FEATURE_OFFSET_LIST="0" LIMIT_SEARCH_LIST="-1" \
bash scripts/a375_feature_sweep.sh
```

Output:
- `/tmp/a375_feature_1000_L001_pair_mc0/m1_fn1_bn1_off0_lsn1/downsampled_1000_L001_pair`

## Why barcodes differ (root cause)
This was not an offset issue:
- `feature_sequences.txt` showed match position 63, and `limit_search=-1` means the full read is searched.
- Logs showed ~906/1000 reads matched to features, only ~72 unmatched.

The major differences were:
1) **Zero-count barcodes in CR raw**: CR raw matrix includes many barcodes with no counts.
2) **UMI clique filtering**: default `min_counts=1` drops singletons in tiny 1k runs.
   - `assignBarcodes.log` showed **471 whitelisted barcodes** but only **“Processed barcodes: 111”**.
   - That drop is mostly singleton UMIs being filtered when `min_counts=1` (UMI needs >1 read).

## Effect of `MIN_COUNTS`
- `min_counts=1`: `process_features` collapsed to ~111 barcodes / 124 counts.
- `min_counts=0`: `process_features` rises to **471 barcodes / 765 counts** (very close to CR nonzero **475 / 772**).

Remaining gap:
- Common barcodes: 464
- Missing from `process_features`: 11 (all count=1 in CR)
- Extra in `process_features`: 7 (counts 1–2)

## What not to do
- Do **not** compare 1k-downsampled STAR runs to full-depth CR matrices.
- Do **not** compare raw CR matrices without first removing zero-count barcodes.

## Next step to close the remaining gap
Sweep CB rescue knobs on the L001 1k pair with `MIN_COUNTS=0` and CR nonzero barcodes forced as the filtered set:
- `barcode_n`: 0–2
- `max_barcode_mismatches`: 0–3

Use:
- `FILTERED_BARCODES=/tmp/a375_crispr_1000_pair/outs/multi/count/raw_feature_bc_matrix_crispr_only_nonzero/barcodes.tsv`
- `FASTQ_DIR=/storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/crispr/downsampled_1000_L001_pair`

This isolates CB rescue effects without UMI/min-count noise.
