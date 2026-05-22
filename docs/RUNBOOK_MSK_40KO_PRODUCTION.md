# MSK 40KO Production Runbook

Date: 2026-05-21
Status: scripts staged; whitelist preflight and command dry-run passed. The
production wrapper supports asynchronous remote downstream submission and
Globus BAM transfer.

## Dataset

Primary FASTQ root: `/mnt/pikachu/scRNAseq_40KO`.

The path supplied as `/mnt/scRNAseq_40KO` is not mounted on this host; the data
are under `/mnt/pikachu/scRNAseq_40KO`.

Observed logical samples:

| Sample | Libraries | FASTQ size |
| --- | --- | ---: |
| `40_KO_ES` | GEX, PolyIII/gRNA, LARRY | 58.1 GiB |
| `40_KO_DE` | GEX, PolyIII/gRNA, LARRY | 62.6 GiB |

The FASTQs are flat, not separated by library directory. The production wrapper
therefore creates per-sample symlink directories under
`<out-root>/samples/<sample>/staged_fastqs/{mRNA,PolyIII,LARRY}` before running
STAR.

Read layout matches the 30KO production surface:

- `I1=10`, `I2=10`, `R1=29`, `R2=89`
- GEX and LARRY use February-2018 TRU barcodes
- PolyIII/gRNA uses February-2018 NXT barcodes

Initial preflight rejected the GEM-X whitelist hypothesis: all six rows called
the February-2018 family. Do not use the 30KO DE-GemX May-2023/GEM-X
whitelists for this dataset.

## CRISPR Methodology

Use the same CR-compatible PolyIII guide semantics as the 30KO DE-GemX run, but
with the February-2018 NXT barcode namespace:

- `library_type=CRISPR Guide Capture`
- `feature_types=CRISPR Guide Capture`
- feature reference:
  `/mnt/pikachu/MSK-whitelists/ref_feature_geneBC_crispr.csv`
- `--crMinUmi 2`
- `--crAssignFeatureOffset 0`
- `--crOutputChemistry TRU`
- NXT guide whitelist for STAR assignment:
  `/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt`

Targeted R2 motif validation on 50,000 representative reads per library supports
this routing:

| Library | CRISPR scaffold, <=2 mismatches | LARRY scaffold, <=2 mismatches |
| --- | ---: | ---: |
| ES GEX | 0.000 | 0.001 |
| ES PolyIII | 0.966 | 0.000 |
| ES LARRY | 0.000 | 0.913 |
| DE GEX | 0.000 | 0.001 |
| DE PolyIII | 0.957 | 0.000 |
| DE LARRY | 0.000 | 0.898 |

The manifest records raw whitelist expectations separately from STAR assignment
whitelists so preflight can validate the observed barcode namespace:

`docs/MSK_40KO_FASTQ_MANIFEST.tsv`

## Preflight

Run the whitelist-family preflight first:

```bash
scripts/run_msk_40ko_fastq_preflight.sh \
  --outdir plans/artifacts/msk_40ko_fastq_preflight_$(date +%Y%m%d_%H%M%S)
```

Expected result: all six rows pass with `feb2018:TRU` for GEX/LARRY and
`feb2018:NXT` for PolyIII/gRNA.

## Dry Run

Generate staged symlinks, `pf_multi_config.csv`, and `RUN_COMMAND.sh` without
starting STAR:

```bash
scripts/run_msk_40ko_pipeline_from_manifest.py \
  --out-root /storage/MSK-perturb-comparison/msk40ko_dryrun_$(date +%Y%m%d_%H%M%S) \
  --dry-run
```

Inspect:

- `<out-root>/samples/40_KO_ES/pf_multi_config.csv`
- `<out-root>/samples/40_KO_DE/pf_multi_config.csv`
- `<out-root>/samples/<sample>/RUN_COMMAND.sh`

## Production

Run one sample at a time unless explicitly scheduling on a dedicated machine:

```bash
scripts/run_msk_40ko_pipeline_from_manifest.py \
  --out-root /storage/MSK-perturb-comparison/msk40ko_prod_$(date +%Y%m%d_%H%M%S) \
  --samples 40_KO_ES \
  --threads 32 \
  --star-bin /mnt/pikachu/STAR-suite/core/legacy/source/STAR \
  --genome-dir /storage/autoindex_110_44/bulk_index \
  --out-samtype "BAM Unsorted" \
  --run-downstream \
  --remote-host 10.159.4.53 \
  --remote-root /tmp/msk40ko_cellbender \
  --downstream-async \
  --run-cellbender \
  --cellbender-cpu-cores 8
```

Then run `40_KO_DE` after `40_KO_ES` completes.

The wrapper uses:

- `--soloFeatures GeneFull Velocyto`
- `--soloCellFilter EmptyDrops_CR`
- `--soloCrMultimapRescue yes`
- `--dynamicThreadInterface 1`
- `--crAssignConsumerThreads -1`
- `--crAssignSearchThreads 1`
- `STAR_VELOCYTO_LOW_MEM=1`
- `STAR_SOLO_BINARY_SPOOL=1`

No Y-removal is enabled, matching the MSK 30KO production policy.

## Globus BAM Transfer

If BAMs should be uploaded after each sample, add:

```bash
  --globus-src-endpoint 07446cad-33b8-11f0-8c0c-0afffb017b7d \
  --globus-dst-endpoint 61fb8b9a-9b52-456e-928c-30c0fb0140bf \
  --globus-dst-root /MSK-40KO-large-files \
  --bam-archive-root /mnt/pikachu/msk40ko-bam-archive/MSK-40KO-large-files \
  --globus-poll-seconds 60
```

Use `--delete-local-bam-after-transfer` only after confirming the destination
task succeeds and the local archive policy is acceptable.

The wrapper waits for each BAM Globus task by default. If
`--delete-local-bam-after-transfer` is supplied, deletion is gated on a
`SUCCEEDED` task status and records `deleted_generated_large_files.tsv`.

## Expected Output Layout

```text
<out-root>/samples/<sample>/
  staged_fastqs/
    mRNA/
    PolyIII/
    LARRY/
  pf_multi_config.csv
  RUN_COMMAND.sh
  RUN_MANIFEST.txt
  run/
    Aligned.out.bam
    outs/raw_feature_bc_matrix/
    outs/filtered_feature_bc_matrix/
    outs/raw_velocyto_feature_bc_matrix/
    outs/filtered_velocyto_feature_bc_matrix/
    outs/crispr_analysis/
    cr_assign/
  downstream_genefull_velocyto_cellbender/
    counts.h5ad
    unfiltered_counts.h5ad
    filtered_counts.h5ad
    default_singlet_filtered_counts.h5ad
    final_counts.h5ad
    adaptive_qc_threshold.json
    feature_libraries/
```
