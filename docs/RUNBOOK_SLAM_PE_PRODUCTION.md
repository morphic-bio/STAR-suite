# SLAM PE Production Runbook

Date: 2026-05-12

This runbook records the current paired-end SLAM production surface for the NW
SLAM-Seq panel. The active run was launched from the isolated checkout
`/mnt/pikachu/STAR-suite-slam-pe-prod`, not from the working repository.

Do not rebuild, relink, or restart STAR from a checkout that is running the
production job.

## Inputs and Outputs

- FASTQs: `/mnt/pikachu/SLAM-Seq`
- Output root: `/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z`
- STAR binary: `/mnt/pikachu/STAR-suite-slam-pe-prod/core/legacy/source/STAR`
- Genome index: `/storage/autoindex_110_44/bulk_index`
- SNP mask: `/mnt/pikachu/slam_blank_artifacts_20260201/mask/snps_from_vcf.bed.gz`
- Globus destination endpoint: `61fb8b9a-9b52-456e-928c-30c0fb0140bf`
- Globus destination root: `SLAM-seq-PE-results`

The local production root contains per-sample manifests, command records,
`SlamQuant.out`, `SlamQuant.grandslam.tsv`, `SlamQuant.cB.tsv`, QC JSON/HTML,
and TranscriptVB/tximport quantification outputs. Y/noY BAM and FASTQ
derivatives are uploaded by Globus and then removed locally when the transfer
succeeds. Original input FASTQs are not uploaded.

## Production Parameters

The production runner is:

```bash
bash scripts/run_slam_prod_set.sh \
  --outdir /mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z \
  --threads 16 \
  --globus-dst-endpoint 61fb8b9a-9b52-456e-928c-30c0fb0140bf \
  --globus-dst-root SLAM-seq-PE-results
```

Important STAR settings emitted by the runner:

- `--outSAMtype BAM SortedByCoordinate`
- `--emitNoYBAM yes --keepBAM no`
- `--emitYNoYFastq yes --emitYNoYFastqCompression gz`
- `--quantMode TranscriptVB --quantTranscriptomeSAMoutput BanSingleEnd`
- `--quantVBgcBias 1 --quantVBgenes 1 --quantVBgenesMode Tximport`
- `--slamQuantMode 1 --slamGrandSlamOut 1`
- `--slamCbOut 1 --slamCbFormat star`
- `--slamMinCallableLength 30`
- PE trims: R1 `8/13`, R2 `19/14`

The source of truth for a completed run is
`manifests/run_config.env`, `manifests/samples.tsv`, and
`manifests/commands.tsv` inside the output root.

## Current Snapshot

Read-only snapshot at `2026-05-12T01:12:09Z`:

- Manifested samples: `143` (`manifests/samples.tsv`, excluding header)
- Samples started: `143`
- STAR final logs present: `142`
- `SlamQuant.out`: `142`
- `SlamQuant.cB.tsv`: `142`
- `star_quant.sf`: `142`
- Globus transfers submitted: `142`
- Globus-cleaned sample roots: `140`
- Active sample at snapshot: `NUP_62-24h-1_S235`

These numbers are not final benchmark results because the production run was
still active at the snapshot time.

## Status Checks

Use read-only checks while the production run is active:

```bash
OUT=/mnt/pikachu/SLAM-Seq-PE-results/prod_full_20260511T103908Z
wc -l "$OUT"/manifests/samples.tsv "$OUT"/manifests/commands.tsv
find "$OUT/samples" -mindepth 1 -maxdepth 1 -type d | wc -l
find "$OUT/samples" -name '*Log.final.out' -type f | wc -l
find "$OUT/samples" -name '*SlamQuant.out' -type f | wc -l
find "$OUT/samples" -name '*SlamQuant.cB.tsv' -type f | wc -l
find "$OUT/samples" -name 'star_quant.sf' -type f | wc -l
find "$OUT/samples" -name 'TRANSFER_SUBMITTED.ok' -type f | wc -l
find "$OUT/samples" -name 'TRANSFER_CLEANED.ok' -type f | wc -l
tail -n 40 "$OUT/runner.nohup.log"
```

Do not infer completion from `Log.out`, `Log.progress.out`, or `pgrep` alone.
Use wrapper-written completion/transfer markers plus `*Log.final.out`.

## DESeq2 Handoff

Primary DESeq2 input is the Salmon-compatible TranscriptVB surface exported
through tximport, not STAR internal `ReadsPerGene.out.tab` and not SLAM cB
counts. The pinned analysis image is:

```text
star-suite/slam-deseq2:bioc3.23-deseq2-1.52.0-tximport-1.40.0
```

Pinned versions:

- R/Bioconductor base: `bioconductor/bioconductor_docker:RELEASE_3_23-r-4.6.0`
- Bioconductor: `3.23`
- DESeq2: `1.52.0`
- tximport: `1.40.0`

See `docs/RUNBOOK_SLAM_PE_DESEQ2_COUNT_SURFACES.md` for the design decision on
full TranscriptVB counts versus callable-filtered sensitivity analysis.
