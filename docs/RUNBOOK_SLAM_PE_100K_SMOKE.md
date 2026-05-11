# SLAM PE 100K SE/PE Smoke

This runbook records the reproducible 100K smoke used before scaling the
SLAM-seq paired-end production run. It compares the same downsampled R1 data as
single-end against R1/R2 paired-end, using the fixed noSU-derived trim windows.

## Script

Use the tracked harness:

```bash
bash scripts/run_slam_100k_se_pe_smoke.sh \
  --sample ARID1A-no4su_S50 \
  --sample ARID1A-6h-1_S43 \
  --threads 16
```

Default outputs go under `/tmp/slam_pe_100k_smoke_<timestamp>/`. The script
writes:

- `manifests/run_config.env`
- `manifests/samples.tsv`
- `manifests/commands.tsv`
- `runs/<sample>/se/RUN_COMMAND.sh`
- `runs/<sample>/pe/RUN_COMMAND.sh`
- `compare/*.log`
- `compare/smoke_metrics.tsv`

Do not launch this while a production or benchmark STAR-SLAM run is active on
the same host.

## Pinned Smoke Parameters

- Downsample: 100000 paired records per sample.
- SE input: R1 only.
- PE input: R1 and R2.
- SE trim: `--slamCompatTrim5p 8 --slamCompatTrim3p 12`.
- PE trim: R1 `8/13`, R2 `19/14`.
- Callable gate: `--slamMinCallableLength 30`.
- TranscriptVB: `--quantTranscriptomeSAMoutput BanSingleEnd`,
  `--quantVBgenesMode Tximport`.
- Count outputs: GrandSLAM plus cB STAR format.
- SNP mask: `/mnt/pikachu/slam_blank_artifacts_20260201/mask/snps_from_vcf.bed.gz`.

## Latest Smoke Results

Latest artifact used for these numbers:
`/tmp/slam_min_callable_smoke_20260511T071911Z`.

Comparison logs were regenerated from that artifact with:

```bash
python3 scripts/compare_slam_summary.py SE/SlamQuant.out PE/SlamQuant.out \
  --allow-gene-set-mismatch \
  --min-read-count 100 \
  --max-abs-delta 1.0
```

For treatment samples, add `--min-ntr-pearson 0.90`. Do not apply the treatment
NTR correlation gate to noSU/no4SU controls; there are too few true conversions
for that correlation to be biologically meaningful.

| Sample | Shared SLAM genes | Genes with ReadCount >=100 in both | NTR Pearson | NTR Spearman | SLAM ReadCount Pearson | SLAM Coverage Pearson | Tximport gene NumReads Pearson |
|--------|------------------:|-----------------------------------:|------------:|-------------:|-----------------------:|----------------------:|-------------------------------:|
| ARID1A_6h_1 | 7250 | 36 | 0.972806 | 0.846589 | 0.902310 | 0.873196 | 0.932203 |
| noSU | 7872 | 37 | 0.238890 | 0.137980 | 0.921661 | 0.900751 | 0.938455 |

STAR mapping summary:

| Sample | Mode | Unique mapped | Too short |
|--------|------|--------------:|----------:|
| ARID1A_6h_1 | SE | 75.58% | 4.76% |
| ARID1A_6h_1 | PE | 62.91% | 29.93% |
| noSU | SE | 76.89% | 3.00% |
| noSU | PE | 66.62% | 24.90% |

Count totals from the same smoke:

| Sample | Surface | SE total | PE total |
|--------|---------|---------:|---------:|
| ARID1A_6h_1 | SLAM ReadCount | 61079.988136 | 43392.509368 |
| ARID1A_6h_1 | SLAM Coverage | 544107.767236 | 585308.307178 |
| ARID1A_6h_1 | Tximport gene NumReads | 75262.001 | 50944.999 |
| noSU | SLAM ReadCount | 63011.014967 | 46469.003879 |
| noSU | SLAM Coverage | 578542.136163 | 644684.797250 |
| noSU | Tximport gene NumReads | 77725.002 | 55004.999 |

Interpretation:

- The treatment smoke passed the high-confidence-gene NTR gate
  (`0.972806 >= 0.90`).
- The noSU NTR Pearson is low and should not be used as a pass/fail gate. Its
  count and coverage correlations are high, but the sparse/background
  conversion surface is dominated by noise.
- PE still has fewer SLAM reads and TranscriptVB counts than SE, while PE SLAM
  coverage is higher. This is consistent with the R2 quality/length effect seen
  in the mapping summaries rather than an EM/VB-only failure.
