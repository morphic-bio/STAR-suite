# Public Bulk PE Fixture and Regression Tests

## Purpose

This document pins a public human bulk RNA-seq PE fixture for the PE bulk
regression checks:

- Y/noY FASTQ emission parity
- TranscriptVB PE auto-detect health
- optional Salmon parity on the same public fixture

The fixture is public GEO/SRA data and can be downloaded locally without using
restricted Morphic/KOLF-related datasets.

## Pinned Public Fixture

Default fixture used by `scripts/download_public_bulk_fixture.sh`:

- GEO series: `GSE88509`
- GEO sample: `GSM2344101`
- SRA run: `SRR4422207`
- organism: `Homo sapiens`
- sample: male esophagus muscularis mucosa
- layout: paired-end

Default local output:

```bash
/tmp/public_bulk_fixture_SRR4422207/
```

Default download behavior:

- `fastq-dump --split-files -X 500000`
- deterministic first `500000` spots, not a random resample

Default STAR mate order for the public PE smoke:

- `PUBLIC_BULK_STAR_MATE_ORDER=R2R1`
- this matches the current TranscriptVB PE documentation for
  mate-order-sensitive libraries

Fixture-specific expectations for the pinned smoke:

- expected auto-detected format: `ISF`
- Salmon parity thresholds:
  - Spearman all `>= 0.95`
  - Spearman expressed `>= 0.99`
  - Pearson expressed `>= 0.99`
- `dropped_incompat` is not forced to zero for this stranded fixture; parity
  with Salmon is the main acceptance criterion

## Reference Paths

Current default human bulk reference used by the smoke:

- STAR index:
  `/storage/autoindex_110_44/bulk_index`
- transcriptome FASTA:
  `/storage/autoindex_110_44/bulk_index/transcriptome.fa`

These paths are wired into `tests/external_fixtures_env.sh`.

## Scripts

Prepare the public fixture:

```bash
scripts/download_public_bulk_fixture.sh
```

Replay archived bad outputs and verify they fail:

```bash
tests/run_pe_bulk_regression_replay.sh
```

Run the public PE smoke:

```bash
tests/run_public_bulk_pe_smoke.sh
```

## Expected Behavior

`tests/run_pe_bulk_regression_replay.sh` should:

- fail on the old broken Y/noY FASTQ emission artifact
- fail on the old broken TranscriptVB auto-detect artifact
- pass on the fixed Y/noY FASTQ emission artifact
- pass on the fixed TranscriptVB auto-detect artifact

`tests/run_public_bulk_pe_smoke.sh` should:

- run STAR integrated bulk PE with:
  - `--trimCutadapt Yes`
  - `--emitNoYBAM yes`
  - `--emitYNoYFastq yes`
  - `--quantMode TranscriptomeSAM TranscriptVB`
  - `--quantVBLibType A`
- validate BAM-to-FASTQ parity for Y/noY emission
- run Salmon on `Aligned.toTranscriptome.out.bam`
- validate TranscriptVB auto-detect log health
- validate Salmon parity thresholds

## Validated Local Run

Validated on `2026-03-10` with output:

- `tests/public_bulk_pe_smoke_output_20260310_005835/`

Observed results:

- Y/noY BAM-to-FASTQ parity: `PASS`
- detected format: `ISF`
- `dropped_missing_mate_fields: 0`
- `dropped_unknown_obs_fmt: 0`
- Salmon parity:
  - Spearman all: `0.9627`
  - Spearman expressed: `0.9926`
  - Pearson expressed: `0.99998`
  - Salmon total reads: `233117`
  - STAR total reads: `233090`

## Restrictions

This public fixture is safe for distribution.

Do not mix it with restricted local Y-derived artifacts from non-public PE
datasets. Public fixture outputs and restricted local outputs should remain
separate in documentation and file locations.
