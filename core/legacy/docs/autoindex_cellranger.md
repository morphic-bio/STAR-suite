# AutoIndex + CellRanger-Style Reference Notes

This page documents STAR-suite AutoIndex behavior for CellRanger-style reference
formatting and the current parity scope.

## Scope and Parity Contract

For CellRanger-style indexing, STAR-suite guarantees parity at the level of the
formatted reference inputs:

- `${genomeDir}/cellranger_ref/genome.fa`
- `${genomeDir}/cellranger_ref/genes.gtf`

For GRCh38 `2024-A`, when using the same downloaded Ensembl/GENCODE source files
and updated filtering mode, these files are byte-identical to 10x downloaded
reference inputs.

Parity is **not** guaranteed for generated STAR index binaries (`Genome`, `SA`,
`SAindex`, etc.), because those depend on STAR index engine/version and index
build parameters.

## Release Mapping and Modes

AutoIndex release mapping:

- `--cellrangerRefRelease 2024-A` -> Ensembl 110 FASTA + GENCODE v44 GTF
- `--cellrangerRefRelease 2020-A` -> Ensembl 98 FASTA + GENCODE v32 GTF

GTF/FASTA filtering mode:

- `--cellrangerLegacyGtfFilter Auto` (default)
  - `2020-A` -> legacy mode
  - `2024-A` -> updated mode
- `--cellrangerLegacyGtfFilter Yes` -> force legacy mode
- `--cellrangerLegacyGtfFilter No` -> force updated mode

## Legacy vs Updated Gene-Set Behavior

For GRCh38 `2024-A`:

- Updated mode includes 47 genes that legacy mode excluded.
  - Cause: support for `protein_coding_LoF` transcripts.
- Updated mode excludes 33 chrY PAR genes that legacy mode kept.
  - Cause: chrY PAR terminal masking/filtering to match current CR reference behavior.

Net result: updated mode aligns with downloaded 10x `2024-A` FASTA/GTF inputs,
while legacy mode preserves historical behavior.

## Cell Ranger Provenance (A375 Benchmark Reference)

From the downloaded `refdata-gex-GRCh38-2024-A/star/genomeParameters.txt`:

- Command header points to `cellranger-8.0.0/lib/bin/STAR`
- `versionGenome` is `2.7.1a`

This is the reference-input target used for AutoIndex parity checks.

Some benchmark bundles also contain an auxiliary `star_2.7.11b/` index directory
for local compatibility testing. That directory is not the original downloaded
Cell Ranger `star/` index and can report a different `versionGenome`.

## What Can Still Differ After Input Parity

After `genome.fa`/`genes.gtf` parity is achieved, remaining differences come from
index engine/build behavior, e.g.:

- STAR `versionGenome` and indexer implementation
- Index build knobs (`genomeSAsparseD`, `sjdbOverhang`, RAM-dependent chunking)

These differences can change index binaries and downstream mapping/count outputs
even when reference inputs are identical.

## Example (2024-A Updated Mode)

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/index \
  --autoIndex Yes \
  --cellrangerStyleIndex Yes \
  --cellrangerRefRelease 2024-A \
  --cellrangerLegacyGtfFilter No \
  --autoCksumUpdate Yes \
  --sjdbOverhang 100 \
  --runThreadN 16
```
