# SPLiT-seq RNA baseline

`--soloSplitseqMode Fixed` decodes fixed split-pool barcode segments and a real
UMI, then uses ordinary STAR alignment and UMI counting. It defaults to
`GeneFull` and the existing `EmptyDrops_CR` caller without assay-specific tuning.
The ordinary `soloUMIdedup` default (`1MM_All`) remains in effect; explicit
standard UMI deduplication/filtering options remain available.

## Invocation

Supply a layout and allowlists for the actual library. Coordinates below are
zero-based `offset_length`, in **round 1, round 2, round 3 cell-ID order**:

```bash
STAR --genomeDir /path/to/index \
  --readFilesIn RNA_R1.fastq.gz barcode_R2.fastq.gz --readFilesCommand zcat \
  --soloSplitseqMode Fixed \
  --soloSplitseqCBposition 78_8 48_8 10_8 \
  --soloSplitseqUMIposition 0_10 \
  --soloSplitseqWhitelist round1.tsv round2.tsv round3.tsv \
  --outSAMtype None --runThreadN 8 --outFileNamePrefix /path/to/output/
```

This example layout is declared by the deposited
[IGVF RNA layout IGVFFI2031CEME](https://igvf-public.s3.amazonaws.com/2026/08/29/7503b73d-1543-4c6d-84fd-1da5fe11602d/IGVFFI2031CEME.yaml.gz).
It is not a universal Evercode or SPLiT-seq version preset. Confirm the read and
allowlist orientations and any well replacement mapping for the target data.
There is no seqspec runtime dependency. Linkers are skipped at the supplied
coordinates; no linker alignment or variable-phase search is performed.

## Barcode contract

- Supply 2–4 nonoverlapping CB segments and one nonoverlapping UMI region.
  Total CB length is at most 31 bases; UMI length is 1–16 bases.
- Allowlist files are plain text, one observed A/C/G/T sequence per line, with
  an optional second column giving its same-length canonical replacement.
  Files have no header. Empty lines are ignored; duplicate observed entries,
  invalid bases, wrong lengths, extra columns and empty files are errors.
- Each round uses an exact match first, then a unique Hamming-distance-one
  match. Several rounds can each be corrected in one read. A segment one
  substitution from several observed entries is rejected as ambiguous, unless
  all of them map to the same canonical sequence. Non-ACGT barcode bases are
  rejected; no abundance-based tie resolution or requirement for a whole-tuple
  exact supporting read.
- Replacement is applied **once after correction**, never chained. Multiple
  observed entries may map to one canonical sequence, allowing primer-barcode
  aliases or explicit well replacement. The caller must supply the correct
  mapping; no sample/well mapping is guessed or downloaded.
- Cell IDs concatenate canonical segments in the supplied order. Only observed
  tuples enter counting; the Cartesian barcode product is never allocated.
  IDs are deterministic and do not depend on allowlist ordering or threads.
- Process one biological sublibrary per invocation; lanes from that same
  sublibrary may share a run. Keep the experiment and sublibrary identity in
  downstream metadata when joining matrices. Distinct sublibraries must not
  be pooled under these sequence-only cell IDs.

The UMI follows the normal STAR path, including rejection of non-ACGT UMIs and
homopolymers, followed by the configured UMI collapse. Split barcode correction
is separate from UMI correction. The mode disables the simple-CB correction and
inline hash paths that assume a global whitelist. Conflicting assay-owned CLI
settings are errors.

## Outputs and scope

The standard `Solo.out/GeneFull/raw/` and `filtered/` MEX directories are emitted;
an enabled caller with zero called cells still produces a valid empty filtered
matrix. `--soloCellFilter None` is available for upstream-only checks.
`Solo.out/Splitseq.metrics.tsv` reports barcode acceptance, rejection reasons,
UMI rejection and per-round correction counts. Every read is counted in
exactly one of `barcode_exact`, `barcode_corrected`, `read_too_short`,
`barcode_non_acgt`, `barcode_no_match`, `barcode_ambiguous`, `umi_non_acgt` and
`umi_homopolymer`, so these sum to `reads_total`. `barcode_exact` and
`barcode_corrected` are the reads passed to alignment and counting; the
`roundN_corrected` rows count those reads again per corrected round. The
standard feature statistics describe the subsequent counting stage.

This first implementation supports RNA with one `Gene` or `GeneFull` feature.
Guide capture, BAM/tag output, automatic kit detection, sample demultiplexing,
Flex, two-pass and restart modes are not integrated. This is a STAR baseline,
not a reproduction of Parse's pipeline or cell-selection rule. Real-data
qualification and cell-calling optimization remain separate steps.

## Synthetic validation

Decoder checks: `core/legacy/test/test_splitseq.cpp`
(`make -C core/legacy/source test_SplitSeq`).
FASTQ-to-MEX checks: `tests/test_splitseq_smoke.py` (small generated genome;
no external data). These cover correction in multiple rounds, ambiguous and
invalid barcodes, replacement aliases, 24-base cell IDs, ordinary UMI collapse,
default EmptyDrops, empty results and thread consistency. They do not establish
agreement with Parse counts or biological cell calls.
