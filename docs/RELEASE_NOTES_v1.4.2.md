# STAR Suite v1.4.2 Release Notes

Date: 2026-06-29

`v1.4.2` is a patch release for TranscriptVB Salmon parity, STAR-generated
index sidecars, PE bulk benchmark reproducibility, and samtools-style BAM sort
parity. The release artifact version is `v1.4.2`, and `STAR --version` reports
`1.4.2`. The upstream STAR base remains available through
`STAR --upstream-version` (`2.7.11b`), and the genome index compatibility string
remains available through `STAR --genome-compat-version` (`2.7.4a`).

## Fixes

- Improve TranscriptVB parity against Salmon alignment mode by tightening
  library-format detection and reusing a validated detection/read-order path.
- Add STAR index integration for transcriptome GC/effective-length sidecars used
  by TranscriptVB GC-bias correction.
- Link the `compute_gc_bias` helper with OpenMP so release and local builds do
  not fail on OpenMP symbols.
- Align the STAR-suite samtools-style BAM sorter tie-break with samtools
  coordinate order (`tid`, `pos`, reverse-strand bit, then input order).
- Add a synthetic SE/PE BAM-sort regression that compares STAR-suite sorted BAMs
  to `samtools sort --no-PG` at the decompressed BAM byte-stream level.
- Clarify that the PE bulk production benchmark uses internal TranscriptVB only
  on the integrated arm; integrated transcriptome BAM and Salmon QC are opt-in
  parity artifacts.

## Validation

- Clean rebuilt STAR from `core/legacy/source`.
- Passed `tests/run_transcriptvb_chr22_parity_smoke.sh`.
- Passed `tests/test_genome_generate_expected_gc.sh`.
- Passed `tests/test_transcriptvb_gc_effective_length.sh`.
- Passed `tests/run_bamsort_samtools_byte_parity.sh`:
  - SE: 14 records, 2 tied coordinates
  - PE: 24 records, 3 tied coordinates
- Passed `tests/run_ychrom_spill_sorter.sh`.
- Completed the matched PPARG no-Y production benchmark:
  STAR-suite `8:54.52`, external control `16:09.61`, `1.81x` faster.
