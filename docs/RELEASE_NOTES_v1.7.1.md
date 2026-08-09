# STAR Suite v1.7.1 Release Notes

Date: 2026-08-09

`v1.7.1` is a bug-fix release for ordinary Fastx STARsolo throughput and
memory use. It restores the established STAR input hot path and avoids an
unused barcode-correction precompute in plain scRNA-seq runs.

The release artifact version is `v1.7.1`, Debian packages use `1.7.1-1`, and
`STAR --version` reports `1.7.1`. The upstream STAR base remains `2.7.11b`;
genome-index compatibility remains `2.7.4a`. Existing indexes do not need to
be rebuilt.

## STARsolo Fastx scaling

- Route ordinary `--readFilesType Fastx` runs through STAR's direct chunk
  reader again. The replaced record-at-a-time adapter parsed, allocated, and
  reformatted records while STAR held its global input mutex, which serialized
  input and starved mapping threads at higher `--runThreadN` values.
- Keep the generic Fastx module as the independently tested input-contract
  implementation rather than installing it on the mapping hot path.
- Preserve plain FASTQ, internal gzip, explicit `zcat`, comma-separated lanes,
  manifest, single-end, and Y/no-Y input behavior.
- No dynamic-thread flag is required for the fix. In particular,
  `--dynamicThreadInterface` and `--dynamicThreadConstMapPermits` did not
  remove the serialized reader bottleneck in affected binaries.

## Plain STARsolo startup and memory

- Construct the expanded `CbCorrector` only for correction modes that request
  it. Ordinary STARsolo continues to use its established whitelist hash and
  `matchCBtoWL()` behavior.
- With the 10x 3M whitelist, the removed precompute generated up to roughly
  177 million one-mismatch probes despite not contributing to the plain-Solo
  matrix path. On the validation host, avoiding it removed 46–48 seconds of
  startup and approximately 8.5 decimal GB of peak RSS.
- This change does not add persistent barcode lookup tables or alter the
  separate spatial/Visium decoder. Those are distinct future optimizations.

## Validation

- Clean-build synthetic tests passed for plain FASTQ, internal gzip, explicit
  decompression, multiple lanes, manifest input, single-end input, Y/no-Y
  handling, GTF-backed genome generation, sorted BAM, and STARsolo matrices.
- Plain STARsolo raw and filtered GeneFull MEX outputs were byte-identical
  before and after the fix and across 8- and 32-thread fixed runs.
- In upstream-compatible mode, fixed STAR Suite raw and filtered GeneFull MEX
  outputs were byte-identical to official upstream STAR 2.7.11b.
- On a serialized 10-million-read local comparison, mapping improved from
  24 to 20 seconds at 8 versus 32 threads in the affected Suite binary to
  22 to 12 seconds after the fix. Official upstream required 21 and 11
  seconds. Fixed peak RSS was within 0.4 GB of upstream at both thread counts.

The full diagnosis, commands, artifact locations, and timing table are in
`docs/SOLO_FASTX_SCALING_FIX_20260809.md`.
