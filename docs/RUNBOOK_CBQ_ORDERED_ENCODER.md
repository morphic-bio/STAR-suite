# Runbook: Ordered C++ CBQ Encoder

## Context

The native CBQ reader and STAR/process_features adapters can consume CBQ without
FASTQ materialization. In A375 perturb-seq parity testing, CBQ and FASTQ totals
matched at mapping scale, but filtered cell and matrix counts differed slightly.
The first mismatch appeared before mapping: CBQs produced with the Rust
`bqtools encode` path did not preserve the original FASTQ record order. Its
parallel batch writer can ingest completed batches in completion order rather
than source order.

For STAR-suite production parity tests, CBQ generation must be order-preserving
unless an explicit downstream test is designed to be order-independent.

## Policy

- The STAR-suite ordered encoder is the reference encoder for production CBQ
  parity testing.
- It must write records in the exact order they are read from FASTQ.
- It must be lossless for `A/C/G/T/N` sequence content, headers, qualities, and
  paired read association.
- It must not sort by read name unless the caller explicitly asks for a separate
  order transform. Sorting by read name is deterministic but not equivalent to
  original FASTQ order.
- Existing optimized FASTQ paths remain unchanged. The encoder is only for
  creating ordered CBQ inputs.

## Implementation Surface

Tool target:

```bash
make -C core/legacy/source cbq-ordered-encoder
```

Binary:

```bash
core/legacy/source/cbq_ordered_encoder
```

Basic paired usage:

```bash
core/legacy/source/cbq_ordered_encoder \
  --readFilesIn R1.fastq.gz R2.fastq.gz \
  --outFile lane.cbq
```

Single-end usage:

```bash
core/legacy/source/cbq_ordered_encoder \
  --readFilesIn reads.fastq.gz \
  --outFile reads.cbq
```

Supported inputs:

- Plain FASTQ
- Gzipped FASTQ (`.gz`)
- Single-end or paired-end
- Headers and qualities are always retained
- Paired read names are validated by default after stripping leading `@`,
  whitespace suffixes, and trailing `/1` or `/2`

Format details:

- CBQ file/header/block layout matches the native reader in
  `core/legacy/source/input/CbqInputModule.cpp`.
- Sequence data is packed as 2-bit `A/C/G/T`.
- `N`/ambiguous positions are stored through a sucds-compatible Elias-Fano
  payload; ambiguous non-ACGT bases are restored as `N`.
- A CBQ index is written so external mmap-based tools can inspect the file.
- The write boundary is sequential; no worker can publish a later batch before
  an earlier batch.

## Smoke Test

```bash
tests/run_cbq_ordered_encoder_smoke.sh
```

The smoke test checks:

- paired FASTQ -> CBQ exact ordered contract TSV parity
- paired FASTQ.gz -> CBQ exact ordered contract TSV parity
- single-end FASTQ -> CBQ exact ordered contract TSV parity
- N-position round trip through the native CBQ reader
- optional `bqtools info` index readability and `bqtools decode -T 1`
  byte-identical FASTQ round trip when `bqtools` is available

Artifacts default to:

```bash
/tmp/star_suite_cbq_ordered_encoder_smoke/
```

## Production Validation Plan

1. Build the encoder and native reader harness.
2. Re-encode the A375 GEX and CRISPR FASTQs with `cbq_ordered_encoder`.
3. Verify CBQ reader order against FASTQ on a prefix and on full input contract
   dumps if feasible.
4. Rerun the full A375 perturb-seq E2E using ordered CBQ inputs.
5. Compare against the same-binary FASTQ baseline:
   - mapping totals
   - raw and filtered GeneFull MEX
   - raw and filtered feature MEX
   - filtered barcode set
6. If differences remain, debug downstream determinism rather than input order.

## Notes

- Sorting an already encoded CBQ by read name is a separate utility and should
  not replace this encoder for source-order parity.
- The encoder intentionally does not introduce a new reader contract or alter
  the production FASTQ path.
- Compression level defaults to `0`, matching current CBQ defaults. Use
  `--compressionLevel N` only when testing storage/performance tradeoffs.
