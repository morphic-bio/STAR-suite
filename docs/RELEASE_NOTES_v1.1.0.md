# STAR Suite v1.1.0 Release Notes

Date: 2026-05-29

`v1.1.0` is a feature release for native CBQ/BINSEQ input support. The release
artifact version is `1.1.0`, and `STAR --version` reports `1.1.0`. The upstream
STAR base remains available through `STAR --upstream-version` (`2.7.11b`), and
the genome index compatibility string remains available through
`STAR --genome-compat-version` (`2.7.4a`).

## Feature Scope

- Add a native C++ CBQ reader and direct STAR read-buffer adapter behind
  `--readFilesType Binseq PE|SE`.
- Add `cbq_ordered_encoder`, an order-preserving FASTQ/FASTQ.gz-to-CBQ encoder
  for parity-sensitive production tests.
- Cover paired, single-end, manifest, multisample, and batch CBQ input through
  STAR mapper smoke tests.
- Cover STARsolo, OCM composite barcode, Flex, SLAM, and process_features CBQ
  adapter workflows with registered smoke/regression cases.
- Add PE-SLAM CBQ divergence coverage so paired-end SLAM compares FASTQ and
  ordered CBQ at input, alignment, diagnostics, transitions, and pre-NTR
  boundaries.
- Keep Y/noY FASTQ emission gated to FASTQ input for this release; CBQ support
  for Y-removal remains a documented follow-up.
- Keep Chromap ATAC CBQ on the compatibility adapter path; native in-memory
  libchromap CBQ input remains a documented follow-up.

## Validation

- Clean rebuilt STAR and CBQ harness binaries from the release candidate.
- Ran the network-enabled CBQ E2E/module regression suite, including native
  reader, ordered encoder, STAR mapper, STARsolo, process_features adapter,
  Chromap adapter, ARC fixture, and Flex tiny public smoke coverage.
- Ran the JAX scRNAseq02 OCM FASTQ-vs-CBQ local-production smoke.
- Ran SLAM FASTQ-vs-CBQ divergence harnesses with exact pre-NTR parity.
- Ran Fastx input contract, Y-removal FASTQ gate, and assignBarcodes contract
  regression checks.
