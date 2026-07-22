# STAR Suite v1.0.3 Release Notes

Date: 2026-05-29

`v1.0.3` is a hotfix release for non-Flex Solo direct bridge determinism. The
release artifact version is `1.0.3`, and `STAR --version` reports `1.0.3`. The
upstream STAR base remains available through `STAR --upstream-version`
(`2.7.11b`), and the genome index compatibility string remains available
through `STAR --genome-compat-version` (`2.7.4a`).

## Hotfix Scope

- Canonicalize worker-local direct bridge hash entries before non-Flex Solo
  UMI collapse so duplicate `(cell barcode, gene, UMI)` observations have a
  deterministic single count surface before `1MM_CR`.
- Aggregate ambiguous cell-barcode likelihood evidence deterministically before
  Bayesian resolution instead of depending on a schedule-sensitive
  representative read context.
- Use the deterministic bridge-resolved assignment for ambiguous-read
  accounting statistics, eliminating residual repeat-run drift in
  `Features.stats`, `Summary.csv`, and `Barcodes.stats`.
- Add a production-scale A375 repeat-run bridge determinism harness that
  compares bridge stages, per-cell-barcode digests, raw/filtered GeneFull MEX,
  and Solo stats.

## Validation

- Clean rebuilt STAR from the hotfix branch.
- Re-ran the full A375 non-Flex Solo bridge determinism harness twice at
  32 threads. Bridge stages, per-cell-barcode bridge digests, raw and filtered
  GeneFull MEX matrices, `Features.stats`, `Summary.csv`, and `Barcodes.stats`
  matched exactly.
- Ran a non-instrumented STAR-only A375 timing pass with the corrected bridge
  path: 242 seconds wall time with 1187 filtered cells, matching the prior
  241-second STAR-only production benchmark within noise.
