# STAR Suite v1.12.0 Release Notes

Date: 2026-10-14

STAR Suite 1.12.0 restricts feature-read file discovery to FASTQ filenames.
Previously, a filename containing the requested read pattern was eligible
even when it ended in another suffix. A receipt such as
`sample_R1_001.fastq.gz.receipt.json` beside the reads could therefore be
opened as read data and cause a FASTQ header error.

Discovery now accepts `.fastq` and `.fq`, optionally followed by `.gz`,
`.bgz` or `.bgzf`, ignoring case. These suffixes represent the plain, gzip
and BGZF inputs handled by the existing zlib reader. Selection uses the
directory-entry name: a FASTQ-named symbolic link remains eligible, while
a receipt-named link is excluded. All seven calls to this module's finder
select FASTQ inputs, including the barcode pattern and both default and
explicit split-read triplet patterns.

The change only selects files. Reading, file order among eligible inputs,
threading, pairing and triplet validation are unchanged. On the existing
module fixture, directories containing only FASTQ inputs produce the same
output bytes as v1.11.0. The CAT-ATAC guide check compares the original
FASTQ-only directory with the same reads and adjacent receipts.

Regression coverage includes each supported suffix, receipt-only
directories, FASTQ files beside receipts, symbolic links, unsupported
suffixes, barcode discovery and default/explicit triplet patterns. Two
test sources already required by the module Makefile are also supplied:
filtered-barcode namespace validation and isolation of independent API
contexts. No additional processing implementation is included.

`STAR --version` reports `1.12.0`. Debian source packaging uses `1.12.0-1`;
Ubuntu packages use `1.12.0-1~ubuntu22.04.1` and `1.12.0-1~ubuntu24.04.1`.
Upstream STAR remains `2.7.11b`, genome-index compatibility remains
`2.7.4a`, and legacy compatibility remains `2.7.1a`. Existing indexes do
not need rebuilding.
