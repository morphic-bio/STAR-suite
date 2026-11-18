# STAR Suite v1.13.1

## Shared GEX multimapper rescue

The resolver could skip a transcript whose exon covered the read when the
preceding transcripts ended earlier. A read mapping to that exon and to an
intron at a second locus could consequently lose its gene count. The corrected
search checks the current transcript before using the index to rule out earlier
transcripts. It works with existing genome indexes.

Gene counts can change in assays that use `--soloCrMultimapRescue yes`. To
reproduce a previous rescue-enabled analysis exactly, retain its STAR Suite
version along with its index and settings.

## Experimental split-pool RNA

`--soloSplitseqMode Fixed` now defaults to `GeneFull`,
`--soloCrMultimapRescue yes` and `--soloMultiMappers Unique`. A single
eligible exonic alignment can rescue a genomic multimapper. In the default
gene-body mode, a single eligible intronic alignment can also win when there
is no eligible exonic alignment. On one alignment where a nested exon also
overlaps a host intron, transcript-concordant exons take priority. Reads
resolved to one gene enter UMI counting.

For the v1.13.0 split-pool counting policy, set `--soloCrMultimapRescue no`.
The split-pool Launchpad recipe exposes the same option. [Split-pool RNA guide](SPLITSEQ_RNA.md).

## Regression coverage

A synthetic FASTQ-to-matrix test checks an exon against an intronic alternative,
intronic fallback, a nested exon, and ambiguity at shared exons and competing
exonic loci. The test also checks that ordinary GEX and split-pool GEX produce
the same counts from the same alignments and canonical cell barcodes. The
historical 100K split-pool regression uses the v1.13.0 setting and its original
expected matrix.
