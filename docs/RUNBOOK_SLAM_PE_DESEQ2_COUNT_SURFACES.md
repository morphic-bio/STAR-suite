# SLAM PE DESeq2 Count Surfaces

Date: 2026-05-11

This note records the count-surface decision for the SLAM paired-end production
run.

## Primary DESeq2 Surface

Use the full TranscriptVB/Salmon-style expression surface for primary DESeq2:

```text
star_quant.sf -> containerized R tximport -> DESeqDataSetFromTximport
```

Do not use STAR `ReadsPerGene.out.tab`, `SlamQuant.out`, or count-binomial `n`
values for the primary DESeq2 analysis.

The primary DESeq2 estimand is expression abundance from the best available
transcript-quantifiable evidence. It should not be restricted to reads that are
also informative for T->C conversion, because that changes the estimand to
expression among SLAM-callable fragments.

## SLAM Callable Gate

`--slamMinCallableLength` is a SLAM evidence gate, not an expression gate. It is
applied after alignment, SLAM trim windows, and paired-end overlap consensus.
Reads or read pairs below the threshold do not contribute to SLAM transition,
NTR, GrandSLAM-style, or count-binomial statistics.

The gate remains justified for SLAM-specific modeling because very short
post-trim callable fragments can add per-read count mass while contributing
little T coverage or conversion information. The gate provides a minimum sample
size for the per-read T->C evidence.

## Why Primary DESeq2 Should Stay Unfiltered

Filtering DESeq2 to only SLAM-callable reads sounds clean, but it risks adding
a different bias:

- genes and transcripts with lower T content would lose more reads;
- short exons, short isoforms, and fragment geometries near trim windows could
  be penalized;
- tximport length correction does not model this extra read-level callable
  filter;
- condition-specific callable fractions would become an additional nuisance
  variable in the expression analysis.

For primary DESeq2, total expression evidence is the correct default. The SLAM
count surface and the expression count surface should be interpreted as related
but distinct measurements.

Important caveat: callable-filtering both the SLAM surface and the DESeq2
surface does not remove the core assumption. If short-read or callable-length
bias is strong enough to distort interpretation, it can affect both T->C
assignment and expression detection. Using the same filtered read population for
both analyses simply assumes that the filtered population is representative,
while potentially sacrificing substantial expression signal. Keeping the full
DESeq2 surface and using the callable-filtered surface only as a sensitivity
analysis makes that assumption testable instead of hiding it in the primary
analysis.

## Callable-Filtered DESeq2 Sensitivity

Callable-filtered DESeq2 is possible as a sensitivity analysis, but it requires
creating a separate matched quantification surface. The current production
`star_quant.sf` files are not filtered by `--slamMinCallableLength`.

Do not try to recover callable-filtered DESeq2 by post-hoc filtering or scaling
the gene-level tximport matrix. The callable gate is read/pair-level and occurs
before transcript/gene assignment weights are finalized.

Acceptable implementation options:

1. Add a STAR/TranscriptVB option that applies the SLAM callable gate before
   adding a read or read pair to TranscriptVB equivalence classes. This is the
   preferred exact integrated approach.
2. Emit a read-name pass list from the SLAM callable gate, then re-run
   transcript quantification on filtered paired FASTQs or BAMs. This is more
   IO-heavy and must keep both mates synchronized.

For the current full panel, use the existing run as the primary DESeq2 input.
If callable-filtered DESeq2 is needed, run a pilot first and then rerun the
full panel or a quantification-only equivalent after one of the explicit
filtered count-surface paths exists.

Recommended output separation:

```text
de/deseq2_primary_full_txvb/
de/deseq2_sensitivity_slam_callable_min30/
```

Required sensitivity QC:

- per-sample callable pass fraction;
- library-size change after filtering;
- gene-level count correlation vs primary DESeq2 input;
- association of count loss with gene length, transcript length, T content, and
  target/timepoint;
- DE contrast concordance for top genes and global log2 fold changes.
