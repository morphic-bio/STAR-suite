# Flex Sample-Tag Resolution Policy

This document fixes the production sample-tag policy used by the fused Flex
FASTQ, CBQ, and BAM paths. The same policy must be used by every input path.

## Lookup policy

STAR extracts exactly eight bases beginning at `--soloSampleProbeOffset`
(default 68). It does not infer a different offset and does not search adjacent
positions.

Sample-tag resolution is a two-tier lookup:

1. H0 is the supplied sample-probe table. An unambiguous exact H0 entry is
   authoritative and returns immediately.
2. Only after an H0 miss, STAR queries a separate H1 cache. The cache is built
   by enumerating all 24 single-base substitutions of every H0 entry. An H1 key
   is assigned only if every table entry that generated it belongs to the same
   sample. A key generated for more than one sample is marked ambiguous during
   construction and is rejected at runtime.

An eight-base window containing an uncalled base cannot address either cache
and is rejected. `--soloSampleTagMismatch 0` or `--soloSampleStrictMatch yes`
disables H1 and leaves H0-only lookup. The default is
`--soloSampleTagMismatch 1`.

`--soloSampleSearchNearby` is retained as a compatibility parameter, but its
only supported value is `no`. Supplying `yes` is a parameter error.

## Rationale and validation boundary

The fixed offset and supplied table define the assay contract. Neighboring
offset searches can assign unrelated sequence and were observed to recover
reads for which the reference pipeline had not emitted a corrected sample tag.
The constructed H1 cache makes ambiguity a property of the complete supplied
table rather than of runtime search order.

This policy does not claim to reproduce an undocumented reference
implementation detail. It is a systematic rule derived from the assay table.
Any future offset or tolerance relaxation requires validation against an
independent public gold standard; a dataset used to tune a rule must not also
serve as its sole evaluation set.
