# Runtime alignment exclusion

`--alignExcludeRegions regions.bed` excludes genomic placements while retaining
the same FASTA, annotation and STAR index across runs. The default is `-` (off).
This opt-in feature is independent of the single-cell assay.

BED coordinates are zero-based, half-open and must use the exact index contig
names. Overlapping/adjacent intervals are merged. Comment, `track` and `browser`
lines are accepted; columns after the first three are ignored. Unknown contigs,
empty files, malformed coordinates and out-of-bounds intervals cause an error.
Resolve contig aliases against the reference assembly report before preparing
the BED; verify accession and length, and retain the mapping in run provenance.
If a public BED includes genuinely absent alternate sequences, prepare a validated
subset and retain the omitted entries. No bundled mask is selected automatically.

Excluded seed placements do not count toward seed or anchor multiplicity limits.
When a longest exact seed matches only excluded regions, STAR searches shorter
seeds in the unchanged suffix array to expose permitted alternatives. Candidate
alignments whose aligned blocks touch an exclusion are rejected before updating
best scores and competing with other alignments. Introns and deletions can span
excluded regions; aligned exon bases cannot. Reverse strands and annotated
junction seeds are checked in genomic coordinates. Run logs record the BED path,
record count, merged interval count and number of excluded bases.

This is placement exclusion, not a claim of bit-for-bit equivalence to rebuilding
an index with those bases replaced by Ns. In particular, boundary-crossing
alignments are rejected rather than clipped to the boundary. Ordinary STAR
seeding, window and score limits still apply. The reference sequence and suffix
array must not be edited in memory: their ordering must remain consistent.

The initial implementation supports ordinary short-read genomic alignment.
Flex, chimeric alignment, long-read builds and transformed/supertranscriptome
indexes are rejected when the option is supplied. Flex half-probe caches do not
perform genomic alignment, so this option cannot alter their probe assignments.

Validation: `tests/test_runtime_alignment_exclusion.py --star PATH --out FRESH_DIR`
creates one synthetic index and checks duplicate resolution, rescue of weaker
permitted matches, repeat limits, both strands, intron versus exon overlap,
boundary behavior, an unrelated mask, redundant junction flanks, input errors
and index immutability. Add `--baseline-star PATH`
to compare mask-off alignments against an earlier binary. Serialize this native
test with other benchmark jobs on the host.
