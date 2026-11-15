# PIPseq RNA BI/IMI counting

`--soloPipseqMode Illumina` enables an initial native RNA implementation for
the newer PIP/Illumina single-cell chemistry with a three-base binning index.
It does not process the older v3 12-base UMI chemistry or guide-capture libraries.

The implementation follows the published
[DRAGEN 4.4 transcript-counting description](https://help.dragen.illumina.com/dragen-v4.4/product-guide/dragen-v4.4/dragen-single-cell-pipeline/dragen-scrna-pipseq)
and [4.6 barcode-processing description](https://help.dragen.illumina.com/dragen-v4.6/product-guides/dragen-v4.6/dragen-single-cell-pipeline/dragen-scrna-illumina).
It does not require DRAGEN software or packaged reference files. Read-level
barcode agreement and count concordance against a matched DRAGEN dataset remain
to be measured; this mode does not claim exact reproduction.

## Run

Supply four uncompressed, single-column allowlists in physical tier order,
with sequence lengths 8, 6, 6 and 8. Duplicate, empty, malformed and non-ACGT
lists are rejected. No lists are bundled. The 96 sequences of each tier can be
taken from the reads themselves: after the phase is fixed by the constant
linkers, they are the 96 most frequent sequences at each tier position. With
the default `--soloPipseqMultiMatch Reject` the output does not depend on the
order of the entries in a list.

```bash
core/legacy/source/STAR \
  --runThreadN 8 \
  --genomeDir /path/to/STAR_index \
  --readFilesIn /path/to/R2.fastq.gz /path/to/R1.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix /path/to/fresh_output/ \
  --outSAMtype None \
  --soloPipseqMode Illumina \
  --soloPipseqWhitelist tier1.txt tier2.txt tier3.txt tier4.txt \
  --soloPipseqIPM 0 \
  --soloPipseqMultiMatch Reject \
  --soloFeatures GeneFull \
  --soloStrand Forward \
  --soloCellFilter None
```

The transcript read is first and the barcode read last, following STAR's
ordinary separate-barcode input convention. Multiple lanes of one biological
library may be supplied using the normal comma-separated FASTQ lists. Process
independent libraries in separate invocations: this implementation does not add
a library namespace to otherwise identical cell barcodes.

Exactly one of `GeneFull` (default) or `Gene` is supported. Only reads with
one genomic alignment and one assigned gene count. Normal alignment and
trimming options remain explicit; no vendor aligner equivalence is implied.
Raw output is the default. Standard `--soloCellFilter` options, including
`EmptyDrops_CR`, can filter the estimated-transcript matrix using STAR's
existing caller. Calling policy is distinct from barcode correction.

Prepare the cDNA mate with assay-appropriate trimming before alignment. The
mode does not automatically remove the leading technical base, 5-prime TSO or
3-prime poly-A. Use paired-aware preprocessing so discarded cDNA reads also
remove their barcode mates. The published PIPseq V kit TSO is
`AAGCAGTGGTATCAACGCAGAGTGAATGGG`; see Appendix II of the
[public kit guide](https://crtp.ccr.cancer.gov/wp-content/uploads/fsgc-protocols/PIPseq-V-T10-3-Single-Cell-RNA-Kit-User-Guide-1.pdf).
Record trimming thresholds explicitly; no vendor trimmer equivalence is implied.

The mode owns the barcode and UMI parameters. It bypasses ordinary UMI
correction, homopolymer filtering, cross-gene UMI filtering and UMI counting.
Conflicting explicit options are rejected. The internal `NoDedup` setting
selects a single matrix column; PIP's estimator supplies its values.
The initial implementation requires `--outSAMtype None` and excludes Flex,
two-pass/restart and pfMulti guide integration. No CB/UB sidecar is emitted.

## Barcode and fragment contract

After a phase of 0–3 bases, barcode intervals are zero-based half-open
`[0,8), [11,17), [20,26), [31,39)`; BI is `[39,42)`.
Each phase is scored by the sum of minimum Hamming distances to the four
tier lists. A unique lowest-scoring phase is required. Each tier is then
matched on its own:

1. A sequence that is in the list is that entry. It is never rewritten to
   another entry, even when the list also holds a sequence one base away.
2. Otherwise, a sequence one substitution from exactly one entry is corrected
   to it.
3. Otherwise, a sequence one substitution from several entries is rejected
   (`barcode_ambiguous`). `--soloPipseqMultiMatch FirstListed` takes the first
   listed entry instead; that result depends on the order of the list.
4. A sequence more than one substitution from every entry is a non-match.

Multiple tiers may each have one substitution. Linker sequences are not an
additional gate.

The published description says barcodes within one substitution "will be
corrected to the first matching barcode in the barcode sequence list". It gives
no separate rule for a sequence that is itself listed and is also one base from
another entry, and the order of the reference list is not published. A tier
list can hold entries one base apart, and then the rule matters: taking the
first entry within one substitution regardless of an exact match gives both
members of such a pair the same label. This mode therefore states its own
rule, above. Agreement of barcode assignments with any other pipeline has not
been measured.

Phase ties are rejected because the published tie rule has not been established.
Non-ACGT bases in a candidate barcode tier invalidate that phase. A valid
barcode with a missing/non-ACGT BI retains its barcode metrics but contributes
no RNA count. All 64 A/C/G/T BIs, including homopolymers, are valid and uncorrected.

The corrected cell identifier is the concatenation of the four tier sequences,
without separators or an added suffix. Lookup tables cover the individual
tiers; only observed cell/gene combinations are stored. The 96^4 Cartesian
product is never allocated.

Within each cell/gene, exact `(BI, IMI)` observations collapse. This version
defines IMI as the global genomic coordinate of the strand-aware aligned
5-prime base, with strand included in the key. Soft clipping is not extrapolated.
Reverse reads use the last aligned base; internal splice structure is not
part of the key. These coordinate conventions require validation against
reference outputs. Distinct fragments are held in per-thread sparse sets and
merged before calibration, so duplicate reads in different threads collapse.
Memory scales with observed distinct fragments; disk spilling is not yet
implemented. Use development downsample sizes until memory scaling is measured.

## Transcript estimator and outputs

Let B be occupied BI bins and I distinct BI/IMI pairs for one cell/gene:

```text
B <= 10: count = B
B > 10:  count = max(B, floor(I / IPM))
```

With `--soloPipseqIPM 0`, calibrate using all assigned cell/gene groups
with 5–32 occupied bins, before cell filtering. The coupon-collector expectation
is the stopping-time mean `E[M|reach B] = sum(j=0..B-1, 64/(64-j))`.
IPM is the arithmetic mean of I/E across eligible groups. This is an explicit
interpretation of the documentation, not a claim that its numerical edge cases
have been verified against DRAGEN.

Effective IPM is at least one, enforcing the documented bound B <= count <= I.
Both the raw estimate and this floor are reported. A finite explicit IPM >=1
overrides calibration. With no eligible groups, high-occupancy groups cause an
error requesting an explicit IPM or more reads; if every group has at most ten
bins, IPM is unused and counts can still be emitted. There is no silent
high-occupancy fallback.

Under `Solo.out/GeneFull/` (or `Gene/`):

- `raw/{matrix.mtx,features.tsv,barcodes.tsv}`: estimated transcript counts.
- `filtered/`: the same units after optional standard cell calling.
- `Pipseq.groups.tsv`: barcode, gene ID, assigned reads, occupied bins,
  unique BI/IMI pairs and estimated transcripts.
- `Pipseq.metrics.tsv`: effective counting policies, calibration population,
  raw/effective IPM, calibration size and barcode/read rejection categories.
  `reads_total` equals `barcode_exact + barcode_corrected + barcode_invalid +
  barcode_nonmatch + barcode_phase_tie + barcode_ambiguous`.
- `Summary.csv`: read, barcode and mapping fractions, distinct BI/IMI
  fragments, estimated transcripts, IPM and, with a cell filter, the called-cell
  summary. It is written from the PIP counters; `Barcodes.stats` and
  `Features.stats` of the ordinary UMI path are not produced in this mode.

The raw barcode axis contains cells with at least one counted RNA observation.
Inputs with no accepted RNA observations yield a valid empty matrix. Matrix counts and offsets are checked
against STAR's 32-bit matrix limits. The sample-level calibration and resulting
counts depend on the reads supplied; a downsample is not expected to reproduce
full-library counts.

## Validation

`core/legacy/test/test_pipseq.cpp` (`make -C core/legacy/source test_PipSeq`)
checks decoding, all 64 BIs, exact-first matching in both list orders, single
and multiple one-substitution candidates under both `--soloPipseqMultiMatch`
values, phase ties, malformed input, fragment merging and estimator edges.
`tests/test_pipseq_smoke.py` builds a small synthetic index and runs complete
FASTQ-to-matrix fixtures, including reverse reads, duplicate fragment starts,
fixed/estimated IPM, cell filtering and mode-off ordinary UMI counting.
`tests/test_pipseq_edges.py` adds EmptyDrops, empty matrices, threaded merging
and calibration-failure cases using the smoke fixture's generated index.

The initial clean build and these synthetic checks passed. The existing
complex-barcode whitelist lookup regression also passed. The external CB/UB
integration script skipped because its sample-whitelist fixture was absent.
No real-data or DRAGEN parity result is implied by these checks.

```bash
python3 tests/test_pipseq_smoke.py \
  --star core/legacy/source/STAR \
  --out tests/pipseq_output_example
```

Use a fresh output directory. Clean-build STAR before investigating a regression.
Guide hook matching, guide-to-RNA barcode transformations, guide read counting,
v3 profiles, BAM annotation and reference-output parity are subsequent milestones.
