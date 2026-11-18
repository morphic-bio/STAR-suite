# SPLiT-seq RNA baseline

`--soloSplitseqMode Fixed` decodes fixed split-pool barcode segments and a real
UMI, then uses ordinary STAR alignment and UMI counting. It defaults to
`GeneFull` and the existing `EmptyDrops_CR` caller without assay-specific tuning.
The ordinary `soloUMIdedup` default (`1MM_All`) remains in effect; explicit
standard UMI deduplication/filtering options remain available.

Fixed mode uses the shared STAR Suite GEX resolver (`--soloCrMultimapRescue yes`).
With its default evidence policy, one countable exonic alignment wins across
genomic loci. With an intron-inclusive feature, one countable intronic alignment
can win when there are no countable exonic alignments. On a single alignment
overlapping several gene bodies, `GeneFull` gives transcript-concordant exons
priority. The default `--soloMultiMappers Unique` counts reads resolved to one
gene before UMI collapse. For the v1.13.0 split-pool behavior, specify
`--soloCrMultimapRescue no`.

## Invocation

Supply a layout and allowlists for the actual library. Coordinates below are
zero-based `offset_length`, in **round 1, round 2, round 3 cell-ID order**:

```bash
STAR --genomeDir /path/to/index \
  --readFilesIn RNA_R1.fastq.gz barcode_R2.fastq.gz --readFilesCommand zcat \
  --soloSplitseqMode Fixed \
  --soloSplitseqCBposition 78_8 48_8 10_8 \
  --soloSplitseqUMIposition 0_10 \
  --soloSplitseqWhitelist round1.tsv round2.tsv round3.tsv \
  --outSAMtype None --runThreadN 8 --outFileNamePrefix /path/to/output/
```

This example layout is declared by the deposited
[IGVF RNA layout IGVFFI2031CEME](https://igvf-public.s3.amazonaws.com/2026/08/29/7503b73d-1543-4c6d-84fd-1da5fe11602d/IGVFFI2031CEME.yaml.gz).
It is not a universal Evercode or SPLiT-seq version preset. Confirm the read and
allowlist orientations and any well replacement mapping for the target data.
There is no seqspec runtime dependency. Linkers are skipped at the supplied
coordinates; no linker alignment or variable-phase search is performed.

### v3 chemistry (short linkers)

The barcode read of the v3 chemistry holds the UMI at bases 0-9, round 3 at
10-17, the linker `ATGAGGGGTCAG`, round 2 at 30-37, the linker
`TCCAACCACCTC` and round 1 at 50-57 (zero-based). The kallisto-bustools
documentation gives the same layout for its `SPLIT-SEQ` v3 technology string
(`1,10,18,1,30,38,1,50,58:1,0,10:0,0,0`).

```bash
L=share/star-suite/barcode_lists/splitseq_v3/48_wells   # 96_wells for a 96-well round 1
STAR --genomeDir /path/to/index \
  --readFilesIn RNA_R1.fastq.gz barcode_R2.fastq.gz --readFilesCommand zcat \
  --soloSplitseqMode Fixed \
  --soloSplitseqCBposition 50_8 30_8 10_8 \
  --soloSplitseqUMIposition 0_10 \
  --soloSplitseqWhitelist $L/round1.tsv $L/round2.tsv $L/round3.tsv \
  --soloStrand Forward --outSAMtype None --runThreadN 16 --outFileNamePrefix /path/to/output/
```

The lists in `share/star-suite/barcode_lists/splitseq_v3/` are written by
`scripts/make_splitseq_lists.py` from the on-list and the replacement list of
the kallisto-bustools documentation (Pachter lab, BSD-2-Clause);
`PROVENANCE.md` there records the source commit, the checksums and what was
checked on data. `round1.tsv` counts the random-hexamer barcode of a well as
the oligo-dT barcode of the same well, so a cell is one barcode. Use the list
that matches the round-1 wells of the library: a list with unused wells lets a
sequencing error be corrected into one of them. The lists have been checked
on reads for round-1 wells 1 to 48 (one library, below); wells 49 to 96 of
`96_wells` come from the source files alone.

Well labels for the cell barcodes of a run:

```bash
python3 scripts/make_splitseq_lists.py \
  --onlist share/star-suite/barcode_lists/splitseq_v3/source/splitseqv3_barcodes.txt \
  --replace share/star-suite/barcode_lists/splitseq_v3/source/splitseqv3_replace.txt \
  --wells 48 --label-barcodes /path/to/output/Solo.out/GeneFull/raw/barcodes.tsv --out-dir labels
```

`labels/barcode_wells.tsv` gives the three wells of every barcode and the
label `%02d_%02d_%02d` of their indices (row-major: A1 = 1, B1 = 13). Sample
assignment from round-1 wells is the experiment's own table and is not
applied here. An output directory can be used again. A run that fails while
it makes or checks its files leaves the directory as it was. A run that gets
further replaces the tool's own files: it removes `summary.json` first and
writes it last. If that step fails, the directory can hold files of both runs
and has no `summary.json`; use a directory only when `summary.json` is there.

## Barcode contract

- Supply 2–4 nonoverlapping CB segments and one nonoverlapping UMI region.
  Total CB length is at most 31 bases; UMI length is 1–16 bases.
- Allowlist files are plain text, one observed A/C/G/T sequence per line, with
  an optional second column giving its same-length canonical replacement.
  Files have no header. Empty lines are ignored; duplicate observed entries,
  invalid bases, wrong lengths, extra columns and empty files are errors.
- Each round uses an exact match first, then a unique Hamming-distance-one
  match. Several rounds can each be corrected in one read. A segment one
  substitution from several observed entries is rejected as ambiguous, unless
  all of them map to the same canonical sequence. One non-ACGT base in a round
  counts as that round's one substitution; two or more reject the read
  (`barcode_non_acgt`). There is no abundance-based tie resolution and no
  requirement for a whole-tuple exact supporting read.
- Replacement is applied **once after correction**, never chained. Multiple
  observed entries may map to one canonical sequence, allowing primer-barcode
  aliases or explicit well replacement. The caller must supply the correct
  mapping; no sample/well mapping is guessed or downloaded.
- Cell IDs concatenate canonical segments in the supplied order. Only observed
  tuples enter counting; the Cartesian barcode product is never allocated.
  IDs are deterministic and do not depend on allowlist ordering or threads.
- Process one biological sublibrary per invocation; lanes from that same
  sublibrary may share a run. Keep the experiment and sublibrary identity in
  downstream metadata when joining matrices. Distinct sublibraries must not
  be pooled under these sequence-only cell IDs.

The UMI follows the normal STAR path, including rejection of non-ACGT UMIs and
homopolymers, followed by the configured UMI collapse. Split barcode correction
is separate from UMI correction. The mode disables the simple-CB correction and
inline hash paths that assume a global whitelist. Conflicting assay-owned CLI
settings are errors.

## Outputs and scope

The standard `Solo.out/GeneFull/raw/` and `filtered/` MEX directories are emitted;
an enabled caller with zero called cells still produces a valid empty filtered
matrix. `--soloCellFilter None` is available for upstream-only checks.
`Solo.out/Splitseq.metrics.tsv` reports barcode acceptance, rejection reasons,
UMI rejection and per-round correction counts. Every read is counted in
exactly one of `barcode_exact`, `barcode_corrected`, `read_too_short`,
`barcode_non_acgt`, `barcode_no_match`, `barcode_ambiguous`, `umi_non_acgt` and
`umi_homopolymer`, so these sum to `reads_total`. `barcode_exact` and
`barcode_corrected` are the reads passed to alignment and counting; the
`roundN_corrected` rows count those reads again per corrected round. The
standard feature statistics describe the subsequent counting stage.

This implementation supports RNA with one feature: `Gene`, `GeneFull`,
`GeneFull_Ex50pAS` or `GeneFull_ExonOverIntron`.
Guide capture, BAM/tag output, automatic kit detection, sample demultiplexing,
Flex, two-pass and restart modes are not integrated. This is a STAR baseline,
not a reproduction of Parse's pipeline or cell-selection rule. The caller is
not tuned for this assay.

## Validation

### Tests in the repository

Decoder checks: `core/legacy/test/test_splitseq.cpp`
(`make -C core/legacy/source test_SplitSeq`).
FASTQ-to-MEX checks: `tests/test_splitseq_smoke.py` (small generated genome;
no external data). These cover correction in multiple rounds, a single
non-ACGT base, ambiguous and invalid barcodes, replacement aliases, 24-base
cell IDs, ordinary UMI collapse, default EmptyDrops, empty results, thread
consistency, nested and shared exons, and genomic multimapper rescue.
Two and four rounds: `tests/test_splitseq_rounds.py` runs both ends of the
range from FASTQ to matrix, the four-round layout with the largest cell
barcode the mode allows (31 bases), and checks that 32 bases and five rounds
are refused.
List tool: `tests/test_make_splitseq_lists.py` (generated lists, well
restriction, labels, errors, an output directory used again, a replacement
of its files that fails part way, and the bundled lists against their source
files).
`tests/run_perturb_addons_contract.sh` runs these Python tests together.
Public reads: `tests/run_perturb_addons_100k_regression.sh` runs the mode on
the first 100,000 read pairs of SRR41469160 (GSE352893) with the bundled
`48_wells` lists and the v1.13.0 rescue-off setting. It compares mapping counts, the raw and filtered matrices
and every row of `Splitseq.metrics.tsv` with
`tests/fixtures/perturb_addons_100k_golden.json`. The fixture is made by
`tests/setup_perturb_addons_100k_fixture.py` and is not tracked.

### Public libraries

On the first 100,000 read pairs of the public IGVF library IGVFDS5482QPWM, all
94,999 accepted reads that carry the submitter's barcode and UMI annotation
agree with it, and the raw matrix equals that of an independent decoder
followed by ordinary STAR counting on the same reads. Both were measured
before the single non-ACGT base rule, which accepts 12 more of those reads and
adds five matrix entries.

One public v3 library at full depth: GEO GSE352893, run SRR41469160
(97,458,666 read pairs, 48 round-1 wells), with the bundled `48_wells` lists,
the v1.13.0 defaults and a GENCODE v44 index with 38,606 genes; 5.5 minutes on
16 threads, 31 GB. The deposited matrix
(`GSE352893_iPSC_K562_mixed_integrated.h5mu`, 16,818 cells, named by well
indices, filtered to fewer than 5,000 counts per cell) comes from the
depositor's own processing. Locations and checksums are in
`tests/ARTIFACTS.md`.

| Check | Result |
| --- | --- |
| Reads with a valid barcode | 64.8% (60.4% exact, 4.4% corrected) |
| Deposited cells found among the barcodes counted, by well label | 16,818 of 16,818 |
| The 16,818 barcodes with the highest counts that are deposited cells | 16,419; the other 399 have 4,103 or more counts here (about 620 match when the well indices are shuffled) |
| Total counts per deposited cell, this mode / deposited | 0.880 (Pearson of the per-cell totals 0.9987) |
| Counts per gene on the 38,597 gene IDs present in both, summed over the deposited cells, this mode / deposited | 0.937 (Pearson over genes 0.908; 0.970 on log(1 + count)) |

- This compares one run with one depositor's processing under another
  annotation. The per-cell figure is about total counts and says nothing
  about individual genes; the per-gene figures show agreement, not
  equivalence. 6.1% of the deposited counts are in gene IDs that the index
  used here does not have (the deposit lists 62,710 gene IDs). No run with a
  matched annotation has been made, so the remaining difference is not
  attributed.
- Why a third of the reads have no valid barcode was looked at on the first
  2,000,000 read pairs. 24.4% of the barcode reads have neither linker within
  three mismatches at its position nor an exact copy of a linker anywhere.
  300,000 of those pairs were aligned: 62.2% mapped uniquely (18.8% to
  several places), and 186,478 of the 186,495 uniquely mapped pairs face each
  other, with a median template length of 275 bases. That is what fragments
  with cDNA at both ends and no barcode give; the other reads of the group
  were not examined, and how such fragments entered the library is not
  established. Of the reads that do have both linkers within two mismatches,
  95.5% decode (measured on the first 300,000 barcode reads).
- Decoding alone was measured, with the same rule and the unrestricted lists,
  on the first 120,000 to 360,000 barcode reads of seven other public
  libraries whose records name the kit; the seven were chosen before any was
  looked at, and nothing was aligned. Four decode at 83 to 85%
  (SRR34942214, SRR34992469, SRR38228195, SRR40643600). Three do not fit
  this layout as supplied and were not followed up: SRR33458973 decodes at
  69%, SRR34861592 at 27% and SRR25750549 at 0%; in the last two, round 3
  matches while rounds 2 and 1 fall away, as for a different barcode-read
  layout. These figures describe the libraries named, not every library of
  the kit.
- Reads whose linker is displaced by one base (about 2% in these libraries)
  are not decoded: positions are fixed.
