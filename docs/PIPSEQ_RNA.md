# PIPseq RNA BI/IMI counting

`--soloPipseqMode Illumina` counts RNA for the newer PIP/Illumina single-cell
chemistry, whose barcode read carries a four-tier cell barcode and a three-base
binning index (BI) and has no UMI. It does not process the older v3 chemistry
with a 12-base UMI, or guide-capture libraries.

The mode is written from public descriptions only:

- the published DRAGEN documentation of
  [transcript counting](https://help.dragen.illumina.com/dragen-v4.4/product-guide/dragen-v4.4/dragen-single-cell-pipeline/dragen-scrna-pipseq)
  and of [barcode processing, read classification and the STAR aligner option](https://help.dragen.illumina.com/dragen-v4.6/product-guides/dragen-v4.6/dragen-single-cell-pipeline/dragen-scrna-illumina);
- Illumina's note on [aligning PIPseq data with STAR](https://developer.illumina.com/news-updates/your-pipseq-workflow-is-consolidating-into-dragen)
  (15 July 2026);
- the documentation of the [default annotation and ignored biotypes](https://help.multiomics.illumina.com/dragen-single-cell-rna/run-set-up-in-bssh/manual-launch)
  and of the [annotation files](https://support.illumina.com/sequencing/sequencing_software/dragen-bio-it-platform/product_files.html);
- the [public kit guide](https://crtp.ccr.cancer.gov/wp-content/uploads/fsgc-protocols/PIPseq-V-T10-3-Single-Cell-RNA-Kit-User-Guide-1.pdf).

No vendor software and no packaged reference file is needed or was used. The
mode states its own rule wherever those descriptions leave one open, and does
not claim to reproduce another pipeline's output. A comparison with a
deposited matrix is in [Validation](#validation).

## Run

```bash
# once per genome: the reference (see Reference below)
python3 scripts/make_pipseq_reference_gtf.py \
  --gtf gencode.v44.annotation.gtf.gz --out-gtf pipseq_ref/genes.gtf --report pipseq_ref/report.json
core/legacy/source/STAR --runMode genomeGenerate --runThreadN 16 \
  --genomeDir pipseq_ref/index --genomeFastaFiles genome.fa \
  --sjdbGTFfile pipseq_ref/genes.gtf --sjdbOverhang 100

# once per kit: the tier lists (see Tier lists below)
python3 scripts/make_pipseq_tier_lists.py --barcode-fastq R1.fastq.gz --out-dir tiers

core/legacy/source/STAR \
  --runThreadN 16 \
  --genomeDir pipseq_ref/index \
  --readFilesIn R2.fastq.gz R1.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix /path/to/fresh_output/ \
  --soloPipseqMode Illumina \
  --soloPipseqWhitelist tiers/readderived_tier1.tsv tiers/readderived_tier2.tsv \
                        tiers/readderived_tier3.tsv tiers/readderived_tier4.tsv \
  --soloFeatures GeneFull_Ex50pAS \
  --soloStrand Forward \
  --outSAMtype None
```

The cDNA read (R2) is first and the barcode read (R1) last. Lanes of one
library go in as comma-separated lists. Process each library in its own run:
cell barcodes carry no library label.

The raw FASTQ pair is the input. The mode trims the cDNA read itself; there is
no separate trimming step.

### Reference

Any STAR index works. The reference meant for this mode follows the gene list
that Illumina documents for the chemistry, less the read-through genes. Build
it from a GTF written by `scripts/make_pipseq_reference_gtf.py`:

- **Annotation.** The documentation names GENCODE v44 as the default for
  GRCh38: the comprehensive annotation on the reference chromosomes
  (`gencode.v44.annotation.gtf.gz`, 62,700 genes), not the primary-assembly
  file, which adds 54 genes on unplaced scaffolds.
- **Genes the documentation leaves out.** It says the biotypes pseudogene,
  shortRNA and rRNA are ignored with this kit. The tool drops every gene type
  whose name contains `pseudogene` (15,233 genes in v44), the short RNA types
  `misc_RNA`, `snRNA`, `miRNA`, `snoRNA`, `scaRNA`, `sRNA` and `scRNA` (6,989)
  and `rRNA` and `Mt_rRNA` (49). 40,429 genes remain, among them lncRNA, TEC
  and the mitochondrial tRNA genes. `--drop pseudogene,shortRNA,rRNA` writes
  this list.
- **Where that assignment comes from.** The documentation names the three
  classes and not the gene types in them. The assignment was read from the
  feature list of a public deposit (GSE295097), which it reproduces exactly;
  `--check-features` repeats that comparison for any features file.
- **Read-through genes, left out by default.** A read-through gene joins
  neighbouring genes and shares their exons, and STAR leaves out a read that
  fits two genes. The tool drops the transcripts GENCODE tags
  `readthrough_transcript` and the genes left without a transcript: 906 genes,
  and 176 transcripts of 65 further genes. 39,523 genes remain. The numbers
  behind this choice are in [Validation](#validation).
- **Genome.** The sequence is not filtered. Unplaced scaffolds stay in the
  index as places a read can align; they carry no gene.

A reference built for another assay differs in the gene list. STAR Suite's
2024-A-style list (38,606 genes from the primary-assembly file of the same
GENCODE release) has no mitochondrial tRNA genes and no TEC genes, and it
holds 236 immunoglobulin and T-cell-receptor pseudogenes and 24 protein-coding
genes on unplaced scaffolds.

### Tier lists

`--soloPipseqWhitelist` takes four one-column lists in read order, with
sequences of 8, 6, 6 and 8 bases. None is bundled.
`scripts/make_pipseq_tier_lists.py` reads them off the barcode reads: among
reads in which the three constant linkers match at exactly one phase, the 96
most frequent sequences at each tier position. Its `summary.json` shows how far
the last listed sequence stands above the next. A library may use fewer than 96
sequences in a tier. With the default `--soloPipseqMultiMatch Reject` the
output does not depend on the order of a list.

### What the mode sets

| Setting | Value | Source |
| --- | --- | --- |
| cDNA trimming (`--soloPipseqTrim Kit`) | one leading base; adapters `AGAGTGAATGGG` (template switch) and `TCAACGCAGAGT` (SMART PCR) at the 5' end; poly-A and poly-G at the 3' end; reads left under 20 bases are not aligned | DRAGEN 4.6 documentation, STAR aligner section; Illumina note |
| `--alignIntronMax` | 25000 | same |
| Alignments counted | reads with one alignment | this mode's rule; the documentation counts the primary alignment above a mapping-quality limit |
| `--soloFeatures` | `GeneFull` unless given | one of `Gene`, `GeneFull`, `GeneFull_Ex50pAS`, `GeneFull_ExonOverIntron` |
| `--soloCellFilter` | `None` unless given | |

A value given on the command line replaces a default. The barcode, UMI and
deduplication parameters belong to the mode; a conflicting value is an error.

**Trimming.** `--soloPipseqTrim Kit` sets `--clipAdapterType PIPseq`,
`--clip5pNbases 1`, `--clip5pAdapterSeq AGAGTGAATGGG,TCAACGCAGAGT` and
`--clip3pAdapterSeq polyA,polyG`. An adapter is found by overlap: at least 12
aligned bases, mismatches up to `--clip5pAdapterMMp` (0.1) of them, no indels;
the read is cut after the adapter. A 3' homopolymer is a window of 12 bases
with at most one other base; the read is cut at the window's start. Give other
sequences with the clip parameters. Use `--soloPipseqTrim None` for reads that
were trimmed beforehand.

**Gene assignment.** The documentation describes reads that are at least half
exonic as assigned first, then reads inside a gene, with antisense and
multi-gene reads left out. `GeneFull_Ex50pAS` is the closest STAR rule.
`GeneFull` drops every read that lies in two overlapping genes, which removes
whole genes that sit inside another gene.

### BAM output

`--outSAMtype BAM Unsorted` or `SortedByCoordinate` may be used. With
`--outSAMattributes ... CR CB UR GX GN sF`, each alignment carries:

- `CR`: the barcode as read (the four tier sequences at the chosen phase);
- `CB`: the corrected barcode, the same string as in `barcodes.tsv`;
- `UR`: the three BI bases;
- `GX`, `GN`, `sF`: the gene and the overlap type of the alignment.

`CR`, `CB` and `UR` are `-` where undefined (no phase, tied phases, rejected
barcode). `UB` is not defined: there is no molecule identifier.

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
rule, above.

Phase ties are rejected: the published description gives no tie rule.
Non-ACGT bases in a candidate barcode tier invalidate that phase. A valid
barcode with a missing/non-ACGT BI retains its barcode metrics but contributes
no RNA count. All 64 A/C/G/T BIs, including homopolymers, are valid and uncorrected.

The corrected cell identifier is the concatenation of the four tier sequences,
without separators or an added suffix. Lookup tables cover the individual
tiers; only observed cell/gene combinations are stored. The 96^4 Cartesian
product is never allocated.

Within each cell/gene, exact `(BI, IMI)` observations collapse. The IMI is the
global genomic coordinate of the strand-aware aligned 5-prime base, with strand
included in the key. Soft clipping is not extrapolated. Reverse reads use the
last aligned base; internal splice structure is not part of the key. These are
this mode's conventions.

Each mapping thread keeps its distinct fragments as flat 24-byte records and
compacts them as it goes; the threads' records are merged after mapping, so
duplicate reads in different threads collapse. Memory grows with the number of
distinct fragments; there is no spilling to disk.

## Transcript estimator and outputs

Let B be occupied BI bins and I distinct BI/IMI pairs for one cell/gene:

```text
B <= 10: count = B
B > 10:  count = max(B, floor(I / IPM))
```

With `--soloPipseqIPM 0`, calibrate using all assigned cell/gene groups
with 5–32 occupied bins, before cell filtering. The coupon-collector expectation
is the stopping-time mean `E[M|reach B] = sum(j=0..B-1, 64/(64-j))`.
IPM is the arithmetic mean of I/E across eligible groups. This is this mode's
reading of the documentation.

Effective IPM is at least one, enforcing the documented bound B <= count <= I.
Both the raw estimate and this floor are reported. A finite explicit IPM >=1
overrides calibration. With no eligible groups, high-occupancy groups cause an
error requesting an explicit IPM or more reads; if every group has at most ten
bins, IPM is unused and counts can still be emitted. There is no silent
high-occupancy fallback.

Under `Solo.out/<feature>/`:

- `raw/{matrix.mtx,features.tsv,barcodes.tsv}`: estimated transcript counts.
- `filtered/`: the same units after optional standard cell calling.
- `Pipseq.groups.tsv`: barcode, gene ID, assigned reads, occupied bins,
  unique BI/IMI pairs and estimated transcripts.
- `Pipseq.metrics.tsv`: effective counting policies (including `cdna_trim`),
  calibration population, raw/effective IPM, calibration size and barcode/read
  rejection categories. `reads_total` equals `barcode_exact +
  barcode_corrected + barcode_invalid + barcode_nonmatch + barcode_phase_tie +
  barcode_ambiguous`.
- `Summary.csv`: read, barcode and mapping fractions, distinct BI/IMI
  fragments, estimated transcripts, IPM and, with a cell filter, the called-cell
  summary. It is written from the PIP counters; `Barcodes.stats` and
  `Features.stats` of the ordinary UMI path are not produced in this mode.

The raw barcode axis contains cells with at least one counted RNA observation.
Inputs with no accepted RNA observations yield a valid empty matrix. Matrix
counts and offsets are checked against STAR's 32-bit matrix limits. The
sample-level calibration and resulting counts depend on the reads supplied; a
downsample is not expected to reproduce full-library counts.

Standard `--soloCellFilter` options, including `EmptyDrops_CR`, can filter the
estimated-transcript matrix with STAR's existing callers. The threshold rule of
the published description is not implemented.

Not in this mode: Flex, two-pass, restart, pfMulti guide integration, the v3
chemistry, guide reads.

## Validation

### Tests in the repository

- `core/legacy/test/test_pipseq.cpp` (`make -C core/legacy/source test_PipSeq`):
  decoding in all phases and bins, exact-first matching in both list orders,
  one and several one-substitution candidates under both
  `--soloPipseqMultiMatch` values, phase ties, malformed input, the fragment
  store (collapse, merging of stores, compaction), estimator edges and the
  trimming functions.
- `tests/test_pipseq_smoke.py`: FASTQ to matrix on a small generated genome,
  including reverse reads, duplicate fragment starts, fixed and estimated IPM,
  cell filtering, the mode's trimming against the same reads without technical
  sequence, BAM tags (unsorted and sorted), exon-priority gene assignment and
  ordinary UMI counting with the mode off.
- `tests/test_pipseq_edges.py`: EmptyDrops, empty matrices, merging across
  threads and calibration failures, on the smoke fixture's index.
- `tests/test_make_pipseq_reference_gtf.py`: the reference tool on a small
  generated GTF: each dropped class, read-through genes and transcripts, the
  documented list, both attribute names and the features-file check.

```bash
python3 tests/test_pipseq_smoke.py --star core/legacy/source/STAR --out tests/pipseq_output_example
```

### A public library with a deposited matrix

GSE295097 / GSM8942285 (run SRR33231853, 180,144,123 read pairs) has a
deposited filtered matrix of 5,675 cells. The GEO record names PIPseeker 3.3.3
and the paper names the DRAGEN Single Cell RNA app on BaseSpace; neither
states a version, an annotation or a gene filter, and neither links code. The
file names follow the published DRAGEN documentation, and the 40,429 features
are exactly the documented default list: GENCODE v44 on the reference
chromosomes without the pseudogene, shortRNA and rRNA classes. Sources and
checksums are in `tests/ARTIFACTS.md`. Only the filtered matrix is deposited:
there are no alignments and no unfiltered matrix, so nothing can be compared
read by read.

The runs used the raw read pair, the mode's defaults, tier lists derived from
reads, 16 threads and the reference of [Reference](#reference) on a GRCh38
primary-assembly genome. Its 39,523 genes are all in the deposit's list; the
906 read-through genes it leaves out hold 0.2% of the deposited counts.

| | `GeneFull` | `GeneFull_Ex50pAS` |
| --- | ---: | ---: |
| Reads with an accepted barcode | 167,966,621 (93.2%) | same |
| Reads counted | 123,273,928 | 122,362,778 |
| Estimated transcripts, all barcodes | 38,704,199 | 38,410,374 |
| IPM | 1.484 | 1.482 |
| Deposited cell labels that have our counts | 5,675 of 5,675 | 5,675 of 5,675 |
| Per-cell total, ours / deposited | 1.068 | 1.063 |
| Per-cell totals, Pearson on log10 | 0.9997 | 0.9997 |
| Per-gene totals, Pearson on log1p | 0.927 | 0.956 |
| Genes with counts only in the deposit | 3,274 | 2,165 |
| (cell, gene) entries in both; equal | 11,034,381; 85.0% | 11,128,511; 84.8% |
| Wall time; peak memory | 480 s; 33.4 GB | 483 s; 33.7 GB |

The index accounts for about 30 GB of that memory. Every deposited label is
among our barcodes, with at least 6,074 reads.

#### The gene list

The same library was counted on three gene lists, each from GENCODE v44 on the
same genome. Counts are totals over the deposited cells with
`GeneFull_Ex50pAS`:

| Gene group | Genes | Deposited | 2024-A-style list | Documented list | Reference |
| --- | ---: | ---: | ---: | ---: | ---: |
| In every list, not touching an added gene | 36,492 | 23,559,202 | 24,858,738 | 24,971,711 | 24,958,906 |
| In every list, sharing a span and strand with an added protein-coding or lncRNA gene | 1,854 | 1,852,322 | 2,156,490 | 1,575,375 | 2,136,393 |
| Added protein-coding and lncRNA genes, 906 of them read-through | 977 | 61,438 | not listed | 69,796 | 2,693 |
| Mitochondrial tRNA | 22 | 484,807 | not listed | 494,499 | 494,301 |
| TEC | 1,056 | 10,561 | not listed | 10,045 | 10,045 |
| Ribozyme, artifact, vault RNA | 28 | 287 | not listed | 179 | 179 |
| Per-cell total, ours / deposited | | | 1.040 | 1.044 | 1.063 |

- **The 2024-A-style list** (38,606 genes) lacks genes that hold 2.1% of the
  deposited counts, most of it in the mitochondrial tRNA genes. Its ratio of
  1.040 is low for that reason.
- **The documented list** (40,429 genes) brings those genes in, and they agree
  with the deposit. Its read-through genes take reads away from their
  neighbours, because a read in a shared exon fits two genes: MIF falls from
  43,807 to 3,455 (deposited 35,684), TUBB3 from 60,852 to 18,914 (53,206) and
  EEF1G from 36,388 to 2,474 (33,677). With `GeneFull` the whole group falls
  to 492,925. The deposit shows no such loss; how it assigns these reads is
  not documented.
- **The reference** (39,523 genes) keeps the added genes and restores the
  neighbours. On it, our counts are 6.3% above the deposit per cell.

#### Where the two differ

The figures below come from an earlier `GeneFull_Ex50pAS` run on the
2024-A-style list and its tagged BAM, except where the reference is named;
every gene named is in all three lists.

- **About 10% more observations per entry in ours.** For an entry with B
  occupied bins in ours, the deposit holds about 0.90 B, at every B up to ten;
  above ten bins the same factor passes through the division by IPM. In
  entries whose bins all have several reads the deposit holds 95–97% as many
  counts as we have bins; in entries whose bins all have one read, 74–87%.
  Reads in entries the deposit lacks align over
  shorter stretches (25% under 50 aligned bases, against 5% in entries that
  agree). These entries follow expression: across genes their number has a
  rank correlation of 0.94 with all our bins, over 27,000 genes. The deposit
  in turn has 856,052 entries we lack; 40% of their counts are in 100 genes.
- **Reads of multi-copy genes.** Ribosomal-protein genes are 19% higher in
  ours on the reference. Genes that contain a same-strand processed pseudogene which is not in
  the deposit's annotation are higher in the deposit (for example PLEKHA7 with
  `RPL36AP37`, SNX6 with `RPS19P3`, OOEP with `RPL39P3` and `RPS6P8`). The 100
  genes that hold most deposit-only counts total 711,440 there and 114,018 in
  ours; our alignments in them are mostly reads with several alignments
  (612,198 distinct cell and BI pairs), which this mode does not count, and
  few are left out by a gene rule (about 38,000). Trimming decides where such
  reads go: with the poly-A left on, STAR places them uniquely on the
  pseudogene copy, which carries a poly-A stretch in the genome (unique sense
  reads in the first four million, untrimmed against trimmed: PLEKHA7 2,313
  against 240, SNX6 1,083 against 199).
- **SLC26A3** has 102,220 deposited counts, in 5,672 of the 5,675 cells, and 2
  in ours. STAR places no read on the gene's strand there in the first four
  million reads; the reads it does place in the gene are poly-T reads that
  begin with reversed barcode sequence and have several alignments.
- **Gene ends.** `GeneFull_Ex50pAS` leaves out reads that run past the
  annotated end of a gene unless half of the read is exonic; `GeneFull` counts
  any overlap. The documentation's rule (a read overlapping a gene by at least
  20% of its length) lies between the two and is not available in STAR.
  MT-CO2 shows it: 118,960 deposited, 62,768 with `GeneFull_Ex50pAS` on the
  reference.
- **Cell set.** 675 barcodes outside the deposit total at least as much as the
  lowest deposited cell. None is one substitution from a deposited label;
  their genes detected and mitochondrial share are close to deposited cells of
  the same size. The documented threshold rule, applied to our totals with the
  documented expected-cell values of 2,000, 10,000 and 100,000, gives 4,298,
  18,101 and 25,946 cells.

**Aligners.** Illumina's note of 15 July 2026 compares its own aligner with
STAR on this chemistry. It reports the largest gene-level differences at
pseudogenes, lncRNAs, HLA and other repetitive loci, attributes them to how
each aligner scores a read with several placements, and says settings cannot
remove them. The multi-copy differences above are of that kind, and this
comparison does not take them further.

Checked and found to change little: trimming with cutadapt instead of the
mode's own, and no intron limit (per-gene correlation 0.929 against 0.930);
merging fragment starts within two bases of each other (1.3% fewer counts
above ten bins; 3.6% of neighbouring starts in a bin are that close).

### Existing modes

The same commands run with the `master` binary (024b72f) and the branch binary
give identical output: complex barcodes on two libraries, whitelist counting
with adapter clipping and BAM barcode tags, plain alignment with clipping, and
the Flex default-route test with `master` as the reference binary.
