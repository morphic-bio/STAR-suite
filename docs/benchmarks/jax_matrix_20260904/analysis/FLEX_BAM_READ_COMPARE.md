# Native Flex BAM/read difference analyzer

`flex_bam_read_compare` streams a Cell Ranger Flex BAM (through `samtools view`)
or SAM and reports where its read evidence falls relative to STAR's filtered
cell calls and sample-aware probe-cache policy. It does not link htslib and does not
modify either input.

Cell identity is always the corrected composite `CB16|TAG`: the two occurrences
of a CB16 under two different sample tags are distinct cells. Supply STAR's
per-tag filtered barcode lists as `TAG=PATH`; Cell Ranger barcode files normally
already contain the 24-base `CB16+TAG8` identity.

Build and test:

```bash
make -C docs/benchmarks/jax_matrix_20260904/analysis
make -C docs/benchmarks/jax_matrix_20260904/analysis test
```

Example for one L004 per-sample BAM:

```bash
analysis/flex_bam_read_compare \
  --input SAMPLE/count/sample_alignments.bam \
  --out-prefix OUTPUT/SAMPLE \
  --tag-map refs/sample_whitelist_full_16.tsv \
  --star-cells BC001=STAR/per_sample/BC001/Gene/filtered/barcodes.tsv \
  --star-cells BC002=STAR/per_sample/BC002/Gene/filtered/barcodes.tsv \
  --cr-cells SAMPLE/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz \
  --hash-cache refs/h01_cache.bin \
  --gene-list refs/star_index/flex_probe_artifacts/probe_list.txt \
  --samtools-threads 8
```

The exact BAM and filtered-matrix filenames should be taken from the run rather
than inferred from this example. The gene list must be the exact ordered
`probe_list.txt` used to generate the cache because cache genes are stored as
1-based indices. `--cr-cells` can be repeated, and accepts `TAG=PATH` when a
barcode file contains CB16 rather than CB16+TAG8.

Outputs:

- `.metrics.tsv`: call-set overlap and stream/cache provenance counts.
- `.cells.tsv`: every exact composite cell and its shared/STAR-only/CR-only class.
- `.summary.tsv`: exhaustive read strata by cell class, cache verdict, CR `xf&8`
  status, and gene agreement.
- `.genes.tsv`: exhaustive gene-pair strata for prioritizing systematic misses.
- `.cigars.tsv`: exact CIGAR and MAPQ-bin strata.
- `.regions.tsv`: CR `RE` crossed with cache `probeRegion`, cache verdict,
  quantified `fx` gene, and the separate genomic `GX`/`GN` evidence.
- `.counted_records.tsv`: descending outcome/gene/probe strata restricted to
  primary records with `xf&8`. The tool does not deduplicate CB+UB; this is a
  counted-record report, not an independently reconstructed molecule table.
- `.molecules.tsv`: true molecule-level strata. The molecule key is packed
  corrected CB24 + corrected `UB` + quantified `fx` gene. `UR` is used only
  when `UB` is absent, is reported as `UR_fallback`, and has separate metrics.
  A key enters the report only when at least one of its primary records has
  `xf&8`, but its best cache evidence is selected across every primary record
  for the key. Evidence priority is same-gene KEEP, different-gene KEEP, MISS,
  then DENY/other. Thus `rescued_same_gene` identifies molecules whose `xf&8`
  representative was not a same-gene KEEP but another read was.
- `.details.tsv`: relevant BAM tags and per-read evidence. By default this is a
  deterministic 100,000-read sample of discordant evidence, selected by QNAME
  hash; use `--detail-max 0` for every eligible record or `--details none` to
  suppress detail rows.

Primary alignments are analyzed by default. Secondary and supplementary records
are counted in `.metrics.tsv` but excluded from all other reports unless
`--include-secondary yes` is requested. Reverse-strand SAM sequences and
qualities are restored to original read orientation before cache replay. Cache
replay mirrors the sample-aware full classifier: it queries physical offsets 0
and +1 (the configured -1 offset is skipped because the probe start is zero),
prefers a sample-specific record, falls back to sample index zero, and applies
the existing conflict/deny precedence. Optional single-`N` replay is off by
default because it is not part of that production classifier.

`fx` is treated as the quantified gene assignment and `pr` as probe evidence.
`GX`, `GN`, `RE`, reference locus, and CIGAR remain separate explanatory
alignment evidence; they are never substituted for `fx`.

The molecule table uses a flat open-addressed table with 2-bit CB24/UMI keys
and interned gene/evidence strings. `.metrics.tsv` reports all valid molecule
keys, unique counted keys, `xf&8` records represented, keys with multiple
`xf&8` records, invalid/missing CB/UMI/`fx`, UR fallbacks, table capacity/bytes,
and a total for every molecule outcome. MAPQ 255 is reported as
`255/unavailable`, not in the `60+` bin.
