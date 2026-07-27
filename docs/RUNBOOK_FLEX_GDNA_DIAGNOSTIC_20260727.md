# Flex gDNA Diagnostic Runbook

## Objective

Add the Cell Ranger Flex genomic-DNA diagnostic to STAR's final
barcode/gene/UMI-family path without making BAM output mandatory. The result
must be computed from filtered molecules after barcode and UMI correction, must
cover both H0/H1 cache keeps and alignment fallback, and must remain inert for
non-Flex analyses.

This work belongs to the private 1.6.0 candidate together with the native
Visium HD Flex feature. After private validation, the complete candidate is
promoted to public `master` and released as public `v1.6.0`.

## Cell Ranger contract

The implementation follows the Cell Ranger 10.0.0 source at commit
`ae7fcd195bd5b0caafa928fd55904663328564bc`:

- `lib/rust/cr_lib/src/gdna_utils.rs`
- `lib/rust/cr_lib/src/fit_piecewise_linear_model.rs`
- `lib/rust/cr_lib/src/stages/get_gdna_metrics.rs`

For included, nondeprecated probes:

1. Count spliced and unspliced probes per gene.
2. Keep control genes having at least one probe of each region.
3. Restrict molecules to final filtered cell barcodes.
4. For a region-classified molecule, add
   `1 / number_of_region_probes_for_gene` to that gene's region count.
5. Apply `log1p` to the normalized spliced and unspliced counts.
6. Sort by log-spliced count and fit every internal pivot to:

   - `y = constant` before the pivot;
   - `y = constant + slope * (x - pivot)` at and after the pivot.

7. Select the minimum-RSS fit.
8. Report:

   - `estimated_gdna_per_probe = exp(constant) - 1`;
   - `estimated_gdna_unspliced_threshold = ceil(estimated_gdna_per_probe)`;
   - `estimated_gdna_fraction =
     min(gene_assigned_filtered_molecules,
         round(estimated_gdna_per_probe * total_unspliced_probes))
     / gene_assigned_filtered_molecules`.

The denominator includes filtered molecules with a resolved gene even when
their probe region is unknown or conflicting. Such molecules do not enter the
spliced/unspliced model axes, but Cell Ranger includes them in total filtered
UMIs. Molecules without a resolved gene are audited separately and excluded
from the denominator.

At least 10 control genes are required. This is Cell Ranger's
`GDNA_GENE_THRESHOLD`.

The diagnostic is probe-region based. It is not the fraction of reads aligning
to genomic chromosomes.

## Runtime design

### Probe metadata

Runtime metadata comes from the filtered probe CSV created by
`FlexProbeIndex`. STAR accepts an explicit `--soloFlexGdnaProbeSet` path and,
in `auto` mode, searches next to `--soloProbeList` and under
`<genomeDir>/flex_probe_artifacts/filtered_probe_set.csv`.

The loader retains:

- probe ID to `spliced`/`unspliced` region;
- gene ID to its spliced and unspliced probe counts;
- the ordered 15-bit gene axis from `probe_list.txt`.

Missing or malformed region metadata disables only the diagnostic. Mapping,
family correction, matrices, and ordinary scRNA-seq behavior remain unchanged.

### Per-family evidence

The existing Flex aggregate hash key remains:

`[CB20][UMI24][GENE15][TAG5]`

The 32-bit value stores a saturating 30-bit read count plus a two-bit region
state:

- `0`: no region evidence;
- `1`: spliced;
- `2`: unspliced;
- `3`: conflicting spliced/unspliced evidence.

This adds no per-molecule allocation. Region state is merged when thread hashes
are combined, when ambiguous barcodes are resolved, and when corrected UMIs
are re-keyed. A final molecule contributes to the model only when its state is
uniquely spliced or unspliced.

### Cache and alignment paths

H0/H1 cache format v3 preserves the 24-byte record size and stores the region
state in previously unused high bits of the 32-bit resolved-gene field. Readers
remain compatible with v1/v2 caches; old caches yield unknown region evidence
and therefore cannot produce a complete gDNA metric.

Alignment fallback obtains the region from the winning probe pseudo-contig's
probe ID. Genomic-only rescue has no exact probe-region identity and remains
unclassified for this diagnostic.

### Filtering and output

The estimator runs after clique UMI correction and after FlexFilter produces
the final passing barcode set for each sample, but before the aggregate hash is
released.

Outputs:

- `<soloFlexOutputPrefix>/flex_gdna_summary.tsv`
- `<soloFlexOutputPrefix>/<sample>/Gene/filtered/gdna_metrics.json`

The TSV contains per-sample rows and a library row with model parameters,
classified-molecule coverage, excluded unknown/conflicting molecule counts,
and a status. JSON exposes the two Cell Ranger-facing metrics plus the same
audit fields.

No metric is emitted as valid unless all of the following hold:

- Flex inline-hash mode is active;
- FlexFilter completed successfully;
- probe-region metadata loaded;
- the active H0/H1 cache, if used, carries v3 region metadata;
- at least 10 control genes and at least one gene-assigned filtered molecule are
  available.

## Gating

The feature is gated by `--flex yes`. It does not alter:

- non-Flex STAR/STARsolo;
- the non-Flex inline bridge;
- GEX or SLAM counting;
- barcode correction policy;
- UMI correction policy;
- gene resolution;
- matrix values.

`--soloFlexGdna auto|yes|no` controls diagnostics only:

- `auto` (default): run when complete metadata are available, otherwise write
  a skipped status;
- `yes`: require complete metadata and final FlexFilter cell calls, and fail
  clearly if either is unavailable;
- `no`: do not load metadata or compute the diagnostic.

## Test plan

### Tier A: deterministic tests

1. Probe CSV parsing, included/deprecated filtering, and path discovery.
2. Two-bit count/evidence packing, saturation, and merge truth table.
3. H0/H1 cache v1/v2 compatibility and v3 region round trip.
4. Evidence survival through thread merge and UMI re-keying.
5. Piecewise fit and metric calculation against fixed Cell Ranger-derived
   fixtures, including threshold rounding and insufficient-data status.
6. Non-Flex parameter and count-path invariance.

The deterministic Flex diagnostic test is part of the standard Flex test
suite.

### Tier B: authorized 100K fixture

Use the canonical 2024-A fixture once:

- FASTQs: `/storage/downsampled_100K/SC2300771`
- STAR index: `/storage/flex_filtered_reference_2024/star_index`
- probe CSV:
  `/storage/flex_filtered_reference_2024/star_index/flex_probe_artifacts/filtered_probe_set.csv`
- probe list:
  `/storage/flex_filtered_reference_2024/star_index/flex_probe_artifacts/probe_list.txt`
- Cell Ranger sanity oracle:
  `/storage/SC2300771_filtered_100K/cellranger/outs`

Generate a fresh v3 H0/H1 cache, run the current 100K Flex recipe with
`--soloMapqMode off`, and record:

- completion and binary SHA-256;
- raw and filtered MEX hashes;
- per-sample and library gDNA metrics;
- classified/unknown/conflicting molecule coverage;
- runtime and maximum RSS;
- comparison to Cell Ranger metrics where the same downsample oracle exists.

Use a fresh output directory. Do not reuse earlier derived STAR outputs.

#### 2026-07-27 acceptance record

The authorized fixture was executed once at:

`/storage/downsampled_100K/SC2300771/results/flex_gdna_100k_20260727_v1/`

The wrapper completion marker records success with a 44.79-second wall time,
38,020,388 KiB maximum RSS, STAR binary SHA-256
`c11d834a0937426b7591d8326ba7f2d4a76e5beb6cccfa71b6fe35f7cdc08358`,
and v3 cache SHA-256
`4b7fa324b1de428b5166e21c634b241192900588a77b10c1bd8416fe1e04d006`.
The cache contains 8,779,068 H0/H1 records and complete two-bit probe-region
metadata.

That single run also exposed and resolved one denominator bug without
re-executing the fixture. The model's molecule numerator and per-probe estimate
were already correct, but the first report divided by region-classified
molecules only. Applying the corrected, unit-tested denominator to the
run-audited fields gives:

| Scope | STAR estimate | Gene-assigned molecules | Corrected gDNA | Cell Ranger sanity value |
| --- | ---: | ---: | ---: | ---: |
| Library | 0.3868791783 per probe | 581,029 | 3.06095% | 3.08% |
| BC007 | 0.3415846800 per probe | 214,014 | 7.33737% | 7.32% |

The other sample rows follow the same deterministic correction:
`classified + unknown + conflicting`. Gene-unassigned molecules remain
separate audit records and do not enter the denominator.

### Release regression

After the 100K gate:

1. clean STAR rebuild;
2. standard parameter-generation regression;
3. Flex inline hash and alignment-fallback tests;
4. `soloSkipProcessing` tests;
5. barcode-without-feature bookkeeping tests;
6. resolver/materializer unit and parity tests;
7. SLAM parity and standard smoke tests;
8. verify the new gDNA test is enumerated by the Flex suite.

Only then merge to the private 1.6.0 lineage, build and hash the private binary,
promote with `git merge --no-ff` to public `master`, create public `v1.6.0`,
and push the public branch and tag.

## Failure handling

- Do not infer gDNA from genomic BAM alignment fraction.
- Do not silently treat v1/v2 cache keeps as spliced or unspliced.
- Do not bypass final filtered-cell selection.
- Do not count a conflicting-region family in either model axis.
- Do not change molecule matrices merely to produce the diagnostic.
- On incomplete metadata, emit an explicit skipped/incomplete status in auto
  mode; fail only when the user selected `--soloFlexGdna yes`.
