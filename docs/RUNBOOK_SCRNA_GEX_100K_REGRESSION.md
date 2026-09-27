# Plain scRNA-seq Regression Gates

## Scope

Plain 10x GEX must be tested independently of Flex and perturb-seq. Production
use of `MultiGeneUMI_CR` does not exercise all conventional Solo UMI methods.
The existing UCSF GEX+guide production case remains a separate gate.

## Fast, Analytical Gate

The existing `tests/run_solo_smoke.sh` now uses three uniquely mapped reads
with two known UMIs. It requires all three primary BAM alignments and exactly
two molecules at the expected gene/barcode coordinate. File existence or a
header-only MEX cannot pass. Use a fresh `SOLO_SMOKE_OUTDIR` to retain a run
outside the default `tests/solo_smoke/` directory.

```bash
STAR_BIN=/path/to/STAR python3 tests/test_scrna_gex_counts.py
```

The synthetic fixture has four intended cells, one low-count barcode, two genes,
exonic and intronic reads, PCR duplicates, and one-mismatch UMI neighborhoods.
It uses paired gzip FASTQs with poly-G tails. Expected counts are specified from
the fixture, not copied from another STAR execution.

It checks raw and filtered Gene/GeneFull matrices for default Solo, `1MM_All`,
`1MM_CR`, both directional methods, `Exact`, `NoDedup`, and the production-style
`1MM_CR` + `MultiGeneUMI_CR` combination. A deterministic CellRanger2.2 cutoff
checks that the four intended cells survive and the low-count barcode does not.
This is not an EmptyDrops sensitivity/specificity test.

This gate runs in Tier A and against independently built portable/static core
binaries in the partial-Make CI matrix.

## Real 100K Gate

Dataset: public 10x PBMC 10K v3, GEX only. The fixture contains exactly 100,000
matched read pairs: the first 50,000 pairs from each of L001 and L002. It does not
use Flex probes, a feature reference, CRISPR, LARRY, or a multi-library config.

- Original FASTQs: `/storage/pbmc10k_stage/pbmc_10k_v3_fastqs/`.
- Prepared fixture: `/storage/downsampled_100K/pbmc10k_gex/`.
- Reference: `/storage/autoindex_110_44/bulk_index`, GRCh38-2024-A vintage.
- Whitelist: `/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt`.
- Golden summaries: `tests/fixtures/scrna_gex_100k_golden.json`.

Recreate only into a new directory:

```bash
python3 tests/setup_scrna_gex_100k_fixture.py \
  --r1 /path/to/pbmc_10k_v3_S1_L001_R1_001.fastq.gz \
  --r2 /path/to/pbmc_10k_v3_S1_L001_R2_001.fastq.gz \
  --r1 /path/to/pbmc_10k_v3_S1_L002_R1_001.fastq.gz \
  --r2 /path/to/pbmc_10k_v3_S1_L002_R2_001.fastq.gz \
  --outdir /new/fixture/directory
```

Preparation verifies mate names and FASTQ structure, then records file checksums
and lane provenance. Runs verify fixture/whitelist hashes and reference metadata
before launching STAR. The large genome index itself is not rehashed per smoke;
exact output digests provide an additional reference-drift check.

```bash
bash tests/run_production_module_regression_suite.sh \
  --preflight --case gex-only-pbmc-100k --strict-preflight

STAR_BIN=/path/to/STAR bash tests/run_scrna_gex_100k_regression.sh \
  --outdir /new/run/directory
```

The default run executes three profiles serially:

1. Vanilla Solo defaults, no BAM.
2. `--defaultCoreScrna Yes`, no BAM (EmptyDrops_CR and CR-style UMI filtering).
3. The same modern profile with sorted BAM and CB/UB tags.

All request Gene and GeneFull and enable CellRanger4/poly-G clipping. Every
profile must read 100,000 pairs and produce nontrivial raw counts. The gate
checks exact biological-key matrix digests, feature/barcode sets, filtered-cell
identity, and integer mapping QC against the reviewed baseline. Modern counts
and cell calls must be identical with and without BAM; `samtools quickcheck`
validates the BAM. Commands, binary SHA/version, reports, and logs are retained.

A shallow 100K sample is a regression fixture, not a biological cell-count
estimate. If a caller legitimately returns zero cells, the gate checks that
explicitly; zero raw counts never pass. The analytical gate separately requires
positive counts and a known positive filtered-cell set.

The 100K case is fixture-backed and is not run on ordinary GitHub-hosted PR
runners. Select it explicitly from the production-module suite on a host with
the fixtures. Serialize it with other production tests and benchmarks.

## Baseline Policy

`--record-golden /new/path.json` is an explicit bootstrap operation. It refuses
to overwrite a file and does not report a regression pass. Review the analytical
test and old/new comparisons before adopting a newly recorded golden; do not
regenerate expected counts just to make a failure disappear.

Initial validation artifacts: `/tmp/star-htslib-build-20260927/gex100k-*`.
The released v1.9.5 vanilla profile is NOT a valid golden: it generated empty
raw matrices despite successful mapping. The modern released profiles are
retained as an independent unchanged-path comparator. The corrected vanilla
profile is established only after the analytical UMI-method tests pass.

Validation on 2026-09-27, 100,000 read pairs per profile:

| Profile | Raw Gene UMIs | Raw GeneFull UMIs | Result |
| --- | ---: | ---: | --- |
| Released v1.9.5 vanilla | 0 | 0 | Rejected; reproduced on a clean rebuild |
| Corrected vanilla | 47,358 | 70,527 | Nonempty analytical-method correction |
| Released modern | 47,203 | 70,270 | Independent unchanged-path baseline |
| Corrected modern | 47,203 | 70,270 | Exact keyed counts/cell-set/QC match |
| Corrected modern + BAM | 47,203 | 70,270 | Exact match to modern without BAM |

All three corrected profiles map 89,527 pairs uniquely. The golden metadata
identifies the released comparator and the corrected vanilla source separately.
The v1.9.5 binary artifact and APT `.deb` installation checks passed on Ubuntu
22.04 and 24.04; successful installation did not detect this independent runtime
counting bug. The public release was not changed by this local correction.

### Separate Observation

The initial tiny fixture also showed that `TopCells 4` selected five barcodes:
the current implementation uses the zero-based element at `topCells` as its
inclusive UMI threshold. That pre-existing behavior is not changed here. The
analytical test uses an explicit CellRanger2.2 threshold instead; a TopCells
boundary/tie-policy change should be reviewed separately.
