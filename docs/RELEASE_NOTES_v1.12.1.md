# STAR Suite v1.12.1 Release Notes

Date: 2026-10-16

This patch makes OrdMag bootstrap draws independent of the execution worker
count and includes the comprehensive-QC read-length fix through `cd3dadc`.

Both bootstrap stages create a separate PCG32 generator for every resample.
For recovery resample `b`, the generator uses the configured nonzero bootstrap
seed (or default 1) and stream `b`. Top-N uses seed plus 10,000 (or default
100) and stream `10,000 + b`. Resample indices start at zero. Existing seed
options and their mapping into the bootstrap configuration remain unchanged.
Each resample produces the same draws regardless of which worker computes it;
results are averaged in resample-index order. Worker count, detected hardware
concurrency and `OMP_NUM_THREADS` control execution, without changing draws.

The UMI/quality/barcode ranking, its optional mitochondrial tie metadata and
the prefix-selection rule are unchanged. The new stream rule makes a one-time
change from v1.12.0 bootstrap estimates and can add or remove boundary cells.
The validation record reports the affected cells and UMI totals. Selecting the
old worker count is not a way to reproduce the previous release’s streams.

Comprehensive SLAM QC initializes mate lengths, ignores workers that processed
no reads and uses the run’s mate count. Empty, single-end and paired-end
inputs have explicit regression coverage, including allocator perturbations.
The fix changes QC length reporting; it does not change mapping or counting.

`STAR --version` reports `1.12.1`. Debian source packaging uses `1.12.1-1`;
Ubuntu packages use `1.12.1-1~ubuntu22.04.1` and `1.12.1-1~ubuntu24.04.1`.
Upstream STAR remains `2.7.11b`, genome-index compatibility remains `2.7.4a`
and legacy compatibility remains `2.7.1a`. Existing indexes need no rebuild.

The public PBMC 100K GEX regression has a corresponding golden-fixture update:
GeneFull retains 10,530 rather than 10,540 cells. The ten removed boundary
barcodes each have one UMI and one detected gene. All raw files and the quality
ranking are unchanged, and both retained sets are exact ranking prefixes.
The Gene profile and all other profile quantities are unchanged. Independent
modern unsorted and sorted-BAM runs match for every profile quantity with the
new streams. The fixture records the affected barcodes and original provenance.

On the bounded CAT-ATAC DMSO1 input (two million GEX reads, seed 777), runs
with 1, 4, 16 and 32 workers retain the same 12,173 barcodes in the same order;
their filtered matrices and caller debug summaries are identical. Relative to
the previous 16-worker calls, 18 boundary cells are removed, each with 29 UMIs;
no cells are added. Raw counts and quality ranking are unchanged. The previous
one-worker calls lose 13 cells, also all with 29 UMIs.

Other shallow release fixtures illustrate the one-time boundary change: the
UCSF GEX fixture adds two one-UMI cells, and the public A375 scheduler fixture
now retains 1,135 cells in every full-depth scheduler arm. Legacy downstream
doublet and ambient-guide scores can change when their input cell set changes.
Restoring each previous cell mask on the new raw counts reproduces the old
seeded downstream values exactly, including all guide q-values.

Validation includes all ten EmptyDrops test executables; the bootstrap matrix
at 1, 2, 4, 16 and 32 workers, automatic worker selection and OMP unset/2/32;
the old implementation rejected by the new regression; 15 feature-module test
executables; 602 server tests; all 25 release-manifest rows with the previously
approved public fixtures; nine partial builds; 14 container smokes; seven
tarball/installer runtime checks; and three Debian install/purge combinations.
The four SLAM stage harnesses pass and 58 mapping/quantification tables are
byte-identical to v1.12.0. Comprehensive QC is not a whole-file identity claim:
the fixed uninitialized lengths change its reports, and variable-length reads
retain the existing sampled-worker length convention.
