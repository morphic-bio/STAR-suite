# STAR Suite v1.13.0 Release Notes

Date: 2026-11-17

STAR Suite 1.13.0 adds two experimental add-ons for STAR-perturb, RNA
counting for split-pool barcoding and for the newer PIP chemistry, and
corrects cell filtering of a raw matrix. Counting and cell calling of the
existing modules are not changed, with one exception: `--soloCellFilter
TopCells` (below). The pinned regressions give the matrices and cell calls
they gave before ([Validation](#validation)).

## Experimental add-ons for STAR-perturb

Two modes count gene expression for single-cell libraries whose cell barcode
is not one whitelist segment followed by a UMI:

- `--soloSplitseqMode Fixed` decodes two to four barcode rounds at fixed
  positions on the barcode read, each against its own allowlist (exact match,
  then one substitution), and a UMI, then uses STAR's alignment and UMI
  counting. Runbook: [SPLITSEQ_RNA.md](SPLITSEQ_RNA.md).
- `--soloPipseqMode Illumina` decodes the four-tier cell barcode and the
  three-base binning index of the newer PIP chemistry, which has no UMI,
  trims the cDNA read, and estimates transcripts per cell and gene from
  binning indices and fragment ends, following the published description of
  the chemistry. Runbook: [PIPSEQ_RNA.md](PIPSEQ_RNA.md).

Both default to `GeneFull` and to the existing `EmptyDrops_CR` caller, and
write `Solo.out/<feature>/raw` and `filtered` matrices with a metrics table.

They are experimental. They count RNA only: guide capture, `--pfMultiConfig`
and the Cell Ranger-compatible outputs are not connected to them. Alignment,
gene assignment and cell calling are STAR's own; neither mode reproduces a
kit vendor's pipeline or its cell-selection rule, and the caller is not tuned
for these assays. Each mode has been checked on few public libraries; the
runbooks give the figures with their sampling and limits. Option names,
defaults and outputs can change in a later release. With both modes off,
nothing changes for other assays.

With the modes come:

- **Barcode lists**, packaged under `share/star-suite/barcode_lists/`:
  split-pool lists for the v3 chemistry (96 and 48 round-1 wells), written
  from the on-list and replacement list of the kallisto-bustools
  documentation (BSD-2-Clause), and a PIP tier set read off the barcode reads
  of one public library, which is not a manufacturer's list. Each folder has a
  `PROVENANCE.md` with the source, the checksums and what was checked.
- **Tools**: `scripts/make_splitseq_lists.py` (allowlists, well crosswalk and
  well labels for a run's barcodes), `scripts/make_pipseq_tier_lists.py`
  (tier lists from the barcode reads of a library) and
  `scripts/make_pipseq_reference_gtf.py` (the annotation meant for the PIP
  mode, from a GENCODE GTF).
- **Launchpad recipes** `star_splitseq_rna` and `star_pipseq_rna`, one STAR
  command each. The barcode layout, the lists and the strand have no default,
  so they are chosen for the library.

## Cell filtering of a raw matrix

- **`--runMode soloCellFiltering` left out the last barcode.** The mode reads
  a raw matrix and calls cells on it. Its loader did not offer the last
  barcode that has counts to the caller. With a rank-based caller one barcode
  fewer can move the thresholds, so more than one call could differ from the
  call made during the alignment run. On the two libraries where this was
  measured, a re-call gave 16,546 cells where the run had called 16,526, and
  5,331 where the run had called 5,421. With this release a re-call of the
  same raw matrix gives the run's cells. Calls made during an alignment run
  never went through this loader and are unchanged.
- **Earlier figures made with re-calls.** The PBMC 10K and iPSC2 cell counts
  in the [v1.9.2 release notes](RELEASE_NOTES_v1.9.2.md) and the counts in
  [MSK_AMBIENT_RESCUE_ISOLATION_20260402.md](MSK_AMBIENT_RESCUE_ISOLATION_20260402.md)
  came from re-calls. Both documents now carry a dated note. The figures were
  not recomputed.
- **`--soloCellFilter TopCells N` reported N + 1 barcodes.** It now reports
  the N barcodes with the highest counts, plus any tied with the N-th. This
  applies during an alignment run and in a re-call. The number has to be
  digits only, from 1 to 4294967295; a value such as `1e3` or `2junk` was
  read by its numeric prefix before and is now a parameter error.
- **EmptyDrops backend parity test.** `tests/run_emptydrops_parity.sh`
  compared the default libscrna backend with itself. It now runs the legacy
  backend in its first arm, checks each arm's log for the backend that ran
  and uses the same knee in both. On the A375 matrix the legacy backend calls
  1,167 cells and the default backend 1,168; one barcode differs, within the
  test's tolerance of one. Neither backend is changed, and the libscrna
  backend stays the default.

## Tests added

- `tests/run_solo_cell_filtering_contract.sh`: cell filtering of generated
  raw matrices (the last barcode, `TopCells` counts and ties, malformed and
  boundary numbers, a matrix with one barcode).
- `tests/run_perturb_addons_contract.sh`: both add-on modes from FASTQ to
  matrices on generated inputs, the PIP edge cases, and the list and
  reference tools, including the packaged split-pool lists against their
  source files.
- `make -C core/legacy/source test_SplitSeq test_PipSeq`, now part of
  `make test` and of the container build.
- `mcp_server/tests/test_star_splitseq_pipseq_recipes.py`: the two recipes,
  rendered and run on generated reads with the packaged lists.

The two contract scripts are rows of the production module regression
manifest and Tier A container smokes.

## Validation

The checks below ran on one host with the binary built from commit `ba470a1`
of `dev-release-v1.13.0`. The commits after it change the link rule of one
test harness and documents.

- **Unit tests.** All 17 targets of `make test` in `core/legacy/source` pass,
  the two new ones among them.
- **Production module regression manifest** (27 rows).
  - 23 rows pass. Among them are the pinned PBMC 100K gene-expression
    regression (raw and filtered `Gene` and `GeneFull` in three run
    configurations), the UCSF 100K perturb-seq smoke, the feature-assignment
    and Flex hash-screen rows, the SLAM pre-NTR parity and both new contract
    rows.
  - Two of the 23 failed at first and pass after a fix made in this release:
    the link rule of a CBQ test harness lists its objects by hand and lacked
    the new `PipSeq.o`. The STAR binary was not involved.
  - The Flex half-probe 100K row reports differences from its preserved
    v1.9.5.a reference: the OrdMag bootstrap diagnostics of four samples and
    one column of `flexfilter_summary.tsv`. Matrices and cell lists equal the
    reference. The released v1.12.2 binary gives the same 121 output checksums
    as this release on that fixture, so 1.13.0 does not change this output;
    the stored reference predates the bootstrap-stream change of v1.12.1 and
    was not refreshed here.
  - Three rows did not run because their inputs are not on the validation
    host: the OCM composite smoke, the SLAM paired-end 100K smoke and the
    paired-end SLAM CBQ divergence harness.
- **EmptyDrops backend parity** on the A375 matrix: pass (1,167 and 1,168
  cells, one barcode apart).
- **Server and recipe tests.** 643 pass and 1 is skipped; the official
  snapshot check passes.
- **Both add-on modes at full depth with the packaged lists.** The split-pool
  library (GEO GSE352893, 97,458,666 read pairs) and the PIP library (GEO
  GSE295097, 180,144,123 read pairs) give raw and filtered matrices equal byte
  for byte to the runs made with the feature-branch build before the merge.
  The split-pool metrics table and the PIP group table are equal as well. A
  re-call of the split-pool raw matrix with `--runMode soloCellFiltering`
  gives the 16,526 cells of the run.
- **Container.** The image builds from this tree, with `test_PipSeq` and
  `test_SplitSeq` in the build. All 16 Tier A smokes pass, the two new ones
  among them, and the image carries `share/star-suite/barcode_lists/`.
- **Partial builds.** The eight CI cases and `host-api` build from fresh
  source exports.

The release workflow checks the tarballs, the installer bundle and the Debian
packages on the tag.

`STAR --version` reports `1.13.0`. Debian source packaging uses `1.13.0-1`;
Ubuntu packages use `1.13.0-1~ubuntu22.04.1` and `1.13.0-1~ubuntu24.04.1`.
Upstream STAR remains `2.7.11b`, genome-index compatibility remains `2.7.4a`
and legacy compatibility remains `2.7.1a`. Existing indexes need no rebuild.
