# STAR libchromap Contract

Status: implemented STAR orchestration hook. Current production
multiome runs use STAR/Chromap for local BAM plus binary sidecar
materialization, then call peaks and build ATAC peak MEX from that sidecar
outside STAR.

The STAR-side Chromap integration uses a direct contract instead of a process
sidecar. The contract lives in:

- `core/features/libchromap_contract/include/star_chromap_contract.h`
- `core/features/libchromap_contract/src/star_chromap_contract.cpp`

Orchestration (default compiled in + runtime opt-in):

- **Build:** Default `make core` builds the Chromap-enabled multiome-capable
  STAR binary and links `libchromap`, `libMACS3`, the STAR libchromap contract,
  and system `libhts`. Use `make core-portable` or `make core WITH_CHROMAP=0`
  only for explicit no-Chromap compatibility builds. The stub
  (`star_chromap_orchestration_stub.cpp`) remains available for that portable
  mode.
- **Runtime:** With a Chromap-enabled binary, `star_chromap_orchestration.cpp`
  calls `star::multiome::runChromapAtac` either after STAR mapping completes
  (`--chromapAtacStartMode postMapping`, default) or in a background thread
  before STAR mapping (`--chromapAtacStartMode concurrent`). Concurrent mode
  joins before `Log.final.out` / “finished successfully”.
- **Preflight:** When `--chromapAtacEnable 1`, required Chromap inputs, batch
  compatibility, thread count, and Tn5 shift mode are checked before STAR starts
  expensive mapping work.
- GEX/mapping behavior is unchanged when `--chromapAtacEnable` is absent or `0`.
- A portable/stub-only binary with `--chromapAtacEnable 1` **fails** with a
  clear rebuild message (no silent skip).
- When `--chromapAtacEnable 1`, Chromap failures call `exitWithError` with
  `EXIT_CODE_RUNTIME` so the run does not report success.
- **Batch mode** (`--batchMode` / SLAM batch) with `--chromapAtacEnable 1` **fails**
  fast (Chromap is not supported in batch in this gate).
- Required `--chromapAtac*` paths must not be `-` / empty / `none`; optional
  fields (`summary`, `temp`, `translate`, `chromapAtacSecondaryFragments`) may be `-`.
- **Input format:** `chromapAtacInputFormat fastq` is the default and preserves
  the existing `chromapAtacRead1` / `chromapAtacRead2` / `chromapAtacBarcode`
  surface. `chromapAtacInputFormat cbq` routes STAR's libchromap contract to
  Chromap-suite's native CBQ reader using `chromapAtacReadPairCbq` and
  `chromapAtacBarcodeCbq`; it does not materialize FASTQ files. The paired-read
  and barcode CBQ lane lists must have the same length and stay record-aligned.
- **Dual ATAC output (BAM/CRAM + fragments):** When `chromapAtacOutputFormat` is
  `BAM` or `CRAM`, set `chromapAtacOutputFragments` to the primary BAM/CRAM path
  and `chromapAtacSecondaryFragments` to the scATAC fragment sidecar. Production
  uses a binary sidecar path such as `atac_fragments.bin` plus the generated
  `atac_fragments.bin.chroms.tsv`. This maps to Chromap `--atac-fragments`.
  **Intended invariants:** retained fragment events match a fragment-only run
  with the same mapping options; BAM record count is **exactly 2 x** fragment
  rows (fragment-level BAM, not the read-level `--BAM` path). Do **not** expect
  dual BAM to match a standalone BAM-only row count.
- **ATAC Y/noY BAM streams:** When `chromapAtacOutputFormat` is `SAM`, `BAM`,
  or `CRAM`, `chromapAtacEmitNoYBam 1` and `chromapAtacEmitYBam 1` request
  additional Chromap streams from the same libchromap mapping run. Explicit
  `chromapAtacNoYOutput` / `chromapAtacYOutput` paths are accepted; if omitted,
  STAR derives Chromap-style `.noY` and `.Y` paths from
  `chromapAtacOutputFragments`.
- **`--chromapAtacTn5ShiftMode`** must be exactly `classical` or `symmetric`
  (case-insensitive); other values error.
- **Threads:** `--chromapAtacThreads` and `--chromapAtacHtsThreads` are passed
  straight through to the contract. They are **not** coordinated with STAR
  `runThreadN` or dynamic permits (fixed thread counts for this integration
  gate), so concurrent runs should split the machine thread budget explicitly.
- **MACS3 FRAG peaks:** Optional and disabled by default in STAR. The current
  production multiome wrapper does not run peak calling inside STAR; it writes
  BAM plus the binary sidecar locally and then invokes
  `core/features/libchromap_contract/star_multiome_atac_peak_mex` to call peaks
  from the sidecar and build the ATAC peak MEX. This avoids the legacy
  file-source spill/re-read path.

CLI / parameters file (all names also work in a parameters file):

| Parameter | Role |
|-----------|------|
| `chromapAtacEnable` | `0` (default) off; `1` run Chromap (non-batch only). |
| `chromapAtacStartMode` | `postMapping` (default) or `concurrent`. |
| `chromapAtacReferenceFasta` | Reference FASTA. |
| `chromapAtacIndex` | Chromap index path. |
| `chromapAtacInputFormat` | `fastq` (default) or `cbq`. |
| `chromapAtacRead1` | Comma-separated ATAC R1 FASTQs (same lane order as R2/barcode). |
| `chromapAtacRead2` | Comma-separated ATAC R2 FASTQs. |
| `chromapAtacBarcode` | Comma-separated ATAC barcode FASTQs. |
| `chromapAtacReadPairCbq` | Comma-separated paired-read ATAC CBQs when `chromapAtacInputFormat=cbq`. |
| `chromapAtacBarcodeCbq` | Comma-separated barcode ATAC CBQs when `chromapAtacInputFormat=cbq`. |
| `chromapAtacReadFormat` | Optional Chromap read-format string. Use `bc:8:23:-` for ARC barcode reads where 1-based bases 9-24 are reverse-complemented. |
| `chromapAtacBarcodeWhitelist` | Whitelist file. |
| `chromapAtacBarcodeTranslate` | Optional translation table; `-` to omit. |
| `chromapAtacBarcodeTranslateFromFirst` | `1` means treat the first column of the translation table as the source barcode namespace. Required for the current ARC ATAC-to-GEX table. |
| `chromapAtacOutputFragments` | Chromap **primary** output path: fragment file when format is BED/TagAlign, or SAM/BAM/CRAM path when format is SAM/BAM/CRAM. |
| `chromapAtacSecondaryFragments` | Optional **secondary** scATAC fragment sidecar path when `chromapAtacOutputFormat` is `BAM` or `CRAM` (dual output in one Chromap pass). Production uses `atac_fragments.bin`, which also writes `atac_fragments.bin.chroms.tsv`. Must differ from `chromapAtacOutputFragments`; `-` to omit. |
| `chromapAtacOutputFormat` | `BED`, `fragments`, `TagAlign`, `SAM`, `BAM`, `CRAM`, or `pairs`. |
| `chromapAtacSummary` | Optional summary path; `-` to omit. |
| `chromapAtacTempDir` | Optional temp dir; `-` for Chromap default. |
| `chromapAtacThreads` | Chromap thread count (default `1`). |
| `chromapAtacHtsThreads` | HTS compression threads for BAM/CRAM output. |
| `chromapAtacSortBam` | Sort BAM/CRAM output when `1`. |
| `chromapAtacWriteIndex` | Write BAM/CRAM index when `1`; requires sorted BAM/CRAM output. |
| `chromapAtacSortBamRam` | Chromap BAM/CRAM sort RAM limit in bytes. |
| `chromapAtacEmitNoYBam` | Emit an additional SAM/BAM/CRAM stream excluding reads with Y-chromosome alignments when `1`. |
| `chromapAtacEmitYBam` | Emit an additional SAM/BAM/CRAM stream containing reads with Y-chromosome alignments when `1`. |
| `chromapAtacNoYOutput` | Optional explicit noY stream path; `-` derives from `chromapAtacOutputFragments`. |
| `chromapAtacYOutput` | Optional explicit Y stream path; `-` derives from `chromapAtacOutputFragments`. |
| `chromapAtacLowMem` | `0` (default) off; `1` enables Chromap low-memory overflow-spill mode. Production JAX multiome uses `1`. |
| `chromapAtacLowMemRam` | RAM threshold in bytes for low-memory spill; `0` uses Chromap defaults. |
| `chromapAtacMacs3FragLowMem` | `1` uses the low-memory libMACS3 fragment workspace for in-process peak calling. Production keeps this enabled in wrappers but performs peak/MEX construction from the sidecar outside STAR. |
| `chromapAtacTn5ShiftMode` | `classical` or `symmetric`. |
| `chromapAtacCallMacs3FragPeaks` | `0` (default) off; `1` run MACS3-compatible FRAG peak calling after Chromap ATAC mapping. |
| `chromapAtacMacs3FragPeaksOutput` | narrowPeak output path, required when peak calling is enabled. |
| `chromapAtacMacs3FragSummitsOutput` | summits output path, required when peak calling is enabled. |
| `chromapAtacMacs3FragPvalue` | MACS3 FRAG p-value threshold; default `1e-5`. |
| `chromapAtacMacs3FragQvalue` | MACS3 FRAG q-value/FDR threshold; default `0` disables q-value mode. A positive value in `(0,1]` uses q-value mode and rejects an explicit `chromapAtacMacs3FragPvalue`. |
| `chromapAtacMacs3FragMinLength` | minimum peak length; default `200`. |
| `chromapAtacMacs3FragMaxGap` | max gap for peak merging; default `30`. |
| `chromapAtacMacs3FragKeepIntermediates` | optional directory to keep intermediate bedGraph-style outputs; `-` to omit. |
| `chromapAtacMacs3FragUint8Counts` | `1` (default) use MACS3 uint8 count semantics; `0` disables. |

The STAR CLI no longer exposes `chromapAtacMacs3FragPeaksSource`. Production
multiome should use the post-materialization boundary: write BAM plus the
binary sidecar locally, then call peaks and build ATAC peak MEX from that
sidecar outside STAR. Do not route production through a gzipped BED/fragments
spill file that will be read again for peak/MEX construction.

Build:

```bash
make star-libchromap-contract          # contract runner + libstar_chromap_contract.a
make core                              # default: STAR + libchromap + system libhts
make core-portable                     # explicit portable STAR (no Chromap link)
make core WITH_CHROMAP=0               # same no-Chromap compatibility mode
```

With the default Chromap-enabled build, the legacy Makefile builds
`libchromap.a` in `Chromap-suite` if needed. Override paths:

- `CHROMAP_SUITE_DIR` — default `/mnt/pikachu/Chromap-suite`
- `CHROMAP_SYS_HTS` — shared `libhts` passed **after** `libchromap.a` on the link
  line (path to `.so`, e.g. `/lib/x86_64-linux-gnu/libhts.so`). Required because
  Chromap may need a newer hts ABI than STAR’s vendored static archive. Override
  on non-Debian layouts.

Chromap-enabled builds are currently supported only for `STAR`, `STARstatic`, and
`gdb` targets. Top-level `make core-long` forces `WITH_CHROMAP=0`; direct
`STARlong`, `STARlongStatic`, `POSIXSHARED`, `gdb-long`, and Mac static targets
fail fast with `WITH_CHROMAP=1` rather than attempting an incomplete link.
Switching between `WITH_CHROMAP=0` and `WITH_CHROMAP=1`
regenerates dependency state automatically.

The contract implementation links `Chromap-suite` and calls
`chromap::RunAtacMapping`. Chromap types do not appear in `Parameters.h`; only
the adapter header is included from orchestration code.

The permit hooks are part of the STAR contract but are intentionally rejected at
runtime until Chromap exposes batch-level permit acquire/release points. This
keeps the ABI shape in place without pretending permits are already enforced.

Validation order:

1. Build `make star-libchromap-contract`.
2. Run `star_libchromap_contract_runner` on the 100K ATAC fixture.
3. Compare its fragments against the existing Chromap CLI/libchromap parity
   outputs with `scripts/compare_fragment_tuples.sh` from the coordination repo.
4. Repeat the contract runner with `--call-macs3-frag-peaks` plus narrowPeak and
   summits paths; compare peaks against a Chromap CLI/libchromap reference.
5. Build `make core` and re-run the same fixture through STAR with
   `--chromapAtacEnable 1` and matching `--chromapAtac*` paths; compare fragments
   again.
6. Build peaks and ATAC peak MEX from the binary sidecar with
   `star_multiome_atac_peak_mex`; this is the production boundary.
7. Repeat the STAR run with `--chromapAtacStartMode concurrent`, low-memory
   Chromap enabled, and a split thread budget to validate same-process
   simultaneous execution.
8. Treat STAR `--chromapAtacCallMacs3FragPeaks 1` as legacy/integration
   validation only. Do not use it as the production multiome path unless the
   runbook is intentionally revised.
9. Run the end-to-end benchmark gate on the 100K multiome fixture, then the full
   10x demo dataset, before adding dynamic permit sharing.

Initial validation:

- Run root: `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/star_libchromap_contract_20260424_045451`
- Comparator: `/mnt/pikachu/multiomic-atac-scrna/scripts/compare_fragment_tuples.sh`
- Expected source: validated Chromap CLI output from `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/libchromap_parity_20260424_041454/chromap_cli/fragments.tsv`
- Expected rows: `320,017`
- Observed rows: `320,017`
- Shared rows: `320,017`
- Recall: `1.000000000000`
- Precision: `1.000000000000`
- Jaccard: `1.000000000000`
- Per-barcode fragment counts: exact match

STAR-integrated spot check (same tuples as CLI reference):

- Example output dir: `/tmp/star_chromap_integration_100k` (local; not committed)
- After a successful run, `fragment_tuple_metrics.tsv` from
  `compare_fragment_tuples.sh` should report `exact_fragment_tuple_match` true and
  `per_barcode_counts_match` true. Summary CSV differences may still appear in
  cache telemetry fields documented in the multiomic runbook.

**Example (dual fragments + sorted BAM, 100K fixture paths illustrative):**

```text
--chromapAtacEnable 1
--chromapAtacOutputFormat BAM
--chromapAtacSortBam 1
--chromapAtacOutputFragments ./atac_possorted.bam
--chromapAtacSecondaryFragments ./atac_fragments.bin
--chromapAtacEmitNoYBam 1
--chromapAtacEmitYBam 1
--chromapAtacReadFormat bc:8:23:-
--chromapAtacLowMem 1
--chromapAtacLowMemRam 0
... (reference, index, read1/read2/barcode CSVs, whitelist, threads, etc.)
```

Contract runner equivalent:

```bash
star_libchromap_contract_runner \
  --ref genome.fa --index genome.index \
  --read1 'R1_L001.fastq.gz,...' --read2 'R3_L001.fastq.gz,...' --barcode 'R2_L001.fastq.gz,...' \
  --barcode-whitelist whitelist.txt \
  --read-format 'bc:8:23:-' \
  --output ./atac_possorted.bam --output-format BAM --sort-bam \
  --atac-fragments ./atac_fragments.bin \
  --emit-noY-bam --noY-output ./atac_possorted.noY.bam \
  --emit-Y-bam --Y-output ./atac_possorted.Y.bam \
  --threads 8
```

CBQ contract-runner equivalent:

```bash
star_libchromap_contract_runner \
  --ref genome.fa --index genome.index \
  --input-format cbq \
  --read-pair-cbq 'ATAC_L001.reads.cbq,...' \
  --barcode-cbq 'ATAC_L001.barcodes.cbq,...' \
  --barcode-whitelist whitelist.txt \
  --read-format 'bc:8:23:-' \
  --output ./atac_possorted.bam --output-format BAM --sort-bam \
  --atac-fragments ./atac_fragments.bin \
  --emit-noY-bam --noY-output ./atac_possorted.noY.bam \
  --emit-Y-bam --Y-output ./atac_possorted.Y.bam \
  --threads 8
```

Peak-calling extension:

```bash
star_multiome_atac_peak_mex \
  --sidecar ./atac_fragments.bin \
  --barcode-translate ./atac2gex.tsv \
  --barcode-translate-from-first \
  --call-peaks-from-sidecar \
  --peaks ./atac_peaks.narrowPeak \
  --summits-out ./atac_summits.bed \
  --out-dir ./atac/peak_mex \
  --metrics-tsv ./atac/atac_metrics.tsv \
  --threads 16 \
  --macs3-frag-qvalue 0.05 \
  --temp-dir ./chromap_tmp
```

Use either `--macs3-frag-pvalue` or `--macs3-frag-qvalue` with
`star_multiome_atac_peak_mex`; the flags are mutually exclusive. The STAR
parameter equivalent is `--chromapAtacMacs3FragQvalue 0.05` for the in-process
Chromap/libMACS3 path or inline peak/MEX path.

Concurrent BAM smoke from a clean Chromap-enabled build (100K multiome fixture):

- Run root: `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/star_chromap_concurrent_bam_final_20260424_162434`
- Build: `make -C /mnt/pikachu/STAR-suite/core/legacy/source STAR WITH_CHROMAP=1`
- Launch mode: `--chromapAtacStartMode concurrent`
- Thread budget: STAR `--runThreadN 8`; Chromap `--chromapAtacThreads 8`;
  Chromap HTS `--chromapAtacHtsThreads 2`
- Wall time: `31.20 s`
- Max RSS: `21,780,764 KB`
- Output checks: GEX and ATAC BAMs pass `samtools quickcheck`
- ATAC records: `640,088`
- ATAC record MD5: `5d9a2f771d821e5f203bf3b2a28ab805`
- ATAC idxstats MD5: `e1072493d751ead3ac0f6808c349a086`
- GEX records: `246,287`
- GEX record MD5: `9d6bd9242cd27c64978619de550d8edb`
- Log ordering confirms Chromap starts before STAR mapping and STAR waits for it
  before final success.
