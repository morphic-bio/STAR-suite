# STAR libchromap Contract

Status: initial implementation + STAR orchestration hook (opt-in)

The STAR-side Chromap integration uses a direct contract instead of a process
sidecar. The contract lives in:

- `core/features/libchromap_contract/include/star_chromap_contract.h`
- `core/features/libchromap_contract/src/star_chromap_contract.cpp`

Orchestration (compile-time optional + runtime opt-in):

- **Build:** Default `make core` uses a stub (`star_chromap_orchestration_stub.cpp`):
  no `libchromap`, no extra `libhts.so` at link time. Chromap integration requires
  `make core WITH_CHROMAP=1` (same flag for `STAR` / `STARstatic` / `gdb` in
  `core/legacy/source`).
- **Runtime:** With a Chromap-enabled binary, `star_chromap_orchestration.cpp`
  calls `star::multiome::runChromapAtac` either after STAR mapping completes
  (`--chromapAtacStartMode postMapping`, default) or in a background thread
  before STAR mapping (`--chromapAtacStartMode concurrent`). Concurrent mode
  joins before `Log.final.out` / “finished successfully”.
- **Preflight:** When `--chromapAtacEnable 1`, required Chromap inputs, batch
  compatibility, thread count, and Tn5 shift mode are checked before STAR starts
  expensive mapping work.
- GEX/mapping behavior is unchanged when `--chromapAtacEnable` is absent or `0`.
- A stub-only binary with `--chromapAtacEnable 1` **fails** with a clear rebuild
  message (no silent skip).
- When `--chromapAtacEnable 1`, Chromap failures call `exitWithError` with
  `EXIT_CODE_RUNTIME` so the run does not report success.
- **Batch mode** (`--batchMode` / SLAM batch) with `--chromapAtacEnable 1` **fails**
  fast (Chromap is not supported in batch in this gate).
- Required `--chromapAtac*` paths must not be `-` / empty / `none`; optional
  fields (`summary`, `temp`, `translate`, `chromapAtacSecondaryFragments`) may be `-`.
- **Dual ATAC output (BAM/CRAM + fragments):** When `chromapAtacOutputFormat` is
  `BAM` or `CRAM`, set `chromapAtacOutputFragments` to the primary BAM/CRAM path
  and `chromapAtacSecondaryFragments` to the fragments file (e.g. `.tsv.gz`).
  This maps to Chromap `--atac-fragments`. **Intended invariants:** retained
  fragment lines match a fragment-only BED run with the same mapping options;
  BAM record count is **exactly 2 ×** fragment rows (fragment-level BAM, not
  the read-level `--BAM` path). Do **not** expect dual BAM to match a standalone
  BAM-only row count.
- **`--chromapAtacTn5ShiftMode`** must be exactly `classical` or `symmetric`
  (case-insensitive); other values error.
- **Threads:** `--chromapAtacThreads` and `--chromapAtacHtsThreads` are passed
  straight through to the contract. They are **not** coordinated with STAR
  `runThreadN` or dynamic permits (fixed thread counts for this integration
  gate), so concurrent runs should split the machine thread budget explicitly.
- **MACS3 FRAG peaks:** Optional and disabled by default. Peak calling uses the
  same libchromap/Chromap C++ path as `chromap_callpeaks`; when enabled, provide
  narrowPeak and summits outputs. File source re-reads
  `chromapAtacSecondaryFragments`; memory source accumulates fragment events
  during the Chromap mapping pass.

CLI / parameters file (all names also work in a parameters file):

| Parameter | Role |
|-----------|------|
| `chromapAtacEnable` | `0` (default) off; `1` run Chromap (non-batch only). |
| `chromapAtacStartMode` | `postMapping` (default) or `concurrent`. |
| `chromapAtacReferenceFasta` | Reference FASTA. |
| `chromapAtacIndex` | Chromap index path. |
| `chromapAtacRead1` | Comma-separated ATAC R1 FASTQs (same lane order as R2/barcode). |
| `chromapAtacRead2` | Comma-separated ATAC R2 FASTQs. |
| `chromapAtacBarcode` | Comma-separated ATAC barcode FASTQs. |
| `chromapAtacReadFormat` | Optional Chromap read-format string. Use `bc:8:23:-` for ARC barcode reads where 1-based bases 9-24 are reverse-complemented. |
| `chromapAtacBarcodeWhitelist` | Whitelist file. |
| `chromapAtacBarcodeTranslate` | Optional translation table; `-` to omit. |
| `chromapAtacOutputFragments` | Chromap **primary** output path: fragment file when format is BED/TagAlign, or SAM/BAM/CRAM path when format is SAM/BAM/CRAM. |
| `chromapAtacSecondaryFragments` | Optional **secondary** scATAC fragments path when `chromapAtacOutputFormat` is `BAM` or `CRAM` (dual output in one Chromap pass). Must differ from `chromapAtacOutputFragments`; `-` to omit. |
| `chromapAtacOutputFormat` | `BED`, `fragments`, `TagAlign`, `SAM`, `BAM`, `CRAM`, or `pairs`. |
| `chromapAtacSummary` | Optional summary path; `-` to omit. |
| `chromapAtacTempDir` | Optional temp dir; `-` for Chromap default. |
| `chromapAtacThreads` | Chromap thread count (default `1`). |
| `chromapAtacHtsThreads` | HTS compression threads for BAM/CRAM output. |
| `chromapAtacSortBam` | Sort BAM/CRAM output when `1`. |
| `chromapAtacWriteIndex` | Write BAM/CRAM index when `1`; requires sorted BAM/CRAM output. |
| `chromapAtacSortBamRam` | Chromap BAM/CRAM sort RAM limit in bytes. |
| `chromapAtacTn5ShiftMode` | `classical` or `symmetric`. |
| `chromapAtacCallMacs3FragPeaks` | `0` (default) off; `1` run MACS3-compatible FRAG peak calling after Chromap ATAC mapping. |
| `chromapAtacMacs3FragPeaksOutput` | narrowPeak output path, required when peak calling is enabled. |
| `chromapAtacMacs3FragSummitsOutput` | summits output path, required when peak calling is enabled. |
| `chromapAtacMacs3FragPvalue` | MACS3 FRAG p-value threshold; default `1e-5`. |
| `chromapAtacMacs3FragMinLength` | minimum peak length; default `200`. |
| `chromapAtacMacs3FragMaxGap` | max gap for peak merging; default `30`. |
| `chromapAtacMacs3FragPeaksSource` | `file` (default) re-read fragments, or `memory` use in-memory fragment events from the mapping pass. |
| `chromapAtacMacs3FragKeepIntermediates` | optional directory to keep intermediate bedGraph-style outputs; `-` to omit. |
| `chromapAtacMacs3FragUint8Counts` | `1` (default) use MACS3 uint8 count semantics; `0` disables. |

Build:

```bash
make star-libchromap-contract          # contract runner + libstar_chromap_contract.a
make core                              # default: portable STAR (no Chromap link)
make core WITH_CHROMAP=1               # STAR + libchromap + system libhts for Chromap symbols
```

With `WITH_CHROMAP=1`, the legacy Makefile builds `libchromap.a` in
`Chromap-suite` if needed. Override paths:

- `CHROMAP_SUITE_DIR` — default `/mnt/pikachu/Chromap-suite`
- `CHROMAP_SYS_HTS` — shared `libhts` passed **after** `libchromap.a` on the link
  line (path to `.so`, e.g. `/lib/x86_64-linux-gnu/libhts.so`). Required because
  Chromap may need a newer hts ABI than STAR’s vendored static archive. Override
  on non-Debian layouts.

Chromap-enabled builds are currently supported only for `STAR`, `STARstatic`, and
`gdb` targets. `STARlong`, `STARlongStatic`, `POSIXSHARED`, `gdb-long`, and Mac
static targets fail fast with `WITH_CHROMAP=1` rather than attempting an
incomplete link. Switching between `WITH_CHROMAP=0` and `WITH_CHROMAP=1`
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
5. Build `make core WITH_CHROMAP=1` and re-run the same fixture through STAR with
   `--chromapAtacEnable 1` and matching `--chromapAtac*` paths; compare fragments
   again.
6. Repeat the STAR run with `--chromapAtacCallMacs3FragPeaks 1` to validate peak
   outputs from STAR orchestration.
7. Repeat the STAR run with `--chromapAtacStartMode concurrent` and a split
   thread budget to validate same-process simultaneous execution.
8. Run the end-to-end benchmark gate on the 100K multiome fixture, then the full
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
--chromapAtacSecondaryFragments ./atac_fragments.tsv.gz
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
  --atac-fragments ./atac_fragments.tsv.gz --threads 8
```

Peak-calling extension:

```bash
star_libchromap_contract_runner \
  --ref genome.fa --index genome.index \
  --read1 'R1_L001.fastq.gz,...' --read2 'R3_L001.fastq.gz,...' --barcode 'R2_L001.fastq.gz,...' \
  --barcode-whitelist whitelist.txt \
  --read-format 'bc:8:23:-' \
  --output ./atac_possorted.bam --output-format BAM --sort-bam \
  --atac-fragments ./atac_fragments.tsv.gz --threads 8 \
  --call-macs3-frag-peaks \
  --macs3-frag-peaks-output ./atac_peaks.narrowPeak \
  --macs3-frag-summits-output ./atac_summits.bed \
  --macs3-frag-peaks-source memory
```

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
