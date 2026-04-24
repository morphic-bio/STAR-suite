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
  calls `star::multiome::runChromapAtac` after STAR mapping completes (after
  optional wiggle output, before `Log.final.out` / “finished successfully”).
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
  fields (`summary`, `temp`, `translate`) may be `-`.
- **`--chromapAtacTn5ShiftMode`** must be exactly `classical` or `symmetric`
  (case-insensitive); other values error.
- **Threads:** `--chromapAtacThreads` is passed straight through to the contract.
  It is **not** coordinated with STAR `runThreadN` or dynamic permits (fixed
  thread count for this integration gate).

CLI / parameters file (all names also work in a parameters file):

| Parameter | Role |
|-----------|------|
| `chromapAtacEnable` | `0` (default) off; `1` run Chromap after mapping (non-batch only). |
| `chromapAtacReferenceFasta` | Reference FASTA. |
| `chromapAtacIndex` | Chromap index path. |
| `chromapAtacRead1` | Comma-separated ATAC R1 FASTQs (same lane order as R2/barcode). |
| `chromapAtacRead2` | Comma-separated ATAC R2 FASTQs. |
| `chromapAtacBarcode` | Comma-separated ATAC barcode FASTQs. |
| `chromapAtacBarcodeWhitelist` | Whitelist file. |
| `chromapAtacBarcodeTranslate` | Optional translation table; `-` to omit. |
| `chromapAtacOutputFragments` | Fragments output path. |
| `chromapAtacSummary` | Optional summary path; `-` to omit. |
| `chromapAtacTempDir` | Optional temp dir; `-` for Chromap default. |
| `chromapAtacThreads` | Chromap thread count (default `1`). |
| `chromapAtacTn5ShiftMode` | `classical` or `symmetric`. |

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
4. Build `make core WITH_CHROMAP=1` and re-run the same fixture through STAR with
   `--chromapAtacEnable 1` and matching `--chromapAtac*` paths; compare fragments
   again.
5. Run the end-to-end benchmark gate on the 100K multiome fixture, then the full
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
