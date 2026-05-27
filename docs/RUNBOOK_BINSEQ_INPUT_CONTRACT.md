# Runbook: BINSEQ Input Contract

Date: 2026-05-27

This runbook sketches the path for making BINSEQ an input option beside the
current FASTA/FASTQ input regime. The first implementation milestone is not
BINSEQ itself; it is a read-input contract that makes the current FASTQ path a
module behind a stable interface. BINSEQ should then plug into the same
contract and be validated with the same harness.

## External Format Notes

ARC Institute `binseq` is a Rust library for binary sequencing files. It
defines three variants:

- `*.bq`: fixed-length records without quality scores.
- `*.vbq`: variable-length records with optional quality scores and headers.
- `*.cbq`: columnar variable-length records with optional quality scores and
  headers.

For STAR-suite integration, assume `*.cbq` is the target variant because it is
the current recommended, lossless, flexible format and natively supports `N`
sequence content. Treat `*.bq` as unsuitable for general STAR input because it
does not carry qualities and is lossy for ambiguous bases. Treat `*.vbq` as
legacy-compatible but lower priority.

`binseq` is a library, not a primary command-line interface. The companion
`bqtools` CLI can encode FASTQ/FASTA/SAM/BAM/CRAM to BINSEQ, decode back to
FASTQ/FASTA/TSV, inspect files with `info --json`, and create named pipes for
legacy tools. Use `bqtools` for fixture creation and oracle checks during early
contract work; do not make the production STAR path depend permanently on
decode-to-FASTQ unless the user explicitly chooses a compatibility mode.

References:

- `https://github.com/arcinstitute/binseq`
- `https://docs.rs/binseq/latest/binseq/`
- `https://github.com/arcinstitute/bqtools`

## Current STAR-suite Input Surface

Current input is centered on `--readFilesType`, `--readFilesIn`,
`--readFilesManifest`, `--readFilesPrefix`, and `--readFilesCommand`.

Relevant files:

- `core/legacy/source/Parameters_readFilesInit.cpp`
- `core/legacy/source/Parameters_openReadsFiles.cpp`
- `core/legacy/source/Parameters_closeReadsFiles.cpp`
- `core/legacy/source/readLoad.cpp`
- `core/legacy/source/ReadAlign_oneRead.cpp`
- `core/legacy/source/mapThreadsSpawn.cpp`

The existing FASTQ path has two separable responsibilities that are currently
interleaved:

- Source discovery and lane/mate normalization: `readFilesInit()` expands
  comma-separated lanes or a manifest into `readFilesNames[mate][lane]`, sets
  `readFilesN`, `readNends`, and read group state.
- Stream materialization: `openReadsFiles()` opens plain FASTX directly or
  creates per-mate FIFOs for `readFilesCommand` or internal gzip streaming.
  `readLoad()` then parses a FASTA/FASTQ-like text stream into per-mate buffers.

The new contract should preserve STAR CLI compatibility but move these
responsibilities behind explicit module boundaries.

## Target Contract

Introduce a format-independent read input contract before adding BINSEQ. The
contract should expose logical read records, not FASTQ lines.

Minimum record fields:

- `read_name`: canonical name after STAR read-name separator policy is applied.
- `read_name_extra`: untouched trailing header payload needed by SAM/BAM input
  and Solo barcode extraction paths.
- `lane_index`: existing `readFilesIndex` semantics for read group lookup,
  batch accounting, and SLAM resume behavior.
- `read_ordinal`: existing `iReadAll` semantics where available.
- `read_filter`: existing filter character, defaulting to `Y`.
- `mate_count`: 1, 2, or current `MAX_N_MATES` upper bound.
- For each mate: sequence, quality, original length, and a flag indicating
  FASTA/no-quality fallback.

Minimum module interface:

```text
InputModule::configure(Parameters)
InputModule::describe_sources() -> lanes, mates, read groups, counts-if-known
InputModule::open()
InputModule::next_record() -> InputRecord | EOF
InputModule::close()
```

The first C++ implementation can use concrete classes rather than a broad
framework, but the call site should stop assuming "input equals `istream`
containing FASTQ text". Keep the old `readLoad()` parser as the FASTX module's
implementation detail.

## CLI Policy

Keep existing CLI behavior working unchanged:

```bash
--readFilesType Fastx
--readFilesIn R1.fastq.gz R2.fastq.gz
--readFilesCommand -
```

Add BINSEQ through the existing `readFilesType` family rather than introducing
a second parallel input namespace:

```bash
--readFilesType Binseq PE
--readFilesIn sample.cbq
```

Planned support order:

1. Paired `*.cbq` file per lane, where one BINSEQ record carries R1/R2.
2. Single-end `*.cbq` file per lane.
3. Split-mate BINSEQ files, only if needed for real datasets.
4. `*.vbq` support after CBQ is stable.
5. `*.bq` only for explicit sequence-only experiments, not general alignment.

`--readFilesManifest` should remain the preferred multi-sample/multi-lane
surface. For BINSEQ, the manifest should accept one input path for paired CBQ:

```text
paired_cbq_file    -    ID:sample1
```

The `-` second column means "not a second external mate path", not "single-end",
when `--readFilesType Binseq PE` is active. The module must validate this
strictly and write a clear error when the manifest and BINSEQ file metadata
disagree.

## Implementation Phases

### Phase 1: Extract The Input Contract

Goal: no behavior change.

Tasks:

- Create a small contract layer under `core/legacy/source/input/` or
  `core/features/input/`.
- Move source normalization out of `Parameters_readFilesInit.cpp` into a
  reusable `InputSourcePlan`.
- Keep `readFilesNames`, `readFilesN`, `readNends`, and `outSAMattrRG` populated
  for compatibility while the rest of STAR still reads those fields directly.
- Wrap the existing FASTX stream path as `FastxInputModule`.
- Preserve internal gzip behavior, `readFilesCommand`, FILE markers, manifest
  behavior, read groups, `readMapNumber`, batch mode, and SLAM resume semantics.

Exit criteria:

- Existing FASTQ, gzip, `readFilesCommand`, manifest, Solo, Flex, SLAM, and
  batch smoke tests pass with no intentional output changes.
- No BINSEQ code is required for this phase.

### Phase 2: Build The Input Harness

Goal: make input behavior testable independent of alignment.

Tasks:

- Add a harness binary or test helper that opens an input module and writes a
  normalized record dump:

```text
lane_index	read_ordinal	read_filter	read_name	mate	mate_len	seq_sha256	qual_sha256
```

- Add a mode that writes small literal FASTQ for debugging.
- Add fixtures covering:
  - plain FASTQ;
  - gzip FASTQ using internal zlib;
  - `--readFilesCommand zcat`;
  - comma-separated lanes;
  - `--readFilesManifest`;
  - paired reads;
  - single-end reads;
  - optional barcode read / 3-mate surface if there is a compact fixture.
- Make the harness compare old FASTQ parser output against the new
  `FastxInputModule`.

Exit criteria:

- Harness proves the contract is byte-equivalent for current FASTQ workflows.
- Harness is fast enough for CI-tier smoke tests.

Implementation status on `feature/binseq-input`:

- Added `core/legacy/source/fastx_input_harness` via the
  `fastx-input-harness` Makefile target.
- The harness writes the normalized TSV dump above and supports `--dump-fastq`
  for literal FASTQ debugging.
- Smoke coverage lives in `tests/run_fastx_input_harness_smoke.sh` and covers
  plain paired FASTQ, gzip/internal zlib, explicit `zcat`, comma-separated
  lanes, manifest lanes, single-end FASTQ, and a compact three-mate FASTQ
  case.
- The smoke compares normalized dumps across plain/gzip/zcat and
  comma-lane/manifest paths. Phase 3 now exercises the same module through
  STAR's mapping chunk path.

### Phase 3: Re-implement FASTQ Behind The Contract

Goal: make FASTQ a normal module, not the implicit core input mechanism.

Tasks:

- Move FASTX parsing, gzip streaming, command FIFO handling, and close/cleanup
  into the FASTX module boundary.
- Keep external behavior and logs stable.
- Keep direct file and FIFO support only inside the module.
- Remove new call-site assumptions that a read source is an `istream`.

Exit criteria:

- All Phase 2 harness cases pass.
- Core smoke tests pass.
- A clean rebuild is performed before investigating any mismatch.

Implementation status on `feature/binseq-input`:

- `FastxInputModule` is now wired into STAR chunk ingestion for
  `--readFilesType Fastx`, while preserving `readFilesNames`, `readFilesN`,
  `readNends`, and read group compatibility fields for existing code paths.
- The runtime path consumes module `InputRecord`s directly and handles EOF
  without falling back to the legacy `istream` parser for Fastx input.
- `tests/run_fastx_contract_star_smoke.sh` covers plain FASTQ, internal gzip,
  explicit `zcat`, comma-separated lanes, manifest input, and the current
  Y/noY FASTQ gate.
- Y/noY FASTQ emission is intentionally FASTQ-only for now: when
  `--emitYNoYFastq yes` is set, every Fastx input path must end in
  `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`. TODO: extend Y-removal/Y-noY
  emission to future non-FASTQ input modules, including BINSEQ, before enabling
  those combinations.

### Phase 4: BINSEQ Probe Module

Goal: validate real BINSEQ data against the contract before STAR integration.

Tasks:

- Install or build `bqtools` in a disposable tool location for fixture work:

```bash
cargo install bqtools
```

- Create synthetic CBQ fixtures from existing tiny FASTQ fixtures:

```bash
bqtools encode R1.fastq.gz R2.fastq.gz --mode cbq -o tiny_pair.cbq -T 4
bqtools info tiny_pair.cbq --json
bqtools decode tiny_pair.cbq --prefix decoded
```

- Use decoded FASTQ as the oracle only for test comparison.
- Add a prototype `BinseqInputModule` as either:
  - a helper process that emits the contract stream from CBQ; or
  - a Rust static library/FFI adapter around the `binseq` crate.

Preferred first implementation: a helper binary built from Rust and launched by
STAR, because it isolates Cargo/Rust dependencies from the core C++ build while
the contract is still settling. Revisit direct FFI only after the harness and
real data prove this is worth the build complexity.

Exit criteria:

- Synthetic CBQ fixture produces the same normalized record dump as the source
  FASTQ fixture.
- Paired-read synchronization, names, qualities, and `N` handling are verified.
- Module failures produce actionable errors for unsupported variants or missing
  qualities.

Implementation status on `feature/binseq-input`:

- Added `BinseqInputProbeModule`, a non-production `InputModule` that shells out
  to `bqtools decode`, writes decoded FASTQ into a caller-provided scratch
  directory, and delegates record iteration to `FastxInputModule`.
- Added the `binseq-probe-harness` Makefile target and
  `core/legacy/source/binseq_probe_harness` for normalized contract dumps from
  paired CBQ files or a BINSEQ manifest row (`CBQ<TAB>-<TAB>ReadGroup`).
- Added `tests/run_binseq_probe_smoke.sh`, which creates paired synthetic FASTQ,
  encodes it to CBQ with `bqtools`, records `bqtools info`, and compares direct
  plus manifest BINSEQ probe dumps against the source FASTQ dump.
- The smoke skips cleanly when `bqtools` is unavailable. Set
  `BQTOOLS=/path/to/bqtools` to run it on hosts without `bqtools` in `PATH`.
- This is intentionally a probe-only stage: STAR CLI routing for
  `--readFilesType Binseq` remains Phase 6 work after real CBQ fixtures pass
  the harness.

### Phase 5: Find Real BINSEQ Files

Goal: test against externally produced BINSEQ, not only files generated from our
own FASTQ fixtures.

Tasks:

- Search ARC/bqtools examples, public releases, internal Morphic staging, and
  collaborator handoffs for `.cbq`, `.vbq`, and `.bq`.
- Record every candidate in `tests/ARTIFACTS.md` if it stays local and is too
  large or private for git.
- For each candidate:
  - run `bqtools info --json`;
  - decode a small prefix for spot-checking;
  - run the input harness;
  - compare record count, mate count, names, sequence hashes, quality hashes,
    and lane/read group handling.

Exit criteria:

- At least one real paired CBQ and one real single-end CBQ pass the contract
  harness.
- Any unsupported BINSEQ variants are documented with explicit errors.

### Phase 6: STAR-suite Integration

Goal: make BINSEQ a production input option.

Tasks:

- Wire `--readFilesType Binseq PE|SE` through `Parameters`.
- Ensure `readFilesManifest`, read groups, batch mode, Solo/Flex barcode paths,
  SLAM, and `readMapNumber` work through the contract.
- Add smoke tests that map a tiny CBQ fixture and compare against FASTQ output:
  - alignment summary;
  - BAM read names/tags where relevant;
  - Solo raw/filtered MEX checks for a tiny scRNA fixture if available.
- Add documentation to `parametersDefault` and regenerate
  `parametersDefault.xxd`.

Exit criteria:

- `make core` passes.
- `make core WITH_CHROMAP=1` passes.
- Input harness tests pass.
- At least one STAR smoke test using `--readFilesType Binseq PE` passes.

## Risk Register

- Rust dependency in STAR core: keep it out of the default C++ build until
  contract and fixture behavior are proven.
- Three-mate workflows: BINSEQ paired records cover common R1/R2 cases; barcode
  reads may require split sources or an extended adapter plan.
- Quality semantics: BQ lacks qualities, so it must not silently substitute
  qualities for production alignment.
- Header semantics: STAR currently stores useful metadata in FASTQ header
  extras and uses `readFilesIndex`; the contract must preserve those semantics
  explicitly.
- Performance: `bqtools decode` is fine for oracle tests, but production should
  avoid materializing full FASTQ unless explicitly requested.

## Recommended First PR

The first PR on `feature/binseq-input` should contain only the input contract
scaffold and FASTX module wrapping:

- no BINSEQ dependency;
- no CLI behavior change;
- harness output for a tiny FASTQ fixture;
- focused tests proving current FASTQ behavior is preserved.

That gives us a stable seam for BINSEQ without coupling the project to a Rust
format dependency before we can measure and validate it.
