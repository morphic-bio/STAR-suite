# CBQ/BINSEQ Input

Status: released in STAR Suite `v1.1.0` for the STAR mapper, STARsolo, OCM,
Flex, SLAM, and process_features adapter surfaces listed below. Current BINSEQ
support validates `.cbq` files against the STAR-suite input contract through
both the original bqtools-backed probe and the native C++ CBQ reader, and it is
wired into the STAR mapper as a native input path:

```bash
STAR --readFilesType Binseq PE --readFilesIn sample.cbq ...
STAR --readFilesType Binseq SE --readFilesIn sample.cbq ...
```

Technical details for the on-disk CBQ subset, ordered encoder, native reader,
and STAR/process_features/Chromap adapter surfaces are documented in
`docs/CBQ_FORMAT_AND_IMPLEMENTATION.md`.

## What Works Now

- Paired `.cbq` files can be decoded through `bqtools` and read by
  `core/legacy/source/binseq_probe_harness`.
- Paired and single-end `.cbq` files can be read directly by the native C++
  prototype harness `core/legacy/source/cbq_reader_harness`.
- `core/legacy/source/cbq_ordered_encoder` creates paired or single-end CBQ
  from FASTQ/FASTQ.gz while preserving source record order. Use this encoder
  for production parity tests where read order can affect downstream
  deterministic boundaries.
- The native reader now emits a borrowed, block-backed `CbqReadBatchView`
  interchange format. The harness consumes that view directly; the older
  `InputRecord` path remains only as a compatibility adapter.
- `core/legacy/source/cbq_star_adapter_harness` verifies that the direct CBQ
  adapter produces byte-identical STAR internal read-buffer state relative to a
  FASTQ/readLoad reference emulation, and has a decoded-batch benchmark mode for
  isolating adapter cost after CBQ block decode.
- `core/legacy/source/cbq_pf_adapter_harness` verifies that CBQ
  barcode/feature records can stream through the process_features decoded input
  gate and existing consumer path without materializing temporary FASTQ.
- `core/legacy/source/cbq_chromap_adapter_harness` materializes synchronized
  R1/R2/barcode FASTQs from paired-read and barcode CBQ files for Chromap's
  current FASTQ-path contract.
- `STAR --readFilesType Binseq PE|SE` routes CBQ batches directly into STAR
  read buffers without materializing synthetic FASTQ streams.
- `tests/run_cbq_star_input_smoke.sh` maps the same synthetic paired and
  single-end reads from FASTQ and CBQ through production STAR and verifies
  byte-identical SAM body output, including manifest-style paired CBQ input,
  comma-separated multi-lane/multi-sample CBQ input, and paired-CBQ
  `--batchMode 1` output routing. Default-compressed and level-0 CBQ inputs
  are both covered.
- `tests/run_cbq_solo_e2e_smoke.sh` runs a synthetic STARsolo fixture through
  FASTQ and CBQ input and requires byte-identical raw Gene MEX output for
  direct CBQ, level-0 CBQ, and manifest-style CBQ input.
- Single-end `.cbq` files are covered by the synthetic conversion smoke.
- process_features CBQ input is covered by a FASTQ-vs-CBQ MEX/count parity
  smoke on a synthetic feature-barcode fixture with valid nucleotide UMIs. This
  currently exercises the harness/API surface; the production PF CLI flag is
  still pending.
- FLEX CBQ input uses the standard STAR CBQ adapter path for general
  genome-backed runs. Fully fused count-only FlexPipeline runs
  (`--flexPipelineNTriage 0 --flexPipelineNSolo 0 --flexNoAlign 1`) use the
  CBQ-native lane producer. `tests/run_cbq_flex_tiny_public_smoke.sh` covers
  FASTQ-vs-CBQ parity on a generated public tiny FLEX fixture; the host-local
  SC2300771 100K FLEX downsample also passed count parity and order-normalized
  BAM payload parity.
- FLEX count-only no-genome production is the first full-size topline CBQ use
  case. On SC2300771, level-0 CBQ no-genome completed in `8:38.52` versus
  FASTQ.gz no-genome in `12:01.95`, with byte-identical `Solo.out/Gene`,
  `Barcodes.stats`, and `per_sample_filtered` outputs.
- Chromap CBQ input is covered by a CBQ-to-Chromap-FASTQ adapter smoke and, on
  hosts with Chromap available, a tiny synthetic mapping run.
- `tests/run_cbq_e2e_module_regression.sh` runs the downsampled CBQ module
  suite for the reader, ordered encoder, STAR mapper, STARsolo,
  process_features, and Chromap adapter surfaces.
- Flex, OCM, and SLAM CBQ workflows are covered by production-shaped smokes in
  `tests/production_module_regression_manifest.tsv`.
- Direct `--readFilesIn sample.cbq` and manifest-style
  `CBQ<TAB>-<TAB>ReadGroup` probe inputs are covered.
- Default-compressed and uncompressed CBQ files are smoke-tested.
- The upstream ARC `data/subset.cbq` fixture is tested as an external paired
  CBQ example.

## Current Limitations

- FLEX has full-size production count-only no-genome parity/timing, but
  genome-backed FLEX alignment and BAM/SAM output surfaces still need separate
  full-size validation.
- Chromap integration currently adapts CBQ to Chromap's existing FASTQ path
  contract; it is not yet an in-memory libchromap reader API.
- SLAM per-file skipping is currently rejected for BINSEQ input.
- Y/noY FASTQ emission remains FASTQ-only until non-FASTQ input modules have an
  explicit emission contract.
- `--readFilesCommand` is not supported with BINSEQ input; pass `.cbq` files
  directly through `--readFilesIn` or `--readFilesManifest`.
- BINSEQ input order is not guaranteed by the generic input contract. Downstream
  code and tests must compare by read identity, sequence, quality, mate, and
  lane/read group unless a module advertises source-order preservation.
- CBQs produced by external parallel encoders may not preserve FASTQ source
  order. Use `cbq_ordered_encoder` when exact FASTQ-vs-CBQ order parity matters.

## How This Differs From FASTQ

FASTQ/FASTX remains the normal production input path. `STAR --readFilesType
Fastx` supports plain FASTQ, `.gz` FASTQ, manifest inputs, helper commands,
Solo, Flex, SLAM, and production alignment outputs.

BINSEQ now has a narrow native mapper path. The first probe module still decodes
CBQ through `bqtools` for oracle checks and normalized record dumps. The native
C++ reader reads CBQ directly, emits a CBQ-native batch view, and fills
STAR-style internal read buffers through the adapter used by the production
mapper path. This avoids decode-to-FASTQ in the STAR run while preserving the
optimized FASTQ path unchanged.

For FLEX count-only production runs, CBQ is now a validated performance path.
With the strict no-genome FLEX surface, both FASTQ.gz and level-0 CBQ avoid STAR
genome loading and produce byte-identical counts. On the full SC2300771
production run, CBQ completed in `8:38.52` while FASTQ.gz completed in
`12:01.95`, with essentially identical peak RSS (`43.4 GB` vs `43.3 GB`). The
remaining speedup is attributable to avoiding gzip FASTQ decode and using the
native CBQ lane reader for the same Flex hash/count workload.

## Production Implementation Direction

The production path should be native C++ unless Rust interop becomes trivial.
The current ARC `binseq` crate does not expose a stable C ABI, generated C
headers, or a simple `staticlib` / `cdylib` target for direct STAR linking.
`bqtools` remains useful for fixture generation and oracle comparison, but the
planned STAR integration should read CBQ directly in C++. The reader is
intentionally narrow: CBQ file/block/index parsing, zstd column decompression,
two-bit sequence decoding, block-level `N` position restoration, qualities,
headers, paired/single-end record iteration, and adapter handoff into native
app input structs.

For CBQ, "production support" means the consumer receives either an in-memory
read-batch view such as `CbqReadBatchView` or a direct adapter into its native
structs. The production STAR view keeps sequence in CBQ's packed two-bit block
storage and carries `N` positions as side metadata; it does not expand the full
block to ASCII sequence text. ASCII materialization is acceptable for oracle
tests, debugging, process_features compatibility, or explicit path-based
adapters, but it is not considered the production STAR CBQ integration surface.

## Required Tool

Install or provide `bqtools`, then set `BQTOOLS` if it is not on `PATH`:

```bash
export BQTOOLS=/path/to/bqtools
```

## Build The Probe Harness

```bash
make -C core/legacy/source binseq-probe-harness
```

To build STAR plus FASTX and BINSEQ comparison harnesses:

```bash
make -C core/legacy/source STAR fastx-input-harness binseq-probe-harness cbq-reader-harness cbq-ordered-encoder cbq-star-adapter-harness cbq-pf-adapter-harness cbq-chromap-adapter-harness
```

## Run The Smokes

Synthetic FASTQ-to-CBQ parity, covering paired and single-end CBQ with default
compression and uncompressed `-l 0`:

```bash
BQTOOLS=/path/to/bqtools tests/run_binseq_probe_smoke.sh
```

Native C++ CBQ reader parity, covering the same paired/single-end and
compression cases plus byte-identical STAR internal adapter dumps:

```bash
BQTOOLS=/path/to/bqtools tests/run_cbq_cpp_reader_smoke.sh
```

Ordered C++ FASTQ/FASTQ.gz-to-CBQ encoder smoke, covering paired, gzipped
paired, single-end, N-position round trip, and exact ordered contract parity:

```bash
tests/run_cbq_ordered_encoder_smoke.sh
```

Production STAR mapper smoke, covering paired and single-end FASTQ-vs-CBQ SAM
body parity plus manifest, multi-sample, paired-CBQ batch inputs, and level-0
CBQ:

```bash
BQTOOLS=/path/to/bqtools tests/run_cbq_star_input_smoke.sh
```

STARsolo E2E smoke, covering raw Gene MEX parity for FASTQ, direct CBQ,
level-0 CBQ, and manifest CBQ:

```bash
BQTOOLS=/path/to/bqtools tests/run_cbq_solo_e2e_smoke.sh
```

process_features adapter smoke, covering gzipped FASTQ reference input versus
native CBQ in-memory records:

```bash
BQTOOLS=/path/to/bqtools tests/run_cbq_pf_adapter_smoke.sh
```

Chromap adapter smoke, covering paired-read CBQ plus barcode CBQ materialized
to Chromap-compatible FASTQs. If `CHROMAP_BIN` exists, the smoke also runs a
tiny Chromap mapping check:

```bash
BQTOOLS=/path/to/bqtools tests/run_cbq_chromap_adapter_smoke.sh
```

Downsampled CBQ E2E/module regression suite:

```bash
BQTOOLS=/path/to/bqtools tests/run_cbq_e2e_module_regression.sh
```

Include network/public-fixture CBQ smokes, including FLEX tiny public:

```bash
RUN_NETWORK=1 BQTOOLS=/path/to/bqtools tests/run_cbq_e2e_module_regression.sh
```

FLEX FASTQ-vs-CBQ public tiny smoke only:

```bash
tests/run_cbq_flex_tiny_public_smoke.sh
```

Upstream ARC paired-CBQ fixture smoke:

```bash
BQTOOLS=/path/to/bqtools tests/run_binseq_upstream_fixture_smoke.sh
```

The upstream fixture smoke uses network access to clone `ArcInstitute/binseq`.

## Probe A CBQ Manually

Direct paired-CBQ probe:

```bash
core/legacy/source/binseq_probe_harness \
  --readFilesIn sample.cbq \
  --mateCount 2 \
  --bqtools "$BQTOOLS" \
  --workDir /tmp/star_suite_binseq_probe_sample \
  > sample.binseq_probe.tsv
```

Manifest-style paired-CBQ probe:

```bash
printf '%s\t-\tID:sample1\n' sample.cbq > binseq_manifest.tsv

core/legacy/source/binseq_probe_harness \
  --readFilesManifest binseq_manifest.tsv \
  --mateCount 2 \
  --bqtools "$BQTOOLS" \
  --workDir /tmp/star_suite_binseq_probe_manifest \
  > sample.manifest_probe.tsv
```

For single-end CBQ, use `--mateCount 1`.

Experimental STAR alignment path:

```bash
STAR \
  --runMode alignReads \
  --genomeDir /path/to/star-index \
  --readFilesType Binseq PE \
  --readFilesIn sample.cbq \
  --outFileNamePrefix out/
```

For manifest input, use one CBQ path in column 1 and `-` in column 2:

```bash
printf '%s\t-\tID:sample1\n' sample.cbq > binseq_manifest.tsv

STAR \
  --runMode alignReads \
  --genomeDir /path/to/star-index \
  --readFilesType Binseq PE \
  --readFilesManifest binseq_manifest.tsv \
  --outFileNamePrefix out/
```

Native C++ reader probe:

```bash
core/legacy/source/cbq_reader_harness \
  --readFilesIn sample.cbq \
  --mateCount 2 \
  > sample.cbq_reader.tsv
```

STAR internal adapter probe:

```bash
core/legacy/source/cbq_star_adapter_harness \
  --readFilesIn sample.cbq \
  --mateCount 2 \
  --mode direct \
  > sample.cbq_star_adapter.direct.bin

core/legacy/source/cbq_star_adapter_harness \
  --readFilesIn sample.cbq \
  --mateCount 2 \
  --mode reference \
  > sample.cbq_star_adapter.reference.bin

cmp sample.cbq_star_adapter.direct.bin sample.cbq_star_adapter.reference.bin
```

Decoded-batch STAR adapter benchmark:

```bash
core/legacy/source/cbq_star_adapter_harness \
  --readFilesIn sample.cbq \
  --mateCount 2 \
  --mode benchmark \
  --maxRecords 1000000 \
  --iterations 5
```

The benchmark preloads decoded `CbqReadBatchView` batches once, then replays the
STAR adapter into reusable STAR read buffers. `preload_seconds` measures CBQ
block decode plus view construction; `adapter_seconds` measures the STAR adapter
only.

process_features adapter probe:

```bash
core/legacy/source/cbq_pf_adapter_harness \
  --mode cbq \
  --readFilesIn feature_reads.cbq \
  --whitelist whitelist.txt \
  --featureRef features.csv \
  --outputDir pf_out \
  --sampleName sample
```

Chromap adapter probe:

```bash
core/legacy/source/cbq_chromap_adapter_harness \
  --readPairCbq atac_reads.cbq \
  --barcodeCbq atac_barcodes.cbq \
  --outputDir chromap_fastqs \
  --sampleName sample
```

## Output

The probe harness writes one normalized TSV row per mate record. This is a
contract validation format, not a STAR alignment output. Production STAR runs
with `--readFilesType Binseq` write the same alignment output surfaces as the
normal mapper path for the modes that have been validated.

## Developer References

- Technical format and implementation reference:
  `docs/CBQ_FORMAT_AND_IMPLEMENTATION.md`
- Implementation runbook: `docs/RUNBOOK_BINSEQ_INPUT_CONTRACT.md`
- Native reader and STAR adapter prototype:
  `docs/RUNBOOK_BINSEQ_CPP_READER_PROTOTYPE.md`
- process_features native CBQ production plan:
  `docs/RUNBOOK_PROCESS_FEATURES_CBQ_NATIVE.md`
- Chromap ATAC in-memory integration plan:
  `docs/RUNBOOK_CHROMAP_ATAC_CBQ_IN_MEMORY.md`
- Synthetic smoke: `tests/run_binseq_probe_smoke.sh`
- Native C++ reader smoke: `tests/run_cbq_cpp_reader_smoke.sh`
- STAR mapper smoke: `tests/run_cbq_star_input_smoke.sh`
- STARsolo E2E smoke: `tests/run_cbq_solo_e2e_smoke.sh`
- process_features adapter smoke: `tests/run_cbq_pf_adapter_smoke.sh`
- Chromap adapter smoke: `tests/run_cbq_chromap_adapter_smoke.sh`
- CBQ module regression suite: `tests/run_cbq_e2e_module_regression.sh`
- Upstream fixture smoke: `tests/run_binseq_upstream_fixture_smoke.sh`
- Artifact policy: `tests/ARTIFACTS.md`
