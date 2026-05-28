# Experimental BINSEQ Input

Status: experimental on `dev-release`. Current BINSEQ support validates `.cbq`
files against the STAR-suite input contract through both the original
bqtools-backed probe and the native C++ CBQ reader, and it is wired into the
STAR mapper as an experimental native input path:

```bash
STAR --readFilesType Binseq PE --readFilesIn sample.cbq ...
STAR --readFilesType Binseq SE --readFilesIn sample.cbq ...
```

## What Works Now

- Paired `.cbq` files can be decoded through `bqtools` and read by
  `core/legacy/source/binseq_probe_harness`.
- Paired and single-end `.cbq` files can be read directly by the native C++
  prototype harness `core/legacy/source/cbq_reader_harness`.
- The native reader now emits a borrowed, block-backed `CbqReadBatchView`
  interchange format. The harness consumes that view directly; the older
  `InputRecord` path remains only as a compatibility adapter.
- `core/legacy/source/cbq_star_adapter_harness` verifies that the direct CBQ
  adapter produces byte-identical STAR internal read-buffer state relative to a
  FASTQ/readLoad reference emulation.
- `core/legacy/source/cbq_pf_adapter_harness` verifies a direct
  process_features in-memory record API for CBQ barcode/feature records without
  materializing temporary FASTQ.
- `core/legacy/source/cbq_chromap_adapter_harness` materializes synchronized
  R1/R2/barcode FASTQs from paired-read and barcode CBQ files for Chromap's
  current FASTQ-path contract.
- `STAR --readFilesType Binseq PE|SE` routes CBQ batches directly into STAR
  read buffers without materializing synthetic FASTQ streams.
- `tests/run_cbq_star_input_smoke.sh` maps the same synthetic paired and
  single-end reads from FASTQ and CBQ through production STAR and verifies
  byte-identical SAM body output, including manifest-style paired CBQ input,
  comma-separated multi-lane/multi-sample CBQ input, and paired-CBQ
  `--batchMode 1` output routing.
- Single-end `.cbq` files are covered by the synthetic conversion smoke.
- process_features CBQ input is covered by a FASTQ-vs-CBQ MEX/count parity
  smoke on a synthetic feature-barcode fixture with valid nucleotide UMIs.
- Chromap CBQ input is covered by a CBQ-to-Chromap-FASTQ adapter smoke and, on
  hosts with Chromap available, a tiny synthetic mapping run.
- Direct `--readFilesIn sample.cbq` and manifest-style
  `CBQ<TAB>-<TAB>ReadGroup` probe inputs are covered.
- Default-compressed and uncompressed CBQ files are smoke-tested.
- The upstream ARC `data/subset.cbq` fixture is tested as an external paired
  CBQ example.

## Experimental Limitations

- Solo, Flex, SLAM, batch mode, and large production alignment outputs have not
  been validated directly from BINSEQ input.
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

## Production Implementation Direction

The production path should be native C++ unless Rust interop becomes trivial.
The current ARC `binseq` crate does not expose a stable C ABI, generated C
headers, or a simple `staticlib` / `cdylib` target for direct STAR linking.
`bqtools` remains useful for fixture generation and oracle comparison, but the
planned STAR integration should read CBQ directly in C++. The reader is
intentionally narrow: CBQ file/block/index parsing, zstd column decompression,
two-bit sequence decoding, block-level `N` position restoration, qualities,
headers, paired/single-end record iteration, and adapter handoff into native
app input structs. Production app integrations should adapt from
`CbqReadBatchView` spans into the target app's native input structs rather than
materializing generic records or synthetic FASTQ.

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
make -C core/legacy/source STAR fastx-input-harness binseq-probe-harness cbq-reader-harness cbq-star-adapter-harness cbq-pf-adapter-harness cbq-chromap-adapter-harness
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

Production STAR mapper smoke, covering paired and single-end FASTQ-vs-CBQ SAM
body parity plus manifest, multi-sample, and paired-CBQ batch inputs:

```bash
BQTOOLS=/path/to/bqtools tests/run_cbq_star_input_smoke.sh
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

- Implementation runbook: `docs/RUNBOOK_BINSEQ_INPUT_CONTRACT.md`
- Native reader and STAR adapter prototype:
  `docs/RUNBOOK_BINSEQ_CPP_READER_PROTOTYPE.md`
- Synthetic smoke: `tests/run_binseq_probe_smoke.sh`
- Native C++ reader smoke: `tests/run_cbq_cpp_reader_smoke.sh`
- STAR mapper smoke: `tests/run_cbq_star_input_smoke.sh`
- process_features adapter smoke: `tests/run_cbq_pf_adapter_smoke.sh`
- Chromap adapter smoke: `tests/run_cbq_chromap_adapter_smoke.sh`
- Upstream fixture smoke: `tests/run_binseq_upstream_fixture_smoke.sh`
- Artifact policy: `tests/ARTIFACTS.md`
