# Experimental BINSEQ Input

Status: experimental on `dev-release`. Current BINSEQ support validates `.cbq`
files against the STAR-suite input contract through both the original
bqtools-backed probe and the native C++ CBQ reader prototype. It is not yet a
production STAR alignment input path.

## What Works Now

- Paired `.cbq` files can be decoded through `bqtools` and read by
  `core/legacy/source/binseq_probe_harness`.
- Paired and single-end `.cbq` files can be read directly by the native C++
  prototype harness `core/legacy/source/cbq_reader_harness`.
- Single-end `.cbq` files are covered by the synthetic conversion smoke.
- Direct `--readFilesIn sample.cbq` and manifest-style
  `CBQ<TAB>-<TAB>ReadGroup` probe inputs are covered.
- Default-compressed and uncompressed CBQ files are smoke-tested.
- The upstream ARC `data/subset.cbq` fixture is tested as an external paired
  CBQ example.

## What Is Not Production-Ready Yet

- Normal STAR mapping with `STAR --readFilesType Binseq PE|SE` is not wired
  through the production CLI yet.
- Solo, Flex, SLAM, batch mode, and full alignment outputs have not been
  validated directly from BINSEQ input.
- Y/noY FASTQ emission remains FASTQ-only until non-FASTQ input modules have an
  explicit emission contract.
- BINSEQ input order is not guaranteed by the generic input contract. Downstream
  code and tests must compare by read identity, sequence, quality, mate, and
  lane/read group unless a module advertises source-order preservation.

## How This Differs From FASTQ

FASTQ/FASTX is the normal production input path. `STAR --readFilesType Fastx`
is wired into the STAR execution path and supports plain FASTQ, `.gz` FASTQ,
manifest inputs, helper commands, Solo, Flex, SLAM, and production alignment
outputs.

BINSEQ is currently a probe/contract path. The first probe module decodes CBQ
through `bqtools`, feeds the decoded records through the same input-contract
interface, and writes normalized record dumps for parity checks. The native C++
prototype reads CBQ directly and emits the same contract stream, proving the
reader can reconstruct the FASTQ-equivalent input state before we route CBQ
through the full mapper.

## Production Implementation Direction

The production path should be native C++ unless Rust interop becomes trivial.
The current ARC `binseq` crate does not expose a stable C ABI, generated C
headers, or a simple `staticlib` / `cdylib` target for direct STAR linking.
`bqtools` remains useful for fixture generation and oracle comparison, but the
planned STAR integration should read CBQ directly in C++. The planned reader is
intentionally narrow: CBQ file/block/index parsing, zstd column decompression,
two-bit sequence decoding, block-level `N` position restoration, qualities,
headers, and paired/single-end record iteration.

## Required Tool

Install or provide `bqtools`, then set `BQTOOLS` if it is not on `PATH`:

```bash
export BQTOOLS=/path/to/bqtools
```

## Build The Probe Harness

```bash
make -C core/legacy/source binseq-probe-harness
```

To build FASTX plus both BINSEQ comparison harnesses:

```bash
make -C core/legacy/source fastx-input-harness binseq-probe-harness cbq-reader-harness
```

## Run The Smokes

Synthetic FASTQ-to-CBQ parity, covering paired and single-end CBQ with default
compression and uncompressed `-l 0`:

```bash
BQTOOLS=/path/to/bqtools tests/run_binseq_probe_smoke.sh
```

Native C++ CBQ reader parity, covering the same paired/single-end and
compression cases:

```bash
BQTOOLS=/path/to/bqtools tests/run_cbq_cpp_reader_smoke.sh
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

Native C++ reader probe:

```bash
core/legacy/source/cbq_reader_harness \
  --readFilesIn sample.cbq \
  --mateCount 2 \
  > sample.cbq_reader.tsv
```

## Output

The probe harness writes one normalized TSV row per mate record. This is a
contract validation format, not a STAR alignment output. It is intended for
parity checks against FASTQ and for validating BINSEQ metadata, mate handling,
sequence, quality, read names, lanes, and read groups.

## Developer References

- Implementation runbook: `docs/RUNBOOK_BINSEQ_INPUT_CONTRACT.md`
- Native reader prototype: `docs/RUNBOOK_BINSEQ_CPP_READER_PROTOTYPE.md`
- Synthetic smoke: `tests/run_binseq_probe_smoke.sh`
- Native C++ reader smoke: `tests/run_cbq_cpp_reader_smoke.sh`
- Upstream fixture smoke: `tests/run_binseq_upstream_fixture_smoke.sh`
- Artifact policy: `tests/ARTIFACTS.md`
