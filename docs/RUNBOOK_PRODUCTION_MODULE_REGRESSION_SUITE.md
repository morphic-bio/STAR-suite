# Production Module Regression Suite

Status: initial manifest and preflight runner are implemented on
`feature/binseq-input`.

## Goal

Maintain one canonical production-file regression case per STAR-suite module,
plus smaller contract cases for interfaces that are not ready for production
scale. The suite is manifest-driven so agents can discover what should run,
preflight local fixtures, and execute only the selected cases.

## Files

- Manifest: `tests/production_module_regression_manifest.tsv`
- Runner: `tests/run_production_module_regression_suite.sh`
- Wrapper logs: `/tmp/star_suite_production_module_regression_<timestamp>/`

## Basic Commands

List the registered cases:

```bash
tests/run_production_module_regression_suite.sh --list
```

Preflight without running expensive jobs:

```bash
tests/run_production_module_regression_suite.sh --preflight
```

Run only contract cases:

```bash
tests/run_production_module_regression_suite.sh --run --scale contract
```

Run optional network-fetched upstream fixtures:

```bash
tests/run_production_module_regression_suite.sh --run --scale contract-network
```

Run one production case:

```bash
tests/run_production_module_regression_suite.sh --run --case ucsf-corrected-100k
```

Run all 100K production-file cases only after confirming the machine is free for
long jobs:

```bash
tests/run_production_module_regression_suite.sh --preflight --scale production-100k
tests/run_production_module_regression_suite.sh --run --scale production-100k --continue-on-fail
```

## Current Case Policy

Contract cases validate the input/module boundary and should be safe for routine
agent use:

- `input-fastx / fastx-contract`
- `input-binseq / binseq-conversion-parity`
- `core-yremove / yremove-fastq-se`
- `feature-tools / assignbarcodes-regression`

Network-fetched contract cases validate external fixtures but are not part of
the offline contract set:

- `input-binseq / binseq-upstream-fixture`

Production-file cases use stable local fixtures and are not part of routine
developer smoke unless requested:

- `cr-compat / ucsf-corrected-100k`
- `feature-tools / pf-dynamic-permit-100k`
- `flex / flex-hash-screen-100k`
- `slam / slam-pe-100k`
- `multiome / chromap-macs3-100k`

The manifest is the source of truth for required fixture paths, tools, and
per-case environment fallbacks. If a path moves, update the manifest rather than
embedding host assumptions in the runner.

## BINSEQ Note

`binseq-conversion-parity` intentionally converts FASTQ to CBQ with `bqtools`
and compares the BINSEQ probe module against the FASTX module. This is a
contract smoke, not a real external BINSEQ fixture. The comparison is
order-independent because the input contract does not require BINSEQ modules to
preserve original FASTQ/source order. The case covers paired and single-end CBQ
with default zstd compression and uncompressed `-l 0` CBQ. Phase 5 for BINSEQ is
not complete until at least one externally produced paired CBQ and one
externally produced single-end CBQ pass the same harness.

`binseq-upstream-fixture` clones `ArcInstitute/binseq` and validates the
repository's `data/subset.cbq` fixture. That file is an externally produced
paired CBQ. The bundled `subset_R1.fastq.gz` and `subset_R2.fastq.gz` contain
the same reads in a different order, so this case checks order-independent
FASTQ parity rather than direct ordered `diff`.

## Adding A Case

1. Add one manifest row with:
   `module`, `case`, `scale`, `script`, `requirements`, `env`, `description`.
2. Keep requirements explicit. Use `file:`, `dir:`, `cmd:`, or
   `any:alt1|alt2` so preflight can explain why a case will skip.
3. Put generated outputs under `/tmp`, `tests/*_output*/`, or another path
   already listed in `tests/ARTIFACTS.md`.
4. Prefer scripts that write their own command manifests and validation
   summaries.
5. Do not add full-depth production launches to this suite unless they are
   gated by a distinct `scale` value and are never selected by default.

## Agent Policy

- Default to `--preflight` when asked for status.
- Run `--scale contract` for quick validation after module-boundary edits.
- Run `--scale contract-network` only when network access is acceptable.
- Ask or wait for explicit instruction before running `production-100k` cases.
- Never run production-file regression cases in parallel with benchmarks or
  other long STAR jobs from the same checkout.
