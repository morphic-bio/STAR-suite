# Perturb (work in progress)

This guide is not finished yet.

You can still use it to see the current command shape and file layout.

## Preview the current version

```bash
bash scripts/codespaces/run_perturb_public_demo.sh --dry-run
```

## Run the current version

```bash
bash scripts/codespaces/run_perturb_public_demo.sh
```

## What works now

- small demo reference setup
- small shared single-cell input
- CR-style config input path

## What is still missing

- a stronger public perturb demo dataset
- a cleaner end-to-end demo result

## Using your own data

For a real perturb run, you will usually need your own GEX FASTQs, guide FASTQs, whitelist, barcode settings, and feature reference.

Also set the strand to match your library. The demo command uses `--soloStrand Unstranded` for its small fixture. Use `Forward` for 3' libraries and `Reverse` for 10x 5' libraries sequenced from read 2 only (for example A375). The options used for the manuscript's Perturb-seq benchmarks are listed in [PAPER_BENCHMARK_METHODOLOGY.md](../PAPER_BENCHMARK_METHODOLOGY.md) Section 1.6.

See [Using your own data](./08_using_your_own_data.md).
