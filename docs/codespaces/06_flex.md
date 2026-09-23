# Flex (work in progress)

This guide is not finished yet.

You can still use it to see the current command shape and file layout.

## Preview the current version

```bash
bash scripts/codespaces/run_flex_public_demo.sh --dry-run
```

## Run the current version

```bash
bash scripts/codespaces/run_flex_public_demo.sh
```

## What works now

- small demo reference setup
- small shared single-cell input
- Flex config and wrapper wiring

## What is still missing

- a stronger public Flex demo dataset
- a cleaner end-to-end demo result

## Using your own data

For a real Flex run, you will usually need your own FASTQs, whitelist, barcode settings, probe set, and a new filtered reference and STAR index. From STAR Suite 1.9.4, Flex runs also need a half-probe cache built once from that index; the demo wrapper builds it for you.

See [Using your own data](./08_using_your_own_data.md).
