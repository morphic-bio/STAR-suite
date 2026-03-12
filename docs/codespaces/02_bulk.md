# Bulk

This is the easiest full demo.

## What you need

- a STAR index
- if you do not have one, you can use your own or start with the small setup guide

## Preview first

```bash
bash scripts/codespaces/run_bulk_public_demo.sh \
  --dry-run \
  --genome-dir /path/to/star_bulk_index
```

## Run it

```bash
bash scripts/codespaces/run_bulk_public_demo.sh \
  --run \
  --genome-dir /path/to/star_bulk_index
```

## What you get

- `.codespaces-demo/runs/bulk_public_demo/RUN_COMMAND.sh`
- `.codespaces-demo/runs/bulk_public_demo/Aligned.sortedByCoord.out.bam`
- `.codespaces-demo/runs/bulk_public_demo/ReadsPerGene.out.tab`

## If you already know STAR

This guide mostly saves time by downloading a small public dataset and writing the full command for you.
