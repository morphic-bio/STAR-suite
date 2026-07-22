# Bulk

This is the easiest full demo.

The demo uses a small public bulk RNA-seq subset with `25000` read pairs.

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

These are output files, not shell commands.

## Verify outputs

Run the verifier:

```bash
bash scripts/codespaces/verify_bulk_public_demo_outputs.sh \
  --outdir .codespaces-demo/runs/bulk_public_demo
```

If you want to inspect the files manually:

```bash
ls -lh \
  .codespaces-demo/runs/bulk_public_demo/Aligned.sortedByCoord.out.bam \
  .codespaces-demo/runs/bulk_public_demo/ReadsPerGene.out.tab

head -n 8 .codespaces-demo/runs/bulk_public_demo/ReadsPerGene.out.tab
```

If you only prepared a dry run, execute the saved command script with:

```bash
bash .codespaces-demo/runs/bulk_public_demo/RUN_COMMAND.sh
```

## If you already know STAR

This guide mostly saves time by downloading a small public dataset and writing the full command for you.

## Using your own data

For a real bulk dataset, replace the demo FASTQs with your own files and point `--genome-dir` to the STAR index built from the reference you want to use.

See [Using your own data](./08_using_your_own_data.md).
