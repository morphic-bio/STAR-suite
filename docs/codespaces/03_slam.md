# SLAM

This is the small public SLAM demo.

## What you need

- a STAR index
- if you do not have one, you can use your own or start with the small setup guide

## Preview first

```bash
bash scripts/codespaces/run_slam_public_demo.sh \
  --dry-run \
  --genome-dir /path/to/star_bulk_index
```

## Run it

```bash
bash scripts/codespaces/run_slam_public_demo.sh \
  --run \
  --genome-dir /path/to/star_bulk_index
```

## What you get

- `.codespaces-demo/runs/slam_public_demo/RUN_COMMAND.sh`
- `.codespaces-demo/runs/slam_public_demo/Aligned.sortedByCoord.out.bam`
- `.codespaces-demo/runs/slam_public_demo/slam_qc/`

## If you already know STAR

The main new pieces here are the SLAM-specific options. The guide still writes the full command so you can inspect or edit it.

## Using your own data

For a real SLAM dataset, replace the demo FASTQ with your own file and point `--genome-dir` to the STAR index built from the reference you want to use.

See [Using your own data](./08_using_your_own_data.md).
