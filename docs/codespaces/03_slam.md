# SLAM Demo

Current status: stable.

This walkthrough uses a small public SLAM-seq fixture and emits the STAR-SLAM command shape. As with the bulk walkthrough, the first step in Codespaces should be a `--dry-run` preview.

Preview:
```bash
bash scripts/codespaces/run_slam_public_demo.sh \
  --dry-run \
  --genome-dir /path/to/star_bulk_index
```

Run:
```bash
bash scripts/codespaces/run_slam_public_demo.sh \
  --run \
  --genome-dir /path/to/star_bulk_index
```

Outputs:
- `.codespaces-demo/runs/slam_public_demo/RUN_COMMAND.sh`
- `.codespaces-demo/runs/slam_public_demo/Aligned.sortedByCoord.out.bam`
- `.codespaces-demo/runs/slam_public_demo/slam_qc/`
