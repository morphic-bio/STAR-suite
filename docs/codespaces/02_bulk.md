# Bulk Demo

Current status: stable.

This walkthrough uses a small public bulk RNA-seq fixture and emits the STAR bulk command shape. The script supports `--run`, but the recommended first step in Codespaces is to preview it with `--dry-run` and then swap in the index you actually want to use.

Preview:
```bash
bash scripts/codespaces/run_bulk_public_demo.sh \
  --dry-run \
  --genome-dir /path/to/star_bulk_index
```

Run:
```bash
bash scripts/codespaces/run_bulk_public_demo.sh \
  --run \
  --genome-dir /path/to/star_bulk_index
```

Outputs:
- `.codespaces-demo/runs/bulk_public_demo/RUN_COMMAND.sh`
- `.codespaces-demo/runs/bulk_public_demo/Aligned.sortedByCoord.out.bam`
- `.codespaces-demo/runs/bulk_public_demo/ReadsPerGene.out.tab`
