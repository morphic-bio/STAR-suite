# Optional setup: small demo reference

Use this guide only if you want to build the small demo reference inside Codespaces.

If you already have a STAR index that you want to use, skip this guide.

## Run

```bash
bash scripts/codespaces/fetch_public_chr22y_reference.sh
bash scripts/codespaces/build_public_chr22y_index.sh --threads 4
```

## What you get

- `.codespaces-demo/data/public_human_chr22y_ref/fasta/genome.fa`
- `.codespaces-demo/data/public_human_chr22y_ref/genes/genes.gtf`
- `.codespaces-demo/indices/public_human_chr22y_star`

## In plain terms

This builds a very small human reference that only contains `chr22` and `chrY`.
It is small enough to use in Codespaces and is shared by the single-cell demo guides.
