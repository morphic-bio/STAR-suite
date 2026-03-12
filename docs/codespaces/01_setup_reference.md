# Setup + Mini Reference

Goal:
- prepare the shared public `chr22 + chrY` source reference
- build the STAR demo index inside Codespaces

This is the common setup surface for the single-cell walkthroughs.

Commands:
```bash
bash scripts/codespaces/fetch_public_chr22y_reference.sh
bash scripts/codespaces/build_public_chr22y_index.sh --threads 4
```

Outputs:
- `.codespaces-demo/data/public_human_chr22y_ref/fasta/genome.fa`
- `.codespaces-demo/data/public_human_chr22y_ref/genes/genes.gtf`
- `.codespaces-demo/indices/public_human_chr22y_star`

Why this is the reference contract:
- the important compatibility surface is the processed `FASTA + GTF`
- we do not introduce `cellranger mkref` into the Codespaces path
- users can replace this source surface with their own full reference later
