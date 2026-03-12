# Single-cell Fixture Derivation

Current status: stable.

This is the shared single-cell preparation path. It derives a bounded public GEX fixture by aligning a limited number of read pairs from a public male 10x dataset to the shared `chr22 + chrY` mini-reference and retaining only read pairs whose cDNA mates map to `chr22` or `chrY`.

Preview:
```bash
bash scripts/codespaces/derive_public_chr22y_gex_fixture.sh \
  --threads 4 \
  --read-limit 200000 \
  --dry-run
```

Run:
```bash
bash scripts/codespaces/derive_public_chr22y_gex_fixture.sh \
  --threads 4 \
  --read-limit 200000
```

Outputs:
- `.codespaces-demo/data/public_10x_male_bcell/raw_fastqs/`
- `.codespaces-demo/data/public_chr22y_gex_fixture/gex/public_chr22y_demo/`
- `.codespaces-demo/data/public_chr22y_gex_fixture/MANIFEST.txt`

Important note:
- barcode geometry is pinned as source metadata for the chosen public assay
- it is not treated as a global STAR-suite constant
