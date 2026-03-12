# Single-cell input for perturb and Flex

Use this guide only if you want the small shared single-cell input used by the perturb and Flex demos.

You do not need this guide for bulk or SLAM.

## Preview first

```bash
bash scripts/codespaces/derive_public_chr22y_gex_fixture.sh \
  --threads 4 \
  --read-limit 200000 \
  --dry-run
```

## Run it

```bash
bash scripts/codespaces/derive_public_chr22y_gex_fixture.sh \
  --threads 4 \
  --read-limit 200000
```

## What you get

- `.codespaces-demo/data/public_10x_male_bcell/raw_fastqs/`
- `.codespaces-demo/data/public_chr22y_gex_fixture/gex/public_chr22y_demo/`
- `.codespaces-demo/data/public_chr22y_gex_fixture/MANIFEST.txt`

## In plain terms

This takes a small amount of public single-cell data and makes an even smaller input set that fits more easily in Codespaces.
