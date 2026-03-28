# Publication Benchmark Scripts

This directory provides publication-facing entrypoints for the benchmark and
parity surfaces summarized in the top-level `README.md`.

The wrappers here do not reimplement the benchmark logic. They pin the
canonical script to run for each README surface so readers have a stable place
to start.

## Public vs private data

- Public CI/release smoke coverage must avoid the protected UCSF and genomic
  JAX surfaces.
- JAX Flex is allowed for publication and public validation because the
  benchmarked surface is probe-based rather than exposing genomic Y content.
- Several full paper benchmarks still depend on local `/storage` fixtures and
  are therefore documented here as private-data wrappers.

## Surface map

| README surface | Publication wrapper | Canonical implementation | Data class |
|---|---|---|---|
| Bulk RNA-seq wall time/parity | `run_bulk_rnaseq_benchmarks.sh` | `scripts/paper/run_pe_bulk_feature_benchmark.sh` | private |
| UCSF GEX-only Solo | `run_ucsf_gexonly_solo_benchmarks.sh` | `scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh` | private |
| Velocyto bridge (`Gene/GeneFull/Velocyto`) | `run_velocyto_bridge_benchmarks.sh` | `scripts/run_star_velocyto_gexonly_canonical.sh` | private |
| A375 perturb benchmark | `run_a375_perturb_benchmark.sh` | `scripts/paper/run_a375_benchmark.sh` | public/private local fixture wrapper |
| UCSF perturb benchmark | `run_ucsf_perturb_benchmark.sh` | `scripts/paper/run_ucsf_ebs2_2_benchmark.sh` | private |
| MSK perturb benchmark | `run_msk_perturb_benchmark.sh` | `scripts/paper/run_msk_30polyko_benchmark.sh` | private |
| Flex benchmark run surface | `run_flex_benchmark.sh` | `scripts/run_flex_cr_config.sh` | private/JAX Flex |
| Flex parity summary | `run_flex_parity.sh` | `scripts/paper/run_flex_parity.sh` | private/JAX Flex |
| SLAM production benchmark surface | `run_slam_benchmark.sh` | `scripts/run_slam_prod_batch_mode.sh` | private |

## Public smoke equivalents

These are the Docker/public validation surfaces used for CI and release smoke
checks rather than the full paper fixtures:

- `tests/run_solo_smoke.sh`
- `tests/slam/test_snp_mask_build_smoke.sh`
- `tests/run_flex_tiny_public_smoke.sh`
