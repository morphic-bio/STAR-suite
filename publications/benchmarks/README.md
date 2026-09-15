# Publication Benchmark Scripts

This directory provides publication-facing entrypoints for the benchmark and
parity surfaces summarized in the top-level `README.md`.

The wrappers here do not reimplement the benchmark logic. They pin the
canonical script to run for each surface so readers have a stable place to
start. The protocol, comparator commands and parity definitions behind the
published STAR Suite 1.9.4 numbers are in
[`docs/PAPER_BENCHMARK_METHODOLOGY.md`](../../docs/PAPER_BENCHMARK_METHODOLOGY.md)
(Section 1), which also lists the corrections applied to these wrappers for
the paper runs.

## Public vs private data

- Public CI/release smoke coverage must avoid the protected UCSF and genomic
  JAX surfaces.
- JAX Flex is allowed for publication and public validation because the
  benchmarked surface is probe-based rather than exposing genomic Y content.
- Of the manuscript datasets, PPARG bulk (GEO GSE288287), 10x PBMC 10K, A375,
  GSE325982, the 10x 320k scFFPE set and the GRAND-SLAM fixture are public; the
  MSK 30-KO and JAX SC2300771 datasets are MorPhiC consortium data.

## Manuscript benchmark set (STAR Suite 1.9.4)

| Surface | Publication wrapper | Canonical implementation | Data class |
|---|---|---|---|
| Bulk RNA-seq, PPARG with and without Y removal | `run_bulk_rnaseq_benchmarks.sh` | `scripts/paper/run_pe_bulk_feature_benchmark.sh --integrated-only`; the external reference pipeline is given in the methodology | public |
| scRNA-seq, 10x PBMC 10K | — | direct STAR invocation (methodology, Section 1.6) | public |
| Perturb-seq, A375 | `run_a375_perturb_benchmark.sh` | `scripts/paper/run_a375_benchmark.sh` | public |
| Perturb-seq, MSK 30-KO ES | `run_msk_perturb_benchmark.sh` | `scripts/paper/run_msk_30polyko_benchmark.sh` (ES staging) | consortium |
| 10x Flex, JAX SC2300771 / GSE325982 / 320k scFFPE | — | direct STAR invocation (methodology, Section 1.6) | consortium / public / public |
| SLAM-seq, GRAND-SLAM 100K fixture | — | `tests/run_slam_fixture_parity.sh` | public |

## Other surfaces (not part of the manuscript)

| Surface | Publication wrapper | Canonical implementation | Data class |
|---|---|---|---|
| UCSF GEX-only Solo | `run_ucsf_gexonly_solo_benchmarks.sh` | `scripts/paper/run_ucsf_gexonly_no_bam_benchmark.sh` | private |
| Velocyto bridge (`Gene/GeneFull/Velocyto`) | `run_velocyto_bridge_benchmarks.sh` | `scripts/run_star_velocyto_gexonly_canonical.sh` | private |
| UCSF perturb benchmark | `run_ucsf_perturb_benchmark.sh` | `scripts/paper/run_ucsf_ebs2_2_benchmark.sh` | private |
| Flex Cell Ranger-config run surface | `run_flex_benchmark.sh` | `scripts/run_flex_cr_config.sh` | private/JAX Flex |
| Flex parity summary | `run_flex_parity.sh` | `scripts/paper/run_flex_parity.sh` | private/JAX Flex |
| SLAM PE smoke / production surface | `run_slam_benchmark.sh` | `scripts/run_slam_prod_set.sh` (`scripts/run_slam_100k_se_pe_smoke.sh` for 100K SE/PE smoke) | private |

## Public smoke equivalents

These are the Docker/public validation surfaces used for CI and release smoke
checks rather than the full paper fixtures:

- `tests/run_solo_smoke.sh`
- `tests/slam/test_snp_mask_build_smoke.sh`
- `scripts/run_slam_100k_se_pe_smoke.sh` (private local FASTQ fixture; use for pre-production PE validation, not public CI)
- `tests/run_flex_tiny_public_smoke.sh`
