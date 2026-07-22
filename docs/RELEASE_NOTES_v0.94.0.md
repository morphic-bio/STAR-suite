# STAR Suite v0.94.0 Release Notes

Release prep date: 2026-05-12

This release note is prepared from `master` after `v0.93.1`. Do not tag or
build release artifacts from this host while the SLAM PE production run is
active; use the GitHub release workflow or wait until the local production run
is complete.

## Highlights

- SLAM paired-end production support:
  - PE transition coordinate handling and overlap consensus before counting.
  - PE noSU-derived trim surface: R1 `8/13`, R2 `19/14`; SE smoke trim `8/12`.
  - `--slamMinCallableLength` gate for SLAM transition/NTR/cB/GrandSLAM evidence.
  - Count-binomial output via `--slamCbOut` with `star` and `ezbakr` formats.
  - TranscriptVB tximport gene-count surface for DESeq2.
  - Production runner with Y/noY BAM/FASTQ emission, Globus transfer, and
    post-transfer cleanup.
- Reproducible SLAM PE smoke:
  - `scripts/run_slam_100k_se_pe_smoke.sh`
  - Treatment 100K PE vs R1-only SE NTR Pearson `0.972806`.
  - Tximport gene NumReads Pearson `0.932203` treatment and `0.938455` noSU.
- Pinned SLAM DESeq2 image:
  - `star-suite/slam-deseq2:bioc3.23-deseq2-1.52.0-tximport-1.40.0`
  - Bioconductor `3.23`, R `4.6.0`, DESeq2 `1.52.0`, tximport `1.40.0`.
- MCP/Launchpad updates:
  - Private workflows for `slam_pe_100k_smoke`, `slam_pe_production`, and
    `slam_deseq2_container`.
  - Launchpad ordering keeps STAR CLI recipes first, then SLAM private recipes
    when "Include test & other recipes" is enabled.
  - Workflow script metadata and trusted roots updated for the local SLAM PE
    dataset and output roots.
- ATAC and Launchpad work since `v0.93.1` remains included:
  - libchromap/Chromap ATAC integration and evidence sidecars.
  - ATAC permit-domain and drain-time controller work.
  - Launchpad recipe builder, server browse/upload, Script Lane, composition
    tab, and editor bridge.

## Release Gates

Run only after active local production/benchmark jobs are clear:

```bash
python3 -m pytest \
  mcp_server/tests/test_workflow_config.py \
  mcp_server/tests/test_workflow_discovery.py \
  mcp_server/tests/test_workflow_render.py \
  mcp_server/tests/test_launchpad.py
```

Before tagging, run the normal clean build and release smoke matrix on an idle
host or in CI. The local documentation/config update intentionally stops before
any compilation.

## Candidate Tag

```text
v0.94.0
```
