#!/usr/bin/env bash
# External fixture paths for STAR-suite tests.
# Source this file to set defaults for datasets that live outside the repo.
# Each variable honors any pre-existing value so you can override per shell.

# CR multi / A375 fixtures
export CR_MULTI_ROOT="${CR_MULTI_ROOT:-/storage/A375}"
export CR_MULTI_FASTQ_ROOT="${CR_MULTI_FASTQ_ROOT:-$CR_MULTI_ROOT/fastqs/1k_CRISPR_5p_gemx_fastqs}"
export CR_MULTI_GEX_DIR="${CR_MULTI_GEX_DIR:-$CR_MULTI_FASTQ_ROOT/gex}"
export CR_MULTI_CRISPR_DIR="${CR_MULTI_CRISPR_DIR:-$CR_MULTI_FASTQ_ROOT/crispr}"
export CR_MULTI_FEATURE_REF="${CR_MULTI_FEATURE_REF:-$CR_MULTI_ROOT/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv}"
export CR_MULTI_WHITELIST_GZ="${CR_MULTI_WHITELIST_GZ:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-5pgex-jan-2023.txt.gz}"
export CR_MULTI_WHITELIST="${CR_MULTI_WHITELIST:-$CR_MULTI_ROOT/3M-5pgex-jan-2023.txt}"
export CR_MULTI_DOWNSAMPLE_SCRIPT="${CR_MULTI_DOWNSAMPLE_SCRIPT:-/mnt/pikachu/process_features/scripts/downsample_fastq_directory.sh}"
export CR_MULTI_TIER_SCRIPT="${CR_MULTI_TIER_SCRIPT:-/mnt/pikachu/STAR-suite/scripts/a375_make_downsample_tiers.sh}"
export CR_GENOME_DIR="${CR_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
export CR_MULTI_GEX_OUTPREFIX="${CR_MULTI_GEX_OUTPREFIX:-/tmp/star_gex_smoke/}"
export CR_MULTI_OUTPREFIX="${CR_MULTI_OUTPREFIX:-/tmp/star_multi_smoke_cpp/}"

# Flex fixtures
export FLEX_INDEX="${FLEX_INDEX:-/storage/flex_filtered_reference/star_index}"
export FLEX_WHITELIST="${FLEX_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
export FLEX_PROBE_LIST="${FLEX_PROBE_LIST:-/storage/flex_filtered_reference/filtered_reference/probe_list.txt}"
export FLEX_SAMPLE_WHITELIST="${FLEX_SAMPLE_WHITELIST:-/storage/SC2300771_filtered_2M/sample_whitelist.tsv}"
export FLEX_SAMPLE_PROBES="${FLEX_SAMPLE_PROBES:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"
export FLEX_ALLOWED_TAGS="${FLEX_ALLOWED_TAGS:-/storage/SC2300771_filtered_2M/sample_whitelist.tsv}"
export FLEX_FASTQ_R1="${FLEX_FASTQ_R1:-/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz}"
export FLEX_FASTQ_R2="${FLEX_FASTQ_R2:-/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz}"
export FLEX_FIXTURE_ROOT="${FLEX_FIXTURE_ROOT:-/mnt/pikachu/STAR-Flex/reference/tests/100K}"
export FLEX_INLINE_TMP_DIR="${FLEX_INLINE_TMP_DIR:-/tmp/flex_inline_test}"

# SLAM fixture (external, not in repo)
export SLAM_FIXTURE_ROOT="${SLAM_FIXTURE_ROOT:-/mnt/pikachu/STAR-Flex/test/tmp_slam_fixture}"
export SLAM_FIXTURE_FASTQ="${SLAM_FIXTURE_FASTQ:-/mnt/pikachu/STAR-Flex/test/tmp_slam_fixture/SRR32576116.fastq.gz}"
export SLAM_FIXTURE_REF_TSV="${SLAM_FIXTURE_REF_TSV:-/mnt/pikachu/STAR-Flex/test/tmp_slam_fixture/fixture_ref_human.tsv.gz}"
export SLAM_FIXTURE_STAR_INDEX="${SLAM_FIXTURE_STAR_INDEX:-/mnt/pikachu/STAR-Flex/test/fixtures/slam/ref/STAR-index}"
export SLAM_FIXTURE_SNPS_BED="${SLAM_FIXTURE_SNPS_BED:-/mnt/pikachu/STAR-Flex/test/tmp_slam_fixture/snps.bed}"

# Export generic SLAM variables used by fixture scripts
export FASTQ="${FASTQ:-$SLAM_FIXTURE_FASTQ}"
export STAR_INDEX="${STAR_INDEX:-$SLAM_FIXTURE_STAR_INDEX}"
export SNPS_BED="${SNPS_BED:-$SLAM_FIXTURE_SNPS_BED}"
export REF_TSV="${REF_TSV:-$SLAM_FIXTURE_REF_TSV}"

# Y-chrom bulk dataset
export YCHROM_BULK_GENOME_DIR="${YCHROM_BULK_GENOME_DIR:-/storage/flex_filtered_reference/star_index}"
export YCHROM_BULK_R1="${YCHROM_BULK_R1:-/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R1_001.fastq.gz}"
export YCHROM_BULK_R2="${YCHROM_BULK_R2:-/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R2_001.fastq.gz}"

# Y-chrom flex dataset
export YCHROM_FLEX_FASTQ_BASE="${YCHROM_FLEX_FASTQ_BASE:-/storage/downsampled/SC2300771}"
