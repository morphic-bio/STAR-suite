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

# A375 5' GEM-X downloadable fixture (scripts/download_a375_5prime_fixture.sh)
export A375_5P_FIXTURE_ROOT="${A375_5P_FIXTURE_ROOT:-/tmp/a375_5prime_fixture}"
export A375_5P_GEX_DIR="${A375_5P_GEX_DIR:-$A375_5P_FIXTURE_ROOT/fastqs/gex}"
export A375_5P_CRISPR_DIR="${A375_5P_CRISPR_DIR:-$A375_5P_FIXTURE_ROOT/fastqs/crispr}"
export A375_5P_FEATURE_REF="${A375_5P_FEATURE_REF:-$A375_5P_FIXTURE_ROOT/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv}"
export A375_5P_DOWNSAMPLE_DIR="${A375_5P_DOWNSAMPLE_DIR:-$A375_5P_FIXTURE_ROOT/downsampled_100000}"

# A375 3' v3 NextGEM downloadable fixture (scripts/download_a375_3prime_fixture.sh)
export A375_3P_FIXTURE_ROOT="${A375_3P_FIXTURE_ROOT:-/tmp/a375_3prime_fixture}"
export A375_3P_FEATURE_REF="${A375_3P_FEATURE_REF:-$A375_3P_FIXTURE_ROOT/SC3_v3_NextGem_DI_CRISPR_10K_feature_ref.csv}"
export A375_3P_DOWNSAMPLE_DIR="${A375_3P_DOWNSAMPLE_DIR:-$A375_3P_FIXTURE_ROOT/downsampled_100000}"

# UCSF 2M perturb fixtures (current canonical runs)
export UCSF_2M_ROOT="${UCSF_2M_ROOT:-/storage/ucsf-2M}"
export UCSF_2M_STAR_RUNS_ROOT="${UCSF_2M_STAR_RUNS_ROOT:-$UCSF_2M_ROOT/star_runs}"
export UCSF_2M_FIXTURE_SEQUENTIAL="${UCSF_2M_FIXTURE_SEQUENTIAL:-$UCSF_2M_STAR_RUNS_ROOT/fixture_ucsf2m_current_sequential}"
export UCSF_2M_FIXTURE_DYNAMIC="${UCSF_2M_FIXTURE_DYNAMIC:-$UCSF_2M_STAR_RUNS_ROOT/fixture_ucsf2m_current_dynamic}"
export UCSF_2M_FIXTURE_COMPARE_REPORT="${UCSF_2M_FIXTURE_COMPARE_REPORT:-/tmp/ucsf2m_seq_vs_fixture_20260224_090754/COMPARE_REPORT.txt}"
export UCSF_100K_PFCONFIG_ROOT="${UCSF_100K_PFCONFIG_ROOT:-$UCSF_2M_ROOT/fixtures/ucsf2m_iPSC2_AALG2_100k_pfconfig}"
export UCSF_100K_GEX_DIR="${UCSF_100K_GEX_DIR:-$UCSF_100K_PFCONFIG_ROOT/GEX/iPSC2_1_AALG2}"
export UCSF_100K_GUIDE_DIR="${UCSF_100K_GUIDE_DIR:-$UCSF_100K_PFCONFIG_ROOT/guides/iPSC2_1_AALG2}"
export UCSF_100K_PFCONFIG="${UCSF_100K_PFCONFIG:-$UCSF_100K_PFCONFIG_ROOT/multi_config_100k.csv}"
export UCSF_100K_FEATURE_REF="${UCSF_100K_FEATURE_REF:-/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv}"
export UCSF_100K_CB_WHITELIST="${UCSF_100K_CB_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}"
export UCSF_100K_GENOME_DIR="${UCSF_100K_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"

# UCSF full-sample perturb (canonical benchmark surface: corrected EBs2_2).
# Velocyto / run_star_velocyto_canonical profile=2m uses UCSF_2M_*; set UCSF_2M_PFCONFIG locally (no default).
# Legacy iPSC2 under /storage/ucsf-2M is not the default here — override UCSF_2M_GEX_DIR if you still use it.
export UCSF_EBS2_2_ROOT="${UCSF_EBS2_2_ROOT:-/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2}"
export UCSF_2M_GEX_DIR="${UCSF_2M_GEX_DIR:-$UCSF_EBS2_2_ROOT/GEX}"
export UCSF_2M_GUIDE_DIR="${UCSF_2M_GUIDE_DIR:-$UCSF_EBS2_2_ROOT/guides}"
export UCSF_2M_FEATURE_REF="${UCSF_2M_FEATURE_REF:-$UCSF_100K_FEATURE_REF}"
export UCSF_2M_CB_WHITELIST="${UCSF_2M_CB_WHITELIST:-$UCSF_100K_CB_WHITELIST}"
# export UCSF_2M_PFCONFIG=/path/to/full_sample_multi_config.csv

# UCSF full-sample perturb references
export UCSF_FULL_ROOT="${UCSF_FULL_ROOT:-/storage/ucsf-full/bench_20260218_dynamic_first}"
export UCSF_FULL_STAR_RUNS_ROOT="${UCSF_FULL_STAR_RUNS_ROOT:-$UCSF_FULL_ROOT/runs}"
export UCSF_FULL_FIXTURE_STAR="${UCSF_FULL_FIXTURE_STAR:-$UCSF_FULL_STAR_RUNS_ROOT/star_full_iPSC2_1_AALG2_forward_rescue_guides_bootstrap_20260221_010635}"
export UCSF_FULL_DYNAMIC_32X32="${UCSF_FULL_DYNAMIC_32X32:-$UCSF_FULL_STAR_RUNS_ROOT/star_full_dynamic_32x32_20260224_092512}"
export UCSF_FULL_CR_REF="${UCSF_FULL_CR_REF:-$UCSF_FULL_ROOT/cellranger_runs/cr_full_iPSC2_1_AALG2_crstar32_20260218_205804}"

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
export BULK_MULTISAMPLE_GENOME_DIR="${BULK_MULTISAMPLE_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
export BULK_MULTISAMPLE_R1="${BULK_MULTISAMPLE_R1:-$YCHROM_BULK_R1}"
export BULK_MULTISAMPLE_R2="${BULK_MULTISAMPLE_R2:-$YCHROM_BULK_R2}"
export BULK_MULTISAMPLE_PAIRS_PER_SAMPLE="${BULK_MULTISAMPLE_PAIRS_PER_SAMPLE:-5000}"
export BULK_MULTISAMPLE_THREADS="${BULK_MULTISAMPLE_THREADS:-1}"

# Y-chrom flex dataset
export YCHROM_FLEX_FASTQ_BASE="${YCHROM_FLEX_FASTQ_BASE:-/storage/downsampled/SC2300771}"
export FLEX_YREMOVE_THREADS="${FLEX_YREMOVE_THREADS:-1}"
export SOLO_YREMOVE_THREADS="${SOLO_YREMOVE_THREADS:-1}"
export PERTURB_YREMOVE_THREADS="${PERTURB_YREMOVE_THREADS:-1}"

# Public bulk PE regression fixture
# Note: Salmon alignment-mode autodetect (-l A) mis-detects this fixture as ISR.
# Keep the public smoke pinned to IU unless specifically testing Salmon autodetect.
export PUBLIC_BULK_FIXTURE_ROOT="${PUBLIC_BULK_FIXTURE_ROOT:-/tmp/public_bulk_fixture_SRR4422207}"
export PUBLIC_BULK_ACCESSION="${PUBLIC_BULK_ACCESSION:-SRR4422207}"
export PUBLIC_BULK_GEO_SERIES="${PUBLIC_BULK_GEO_SERIES:-GSE88509}"
export PUBLIC_BULK_GEO_SAMPLE="${PUBLIC_BULK_GEO_SAMPLE:-GSM2344101}"
export PUBLIC_BULK_R1="${PUBLIC_BULK_R1:-$PUBLIC_BULK_FIXTURE_ROOT/${PUBLIC_BULK_ACCESSION}_1.fastq.gz}"
export PUBLIC_BULK_R2="${PUBLIC_BULK_R2:-$PUBLIC_BULK_FIXTURE_ROOT/${PUBLIC_BULK_ACCESSION}_2.fastq.gz}"
export PUBLIC_BULK_GENOME_DIR="${PUBLIC_BULK_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
export PUBLIC_BULK_TRANSCRIPTOME="${PUBLIC_BULK_TRANSCRIPTOME:-/storage/autoindex_110_44/bulk_index/transcriptome.fa}"
export PUBLIC_BULK_STAR_MATE_ORDER="${PUBLIC_BULK_STAR_MATE_ORDER:-R2R1}"
export PUBLIC_BULK_SALMON_LIBTYPE="${PUBLIC_BULK_SALMON_LIBTYPE:-IU}"
# With Salmon pinned to IU on this fixture, keep smoke thresholds slightly looser
# than the old auto-detect-based defaults.
export PUBLIC_BULK_EXPECTED_FORMAT="${PUBLIC_BULK_EXPECTED_FORMAT:-ISF}"
export PUBLIC_BULK_MAX_DROPPED_INCOMPAT="${PUBLIC_BULK_MAX_DROPPED_INCOMPAT:--1}"
