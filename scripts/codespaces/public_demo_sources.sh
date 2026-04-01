#!/usr/bin/env bash
# Public source pins for Codespaces demos.

export CODESPACES_PUBLIC_BULK_SRR="SRR4422207"
export CODESPACES_PUBLIC_BULK_GSE="GSE88509"
export CODESPACES_PUBLIC_BULK_GSM="GSM2344101"
export CODESPACES_PUBLIC_BULK_MAX_SPOTS="25000"

export CODESPACES_PUBLIC_SLAM_SRR="SRR32576116"

# Public human scRNA source used for the chromosome-limited Codespaces demo
# fixture. This 10x dataset is explicitly a healthy male donor, so chrY is
# meaningful and reproducible for the demo surface.
export CODESPACES_PUBLIC_10X_MALE_DATASET_URL="https://www.10xgenomics.com/datasets/human-b-cells-from-a-healthy-donor-1-k-cells-2-standard-6-0-0"
export CODESPACES_PUBLIC_10X_MALE_FASTQS_TAR_URL="https://cf.10xgenomics.com/samples/cell-vdj/6.0.0/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar"
export CODESPACES_PUBLIC_10X_MALE_FASTQS_TAR_BYTES="3401635840"
# Assay-specific barcode geometry for this public source. These are demo-source
# defaults, not global STAR-suite assumptions.
export CODESPACES_PUBLIC_10X_MALE_SOLO_CB_START="1"
export CODESPACES_PUBLIC_10X_MALE_SOLO_CB_LEN="16"
export CODESPACES_PUBLIC_10X_MALE_SOLO_UMI_START="17"
export CODESPACES_PUBLIC_10X_MALE_SOLO_UMI_LEN="10"

# Public human reference source for the demo mini-reference.
export CODESPACES_PUBLIC_HUMAN_ENSEMBL_RELEASE="110"
export CODESPACES_PUBLIC_HUMAN_CHR22_FASTA_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"
export CODESPACES_PUBLIC_HUMAN_CHRY_FASTA_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz"
export CODESPACES_PUBLIC_HUMAN_GTF_URL="https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"

# Public assay-specific replacements for the current synthetic single-cell
# overlays. These are official 10x dataset assets and configs. The raw FASTQ
# tars are large, but they are uncompressed tar archives, so a future Codespaces
# path can stream selected entries without downloading the full archive.

# A375 5' GEM-X CRISPR (1K cells, 11 guides, ~4.7 GB tar, CC BY 4.0)
export CODESPACES_PUBLIC_A375_5P_DATASET_URL="https://www.10xgenomics.com/datasets/1k-CRISPR-5p-gemx"
export CODESPACES_PUBLIC_A375_5P_FASTQS_TAR_URL="https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/1k_CRISPR_5p_gemx_Multiplex/1k_CRISPR_5p_gemx_Multiplex_fastqs.tar"
export CODESPACES_PUBLIC_A375_5P_FEATURE_REF_URL="https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/1k_CRISPR_5p_gemx_Multiplex/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv"

# A375 3' v3 NextGEM CRISPR (10K cells, 2 guides, ~59 GB tar, CC BY 4.0)
export CODESPACES_PUBLIC_A375_3P_DATASET_URL="https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/SC3_v3_NextGem_DI_CRISPR_10K"
export CODESPACES_PUBLIC_A375_3P_FASTQS_TAR_URL="https://cg.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_DI_CRISPR_10K/SC3_v3_NextGem_DI_CRISPR_10K_fastqs.tar"
export CODESPACES_PUBLIC_A375_3P_FEATURE_REF_URL="https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_DI_CRISPR_10K/SC3_v3_NextGem_DI_CRISPR_10K_feature_ref.csv"

# Perturb / CRISPR Guide Capture
export CODESPACES_PUBLIC_PERTURB_DATASET_URL="https://www.10xgenomics.com/datasets/10k-k562-transduced-with-small-guide-library-5-ht-v2-0-chromium-x-2-standard"
export CODESPACES_PUBLIC_PERTURB_FASTQS_TAR_URL="https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/10k_K562sl_5pv2_nextgem_CRISPR_Chromium_X_Multiplex/10k_K562sl_5pv2_nextgem_CRISPR_Chromium_X_Multiplex_fastqs.tar"
export CODESPACES_PUBLIC_PERTURB_FASTQS_TAR_BYTES="32199505920"
export CODESPACES_PUBLIC_PERTURB_CONFIG_URL="https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/10k_K562sl_5pv2_nextgem_CRISPR_Chromium_X_Multiplex/10k_K562sl_5pv2_nextgem_CRISPR_Chromium_X_Multiplex_config.csv"
export CODESPACES_PUBLIC_PERTURB_FEATURE_REF_URL="https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/10k_K562sl_5pv2_nextgem_CRISPR_Chromium_X_10k_K562sl_5pv2_nextgem_CRISPR_Chromium_X/10k_K562sl_5pv2_nextgem_CRISPR_Chromium_X_10k_K562sl_5pv2_nextgem_CRISPR_Chromium_X_count_feature_reference.csv"

# Flex / Fixed RNA Profiling
export CODESPACES_PUBLIC_FLEX_DATASET_URL="https://www.10xgenomics.com/datasets/10k-human-k562-r-cells-singleplex-sample-1-standard"
export CODESPACES_PUBLIC_FLEX_FASTQS_TAR_URL="https://cf.10xgenomics.com/samples/cell-exp/7.0.0/10k_K562_singleplex_Multiplex/10k_K562_singleplex_Multiplex_fastqs.tar"
export CODESPACES_PUBLIC_FLEX_FASTQS_TAR_BYTES="41699727360"
export CODESPACES_PUBLIC_FLEX_CONFIG_URL="https://cf.10xgenomics.com/samples/cell-exp/7.0.0/10k_K562_singleplex_Multiplex/10k_K562_singleplex_Multiplex_config.csv"
export CODESPACES_PUBLIC_FLEX_PROBE_SET_URL="https://cf.10xgenomics.com/samples/cell-exp/7.0.0/10k_K562_singleplex_10k_K562_singleplex/10k_K562_singleplex_10k_K562_singleplex_count_probe_set.csv"
