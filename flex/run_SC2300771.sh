#!/bin/bash
# Production script for STAR-Flex pipeline
# Uses --flex yes to enable the inline hash pipeline with sample detection and FlexFilter

# Optional debug/trace environment variables:
#export STAR_INLINE_REJECT_LOG=/storage/trace/inline_resolver_BC004.tsv
#export STAR_INLINE_TRACE_QNAME=1    # optional, for qnames
#export STAR_DEBUG_TAG=1             # optional, for UMI correction debug

OMP_NUM_THREADS=24
OUTPUT_DIR="/storage/production/SC2300771"
TMP_DIR="/storage/production/tmp/SC2300771"
fastq_dir="/storage/JAX_sequences"
rm -rf "$TMP_DIR"
mkdir -p "$OUTPUT_DIR"

/mnt/pikachu/STAR-Flex/source/STAR \
  --runThreadN 24 \
  --outTmpDir "$TMP_DIR" \
  --genomeDir /storage/flex-index-44 \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
  --soloCBwhitelist /storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt \
  --flex yes \
  --soloFlexExpectedCellsPerTag 3000 \
  --soloSampleWhitelist /storage/SC2300771_whitelist.tsv \
  --soloSampleProbes /mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt \
  --soloSampleProbeOffset 68 \
  --soloFlexOutputPrefix "$OUTPUT_DIR/per_sample" \
  --limitIObufferSize 50000000 50000000 \
  --outSJtype None \
  --outBAMcompression 6 \
  --soloMultiMappers Rescue \
  --alignIntronMax 500000 \
  --outFilterMismatchNmax 6 \
  --outFilterMismatchNoverReadLmax 1.0 \
  --outFilterMatchNmin 25 \
  --outSAMunmapped None \
  --outFilterMatchNminOverLread 0 \
  --outFilterMultimapNmax 10000 \
  --outFilterMultimapScoreRange 4 \
  --outSAMmultNmax 10000 \
  --winAnchorMultimapNmax 200 \
  --outSAMprimaryFlag AllBestScore \
  --outFilterScoreMin 0 \
  --outFilterScoreMinOverLread 0 \
  --outSAMattributes NH HI AS nM NM GX GN \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloCellFilter None \
  --clipAdapterType CellRanger4 \
  --soloFeatures Gene \
  --alignEndsType Local \
  --soloStrand Unstranded \
  --chimSegmentMin 1000000 \
  --soloSampleSearchNearby no \
  --outSAMtype BAM Unsorted \
  --readFilesCommand zcat \
  --readFilesIn $fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R2_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R2_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R2_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R2_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R2_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R2_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R2_001.fastq.gz $fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R1_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R1_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R1_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R1_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R1_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R1_001.fastq.gz,$fastq_dir/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R1_001.fastq.gz \
  --outFileNamePrefix "$OUTPUT_DIR/"

echo ""
echo "=== Run complete ==="
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Check logs:"
echo "  $OUTPUT_DIR/Log.out"
echo "  $OUTPUT_DIR/Log.final.out"
echo ""
echo "Solo output:"
echo "  $OUTPUT_DIR/Solo.out/Gene/raw/"
echo ""
echo "FlexFilter output (per-sample):"
echo "  $OUTPUT_DIR/per_sample/"

