./demux_bam \
    --bam /storage/scRNAseq_output/Alignments/SC2300771/star/Aligned.sortedByCoord.out.bam \
    --outdir test_output              \
    --sample_probes tables/probe-barcodes-fixed-rna-profiling-rna.txt \
    --probe_offset 68                   \
    --threads 4                         \
    --hts_threads 4                    \
    --min_mapq 20                       \
    --ub_tag UB                         \
    -v
