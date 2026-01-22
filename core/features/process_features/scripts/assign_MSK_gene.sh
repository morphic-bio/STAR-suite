#!/bin/bash
./assignBarcodes -v  -w  /storage/scRNAseq_output/whitelists/3M-february-2018_NXT.txt -o 20 -m 1 -s 1 -i 0 -u 12 --limit_search 2 -c 4   --max_barcode_mismatches 1 -f /mnt/pikachu/ref_feature_geneBC.csv --translate_NXT -d /storage/gene_features  \
/storage/gene/30_KO_ES/ \
/storage/gene/30_KO_S6_1 \
/storage/gene/30_KO_S5_2 \
/storage/gene/30_KO_S5_1 \
/storage/gene/30_KO_PP1 \
/storage/gene/30_KO_S6_2 \
/storage/gene/30_KO_DE_XM \
/storage/gene/30_KO_PP2 