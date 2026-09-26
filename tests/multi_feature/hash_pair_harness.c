/* Synthetic counts exercise the native output function without FASTQ/dedup noise. */
#ifndef HASH_IMPL
#define HASH_IMPL "../../core/features/process_features/src/hash_demux.c"
#endif
#include HASH_IMPL

int main(int argc, char **argv) {
    if (argc != 4) return 2; /* outdir, method, pair table */
    char *ids[] = {"B", "A", "C", "D"}; /* tie ordering differs from ID ordering */
    char *types[] = {"Multiplexing Capture", "Multiplexing Capture", "Multiplexing Capture", "Multiplexing Capture"};
    char *barcodes[] = {"AAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCC", "GGGGGGGGGGGGGGGG", "TTTTTTTTTTTTTTTT", "ACGTACGTACGTACGT", "TGCATGCATGCATGCA"};
    unsigned int values[][4] = {{40,40,0,0}, {90,10,4,0}, {10,8,6,0}, {0,0,0,0}, {10,0,9,0}, {20,10,5,0}};
    feature_arrays features = {0};
    features.number_of_features = 4;
    features.feature_names = ids; features.feature_ids = ids; features.feature_types = types;
    pf_counts_result counts = {0};
    counts.barcode_to_deduped_hash = kh_init(u32ptr);
    initseq2Code();
    barcode_length = 16;
    char path[4096]; snprintf(path, sizeof(path), "%s/barcodes.txt", argv[1]);
    FILE *fp = fopen(path,"w"); if (!fp) return 2;
    for (int i=0; i<6; ++i) {
        fprintf(fp,"%s\n",barcodes[i]);
        unsigned char code[8] = {0}; string2code(barcodes[i],16,code);
        uint32_t key; memcpy(&key,code,4);
        int absent;
        khint_t k = kh_put(u32ptr, counts.barcode_to_deduped_hash, key, &absent);
        khash_t(u32u32) *h = kh_init(u32u32);
        kh_val(counts.barcode_to_deduped_hash,k) = h;
        for (int j=0;j<4;++j) {
            khint_t f=kh_put(u32u32,h,j+1,&absent); kh_val(h,f)=values[i][j];
        }
    }
    fclose(fp);
    pf_hash_mex_config cfg = {0};
    cfg.base.assign_output_dir=argv[1]; cfg.base.mex_output_dir=argv[1];
    cfg.base.features=&features; cfg.base.counts=&counts;
    cfg.hash_demux_method=argv[2]; cfg.hash_min_total=3; cfg.hash_min_top=3; cfg.hash_min_ratio=2.0;
    cfg.hash_sample_table=argv[3]; cfg.hash_min_pair_ratio=2.0;
    unsigned char mask[]={1,1,1,1}; pf_hash_demux_stats stats;
    int result = write_hash_demux_outputs(&cfg,mask,4,&stats);
    for (khint_t k=kh_begin(counts.barcode_to_deduped_hash);k!=kh_end(counts.barcode_to_deduped_hash);++k)
        if(kh_exist(counts.barcode_to_deduped_hash,k)) kh_destroy(u32u32,(khash_t(u32u32)*)kh_val(counts.barcode_to_deduped_hash,k));
    kh_destroy(u32ptr,counts.barcode_to_deduped_hash);
    return result == 0 ? 0 : 1;
}
