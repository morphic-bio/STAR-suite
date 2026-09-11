/* Differential integration fixture for assignment counts, gather, connected
 * components and tie handling. Does not inspect the per-UMI representation. */
#include "memory.h"
#include "barcode_match.h"
#include "pf_counts.h"
#include "prototypes.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

static unsigned rng;
static unsigned next_random(void) { rng = rng * 1664525u + 1013904223u; return rng; }
static void add(data_structures *h, memory_pool_collection *p, unsigned bc, unsigned umi, unsigned feature) {
    unsigned char code[4];
    memcpy(code, &bc, 4);
    char seq[13];
    unsigned char packed[3] = {umi & 255u, (umi >> 8) & 255u, (umi >> 16) & 255u};
    code2string(packed, seq, 3);
    update_feature_counts_from_code(code, seq, feature, h, p);
}
int main(void) {
    barcode_match_init(); initialize_complement();
    barcode_length = 16; umi_length = 12; number_of_features = 64; maximum_feature_length = 40;
    initialize_unit_sizes();
    for (unsigned scenario = 0; scenario < 6; ++scenario) {
        data_structures hashes[4] = {0};
        statistics stats[4] = {0};
        memory_pool_collection *pools[4];
        for (int i=0;i<4;++i) { initialize_data_structures(&hashes[i]); pools[i]=initialize_memory_pool_collection(); }
        rng = 9371;
        for (unsigned i=0;i<32000;++i) {
            unsigned n=next_random(), bc=(n>>16)%128, u=next_random() & 0xffffffu;
            unsigned feature=1+(bc%64), thread=i%4;
            add(&hashes[thread],pools[thread],bc,u,feature);
            if (i%3==0) add(&hashes[(thread+1)%4],pools[(thread+1)%4],bc,u,feature);
            if (i%5==0) add(&hashes[(thread+2)%4],pools[(thread+2)%4],bc,u^1u,feature);
            if (i%7==0) add(&hashes[thread],pools[thread],bc,u,1+(feature%64));
        }
        // Many competing features on the same UMI; ties and shared counters must survive gather.
        for (unsigned f=1;f<=64;++f) for(unsigned n=0;n<3;++n)
            add(&hashes[f%4],pools[f%4],129,0xaaaaaa,f);
        for (int i=1;i<4;++i) {
            merge_process_feature_thread_data(&hashes[0],pools[0],&stats[0],&hashes[i],pools[i],&stats[i]);
            destroy_data_structures(&hashes[i]);free_memory_pool_collection(pools[i]);
        }
        unsigned stringency=(scenario%3)*500, minimum=scenario/3*2;
        pf_counts_result *result=pf_build_deduped_counts(&hashes[0],64,stringency,minimum);assert(result);
        for(unsigned bc=0;bc<=129;++bc) {
            khint_t k=kh_get(u32ptr,result->barcode_to_deduped_hash,bc);
            if(k==kh_end(result->barcode_to_deduped_hash))continue;
            khash_t(u32u32)*counts=kh_val(result->barcode_to_deduped_hash,k);
            for(unsigned f=1;f<=64;++f){
                khint_t j=kh_get(u32u32,counts,f);
                if(j!=kh_end(counts))printf("%u\t%u\t%u\t%u\n",scenario,bc,f,kh_val(counts,j));
            }
        }
        pf_counts_result_free(result);destroy_data_structures(&hashes[0]);free_memory_pool_collection(pools[0]);
    }
    return 0;
}
