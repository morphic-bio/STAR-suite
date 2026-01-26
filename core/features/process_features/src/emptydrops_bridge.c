/**
 * @file emptydrops_bridge.c
 * @brief Bridge between pf_counts_result and EmptyDrops API
 */

#include "../include/emptydrops_bridge.h"
#include "../include/pf_api.h"
#include "../include/globals.h"
#include "../include/prototypes.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ============================================================================
 * EmptyDrops Input Construction
 * ============================================================================ */

emptydrops_input* emptydrops_input_from_counts(
    pf_counts_result *counts,
    char **feature_names,
    int do_translate_nxt
) {
    if (!counts || !counts->barcode_to_deduped_hash) {
        return NULL;
    }
    
    khash_t(u32ptr) *hash = counts->barcode_to_deduped_hash;
    uint32_t n_barcodes = kh_size(hash);
    
    if (n_barcodes == 0) {
        fprintf(stderr, "[emptydrops_bridge] No barcodes in counts\n");
        return NULL;
    }
    
    /* Allocate result structure */
    emptydrops_input *input = calloc(1, sizeof(emptydrops_input));
    if (!input) return NULL;
    
    input->n_barcodes = n_barcodes;
    input->n_features = counts->n_features;
    
    /* Allocate barcode and UMI arrays */
    input->barcodes = calloc(n_barcodes, sizeof(char*));
    input->umi_counts = calloc(n_barcodes, sizeof(uint32_t));
    input->sparse_cell_index = calloc(n_barcodes + 1, sizeof(uint32_t));
    input->n_genes_per_cell = calloc(n_barcodes, sizeof(uint32_t));
    
    if (!input->barcodes || !input->umi_counts || 
        !input->sparse_cell_index || !input->n_genes_per_cell) {
        emptydrops_input_free(input);
        return NULL;
    }
    
    /* First pass: count total non-zero entries and fill per-cell arrays */
    size_t nnz = 0;
    uint32_t bc_idx = 0;
    khint_t k;
    
    for (k = kh_begin(hash); k != kh_end(hash); ++k) {
        if (!kh_exist(hash, k)) continue;
        
        khash_t(u32u32) *deduped = (khash_t(u32u32)*)kh_val(hash, k);
        uint32_t n_genes = kh_size(deduped);
        
        /* Count UMIs for this barcode */
        uint32_t umi_sum = 0;
        khint_t kf;
        for (kf = kh_begin(deduped); kf != kh_end(deduped); ++kf) {
            if (!kh_exist(deduped, kf)) continue;
            umi_sum += kh_val(deduped, kf);
        }
        
        input->umi_counts[bc_idx] = umi_sum;
        input->n_genes_per_cell[bc_idx] = n_genes;
        input->sparse_cell_index[bc_idx] = nnz;
        nnz += n_genes;
        bc_idx++;
    }
    input->sparse_cell_index[n_barcodes] = nnz;
    input->sparse_nnz = nnz;
    
    /* Allocate sparse arrays */
    if (nnz > 0) {
        input->sparse_gene_ids = calloc(nnz, sizeof(uint32_t));
        input->sparse_counts = calloc(nnz, sizeof(uint32_t));
        if (!input->sparse_gene_ids || !input->sparse_counts) {
            emptydrops_input_free(input);
            return NULL;
        }
    }
    
    /* Second pass: fill barcode strings and sparse data */
    bc_idx = 0;
    size_t sparse_idx = 0;
    
    for (k = kh_begin(hash); k != kh_end(hash); ++k) {
        if (!kh_exist(hash, k)) continue;
        
        uint32_t barcode_key = kh_key(hash, k);
        khash_t(u32u32) *deduped = (khash_t(u32u32)*)kh_val(hash, k);
        
        /* Decode barcode to string */
        char *barcode_str = calloc(barcode_length + 1, 1);
        if (!barcode_str) {
            emptydrops_input_free(input);
            return NULL;
        }
        code2string((unsigned char*)&barcode_key, barcode_str, barcode_code_length);
        
        /* Apply NXT translation if requested */
        if (do_translate_nxt && translate_NXT) {
            translate_nxt_inplace(barcode_str, barcode_length);
        }
        
        input->barcodes[bc_idx] = barcode_str;
        
        /* Fill sparse entries for this barcode */
        khint_t kf;
        for (kf = kh_begin(deduped); kf != kh_end(deduped); ++kf) {
            if (!kh_exist(deduped, kf)) continue;
            
            uint32_t feature_id = kh_key(deduped, kf);
            uint32_t count = kh_val(deduped, kf);
            
            /* Convert 1-based feature_id to 0-based gene index */
            input->sparse_gene_ids[sparse_idx] = feature_id - 1;
            input->sparse_counts[sparse_idx] = count;
            sparse_idx++;
        }
        
        bc_idx++;
    }
    
    /* Copy feature names (shallow copy - we don't own these) */
    input->features = feature_names;  /* Borrowed reference */
    
    return input;
}

void emptydrops_input_free(emptydrops_input *input) {
    if (!input) return;
    
    /* Free barcode strings (we allocated these) */
    if (input->barcodes) {
        for (uint32_t i = 0; i < input->n_barcodes; i++) {
            free(input->barcodes[i]);
        }
        free(input->barcodes);
    }
    
    /* Note: features is a borrowed reference, don't free */
    
    free(input->umi_counts);
    free(input->sparse_gene_ids);
    free(input->sparse_counts);
    free(input->sparse_cell_index);
    free(input->n_genes_per_cell);
    free(input);
}

/* ============================================================================
 * Main EmptyDrops Entry Point
 * ============================================================================ */

khash_t(strptr)* run_emptydrops_on_counts(
    pf_counts_result *counts,
    char **feature_names,
    int n_features,
    const char *output_dir,
    int n_expected_cells,
    int use_fdr_gate,
    int do_translate_nxt
) {
    if (!counts || !output_dir) {
        fprintf(stderr, "[emptydrops_bridge] Invalid arguments\n");
        return NULL;
    }
    
    /* Validate n_features matches counts */
    if (n_features != counts->n_features) {
        fprintf(stderr, "[emptydrops_bridge] WARNING: n_features mismatch (config=%d, counts=%d), using counts->n_features\n",
                n_features, counts->n_features);
    }
    
    /* Convert counts to EmptyDrops input format */
    fprintf(stderr, "[emptydrops_bridge] Converting counts to EmptyDrops format...\n");
    emptydrops_input *ed_input = emptydrops_input_from_counts(
        counts, feature_names, do_translate_nxt);
    
    if (!ed_input) {
        fprintf(stderr, "[emptydrops_bridge] Failed to convert counts\n");
        return NULL;
    }
    
    fprintf(stderr, "[emptydrops_bridge] Converted %u barcodes, %zu non-zero entries\n",
            ed_input->n_barcodes, ed_input->sparse_nnz);
    
    /* Call EmptyDrops */
    char **filtered_barcodes = NULL;
    uint32_t n_filtered = 0;
    
    pf_error err = pf_run_emptydrops_premex(
        ed_input->umi_counts,
        (const char**)ed_input->barcodes,
        ed_input->n_barcodes,
        (const char**)ed_input->features,
        ed_input->n_features,
        ed_input->sparse_gene_ids,
        ed_input->sparse_counts,
        ed_input->sparse_cell_index,
        ed_input->n_genes_per_cell,
        output_dir,
        n_expected_cells,
        use_fdr_gate,
        &filtered_barcodes,
        &n_filtered
    );
    
    emptydrops_input_free(ed_input);
    
    if (err != PF_OK) {
        fprintf(stderr, "[emptydrops_bridge] EmptyDrops failed with error %d\n", err);
        return NULL;
    }
    
    fprintf(stderr, "[emptydrops_bridge] EmptyDrops returned %u filtered barcodes\n", n_filtered);
    
    /* Convert filtered_barcodes array to khash_t(strptr)* */
    khash_t(strptr) *result_hash = kh_init(strptr);
    if (!result_hash) {
        /* Free filtered_barcodes */
        for (uint32_t i = 0; i < n_filtered; i++) {
            free(filtered_barcodes[i]);
        }
        free(filtered_barcodes);
        return NULL;
    }
    
    int ret;
    for (uint32_t i = 0; i < n_filtered; i++) {
        /* The hash takes ownership of the barcode string */
        khint_t k = kh_put(strptr, result_hash, filtered_barcodes[i], &ret);
        kh_val(result_hash, k) = NULL;  /* No associated value needed */
    }
    
    /* Free the array but not the strings (now owned by hash) */
    free(filtered_barcodes);
    
    return result_hash;
}
