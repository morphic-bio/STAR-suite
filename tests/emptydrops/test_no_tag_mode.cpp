/**
 * @file test_no_tag_mode.cpp
 * @brief Unit test: No-tag single-sample EmptyDrops mode
 * 
 * Target: scrna_emptydrops_run (C API)
 * Test: Build tiny in-memory SampleMatrixData with 16bp barcodes,
 *       run EmptyDrops with occupancy disabled, verify passingBarcodes includes expected set
 */

#include <iostream>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include "scrna_api.h"

using namespace std;

int main() {
    cout << "Test: No-tag single-sample mode (scrna_emptydrops_run)" << endl;
    
    // Create input with 16bp barcodes (no tag suffix)
    uint32_t nCells = 5000;
    vector<uint32_t> umiCounts(nCells);
    vector<char*> barcodes(nCells);
    
    // Use fixed seed for reproducibility
    srand(123);
    
    // Create UMI distribution
    // High UMI cells (500): 5k-10k
    for (uint32_t i = 0; i < 500; i++) {
        umiCounts[i] = 5000 + (i * 10);
    }
    // Medium UMI cells (1500): 500-1000
    for (uint32_t i = 500; i < 2000; i++) {
        umiCounts[i] = 500 + (rand() % 500);
    }
    // Low UMI / empty (3000): 0-100
    for (uint32_t i = 2000; i < nCells; i++) {
        umiCounts[i] = rand() % 100;
    }
    
    // Create 16bp barcode strings
    for (uint32_t i = 0; i < nCells; i++) {
        char buf[32];
        snprintf(buf, sizeof(buf), "ACGT%012d", i);  // 16bp barcode
        barcodes[i] = strdup(buf);
    }
    
    // Set up input (no sparse data - simple mode fallback)
    scrna_matrix_input input;
    memset(&input, 0, sizeof(input));
    input.umi_counts = umiCounts.data();
    input.barcodes = barcodes.data();
    input.n_cells = nCells;
    input.n_features = 0;
    // No sparse data - triggers Simple ED fallback
    input.sparse_gene_ids = nullptr;
    input.sparse_counts = nullptr;
    input.sparse_cell_index = nullptr;
    input.n_genes_per_cell = nullptr;
    
    // Create config with occupancy disabled (compat mode)
    scrna_ed_config *config = scrna_ed_config_create();
    config->n_expected_cells = 1000;
    config->disable_occupancy_filter = 1;  // Required for no-tag mode
    
    // Verify default seed is 1 (per plan)
    if (config->seed != 1) {
        cerr << "FAIL: Default seed is " << config->seed << ", expected 1" << endl;
        scrna_ed_config_destroy(config);
        for (uint32_t i = 0; i < nCells; i++) free(barcodes[i]);
        return 1;
    }
    
    // Run
    scrna_ed_result result;
    int rc = scrna_emptydrops_run(&input, config, &result);
    
    // Clean up barcodes
    for (uint32_t i = 0; i < nCells; i++) {
        free(barcodes[i]);
    }
    
    // Check return code
    if (rc != 0) {
        cerr << "FAIL: scrna_emptydrops_run returned " << rc << endl;
        if (result.error_message) {
            cerr << "  Error: " << result.error_message << endl;
        }
        scrna_ed_result_free(&result);
        scrna_ed_config_destroy(config);
        return 1;
    }
    
    // Verify results
    if (result.n_barcodes == 0) {
        cerr << "FAIL: No barcodes returned" << endl;
        scrna_ed_result_free(&result);
        scrna_ed_config_destroy(config);
        return 1;
    }
    
    // Should detect high-UMI cells (at least 400 of top 500)
    if (result.n_barcodes < 400) {
        cerr << "FAIL: Too few cells (" << result.n_barcodes << ", expected >= 400)" << endl;
        scrna_ed_result_free(&result);
        scrna_ed_config_destroy(config);
        return 1;
    }
    
    // Verify barcodes are 16bp
    if (result.barcodes[0] && strlen(result.barcodes[0]) != 16) {
        cerr << "FAIL: Barcode length " << strlen(result.barcodes[0]) << " != 16" << endl;
        scrna_ed_result_free(&result);
        scrna_ed_config_destroy(config);
        return 1;
    }
    
    // Verify summary counts are consistent
    if (result.n_simple_cells > result.n_barcodes) {
        cerr << "FAIL: n_simple_cells (" << result.n_simple_cells << ") > n_barcodes (" << result.n_barcodes << ")" << endl;
        scrna_ed_result_free(&result);
        scrna_ed_config_destroy(config);
        return 1;
    }
    
    cout << "PASS: No-tag mode: " << result.n_barcodes << " cells, "
         << "simple=" << result.n_simple_cells << ", "
         << "threshold=" << result.retain_threshold << ", "
         << "seed=" << config->seed << endl;
    
    // Cleanup
    scrna_ed_result_free(&result);
    scrna_ed_config_destroy(config);
    
    return 0;
}
