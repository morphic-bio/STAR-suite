/**
 * @file scrna_api.cpp
 * @brief C ABI wrapper implementation for libscrna
 * 
 * EmptyDrops filtering for no-tag/single-sample mode:
 * - Default seed=1 (per emptydrops_refactor_plan.md)
 * - Simple ED is FALLBACK only (runs when no sparse data or explicit --use-simple-empty-drops)
 * - Full Monte Carlo EmptyDrops is the primary filtering method
 * - Occupancy filter disabled (Flex-only feature)
 * - Outputs: filtered_barcodes.txt + EmptyDrops/emptydrops_results.tsv
 */

#include "scrna_api.h"
#include "SampleMatrixData.h"
#include "OrdMagStage.h"
#include "EmptyDropsMultinomial.h"
#include "OccupancyGuard.h"
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>
#include <utility>

using namespace std;

// ============================================================================
// Helper functions
// ============================================================================

static char* strdup_safe(const char *s) {
    if (!s) return nullptr;
    size_t len = strlen(s) + 1;
    char *dup = (char*)malloc(len);
    if (dup) memcpy(dup, s, len);
    return dup;
}

static void mkdir_p(const char *path) {
    char tmp[4096];
    char *p = nullptr;
    size_t len;
    
    snprintf(tmp, sizeof(tmp), "%s", path);
    len = strlen(tmp);
    if (tmp[len - 1] == '/') {
        tmp[len - 1] = 0;
    }
    
    for (p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = 0;
            mkdir(tmp, 0755);
            *p = '/';
        }
    }
    mkdir(tmp, 0755);
}

// ============================================================================
// Config functions
// ============================================================================

extern "C" scrna_ed_config* scrna_ed_config_create(void) {
    scrna_ed_config *config = (scrna_ed_config*)calloc(1, sizeof(scrna_ed_config));
    if (!config) return nullptr;
    
    // Set defaults matching CR/Flex (from FlexFilter::getExternalEmptyDropsDefaults)
    config->n_expected_cells = 3000;
    config->max_percentile = 0.99;
    config->max_min_ratio = 10.0;
    config->umi_min = 500;
    config->umi_min_frac_median = 0.01;
    config->cand_max_n = 20000;
    config->ind_min = 45000;
    config->ind_max = 90000;
    
    config->sim_n = 10000;
    config->fdr = 0.01;
    config->raw_pvalue_threshold = 0.05;
    config->seed = 1;                    // Default seed=1 per emptydrops_refactor_plan.md
    config->lower_testing_bound = 500;   // R's umi.min: cells with UMI <= 500 excluded from testing
    config->ambient_umi_max = 100;       // R's lower: cells with UMI <= 100 used for ambient
    config->mc_threads = 0;              // Single-threaded by default
    
    config->disable_occupancy_filter = 1; // Disabled for compat mode (occupancy is Flex-only)
    config->ed_retain_count = 120000;     // Default retain window
    config->use_fdr_gate = 0;             // Use raw p-value by default (Flex behavior)
    config->use_simple_emptydrops = 0;    // Simple ED is fallback only by default
    
    return config;
}

extern "C" void scrna_ed_config_destroy(scrna_ed_config *config) {
    free(config);
}

// ============================================================================
// Result functions
// ============================================================================

extern "C" void scrna_ed_result_free(scrna_ed_result *result) {
    if (!result) return;
    
    // Free barcodes
    if (result->barcodes) {
        for (size_t i = 0; i < result->n_barcodes; i++) {
            free(result->barcodes[i]);
        }
        free(result->barcodes);
    }
    
    // Free candidates
    if (result->candidates) {
        for (size_t i = 0; i < result->n_candidates; i++) {
            free(result->candidates[i].barcode);
        }
        free(result->candidates);
    }
    
    // Free error message
    free(result->error_message);
    
    // Zero out struct
    memset(result, 0, sizeof(scrna_ed_result));
}

// ============================================================================
// Main EmptyDrops function - Aligned with Flex pipeline
// ============================================================================

extern "C" int scrna_emptydrops_run(
    const scrna_matrix_input *input,
    const scrna_ed_config *config,
    scrna_ed_result *result
) {
    if (!input || !config || !result) {
        return -1;
    }
    
    memset(result, 0, sizeof(scrna_ed_result));
    
    if (input->n_cells == 0) {
        result->error_message = strdup_safe("No cells in input");
        return -1;
    }
    
    cerr << "[scrna_api] Input: " << input->n_cells << " cells, " << input->n_features << " features" << endl;
    
    // Convert input to C++ vectors
    vector<uint32_t> nUMIperCB(input->umi_counts, input->umi_counts + input->n_cells);
    vector<string> barcodes(input->n_cells);
    for (uint32_t i = 0; i < input->n_cells; i++) {
        barcodes[i] = input->barcodes[i] ? input->barcodes[i] : "";
    }
    
    // ========================================================================
    // Step 1: Sort cells by UMI descending and apply retain window
    // ========================================================================
    vector<pair<uint32_t, uint32_t>> umiIdx;  // (umi, original_index)
    umiIdx.reserve(input->n_cells);
    for (uint32_t i = 0; i < input->n_cells; i++) {
        umiIdx.push_back({nUMIperCB[i], i});
    }
    stable_sort(umiIdx.begin(), umiIdx.end(), [](const pair<uint32_t,uint32_t>& a, const pair<uint32_t,uint32_t>& b) {
        return a.first > b.first;
    });
    
    // Apply retain window (top N by UMI, capped by ed_retain_count)
    uint32_t retainCount = (config->ed_retain_count > 0) 
        ? min(config->ed_retain_count, input->n_cells) 
        : input->n_cells;
    
    vector<uint32_t> retainIndices;
    retainIndices.reserve(retainCount);
    for (uint32_t i = 0; i < retainCount && i < umiIdx.size(); i++) {
        retainIndices.push_back(umiIdx[i].second);
    }
    
    cerr << "[scrna_api] Retain window: " << retainIndices.size() << " cells (of " << input->n_cells << " total)" << endl;
    
    // Build mapping from original index to retain index
    vector<uint32_t> origToRetainIdx(input->n_cells, UINT32_MAX);
    for (size_t ri = 0; ri < retainIndices.size(); ri++) {
        origToRetainIdx[retainIndices[ri]] = ri;
    }
    
    // ========================================================================
    // Step 2: Build retainUMI for threshold calculation
    // ========================================================================
    vector<uint32_t> retainUMI;
    retainUMI.reserve(retainIndices.size());
    for (uint32_t idx : retainIndices) {
        retainUMI.push_back(nUMIperCB[idx]);
    }
    
    // ========================================================================
    // Step 3: Check if we should use Simple ED (fallback or explicit request)
    // Simple ED runs ONLY as fallback (no sparse data) or if explicitly requested
    // ========================================================================
    bool useSimpleED = config->use_simple_emptydrops || 
                       !input->sparse_gene_ids || 
                       !input->sparse_counts || 
                       !input->sparse_cell_index;
    
    if (useSimpleED) {
        cerr << "[scrna_api] Using Simple EmptyDrops (";
        if (config->use_simple_emptydrops) {
            cerr << "explicit --use-simple-empty-drops";
        } else {
            cerr << "fallback: no sparse matrix data";
        }
        cerr << ")" << endl;
        
        SimpleEmptyDropsParams simpleParams;
        simpleParams.nExpectedCells = config->n_expected_cells;
        simpleParams.maxPercentile = config->max_percentile;
        simpleParams.maxMinRatio = config->max_min_ratio;
        simpleParams.umiMin = config->umi_min;
        simpleParams.umiMinFracMedian = config->umi_min_frac_median;
        simpleParams.candMaxN = config->cand_max_n;
        simpleParams.indMin = config->ind_min;
        simpleParams.indMax = config->ind_max;
        simpleParams.useBootstrap = false;
        
        SimpleEmptyDropsResult simpleResult = SimpleEmptyDropsStage::runCRSimpleFilter(
            retainUMI, retainIndices.size(), simpleParams);
        
        result->n_simple_cells = simpleResult.nCellsSimple;
        result->retain_threshold = simpleResult.retainThreshold;
        result->min_umi = simpleResult.minUMI;
        
        cerr << "[scrna_api] Simple filter: " << simpleResult.nCellsSimple << " cells, threshold=" << simpleResult.retainThreshold << endl;
        
        // Map simple filter passing indices back to original indices
        result->n_barcodes = simpleResult.passingIndices.size();
        result->barcodes = (char**)malloc(result->n_barcodes * sizeof(char*));
        if (!result->barcodes) {
            result->error_message = strdup_safe("Memory allocation failed");
            return -1;
        }
        
        for (size_t i = 0; i < result->n_barcodes; i++) {
            uint32_t retainIdx = simpleResult.passingIndices[i];
            uint32_t origIdx = retainIndices[retainIdx];
            result->barcodes[i] = strdup_safe(barcodes[origIdx].c_str());
        }
        
        return 0;
    }
    
    // ========================================================================
    // Step 4 onwards: Full Monte Carlo EmptyDrops (primary method)
    // ========================================================================
    cerr << "[scrna_api] Running full Monte Carlo EmptyDrops (seed=" << config->seed << ")" << endl;
    
    // Compute retain threshold for Monte Carlo mode (used for "simple cells" auto-pass)
    // Sort retainUMI to find percentile-based threshold
    vector<uint32_t> sortedUMI = retainUMI;
    sort(sortedUMI.begin(), sortedUMI.end(), greater<uint32_t>());
    
    uint32_t nExpected = config->n_expected_cells;
    uint32_t baselineIdx = min(nExpected, (uint32_t)sortedUMI.size()) - 1;
    uint32_t retainThreshold = baselineIdx < sortedUMI.size() ? 
        sortedUMI[baselineIdx] / 10 : config->umi_min;  // OrdMag style: baseline / 10
    retainThreshold = max(retainThreshold, config->umi_min);
    
    result->retain_threshold = retainThreshold;
    result->min_umi = config->umi_min;
    
    // ========================================================================
    // Build ambient profile from cells with UMI <= ambientUmiMax WITHIN retain window
    // ========================================================================
    vector<uint32_t> ambCount(input->n_features, 0);
    size_t ambientCellsUsed = 0;
    uint32_t ambientUmiMax = config->ambient_umi_max;
    
    for (size_t ri = 0; ri < retainIndices.size(); ri++) {
        uint32_t origIdx = retainIndices[ri];
        uint32_t umi = nUMIperCB[origIdx];
        if (umi > ambientUmiMax) continue;  // Only cells with UMI <= ambientUmiMax
        
        uint32_t start = input->sparse_cell_index[origIdx];
        uint32_t nGenes = input->n_genes_per_cell[origIdx];
        for (uint32_t g = 0; g < nGenes; g++) {
            uint32_t geneId = input->sparse_gene_ids[start + g];
            uint32_t count = input->sparse_counts[start + g];
            if (geneId < input->n_features) {
                ambCount[geneId] += count;
            }
        }
        ambientCellsUsed++;
    }
    
    cerr << "[scrna_api] Ambient profile: " << ambientCellsUsed << " cells (UMI <= " << ambientUmiMax << " within retain)" << endl;
    
    // Find detected features
    vector<uint32_t> featDetVec;
    for (uint32_t i = 0; i < input->n_features; i++) {
        if (ambCount[i] > 0) featDetVec.push_back(i);
    }
    
    AmbientProfile ambProfile = EmptyDropsMultinomial::computeAmbientProfile(
        ambCount, input->n_features, featDetVec, featDetVec.size());
    
    cerr << "[scrna_api] Ambient profile: " << ambProfile.featuresDetected << " features detected" << endl;
    
    // ========================================================================
    // Build candidate list (UMI > lowerTestingBound within retain window)
    // ========================================================================
    uint32_t lowerTestingBound = config->lower_testing_bound;
    
    // Build sorted index within retain window (by UMI descending)
    vector<pair<uint32_t, uint32_t>> retainUmiIdx;  // (umi, retain_index)
    for (size_t ri = 0; ri < retainIndices.size(); ri++) {
        retainUmiIdx.push_back({retainUMI[ri], ri});
    }
    stable_sort(retainUmiIdx.begin(), retainUmiIdx.end(), [](const pair<uint32_t,uint32_t>& a, const pair<uint32_t,uint32_t>& b) {
        return a.first > b.first;
    });
    
    // Candidates: cells with UMI > lowerTestingBound, capped by candMaxN
    vector<uint32_t> candidateRetainIndices;  // indices into retain window
    vector<uint32_t> candidateCounts;
    uint32_t nSimpleCells = 0;
    
    for (size_t i = 0; i < retainUmiIdx.size(); i++) {
        uint32_t umi = retainUmiIdx[i].first;
        uint32_t retainIdx = retainUmiIdx[i].second;
        
        if (umi <= lowerTestingBound) break;  // No more candidates (sorted descending)
        
        candidateRetainIndices.push_back(retainIdx);
        candidateCounts.push_back(umi);
        
        // Count simple cells (above retain threshold)
        if (umi >= result->retain_threshold) {
            nSimpleCells++;
        }
        
        if (candidateRetainIndices.size() >= config->cand_max_n) break;
    }
    
    cerr << "[scrna_api] Candidates: " << candidateRetainIndices.size() 
         << " (UMI > " << lowerTestingBound << ", cap=" << config->cand_max_n << ")" << endl;
    cerr << "[scrna_api] Simple cells (UMI >= " << result->retain_threshold << "): " << nSimpleCells << endl;
    
    // ========================================================================
    // Build sparse matrix for EmptyDrops
    // ========================================================================
    vector<uint32_t> countCellGeneUMI;
    vector<uint32_t> countCellGeneUMIindex(input->n_cells);
    vector<uint32_t> nGenePerCB(input->n_cells);
    
    for (uint32_t c = 0; c < input->n_cells; c++) {
        countCellGeneUMIindex[c] = countCellGeneUMI.size();
        uint32_t start = input->sparse_cell_index[c];
        uint32_t nGenes = input->n_genes_per_cell[c];
        nGenePerCB[c] = nGenes;
        for (uint32_t g = 0; g < nGenes; g++) {
            countCellGeneUMI.push_back(input->sparse_gene_ids[start + g]);  // gene ID
            countCellGeneUMI.push_back(input->sparse_counts[start + g]);    // count
        }
    }
    
    // Map candidate retain indices back to original indices for ED
    vector<uint32_t> candidateOrigIndices;
    for (uint32_t retainIdx : candidateRetainIndices) {
        candidateOrigIndices.push_back(retainIndices[retainIdx]);
    }
    
    // Set up EmptyDrops parameters
    EmptyDropsParams edParams;
    edParams.indMin = config->ind_min;
    edParams.indMax = config->ind_max;
    edParams.umiMin = config->umi_min;
    edParams.umiMinFracMedian = config->umi_min_frac_median;
    edParams.candMaxN = config->cand_max_n;
    edParams.FDR = config->fdr;
    edParams.rawPvalueThreshold = config->raw_pvalue_threshold;
    edParams.simN = config->sim_n;
    edParams.seed = config->seed;
    edParams.lowerTestingBound = config->lower_testing_bound;
    edParams.ambientUmiMax = config->ambient_umi_max;
    edParams.mcThreads = config->mc_threads;
    
    // ========================================================================
    // Run EmptyDrops Monte Carlo
    // ========================================================================
    cerr << "[scrna_api] Running EmptyDrops MC with " << edParams.simN << " simulations..." << endl;
    
    vector<EmptyDropsResult> edResults = EmptyDropsMultinomial::computePValues(
        ambProfile,
        candidateOrigIndices,
        candidateCounts,
        countCellGeneUMI,
        countCellGeneUMIindex,
        nGenePerCB,
        2,  // stride (gene, count pairs)
        1,  // count offset
        edParams,
        nSimpleCells,  // Simple cells auto-pass
        vector<string>(),  // featureNames (unused)
        input->n_cells,
        "",  // debugOutputDir
        "",  // tagName
        false  // enableInvariantChecks
    );
    
    // ========================================================================
    // Collect passing barcodes (use passesRawP unless use_fdr_gate is set)
    // ========================================================================
    vector<string> passingBarcodes;
    uint32_t nEdPassers = 0;
    
    for (size_t i = 0; i < edResults.size(); i++) {
        bool passes = config->use_fdr_gate ? edResults[i].passesFDR : edResults[i].passesRawP;
        if (passes) {
            uint32_t cellIdx = edResults[i].cellIndex;
            if (cellIdx < barcodes.size()) {
                passingBarcodes.push_back(barcodes[cellIdx]);
            }
            if (i >= nSimpleCells) {
                nEdPassers++;
            }
        }
    }
    
    cerr << "[scrna_api] EmptyDrops complete: " << passingBarcodes.size() << " cells pass" << endl;
    cerr << "[scrna_api]   Simple passers: " << nSimpleCells << endl;
    cerr << "[scrna_api]   ED rescues: " << nEdPassers << endl;
    
    // ========================================================================
    // Store results
    // ========================================================================
    result->n_barcodes = passingBarcodes.size();
    result->barcodes = (char**)malloc(result->n_barcodes * sizeof(char*));
    if (!result->barcodes) {
        result->error_message = strdup_safe("Memory allocation failed");
        return -1;
    }
    
    for (size_t i = 0; i < result->n_barcodes; i++) {
        result->barcodes[i] = strdup_safe(passingBarcodes[i].c_str());
    }
    
    // Store candidate results with barcodes
    result->n_candidates = edResults.size();
    result->candidates = (scrna_ed_candidate*)malloc(result->n_candidates * sizeof(scrna_ed_candidate));
    if (result->candidates) {
        for (size_t i = 0; i < result->n_candidates; i++) {
            uint32_t cellIdx = edResults[i].cellIndex;
            result->candidates[i].cell_index = cellIdx;
            result->candidates[i].barcode = (cellIdx < barcodes.size()) 
                ? strdup_safe(barcodes[cellIdx].c_str()) : nullptr;
            result->candidates[i].umi_count = (cellIdx < nUMIperCB.size()) ? nUMIperCB[cellIdx] : 0;
            result->candidates[i].p_value = edResults[i].pValue;
            result->candidates[i].p_adjusted = edResults[i].pAdjusted;
            result->candidates[i].passes_raw_p = edResults[i].passesRawP ? 1 : 0;
            result->candidates[i].passes_fdr = edResults[i].passesFDR ? 1 : 0;
            result->candidates[i].obs_log_prob = edResults[i].obsLogProb;
            result->candidates[i].is_simple_cell = (i < nSimpleCells) ? 1 : 0;
        }
    }
    
    // Store statistics
    result->n_simple_cells = nSimpleCells;
    result->n_tail_cells = edResults.size() > nSimpleCells ? edResults.size() - nSimpleCells : 0;
    result->n_ed_passers = nEdPassers;
    
    return 0;
}

// ============================================================================
// Output functions
// ============================================================================

extern "C" int scrna_write_filtered_barcodes(
    const scrna_ed_result *result,
    const char *filepath
) {
    if (!result || !filepath) return -1;
    
    ofstream out(filepath);
    if (!out.is_open()) return -1;
    
    for (size_t i = 0; i < result->n_barcodes; i++) {
        if (result->barcodes[i]) {
            out << result->barcodes[i] << "\n";
        }
    }
    
    out.close();
    return 0;
}

extern "C" int scrna_emptydrops_write_outputs(
    const scrna_ed_result *result,
    const char *output_dir
) {
    if (!result || !output_dir) return -1;
    
    // Create output directories
    mkdir_p(output_dir);
    string edDir = string(output_dir) + "/EmptyDrops";
    mkdir_p(edDir.c_str());
    
    // Write filtered_barcodes.txt
    string barcodesPath = string(output_dir) + "/filtered_barcodes.txt";
    if (scrna_write_filtered_barcodes(result, barcodesPath.c_str()) != 0) {
        return -1;
    }
    
    // Write emptydrops_results.tsv (Flex-compatible schema with barcodes)
    string resultsPath = edDir + "/emptydrops_results.tsv";
    ofstream tsv(resultsPath);
    if (!tsv.is_open()) return -1;
    
    // Header (matches Flex output)
    tsv << "barcode\tcell_index\tumi_count\tp_value\tp_adjusted\tpasses_raw_p\tpasses_fdr\tobs_log_prob\tis_simple_cell\n";
    
    // Data
    if (result->candidates) {
        for (size_t i = 0; i < result->n_candidates; i++) {
            tsv << (result->candidates[i].barcode ? result->candidates[i].barcode : "") << "\t"
                << result->candidates[i].cell_index << "\t"
                << result->candidates[i].umi_count << "\t"
                << result->candidates[i].p_value << "\t"
                << result->candidates[i].p_adjusted << "\t"
                << result->candidates[i].passes_raw_p << "\t"
                << result->candidates[i].passes_fdr << "\t"
                << result->candidates[i].obs_log_prob << "\t"
                << result->candidates[i].is_simple_cell << "\n";
        }
    }
    
    tsv.close();
    
    // Write summary JSON
    string summaryPath = edDir + "/summary.json";
    ofstream json(summaryPath);
    if (json.is_open()) {
        json << "{\n";
        json << "  \"n_cells_passing\": " << result->n_barcodes << ",\n";
        json << "  \"n_candidates_tested\": " << result->n_candidates << ",\n";
        json << "  \"n_simple_cells\": " << result->n_simple_cells << ",\n";
        json << "  \"n_tail_cells\": " << result->n_tail_cells << ",\n";
        json << "  \"n_ed_passers\": " << result->n_ed_passers << ",\n";
        json << "  \"retain_threshold\": " << result->retain_threshold << ",\n";
        json << "  \"min_umi\": " << result->min_umi << "\n";
        json << "}\n";
        json.close();
    }
    
    return 0;
}
