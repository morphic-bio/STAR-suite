#include "SoloFeature.h"
#include "serviceFuns.cpp"
#include "scrna_api.h"

#include <cstring>
#include <unordered_map>

void SoloFeature::emptyDrops_libscrna()
{
    const bool forceUnionMode = (pSolo.emptyDropsMode == ParametersSolo::EmptyDropsUnion);
    const bool bootstrapEnabled = !pSolo.emptyDropsLegacyKnee;
    if (nCB <= pSolo.cellFilter.eDcr.indMin && !forceUnionMode && !bootstrapEnabled) {
        P.inOut->logMain << "emptyDrops_CR (libscrna) filtering: total number of cells: nCB=" << nCB
                         << " is smaller than emptyCellMinIndex=" << pSolo.cellFilter.eDcr.indMin
                         << ", which is the starting index for the *true empty* cells."
                         << " The additional non-empty cells will not be detected."
                         << " (Use --soloEmptyDropsLegacyKnee no to enable bootstrap OrdMag on shallow data.)\n";
        return;
    }
    if (nCB <= pSolo.cellFilter.eDcr.indMin && forceUnionMode) {
        P.inOut->logMain << "emptyDrops_CR (libscrna): forcing union mode despite nCB=" << nCB
                         << " <= emptyCellMinIndex=" << pSolo.cellFilter.eDcr.indMin
                         << "; attempting EmptyDrops rescue anyway.\n";
    }
    if (nCB <= pSolo.cellFilter.eDcr.indMin && bootstrapEnabled && !forceUnionMode) {
        P.inOut->logMain << "emptyDrops_CR (libscrna): nCB=" << nCB
                         << " < emptyCellMinIndex=" << pSolo.cellFilter.eDcr.indMin
                         << "; bootstrap mode active — MC tail rescue will proceed with ambient fallback.\n";
    }

    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... starting emptyDrops_CR filtering (libscrna)" << endl;

    // Build barcode list for detected CBs
    vector<string> barcodes;
    barcodes.reserve(nCB);
    for (uint32 icb = 0; icb < nCB; icb++) {
        barcodes.push_back(pSolo.cbWLstr[indCB[icb]]);
    }
    vector<char*> barcode_ptrs;
    barcode_ptrs.reserve(barcodes.size());
    for (auto &bc : barcodes) {
        barcode_ptrs.push_back(const_cast<char*>(bc.c_str()));
    }

    // Build sparse matrix (main count only)
    vector<uint32_t> sparse_gene_ids;
    vector<uint32_t> sparse_counts;
    vector<uint32_t> sparse_cell_index(nCB + 1, 0);
    vector<uint32_t> n_genes_per_cell(nCB, 0);
    size_t nnz = 0;

    for (uint32 icb = 0; icb < nCB; icb++) {
        sparse_cell_index[icb] = static_cast<uint32_t>(nnz);
        uint32 nGenes = nGenePerCB[icb];
        for (uint32 ig = 0; ig < nGenes; ig++) {
            uint32 irec = countCellGeneUMIindex[icb] + ig * countMatStride;
            if (irec + pSolo.umiDedup.countInd.main >= countCellGeneUMI.size()) {
                continue;
            }
            uint32 geneId = countCellGeneUMI[irec];
            uint32 count = countCellGeneUMI[irec + pSolo.umiDedup.countInd.main];
            if (count == 0) {
                continue;
            }
            sparse_gene_ids.push_back(geneId);
            sparse_counts.push_back(count);
            n_genes_per_cell[icb]++;
            nnz++;
        }
    }
    sparse_cell_index[nCB] = static_cast<uint32_t>(nnz);

    // Prepare input for libscrna
    scrna_matrix_input input;
    memset(&input, 0, sizeof(input));
    input.umi_counts = nUMIperCB.data();
    input.barcodes = barcode_ptrs.data();
    input.n_cells = nCB;
    input.n_features = featuresNumber;
    input.features = nullptr;
    input.sparse_gene_ids = sparse_gene_ids.data();
    input.sparse_counts = sparse_counts.data();
    input.sparse_cell_index = sparse_cell_index.data();
    input.n_genes_per_cell = n_genes_per_cell.data();
    input.sparse_nnz = nnz;

    scrna_ed_config *config = scrna_ed_config_create();
    if (config == nullptr) {
        P.inOut->logMain << "emptyDrops_CR (libscrna) failed: could not allocate config\n";
        return;
    }

    // Map STAR Solo EmptyDrops_CR parameters (scRNA-seq defaults, not Flex)
    config->max_percentile = pSolo.cellFilter.knee.maxPercentile;
    config->max_min_ratio = pSolo.cellFilter.knee.maxMinRatio;
    config->ind_min = pSolo.cellFilter.eDcr.indMin;
    config->ind_max = pSolo.cellFilter.eDcr.indMax;
    config->umi_min = 100; // scRNA-seq: low floor to allow EmptyDrops rescue below the OrdMag knee
    config->umi_min_frac_median = pSolo.cellFilter.eDcr.umiMinFracMedian;
    config->cand_max_n = pSolo.cellFilter.eDcr.candMaxN;
    config->fdr = pSolo.cellFilter.eDcr.FDR;
    config->sim_n = 100000; // match CR9: 100K MC simulations for better p-value resolution
    config->use_fdr_gate = 1;
    config->apply_bh_correction = 1; // scRNA-seq: proper BH-corrected FDR (matches CR9)
    config->mc_threads = 0;
    config->disable_occupancy_filter = 1;

    if (bootstrapEnabled) {
        config->n_expected_cells = 0;  // triggers bootstrap estimation (CR9 style)
        config->use_bootstrap = 1;
    } else {
        config->n_expected_cells = pSolo.cellFilter.knee.nExpectedCells;
        config->use_bootstrap = 0;
    }

    P.inOut->logMain << "emptyDrops_CR (libscrna) params: indMin=" << config->ind_min
                     << " indMax=" << config->ind_max
                     << " umiMin=" << config->umi_min
                     << " umiMinFracMedian=" << config->umi_min_frac_median
                     << " candMaxN=" << config->cand_max_n
                     << " FDR=" << config->fdr
                     << " simN=" << config->sim_n
                     << " bootstrap=" << (config->use_bootstrap ? "yes" : "no")
                     << " legacyKnee=" << (pSolo.emptyDropsLegacyKnee ? "yes" : "no")
                     << " mode=" << (forceUnionMode ? "union" : "auto") << "\n";

    scrna_ed_result result;
    memset(&result, 0, sizeof(result));
    int rc = scrna_emptydrops_run(&input, config, &result);
    if (rc != 0) {
        P.inOut->logMain << "emptyDrops_CR (libscrna) failed: " 
                         << (result.error_message ? result.error_message : "unknown error") << "\n";
        scrna_ed_result_free(&result);
        scrna_ed_config_destroy(config);
        return;
    }

    // Overwrite filteredCells based on libscrna results
    filteredCells.filtVecBool.assign(nCB, false);
    filteredCells.nCellsSimple = result.n_simple_cells;

    unordered_map<string, uint32_t> bcIndex;
    bcIndex.reserve(nCB * 1.3);
    for (uint32_t i = 0; i < nCB; i++) {
        bcIndex.emplace(barcodes[i], i);
    }

    uint32_t passCount = 0;
    for (size_t i = 0; i < result.n_barcodes; i++) {
        const char *bc = result.barcodes[i];
        if (!bc) {
            continue;
        }
        auto it = bcIndex.find(bc);
        if (it != bcIndex.end()) {
            filteredCells.filtVecBool[it->second] = true;
            passCount++;
        }
    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ... finished emptyDrops_CR filtering (libscrna): total pass=" << passCount
                     << " (simple=" << result.n_simple_cells
                     << ", ED rescues=" << result.n_ed_passers << ")\n";

    // Write EmptyDrops outputs for audit
    string edOutDir = outputPrefixFiltered + "EmptyDrops";
    if (scrna_emptydrops_write_outputs(&result, edOutDir.c_str()) != 0) {
        P.inOut->logMain << "emptyDrops_CR (libscrna) WARNING: failed to write outputs to " << edOutDir << "\n";
    }

    scrna_ed_result_free(&result);
    scrna_ed_config_destroy(config);
}
