#include "SoloFeature.h"
#include "serviceFuns.cpp"
#include "libscrna/OrdMagStage.h"
#include "scrna_api.h"
#include "ScrnaTrace.h"
#include "MitochondrialRankMask.h"
#include "OrdMagRank.h"
#include "ErrorWarning.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>

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

    SparseCountView matrix;
    matrix.genes = countCellGeneUMI.data();
    matrix.geneWords = countCellGeneUMI.size();
    const size_t shift = pSolo.umiDedup.countInd.main;
    if (shift < countCellGeneUMI.size()) {
        matrix.counts = countCellGeneUMI.data() + shift;
        matrix.countWords = countCellGeneUMI.size() - shift;
    }
    matrix.stride = countMatStride;
    matrix.offsets = countCellGeneUMIindex.data();
    matrix.entries = nGenePerCB.data();
    matrix.cells = nCB;
    P.inOut->logMain << "emptyDrops matrix: borrowed " << countCellGeneUMI.size() * sizeof(uint32_t)
                    << " bytes; adapter and interleave copies=0\n";

    // Prepare input for libscrna
    scrna_matrix_input input;
    memset(&input, 0, sizeof(input));
    input.umi_counts = nUMIperCB.data();
    input.barcodes = barcode_ptrs.data();
    input.n_cells = nCB;
    input.n_features = featuresNumber;
    input.features = nullptr;
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
    // EmptyDrops candidates start at the declared umiMin (500 by default).
    // OrdMag primaries get no extra floor: low-depth libraries call cells
    // well below 500, so the declared value must not trim them.
    config->umi_min = pSolo.cellFilter.eDcr.umiMin;
    config->primary_umi_min = 1;
    config->umi_min_frac_median = pSolo.cellFilter.eDcr.umiMinFracMedian;
    config->cand_max_n = pSolo.cellFilter.eDcr.candMaxN;
    config->fdr = pSolo.cellFilter.eDcr.FDR;
    config->sim_n = 100000; // match CR9: 100K MC simulations for better p-value resolution
    config->use_fdr_gate = 1;
    config->apply_bh_correction = 1; // scRNA-seq: proper BH-corrected FDR (matches CR9)
    config->mc_threads = static_cast<uint32_t>(std::max(1, P.runThreadN));
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
                     << " primaryUmiMin=" << config->primary_umi_min
                     << " umiMinFracMedian=" << config->umi_min_frac_median
                     << " candMaxN=" << config->cand_max_n
                     << " FDR=" << config->fdr
                     << " simN=" << config->sim_n
                     << " mcThreads=" << config->mc_threads
                     << " bootstrap=" << (config->use_bootstrap ? "yes" : "no")
                     << " legacyKnee=" << (pSolo.emptyDropsLegacyKnee ? "yes" : "no")
                     << " mode=" << (forceUnionMode ? "union" : "auto") << "\n";

    scrna_ed_result result;
    memset(&result, 0, sizeof(result));
    vector<uint8_t> mitochondrialMask;
    try {
        if (!pSolo.cellFilterMitochondrialGenes.empty()) {
            if (Trans.geID.size() != featuresNumber)
                throw std::runtime_error("MT rank annotation does not match this feature matrix");
            mitochondrialMask = loadMitochondrialRankMask(Trans.geID, pSolo.cellFilterMitochondrialGenes);
        }
    } catch (const std::exception& e) {
        scrna_ed_config_destroy(config);
        exitWithError(string("ERROR: ") + e.what(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    const uint32_t bootstrapThreads = pSolo.cellFilterBootstrapThreads > 0
        ? pSolo.cellFilterBootstrapThreads : static_cast<uint32_t>(std::max(1, P.runThreadN));
    ScrnaTrace trace;
    int rc = scrnaEmptyDropsTrace(&input, config,
        mitochondrialMask.empty() ? nullptr : mitochondrialMask.data(), bootstrapThreads, &result, &trace, &matrix);
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

    // Emit branch-only audit files so libscrna can be compared against legacy on the same ranked cells.
    {
        const auto& retainIndices = trace.retainIndices;
        const auto& simpleResult = trace.ordmag;
        vector<uint32_t> origToRetainRank(nCB, UINT32_MAX);
        for (size_t rank = 0; rank < retainIndices.size(); ++rank)
            origToRetainRank[retainIndices[rank]] = rank;

        unordered_set<uint32_t> ambientOrigIndices;
        ambientOrigIndices.reserve(simpleResult.ambientIndices.size() * 2);
        for (uint32_t retainIdx : simpleResult.ambientIndices) {
            if (retainIdx < retainIndices.size()) {
                ambientOrigIndices.insert(retainIndices[retainIdx]);
            }
        }

        mkdir(edOutDir.c_str(), P.runDirPerm);

        ofstream summaryOut(edOutDir + "/backend_debug_summary.json");
        if (summaryOut.is_open()) {
            summaryOut << "{\n";
            summaryOut << "  \"backend\": \"libscrna\",\n";
            summaryOut << "  \"n_total_cells\": " << nCB << ",\n";
            summaryOut << "  \"n_retain_cells\": " << retainIndices.size() << ",\n";
            summaryOut << "  \"n_simple_cells\": " << result.n_simple_cells << ",\n";
            summaryOut << "  \"n_candidates_total\": " << result.n_candidates << ",\n";
            summaryOut << "  \"n_tail_candidates\": " << result.n_tail_cells << ",\n";
            summaryOut << "  \"n_ambient_cells\": " << simpleResult.ambientIndices.size() << ",\n";
            summaryOut << "  \"ambient_start_rank\": " << simpleResult.ambientRange.first << ",\n";
            summaryOut << "  \"ambient_end_rank_exclusive\": " << simpleResult.ambientRange.second << ",\n";
            summaryOut << "  \"retain_threshold\": " << result.retain_threshold << ",\n";
            summaryOut << "  \"min_umi\": " << result.min_umi << ",\n";
            summaryOut << "  \"n_ed_passers\": " << result.n_ed_passers << ",\n";
            summaryOut << "  \"use_fdr_gate\": " << (config->use_fdr_gate ? 1 : 0) << ",\n";
            summaryOut << "  \"apply_bh_correction\": " << (config->apply_bh_correction ? 1 : 0) << ",\n";
            summaryOut << "  \"use_bootstrap\": " << (config->use_bootstrap ? 1 : 0) << "\n";
            summaryOut << "}\n";
        }

        ofstream candidateOut(edOutDir + "/backend_debug_candidates.tsv");
        if (candidateOut.is_open()) {
            candidateOut << "backend\trank_desc\tbarcode\tcell_index\tumi_count\tis_simple_cell\tis_tail_candidate\tis_ambient_cell\tp_value\tp_adjusted\tpasses_raw_p\tpasses_fdr\tobs_log_prob\n";
            for (size_t i = 0; i < result.n_candidates; i++) {
                const scrna_ed_candidate &cand = result.candidates[i];
                uint32_t retainRank = cand.cell_index < origToRetainRank.size() ? origToRetainRank[cand.cell_index] : (uint32_t)-1;
                bool isAmbient = ambientOrigIndices.count(cand.cell_index) > 0;
                candidateOut << "libscrna\t"
                             << (retainRank == (uint32_t)-1 ? string("") : to_string(retainRank)) << '\t'
                             << (cand.barcode ? cand.barcode : "") << '\t'
                             << cand.cell_index << '\t'
                             << cand.umi_count << '\t'
                             << cand.is_simple_cell << '\t'
                             << (cand.is_simple_cell ? 0 : 1) << '\t'
                             << (isAmbient ? 1 : 0) << '\t'
                             << cand.p_value << '\t'
                             << cand.p_adjusted << '\t'
                             << cand.passes_raw_p << '\t'
                             << cand.passes_fdr << '\t'
                             << cand.obs_log_prob << '\n';
            }
        }

        ofstream ambientOut(edOutDir + "/backend_debug_ambient.tsv");
        if (ambientOut.is_open()) {
            ambientOut << "backend\trank_desc\tbarcode\tcell_index\tumi_count\n";
            for (uint32_t retainIdx : simpleResult.ambientIndices) {
                if (retainIdx >= retainIndices.size()) {
                    continue;
                }
                uint32_t origIdx = retainIndices[retainIdx];
                ambientOut << "libscrna\t"
                           << retainIdx << '\t'
                           << barcodes[origIdx] << '\t'
                           << origIdx << '\t'
                           << nUMIperCB[origIdx] << '\n';
            }
        }
    }

    scrna_ed_result_free(&result);
    scrna_ed_config_destroy(config);
}
