#include "FlexFilter.h"
#include "FlexFilterIO.h"
#include "SimpleEDCaller.h"
#include "MitochondrialRankMask.h"
#include "FlexTagGroup.h"
#include "ObservedTagOccupancy.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>

int FlexFilter::runTagAware(const SampleMatrixData& matrix,
    const std::vector<std::string>& labels, const std::vector<std::string>& tags,
    Outputs* outputs, const Config& config)
{
    if (!outputs) return 1;
    outputs->tagResults.clear();
    try {
        if (labels.size() != tags.size() || labels.empty())
            throw std::runtime_error("Tag-aware caller needs aligned sample labels and TAG8 sequences");
        if (matrix.features.size() != matrix.nGenes || matrix.barcodes.size() != matrix.nCells ||
            matrix.nUMIperCB.size() != matrix.nCells || matrix.nGenePerCB.size() != matrix.nCells ||
            matrix.countCellGeneUMIindex.size() < matrix.nCells || matrix.countMatStride < 2)
            throw std::runtime_error("Misaligned matrix metadata in tag-aware caller");

        std::vector<std::string> groups;
        std::map<std::string, std::vector<std::string>> groupTags;
        std::set<std::string> seenTags;
        for (size_t i = 0; i < labels.size(); ++i) {
            if (labels[i].empty() || labels[i] == "." || labels[i] == ".." ||
                labels[i].find_first_of("/\\\t\r\n") != std::string::npos)
                throw std::runtime_error("Invalid sample label: " + labels[i]);
            if (tags[i].size() != 8 || tags[i].find_first_not_of("ACGT") != std::string::npos ||
                !seenTags.insert(tags[i]).second)
                throw std::runtime_error("Invalid or duplicate TAG8 in sample whitelist: " + tags[i]);
            if (!groupTags.count(labels[i])) groups.push_back(labels[i]);
            groupTags[labels[i]].push_back(tags[i]);
        }
        const std::vector<uint8_t> mt = loadMitochondrialRankMask(matrix.features, config.mitochondrialGenesPath);
        if (!config.debugOutputDir.empty()) {
            if (!FlexFilterIO::defaultCreateDirectory(config.debugOutputDir, 0775))
                throw std::runtime_error("Cannot create Flex caller diagnostics directory");
            std::ofstream annotations(config.debugOutputDir + "/feature_rank_mask.tsv");
            annotations << "feature_index\tfeature_id\tis_mitochondrial\n";
            for (uint32_t g = 0; g < matrix.nGenes; ++g)
                annotations << g << '\t' << matrix.features[g] << '\t' << unsigned(mt.empty() ? 0 : mt[g]) << '\n';
            annotations.close();
            if (!annotations) throw std::runtime_error("Cannot write feature rank mask");
        }

        // Groups run serially. Bootstrap receives the declared worker budget;
        // EmptyDrops uses the same eight-worker default as the external caller.
        for (const std::string& label : groups) {
            const auto& sampleTags = groupTags.at(label);
            std::cerr << "[Flex tag-aware] sample=" << label << " tags=" << sampleTags.size() << "\n";
            if (sampleTags.size() > std::numeric_limits<uint32_t>::max() / 90000u)
                throw std::runtime_error("Too many tags for caller rank limits");
            const uint32_t nTags = static_cast<uint32_t>(sampleTags.size());
            std::vector<std::string> barcodes;
            std::vector<uint32_t> umi, genes, counts, starts, nGenes;
            for (uint32_t cell = 0; cell < matrix.nCells; ++cell) {
                if (!matrix.nUMIperCB[cell] || !barcodeHasAnyFlexTag(matrix.barcodes[cell], sampleTags)) continue;
                if (counts.size() > std::numeric_limits<uint32_t>::max())
                    throw std::runtime_error("Sparse group matrix exceeds 32-bit index range");
                barcodes.push_back(matrix.barcodes[cell]);
                umi.push_back(matrix.nUMIperCB[cell]);
                starts.push_back(static_cast<uint32_t>(counts.size()));
                nGenes.push_back(matrix.nGenePerCB[cell]);
                uint64_t summed = 0;
                for (uint32_t j = 0; j < matrix.nGenePerCB[cell]; ++j) {
                    const size_t pos = static_cast<size_t>(matrix.countCellGeneUMIindex[cell])
                        + static_cast<size_t>(j) * matrix.countMatStride;
                    if (pos + 1 >= matrix.countCellGeneUMI.size() || matrix.countCellGeneUMI[pos] >= matrix.nGenes)
                        throw std::runtime_error("Invalid sparse coordinate in group matrix");
                    genes.push_back(matrix.countCellGeneUMI[pos]);
                    counts.push_back(matrix.countCellGeneUMI[pos + 1]);
                    summed += matrix.countCellGeneUMI[pos + 1];
                }
                if (config.enableInvariantChecks && summed != umi.back())
                    throw std::runtime_error("Sparse counts do not match total UMIs for " + barcodes.back());
            }
            Outputs::TagResults output;
            output.sampleLabel = label;
            output.tag = sampleTags.front();
            output.expectedCells = 0;
            output.tagBarcodes = barcodes;
            output.retainBarcodes = barcodes; // ED cell indices refer to this group matrix.
            if (barcodes.empty()) {
                std::cerr << "[Flex tag-aware] " << label << ": no nonzero barcodes\n";
                outputs->tagResults.push_back(std::move(output));
                continue;
            }
            std::unique_ptr<scrna_ed_config, decltype(&scrna_ed_config_destroy)> ed(
                scrna_ed_config_create(), &scrna_ed_config_destroy);
            if (!ed) throw std::bad_alloc();
            ed->use_bootstrap = 1;
            ed->n_expected_cells = 0;
            ed->ind_min = 45000u * nTags;
            ed->ind_max = 90000u * nTags;
            ed->ed_retain_count = ed->ind_max;
            ed->cand_max_n = 100000;
            ed->sim_n = config.emptydropsParams.simN ? config.emptydropsParams.simN : 10000;
            ed->fdr = config.emptydropsParams.FDR > 0 ? config.emptydropsParams.FDR : .01;
            ed->mc_threads = config.emptydropsParams.mcThreads ? config.emptydropsParams.mcThreads : 8;
            ed->use_fdr_gate = 1;
            ed->apply_bh_correction = 1;
            SimpleEDOptions options;
            options.ordmagRetainCount = 90000u * nTags;
            options.maxExpectedCells = 22500u * nTags;
            options.ambientUmiTarget = 0; // Fixed tag-scaled ambient ranks.
            options.bootstrapThreads = config.simpleEmptyDropsParams.maxThreads;
            options.invariantChecks = config.enableInvariantChecks;
            if (!config.debugOutputDir.empty()) options.diagnosticsDir = config.debugOutputDir + "/" + label;
            scrna_ed_result result = {};
            SimpleEDRunInfo info;
            try {
                const int rc = runSimpleEDWithAmbient(barcodes, umi, genes, counts, starts, nGenes,
                    matrix.nGenes, ed.get(), options, mt, &result, &info);
                if (rc) throw std::runtime_error(result.error_message ? result.error_message : "Shared caller failed");
                output.expectedCells = info.bootstrap.recoveredCells;
                output.nRetainWindow = info.retainCount;
                output.nSimpleCells = result.n_simple_cells;
                output.nSimplePassers = result.n_simple_cells;
                output.nTailTested = result.n_tail_cells;
                output.nTailPassers = result.n_ed_passers;
                for (size_t i = 0; i < result.n_barcodes; ++i) output.passingBarcodes.emplace_back(result.barcodes[i]);
                output.filteredBarcodes = output.passingBarcodes;
                output.edPasserBarcodes = output.passingBarcodes;
                for (size_t i = 0; i < result.n_candidates; ++i) {
                    const auto& c = result.candidates[i];
                    EmptyDropsResult e = {};
                    e.cellIndex = c.cell_index; e.pValue = c.p_value; e.pAdjusted = c.p_adjusted;
                    e.passesRawP = c.passes_raw_p; e.passesFDR = c.passes_fdr; e.obsLogProb = c.obs_log_prob;
                    output.emptydropsResults.push_back(e);
                }
                if (!options.diagnosticsDir.empty()) {
                    if (scrna_emptydrops_write_outputs(&result, options.diagnosticsDir.c_str()) != 0)
                        throw std::runtime_error("Cannot write EmptyDrops diagnostic results");
                    std::ofstream membership(options.diagnosticsDir + "/sample_tags.tsv");
                    for (const auto& tag : sampleTags) membership << label << '\t' << tag << '\n';
                    membership.close();
                    if (!membership) throw std::runtime_error("Cannot write sample tag membership");
                }
                scrna_ed_result_free(&result);
            } catch (...) { scrna_ed_result_free(&result); throw; }
            outputs->tagResults.push_back(std::move(output));
        }
        if (!config.disableOccupancyFilter) {
            std::vector<std::string> calls;
            for (const auto& group : outputs->tagResults)
                calls.insert(calls.end(), group.passingBarcodes.begin(), group.passingBarcodes.end());
            const auto occupancy = fitObservedTagOccupancy(calls, config.occupancyPercentile);
            std::cerr << "[Flex tag-aware occupancy] occupiedGems=" << occupancy.occupiedGems
                << " mean=" << occupancy.occupiedMean << " lambda=" << occupancy.lambda
                << " cutoff=" << occupancy.cutoff << " removedGems=" << occupancy.rejectedGems.size() << "\n";
            for (auto& group : outputs->tagResults) {
                const size_t before = group.passingBarcodes.size();
                auto& kept = group.passingBarcodes;
                kept.erase(std::remove_if(kept.begin(), kept.end(), [&](const std::string& bc) {
                    return occupancy.rejectedGems.count(bc.substr(0, 16)) != 0;
                }), kept.end());
                group.filteredBarcodes = kept;
                group.occupancyRemoved = before - kept.size();
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "[Flex tag-aware] ERROR: " << e.what() << "\n";
        outputs->tagResults.clear();
        return 1;
    }
    return 0;
}
