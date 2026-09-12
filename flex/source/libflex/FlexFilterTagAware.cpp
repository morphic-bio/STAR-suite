#include "FlexFilter.h"
#include "FlexFilterIO.h"
#include "SimpleEDCaller.h"
#include "MitochondrialRankMask.h"
#include "FlexTagGroup.h"
#include "ObservedTagOccupancy.h"
#include "ParallelTasks.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <omp.h>

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

        const uint32_t totalThreads = config.totalThreads ? config.totalThreads
            : static_cast<uint32_t>(std::max(1, omp_get_max_threads()));
        const size_t groupWorkers = std::min<size_t>(groups.size(), totalThreads);
        const uint32_t threadsPerGroup = totalThreads / groupWorkers;
        const uint32_t extraThreads = totalThreads % groupWorkers;
        // Keep paired tags in one sample model. Workers share a bounded CPU
        // budget and write disjoint, ordered result slots.
        // Bootstrap RNG partitions remain unchanged when execution is capped.
        outputs->tagResults.resize(groups.size());
        struct GroupTiming { double prepare = 0, call = 0; uint32_t workers = 0, mc = 0; };
        std::vector<GroupTiming> timings(groups.size());
        const auto started = std::chrono::steady_clock::now();
        std::cerr << "[Flex tag-aware parallel] groups=" << groups.size()
            << " group_workers=" << groupWorkers << " total_threads=" << totalThreads
            << " scheduling=" << (config.useThreadPermits ? "permits" : "fixed")
            << " bootstrap_streams=" << config.simpleEmptyDropsParams.maxThreads << "\n";
        auto callGroup = [&](size_t group, size_t worker) {
            const auto groupStarted = std::chrono::steady_clock::now();
            const std::string& label = groups[group];
            const uint32_t workerBudget = config.useThreadPermits ? totalThreads
                : threadsPerGroup + (worker < extraThreads ? 1 : 0);
            timings[group].workers = workerBudget;
            const auto& sampleTags = groupTags.at(label);
            if (sampleTags.size() > std::numeric_limits<uint32_t>::max() / 90000u)
                throw std::runtime_error("Too many tags for caller rank limits");
            const uint32_t nTags = static_cast<uint32_t>(sampleTags.size());
            auto& output = outputs->tagResults[group];
            auto& barcodes = output.tagBarcodes;
            std::vector<uint32_t> umi, starts, nGenes;
            for (uint32_t cell = 0; cell < matrix.nCells; ++cell) {
                if (!matrix.nUMIperCB[cell] || !barcodeHasAnyFlexTag(matrix.barcodes[cell], sampleTags)) continue;
                barcodes.push_back(matrix.barcodes[cell]);
                umi.push_back(matrix.nUMIperCB[cell]);
                starts.push_back(matrix.countCellGeneUMIindex[cell]);
                nGenes.push_back(matrix.nGenePerCB[cell]);
                uint64_t summed = 0;
                for (uint32_t j = 0; j < matrix.nGenePerCB[cell]; ++j) {
                    const size_t pos = static_cast<size_t>(matrix.countCellGeneUMIindex[cell])
                        + static_cast<size_t>(j) * matrix.countMatStride;
                    if (pos + 1 >= matrix.countCellGeneUMI.size() || matrix.countCellGeneUMI[pos] >= matrix.nGenes)
                        throw std::runtime_error("Invalid sparse coordinate in group matrix");
                    summed += matrix.countCellGeneUMI[pos + 1];
                }
                if (config.enableInvariantChecks && summed != umi.back())
                    throw std::runtime_error("Sparse counts do not match total UMIs for " + barcodes.back());
            }
            output.sampleLabel = label;
            output.tag = sampleTags.front();
            output.expectedCells = 0;
            output.retainBarcodes = barcodes; // ED cell indices refer to this group matrix.
            if (barcodes.empty()) {
                return;
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
            ed->mc_threads = std::min(workerBudget, config.emptydropsParams.mcThreads
                ? config.emptydropsParams.mcThreads : workerBudget);
            timings[group].mc = ed->mc_threads;
            ed->use_fdr_gate = 1;
            ed->apply_bh_correction = 1;
            SimpleEDOptions options;
            options.ordmagRetainCount = 90000u * nTags;
            options.maxExpectedCells = 22500u * nTags;
            options.ambientUmiTarget = 0; // Fixed tag-scaled ambient ranks.
            options.bootstrapThreads = config.simpleEmptyDropsParams.maxThreads;
            options.bootstrapWorkers = workerBudget;
            options.invariantChecks = config.enableInvariantChecks;
            if (!config.debugOutputDir.empty()) options.diagnosticsDir = config.debugOutputDir + "/" + label;
            scrna_ed_result result = {};
            SimpleEDRunInfo info;
            const auto callStarted = std::chrono::steady_clock::now();
            timings[group].prepare = std::chrono::duration<double>(callStarted - groupStarted).count();
            SparseCountView view;
            view.genes = matrix.countCellGeneUMI.data();
            view.counts = matrix.countCellGeneUMI.data() + 1;
            view.geneWords = matrix.countCellGeneUMI.size();
            view.countWords = matrix.countCellGeneUMI.size() - 1;
            view.stride = matrix.countMatStride;
            view.offsets = starts.data();
            view.entries = nGenes.data();
            view.cells = barcodes.size();
            try {
                const int rc = runSimpleEDWithAmbientView(barcodes, umi, view,
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
            timings[group].call = std::chrono::duration<double>(std::chrono::steady_clock::now() - callStarted).count();
        };
        if (config.useThreadPermits) {
            scrna::ThreadPermitPool permits(totalThreads);
            // Reserve one base worker per coordinator before any group can
            // borrow helpers. Unstarted/error-path reservations release by RAII.
            std::vector<std::unique_ptr<scrna::ThreadPermit>> bases;
            for (size_t i = 0; i < groupWorkers; ++i)
                bases.emplace_back(new scrna::ThreadPermit(permits));
            std::atomic<size_t> nextGroup(0);
            scrna::parallelFor(groupWorkers, groupWorkers, [&](size_t slot, size_t worker) {
                auto base = std::move(bases[slot]);
                scrna::ScopedThreadPermits context(&permits);
                for (;;) {
                    const size_t group = nextGroup.fetch_add(1);
                    if (group >= groups.size()) break;
                    callGroup(group, worker);
                }
            });
            bases.clear();
            std::cerr << "[Flex caller permits] capacity=" << totalThreads
                << " peak_reserved=" << permits.peakUsed() << " acquisitions=" << permits.acquisitions()
                << " returned=" << permits.available() << "\n";
        } else {
            scrna::parallelFor(groups.size(), groupWorkers, callGroup);
        }
        for (size_t group = 0; group < groups.size(); ++group) {
            const auto& timing = timings[group];
            std::cerr << "[Flex tag-aware timing] sample=" << groups[group]
                << " tags=" << groupTags.at(groups[group]).size() << " worker_limit=" << timing.workers
                << " mc_threads=" << timing.mc << " prepare_seconds=" << timing.prepare
                << " caller_seconds=" << timing.call << "\n";
        }
        std::cerr << "[Flex tag-aware timing] groups_wall_seconds="
            << std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count() << "\n";
        // Occupancy is a joint fit: wait for every independent sample model.
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
