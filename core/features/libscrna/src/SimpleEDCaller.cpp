#include "SimpleEDCaller.h"
#include "OrdMagRank.h"
#include "EmptyDropsMultinomial.h"
#include "AdaptiveAmbientWindow.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <sys/stat.h>
#include <cerrno>

using namespace std;

namespace {
char* strdup_cpp(const string& value) {
    char* copy = static_cast<char*>(std::malloc(value.size() + 1));
    if (copy) std::memcpy(copy, value.c_str(), value.size() + 1);
    return copy;
}

void ensureDiagnosticsDirectory(const string& path) {
    string current;
    for (size_t i = 0; i <= path.size(); ++i) {
        if (i == path.size() || path[i] == '/') {
            if (!current.empty() && mkdir(current.c_str(), 0775) != 0 && errno != EEXIST)
                throw runtime_error("Cannot create caller diagnostics directory: " + current);
        }
        if (i < path.size()) current += path[i];
    }
}

void writeSimpleEDDiagnostics(const SimpleEDOptions& options, const scrna_ed_config& config,
    const vector<string>& barcodes, const vector<uint32_t>& retained,
    const vector<uint32_t>& genes, const vector<uint64_t>& nonMT,
    const vector<uint32_t>& umi, const vector<string>& simpleBarcodes,
    const SimpleEmptyDropsResult& simple,
    const vector<uint32_t>& ambient, const vector<uint32_t>& ambientCounts,
    const AmbientProfile& profile, const OrdMagBootstrapTrace& trace,
    const vector<uint8_t>& mitochondrial, const scrna_ed_result& result)
{
    ensureDiagnosticsDirectory(options.diagnosticsDir);
    ofstream ranks(options.diagnosticsDir + "/ordmag_rank.tsv");
    ofstream counts(options.diagnosticsDir + "/ambient_profile.tsv");
    ofstream summary(options.diagnosticsDir + "/caller_diagnostics.json");
    if (!ranks || !counts || !summary) throw runtime_error("Cannot write caller diagnostics");
    vector<uint32_t> order(umi.size());
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
        return ordMagRankBefore(a, b, umi, &genes, &simpleBarcodes,
                               nonMT.empty() ? nullptr : &nonMT);
    });
    vector<uint32_t> qualityRank(umi.size());
    for (uint32_t r = 0; r < order.size(); ++r) qualityRank[order[r]] = r + 1;
    vector<uint8_t> selected(retained.size(), 0), candidate(retained.size(), 0), background(retained.size(), 0);
    for (uint32_t r : simple.passingIndices) selected.at(r) = 1;
    for (uint32_t r : simple.candidateIndices) candidate.at(r) = 1;
    for (uint32_t r : ambient) background.at(r) = 1;
    ranks << "barcode\tretain_rank\tquality_rank\ttotal_umi\tnon_mt_umi\tdetected_genes\tis_simple\tis_candidate\tis_ambient\n";
    for (size_t r = 0; r < retained.size(); ++r) {
        ranks << barcodes.at(retained[r]) << '\t' << r + 1 << '\t';
        if (r < umi.size()) {
            ranks << qualityRank[r] << '\t' << umi[r] << '\t'
                  << (nonMT.empty() ? umi[r] : nonMT[r]) << '\t' << genes[r];
        } else {
            ranks << "NA\tNA\tNA\tNA";
        }
        ranks << '\t' << unsigned(selected[r]) << '\t' << unsigned(candidate[r])
              << '\t' << unsigned(background[r]) << '\n';
    }
    counts << "feature_index\tambient_umi\tlog_probability\tis_mitochondrial\n" << setprecision(17);
    uint64_t ambientMass = 0;
    for (uint32_t g = 0; g < ambientCounts.size(); ++g) {
        ambientMass += ambientCounts[g];
        counts << g << '\t' << ambientCounts[g] << '\t';
        if (g < profile.ambProfileLogP.size()) counts << profile.ambProfileLogP[g];
        else counts << "NA";
        counts << '\t' << unsigned(mitochondrial.empty() ? 0 : mitochondrial[g]) << '\n';
    }
    summary << setprecision(17) << "{\n"
        << "  \"protocol\": \"tag-aware-quality-ordmag-v1\",\n"
        << "  \"input_cells\": " << barcodes.size() << ",\n"
        << "  \"retain_count\": " << retained.size() << ",\n"
        << "  \"ordmag_input_cells\": " << umi.size() << ",\n"
        << "  \"ambient_start_parameter\": " << config.ind_min << ",\n"
        << "  \"ambient_base_end_parameter\": " << config.ind_max << ",\n"
        << "  \"ambient_cells\": " << ambient.size() << ",\n"
        << "  \"ambient_umi\": " << ambientMass << ",\n"
        << "  \"ambient_umi_target\": " << options.ambientUmiTarget << ",\n"
        << "  \"max_expected_cells\": " << options.maxExpectedCells << ",\n"
        << "  \"recovered_cells\": " << trace.recoveredCells << ",\n"
        << "  \"bootstrap_threads\": " << trace.bootstrapThreads << ",\n"
        << "  \"bootstrap_mean\": " << trace.meanRetained << ",\n"
        << "  \"bootstrap_sd\": " << trace.sdRetained << ",\n"
        << "  \"simple_cells\": " << simple.nCellsSimple << ",\n"
        << "  \"retain_umi\": " << simple.retainThreshold << ",\n"
        << "  \"candidate_umi_floor\": " << simple.minUMI << ",\n"
        << "  \"candidate_cap\": " << config.cand_max_n << ",\n"
        << "  \"candidates\": " << result.n_candidates << ",\n"
        << "  \"tail_candidates\": " << result.n_tail_cells << ",\n"
        << "  \"tail_passers\": " << result.n_ed_passers << ",\n"
        << "  \"final_cells\": " << result.n_barcodes << ",\n"
        << "  \"mc_threads\": " << config.mc_threads << ",\n"
        << "  \"simulations\": " << config.sim_n << ",\n"
        << "  \"fdr\": " << config.fdr << ",\n"
        << "  \"apply_bh\": " << config.apply_bh_correction << ",\n"
        << "  \"gate_on_fdr\": " << config.use_fdr_gate << ",\n"
        << "  \"mt_feature_rows\": " << count(mitochondrial.begin(), mitochondrial.end(), uint8_t(1)) << "\n}\n";
    ranks.close(); counts.close(); summary.close();
    if (!ranks || !counts || !summary) throw runtime_error("Incomplete caller diagnostics write");
}


} // namespace

int runSimpleEDWithAmbient(const vector<string>& barcodes,
                           const vector<uint32_t>& umi_counts,
                           const vector<uint32_t>& sparse_gene_ids,
                           const vector<uint32_t>& sparse_counts,
                           const vector<uint32_t>& sparse_cell_index,
                           const vector<uint32_t>& n_genes_per_cell,
                           uint32_t n_features,
                           const scrna_ed_config* config,
                           const SimpleEDOptions& options,
                           const vector<uint8_t>& mitochondrial_features,
                           scrna_ed_result* result,
                           SimpleEDRunInfo* info) {
    const bool use_legacy_rank_ambient = options.legacyRankAmbient;
    const bool use_guarded_rank_ambient = options.guardedRankAmbient;
    const uint32_t ambient_fallback_min_abs = options.ambientFallbackMinAbs;
    const double ambient_fallback_min_frac = options.ambientFallbackMinFrac;
    const uint32_t max_expected_cells = options.maxExpectedCells;
    const uint32_t ordmag_retain_count = options.ordmagRetainCount;
    const uint64_t ambient_umi_target = options.ambientUmiTarget;
    if (!config || !result) {
        return -1;
    }

    std::vector<std::pair<uint32_t, uint32_t>> umi_idx;
    umi_idx.reserve(umi_counts.size());
    for (uint32_t i = 0; i < umi_counts.size(); i++) {
        umi_idx.push_back({umi_counts[i], i});
    }
    std::stable_sort(umi_idx.begin(), umi_idx.end(),
                     [&barcodes](const std::pair<uint32_t, uint32_t>& a,
                        const std::pair<uint32_t, uint32_t>& b) {
        if (a.first != b.first) return a.first > b.first;
        return barcodes[a.second] < barcodes[b.second];
    });

    uint32_t retain_count = (config->ed_retain_count > 0)
        ? std::min<uint32_t>(config->ed_retain_count, umi_counts.size())
        : static_cast<uint32_t>(umi_counts.size());

    AdaptiveAmbientWindow adaptiveWindow;
    adaptiveWindow.start = std::min<uint32_t>(config->ind_min, umi_counts.size());
    adaptiveWindow.end = retain_count;
    if (ambient_umi_target > 0) {
        adaptiveWindow = selectAdaptiveAmbientWindow(
            umi_idx, config->ind_min, retain_count, ambient_umi_target);
        retain_count = adaptiveWindow.end;
        std::cerr << "[scrna_simpleed] Adaptive ambient window: ["
                  << adaptiveWindow.start << ", " << adaptiveWindow.end
                  << "), mass=" << adaptiveWindow.umiMass
                  << ", target=" << ambient_umi_target << "\n";
    }

    std::vector<uint32_t> retain_indices;
    std::vector<uint32_t> retain_umi;
    retain_indices.reserve(retain_count);
    retain_umi.reserve(retain_count);
    for (uint32_t i = 0; i < retain_count; i++) {
        uint32_t orig_idx = umi_idx[i].second;
        retain_indices.push_back(orig_idx);
        retain_umi.push_back(umi_counts[orig_idx]);
    }

    const uint32_t simple_count = ordmag_retain_count > 0
        ? std::min<uint32_t>(ordmag_retain_count, retain_umi.size())
        : static_cast<uint32_t>(retain_umi.size());
    vector<uint32_t> simple_umi(retain_umi.begin(), retain_umi.begin() + simple_count);
    vector<uint32_t> simple_genes;
    vector<string> simple_barcodes;
    vector<uint64_t> simple_non_mito_umis;
    simple_genes.reserve(simple_count);
    simple_barcodes.reserve(simple_count);
    if (!mitochondrial_features.empty()) simple_non_mito_umis.reserve(simple_count);
    vector<uint32_t> seen_genes(n_features, std::numeric_limits<uint32_t>::max());
    for (uint32_t rank = 0; rank < simple_count; ++rank) {
        const uint32_t cell = retain_indices[rank];
        const uint32_t start = sparse_cell_index[cell];
        const OrdMagCellQuality quality = ordMagCellQuality(
            sparse_gene_ids.data() + start, sparse_counts.data() + start,
            n_genes_per_cell[cell], seen_genes, rank,
            mitochondrial_features.empty() ? nullptr : mitochondrial_features.data());
        simple_genes.push_back(quality.detectedGenes);
        if (!mitochondrial_features.empty()) simple_non_mito_umis.push_back(quality.nonMitoUMIs);
        simple_barcodes.push_back(barcodes[cell]);
    }
    std::cerr << "[scrna_simpleed] OrdMag retain window: " << simple_count
              << "; ambient-accessible retain window: " << retain_count << "\n";

    SimpleEmptyDropsParams simple_params;
    simple_params.nExpectedCells = config->n_expected_cells;
    simple_params.maxPercentile = config->max_percentile;
    simple_params.maxMinRatio = config->max_min_ratio;
    simple_params.umiMin = config->umi_min;
    simple_params.primaryUmiMin = config->primary_umi_min;
    simple_params.umiMinFracMedian = config->umi_min_frac_median;
    simple_params.candMaxN = config->cand_max_n;
    simple_params.indMin = config->ind_min;
    simple_params.indMax = retain_count;

    simple_params.maxThreads = options.bootstrapThreads;
    simple_params.maxConcurrentThreads = options.bootstrapWorkers;
    OrdMagBootstrapTrace bootstrap_trace;
    SimpleEmptyDropsResult simple_result;
    if (config->use_bootstrap) {
        simple_params.useBootstrap = true;
        simple_params.nExpectedCells = 0;
        simple_params.maxExpectedCells = max_expected_cells > 0
            ? max_expected_cells
            : std::min(config->ind_min / 2, static_cast<uint32_t>(262144));
        if (simple_params.maxExpectedCells < 1000) {
            simple_params.maxExpectedCells = 90000;
        }
        simple_result = SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(
            simple_umi, simple_umi.size(), simple_params, simple_genes, simple_barcodes, simple_non_mito_umis, &bootstrap_trace);
    } else {
        simple_result = SimpleEmptyDropsStage::runCRSimpleFilter(
            simple_umi, simple_umi.size(), simple_params);
    }

    std::cerr << "[scrna_simpleed] Custom ambient simple filter: "
              << simple_result.nCellsSimple << " cells, threshold="
              << simple_result.retainThreshold << "\n";

    std::vector<uint32_t> ambient_retain_indices;
    if (use_legacy_rank_ambient) {
        uint32_t ambient_start = std::min<uint32_t>(config->ind_min, retain_indices.size());
        uint32_t ambient_end = ambient_umi_target > 0
            ? std::min<uint32_t>(adaptiveWindow.end, retain_indices.size())
            : std::min<uint32_t>(config->ind_max, retain_indices.size());
        ambient_retain_indices.reserve(ambient_end > ambient_start ? ambient_end - ambient_start : 0);
        for (uint32_t rank = ambient_start; rank < ambient_end; rank++) {
            ambient_retain_indices.push_back(rank);
        }
        std::cerr << "[scrna_simpleed] Ambient source: legacy rank window ["
                  << ambient_start << ", " << ambient_end << ")\n";
    } else if (use_guarded_rank_ambient) {
        uint32_t ambient_start = std::min<uint32_t>(config->ind_min, retain_indices.size());
        uint32_t ambient_end = std::min<uint32_t>(config->ind_max, retain_indices.size());
        uint32_t ambient_window_size = (ambient_end > ambient_start) ? (ambient_end - ambient_start) : 0;
        uint32_t min_from_frac = 0;
        if (ambient_fallback_min_frac > 0.0) {
            min_from_frac = static_cast<uint32_t>(
                ambient_fallback_min_frac * static_cast<double>(retain_indices.size()));
        }
        uint32_t required_ambient = std::max(ambient_fallback_min_abs, min_from_frac);
        if (required_ambient == 0) {
            required_ambient = std::min<uint32_t>(100, retain_indices.size());
        }

        if (ambient_end <= ambient_start || ambient_window_size < required_ambient) {
            uint32_t fallback_size = std::min<uint32_t>(required_ambient, retain_indices.size());
            uint32_t fallback_start = (retain_indices.size() >= fallback_size)
                ? (static_cast<uint32_t>(retain_indices.size()) - fallback_size) : 0;
            ambient_retain_indices.reserve(fallback_size);
            for (uint32_t rank = fallback_start; rank < retain_indices.size(); rank++) {
                ambient_retain_indices.push_back(rank);
            }
            std::cerr << "[scrna_simpleed] Ambient source: guarded rank fallback bottom "
                      << fallback_size << " cells (required=" << required_ambient
                      << ", frac=" << ambient_fallback_min_frac
                      << ", abs=" << ambient_fallback_min_abs << ")\n";
        } else {
            ambient_retain_indices.reserve(ambient_window_size);
            for (uint32_t rank = ambient_start; rank < ambient_end; rank++) {
                ambient_retain_indices.push_back(rank);
            }
            std::cerr << "[scrna_simpleed] Ambient source: guarded legacy rank window ["
                      << ambient_start << ", " << ambient_end << ")"
                      << " (required=" << required_ambient
                      << ", frac=" << ambient_fallback_min_frac
                      << ", abs=" << ambient_fallback_min_abs << ")\n";
        }
    } else {
        ambient_retain_indices = simple_result.ambientIndices;
        std::cerr << "[scrna_simpleed] Ambient source: SimpleED ambient set ("
                  << ambient_retain_indices.size() << " cells)\n";
    }

    std::vector<uint32_t> amb_count(n_features, 0);
    uint32_t ambient_cells_used = 0;
    for (uint32_t retain_idx : ambient_retain_indices) {
        if (retain_idx >= retain_indices.size()) {
            continue;
        }
        uint32_t orig_idx = retain_indices[retain_idx];
        uint32_t start = sparse_cell_index[orig_idx];
        uint32_t n_genes = n_genes_per_cell[orig_idx];
        for (uint32_t g = 0; g < n_genes; g++) {
            size_t pos = static_cast<size_t>(start + g);
            uint32_t gene_id = sparse_gene_ids[pos];
            uint32_t count = sparse_counts[pos];
            if (gene_id < n_features) {
                amb_count[gene_id] += count;
            }
        }
        ambient_cells_used++;
    }

    std::cerr << "[scrna_simpleed] Ambient cells used: " << ambient_cells_used << "\n";

    std::vector<uint32_t> feat_det_vec;
    feat_det_vec.reserve(n_features);
    for (uint32_t i = 0; i < n_features; i++) {
        if (amb_count[i] > 0) {
            feat_det_vec.push_back(i);
        }
    }

    AmbientProfile amb_profile = EmptyDropsMultinomial::computeAmbientProfile(
        amb_count, n_features, feat_det_vec, feat_det_vec.size());

    std::vector<uint32_t> count_cell_gene_umi;
    std::vector<uint32_t> count_cell_gene_umiindex(umi_counts.size(), 0);
    std::vector<uint32_t> n_gene_per_cb(umi_counts.size(), 0);
    count_cell_gene_umi.reserve(sparse_counts.size() * 2);
    for (uint32_t cell_idx = 0; cell_idx < umi_counts.size(); cell_idx++) {
        count_cell_gene_umiindex[cell_idx] = static_cast<uint32_t>(count_cell_gene_umi.size());
        uint32_t start = sparse_cell_index[cell_idx];
        uint32_t n_genes = n_genes_per_cell[cell_idx];
        n_gene_per_cb[cell_idx] = n_genes;
        for (uint32_t g = 0; g < n_genes; g++) {
            size_t pos = static_cast<size_t>(start + g);
            count_cell_gene_umi.push_back(sparse_gene_ids[pos]);
            count_cell_gene_umi.push_back(sparse_counts[pos]);
        }
    }

    std::vector<uint32_t> candidate_orig_indices;
    std::vector<uint32_t> candidate_counts;
    candidate_orig_indices.reserve(simple_result.candidateIndices.size());
    candidate_counts.reserve(simple_result.candidateIndices.size());
    for (uint32_t retain_idx : simple_result.candidateIndices) {
        if (retain_idx >= retain_indices.size()) {
            continue;
        }
        candidate_orig_indices.push_back(retain_indices[retain_idx]);
        candidate_counts.push_back(retain_umi[retain_idx]);
    }

    EmptyDropsParams ed_params;
    ed_params.indMin = config->ind_min;
    ed_params.indMax = retain_count;
    ed_params.umiMin = config->umi_min;
    ed_params.umiMinFracMedian = config->umi_min_frac_median;
    ed_params.candMaxN = config->cand_max_n;
    ed_params.FDR = config->fdr;
    ed_params.rawPvalueThreshold = config->raw_pvalue_threshold;
    ed_params.simN = config->sim_n;
    ed_params.seed = config->seed;
    ed_params.lowerTestingBound = config->lower_testing_bound;
    ed_params.ambientUmiMax = config->ambient_umi_max;
    ed_params.mcThreads = config->mc_threads;
    ed_params.applyBHCorrection = (config->apply_bh_correction != 0);

    std::vector<EmptyDropsResult> ed_results = EmptyDropsMultinomial::computePValues(
        amb_profile,
        candidate_orig_indices,
        candidate_counts,
        count_cell_gene_umi,
        count_cell_gene_umiindex,
        n_gene_per_cb,
        2,
        1,
        ed_params,
        simple_result.nCellsSimple,
        std::vector<string>(),
        static_cast<uint32_t>(umi_counts.size()),
        "",
        "",
        options.invariantChecks
    );

    std::vector<string> passing_barcodes;
    passing_barcodes.reserve(simple_result.passingIndices.size() + ed_results.size());
    std::vector<uint8_t> simple_flags(umi_counts.size(), 0);
    for (uint32_t retain_idx : simple_result.passingIndices) {
        if (retain_idx >= retain_indices.size()) {
            continue;
        }
        uint32_t orig_idx = retain_indices[retain_idx];
        simple_flags[orig_idx] = 1;
        passing_barcodes.push_back(barcodes[orig_idx]);
    }

    uint32_t n_ed_passers = 0;
    for (const auto& ed_result : ed_results) {
        bool passes = config->use_fdr_gate ? ed_result.passesFDR : ed_result.passesRawP;
        if (!passes) {
            continue;
        }
        uint32_t orig_idx = ed_result.cellIndex;
        if (orig_idx >= barcodes.size() || simple_flags[orig_idx]) {
            continue;
        }
        passing_barcodes.push_back(barcodes[orig_idx]);
        n_ed_passers++;
    }

    result->n_barcodes = passing_barcodes.size();
    result->barcodes = static_cast<char**>(std::malloc(result->n_barcodes * sizeof(char*)));
    if (!result->barcodes && result->n_barcodes > 0) {
        result->error_message = strdup_cpp("Memory allocation failed");
        return -1;
    }
    for (size_t i = 0; i < passing_barcodes.size(); i++) {
        result->barcodes[i] = strdup_cpp(passing_barcodes[i]);
    }

    result->n_candidates = ed_results.size();
    result->candidates = static_cast<scrna_ed_candidate*>(
        std::malloc(result->n_candidates * sizeof(scrna_ed_candidate)));
    if (!result->candidates && result->n_candidates > 0) {
        result->error_message = strdup_cpp("Memory allocation failed");
        return -1;
    }
    for (size_t i = 0; i < ed_results.size(); i++) {
        uint32_t orig_idx = ed_results[i].cellIndex;
        result->candidates[i].cell_index = orig_idx;
        result->candidates[i].barcode = (orig_idx < barcodes.size())
            ? strdup_cpp(barcodes[orig_idx]) : nullptr;
        result->candidates[i].umi_count = (orig_idx < umi_counts.size()) ? umi_counts[orig_idx] : 0;
        result->candidates[i].p_value = ed_results[i].pValue;
        result->candidates[i].p_adjusted = ed_results[i].pAdjusted;
        result->candidates[i].passes_raw_p = ed_results[i].passesRawP ? 1 : 0;
        result->candidates[i].passes_fdr = ed_results[i].passesFDR ? 1 : 0;
        result->candidates[i].obs_log_prob = ed_results[i].obsLogProb;
        result->candidates[i].is_simple_cell = (orig_idx < simple_flags.size() && simple_flags[orig_idx]) ? 1 : 0;
    }

    result->n_simple_cells = simple_result.nCellsSimple;
    result->n_tail_cells = ed_results.size() >= simple_result.nCellsSimple
        ? ed_results.size() - simple_result.nCellsSimple : 0;
    result->n_ed_passers = n_ed_passers;
    result->retain_threshold = simple_result.retainThreshold;
    result->min_umi = simple_result.minUMI;

    std::cerr << "[scrna_simpleed] Custom ambient passers: " << result->n_barcodes << "\n";
    if (info) {
        info->retainCount = retain_count;
        info->ordmagCount = simple_count;
        info->ambientCells = ambient_cells_used;
        info->bootstrap = bootstrap_trace;
    }
    if (!options.diagnosticsDir.empty()) {
        writeSimpleEDDiagnostics(options, *config, barcodes, retain_indices, simple_genes,
            simple_non_mito_umis, simple_umi, simple_barcodes, simple_result,
            ambient_retain_indices, amb_count, amb_profile, bootstrap_trace,
            mitochondrial_features, *result);
    }
    return 0;
}
