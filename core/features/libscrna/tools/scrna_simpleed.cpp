#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "EmptyDropsMultinomial.h"
#include "OrdMagStage.h"
#include "scrna_api.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;

namespace {

struct Args {
    string matrix_path;
    string barcodes_path;
    string out_barcodes;
    string out_dir;
    string mode = "simple";
    bool include_zero_umis = false;
    uint32_t expected_cells = 0;
    uint32_t umi_min = 0;
    uint32_t sim_n = 0;
    uint32_t ed_retain_count = 0;
    uint32_t lower_testing_bound = 0;
    uint32_t ambient_umi_max = 0;
    double fdr = 0.0;
    double raw_pvalue = 0.0;
    bool use_fdr_gate = false;
    bool apply_bh_correction = false;
    bool use_bootstrap = false;
    bool direct_ed_surface = false;
    bool use_legacy_rank_ambient = false;
    bool use_guarded_rank_ambient = false;
    uint32_t ambient_fallback_min_abs = 0;
    double ambient_fallback_min_frac = 0.0;
};

void usage(const char* prog) {
    cerr << "Usage: " << prog << " --matrix matrix.mtx --barcodes barcodes.tsv --out barcodes.tsv"
         << " [--mode simple|full] [--out-dir DIR] [--expected-cells N] [--umi-min N]"
         << " [--sim-n N] [--ed-retain-count N] [--lower-testing-bound N]"
         << " [--ambient-umi-max N] [--fdr X] [--raw-pvalue X] [--use-fdr-gate]"
         << " [--apply-bh-correction] [--use-bootstrap] [--direct-ed-surface]"
         << " [--use-legacy-rank-ambient] [--use-guarded-rank-ambient]"
         << " [--ambient-fallback-min-abs N] [--ambient-fallback-min-frac X]"
         << " [--include-zero-umis]\n";
}

bool parse_uint32(const string& s, uint32_t* out) {
    if (!out) return false;
    char* end = nullptr;
    errno = 0;
    unsigned long v = std::strtoul(s.c_str(), &end, 10);
    if (errno != 0 || end == s.c_str() || *end != '\0' || v > std::numeric_limits<uint32_t>::max()) {
        return false;
    }
    *out = static_cast<uint32_t>(v);
    return true;
}

bool parse_double(const string& s, double* out) {
    if (!out) return false;
    char* end = nullptr;
    errno = 0;
    double v = std::strtod(s.c_str(), &end);
    if (errno != 0 || end == s.c_str() || *end != '\0') {
        return false;
    }
    *out = v;
    return true;
}

bool read_barcodes(const string& path, vector<string>* out) {
    if (!out) return false;
    ifstream in(path);
    if (!in.is_open()) {
        return false;
    }
    string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        out->push_back(line);
    }
    return true;
}

bool read_matrix_counts(const string& path, uint32_t* out_rows, uint32_t* out_cols, vector<uint64_t>* out_counts) {
    if (!out_rows || !out_cols || !out_counts) return false;
    ifstream in(path);
    if (!in.is_open()) {
        return false;
    }

    string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '%') continue;
        std::stringstream header(line);
        uint64_t rows = 0, cols = 0, nnz = 0;
        if (!(header >> rows >> cols >> nnz)) {
            return false;
        }
        if (rows == 0 || cols == 0) {
            return false;
        }
        if (rows > std::numeric_limits<uint32_t>::max() || cols > std::numeric_limits<uint32_t>::max()) {
            return false;
        }
        *out_rows = static_cast<uint32_t>(rows);
        *out_cols = static_cast<uint32_t>(cols);
        out_counts->assign(*out_cols, 0);

        uint64_t row = 0, col = 0;
        double val = 0.0;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            if (!(ss >> row >> col >> val)) {
                continue;
            }
            if (col == 0 || col > *out_cols) {
                continue;
            }
            uint64_t add = (val < 0.0) ? 0 : static_cast<uint64_t>(val + 0.5);
            (*out_counts)[static_cast<size_t>(col - 1)] += add;
        }
        return true;
    }
    return false;
}

bool read_matrix_sparse(const string& path,
                        uint32_t* out_rows,
                        uint32_t* out_cols,
                        vector<uint32_t>* umi_counts,
                        vector<uint32_t>* sparse_gene_ids,
                        vector<uint32_t>* sparse_counts,
                        vector<uint32_t>* sparse_cell_index,
                        vector<uint32_t>* n_genes_per_cell,
                        size_t* out_nnz) {
    if (!out_rows || !out_cols || !umi_counts || !sparse_gene_ids || !sparse_counts || !sparse_cell_index || !n_genes_per_cell || !out_nnz) {
        return false;
    }

    FILE* f = fopen(path.c_str(), "r");
    if (!f) {
        return false;
    }

    char line[1024];
    // Skip comments
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '%') break;
    }

    uint64_t nrows = 0, ncols = 0, nnz = 0;
    if (sscanf(line, "%lu %lu %lu", &nrows, &ncols, &nnz) != 3) {
        fclose(f);
        return false;
    }
    if (nrows == 0 || ncols == 0) {
        fclose(f);
        return false;
    }
    if (nrows > std::numeric_limits<uint32_t>::max() || ncols > std::numeric_limits<uint32_t>::max()) {
        fclose(f);
        return false;
    }

    *out_rows = static_cast<uint32_t>(nrows);
    *out_cols = static_cast<uint32_t>(ncols);
    *out_nnz = static_cast<size_t>(nnz);

    sparse_gene_ids->assign(*out_nnz, 0);
    sparse_counts->assign(*out_nnz, 0);

    vector<uint32_t> cell_starts(*out_cols + 1, 0);
    long data_pos = ftell(f);

    // First pass: count per cell
    while (fgets(line, sizeof(line), f)) {
        uint32_t row = 0, col = 0, val = 0;
        if (sscanf(line, "%u %u %u", &row, &col, &val) == 3) {
            if (col >= 1 && col <= *out_cols) {
                cell_starts[col]++;
            }
        }
    }

    sparse_cell_index->assign(*out_cols + 1, 0);
    n_genes_per_cell->assign(*out_cols, 0);
    (*sparse_cell_index)[0] = 0;
    for (uint32_t i = 0; i < *out_cols; i++) {
        (*n_genes_per_cell)[i] = cell_starts[i + 1];
        (*sparse_cell_index)[i + 1] = (*sparse_cell_index)[i] + cell_starts[i + 1];
    }

    std::fill(cell_starts.begin(), cell_starts.end(), 0);
    fseek(f, data_pos, SEEK_SET);
    umi_counts->assign(*out_cols, 0);

    while (fgets(line, sizeof(line), f)) {
        uint32_t row = 0, col = 0, val = 0;
        if (sscanf(line, "%u %u %u", &row, &col, &val) == 3) {
            if (col >= 1 && col <= *out_cols && row >= 1 && row <= *out_rows) {
                uint32_t cell_idx = col - 1;
                uint32_t gene_idx = row - 1;
                uint32_t write_pos = (*sparse_cell_index)[cell_idx] + cell_starts[col];
                if (write_pos < sparse_gene_ids->size()) {
                    (*sparse_gene_ids)[write_pos] = gene_idx;
                    (*sparse_counts)[write_pos] = val;
                }
                cell_starts[col]++;
                (*umi_counts)[cell_idx] += val;
            }
        }
    }

    fclose(f);
    return true;
}

char* strdup_cpp(const string& s) {
    char* out = static_cast<char*>(std::malloc(s.size() + 1));
    if (!out) {
        return nullptr;
    }
    std::memcpy(out, s.c_str(), s.size() + 1);
    return out;
}

int run_direct_ed_surface(const vector<string>& barcodes,
                          const vector<uint32_t>& umi_counts,
                          const vector<uint32_t>& sparse_gene_ids,
                          const vector<uint32_t>& sparse_counts,
                          const vector<uint32_t>& sparse_cell_index,
                          const vector<uint32_t>& n_genes_per_cell,
                          uint32_t n_features,
                          const scrna_ed_config* config,
                          scrna_ed_result* result) {
    if (!config || !result) {
        return -1;
    }

    std::vector<std::pair<uint32_t, uint32_t>> umi_idx;
    umi_idx.reserve(umi_counts.size());
    for (uint32_t i = 0; i < umi_counts.size(); i++) {
        umi_idx.push_back({umi_counts[i], i});
    }
    std::stable_sort(umi_idx.begin(), umi_idx.end(),
                     [](const std::pair<uint32_t, uint32_t>& a,
                        const std::pair<uint32_t, uint32_t>& b) {
        return a.first > b.first;
    });

    uint32_t retain_count = (config->ed_retain_count > 0)
        ? std::min<uint32_t>(config->ed_retain_count, umi_counts.size())
        : static_cast<uint32_t>(umi_counts.size());

    std::vector<uint32_t> retain_indices;
    retain_indices.reserve(retain_count);
    std::vector<string> retain_barcodes;
    retain_barcodes.reserve(retain_count);
    std::vector<uint32_t> retain_umi;
    retain_umi.reserve(retain_count);

    for (uint32_t i = 0; i < retain_count; i++) {
        uint32_t orig_idx = umi_idx[i].second;
        retain_indices.push_back(orig_idx);
        retain_barcodes.push_back(barcodes[orig_idx]);
        retain_umi.push_back(umi_counts[orig_idx]);
    }

    std::cerr << "[scrna_simpleed] Direct ED retain window: " << retain_indices.size()
              << " cells\n";

    std::vector<uint32_t> count_cell_gene_umi;
    std::vector<uint32_t> count_cell_gene_umiindex(retain_indices.size(), 0);
    std::vector<uint32_t> n_gene_per_cb(retain_indices.size(), 0);
    count_cell_gene_umi.reserve(sparse_counts.size() * 2);

    for (size_t local_idx = 0; local_idx < retain_indices.size(); local_idx++) {
        uint32_t orig_idx = retain_indices[local_idx];
        count_cell_gene_umiindex[local_idx] = static_cast<uint32_t>(count_cell_gene_umi.size());
        uint32_t start = sparse_cell_index[orig_idx];
        uint32_t n_genes = n_genes_per_cell[orig_idx];
        n_gene_per_cb[local_idx] = n_genes;
        for (uint32_t g = 0; g < n_genes; g++) {
            size_t pos = static_cast<size_t>(start + g);
            count_cell_gene_umi.push_back(sparse_gene_ids[pos]);
            count_cell_gene_umi.push_back(sparse_counts[pos]);
        }
    }

    std::vector<uint32_t> amb_count(n_features, 0);
    uint32_t ambient_cells_used = 0;
    for (size_t local_idx = 0; local_idx < retain_indices.size(); local_idx++) {
        if (retain_umi[local_idx] > config->ambient_umi_max) {
            continue;
        }
        uint32_t start = count_cell_gene_umiindex[local_idx];
        uint32_t n_genes = n_gene_per_cb[local_idx];
        for (uint32_t g = 0; g < n_genes; g++) {
            uint32_t gene_id = count_cell_gene_umi[start + g * 2];
            uint32_t count = count_cell_gene_umi[start + g * 2 + 1];
            if (gene_id < n_features) {
                amb_count[gene_id] += count;
            }
        }
        ambient_cells_used++;
    }

    std::cerr << "[scrna_simpleed] Direct ED ambient cells (UMI <= "
              << config->ambient_umi_max << "): " << ambient_cells_used << "\n";

    AmbientProfile amb_profile = EmptyDropsMultinomial::computeAmbientProfile(
        amb_count, n_features, std::vector<uint32_t>(), 0);

    std::vector<uint32_t> candidate_indices;
    std::vector<uint32_t> candidate_counts;
    candidate_indices.reserve(retain_indices.size());
    candidate_counts.reserve(retain_indices.size());
    for (size_t local_idx = 0; local_idx < retain_indices.size(); local_idx++) {
        if (retain_umi[local_idx] > config->lower_testing_bound) {
            candidate_indices.push_back(static_cast<uint32_t>(local_idx));
            candidate_counts.push_back(retain_umi[local_idx]);
        }
    }

    std::cerr << "[scrna_simpleed] Direct ED candidates (UMI > "
              << config->lower_testing_bound << "): " << candidate_indices.size() << "\n";

    EmptyDropsParams ed_params;
    ed_params.indMin = config->ind_min;
    ed_params.indMax = config->ind_max;
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
        candidate_indices,
        candidate_counts,
        count_cell_gene_umi,
        count_cell_gene_umiindex,
        n_gene_per_cb,
        2,
        1,
        ed_params,
        0,
        std::vector<string>(),
        static_cast<uint32_t>(retain_indices.size()),
        "",
        "",
        false
    );

    std::vector<string> passing_barcodes;
    passing_barcodes.reserve(ed_results.size());
    uint32_t n_ed_passers = 0;
    for (const auto& ed_result : ed_results) {
        bool passes = config->use_fdr_gate ? ed_result.passesFDR : ed_result.passesRawP;
        if (!passes) {
            continue;
        }
        uint32_t local_idx = ed_result.cellIndex;
        if (local_idx >= retain_barcodes.size()) {
            continue;
        }
        passing_barcodes.push_back(retain_barcodes[local_idx]);
        n_ed_passers++;
    }

    result->n_barcodes = passing_barcodes.size();
    result->barcodes = static_cast<char**>(std::malloc(result->n_barcodes * sizeof(char*)));
    if (!result->barcodes) {
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
        uint32_t local_idx = ed_results[i].cellIndex;
        uint32_t orig_idx = (local_idx < retain_indices.size()) ? retain_indices[local_idx] : UINT32_MAX;
        result->candidates[i].cell_index = orig_idx;
        result->candidates[i].barcode = (orig_idx < barcodes.size())
            ? strdup_cpp(barcodes[orig_idx]) : nullptr;
        result->candidates[i].umi_count = (orig_idx < umi_counts.size()) ? umi_counts[orig_idx] : 0;
        result->candidates[i].p_value = ed_results[i].pValue;
        result->candidates[i].p_adjusted = ed_results[i].pAdjusted;
        result->candidates[i].passes_raw_p = ed_results[i].passesRawP ? 1 : 0;
        result->candidates[i].passes_fdr = ed_results[i].passesFDR ? 1 : 0;
        result->candidates[i].obs_log_prob = ed_results[i].obsLogProb;
        result->candidates[i].is_simple_cell = 0;
    }

    result->n_simple_cells = 0;
    result->n_tail_cells = ed_results.size();
    result->n_ed_passers = n_ed_passers;
    result->retain_threshold = 0;
    result->min_umi = config->lower_testing_bound;

    std::cerr << "[scrna_simpleed] Direct ED passers: " << result->n_barcodes << "\n";
    return 0;
}

int run_simpleed_custom_ambient(const vector<string>& barcodes,
                                const vector<uint32_t>& umi_counts,
                                const vector<uint32_t>& sparse_gene_ids,
                                const vector<uint32_t>& sparse_counts,
                                const vector<uint32_t>& sparse_cell_index,
                                const vector<uint32_t>& n_genes_per_cell,
                                uint32_t n_features,
                                const scrna_ed_config* config,
                                bool use_legacy_rank_ambient,
                                bool use_guarded_rank_ambient,
                                uint32_t ambient_fallback_min_abs,
                                double ambient_fallback_min_frac,
                                scrna_ed_result* result) {
    if (!config || !result) {
        return -1;
    }

    std::vector<std::pair<uint32_t, uint32_t>> umi_idx;
    umi_idx.reserve(umi_counts.size());
    for (uint32_t i = 0; i < umi_counts.size(); i++) {
        umi_idx.push_back({umi_counts[i], i});
    }
    std::stable_sort(umi_idx.begin(), umi_idx.end(),
                     [](const std::pair<uint32_t, uint32_t>& a,
                        const std::pair<uint32_t, uint32_t>& b) {
        return a.first > b.first;
    });

    uint32_t retain_count = (config->ed_retain_count > 0)
        ? std::min<uint32_t>(config->ed_retain_count, umi_counts.size())
        : static_cast<uint32_t>(umi_counts.size());

    std::vector<uint32_t> retain_indices;
    std::vector<uint32_t> retain_umi;
    retain_indices.reserve(retain_count);
    retain_umi.reserve(retain_count);
    for (uint32_t i = 0; i < retain_count; i++) {
        uint32_t orig_idx = umi_idx[i].second;
        retain_indices.push_back(orig_idx);
        retain_umi.push_back(umi_counts[orig_idx]);
    }

    SimpleEmptyDropsParams simple_params;
    simple_params.nExpectedCells = config->n_expected_cells;
    simple_params.maxPercentile = config->max_percentile;
    simple_params.maxMinRatio = config->max_min_ratio;
    simple_params.umiMin = config->umi_min;
    simple_params.umiMinFracMedian = config->umi_min_frac_median;
    simple_params.candMaxN = config->cand_max_n;
    simple_params.indMin = config->ind_min;
    simple_params.indMax = config->ind_max;

    SimpleEmptyDropsResult simple_result;
    if (config->use_bootstrap) {
        simple_params.useBootstrap = true;
        simple_params.nExpectedCells = 0;
        simple_params.maxExpectedCells = std::min(config->ind_min / 2, static_cast<uint32_t>(262144));
        if (simple_params.maxExpectedCells < 1000) {
            simple_params.maxExpectedCells = 90000;
        }
        simple_result = SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(
            retain_umi, retain_indices.size(), simple_params);
    } else {
        simple_result = SimpleEmptyDropsStage::runCRSimpleFilter(
            retain_umi, retain_indices.size(), simple_params);
    }

    std::cerr << "[scrna_simpleed] Custom ambient simple filter: "
              << simple_result.nCellsSimple << " cells, threshold="
              << simple_result.retainThreshold << "\n";

    std::vector<uint32_t> ambient_retain_indices;
    if (use_legacy_rank_ambient) {
        uint32_t ambient_start = std::min<uint32_t>(config->ind_min, retain_indices.size());
        uint32_t ambient_end = std::min<uint32_t>(config->ind_max, retain_indices.size());
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
    ed_params.indMax = config->ind_max;
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
        false
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
        result->candidates[i].is_simple_cell = 0;
    }

    result->n_simple_cells = simple_result.nCellsSimple;
    result->n_tail_cells = ed_results.size();
    result->n_ed_passers = n_ed_passers;
    result->retain_threshold = simple_result.retainThreshold;
    result->min_umi = simple_result.minUMI;

    std::cerr << "[scrna_simpleed] Custom ambient passers: " << result->n_barcodes << "\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; i++) {
        string key = argv[i];
        if (key == "--matrix" && i + 1 < argc) {
            args.matrix_path = argv[++i];
        } else if (key == "--barcodes" && i + 1 < argc) {
            args.barcodes_path = argv[++i];
        } else if (key == "--out" && i + 1 < argc) {
            args.out_barcodes = argv[++i];
        } else if (key == "--out-dir" && i + 1 < argc) {
            args.out_dir = argv[++i];
        } else if (key == "--mode" && i + 1 < argc) {
            args.mode = argv[++i];
        } else if (key == "--include-zero-umis") {
            args.include_zero_umis = true;
        } else if (key == "--expected-cells" && i + 1 < argc) {
            if (!parse_uint32(argv[++i], &args.expected_cells)) {
                cerr << "Invalid --expected-cells value\n";
                return 2;
            }
        } else if (key == "--umi-min" && i + 1 < argc) {
            if (!parse_uint32(argv[++i], &args.umi_min)) {
                cerr << "Invalid --umi-min value\n";
                return 2;
            }
        } else if (key == "--sim-n" && i + 1 < argc) {
            if (!parse_uint32(argv[++i], &args.sim_n)) {
                cerr << "Invalid --sim-n value\n";
                return 2;
            }
        } else if (key == "--ed-retain-count" && i + 1 < argc) {
            if (!parse_uint32(argv[++i], &args.ed_retain_count)) {
                cerr << "Invalid --ed-retain-count value\n";
                return 2;
            }
        } else if (key == "--lower-testing-bound" && i + 1 < argc) {
            if (!parse_uint32(argv[++i], &args.lower_testing_bound)) {
                cerr << "Invalid --lower-testing-bound value\n";
                return 2;
            }
        } else if (key == "--ambient-umi-max" && i + 1 < argc) {
            if (!parse_uint32(argv[++i], &args.ambient_umi_max)) {
                cerr << "Invalid --ambient-umi-max value\n";
                return 2;
            }
        } else if (key == "--ambient-fallback-min-abs" && i + 1 < argc) {
            if (!parse_uint32(argv[++i], &args.ambient_fallback_min_abs)) {
                cerr << "Invalid --ambient-fallback-min-abs value\n";
                return 2;
            }
        } else if (key == "--ambient-fallback-min-frac" && i + 1 < argc) {
            if (!parse_double(argv[++i], &args.ambient_fallback_min_frac)) {
                cerr << "Invalid --ambient-fallback-min-frac value\n";
                return 2;
            }
        } else if (key == "--fdr" && i + 1 < argc) {
            args.fdr = std::atof(argv[++i]);
        } else if (key == "--raw-pvalue" && i + 1 < argc) {
            args.raw_pvalue = std::atof(argv[++i]);
        } else if (key == "--use-fdr-gate") {
            args.use_fdr_gate = true;
        } else if (key == "--apply-bh-correction") {
            args.apply_bh_correction = true;
        } else if (key == "--use-bootstrap") {
            args.use_bootstrap = true;
        } else if (key == "--direct-ed-surface") {
            args.direct_ed_surface = true;
        } else if (key == "--use-legacy-rank-ambient") {
            args.use_legacy_rank_ambient = true;
        } else if (key == "--use-guarded-rank-ambient") {
            args.use_guarded_rank_ambient = true;
        } else if (key == "-h" || key == "--help") {
            usage(argv[0]);
            return 0;
        } else {
            cerr << "Unknown argument: " << key << "\n";
            usage(argv[0]);
            return 2;
        }
    }

    if (args.matrix_path.empty() || args.barcodes_path.empty() || args.out_barcodes.empty()) {
        usage(argv[0]);
        return 2;
    }

    vector<string> barcodes;
    if (!read_barcodes(args.barcodes_path, &barcodes) || barcodes.empty()) {
        cerr << "Failed to read barcodes from " << args.barcodes_path << "\n";
        return 1;
    }

    const bool use_full = (args.mode == "full");

    uint32_t n_rows = 0;
    uint32_t n_cols = 0;
    vector<uint64_t> counts64;
    vector<uint32_t> counts32;
    vector<uint32_t> sparse_gene_ids;
    vector<uint32_t> sparse_counts;
    vector<uint32_t> sparse_cell_index;
    vector<uint32_t> n_genes_per_cell;
    size_t nnz = 0;

    if (!use_full) {
        if (!read_matrix_counts(args.matrix_path, &n_rows, &n_cols, &counts64)) {
            cerr << "Failed to read matrix from " << args.matrix_path << "\n";
            return 1;
        }

        if (n_cols != barcodes.size()) {
            cerr << "Barcode count (" << barcodes.size() << ") does not match matrix columns (" << n_cols << ")\n";
            return 1;
        }
    }

    if (args.direct_ed_surface && !use_full) {
        cerr << "--direct-ed-surface requires --mode full\n";
        return 2;
    }
    if (args.use_legacy_rank_ambient && !use_full) {
        cerr << "--use-legacy-rank-ambient requires --mode full\n";
        return 2;
    }
    if (args.use_guarded_rank_ambient && !use_full) {
        cerr << "--use-guarded-rank-ambient requires --mode full\n";
        return 2;
    }
    if (args.use_legacy_rank_ambient && args.use_guarded_rank_ambient) {
        cerr << "--use-legacy-rank-ambient and --use-guarded-rank-ambient are mutually exclusive\n";
        return 2;
    }

    if (use_full) {
        if (!read_matrix_sparse(args.matrix_path, &n_rows, &n_cols, &counts32,
                                &sparse_gene_ids, &sparse_counts, &sparse_cell_index,
                                &n_genes_per_cell, &nnz)) {
            cerr << "Failed to read sparse matrix for full mode\n";
            return 1;
        }
        if (n_cols != barcodes.size()) {
            cerr << "Barcode count (" << barcodes.size() << ") does not match matrix columns (" << n_cols << ")\n";
            return 1;
        }
        uint64_t total_genes = 0;
        uint64_t total_umis = 0;
        for (uint32_t v : n_genes_per_cell) total_genes += v;
        for (uint32_t v : counts32) total_umis += v;
        cerr << "[scrna_simpleed] full mode: nnz=" << nnz
             << " sum_genes_per_cell=" << total_genes
             << " sum_umis=" << total_umis << "\n";
    } else {
        counts32.assign(n_cols, 0);
        for (size_t i = 0; i < counts64.size(); i++) {
            if (counts64[i] > std::numeric_limits<uint32_t>::max()) {
                counts32[i] = std::numeric_limits<uint32_t>::max();
            } else {
                counts32[i] = static_cast<uint32_t>(counts64[i]);
            }
        }
    }

    if (!args.include_zero_umis) {
        const uint32_t old_cols = n_cols;
        vector<uint32_t> counts_filtered;
        counts_filtered.reserve(n_cols);
        vector<string> barcodes_filtered;
        barcodes_filtered.reserve(barcodes.size());
        vector<uint32_t> old_to_new(old_cols, std::numeric_limits<uint32_t>::max());

        for (uint32_t i = 0; i < old_cols; i++) {
            if (counts32[i] > 0) {
                old_to_new[i] = static_cast<uint32_t>(barcodes_filtered.size());
                barcodes_filtered.push_back(barcodes[i]);
                counts_filtered.push_back(counts32[i]);
            }
        }

        if (counts_filtered.empty()) {
            cerr << "All cells have zero UMIs after filtering; cannot run EmptyDrops\n";
            return 1;
        }

        const uint32_t dropped = old_cols - static_cast<uint32_t>(counts_filtered.size());
        cerr << "[scrna_simpleed] Dropped " << dropped << " zero-UMI cells; kept "
             << counts_filtered.size() << "\n";

        barcodes.swap(barcodes_filtered);
        counts32.swap(counts_filtered);
        n_cols = static_cast<uint32_t>(barcodes.size());

        if (use_full) {
            vector<uint32_t> sparse_gene_ids_f;
            vector<uint32_t> sparse_counts_f;
            vector<uint32_t> sparse_cell_index_f(n_cols + 1, 0);
            vector<uint32_t> n_genes_per_cell_f(n_cols, 0);
            size_t out_pos = 0;

            for (uint32_t old_idx = 0; old_idx < old_cols; old_idx++) {
                uint32_t new_idx = old_to_new[old_idx];
                if (new_idx == std::numeric_limits<uint32_t>::max()) {
                    continue;
                }
                uint32_t n_genes = n_genes_per_cell[old_idx];
                n_genes_per_cell_f[new_idx] = n_genes;
                sparse_cell_index_f[new_idx] = static_cast<uint32_t>(out_pos);
                for (uint32_t k = 0; k < n_genes; k++) {
                    size_t pos = static_cast<size_t>(sparse_cell_index[old_idx] + k);
                    if (pos >= sparse_gene_ids.size()) {
                        continue;
                    }
                    sparse_gene_ids_f.push_back(sparse_gene_ids[pos]);
                    sparse_counts_f.push_back(sparse_counts[pos]);
                    out_pos++;
                }
            }
            sparse_cell_index_f[n_cols] = static_cast<uint32_t>(out_pos);
            nnz = out_pos;
            sparse_gene_ids.swap(sparse_gene_ids_f);
            sparse_counts.swap(sparse_counts_f);
            sparse_cell_index.swap(sparse_cell_index_f);
            n_genes_per_cell.swap(n_genes_per_cell_f);
        }
    }

    vector<char*> barcode_ptrs;
    barcode_ptrs.reserve(barcodes.size());
    for (auto& bc : barcodes) {
        barcode_ptrs.push_back(const_cast<char*>(bc.c_str()));
    }

    scrna_matrix_input input;
    std::memset(&input, 0, sizeof(input));
    input.umi_counts = counts32.data();
    input.barcodes = barcode_ptrs.data();
    input.n_cells = n_cols;
    input.n_features = n_rows;
    input.features = nullptr;
    if (use_full) {
        input.sparse_gene_ids = sparse_gene_ids.data();
        input.sparse_counts = sparse_counts.data();
        input.sparse_cell_index = sparse_cell_index.data();
        input.n_genes_per_cell = n_genes_per_cell.data();
        input.sparse_nnz = nnz;
    } else {
        input.sparse_gene_ids = nullptr;
        input.sparse_counts = nullptr;
        input.sparse_cell_index = nullptr;
        input.n_genes_per_cell = nullptr;
        input.sparse_nnz = 0;
    }

    scrna_ed_config* config = scrna_ed_config_create();
    if (!config) {
        cerr << "Failed to create scrna_ed_config\n";
        return 1;
    }
    if (args.expected_cells > 0) {
        config->n_expected_cells = args.expected_cells;
    }
    if (args.umi_min > 0) {
        config->umi_min = args.umi_min;
    }
    if (args.sim_n > 0) {
        config->sim_n = args.sim_n;
    }
    if (args.ed_retain_count > 0) {
        config->ed_retain_count = args.ed_retain_count;
    }
    if (args.lower_testing_bound > 0) {
        config->lower_testing_bound = args.lower_testing_bound;
    }
    if (args.ambient_umi_max > 0) {
        config->ambient_umi_max = args.ambient_umi_max;
    }
    if (args.fdr > 0.0) {
        config->fdr = args.fdr;
    }
    if (args.raw_pvalue > 0.0) {
        config->raw_pvalue_threshold = args.raw_pvalue;
    }
    if (args.use_fdr_gate) {
        config->use_fdr_gate = 1;
    }
    if (args.apply_bh_correction) {
        config->apply_bh_correction = 1;
    }
    if (args.use_bootstrap) {
        config->use_bootstrap = 1;
    }
    config->disable_occupancy_filter = 1;

    scrna_ed_result result;
    std::memset(&result, 0, sizeof(result));
    int rc = 0;
    if (args.direct_ed_surface) {
        rc = run_direct_ed_surface(barcodes, counts32, sparse_gene_ids, sparse_counts,
                                   sparse_cell_index, n_genes_per_cell, n_rows, config, &result);
    } else if (args.use_legacy_rank_ambient || args.use_guarded_rank_ambient) {
        rc = run_simpleed_custom_ambient(barcodes, counts32, sparse_gene_ids, sparse_counts,
                                         sparse_cell_index, n_genes_per_cell, n_rows, config,
                                         args.use_legacy_rank_ambient,
                                         args.use_guarded_rank_ambient,
                                         args.ambient_fallback_min_abs,
                                         args.ambient_fallback_min_frac,
                                         &result);
    } else {
        rc = scrna_emptydrops_run(&input, config, &result);
    }
    if (rc != 0) {
        cerr << "scrna_emptydrops_run failed\n";
        if (result.error_message) {
            cerr << "  error: " << result.error_message << "\n";
        }
        scrna_ed_result_free(&result);
        scrna_ed_config_destroy(config);
        return 1;
    }

    if (scrna_write_filtered_barcodes(&result, args.out_barcodes.c_str()) != 0) {
        cerr << "Failed to write barcodes to " << args.out_barcodes << "\n";
        scrna_ed_result_free(&result);
        scrna_ed_config_destroy(config);
        return 1;
    }

    if (!args.out_dir.empty()) {
        if (scrna_emptydrops_write_outputs(&result, args.out_dir.c_str()) != 0) {
            cerr << "Warning: failed to write detailed outputs to " << args.out_dir << "\n";
        }
    }

    cout << "SimpleED completed: " << result.n_barcodes << " barcodes written\n";

    scrna_ed_result_free(&result);
    scrna_ed_config_destroy(config);
    return 0;
}
