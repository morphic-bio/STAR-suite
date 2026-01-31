#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

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
    double fdr = 0.0;
    double raw_pvalue = 0.0;
    bool use_fdr_gate = false;
};

void usage(const char* prog) {
    cerr << "Usage: " << prog << " --matrix matrix.mtx --barcodes barcodes.tsv --out barcodes.tsv"
         << " [--mode simple|full] [--out-dir DIR] [--expected-cells N] [--umi-min N]"
         << " [--sim-n N] [--ed-retain-count N] [--fdr X] [--raw-pvalue X] [--use-fdr-gate]"
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
        } else if (key == "--fdr" && i + 1 < argc) {
            args.fdr = std::atof(argv[++i]);
        } else if (key == "--raw-pvalue" && i + 1 < argc) {
            args.raw_pvalue = std::atof(argv[++i]);
        } else if (key == "--use-fdr-gate") {
            args.use_fdr_gate = true;
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
    if (args.fdr > 0.0) {
        config->fdr = args.fdr;
    }
    if (args.raw_pvalue > 0.0) {
        config->raw_pvalue_threshold = args.raw_pvalue;
    }
    if (args.use_fdr_gate) {
        config->use_fdr_gate = 1;
    }
    config->disable_occupancy_filter = 1;

    scrna_ed_result result;
    int rc = scrna_emptydrops_run(&input, config, &result);
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
