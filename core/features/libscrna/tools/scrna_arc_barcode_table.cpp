#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::unordered_map;
using std::vector;

namespace {

struct Args {
    string per_barcode_metrics_path;
    string out_path;
};

void usage(const char* prog) {
    cerr << "Usage: " << prog
         << " --per-barcode-metrics outs/per_barcode_metrics.csv --out arc_barcode_table.tsv\n";
}

void trim_cr(string* value) {
    if (value != nullptr && !value->empty() && value->back() == '\r') {
        value->pop_back();
    }
}

vector<string> parse_csv_line(const string& line) {
    vector<string> fields;
    string current;
    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (in_quotes) {
            if (c == '"') {
                if (i + 1 < line.size() && line[i + 1] == '"') {
                    current.push_back('"');
                    ++i;
                } else {
                    in_quotes = false;
                }
            } else {
                current.push_back(c);
            }
        } else if (c == '"') {
            in_quotes = true;
        } else if (c == ',') {
            fields.push_back(current);
            current.clear();
        } else {
            current.push_back(c);
        }
    }

    fields.push_back(current);
    return fields;
}

bool parse_args(int argc, char** argv, Args* args) {
    if (args == nullptr) {
        return false;
    }

    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];
        if (arg == "--per-barcode-metrics" && i + 1 < argc) {
            args->per_barcode_metrics_path = argv[++i];
        } else if (arg == "--out" && i + 1 < argc) {
            args->out_path = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return false;
        } else {
            cerr << "Unknown or incomplete argument: " << arg << endl;
            usage(argv[0]);
            return false;
        }
    }

    if (args->per_barcode_metrics_path.empty() || args->out_path.empty()) {
        usage(argv[0]);
        return false;
    }

    return true;
}

bool require_column(const unordered_map<string, size_t>& header_index,
                    const string& column_name,
                    size_t* index_out) {
    if (index_out == nullptr) {
        return false;
    }

    const auto it = header_index.find(column_name);
    if (it == header_index.end()) {
        cerr << "Missing required column: " << column_name << endl;
        return false;
    }

    *index_out = it->second;
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, &args)) {
        return 1;
    }

    ifstream input(args.per_barcode_metrics_path.c_str());
    if (!input.is_open()) {
        cerr << "Failed to open input: " << args.per_barcode_metrics_path << endl;
        return 1;
    }

    ofstream output(args.out_path.c_str());
    if (!output.is_open()) {
        cerr << "Failed to open output: " << args.out_path << endl;
        return 1;
    }

    string line;
    if (!std::getline(input, line)) {
        cerr << "Input file is empty: " << args.per_barcode_metrics_path << endl;
        return 1;
    }
    trim_cr(&line);

    const vector<string> header = parse_csv_line(line);
    unordered_map<string, size_t> header_index;
    for (size_t i = 0; i < header.size(); ++i) {
        header_index[header[i]] = i;
    }

    size_t barcode_idx = 0;
    size_t gex_barcode_idx = 0;
    size_t atac_barcode_idx = 0;
    size_t gex_umis_count_idx = 0;
    size_t atac_peak_region_cutsites_idx = 0;
    size_t atac_peak_region_fragments_idx = 0;
    size_t atac_fragments_idx = 0;
    size_t excluded_reason_idx = 0;
    size_t is_cell_idx = 0;

    if (!require_column(header_index, "barcode", &barcode_idx) ||
        !require_column(header_index, "gex_barcode", &gex_barcode_idx) ||
        !require_column(header_index, "atac_barcode", &atac_barcode_idx) ||
        !require_column(header_index, "gex_umis_count", &gex_umis_count_idx) ||
        !require_column(header_index, "atac_peak_region_cutsites", &atac_peak_region_cutsites_idx) ||
        !require_column(header_index, "atac_peak_region_fragments", &atac_peak_region_fragments_idx) ||
        !require_column(header_index, "atac_fragments", &atac_fragments_idx) ||
        !require_column(header_index, "excluded_reason", &excluded_reason_idx) ||
        !require_column(header_index, "is_cell", &is_cell_idx)) {
        return 1;
    }

    output << "barcode\tgex_barcode\tatac_barcode\tgex_umis_count\tatac_peak_region_cutsites\t"
              "atac_peak_region_fragments\tatac_fragments\texcluded_reason\tis_cell\n";

    uint64_t row_count = 0;
    uint64_t cell_count = 0;
    uint64_t line_number = 1;

    while (std::getline(input, line)) {
        ++line_number;
        trim_cr(&line);
        if (line.empty()) {
            continue;
        }

        const vector<string> fields = parse_csv_line(line);
        const size_t max_required_index = std::max(
            std::max(std::max(barcode_idx, gex_barcode_idx), std::max(atac_barcode_idx, gex_umis_count_idx)),
            std::max(std::max(atac_peak_region_cutsites_idx, atac_peak_region_fragments_idx),
                     std::max(atac_fragments_idx, std::max(excluded_reason_idx, is_cell_idx))));

        if (fields.size() <= max_required_index) {
            cerr << "Line " << line_number << " has too few columns: expected at least "
                 << (max_required_index + 1) << ", observed " << fields.size() << endl;
            return 1;
        }

        output << fields[barcode_idx] << '\t'
               << fields[gex_barcode_idx] << '\t'
               << fields[atac_barcode_idx] << '\t'
               << fields[gex_umis_count_idx] << '\t'
               << fields[atac_peak_region_cutsites_idx] << '\t'
               << fields[atac_peak_region_fragments_idx] << '\t'
               << fields[atac_fragments_idx] << '\t'
               << fields[excluded_reason_idx] << '\t'
               << fields[is_cell_idx] << '\n';

        ++row_count;
        if (fields[is_cell_idx] == "1") {
            ++cell_count;
        }
    }

    cerr << "Wrote " << row_count << " barcode rows to " << args.out_path << endl;
    cerr << "Rows with is_cell=1: " << cell_count << endl;
    return 0;
}
