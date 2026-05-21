#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

namespace {

struct Args {
    string arc_table_path;
    string filtered_barcodes_path;
    string emptydrops_results_path;
    string out_path;

    string arc_barcode_col = "barcode";
    string arc_umis_col = "gex_umis_count";
    string barcode_suffix;
    string rescue_mode = "candidate_non_call";
};

struct ArcRow {
    string barcode;
    string gex_umis_count;
};

struct EdEntry {
    int in_results = 0;
    int passes_raw_p = 0;
    int passes_fdr = 0;
    int is_simple_cell = 0;
};

void usage(const char* prog) {
    cerr << "Usage: " << prog
         << " --arc-table arc_barcode_table.tsv"
         << " --filtered-barcodes filtered_barcodes.txt"
         << " [--emptydrops-results EmptyDrops/emptydrops_results.tsv]"
         << " --out gex_evidence.tsv"
         << " [--arc-barcode-col barcode] [--arc-umis-col gex_umis_count]"
         << " [--barcode-suffix -1]"
         << " [--rescue-mode candidate_non_call|raw_p_non_call|fdr_non_call|none]\n";
}

void trim_cr(string* value) {
    if (value != nullptr && !value->empty() && value->back() == '\r') {
        value->pop_back();
    }
}

vector<string> parse_tsv_line(const string& line) {
    vector<string> fields;
    string current;
    for (size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '\t') {
            fields.push_back(current);
            current.clear();
        } else {
            current.push_back(line[i]);
        }
    }
    fields.push_back(current);
    return fields;
}

string normalize_barcode(const string& barcode, const string& suffix) {
    if (suffix.empty()) {
        return barcode;
    }
    if (barcode.size() >= suffix.size() &&
        barcode.compare(barcode.size() - suffix.size(), suffix.size(), suffix) == 0) {
        return barcode;
    }
    return barcode + suffix;
}

bool parse_binary_call(const string& raw_value, int* out) {
    if (out == nullptr) {
        return false;
    }
    if (raw_value == "1" || raw_value == "true" || raw_value == "TRUE" || raw_value == "yes") {
        *out = 1;
        return true;
    }
    if (raw_value == "0" || raw_value == "false" || raw_value == "FALSE" || raw_value == "no") {
        *out = 0;
        return true;
    }
    return false;
}

bool parse_args(int argc, char** argv, Args* args) {
    if (args == nullptr) {
        return false;
    }

    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];
        if (arg == "--arc-table" && i + 1 < argc) {
            args->arc_table_path = argv[++i];
        } else if (arg == "--filtered-barcodes" && i + 1 < argc) {
            args->filtered_barcodes_path = argv[++i];
        } else if (arg == "--emptydrops-results" && i + 1 < argc) {
            args->emptydrops_results_path = argv[++i];
        } else if (arg == "--out" && i + 1 < argc) {
            args->out_path = argv[++i];
        } else if (arg == "--arc-barcode-col" && i + 1 < argc) {
            args->arc_barcode_col = argv[++i];
        } else if (arg == "--arc-umis-col" && i + 1 < argc) {
            args->arc_umis_col = argv[++i];
        } else if (arg == "--barcode-suffix" && i + 1 < argc) {
            args->barcode_suffix = argv[++i];
        } else if (arg == "--rescue-mode" && i + 1 < argc) {
            args->rescue_mode = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return false;
        } else {
            cerr << "Unknown or incomplete argument: " << arg << endl;
            usage(argv[0]);
            return false;
        }
    }

    if (args->arc_table_path.empty() || args->filtered_barcodes_path.empty() || args->out_path.empty()) {
        usage(argv[0]);
        return false;
    }

    if (args->rescue_mode != "candidate_non_call" &&
        args->rescue_mode != "raw_p_non_call" &&
        args->rescue_mode != "fdr_non_call" &&
        args->rescue_mode != "none") {
        cerr << "Unsupported rescue mode: " << args->rescue_mode << endl;
        usage(argv[0]);
        return false;
    }

    return true;
}

bool find_required_column(const unordered_map<string, size_t>& header_index,
                          const string& name,
                          size_t* index_out) {
    if (index_out == nullptr) {
        return false;
    }
    const auto it = header_index.find(name);
    if (it == header_index.end()) {
        cerr << "Missing required column: " << name << endl;
        return false;
    }
    *index_out = it->second;
    return true;
}

bool load_arc_rows(const string& path,
                   const string& barcode_col,
                   const string& umis_col,
                   vector<ArcRow>* rows) {
    if (rows == nullptr) {
        return false;
    }

    ifstream input(path.c_str());
    if (!input.is_open()) {
        cerr << "Failed to open ARC table: " << path << endl;
        return false;
    }

    string line;
    if (!std::getline(input, line)) {
        cerr << "ARC table is empty: " << path << endl;
        return false;
    }
    trim_cr(&line);

    const vector<string> header = parse_tsv_line(line);
    unordered_map<string, size_t> header_index;
    for (size_t i = 0; i < header.size(); ++i) {
        header_index[header[i]] = i;
    }

    size_t barcode_idx = 0;
    size_t umis_idx = 0;
    if (!find_required_column(header_index, barcode_col, &barcode_idx) ||
        !find_required_column(header_index, umis_col, &umis_idx)) {
        return false;
    }

    unordered_set<string> seen;
    uint64_t line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        trim_cr(&line);
        if (line.empty()) {
            continue;
        }

        const vector<string> fields = parse_tsv_line(line);
        const size_t max_required_idx = barcode_idx > umis_idx ? barcode_idx : umis_idx;
        if (fields.size() <= max_required_idx) {
            cerr << "ARC table line " << line_number << " has too few columns" << endl;
            return false;
        }

        const string& barcode = fields[barcode_idx];
        if (seen.count(barcode) != 0) {
            cerr << "Duplicate barcode in ARC table: " << barcode << endl;
            return false;
        }
        seen.insert(barcode);
        rows->push_back(ArcRow{barcode, fields[umis_idx]});
    }

    return true;
}

bool load_filtered_barcodes(const string& path,
                            const string& suffix,
                            unordered_set<string>* filtered) {
    if (filtered == nullptr) {
        return false;
    }

    ifstream input(path.c_str());
    if (!input.is_open()) {
        cerr << "Failed to open filtered_barcodes.txt: " << path << endl;
        return false;
    }

    string line;
    while (std::getline(input, line)) {
        trim_cr(&line);
        if (line.empty()) {
            continue;
        }
        const string barcode = normalize_barcode(line, suffix);
        if (filtered->count(barcode) != 0) {
            cerr << "Duplicate barcode in filtered_barcodes.txt: " << barcode << endl;
            return false;
        }
        filtered->insert(barcode);
    }

    return true;
}

bool load_emptydrops_results(const string& path,
                             const string& suffix,
                             unordered_map<string, EdEntry>* by_barcode) {
    if (by_barcode == nullptr) {
        return false;
    }

    ifstream input(path.c_str());
    if (!input.is_open()) {
        cerr << "Failed to open emptydrops_results.tsv: " << path << endl;
        return false;
    }

    string line;
    if (!std::getline(input, line)) {
        cerr << "emptydrops_results.tsv is empty: " << path << endl;
        return false;
    }
    trim_cr(&line);

    const vector<string> header = parse_tsv_line(line);
    unordered_map<string, size_t> header_index;
    for (size_t i = 0; i < header.size(); ++i) {
        header_index[header[i]] = i;
    }

    size_t barcode_idx = 0;
    size_t passes_raw_idx = 0;
    size_t passes_fdr_idx = 0;
    size_t is_simple_idx = 0;
    if (!find_required_column(header_index, "barcode", &barcode_idx) ||
        !find_required_column(header_index, "passes_raw_p", &passes_raw_idx) ||
        !find_required_column(header_index, "passes_fdr", &passes_fdr_idx) ||
        !find_required_column(header_index, "is_simple_cell", &is_simple_idx)) {
        return false;
    }

    uint64_t line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        trim_cr(&line);
        if (line.empty()) {
            continue;
        }

        const vector<string> fields = parse_tsv_line(line);
        size_t max_required_idx = barcode_idx;
        if (passes_raw_idx > max_required_idx) max_required_idx = passes_raw_idx;
        if (passes_fdr_idx > max_required_idx) max_required_idx = passes_fdr_idx;
        if (is_simple_idx > max_required_idx) max_required_idx = is_simple_idx;
        if (fields.size() <= max_required_idx) {
            cerr << "emptydrops_results.tsv line " << line_number << " has too few columns" << endl;
            return false;
        }

        EdEntry entry;
        entry.in_results = 1;
        if (!parse_binary_call(fields[passes_raw_idx], &entry.passes_raw_p) ||
            !parse_binary_call(fields[passes_fdr_idx], &entry.passes_fdr) ||
            !parse_binary_call(fields[is_simple_idx], &entry.is_simple_cell)) {
            cerr << "emptydrops_results.tsv line " << line_number << " has invalid binary fields" << endl;
            return false;
        }

        const string barcode = normalize_barcode(fields[barcode_idx], suffix);
        if (by_barcode->count(barcode) != 0) {
            cerr << "Duplicate barcode in emptydrops_results.tsv: " << barcode << endl;
            return false;
        }
        (*by_barcode)[barcode] = entry;
    }

    return true;
}

int derive_rescue_eligible(const string& mode,
                           int in_filtered,
                           const EdEntry* ed_entry) {
    if (mode == "none" || ed_entry == nullptr || !ed_entry->in_results || in_filtered) {
        return 0;
    }
    if (mode == "candidate_non_call") {
        return 1;
    }
    if (mode == "raw_p_non_call") {
        return ed_entry->passes_raw_p;
    }
    if (mode == "fdr_non_call") {
        return ed_entry->passes_fdr;
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, &args)) {
        return 1;
    }

    vector<ArcRow> arc_rows;
    if (!load_arc_rows(args.arc_table_path, args.arc_barcode_col, args.arc_umis_col, &arc_rows)) {
        return 1;
    }

    unordered_set<string> filtered_barcodes;
    if (!load_filtered_barcodes(args.filtered_barcodes_path, args.barcode_suffix, &filtered_barcodes)) {
        return 1;
    }

    unordered_map<string, EdEntry> ed_by_barcode;
    if (!args.emptydrops_results_path.empty() &&
        !load_emptydrops_results(args.emptydrops_results_path, args.barcode_suffix, &ed_by_barcode)) {
        return 1;
    }

    ofstream output(args.out_path.c_str());
    if (!output.is_open()) {
        cerr << "Failed to open output: " << args.out_path << endl;
        return 1;
    }

    const string source_value = args.emptydrops_results_path.empty()
        ? "filtered_barcodes_only"
        : "filtered_barcodes_plus_emptydrops_results";

    output << "barcode\tgex_umis_count\tgex_in_filtered_barcodes\tgex_in_emptydrops_results\t"
              "gex_passes_raw_p\tgex_passes_fdr\tgex_is_simple_cell\tgex_module_call\t"
              "gex_rescue_eligible\tgex_source\n";

    uint64_t rows_written = 0;
    uint64_t module_call_count = 0;
    uint64_t rescue_eligible_count = 0;
    uint64_t in_ed_results_count = 0;

    for (size_t i = 0; i < arc_rows.size(); ++i) {
        const ArcRow& row = arc_rows[i];
        const int in_filtered = filtered_barcodes.count(row.barcode) ? 1 : 0;
        const auto ed_it = ed_by_barcode.find(row.barcode);
        const EdEntry* ed_entry = (ed_it != ed_by_barcode.end()) ? &ed_it->second : nullptr;

        const int in_results = (ed_entry != nullptr) ? ed_entry->in_results : 0;
        const int passes_raw_p = (ed_entry != nullptr) ? ed_entry->passes_raw_p : 0;
        const int passes_fdr = (ed_entry != nullptr) ? ed_entry->passes_fdr : 0;
        const int is_simple_cell = (ed_entry != nullptr) ? ed_entry->is_simple_cell : 0;
        const int gex_module_call = in_filtered;
        const int gex_rescue_eligible = derive_rescue_eligible(args.rescue_mode, in_filtered, ed_entry);

        output << row.barcode << '\t'
               << row.gex_umis_count << '\t'
               << in_filtered << '\t'
               << in_results << '\t'
               << passes_raw_p << '\t'
               << passes_fdr << '\t'
               << is_simple_cell << '\t'
               << gex_module_call << '\t'
               << gex_rescue_eligible << '\t'
               << source_value << '\n';

        ++rows_written;
        if (gex_module_call) {
            ++module_call_count;
        }
        if (gex_rescue_eligible) {
            ++rescue_eligible_count;
        }
        if (in_results) {
            ++in_ed_results_count;
        }
    }

    cerr << "Wrote " << rows_written << " GEX evidence rows to " << args.out_path << endl;
    cerr << "gex_module_call=1 rows: " << module_call_count << endl;
    cerr << "gex_rescue_eligible=1 rows: " << rescue_eligible_count << endl;
    cerr << "Rows present in emptydrops_results.tsv: " << in_ed_results_count << endl;
    return 0;
}
