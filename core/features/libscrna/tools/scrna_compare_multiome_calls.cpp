#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
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
    string calls_path;
    string out_path;
    string summary_path;
    string arc_call_col = "is_cell";
    string call_col = "final_is_cell";
    string barcode_col = "barcode";
};

struct ArcEntry {
    string barcode;
    int arc_is_cell = 0;

    ArcEntry() {}

    ArcEntry(const string& barcode_in, int arc_is_cell_in)
        : barcode(barcode_in), arc_is_cell(arc_is_cell_in) {}
};

struct ExternalEntry {
    int final_is_cell = 0;

    ExternalEntry() {}

    explicit ExternalEntry(int final_is_cell_in)
        : final_is_cell(final_is_cell_in) {}
};

void usage(const char* prog) {
    cerr << "Usage: " << prog
         << " --arc-table arc_barcode_table.tsv --calls multiome_calls.tsv"
         << " --out comparison.tsv --summary summary.tsv"
         << " [--arc-call-col is_cell] [--call-col final_is_cell] [--barcode-col barcode]\n";
}

void trim_cr(string* value) {
    if (value != nullptr && !value->empty() && value->back() == '\r') {
        value->pop_back();
    }
}

string to_lower_ascii(string value) {
    for (size_t i = 0; i < value.size(); ++i) {
        value[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(value[i])));
    }
    return value;
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

bool parse_args(int argc, char** argv, Args* args) {
    if (args == nullptr) {
        return false;
    }

    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];
        if (arg == "--arc-table" && i + 1 < argc) {
            args->arc_table_path = argv[++i];
        } else if (arg == "--calls" && i + 1 < argc) {
            args->calls_path = argv[++i];
        } else if (arg == "--out" && i + 1 < argc) {
            args->out_path = argv[++i];
        } else if (arg == "--summary" && i + 1 < argc) {
            args->summary_path = argv[++i];
        } else if (arg == "--arc-call-col" && i + 1 < argc) {
            args->arc_call_col = argv[++i];
        } else if (arg == "--call-col" && i + 1 < argc) {
            args->call_col = argv[++i];
        } else if (arg == "--barcode-col" && i + 1 < argc) {
            args->barcode_col = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return false;
        } else {
            cerr << "Unknown or incomplete argument: " << arg << endl;
            usage(argv[0]);
            return false;
        }
    }

    if (args->arc_table_path.empty() || args->calls_path.empty() ||
        args->out_path.empty() || args->summary_path.empty()) {
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

bool parse_binary_call(const string& raw_value, int* out) {
    if (out == nullptr) {
        return false;
    }

    const string value = to_lower_ascii(raw_value);
    if (value == "1" || value == "true" || value == "yes") {
        *out = 1;
        return true;
    }
    if (value == "0" || value == "false" || value == "no") {
        *out = 0;
        return true;
    }
    return false;
}

bool load_arc_table(const string& path,
                    const string& barcode_col,
                    const string& call_col,
                    vector<ArcEntry>* entries,
                    unordered_map<string, int>* by_barcode) {
    if (entries == nullptr || by_barcode == nullptr) {
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
    size_t call_idx = 0;
    if (!find_required_column(header_index, barcode_col, &barcode_idx) ||
        !find_required_column(header_index, call_col, &call_idx)) {
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
        const size_t max_required = (barcode_idx > call_idx) ? barcode_idx : call_idx;
        if (fields.size() <= max_required) {
            cerr << "ARC table line " << line_number << " has too few columns" << endl;
            return false;
        }

        int call_value = 0;
        if (!parse_binary_call(fields[call_idx], &call_value)) {
            cerr << "ARC table line " << line_number << " has invalid call value: "
                 << fields[call_idx] << endl;
            return false;
        }

        const string& barcode = fields[barcode_idx];
        if (by_barcode->count(barcode) != 0) {
            cerr << "Duplicate barcode in ARC table: " << barcode << endl;
            return false;
        }

        entries->push_back(ArcEntry{barcode, call_value});
        (*by_barcode)[barcode] = call_value;
    }

    return true;
}

bool load_external_calls(const string& path,
                         const string& barcode_col,
                         const string& call_col,
                         unordered_map<string, ExternalEntry>* by_barcode) {
    if (by_barcode == nullptr) {
        return false;
    }

    ifstream input(path.c_str());
    if (!input.is_open()) {
        cerr << "Failed to open external calls table: " << path << endl;
        return false;
    }

    string line;
    if (!std::getline(input, line)) {
        cerr << "External calls table is empty: " << path << endl;
        return false;
    }
    trim_cr(&line);

    const vector<string> header = parse_tsv_line(line);
    unordered_map<string, size_t> header_index;
    for (size_t i = 0; i < header.size(); ++i) {
        header_index[header[i]] = i;
    }

    size_t barcode_idx = 0;
    size_t call_idx = 0;
    if (!find_required_column(header_index, barcode_col, &barcode_idx) ||
        !find_required_column(header_index, call_col, &call_idx)) {
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
        const size_t max_required = (barcode_idx > call_idx) ? barcode_idx : call_idx;
        if (fields.size() <= max_required) {
            cerr << "External calls line " << line_number << " has too few columns" << endl;
            return false;
        }

        int call_value = 0;
        if (!parse_binary_call(fields[call_idx], &call_value)) {
            cerr << "External calls line " << line_number << " has invalid call value: "
                 << fields[call_idx] << endl;
            return false;
        }

        const string& barcode = fields[barcode_idx];
        if (by_barcode->count(barcode) != 0) {
            cerr << "Duplicate barcode in external calls table: " << barcode << endl;
            return false;
        }

        (*by_barcode)[barcode] = ExternalEntry{call_value};
    }

    return true;
}

string format_rate(uint64_t numerator, uint64_t denominator) {
    if (denominator == 0) {
        return "NA";
    }

    std::ostringstream out;
    out << std::fixed << std::setprecision(6)
        << (static_cast<double>(numerator) / static_cast<double>(denominator));
    return out.str();
}

}  // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, &args)) {
        return 1;
    }

    vector<ArcEntry> arc_entries;
    unordered_map<string, int> arc_calls_by_barcode;
    if (!load_arc_table(args.arc_table_path, args.barcode_col, args.arc_call_col,
                        &arc_entries, &arc_calls_by_barcode)) {
        return 1;
    }

    unordered_map<string, ExternalEntry> external_calls_by_barcode;
    if (!load_external_calls(args.calls_path, args.barcode_col, args.call_col,
                             &external_calls_by_barcode)) {
        return 1;
    }

    ofstream out(args.out_path.c_str());
    if (!out.is_open()) {
        cerr << "Failed to open comparison output: " << args.out_path << endl;
        return 1;
    }

    ofstream summary(args.summary_path.c_str());
    if (!summary.is_open()) {
        cerr << "Failed to open summary output: " << args.summary_path << endl;
        return 1;
    }

    out << "barcode\tarc_is_cell\tfinal_is_cell\thas_external_call\tagreement\tcomparison_category\n";

    uint64_t arc_rows = 0;
    uint64_t arc_positive = 0;
    uint64_t arc_negative = 0;
    uint64_t external_rows_total = external_calls_by_barcode.size();
    uint64_t external_rows_matched = 0;
    uint64_t external_rows_unmatched = 0;
    uint64_t external_positive_matched = 0;
    uint64_t external_negative_matched = 0;
    uint64_t rows_missing_external_call = 0;
    uint64_t agreement_rows = 0;
    uint64_t disagreement_rows = 0;
    uint64_t true_positive = 0;
    uint64_t true_negative = 0;
    uint64_t false_positive = 0;
    uint64_t false_negative = 0;

    unordered_set<string> matched_barcodes;
    matched_barcodes.reserve(arc_entries.size());

    for (size_t i = 0; i < arc_entries.size(); ++i) {
        const ArcEntry& entry = arc_entries[i];
        ++arc_rows;
        if (entry.arc_is_cell == 1) {
            ++arc_positive;
        } else {
            ++arc_negative;
        }

        const auto external_it = external_calls_by_barcode.find(entry.barcode);
        if (external_it == external_calls_by_barcode.end()) {
            ++rows_missing_external_call;
            out << entry.barcode << '\t'
                << entry.arc_is_cell << '\t'
                << "" << '\t'
                << 0 << '\t'
                << "" << '\t'
                << "missing_external_call\n";
            continue;
        }

        matched_barcodes.insert(entry.barcode);
        ++external_rows_matched;

        const int external_call = external_it->second.final_is_cell;
        if (external_call == 1) {
            ++external_positive_matched;
        } else {
            ++external_negative_matched;
        }

        if (entry.arc_is_cell == external_call) {
            ++agreement_rows;
            if (entry.arc_is_cell == 1) {
                ++true_positive;
                out << entry.barcode << '\t' << 1 << '\t' << 1 << '\t' << 1 << '\t' << 1
                    << '\t' << "true_positive\n";
            } else {
                ++true_negative;
                out << entry.barcode << '\t' << 0 << '\t' << 0 << '\t' << 1 << '\t' << 1
                    << '\t' << "true_negative\n";
            }
        } else {
            ++disagreement_rows;
            if (entry.arc_is_cell == 0 && external_call == 1) {
                ++false_positive;
                out << entry.barcode << '\t' << 0 << '\t' << 1 << '\t' << 1 << '\t' << 0
                    << '\t' << "false_positive\n";
            } else {
                ++false_negative;
                out << entry.barcode << '\t' << 1 << '\t' << 0 << '\t' << 1 << '\t' << 0
                    << '\t' << "false_negative\n";
            }
        }
    }

    for (unordered_map<string, ExternalEntry>::const_iterator it = external_calls_by_barcode.begin();
         it != external_calls_by_barcode.end(); ++it) {
        if (matched_barcodes.count(it->first) == 0) {
            ++external_rows_unmatched;
        }
    }

    summary << "metric\tvalue\n";
    summary << "arc_rows\t" << arc_rows << '\n';
    summary << "arc_positive\t" << arc_positive << '\n';
    summary << "arc_negative\t" << arc_negative << '\n';
    summary << "external_rows_total\t" << external_rows_total << '\n';
    summary << "external_rows_matched\t" << external_rows_matched << '\n';
    summary << "external_rows_unmatched\t" << external_rows_unmatched << '\n';
    summary << "external_positive_matched\t" << external_positive_matched << '\n';
    summary << "external_negative_matched\t" << external_negative_matched << '\n';
    summary << "rows_missing_external_call\t" << rows_missing_external_call << '\n';
    summary << "agreement_rows\t" << agreement_rows << '\n';
    summary << "disagreement_rows\t" << disagreement_rows << '\n';
    summary << "true_positive\t" << true_positive << '\n';
    summary << "true_negative\t" << true_negative << '\n';
    summary << "false_positive\t" << false_positive << '\n';
    summary << "false_negative\t" << false_negative << '\n';
    summary << "accuracy\t" << format_rate(agreement_rows, external_rows_matched) << '\n';
    summary << "sensitivity\t" << format_rate(true_positive, true_positive + false_negative) << '\n';
    summary << "specificity\t" << format_rate(true_negative, true_negative + false_positive) << '\n';
    summary << "precision\t" << format_rate(true_positive, true_positive + false_positive) << '\n';

    cerr << "Compared " << external_rows_matched << " matched barcode rows against ARC" << endl;
    cerr << "Missing external calls for " << rows_missing_external_call << " ARC barcodes" << endl;
    cerr << "Wrote comparison rows to " << args.out_path << endl;
    cerr << "Wrote summary metrics to " << args.summary_path << endl;
    return 0;
}
