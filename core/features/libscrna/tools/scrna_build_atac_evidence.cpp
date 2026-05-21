#include <cstdint>
#include <cstdlib>
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
    string out_path;

    string barcode_col = "barcode";
    string peak_cutsites_col = "atac_peak_region_cutsites";
    string peak_fragments_col = "atac_peak_region_fragments";
    string atac_fragments_col = "atac_fragments";
    string excluded_reason_col = "excluded_reason";

    uint64_t min_peak_region_cutsites = 1;
    uint64_t min_peak_region_fragments = 0;
    uint64_t min_atac_fragments = 0;
    double min_peak_fraction = 0.0;
};

struct ArcAtacRow {
    string barcode;
    string peak_cutsites_raw;
    string peak_fragments_raw;
    string atac_fragments_raw;
    string excluded_reason;
    uint64_t peak_cutsites = 0;
    uint64_t peak_fragments = 0;
    uint64_t atac_fragments = 0;
};

void usage(const char* prog) {
    cerr << "Usage: " << prog
         << " --arc-table arc_barcode_table.tsv --out atac_evidence.tsv"
         << " [--barcode-col barcode]"
         << " [--peak-cutsites-col atac_peak_region_cutsites]"
         << " [--peak-fragments-col atac_peak_region_fragments]"
         << " [--atac-fragments-col atac_fragments]"
         << " [--excluded-reason-col excluded_reason]"
         << " [--min-peak-region-cutsites N]"
         << " [--min-peak-region-fragments N]"
         << " [--min-atac-fragments N]"
         << " [--min-peak-fraction X]\n";
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

bool parse_uint64(const string& value, uint64_t* out) {
    if (out == nullptr) {
        return false;
    }
    char* end = nullptr;
    const unsigned long long parsed = std::strtoull(value.c_str(), &end, 10);
    if (end == value.c_str() || *end != '\0') {
        return false;
    }
    *out = static_cast<uint64_t>(parsed);
    return true;
}

bool parse_double(const string& value, double* out) {
    if (out == nullptr) {
        return false;
    }
    char* end = nullptr;
    const double parsed = std::strtod(value.c_str(), &end);
    if (end == value.c_str() || *end != '\0') {
        return false;
    }
    *out = parsed;
    return true;
}

bool parse_args(int argc, char** argv, Args* args) {
    if (args == nullptr) {
        return false;
    }

    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];
        if (arg == "--arc-table" && i + 1 < argc) {
            args->arc_table_path = argv[++i];
        } else if (arg == "--out" && i + 1 < argc) {
            args->out_path = argv[++i];
        } else if (arg == "--barcode-col" && i + 1 < argc) {
            args->barcode_col = argv[++i];
        } else if (arg == "--peak-cutsites-col" && i + 1 < argc) {
            args->peak_cutsites_col = argv[++i];
        } else if (arg == "--peak-fragments-col" && i + 1 < argc) {
            args->peak_fragments_col = argv[++i];
        } else if (arg == "--atac-fragments-col" && i + 1 < argc) {
            args->atac_fragments_col = argv[++i];
        } else if (arg == "--excluded-reason-col" && i + 1 < argc) {
            args->excluded_reason_col = argv[++i];
        } else if (arg == "--min-peak-region-cutsites" && i + 1 < argc) {
            if (!parse_uint64(argv[++i], &args->min_peak_region_cutsites)) {
                cerr << "Invalid --min-peak-region-cutsites value\n";
                return false;
            }
        } else if (arg == "--min-peak-region-fragments" && i + 1 < argc) {
            if (!parse_uint64(argv[++i], &args->min_peak_region_fragments)) {
                cerr << "Invalid --min-peak-region-fragments value\n";
                return false;
            }
        } else if (arg == "--min-atac-fragments" && i + 1 < argc) {
            if (!parse_uint64(argv[++i], &args->min_atac_fragments)) {
                cerr << "Invalid --min-atac-fragments value\n";
                return false;
            }
        } else if (arg == "--min-peak-fraction" && i + 1 < argc) {
            if (!parse_double(argv[++i], &args->min_peak_fraction)) {
                cerr << "Invalid --min-peak-fraction value\n";
                return false;
            }
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return false;
        } else {
            cerr << "Unknown or incomplete argument: " << arg << endl;
            usage(argv[0]);
            return false;
        }
    }

    if (args->arc_table_path.empty() || args->out_path.empty()) {
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

bool load_arc_atac_rows(const Args& args, vector<ArcAtacRow>* rows) {
    if (rows == nullptr) {
        return false;
    }

    ifstream input(args.arc_table_path.c_str());
    if (!input.is_open()) {
        cerr << "Failed to open ARC table: " << args.arc_table_path << endl;
        return false;
    }

    string line;
    if (!std::getline(input, line)) {
        cerr << "ARC table is empty: " << args.arc_table_path << endl;
        return false;
    }
    trim_cr(&line);

    const vector<string> header = parse_tsv_line(line);
    unordered_map<string, size_t> header_index;
    for (size_t i = 0; i < header.size(); ++i) {
        header_index[header[i]] = i;
    }

    size_t barcode_idx = 0;
    size_t peak_cutsites_idx = 0;
    size_t peak_fragments_idx = 0;
    size_t atac_fragments_idx = 0;
    size_t excluded_reason_idx = 0;
    if (!find_required_column(header_index, args.barcode_col, &barcode_idx) ||
        !find_required_column(header_index, args.peak_cutsites_col, &peak_cutsites_idx) ||
        !find_required_column(header_index, args.peak_fragments_col, &peak_fragments_idx) ||
        !find_required_column(header_index, args.atac_fragments_col, &atac_fragments_idx) ||
        !find_required_column(header_index, args.excluded_reason_col, &excluded_reason_idx)) {
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
        size_t max_required_idx = barcode_idx;
        if (peak_cutsites_idx > max_required_idx) max_required_idx = peak_cutsites_idx;
        if (peak_fragments_idx > max_required_idx) max_required_idx = peak_fragments_idx;
        if (atac_fragments_idx > max_required_idx) max_required_idx = atac_fragments_idx;
        if (excluded_reason_idx > max_required_idx) max_required_idx = excluded_reason_idx;
        if (fields.size() <= max_required_idx) {
            cerr << "ARC table line " << line_number << " has too few columns" << endl;
            return false;
        }

        ArcAtacRow row;
        row.barcode = fields[barcode_idx];
        row.peak_cutsites_raw = fields[peak_cutsites_idx];
        row.peak_fragments_raw = fields[peak_fragments_idx];
        row.atac_fragments_raw = fields[atac_fragments_idx];
        row.excluded_reason = fields[excluded_reason_idx];

        if (seen.count(row.barcode) != 0) {
            cerr << "Duplicate barcode in ARC table: " << row.barcode << endl;
            return false;
        }
        seen.insert(row.barcode);

        if (!parse_uint64(row.peak_cutsites_raw, &row.peak_cutsites) ||
            !parse_uint64(row.peak_fragments_raw, &row.peak_fragments) ||
            !parse_uint64(row.atac_fragments_raw, &row.atac_fragments)) {
            cerr << "ARC table line " << line_number << " has invalid ATAC counts" << endl;
            return false;
        }

        rows->push_back(row);
    }

    return true;
}

string source_string(const Args& args) {
    std::ostringstream out;
    out << "arc_table_thresholds"
        << "(cutsites>=" << args.min_peak_region_cutsites
        << ",peak_fragments>=" << args.min_peak_region_fragments
        << ",atac_fragments>=" << args.min_atac_fragments
        << ",peak_fraction>=" << std::fixed << std::setprecision(6) << args.min_peak_fraction
        << ")";
    return out.str();
}

}  // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, &args)) {
        return 1;
    }

    vector<ArcAtacRow> rows;
    if (!load_arc_atac_rows(args, &rows)) {
        return 1;
    }

    ofstream output(args.out_path.c_str());
    if (!output.is_open()) {
        cerr << "Failed to open output: " << args.out_path << endl;
        return 1;
    }

    const string source = source_string(args);
    output << "barcode\tatac_peak_region_cutsites\tatac_peak_region_fragments\tatac_fragments\t"
              "atac_peak_fraction\texcluded_reason\tatac_module_call\tatac_low_targeting\t"
              "atac_source\n";

    uint64_t rows_written = 0;
    uint64_t module_call_count = 0;
    uint64_t low_targeting_count = 0;

    for (size_t i = 0; i < rows.size(); ++i) {
        const ArcAtacRow& row = rows[i];
        const double peak_fraction = (row.atac_fragments > 0)
            ? static_cast<double>(row.peak_fragments) / static_cast<double>(row.atac_fragments)
            : 0.0;
        const int atac_module_call =
            (row.peak_cutsites >= args.min_peak_region_cutsites &&
             row.peak_fragments >= args.min_peak_region_fragments &&
             row.atac_fragments >= args.min_atac_fragments) ? 1 : 0;
        const int atac_low_targeting =
            (row.atac_fragments > 0 && peak_fraction < args.min_peak_fraction) ? 1 : 0;

        output << row.barcode << '\t'
               << row.peak_cutsites_raw << '\t'
               << row.peak_fragments_raw << '\t'
               << row.atac_fragments_raw << '\t'
               << std::fixed << std::setprecision(6) << peak_fraction << '\t'
               << row.excluded_reason << '\t'
               << atac_module_call << '\t'
               << atac_low_targeting << '\t'
               << source << '\n';

        ++rows_written;
        if (atac_module_call) {
            ++module_call_count;
        }
        if (atac_low_targeting) {
            ++low_targeting_count;
        }
    }

    cerr << "Wrote " << rows_written << " ATAC evidence rows to " << args.out_path << endl;
    cerr << "atac_module_call=1 rows: " << module_call_count << endl;
    cerr << "atac_low_targeting=1 rows: " << low_targeting_count << endl;
    return 0;
}
