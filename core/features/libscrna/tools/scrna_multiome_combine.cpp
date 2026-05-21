#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using std::cerr;
using std::endl;
using std::ifstream;
using std::map;
using std::ofstream;
using std::string;
using std::unordered_map;
using std::vector;

namespace {

struct Args {
    string arc_table_path;
    string gex_evidence_path;
    string atac_evidence_path;
    string out_path;

    string mode = "gex_atac_rescue";

    string arc_barcode_col = "barcode";
    string gex_barcode_col = "barcode";
    string atac_barcode_col = "barcode";

    string gex_call_col = "gex_module_call";
    string gex_rescue_col;
    string atac_call_col = "atac_module_call";
    string atac_low_targeting_col;
};

struct GexEntry {
    int gex_module_call = 0;
    int has_gex_module_call = 0;
    int gex_rescue_eligible = 0;
    int has_gex_rescue_eligible = 0;
};

struct AtacEntry {
    int atac_module_call = 0;
    int has_atac_module_call = 0;
    int atac_low_targeting = 0;
    int has_atac_low_targeting = 0;
};

struct Decision {
    int final_is_cell = 0;
    int effective_atac_module_call = 0;
    string call_provenance;
};

void usage(const char* prog) {
    cerr << "Usage: " << prog
         << " --arc-table arc_barcode_table.tsv"
         << " --gex-evidence gex_evidence.tsv"
         << " --atac-evidence atac_evidence.tsv"
         << " --out multiome_calls.tsv"
         << " [--mode gex_atac_rescue|gex_only|union|intersection]"
         << " [--arc-barcode-col barcode]"
         << " [--gex-barcode-col barcode] [--gex-call-col gex_module_call]"
         << " [--gex-rescue-col gex_rescue_eligible]"
         << " [--atac-barcode-col barcode] [--atac-call-col atac_module_call]"
         << " [--atac-low-targeting-col atac_low_targeting]\n";
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

string to_lower_ascii(string value) {
    for (size_t i = 0; i < value.size(); ++i) {
        value[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(value[i])));
    }
    return value;
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

bool parse_args(int argc, char** argv, Args* args) {
    if (args == nullptr) {
        return false;
    }

    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];
        if (arg == "--arc-table" && i + 1 < argc) {
            args->arc_table_path = argv[++i];
        } else if (arg == "--gex-evidence" && i + 1 < argc) {
            args->gex_evidence_path = argv[++i];
        } else if (arg == "--atac-evidence" && i + 1 < argc) {
            args->atac_evidence_path = argv[++i];
        } else if (arg == "--out" && i + 1 < argc) {
            args->out_path = argv[++i];
        } else if (arg == "--mode" && i + 1 < argc) {
            args->mode = argv[++i];
        } else if (arg == "--arc-barcode-col" && i + 1 < argc) {
            args->arc_barcode_col = argv[++i];
        } else if (arg == "--gex-barcode-col" && i + 1 < argc) {
            args->gex_barcode_col = argv[++i];
        } else if (arg == "--gex-call-col" && i + 1 < argc) {
            args->gex_call_col = argv[++i];
        } else if (arg == "--gex-rescue-col" && i + 1 < argc) {
            args->gex_rescue_col = argv[++i];
        } else if (arg == "--atac-barcode-col" && i + 1 < argc) {
            args->atac_barcode_col = argv[++i];
        } else if (arg == "--atac-call-col" && i + 1 < argc) {
            args->atac_call_col = argv[++i];
        } else if (arg == "--atac-low-targeting-col" && i + 1 < argc) {
            args->atac_low_targeting_col = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return false;
        } else {
            cerr << "Unknown or incomplete argument: " << arg << endl;
            usage(argv[0]);
            return false;
        }
    }

    if (args->arc_table_path.empty() || args->gex_evidence_path.empty() ||
        args->atac_evidence_path.empty() || args->out_path.empty()) {
        usage(argv[0]);
        return false;
    }

    if (args->mode != "gex_atac_rescue" &&
        args->mode != "gex_only" &&
        args->mode != "union" &&
        args->mode != "intersection") {
        cerr << "Unsupported mode: " << args->mode << endl;
        usage(argv[0]);
        return false;
    }

    if (args->mode == "gex_atac_rescue" && args->gex_rescue_col.empty()) {
        cerr << "--gex-rescue-col is required for mode gex_atac_rescue" << endl;
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

bool find_optional_column(const unordered_map<string, size_t>& header_index,
                          const string& name,
                          bool* found_out,
                          size_t* index_out) {
    if (found_out == nullptr || index_out == nullptr) {
        return false;
    }

    if (name.empty()) {
        *found_out = false;
        *index_out = 0;
        return true;
    }

    const auto it = header_index.find(name);
    if (it == header_index.end()) {
        cerr << "Missing optional column requested on command line: " << name << endl;
        return false;
    }

    *found_out = true;
    *index_out = it->second;
    return true;
}

bool load_arc_barcodes(const string& path,
                       const string& barcode_col,
                       vector<string>* barcodes) {
    if (barcodes == nullptr) {
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
    if (!find_required_column(header_index, barcode_col, &barcode_idx)) {
        return false;
    }

    unordered_map<string, int> seen;
    uint64_t line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        trim_cr(&line);
        if (line.empty()) {
            continue;
        }

        const vector<string> fields = parse_tsv_line(line);
        if (fields.size() <= barcode_idx) {
            cerr << "ARC table line " << line_number << " has too few columns" << endl;
            return false;
        }

        const string& barcode = fields[barcode_idx];
        if (seen.count(barcode) != 0) {
            cerr << "Duplicate barcode in ARC table: " << barcode << endl;
            return false;
        }
        seen[barcode] = 1;
        barcodes->push_back(barcode);
    }

    return true;
}

bool load_gex_evidence(const string& path,
                       const string& barcode_col,
                       const string& gex_call_col,
                       const string& gex_rescue_col,
                       unordered_map<string, GexEntry>* by_barcode) {
    if (by_barcode == nullptr) {
        return false;
    }

    ifstream input(path.c_str());
    if (!input.is_open()) {
        cerr << "Failed to open GEX evidence: " << path << endl;
        return false;
    }

    string line;
    if (!std::getline(input, line)) {
        cerr << "GEX evidence is empty: " << path << endl;
        return false;
    }
    trim_cr(&line);

    const vector<string> header = parse_tsv_line(line);
    unordered_map<string, size_t> header_index;
    for (size_t i = 0; i < header.size(); ++i) {
        header_index[header[i]] = i;
    }

    size_t barcode_idx = 0;
    size_t gex_call_idx = 0;
    size_t gex_rescue_idx = 0;
    bool has_gex_rescue_col = false;
    if (!find_required_column(header_index, barcode_col, &barcode_idx) ||
        !find_required_column(header_index, gex_call_col, &gex_call_idx) ||
        !find_optional_column(header_index, gex_rescue_col, &has_gex_rescue_col, &gex_rescue_idx)) {
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
        size_t max_required_idx = barcode_idx > gex_call_idx ? barcode_idx : gex_call_idx;
        if (has_gex_rescue_col && gex_rescue_idx > max_required_idx) {
            max_required_idx = gex_rescue_idx;
        }
        if (fields.size() <= max_required_idx) {
            cerr << "GEX evidence line " << line_number << " has too few columns" << endl;
            return false;
        }

        GexEntry entry;
        if (!parse_binary_call(fields[gex_call_idx], &entry.gex_module_call)) {
            cerr << "GEX evidence line " << line_number << " has invalid "
                 << gex_call_col << " value: " << fields[gex_call_idx] << endl;
            return false;
        }
        entry.has_gex_module_call = 1;

        if (has_gex_rescue_col) {
            if (!parse_binary_call(fields[gex_rescue_idx], &entry.gex_rescue_eligible)) {
                cerr << "GEX evidence line " << line_number << " has invalid "
                     << gex_rescue_col << " value: " << fields[gex_rescue_idx] << endl;
                return false;
            }
            entry.has_gex_rescue_eligible = 1;
        }

        const string& barcode = fields[barcode_idx];
        if (by_barcode->count(barcode) != 0) {
            cerr << "Duplicate barcode in GEX evidence: " << barcode << endl;
            return false;
        }

        (*by_barcode)[barcode] = entry;
    }

    return true;
}

bool load_atac_evidence(const string& path,
                        const string& barcode_col,
                        const string& atac_call_col,
                        const string& atac_low_targeting_col,
                        unordered_map<string, AtacEntry>* by_barcode) {
    if (by_barcode == nullptr) {
        return false;
    }

    ifstream input(path.c_str());
    if (!input.is_open()) {
        cerr << "Failed to open ATAC evidence: " << path << endl;
        return false;
    }

    string line;
    if (!std::getline(input, line)) {
        cerr << "ATAC evidence is empty: " << path << endl;
        return false;
    }
    trim_cr(&line);

    const vector<string> header = parse_tsv_line(line);
    unordered_map<string, size_t> header_index;
    for (size_t i = 0; i < header.size(); ++i) {
        header_index[header[i]] = i;
    }

    size_t barcode_idx = 0;
    size_t atac_call_idx = 0;
    size_t atac_low_targeting_idx = 0;
    bool has_low_targeting_col = false;
    if (!find_required_column(header_index, barcode_col, &barcode_idx) ||
        !find_required_column(header_index, atac_call_col, &atac_call_idx) ||
        !find_optional_column(header_index, atac_low_targeting_col,
                              &has_low_targeting_col, &atac_low_targeting_idx)) {
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
        size_t max_required_idx = barcode_idx > atac_call_idx ? barcode_idx : atac_call_idx;
        if (has_low_targeting_col && atac_low_targeting_idx > max_required_idx) {
            max_required_idx = atac_low_targeting_idx;
        }
        if (fields.size() <= max_required_idx) {
            cerr << "ATAC evidence line " << line_number << " has too few columns" << endl;
            return false;
        }

        AtacEntry entry;
        if (!parse_binary_call(fields[atac_call_idx], &entry.atac_module_call)) {
            cerr << "ATAC evidence line " << line_number << " has invalid "
                 << atac_call_col << " value: " << fields[atac_call_idx] << endl;
            return false;
        }
        entry.has_atac_module_call = 1;

        if (has_low_targeting_col) {
            if (!parse_binary_call(fields[atac_low_targeting_idx], &entry.atac_low_targeting)) {
                cerr << "ATAC evidence line " << line_number << " has invalid "
                     << atac_low_targeting_col << " value: "
                     << fields[atac_low_targeting_idx] << endl;
                return false;
            }
            entry.has_atac_low_targeting = 1;
        }

        const string& barcode = fields[barcode_idx];
        if (by_barcode->count(barcode) != 0) {
            cerr << "Duplicate barcode in ATAC evidence: " << barcode << endl;
            return false;
        }

        (*by_barcode)[barcode] = entry;
    }

    return true;
}

Decision combine_calls(const Args& args,
                       int gex_call,
                       int gex_rescue_eligible,
                       int atac_call,
                       int atac_low_targeting) {
    Decision decision;
    const int effective_atac_call = atac_low_targeting ? 0 : atac_call;
    decision.effective_atac_module_call = effective_atac_call;

    if (args.mode == "gex_only") {
        decision.final_is_cell = gex_call;
        decision.call_provenance = gex_call ? "gex_module_call" : "no_call";
        return decision;
    }

    if (args.mode == "union") {
        if (gex_call) {
            decision.final_is_cell = 1;
            decision.call_provenance = "gex_module_call";
        } else if (effective_atac_call) {
            decision.final_is_cell = 1;
            decision.call_provenance = "atac_module_call";
        } else if (atac_call && atac_low_targeting) {
            decision.final_is_cell = 0;
            decision.call_provenance = "atac_low_targeting_block";
        } else {
            decision.final_is_cell = 0;
            decision.call_provenance = "no_call";
        }
        return decision;
    }

    if (args.mode == "intersection") {
        if (gex_call && effective_atac_call) {
            decision.final_is_cell = 1;
            decision.call_provenance = "gex_and_atac_call";
        } else if (atac_call && atac_low_targeting) {
            decision.final_is_cell = 0;
            decision.call_provenance = "atac_low_targeting_block";
        } else {
            decision.final_is_cell = 0;
            decision.call_provenance = "no_call";
        }
        return decision;
    }

    if (gex_call) {
        decision.final_is_cell = 1;
        decision.call_provenance = "gex_module_call";
    } else if (atac_call && atac_low_targeting) {
        decision.final_is_cell = 0;
        decision.call_provenance = "atac_low_targeting_block";
    } else if (effective_atac_call && gex_rescue_eligible) {
        decision.final_is_cell = 1;
        decision.call_provenance = "gex_atac_rescue";
    } else if (effective_atac_call) {
        decision.final_is_cell = 0;
        decision.call_provenance = "atac_call_not_rescue_eligible";
    } else {
        decision.final_is_cell = 0;
        decision.call_provenance = "no_call";
    }

    return decision;
}

}  // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, &args)) {
        return 1;
    }

    vector<string> arc_barcodes;
    if (!load_arc_barcodes(args.arc_table_path, args.arc_barcode_col, &arc_barcodes)) {
        return 1;
    }

    unordered_map<string, GexEntry> gex_by_barcode;
    if (!load_gex_evidence(args.gex_evidence_path, args.gex_barcode_col,
                           args.gex_call_col, args.gex_rescue_col, &gex_by_barcode)) {
        return 1;
    }

    unordered_map<string, AtacEntry> atac_by_barcode;
    if (!load_atac_evidence(args.atac_evidence_path, args.atac_barcode_col,
                            args.atac_call_col, args.atac_low_targeting_col, &atac_by_barcode)) {
        return 1;
    }

    ofstream out(args.out_path.c_str());
    if (!out.is_open()) {
        cerr << "Failed to open output: " << args.out_path << endl;
        return 1;
    }

    out << "barcode\tgex_module_call\tatac_module_call\teffective_atac_module_call\t"
           "gex_rescue_eligible\tatac_low_targeting\tfinal_is_cell\tcall_provenance\t"
           "has_gex_evidence\thas_atac_evidence\n";

    uint64_t final_is_cell_count = 0;
    uint64_t missing_gex = 0;
    uint64_t missing_atac = 0;
    uint64_t rows_written = 0;
    map<string, uint64_t> provenance_counts;

    for (size_t i = 0; i < arc_barcodes.size(); ++i) {
        const string& barcode = arc_barcodes[i];

        const auto gex_it = gex_by_barcode.find(barcode);
        const auto atac_it = atac_by_barcode.find(barcode);

        const int has_gex_evidence = (gex_it != gex_by_barcode.end()) ? 1 : 0;
        const int has_atac_evidence = (atac_it != atac_by_barcode.end()) ? 1 : 0;
        if (!has_gex_evidence) {
            ++missing_gex;
        }
        if (!has_atac_evidence) {
            ++missing_atac;
        }

        const int gex_call = has_gex_evidence ? gex_it->second.gex_module_call : 0;
        const int gex_rescue_eligible = (has_gex_evidence && gex_it->second.has_gex_rescue_eligible)
            ? gex_it->second.gex_rescue_eligible : 0;
        const int atac_call = has_atac_evidence ? atac_it->second.atac_module_call : 0;
        const int atac_low_targeting = (has_atac_evidence && atac_it->second.has_atac_low_targeting)
            ? atac_it->second.atac_low_targeting : 0;

        const Decision decision = combine_calls(args, gex_call, gex_rescue_eligible,
                                                atac_call, atac_low_targeting);
        if (decision.final_is_cell == 1) {
            ++final_is_cell_count;
        }
        provenance_counts[decision.call_provenance]++;

        out << barcode << '\t'
            << gex_call << '\t'
            << atac_call << '\t'
            << decision.effective_atac_module_call << '\t'
            << gex_rescue_eligible << '\t'
            << atac_low_targeting << '\t'
            << decision.final_is_cell << '\t'
            << decision.call_provenance << '\t'
            << has_gex_evidence << '\t'
            << has_atac_evidence << '\n';

        ++rows_written;
    }

    cerr << "Wrote " << rows_written << " multiome barcode rows to " << args.out_path << endl;
    cerr << "Final is_cell count: " << final_is_cell_count << endl;
    cerr << "Missing GEX evidence rows: " << missing_gex << endl;
    cerr << "Missing ATAC evidence rows: " << missing_atac << endl;
    for (map<string, uint64_t>::const_iterator it = provenance_counts.begin();
         it != provenance_counts.end(); ++it) {
        cerr << "call_provenance[" << it->first << "]=" << it->second << endl;
    }

    return 0;
}
