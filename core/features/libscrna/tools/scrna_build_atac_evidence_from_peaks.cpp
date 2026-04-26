#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::cerr;
using std::endl;
using std::ifstream;
using std::map;
using std::ofstream;
using std::pair;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

namespace {

struct Args {
    string fragments_path;
    string peaks_path;
    string out_path;
    string whitelist_path;
    string barcode_suffix;
};

struct PeakInterval {
    int64_t start = 0;
    int64_t end = 0;
};

struct BarcodeCounts {
    uint64_t atac_fragments = 0;
    uint64_t peak_fragments = 0;
    uint64_t peak_cutsites = 0;
};

void usage(const char* prog) {
    cerr << "Usage: " << prog
         << " --fragments fragments.tsv[.gz]"
         << " --peaks peaks.narrowPeak"
         << " --out atac_evidence.tsv"
         << " [--barcode-whitelist whitelist.txt]"
         << " [--barcode-suffix -1]\n"
         << "Reads a Chromap/ARC 5-col fragments TSV (chrom start end barcode count)\n"
         << "and a narrowPeak file (BED-like; first 3 columns used) and emits per-barcode\n"
         << "ATAC evidence with the following columns:\n"
         << "  barcode, atac_peak_region_cutsites, atac_peak_region_fragments,\n"
         << "  atac_fragments, atac_peak_fraction\n";
}

void trim_cr(string* value) {
    if (value != nullptr && !value->empty() && value->back() == '\r') {
        value->pop_back();
    }
}

vector<string> split(const string& line, char sep) {
    vector<string> fields;
    string current;
    for (size_t i = 0; i < line.size(); ++i) {
        if (line[i] == sep) {
            fields.push_back(current);
            current.clear();
        } else {
            current.push_back(line[i]);
        }
    }
    fields.push_back(current);
    return fields;
}

bool parse_int64(const string& value, int64_t* out) {
    if (out == nullptr || value.empty()) {
        return false;
    }
    char* end = nullptr;
    const long long parsed = std::strtoll(value.c_str(), &end, 10);
    if (end == value.c_str() || *end != '\0') {
        return false;
    }
    *out = static_cast<int64_t>(parsed);
    return true;
}

bool parse_args(int argc, char** argv, Args* args) {
    if (args == nullptr) {
        return false;
    }
    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];
        if (arg == "--fragments" && i + 1 < argc) {
            args->fragments_path = argv[++i];
        } else if (arg == "--peaks" && i + 1 < argc) {
            args->peaks_path = argv[++i];
        } else if (arg == "--out" && i + 1 < argc) {
            args->out_path = argv[++i];
        } else if (arg == "--barcode-whitelist" && i + 1 < argc) {
            args->whitelist_path = argv[++i];
        } else if (arg == "--barcode-suffix" && i + 1 < argc) {
            args->barcode_suffix = argv[++i];
        } else if (arg == "--help" || arg == "-h") {
            usage(argv[0]);
            return false;
        } else {
            cerr << "Unknown or incomplete argument: " << arg << endl;
            usage(argv[0]);
            return false;
        }
    }
    if (args->fragments_path.empty() || args->peaks_path.empty() || args->out_path.empty()) {
        usage(argv[0]);
        return false;
    }
    return true;
}

bool load_whitelist(const string& path, const string& suffix,
                    unordered_set<string>* whitelist) {
    if (path.empty()) {
        return true;
    }
    ifstream in(path.c_str());
    if (!in.is_open()) {
        cerr << "Failed to open barcode whitelist: " << path << endl;
        return false;
    }
    string line;
    while (std::getline(in, line)) {
        trim_cr(&line);
        if (line.empty()) {
            continue;
        }
        whitelist->insert(line + suffix);
    }
    return true;
}

bool load_peaks(const string& path,
                unordered_map<string, vector<PeakInterval>>* peaks_by_chrom) {
    ifstream in(path.c_str());
    if (!in.is_open()) {
        cerr << "Failed to open peaks file: " << path << endl;
        return false;
    }
    string line;
    uint64_t line_no = 0;
    uint64_t kept = 0;
    while (std::getline(in, line)) {
        ++line_no;
        trim_cr(&line);
        if (line.empty() || line[0] == '#') {
            continue;
        }
        if (line.size() >= 5 && line.compare(0, 5, "track") == 0) {
            continue;
        }
        if (line.size() >= 7 && line.compare(0, 7, "browser") == 0) {
            continue;
        }
        const vector<string> fields = split(line, '\t');
        if (fields.size() < 3) {
            cerr << "Peak line " << line_no << " has fewer than 3 columns" << endl;
            return false;
        }
        PeakInterval iv;
        if (!parse_int64(fields[1], &iv.start) || !parse_int64(fields[2], &iv.end)) {
            cerr << "Peak line " << line_no << " has non-numeric start/end" << endl;
            return false;
        }
        if (iv.end <= iv.start) {
            continue;
        }
        (*peaks_by_chrom)[fields[0]].push_back(iv);
        ++kept;
    }
    for (auto& kv : *peaks_by_chrom) {
        std::sort(kv.second.begin(), kv.second.end(),
                  [](const PeakInterval& a, const PeakInterval& b) {
                      if (a.start != b.start) return a.start < b.start;
                      return a.end < b.end;
                  });
        // Merge overlapping/adjacent peaks so overlap queries are O(log n).
        vector<PeakInterval>& v = kv.second;
        size_t w = 0;
        for (size_t r = 0; r < v.size(); ++r) {
            if (w > 0 && v[r].start <= v[w - 1].end) {
                if (v[r].end > v[w - 1].end) {
                    v[w - 1].end = v[r].end;
                }
            } else {
                v[w++] = v[r];
            }
        }
        v.resize(w);
    }
    cerr << "Loaded " << kept << " peak intervals across "
         << peaks_by_chrom->size() << " chromosomes" << endl;
    return true;
}

bool fragment_overlaps_any_peak(const vector<PeakInterval>& peaks,
                                int64_t start, int64_t end) {
    if (peaks.empty() || end <= start) {
        return false;
    }
    // Binary search for first peak whose end > fragment.start.
    size_t lo = 0;
    size_t hi = peaks.size();
    while (lo < hi) {
        const size_t mid = lo + (hi - lo) / 2;
        if (peaks[mid].end <= start) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if (lo >= peaks.size()) {
        return false;
    }
    return peaks[lo].start < end;
}

bool position_in_peak(const vector<PeakInterval>& peaks, int64_t pos) {
    if (peaks.empty()) {
        return false;
    }
    size_t lo = 0;
    size_t hi = peaks.size();
    while (lo < hi) {
        const size_t mid = lo + (hi - lo) / 2;
        if (peaks[mid].end <= pos) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if (lo >= peaks.size()) {
        return false;
    }
    return peaks[lo].start <= pos;
}

class FragmentReader {
public:
    explicit FragmentReader(const string& path)
        : path_(path), pipe_(nullptr), file_(nullptr), finished_(false), exit_status_(0) {
        const bool gz = path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
        if (gz) {
            // popen with /bin/sh -c; quoting with single quotes is sufficient for
            // typical fixture paths (no embedded single quotes). Caller should
            // pass real filesystem paths.
            const string cmd = "gzip -cd '" + path + "'";
            pipe_ = popen(cmd.c_str(), "r");
            if (pipe_ == nullptr) {
                cerr << "Failed to popen: " << cmd << endl;
            }
        } else {
            file_ = std::fopen(path.c_str(), "r");
            if (file_ == nullptr) {
                cerr << "Failed to open fragments file: " << path << endl;
            }
        }
    }

    ~FragmentReader() {
        // Best-effort cleanup; status from finish() is what callers should check.
        if (!finished_) {
            (void)finish();
        }
    }

    bool ok() const { return pipe_ != nullptr || file_ != nullptr; }

    bool getline(string* out) {
        if (out == nullptr) {
            return false;
        }
        out->clear();
        FILE* src = (pipe_ != nullptr) ? pipe_ : file_;
        if (src == nullptr) {
            return false;
        }
        char buf[8192];
        bool any = false;
        while (std::fgets(buf, sizeof(buf), src) != nullptr) {
            any = true;
            out->append(buf);
            if (!out->empty() && out->back() == '\n') {
                out->pop_back();
                trim_cr(out);
                return true;
            }
        }
        if (any) {
            trim_cr(out);
            return true;
        }
        return false;
    }

    // Closes the underlying stream and returns true iff EOF was reached cleanly
    // and (for gzip pipes) the child process exited with status 0. Idempotent.
    bool finish() {
        if (finished_) {
            return exit_status_ == 0 && !read_error_;
        }
        finished_ = true;
        FILE* src = (pipe_ != nullptr) ? pipe_ : file_;
        if (src != nullptr && std::ferror(src)) {
            cerr << "Read error on fragments stream: " << path_ << endl;
            read_error_ = true;
        }
        if (pipe_ != nullptr) {
            const int rc = pclose(pipe_);
            pipe_ = nullptr;
            if (rc == -1) {
                cerr << "pclose failed for fragments pipe: " << path_ << endl;
                exit_status_ = -1;
            } else if (rc != 0) {
                cerr << "Fragments decompression pipeline exited with status "
                     << rc << " for: " << path_ << endl;
                exit_status_ = rc;
            }
        }
        if (file_ != nullptr) {
            std::fclose(file_);
            file_ = nullptr;
        }
        return exit_status_ == 0 && !read_error_;
    }

private:
    string path_;
    FILE* pipe_;
    FILE* file_;
    bool finished_;
    bool read_error_ = false;
    int exit_status_;
};

}  // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, &args)) {
        return 1;
    }

    unordered_set<string> whitelist;
    if (!load_whitelist(args.whitelist_path, args.barcode_suffix, &whitelist)) {
        return 1;
    }
    const bool have_whitelist = !args.whitelist_path.empty();

    unordered_map<string, vector<PeakInterval>> peaks_by_chrom;
    if (!load_peaks(args.peaks_path, &peaks_by_chrom)) {
        return 1;
    }

    FragmentReader reader(args.fragments_path);
    if (!reader.ok()) {
        return 1;
    }

    unordered_map<string, BarcodeCounts> per_barcode;
    string line;
    uint64_t lineno = 0;
    uint64_t total_rows = 0;
    uint64_t total_fragments = 0;
    uint64_t total_peak_fragments = 0;
    uint64_t total_peak_cutsites = 0;
    while (reader.getline(&line)) {
        ++lineno;
        if (line.empty() || line[0] == '#') {
            continue;
        }
        const vector<string> fields = split(line, '\t');
        if (fields.size() < 4) {
            cerr << "Fragment line " << lineno << " has fewer than 4 columns" << endl;
            return 1;
        }
        int64_t start = 0;
        int64_t end = 0;
        if (!parse_int64(fields[1], &start) || !parse_int64(fields[2], &end)) {
            cerr << "Fragment line " << lineno << " has non-numeric start/end" << endl;
            return 1;
        }
        if (end <= start) {
            continue;
        }
        // Column 5 is the fragment-support count (Chromap/ARC convention).
        // ARC parity treats this as a per-row weight: a row with count=N
        // contributes N to atac_fragments and (when overlapping a peak) N to
        // atac_peak_region_fragments; each cut-site endpoint that falls in a
        // peak contributes N to atac_peak_region_cutsites.
        uint64_t count = 1;
        if (fields.size() >= 5 && !fields[4].empty()) {
            int64_t parsed = 0;
            if (!parse_int64(fields[4], &parsed) || parsed < 0) {
                cerr << "Fragment line " << lineno
                     << " has invalid count column: '" << fields[4] << "'" << endl;
                return 1;
            }
            count = static_cast<uint64_t>(parsed);
            if (count == 0) {
                continue;
            }
        }
        const string& chrom = fields[0];
        const string& barcode = fields[3];
        if (have_whitelist && whitelist.count(barcode) == 0) {
            continue;
        }
        BarcodeCounts& counts = per_barcode[barcode];
        counts.atac_fragments += count;
        total_fragments += count;
        ++total_rows;

        const auto it = peaks_by_chrom.find(chrom);
        if (it == peaks_by_chrom.end()) {
            continue;
        }
        const vector<PeakInterval>& peaks = it->second;
        if (fragment_overlaps_any_peak(peaks, start, end)) {
            counts.peak_fragments += count;
            total_peak_fragments += count;
        }
        // Cut sites: Tn5 nicks at the two fragment endpoints. Use [start, end-1]
        // (BED half-open: end-1 is the last covered base). Each endpoint that
        // lands in a peak contributes `count` to atac_peak_region_cutsites.
        if (position_in_peak(peaks, start)) {
            counts.peak_cutsites += count;
            total_peak_cutsites += count;
        }
        if (position_in_peak(peaks, end - 1)) {
            counts.peak_cutsites += count;
            total_peak_cutsites += count;
        }
    }
    if (!reader.finish()) {
        cerr << "Fragments stream ended with error; refusing to write evidence." << endl;
        return 1;
    }

    map<string, BarcodeCounts> sorted(per_barcode.begin(), per_barcode.end());

    ofstream output(args.out_path.c_str());
    if (!output.is_open()) {
        cerr << "Failed to open output: " << args.out_path << endl;
        return 1;
    }
    output << "barcode\tatac_peak_region_cutsites\tatac_peak_region_fragments\t"
              "atac_fragments\tatac_peak_fraction\n";
    output << std::fixed << std::setprecision(6);
    for (const auto& kv : sorted) {
        const BarcodeCounts& c = kv.second;
        const double frac = (c.atac_fragments > 0)
            ? static_cast<double>(c.peak_fragments) /
              static_cast<double>(c.atac_fragments)
            : 0.0;
        output << kv.first << '\t'
               << c.peak_cutsites << '\t'
               << c.peak_fragments << '\t'
               << c.atac_fragments << '\t'
               << frac << '\n';
    }

    cerr << "Wrote evidence for " << sorted.size() << " barcodes to "
         << args.out_path << endl;
    cerr << "Fragment rows kept: " << total_rows << endl;
    cerr << "Sum of count column (atac_fragments total): " << total_fragments << endl;
    cerr << "Sum of count column for fragments overlapping peaks: "
         << total_peak_fragments << endl;
    cerr << "Sum of count column over peak-resident cut-site endpoints: "
         << total_peak_cutsites << endl;
    return 0;
}
