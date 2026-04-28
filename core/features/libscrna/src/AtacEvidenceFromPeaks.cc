#include "libscrna/AtacEvidenceFromPeaks.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace libscrna {
namespace atac {

namespace {

struct PeakInterval {
    int64_t start = 0;
    int64_t end = 0;
};

struct BarcodeCounts {
    uint64_t atac_fragments = 0;
    uint64_t peak_fragments = 0;
    uint64_t peak_cutsites = 0;
};

void trim_cr(std::string* value) {
    if (value != nullptr && !value->empty() && value->back() == '\r') {
        value->pop_back();
    }
}

std::vector<std::string> split(const std::string& line, char sep) {
    std::vector<std::string> fields;
    std::string current;
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

bool parse_int64(const std::string& value, int64_t* out) {
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

bool load_whitelist(const std::string& path, const std::string& suffix,
                    std::unordered_set<std::string>* whitelist,
                    std::ostream& err) {
    if (path.empty()) {
        return true;
    }
    std::ifstream in(path.c_str());
    if (!in.is_open()) {
        err << "Failed to open barcode whitelist: " << path << "\n";
        return false;
    }
    std::string line;
    while (std::getline(in, line)) {
        trim_cr(&line);
        if (line.empty()) {
            continue;
        }
        whitelist->insert(line + suffix);
    }
    return true;
}

bool load_peaks(const std::string& path,
                std::unordered_map<std::string, std::vector<PeakInterval>>*
                    peaks_by_chrom,
                std::ostream& err) {
    std::ifstream in(path.c_str());
    if (!in.is_open()) {
        err << "Failed to open peaks file: " << path << "\n";
        return false;
    }
    std::string line;
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
        const std::vector<std::string> fields = split(line, '\t');
        if (fields.size() < 3) {
            err << "Peak line " << line_no << " has fewer than 3 columns\n";
            return false;
        }
        PeakInterval iv;
        if (!parse_int64(fields[1], &iv.start) ||
            !parse_int64(fields[2], &iv.end)) {
            err << "Peak line " << line_no << " has non-numeric start/end\n";
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
        std::vector<PeakInterval>& v = kv.second;
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
    err << "Loaded " << kept << " peak intervals across "
        << peaks_by_chrom->size() << " chromosomes\n";
    return true;
}

bool fragment_overlaps_any_peak(const std::vector<PeakInterval>& peaks,
                                int64_t start, int64_t end) {
    if (peaks.empty() || end <= start) {
        return false;
    }
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

bool position_in_peak(const std::vector<PeakInterval>& peaks, int64_t pos) {
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
    FragmentReader(const std::string& path, std::ostream& err)
        : path_(path), err_(err), pipe_(nullptr), file_(nullptr),
          finished_(false), exit_status_(0) {
        const bool gz = path.size() >= 3 &&
                        path.compare(path.size() - 3, 3, ".gz") == 0;
        if (gz) {
            const std::string cmd = "gzip -cd '" + path + "'";
            pipe_ = popen(cmd.c_str(), "r");
            if (pipe_ == nullptr) {
                err_ << "Failed to popen: " << cmd << "\n";
            }
        } else {
            file_ = std::fopen(path.c_str(), "r");
            if (file_ == nullptr) {
                err_ << "Failed to open fragments file: " << path << "\n";
            }
        }
    }

    ~FragmentReader() {
        if (!finished_) {
            (void)finish();
        }
    }

    bool ok() const { return pipe_ != nullptr || file_ != nullptr; }

    bool getline(std::string* out) {
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

    bool finish() {
        if (finished_) {
            return exit_status_ == 0 && !read_error_;
        }
        finished_ = true;
        FILE* src = (pipe_ != nullptr) ? pipe_ : file_;
        if (src != nullptr && std::ferror(src)) {
            err_ << "Read error on fragments stream: " << path_ << "\n";
            read_error_ = true;
        }
        if (pipe_ != nullptr) {
            const int rc = pclose(pipe_);
            pipe_ = nullptr;
            if (rc == -1) {
                err_ << "pclose failed for fragments pipe: " << path_ << "\n";
                exit_status_ = -1;
            } else if (rc != 0) {
                err_ << "Fragments decompression pipeline exited with status "
                     << rc << " for: " << path_ << "\n";
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
    std::string path_;
    std::ostream& err_;
    FILE* pipe_;
    FILE* file_;
    bool finished_;
    bool read_error_ = false;
    int exit_status_;
};

}  // namespace

int RunAtacEvidenceFromPeaks(const AtacEvidenceFromPeaksOptions& opts,
                             std::ostream* err_in) {
    std::ostream& err = (err_in != nullptr) ? *err_in : std::cerr;

    if (opts.fragments_path.empty() || opts.peaks_path.empty() ||
        opts.out_path.empty()) {
        err << "RunAtacEvidenceFromPeaks: fragments_path, peaks_path and "
               "out_path are all required.\n";
        return 1;
    }

    std::unordered_set<std::string> whitelist;
    if (!load_whitelist(opts.whitelist_path, opts.barcode_suffix, &whitelist,
                         err)) {
        return 1;
    }
    const bool have_whitelist = !opts.whitelist_path.empty();

    std::unordered_map<std::string, std::vector<PeakInterval>> peaks_by_chrom;
    if (!load_peaks(opts.peaks_path, &peaks_by_chrom, err)) {
        return 1;
    }

    FragmentReader reader(opts.fragments_path, err);
    if (!reader.ok()) {
        return 1;
    }

    std::unordered_map<std::string, BarcodeCounts> per_barcode;
    std::string line;
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
        const std::vector<std::string> fields = split(line, '\t');
        if (fields.size() < 4) {
            err << "Fragment line " << lineno
                << " has fewer than 4 columns\n";
            return 1;
        }
        int64_t start = 0;
        int64_t end = 0;
        if (!parse_int64(fields[1], &start) ||
            !parse_int64(fields[2], &end)) {
            err << "Fragment line " << lineno
                << " has non-numeric start/end\n";
            return 1;
        }
        if (end <= start) {
            continue;
        }
        // Column 5 is the fragment-support count (Chromap/ARC convention).
        // ARC parity treats this as a per-row weight.
        uint64_t count = 1;
        if (fields.size() >= 5 && !fields[4].empty()) {
            int64_t parsed = 0;
            if (!parse_int64(fields[4], &parsed) || parsed < 0) {
                err << "Fragment line " << lineno
                    << " has invalid count column: '" << fields[4] << "'\n";
                return 1;
            }
            count = static_cast<uint64_t>(parsed);
            if (count == 0) {
                continue;
            }
        }
        const std::string& chrom = fields[0];
        const std::string& barcode = fields[3];
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
        const std::vector<PeakInterval>& peaks = it->second;
        if (fragment_overlaps_any_peak(peaks, start, end)) {
            counts.peak_fragments += count;
            total_peak_fragments += count;
        }
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
        err << "Fragments stream ended with error; refusing to write evidence.\n";
        return 1;
    }

    std::map<std::string, BarcodeCounts> sorted(per_barcode.begin(),
                                                  per_barcode.end());

    std::ofstream output(opts.out_path.c_str());
    if (!output.is_open()) {
        err << "Failed to open output: " << opts.out_path << "\n";
        return 1;
    }
    output << "barcode\tatac_peak_region_cutsites\tatac_peak_region_fragments\t"
              "atac_fragments\tatac_peak_fraction\n";
    output << std::fixed << std::setprecision(6);
    for (const auto& kv : sorted) {
        const BarcodeCounts& c = kv.second;
        const double frac =
            (c.atac_fragments > 0)
                ? static_cast<double>(c.peak_fragments) /
                      static_cast<double>(c.atac_fragments)
                : 0.0;
        output << kv.first << '\t' << c.peak_cutsites << '\t'
               << c.peak_fragments << '\t' << c.atac_fragments << '\t' << frac
               << '\n';
    }

    err << "Wrote evidence for " << sorted.size() << " barcodes to "
        << opts.out_path << "\n"
        << "Fragment rows kept: " << total_rows << "\n"
        << "Sum of count column (atac_fragments total): " << total_fragments
        << "\n"
        << "Sum of count column for fragments overlapping peaks: "
        << total_peak_fragments << "\n"
        << "Sum of count column over peak-resident cut-site endpoints: "
        << total_peak_cutsites << "\n";
    return 0;
}

}  // namespace atac
}  // namespace libscrna
