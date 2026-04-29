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
#include <sys/types.h>
#include <unordered_map>
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

struct AtacEvidenceBinaryHeader {
    char magic[4];
    uint32_t format_version;
    uint32_t record_size;
    uint32_t barcode_length;
    uint32_t num_chroms;
    uint32_t flags;
    uint64_t num_records;
};

struct AtacEvidenceBinaryRecord {
    int32_t chrom_id;
    int32_t start;
    int32_t end;
    uint32_t count;
    uint64_t barcode_key;
};

static_assert(sizeof(AtacEvidenceBinaryHeader) == 32,
              "ATAC evidence binary header must be 32 bytes");
static_assert(sizeof(AtacEvidenceBinaryRecord) == 24,
              "ATAC evidence binary record must be 24 bytes");

bool load_chroms(const std::string& path, std::vector<std::string>* chrom_names,
                 std::ostream& err) {
    std::ifstream in(path.c_str());
    if (!in.is_open()) {
        err << "Failed to open chrom metadata: " << path << "\n";
        return false;
    }
    std::string line;
    uint64_t line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        trim_cr(&line);
        if (line.empty()) {
            continue;
        }
        const std::vector<std::string> fields = split(line, '\t');
        if (fields.size() < 2) {
            err << "Chrom metadata line " << line_no
                << " has fewer than 2 columns\n";
            return false;
        }
        int64_t chrom_id = -1;
        if (!parse_int64(fields[0], &chrom_id) || chrom_id < 0) {
            err << "Chrom metadata line " << line_no
                << " has invalid chrom_id\n";
            return false;
        }
        if (static_cast<size_t>(chrom_id) >= chrom_names->size()) {
            chrom_names->resize(static_cast<size_t>(chrom_id) + 1);
        }
        (*chrom_names)[static_cast<size_t>(chrom_id)] = fields[1];
    }
    return true;
}

int write_evidence_output(const std::map<std::string, BarcodeCounts>& sorted,
                          const std::string& out_path,
                          uint64_t total_rows,
                          uint64_t total_fragments,
                          uint64_t total_peak_fragments,
                          uint64_t total_peak_cutsites,
                          std::ostream& err) {
    std::ofstream output(out_path.c_str());
    if (!output.is_open()) {
        err << "Failed to open output: " << out_path << "\n";
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
        << out_path << "\n"
        << "Fragment rows kept: " << total_rows << "\n"
        << "Sum of count column (atac_fragments total): " << total_fragments
        << "\n"
        << "Sum of count column for fragments overlapping peaks: "
        << total_peak_fragments << "\n"
        << "Sum of count column over peak-resident cut-site endpoints: "
        << total_peak_cutsites << "\n";
    return 0;
}

}  // namespace

int RunAtacEvidenceFromBinary(
    const std::string& sidecar_path,
    const std::string& peaks_path,
    const std::string& out_path,
    AtacEvidenceBarcodeDecoderFn decoder,
    void* decoder_ctx,
    std::ostream* err_in) {
    std::ostream& err = (err_in != nullptr) ? *err_in : std::cerr;

    if (sidecar_path.empty() || peaks_path.empty() || out_path.empty()) {
        err << "RunAtacEvidenceFromBinary: sidecar_path, peaks_path and out_path are required.\n";
        return 1;
    }
    if (decoder == nullptr) {
        err << "RunAtacEvidenceFromBinary: barcode decoder is required.\n";
        return 1;
    }

    FILE* fp = std::fopen(sidecar_path.c_str(), "rb");
    if (fp == nullptr) {
        err << "Failed to open binary sidecar: " << sidecar_path << "\n";
        return 1;
    }
    AtacEvidenceBinaryHeader header;
    if (std::fread(&header, sizeof(header), 1, fp) != 1) {
        err << "Failed to read binary sidecar header: " << sidecar_path << "\n";
        std::fclose(fp);
        return 1;
    }
    if (header.magic[0] != 'A' || header.magic[1] != 'E' ||
        header.magic[2] != 'V' || header.magic[3] != '1' ||
        header.format_version != 1 ||
        header.record_size != sizeof(AtacEvidenceBinaryRecord) ||
        header.flags != 0) {
        err << "Invalid ATAC evidence binary header: " << sidecar_path << "\n";
        std::fclose(fp);
        return 1;
    }
    if (fseeko(fp, 0, SEEK_END) != 0) {
        err << "Failed to seek binary sidecar: " << sidecar_path << "\n";
        std::fclose(fp);
        return 1;
    }
    const off_t file_size = ftello(fp);
    if (file_size < static_cast<off_t>(sizeof(AtacEvidenceBinaryHeader)) ||
        (static_cast<uint64_t>(file_size) - sizeof(AtacEvidenceBinaryHeader)) %
            sizeof(AtacEvidenceBinaryRecord) != 0) {
        err << "Binary sidecar size is not header + N fixed records: "
            << sidecar_path << "\n";
        std::fclose(fp);
        return 1;
    }
    const uint64_t records_by_size =
        (static_cast<uint64_t>(file_size) - sizeof(AtacEvidenceBinaryHeader)) /
        sizeof(AtacEvidenceBinaryRecord);
    const uint64_t records_to_read =
        (header.num_records == 0) ? records_by_size : header.num_records;
    if (records_to_read != records_by_size) {
        err << "Binary sidecar header num_records (" << header.num_records
            << ") does not match file size records (" << records_by_size
            << "): " << sidecar_path << "\n";
        std::fclose(fp);
        return 1;
    }
    if (fseeko(fp, sizeof(AtacEvidenceBinaryHeader), SEEK_SET) != 0) {
        err << "Failed to seek to binary sidecar records: " << sidecar_path
            << "\n";
        std::fclose(fp);
        return 1;
    }

    std::vector<std::string> chrom_names;
    if (!load_chroms(sidecar_path + ".chroms.tsv", &chrom_names, err)) {
        std::fclose(fp);
        return 1;
    }
    if (chrom_names.size() != header.num_chroms) {
        err << "Chrom metadata row count (" << chrom_names.size()
            << ") does not match binary header num_chroms ("
            << header.num_chroms << ")\n";
        std::fclose(fp);
        return 1;
    }

    std::unordered_map<std::string, std::vector<PeakInterval>> peaks_by_chrom;
    if (!load_peaks(peaks_path, &peaks_by_chrom, err)) {
        std::fclose(fp);
        return 1;
    }
    std::vector<const std::vector<PeakInterval>*> peaks_by_chrom_id(
        chrom_names.size(), nullptr);
    for (size_t chrom_id = 0; chrom_id < chrom_names.size(); ++chrom_id) {
        const auto it = peaks_by_chrom.find(chrom_names[chrom_id]);
        if (it != peaks_by_chrom.end()) {
            peaks_by_chrom_id[chrom_id] = &it->second;
        }
    }

    std::unordered_map<uint64_t, BarcodeCounts> per_barcode_key;
    uint64_t total_rows = 0;
    uint64_t total_fragments = 0;
    uint64_t total_peak_fragments = 0;
    uint64_t total_peak_cutsites = 0;

    const size_t chunk_records = 1 << 14;
    std::vector<AtacEvidenceBinaryRecord> records(chunk_records);
    uint64_t records_read = 0;
    while (records_read < records_to_read) {
        const size_t want = static_cast<size_t>(
            std::min<uint64_t>(chunk_records, records_to_read - records_read));
        const size_t got =
            std::fread(records.data(), sizeof(AtacEvidenceBinaryRecord), want, fp);
        if (got != want) {
            err << "Short read in binary sidecar records: " << sidecar_path
                << "\n";
            std::fclose(fp);
            return 1;
        }
        records_read += got;
        for (size_t i = 0; i < got; ++i) {
            const AtacEvidenceBinaryRecord& row = records[i];
            if (row.chrom_id < 0 ||
                static_cast<size_t>(row.chrom_id) >= chrom_names.size()) {
                err << "RunAtacEvidenceFromBinary: chrom_id out of range: "
                    << row.chrom_id << " (chrom_names=" << chrom_names.size()
                    << ")\n";
                std::fclose(fp);
                return 1;
            }
            const int64_t start = row.start;
            const int64_t end = row.end;
            const uint64_t count = row.count;
            if (end <= start || count == 0) {
                continue;
            }
            BarcodeCounts& bc_counts = per_barcode_key[row.barcode_key];
            bc_counts.atac_fragments += count;
            total_fragments += count;
            ++total_rows;
            const std::vector<PeakInterval>* peaks =
                peaks_by_chrom_id[static_cast<size_t>(row.chrom_id)];
            if (peaks == nullptr) {
                continue;
            }
            if (fragment_overlaps_any_peak(*peaks, start, end)) {
                bc_counts.peak_fragments += count;
                total_peak_fragments += count;
            }
            if (position_in_peak(*peaks, start)) {
                bc_counts.peak_cutsites += count;
                total_peak_cutsites += count;
            }
            if (position_in_peak(*peaks, end - 1)) {
                bc_counts.peak_cutsites += count;
                total_peak_cutsites += count;
            }
        }
    }
    if (std::fclose(fp) != 0) {
        err << "Failed to close binary sidecar: " << sidecar_path << "\n";
        return 1;
    }

    std::map<std::string, BarcodeCounts> sorted;
    for (const auto& kv : per_barcode_key) {
        std::string barcode;
        if (!decoder(kv.first, header.barcode_length, decoder_ctx, &barcode)) {
            err << "Failed to decode barcode key: " << kv.first << "\n";
            return 1;
        }
        BarcodeCounts& out_counts = sorted[barcode];
        out_counts.atac_fragments += kv.second.atac_fragments;
        out_counts.peak_fragments += kv.second.peak_fragments;
        out_counts.peak_cutsites += kv.second.peak_cutsites;
    }

    return write_evidence_output(sorted, out_path, total_rows, total_fragments,
                                 total_peak_fragments, total_peak_cutsites,
                                 err);
}

}  // namespace atac
}  // namespace libscrna
