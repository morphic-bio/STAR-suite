#include "multiome_atac_peak_mex.h"

#include <zlib.h>

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "libmacs3/macs3_frag_peak_pipeline.h"

namespace star {
namespace multiome {
namespace {

struct Args {
  std::string fragments;
  std::string sidecar;
  std::string peaks;
  std::string summits_out;
  std::string out_dir;
  std::string metrics_tsv;
  std::string barcode_translate;
  bool barcode_translate_from_first = false;
  bool call_peaks_from_sidecar = false;
  bool force = false;
  uint64_t max_barcodes = 0;
  int threads = 1;
  double macs3_pvalue = 1e-5;
  int macs3_min_length = 200;
  int macs3_max_gap = 30;
  bool macs3_uint8_counts = true;
  std::string temp_dir;
  std::string keep_intermediates_dir;
};

struct Peak {
  std::string chrom;
  int64_t start = 0;
  int64_t end = 0;
  std::string id;
};

struct PeakRef {
  int64_t start = 0;
  int64_t end = 0;
  uint32_t index = 0;
};

struct ChromPeakIndex {
  std::vector<PeakRef> peaks;
  int64_t max_len = 0;
};

struct Fragment {
  std::string chrom;
  int64_t start = 0;
  int64_t end = 0;
  std::string barcode;
  uint64_t count = 1;
};

struct BarcodeMetrics {
  uint64_t total_fragments = 0;
  uint64_t peak_fragments = 0;
};

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

bool parse_int64_range(const char *begin, const char *end, int64_t *out) {
  if (out == nullptr || begin == end) {
    return false;
  }
  errno = 0;
  char *parse_end = nullptr;
  const long long parsed = std::strtoll(begin, &parse_end, 10);
  if (errno != 0 || parse_end != end) {
    return false;
  }
  *out = static_cast<int64_t>(parsed);
  return true;
}

bool parse_uint64_range(const char *begin, const char *end, uint64_t *out) {
  if (out == nullptr || begin == end) {
    return false;
  }
  errno = 0;
  char *parse_end = nullptr;
  const unsigned long long parsed = std::strtoull(begin, &parse_end, 10);
  if (errno != 0 || parse_end != end) {
    return false;
  }
  *out = static_cast<uint64_t>(parsed);
  return true;
}

bool file_exists(const std::string &path) {
  struct stat st;
  return !path.empty() && stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

bool is_gzip_file(const std::string &path) {
  std::ifstream input(path.c_str(), std::ios::binary);
  if (!input.is_open()) {
    return false;
  }
  unsigned char magic[2] = {0, 0};
  input.read(reinterpret_cast<char *>(magic), 2);
  return input.gcount() == 2 && magic[0] == 0x1f && magic[1] == 0x8b;
}

class LineReader {
 public:
  explicit LineReader(const std::string &path)
      : path_(path), gzip_(is_gzip_file(path)) {}

  bool open(std::ostream &err) {
    if (gzip_) {
      gz_ = gzopen(path_.c_str(), "rb");
      if (gz_ == nullptr) {
        err << "Failed to open gzip file: " << path_ << "\n";
        return false;
      }
    } else {
      in_.open(path_.c_str());
      if (!in_.is_open()) {
        err << "Failed to open file: " << path_ << "\n";
        return false;
      }
    }
    return true;
  }

  bool getline(std::string *line) {
    if (line == nullptr) {
      return false;
    }
    line->clear();
    if (!gzip_) {
      return static_cast<bool>(std::getline(in_, *line));
    }

    char buffer[65536];
    while (true) {
      char *got = gzgets(gz_, buffer, static_cast<int>(sizeof(buffer)));
      if (got == nullptr) {
        return !line->empty();
      }
      line->append(buffer);
      if (!line->empty() && line->back() == '\n') {
        break;
      }
    }
    if (!line->empty() && line->back() == '\n') {
      line->pop_back();
    }
    return true;
  }

  ~LineReader() {
    if (gz_ != nullptr) {
      gzclose(gz_);
    }
  }

 private:
  std::string path_;
  bool gzip_ = false;
  gzFile gz_ = nullptr;
  std::ifstream in_;
};

void trim_cr(std::string *line) {
  if (line != nullptr && !line->empty() && line->back() == '\r') {
    line->pop_back();
  }
}

std::vector<std::string> split_tab(const std::string &line) {
  std::vector<std::string> fields;
  size_t begin = 0;
  while (begin <= line.size()) {
    const size_t pos = line.find('\t', begin);
    if (pos == std::string::npos) {
      fields.push_back(line.substr(begin));
      break;
    }
    fields.push_back(line.substr(begin, pos - begin));
    begin = pos + 1;
  }
  return fields;
}

std::vector<std::string> split_translate_line(const std::string &line) {
  const size_t tab = line.find('\t');
  const size_t comma = line.find(',');
  size_t sep = std::string::npos;
  if (tab == std::string::npos) {
    sep = comma;
  } else if (comma == std::string::npos) {
    sep = tab;
  } else {
    sep = std::min(tab, comma);
  }
  if (sep == std::string::npos) {
    return std::vector<std::string>();
  }
  std::vector<std::string> fields;
  fields.push_back(line.substr(0, sep));
  fields.push_back(line.substr(sep + 1));
  return fields;
}

bool parse_fragment(const std::string &line, Fragment *fragment,
                    std::ostream &err, uint64_t line_no) {
  if (fragment == nullptr) {
    return false;
  }
  const size_t t1 = line.find('\t');
  const size_t t2 = (t1 == std::string::npos) ? std::string::npos
                                               : line.find('\t', t1 + 1);
  const size_t t3 = (t2 == std::string::npos) ? std::string::npos
                                               : line.find('\t', t2 + 1);
  if (t1 == std::string::npos || t2 == std::string::npos ||
      t3 == std::string::npos) {
    err << "Fragment line " << line_no << " has fewer than 4 columns\n";
    return false;
  }
  const size_t t4 = line.find('\t', t3 + 1);

  int64_t start = 0;
  int64_t end = 0;
  if (!parse_int64_range(line.data() + t1 + 1, line.data() + t2, &start) ||
      !parse_int64_range(line.data() + t2 + 1, line.data() + t3, &end)) {
    err << "Fragment line " << line_no << " has non-numeric start/end\n";
    return false;
  }

  uint64_t count = 1;
  if (t4 != std::string::npos && t4 + 1 < line.size()) {
    const size_t t5 = line.find('\t', t4 + 1);
    const char *count_begin = line.data() + t4 + 1;
    const char *count_end = (t5 == std::string::npos) ? line.data() + line.size()
                                                       : line.data() + t5;
    if (!parse_uint64_range(count_begin, count_end, &count)) {
      err << "Fragment line " << line_no << " has non-numeric count\n";
      return false;
    }
  }

  fragment->chrom.assign(line.data(), t1);
  fragment->start = start;
  fragment->end = end;
  fragment->barcode =
      (t4 == std::string::npos) ? line.substr(t3 + 1)
                                : line.substr(t3 + 1, t4 - t3 - 1);
  fragment->count = count;
  return true;
}

bool load_peaks(const std::string &path, std::vector<Peak> *peaks,
                std::unordered_map<std::string, ChromPeakIndex> *index,
                std::ostream &err) {
  if (peaks == nullptr || index == nullptr) {
    return false;
  }
  LineReader reader(path);
  if (!reader.open(err)) {
    return false;
  }

  std::string line;
  uint64_t line_no = 0;
  while (reader.getline(&line)) {
    ++line_no;
    trim_cr(&line);
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (line.compare(0, 5, "track") == 0 ||
        line.compare(0, 7, "browser") == 0) {
      continue;
    }
    const std::vector<std::string> fields = split_tab(line);
    if (fields.size() < 3) {
      err << "Peak line " << line_no << " has fewer than 3 columns\n";
      return false;
    }

    int64_t start = 0;
    int64_t end = 0;
    if (!parse_int64_range(fields[1].data(), fields[1].data() + fields[1].size(),
                           &start) ||
        !parse_int64_range(fields[2].data(), fields[2].data() + fields[2].size(),
                           &end)) {
      err << "Peak line " << line_no << " has non-numeric start/end\n";
      return false;
    }
    if (end <= start) {
      continue;
    }

    Peak peak;
    peak.chrom = fields[0];
    peak.start = start;
    peak.end = end;
    std::ostringstream peak_id;
    peak_id << peak.chrom << ":" << peak.start << "-" << peak.end;
    peak.id = peak_id.str();
    const uint32_t peak_index = static_cast<uint32_t>(peaks->size());
    peaks->push_back(peak);

    PeakRef ref;
    ref.start = start;
    ref.end = end;
    ref.index = peak_index;
    ChromPeakIndex &chrom_index = (*index)[peak.chrom];
    chrom_index.peaks.push_back(ref);
    const int64_t len = end - start;
    if (len > chrom_index.max_len) {
      chrom_index.max_len = len;
    }
  }

  for (auto &kv : *index) {
    std::sort(kv.second.peaks.begin(), kv.second.peaks.end(),
              [](const PeakRef &a, const PeakRef &b) {
                if (a.start != b.start) {
                  return a.start < b.start;
                }
                if (a.end != b.end) {
                  return a.end < b.end;
                }
                return a.index < b.index;
              });
  }

  if (peaks->empty()) {
    err << "No peaks found in " << path << "\n";
    return false;
  }
  return true;
}

bool validate_sidecar_header(FILE *fp, const std::string &path,
                             AtacEvidenceBinaryHeader *header,
                             uint64_t *records_to_read,
                             std::ostream &err) {
  if (fp == nullptr || header == nullptr || records_to_read == nullptr) {
    return false;
  }
  if (std::fread(header, sizeof(*header), 1, fp) != 1) {
    err << "Failed to read binary sidecar header: " << path << "\n";
    return false;
  }
  if (header->magic[0] != 'A' || header->magic[1] != 'E' ||
      header->magic[2] != 'V' || header->magic[3] != '1' ||
      header->format_version != 1 ||
      header->record_size != sizeof(AtacEvidenceBinaryRecord) ||
      header->flags != 0 || header->barcode_length == 0 ||
      header->barcode_length > 32) {
    err << "Invalid ATAC evidence binary header: " << path << "\n";
    return false;
  }
  if (fseeko(fp, 0, SEEK_END) != 0) {
    err << "Failed to seek binary sidecar: " << path << "\n";
    return false;
  }
  const off_t file_size = ftello(fp);
  if (file_size < static_cast<off_t>(sizeof(AtacEvidenceBinaryHeader)) ||
      (static_cast<uint64_t>(file_size) - sizeof(AtacEvidenceBinaryHeader)) %
          sizeof(AtacEvidenceBinaryRecord) != 0) {
    err << "Binary sidecar size is not header + N fixed records: " << path
        << "\n";
    return false;
  }
  const uint64_t records_by_size =
      (static_cast<uint64_t>(file_size) - sizeof(AtacEvidenceBinaryHeader)) /
      sizeof(AtacEvidenceBinaryRecord);
  *records_to_read =
      (header->num_records == 0) ? records_by_size : header->num_records;
  if (*records_to_read != records_by_size) {
    err << "Binary sidecar header num_records (" << header->num_records
        << ") does not match file size records (" << records_by_size
        << "): " << path << "\n";
    return false;
  }
  if (fseeko(fp, sizeof(AtacEvidenceBinaryHeader), SEEK_SET) != 0) {
    err << "Failed to seek to binary sidecar records: " << path << "\n";
    return false;
  }
  return true;
}

bool open_sidecar_records(const std::string &path, FILE **fp_out,
                          AtacEvidenceBinaryHeader *header,
                          uint64_t *records_to_read, std::ostream &err) {
  FILE *fp = std::fopen(path.c_str(), "rb");
  if (fp == nullptr) {
    err << "Failed to open binary sidecar: " << path << "\n";
    return false;
  }
  if (!validate_sidecar_header(fp, path, header, records_to_read, err)) {
    std::fclose(fp);
    return false;
  }
  *fp_out = fp;
  return true;
}

bool parse_chrom_id(const std::string &value, int64_t *out) {
  if (value.empty()) {
    return false;
  }
  char *end = nullptr;
  const long long parsed = std::strtoll(value.c_str(), &end, 10);
  if (end == value.c_str() || *end != '\0') {
    return false;
  }
  *out = static_cast<int64_t>(parsed);
  return true;
}

bool load_chroms(const std::string &path, std::vector<std::string> *chrom_names,
                 std::ostream &err) {
  if (chrom_names == nullptr) {
    return false;
  }
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
    const std::vector<std::string> fields = split_tab(line);
    if (fields.size() < 2) {
      err << "Chrom metadata line " << line_no
          << " has fewer than 2 columns\n";
      return false;
    }
    int64_t chrom_id = -1;
    if (!parse_chrom_id(fields[0], &chrom_id) || chrom_id < 0) {
      err << "Chrom metadata line " << line_no << " has invalid chrom_id\n";
      return false;
    }
    if (static_cast<size_t>(chrom_id) >= chrom_names->size()) {
      chrom_names->resize(static_cast<size_t>(chrom_id) + 1);
    }
    (*chrom_names)[static_cast<size_t>(chrom_id)] = fields[1];
  }
  return true;
}

bool load_sidecar_info(const std::string &path,
                       AtacEvidenceBinaryHeader *header,
                       uint64_t *records_to_read,
                       std::vector<std::string> *chrom_names,
                       std::ostream &err) {
  FILE *fp = nullptr;
  if (!open_sidecar_records(path, &fp, header, records_to_read, err)) {
    return false;
  }
  std::fclose(fp);
  if (!load_chroms(path + ".chroms.tsv", chrom_names, err)) {
    return false;
  }
  if (chrom_names->size() != header->num_chroms) {
    err << "Chrom metadata row count (" << chrom_names->size()
        << ") does not match binary header num_chroms ("
        << header->num_chroms << ")\n";
    return false;
  }
  return true;
}

class SidecarFragmentIterator final : public macs3::FragmentIterator {
 public:
  explicit SidecarFragmentIterator(const std::string &path) : path_(path) {
    std::ostringstream err;
    if (!open_sidecar_records(path_, &fp_, &header_, &records_to_read_, err)) {
      error_ = err.str();
      return;
    }
    if (!load_chroms(path_ + ".chroms.tsv", &chrom_names_, err)) {
      error_ = err.str();
      return;
    }
    if (chrom_names_.size() != header_.num_chroms) {
      std::ostringstream msg;
      msg << "Chrom metadata row count (" << chrom_names_.size()
          << ") does not match binary header num_chroms ("
          << header_.num_chroms << ")";
      error_ = msg.str();
      return;
    }
  }

  ~SidecarFragmentIterator() override {
    if (fp_ != nullptr) {
      std::fclose(fp_);
    }
  }

  bool Next(macs3::FragmentRecord *out) override {
    if (out == nullptr || fp_ == nullptr || !error_.empty() ||
        records_read_ >= records_to_read_) {
      return false;
    }
    AtacEvidenceBinaryRecord row;
    if (std::fread(&row, sizeof(row), 1, fp_) != 1) {
      error_ = "Short read in binary sidecar records: " + path_;
      return false;
    }
    ++records_read_;
    if (row.chrom_id < 0 ||
        static_cast<size_t>(row.chrom_id) >= chrom_names_.size()) {
      std::ostringstream msg;
      msg << "Sidecar chrom_id out of range: " << row.chrom_id;
      error_ = msg.str();
      return false;
    }
    if (row.end <= row.start || row.count == 0) {
      return Next(out);
    }
    out->chrom_id = row.chrom_id;
    out->start = row.start;
    out->end = row.end;
    out->count = row.count;
    return true;
  }

  const std::vector<std::string> &ChromNames() const override {
    return chrom_names_;
  }

  const std::string &Error() const override { return error_; }

 private:
  std::string path_;
  FILE *fp_ = nullptr;
  AtacEvidenceBinaryHeader header_{};
  uint64_t records_to_read_ = 0;
  uint64_t records_read_ = 0;
  std::vector<std::string> chrom_names_;
  std::string error_;
};

std::string decode_barcode_key(uint64_t key, uint32_t length) {
  static const char bases[4] = {'A', 'C', 'G', 'T'};
  std::string out;
  out.reserve(length);
  for (uint32_t i = 0; i < length; ++i) {
    const uint32_t shift = (length - 1 - i) * 2;
    out.push_back(bases[(key >> shift) & 0x3u]);
  }
  return out;
}

bool load_translate_table(
    const std::string &path, bool from_first,
    std::unordered_map<std::string, std::string> *translate,
    std::ostream &err) {
  if (translate == nullptr || path.empty()) {
    return true;
  }
  LineReader reader(path);
  if (!reader.open(err)) {
    return false;
  }
  std::string line;
  uint64_t line_no = 0;
  while (reader.getline(&line)) {
    ++line_no;
    trim_cr(&line);
    if (line.empty() || line[0] == '#') {
      continue;
    }
    const std::vector<std::string> fields = split_translate_line(line);
    if (fields.size() < 2) {
      continue;
    }
    const std::string source = from_first ? fields[0] : fields[1];
    const std::string dest = from_first ? fields[1] : fields[0];
    if (source.empty() || dest.empty()) {
      err << "Barcode translate line " << line_no
          << " has an empty source or destination barcode\n";
      return false;
    }
    (*translate)[source] = dest;
  }
  return true;
}

bool resolve_barcode_name(
    uint64_t key, uint32_t barcode_length,
    const std::unordered_map<std::string, std::string> &translate,
    std::string *out, std::ostream &err) {
  if (out == nullptr) {
    return false;
  }
  const std::string source = decode_barcode_key(key, barcode_length);
  if (translate.empty()) {
    *out = source;
    return true;
  }
  const auto it = translate.find(source);
  if (it == translate.end()) {
    err << "Barcode " << source
        << " was not found in the translation table\n";
    return false;
  }
  *out = it->second;
  return true;
}

std::string parent_dir(const std::string &path);
bool mkdir_p(const std::string &path, std::ostream &err);

bool call_peaks_from_sidecar_if_needed(const Args &args, std::ostream &err) {
  if (!args.call_peaks_from_sidecar) {
    return true;
  }
  if (!args.force && file_exists(args.peaks) && file_exists(args.summits_out)) {
    err << "Reusing existing sidecar-derived peaks: " << args.peaks << "\n";
    return true;
  }

  SidecarFragmentIterator iter(args.sidecar);
  if (!iter.Error().empty()) {
    err << iter.Error() << "\n";
    return false;
  }

  if (!mkdir_p(parent_dir(args.peaks), err) ||
      !mkdir_p(parent_dir(args.summits_out), err)) {
    return false;
  }

  chromap::peaks::Macs3FragPeakPipelineParams pr;
  pr.bdgpeakcall_cutoff =
      chromap::peaks::BdgPeakCallCutoffFromPValue(args.macs3_pvalue);
  if (pr.bdgpeakcall_cutoff <= 0.f) {
    err << "Invalid MACS3 FRAG p-value for bdgpeakcall cutoff\n";
    return false;
  }
  pr.min_length = args.macs3_min_length;
  pr.max_gap = args.macs3_max_gap;
  pr.macs3_uint8_counts = args.macs3_uint8_counts;
  pr.peak_caller_threads = args.threads;

  std::string work_used;
  std::string peak_error;
  if (!chromap::peaks::RunMacs3FragPeakPipelineFromSortedIterator(
          iter, pr, chromap::peaks::Macs3FragPeakPipelinePaths(), args.peaks,
          args.summits_out, args.keep_intermediates_dir, args.temp_dir,
          &work_used, &peak_error)) {
    err << "MACS3 FRAG peaks from binary sidecar failed: " << peak_error
        << "\n";
    return false;
  }
  err << "Wrote sidecar-derived peaks: " << args.peaks << "\n";
  err << "Wrote sidecar-derived summits: " << args.summits_out << "\n";
  return true;
}

bool collect_total_fragments_text(
    const std::string &fragments_path,
    std::unordered_map<std::string, BarcodeMetrics> *metrics,
    std::ostream &err) {
  if (metrics == nullptr) {
    return false;
  }
  LineReader reader(fragments_path);
  if (!reader.open(err)) {
    return false;
  }
  std::string line;
  uint64_t line_no = 0;
  Fragment fragment;
  while (reader.getline(&line)) {
    ++line_no;
    trim_cr(&line);
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (!parse_fragment(line, &fragment, err, line_no)) {
      return false;
    }
    (*metrics)[fragment.barcode].total_fragments += fragment.count;
  }
  return true;
}

bool collect_total_fragments_sidecar(
    const std::string &sidecar_path,
    std::unordered_map<uint64_t, BarcodeMetrics> *metrics,
    std::ostream &err) {
  if (metrics == nullptr) {
    return false;
  }
  FILE *fp = nullptr;
  AtacEvidenceBinaryHeader header{};
  uint64_t records_to_read = 0;
  if (!open_sidecar_records(sidecar_path, &fp, &header, &records_to_read, err)) {
    return false;
  }
  const size_t chunk_records = 1 << 14;
  std::vector<AtacEvidenceBinaryRecord> records(chunk_records);
  uint64_t records_read = 0;
  while (records_read < records_to_read) {
    const size_t want = static_cast<size_t>(
        std::min<uint64_t>(chunk_records, records_to_read - records_read));
    const size_t got =
        std::fread(records.data(), sizeof(AtacEvidenceBinaryRecord), want, fp);
    if (got != want) {
      err << "Short read in binary sidecar records: " << sidecar_path << "\n";
      std::fclose(fp);
      return false;
    }
    records_read += got;
    for (size_t i = 0; i < got; ++i) {
      const AtacEvidenceBinaryRecord &row = records[i];
      if (row.end <= row.start || row.count == 0) {
        continue;
      }
      (*metrics)[row.barcode_key].total_fragments += row.count;
    }
  }
  std::fclose(fp);
  return true;
}

std::unordered_set<std::string> select_barcodes_text(
    const std::unordered_map<std::string, BarcodeMetrics> &metrics,
    uint64_t max_barcodes) {
  std::unordered_set<std::string> keep;
  if (max_barcodes == 0 || metrics.size() <= max_barcodes) {
    keep.reserve(metrics.size());
    for (const auto &kv : metrics) {
      keep.insert(kv.first);
    }
    return keep;
  }

  std::vector<std::pair<std::string, uint64_t> > ranked;
  ranked.reserve(metrics.size());
  for (const auto &kv : metrics) {
    ranked.push_back(std::make_pair(kv.first, kv.second.total_fragments));
  }
  std::sort(ranked.begin(), ranked.end(),
            [](const std::pair<std::string, uint64_t> &a,
               const std::pair<std::string, uint64_t> &b) {
              if (a.second != b.second) {
                return a.second > b.second;
              }
              return a.first < b.first;
            });
  keep.reserve(static_cast<size_t>(max_barcodes));
  for (uint64_t i = 0; i < max_barcodes && i < ranked.size(); ++i) {
    keep.insert(ranked[static_cast<size_t>(i)].first);
  }
  return keep;
}

std::unordered_set<uint64_t> select_barcodes_sidecar(
    const std::unordered_map<uint64_t, BarcodeMetrics> &metrics,
    uint64_t max_barcodes) {
  std::unordered_set<uint64_t> keep;
  if (max_barcodes == 0 || metrics.size() <= max_barcodes) {
    keep.reserve(metrics.size());
    for (const auto &kv : metrics) {
      keep.insert(kv.first);
    }
    return keep;
  }

  std::vector<std::pair<uint64_t, uint64_t> > ranked;
  ranked.reserve(metrics.size());
  for (const auto &kv : metrics) {
    ranked.push_back(std::make_pair(kv.first, kv.second.total_fragments));
  }
  std::sort(ranked.begin(), ranked.end(),
            [](const std::pair<uint64_t, uint64_t> &a,
               const std::pair<uint64_t, uint64_t> &b) {
              if (a.second != b.second) {
                return a.second > b.second;
              }
              return a.first < b.first;
            });
  keep.reserve(static_cast<size_t>(max_barcodes));
  for (uint64_t i = 0; i < max_barcodes && i < ranked.size(); ++i) {
    keep.insert(ranked[static_cast<size_t>(i)].first);
  }
  return keep;
}

template <typename Fn>
void for_each_overlap(const ChromPeakIndex &index, int64_t start, int64_t end,
                      Fn fn) {
  if (end <= start || index.peaks.empty()) {
    return;
  }
  const int64_t lower_start = start - index.max_len;
  const auto lower =
      std::lower_bound(index.peaks.begin(), index.peaks.end(), lower_start,
                       [](const PeakRef &peak, int64_t value) {
                         return peak.start < value;
                       });
  for (auto it = lower; it != index.peaks.end() && it->start < end; ++it) {
    if (it->end > start) {
      fn(*it);
    }
  }
}

uint32_t barcode_index_for_text(
    const std::string &barcode, std::unordered_map<std::string, uint32_t> *lookup,
    std::vector<std::string> *barcodes) {
  auto found = lookup->find(barcode);
  if (found != lookup->end()) {
    return found->second;
  }
  const uint32_t index = static_cast<uint32_t>(barcodes->size());
  (*lookup)[barcode] = index;
  barcodes->push_back(barcode);
  return index;
}

uint32_t barcode_index_for_key(
    uint64_t barcode_key, std::unordered_map<uint64_t, uint32_t> *lookup,
    std::vector<uint64_t> *barcode_keys) {
  auto found = lookup->find(barcode_key);
  if (found != lookup->end()) {
    return found->second;
  }
  const uint32_t index = static_cast<uint32_t>(barcode_keys->size());
  (*lookup)[barcode_key] = index;
  barcode_keys->push_back(barcode_key);
  return index;
}

bool build_matrix_text(
    const std::string &fragments_path,
    const std::unordered_map<std::string, ChromPeakIndex> &peaks_by_chrom,
    const std::unordered_set<std::string> &keep_barcodes,
    std::unordered_map<std::string, BarcodeMetrics> *metrics,
    std::unordered_map<uint64_t, uint64_t> *entries,
    std::vector<std::string> *barcodes, std::ostream &err) {
  if (metrics == nullptr || entries == nullptr || barcodes == nullptr) {
    return false;
  }
  LineReader reader(fragments_path);
  if (!reader.open(err)) {
    return false;
  }

  std::unordered_map<std::string, uint32_t> barcode_to_index;
  std::string line;
  uint64_t line_no = 0;
  Fragment fragment;
  while (reader.getline(&line)) {
    ++line_no;
    trim_cr(&line);
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (!parse_fragment(line, &fragment, err, line_no)) {
      return false;
    }
    if (keep_barcodes.find(fragment.barcode) == keep_barcodes.end()) {
      continue;
    }
    const auto chrom_it = peaks_by_chrom.find(fragment.chrom);
    if (chrom_it == peaks_by_chrom.end()) {
      continue;
    }

    bool overlapped_any = false;
    uint32_t barcode_index = 0;
    for_each_overlap(chrom_it->second, fragment.start, fragment.end,
                     [&](const PeakRef &peak) {
                       if (!overlapped_any) {
                         barcode_index = barcode_index_for_text(
                             fragment.barcode, &barcode_to_index, barcodes);
                         overlapped_any = true;
                       }
                       const uint64_t key =
                           (static_cast<uint64_t>(peak.index) << 32) |
                           static_cast<uint64_t>(barcode_index);
                       (*entries)[key] += fragment.count;
                     });

    if (overlapped_any) {
      (*metrics)[fragment.barcode].peak_fragments += fragment.count;
    }
  }
  return true;
}

bool build_matrix_sidecar(
    const std::string &sidecar_path,
    const std::vector<std::string> &chrom_names,
    const std::unordered_map<std::string, ChromPeakIndex> &peaks_by_chrom,
    const std::unordered_set<uint64_t> &keep_barcodes,
    std::unordered_map<uint64_t, BarcodeMetrics> *metrics,
    std::unordered_map<uint64_t, uint64_t> *entries,
    std::vector<uint64_t> *barcode_keys, std::ostream &err) {
  if (metrics == nullptr || entries == nullptr || barcode_keys == nullptr) {
    return false;
  }
  std::vector<const ChromPeakIndex *> peaks_by_chrom_id(chrom_names.size(),
                                                        nullptr);
  for (size_t i = 0; i < chrom_names.size(); ++i) {
    const auto it = peaks_by_chrom.find(chrom_names[i]);
    if (it != peaks_by_chrom.end()) {
      peaks_by_chrom_id[i] = &it->second;
    }
  }

  FILE *fp = nullptr;
  AtacEvidenceBinaryHeader header{};
  uint64_t records_to_read = 0;
  if (!open_sidecar_records(sidecar_path, &fp, &header, &records_to_read, err)) {
    return false;
  }
  std::unordered_map<uint64_t, uint32_t> barcode_to_index;
  const size_t chunk_records = 1 << 14;
  std::vector<AtacEvidenceBinaryRecord> records(chunk_records);
  uint64_t records_read = 0;
  while (records_read < records_to_read) {
    const size_t want = static_cast<size_t>(
        std::min<uint64_t>(chunk_records, records_to_read - records_read));
    const size_t got =
        std::fread(records.data(), sizeof(AtacEvidenceBinaryRecord), want, fp);
    if (got != want) {
      err << "Short read in binary sidecar records: " << sidecar_path << "\n";
      std::fclose(fp);
      return false;
    }
    records_read += got;
    for (size_t i = 0; i < got; ++i) {
      const AtacEvidenceBinaryRecord &row = records[i];
      if (row.end <= row.start || row.count == 0) {
        continue;
      }
      if (keep_barcodes.find(row.barcode_key) == keep_barcodes.end()) {
        continue;
      }
      if (row.chrom_id < 0 ||
          static_cast<size_t>(row.chrom_id) >= peaks_by_chrom_id.size()) {
        err << "Sidecar chrom_id out of range: " << row.chrom_id << "\n";
        std::fclose(fp);
        return false;
      }
      const ChromPeakIndex *chrom_index =
          peaks_by_chrom_id[static_cast<size_t>(row.chrom_id)];
      if (chrom_index == nullptr) {
        continue;
      }

      bool overlapped_any = false;
      uint32_t barcode_index = 0;
      for_each_overlap(*chrom_index, row.start, row.end,
                       [&](const PeakRef &peak) {
                         if (!overlapped_any) {
                           barcode_index = barcode_index_for_key(
                               row.barcode_key, &barcode_to_index,
                               barcode_keys);
                           overlapped_any = true;
                         }
                         const uint64_t key =
                             (static_cast<uint64_t>(peak.index) << 32) |
                             static_cast<uint64_t>(barcode_index);
                         (*entries)[key] += row.count;
                       });
      if (overlapped_any) {
        (*metrics)[row.barcode_key].peak_fragments += row.count;
      }
    }
  }
  std::fclose(fp);
  return true;
}

std::string parent_dir(const std::string &path) {
  const size_t pos = path.find_last_of('/');
  if (pos == std::string::npos) {
    return ".";
  }
  if (pos == 0) {
    return "/";
  }
  return path.substr(0, pos);
}

bool mkdir_p(const std::string &path, std::ostream &err) {
  if (path.empty()) {
    return true;
  }
  std::string current;
  size_t i = 0;
  if (path[0] == '/') {
    current = "/";
    i = 1;
  }
  while (i <= path.size()) {
    const size_t slash = path.find('/', i);
    const std::string part =
        path.substr(i, slash == std::string::npos ? std::string::npos
                                                  : slash - i);
    if (!part.empty()) {
      if (!current.empty() && current.back() != '/') {
        current.push_back('/');
      }
      current += part;
      if (mkdir(current.c_str(), 0775) != 0 && errno != EEXIST) {
        err << "Failed to create directory " << current << ": "
            << std::strerror(errno) << "\n";
        return false;
      }
    }
    if (slash == std::string::npos) {
      break;
    }
    i = slash + 1;
  }
  return true;
}

bool gz_write(gzFile out, const std::string &text) {
  return gzwrite(out, text.data(), static_cast<unsigned int>(text.size())) ==
         static_cast<int>(text.size());
}

bool write_gzip_lines(const std::string &path,
                      const std::vector<std::string> &lines,
                      std::ostream &err) {
  gzFile out = gzopen(path.c_str(), "wb");
  if (out == nullptr) {
    err << "Failed to open gzip output: " << path << "\n";
    return false;
  }
  for (const std::string &line : lines) {
    if (!gz_write(out, line)) {
      err << "Failed to write gzip output: " << path << "\n";
      gzclose(out);
      return false;
    }
  }
  gzclose(out);
  return true;
}

bool write_features(const std::string &path, const std::vector<Peak> &peaks,
                    std::ostream &err) {
  std::vector<std::string> lines;
  lines.reserve(peaks.size());
  for (const Peak &peak : peaks) {
    std::ostringstream out;
    out << peak.id << '\t' << peak.id << "\tPeaks\t" << peak.chrom << '\t'
        << peak.start << '\t' << peak.end << '\n';
    lines.push_back(out.str());
  }
  return write_gzip_lines(path, lines, err);
}

bool write_barcodes(const std::string &path,
                    const std::vector<std::string> &barcodes,
                    std::ostream &err) {
  std::vector<std::string> lines;
  lines.reserve(barcodes.size());
  for (const std::string &barcode : barcodes) {
    lines.push_back(barcode + "\n");
  }
  return write_gzip_lines(path, lines, err);
}

bool write_matrix(const std::string &path, size_t num_peaks,
                  size_t num_barcodes,
                  const std::unordered_map<uint64_t, uint64_t> &entries,
                  std::ostream &err) {
  gzFile out = gzopen(path.c_str(), "wb");
  if (out == nullptr) {
    err << "Failed to open gzip output: " << path << "\n";
    return false;
  }
  std::ostringstream header;
  header << "%%MatrixMarket matrix coordinate integer general\n"
         << "% generated by star_multiome_atac_peak_mex\n"
         << num_peaks << ' ' << num_barcodes << ' ' << entries.size() << '\n';
  if (!gz_write(out, header.str())) {
    err << "Failed to write matrix header\n";
    gzclose(out);
    return false;
  }

  std::vector<uint64_t> keys;
  keys.reserve(entries.size());
  for (const auto &kv : entries) {
    keys.push_back(kv.first);
  }
  std::sort(keys.begin(), keys.end());
  for (uint64_t key : keys) {
    const uint32_t peak_index = static_cast<uint32_t>(key >> 32);
    const uint32_t barcode_index = static_cast<uint32_t>(key & 0xffffffffULL);
    std::ostringstream row;
    row << static_cast<uint64_t>(peak_index) + 1 << ' '
        << static_cast<uint64_t>(barcode_index) + 1 << ' ' << entries.at(key)
        << '\n';
    if (!gz_write(out, row.str())) {
      err << "Failed to write matrix row\n";
      gzclose(out);
      return false;
    }
  }
  gzclose(out);
  return true;
}

bool write_metrics_text(
    const std::string &path, const std::vector<std::string> &barcodes,
    const std::unordered_map<std::string, BarcodeMetrics> &metrics,
    std::ostream &err) {
  if (!mkdir_p(parent_dir(path), err)) {
    return false;
  }
  std::ofstream out(path.c_str());
  if (!out.is_open()) {
    err << "Failed to open metrics output: " << path << "\n";
    return false;
  }
  out << "barcode\tatac_fragments\tatac_peak_region_fragments\t"
         "atac_peak_region_cutsites\tatac_peak_fraction\n";
  out << std::setprecision(8);
  for (const std::string &barcode : barcodes) {
    const auto it = metrics.find(barcode);
    const BarcodeMetrics m =
        (it == metrics.end()) ? BarcodeMetrics() : it->second;
    const double frac =
        (m.total_fragments == 0)
            ? 0.0
            : static_cast<double>(m.peak_fragments) /
                  static_cast<double>(m.total_fragments);
    out << barcode << '\t' << m.total_fragments << '\t' << m.peak_fragments
        << '\t' << (m.peak_fragments * 2) << '\t' << frac << '\n';
  }
  return true;
}

bool write_metrics_sidecar(
    const std::string &path, const std::vector<uint64_t> &barcode_keys,
    uint32_t barcode_length,
    const std::unordered_map<uint64_t, BarcodeMetrics> &metrics,
    const std::unordered_map<std::string, std::string> &translate,
    const std::vector<std::string> *barcode_names, std::ostream &err) {
  if (barcode_names == nullptr || barcode_names->size() != barcode_keys.size()) {
    return false;
  }
  (void)barcode_length;
  (void)translate;
  if (!mkdir_p(parent_dir(path), err)) {
    return false;
  }
  std::ofstream out(path.c_str());
  if (!out.is_open()) {
    err << "Failed to open metrics output: " << path << "\n";
    return false;
  }
  out << "barcode\tatac_fragments\tatac_peak_region_fragments\t"
         "atac_peak_region_cutsites\tatac_peak_fraction\n";
  out << std::setprecision(8);
  for (size_t i = 0; i < barcode_keys.size(); ++i) {
    const auto it = metrics.find(barcode_keys[i]);
    const BarcodeMetrics m =
        (it == metrics.end()) ? BarcodeMetrics() : it->second;
    const double frac =
        (m.total_fragments == 0)
            ? 0.0
            : static_cast<double>(m.peak_fragments) /
                  static_cast<double>(m.total_fragments);
    out << (*barcode_names)[i] << '\t' << m.total_fragments << '\t'
        << m.peak_fragments << '\t' << (m.peak_fragments * 2) << '\t'
        << frac << '\n';
  }
  return true;
}

bool write_core_outputs(const Args &args, const std::vector<Peak> &peaks,
                        const std::vector<std::string> &barcodes,
                        const std::unordered_map<uint64_t, uint64_t> &entries,
                        std::ostream &err) {
  if (!mkdir_p(args.out_dir, err)) {
    return false;
  }
  if (!write_features(args.out_dir + "/features.tsv.gz", peaks, err)) {
    return false;
  }
  if (!write_barcodes(args.out_dir + "/barcodes.tsv.gz", barcodes, err)) {
    return false;
  }
  return write_matrix(args.out_dir + "/matrix.mtx.gz", peaks.size(),
                      barcodes.size(), entries, err);
}

bool run_text_path(const Args &args, const std::vector<Peak> &peaks,
                   const std::unordered_map<std::string, ChromPeakIndex> &peaks_by_chrom,
                   std::ostream &err) {
  std::unordered_map<std::string, BarcodeMetrics> metrics;
  if (!collect_total_fragments_text(args.fragments, &metrics, err)) {
    return false;
  }
  const std::unordered_set<std::string> keep =
      select_barcodes_text(metrics, args.max_barcodes);

  std::unordered_map<uint64_t, uint64_t> entries;
  std::vector<std::string> barcodes;
  if (!build_matrix_text(args.fragments, peaks_by_chrom, keep, &metrics,
                         &entries, &barcodes, err)) {
    return false;
  }
  if (barcodes.empty()) {
    err << "No peak-overlapping barcodes were found\n";
    return false;
  }
  if (!write_core_outputs(args, peaks, barcodes, entries, err)) {
    return false;
  }
  return write_metrics_text(args.metrics_tsv, barcodes, metrics, err);
}

bool run_sidecar_path(
    const Args &args, const std::vector<Peak> &peaks,
    const std::unordered_map<std::string, ChromPeakIndex> &peaks_by_chrom,
    std::ostream &err) {
  AtacEvidenceBinaryHeader header{};
  uint64_t records_to_read = 0;
  std::vector<std::string> chrom_names;
  if (!load_sidecar_info(args.sidecar, &header, &records_to_read, &chrom_names,
                         err)) {
    return false;
  }
  (void)records_to_read;

  std::unordered_map<uint64_t, BarcodeMetrics> metrics;
  if (!collect_total_fragments_sidecar(args.sidecar, &metrics, err)) {
    return false;
  }
  const std::unordered_set<uint64_t> keep =
      select_barcodes_sidecar(metrics, args.max_barcodes);

  std::unordered_map<uint64_t, uint64_t> entries;
  std::vector<uint64_t> barcode_keys;
  if (!build_matrix_sidecar(args.sidecar, chrom_names, peaks_by_chrom, keep,
                            &metrics, &entries, &barcode_keys, err)) {
    return false;
  }
  if (barcode_keys.empty()) {
    err << "No peak-overlapping barcodes were found\n";
    return false;
  }

  std::unordered_map<std::string, std::string> translate;
  if (!args.barcode_translate.empty() &&
      !load_translate_table(args.barcode_translate,
                            args.barcode_translate_from_first, &translate,
                            err)) {
    return false;
  }

  std::vector<std::string> barcodes;
  barcodes.reserve(barcode_keys.size());
  for (uint64_t key : barcode_keys) {
    std::string barcode;
    if (!resolve_barcode_name(key, header.barcode_length, translate, &barcode,
                              err)) {
      return false;
    }
    barcodes.push_back(barcode);
  }

  if (!write_core_outputs(args, peaks, barcodes, entries, err)) {
    return false;
  }
  return write_metrics_sidecar(args.metrics_tsv, barcode_keys,
                               header.barcode_length, metrics, translate,
                               &barcodes, err);
}

Args to_internal_args(const MultiomeAtacPeakMexArgs &input) {
  Args args;
  args.fragments = input.fragments;
  args.sidecar = input.sidecar;
  args.peaks = input.peaks;
  args.summits_out = input.summits_out;
  args.out_dir = input.out_dir;
  args.metrics_tsv = input.metrics_tsv;
  args.barcode_translate = input.barcode_translate;
  args.barcode_translate_from_first = input.barcode_translate_from_first;
  args.call_peaks_from_sidecar = input.call_peaks_from_sidecar;
  args.force = input.force;
  args.max_barcodes = input.max_barcodes;
  args.threads = input.threads;
  args.macs3_pvalue = input.macs3_pvalue;
  args.macs3_min_length = input.macs3_min_length;
  args.macs3_max_gap = input.macs3_max_gap;
  args.macs3_uint8_counts = input.macs3_uint8_counts;
  args.temp_dir = input.temp_dir;
  args.keep_intermediates_dir = input.keep_intermediates_dir;
  return args;
}

bool validate_args(const Args &args, std::ostream &err) {
  const bool have_fragments = !args.fragments.empty();
  const bool have_sidecar = !args.sidecar.empty();
  if (have_fragments == have_sidecar) {
    err << "Exactly one of fragments or sidecar is required\n";
    return false;
  }
  if (args.peaks.empty()) {
    err << "Peak path is required\n";
    return false;
  }
  if (args.out_dir.empty()) {
    err << "Peak MEX output directory is required\n";
    return false;
  }
  if (args.metrics_tsv.empty()) {
    err << "ATAC metrics TSV path is required\n";
    return false;
  }
  if (args.threads < 1) {
    err << "Thread count must be >= 1\n";
    return false;
  }
  if (args.macs3_pvalue <= 0.0 || args.macs3_pvalue > 1.0) {
    err << "MACS3 FRAG p-value must be in (0, 1]\n";
    return false;
  }
  if (args.macs3_min_length < 1) {
    err << "MACS3 FRAG min length must be >= 1\n";
    return false;
  }
  if (args.macs3_max_gap < 0) {
    err << "MACS3 FRAG max gap must be >= 0\n";
    return false;
  }
  if (args.call_peaks_from_sidecar) {
    if (!have_sidecar) {
      err << "--call-peaks-from-sidecar requires --sidecar\n";
      return false;
    }
    if (args.summits_out.empty()) {
      err << "--call-peaks-from-sidecar requires --summits-out\n";
      return false;
    }
  }
  return true;
}

}  // namespace

int RunMultiomeAtacPeakMex(const MultiomeAtacPeakMexArgs &input) {
  const Args args = to_internal_args(input);
  if (!validate_args(args, std::cerr)) {
    return 2;
  }

  if (!call_peaks_from_sidecar_if_needed(args, std::cerr)) {
    return 1;
  }

  std::vector<Peak> peaks;
  std::unordered_map<std::string, ChromPeakIndex> peaks_by_chrom;
  if (!load_peaks(args.peaks, &peaks, &peaks_by_chrom, std::cerr)) {
    return 1;
  }

  const bool ok =
      !args.sidecar.empty()
          ? run_sidecar_path(args, peaks, peaks_by_chrom, std::cerr)
          : run_text_path(args, peaks, peaks_by_chrom, std::cerr);
  if (!ok) {
    return 1;
  }

  std::cout << "Wrote " << args.out_dir << "\n";
  std::cout << "Wrote " << args.metrics_tsv << "\n";
  std::cout << "peaks=" << peaks.size() << "\n";
  return 0;
}

}  // namespace multiome
}  // namespace star
