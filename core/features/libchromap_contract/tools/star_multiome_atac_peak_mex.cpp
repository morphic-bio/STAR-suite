#include <zlib.h>

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

struct Args {
  std::string fragments;
  std::string peaks;
  std::string out_dir;
  std::string metrics_tsv;
  uint64_t max_barcodes = 0;
  bool help = false;
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

void usage(const char *prog) {
  std::cerr
      << "Usage: " << prog
      << " --fragments fragments.tsv[.gz] --peaks peaks.narrowPeak "
         "--out-dir peak_mex --metrics-tsv atac_metrics.tsv [options]\n"
      << "\nOptions:\n"
      << "  --max-barcodes N   Keep only the top N barcodes by total fragments; "
         "0 keeps all (default 0)\n";
}

bool parse_uint64(const std::string &value, uint64_t *out) {
  if (out == nullptr || value.empty()) {
    return false;
  }
  char *end = nullptr;
  const unsigned long long parsed = std::strtoull(value.c_str(), &end, 10);
  if (end == value.c_str() || *end != '\0') {
    return false;
  }
  *out = static_cast<uint64_t>(parsed);
  return true;
}

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

bool parse_args(int argc, char **argv, Args *args) {
  if (args == nullptr) {
    return false;
  }
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      usage(argv[0]);
      args->help = true;
      return true;
    } else if (arg == "--fragments" && i + 1 < argc) {
      args->fragments = argv[++i];
    } else if (arg == "--peaks" && i + 1 < argc) {
      args->peaks = argv[++i];
    } else if (arg == "--out-dir" && i + 1 < argc) {
      args->out_dir = argv[++i];
    } else if (arg == "--metrics-tsv" && i + 1 < argc) {
      args->metrics_tsv = argv[++i];
    } else if (arg == "--max-barcodes" && i + 1 < argc) {
      if (!parse_uint64(argv[++i], &args->max_barcodes)) {
        std::cerr << "Invalid --max-barcodes value\n";
        return false;
      }
    } else {
      std::cerr << "Unknown or incomplete argument: " << arg << "\n";
      usage(argv[0]);
      return false;
    }
  }
  if (args->fragments.empty() || args->peaks.empty() || args->out_dir.empty() ||
      args->metrics_tsv.empty()) {
    usage(argv[0]);
    return false;
  }
  return true;
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

bool collect_total_fragments(
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

std::unordered_set<std::string> select_barcodes(
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

uint32_t barcode_index_for(
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

bool build_matrix(
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
                         barcode_index = barcode_index_for(
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

bool write_metrics(const std::string &path,
                   const std::vector<std::string> &barcodes,
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
    const BarcodeMetrics m = (it == metrics.end()) ? BarcodeMetrics() : it->second;
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

bool write_outputs(const Args &args, const std::vector<Peak> &peaks,
                   const std::vector<std::string> &barcodes,
                   const std::unordered_map<uint64_t, uint64_t> &entries,
                   const std::unordered_map<std::string, BarcodeMetrics> &metrics,
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
  if (!write_matrix(args.out_dir + "/matrix.mtx.gz", peaks.size(),
                    barcodes.size(), entries, err)) {
    return false;
  }
  return write_metrics(args.metrics_tsv, barcodes, metrics, err);
}

}  // namespace

int main(int argc, char **argv) {
  Args args;
  if (!parse_args(argc, argv, &args)) {
    return 2;
  }
  if (args.help) {
    return 0;
  }

  std::vector<Peak> peaks;
  std::unordered_map<std::string, ChromPeakIndex> peaks_by_chrom;
  if (!load_peaks(args.peaks, &peaks, &peaks_by_chrom, std::cerr)) {
    return 1;
  }

  std::unordered_map<std::string, BarcodeMetrics> metrics;
  if (!collect_total_fragments(args.fragments, &metrics, std::cerr)) {
    return 1;
  }
  const std::unordered_set<std::string> keep =
      select_barcodes(metrics, args.max_barcodes);

  std::unordered_map<uint64_t, uint64_t> entries;
  std::vector<std::string> barcodes;
  if (!build_matrix(args.fragments, peaks_by_chrom, keep, &metrics, &entries,
                    &barcodes, std::cerr)) {
    return 1;
  }
  if (barcodes.empty()) {
    std::cerr << "No peak-overlapping barcodes were found\n";
    return 1;
  }

  if (!write_outputs(args, peaks, barcodes, entries, metrics, std::cerr)) {
    return 1;
  }

  std::cout << "Wrote " << args.out_dir << "\n";
  std::cout << "Wrote " << args.metrics_tsv << "\n";
  std::cout << "peaks=" << peaks.size() << "\n";
  std::cout << "barcodes=" << barcodes.size() << "\n";
  std::cout << "matrix_nnz=" << entries.size() << "\n";
  return 0;
}
