// Stream a Cell Ranger Flex BAM/SAM and stratify its primary read records by
// STAR/Cell Ranger cell membership and by the sample-aware STAR hash verdict.
//
// BAM decoding is delegated to samtools so this diagnostic has no htslib
// dependency.  The FH01SEQ1 cache is mmap'ed read-only and binary-searched in
// place; a multi-gigabyte cache is therefore not copied into an unordered map.

#include <algorithm>
#include <cerrno>
#include <charconv>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

bool ends_with(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string shell_quote(const std::string& value) {
    std::string out("'");
    for (char c : value) {
        if (c == '\'') out += "'\\''";
        else out.push_back(c);
    }
    out.push_back('\'');
    return out;
}

std::string trim(std::string value) {
    const std::string spaces = " \t\r\n";
    const std::size_t first = value.find_first_not_of(spaces);
    if (first == std::string::npos) return {};
    const std::size_t last = value.find_last_not_of(spaces);
    return value.substr(first, last - first + 1);
}

std::vector<std::string_view> split_tabs(const std::string& line) {
    std::vector<std::string_view> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t end = line.find('\t', begin);
        fields.emplace_back(line.data() + begin,
                            (end == std::string::npos ? line.size() : end) - begin);
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return fields;
}

template <typename Integer>
bool parse_integer(std::string_view text, Integer& value) {
    if (text.empty()) return false;
    const char* begin = text.data();
    const char* end = begin + text.size();
    const auto result = std::from_chars(begin, end, value);
    return result.ec == std::errc() && result.ptr == end;
}

class LineInput {
public:
    LineInput(const std::string& path, bool gzip, const std::string& command = {}) {
        if (!command.empty()) {
            file_ = ::popen(command.c_str(), "r");
            is_pipe_ = true;
            description_ = command;
        } else if (path == "-") {
            file_ = stdin;
            description_ = "stdin";
        } else if (gzip) {
            description_ = path;
            const std::string gzip_command = "gzip -dc -- " + shell_quote(path);
            file_ = ::popen(gzip_command.c_str(), "r");
            is_pipe_ = true;
        } else {
            description_ = path;
            file_ = std::fopen(path.c_str(), "rb");
        }
        if (file_ == nullptr) {
            fail("cannot open " + description_ + ": " + std::strerror(errno));
        }
    }

    LineInput(const LineInput&) = delete;
    LineInput& operator=(const LineInput&) = delete;

    ~LineInput() {
        if (buffer_ != nullptr) std::free(buffer_);
        close_noexcept();
    }

    bool getline(std::string& line) {
        errno = 0;
        const ssize_t count = ::getline(&buffer_, &capacity_, file_);
        if (count < 0) {
            if (std::ferror(file_)) {
                fail("read failed for " + description_ + ": " + std::strerror(errno));
            }
            return false;
        }
        std::size_t size = static_cast<std::size_t>(count);
        while (size > 0 && (buffer_[size - 1] == '\n' || buffer_[size - 1] == '\r')) --size;
        line.assign(buffer_, size);
        return true;
    }

    void close_checked() {
        if (closed_) return;
        int status = 0;
        if (file_ != stdin) status = is_pipe_ ? ::pclose(file_) : std::fclose(file_);
        file_ = nullptr;
        closed_ = true;
        if (status != 0) fail("input command failed for " + description_);
    }

private:
    void close_noexcept() {
        if (closed_ || file_ == nullptr) return;
        if (file_ != stdin) {
            if (is_pipe_) ::pclose(file_);
            else std::fclose(file_);
        }
        file_ = nullptr;
        closed_ = true;
    }

    FILE* file_ = nullptr;
    char* buffer_ = nullptr;
    std::size_t capacity_ = 0;
    bool is_pipe_ = false;
    bool closed_ = false;
    std::string description_;
};

struct TaggedCellFile {
    std::string tag;
    std::string path;
};

struct Options {
    std::string input;
    std::string input_kind = "auto";
    std::string output_prefix;
    std::string samtools = "samtools";
    unsigned samtools_threads = 4;
    std::string region;
    std::string tag_map;
    std::string sample_tag;
    std::vector<TaggedCellFile> star_cells;
    std::vector<TaggedCellFile> cr_cells;
    std::string cache;
    std::string gene_list;
    bool single_n = false;
    bool include_secondary = false;
    std::string details = "discordant";
    std::uint64_t detail_max = 100000;
    std::uint64_t limit = 0;
};

void usage(std::ostream& out) {
    out <<
        "Usage: flex_bam_read_compare --input FILE.bam|FILE.sam|- --out-prefix PATH\n"
        "  --tag-map FILE                         BC001<TAB>TAG8 mapping\n"
        "  --star-cells [TAG=]barcodes.tsv[.gz]   repeat for each STAR tag\n"
        "  --cr-cells [TAG=]barcodes.tsv[.gz]     repeat; CR files normally contain CB16+TAG8\n"
        "  --hash-cache h01_cache.bin --gene-list probe_list.txt\n"
        "Options:\n"
        "  --input-kind auto|bam|sam              default auto\n"
        "  --sample-tag BC001|TAG8                fallback when a BAM barcode is only CB16\n"
        "  --samtools PATH --samtools-threads N   default samtools, 4\n"
        "  --region REGION                        passed to samtools view (BAM only)\n"
        "  --single-n yes|no                      optional single-N replay (default no; full classifier parity)\n"
        "  --include-secondary yes|no             default no\n"
        "  --details discordant|all|none          default discordant\n"
        "  --detail-max N                         deterministic sample size; 0 means unbounded\n"
        "  --limit N                              stop after N included primary records\n";
}

bool parse_bool(const std::string& value, const std::string& option) {
    if (value == "yes" || value == "true" || value == "1") return true;
    if (value == "no" || value == "false" || value == "0") return false;
    fail(option + " expects yes or no");
}

TaggedCellFile parse_cell_spec(const std::string& spec, bool tag_required) {
    const std::size_t equal = spec.find('=');
    if (equal == std::string::npos) {
        if (tag_required) fail("STAR cell input must be TAG=PATH: " + spec);
        return {"", spec};
    }
    if (equal == 0 || equal + 1 == spec.size()) fail("invalid TAG=PATH cell input: " + spec);
    return {spec.substr(0, equal), spec.substr(equal + 1)};
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&]() -> std::string {
            if (++i >= argc) fail(arg + " requires a value");
            return argv[i];
        };
        if (arg == "--input") options.input = value();
        else if (arg == "--input-kind") options.input_kind = value();
        else if (arg == "--out-prefix") options.output_prefix = value();
        else if (arg == "--samtools") options.samtools = value();
        else if (arg == "--samtools-threads") options.samtools_threads = std::stoul(value());
        else if (arg == "--region") options.region = value();
        else if (arg == "--tag-map") options.tag_map = value();
        else if (arg == "--sample-tag") options.sample_tag = value();
        else if (arg == "--star-cells") options.star_cells.push_back(parse_cell_spec(value(), false));
        else if (arg == "--cr-cells") options.cr_cells.push_back(parse_cell_spec(value(), false));
        else if (arg == "--hash-cache") options.cache = value();
        else if (arg == "--gene-list") options.gene_list = value();
        else if (arg == "--single-n") options.single_n = parse_bool(value(), arg);
        else if (arg == "--include-secondary") options.include_secondary = parse_bool(value(), arg);
        else if (arg == "--details") options.details = value();
        else if (arg == "--detail-max") options.detail_max = std::stoull(value());
        else if (arg == "--limit") options.limit = std::stoull(value());
        else if (arg == "--help" || arg == "-h") {
            usage(std::cout);
            std::exit(0);
        } else fail("unknown option: " + arg);
    }
    if (options.input.empty() || options.output_prefix.empty()) {
        usage(std::cerr);
        fail("--input and --out-prefix are required");
    }
    if (options.input_kind == "auto") {
        options.input_kind = ends_with(options.input, ".bam") ? "bam" : "sam";
    }
    if (options.input_kind != "bam" && options.input_kind != "sam") {
        fail("--input-kind must be auto, bam, or sam");
    }
    if (options.details != "discordant" && options.details != "all" && options.details != "none") {
        fail("--details must be discordant, all, or none");
    }
    if (options.star_cells.empty() || options.cr_cells.empty()) {
        fail("at least one --star-cells and one --cr-cells input are required");
    }
    if (options.cache.empty() != options.gene_list.empty()) {
        fail("--hash-cache and --gene-list must be supplied together");
    }
    if (options.samtools_threads == 0) fail("--samtools-threads must be positive");
    if (options.input_kind == "sam" && !options.region.empty()) {
        fail("--region is supported only for BAM input");
    }
    return options;
}

struct TagMap {
    std::unordered_map<std::string, std::string> id_to_sequence;
    std::unordered_map<std::string, std::string> sequence_to_id;
    std::unordered_map<std::string, std::uint16_t> id_to_index;
    std::unordered_map<std::string, std::uint16_t> sequence_to_index;

    std::string id_for(std::string tag) const {
        const auto seq_it = sequence_to_id.find(tag);
        if (seq_it != sequence_to_id.end()) return seq_it->second;
        return tag;
    }

    std::string sequence_for(const std::string& tag) const {
        const auto id_it = id_to_sequence.find(tag);
        if (id_it != id_to_sequence.end()) return id_it->second;
        return tag;
    }

    std::uint16_t index_for(const std::string& tag) const {
        const auto by_id = id_to_index.find(tag);
        if (by_id != id_to_index.end()) return by_id->second;
        const auto by_sequence = sequence_to_index.find(tag);
        return by_sequence == sequence_to_index.end() ? 0 : by_sequence->second;
    }
};

TagMap load_tag_map(const std::string& path) {
    TagMap result;
    if (path.empty()) return result;
    LineInput input(path, ends_with(path, ".gz"));
    std::string line;
    std::uint64_t line_number = 0;
    std::uint32_t sequential_index = 0;
    while (input.getline(line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::string id, sequence;
        if (!(fields >> id >> sequence)) fail(path + ": invalid tag map line " + std::to_string(line_number));
        if (id.size() == 8 && sequence.rfind("BC", 0) == 0) std::swap(id, sequence);
        if (sequence.size() != 8) fail(path + ": tag sequence is not eight bases at line " + std::to_string(line_number));
        if (++sequential_index > std::numeric_limits<std::uint16_t>::max()) {
            fail(path + ": too many sample tags");
        }
        const auto old_id = result.id_to_sequence.emplace(id, sequence);
        const auto old_seq = result.sequence_to_id.emplace(sequence, id);
        if ((!old_id.second && old_id.first->second != sequence) ||
            (!old_seq.second && old_seq.first->second != id)) {
            fail(path + ": conflicting tag mapping at line " + std::to_string(line_number));
        }
        result.id_to_index.emplace(id, static_cast<std::uint16_t>(sequential_index));
        result.sequence_to_index.emplace(sequence, static_cast<std::uint16_t>(sequential_index));
    }
    input.close_checked();
    return result;
}

std::string strip_dash_one(std::string barcode) {
    barcode = trim(std::move(barcode));
    if (ends_with(barcode, "-1")) barcode.resize(barcode.size() - 2);
    return barcode;
}

std::string cell_key(std::string barcode, const std::string& forced_tag, const TagMap& tag_map) {
    barcode = strip_dash_one(std::move(barcode));
    if (barcode.empty()) return {};
    std::string cb;
    std::string tag;
    const std::size_t pipe = barcode.find('|');
    if (pipe != std::string::npos) {
        cb = barcode.substr(0, pipe);
        tag = barcode.substr(pipe + 1);
    } else {
        if (barcode.size() < 16) return {};
        cb = barcode.substr(0, 16);
        if (barcode.size() >= 24) tag = barcode.substr(16, 8);
    }
    if (tag.empty()) tag = forced_tag;
    if (tag.empty()) return cb;
    return cb + "|" + tag_map.id_for(tag);
}

std::unordered_set<std::string> load_cells(const std::vector<TaggedCellFile>& files,
                                           const TagMap& tag_map,
                                           const std::string& label) {
    std::unordered_set<std::string> result;
    for (const auto& input_spec : files) {
        LineInput input(input_spec.path, ends_with(input_spec.path, ".gz"));
        std::string line;
        std::uint64_t line_number = 0;
        while (input.getline(line)) {
            ++line_number;
            if (line.empty() || line[0] == '#') continue;
            const std::size_t delimiter = line.find_first_of("\t ,");
            const std::string raw = line.substr(0, delimiter);
            const std::string key = cell_key(raw, input_spec.tag, tag_map);
            if (key.empty() || key.size() < 16) {
                fail(input_spec.path + ": invalid " + label + " barcode at line " +
                     std::to_string(line_number));
            }
            result.insert(key);
        }
        input.close_checked();
    }
    return result;
}

#pragma pack(push, 1)
struct CacheHeaderRaw {
    char magic[8];
    std::uint16_t version;
    std::uint16_t kmer_length;
    std::uint32_t record_size;
    std::uint64_t record_count;
};

struct CacheRecordRaw {
    std::uint64_t seq_lo;
    std::uint64_t seq_hi;
    std::uint32_t resolved_gene;
    std::uint8_t cache_class;
    std::uint8_t negative_code;
    std::uint16_t sample_index;
};
#pragma pack(pop)

static_assert(sizeof(CacheHeaderRaw) == 24, "unexpected cache header layout");
static_assert(sizeof(CacheRecordRaw) == 24, "unexpected cache record layout");

struct HashDecision {
    std::string verdict = "DISABLED";
    std::uint16_t gene_index = 0;
    std::uint8_t cache_class = 0;
    std::uint8_t negative_code = 0;
    std::uint8_t probe_region = 0;
    int offset = 0;
};

class HashCacheView {
public:
    HashCacheView() = default;
    HashCacheView(const HashCacheView&) = delete;
    HashCacheView& operator=(const HashCacheView&) = delete;

    ~HashCacheView() {
        if (mapping_ != MAP_FAILED) ::munmap(mapping_, mapping_size_);
        if (fd_ >= 0) ::close(fd_);
    }

    void open(const std::string& path) {
        const std::uint16_t endian_test = 1;
        if (*reinterpret_cast<const std::uint8_t*>(&endian_test) != 1) {
            fail("FH01SEQ1 cache reader requires a little-endian host");
        }
        fd_ = ::open(path.c_str(), O_RDONLY);
        if (fd_ < 0) fail("cannot open hash cache " + path + ": " + std::strerror(errno));
        struct stat info {};
        if (::fstat(fd_, &info) != 0 || info.st_size < static_cast<off_t>(sizeof(CacheHeaderRaw))) {
            fail("cannot stat or cache is too short: " + path);
        }
        mapping_size_ = static_cast<std::size_t>(info.st_size);
        mapping_ = ::mmap(nullptr, mapping_size_, PROT_READ, MAP_PRIVATE, fd_, 0);
        if (mapping_ == MAP_FAILED) fail("cannot mmap hash cache " + path + ": " + std::strerror(errno));
        header_ = static_cast<const CacheHeaderRaw*>(mapping_);
        const char expected[8] = {'F','H','0','1','S','E','Q','1'};
        if (std::memcmp(header_->magic, expected, sizeof(expected)) != 0 ||
            header_->version < 1 || header_->version > 3 ||
            header_->kmer_length != 50 || header_->record_size != sizeof(CacheRecordRaw)) {
            fail("unsupported FH01SEQ1 cache header: " + path);
        }
        if (header_->record_count >
            (mapping_size_ - sizeof(CacheHeaderRaw)) / sizeof(CacheRecordRaw)) {
            fail("truncated FH01SEQ1 cache: " + path);
        }
        const std::size_t expected_size = sizeof(CacheHeaderRaw) +
            static_cast<std::size_t>(header_->record_count) * sizeof(CacheRecordRaw);
        if (expected_size != mapping_size_) fail("FH01SEQ1 cache has trailing or missing bytes: " + path);
        records_ = reinterpret_cast<const CacheRecordRaw*>(
            static_cast<const std::uint8_t*>(mapping_) + sizeof(CacheHeaderRaw));
        enabled_ = true;
    }

    bool enabled() const { return enabled_; }
    std::uint64_t record_count() const { return enabled_ ? header_->record_count : 0; }
    std::uint16_t version() const { return enabled_ ? header_->version : 0; }

    HashDecision classify(const std::string& sequence, std::uint16_t sample_index,
                          bool single_n) const {
        if (!enabled_) return {};
        if (sequence.size() < 50) return {"SHORT", 0, 0, 0, 0, 0};
        HashDecision direct = classify_full(sequence, sample_index);
        if (!single_n || direct.verdict != "MISS") return direct;
        int invalid_position = -1;
        const int scan_length = static_cast<int>(std::min<std::size_t>(sequence.size(), 51));
        for (int i = 0; i < scan_length; ++i) {
            const char c = sequence[static_cast<std::size_t>(i)];
            if (c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
                c == 'a' || c == 'c' || c == 'g' || c == 't') continue;
            if (invalid_position >= 0) return {"N_MISS", 0, 0, 0, 0, 0};
            invalid_position = i;
        }
        if (invalid_position < 0) return direct;
        std::string repaired = sequence;
        bool have_keep = false;
        bool conflict = false;
        bool saw_deny = false;
        HashDecision keep;
        for (char base : std::string("ACGT")) {
            repaired[static_cast<std::size_t>(invalid_position)] = base;
            const HashDecision decision = classify_full(repaired, sample_index);
            if (decision.verdict.rfind("DENY", 0) == 0) saw_deny = true;
            if (decision.verdict.rfind("KEEP", 0) == 0) {
                if (!have_keep) {
                    keep = decision;
                    have_keep = true;
                } else if (keep.gene_index != decision.gene_index) {
                    conflict = true;
                }
            }
        }
        if (have_keep && !conflict && !saw_deny) {
            keep.verdict = "N_KEEP";
            keep.cache_class = 1;
            return keep;
        }
        return {"N_MISS", 0, 0, 0, 0, 0};
    }

private:
    static bool encode(const char* sequence, std::uint64_t& lo, std::uint64_t& hi) {
        lo = 0;
        hi = 0;
        for (unsigned i = 0; i < 50; ++i) {
            std::uint64_t code = 0;
            switch (sequence[i]) {
                case 'C': case 'c': code = 1; break;
                case 'G': case 'g': code = 2; break;
                case 'T': case 't': code = 3; break;
                case 'A': case 'a': break;
                default: return false;
            }
            hi = (hi << 2) | (lo >> 62);
            lo = (lo << 2) | code;
        }
        return true;
    }

    bool key_less(std::size_t index, std::uint64_t hi, std::uint64_t lo) const {
        return records_[index].seq_hi < hi ||
               (records_[index].seq_hi == hi && records_[index].seq_lo < lo);
    }

    std::pair<std::size_t, std::size_t> equal_range(std::uint64_t lo, std::uint64_t hi) const {
        std::size_t first = 0;
        std::size_t last = static_cast<std::size_t>(header_->record_count);
        while (first < last) {
            const std::size_t middle = first + (last - first) / 2;
            if (key_less(middle, hi, lo)) first = middle + 1;
            else last = middle;
        }
        const std::size_t begin = first;
        while (first < static_cast<std::size_t>(header_->record_count) &&
               records_[first].seq_hi == hi && records_[first].seq_lo == lo) ++first;
        return {begin, first};
    }

    bool find_record(std::uint64_t lo, std::uint64_t hi, std::uint16_t sample_index,
                     CacheRecordRaw& output) const {
        const auto range = equal_range(lo, hi);
        for (std::size_t i = range.first; i < range.second; ++i) {
            if (records_[i].sample_index == sample_index) {
                output = records_[i];
                return true;
            }
        }
        if (sample_index == 0) return false;
        for (std::size_t i = range.first; i < range.second; ++i) {
            if (records_[i].sample_index == 0) {
                output = records_[i];
                return true;
            }
        }
        return false;
    }

    std::uint8_t region(const CacheRecordRaw& record) const {
        return header_->version >= 3
            ? static_cast<std::uint8_t>((record.resolved_gene >> 30) & 0x3u)
            : 0;
    }

    static std::uint8_t merge_region(std::uint8_t lhs, std::uint8_t rhs) {
        if (lhs == 3 || rhs == 3) return 3;
        if (lhs == 0) return rhs;
        if (rhs == 0) return lhs;
        return lhs == rhs ? lhs : 3;
    }

    HashDecision classify_hits(const CacheRecordRaw* const hits[3],
                               const int relative_offsets[3],
                               std::uint16_t runtime_sample) const {
        bool saw_ambiguous = false;
        bool saw_sample_mismatch = false;
        int sample_mismatch_offset = 0;
        bool saw_keep = false;
        std::uint16_t keep_gene = 0;
        std::uint16_t keep_sample = 0;
        std::uint8_t keep_class = 0;
        std::uint8_t keep_region = 0;
        int keep_offset = 0;
        bool saw_gene_conflict = false;
        int conflict_offset = 0;
        std::uint8_t negative_code = 0;
        int negative_offset = 0;

        for (int i = 0; i < 3; ++i) {
            if (hits[i] == nullptr) continue;
            const CacheRecordRaw& record = *hits[i];
            const std::uint16_t gene = static_cast<std::uint16_t>(record.resolved_gene & 0x7fffu);
            const bool sample_matched = record.sample_index != 0 && record.sample_index == runtime_sample;
            const bool sample_mismatch = record.sample_index != 0 && record.sample_index != runtime_sample;
            if (gene == 0 || record.cache_class == 2) {
                saw_ambiguous = true;
                negative_code = record.negative_code == 0 ? 1 : record.negative_code;
                negative_offset = relative_offsets[i];
                continue;
            }
            if (record.cache_class == 0 && sample_matched) {
                return {"KEEP_H0", gene, record.cache_class, 0, region(record), relative_offsets[i]};
            }
            if ((record.cache_class == 0 || record.cache_class == 1 || record.cache_class == 3) &&
                sample_mismatch) {
                if (!saw_sample_mismatch) {
                    saw_sample_mismatch = true;
                    sample_mismatch_offset = relative_offsets[i];
                }
                continue;
            }
            if (record.cache_class == 0 || record.cache_class == 1 || record.cache_class == 3) {
                const std::uint16_t sample_key = sample_matched ? runtime_sample : 0;
                if (!saw_keep) {
                    saw_keep = true;
                    keep_gene = gene;
                    keep_sample = sample_key;
                    keep_class = record.cache_class;
                    keep_region = region(record);
                    keep_offset = relative_offsets[i];
                } else if (keep_gene != gene || keep_sample != sample_key) {
                    saw_gene_conflict = true;
                    conflict_offset = relative_offsets[i];
                } else {
                    keep_region = merge_region(keep_region, region(record));
                }
            }
        }
        if (saw_ambiguous) return {"DENY_CACHE", 0, 0, negative_code, 0, negative_offset};
        if (saw_gene_conflict) return {"DENY_GENE_CONFLICT", 0, 0, 1, 0, conflict_offset};
        if (saw_keep) {
            const std::string verdict = keep_class == 0 ? "KEEP_H0_GLOBAL" :
                                        keep_class == 1 ? "KEEP_H1" : "KEEP_H2";
            return {verdict, keep_gene, keep_class, 0, keep_region, keep_offset};
        }
        if (saw_sample_mismatch) return {"DENY_SAMPLE_MISMATCH", 0, 0, 1, 0,
                                         sample_mismatch_offset};
        return {"MISS", 0, 0, 0, 0, 0};
    }

    HashDecision classify_full(const std::string& sequence, std::uint16_t sample_index) const {
        const int relative_offsets[3] = {0, 1, -1};
        CacheRecordRaw hit_storage[3] {};
        const CacheRecordRaw* hits[3] = {nullptr, nullptr, nullptr};
        for (int i = 0; i < 3; ++i) {
            const int start = relative_offsets[i];
            if (start < 0 || static_cast<std::size_t>(start + 50) > sequence.size()) continue;
            std::uint64_t lo = 0, hi = 0;
            if (!encode(sequence.data() + start, lo, hi)) continue;
            if (find_record(lo, hi, sample_index, hit_storage[i])) hits[i] = &hit_storage[i];
        }
        return classify_hits(hits, relative_offsets, sample_index);
    }

    int fd_ = -1;
    void* mapping_ = MAP_FAILED;
    std::size_t mapping_size_ = 0;
    const CacheHeaderRaw* header_ = nullptr;
    const CacheRecordRaw* records_ = nullptr;
    bool enabled_ = false;
};

std::vector<std::string> load_gene_list(const std::string& path) {
    std::vector<std::string> genes;
    if (path.empty()) return genes;
    LineInput input(path, ends_with(path, ".gz"));
    std::string line;
    while (input.getline(line)) {
        line = trim(std::move(line));
        if (!line.empty() && line[0] != '#') genes.push_back(line);
    }
    input.close_checked();
    return genes;
}

std::string gene_id(std::uint16_t index, const std::vector<std::string>& genes) {
    if (index == 0) return {};
    if (index <= genes.size()) return genes[index - 1];
    return "#" + std::to_string(index);
}

std::string reverse_complement(std::string_view sequence) {
    std::string result(sequence.size(), 'N');
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        switch (sequence[sequence.size() - 1 - i]) {
            case 'A': case 'a': result[i] = 'T'; break;
            case 'C': case 'c': result[i] = 'G'; break;
            case 'G': case 'g': result[i] = 'C'; break;
            case 'T': case 't': result[i] = 'A'; break;
            default: result[i] = 'N'; break;
        }
    }
    return result;
}

std::string reverse_string(std::string_view value) {
    return std::string(value.rbegin(), value.rend());
}

std::string get_tag_value(const std::vector<std::string_view>& fields, std::string_view tag) {
    for (std::size_t i = 11; i < fields.size(); ++i) {
        const std::string_view field = fields[i];
        if (field.size() >= 5 && field.substr(0, 2) == tag && field[2] == ':' && field[4] == ':') {
            return std::string(field.substr(5));
        }
    }
    return {};
}

std::string all_tags(const std::vector<std::string_view>& fields) {
    std::string result;
    for (std::size_t i = 11; i < fields.size(); ++i) {
        if (!result.empty()) result.push_back(';');
        for (char c : fields[i]) result.push_back(c == '\t' ? ' ' : c);
    }
    return result;
}

std::vector<std::string> gene_set(std::string value) {
    std::vector<std::string> genes;
    std::size_t begin = 0;
    while (begin <= value.size()) {
        std::size_t end = value.find_first_of(";,", begin);
        if (end == std::string::npos) end = value.size();
        std::string gene = trim(value.substr(begin, end - begin));
        if (!gene.empty()) {
            const std::size_t dot = gene.rfind('.');
            if (dot != std::string::npos && dot + 1 < gene.size() &&
                std::all_of(gene.begin() + static_cast<std::ptrdiff_t>(dot + 1), gene.end(),
                            [](unsigned char c) { return std::isdigit(c) != 0; })) {
                gene.resize(dot);
            }
            genes.push_back(std::move(gene));
        }
        if (end == value.size()) break;
        begin = end + 1;
    }
    std::sort(genes.begin(), genes.end());
    genes.erase(std::unique(genes.begin(), genes.end()), genes.end());
    return genes;
}

std::string join(const std::vector<std::string>& values, char delimiter) {
    std::string result;
    for (const auto& value : values) {
        if (!result.empty()) result.push_back(delimiter);
        result += value;
    }
    return result;
}

std::string gene_relation(const std::string& hash_gene, const std::vector<std::string>& cr_genes) {
    if (hash_gene.empty() && cr_genes.empty()) return "neither";
    if (hash_gene.empty()) return cr_genes.size() == 1 ? "cr_only" : "cr_multi_only";
    if (cr_genes.empty()) return "hash_only";
    const bool contains = std::binary_search(cr_genes.begin(), cr_genes.end(), hash_gene);
    if (cr_genes.size() == 1 && contains) return "same";
    if (contains) return "cr_multi_contains_hash";
    return cr_genes.size() == 1 ? "conflict" : "cr_multi_conflict";
}

bool hash_keep(const HashDecision& decision) {
    return decision.verdict == "N_KEEP" || decision.verdict.rfind("KEEP_", 0) == 0;
}

std::string probe_region_name(std::uint8_t region) {
    switch (region) {
        case 1: return "spliced";
        case 2: return "unspliced";
        case 3: return "conflicting";
        default: return "unknown";
    }
}

std::string counted_record_outcome(const HashDecision& hash, const std::string& relation) {
    if (hash_keep(hash)) {
        if (relation == "same") return "same_gene_accept";
        if (relation == "cr_multi_contains_hash") return "multi_gene_contains_accept";
        if (relation == "conflict" || relation == "cr_multi_conflict") return "different_gene_conflict";
        return "hash_accept_" + relation;
    }
    if (hash.verdict == "DENY_CACHE") return "cache_deny";
    if (hash.verdict == "DENY_GENE_CONFLICT") return "cache_nearby_gene_conflict";
    if (hash.verdict == "DENY_SAMPLE_MISMATCH") return "cache_sample_mismatch";
    if (hash.verdict == "MISS") return "cache_miss";
    return hash.verdict;
}

std::string cell_class(const std::string& key,
                       const std::unordered_set<std::string>& star,
                       const std::unordered_set<std::string>& cr) {
    if (key.empty()) return "no_barcode";
    const bool in_star = star.count(key) != 0;
    const bool in_cr = cr.count(key) != 0;
    if (in_star && in_cr) return "shared";
    if (in_star) return "star_only";
    if (in_cr) return "cr_only";
    return "neither";
}

double mean_quality(const std::string& quality, std::size_t limit = std::numeric_limits<std::size_t>::max()) {
    if (quality.empty() || quality == "*") return std::numeric_limits<double>::quiet_NaN();
    const std::size_t size = std::min(quality.size(), limit);
    if (size == 0) return std::numeric_limits<double>::quiet_NaN();
    std::uint64_t sum = 0;
    for (std::size_t i = 0; i < size; ++i) sum += static_cast<unsigned char>(quality[i]) - 33u;
    return static_cast<double>(sum) / static_cast<double>(size);
}

int min_quality(const std::string& quality) {
    if (quality.empty() || quality == "*") return -1;
    int result = 255;
    for (unsigned char c : quality) result = std::min(result, static_cast<int>(c) - 33);
    return result;
}

std::string mapq_bin(unsigned mapq) {
    if (mapq == 255) return "255/unavailable";
    if (mapq == 0) return "0";
    if (mapq < 10) return "1-9";
    if (mapq < 30) return "10-29";
    if (mapq < 60) return "30-59";
    return "60+";
}

struct SummaryKey {
    std::string cell;
    std::string hash;
    bool counted = false;
    std::string gene;

    bool operator<(const SummaryKey& other) const {
        return std::tie(cell, hash, counted, gene) <
               std::tie(other.cell, other.hash, other.counted, other.gene);
    }
};

struct Stats {
    std::uint64_t records = 0;
    std::uint64_t mapped = 0;
    std::uint64_t duplicate = 0;
    std::uint64_t mapq_sum = 0;
    long double probe_quality_sum = 0;
    std::uint64_t probe_quality_n = 0;
    std::uint64_t umi_q_lt20 = 0;
    std::uint64_t umi_quality_n = 0;

    void add(bool is_mapped, bool is_duplicate, unsigned mapq, double probe_q, int umi_min_q) {
        ++records;
        mapped += is_mapped;
        duplicate += is_duplicate;
        mapq_sum += mapq;
        if (std::isfinite(probe_q)) {
            probe_quality_sum += probe_q;
            ++probe_quality_n;
        }
        if (umi_min_q >= 0) {
            ++umi_quality_n;
            umi_q_lt20 += umi_min_q < 20;
        }
    }
};

std::uint64_t fnv1a(std::string_view value) {
    std::uint64_t hash = UINT64_C(14695981039346656037);
    for (unsigned char byte : value) {
        hash ^= byte;
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

struct SampledDetail {
    std::uint64_t score;
    std::string line;
    bool operator<(const SampledDetail& other) const {
        if (score != other.score) return score < other.score;
        return line < other.line;
    }
};

class DetailWriter {
public:
    DetailWriter(const std::string& path, const std::string& mode, std::uint64_t maximum)
        : mode_(mode), maximum_(maximum), output_(path) {
        if (!output_) fail("cannot write " + path);
        output_ << "qname\tcell_key\tcell_class\tcr_counted\thash_verdict\tcache_class"
                   "\tnegative_code\thash_offset\tsample_index\tprobe_region\tcr_region"
                   "\thash_gene\tquant_gene_ids\tgenomic_gene_ids\tgenomic_gene_names"
                   "\tgene_relation\tflag\trname\tpos"
                   "\tmapq\tmapq_bin\tcigar\txf\tCB\tCR\tUB\tUR\tprobe_tag\tprobe_mean_q"
                   "\tcb_mean_q\tumi_mean_q\tumi_min_q\tall_tags\n";
    }

    bool wants(const std::string& cell, const HashDecision& hash, bool counted,
               const std::string& relation) const {
        if (mode_ == "none") return false;
        if (mode_ == "all") return true;
        return cell == "star_only" || cell == "cr_only" ||
               (counted && !hash_keep(hash)) || (hash_keep(hash) && relation != "same");
    }

    void add(std::string_view qname, std::string line) {
        ++eligible_;
        if (maximum_ == 0) {
            output_ << line << '\n';
            ++written_;
            return;
        }
        SampledDetail item{fnv1a(qname), std::move(line)};
        if (sample_.size() < maximum_) sample_.push(std::move(item));
        else if (item < sample_.top()) {
            sample_.pop();
            sample_.push(std::move(item));
        }
    }

    void finish() {
        if (maximum_ == 0) return;
        std::vector<SampledDetail> rows;
        rows.reserve(sample_.size());
        while (!sample_.empty()) {
            // priority_queue::top is const; copy before pop so the heap remains
            // valid until pop_heap has completed.
            rows.push_back(sample_.top());
            sample_.pop();
        }
        std::sort(rows.begin(), rows.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.score < rhs.score || (lhs.score == rhs.score && lhs.line < rhs.line);
        });
        for (const auto& row : rows) output_ << row.line << '\n';
        written_ = rows.size();
    }

    std::uint64_t eligible() const { return eligible_; }
    std::uint64_t written() const { return written_; }

private:
    std::string mode_;
    std::uint64_t maximum_;
    std::ofstream output_;
    std::priority_queue<SampledDetail> sample_;
    std::uint64_t eligible_ = 0;
    std::uint64_t written_ = 0;
};

struct GeneKey {
    std::string cell;
    std::string hash_gene;
    std::string cr_genes;
    std::string relation;
    bool counted = false;
    bool operator<(const GeneKey& other) const {
        return std::tie(cell, hash_gene, cr_genes, relation, counted) <
               std::tie(other.cell, other.hash_gene, other.cr_genes, other.relation, other.counted);
    }
};

struct CigarKey {
    std::string cell;
    std::string hash;
    bool counted = false;
    std::string relation;
    std::string cigar;
    std::string mapq;
    bool operator<(const CigarKey& other) const {
        return std::tie(cell, hash, counted, relation, cigar, mapq) <
               std::tie(other.cell, other.hash, other.counted, other.relation, other.cigar, other.mapq);
    }
};

struct RegionKey {
    std::string cell;
    std::string cr_region;
    std::string probe_region;
    std::string hash;
    std::string hash_gene;
    std::string quant_genes;
    std::string genomic_genes;
    std::string genomic_names;
    std::string relation;
    bool counted = false;
    bool operator<(const RegionKey& other) const {
        return std::tie(cell, cr_region, probe_region, hash, hash_gene, quant_genes,
                        genomic_genes, genomic_names, relation, counted) <
               std::tie(other.cell, other.cr_region, other.probe_region, other.hash,
                        other.hash_gene, other.quant_genes, other.genomic_genes,
                        other.genomic_names, other.relation, other.counted);
    }
};

struct CountedRecordKey {
    std::string outcome;
    std::string cell;
    std::string quant_gene_ids;
    std::string probe;
    std::string hash_gene;
    std::string genomic_gene_ids;
    std::string genomic_gene_names;
    std::string cr_region;
    std::string probe_region;
    bool operator<(const CountedRecordKey& other) const {
        return std::tie(outcome, cell, quant_gene_ids, probe, hash_gene, genomic_gene_ids,
                        genomic_gene_names, cr_region, probe_region) <
               std::tie(other.outcome, other.cell, other.quant_gene_ids, other.probe,
                        other.hash_gene, other.genomic_gene_ids, other.genomic_gene_names,
                        other.cr_region, other.probe_region);
    }
};

class StringInterner {
public:
    StringInterner() {
        values_.emplace_back();
        ids_.emplace(std::string(), 0);
    }

    std::uint32_t intern(const std::string& value) {
        const auto found = ids_.find(value);
        if (found != ids_.end()) return found->second;
        if (values_.size() > std::numeric_limits<std::uint32_t>::max()) {
            fail("too many interned evidence strings");
        }
        const std::uint32_t id = static_cast<std::uint32_t>(values_.size());
        values_.push_back(value);
        ids_.emplace(values_.back(), id);
        return id;
    }

    const std::string& get(std::uint32_t id) const {
        if (id >= values_.size()) fail("invalid interned string id");
        return values_[id];
    }

private:
    std::vector<std::string> values_;
    std::unordered_map<std::string, std::uint32_t> ids_;
};

bool pack_acgt(std::string_view sequence, std::uint64_t& packed) {
    if (sequence.empty() || sequence.size() > 32) return false;
    packed = 0;
    for (char base : sequence) {
        std::uint64_t code = 0;
        switch (base) {
            case 'A': case 'a': break;
            case 'C': case 'c': code = 1; break;
            case 'G': case 'g': code = 2; break;
            case 'T': case 't': code = 3; break;
            default: return false;
        }
        packed = (packed << 2) | code;
    }
    return true;
}

bool pack_composite_cell(const std::string& key, const TagMap& tag_map,
                         std::uint64_t& packed) {
    const std::size_t pipe = key.find('|');
    if (pipe != 16 || pipe + 1 >= key.size()) return false;
    const std::string tag_sequence = tag_map.sequence_for(key.substr(pipe + 1));
    if (tag_sequence.size() != 8) return false;
    return pack_acgt(key.substr(0, 16) + tag_sequence, packed);
}

std::uint64_t mix64(std::uint64_t value) {
    value ^= value >> 30;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27;
    value *= UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

struct PackedMoleculeKey {
    std::uint64_t cell24 = 0;
    std::uint64_t umi = 0;
    std::uint32_t quant_gene = 0;
    std::uint8_t umi_length = 0;

    bool operator==(const PackedMoleculeKey& other) const {
        return cell24 == other.cell24 && umi == other.umi &&
               quant_gene == other.quant_gene && umi_length == other.umi_length;
    }
};

std::uint64_t molecule_hash(const PackedMoleculeKey& key) {
    std::uint64_t hash = mix64(key.cell24 + UINT64_C(0x9e3779b97f4a7c15));
    hash ^= mix64(key.umi + UINT64_C(0x243f6a8885a308d3));
    hash ^= mix64((static_cast<std::uint64_t>(key.quant_gene) << 8) | key.umi_length);
    return hash == 0 ? 1 : hash; // zero is the empty-slot sentinel
}

struct MoleculeEvidence {
    std::uint32_t hash_verdict = 0;
    std::uint32_t hash_gene = 0;
    std::uint32_t probe = 0;
    std::uint32_t genomic_gene_ids = 0;
    std::uint32_t genomic_gene_names = 0;
    std::uint32_t cr_region = 0;
    std::uint32_t rname = 0;
    std::uint32_t cigar = 0;
    std::uint32_t position = 0;
    std::uint16_t selection_score = 0;
    std::uint16_t mapq = 0;
    std::int8_t hash_offset = 0;
    std::uint8_t probe_region = 0;
    std::uint8_t priority = 0;
};

struct MoleculeState {
    MoleculeEvidence best;
    std::uint32_t records = 0;
    std::uint32_t xf8_records = 0;
    std::uint32_t cell_class = 0;
    std::uint8_t best_xf8_priority = 0;
    bool used_ur_fallback = false;
};

struct MoleculeSlot {
    std::uint64_t hash = 0;
    PackedMoleculeKey key;
    MoleculeState state;
};

class MoleculeTable {
public:
    MoleculeTable() : slots_(1024) {}

    MoleculeState& get_or_insert(const PackedMoleculeKey& key) {
        if ((size_ + 1) * 10 >= slots_.size() * 7) rehash(slots_.size() * 2);
        return find_or_insert(key, molecule_hash(key), slots_, size_);
    }

    const std::vector<MoleculeSlot>& slots() const { return slots_; }
    std::size_t size() const { return size_; }
    std::size_t capacity() const { return slots_.size(); }

private:
    static MoleculeState& find_or_insert(const PackedMoleculeKey& key, std::uint64_t hash,
                                         std::vector<MoleculeSlot>& slots, std::size_t& size) {
        const std::size_t mask = slots.size() - 1;
        std::size_t index = static_cast<std::size_t>(hash) & mask;
        while (true) {
            MoleculeSlot& slot = slots[index];
            if (slot.hash == 0) {
                slot.hash = hash;
                slot.key = key;
                slot.state = MoleculeState{};
                ++size;
                return slot.state;
            }
            if (slot.hash == hash && slot.key == key) return slot.state;
            index = (index + 1) & mask;
        }
    }

    void rehash(std::size_t capacity) {
        std::vector<MoleculeSlot> replacement(capacity);
        std::size_t replacement_size = 0;
        for (const MoleculeSlot& old : slots_) {
            if (old.hash == 0) continue;
            MoleculeState& state = find_or_insert(old.key, old.hash, replacement,
                                                   replacement_size);
            state = old.state;
        }
        slots_.swap(replacement);
        size_ = replacement_size;
    }

    std::vector<MoleculeSlot> slots_;
    std::size_t size_ = 0;
};

std::uint8_t evidence_priority(const HashDecision& hash, const std::string& relation) {
    if (hash_keep(hash)) {
        return relation == "same" || relation == "cr_multi_contains_hash" ? 4 : 3;
    }
    return hash.verdict == "MISS" ? 2 : 1;
}

std::uint16_t evidence_score(const HashDecision& hash, std::uint8_t priority, unsigned mapq) {
    unsigned verdict_rank = 0;
    if (hash.verdict == "KEEP_H0") verdict_rank = 6;
    else if (hash.verdict == "KEEP_H0_GLOBAL") verdict_rank = 5;
    else if (hash.verdict == "KEEP_H1") verdict_rank = 4;
    else if (hash.verdict == "KEEP_H2") verdict_rank = 3;
    else if (hash.verdict == "N_KEEP") verdict_rank = 2;
    else if (hash.verdict == "MISS") verdict_rank = 1;
    const unsigned offset_rank = hash.offset == 0 ? 1 : 0;
    return static_cast<std::uint16_t>(priority * 10000 + verdict_rank * 1000 +
                                      offset_rank * 100 + std::min(mapq, 99u));
}

struct MoleculeReportKey {
    std::string outcome;
    std::uint32_t cell_class = 0;
    std::uint32_t quant_gene = 0;
    std::uint32_t probe = 0;
    std::uint32_t hash_verdict = 0;
    std::uint32_t hash_gene = 0;
    std::uint32_t genomic_gene_ids = 0;
    std::uint32_t genomic_gene_names = 0;
    std::uint32_t cr_region = 0;
    std::uint8_t probe_region = 0;
    std::uint32_t rname = 0;
    std::uint32_t position_mb = 0;
    std::uint32_t cigar = 0;
    std::uint16_t mapq = 0;
    std::int8_t hash_offset = 0;
    bool used_ur = false;

    bool operator<(const MoleculeReportKey& other) const {
        return std::tie(outcome, cell_class, quant_gene, probe, hash_verdict, hash_gene,
                        genomic_gene_ids, genomic_gene_names, cr_region, probe_region,
                        rname, position_mb, cigar, mapq, hash_offset, used_ur) <
               std::tie(other.outcome, other.cell_class, other.quant_gene, other.probe,
                        other.hash_verdict, other.hash_gene, other.genomic_gene_ids,
                        other.genomic_gene_names, other.cr_region, other.probe_region,
                        other.rname, other.position_mb, other.cigar, other.mapq,
                        other.hash_offset, other.used_ur);
    }
};

struct MoleculeReportStats {
    std::uint64_t molecules = 0;
    std::uint64_t records = 0;
    std::uint64_t xf8_records = 0;
};

std::string molecule_state_outcome(const MoleculeState& state,
                                   const StringInterner& strings) {
    if (state.best.priority == 4) {
        return state.best_xf8_priority == 4 ? "same_gene_accept" : "rescued_same_gene";
    }
    if (state.best.priority == 3) {
        return state.best_xf8_priority == 3 ? "different_gene_conflict" :
                                              "rescued_different_gene_conflict";
    }
    if (state.best.priority == 2) return "cache_miss";
    const std::string& verdict = strings.get(state.best.hash_verdict);
    if (verdict == "DENY_CACHE") return "cache_deny";
    if (verdict == "DENY_GENE_CONFLICT") return "cache_nearby_gene_conflict";
    if (verdict == "DENY_SAMPLE_MISMATCH") return "cache_sample_mismatch";
    return "cache_other_" + verdict;
}

std::string tsv_number(double value) {
    if (!std::isfinite(value)) return "NA";
    std::ostringstream out;
    out << std::fixed << std::setprecision(3) << value;
    return out.str();
}

void write_cell_sets(const std::string& path,
                     const std::unordered_set<std::string>& star,
                     const std::unordered_set<std::string>& cr) {
    std::vector<std::string> cells;
    cells.reserve(star.size() + cr.size());
    cells.insert(cells.end(), star.begin(), star.end());
    for (const auto& cell : cr) if (star.count(cell) == 0) cells.push_back(cell);
    std::sort(cells.begin(), cells.end());
    std::ofstream out(path);
    if (!out) fail("cannot write " + path);
    out << "cell_key\tcell_class\tstar_called\tcr_called\n";
    for (const auto& cell : cells) {
        const bool in_star = star.count(cell) != 0;
        const bool in_cr = cr.count(cell) != 0;
        out << cell << '\t' << (in_star && in_cr ? "shared" : in_star ? "star_only" : "cr_only")
            << '\t' << in_star << '\t' << in_cr << '\n';
    }
}

int run(const Options& options) {
    const TagMap tag_map = load_tag_map(options.tag_map);
    const auto star_cells = load_cells(options.star_cells, tag_map, "STAR");
    const auto cr_cells = load_cells(options.cr_cells, tag_map, "Cell Ranger");
    write_cell_sets(options.output_prefix + ".cells.tsv", star_cells, cr_cells);

    HashCacheView cache;
    if (!options.cache.empty()) cache.open(options.cache);
    const auto genes = load_gene_list(options.gene_list);
    if (cache.enabled() && genes.empty()) fail("gene list is empty");

    std::string fallback_tag = options.sample_tag;
    if (!fallback_tag.empty()) fallback_tag = tag_map.id_for(fallback_tag);

    std::unique_ptr<LineInput> input;
    if (options.input_kind == "bam") {
        std::string command = shell_quote(options.samtools) + " view -h -@ " +
            std::to_string(options.samtools_threads) + " " + shell_quote(options.input);
        if (!options.region.empty()) command += " " + shell_quote(options.region);
        input.reset(new LineInput("", false, command));
    } else {
        input.reset(new LineInput(options.input, ends_with(options.input, ".gz")));
    }

    std::map<SummaryKey, Stats> summary;
    std::map<GeneKey, std::uint64_t> gene_summary;
    std::map<CigarKey, std::uint64_t> cigar_summary;
    std::map<RegionKey, std::uint64_t> region_summary;
    std::map<CountedRecordKey, std::uint64_t> counted_record_summary;
    StringInterner molecule_strings;
    MoleculeTable molecules;
    DetailWriter details(options.output_prefix + ".details.tsv", options.details, options.detail_max);

    std::uint64_t sam_lines = 0;
    std::uint64_t included = 0;
    std::uint64_t skipped_secondary = 0;
    std::uint64_t malformed = 0;
    std::uint64_t cr_xf8_records = 0;
    std::uint64_t molecule_invalid_cell_records = 0;
    std::uint64_t molecule_missing_ub_records = 0;
    std::uint64_t molecule_missing_umi_records = 0;
    std::uint64_t molecule_invalid_umi_records = 0;
    std::uint64_t molecule_ur_fallback_records = 0;
    std::uint64_t molecule_missing_fx_records = 0;
    std::string line;
    while (input->getline(line)) {
        if (line.empty() || line[0] == '@') continue;
        ++sam_lines;
        const auto fields = split_tabs(line);
        if (fields.size() < 11) {
            ++malformed;
            continue;
        }
        unsigned flag = 0, position = 0, mapq = 0;
        if (!parse_integer(fields[1], flag) || !parse_integer(fields[3], position) ||
            !parse_integer(fields[4], mapq)) {
            ++malformed;
            continue;
        }
        const bool primary = (flag & (0x100u | 0x800u)) == 0;
        if (!options.include_secondary && !primary) {
            ++skipped_secondary;
            continue;
        }
        ++included;

        std::string sequence(fields[9]);
        std::string quality(fields[10]);
        if (flag & 0x10u) {
            sequence = reverse_complement(sequence);
            quality = reverse_string(quality);
        }
        const std::string cb = get_tag_value(fields, "CB");
        const std::string raw_cb = get_tag_value(fields, "CR");
        const std::string key = cell_key(cb.empty() ? raw_cb : cb, fallback_tag, tag_map);
        std::string record_tag = fallback_tag;
        const std::size_t key_pipe = key.find('|');
        if (key_pipe != std::string::npos) record_tag = key.substr(key_pipe + 1);
        const std::uint16_t sample_index = tag_map.index_for(record_tag);
        const HashDecision hash = cache.classify(sequence, sample_index, options.single_n);
        const std::string hash_gene = gene_id(hash.gene_index, genes);
        const std::string membership = cell_class(key, star_cells, cr_cells);
        const std::string gx = get_tag_value(fields, "GX");
        const std::string gn = get_tag_value(fields, "GN");
        const std::string fx_gene = get_tag_value(fields, "fx");
        const auto cr_gene_set = gene_set(fx_gene);
        const std::string cr_genes = join(cr_gene_set, ';');
        const std::string relation = gene_relation(hash_gene, cr_gene_set);
        std::string cr_region = get_tag_value(fields, "RE");
        if (cr_region.empty()) cr_region = "other";
        const std::string probe_region = probe_region_name(hash.probe_region);

        const std::string xf_text = get_tag_value(fields, "xf");
        unsigned xf = 0;
        if (!xf_text.empty() && !parse_integer(std::string_view(xf_text), xf)) xf = 0;
        const bool counted = (xf & 8u) != 0;
        cr_xf8_records += counted && primary;
        const bool mapped = (flag & 0x4u) == 0;
        const bool duplicate = (flag & 0x400u) != 0;
        const double probe_q = mean_quality(quality, 50);
        const std::string umi_quality = get_tag_value(fields, "UY");
        const int umi_min_q = min_quality(umi_quality);

        SummaryKey summary_key{membership, hash.verdict, counted, relation};
        summary[summary_key].add(mapped, duplicate, mapq, probe_q, umi_min_q);
        ++gene_summary[{membership, hash_gene, cr_genes, relation, counted}];
        ++cigar_summary[{membership, hash.verdict, counted, relation,
                         std::string(fields[5]), mapq_bin(mapq)}];
        ++region_summary[{membership, cr_region, probe_region, hash.verdict, hash_gene,
                          cr_genes, gx, gn, relation, counted}];

        std::string probe = get_tag_value(fields, "pr");
        if (probe.empty()) probe = get_tag_value(fields, "PR");
        if (probe.empty()) probe = get_tag_value(fields, "PB");
        const std::string ub = get_tag_value(fields, "UB");
        const std::string ur = get_tag_value(fields, "UR");
        if (counted) {
            ++counted_record_summary[{counted_record_outcome(hash, relation), membership,
                                      cr_genes, probe, hash_gene, gx, gn, cr_region,
                                      probe_region}];
        }

        if (primary) {
            const bool missing_ub = ub.empty() || ub == "*";
            molecule_missing_ub_records += missing_ub;
            std::uint64_t packed_cell = 0;
            if (!pack_composite_cell(key, tag_map, packed_cell)) {
                ++molecule_invalid_cell_records;
            } else if (cr_genes.empty()) {
                ++molecule_missing_fx_records;
            } else {
                const bool used_ur = missing_ub;
                const std::string& molecule_umi = used_ur ? ur : ub;
                if (molecule_umi.empty() || molecule_umi == "*") {
                    ++molecule_missing_umi_records;
                } else {
                    std::uint64_t packed_umi = 0;
                    if (!pack_acgt(molecule_umi, packed_umi)) {
                        ++molecule_invalid_umi_records;
                    } else {
                        molecule_ur_fallback_records += used_ur;
                        const PackedMoleculeKey molecule_key{
                            packed_cell, packed_umi, molecule_strings.intern(cr_genes),
                            static_cast<std::uint8_t>(molecule_umi.size())};
                        MoleculeState& state = molecules.get_or_insert(molecule_key);
                        const bool first_record = state.records == 0;
                        ++state.records;
                        state.xf8_records += counted;
                        state.used_ur_fallback = state.used_ur_fallback || used_ur;
                        if (first_record) state.cell_class = molecule_strings.intern(membership);

                        MoleculeEvidence evidence;
                        evidence.hash_verdict = molecule_strings.intern(hash.verdict);
                        evidence.hash_gene = molecule_strings.intern(hash_gene);
                        evidence.probe = molecule_strings.intern(probe);
                        evidence.genomic_gene_ids = molecule_strings.intern(gx);
                        evidence.genomic_gene_names = molecule_strings.intern(gn);
                        evidence.cr_region = molecule_strings.intern(cr_region);
                        evidence.rname = molecule_strings.intern(std::string(fields[2]));
                        evidence.cigar = molecule_strings.intern(std::string(fields[5]));
                        evidence.position = position;
                        evidence.mapq = static_cast<std::uint16_t>(std::min(mapq, 65535u));
                        evidence.hash_offset = static_cast<std::int8_t>(hash.offset);
                        evidence.probe_region = hash.probe_region;
                        evidence.priority = evidence_priority(hash, relation);
                        evidence.selection_score = evidence_score(
                            hash, evidence.priority, mapq == 255 ? 0 : mapq);
                        if (first_record || evidence.selection_score > state.best.selection_score) {
                            state.best = evidence;
                        }
                        if (counted) {
                            state.best_xf8_priority = std::max(state.best_xf8_priority,
                                                               evidence.priority);
                        }
                    }
                }
            }
        }

        if (details.wants(membership, hash, counted, relation)) {
            const std::string cy = get_tag_value(fields, "CY");
            std::ostringstream row;
            row << fields[0] << '\t' << key << '\t' << membership << '\t' << counted
                << '\t' << hash.verdict << '\t' << static_cast<unsigned>(hash.cache_class)
                << '\t' << static_cast<unsigned>(hash.negative_code) << '\t' << hash.offset
                << '\t' << sample_index << '\t' << probe_region << '\t' << cr_region
                << '\t' << hash_gene << '\t' << cr_genes << '\t' << gx << '\t' << gn
                << '\t' << relation
                << '\t' << flag << '\t' << fields[2]
                << '\t' << position << '\t' << mapq << '\t' << mapq_bin(mapq) << '\t' << fields[5]
                << '\t' << xf_text << '\t' << cb << '\t' << raw_cb << '\t' << ub << '\t' << ur
                << '\t' << probe << '\t' << tsv_number(probe_q)
                << '\t' << tsv_number(mean_quality(cy))
                << '\t' << tsv_number(mean_quality(umi_quality)) << '\t' << umi_min_q
                << '\t' << all_tags(fields);
            details.add(fields[0], row.str());
        }
        if (options.limit != 0 && included >= options.limit) break;
    }
    input->close_checked();
    details.finish();

    std::map<MoleculeReportKey, MoleculeReportStats> molecule_report;
    std::map<std::string, std::uint64_t> molecule_outcome_totals;
    std::uint64_t unique_counted_molecules = 0;
    std::uint64_t counted_molecules_multiple_xf8 = 0;
    std::uint64_t counted_molecules_ur_fallback = 0;
    std::uint64_t valid_xf8_records_in_molecules = 0;
    for (const MoleculeSlot& slot : molecules.slots()) {
        if (slot.hash == 0 || slot.state.xf8_records == 0) continue;
        const MoleculeState& state = slot.state;
        ++unique_counted_molecules;
        counted_molecules_multiple_xf8 += state.xf8_records > 1;
        counted_molecules_ur_fallback += state.used_ur_fallback;
        valid_xf8_records_in_molecules += state.xf8_records;
        const std::string outcome = molecule_state_outcome(state, molecule_strings);
        ++molecule_outcome_totals[outcome];
        const std::uint32_t position_mb = state.best.position == 0
            ? std::numeric_limits<std::uint32_t>::max()
            : (state.best.position - 1) / 1000000;
        const MoleculeReportKey report_key{
            outcome, state.cell_class, slot.key.quant_gene, state.best.probe,
            state.best.hash_verdict, state.best.hash_gene, state.best.genomic_gene_ids,
            state.best.genomic_gene_names, state.best.cr_region, state.best.probe_region,
            state.best.rname, position_mb, state.best.cigar, state.best.mapq,
            state.best.hash_offset, state.used_ur_fallback};
        MoleculeReportStats& report = molecule_report[report_key];
        ++report.molecules;
        report.records += state.records;
        report.xf8_records += state.xf8_records;
    }

    {
        const std::string path = options.output_prefix + ".summary.tsv";
        std::ofstream out(path);
        if (!out) fail("cannot write " + path);
        out << "cell_class\thash_verdict\tcr_counted\tgene_relation\trecords\tmapped"
               "\tduplicate\tmean_mapq\tmean_probe_q\tumi_q_lt20\tumi_q_observed\n";
        for (const auto& item : summary) {
            const auto& key = item.first;
            const auto& stats = item.second;
            out << key.cell << '\t' << key.hash << '\t' << key.counted << '\t' << key.gene
                << '\t' << stats.records << '\t' << stats.mapped << '\t' << stats.duplicate
                << '\t' << tsv_number(stats.records == 0 ? NAN :
                                      static_cast<double>(stats.mapq_sum) / stats.records)
                << '\t' << tsv_number(stats.probe_quality_n == 0 ? NAN :
                                      static_cast<double>(stats.probe_quality_sum / stats.probe_quality_n))
                << '\t' << stats.umi_q_lt20 << '\t' << stats.umi_quality_n << '\n';
        }
    }
    {
        const std::string path = options.output_prefix + ".genes.tsv";
        std::ofstream out(path);
        if (!out) fail("cannot write " + path);
        out << "cell_class\thash_gene\tquant_gene_ids\tgene_relation\tcr_counted\trecords\n";
        for (const auto& item : gene_summary) {
            out << item.first.cell << '\t' << item.first.hash_gene << '\t' << item.first.cr_genes
                << '\t' << item.first.relation << '\t' << item.first.counted << '\t' << item.second << '\n';
        }
    }
    {
        const std::string path = options.output_prefix + ".cigars.tsv";
        std::ofstream out(path);
        if (!out) fail("cannot write " + path);
        out << "cell_class\thash_verdict\tcr_counted\tgene_relation\tcigar\tmapq_bin\trecords\n";
        for (const auto& item : cigar_summary) {
            out << item.first.cell << '\t' << item.first.hash << '\t' << item.first.counted
                << '\t' << item.first.relation << '\t' << item.first.cigar << '\t'
                << item.first.mapq << '\t' << item.second << '\n';
        }
    }
    {
        const std::string path = options.output_prefix + ".regions.tsv";
        std::ofstream out(path);
        if (!out) fail("cannot write " + path);
        out << "cell_class\tcr_region\tprobe_region\thash_verdict\thash_gene"
               "\tquant_gene_ids\tgenomic_gene_ids\tgenomic_gene_names"
               "\tgene_relation\tcr_counted\trecords\n";
        for (const auto& item : region_summary) {
            const auto& key = item.first;
            out << key.cell << '\t' << key.cr_region << '\t' << key.probe_region << '\t'
                << key.hash << '\t' << key.hash_gene << '\t' << key.quant_genes << '\t'
                << key.genomic_genes << '\t' << key.genomic_names << '\t' << key.relation
                << '\t' << key.counted << '\t' << item.second << '\n';
        }
    }
    {
        std::vector<std::pair<CountedRecordKey, std::uint64_t>> rows(
            counted_record_summary.begin(), counted_record_summary.end());
        std::sort(rows.begin(), rows.end(), [](const auto& lhs, const auto& rhs) {
            if (lhs.second != rhs.second) return lhs.second > rhs.second;
            return lhs.first < rhs.first;
        });
        const std::string path = options.output_prefix + ".counted_records.tsv";
        std::ofstream out(path);
        if (!out) fail("cannot write " + path);
        out << "rank\toutcome\tcell_class\tquant_gene_ids\tprobe_tag\thash_gene"
               "\tgenomic_gene_ids\tgenomic_gene_names\tcr_region\tprobe_region\trecords\n";
        std::uint64_t rank = 0;
        for (const auto& item : rows) {
            const auto& key = item.first;
            out << ++rank << '\t' << key.outcome << '\t' << key.cell << '\t'
                << key.quant_gene_ids << '\t' << key.probe << '\t' << key.hash_gene << '\t'
                << key.genomic_gene_ids << '\t' << key.genomic_gene_names << '\t'
                << key.cr_region << '\t' << key.probe_region << '\t'
                << item.second << '\n';
        }
    }
    {
        std::vector<std::pair<MoleculeReportKey, MoleculeReportStats>> rows(
            molecule_report.begin(), molecule_report.end());
        std::sort(rows.begin(), rows.end(), [](const auto& lhs, const auto& rhs) {
            if (lhs.second.molecules != rhs.second.molecules) {
                return lhs.second.molecules > rhs.second.molecules;
            }
            return lhs.first < rhs.first;
        });
        const std::string path = options.output_prefix + ".molecules.tsv";
        std::ofstream out(path);
        if (!out) fail("cannot write " + path);
        out << "rank\toutcome\tcell_class\tquant_gene_ids\tprobe_tag\thash_verdict"
               "\thash_gene\tgenomic_gene_ids\tgenomic_gene_names\tcr_region"
               "\tprobe_region\trname\tposition_1mb_bin\tcigar"
               "\tmapq\tmapq_bin\thash_offset\tumi_source\tmolecules"
               "\tinput_records\txf8_records\n";
        std::uint64_t rank = 0;
        for (const auto& item : rows) {
            const MoleculeReportKey& key = item.first;
            const MoleculeReportStats& stats = item.second;
            std::string position_bin = "unmapped";
            if (key.position_mb != std::numeric_limits<std::uint32_t>::max()) {
                const std::uint64_t start = static_cast<std::uint64_t>(key.position_mb) *
                                                UINT64_C(1000000) + 1;
                position_bin = std::to_string(start) + "-" +
                               std::to_string(start + UINT64_C(999999));
            }
            out << ++rank << '\t' << key.outcome << '\t'
                << molecule_strings.get(key.cell_class) << '\t'
                << molecule_strings.get(key.quant_gene) << '\t'
                << molecule_strings.get(key.probe) << '\t'
                << molecule_strings.get(key.hash_verdict) << '\t'
                << molecule_strings.get(key.hash_gene) << '\t'
                << molecule_strings.get(key.genomic_gene_ids) << '\t'
                << molecule_strings.get(key.genomic_gene_names) << '\t'
                << molecule_strings.get(key.cr_region) << '\t'
                << probe_region_name(key.probe_region) << '\t'
                << molecule_strings.get(key.rname) << '\t' << position_bin << '\t'
                << molecule_strings.get(key.cigar) << '\t' << key.mapq << '\t'
                << mapq_bin(key.mapq) << '\t' << static_cast<int>(key.hash_offset) << '\t'
                << (key.used_ur ? "UR_fallback" : "UB") << '\t'
                << stats.molecules << '\t' << stats.records << '\t'
                << stats.xf8_records << '\n';
        }
    }
    {
        std::uint64_t shared = 0;
        for (const auto& cell : star_cells) shared += cr_cells.count(cell) != 0;
        const std::uint64_t star_only = star_cells.size() - shared;
        const std::uint64_t cr_only = cr_cells.size() - shared;
        const std::uint64_t union_size = shared + star_only + cr_only;
        const std::string path = options.output_prefix + ".metrics.tsv";
        std::ofstream out(path);
        if (!out) fail("cannot write " + path);
        out << "metric\tvalue\n"
            << "star_cells\t" << star_cells.size() << '\n'
            << "cr_cells\t" << cr_cells.size() << '\n'
            << "shared_cells\t" << shared << '\n'
            << "star_only_cells\t" << star_only << '\n'
            << "cr_only_cells\t" << cr_only << '\n'
            << "cell_precision\t" << std::setprecision(12)
            << (star_cells.empty() ? 0.0 : static_cast<double>(shared) / star_cells.size()) << '\n'
            << "cell_recall\t" << (cr_cells.empty() ? 0.0 : static_cast<double>(shared) / cr_cells.size()) << '\n'
            << "cell_jaccard\t" << (union_size == 0 ? 0.0 : static_cast<double>(shared) / union_size) << '\n'
            << "sam_alignment_lines\t" << sam_lines << '\n'
            << "included_primary_records\t" << included << '\n'
            << "cr_xf8_primary_records\t" << cr_xf8_records << '\n'
            << "molecule_keys_all_valid_records\t" << molecules.size() << '\n'
            << "unique_counted_molecule_keys\t" << unique_counted_molecules << '\n'
            << "valid_xf8_records_in_molecule_keys\t" << valid_xf8_records_in_molecules << '\n'
            << "counted_molecule_keys_multiple_xf8\t" << counted_molecules_multiple_xf8 << '\n'
            << "counted_molecule_keys_ur_fallback\t" << counted_molecules_ur_fallback << '\n'
            << "molecule_invalid_cell_records\t" << molecule_invalid_cell_records << '\n'
            << "molecule_missing_ub_records\t" << molecule_missing_ub_records << '\n'
            << "molecule_missing_umi_records\t" << molecule_missing_umi_records << '\n'
            << "molecule_invalid_umi_records\t" << molecule_invalid_umi_records << '\n'
            << "molecule_ur_fallback_records\t" << molecule_ur_fallback_records << '\n'
            << "molecule_missing_fx_records\t" << molecule_missing_fx_records << '\n'
            << "molecule_table_capacity\t" << molecules.capacity() << '\n'
            << "molecule_table_bytes\t" << molecules.capacity() * sizeof(MoleculeSlot) << '\n'
            << "skipped_secondary_or_supplementary\t" << skipped_secondary << '\n'
            << "malformed_sam_lines\t" << malformed << '\n'
            << "cache_version\t" << cache.version() << '\n'
            << "cache_records\t" << cache.record_count() << '\n'
            << "detail_eligible\t" << details.eligible() << '\n'
            << "detail_written\t" << details.written() << '\n';
        for (const auto& outcome : molecule_outcome_totals) {
            out << "molecule_outcome_" << outcome.first << '\t' << outcome.second << '\n';
        }
    }
    return malformed == 0 ? 0 : 2;
}

} // namespace

int main(int argc, char** argv) {
    try {
        return run(parse_options(argc, argv));
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
