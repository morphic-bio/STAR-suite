// Deterministic read-level comparison of STAR's authoritative fused-Flex
// decision ledger with Cell Ranger BAM evidence. BAM decoding is delegated to
// samtools; all grouping, inclusion policy, and comparison logic is C++17.

#include <algorithm>
#include <cerrno>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

[[noreturn]] void fail(const std::string &message) {
    throw std::runtime_error(message);
}

std::string shell_quote(const std::string &value) {
    std::string out("'");
    for (char c : value) out += c == '\'' ? "'\\''" : std::string(1, c);
    out += "'";
    return out;
}

std::vector<std::string_view> split_tabs(const std::string &line) {
    std::vector<std::string_view> out;
    std::size_t begin = 0;
    while (true) {
        const std::size_t end = line.find('\t', begin);
        out.emplace_back(line.data() + begin,
                         (end == std::string::npos ? line.size() : end) - begin);
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return out;
}

std::vector<std::string> split_csv(const std::string &line) {
    std::vector<std::string> out;
    std::string field;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                field.push_back('"');
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (c == ',' && !quoted) {
            out.push_back(field);
            field.clear();
        } else {
            field.push_back(c);
        }
    }
    out.push_back(field);
    return out;
}

template <typename T> T parse_number(std::string_view text, const char *what) {
    T value{};
    const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
    if (result.ec != std::errc() || result.ptr != text.data() + text.size())
        fail(std::string("invalid ") + what + ": " + std::string(text));
    return value;
}

std::string value_or_dot(const std::string &s) { return s.empty() ? "." : s; }

class LineInput {
public:
    LineInput(const std::string &path, const std::string &command = {}) {
        if (!command.empty()) {
            fp_ = ::popen(command.c_str(), "r");
            pipe_ = true;
            description_ = command;
        } else {
            fp_ = std::fopen(path.c_str(), "rb");
            description_ = path;
        }
        if (fp_ == nullptr) fail("cannot open " + description_ + ": " + std::strerror(errno));
    }
    ~LineInput() { close_noexcept(); std::free(buffer_); }
    bool getline(std::string &line) {
        const ssize_t n = ::getline(&buffer_, &capacity_, fp_);
        if (n < 0) {
            if (std::ferror(fp_)) fail("read error in " + description_);
            return false;
        }
        std::size_t size = static_cast<std::size_t>(n);
        while (size && (buffer_[size - 1] == '\n' || buffer_[size - 1] == '\r')) --size;
        line.assign(buffer_, size);
        return true;
    }
    void close_checked() {
        if (closed_) return;
        const int status = pipe_ ? ::pclose(fp_) : std::fclose(fp_);
        fp_ = nullptr;
        closed_ = true;
        if (status != 0) fail("input command failed: " + description_);
    }
private:
    void close_noexcept() {
        if (closed_ || fp_ == nullptr) return;
        if (pipe_) ::pclose(fp_); else std::fclose(fp_);
        fp_ = nullptr;
        closed_ = true;
    }
    FILE *fp_ = nullptr;
    char *buffer_ = nullptr;
    std::size_t capacity_ = 0;
    bool pipe_ = false;
    bool closed_ = false;
    std::string description_;
};

struct NamedBam { std::string label, path; };
struct Options {
    std::string ledger, star_bam, probe_set, probe_list, tag_map, reference, out_dir;
    std::string samtools = "samtools";
    unsigned threads = 4;
    std::uint64_t limit = 0;
    std::vector<NamedBam> cr_bams;
};

void usage() {
    std::cerr <<
      "Usage: flex_decision_ledger_compare --star-ledger FILE [--star-bam FILE] --cr-bam LABEL=FILE ...\n"
      "  --probe-set FILE.csv --probe-list FILE --tag-map FILE --reference genome.fa --out-dir DIR\n"
      "  [--samtools PATH] [--samtools-threads N] [--limit N]\n";
}

Options parse_options(int argc, char **argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&]() -> std::string {
            if (++i >= argc) fail("missing value after " + arg);
            return argv[i];
        };
        if (arg == "--star-ledger") o.ledger = value();
        else if (arg == "--star-bam") o.star_bam = value();
        else if (arg == "--probe-set") o.probe_set = value();
        else if (arg == "--probe-list") o.probe_list = value();
        else if (arg == "--tag-map") o.tag_map = value();
        else if (arg == "--reference") o.reference = value();
        else if (arg == "--out-dir") o.out_dir = value();
        else if (arg == "--samtools") o.samtools = value();
        else if (arg == "--samtools-threads") o.threads = std::stoul(value());
        else if (arg == "--limit") o.limit = std::stoull(value());
        else if (arg == "--cr-bam") {
            const std::string spec = value();
            const std::size_t pos = spec.find('=');
            if (pos == 0 || pos == std::string::npos || pos + 1 == spec.size())
                fail("--cr-bam requires LABEL=FILE");
            o.cr_bams.push_back({spec.substr(0, pos), spec.substr(pos + 1)});
        } else if (arg == "--help" || arg == "-h") { usage(); std::exit(0); }
        else fail("unknown option: " + arg);
    }
    if (o.ledger.empty() || o.cr_bams.empty() ||
        o.probe_set.empty() || o.probe_list.empty() || o.tag_map.empty() ||
        o.reference.empty() ||
        o.out_dir.empty()) fail("missing required option (use --help)");
    if (o.threads == 0) fail("--samtools-threads must be positive");
    return o;
}

struct ProbeInfo {
    std::string gene;
    bool included = false;
    bool deprecated = false;
};

std::unordered_map<std::string, ProbeInfo> load_probes(const std::string &path) {
    LineInput in(path);
    std::string line;
    while (in.getline(line) && (line.empty() || line[0] == '#')) {}
    if (line.empty() || line[0] == '#') fail("probe CSV has no header: " + path);
    const auto header = split_csv(line);
    std::unordered_map<std::string, std::size_t> col;
    for (std::size_t i = 0; i < header.size(); ++i) col[header[i]] = i;
    for (const char *required : {"gene_id", "probe_id", "included"})
        if (!col.count(required)) fail(std::string("probe CSV missing column ") + required);
    std::unordered_map<std::string, ProbeInfo> probes;
    while (in.getline(line)) {
        const auto f = split_csv(line);
        const std::size_t need = std::max({col["gene_id"], col["probe_id"], col["included"]});
        if (f.size() <= need) continue;
        ProbeInfo p;
        p.gene = f[col["gene_id"]];
        std::string yes = f[col["included"]];
        std::transform(yes.begin(), yes.end(), yes.begin(), ::toupper);
        p.included = yes == "TRUE" || yes == "1" || yes == "YES";
        p.deprecated = p.gene.rfind("DEPRECATED_", 0) == 0;
        probes[f[col["probe_id"]]] = std::move(p);
    }
    in.close_checked();
    return probes;
}

std::vector<std::string> load_lines(const std::string &path) {
    LineInput in(path);
    std::vector<std::string> out;
    std::string line;
    while (in.getline(line)) if (!line.empty() && line[0] != '#') out.push_back(line);
    in.close_checked();
    return out;
}

std::vector<std::string> load_tag_sequences(const std::string &path) {
    LineInput in(path);
    std::vector<std::string> result(1);
    std::string line;
    while (in.getline(line)) {
        const auto f = split_tabs(line);
        if (f.size() < 2) continue;
        std::string label(f[0]), seq(f[1]);
        std::uint32_t index = 0;
        if (label.rfind("BC", 0) == 0) {
            try { index = static_cast<std::uint32_t>(std::stoul(label.substr(2))); }
            catch (...) { index = 0; }
        }
        // Production Flex whitelists normally use a biological sample label
        // in column 1, repeated for its tags. In that form the runtime token
        // is the 1-based row number.
        if (index == 0) index = static_cast<std::uint32_t>(result.size());
        if (result.size() <= index) result.resize(index + 1);
        result[index] = seq;
    }
    in.close_checked();
    return result;
}

struct Resolution { bool resolved = false; std::uint32_t cb_idx = 0; std::string cb; };
struct StarDecision {
    std::uint64_t iread = 0, lane_ordinal = 0, molecule_key = 0, ambiguous_key = 0;
    std::uint32_t lane = 0, cache_gene = 0, gene_idx = 0, umi_packed = 0;
    int cache_class = -1, negative_code = -1, hash_offset = 0, sample_ok = -1,
        tag_token = -1, cb_match = -999, umi_valid = -1,
        probe_hamming_distance = -1;
    bool have_molecule_key = false, have_ambiguous_key = false;
    unsigned triage_count = 0, final_count = 0;
    std::string qname, triage, cb_raw, cb_corrected, cb_candidates, umi_raw,
        source, state, reason, resolver_reason, detail;
};

std::uint64_t parse_hex(std::string_view text) {
    if (text == "." || text.empty()) return 0;
    if (text.size() > 2 && text.substr(0, 2) == "0x") text.remove_prefix(2);
    std::uint64_t value = 0;
    const auto r = std::from_chars(text.data(), text.data() + text.size(), value, 16);
    if (r.ec != std::errc() || r.ptr != text.data() + text.size()) fail("bad hex value");
    return value;
}

std::string detail_value(const std::string &detail, const std::string &key) {
    const std::string needle = key + "=";
    std::size_t pos = detail.find(needle);
    if (pos == std::string::npos) return {};
    pos += needle.size();
    const std::size_t end = detail.find(';', pos);
    return detail.substr(pos, end == std::string::npos ? std::string::npos : end - pos);
}

struct LedgerData {
    std::unordered_map<std::uint64_t, StarDecision> by_id;
    std::unordered_map<std::uint64_t, Resolution> resolutions;
};

LedgerData load_ledger(const std::string &path) {
    LineInput in(path);
    std::string line;
    if (!in.getline(line)) fail("empty STAR decision ledger");
    const auto header = split_tabs(line);
    std::unordered_map<std::string, std::size_t> c;
    for (std::size_t i = 0; i < header.size(); ++i) c[std::string(header[i])] = i;
    auto require = [&](const char *name) { if (!c.count(name)) fail(std::string("ledger missing ") + name); };
    for (const char *name : {"event","iread","qname","lane","lane_ordinal","triage_action",
         "cache_class","cache_gene_idx","negative_code","hash_offset","sample_ok","tag_token",
         "cb_raw","cb_corrected","cb_match","cb_candidate_indices","umi_raw","umi_packed",
         "umi_valid","chosen_source","chosen_gene_idx","final_state","final_reason",
         "molecule_key_hex","ambiguous_key_hex","detail","probe_hamming_distance"}) require(name);
    LedgerData data;
    std::uint64_t rows = 0;
    while (in.getline(line)) {
        const auto f = split_tabs(line);
        if (f.size() != header.size()) fail("malformed ledger row " + std::to_string(rows + 2));
        ++rows;
        auto s = [&](const char *name) { return std::string(f[c[name]]); };
        const std::string event = s("event");
        if (event == "AMBIGUOUS_RESOLUTION") {
            const std::uint64_t key = parse_hex(f[c["ambiguous_key_hex"]]);
            Resolution r;
            const std::string resolution_state = s("final_state");
            r.resolved = resolution_state == "KEEP" || resolution_state == "RESOLVED";
            const std::string idx = detail_value(s("detail"), "corrected_cb_index");
            if (!idx.empty()) r.cb_idx = static_cast<std::uint32_t>(std::stoul(idx));
            r.cb = detail_value(s("detail"), "corrected_cb");
            data.resolutions[key] = std::move(r);
            continue;
        }
        const std::uint64_t id = parse_number<std::uint64_t>(f[c["iread"]], "iread");
        StarDecision &d = data.by_id[id];
        d.iread = id;
        if (event == "TRIAGE") {
            ++d.triage_count;
            d.qname = s("qname");
            // readLoad retains the FASTQ marker in the reader-side name used
            // by the diagnostic hook; SAM/BAM and raw FASTQ comparison keys
            // conventionally omit it.
            if (!d.qname.empty() && d.qname[0] == '@') d.qname.erase(0, 1);
            d.lane = static_cast<std::uint32_t>(std::stoul(s("lane")));
            d.lane_ordinal = std::stoull(s("lane_ordinal"));
            d.triage = s("triage_action");
            d.cache_class = std::stoi(s("cache_class"));
            d.cache_gene = static_cast<std::uint32_t>(std::stoul(s("cache_gene_idx")));
            d.negative_code = std::stoi(s("negative_code"));
            d.hash_offset = std::stoi(s("hash_offset"));
            d.sample_ok = std::stoi(s("sample_ok"));
            d.tag_token = std::stoi(s("tag_token"));
            const std::string probe_distance = s("probe_hamming_distance");
            d.probe_hamming_distance = probe_distance == "."
                ? -1 : std::stoi(probe_distance);
        } else {
            if (event == "RESOLVER") d.resolver_reason = s("final_reason");
            if (event == "FINAL") {
                ++d.final_count;
                d.cb_raw = s("cb_raw"); d.cb_corrected = s("cb_corrected");
                d.cb_match = std::stoi(s("cb_match"));
                d.cb_candidates = s("cb_candidate_indices"); d.umi_raw = s("umi_raw");
                d.umi_packed = static_cast<std::uint32_t>(std::stoul(s("umi_packed")));
                d.umi_valid = std::stoi(s("umi_valid")); d.source = s("chosen_source");
                d.gene_idx = static_cast<std::uint32_t>(std::stoul(s("chosen_gene_idx")));
                d.state = s("final_state"); d.reason = s("final_reason"); d.detail = s("detail");
                if (s("molecule_key_hex") != ".") {
                    d.molecule_key = parse_hex(f[c["molecule_key_hex"]]); d.have_molecule_key = true;
                }
                if (s("ambiguous_key_hex") != ".") {
                    d.ambiguous_key = parse_hex(f[c["ambiguous_key_hex"]]); d.have_ambiguous_key = true;
                }
            }
        }
    }
    in.close_checked();
    return data;
}

void add_unique(std::vector<std::string> &values, const std::string &value) {
    if (value.empty()) return;
    if (std::find(values.begin(), values.end(), value) == values.end()) values.push_back(value);
}

std::string join(std::vector<std::string> values, char sep = ',') {
    if (values.empty()) return ".";
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    std::string out;
    for (const auto &v : values) { if (!out.empty()) out.push_back(sep); out += v; }
    return out;
}

struct CigarStats {
    std::uint32_t query = 0, reference = 0, soft5 = 0, soft3 = 0;
    bool indel = false, splice = false;
    std::string shape = "none";
};

CigarStats cigar_stats(const std::string &cigar) {
    CigarStats out;
    if (cigar.empty() || cigar == "*") return out;
    std::vector<std::pair<std::uint32_t,char>> ops;
    std::uint32_t n = 0;
    for (char c : cigar) {
        if (c >= '0' && c <= '9') { n = n * 10 + static_cast<unsigned>(c - '0'); continue; }
        if (n == 0) fail("invalid CIGAR: " + cigar);
        ops.push_back({n, c});
        if (c == 'M' || c == 'I' || c == '=' || c == 'X') out.query += n;
        if (c == 'M' || c == 'D' || c == 'N' || c == '=' || c == 'X') out.reference += n;
        if (c == 'I' || c == 'D') out.indel = true;
        if (c == 'N') out.splice = true;
        n = 0;
    }
    if (n != 0) fail("invalid CIGAR tail: " + cigar);
    if (!ops.empty() && ops.front().second == 'S') out.soft5 = ops.front().first;
    if (!ops.empty() && ops.back().second == 'S') out.soft3 = ops.back().first;
    if (out.splice) out.shape = "spliced";
    else if (out.indel) out.shape = "indel";
    else if (out.soft5 || out.soft3) out.shape = "softclip";
    else out.shape = "contiguous";
    return out;
}

struct SamRepresentative {
    bool set = false;
    int flag = 0, mapq = -1, nm = -1, nh = -1;
    std::string rname = ".", pos = ".", cigar = ".", md = ".";
    CigarStats geometry;
};

struct BamAggregate {
    std::uint32_t records = 0, primary = 0, secondary = 0, supplementary = 0,
                  mapped = 0, unmapped = 0;
    SamRepresentative rep, genomic_rep;
    std::vector<std::string> flags, cb, cr, gr, ub, fx, xf, pr, sources, mapqs, cigars, loci;
};

int representative_rank(int flag, int mapq, const std::string &rname,
                        const std::string &pos, const std::string &cigar) {
    int score = 0;
    if (!(flag & 0x100) && !(flag & 0x800)) score += 1000000;
    if (!(flag & 0x4)) score += 100000;
    score += std::max(0, mapq) * 100;
    score += cigar != "*" ? 1 : 0;
    (void)rname; (void)pos;
    return score;
}

std::string tag_value(const std::vector<std::string_view> &f, std::string_view name) {
    for (std::size_t i = 11; i < f.size(); ++i) {
        if (f[i].size() >= 5 && f[i].substr(0, 2) == name && f[i][2] == ':') {
            const std::size_t second = f[i].find(':', 3);
            if (second != std::string_view::npos) return std::string(f[i].substr(second + 1));
        }
    }
    return {};
}

using BamMap = std::unordered_map<std::string, BamAggregate>;

void load_bam(const std::string &path, const std::string &source, const Options &o,
              BamMap &result) {
    // Add MD tags in a read-only stream so physical-half mismatch burdens can
    // be reconstructed exactly; a global NM count cannot be split safely.
    const std::string command = shell_quote(o.samtools) + " calmd -@ " +
        std::to_string(o.threads) + " -u " + shell_quote(path) + " " +
        shell_quote(o.reference) + " | " + shell_quote(o.samtools) +
        " view -@ " + std::to_string(o.threads) + " -";
    LineInput in("", command);
    std::string line;
    std::uint64_t accepted = 0;
    while (in.getline(line)) {
        if (line.empty() || line[0] == '@') continue;
        const auto f = split_tabs(line);
        if (f.size() < 11) fail("malformed SAM row from " + source);
        const std::string qname(f[0]);
        BamAggregate &a = result[qname];
        const int flag = parse_number<int>(f[1], "SAM flag");
        const int mapq = parse_number<int>(f[4], "SAM MAPQ");
        const bool secondary = flag & 0x100, supplementary = flag & 0x800;
        ++a.records;
        if (!secondary && !supplementary) ++a.primary;
        if (secondary) ++a.secondary;
        if (supplementary) ++a.supplementary;
        if (flag & 0x4) ++a.unmapped; else ++a.mapped;
        add_unique(a.flags, std::to_string(flag)); add_unique(a.mapqs, std::to_string(mapq));
        add_unique(a.cigars, std::string(f[5])); add_unique(a.sources, source);
        if (!(flag & 0x4)) add_unique(a.loci, std::string(f[2]) + ":" + std::string(f[3]));
        add_unique(a.cb, tag_value(f, "CB")); add_unique(a.cr, tag_value(f, "CR"));
        add_unique(a.gr, tag_value(f, "GR")); add_unique(a.ub, tag_value(f, "UB"));
        add_unique(a.fx, tag_value(f, "fx")); add_unique(a.xf, tag_value(f, "xf"));
        add_unique(a.pr, tag_value(f, "pr"));
        auto update_representative = [&](SamRepresentative &rep) {
            const int rank = representative_rank(flag, mapq, std::string(f[2]), std::string(f[3]), std::string(f[5]));
            const int old_rank = rep.set
                ? representative_rank(rep.flag, rep.mapq, rep.rname, rep.pos, rep.cigar) : -1;
            const std::string tie = std::string(f[2]) + "\t" + std::string(f[3]) + "\t" + std::string(f[5]);
            const std::string old_tie = rep.rname + "\t" + rep.pos + "\t" + rep.cigar;
            if (rep.set && (rank < old_rank || (rank == old_rank && tie >= old_tie))) return;
            rep.set = true; rep.flag = flag; rep.mapq = mapq;
            rep.rname = std::string(f[2]); rep.pos = std::string(f[3]);
            rep.cigar = std::string(f[5]); rep.geometry = cigar_stats(rep.cigar);
            const std::string nm = tag_value(f, "NM"), nh = tag_value(f, "NH");
            rep.nm = nm.empty() ? -1 : std::stoi(nm);
            rep.nh = nh.empty() ? -1 : std::stoi(nh);
            const std::string md = tag_value(f, "MD");
            rep.md = md.empty() ? "." : md;
        };
        update_representative(a.rep);
        // Probe pseudo-contigs use the probe ID (and therefore contain '|').
        // Keep a separate genomic representative for the physical-half edit
        // comparison rather than silently treating a probe-contig CIGAR as a
        // competing genomic alignment.
        if (!(flag & 0x4) && std::string_view(f[2]).find('|') == std::string_view::npos)
            update_representative(a.genomic_rep);
        if (o.limit && ++accepted >= o.limit) break;
    }
    in.close_checked();
}

struct CrAssignment {
    bool included = false, excluded = false, deprecated = false, unknown = false;
    bool assigned = false, counted = false;
    std::vector<std::string> genes;
    std::string status = "none";
};

std::vector<std::string> split_list(const std::string &text) {
    std::vector<std::string> out;
    std::string item;
    for (char c : text) {
        if (c == ',' || c == ';') {
            if (!item.empty()) { out.push_back(item); item.clear(); }
        } else item.push_back(c);
    }
    if (!item.empty()) out.push_back(item);
    return out;
}

CrAssignment classify_cr(const BamAggregate *a,
                         const std::unordered_map<std::string, ProbeInfo> &probes,
                         const std::set<std::string> &included_genes) {
    CrAssignment r;
    if (a == nullptr) { r.status = "missing"; return r; }
    for (const std::string &value : a->pr) {
        for (const std::string &probe : split_list(value)) {
            const auto it = probes.find(probe);
            if (it == probes.end()) { r.unknown = true; continue; }
            if (it->second.deprecated) r.deprecated = true;
            else if (it->second.included) { r.included = true; add_unique(r.genes, it->second.gene); }
            else r.excluded = true;
        }
    }
    const bool had_probe_evidence = !a->pr.empty();
    if (!had_probe_evidence) {
        for (const auto &fx : a->fx) {
            for (const auto &gene : split_list(fx)) {
                if (included_genes.count(gene)) { r.included = true; add_unique(r.genes, gene); }
                else r.unknown = true;
            }
        }
    }
    // fx is feature/gene evidence. It is authoritative only after pr has shown
    // that the source probe is included in the 2024 set.
    if (r.included && !a->fx.empty()) {
        r.genes.clear();
        for (const auto &fx : a->fx) for (const auto &gene : split_list(fx)) add_unique(r.genes, gene);
    }
    for (const auto &xf : a->xf) {
        try { if ((std::stoi(xf) & 8) != 0) r.counted = true; } catch (...) {}
    }
    // fx is the read-level feature assignment. xf bit 8 is narrower: it marks
    // the record selected to count after duplicate/UMI policy. Preserve both;
    // STAR's pre-collapse count-record insertion is compared to fx assignment,
    // not to the single xf8 molecule representative.
    r.assigned = r.included && !a->fx.empty();
    r.counted = r.counted && r.assigned;
    if (r.included && (r.excluded || r.deprecated || r.unknown)) r.status = "mixed_included";
    else if (r.included) r.status = "included";
    else if (r.deprecated) r.status = "excluded_deprecated";
    else if (r.excluded) r.status = "excluded";
    else if (r.unknown) r.status = "unknown_probe";
    else r.status = "none";
    return r;
}

std::string gene_id(std::uint32_t idx, const std::vector<std::string> &genes) {
    if (idx == 0) return ".";
    return idx <= genes.size() ? genes[idx - 1] : "#" + std::to_string(idx);
}

std::string compare_category(const StarDecision &s, const CrAssignment &cr) {
    const bool star = s.state == "KEEP";
    const bool cra = cr.assigned;
    if (star && cra) {
        const std::string sg = std::to_string(s.gene_idx);
        (void)sg;
        return "both_assigned"; // split into same/conflict after gene ID mapping
    }
    if (star) {
        if (cr.deprecated || cr.excluded) return "star_only_cr_excluded";
        if (cr.unknown) return "star_only_cr_unknown_probe";
        return "star_only_cr_unassigned";
    }
    if (cra) {
        if (s.reason == "NO_CB_MATCH") return "cr_only_star_no_cb";
        if (s.source == "CACHE_DENY") return "cr_only_star_cache_deny";
        if (s.reason == "RESOLVER_REJECT" || !s.resolver_reason.empty()) return "cr_only_star_resolver_reject";
        if (s.reason == "UNMATCHED_TAG") return "cr_only_star_unmatched_tag";
        return "cr_only_star_other_reject";
    }
    return "neither_assigned";
}

std::string routing_class(const StarDecision &s) {
    if (s.source == "H0" || (s.triage == "KEEP" && s.cache_class == 0)) return "H0";
    if (s.source == "H1") return "H1";
    if (s.source == "H1X2" || (s.triage == "KEEP" && s.cache_class == 4))
        return s.probe_hamming_distance < 0 ? "H1X2_DIRECT" : "H1X2_EXTEND";
    if (s.triage == "KEEP" && s.cache_class != 0) return "H1";
    if (s.source == "CACHE_DENY" || s.triage == "DENY") return "CACHE_DENY";
    if (s.triage == "MISS") return "RESIDUAL_ALIGNMENT";
    return "OTHER";
}

std::uint32_t route_bit(const std::string &route) {
    if (route == "H0") return 1u;
    if (route == "H1") return 2u;
    if (route == "H1X2_DIRECT") return 4u;
    if (route == "H1X2_EXTEND") return 8u;
    if (route == "CACHE_DENY") return 16u;
    if (route == "RESIDUAL_ALIGNMENT") return 32u;
    return 64u;
}

std::string route_mask_name(std::uint32_t mask) {
    std::string out;
    for (const auto &entry : std::vector<std::pair<std::uint32_t,const char*>>{
             {1u,"H0"},{2u,"H1"},{4u,"H1X2_DIRECT"},{8u,"H1X2_EXTEND"},
             {16u,"CACHE_DENY"},{32u,"RESIDUAL_ALIGNMENT"},{64u,"OTHER"}}) {
        if (mask & entry.first) { if (!out.empty()) out += '+'; out += entry.second; }
    }
    return out.empty() ? "NONE" : out;
}

void ensure_dir(const std::string &path) {
    const std::string command = "mkdir -p -- " + shell_quote(path);
    if (std::system(command.c_str()) != 0) fail("cannot create output directory");
}

template <typename Map>
void write_summary(const std::string &path, const std::string &header, const Map &counts) {
    std::ofstream out(path);
    if (!out) fail("cannot write " + path);
    out << header << "\treads\n";
    for (const auto &kv : counts) out << kv.first << '\t' << kv.second << '\n';
}

std::string hex_key(std::uint64_t key, bool present) {
    if (!present) return ".";
    std::ostringstream out;
    out << "0x" << std::hex << std::setw(16) << std::setfill('0') << key;
    return out.str();
}

int run(const Options &o) {
    ensure_dir(o.out_dir);
    auto probes = load_probes(o.probe_set);
    std::set<std::string> included_genes;
    for (const auto &kv : probes) {
        if (kv.second.included && !kv.second.deprecated) included_genes.insert(kv.second.gene);
    }
    auto genes = load_lines(o.probe_list);
    auto tag_sequences = load_tag_sequences(o.tag_map);
    auto ledger = load_ledger(o.ledger);
    for (const auto &kv : ledger.by_id) {
        if (kv.second.triage_count != 1 || kv.second.final_count != 1) {
            fail("STAR ledger does not have exactly one TRIAGE and FINAL for iread=" +
                 std::to_string(kv.first) + " (triage=" +
                 std::to_string(kv.second.triage_count) + ", final=" +
                 std::to_string(kv.second.final_count) + ")");
        }
    }

    // Apply the actual Bayesian ambiguous-CB decisions and reconstruct the
    // exact key that was inserted at resolution time.
    for (auto &kv : ledger.by_id) {
        StarDecision &s = kv.second;
        if (s.state != "PENDING" || !s.have_ambiguous_key) continue;
        const auto it = ledger.resolutions.find(s.ambiguous_key);
        if (it == ledger.resolutions.end() || !it->second.resolved) {
            s.state = "REJECT"; s.reason = "BAYES_UNRESOLVED"; continue;
        }
        const std::string tag_text = detail_value(s.detail, "tag_idx");
        const std::uint32_t tag = tag_text.empty()
            ? (s.tag_token > 0 && s.tag_token < 32
                   ? static_cast<std::uint32_t>(s.tag_token) : 0u)
            : static_cast<std::uint32_t>(std::stoul(tag_text));
        s.molecule_key = (static_cast<std::uint64_t>(it->second.cb_idx & 0xFFFFFu) << 44) |
                         (static_cast<std::uint64_t>(s.umi_packed & 0xFFFFFFu) << 20) |
                         (static_cast<std::uint64_t>(s.gene_idx & 0x7FFFu) << 5) |
                         (tag & 0x1Fu);
        s.have_molecule_key = true; s.state = "KEEP"; s.reason = "AMBIG_CB_RESOLVED";
        s.cb_corrected = it->second.cb;
    }

    BamMap star_bam, cr_bam;
    if (!o.star_bam.empty()) load_bam(o.star_bam, "STAR", o, star_bam);
    for (const auto &bam : o.cr_bams) load_bam(bam.path, bam.label, o, cr_bam);

    std::vector<StarDecision*> ordered;
    ordered.reserve(ledger.by_id.size());
    for (auto &kv : ledger.by_id) {
        if (!kv.second.qname.empty() && kv.second.qname != ".") ordered.push_back(&kv.second);
    }
    std::sort(ordered.begin(), ordered.end(), [](const auto *a, const auto *b) {
        return a->qname < b->qname || (a->qname == b->qname && a->iread < b->iread);
    });

    std::ofstream rows(o.out_dir + "/read_ledger.tsv");
    if (!rows) fail("cannot write read_ledger.tsv");
    rows << "qname\tiread\tlane\tlane_ordinal\ttriage\tcache_class\tcache_gene_idx"
         << "\tnegative_code\thash_offset\tprobe_hamming_distance\tsample_ok\ttag_token\tstar_source\tstar_state"
         << "\tstar_route\tstar_tag_index\tstar_tag_sequence\tstar_composite_cb"
         << "\tstar_reason\tstar_resolver_reason\tstar_gene_idx\tstar_gene_id\tstar_cb_raw"
         << "\tstar_cb_corrected\tstar_cb_candidates\tstar_umi_raw\tstar_umi_corrected"
         << "\tstar_umi_packed\tstar_umi_valid\tmolecule_key\tambiguous_key"
         << "\tstar_bam_present\tstar_bam_records\tstar_flags\tstar_cb_bam\tstar_ub_bam"
         << "\tstar_rep_flag\tstar_mapq\tstar_bam_geometry_scope\tstar_cigar\tstar_md"
         << "\tstar_aligned_query_bases\tstar_aligned_reference_bases"
         << "\tstar_softclip_5p\tstar_softclip_3p\tstar_indel\tstar_splice\tstar_nh\tstar_nm"
         << "\tstar_genomic_flag\tstar_genomic_mapq\tstar_genomic_cigar\tstar_genomic_md\tstar_genomic_nm"
         << "\tcr_bam_present\tcr_sources\tcr_records\tcr_flags\tcr_cb_corrected\tcr_cb_raw\tcr_ub\tcr_fx\tcr_xf\tcr_pr"
         << "\tcr_feature_status\tcr_assigned_included\tcr_counted_included\tcr_gene_ids"
         << "\tcr_rep_flag\tcr_mapq\tcr_cigar\tcr_md"
         << "\tcr_aligned_query_bases\tcr_aligned_reference_bases\tcr_softclip_5p\tcr_softclip_3p"
         << "\tcr_indel\tcr_splice\tcr_nh\tcr_nm"
         << "\tcr_genomic_flag\tcr_genomic_mapq\tcr_genomic_cigar\tcr_genomic_md\tcr_genomic_nm\tcategory\n";

    std::map<std::string,std::uint64_t> categories, gene_summary, mapq_summary,
        cigar_summary, feature_summary, routing_summary, xf_summary, balance_summary;
    std::unordered_map<std::uint64_t,std::uint32_t> raw_key_routes;
    std::unordered_map<std::string,std::uint32_t> qname_count;
    for (const StarDecision *sp : ordered) {
        const StarDecision &s = *sp;
        ++qname_count[s.qname];
        const auto sit = star_bam.find(s.qname), cit = cr_bam.find(s.qname);
        const BamAggregate *sb = sit == star_bam.end() ? nullptr : &sit->second;
        const BamAggregate *cb = cit == cr_bam.end() ? nullptr : &cit->second;
        CrAssignment cr = classify_cr(cb, probes, included_genes);
        const std::string sgene = gene_id(s.gene_idx, genes);
        const std::string cgenes = join(cr.genes);
        const std::string route = routing_class(s);
        std::uint32_t tag_idx = s.have_molecule_key ? static_cast<std::uint32_t>(s.molecule_key & 0x1Fu) : 0;
        const std::string detail_tag = detail_value(s.detail, "tag_idx");
        if (!s.have_molecule_key && !detail_tag.empty()) tag_idx = static_cast<std::uint32_t>(std::stoul(detail_tag));
        if (tag_idx == 0 && s.tag_token > 0 && s.tag_token < 32) tag_idx = static_cast<std::uint32_t>(s.tag_token);
        const std::string tag_seq = tag_idx < tag_sequences.size() ? tag_sequences[tag_idx] : std::string();
        const std::string composite_cb = s.cb_corrected != "." && !s.cb_corrected.empty() && !tag_seq.empty()
            ? s.cb_corrected + tag_seq : std::string();
        std::string category = compare_category(s, cr);
        if (category == "both_assigned") {
            category = std::find(cr.genes.begin(), cr.genes.end(), sgene) != cr.genes.end()
                ? "both_assigned_same_gene" : "both_assigned_gene_conflict";
        }
        ++categories[category];
        ++gene_summary[category + "\t" + sgene + "\t" + cgenes];
        const std::string smapq = sb ? std::to_string(sb->rep.mapq) : ".";
        const std::string cmapq = cb ? std::to_string(cb->rep.mapq) : ".";
        ++mapq_summary[category + "\t" + smapq + "\t" + cmapq];
        const std::string sshape = sb ? sb->rep.geometry.shape : "missing";
        const std::string cshape = cb ? cb->rep.geometry.shape : "missing";
        ++cigar_summary[category + "\t" + sshape + "\t" + cshape];
        ++feature_summary[category + "\t" + cr.status];
        ++routing_summary[route + "\t" + category];
        ++xf_summary[category + "\t" + std::to_string(cr.assigned)
                     + "\t" + std::to_string(cr.counted) + "\t"
                     + (cb ? join(cb->xf) : ".")];
        if (s.state == "KEEP") {
            ++balance_summary[route + "\tSTAR_KEPT_READ_PRE_DEDUP"];
            ++balance_summary["TOTAL\tSTAR_KEPT_READ_PRE_DEDUP"];
            if (s.have_molecule_key) raw_key_routes[s.molecule_key] |= route_bit(route);
            else ++balance_summary[route + "\tSTAR_KEEP_WITHOUT_RAW_KEY"];
        }
        if (cr.assigned) {
            ++balance_summary[route + "\tCR_INCLUDED_FX_ASSIGNED_READ"];
            ++balance_summary["TOTAL\tCR_INCLUDED_FX_ASSIGNED_READ"];
        }
        if (cr.counted) {
            ++balance_summary[route + "\tCR_INCLUDED_XF8_COUNTED_RECORD"];
            ++balance_summary["TOTAL\tCR_INCLUDED_XF8_COUNTED_RECORD"];
        }
        if (category.rfind("star_only",0)==0)
            ++balance_summary[route + "\tSTAR_ONLY_INCLUDED_ASSIGNMENT_READ"];
        if (category.rfind("cr_only",0)==0)
            ++balance_summary[route + "\tCR_ONLY_INCLUDED_ASSIGNMENT_READ"];
        const std::string geometry_scope = sb == nullptr
            ? "NO_STAR_BAM_LEDGER_ONLY"
            : (route == "H1" ? "RUNTIME_BAM_ONLY_NOT_H1_VERIFIER"
                              : (route == "H0" || route.rfind("H1X2_", 0) == 0
                                     ? "RUNTIME_BAM_ONLY_NOT_CACHE_VERIFIER"
                                     : "RUNTIME_RESIDUAL_ALIGNMENT"));

        const CigarStats empty;
        const SamRepresentative sr = sb ? sb->rep : SamRepresentative{};
        const SamRepresentative crp = cb ? cb->rep : SamRepresentative{};
        const SamRepresentative sgr = sb ? sb->genomic_rep : SamRepresentative{};
        const SamRepresentative cgr = cb ? cb->genomic_rep : SamRepresentative{};
        rows << s.qname << '\t' << s.iread << '\t' << s.lane << '\t' << s.lane_ordinal
             << '\t' << value_or_dot(s.triage) << '\t' << s.cache_class << '\t' << s.cache_gene
             << '\t' << s.negative_code << '\t' << s.hash_offset << '\t' << s.probe_hamming_distance
             << '\t' << s.sample_ok
             << '\t' << s.tag_token << '\t' << value_or_dot(s.source) << '\t' << value_or_dot(s.state)
             << '\t' << route << '\t' << tag_idx << '\t' << value_or_dot(tag_seq) << '\t' << value_or_dot(composite_cb)
             << '\t' << value_or_dot(s.reason) << '\t' << value_or_dot(s.resolver_reason)
             << '\t' << s.gene_idx << '\t' << sgene << '\t' << value_or_dot(s.cb_raw)
             << '\t' << value_or_dot(s.cb_corrected) << '\t' << value_or_dot(s.cb_candidates)
             << '\t' << value_or_dot(s.umi_raw) << "\t.\t" << s.umi_packed << '\t' << s.umi_valid
             << '\t' << hex_key(s.molecule_key, s.have_molecule_key)
             << '\t' << hex_key(s.ambiguous_key, s.have_ambiguous_key)
             << '\t' << (sb ? 1 : 0) << '\t' << (sb ? sb->records : 0)
             << '\t' << (sb ? join(sb->flags) : ".") << '\t' << (sb ? join(sb->cb) : ".")
             << '\t' << (sb ? join(sb->ub) : ".") << '\t' << (sb ? sr.flag : -1)
             << '\t' << (sb ? sr.mapq : -1)
             << '\t' << geometry_scope
             << '\t' << (sb ? sr.cigar : ".") << '\t' << (sb ? sr.md : ".")
             << '\t' << sr.geometry.query
             << '\t' << sr.geometry.reference << '\t' << sr.geometry.soft5 << '\t' << sr.geometry.soft3
             << '\t' << sr.geometry.indel << '\t' << sr.geometry.splice << '\t' << sr.nh << '\t' << sr.nm
             << '\t' << (sgr.set ? sgr.flag : -1) << '\t' << (sgr.set ? sgr.mapq : -1)
             << '\t' << (sgr.set ? sgr.cigar : ".") << '\t' << (sgr.set ? sgr.md : ".")
             << '\t' << (sgr.set ? sgr.nm : -1)
             << '\t' << (cb ? 1 : 0) << '\t' << (cb ? join(cb->sources) : ".")
             << '\t' << (cb ? cb->records : 0) << '\t' << (cb ? join(cb->flags) : ".")
             << '\t' << (cb ? (!cb->cb.empty() ? join(cb->cb) : join(cb->gr)) : ".")
             << '\t' << (cb ? join(cb->cr) : ".") << '\t' << (cb ? join(cb->ub) : ".")
             << '\t' << (cb ? join(cb->fx) : ".") << '\t' << (cb ? join(cb->xf) : ".")
             << '\t' << (cb ? join(cb->pr) : ".") << '\t' << cr.status << '\t' << cr.assigned
             << '\t' << cr.counted << '\t' << cgenes << '\t' << (cb ? crp.flag : -1)
             << '\t' << (cb ? crp.mapq : -1) << '\t' << (cb ? crp.cigar : ".")
             << '\t' << (cb ? crp.md : ".")
             << '\t' << crp.geometry.query << '\t' << crp.geometry.reference
             << '\t' << crp.geometry.soft5 << '\t' << crp.geometry.soft3 << '\t' << crp.geometry.indel
             << '\t' << crp.geometry.splice << '\t' << crp.nh << '\t' << crp.nm
             << '\t' << (cgr.set ? cgr.flag : -1) << '\t' << (cgr.set ? cgr.mapq : -1)
             << '\t' << (cgr.set ? cgr.cigar : ".") << '\t' << (cgr.set ? cgr.md : ".")
             << '\t' << (cgr.set ? cgr.nm : -1)
             << '\t' << category << '\n';
    }
    rows.close();

    std::uint64_t duplicate_qnames = 0;
    for (const auto &kv : qname_count) if (kv.second > 1) ++duplicate_qnames;
    write_summary(o.out_dir + "/summary.categories.tsv", "category", categories);
    write_summary(o.out_dir + "/summary.genes.tsv", "category\tstar_gene\tcr_gene", gene_summary);
    write_summary(o.out_dir + "/summary.mapq.tsv", "category\tstar_mapq\tcr_mapq", mapq_summary);
    write_summary(o.out_dir + "/summary.cigar.tsv", "category\tstar_shape\tcr_shape", cigar_summary);
    write_summary(o.out_dir + "/summary.feature_status.tsv", "category\tcr_feature_status", feature_summary);
    write_summary(o.out_dir + "/summary.routing.tsv", "star_route\tcategory", routing_summary);
    write_summary(o.out_dir + "/summary.xf.tsv",
                  "category\tcr_assigned_included\tcr_counted_included\tcr_xf",
                  xf_summary);
    for (const auto &kv : raw_key_routes) {
        ++balance_summary[route_mask_name(kv.second) + "\tSTAR_UNIQUE_RAW_UMI_KEY_PRE_CORRECTION"];
        ++balance_summary["TOTAL\tSTAR_UNIQUE_RAW_UMI_KEY_PRE_CORRECTION"];
    }
    write_summary(o.out_dir + "/summary.balance.tsv", "star_route_or_route_set\tmetric",
                  balance_summary);

    std::ofstream sanity(o.out_dir + "/sanity.tsv");
    sanity << "metric\tvalue\nledger_ids\t" << ledger.by_id.size()
           << "\nledger_qname_rows\t" << ordered.size()
           << "\nduplicate_qnames\t" << duplicate_qnames
           << "\nambiguous_resolution_keys\t" << ledger.resolutions.size()
           << "\nstar_bam_qnames\t" << star_bam.size()
           << "\ncr_bam_qnames\t" << cr_bam.size() << '\n';
    sanity.close();
    if (duplicate_qnames != 0) fail("duplicate QNAMEs in STAR decision ledger");
    std::cout << "ledger_qnames=" << ordered.size() << " star_bam_qnames=" << star_bam.size()
              << " cr_bam_qnames=" << cr_bam.size() << '\n';
    for (const auto &kv : categories) std::cout << kv.first << '\t' << kv.second << '\n';
    return 0;
}

} // namespace

int main(int argc, char **argv) {
    try { return run(parse_options(argc, argv)); }
    catch (const std::exception &e) { std::cerr << "ERROR: " << e.what() << '\n'; return 1; }
}
