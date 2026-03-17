#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

#include "Stats.h"
#include "TrimQc.h"
#include "TrimQcOutput.h"
#include "htslib/kseq.h"

KSEQ_INIT(gzFile, gzread)

using std::cerr;
using std::cout;
using std::runtime_error;
using std::string;
using std::vector;

namespace {

struct Config {
    vector<string> inputPaths;
    string reportPrefix;
    string stageLabel = "fastq_replay";
    uint64 maxReads = 0;
    uint32 mateCount = 1;
    uint32 mateIndex = 1;
    uint32 qualityBase = 33;
};

void usage(const char* argv0) {
    cerr
        << "Usage: " << argv0 << " --input <fastq[.gz]> [--input <fastq2> ...] --report <prefix> [options]\n"
        << "\n"
        << "Options:\n"
        << "  --stage <label>       Stage label stored in JSON/HTML (default: fastq_replay)\n"
        << "  --max-reads <N>       Limit reads collected for QC (default: 0 = all)\n"
        << "  --mate-count <1|2>    Total mate count reported in JSON (default: 1)\n"
        << "  --mate-index <1|2>    Mate slot to populate (default: 1)\n"
        << "  --quality-base <N>    Phred quality ASCII base (default: 33)\n"
        << "  -h, --help            Show this help\n";
}

void appendCommaSeparated(const string& value, vector<string>& out) {
    size_t start = 0;
    while (start <= value.size()) {
        const size_t comma = value.find(',', start);
        const string token = value.substr(start, comma == string::npos ? string::npos : comma - start);
        if (!token.empty()) {
            out.push_back(token);
        }
        if (comma == string::npos) {
            break;
        }
        start = comma + 1;
    }
}

char toTrimQcBase(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
        case 'U':
            return 3;
        default:
            return 4;
    }
}

void processFastq(const string& inputPath, const Config& cfg, TrimQcCollector& qc,
                  uint64& readsSeen) {
    gzFile file = gzopen(inputPath.c_str(), "rb");
    if (file == NULL) {
        throw runtime_error("could not open FASTQ input: " + inputPath);
    }

    kseq_t* seq = kseq_init(file);
    if (seq == NULL) {
        gzclose(file);
        throw runtime_error("could not initialize FASTQ reader for: " + inputPath);
    }

    while (true) {
        if (cfg.maxReads > 0 && qc.totalReads() >= cfg.maxReads) {
            break;
        }

        const int status = kseq_read(seq);
        if (status < 0) {
            if (status != -1) {
                kseq_destroy(seq);
                gzclose(file);
                throw runtime_error("malformed FASTQ record in: " + inputPath);
            }
            break;
        }

        ++readsSeen;
        if (seq->qual.l == 0) {
            kseq_destroy(seq);
            gzclose(file);
            throw runtime_error("missing quality string in FASTQ input: " + inputPath);
        }
        if (seq->seq.l != seq->qual.l) {
            kseq_destroy(seq);
            gzclose(file);
            throw runtime_error("sequence/quality length mismatch in FASTQ input: " + inputPath);
        }

        vector<char> seqNum(static_cast<size_t>(seq->seq.l));
        for (size_t i = 0; i < static_cast<size_t>(seq->seq.l); ++i) {
            seqNum[i] = toTrimQcBase(seq->seq.s[i]);
        }

        qc.addRead(seqNum.data(), seq->qual.s, static_cast<uint32>(seq->seq.l), cfg.mateIndex - 1, 0);
    }

    kseq_destroy(seq);
    gzclose(file);
}

} // namespace

int main(int argc, char* argv[]) {
    Config cfg;

    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];
        if (arg == "--input" && i + 1 < argc) {
            appendCommaSeparated(argv[++i], cfg.inputPaths);
        } else if (arg == "--report" && i + 1 < argc) {
            cfg.reportPrefix = argv[++i];
        } else if (arg == "--stage" && i + 1 < argc) {
            cfg.stageLabel = argv[++i];
        } else if (arg == "--max-reads" && i + 1 < argc) {
            cfg.maxReads = strtoull(argv[++i], NULL, 10);
        } else if (arg == "--mate-count" && i + 1 < argc) {
            cfg.mateCount = static_cast<uint32>(strtoul(argv[++i], NULL, 10));
        } else if (arg == "--mate-index" && i + 1 < argc) {
            cfg.mateIndex = static_cast<uint32>(strtoul(argv[++i], NULL, 10));
        } else if (arg == "--quality-base" && i + 1 < argc) {
            cfg.qualityBase = static_cast<uint32>(strtoul(argv[++i], NULL, 10));
        } else if (arg == "-h" || arg == "--help") {
            usage(argv[0]);
            return 0;
        } else {
            cerr << "ERROR: unknown argument: " << arg << "\n";
            usage(argv[0]);
            return 2;
        }
    }

    if (cfg.inputPaths.empty() || cfg.reportPrefix.empty()) {
        usage(argv[0]);
        return 2;
    }
    if (cfg.mateCount < 1 || cfg.mateCount > 2) {
        cerr << "ERROR: --mate-count must be 1 or 2\n";
        return 2;
    }
    if (cfg.mateIndex < 1 || cfg.mateIndex > cfg.mateCount) {
        cerr << "ERROR: --mate-index must be between 1 and --mate-count\n";
        return 2;
    }

    try {
        TrimQcCollector qc;
        qc.init(cfg.mateCount, cfg.maxReads, cfg.qualityBase);
        Stats stats;
        uint64 readsSeen = 0;

        for (const string& inputPath : cfg.inputPaths) {
            if (cfg.maxReads > 0 && qc.totalReads() >= cfg.maxReads) {
                break;
            }
            processFastq(inputPath, cfg, qc, readsSeen);
        }

        const string jsonPath = cfg.reportPrefix + ".trim_qc.json";
        const string htmlPath = cfg.reportPrefix + ".trim_qc.html";
        if (!writeTrimQcJson(qc, stats, jsonPath, cfg.stageLabel)) {
            throw runtime_error("failed to write trim QC JSON: " + jsonPath);
        }
        if (!writeTrimQcHtml(jsonPath, htmlPath)) {
            throw runtime_error("failed to write trim QC HTML: " + htmlPath);
        }

        cout
            << "inputs\t" << cfg.inputPaths.size() << "\n"
            << "report_prefix\t" << cfg.reportPrefix << "\n"
            << "stage\t" << cfg.stageLabel << "\n"
            << "mate_count\t" << cfg.mateCount << "\n"
            << "mate_index\t" << cfg.mateIndex << "\n"
            << "reads_seen\t" << readsSeen << "\n"
            << "reads_added\t" << qc.totalReads() << "\n"
            << "json\t" << jsonPath << "\n"
            << "html\t" << htmlPath << "\n";
        return 0;
    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
