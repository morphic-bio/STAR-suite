#include "IncludeDefine.h"
#include SAMTOOLS_SAM_H

#include "Stats.h"
#include "TrimQc.h"
#include "TrimQcOutput.h"

namespace {

struct Config {
    string inputPath;
    string reportPrefix;
    string stageLabel = "alignment_replay";
    uint64 maxReads = 0;
    int mateCount = 0; // 0 = auto
    bool primaryOnly = true;
};

void usage(const char* argv0) {
    cerr
        << "Usage: " << argv0 << " --input <alignment.sam|bam|cram> --report <prefix> [options]\n"
        << "\n"
        << "Options:\n"
        << "  --stage <label>       Stage label stored in JSON/HTML (default: alignment_replay)\n"
        << "  --max-reads <N>       Limit reads collected for QC (default: 0 = all)\n"
        << "  --mate-count <1|2|auto>\n"
        << "                        Override mate count; default is auto-detect from first primary record\n"
        << "  --all-records         Include secondary and supplementary alignments\n"
        << "  -h, --help            Show this help\n";
}

bool shouldSkipRecord(const bam1_t* record, bool primaryOnly) {
    if (!primaryOnly) {
        return false;
    }
    return (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) != 0;
}

uint32_t inferMateCount(const Config& cfg) {
    samFile* in = sam_open(cfg.inputPath.c_str(), "r");
    if (in == NULL) {
        throw runtime_error("could not open alignment input: " + cfg.inputPath);
    }

    bam_hdr_t* hdr = sam_hdr_read(in);
    bam1_t* record = bam_init1();
    if (hdr == NULL || record == NULL) {
        if (record != NULL) {
            bam_destroy1(record);
        }
        if (hdr != NULL) {
            bam_hdr_destroy(hdr);
        }
        sam_close(in);
        throw runtime_error("could not initialize SAM/BAM reader for: " + cfg.inputPath);
    }

    uint32_t mateCount = 1;
    while (sam_read1(in, hdr, record) >= 0) {
        if (shouldSkipRecord(record, cfg.primaryOnly)) {
            continue;
        }
        if ((record->core.flag & BAM_FPAIRED) != 0) {
            mateCount = 2;
        }
        break;
    }

    bam_destroy1(record);
    bam_hdr_destroy(hdr);
    sam_close(in);
    return mateCount;
}

uint8_t toTrimQcBase(uint8_t bamBase) {
    static const char* bamNt16 = seq_nt16_str;
    switch (bamNt16[bamBase & 0x0f]) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return 4;
    }
}

char complementBaseCode(char baseCode) {
    switch (baseCode) {
        case 0:
            return 3; // A -> T
        case 1:
            return 2; // C -> G
        case 2:
            return 1; // G -> C
        case 3:
            return 0; // T -> A
        default:
            return 4; // N
    }
}

} // namespace

int main(int argc, char* argv[]) {
    Config cfg;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--input" && i + 1 < argc) {
            cfg.inputPath = argv[++i];
        } else if (arg == "--report" && i + 1 < argc) {
            cfg.reportPrefix = argv[++i];
        } else if (arg == "--stage" && i + 1 < argc) {
            cfg.stageLabel = argv[++i];
        } else if (arg == "--max-reads" && i + 1 < argc) {
            cfg.maxReads = strtoull(argv[++i], NULL, 10);
        } else if (arg == "--mate-count" && i + 1 < argc) {
            string value = argv[++i];
            if (value == "auto") {
                cfg.mateCount = 0;
            } else if (value == "1" || value == "2") {
                cfg.mateCount = atoi(value.c_str());
            } else {
                cerr << "ERROR: --mate-count must be 1, 2, or auto\n";
                return 2;
            }
        } else if (arg == "--all-records") {
            cfg.primaryOnly = false;
        } else if (arg == "-h" || arg == "--help") {
            usage(argv[0]);
            return 0;
        } else {
            cerr << "ERROR: unknown argument: " << arg << "\n";
            usage(argv[0]);
            return 2;
        }
    }

    if (cfg.inputPath.empty() || cfg.reportPrefix.empty()) {
        usage(argv[0]);
        return 2;
    }

    try {
        uint32_t mateCount = cfg.mateCount == 0 ? inferMateCount(cfg) : static_cast<uint32_t>(cfg.mateCount);

        TrimQcCollector qc;
        qc.init(mateCount, cfg.maxReads, 33);
        Stats stats;

        samFile* in = sam_open(cfg.inputPath.c_str(), "r");
        if (in == NULL) {
            throw runtime_error("could not open alignment input: " + cfg.inputPath);
        }

        bam_hdr_t* hdr = sam_hdr_read(in);
        bam1_t* record = bam_init1();
        if (hdr == NULL || record == NULL) {
            if (record != NULL) {
                bam_destroy1(record);
            }
            if (hdr != NULL) {
                bam_hdr_destroy(hdr);
            }
            sam_close(in);
            throw runtime_error("could not initialize SAM/BAM reader for: " + cfg.inputPath);
        }

        uint64_t recordsSeen = 0;
        uint64_t recordsAdded = 0;
        uint64_t recordsSkippedSecondary = 0;

        while (sam_read1(in, hdr, record) >= 0) {
            ++recordsSeen;
            if (shouldSkipRecord(record, cfg.primaryOnly)) {
                ++recordsSkippedSecondary;
                continue;
            }

            const int32_t readLength = record->core.l_qseq;
            if (readLength <= 0) {
                continue;
            }

            vector<char> seqNum(static_cast<size_t>(readLength));
            string qualAscii(static_cast<size_t>(readLength), '!');

            const uint8_t* bamSeq = bam_get_seq(record);
            const uint8_t* bamQual = bam_get_qual(record);
            for (int32_t i = 0; i < readLength; ++i) {
                seqNum[static_cast<size_t>(i)] = static_cast<char>(toTrimQcBase(bam_seqi(bamSeq, i)));
                if (bamQual[i] != 255) {
                    qualAscii[static_cast<size_t>(i)] = static_cast<char>(bamQual[i] + 33);
                }
            }
            if ((record->core.flag & BAM_FREVERSE) != 0) {
                reverse(seqNum.begin(), seqNum.end());
                reverse(qualAscii.begin(), qualAscii.end());
                for (char& baseCode : seqNum) {
                    baseCode = complementBaseCode(baseCode);
                }
            }

            const uint32_t mateIndex = ((record->core.flag & BAM_FREAD2) != 0 && mateCount > 1) ? 1u : 0u;
            qc.addRead(seqNum.data(), qualAscii.data(), static_cast<uint32_t>(readLength), mateIndex, 0);
            if (qc.totalReads() > recordsAdded) {
                ++recordsAdded;
            }
        }

        bam_destroy1(record);
        bam_hdr_destroy(hdr);
        sam_close(in);

        const string jsonPath = cfg.reportPrefix + ".trim_qc.json";
        const string htmlPath = cfg.reportPrefix + ".trim_qc.html";
        if (!writeTrimQcJson(qc, stats, jsonPath, cfg.stageLabel)) {
            throw runtime_error("failed to write trim QC JSON: " + jsonPath);
        }
        if (!writeTrimQcHtml(jsonPath, htmlPath)) {
            throw runtime_error("failed to write trim QC HTML: " + htmlPath);
        }

        cout
            << "input\t" << cfg.inputPath << "\n"
            << "report_prefix\t" << cfg.reportPrefix << "\n"
            << "stage\t" << cfg.stageLabel << "\n"
            << "mate_count\t" << mateCount << "\n"
            << "records_seen\t" << recordsSeen << "\n"
            << "records_added\t" << recordsAdded << "\n"
            << "records_skipped_secondary\t" << recordsSkippedSecondary << "\n"
            << "json\t" << jsonPath << "\n"
            << "html\t" << htmlPath << "\n";
        return 0;
    } catch (const exception& e) {
        cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
