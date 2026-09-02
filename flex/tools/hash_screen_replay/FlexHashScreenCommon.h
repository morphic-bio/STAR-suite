#ifndef H_FlexHashScreenCommon
#define H_FlexHashScreenCommon

// Standalone hash screen types and utilities, extracted from
// core/legacy/source/FlexHashScreen.{h,cpp} with no STAR dependencies.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

// ── Decision ────────────────────────────────────────────────────────────────

struct FlexHashScreenDecision {
    enum Action : uint8_t {
        Disabled = 0,
        Pass     = 1,
        Keep     = 2,
        Deny     = 3
    };

    Action   action     = Disabled;
    uint16_t geneIdx15  = 0;
    uint8_t  cacheClass = 0;
    uint8_t  negativeCode = 0;
    int8_t   offset     = 0;

    bool operator==(const FlexHashScreenDecision& o) const {
        return action == o.action && geneIdx15 == o.geneIdx15 &&
               cacheClass == o.cacheClass && negativeCode == o.negativeCode &&
               offset == o.offset;
    }
    bool operator!=(const FlexHashScreenDecision& o) const { return !(*this == o); }
};

enum FlexHashScreenNegativeCode : uint8_t {
    FlexHashNegNone      = 0,
    FlexHashNegProbeAmbig = 1
};

// ── Record ──────────────────────────────────────────────────────────────────

struct Record {
    uint64_t seqLo            = 0;
    uint64_t seqHi            = 0;
    uint32_t resolvedGeneIdx15 = 0;
    uint8_t  cacheClass       = 0;
    uint8_t  negativeCode     = 0;
    uint16_t sampleIdx        = 0;
};

inline bool recordLess(const Record& lhs, const Record& rhs) {
    if (lhs.seqHi != rhs.seqHi) return lhs.seqHi < rhs.seqHi;
    if (lhs.seqLo != rhs.seqLo) return lhs.seqLo < rhs.seqLo;
    return lhs.sampleIdx < rhs.sampleIdx;
}

// ── Constants ───────────────────────────────────────────────────────────────

static const uint16_t kCacheKmerLength = 50;
static const int8_t   kRelativeProbeOffsets[] = {0, 1, -1};
static const int32_t  kProbeStartOffset = 0;
static const size_t   kNumOffsets = sizeof(kRelativeProbeOffsets) / sizeof(kRelativeProbeOffsets[0]);

// ── encodeWindow ────────────────────────────────────────────────────────────

inline bool encodeWindow(const char* readSeq, uint32_t offset,
                         uint64_t& seqLo, uint64_t& seqHi) {
    uint64_t packedLo = 0, packedHi = 0;
    for (uint32_t i = 0; i < kCacheKmerLength; ++i) {
        char base = readSeq[offset + i];
        uint64_t code;
        switch (base) {
            case 'A': case 'a': code = 0; break;
            case 'C': case 'c': code = 1; break;
            case 'G': case 'g': code = 2; break;
            case 'T': case 't': code = 3; break;
            default: return false;
        }
        packedHi = (packedHi << 2) | (packedLo >> 62);
        packedLo = (packedLo << 2) | code;
    }
    seqLo = packedLo;
    seqHi = packedHi;
    return true;
}

// ── Cache file loader (FH01SEQ1 binary format) ─────────────────────────────

struct CacheHeaderRaw {
    char     magic[8];
    uint16_t version;
    uint16_t kmerLength;
    uint32_t recordSize;
    uint64_t recordCount;
};

struct CacheRecordRaw {
    uint64_t seqLo;
    uint64_t seqHi;
    uint32_t resolvedGeneIdx15;
    uint8_t  cacheClass;
    uint8_t  negativeCode;
    uint16_t reserved;
};

static const char     kCacheMagic[] = {'F','H','0','1','S','E','Q','1'};
static const uint16_t kCacheVersionSampleAware = 2;
static const uint16_t kCacheVersionProbeRegion = 3;
static const uint16_t kCacheRecordSize = 24;

inline bool loadCacheRecords(const char* path, std::vector<Record>& out,
                             std::string* errorOut = nullptr) {
    FILE* f = std::fopen(path, "rb");
    if (!f) {
        if (errorOut) *errorOut = "cannot open cache";
        return false;
    }

    CacheHeaderRaw header {};
    if (std::fread(&header, sizeof(header), 1, f) != 1) {
        if (errorOut) *errorOut = "cannot read cache header";
        std::fclose(f);
        return false;
    }

    if (!std::equal(header.magic, header.magic + 8, kCacheMagic)) {
        if (errorOut) *errorOut = "cache magic mismatch";
        std::fclose(f);
        return false;
    }
    if ((header.version != 1 && header.version != kCacheVersionSampleAware &&
         header.version != kCacheVersionProbeRegion) ||
        header.kmerLength != kCacheKmerLength || header.recordSize != kCacheRecordSize) {
        if (errorOut) *errorOut = "cache format mismatch";
        std::fclose(f);
        return false;
    }

    out.resize(static_cast<size_t>(header.recordCount));
    for (size_t i = 0; i < out.size(); ++i) {
        CacheRecordRaw raw {};
        if (std::fread(&raw, sizeof(raw), 1, f) != 1) {
            out.clear();
            if (errorOut) *errorOut = "cache truncated";
            std::fclose(f);
            return false;
        }
        Record& rec = out[i];
        rec.seqLo            = raw.seqLo;
        rec.seqHi            = raw.seqHi;
        // v3 stores the two-bit probe region in bits 30-31. The standalone
        // replay tools classify genes only, so discard the metadata while
        // preserving v1/v2 decision parity.
        rec.resolvedGeneIdx15 = raw.resolvedGeneIdx15 & 0x7FFFu;
        rec.cacheClass       = raw.cacheClass;
        rec.negativeCode     = raw.negativeCode;
        rec.sampleIdx        = (header.version >= kCacheVersionSampleAware) ? raw.reserved : 0;
    }
    std::fclose(f);

    std::sort(out.begin(), out.end(), recordLess);
    return true;
}

// ── Dump file reader (HSCRN001 binary format) ──────────────────────────────

struct DumpRecord {
    std::string readSeq;
    uint32_t    readLen   = 0;
    uint16_t    sampleIdx = 0;
    FlexHashScreenDecision truth;
};

inline bool loadDumpRecords(const char* path, std::vector<DumpRecord>& out,
                            std::string* errorOut = nullptr) {
    FILE* f = std::fopen(path, "rb");
    if (!f) {
        if (errorOut) *errorOut = "cannot open dump";
        return false;
    }

    char magic[8];
    if (std::fread(magic, 1, 8, f) != 8 ||
        std::memcmp(magic, "HSCRN001", 8) != 0) {
        if (errorOut) *errorOut = "dump magic mismatch";
        std::fclose(f);
        return false;
    }

    uint64_t nReads = 0;
    if (std::fread(&nReads, 8, 1, f) != 1) {
        if (errorOut) *errorOut = "cannot read dump header";
        std::fclose(f);
        return false;
    }

    out.reserve(static_cast<size_t>(nReads));
    for (uint64_t i = 0; i < nReads; ++i) {
        DumpRecord dr;
        uint8_t action, cacheClass;
        int8_t offset;

        if (std::fread(&dr.readLen, 4, 1, f) != 1) break;
        if (std::fread(&dr.sampleIdx, 2, 1, f) != 1) break;
        if (std::fread(&action, 1, 1, f) != 1) break;
        if (std::fread(&cacheClass, 1, 1, f) != 1) break;
        if (std::fread(&dr.truth.geneIdx15, 2, 1, f) != 1) break;
        if (std::fread(&offset, 1, 1, f) != 1) break;
        uint8_t negativeCode;
        if (std::fread(&negativeCode, 1, 1, f) != 1) break;

        dr.truth.action       = static_cast<FlexHashScreenDecision::Action>(action);
        dr.truth.cacheClass   = cacheClass;
        dr.truth.offset       = offset;
        dr.truth.negativeCode = negativeCode;

        dr.readSeq.resize(dr.readLen);
        if (dr.readLen > 0 &&
            std::fread(&dr.readSeq[0], 1, dr.readLen, f) != dr.readLen) break;

        out.push_back(std::move(dr));
    }
    std::fclose(f);

    if (out.size() != static_cast<size_t>(nReads)) {
        if (errorOut) *errorOut = "dump truncated";
        return false;
    }
    return true;
}

#endif // H_FlexHashScreenCommon
