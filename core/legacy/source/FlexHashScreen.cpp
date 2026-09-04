#include "FlexHashScreen.h"
#include "ParametersSolo.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>

namespace {

const char kCacheMagic[] = {'F', 'H', '0', '1', 'S', 'E', 'Q', '1'};
const uint16_t kCacheVersionSampleAware = 2;
const uint16_t kCacheVersionProbeRegion = 3;
const uint16_t kCacheKmerLength = 50;
const uint16_t kCacheRecordSize = 24;
const int8_t kRelativeProbeOffsets[] = {0, 1, -1};
const int32_t kProbeStartOffset = 0;

struct CacheHeaderRaw {
    char magic[8];
    uint16_t version;
    uint16_t kmerLength;
    uint32_t recordSize;
    uint64_t recordCount;
};

struct CacheRecordRaw {
    uint64_t seqLo;
    uint64_t seqHi;
    uint32_t resolvedGeneIdx15;
    uint8_t cacheClass;
    uint8_t negativeCode;
    uint16_t reserved;
};

inline bool recordLess(const FlexHashScreenCache::Record& lhs, const FlexHashScreenCache::Record& rhs) {
    if (lhs.seqHi != rhs.seqHi) {
        return lhs.seqHi < rhs.seqHi;
    }
    if (lhs.seqLo != rhs.seqLo) {
        return lhs.seqLo < rhs.seqLo;
    }
    return lhs.sampleIdx < rhs.sampleIdx;
}

} // namespace

FlexHashScreenCache& FlexHashScreenCache::instance() {
    static FlexHashScreenCache cache;
    return cache;
}

bool FlexHashScreenCache::ensureLoaded(const ParametersSolo& pSolo, std::string* errorOut) {
    if (initialized_) {
        if (!enabled_ && errorOut != nullptr && !loadedPath_.empty()) {
            *errorOut = loadedPath_;
        }
        return enabled_;
    }

    initialized_ = true;
    enabled_ = false;
    loadedPath_.clear();
    records_.clear();
    cacheVersion_ = 0;
    regionMetadataComplete_ = false;

    if (!pSolo.hashScreenEnabled || pSolo.hashScreenFile.empty()) {
        return false;
    }

    enabled_ = loadFile(pSolo.hashScreenFile, errorOut);
    if (enabled_) {
        loadedPath_ = pSolo.hashScreenFile;
    }
    return enabled_;
}

bool FlexHashScreenCache::loadFile(const std::string& path, std::string* errorOut) {
    std::ifstream in(path.c_str(), std::ios::binary);
    if (!in.good()) {
        if (errorOut != nullptr) {
            *errorOut = "cannot open cache";
        }
        return false;
    }

    CacheHeaderRaw header {};
    in.read(reinterpret_cast<char*>(&header), sizeof(header));
    if (!in.good()) {
        if (errorOut != nullptr) {
            *errorOut = "cannot read cache header";
        }
        return false;
    }

    if (!std::equal(header.magic, header.magic + 8, kCacheMagic)) {
        if (errorOut != nullptr) {
            *errorOut = "cache magic mismatch";
        }
        return false;
    }
    if ((header.version != 1 && header.version != kCacheVersionSampleAware
         && header.version != kCacheVersionProbeRegion) ||
        header.kmerLength != kCacheKmerLength || header.recordSize != kCacheRecordSize) {
        if (errorOut != nullptr) {
            *errorOut = "cache format mismatch";
        }
        return false;
    }

    records_.resize(static_cast<size_t>(header.recordCount));
    for (size_t i = 0; i < records_.size(); ++i) {
        CacheRecordRaw raw {};
        in.read(reinterpret_cast<char*>(&raw), sizeof(raw));
        if (!in.good()) {
            records_.clear();
            if (errorOut != nullptr) {
                *errorOut = "cache truncated";
            }
            return false;
        }
        Record rec;
        rec.seqLo = raw.seqLo;
        rec.seqHi = raw.seqHi;
        rec.resolvedGeneIdx15 = raw.resolvedGeneIdx15 & 0x7FFFu;
        rec.probeRegion = header.version >= kCacheVersionProbeRegion
            ? static_cast<FlexGdnaRegion>((raw.resolvedGeneIdx15 >> 30) & 0x3u)
            : FlexGdnaUnknown;
        rec.cacheClass = raw.cacheClass;
        rec.negativeCode = raw.negativeCode;
        rec.sampleIdx = (header.version >= kCacheVersionSampleAware) ? raw.reserved : 0;
        records_[i] = rec;
    }

    cacheVersion_ = header.version;
    regionMetadataComplete_ = header.version >= kCacheVersionProbeRegion;
    std::sort(records_.begin(), records_.end(), recordLess);
    buildTieredVectors();
    return true;
}

bool FlexHashScreenCache::encodeWindow(const char* readSeq, uint32_t offset, uint64_t& seqLo, uint64_t& seqHi) const {
    uint64_t packedLo = 0;
    uint64_t packedHi = 0;
    for (uint32_t i = 0; i < kCacheKmerLength; ++i) {
        char base = readSeq[offset + i];
        uint64_t code = 0;
        switch (base) {
            case 'A':
            case 'a':
                code = 0;
                break;
            case 'C':
            case 'c':
                code = 1;
                break;
            case 'G':
            case 'g':
                code = 2;
                break;
            case 'T':
            case 't':
                code = 3;
                break;
            default:
                return false;
        }
        packedHi = (packedHi << 2) | (packedLo >> 62);
        packedLo = (packedLo << 2) | code;
    }
    seqLo = packedLo;
    seqHi = packedHi;
    return true;
}

bool FlexHashScreenCache::findRecord(uint64_t seqLo, uint64_t seqHi, uint16_t sampleIdx, Record& out) const {
    Record needle;
    needle.seqLo = seqLo;
    needle.seqHi = seqHi;
    needle.sampleIdx = sampleIdx;
    auto it = std::lower_bound(records_.begin(), records_.end(), needle, recordLess);
    if (it != records_.end() && it->seqLo == seqLo && it->seqHi == seqHi && it->sampleIdx == sampleIdx) {
        out = *it;
        return true;
    }

    if (sampleIdx == 0) {
        return false;
    }

    needle.sampleIdx = 0;
    it = std::lower_bound(records_.begin(), records_.end(), needle, recordLess);
    if (it == records_.end() || it->seqLo != seqLo || it->seqHi != seqHi || it->sampleIdx != 0) {
        return false;
    }
    out = *it;
    return true;
}

FlexHashScreenDecision FlexHashScreenCache::classifyHits(const Record* const* hits, const int8_t* relativeOffsets, size_t nHits, uint16_t runtimeSampleIdx) const {
    FlexHashScreenDecision out;
    out.action = FlexHashScreenDecision::Pass;
    bool sawAmbig = false;
    bool sawSampleMismatch = false;
    int8_t sampleMismatchOffset = 0;
    bool sawNonExactKeep = false;
    uint16_t nonExactGene = 0;
    uint16_t nonExactSample = 0;
    uint8_t nonExactClass = 0;
    FlexGdnaRegion nonExactRegion = FlexGdnaUnknown;
    int8_t nonExactOffset = 0;
    bool sawGeneConflict = false;
    int8_t geneConflictOffset = 0;

    for (size_t idx = 0; idx < nHits; ++idx) {
        const Record* rec = hits[idx];
        if (rec == nullptr) {
            continue;
        }

        const bool sampleMatched = (rec->sampleIdx != 0 && rec->sampleIdx == runtimeSampleIdx);
        const bool sampleSpecifiedMismatch = (rec->sampleIdx != 0 && rec->sampleIdx != runtimeSampleIdx);

        if (rec->resolvedGeneIdx15 == 0 || rec->cacheClass == 2) {
            sawAmbig = true;
            out.negativeCode = rec->negativeCode != FlexHashNegNone
                ? rec->negativeCode : FlexHashNegProbeAmbig;
            out.offset = relativeOffsets[idx];
            continue;
        }

        if (rec->cacheClass == 0 && sampleMatched) {
            out.action = FlexHashScreenDecision::Keep;
            out.geneIdx15 = static_cast<uint16_t>(rec->resolvedGeneIdx15);
            out.cacheClass = rec->cacheClass;
            out.probeRegion = rec->probeRegion;
            out.offset = relativeOffsets[idx];
            return out;
        }

        if ((rec->cacheClass == 0 || rec->cacheClass == 1 || rec->cacheClass == 3) &&
            sampleSpecifiedMismatch) {
            if (!sawSampleMismatch) {
                sawSampleMismatch = true;
                sampleMismatchOffset = relativeOffsets[idx];
            }
            continue;
        }

        if (rec->cacheClass == 0 || rec->cacheClass == 1 || rec->cacheClass == 3) {
            const uint16_t geneIdx15 = static_cast<uint16_t>(rec->resolvedGeneIdx15);
            const uint16_t sampleKey = sampleMatched ? runtimeSampleIdx : 0;
            if (!sawNonExactKeep) {
                sawNonExactKeep = true;
                nonExactGene = geneIdx15;
                nonExactSample = sampleKey;
                nonExactClass = rec->cacheClass;
                nonExactRegion = rec->probeRegion;
                nonExactOffset = relativeOffsets[idx];
            } else if (nonExactGene != geneIdx15 || nonExactSample != sampleKey) {
                sawGeneConflict = true;
                geneConflictOffset = relativeOffsets[idx];
            } else {
                nonExactRegion = flexGdnaMergeRegion(nonExactRegion, rec->probeRegion);
            }
        }
    }

    if (sawAmbig) {
        out.action = FlexHashScreenDecision::Deny;
        return out;
    }

    if (sawGeneConflict) {
        out.action = FlexHashScreenDecision::Deny;
        out.geneIdx15 = 0;
        out.cacheClass = 0;
        out.negativeCode = FlexHashNegProbeAmbig;
        out.offset = geneConflictOffset;
        return out;
    }

    if (sawNonExactKeep) {
        out.action = FlexHashScreenDecision::Keep;
        out.geneIdx15 = nonExactGene;
        out.cacheClass = nonExactClass;
        out.probeRegion = nonExactRegion;
        out.offset = nonExactOffset;
        return out;
    }

    if (sawSampleMismatch) {
        out.action = FlexHashScreenDecision::Deny;
        out.geneIdx15 = 0;
        out.cacheClass = 0;
        out.negativeCode = FlexHashNegProbeAmbig;
        out.offset = sampleMismatchOffset;
        return out;
    }

    return out;
}

void FlexHashScreenCache::buildTieredVectors() {
    h0Records_.clear();
    h1DenyRecords_.clear();
    for (const Record& rec : records_) {
        if (rec.cacheClass == 0 && rec.resolvedGeneIdx15 > 0) {
            h0Records_.push_back(rec);
        } else if ((rec.cacheClass == 1 && rec.resolvedGeneIdx15 > 0) ||
                   rec.resolvedGeneIdx15 == 0) {
            h1DenyRecords_.push_back(rec);
        }
    }
    std::sort(h0Records_.begin(), h0Records_.end(), recordLess);
    std::sort(h1DenyRecords_.begin(), h1DenyRecords_.end(), recordLess);
    buildH0NoSampleMap();
    buildH1DenyNoSampleMap();
}

// ── LUT for branchless base-to-2bit conversion ──────────────────────────────

uint8_t FlexHashScreenCache::baseLUT_[256] = {};
bool FlexHashScreenCache::lutInitialized_ = false;

static void initBaseLUT(uint8_t table[256]) {
    std::memset(table, 0xFF, 256);
    table['A'] = table['a'] = 0;
    table['C'] = table['c'] = 1;
    table['G'] = table['g'] = 2;
    table['T'] = table['t'] = 3;
}

bool FlexHashScreenCache::encodeWindowLUT(const char* readSeq, uint32_t offset,
                                           uint64_t& seqLo, uint64_t& seqHi) const {
    if (!lutInitialized_) {
        initBaseLUT(baseLUT_);
        lutInitialized_ = true;
    }
    uint64_t lo = 0, hi = 0;
    for (uint32_t i = 0; i < kCacheKmerLength; ++i) {
        uint8_t code = baseLUT_[static_cast<uint8_t>(readSeq[offset + i])];
        if (code == 0xFF) return false;
        hi = (hi << 2) | (lo >> 62);
        lo = (lo << 2) | code;
    }
    seqLo = lo;
    seqHi = hi;
    return true;
}

void FlexHashScreenCache::buildH0NoSampleMap() {
    h0NoSampleMap_.clear();
    h0NoSampleMap_.reserve(h0Records_.size() * 2);
    for (const Record& r : h0Records_) {
        SeqKeyNoSample key{r.seqLo, r.seqHi};
        h0NoSampleMap_.emplace(key, r);
    }
}

void FlexHashScreenCache::buildH1DenyNoSampleMap() {
    h1DenyNoSampleMap_.clear();
    h1DenyNoSampleMap_.reserve(h1DenyRecords_.size() * 2);
    for (const Record& r : h1DenyRecords_) {
        SeqKeyNoSample key{r.seqLo, r.seqHi};
        h1DenyNoSampleMap_.emplace(key, r);
    }
}

FlexHashScreenDecision FlexHashScreenCache::classifyReadH0Offset0(const char* readSeq, uint32_t readLen) const {
    if (!initialized_ || !enabled_) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Disabled;
        return out;
    }

    if (readSeq == nullptr || readLen < kCacheKmerLength) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Pass;
        return out;
    }

    const uint32_t offset = static_cast<uint32_t>(kProbeStartOffset);
    if (offset + kCacheKmerLength > readLen) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Pass;
        return out;
    }

    uint64_t seqLo = 0, seqHi = 0;
    if (!encodeWindowLUT(readSeq, offset, seqLo, seqHi)) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Pass;
        return out;
    }

    SeqKeyNoSample key{seqLo, seqHi};
    auto it = h0NoSampleMap_.find(key);
    if (it == h0NoSampleMap_.end()) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Pass;
        return out;
    }

    const Record& rec = it->second;
    FlexHashScreenDecision out;
    if (rec.negativeCode == FlexHashNegProbeAmbig) {
        out.action = FlexHashScreenDecision::Deny;
        out.geneIdx15 = 0;
        out.cacheClass = rec.cacheClass;
        out.probeRegion = rec.probeRegion;
        out.negativeCode = rec.negativeCode;
        out.offset = 0;
    } else {
        out.action = FlexHashScreenDecision::Keep;
        out.geneIdx15 = static_cast<uint16_t>(rec.resolvedGeneIdx15);
        out.cacheClass = rec.cacheClass;
        out.probeRegion = rec.probeRegion;
        out.negativeCode = 0;
        out.offset = 0;
    }
    return out;
}

FlexHashScreenDecision FlexHashScreenCache::classifyReadH0H1Offset0(const char* readSeq, uint32_t readLen) const {
    if (!initialized_ || !enabled_) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Disabled;
        return out;
    }

    if (readSeq == nullptr || readLen < kCacheKmerLength) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Pass;
        return out;
    }

    const uint32_t offset = static_cast<uint32_t>(kProbeStartOffset);
    if (offset + kCacheKmerLength > readLen) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Pass;
        return out;
    }

    uint64_t seqLo = 0, seqHi = 0;
    if (!encodeWindowLUT(readSeq, offset, seqLo, seqHi)) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Pass;
        return out;
    }

    SeqKeyNoSample key{seqLo, seqHi};

    // H0 check
    auto it = h0NoSampleMap_.find(key);
    if (it != h0NoSampleMap_.end()) {
        const Record& rec = it->second;
        FlexHashScreenDecision out;
        if (rec.negativeCode == FlexHashNegProbeAmbig) {
            out.action = FlexHashScreenDecision::Deny;
            out.geneIdx15 = 0;
            out.cacheClass = rec.cacheClass;
            out.probeRegion = rec.probeRegion;
            out.negativeCode = rec.negativeCode;
        } else {
            out.action = FlexHashScreenDecision::Keep;
            out.geneIdx15 = static_cast<uint16_t>(rec.resolvedGeneIdx15);
            out.probeRegion = rec.probeRegion;
            out.cacheClass = rec.cacheClass;
        }
        out.offset = 0;
        return out;
    }

    // H1+Deny check (same key, different tier)
    it = h1DenyNoSampleMap_.find(key);
    if (it != h1DenyNoSampleMap_.end()) {
        const Record& rec = it->second;
        FlexHashScreenDecision out;
        if (rec.cacheClass == 2 && rec.negativeCode == FlexHashNegProbeAmbig) {
            out.action = FlexHashScreenDecision::Deny;
            out.geneIdx15 = 0;
            out.negativeCode = rec.negativeCode;
        } else if (rec.resolvedGeneIdx15 > 0) {
            out.action = FlexHashScreenDecision::Keep;
            out.geneIdx15 = static_cast<uint16_t>(rec.resolvedGeneIdx15);
            out.probeRegion = rec.probeRegion;
        } else {
            // A cache hit is resolved evidence. If it is not a valid KEEP,
            // it is a certified negative rather than an unknown sequence;
            // only a cache miss is eligible for residual genome alignment.
            out.action = FlexHashScreenDecision::Deny;
            out.negativeCode = rec.negativeCode != FlexHashNegNone
                ? rec.negativeCode : FlexHashNegProbeAmbig;
        }
        out.cacheClass = rec.cacheClass;
        out.offset = 0;
        return out;
    }

    FlexHashScreenDecision out;
    out.action = FlexHashScreenDecision::Pass;
    return out;
}

bool FlexHashScreenCache::findRecordInVec(const std::vector<Record>& vec, uint64_t seqLo, uint64_t seqHi, uint16_t sampleIdx, Record& out) const {
    Record needle;
    needle.seqLo = seqLo;
    needle.seqHi = seqHi;
    needle.sampleIdx = sampleIdx;
    auto it = std::lower_bound(vec.begin(), vec.end(), needle, recordLess);
    if (it != vec.end() && it->seqLo == seqLo && it->seqHi == seqHi && it->sampleIdx == sampleIdx) {
        out = *it;
        return true;
    }
    if (sampleIdx == 0) {
        return false;
    }
    needle.sampleIdx = 0;
    it = std::lower_bound(vec.begin(), vec.end(), needle, recordLess);
    if (it == vec.end() || it->seqLo != seqLo || it->seqHi != seqHi || it->sampleIdx != 0) {
        return false;
    }
    out = *it;
    return true;
}

FlexHashScreenDecision FlexHashScreenCache::classifyReadH0Only(const char* readSeq, uint32_t readLen, uint16_t sampleIdx) const {
    if (!initialized_ || !enabled_) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Disabled;
        return out;
    }

    if (readSeq == nullptr || readLen < kCacheKmerLength) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Pass;
        return out;
    }

    Record hits[sizeof(kRelativeProbeOffsets) / sizeof(kRelativeProbeOffsets[0])];
    const Record* hitPtr[sizeof(kRelativeProbeOffsets) / sizeof(kRelativeProbeOffsets[0])] = {nullptr, nullptr, nullptr};
    for (size_t idx = 0; idx < sizeof(kRelativeProbeOffsets) / sizeof(kRelativeProbeOffsets[0]); ++idx) {
        const int32_t start = kProbeStartOffset + kRelativeProbeOffsets[idx];
        if (start < 0) {
            continue;
        }
        const uint32_t offset = static_cast<uint32_t>(start);
        if (offset + kCacheKmerLength > readLen) {
            continue;
        }
        uint64_t seqLo = 0;
        uint64_t seqHi = 0;
        if (!encodeWindow(readSeq, offset, seqLo, seqHi)) {
            continue;
        }
        if (findRecordInVec(h0Records_, seqLo, seqHi, sampleIdx, hits[idx])) {
            hitPtr[idx] = &hits[idx];
        }
    }

    return classifyHits(hitPtr, kRelativeProbeOffsets, sizeof(hitPtr) / sizeof(hitPtr[0]), sampleIdx);
}

bool FlexHashScreenCache::encodeProbeWindow(const char* readSeq, uint32_t offset, uint64_t& seqLo, uint64_t& seqHi) {
    uint64_t packedLo = 0;
    uint64_t packedHi = 0;
    for (uint32_t i = 0; i < kCacheKmerLength; ++i) {
        char base = readSeq[offset + i];
        uint64_t code = 0;
        switch (base) {
            case 'A':
            case 'a':
                code = 0;
                break;
            case 'C':
            case 'c':
                code = 1;
                break;
            case 'G':
            case 'g':
                code = 2;
                break;
            case 'T':
            case 't':
                code = 3;
                break;
            default:
                return false;
        }
        packedHi = (packedHi << 2) | (packedLo >> 62);
        packedLo = (packedLo << 2) | code;
    }
    seqLo = packedLo;
    seqHi = packedHi;
    return true;
}

bool FlexHashScreenCache::writeHashCacheFile(const std::string& path, std::vector<Record>& records,
                                             std::string* errorOut, bool regionMetadataComplete) {
    std::sort(records.begin(), records.end(), recordLess);
    std::ofstream out(path.c_str(), std::ios::binary);
    if (!out.good()) {
        if (errorOut != nullptr) {
            *errorOut = "cannot open hash cache output for write";
        }
        return false;
    }
    CacheHeaderRaw header {};
    std::memcpy(header.magic, kCacheMagic, 8);
    header.version = regionMetadataComplete ? kCacheVersionProbeRegion : kCacheVersionSampleAware;
    header.kmerLength = kCacheKmerLength;
    header.recordSize = kCacheRecordSize;
    header.recordCount = static_cast<uint64_t>(records.size());
    out.write(reinterpret_cast<const char*>(&header), sizeof(header));
    if (!out.good()) {
        if (errorOut != nullptr) {
            *errorOut = "cannot write hash cache header";
        }
        return false;
    }
    for (const Record& rec : records) {
        CacheRecordRaw raw {};
        raw.seqLo = rec.seqLo;
        raw.seqHi = rec.seqHi;
        raw.resolvedGeneIdx15 = rec.resolvedGeneIdx15 & 0x7FFFu;
        if (regionMetadataComplete) {
            raw.resolvedGeneIdx15 |=
                (static_cast<uint32_t>(rec.probeRegion) & 0x3u) << 30;
        }
        raw.cacheClass = rec.cacheClass;
        raw.negativeCode = rec.negativeCode;
        raw.reserved = rec.sampleIdx;
        out.write(reinterpret_cast<const char*>(&raw), sizeof(raw));
        if (!out.good()) {
            if (errorOut != nullptr) {
                *errorOut = "cannot write hash cache record";
            }
            return false;
        }
    }
    return true;
}

FlexHashScreenDecision FlexHashScreenCache::classifyRead(const char* readSeq, uint32_t readLen, uint16_t sampleIdx) const {
    if (!initialized_ || !enabled_) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Disabled;
        return out;
    }

    if (readSeq == nullptr || readLen < kCacheKmerLength) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Pass;
        return out;
    }

    Record hits[sizeof(kRelativeProbeOffsets) / sizeof(kRelativeProbeOffsets[0])];
    const Record* hitPtr[sizeof(kRelativeProbeOffsets) / sizeof(kRelativeProbeOffsets[0])] = {nullptr, nullptr, nullptr};
    for (size_t idx = 0; idx < sizeof(kRelativeProbeOffsets) / sizeof(kRelativeProbeOffsets[0]); ++idx) {
        const int32_t start = kProbeStartOffset + kRelativeProbeOffsets[idx];
        if (start < 0) {
            continue;
        }
        const uint32_t offset = static_cast<uint32_t>(start);
        if (offset + kCacheKmerLength > readLen) {
            continue;
        }
        uint64_t seqLo = 0;
        uint64_t seqHi = 0;
        if (!encodeWindow(readSeq, offset, seqLo, seqHi)) {
            continue;
        }
        if (findRecord(seqLo, seqHi, sampleIdx, hits[idx])) {
            hitPtr[idx] = &hits[idx];
        }
    }

    return classifyHits(hitPtr, kRelativeProbeOffsets, sizeof(hitPtr) / sizeof(hitPtr[0]), sampleIdx);
}
