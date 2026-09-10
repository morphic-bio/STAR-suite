#include "FlexHashScreen.h"
#include "Parameters.h"
#include "ParametersSolo.h"

#include <algorithm>
#include <cstring>

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
    cacheVersion_ = 0;
    regionMetadataComplete_ = false;
    hasH1X2_ = false;

    if (!pSolo.hashScreenEnabled || pSolo.hashScreenFile.empty()) {
        return false;
    }

    enabled_ = loadFile(pSolo.hashScreenFile, errorOut);
    if (enabled_ && hasH1X2_ && !pair_.ready() && !buildH1X2ProbeIndex(errorOut)) {
        enabled_ = false;
    }
    if (enabled_) {
        loadedPath_ = pSolo.hashScreenFile;
    }
    return enabled_;
}

bool FlexHashScreenCache::loadFile(const std::string& path, std::string* errorOut) {
    if (FlexProbePairKhash::isSnapshot(path)) {
        if (!pair_.open(path,errorOut)) return false;
        cacheVersion_ = 3; regionMetadataComplete_ = true; hasH1X2_ = true;
        return true;
    }
    if (!storage_.open(path, errorOut)) return false;
    cacheVersion_ = storage_.sourceVersion();
    regionMetadataComplete_ = cacheVersion_ >= kCacheVersionProbeRegion;
    hasH1X2_ = storage_.hasH1X2();
    return true;
}

FlexHashScreenCache::Record FlexHashScreenCache::decodeStorageRecord(const FlexProbeRecord& raw) const {
    Record rec;
    rec.seqLo = raw.key.lo; rec.seqHi = raw.key.hi;
    rec.resolvedGeneIdx15 = raw.value.geneAndRegion & 0x7fffu;
    rec.probeRegion = regionMetadataComplete_
        ? static_cast<FlexGdnaRegion>((raw.value.geneAndRegion >> 30) & 3u) : FlexGdnaUnknown;
    rec.cacheClass = raw.value.cacheClass; rec.negativeCode = raw.value.negativeCode;
    rec.sampleIdx = cacheVersion_ >= kCacheVersionSampleAware ? raw.value.sample : 0;
    return rec;
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

bool FlexHashScreenCache::findRecord(uint64_t seqLo, uint64_t seqHi, uint16_t sampleIdx,
                                      Record& out, bool h0Only) const {
    FlexProbeRecord raw;
    if (pair_.ready()) {
        if (!pair_.find(SeqKeyNoSample{seqLo,seqHi}, sampleIdx, h0Only, raw)) return false;
    } else if (!storage_.find(SeqKeyNoSample{seqLo, seqHi}, sampleIdx, h0Only, raw)) return false;
    out = decodeStorageRecord(raw);
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

        if ((rec->cacheClass == FlexHashCacheH0 || rec->cacheClass == FlexHashCacheH1 ||
             rec->cacheClass == FlexHashCacheH2 || rec->cacheClass == FlexHashCacheH1X2) &&
            sampleSpecifiedMismatch) {
            if (!sawSampleMismatch) {
                sawSampleMismatch = true;
                sampleMismatchOffset = relativeOffsets[idx];
            }
            continue;
        }

        if (rec->cacheClass == FlexHashCacheH0 || rec->cacheClass == FlexHashCacheH1 ||
            rec->cacheClass == FlexHashCacheH2 || rec->cacheClass == FlexHashCacheH1X2) {
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

std::string FlexHashScreenCache::decodeCacheSequence(uint64_t seqLo,
                                                      uint64_t seqHi)
{
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::string sequence(kCacheKmerLength, 'N');
    for (unsigned position = 0; position < kCacheKmerLength; ++position) {
        const uint64_t packed = position < 18 ? seqHi : seqLo;
        const unsigned shift = position < 18
            ? 2U * (17U - position)
            : 2U * (49U - position);
        sequence[position] = bases[(packed >> shift) & UINT64_C(3)];
    }
    return sequence;
}

bool FlexHashScreenCache::buildH1X2ProbeIndex(std::string* errorOut)
{
    if (!hasH1X2_ || probeSeedIndex_.ready())
        return true;

    // H0 entries are the authoritative active probe parents already embedded
    // in the cache. Reconstructing from them keeps count-only/no-genome and
    // alignment runs on the identical probe set and removes a genome-loading
    // dependency from this probe-only classifier. H0 is sample-specific, so
    // adjacent copies of the same packed probe are collapsed here.
    std::vector<FlexProbeHalfIndex::Probe> probes;
    probes.reserve(storage_.h0Count());
    uint64_t previousLo = 0;
    uint64_t previousHi = 0;
    uint16_t previousGene = 0;
    bool havePrevious = false;
    for (uint64_t i = 0; i < storage_.h0Count(); ++i) {
        const Record record = decodeStorageRecord(storage_.h0Record(i));
        if (havePrevious && record.seqLo == previousLo &&
            record.seqHi == previousHi) {
            if (record.resolvedGeneIdx15 != previousGene) {
                if (errorOut != nullptr)
                    *errorOut = "H1X2 cache has one H0 probe assigned to multiple genes";
                return false;
            }
            continue;
        }
        probes.push_back(FlexProbeHalfIndex::Probe(
            decodeCacheSequence(record.seqLo, record.seqHi),
            static_cast<uint16_t>(record.resolvedGeneIdx15)));
        previousLo = record.seqLo;
        previousHi = record.seqHi;
        previousGene = static_cast<uint16_t>(record.resolvedGeneIdx15);
        havePrevious = true;
    }

    if (probes.empty()) {
        if (errorOut != nullptr)
            *errorOut = "H1X2 seed-and-extend requires H0 parent records in the same cache";
        return false;
    }
    if (!probeSeedIndex_.build(probes, errorOut))
        return false;
    return true;
}

namespace {

FlexHashScreenDecision probeSeedDecision(
    const FlexProbeHalfIndex::Result& anchor)
{
    FlexHashScreenDecision decision;
    if (anchor.status == FlexProbeHalfIndex::Unique) {
        // The fixed-position 25-base seed plus 50-base Hamming extension is
        // the complete H1X2 rescue. Do not send this read through genomic
        // alignment: Cell Ranger's Flex gene assignment is probe-based, and
        // its genomic alignment is output annotation only.
        decision.action = FlexHashScreenDecision::Keep;
        decision.geneIdx15 = anchor.geneIdx15;
        decision.cacheClass = FlexHashCacheH1X2;
        decision.probeHammingDistance = anchor.hammingDistance;
    } else {
        decision.action = FlexHashScreenDecision::Deny;
        if (anchor.status == FlexProbeHalfIndex::ScoreFail)
            decision.probeHammingDistance = anchor.hammingDistance;
        if (anchor.status == FlexProbeHalfIndex::Ambiguous) {
            decision.negativeCode = FlexHashNegHalfGeneAmbig;
        } else if (anchor.status == FlexProbeHalfIndex::ScoreFail) {
            decision.negativeCode = FlexHashNegHalfScoreFail;
        } else if (anchor.status == FlexProbeHalfIndex::SplitProbe) {
            decision.negativeCode = FlexHashNegHalfSplitProbe;
        } else {
            decision.negativeCode = FlexHashNegHalfNoAnchor;
        }
    }
    return decision;
}

} // namespace

FlexHashScreenDecision FlexHashScreenCache::classifyReadH1X2SeedExtend(
    const char* readSeq, uint32_t readLen) const
{
    if (pair_.ready()) {
        FlexHashScreenDecision pass; pass.action=FlexHashScreenDecision::Pass;
        if (!readSeq || readLen<50) { pass.action=FlexHashScreenDecision::Deny; pass.negativeCode=FlexHashNegHalfNoAnchor; return pass; }
        uint64_t lo=0,hi=0,nmask=0;
        for (unsigned i=0;i<50;++i) {
            const char c=readSeq[i]; unsigned b=0;
            switch(c) { case 'A':case 'a':break;case 'C':case 'c':b=1;break;case 'G':case 'g':b=2;break;case 'T':case 't':b=3;break;default:nmask|=UINT64_C(1)<<i; }
            if(i<32)lo|=uint64_t(b)<<(2*i);else hi|=uint64_t(b)<<(2*(i-32));
        }
        return pair_.seedOnly({lo,hi},nmask).decision;
    }
    if (!hasH1X2_ || !probeSeedIndex_.ready()) {
        FlexHashScreenDecision decision;
        decision.action = FlexHashScreenDecision::Pass;
        return decision;
    }
    return probeSeedDecision(probeSeedIndex_.classify(readSeq, readLen));
}

FlexHashScreenDecision FlexHashScreenCache::classifyCbqH1X2SeedExtend(
    uint64_t seqLo, uint64_t seqHi, uint64_t nMask) const
{
    if (pair_.ready()) return pair_.seedOnly({seqLo,seqHi},nMask).decision;
    if (!hasH1X2_ || !probeSeedIndex_.ready()) {
        FlexHashScreenDecision decision;
        decision.action = FlexHashScreenDecision::Pass;
        return decision;
    }
    const SeqKeyNoSample cacheKey = cbqKeyToCacheKey(seqLo, seqHi);
    return probeSeedDecision(probeSeedIndex_.classifyCachePacked(
        cacheKey.lo, cacheKey.hi, nMask));
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

static inline uint64_t reverseTwoBitGroups(uint64_t value) {
    value = ((value & UINT64_C(0x3333333333333333)) << 2) |
            ((value >> 2) & UINT64_C(0x3333333333333333));
    value = ((value & UINT64_C(0x0F0F0F0F0F0F0F0F)) << 4) |
            ((value >> 4) & UINT64_C(0x0F0F0F0F0F0F0F0F));
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_bswap64(value);
#else
    value = ((value & UINT64_C(0x00FF00FF00FF00FF)) << 8) |
            ((value >> 8) & UINT64_C(0x00FF00FF00FF00FF));
    value = ((value & UINT64_C(0x0000FFFF0000FFFF)) << 16) |
            ((value >> 16) & UINT64_C(0x0000FFFF0000FFFF));
    return (value << 32) | (value >> 32);
#endif
}

FlexHashScreenCache::SeqKeyNoSample FlexHashScreenCache::cacheKeyToCbqKey(
    uint64_t seqLo, uint64_t seqHi) {
    return FlexHashCacheStorage::cbqKey(SeqKeyNoSample{seqLo, seqHi});
}

FlexHashScreenCache::SeqKeyNoSample FlexHashScreenCache::cbqKeyToCacheKey(
    uint64_t seqLo, uint64_t seqHi) {
    const uint64_t last32 = (seqLo >> 36) |
        ((seqHi & UINT64_C(0xFFFFFFFFF)) << 28);
    const uint64_t first18 = seqLo & UINT64_C(0xFFFFFFFFF);
    SeqKeyNoSample out;
    out.lo = reverseTwoBitGroups(last32);
    out.hi = reverseTwoBitGroups(first18 << 28);
    return out;
}

FlexHashScreenCache::SeqKeyNoSample FlexHashScreenCache::offset0MapKeyFromCacheKey(
    uint64_t seqLo, uint64_t seqHi) const {
    return cacheKeyToCbqKey(seqLo, seqHi);
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

    return classifyH0Offset0MapKey(offset0MapKeyFromCacheKey(seqLo, seqHi));
}

FlexHashScreenDecision FlexHashScreenCache::classifyCbqH0Offset0(
    uint64_t seqLo, uint64_t seqHi) const {
    const SeqKeyNoSample key{seqLo, seqHi};
    return classifyH0Offset0MapKey(key);
}

FlexHashScreenDecision FlexHashScreenCache::classifyH0Offset0MapKey(
    const SeqKeyNoSample& key) const {
    if (!initialized_ || !enabled_) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Disabled;
        return out;
    }

    const auto* value = pair_.ready() ? pair_.lookupH0(key) : storage_.lookup(0, key);
    if (!value) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Pass;
        return out;
    }

    const Record rec = decodeStorageRecord(FlexProbeRecord{key, *value});
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

    return classifyH0H1Offset0MapKey(offset0MapKeyFromCacheKey(seqLo, seqHi));
}

namespace {
// Merge the four single-N substitutions into one verdict.
struct SingleNMerge {
    FlexHashScreenDecision keep;
    bool haveKeep = false, ambiguous = false, sawDeny = false;
    uint8_t matchedClass = 0xFF;
    void add(const FlexHashScreenDecision& d) {
        if (d.action == FlexHashScreenDecision::Deny) { sawDeny = true; return; }
        if (d.action != FlexHashScreenDecision::Keep) return;
        if (!haveKeep) { keep = d; haveKeep = true; matchedClass = d.cacheClass; }
        else if (d.geneIdx15 != keep.geneIdx15) ambiguous = true;
        else if (d.cacheClass != matchedClass) matchedClass = 0xFE;
    }
    FlexHashScreenDecision result() const {
        FlexHashScreenDecision out;
        // The true base is unknown, so a disagreement between the substitutions (or a
        // DENY record among them) is not resolved evidence: report a miss, which keeps
        // the read eligible for residual alignment, rather than a certified negative.
        if (haveKeep && !ambiguous && !sawDeny) {
            out = keep;
            out.singleN = true;
            out.singleNCacheClass = matchedClass;
            out.cacheClass = FlexHashCacheH1;
            return out;
        }
        out.singleN = true;
        out.action = FlexHashScreenDecision::Pass; return out;
    }
};
} // namespace

FlexHashScreenDecision FlexHashScreenCache::classifyCbqH0H1Offset0SingleN(
    uint64_t seqLo, uint64_t seqHi, uint64_t nMask) const {
    if (nMask == 0 || (nMask & (nMask - 1)) != 0) {
        FlexHashScreenDecision out; out.action = FlexHashScreenDecision::Pass; return out;
    }
    const unsigned pos = static_cast<unsigned>(__builtin_ctzll(nMask));
    SingleNMerge m;
    for (uint64_t b = 0; b < 4; ++b) {
        uint64_t lo = seqLo, hi = seqHi;
        if (pos < 32) lo = (lo & ~(UINT64_C(3) << (2 * pos))) | (b << (2 * pos));
        else          hi = (hi & ~(UINT64_C(3) << (2 * (pos - 32)))) | (b << (2 * (pos - 32)));
        m.add(classifyCbqH0H1Offset0(lo, hi));
    }
    return m.result();
}

FlexHashScreenDecision FlexHashScreenCache::classifyReadH0H1Offset0SingleN(const char* readSeq, uint32_t readLen) const {
    FlexHashScreenDecision pass; pass.action = FlexHashScreenDecision::Pass;
    if (readSeq == nullptr || readLen < kCacheKmerLength) return pass;
    int nPos = -1;
    for (uint32_t i = 0; i < kCacheKmerLength; ++i) {
        const char c = readSeq[i];
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'a' || c == 'c' || c == 'g' || c == 't') continue;
        if (nPos >= 0) return pass;   // two or more Ns
        nPos = static_cast<int>(i);
    }
    if (nPos < 0) return pass;
    char buf[kCacheKmerLength];
    std::memcpy(buf, readSeq, kCacheKmerLength);
    SingleNMerge m;
    for (const char* p = "ACGT"; *p; ++p) { buf[nPos] = *p; m.add(classifyReadH0H1Offset0(buf, kCacheKmerLength)); }
    return m.result();
}

FlexHashScreenDecision FlexHashScreenCache::classifyCbqH0H1Offset0(
    uint64_t seqLo, uint64_t seqHi) const {
    const SeqKeyNoSample key{seqLo, seqHi};
    return classifyH0H1Offset0MapKey(key);
}

FlexHashScreenDecision FlexHashScreenCache::classifyH0H1Offset0MapKey(
    const SeqKeyNoSample& key) const {
    if (!initialized_ || !enabled_) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Disabled;
        return out;
    }

    if (pair_.ready()) return pair_.primaryOnly(key).decision;

    // H0 check
    const auto* value = storage_.lookup(0, key);
    if (value) {
        const Record rec = decodeStorageRecord(FlexProbeRecord{key, *value});
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
    value = storage_.lookup(1, key);
    if (value) {
        const Record rec = decodeStorageRecord(FlexProbeRecord{key, *value});
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

uint32_t FlexHashScreenCache::probeWindowLength() {
    return kCacheKmerLength;
}

const char* flexHashScreenDenyReason(uint8_t negativeCode)
{
    switch (negativeCode) {
        case FlexHashNegHalfNoAnchor:
            return "H1X2_HALF_NO_ANCHOR";
        case FlexHashNegHalfGeneAmbig:
            return "H1X2_HALF_PROBE_AMBIG";
        case FlexHashNegHalfScoreFail:
            return "H1X2_PROBE_SCORE_FAIL";
        case FlexHashNegHalfSplitProbe:
            return "H1X2_SPLIT_PROBE";
        default:
            return "NEG_PROBE_AMBIG";
    }
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
        if (findRecord(seqLo, seqHi, sampleIdx, hits[idx], true)) {
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

FlexHashScreenDecision FlexHashScreenCache::classifyReadComplete(const char* read,uint32_t len,bool singleN) const {
    if (pair_.ready() && read && len>=50) return pair_.classifyRead(read,singleN).decision;
    auto d=classifyReadH0H1Offset0(read,len);
    if(d.action==FlexHashScreenDecision::Pass && singleN) d=classifyReadH0H1Offset0SingleN(read,len);
    if(d.action==FlexHashScreenDecision::Pass && h1x2ProbeIndexReady())d=classifyReadH1X2SeedExtend(read,len);
    return d;
}
FlexHashScreenDecision FlexHashScreenCache::classifyCbqComplete(uint64_t lo,uint64_t hi,uint64_t nmask,bool singleN) const {
    if(pair_.ready())return pair_.classify({lo,hi},nmask,singleN).decision;
    FlexHashScreenDecision d;d.action=FlexHashScreenDecision::Pass;
    if(!nmask)d=classifyCbqH0H1Offset0(lo,hi);
    else if(singleN)d=classifyCbqH0H1Offset0SingleN(lo,hi,nmask);
    if(d.action==FlexHashScreenDecision::Pass && h1x2ProbeIndexReady())d=classifyCbqH1X2SeedExtend(lo,hi,nmask);
    return d;
}
