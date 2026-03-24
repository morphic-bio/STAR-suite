// Reference implementation: flat binary search through the full sorted record
// array. This is a standalone copy of FlexHashScreen.cpp::classifyRead and
// must produce identical output to STAR's in-process decisions.

#include "FlexHashScreenFlat.h"
#include <algorithm>

void FlatCache::init(const std::vector<Record>& allRecords) {
    records_ = allRecords;
    std::sort(records_.begin(), records_.end(), recordLess);
}

bool FlatCache::findRecord(uint64_t seqLo, uint64_t seqHi,
                           uint16_t sampleIdx, Record& out) const {
    Record needle;
    needle.seqLo = seqLo;
    needle.seqHi = seqHi;
    needle.sampleIdx = sampleIdx;

    auto it = std::lower_bound(records_.begin(), records_.end(), needle, recordLess);
    if (it != records_.end() && it->seqLo == seqLo && it->seqHi == seqHi &&
        it->sampleIdx == sampleIdx) {
        out = *it;
        return true;
    }

    if (sampleIdx == 0)
        return false;

    needle.sampleIdx = 0;
    it = std::lower_bound(records_.begin(), records_.end(), needle, recordLess);
    if (it == records_.end() || it->seqLo != seqLo || it->seqHi != seqHi ||
        it->sampleIdx != 0)
        return false;

    out = *it;
    return true;
}

FlexHashScreenDecision FlatCache::classifyHits(
        const Record* const* hits, const int8_t* relativeOffsets,
        size_t nHits, uint16_t runtimeSampleIdx) const {
    FlexHashScreenDecision out;
    out.action = FlexHashScreenDecision::Pass;

    bool sawAmbig = false;
    bool sawSampleMismatch = false;
    int8_t sampleMismatchOffset = 0;
    bool sawNonExactKeep = false;
    uint16_t nonExactGene = 0;
    uint16_t nonExactSample = 0;
    uint8_t nonExactClass = 0;
    int8_t nonExactOffset = 0;
    bool sawGeneConflict = false;
    int8_t geneConflictOffset = 0;

    for (size_t idx = 0; idx < nHits; ++idx) {
        const Record* rec = hits[idx];
        if (!rec) continue;

        const bool sampleMatched = (rec->sampleIdx != 0 && rec->sampleIdx == runtimeSampleIdx);
        const bool sampleSpecifiedMismatch = (rec->sampleIdx != 0 && rec->sampleIdx != runtimeSampleIdx);

        if (rec->cacheClass == 2 && rec->negativeCode == FlexHashNegProbeAmbig) {
            sawAmbig = true;
            out.negativeCode = rec->negativeCode;
            out.offset = relativeOffsets[idx];
            continue;
        }

        if (rec->resolvedGeneIdx15 == 0)
            continue;

        if (rec->cacheClass == 0 && sampleMatched) {
            out.action = FlexHashScreenDecision::Keep;
            out.geneIdx15 = static_cast<uint16_t>(rec->resolvedGeneIdx15);
            out.cacheClass = rec->cacheClass;
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
                nonExactOffset = relativeOffsets[idx];
            } else if (nonExactGene != geneIdx15 || nonExactSample != sampleKey) {
                sawGeneConflict = true;
                geneConflictOffset = relativeOffsets[idx];
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

FlexHashScreenDecision FlatCache::classifyRead(
        const char* readSeq, uint32_t readLen, uint16_t sampleIdx) const {
    if (records_.empty()) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Pass;
        return out;
    }

    if (!readSeq || readLen < kCacheKmerLength) {
        FlexHashScreenDecision out;
        out.action = FlexHashScreenDecision::Pass;
        return out;
    }

    Record hits[kNumOffsets];
    const Record* hitPtr[kNumOffsets] = {};

    for (size_t idx = 0; idx < kNumOffsets; ++idx) {
        const int32_t start = kProbeStartOffset + kRelativeProbeOffsets[idx];
        if (start < 0) continue;
        const uint32_t offset = static_cast<uint32_t>(start);
        if (offset + kCacheKmerLength > readLen) continue;

        uint64_t seqLo = 0, seqHi = 0;
        if (!encodeWindow(readSeq, offset, seqLo, seqHi)) continue;
        if (findRecord(seqLo, seqHi, sampleIdx, hits[idx]))
            hitPtr[idx] = &hits[idx];
    }

    return classifyHits(hitPtr, kRelativeProbeOffsets, kNumOffsets, sampleIdx);
}
