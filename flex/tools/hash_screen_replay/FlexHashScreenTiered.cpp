// Two-tier hash screen lookup.
//
// Tier 1 (fast): sample-specific H0 records — exact-match-only search at all
//   3 offsets. On hit with gene>0 → immediate KEEP. ~78% of reads resolved
//   with a single binary search through 53K records (1.3 MB, L2-hot).
//
// Tier 2 (slow): everything else (H0 globals + H1 + deny) as a single flat
//   array. Searched with full findRecord (sampleIdx=0 fallback) + classifyHits.
//   Same cache locality as the pure flat implementation.

#include "FlexHashScreenTiered.h"
#include <algorithm>
#include <unordered_set>

struct SeqSampleKey {
    uint64_t seqHi, seqLo;
    uint16_t sampleIdx;
    bool operator==(const SeqSampleKey& o) const {
        return seqHi == o.seqHi && seqLo == o.seqLo && sampleIdx == o.sampleIdx;
    }
};

struct SeqSampleKeyHash {
    size_t operator()(const SeqSampleKey& k) const {
        // seqLo/seqHi are already well-distributed packed bases; mix in sampleIdx
        return k.seqLo ^ (k.seqHi * 0x9e3779b97f4a7c15ULL) ^ (uint64_t(k.sampleIdx) << 48);
    }
};

void TieredCache::init(const std::vector<Record>& allRecords) {
    h0_.clear();
    rest_.clear();
    h0Count_ = 0;
    h1Count_ = 0;
    h2Count_ = 0;
    denyCount_ = 0;
    dropped_ = 0;
    crossTierDups_ = 0;

    for (const Record& r : allRecords) {
        if (r.cacheClass == 0 && r.sampleIdx != 0 && r.resolvedGeneIdx15 > 0) {
            h0_.push_back(r);
            ++h0Count_;
        } else {
            rest_.push_back(r);
            if (r.resolvedGeneIdx15 == 0 || r.cacheClass == 2)
                ++denyCount_;
            else if (r.cacheClass == 0)
                ++h0Count_;   // H0 globals go in rest_ but count for stats
            else if (r.cacheClass == 1)
                ++h1Count_;
            else if (r.cacheClass == 3)
                ++h2Count_;
            else
                ++dropped_;
        }
    }

    std::sort(h0_.begin(), h0_.end(), recordLess);
    std::sort(rest_.begin(), rest_.end(), recordLess);

    // Cross-tier duplicate detection — seqLo/seqHi are already packed hashes
    using KeySet = std::unordered_set<SeqSampleKey, SeqSampleKeyHash>;
    KeySet h0Keys, h1Keys, h2Keys;
    for (const auto& r : allRecords) {
        if (r.cacheClass == 0)
            h0Keys.insert({r.seqHi, r.seqLo, r.sampleIdx});
        else if (r.cacheClass == 1)
            h1Keys.insert({r.seqHi, r.seqLo, r.sampleIdx});
        else if (r.cacheClass == 3)
            h2Keys.insert({r.seqHi, r.seqLo, r.sampleIdx});
    }
    for (const auto& r : allRecords) {
        if (r.cacheClass == 1) {
            SeqSampleKey k = {r.seqHi, r.seqLo, r.sampleIdx};
            if (h0Keys.count(k)) ++crossTierDups_;
        } else if (r.cacheClass == 3) {
            SeqSampleKey k = {r.seqHi, r.seqLo, r.sampleIdx};
            if (h0Keys.count(k)) ++crossTierDups_;
            if (h1Keys.count(k)) ++crossTierDups_;
        } else if (r.cacheClass == 2 && r.negativeCode == FlexHashNegProbeAmbig) {
            SeqSampleKey k = {r.seqHi, r.seqLo, r.sampleIdx};
            if (h0Keys.count(k)) ++crossTierDups_;
            if (h1Keys.count(k)) ++crossTierDups_;
            if (h2Keys.count(k)) ++crossTierDups_;
        }
    }
}

// Exact-match only: no sampleIdx=0 fallback. Used for the fast H0 tier
// where all records have sampleIdx>0.
bool TieredCache::findExact(const std::vector<Record>& arr,
                            uint64_t seqLo, uint64_t seqHi,
                            uint16_t sampleIdx, Record& out) const {
    Record needle;
    needle.seqLo = seqLo;
    needle.seqHi = seqHi;
    needle.sampleIdx = sampleIdx;

    auto it = std::lower_bound(arr.begin(), arr.end(), needle, recordLess);
    if (it != arr.end() && it->seqLo == seqLo && it->seqHi == seqHi &&
        it->sampleIdx == sampleIdx) {
        out = *it;
        return true;
    }
    return false;
}

// Full findRecord with sampleIdx=0 fallback. Same as FlatCache::findRecord.
bool TieredCache::findRecord(const std::vector<Record>& arr,
                             uint64_t seqLo, uint64_t seqHi,
                             uint16_t sampleIdx, Record& out) const {
    Record needle;
    needle.seqLo = seqLo;
    needle.seqHi = seqHi;
    needle.sampleIdx = sampleIdx;

    auto it = std::lower_bound(arr.begin(), arr.end(), needle, recordLess);
    if (it != arr.end() && it->seqLo == seqLo && it->seqHi == seqHi &&
        it->sampleIdx == sampleIdx) {
        out = *it;
        return true;
    }

    if (sampleIdx == 0)
        return false;

    needle.sampleIdx = 0;
    it = std::lower_bound(arr.begin(), arr.end(), needle, recordLess);
    if (it == arr.end() || it->seqLo != seqLo || it->seqHi != seqHi ||
        it->sampleIdx != 0)
        return false;

    out = *it;
    return true;
}

FlexHashScreenDecision TieredCache::classifyHits(
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

        if (rec->resolvedGeneIdx15 == 0 || rec->cacheClass == 2) {
            sawAmbig = true;
            out.negativeCode = rec->negativeCode != FlexHashNegNone
                ? rec->negativeCode : static_cast<uint8_t>(FlexHashNegProbeAmbig);
            out.offset = relativeOffsets[idx];
            continue;
        }

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

FlexHashScreenDecision TieredCache::classifyRead(
        const char* readSeq, uint32_t readLen, uint16_t sampleIdx) const {
    FlexHashScreenDecision out;
    out.action = FlexHashScreenDecision::Pass;

    if (h0_.empty() && rest_.empty())
        return out;
    if (!readSeq || readLen < kCacheKmerLength)
        return out;

    // Pre-encode all valid offset windows once.
    uint64_t wLo[kNumOffsets] = {}, wHi[kNumOffsets] = {};
    bool wOk[kNumOffsets] = {};

    for (size_t oi = 0; oi < kNumOffsets; ++oi) {
        const int32_t start = kProbeStartOffset + kRelativeProbeOffsets[oi];
        if (start < 0) continue;
        const uint32_t off = static_cast<uint32_t>(start);
        if (off + kCacheKmerLength > readLen) continue;
        wOk[oi] = encodeWindow(readSeq, off, wLo[oi], wHi[oi]);
    }

    // ── Fast tier: H0 exact-match at all 3 offsets ──────────────────────
    for (size_t oi = 0; oi < kNumOffsets; ++oi) {
        if (!wOk[oi]) continue;
        Record h0Hit;
        if (findExact(h0_, wLo[oi], wHi[oi], sampleIdx, h0Hit) &&
            h0Hit.resolvedGeneIdx15 > 0) {
            out.action = FlexHashScreenDecision::Keep;
            out.geneIdx15 = static_cast<uint16_t>(h0Hit.resolvedGeneIdx15);
            out.cacheClass = h0Hit.cacheClass;
            out.offset = kRelativeProbeOffsets[oi];
            return out;
        }
    }

    // ── Slow tier: rest_ (H0 globals + H1 + deny) at all 3 offsets ─────
    Record hits[kNumOffsets];
    const Record* hitPtr[kNumOffsets] = {};

    for (size_t oi = 0; oi < kNumOffsets; ++oi) {
        if (!wOk[oi]) continue;
        if (findRecord(rest_, wLo[oi], wHi[oi], sampleIdx, hits[oi]))
            hitPtr[oi] = &hits[oi];
    }

    return classifyHits(hitPtr, kRelativeProbeOffsets, kNumOffsets, sampleIdx);
}
