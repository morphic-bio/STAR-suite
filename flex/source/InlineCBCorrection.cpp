#include "InlineCBCorrection.h"
#include "SequenceFuns.h"
#include <algorithm>
#include <cstring>
#include <thread>
#include <pthread.h>
#include <atomic>

// Static members
khash_t(cbwl) *InlineCBCorrection::exactHash_ = nullptr;
khash_t(cbwl) *InlineCBCorrection::variantHash_ = nullptr;
std::unordered_set<uint64_t> InlineCBCorrection::whitelistHash_;
std::unordered_map<uint64_t, std::vector<uint32_t>> InlineCBCorrection::variantCollisions_;
const std::vector<std::string>* InlineCBCorrection::wlStrings_ = nullptr;
bool InlineCBCorrection::whitelistInitialized_ = false;
std::vector<std::vector<uint64_t>*> InlineCBCorrection::evidenceShards_;
std::mutex InlineCBCorrection::evidenceMutex_;
size_t InlineCBCorrection::evidenceSize_ = 0;
static std::atomic<uint64_t> gEvidenceHits{0};
static std::atomic<uint64_t> gAmbigResolved{0};
std::vector<InlineCBCorrection::AmbigShard*> InlineCBCorrection::ambigShards_;
std::mutex InlineCBCorrection::ambigMutex_;
uint64_t InlineCBCorrection::ambigTotal_ = 0;

size_t InlineCBCorrection::exactMapSize() { return exactHash_ ? kh_size(exactHash_) : 0; }
size_t InlineCBCorrection::variantMapSize() { return variantHash_ ? kh_size(variantHash_) : 0; }
size_t InlineCBCorrection::variantCollisionSize() { return variantCollisions_.size(); }
size_t InlineCBCorrection::variantCollisionMaxFanout() {
    size_t maxFan = 0;
    for (const auto &kv : variantCollisions_) {
        if (kv.second.size() > maxFan) {
            maxFan = kv.second.size();
        }
    }
    return maxFan;
}
uint64_t InlineCBCorrection::parentEvidenceTotal() {
    uint64_t total = 0;
    for (auto *shard : evidenceShards_) {
        if (!shard) continue;
        for (uint64_t v : *shard) {
            total += v;
        }
    }
    // include hits for sanity visibility
    if (total == 0) {
        total = gEvidenceHits.load();
    }
    return total;
}
uint64_t InlineCBCorrection::ambigCapturedTotal() { return ambigTotal_; }
size_t InlineCBCorrection::ambigUniqueVariants() {
    size_t total = 0;
    for (auto *shard : ambigShards_) {
        if (!shard) continue;
        total += shard->records.size();
    }
    return total;
}
size_t InlineCBCorrection::ambigMaxParents() {
    size_t maxP = 0;
    for (auto *shard : ambigShards_) {
        if (!shard) continue;
        for (auto &rec : shard->records) {
            if (rec.parents.size() > maxP) maxP = rec.parents.size();
        }
    }
    return maxP;
}
uint64_t InlineCBCorrection::getResolvedAmbigCount() { return gAmbigResolved.load(); }

uint64_t InlineCBCorrection::packCBForLookup(const std::string &cbSeq) {
    return encodeCB(cbSeq);
}

bool InlineCBCorrection::isAmbiguousVariant(uint64_t packedVariant) {
    if (!variantHash_) return false;
    khiter_t itVar = kh_get(cbwl, variantHash_, packedVariant);
    return itVar != kh_end(variantHash_) && kh_val(variantHash_, itVar) == 0;
}

uint32_t InlineCBCorrection::exactIndex(uint64_t packed) {
    if (!exactHash_) return 0;
    khiter_t it = kh_get(cbwl, exactHash_, packed);
    return (it != kh_end(exactHash_)) ? kh_val(exactHash_, it) : 0;
}

uint32_t InlineCBCorrection::exactIndex(const std::string &cbSeq) {
    return exactIndex(packCBForLookup(cbSeq));
}

void InlineCBCorrection::recordParentEvidence(uint32_t wlIndex1) {
    if (wlIndex1 == 0 || evidenceSize_ == 0) return;
    // Thread-local pointer to a heap-owned shard to survive thread exit
    thread_local std::vector<uint64_t>* localShard = nullptr;
    if (!localShard) {
        std::lock_guard<std::mutex> lock(evidenceMutex_);
        localShard = new std::vector<uint64_t>(evidenceSize_, 0);
        evidenceShards_.push_back(localShard);
    }
    if (wlIndex1 <= localShard->size()) {
        (*localShard)[wlIndex1 - 1] += 1;
        ++gEvidenceHits;
    }
}

void InlineCBCorrection::recordAmbiguousCB(uint64_t packedVariant, const std::string &cbSeq, const std::string &cbQual) {
    if (packedVariant == UINT64_MAX || !whitelistInitialized_) return;
    auto parentIt = variantCollisions_.find(packedVariant);
    if (parentIt == variantCollisions_.end()) return;
    // Thread-local shard pointer
    thread_local AmbigShard* shard = nullptr;
    if (!shard) {
        std::lock_guard<std::mutex> lock(ambigMutex_);
        shard = new AmbigShard();
        ambigShards_.push_back(shard);
    }
    auto idxIt = shard->index.find(packedVariant);
    if (idxIt == shard->index.end()) {
        AmbigRecord rec;
        rec.packedVariant = packedVariant;
        rec.cbSeq = cbSeq;
        rec.cbQual = cbQual;
        rec.parents = parentIt->second;
        rec.count = 1;
        shard->records.push_back(std::move(rec));
        shard->index[packedVariant] = shard->records.size() - 1;
    } else {
        shard->records[idxIt->second].count += 1;
    }
    ++ambigTotal_;
}

uint32_t InlineCBCorrection::resolveAmbiguousVariant(uint64_t packedVariant,
                                                    const std::string &cbSeq,
                                                    const std::string &cbQual,
                                                    const ParametersSolo &pSolo) {
    if (packedVariant == UINT64_MAX || !whitelistInitialized_ || !pSolo.cbWLyes || pSolo.cbWLstr.empty()) return 0;
    auto parentIt = variantCollisions_.find(packedVariant);
    if (parentIt == variantCollisions_.end()) return 0;
    const auto &parents = parentIt->second;
    if (parents.empty()) return 0;

    CbBayesianResolver resolver(pSolo.cbWLstr.size(), &pSolo.cbWLstr);
    std::vector<Candidate> candidates;
    candidates.reserve(parents.size());
    for (uint32_t idx1 : parents) {
        if (idx1 >= 1 && idx1 <= pSolo.cbWLstr.size()) {
            candidates.emplace_back(idx1, pSolo.cbWLstr[idx1 - 1]);
        }
    }
    if (candidates.empty()) return 0;

    std::unordered_map<uint32_t, uint32_t> umiCounts;
    umiCounts[0] = 1; // single-read evidence

    CBContext ctx(cbSeq, cbQual);
    BayesianResult result = resolver.resolve(ctx, candidates, umiCounts);
    if (result.status == BayesianResult::Resolved && result.bestIdx != 0) {
        ++gAmbigResolved;
        return result.bestIdx;
    }
    return 0;
}

void InlineCBCorrection::clearEvidence() {
    for (auto *shard : evidenceShards_) {
        delete shard;
    }
    evidenceShards_.clear();
    evidenceSize_ = 0;
}

void InlineCBCorrection::clearAmbiguous() {
    for (auto *shard : ambigShards_) {
        delete shard;
    }
    ambigShards_.clear();
    ambigTotal_ = 0;
}

void InlineCBCorrection::mergeAmbiguousShards(std::unordered_map<uint64_t, MergedAmbigEntry> &out) {
    out.clear();
    for (auto *shard : ambigShards_) {
        if (!shard) continue;
        for (const auto &rec : shard->records) {
            auto it = out.find(rec.packedVariant);
            if (it == out.end()) {
                MergedAmbigEntry e;
                e.cbSeq = rec.cbSeq;
                e.cbQual = rec.cbQual;
                e.parents = rec.parents;
                e.count = rec.count;
                out[rec.packedVariant] = std::move(e);
            } else {
                it->second.count += rec.count;
            }
        }
    }
}

InlineCBCorrection::AmbigResolveStats InlineCBCorrection::resolveAmbiguousMerged(
    const std::unordered_map<uint64_t, MergedAmbigEntry> &merged,
    const ParametersSolo &pSolo,
    std::vector<uint32_t> &resolvedIdx) {

    AmbigResolveStats stats;
    resolvedIdx.clear();
    resolvedIdx.reserve(merged.size());

    if (merged.empty() || !pSolo.cbWLyes || pSolo.cbWLstr.empty()) {
        return stats;
    }

    CbBayesianResolver resolver(pSolo.cbWLstr.size(), &pSolo.cbWLstr);

    for (const auto &kv : merged) {
        const auto &entry = kv.second;
        stats.total += entry.count;

        // Build candidates from parents
        std::vector<Candidate> candidates;
        candidates.reserve(entry.parents.size());
        for (uint32_t idx1 : entry.parents) {
            if (idx1 >= 1 && idx1 <= pSolo.cbWLstr.size()) {
                candidates.emplace_back(idx1, pSolo.cbWLstr[idx1 - 1]);
            }
        }
        if (candidates.empty()) {
            stats.unresolved += entry.count;
            resolvedIdx.push_back(0);
            continue;
        }

        // Evidence: use total count as a single UMI
        std::unordered_map<uint32_t, uint32_t> umiCounts;
        umiCounts[0] = entry.count;

        CBContext ctx(entry.cbSeq, entry.cbQual);
        BayesianResult result = resolver.resolve(ctx, candidates, umiCounts);
        if (result.status == BayesianResult::Resolved && result.bestIdx != 0) {
            resolvedIdx.push_back(result.bestIdx);
            stats.resolved += entry.count;
        } else if (result.status == BayesianResult::Ambiguous) {
            resolvedIdx.push_back(0);
            stats.ambiguous += entry.count;
        } else {
            resolvedIdx.push_back(0);
            stats.unresolved += entry.count;
        }
    }

    return stats;
}

// Encode CB string to packed uint64 (2 bits per base, LSB-first)
// Same encoding as encodeUMI12 but for variable-length CB
uint64_t InlineCBCorrection::encodeCB(const std::string &cbSeq) {
    if (cbSeq.empty() || cbSeq.length() > 32) {
        return UINT64_MAX;  // Max 32bp fits in 64 bits
    }
    
    uint64_t packed = 0;
    for (size_t i = 0; i < cbSeq.length(); i++) {
        char base = cbSeq[i];
        uint32_t code = 0;
        if (base == 'A' || base == 'a') code = 0;
        else if (base == 'C' || base == 'c') code = 1;
        else if (base == 'G' || base == 'g') code = 2;
        else if (base == 'T' || base == 't') code = 3;
        else return UINT64_MAX; // Invalid base
        
        packed |= ((uint64_t)code << (i * 2));  // LSB-first: first base in lowest bits
    }
    return packed;
}

// Decode packed uint64 back to DNA string
std::string InlineCBCorrection::decodeCB(uint64_t packed, int length) {
    const char bases[] = "ACGT";
    std::string cb(length, 'N');
    for (int i = 0; i < length; i++) {
        uint32_t code = (packed >> (i * 2)) & 3;
        cb[i] = bases[code];
    }
    return cb;
}

void InlineCBCorrection::initializeWhitelist(const ParametersSolo &pSolo) {
    // Cleanup prior state
    if (exactHash_) kh_destroy(cbwl, exactHash_);
    if (variantHash_) kh_destroy(cbwl, variantHash_);
    exactHash_ = kh_init(cbwl);
    variantHash_ = kh_init(cbwl);
    whitelistHash_.clear();
    variantCollisions_.clear();
    // reset evidence shards
    evidenceShards_.clear();
    evidenceSize_ = pSolo.cbWLstr.size();
    wlStrings_ = &pSolo.cbWLstr;

    // Reserve space to reduce rehashing
    if (exactHash_) kh_resize(cbwl, exactHash_, pSolo.cbWLstr.size() * 1.3);
    // Variant upper bound: len*3 per entry
    size_t maxVariants = pSolo.cbWLstr.empty() ? 0 : pSolo.cbWLstr[0].size() * 3ull * pSolo.cbWLstr.size();
    if (variantHash_) kh_resize(cbwl, variantHash_, maxVariants);

    // Build exact map (1-based index)
    for (size_t i = 0; i < pSolo.cbWLstr.size(); i++) {
        const std::string &cbStr = pSolo.cbWLstr[i];
        uint64_t packed = encodeCB(cbStr);
        if (packed == UINT64_MAX) continue;
        int absent;
        khiter_t it = kh_put(cbwl, exactHash_, packed, &absent);
        kh_val(exactHash_, it) = static_cast<uint32_t>(i + 1);
        whitelistHash_.insert(packed);
    }

    // Build 1MM variant map (unique -> index, ambiguous -> 0) and collision tracking
    if (wlStrings_ != nullptr) {
        for (khiter_t it = kh_begin(exactHash_); it != kh_end(exactHash_); ++it) {
            if (!kh_exist(exactHash_, it)) continue;
            uint64_t packed = kh_key(exactHash_, it);
            uint32_t idx1 = kh_val(exactHash_, it);
            const std::string &cbStr = (*wlStrings_)[idx1 - 1];
            int len = static_cast<int>(cbStr.length());
            for (int pos = 0; pos < len; pos++) {
                uint32_t baseCode = (packed >> (pos * 2)) & 3;
                for (uint32_t altCode = 0; altCode < 4; altCode++) {
                    if (altCode == baseCode) continue;
                    uint64_t variant = packed;
                    variant &= ~((uint64_t)3 << (pos * 2));
                    variant |= ((uint64_t)altCode << (pos * 2));
                    // Skip if this variant is an exact whitelist barcode
                    if (kh_get(cbwl, exactHash_, variant) != kh_end(exactHash_)) {
                        continue;
                    }
                    int absentVar;
                    khiter_t itVar = kh_put(cbwl, variantHash_, variant, &absentVar);
                    if (absentVar) {
                        // First claim: store as unique, no collision entry yet
                        kh_val(variantHash_, itVar) = idx1;
                    } else {
                        uint32_t prev = kh_val(variantHash_, itVar);
                        if (prev != 0) {
                            // First collision: seed collision list with previous parent
                            variantCollisions_[variant].push_back(prev);
                            kh_val(variantHash_, itVar) = 0; // mark ambiguous
                        }
                        // Record current parent for this colliding variant
                        variantCollisions_[variant].push_back(idx1);
                    }
                }
            }
        }
    }

    whitelistInitialized_ = true;
}

// Find variant matches at a specific position
int InlineCBCorrection::findVariantMatch(uint64_t packedCode, int position, std::vector<uint64_t> &variants) {
    variants.clear();
    
    // Extract the base at this position
    uint32_t baseCode = (packedCode >> (position * 2)) & 3;
    
    // Try the other 3 bases
    for (uint32_t altCode = 0; altCode < 4; altCode++) {
        if (altCode == baseCode) continue;
        
        // Create variant: clear current base, set new base
        uint64_t variant = packedCode;
        variant &= ~((uint64_t)3 << (position * 2));  // Clear current base
        variant |= ((uint64_t)altCode << (position * 2));  // Set new base
        
        // Check if variant is in whitelist
        if (whitelistHash_.find(variant) != whitelistHash_.end()) {
            variants.push_back(variant);
        }
    }
    
    return variants.size();
}

// Find closest barcodes with single mismatch
int InlineCBCorrection::findClosestBarcodes(const std::string &cbSeq, std::vector<std::string> &matches) {
    matches.clear();
    
    if (!whitelistInitialized_) {
        return 0;
    }
    
    uint64_t packedCode = encodeCB(cbSeq);
    if (packedCode == UINT64_MAX) {
        return 0;
    }
    
    // Check exact match first
    khiter_t itExact = exactHash_ ? kh_get(cbwl, exactHash_, packedCode) : kh_end(exactHash_);
    if (itExact != kh_end(exactHash_)) {
        uint32_t idx1 = kh_val(exactHash_, itExact);
        if (wlStrings_ != nullptr && idx1 >= 1 && idx1 <= wlStrings_->size()) {
            matches.push_back((*wlStrings_)[idx1 - 1]);
        } else {
            matches.push_back(cbSeq);
        }
        return 1;
    }
    
    // Single-mismatch lookup via precomputed variant map
    khiter_t itVar = variantHash_ ? kh_get(cbwl, variantHash_, packedCode) : kh_end(variantHash_);
    if (itVar == kh_end(variantHash_)) {
        return 0; // no 1MM neighbor
    }
    uint32_t idx1 = kh_val(variantHash_, itVar);
    if (idx1 == 0) {
        return 2; // ambiguous neighbors exist
    }
    // Unique 1MM neighbor
    if (wlStrings_ != nullptr && idx1 >= 1 && idx1 <= wlStrings_->size()) {
        matches.push_back((*wlStrings_)[idx1 - 1]);
    } else {
        matches.push_back(cbSeq);
    }
    return 1;
}

// Fast-path correction: exact match or single variant
int InlineCBCorrection::fastPathCorrection(const std::string &cbSeq, std::string &correctedCB) {
    correctedCB.clear();
    
    if (!whitelistInitialized_) {
        return -1;
    }
    
    // Check for Ns in sequence (Phase 1: reject)
    for (char c : cbSeq) {
        if (c == 'N' || c == 'n') {
            return -1;  // Reject sequences with N in Phase 1
        }
    }
    
    std::vector<std::string> matches;
    int nMatches = findClosestBarcodes(cbSeq, matches);
    
    if (nMatches == 0) {
        return -1;  // No match
    } else if (nMatches == 1) {
        correctedCB = matches[0];
        // Check if it's exact match or correction
        if (matches[0] == cbSeq) {
            return 0;  // Exact match
        } else {
            return 1;  // Single variant corrected
        }
    } else {
        return -1;  // Multiple matches - ambiguous (will be handled in Phase 3)
    }
}

// Check sequence for Ns and correct if possible (Phase 2)
int InlineCBCorrection::checkSequenceAndCorrectForN(const std::string &seq, int maxN, std::string &correctedSeq) {
    correctedSeq.clear();
    
    if (!whitelistInitialized_) {
        return 0;
    }
    
    // Count Ns
    int nCount = 0;
    std::vector<int> nPositions;
    for (size_t i = 0; i < seq.length(); i++) {
        if (seq[i] == 'N' || seq[i] == 'n') {
            nCount++;
            nPositions.push_back(i);
        }
    }
    
    if (nCount == 0) {
        // No Ns - use fast-path
        int result = fastPathCorrection(seq, correctedSeq);
        return (result >= 0) ? 1 : 0;
    }
    
    if (nCount > maxN) {
        return 0;  // Too many Ns
    }
    
    // Enumerate all possible sequences by replacing Ns with A/C/G/T
    int nCandidates = 1 << (nCount * 2);  // 4^nCount possibilities
    std::vector<std::string> validCandidates;
    
    for (int candidateIdx = 0; candidateIdx < nCandidates; candidateIdx++) {
        std::string candidate = seq;
        int tempIdx = candidateIdx;
        
        // Replace each N with a base based on candidateIdx
        for (int nIdx = 0; nIdx < nCount; nIdx++) {
            int baseCode = tempIdx & 3;
            tempIdx >>= 2;
            
            const char bases[] = "ACGT";
            candidate[nPositions[nIdx]] = bases[baseCode];
        }
        
        // Check if this candidate matches (exact or single variant)
        std::string corrected;
        int result = fastPathCorrection(candidate, corrected);
        if (result >= 0) {
            validCandidates.push_back(corrected);
            
            // If we find more than one valid candidate, it's ambiguous
            if (validCandidates.size() > 1) {
                return validCandidates.size();
            }
        }
    }
    
    if (validCandidates.size() == 1) {
        correctedSeq = validCandidates[0];
        return 1;
    }
    
    return 0;
}

// Clear whitelist hash
void InlineCBCorrection::clearWhitelist() {
    whitelistHash_.clear();
    if (exactHash_) { kh_destroy(cbwl, exactHash_); exactHash_ = nullptr; }
    if (variantHash_) { kh_destroy(cbwl, variantHash_); variantHash_ = nullptr; }
    variantCollisions_.clear();
    evidenceShards_.clear();
    evidenceSize_ = 0;
    wlStrings_ = nullptr;
    whitelistInitialized_ = false;
}
