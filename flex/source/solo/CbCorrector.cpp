#include "solo/CbCorrector.h"
#include "SequenceFuns.h"
#include <algorithm>
#include <cstring>
#include <functional>
#include <limits>
#include <stdexcept>
#include <new>

CbCorrector::CbCorrector(const std::vector<std::string> &whitelist, int maxHamming,
                         bool cbqNativeOrder)
    : maxHamming_(maxHamming), cbqNativeOrder_(cbqNativeOrder) {
    
    if (whitelist.size() > std::numeric_limits<uint32_t>::max())
        throw std::length_error("CbCorrector whitelist exceeds 32-bit indices");
    whitelistSize_ = whitelist.size();
    cbLength_ = whitelist.empty() ? 0 : whitelist[0].length();
    exactMap_.reset(kh_init(cbCorrectorLookup));
    variantMap_.reset(kh_init(cbCorrectorLookup));
    ambiguousRanges_.reset(kh_init(cbCorrectorRanges));
    if (!exactMap_ || !variantMap_ || !ambiguousRanges_) throw std::bad_alloc();

    for (size_t i = 0; i < whitelist.size(); ++i) {
        PackedCB packed;
        if (!encodeCB(whitelist[i], packed)) continue;
        int absent;
        const khint_t k = kh_put(cbCorrectorLookup, exactMap_.get(), packed.key, &absent);
        if (absent < 0) throw std::bad_alloc();
        kh_val(exactMap_.get(), k) = static_cast<uint32_t>(i);
    }
    if (maxHamming_ < 1) return;

    // Pass 1: preserve unique-hit/sentinel semantics and count ambiguity sizes.
    // No per-key vectors are allocated, including during construction.
    forEachVariant(whitelist, [&](uint32_t key, uint32_t index) {
        int absent;
        const khint_t k = kh_put(cbCorrectorLookup, variantMap_.get(), key, &absent);
        if (absent < 0) throw std::bad_alloc();
        uint32_t& hit = kh_val(variantMap_.get(), k);
        if (absent) {
            hit = to1Based(index);
        } else {
            int rangeAbsent;
            const khint_t a = kh_put(cbCorrectorRanges, ambiguousRanges_.get(), key, &rangeAbsent);
            if (rangeAbsent < 0) throw std::bad_alloc();
            auto& range = kh_val(ambiguousRanges_.get(), a);
            if (rangeAbsent) {
                range.offset = 0;
                range.count = 2;
                hit = 0;
            } else {
                if (range.count == std::numeric_limits<uint32_t>::max())
                    throw std::length_error("CbCorrector candidate count overflow");
                ++range.count;
            }
        }
    });

    size_t total = 0;
    auto* ranges = ambiguousRanges_.get();
    for (khint_t k = kh_begin(ranges); k != kh_end(ranges); ++k) {
        if (!kh_exist(ranges, k)) continue;
        auto& range = kh_val(ranges, k);
        if (range.count > candidateIndices_.max_size() - total)
            throw std::length_error("CbCorrector candidate array overflow");
        range.offset = total;
        total += range.count;
        range.count = 0; // Reused as the fill cursor, then the final list length.
    }
    candidateIndices_.resize(total);

    // Pass 2: replay the original whitelist/base traversal so each candidate
    // list has exactly the former vector's order (including duplicate hits).
    forEachVariant(whitelist, [&](uint32_t key, uint32_t index) {
        const khint_t k = kh_get(cbCorrectorRanges, ranges, key);
        if (k != kh_end(ranges)) {
            auto& range = kh_val(ranges, k);
            candidateIndices_[range.offset + range.count++] = index;
        }
    });
}

uint32_t CbCorrector::packedShift(size_t pos, size_t cbLength) const {
    return static_cast<uint32_t>(2 * (cbqNativeOrder_ ? pos : cbLength - 1 - pos));
}

// Encode CB string to packed key + N mask
bool CbCorrector::encodeCB(const std::string &cb, PackedCB &packed) const {
    packed.key = 0;
    packed.nMask = 0;
    
    if (cb.length() > 16) {
        return false; // CB too long for 32-bit packing
    }
    
    for (size_t i = 0; i < cb.length(); i++) {
        char base = cb[i];
        uint32_t code = 0;
        
        if (base == 'A' || base == 'a') code = 0;
        else if (base == 'C' || base == 'c') code = 1;
        else if (base == 'G' || base == 'g') code = 2;
        else if (base == 'T' || base == 't') code = 3;
        else if (base == 'N' || base == 'n') {
            // Mark N position in mask
            packed.nMask |= (1u << i);
            code = 0; // Use A as placeholder
        } else {
            return false; // Invalid base
        }
        
        if (cbqNativeOrder_) {
            packed.key |= code << packedShift(i, cb.length());
        } else {
            // Legacy STAR order: the first base occupies the high used bits.
            packed.key = (packed.key << 2) | code;
        }
    }
    
    return true;
}

// Decode packed key back to CB string
std::string CbCorrector::decodeCB(const PackedCB &packed, size_t cbLength) const {
    const char bases[] = "ACGT";
    std::string cb(cbLength, 'N');
    
    for (size_t pos = 0; pos < cbLength; pos++) {
        uint32_t code = (packed.key >> packedShift(pos, cbLength)) & 3u;

        // Check if this position has N
        if (packed.nMask & (1u << pos)) {
            cb[pos] = 'N';
        } else {
            cb[pos] = bases[code];
        }
    }
    
    return cb;
}

// Generate all 1-hamming variants of a packed CB
void CbCorrector::generateVariants(const PackedCB &packed, std::vector<PackedCB> &variants, size_t cbLength) const {
    variants.clear();
    
    // For each position, try replacing with the other 3 nucleotides
    for (size_t pos = 0; pos < cbLength; pos++) {
        // Extract current base at this position
        uint32_t shift = packedShift(pos, cbLength);
        uint32_t currentBase = (packed.key >> shift) & 3;
        
        // Try the other 3 bases
        for (uint32_t altBase = 0; altBase < 4; altBase++) {
            if (altBase == currentBase) continue;
            
            PackedCB variant = packed;
            // Clear current base and set new base
            variant.key &= ~(3u << shift);
            variant.key |= (altBase << shift);
            variants.push_back(variant);
        }
    }
}

bool CbCorrector::hasN(const std::string &cb) const {
    for (char c : cb) {
        if (c == 'N' || c == 'n') {
            return true;
        }
    }
    return false;
}

bool CbCorrector::expandN(const PackedCB &packed, size_t cbLength, PackedCB &corrected, int &hammingDist, 
                           std::vector<PackedCB> &ambiguousSeqs) const {
    corrected.key = 0;
    corrected.nMask = 0;
    hammingDist = 255;
    ambiguousSeqs.clear();
    
    if (!packed.hasN()) {
        return false; // No Ns to expand
    }
    
    // Extract N positions from mask
    std::vector<size_t> nPositions;
    for (size_t i = 0; i < cbLength; i++) {
        if (packed.nMask & (1u << i)) {
            nPositions.push_back(i);
        }
    }
    
    // Iterative expansion: use flat array like process_features
    // Preallocate array for expanded packed keys
    const size_t MAX_EXPANDED_SEQS = 256; // Reasonable limit: 4^4 = 256
    PackedCB corrected_keys[MAX_EXPANDED_SEQS];
    size_t corrected_seq_count = 0;
    
    // Copy original packed CB into first slot
    corrected_keys[0] = packed;
    corrected_seq_count = 1;
    
    // For each N position, expand all current sequences by replacing N with A/C/G/T
    for (size_t nPos : nPositions) {
        size_t current_count = corrected_seq_count;
        uint32_t shift = packedShift(nPos, cbLength);
        
        // For each existing sequence, create 4 new sequences (one for each base)
        for (size_t i = 0; i < current_count; i++) {
            PackedCB original = corrected_keys[i];
            
            // Create 4 sequences: one for each base (A, C, G, T)
            for (uint32_t baseIdx = 0; baseIdx < 4; baseIdx++) {
                if (baseIdx == 0) {
                    // First base (A): overwrite current slot
                    corrected_keys[i].key &= ~(3u << shift);
                    corrected_keys[i].key |= (baseIdx << shift);
                    corrected_keys[i].nMask &= ~(1u << nPos); // Clear N bit
                } else {
                    // Additional bases (C, G, T): create new slots
                    if (corrected_seq_count >= MAX_EXPANDED_SEQS) {
                        return false; // Too many sequences
                    }
                    
                    // Copy original packed CB to new slot, then replace N
                    corrected_keys[corrected_seq_count] = original;
                    corrected_keys[corrected_seq_count].key &= ~(3u << shift);
                    corrected_keys[corrected_seq_count].key |= (baseIdx << shift);
                    corrected_keys[corrected_seq_count].nMask &= ~(1u << nPos); // Clear N bit
                    corrected_seq_count++;
                }
            }
        }
    }
    
    // Now check each expanded sequence against whitelist using packed keys
    // Strategy: 
    // 1. First pass: check ALL expanded sequences for exact matches (0 hash)
    // 2. Second pass: ONLY if no exact matches, check ALL for 1-hamming variants (1 hash)
    // 3. Do NOT early return - check all sequences first
    // 4. Track ambiguous hits (up to 5)
    const size_t MAX_AMBIGUOUS_HITS = 5;
    std::vector<std::pair<PackedCB, int>> bestHits; // (packedCB, hammingDist)
    
    // First pass: check ALL expanded sequences for exact matches (0 hash)
    std::vector<PackedCB> exactMatches;
    for (size_t i = 0; i < corrected_seq_count; i++) {
        const PackedCB &candidate = corrected_keys[i];
        auto exactIt = findLookup(exactMap_, candidate.key);
        if (exactIt != nullptr) {
            exactMatches.push_back(candidate);
        }
    }
    
    // If we have exact matches, use those (don't check 1-hamming hash)
    if (!exactMatches.empty()) {
        if (exactMatches.size() == 1) {
            // Exactly one exact match - return it
            corrected = exactMatches[0];
            hammingDist = 0;
            return true;
        } else {
            // Multiple exact matches - ambiguous
            for (const PackedCB &seq : exactMatches) {
                if (bestHits.size() < MAX_AMBIGUOUS_HITS) {
                    bestHits.push_back(std::make_pair(seq, 0));
                }
            }
            // Return ambiguous result
            corrected = bestHits[0].first;
            hammingDist = 0;
            ambiguousSeqs.clear();
            for (const auto &hit : bestHits) {
                if (ambiguousSeqs.size() < MAX_AMBIGUOUS_HITS) {
                    ambiguousSeqs.push_back(hit.first);
                }
            }
            return false; // Ambiguous
        }
    }
    
    // Second pass: ONLY if no exact matches, check ALL for 1-hamming variants (1 hash)
    if (maxHamming_ >= 1) {
        std::vector<PackedCB> variantMatches;
        std::vector<uint32_t> ambiguousIndices; // Collect all ambiguous indices
        
        for (size_t i = 0; i < corrected_seq_count; i++) {
            const PackedCB &candidate = corrected_keys[i];
            auto variantIt = findLookup(variantMap_, candidate.key);
            if (variantIt != nullptr) {
                if (*variantIt == 0) {
                    // Ambiguous variant - look up ambiguous hash
                    auto ambIt = ambiguousCandidates(candidate.key);
                    if (ambIt.size() != 0) {
                        // Collect all ambiguous indices
                        for (uint32_t idx : ambIt) {
                            ambiguousIndices.push_back(idx);
                        }
                    }
                } else {
                    // Non-ambiguous variant - track for uniqueness check
                    variantMatches.push_back(candidate);
                }
            }
        }
        
        // Count total ambiguous matches
        size_t totalAmbiguous = ambiguousIndices.size() + variantMatches.size();
        
        // If more than 5 ambiguous matches, return no match
        if (totalAmbiguous > 5) {
            return false; // No match
        }
        
        if (variantMatches.size() == 1 && ambiguousIndices.empty()) {
            // Exactly one non-ambiguous variant - return it
            corrected = variantMatches[0];
            hammingDist = 1;
            return true;
        } else if (totalAmbiguous > 0) {
            // Ambiguous: collect sequences for return
            // First add non-ambiguous variants
            for (const PackedCB &seq : variantMatches) {
                if (bestHits.size() < MAX_AMBIGUOUS_HITS) {
                    bestHits.push_back(std::make_pair(seq, 1));
                }
            }
            // Then add ambiguous variants (we already have indices, but need sequences)
            // For ambiguous variants, we track the sequences that generated them
            for (size_t i = 0; i < corrected_seq_count && bestHits.size() < MAX_AMBIGUOUS_HITS; i++) {
                const PackedCB &candidate = corrected_keys[i];
                auto variantIt = findLookup(variantMap_, candidate.key);
                if (variantIt != nullptr && *variantIt == 0) {
                    bestHits.push_back(std::make_pair(candidate, 1));
                }
            }
            // Return ambiguous result
            corrected = bestHits[0].first;
            hammingDist = 1;
            ambiguousSeqs.clear();
            for (const auto &hit : bestHits) {
                if (ambiguousSeqs.size() < MAX_AMBIGUOUS_HITS) {
                    ambiguousSeqs.push_back(hit.first);
                }
            }
            return false; // Ambiguous
        }
    }
    
    // No match found
    return false;
}

CbMatch CbCorrector::correct(const std::string &cb) const {
    CbMatch result;
    
    // Normalize input: convert to uppercase
    std::string normalized = cb;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), ::toupper);
    
    // Encode to packed key
    PackedCB packed;
    if (!encodeCB(normalized, packed)) {
        return result; // Invalid CB, no match
    }
    
    // Check for Ns
    if (packed.hasN()) {
        // Try N expansion using packed keys
        PackedCB expanded;
        int hammingDist;
        std::vector<PackedCB> ambiguousSeqs;
        if (expandN(packed, cbLength_, expanded, hammingDist, ambiguousSeqs)) {
            // Successfully expanded Ns to unique match (early return from expandN)
            // Look up the whitelist index using packed key
            if (hammingDist == 0) {
                // Exact match
                auto exactIt = findLookup(exactMap_, expanded.key);
                if (exactIt != nullptr) {
                    result.whitelistIdx = to1Based(*exactIt);
                    result.hammingDist = 0;
                    result.ambiguous = false;
                    return result;
                }
            } else if (hammingDist == 1) {
                // 1-hamming variant
                auto variantIt = findLookup(variantMap_, expanded.key);
                if (variantIt != nullptr && *variantIt != 0) {
                    result.whitelistIdx = *variantIt;
                    result.hammingDist = 1;
                    result.ambiguous = false;
                    return result;
                }
            }
        } else {
            // expandN returned false - could be ambiguous or no match
            // Check if we got ambiguous sequences
            if (!ambiguousSeqs.empty()) {
                // Collect all whitelist indices for ambiguous sequences
                std::vector<uint32_t> allAmbiguousIndices;
                
                for (const PackedCB &seq : ambiguousSeqs) {
                    // Check exact match first
                    auto exactIt = findLookup(exactMap_, seq.key);
                    if (exactIt != nullptr) {
                        allAmbiguousIndices.push_back(*exactIt);
                        continue;
                    }
                    
                    // Check variant match
                    if (maxHamming_ >= 1 && hammingDist == 1) {
                        auto variantIt = findLookup(variantMap_, seq.key);
                        if (variantIt != nullptr) {
                            if (*variantIt == 0) {
                                // Get the ordered range in the flat candidate array
                                auto ambIt = ambiguousCandidates(seq.key);
                                if (ambIt.size() != 0) {
                                    for (uint32_t idx : ambIt) {
                                        allAmbiguousIndices.push_back(idx);
                                    }
                                }
                            } else {
                                allAmbiguousIndices.push_back(*variantIt - 1); // Convert 1-based to 0-based
                            }
                        }
                    }
                }
                
                // If more than 5 ambiguous matches, return no match
                if (allAmbiguousIndices.size() > 5) {
                    return result; // No match
                }
                
                // Return ambiguous result (up to 5 indices)
                if (!allAmbiguousIndices.empty()) {
                    result.ambiguous = true;
                    result.hammingDist = static_cast<uint8_t>(hammingDist);
                    for (uint32_t idx : allAmbiguousIndices) {
                        result.ambiguousIdx.push_back(to1Based(idx));
                    }
                    return result;
                }
            }
        }
        // If expansion failed or no match, result stays as no match
        return result;
    }
    
    // Check exact match first using packed key
    auto exactIt = findLookup(exactMap_, packed.key);
    if (exactIt != nullptr) {
        result.whitelistIdx = to1Based(*exactIt);
        result.hammingDist = 0;
        result.ambiguous = false;
        return result;
    }
    
    // If maxHamming is 0, only exact matches allowed
    if (maxHamming_ < 1) {
        return result; // No match
    }
    
    // Check 1-hamming variants using packed key
    auto variantIt = findLookup(variantMap_, packed.key);
    if (variantIt != nullptr) {
        if (*variantIt == 0) {
            // Ambiguous: multiple WL entries match - look up ambiguous hash
            auto ambIt = ambiguousCandidates(packed.key);
            if (ambIt.size() != 0) {
                // If more than 5 ambiguous matches, return no match
                if (ambIt.size() > 5) {
                    return result; // No match
                }
                
                // Return ambiguous result (up to 5 indices)
                result.ambiguous = true;
                result.hammingDist = 1;
                for (uint32_t idx : ambIt) {
                    result.ambiguousIdx.push_back(to1Based(idx));
                }
            }
        } else {
            // Unique 1-hamming match
            result.whitelistIdx = *variantIt;
            result.hammingDist = 1;
            result.ambiguous = false;
        }
        return result;
    }
    
    // No match found
    return result;
}

bool CbCorrector::correctPackedCbq(uint32_t packedKey, uint32_t &whitelistIdx,
                                   uint8_t &hammingDist) const {
    whitelistIdx = 0;
    hammingDist = 255;
    if (!cbqNativeOrder_) {
        return false;
    }

    auto exactIt = findLookup(exactMap_, packedKey);
    if (exactIt != nullptr) {
        whitelistIdx = to1Based(*exactIt);
        hammingDist = 0;
        return true;
    }

    if (maxHamming_ < 1) {
        return false;
    }
    auto variantIt = findLookup(variantMap_, packedKey);
    if (variantIt == nullptr || *variantIt == 0) {
        return false;
    }
    whitelistIdx = *variantIt;
    hammingDist = 1;
    return true;
}

// Decode packed key to CB string (public helper)
std::string CbCorrector::decodePackedKey(uint32_t packedKey, size_t cbLength) const {
    const char bases[] = "ACGT";
    std::string cb(cbLength, 'N');
    
    for (size_t pos = 0; pos < cbLength; pos++) {
        uint32_t code = (packedKey >> packedShift(pos, cbLength)) & 3u;
        cb[pos] = bases[code];
    }
    
    return cb;
}
