#include "SampleDetector.h"
#include "ParametersSolo.h"
#include "GlobalVariables.h"
#include <fstream>
#include <sstream>
#include "htslib/sam.h"

// Global canonical tag table (1-based; index 0 unused)
std::vector<std::string> gCanonicalTags;

std::array<uint16_t, 32> SampleDetector::tokenToSampleIdx_ = {};
std::mutex SampleDetector::tokenLUTMutex_;
std::vector<std::string> SampleDetector::canonicalByIdx_;
std::vector<std::string> SampleDetector::labelsByIdx_;

SampleDetector::SampleDetector(const ParametersSolo &p) : p_(p) {}

namespace {
int parseWhitelistIndex(const std::string &id) {
    if (id.empty()) return 0;
    if (id.size() >= 3 && (id[0]=='B' || id[0]=='b') && (id[1]=='C' || id[1]=='c')) {
        return std::atoi(id.c_str() + 2);
    }
    return std::atoi(id.c_str());
}
}

void SampleDetector::registerSampleToken(uint8_t token, uint16_t sampleIdx) {
    if (token == 0 || token >= tokenToSampleIdx_.size() || sampleIdx == 0) {
        return;
    }
    std::lock_guard<std::mutex> lock(tokenLUTMutex_);
    uint16_t &slot = tokenToSampleIdx_[token];
    if (slot == 0) {
        slot = sampleIdx;
    }
}

uint16_t SampleDetector::sampleIndexForToken(uint8_t token) {
    if (token >= tokenToSampleIdx_.size()) {
        return 0;
    }
    return tokenToSampleIdx_[token];
}

bool SampleDetector::loadWhitelist(const std::string &path) {
    canonicalToIndex_.clear();
    canonicalToWhitelistIndex_.clear();
    whitelistIndexToSequential_.clear();
    whitelistIndexToCanonical_.clear();
    whitelistIndexToLabel_.clear();
    indexToLabel_.clear();
    indexToCanonical_.clear();
    sampleCodes_.clear();
    variantCountsPerSample_.clear();
    if (path.empty() || path=="-") return false;
    std::ifstream in(path.c_str());
    if (!in.is_open()) return false;
    std::string line;
    uint32_t idx = 0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream is(line);
        std::string id, canonical;
        if (!(is >> id >> canonical)) continue;
        if (!isACGT8(canonical)) continue;
        // 1-based index
        ++idx;
        canonicalToIndex_[canonical] = idx;
        indexToLabel_.push_back(id);
        indexToCanonical_.push_back(canonical);
        int whitelistIdx = parseWhitelistIndex(id);
        if (whitelistIdx <= 0) {
            whitelistIdx = static_cast<int>(idx);
        }
        canonicalToWhitelistIndex_[canonical] = static_cast<uint32_t>(whitelistIdx);
        whitelistIndexToSequential_[static_cast<uint32_t>(whitelistIdx)] = idx;
        whitelistIndexToCanonical_[static_cast<uint32_t>(whitelistIdx)] = canonical;
        whitelistIndexToLabel_[static_cast<uint32_t>(whitelistIdx)] = id;
    }
    size_t nSamples = indexToCanonical_.size();
    sampleCodes_.assign(nSamples * 8, 0);
    variantCountsPerSample_.assign(nSamples, 0);
    // Populate global canonical and label tables (1-based)
    canonicalByIdx_ = indexToCanonical_;
    labelsByIdx_ = indexToLabel_;
    gCanonicalTags.assign(nSamples + 1, std::string()); // index 0 unused
    for (size_t i = 0; i < indexToCanonical_.size(); ++i) {
        gCanonicalTags[i + 1] = indexToCanonical_[i];
    }
    return idx>0;
}

bool SampleDetector::loadProbes(const std::string &path) {
    if (indexToCanonical_.empty()) {
        return false;
    }
    // Seed canonical sequence for each sample
    for (size_t i = 0; i < indexToCanonical_.size(); ++i) {
        uint32_t code = 0;
        if (encodeStringToCode(indexToCanonical_[i], code)) {
            sampleCodes_[i * 8] = code;
            variantCountsPerSample_[i] = 1;
        } else {
            sampleCodes_[i * 8] = 0;
            variantCountsPerSample_[i] = 0;
        }
    }
    if (path.empty() || path=="-") return false;
    std::ifstream in(path.c_str());
    if (!in.is_open()) return false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream is(line);
        std::string variant, canonical, barcodeId;
        if (!(is >> variant >> canonical >> barcodeId)) continue;
        if (!isACGT8(variant) || !isACGT8(canonical)) continue;
        uint32_t variantCode = 0;
        if (!encodeStringToCode(variant, variantCode)) {
            continue;
        }
        auto it = canonicalToIndex_.find(canonical);
        if (it == canonicalToIndex_.end()) {
            continue;
        }
        uint32_t sampleIdx = it->second;
        if (sampleIdx == 0) continue;
        size_t sampleZeroIdx = static_cast<size_t>(sampleIdx - 1);
        if (sampleZeroIdx >= variantCountsPerSample_.size()) continue;
        uint8_t &count = variantCountsPerSample_[sampleZeroIdx];
        // Skip if variant already stored for this sample
        bool duplicate = false;
        for (uint8_t k = 0; k < count; ++k) {
            if (sampleCodes_[sampleZeroIdx * 8 + k] == variantCode) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) continue;
        if (count >= 8) continue;
        size_t offset = sampleZeroIdx * 8 + count;
        sampleCodes_[offset] = variantCode;
        ++count;
    }
    return true;
}

std::string SampleDetector::canonicalForIndexStatic(uint32_t sampleIdx) {
    if (sampleIdx == 0) return std::string();
    if (sampleIdx <= canonicalByIdx_.size()) {
        return canonicalByIdx_[sampleIdx-1];
    }
    return std::string();
}

void SampleDetector::setCanonicalTable(const std::vector<std::string>& canon) {
    canonicalByIdx_ = canon;
}

std::string SampleDetector::labelForIndexStatic(uint32_t sampleIdx) {
    if (sampleIdx == 0) return std::string();
    if (sampleIdx <= labelsByIdx_.size()) {
        return labelsByIdx_[sampleIdx-1];
    }
    return std::string();
}

void SampleDetector::setLabelTable(const std::vector<std::string>& labels) {
    labelsByIdx_ = labels;
}

uint32_t SampleDetector::detectSampleIndex(const uint8_t *seqData, int32_t readLength, bool reverseStrand) const {
    if (seqData == nullptr || readLength < 8) return 0u;

    auto tryProbeAt = [&](int pos)->uint32_t {
        if (pos < 0 || pos + 8 > readLength) return 0u;
        uint32_t code = 0;
        if (!reverseStrand) {
            for (int i = 0; i < 8; ++i) {
                uint32_t nib = static_cast<uint32_t>(bam_seqi(seqData, pos + i) & 0xFu);
                if (!isACGTNib(nib)) {
                    return 0u;
                }
                code = (code << 4) | nib;
            }
        } else {
            for (int i = 0; i < 8; ++i) {
                int idx = pos + (7 - i);
                uint32_t nib = static_cast<uint32_t>(bam_seqi(seqData, idx) & 0xFu);
                if (!isACGTNib(nib)) {
                    return 0u;
                }
                uint32_t comp = complementNib(nib);
                if (comp == 0u) {
                    return 0u;
                }
                code = (code << 4) | comp;
            }
        }
        for (size_t sample = 0; sample < variantCountsPerSample_.size(); ++sample) {
            size_t base = sample * 8;
            uint8_t count = variantCountsPerSample_[sample];
            for (uint8_t k = 0; k < count; ++k) {
                if (sampleCodes_[base + k] == code) {
                    return static_cast<uint32_t>(sample + 1);
                }
            }
        }
        return 0u;
    };

    int primary = static_cast<int>(p_.sampleProbeOffset);
    uint32_t idx = tryProbeAt(primary);
    if (idx) {
        registerSampleToken(static_cast<uint8_t>(idx & 0x1Fu), static_cast<uint16_t>(idx));
        return idx;
    }

    if (p_.sampleStrictMatch) return 0u;

    if (p_.sampleSearchNearby) {
        static const int deltas[4] = {-1, 1, -2, 2};
        for (int delta : deltas) {
            idx = tryProbeAt(primary + delta);
            if (idx) {
                registerSampleToken(static_cast<uint8_t>(idx & 0x1Fu), static_cast<uint16_t>(idx));
                return idx;
            }
        }
    }

    return 0u;
}

uint32_t SampleDetector::detectSampleIndexFromBam(const char *bamRecord, uint32_t bamSize) const {
    if (bamRecord == nullptr || bamSize < 32u || !ready()) {
        return 0u;
    }

    const uint8_t *p = reinterpret_cast<const uint8_t*>(bamRecord);
    const uint8_t *ptr = p + 4; // skip block_size
    const uint8_t *end = p + bamSize;
    if (ptr >= end) {
        return 0u;
    }

    if (end - ptr < 32) {
        return 0u;
    }

    ptr += 4; // refID
    ptr += 4; // pos
    if (ptr >= end) return 0u;

    uint8_t l_read_name = *ptr++;
    ptr += 1; // mapq
    ptr += 2; // bin
    if (ptr + 2 > end) return 0u;
    uint16_t n_cigar_op = *reinterpret_cast<const uint16_t*>(ptr); ptr += 2;
    if (ptr + 2 > end) return 0u;
    uint16_t flag = *reinterpret_cast<const uint16_t*>(ptr); ptr += 2;
    if (ptr + 4 > end) return 0u;
    int32_t l_seq = *reinterpret_cast<const int32_t*>(ptr); ptr += 4;
    if (ptr + 12 > end) return 0u;
    ptr += 4; // next_refID
    ptr += 4; // next_pos
    ptr += 4; // tlen

    if (ptr + l_read_name > end) {
        return 0u;
    }
    ptr += l_read_name; // skip qname

    size_t cigarBytes = static_cast<size_t>(n_cigar_op) * 4u;
    if (ptr + cigarBytes > end) return 0u;
    ptr += cigarBytes;

    size_t seqBytes = static_cast<size_t>((l_seq + 1) / 2);
    if (ptr + seqBytes > end) return 0u;
    const uint8_t *seqPtr = ptr;
    ptr += seqBytes;

    if (ptr + l_seq > end) return 0u; // qualities

    bool reverseStrand = (flag & 0x10u) != 0;
    return detectSampleIndex(seqPtr, l_seq, reverseStrand);
}

std::string SampleDetector::labelFor(uint32_t sampleIdx) const {
    if (sampleIdx==0) return std::string();
    if (sampleIdx > indexToLabel_.size()) return std::string();
    return indexToLabel_[sampleIdx-1];
}

std::string SampleDetector::labelForWhitelistIndex(uint32_t sampleIdx) const {
    auto it = whitelistIndexToLabel_.find(sampleIdx);
    if (it == whitelistIndexToLabel_.end()) return std::string();
    return it->second;
}

std::string SampleDetector::canonicalFor(uint32_t sampleIdx) const {
    if (sampleIdx==0) return std::string();
    if (sampleIdx > indexToCanonical_.size()) return std::string();
    return indexToCanonical_[sampleIdx-1];
}

std::string SampleDetector::canonicalForWhitelistIndex(uint32_t sampleIdx) const {
    auto it = whitelistIndexToCanonical_.find(sampleIdx);
    if (it == whitelistIndexToCanonical_.end()) return std::string();
    return it->second;
}

uint32_t SampleDetector::whitelistIndexForCanonical(uint32_t sampleIdx) const {
    if (sampleIdx == 0) return 0;
    std::string canonical = canonicalFor(sampleIdx);
    if (canonical.empty()) return 0;
    auto it = canonicalToWhitelistIndex_.find(canonical);
    if (it == canonicalToWhitelistIndex_.end()) return 0;
    return it->second;
}

uint32_t SampleDetector::sequentialIndexForWhitelist(uint32_t sampleIdx) const {
    auto it = whitelistIndexToSequential_.find(sampleIdx);
    if (it == whitelistIndexToSequential_.end()) return 0u;
    return it->second;
}

bool SampleDetector::encodeStringToCode(const std::string &s, uint32_t &out) {
    out = 0;
    if (s.size() != 8) return false;
    for (char c : s) {
        int nib = seq_nt16_table[static_cast<uint8_t>(c)];
        if (!isACGTNib(static_cast<uint32_t>(nib & 0xFu))) return false;
        out = (out << 4) | static_cast<uint32_t>(nib & 0xFu);
    }
    return true;
}
