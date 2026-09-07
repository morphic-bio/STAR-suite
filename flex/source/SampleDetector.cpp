#include "SampleDetector.h"
#include "ParametersSolo.h"
#include "GlobalVariables.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <utility>
#include "htslib/sam.h"

// Global canonical tag table (1-based; index 0 unused)
std::vector<std::string> gCanonicalTags;

std::array<std::atomic<uint16_t>, 32> SampleDetector::tokenToSampleIdx_{};
std::mutex SampleDetector::tokenLUTMutex_;
std::once_flag SampleDetector::canonicalByIdxOnce_;
std::once_flag SampleDetector::labelsByIdxOnce_;
std::once_flag SampleDetector::canonicalTagsOnce_;
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
    // This is called from detectSampleFromPackedTag on every read that carries
    // a sample tag, by every fused thread. The table is written at most once
    // per token, so check without the lock first: taking the process-wide
    // mutex here cost ~30% of all CPU on a 32-thread Flex run (kernel futex
    // spinlock) for a store that had already happened.
    std::atomic<uint16_t> &slot = tokenToSampleIdx_[token];
    if (slot.load(std::memory_order_relaxed) != 0) {
        return;
    }
    std::lock_guard<std::mutex> lock(tokenLUTMutex_);
    if (slot.load(std::memory_order_relaxed) == 0) {
        slot.store(sampleIdx, std::memory_order_relaxed);
    }
}

uint16_t SampleDetector::sampleIndexForToken(uint8_t token) {
    if (token >= tokenToSampleIdx_.size()) {
        return 0;
    }
    return tokenToSampleIdx_[token].load(std::memory_order_relaxed);
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
    // All worker-local detectors load the same files concurrently. Publish
    // each process-wide table once; concurrent vector<string> assignment here
    // corrupts allocator state before read processing begins.
    std::call_once(canonicalByIdxOnce_, [this]() { canonicalByIdx_ = indexToCanonical_; });
    std::call_once(labelsByIdxOnce_, [this]() { labelsByIdx_ = indexToLabel_; });
    std::call_once(canonicalTagsOnce_, [this, nSamples]() {
        gCanonicalTags.assign(nSamples + 1, std::string()); // index 0 unused
        for (size_t i = 0; i < indexToCanonical_.size(); ++i) {
            gCanonicalTags[i + 1] = indexToCanonical_[i];
        }
    });
    for (size_t i = 0; i < indexToCanonical_.size(); ++i) {
        uint32_t code = 0;
        if (encodeStringToCode(indexToCanonical_[i], code)) { sampleCodes_[i * 8] = code; variantCountsPerSample_[i] = 1; }
    }
    buildTagTable();
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
    buildTagTable();
    return true;
}

void SampleDetector::buildTagTable() {
    tagCodes_.clear();
    tagOwners_.clear();
    tagStatExact_ = tagStatMismatch_ = tagStatAmbiguous_ = 0;
    tagMaxOccupancy_ = 0;
    const size_t nSamples = variantCountsPerSample_.size();
    for (size_t s = 0; s < nSamples; ++s) {
        for (uint8_t k = 0; k < variantCountsPerSample_[s]; ++k) {
            tagCodes_.push_back(codeToIndex(sampleCodes_[s * 8 + k]));
            tagOwners_.push_back(static_cast<uint16_t>(s + 1));
        }
    }
    tagStatExact_ = static_cast<uint32_t>(tagCodes_.size());

    static const uint16_t slotCounts[3] = {1024u, 1024u, 4096u};
    for (unsigned tableIndex = 0; tableIndex < tagPieceTables_.size(); ++tableIndex) {
        TagPieceTable &table = tagPieceTables_[tableIndex];
        const uint16_t slots = slotCounts[tableIndex];
        table.offsets.assign(static_cast<size_t>(slots) + 1u, 0u);
        for (size_t id = 0; id < tagCodes_.size(); ++id) {
            ++table.offsets[pieceKey(tagCodes_[id], tableIndex) + 1u];
        }
        for (size_t slot = 1; slot < table.offsets.size(); ++slot) {
            table.offsets[slot] = static_cast<uint16_t>(table.offsets[slot] + table.offsets[slot - 1u]);
        }
        table.ids.assign(tagCodes_.size(), 0u);
        std::vector<uint16_t> cursor(table.offsets.begin(), table.offsets.end() - 1);
        for (uint16_t id = 0; id < tagCodes_.size(); ++id) {
            const uint16_t key = pieceKey(tagCodes_[id], tableIndex);
            table.ids[cursor[key]++] = id;
        }
        for (size_t slot = 0; slot < slots; ++slot) {
            const uint32_t occupancy = table.offsets[slot + 1u] - table.offsets[slot];
            tagMaxOccupancy_ = std::max(tagMaxOccupancy_, occupancy);
        }
    }

    // Preserve the diagnostic counts without retaining a direct-address
    // mismatch table. Sorting the build-time neighbour list also makes
    // duplicate variants from one owner harmless.
    const int allow = p_.sampleStrictMatch ? 0 : p_.sampleTagMismatch;
    if (allow >= 1) {
        std::vector<uint16_t> exactCodes(tagCodes_);
        std::sort(exactCodes.begin(), exactCodes.end());
        exactCodes.erase(std::unique(exactCodes.begin(), exactCodes.end()), exactCodes.end());
        std::vector<std::pair<uint16_t, uint16_t> > neighbours;
        neighbours.reserve(tagCodes_.size() * 24u);
        for (size_t id = 0; id < tagCodes_.size(); ++id) {
            const uint16_t index = tagCodes_[id];
            for (unsigned pos = 0; pos < 8; ++pos) {
                const uint16_t have = static_cast<uint16_t>((index >> (2u * pos)) & 0x3u);
                for (uint16_t base = 0; base < 4; ++base) {
                    if (base == have) continue;
                    const uint16_t neighbour = static_cast<uint16_t>(
                        (index & ~(0x3u << (2u * pos))) | (base << (2u * pos)));
                    if (!std::binary_search(exactCodes.begin(), exactCodes.end(), neighbour)) {
                        neighbours.push_back(std::make_pair(neighbour, tagOwners_[id]));
                    }
                }
            }
        }
        std::sort(neighbours.begin(), neighbours.end());
        for (size_t begin = 0; begin < neighbours.size();) {
            size_t end = begin + 1u;
            uint16_t owner = neighbours[begin].second;
            bool multipleOwners = false;
            while (end < neighbours.size() && neighbours[end].first == neighbours[begin].first) {
                if (neighbours[end].second != owner) multipleOwners = true;
                ++end;
            }
            if (multipleOwners) ++tagStatAmbiguous_;
            else ++tagStatMismatch_;
            begin = end;
        }
    }
    tagTableBuilt_ = true;
}

void SampleDetector::tagTableStats(uint32_t &exact, uint32_t &mismatch1, uint32_t &ambiguous) const {
    exact = tagStatExact_; mismatch1 = tagStatMismatch_; ambiguous = tagStatAmbiguous_;
}

uint32_t SampleDetector::lookupIndex(uint16_t index, bool exactOnly, bool *ambiguous) const {
    if (ambiguous) *ambiguous = false;
    const bool allowMismatch = !exactOnly && !p_.sampleStrictMatch && p_.sampleTagMismatch >= 1;
    uint32_t ownerMask = 0u;
    for (unsigned tableIndex = 0; tableIndex < tagPieceTables_.size(); ++tableIndex) {
        const TagPieceTable &table = tagPieceTables_[tableIndex];
        const uint16_t key = pieceKey(index, tableIndex);
        for (uint16_t at = table.offsets[key]; at < table.offsets[key + 1u]; ++at) {
            const uint16_t id = table.ids[at];
            const uint16_t x = static_cast<uint16_t>(index ^ tagCodes_[id]);
            const uint16_t differingBases = static_cast<uint16_t>((x | (x >> 1u)) & 0x5555u);
            const unsigned distance = static_cast<unsigned>(__builtin_popcount(static_cast<unsigned>(differingBases)));
            const uint16_t owner = tagOwners_[id];
            if (distance == 0u) {
                registerSampleToken(static_cast<uint8_t>(owner & 0x1Fu), owner);
                return owner; // listed sequences always decide the exact stage
            }
            if (allowMismatch && distance == 1u && owner > 0u && owner <= 32u) {
                ownerMask |= 1u << (owner - 1u);
            }
        }
    }
    if (!allowMismatch || ownerMask == 0u) return 0u;
    if ((ownerMask & (ownerMask - 1u)) != 0u) {
        if (ambiguous) *ambiguous = true;
        return 0u;
    }
    const uint16_t owner = static_cast<uint16_t>(__builtin_ctz(ownerMask) + 1u);
    registerSampleToken(static_cast<uint8_t>(owner & 0x1Fu), owner);
    return owner;
}

uint32_t SampleDetector::lookupCode(uint32_t code, bool exactOnly, bool *ambiguous) const {
    if (ambiguous) *ambiguous = false;
    if (tagTableBuilt_) return lookupIndex(codeToIndex(code), exactOnly, ambiguous);
    for (size_t sample = 0; sample < variantCountsPerSample_.size(); ++sample) {
        const size_t base = sample * 8;
        for (uint8_t k = 0; k < variantCountsPerSample_[sample]; ++k) {
            if (sampleCodes_[base + k] == code) {
                const uint32_t idx = static_cast<uint32_t>(sample + 1);
                registerSampleToken(static_cast<uint8_t>(idx & 0x1Fu), static_cast<uint16_t>(idx));
                return idx;
            }
        }
    }
    return 0u;
}

std::string SampleDetector::canonicalForIndexStatic(uint32_t sampleIdx) {
    if (sampleIdx == 0) return std::string();
    if (sampleIdx <= canonicalByIdx_.size()) {
        return canonicalByIdx_[sampleIdx-1];
    }
    return std::string();
}

void SampleDetector::setCanonicalTable(const std::vector<std::string>& canon) {
    std::call_once(canonicalByIdxOnce_, [&canon]() { canonicalByIdx_ = canon; });
}

std::string SampleDetector::labelForIndexStatic(uint32_t sampleIdx) {
    if (sampleIdx == 0) return std::string();
    if (sampleIdx <= labelsByIdx_.size()) {
        return labelsByIdx_[sampleIdx-1];
    }
    return std::string();
}

void SampleDetector::setLabelTable(const std::vector<std::string>& labels) {
    std::call_once(labelsByIdxOnce_, [&labels]() { labelsByIdx_ = labels; });
}

uint32_t SampleDetector::detectSampleIndex(const uint8_t *seqData, int32_t readLength, bool reverseStrand) const {
    if (seqData == nullptr || readLength < 8) return 0u;

    bool primaryAmbiguous = false;
    auto tryProbeAt = [&](int pos, bool exactOnly)->uint32_t {
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
        return lookupCode(code, exactOnly, exactOnly ? nullptr : &primaryAmbiguous);
    };

    int primary = static_cast<int>(p_.sampleProbeOffset);
    uint32_t idx = tryProbeAt(primary, false);
    if (idx) {
        registerSampleToken(static_cast<uint8_t>(idx & 0x1Fu), static_cast<uint16_t>(idx));
        return idx;
    }

    if (p_.sampleStrictMatch || primaryAmbiguous) return 0u;   // a tie at the tag position is final

    if (p_.sampleSearchNearby) {
        static const int deltas[4] = {-1, 1, -2, 2};
        for (int delta : deltas) {
            idx = tryProbeAt(primary + delta, true);   // a shifted tag must hit the listed table exactly
            if (idx) {
                registerSampleToken(static_cast<uint8_t>(idx & 0x1Fu), static_cast<uint16_t>(idx));
                return idx;
            }
        }
    }

    return 0u;
}

uint32_t SampleDetector::detectSampleFromPackedTag(const uint8_t *packedTag8, bool exactOnly, bool *ambiguous) const {
    if (ambiguous) *ambiguous = false;
    if (packedTag8 == nullptr || !ready()) return 0u;

    uint32_t code = 0;
    for (int i = 0; i < 8; ++i) {
        uint32_t nib = static_cast<uint32_t>(bam_seqi(packedTag8, i) & 0xFu);
        if (!isACGTNib(nib)) return 0u;
        code = (code << 4) | nib;
    }
    return lookupCode(code, exactOnly, ambiguous);
}

uint32_t SampleDetector::detectSampleFromTwoBitTag(uint16_t packedTag8, bool exactOnly, bool *ambiguous) const {
    if (ambiguous) *ambiguous = false;
    if (!ready()) return 0u;
    if (tagTableBuilt_) return lookupIndex(packedTag8, exactOnly, ambiguous);
    static const uint32_t bamNib[4] = {1u, 2u, 4u, 8u};
    uint32_t code = 0;
    for (unsigned i = 0; i < 8; ++i) {
        code = (code << 4) | bamNib[(packedTag8 >> (2U * i)) & 0x3U];
    }
    return detectSampleCode(code);
}

uint32_t SampleDetector::detectSampleCode(uint32_t code) const {
    return lookupCode(code, false, nullptr);
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
