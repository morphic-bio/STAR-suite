#include "SplitSeq.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace splitseq {
namespace {
Region region(const std::string& value) {
    const auto separator = value.find('_');
    if (separator == std::string::npos || separator == 0 || separator + 1 == value.size() ||
        value.find_first_not_of("0123456789", 0) != separator ||
        value.find_first_not_of("0123456789", separator + 1) != std::string::npos)
        throw std::runtime_error("positions must be zero-based offset_length: " + value);
    const auto offset = std::stoull(value.substr(0, separator));
    const auto length = std::stoull(value.substr(separator + 1));
    if (offset > 1000000 || length == 0 || length > 31)
        throw std::runtime_error("invalid split-seq position: " + value);
    return {static_cast<unsigned>(offset), static_cast<unsigned>(length)};
}

bool pack(const std::string& sequence, size_t start, size_t length, uint64_t& value) {
    value = 0;
    for (size_t i = start; i < start + length; ++i) {
        uint64_t digit;
        switch (sequence[i]) {
            case 'A': digit = 0; break;
            case 'C': digit = 1; break;
            case 'G': digit = 2; break;
            case 'T': digit = 3; break;
            default: return false;
        }
        value = (value << 2) | digit;
    }
    return true;
}
bool pack(const std::string& sequence, uint64_t& value) { return pack(sequence, 0, sequence.size(), value); }

std::string unpack(uint64_t value, unsigned length) {
    std::string result(length, 'A');
    for (unsigned i = length; i > 0; --i) { result[i-1] = "ACGT"[value & 3]; value >>= 2; }
    return result;
}
} // namespace

Decoder::Decoder(const std::vector<std::string>& positions, const std::string& umiPosition,
                 const std::vector<std::string>& whitelists) {
    if (positions.size() < 2 || positions.size() > 4 || whitelists.size() != positions.size())
        throw std::runtime_error("requires 2-4 CB positions and one allowlist per position, in cell-ID round order");
    umi_ = region(umiPosition);
    if (umi_.length > 16) throw std::runtime_error("split-seq UMI length must be 1-16 bases");
    minimumReadLength_ = umi_.offset + umi_.length;
    std::vector<Region> occupied{umi_};
    for (const auto& position : positions) {
        const auto r = region(position);
        for (const auto& other : occupied)
            if (r.offset < other.offset + other.length && other.offset < r.offset + r.length)
                throw std::runtime_error("split-seq CB/UMI positions must not overlap");
        occupied.push_back(r);
        regions_.push_back(r);
        cellLength_ += r.length;
        minimumReadLength_ = std::max(minimumReadLength_, r.offset + r.length);
    }
    if (cellLength_ > 31) throw std::runtime_error("total split-seq cell barcode length must be <=31 bases");

    for (size_t i = 0; i < whitelists.size(); ++i) {
        std::ifstream input(whitelists[i]);
        if (!input) throw std::runtime_error("cannot open split-seq allowlist: " + whitelists[i]);
        Tier tier;
        std::string line;
        size_t lineNumber = 0;
        while (std::getline(input, line)) {
            ++lineNumber;
            std::istringstream row(line);
            std::string observed, canonical, extra;
            if (!(row >> observed)) continue;
            if (!(row >> canonical)) canonical = observed;
            uint64_t key, value;
            if ((row >> extra) || observed.size() != regions_[i].length || canonical.size() != regions_[i].length ||
                !pack(observed, key) || !pack(canonical, value))
                throw std::runtime_error("invalid split-seq allowlist row " + whitelists[i] + ":" + std::to_string(lineNumber) +
                                         " (expected ACGT sequence, optionally followed by same-length canonical sequence)");
            if (!tier.exact.emplace(key, value).second)
                throw std::runtime_error("duplicate observed split-seq barcode: " + observed);
        }
        if (input.bad() || tier.exact.empty()) throw std::runtime_error("empty or unreadable split-seq allowlist: " + whitelists[i]);
        for (const auto& entry : tier.exact) {
            for (unsigned base = 0; base < regions_[i].length; ++base) {
                for (uint64_t change = 1; change < 4; ++change) {
                    const uint64_t neighbor = entry.first ^ (change << (2 * base));
                    auto inserted = tier.neighbors.emplace(neighbor, entry.first);
                    // Two observed entries one substitution away are only a
                    // conflict when they stand for different canonical
                    // sequences; aliases of one sequence resolve to it.
                    if (!inserted.second && inserted.first->second != UINT64_MAX &&
                        tier.exact.at(inserted.first->second) != entry.second)
                        inserted.first->second = UINT64_MAX;
                }
            }
        }
        tiers_.push_back(std::move(tier));
    }
}

Decoded Decoder::decode(const std::string& read) const {
    Decoded result;
    if (read.size() < minimumReadLength_) return result;
    result.umi = read.substr(umi_.offset, umi_.length);
    for (size_t i = 0; i < regions_.size(); ++i) {
        const auto& r = regions_[i];
        result.rawCell.append(read, r.offset, r.length);
        uint64_t observed;
        if (!pack(read, r.offset, r.length, observed)) { result.status = Status::InvalidBarcode; return result; }
        auto found = tiers_[i].exact.find(observed);
        if (found == tiers_[i].exact.end()) {
            const auto neighbor = tiers_[i].neighbors.find(observed);
            if (neighbor == tiers_[i].neighbors.end()) { result.status = Status::NoMatch; return result; }
            if (neighbor->second == UINT64_MAX) { result.status = Status::Ambiguous; return result; }
            found = tiers_[i].exact.find(neighbor->second);
            result.correctedRounds |= 1U << i;
        }
        // Apply the reference's observed->canonical replacement once, after correction.
        result.cell = (result.cell << (2 * r.length)) | found->second;
        result.canonicalCell += unpack(found->second, r.length);
    }
    result.status = Status::Accepted;
    return result;
}
} // namespace splitseq
