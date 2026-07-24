#include "SpatialGex.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>
#include <vector>

namespace spatial_gex {
namespace {

bool checkedAdd(std::uint64_t left, std::uint64_t right, std::uint64_t &out)
{
    if (left > std::numeric_limits<std::uint64_t>::max() - right) return false;
    out = left + right;
    return true;
}

bool checkedMul(std::uint64_t left, std::uint64_t right, std::uint64_t &out)
{
    if (left != 0 && right > std::numeric_limits<std::uint64_t>::max() / left) return false;
    out = left * right;
    return true;
}

std::string lower(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char value) {
        return static_cast<char>(std::tolower(value));
    });
    return value;
}

std::vector<std::string> commaValues(const std::string &value)
{
    std::vector<std::string> output;
    std::string current;
    std::istringstream input(value);
    while (std::getline(input, current, ',')) {
        const std::size_t begin = current.find_first_not_of(" \t");
        const std::size_t end = current.find_last_not_of(" \t");
        output.push_back(begin == std::string::npos
            ? std::string() : lower(current.substr(begin, end - begin + 1)));
    }
    return output;
}

bool bytesFor(std::uint64_t count, std::uint64_t size, std::uint64_t &out,
              const char *label, std::string &error)
{
    if (checkedMul(count, size, out)) return true;
    error = std::string("spatial GEX memory arithmetic overflow for ") + label;
    return false;
}

bool sumBytes(const std::vector<std::uint64_t> &values, std::uint64_t &out,
              const char *label, std::string &error)
{
    out = 0;
    for (std::uint64_t value : values) {
        if (!checkedAdd(out, value, out)) {
            error = std::string("spatial GEX memory arithmetic overflow for ") + label;
            return false;
        }
    }
    return true;
}

} // namespace

bool parseProducts(const std::string &value, std::uint8_t &mask, std::string &error)
{
    mask = 0;
    const std::vector<std::string> values = commaValues(value);
    if (values.size() == 1 && values[0] == "all") {
        mask = ProductAll;
        return true;
    }
    std::set<std::string> seen;
    for (const std::string &item : values) {
        if (item.empty() || !seen.insert(item).second) {
            error = "spatial GEX assignment products contain an empty or duplicate value";
            return false;
        }
        if (item == "strict") mask |= ProductStrict;
        else if (item == "soft_expected") mask |= ProductSoftExpected;
        else if (item == "hard") mask |= ProductHard;
        else if (item == "gated_hard") mask |= ProductGatedHard;
        else {
            error = "unknown spatial GEX assignment product: " + item;
            return false;
        }
    }
    if (mask == 0) {
        error = "spatial GEX assignment products are empty";
        return false;
    }
    return true;
}

bool parseScales(const std::string &value, std::uint8_t &mask, std::string &error)
{
    mask = 0;
    std::set<std::string> seen;
    for (const std::string &item : commaValues(value)) {
        if (item.empty() || !seen.insert(item).second) {
            error = "spatial GEX bin sizes contain an empty or duplicate value";
            return false;
        }
        if (item == "2") mask |= Scale2um;
        else if (item == "8") mask |= Scale8um;
        else if (item == "16") mask |= Scale16um;
        else {
            error = "unsupported spatial GEX bin size: " + item;
            return false;
        }
    }
    if (mask == 0) {
        error = "spatial GEX bin sizes are empty";
        return false;
    }
    return true;
}

bool parseOverflowPolicy(const std::string &value, OverflowPolicy &policy,
                         std::string &error)
{
    const std::string normalized = lower(value);
    if (normalized == "fail") {
        policy = OverflowPolicy::Fail;
        return true;
    }
    if (normalized == "spill") {
        policy = OverflowPolicy::Spill;
        return true;
    }
    error = "unknown spatial GEX overflow policy: " + value;
    return false;
}

bool estimateMemory(const Capacity &capacity, MemoryModel &model, std::string &error)
{
    model = MemoryModel();
    if (capacity.reads == 0 || capacity.candidates == 0 || capacity.threads == 0) {
        error = "spatial GEX memory model requires positive reads, candidates, and threads";
        return false;
    }
    if (capacity.reads > std::numeric_limits<std::uint32_t>::max()
        || capacity.candidates > std::numeric_limits<std::uint32_t>::max()) {
        error = "spatial GEX compact indices require reads and candidates to fit uint32";
        return false;
    }

    std::uint64_t reads = 0, candidates = 0, readOrder = 0, supports = 0;
    std::uint64_t corrected = 0, finalMolecules = 0, matrixEntries = 0;
    std::uint64_t correctionCandidates = 0, correctionReads = 0;
    std::uint64_t threadScratch = 0, h0Counts = 0;
    if (!bytesFor(capacity.reads, sizeof(ReadEvidence), reads, "read evidence", error)
        || !bytesFor(capacity.candidates, sizeof(CandidateEvidence), candidates,
                     "candidate evidence", error)
        || !bytesFor(capacity.reads, sizeof(std::uint64_t), readOrder,
                     "read grouping indices", error)
        || !bytesFor(capacity.candidates, sizeof(PolicySupport), supports,
                     "policy support", error)
        || !bytesFor(capacity.candidates, sizeof(CorrectedSupport), corrected,
                     "corrected support", error)
        || !bytesFor(capacity.reads, sizeof(FinalMolecule), finalMolecules,
                     "final molecules", error)
        || !bytesFor(capacity.reads, 104ULL, matrixEntries,
                     "policy materialization workspace", error)
        // The current direct resolver temporarily co-owns flat raw/provisional
        // vectors and std::map nodes during soft occupancy resolution. These
        // conservative bounds include clique inputs, allocator/node overhead,
        // prior policy molecules, and the result vector. Keep them explicit
        // until those later phases acquire their own compact spill boundary.
        || !bytesFor(capacity.candidates, 224ULL, correctionCandidates,
                     "correction candidate workspace", error)
        || !bytesFor(capacity.reads, 96ULL, correctionReads,
                     "correction read workspace", error)
        || !bytesFor(capacity.threads, 8ULL * 1024ULL * 1024ULL, threadScratch,
                     "decoder thread scratch", error)
        || !bytesFor(capacity.threads, 6518ULL * sizeof(std::uint64_t), h0Counts,
                     "thread-local H0 counts", error)) {
        return false;
    }
    const std::uint64_t fixedDecoder = 64ULL * 1024ULL * 1024ULL;
    if (!sumBytes({fixedDecoder, threadScratch, h0Counts},
                  model.accumulationFixedBytes, "fixed accumulation state", error)
        || !sumBytes({model.accumulationFixedBytes, reads, candidates},
                  model.accumulationBytes, "accumulation phase", error)
        || !sumBytes({reads, candidates, readOrder, supports},
                     model.cliqueBytes, "clique phase", error)
        || !sumBytes({correctionCandidates, correctionReads, threadScratch},
                     model.correctionBytes, "correction phase", error)
        || !sumBytes({corrected, finalMolecules, readOrder},
                     model.reconciliationBytes, "reconciliation phase", error)
        || !sumBytes({finalMolecules, matrixEntries},
                     model.materializationBytes, "materialization phase", error)) {
        return false;
    }
    model.peakBytes = std::max({model.accumulationBytes, model.cliqueBytes,
                                model.correctionBytes, model.reconciliationBytes,
                                model.materializationBytes});
    return true;
}

bool memoryFits(const MemoryModel &model, std::uint64_t availableBytes,
                double fraction, std::uint64_t &budgetBytes, std::string &error)
{
    budgetBytes = 0;
    if (availableBytes == 0 || !std::isfinite(fraction) || fraction <= 0.0
        || fraction > 1.0) {
        error = "spatial GEX memory fraction must be finite in (0,1] and memory must be positive";
        return false;
    }
    const long double budget = static_cast<long double>(availableBytes) * fraction;
    if (budget > std::numeric_limits<std::uint64_t>::max()) {
        error = "spatial GEX memory budget overflow";
        return false;
    }
    budgetBytes = static_cast<std::uint64_t>(budget);
    if (model.peakBytes > budgetBytes) {
        std::ostringstream message;
        message << "spatial GEX estimated peak " << model.peakBytes
                << " exceeds configured budget " << budgetBytes;
        error = message.str();
        return false;
    }
    return true;
}

bool spillBudgetFits(const MemoryModel &model, std::uint64_t budgetBytes,
                     std::string &error)
{
    const std::uint64_t downstreamPeak = std::max(
        std::max(model.cliqueBytes, model.correctionBytes),
        std::max(model.reconciliationBytes, model.materializationBytes));
    if (downstreamPeak > budgetBytes) {
        std::ostringstream message;
        message << "integrated spatial spill can reduce read/candidate "
                << "accumulation only; post-accumulation estimate "
                << downstreamPeak << " exceeds configured budget "
                << budgetBytes;
        error = message.str();
        return false;
    }
    if (model.accumulationFixedBytes >= budgetBytes) {
        std::ostringstream message;
        message << "integrated spatial fixed accumulation state "
                << model.accumulationFixedBytes
                << " leaves no configured budget for read/candidate evidence (budget "
                << budgetBytes << ')';
        error = message.str();
        return false;
    }
    return true;
}

} // namespace spatial_gex
