#include "SpatialGex.h"
#include "SpatialGexSpill.h"
#include "MultiGeneUmiCr.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <mutex>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>
#include <tuple>
#include <unistd.h>
#include <unordered_map>

namespace spatial_gex {
namespace {

const std::uint32_t kGridRows = 3350;
const std::uint32_t kGridColumns = 3350;
const double kGateMinimumPosterior = 0.95;
const double kGateMinimumMargin = 0.90;
std::atomic<std::uint64_t> gPipelineGeneration(1);

bool nearlyEqual(double left, double right)
{
    const double absoluteTolerance = 1e-12;
    const double relativeTolerance = 1e-10;
    return std::fabs(left - right) <= absoluteTolerance
        + relativeTolerance * std::max(std::fabs(left), std::fabs(right));
}

bool hammingOneUmi(std::uint32_t left, std::uint32_t right)
{
    std::uint32_t mismatches = 0;
    for (std::uint32_t base = 0; base < 9; ++base) {
        if ((left & 3u) != (right & 3u) && ++mismatches > 1) return false;
        left >>= 2;
        right >>= 2;
    }
    return mismatches == 1;
}

std::uint32_t coordinateRow(std::uint32_t coordinate)
{
    return coordinate / kGridColumns;
}

std::uint32_t coordinateColumn(std::uint32_t coordinate)
{
    return coordinate % kGridColumns;
}

bool decimalComponentLess(std::uint32_t left, std::uint32_t right,
                          bool followedByUnderscore)
{
    char leftBuffer[11], rightBuffer[11];
    char *leftBegin = leftBuffer + sizeof(leftBuffer);
    char *rightBegin = rightBuffer + sizeof(rightBuffer);
    do {
        *--leftBegin = static_cast<char>('0' + left % 10);
        left /= 10;
    } while (left != 0);
    do {
        *--rightBegin = static_cast<char>('0' + right % 10);
        right /= 10;
    } while (right != 0);
    const std::size_t leftSize = static_cast<std::size_t>(
        leftBuffer + sizeof(leftBuffer) - leftBegin);
    const std::size_t rightSize = static_cast<std::size_t>(
        rightBuffer + sizeof(rightBuffer) - rightBegin);
    const std::size_t shared = std::min(leftSize, rightSize);
    const int comparison = std::memcmp(leftBegin, rightBegin, shared);
    if (comparison != 0) return comparison < 0;
    if (leftSize == rightSize) return false;
    return followedByUnderscore ? leftSize > rightSize : leftSize < rightSize;
}

bool coordinateLess(std::uint32_t left, std::uint32_t right)
{
    const std::uint32_t leftRow = coordinateRow(left);
    const std::uint32_t rightRow = coordinateRow(right);
    if (leftRow != rightRow) {
        return decimalComponentLess(leftRow, rightRow, true);
    }
    return decimalComponentLess(coordinateColumn(left), coordinateColumn(right), false);
}

bool coordinateVectorLess(const std::vector<std::uint32_t> &left,
                          const std::vector<std::uint32_t> &right)
{
    return std::lexicographical_compare(left.begin(), left.end(), right.begin(),
                                        right.end(), coordinateLess);
}

std::vector<std::uint32_t> coordinateIntersection(
    const std::vector<std::uint32_t> &left,
    const std::vector<std::uint32_t> &right)
{
    std::vector<std::uint32_t> result;
    std::size_t li = 0, ri = 0;
    while (li < left.size() && ri < right.size()) {
        if (left[li] == right[ri]) {
            result.push_back(left[li]);
            ++li;
            ++ri;
        } else if (coordinateLess(left[li], right[ri])) {
            ++li;
        } else {
            ++ri;
        }
    }
    return result;
}

void makeDirectory(const std::string &path)
{
    struct stat info;
    if (::stat(path.c_str(), &info) == 0) {
        if (!S_ISDIR(info.st_mode)) {
            throw std::runtime_error("path is not a directory: " + path);
        }
        return;
    }
    if (errno != ENOENT || ::mkdir(path.c_str(), 0755) != 0) {
        throw std::runtime_error("cannot create directory " + path + ": "
                                 + std::strerror(errno));
    }
}

void makeDirectories(const std::string &path)
{
    if (path.empty()) return;
    std::size_t begin = path[0] == '/' ? 1 : 0;
    for (;;) {
        const std::size_t end = path.find('/', begin);
        const std::string prefix = path.substr(0, end);
        if (!prefix.empty()) makeDirectory(prefix);
        if (end == std::string::npos) break;
        begin = end + 1;
    }
}

bool directoryEmpty(const std::string &path)
{
    DIR *directory = ::opendir(path.c_str());
    if (directory == NULL) return false;
    bool empty = true;
    for (dirent *entry = ::readdir(directory); entry != NULL;
         entry = ::readdir(directory)) {
        if (std::strcmp(entry->d_name, ".") != 0
            && std::strcmp(entry->d_name, "..") != 0) {
            empty = false;
            break;
        }
    }
    ::closedir(directory);
    return empty;
}

class AtomicOutput {
  public:
    explicit AtomicOutput(const std::string &path)
        : finalPath_(path), temporaryPath_(path + ".tmp"),
          stream_(temporaryPath_.c_str(), std::ios::binary)
    {
        if (!stream_) throw std::runtime_error("cannot create output: " + temporaryPath_);
    }

    std::ostream &stream() { return stream_; }

    void commit()
    {
        stream_.close();
        if (!stream_ || ::rename(temporaryPath_.c_str(), finalPath_.c_str()) != 0) {
            throw std::runtime_error("cannot commit output: " + finalPath_);
        }
        committed_ = true;
    }

    ~AtomicOutput()
    {
        if (!committed_) {
            stream_.close();
            ::remove(temporaryPath_.c_str());
        }
    }

  private:
    std::string finalPath_;
    std::string temporaryPath_;
    std::ofstream stream_;
    bool committed_ = false;
};

void atomicLinkOrCopy(const std::string &source, const std::string &target)
{
    const std::string temporary = target + ".tmp";
    ::remove(temporary.c_str());
    if (::link(source.c_str(), temporary.c_str()) != 0) {
        std::ifstream input(source.c_str(), std::ios::binary);
        std::ofstream output(temporary.c_str(), std::ios::binary);
        if (!input || !output) {
            throw std::runtime_error("cannot copy axis " + source + " to " + target);
        }
        output << input.rdbuf();
        output.close();
        if (!output) throw std::runtime_error("cannot finish axis copy: " + target);
    }
    if (::rename(temporary.c_str(), target.c_str()) != 0) {
        ::remove(temporary.c_str());
        throw std::runtime_error("cannot commit linked axis: " + target);
    }
}

bool parseUnsigned(const std::string &text, std::uint64_t &value)
{
    if (text.empty()) return false;
    char *end = NULL;
    errno = 0;
    const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
    if (errno != 0 || end == text.c_str() || *end != '\0') return false;
    value = static_cast<std::uint64_t>(parsed);
    return true;
}

bool validateBarcodeContract(const std::string &directory, std::string &error)
{
    const std::string path = directory + "/barcode_coords.tsv";
    std::ifstream input(path.c_str());
    if (!input) {
        error = "cannot open spatial barcode contract: " + path;
        return false;
    }
    std::string line;
    if (!std::getline(input, line) || line != "barcode\trow2\tcol2") {
        error = "unexpected spatial barcode contract header: " + path;
        return false;
    }
    std::uint64_t index = 0;
    while (std::getline(input, line)) {
        if (line.empty()) {
            error = "empty row in spatial barcode contract: " + path;
            return false;
        }
        const std::size_t last = line.rfind('\t');
        const std::size_t second = last == std::string::npos
            ? std::string::npos : line.rfind('\t', last - 1);
        if (second == std::string::npos || last == std::string::npos) {
            error = "malformed row in spatial barcode contract: " + path;
            return false;
        }
        std::uint64_t row = 0, column = 0;
        if (!parseUnsigned(line.substr(second + 1, last - second - 1), row)
            || !parseUnsigned(line.substr(last + 1), column)
            || row != index / kGridColumns || column != index % kGridColumns) {
            std::ostringstream message;
            message << "spatial barcode contract is not the complete row-major 3350x3350 "
                    << "universe at data row " << index;
            error = message.str();
            return false;
        }
        ++index;
    }
    const std::uint64_t expected =
        static_cast<std::uint64_t>(kGridRows) * kGridColumns;
    if (index != expected) {
        std::ostringstream message;
        message << "spatial barcode contract has " << index << " rows; expected "
                << expected;
        error = message.str();
        return false;
    }
    return true;
}

struct Bundle {
    std::uint32_t gene = 0;
    std::uint32_t rawUmi = 0;
    std::uint32_t sourceOrdinal = 0;
    std::vector<CandidateEvidence> candidates;
};

struct CliqueCandidate {
    std::uint32_t coordinate = 0;
    double posterior = 0.0;
};

struct Clique {
    std::uint32_t gene = 0;
    std::uint32_t rawUmi = 0;
    std::uint32_t candidateBegin = 0;
    std::uint32_t memberCount = 0;
    std::uint16_t candidateCount = 0;
};

struct Assignment {
    std::uint32_t gene = 0;
    std::uint32_t coordinate = 0;
    std::uint32_t rawUmi = 0;
    std::uint32_t count = 0;
};

struct IntegerRawSupport {
    std::uint32_t gene = 0;
    std::uint32_t coordinate = 0;
    std::uint32_t rawUmi = 0;
    std::uint32_t correctedUmi = 0;
    std::uint64_t count = 0;
};

struct IntegerProvisional {
    std::uint32_t coordinate = 0;
    std::uint32_t correctedUmi = 0;
    std::uint32_t gene = 0;
    std::uint64_t correctedCount = 0;
    std::uint64_t originalCount = 0;
};

struct SoftRawSupport {
    std::uint32_t gene = 0;
    std::uint32_t coordinate = 0;
    std::uint32_t rawUmi = 0;
    std::uint32_t correctedUmi = 0;
    double support = 0.0;
};

struct SoftProvisional {
    std::uint32_t coordinate = 0;
    std::uint32_t correctedUmi = 0;
    std::uint32_t gene = 0;
    double absent = 1.0;
    double originalAbsent = 1.0;
};

struct MatrixEntry {
    std::uint64_t key = 0;
    double value = 0.0;
};

bool assignmentLess(const Assignment &left, const Assignment &right)
{
    return std::tie(left.gene, left.coordinate, left.rawUmi)
        < std::tie(right.gene, right.coordinate, right.rawUmi);
}

bool integerRawLess(const IntegerRawSupport &left,
                    const IntegerRawSupport &right)
{
    return std::tie(left.gene, left.coordinate, left.rawUmi)
        < std::tie(right.gene, right.coordinate, right.rawUmi);
}

bool softRawLess(const SoftRawSupport &left, const SoftRawSupport &right)
{
    return std::tie(left.gene, left.coordinate, left.rawUmi)
        < std::tie(right.gene, right.coordinate, right.rawUmi);
}

std::uint64_t matrixKey(std::uint32_t column, std::uint32_t feature)
{
    return (static_cast<std::uint64_t>(column) << 32) | feature;
}

std::uint32_t keyColumn(std::uint64_t key)
{
    return static_cast<std::uint32_t>(key >> 32);
}

std::uint32_t keyFeature(std::uint64_t key)
{
    return static_cast<std::uint32_t>(key);
}

std::vector<MatrixEntry> aggregateMatrixEntries(std::vector<MatrixEntry> values)
{
    std::sort(values.begin(), values.end(), [](const MatrixEntry &left,
                                               const MatrixEntry &right) {
        return left.key < right.key;
    });
    std::vector<MatrixEntry> result;
    result.reserve(values.size());
    for (const MatrixEntry &value : values) {
        if (!result.empty() && result.back().key == value.key) {
            result.back().value += value.value;
        } else {
            result.push_back(value);
        }
    }
    return result;
}

std::uint64_t readUnsignedFile(const char *path)
{
    std::ifstream input(path);
    std::string value;
    if (!input || !std::getline(input, value)) return 0;
    if (value == "max") return std::numeric_limits<std::uint64_t>::max();
    std::uint64_t parsed = 0;
    return parseUnsigned(value, parsed) ? parsed : 0;
}

std::uint64_t hostAvailableMemory()
{
    std::ifstream input("/proc/meminfo");
    std::string key, unit;
    std::uint64_t value = 0;
    while (input >> key >> value >> unit) {
        if (key == "MemAvailable:") {
            return value <= std::numeric_limits<std::uint64_t>::max() / 1024
                ? value * 1024 : 0;
        }
    }
    return 0;
}

} // namespace

std::uint64_t availableMemoryBytes()
{
    std::uint64_t available = hostAvailableMemory();
    const char *limitPaths[] = {
        "/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"
    };
    const char *usagePaths[] = {
        "/sys/fs/cgroup/memory.current", "/sys/fs/cgroup/memory/memory.usage_in_bytes"
    };
    for (std::size_t index = 0; index < 2; ++index) {
        const std::uint64_t limit = readUnsignedFile(limitPaths[index]);
        const std::uint64_t usage = readUnsignedFile(usagePaths[index]);
        if (limit == 0 || limit == std::numeric_limits<std::uint64_t>::max()) continue;
        const std::uint64_t remaining = usage < limit ? limit - usage : 0;
        available = available == 0 ? remaining : std::min(available, remaining);
        break;
    }
    return available;
}

struct Pipeline::Impl {
    struct ThreadState {
        std::vector<ReadEvidence> reads;
        std::vector<CandidateEvidence> candidates;
        spatial_r1_decoder::ExactH0Counts h0;
        spatial_r1_decoder::Result currentDecoded;
        bool currentDecodedReady = false;
        std::vector<spill::Run> spillRuns;
        std::uint64_t nextSpillRun = 0;
        std::uint64_t decoded = 0;
        std::uint64_t candidateReads = 0;
        std::uint64_t uniqueGeneReads = 0;
        std::uint64_t barcodeReadsWithN = 0;
        std::uint64_t barcodeNBases = 0;
        std::uint64_t barcodeDpRecoveredReads = 0;
        std::uint64_t barcodeDpAmbiguousReads = 0;
        std::uint64_t barcodeDpUnassignedReads = 0;
        std::uint64_t barcodeUnsupportedReads = 0;
        std::uint64_t umiReadsWithN = 0;
        std::uint64_t umiReadsWithInvalidBase = 0;
    };

    PipelineConfig config;
    MemoryModel memory;
    std::uint64_t budgetBytes = 0;
    std::unique_ptr<spatial_r1_decoder::Decoder> decoder;
    std::vector<ThreadState> threads;
    PipelineSummary summary;
    std::atomic<std::uint64_t> totalJoinedReads;
    std::atomic<std::uint64_t> totalCandidateRows;
    std::atomic<std::uint64_t> residentReads;
    std::atomic<std::uint64_t> residentCandidates;
    std::atomic<std::uint64_t> peakResidentReads;
    std::atomic<std::uint64_t> peakResidentCandidates;
    std::atomic<std::uint32_t> nextThreadBinding;
    std::uint64_t generation = 0;
    std::uint64_t spillCandidateHighWater = 0;
    std::string spillDirectory;
    std::once_flag spillDirectoryOnce;
    bool finalized = false;

    explicit Impl(const PipelineConfig &value)
        : config(value), totalJoinedReads(0), totalCandidateRows(0),
          residentReads(0), residentCandidates(0), peakResidentReads(0),
          peakResidentCandidates(0), nextThreadBinding(0),
          generation(gPipelineGeneration.fetch_add(1, std::memory_order_relaxed)) {}
};

namespace {

void updatePeak(std::atomic<std::uint64_t> &peak, std::uint64_t value)
{
    std::uint64_t current = peak.load(std::memory_order_relaxed);
    while (current < value
           && !peak.compare_exchange_weak(current, value,
                                          std::memory_order_relaxed,
                                          std::memory_order_relaxed)) {
    }
}

void sortBundleCandidates(Bundle &bundle)
{
    std::sort(bundle.candidates.begin(), bundle.candidates.end(),
              [](const CandidateEvidence &left, const CandidateEvidence &right) {
        return coordinateLess(left.coordinateIndex, right.coordinateIndex);
    });
}

Bundle bundleFromMemory(const Pipeline::Impl::ThreadState &thread,
                        const ReadEvidence &read)
{
    if (static_cast<std::uint64_t>(read.candidateBegin) + read.candidateCount
        > thread.candidates.size()) {
        throw std::runtime_error("spatial read candidate span is corrupt");
    }
    Bundle bundle;
    bundle.gene = read.geneIndex;
    bundle.rawUmi = read.rawUmi;
    bundle.sourceOrdinal = read.sourceOrdinal;
    bundle.candidates.assign(
        thread.candidates.begin() + read.candidateBegin,
        thread.candidates.begin() + read.candidateBegin + read.candidateCount);
    sortBundleCandidates(bundle);
    return bundle;
}

void spillThreadSegment(Pipeline::Impl &impl, std::uint32_t threadIndex)
{
    Pipeline::Impl::ThreadState &thread = impl.threads.at(threadIndex);
    if (thread.reads.empty()) return;
    std::call_once(impl.spillDirectoryOnce, [&]() {
        makeDirectories(impl.spillDirectory);
    });
    spill::Run run;
    std::string error;
    if (!spill::writeRun(impl.spillDirectory, impl.config.sourceRevision,
                         threadIndex, thread.nextSpillRun++, thread.reads,
                         thread.candidates, run, error)) {
        throw std::runtime_error(error);
    }
    thread.spillRuns.push_back(run);
    impl.residentReads.fetch_sub(thread.reads.size(), std::memory_order_relaxed);
    impl.residentCandidates.fetch_sub(thread.candidates.size(),
                                      std::memory_order_relaxed);
    std::vector<ReadEvidence>().swap(thread.reads);
    std::vector<CandidateEvidence>().swap(thread.candidates);
    thread.reads.reserve(static_cast<std::size_t>(impl.spillCandidateHighWater));
    thread.candidates.reserve(static_cast<std::size_t>(impl.spillCandidateHighWater));
}

void cleanupSpills(Pipeline::Impl &impl)
{
    for (Pipeline::Impl::ThreadState &thread : impl.threads) {
        for (const spill::Run &run : thread.spillRuns) {
            if (std::remove(run.path.c_str()) != 0 && errno != ENOENT) {
                throw std::runtime_error("cannot remove completed spatial spill "
                                         + run.path + ": " + std::strerror(errno));
            }
        }
        thread.spillRuns.clear();
    }
    if (!impl.spillDirectory.empty()
        && ::rmdir(impl.spillDirectory.c_str()) != 0 && errno != ENOENT) {
        throw std::runtime_error("cannot remove spatial spill directory "
                                 + impl.spillDirectory + ": " + std::strerror(errno));
    }
}

double candidateLikelihood(const Bundle &read, std::uint32_t coordinate)
{
    for (const CandidateEvidence &candidate : read.candidates) {
        if (candidate.coordinateIndex == coordinate) {
            return candidate.logSequenceLikelihood;
        }
    }
    throw std::runtime_error("clique intersection lost a candidate likelihood");
}

void appendCliquesForGroup(
    const std::vector<Bundle> &reads, std::size_t begin, std::size_t end,
    const spatial_r1_decoder::ExactH0Counts &h0, std::vector<Clique> &cliques,
    std::vector<CliqueCandidate> &candidates, std::vector<std::uint64_t> &coordinateBits,
    std::vector<std::uint8_t> &featureUsed)
{
    struct Partition {
        std::vector<std::uint32_t> candidates;
        std::vector<std::size_t> members;
    };
    std::vector<Partition> partitions;
    for (std::size_t readIndex = begin; readIndex < end; ++readIndex) {
        std::vector<std::uint32_t> names;
        names.reserve(reads[readIndex].candidates.size());
        for (const CandidateEvidence &candidate : reads[readIndex].candidates) {
            names.push_back(candidate.coordinateIndex);
        }
        bool haveChoice = false;
        std::vector<std::uint32_t> bestIntersection;
        std::size_t bestPartition = 0;
        for (std::size_t partitionIndex = 0; partitionIndex < partitions.size();
             ++partitionIndex) {
            std::vector<std::uint32_t> overlap = coordinateIntersection(
                partitions[partitionIndex].candidates, names);
            if (overlap.empty()) continue;
            if (!haveChoice || coordinateVectorLess(overlap, bestIntersection)
                || (!coordinateVectorLess(bestIntersection, overlap)
                    && !coordinateVectorLess(overlap, bestIntersection)
                    && partitionIndex < bestPartition)) {
                haveChoice = true;
                bestIntersection.swap(overlap);
                bestPartition = partitionIndex;
            }
        }
        if (haveChoice) {
            partitions[bestPartition].candidates.swap(bestIntersection);
            partitions[bestPartition].members.push_back(readIndex);
        } else {
            Partition partition;
            partition.candidates.swap(names);
            partition.members.push_back(readIndex);
            partitions.push_back(std::move(partition));
        }
    }

    for (const Partition &partition : partitions) {
        if (partition.candidates.empty()
            || partition.candidates.size() > std::numeric_limits<std::uint16_t>::max()
            || partition.members.size() > std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("invalid spatial read clique cardinality");
        }
        Clique clique;
        clique.gene = reads[begin].gene;
        clique.rawUmi = reads[begin].rawUmi;
        clique.candidateBegin = static_cast<std::uint32_t>(candidates.size());
        clique.candidateCount = static_cast<std::uint16_t>(partition.candidates.size());
        clique.memberCount = static_cast<std::uint32_t>(partition.members.size());
        std::vector<double> evidence;
        evidence.reserve(partition.candidates.size());
        for (std::uint32_t coordinate : partition.candidates) {
            double likelihood = 0.0;
            for (std::size_t member : partition.members) {
                likelihood += candidateLikelihood(reads[member], coordinate);
            }
            const std::uint32_t row = coordinateRow(coordinate);
            const std::uint32_t column = coordinateColumn(coordinate);
            if (row >= h0.bc2.size() || column >= h0.bc1.size()) {
                throw std::runtime_error("clique coordinate lies outside H0 prior axes");
            }
            const std::uint64_t left = h0.bc1[column] + 1;
            const std::uint64_t right = h0.bc2[row] + 1;
            if (right != 0 && left > std::numeric_limits<std::uint64_t>::max() / right) {
                throw std::runtime_error("spatial H0 prior multiplication overflow");
            }
            const std::uint64_t exactReadCount = left * right - 1;
            evidence.push_back(likelihood + std::log(static_cast<double>(exactReadCount) + 1.0));
        }
        const double maximum = *std::max_element(evidence.begin(), evidence.end());
        double total = 0.0;
        for (double &value : evidence) {
            value = std::exp(value - maximum);
            total += value;
        }
        if (!std::isfinite(total) || total <= 0.0) {
            throw std::runtime_error("spatial posterior normalization failed");
        }
        for (std::size_t index = 0; index < partition.candidates.size(); ++index) {
            CliqueCandidate candidate;
            candidate.coordinate = partition.candidates[index];
            candidate.posterior = evidence[index] / total;
            candidates.push_back(candidate);
            coordinateBits[candidate.coordinate / 64] |=
                1ULL << (candidate.coordinate % 64);
        }
        featureUsed[clique.gene] = 1;
        cliques.push_back(clique);
    }
}

void appendCompleteGroup(const std::vector<Bundle> &group,
                         const spatial_r1_decoder::ExactH0Counts &h0,
                         std::vector<Clique> &cliques,
                         std::vector<CliqueCandidate> &candidates,
                         std::vector<std::uint64_t> &coordinateBits,
                         std::vector<std::uint8_t> &featureUsed)
{
    if (group.empty()) return;
    if (group.front().gene >= featureUsed.size()) {
        throw std::runtime_error("spatial GeneFull index lies outside transcriptome axis");
    }
    appendCliquesForGroup(group, 0, group.size(), h0, cliques, candidates,
                          coordinateBits, featureUsed);
}

void buildCliquesMemory(Pipeline::Impl &impl,
                        const spatial_r1_decoder::ExactH0Counts &h0,
                        std::vector<Clique> &cliques,
                        std::vector<CliqueCandidate> &candidates,
                        std::vector<std::uint64_t> &coordinateBits,
                        std::vector<std::uint8_t> &featureUsed)
{
    struct ReadReference {
        std::uint32_t thread;
        std::uint32_t index;
    };
    std::vector<ReadReference> order;
    order.reserve(static_cast<std::size_t>(impl.summary.joinedReads));
    for (std::size_t threadIndex = 0; threadIndex < impl.threads.size(); ++threadIndex) {
        const Pipeline::Impl::ThreadState &thread = impl.threads[threadIndex];
        if (threadIndex > std::numeric_limits<std::uint32_t>::max()
            || thread.reads.size() > std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("spatial read reference exceeds compact uint32 range");
        }
        for (std::size_t index = 0; index < thread.reads.size(); ++index) {
            ReadReference reference;
            reference.thread = static_cast<std::uint32_t>(threadIndex);
            reference.index = static_cast<std::uint32_t>(index);
            order.push_back(reference);
        }
    }
    const auto readAt = [&](const ReadReference &reference) -> const ReadEvidence & {
        return impl.threads[reference.thread].reads[reference.index];
    };
    std::sort(order.begin(), order.end(), [&](const ReadReference &left,
                                               const ReadReference &right) {
        const ReadEvidence &a = readAt(left);
        const ReadEvidence &b = readAt(right);
        if (std::tie(a.geneIndex, a.rawUmi, a.sourceOrdinal)
            != std::tie(b.geneIndex, b.rawUmi, b.sourceOrdinal)) {
            return std::tie(a.geneIndex, a.rawUmi, a.sourceOrdinal)
                < std::tie(b.geneIndex, b.rawUmi, b.sourceOrdinal);
        }
        return std::tie(left.thread, left.index) < std::tie(right.thread, right.index);
    });

    cliques.reserve(order.size());
    candidates.reserve(static_cast<std::size_t>(impl.summary.candidateRows));
    std::size_t begin = 0;
    while (begin < order.size()) {
        const ReadEvidence &first = readAt(order[begin]);
        std::size_t end = begin + 1;
        while (end < order.size()) {
            const ReadEvidence &next = readAt(order[end]);
            if (next.geneIndex != first.geneIndex || next.rawUmi != first.rawUmi) break;
            ++end;
        }
        std::vector<Bundle> group;
        group.reserve(end - begin);
        for (std::size_t index = begin; index < end; ++index) {
            const ReadReference &reference = order[index];
            group.push_back(bundleFromMemory(impl.threads[reference.thread],
                                             readAt(reference)));
        }
        appendCompleteGroup(group, h0, cliques, candidates,
                            coordinateBits, featureUsed);
        begin = end;
    }
    impl.summary.readCliques = cliques.size();
}

void buildCliquesSpill(Pipeline::Impl &impl,
                       const spatial_r1_decoder::ExactH0Counts &h0,
                       std::vector<Clique> &cliques,
                       std::vector<CliqueCandidate> &candidates,
                       std::vector<std::uint64_t> &coordinateBits,
                       std::vector<std::uint8_t> &featureUsed)
{
    for (std::uint32_t thread = 0; thread < impl.threads.size(); ++thread) {
        spillThreadSegment(impl, thread);
    }
    struct CursorState {
        std::unique_ptr<spill::Cursor> cursor;
        ReadEvidence read;
        std::vector<CandidateEvidence> candidates;
        bool active = false;
    };
    std::vector<CursorState> cursors;
    for (const Pipeline::Impl::ThreadState &thread : impl.threads) {
        for (const spill::Run &run : thread.spillRuns) {
            std::string error;
            CursorState state;
            state.cursor = spill::Cursor::open(run, impl.config.sourceRevision, error);
            if (!state.cursor) throw std::runtime_error(error);
            state.active = state.cursor->next(state.read, state.candidates, error);
            if (!state.active && !error.empty()) throw std::runtime_error(error);
            cursors.push_back(std::move(state));
        }
    }
    struct Node { std::size_t cursor; };
    struct Later {
        const std::vector<CursorState> *cursors;
        bool operator()(const Node &left, const Node &right) const
        {
            const ReadEvidence &a = (*cursors)[left.cursor].read;
            const ReadEvidence &b = (*cursors)[right.cursor].read;
            if (std::tie(a.geneIndex, a.rawUmi, a.sourceOrdinal)
                != std::tie(b.geneIndex, b.rawUmi, b.sourceOrdinal)) {
                return std::tie(a.geneIndex, a.rawUmi, a.sourceOrdinal)
                    > std::tie(b.geneIndex, b.rawUmi, b.sourceOrdinal);
            }
            return left.cursor > right.cursor;
        }
    };
    Later later;
    later.cursors = &cursors;
    std::priority_queue<Node, std::vector<Node>, Later> queue(later);
    for (std::size_t index = 0; index < cursors.size(); ++index) {
        if (cursors[index].active) {
            Node node;
            node.cursor = index;
            queue.push(node);
        }
    }
    cliques.reserve(static_cast<std::size_t>(impl.summary.joinedReads));
    candidates.reserve(static_cast<std::size_t>(impl.summary.candidateRows));
    while (!queue.empty()) {
        const ReadEvidence first = cursors[queue.top().cursor].read;
        std::vector<Bundle> group;
        while (!queue.empty()) {
            const Node node = queue.top();
            CursorState &state = cursors[node.cursor];
            if (state.read.geneIndex != first.geneIndex
                || state.read.rawUmi != first.rawUmi) break;
            queue.pop();
            Bundle bundle;
            bundle.gene = state.read.geneIndex;
            bundle.rawUmi = state.read.rawUmi;
            bundle.sourceOrdinal = state.read.sourceOrdinal;
            bundle.candidates.swap(state.candidates);
            sortBundleCandidates(bundle);
            group.push_back(std::move(bundle));
            std::string error;
            state.active = state.cursor->next(state.read, state.candidates, error);
            if (!state.active && !error.empty()) throw std::runtime_error(error);
            if (state.active) queue.push(node);
        }
        appendCompleteGroup(group, h0, cliques, candidates,
                            coordinateBits, featureUsed);
    }
    impl.summary.readCliques = cliques.size();
}

std::uint32_t topCandidate(const Clique &clique,
                           const std::vector<CliqueCandidate> &candidates,
                           double *topPosterior, double *margin)
{
    const std::size_t begin = clique.candidateBegin;
    std::size_t best = begin;
    for (std::size_t index = begin + 1; index < begin + clique.candidateCount; ++index) {
        const double value = candidates[index].posterior;
        const double current = candidates[best].posterior;
        if ((!nearlyEqual(value, current) && value > current)
            || (nearlyEqual(value, current)
                && coordinateLess(candidates[index].coordinate,
                                  candidates[best].coordinate))) {
            best = index;
        }
    }
    double second = 0.0;
    for (std::size_t index = begin; index < begin + clique.candidateCount; ++index) {
        if (index != best && candidates[index].posterior > second) {
            second = candidates[index].posterior;
        }
    }
    if (topPosterior != NULL) *topPosterior = candidates[best].posterior;
    if (margin != NULL) {
        *margin = nearlyEqual(candidates[best].posterior, second)
            ? 0.0 : candidates[best].posterior - second;
    }
    return candidates[best].coordinate;
}

std::vector<Assignment> integerAssignments(
    const std::vector<Clique> &cliques,
    const std::vector<CliqueCandidate> &candidates, std::uint32_t policy)
{
    std::vector<Assignment> result;
    result.reserve(cliques.size());
    for (const Clique &clique : cliques) {
        std::uint32_t coordinate = 0;
        if (policy == ProductStrict) {
            if (clique.candidateCount != 1) continue;
            coordinate = candidates[clique.candidateBegin].coordinate;
        } else {
            double posterior = 0.0, margin = 0.0;
            coordinate = topCandidate(clique, candidates, &posterior, &margin);
            if (policy == ProductGatedHard
                && (posterior < kGateMinimumPosterior || margin < kGateMinimumMargin)) {
                continue;
            }
        }
        Assignment assignment;
        assignment.gene = clique.gene;
        assignment.coordinate = coordinate;
        assignment.rawUmi = clique.rawUmi;
        assignment.count = clique.memberCount;
        result.push_back(assignment);
    }
    std::sort(result.begin(), result.end(), assignmentLess);
    return result;
}

std::vector<IntegerRawSupport> collapseIntegerSupport(
    const std::vector<Assignment> &assignments)
{
    std::vector<IntegerRawSupport> raw;
    raw.reserve(assignments.size());
    for (const Assignment &assignment : assignments) {
        if (!raw.empty() && raw.back().gene == assignment.gene
            && raw.back().coordinate == assignment.coordinate
            && raw.back().rawUmi == assignment.rawUmi) {
            raw.back().count += assignment.count;
        } else {
            IntegerRawSupport support;
            support.gene = assignment.gene;
            support.coordinate = assignment.coordinate;
            support.rawUmi = assignment.rawUmi;
            support.correctedUmi = assignment.rawUmi;
            support.count = assignment.count;
            raw.push_back(support);
        }
    }

    std::size_t begin = 0;
    while (begin < raw.size()) {
        std::size_t end = begin + 1;
        while (end < raw.size() && raw[end].gene == raw[begin].gene
               && raw[end].coordinate == raw[begin].coordinate) ++end;
        std::vector<std::size_t> order;
        order.reserve(end - begin);
        for (std::size_t index = begin; index < end; ++index) order.push_back(index);
        std::sort(order.begin(), order.end(), [&](std::size_t left, std::size_t right) {
            if (raw[left].count != raw[right].count) return raw[left].count > raw[right].count;
            return raw[left].rawUmi < raw[right].rawUmi;
        });
        for (std::size_t oi = 0; oi < order.size(); ++oi) {
            IntegerRawSupport &child = raw[order[oi]];
            for (std::size_t pi = 0; pi < oi; ++pi) {
                const IntegerRawSupport &parent = raw[order[pi]];
                const bool directional = child.count <= parent.count / 2 + parent.count % 2;
                if (directional && hammingOneUmi(child.rawUmi, parent.rawUmi)) {
                    child.correctedUmi = parent.correctedUmi;
                    break;
                }
            }
        }
        begin = end;
    }
    return raw;
}

std::vector<FinalMolecule> resolveInteger(
    const std::vector<Clique> &cliques,
    const std::vector<CliqueCandidate> &candidates, std::uint32_t policy)
{
    const std::vector<Assignment> assignments = integerAssignments(cliques, candidates, policy);
    const std::vector<IntegerRawSupport> raw = collapseIntegerSupport(assignments);
    std::vector<IntegerProvisional> provisional;
    provisional.reserve(raw.size());
    for (const IntegerRawSupport &support : raw) {
        IntegerProvisional value;
        value.coordinate = support.coordinate;
        value.correctedUmi = support.correctedUmi;
        value.gene = support.gene;
        value.correctedCount = support.count;
        value.originalCount = support.rawUmi == support.correctedUmi ? support.count : 0;
        provisional.push_back(value);
    }
    std::sort(provisional.begin(), provisional.end(), [](const IntegerProvisional &left,
                                                         const IntegerProvisional &right) {
        return std::tie(left.coordinate, left.correctedUmi, left.gene)
            < std::tie(right.coordinate, right.correctedUmi, right.gene);
    });
    std::vector<IntegerProvisional> aggregated;
    aggregated.reserve(provisional.size());
    for (const IntegerProvisional &value : provisional) {
        if (!aggregated.empty() && aggregated.back().coordinate == value.coordinate
            && aggregated.back().correctedUmi == value.correctedUmi
            && aggregated.back().gene == value.gene) {
            aggregated.back().correctedCount += value.correctedCount;
            aggregated.back().originalCount += value.originalCount;
        } else {
            aggregated.push_back(value);
        }
    }

    std::vector<FinalMolecule> result;
    std::size_t begin = 0;
    while (begin < aggregated.size()) {
        std::size_t end = begin + 1;
        while (end < aggregated.size()
               && aggregated[end].coordinate == aggregated[begin].coordinate
               && aggregated[end].correctedUmi == aggregated[begin].correctedUmi) ++end;
        std::vector<multi_gene_umi_cr::GeneSupport> supports;
        supports.reserve(end - begin);
        for (std::size_t index = begin; index < end; ++index) {
            supports.push_back(multi_gene_umi_cr::GeneSupport(
                aggregated[index].gene, aggregated[index].correctedCount,
                aggregated[index].originalCount));
        }
        const multi_gene_umi_cr::Result resolution = multi_gene_umi_cr::resolve(supports);
        if (resolution.accepted) {
            FinalMolecule molecule;
            molecule.geneIndex = resolution.gene;
            molecule.coordinateIndex = aggregated[begin].coordinate;
            molecule.correctedUmi = aggregated[begin].correctedUmi;
            molecule.policy = policy;
            molecule.weight = 1.0;
            result.push_back(molecule);
        }
        begin = end;
    }
    return result;
}

std::vector<SoftRawSupport> buildSoftRawSupport(
    const std::vector<Clique> &cliques,
    const std::vector<CliqueCandidate> &candidates)
{
    std::vector<SoftRawSupport> raw;
    raw.reserve(candidates.size());
    for (const Clique &clique : cliques) {
        for (std::size_t index = clique.candidateBegin;
             index < clique.candidateBegin + clique.candidateCount; ++index) {
            SoftRawSupport value;
            value.gene = clique.gene;
            value.coordinate = candidates[index].coordinate;
            value.rawUmi = clique.rawUmi;
            value.correctedUmi = clique.rawUmi;
            value.support = candidates[index].posterior;
            raw.push_back(value);
        }
    }
    std::sort(raw.begin(), raw.end(), softRawLess);
    std::vector<SoftRawSupport> aggregated;
    aggregated.reserve(raw.size());
    for (const SoftRawSupport &value : raw) {
        if (!aggregated.empty() && aggregated.back().gene == value.gene
            && aggregated.back().coordinate == value.coordinate
            && aggregated.back().rawUmi == value.rawUmi) {
            aggregated.back().support += value.support;
        } else {
            aggregated.push_back(value);
        }
    }
    std::size_t begin = 0;
    while (begin < aggregated.size()) {
        std::size_t end = begin + 1;
        while (end < aggregated.size() && aggregated[end].gene == aggregated[begin].gene
               && aggregated[end].coordinate == aggregated[begin].coordinate) ++end;
        std::vector<std::size_t> order;
        for (std::size_t index = begin; index < end; ++index) order.push_back(index);
        std::sort(order.begin(), order.end(), [&](std::size_t left, std::size_t right) {
            if (!nearlyEqual(aggregated[left].support, aggregated[right].support)) {
                return aggregated[left].support > aggregated[right].support;
            }
            return aggregated[left].rawUmi < aggregated[right].rawUmi;
        });
        for (std::size_t oi = 0; oi < order.size(); ++oi) {
            SoftRawSupport &child = aggregated[order[oi]];
            for (std::size_t pi = 0; pi < oi; ++pi) {
                const SoftRawSupport &parent = aggregated[order[pi]];
                const bool directional = parent.support > 2.0 * child.support - 1.0
                    || nearlyEqual(parent.support, 2.0 * child.support - 1.0);
                if (directional && hammingOneUmi(child.rawUmi, parent.rawUmi)) {
                    child.correctedUmi = parent.correctedUmi;
                    break;
                }
            }
        }
        begin = end;
    }
    std::sort(aggregated.begin(), aggregated.end(), softRawLess);
    return aggregated;
}

const SoftRawSupport &findSoftSupport(const std::vector<SoftRawSupport> &raw,
                                     std::uint32_t gene, std::uint32_t coordinate,
                                     std::uint32_t umi)
{
    SoftRawSupport key;
    key.gene = gene;
    key.coordinate = coordinate;
    key.rawUmi = umi;
    const std::vector<SoftRawSupport>::const_iterator found = std::lower_bound(
        raw.begin(), raw.end(), key, softRawLess);
    if (found == raw.end() || found->gene != gene || found->coordinate != coordinate
        || found->rawUmi != umi) {
        throw std::runtime_error("soft UMI correction lookup failed");
    }
    return *found;
}

std::vector<FinalMolecule> resolveSoft(
    const std::vector<Clique> &cliques,
    const std::vector<CliqueCandidate> &candidates)
{
    const std::vector<SoftRawSupport> raw = buildSoftRawSupport(cliques, candidates);
    std::map<std::tuple<std::uint32_t, std::uint32_t, std::uint32_t>, SoftProvisional> values;
    for (const Clique &clique : cliques) {
        for (std::size_t index = clique.candidateBegin;
             index < clique.candidateBegin + clique.candidateCount; ++index) {
            const CliqueCandidate &candidate = candidates[index];
            const SoftRawSupport &support = findSoftSupport(
                raw, clique.gene, candidate.coordinate, clique.rawUmi);
            const std::tuple<std::uint32_t, std::uint32_t, std::uint32_t> key(
                candidate.coordinate, support.correctedUmi, clique.gene);
            SoftProvisional &value = values[key];
            value.coordinate = candidate.coordinate;
            value.correctedUmi = support.correctedUmi;
            value.gene = clique.gene;
            value.absent *= 1.0 - candidate.posterior;
            if (clique.rawUmi == support.correctedUmi) {
                value.originalAbsent *= 1.0 - candidate.posterior;
            }
        }
    }
    std::vector<SoftProvisional> provisional;
    provisional.reserve(values.size());
    for (const std::pair<const std::tuple<std::uint32_t, std::uint32_t, std::uint32_t>,
                         SoftProvisional> &entry : values) {
        provisional.push_back(entry.second);
    }
    std::vector<FinalMolecule> result;
    std::size_t begin = 0;
    while (begin < provisional.size()) {
        std::size_t end = begin + 1;
        while (end < provisional.size()
               && provisional[end].coordinate == provisional[begin].coordinate
               && provisional[end].correctedUmi == provisional[begin].correctedUmi) ++end;
        std::size_t winner = begin;
        bool tied = false;
        for (std::size_t index = begin + 1; index < end; ++index) {
            const double value = 1.0 - provisional[index].absent;
            const double best = 1.0 - provisional[winner].absent;
            if (!nearlyEqual(value, best) && value > best) {
                winner = index;
                tied = false;
            } else if (nearlyEqual(value, best)) {
                tied = true;
            }
        }
        if (!tied) {
            const double winnerOriginal = 1.0 - provisional[winner].originalAbsent;
            bool originalDominance = false;
            for (std::size_t index = begin; index < end; ++index) {
                const double original = 1.0 - provisional[index].originalAbsent;
                if (!nearlyEqual(original, winnerOriginal) && original > winnerOriginal) {
                    originalDominance = true;
                    break;
                }
            }
            if (!originalDominance) {
                FinalMolecule molecule;
                molecule.geneIndex = provisional[winner].gene;
                molecule.coordinateIndex = provisional[winner].coordinate;
                molecule.correctedUmi = provisional[winner].correctedUmi;
                molecule.policy = ProductSoftExpected;
                molecule.weight = 1.0 - provisional[winner].absent;
                result.push_back(molecule);
            }
        }
        begin = end;
    }
    return result;
}

struct Axis {
    std::string name;
    std::uint32_t factor = 1;
    std::vector<std::uint32_t> coordinates;
    std::unordered_map<std::uint32_t, std::uint32_t> oneBasedIndex;
};

Axis buildAxis(const std::vector<std::uint32_t> &fineCoordinates,
               std::uint32_t factor, const std::string &name)
{
    std::set<std::uint32_t> unique;
    for (std::uint32_t coordinate : fineCoordinates) {
        const std::uint32_t row = coordinateRow(coordinate) / factor;
        const std::uint32_t column = coordinateColumn(coordinate) / factor;
        unique.insert(row * kGridColumns + column);
    }
    Axis axis;
    axis.name = name;
    axis.factor = factor;
    axis.coordinates.assign(unique.begin(), unique.end());
    std::sort(axis.coordinates.begin(), axis.coordinates.end(), coordinateLess);
    axis.oneBasedIndex.reserve(axis.coordinates.size() * 2 + 1);
    for (std::size_t index = 0; index < axis.coordinates.size(); ++index) {
        axis.oneBasedIndex[axis.coordinates[index]] =
            static_cast<std::uint32_t>(index + 1);
    }
    return axis;
}

void writeFeatureAxis(const std::string &path,
                      const std::vector<std::string> &features)
{
    AtomicOutput output(path);
    for (const std::string &feature : features) {
        output.stream() << feature << '\t' << feature << "\tGene Expression\n";
    }
    output.commit();
}

void writeBarcodeAxis(const std::string &path, const Axis &axis)
{
    AtomicOutput output(path);
    const std::uint32_t microns = axis.factor * 2;
    for (std::uint32_t coordinate : axis.coordinates) {
        output.stream() << "s_" << std::setw(3) << std::setfill('0') << microns
                        << "um_" << std::setfill(' ') << coordinateRow(coordinate)
                        << '_' << coordinateColumn(coordinate) << "-1\n";
    }
    output.commit();
}

std::vector<MatrixEntry> matrixAtScale(
    const std::vector<MatrixEntry> &fine, const Axis &fineAxis,
    const Axis &targetAxis)
{
    if (targetAxis.factor == 1) return fine;
    std::vector<MatrixEntry> coarse;
    coarse.reserve(fine.size());
    for (const MatrixEntry &entry : fine) {
        const std::uint32_t fineColumn = keyColumn(entry.key);
        if (fineColumn == 0 || fineColumn > fineAxis.coordinates.size()) {
            throw std::runtime_error("fine matrix column lies outside spatial axis");
        }
        const std::uint32_t coordinate = fineAxis.coordinates[fineColumn - 1];
        const std::uint32_t parent =
            (coordinateRow(coordinate) / targetAxis.factor) * kGridColumns
            + coordinateColumn(coordinate) / targetAxis.factor;
        const std::unordered_map<std::uint32_t, std::uint32_t>::const_iterator found =
            targetAxis.oneBasedIndex.find(parent);
        if (found == targetAxis.oneBasedIndex.end()) {
            throw std::runtime_error("spatial parent axis is incomplete");
        }
        MatrixEntry value;
        value.key = matrixKey(found->second, keyFeature(entry.key));
        value.value = entry.value;
        coarse.push_back(value);
    }
    return aggregateMatrixEntries(std::move(coarse));
}

double writeMatrix(const std::string &path, const std::vector<MatrixEntry> &entries,
                   std::size_t features, std::size_t barcodes, bool real)
{
    AtomicOutput output(path);
    output.stream() << "%%MatrixMarket matrix coordinate "
                    << (real ? "real" : "integer") << " general\n"
                    << "% STAR Suite molecule-first post-collapse policy matrix\n"
                    << features << ' ' << barcodes << ' ' << entries.size() << '\n'
                    << std::setprecision(17);
    double mass = 0.0;
    for (const MatrixEntry &entry : entries) {
        output.stream() << keyFeature(entry.key) << ' ' << keyColumn(entry.key) << ' ';
        if (real) {
            output.stream() << entry.value;
        } else {
            const long long rounded = std::llround(entry.value);
            if (rounded < 0 || !nearlyEqual(entry.value, static_cast<double>(rounded))) {
                throw std::runtime_error("integer spatial matrix contains a non-integer value");
            }
            output.stream() << static_cast<std::uint64_t>(rounded);
        }
        output.stream() << '\n';
        mass += entry.value;
    }
    output.commit();
    return mass;
}

struct ProductDefinition {
    std::uint32_t bit;
    const char *name;
    bool real;
};

void materialize(Pipeline::Impl &impl, const std::vector<std::string> &geneIds,
                 const std::vector<std::uint8_t> &featureUsed,
                 const std::vector<std::uint64_t> &coordinateBits,
                 const std::vector<FinalMolecule> &molecules)
{
    std::vector<std::pair<std::string, std::uint32_t> > featurePairs;
    for (std::size_t gene = 0; gene < featureUsed.size(); ++gene) {
        if (featureUsed[gene]) featurePairs.push_back(
            std::make_pair(geneIds[gene], static_cast<std::uint32_t>(gene)));
    }
    std::sort(featurePairs.begin(), featurePairs.end());
    std::vector<std::string> features;
    std::vector<std::uint32_t> featureIndex(geneIds.size(), 0);
    features.reserve(featurePairs.size());
    for (std::size_t index = 0; index < featurePairs.size(); ++index) {
        if (index > 0 && featurePairs[index - 1].first == featurePairs[index].first) {
            throw std::runtime_error("duplicate canonical gene ID in spatial feature axis");
        }
        features.push_back(featurePairs[index].first);
        featureIndex[featurePairs[index].second] = static_cast<std::uint32_t>(index + 1);
    }

    std::vector<std::uint32_t> fineCoordinates;
    for (std::size_t word = 0; word < coordinateBits.size(); ++word) {
        std::uint64_t bits = coordinateBits[word];
        while (bits != 0) {
            const unsigned offset = static_cast<unsigned>(__builtin_ctzll(bits));
            fineCoordinates.push_back(static_cast<std::uint32_t>(word * 64 + offset));
            bits &= bits - 1;
        }
    }
    std::sort(fineCoordinates.begin(), fineCoordinates.end(), coordinateLess);
    if (features.empty() || fineCoordinates.empty()) {
        throw std::runtime_error("integrated spatial read-clique universe is empty");
    }
    const Axis fine = buildAxis(fineCoordinates, 1, "square_002um");
    const Axis eight = buildAxis(fineCoordinates, 4, "square_008um");
    const Axis sixteen = buildAxis(fineCoordinates, 8, "square_016um");
    const Axis *axes[] = {&fine, &eight, &sixteen};
    const std::uint8_t scaleBits[] = {Scale2um, Scale8um, Scale16um};

    const ProductDefinition products[] = {
        {ProductStrict, "strict", false},
        {ProductSoftExpected, "soft_expected", true},
        {ProductHard, "hard", false},
        {ProductGatedHard, "gated_hard", false}
    };

    std::vector<const ProductDefinition *> selectedProducts;
    for (std::size_t product = 0; product < 4; ++product) {
        if ((impl.config.products & products[product].bit) != 0) {
            selectedProducts.push_back(&products[product]);
        }
    }
    if (selectedProducts.empty()) throw std::runtime_error("no spatial products selected");

    const std::string &root = impl.config.outputDirectory;
    struct stat info;
    if (::stat(root.c_str(), &info) == 0) {
        if (!S_ISDIR(info.st_mode) || !directoryEmpty(root)) {
            throw std::runtime_error("spatial output directory must be empty: " + root);
        }
    } else {
        makeDirectories(root);
    }
    for (const ProductDefinition *product : selectedProducts) {
        for (std::size_t scale = 0; scale < 3; ++scale) {
            if ((impl.config.scales & scaleBits[scale]) != 0) {
                makeDirectories(root + "/" + product->name + "/" + axes[scale]->name);
            }
        }
    }

    std::size_t canonicalScale = 0;
    while (canonicalScale < 3 && (impl.config.scales & scaleBits[canonicalScale]) == 0) {
        ++canonicalScale;
    }
    const std::string canonicalDirectory = root + "/" + selectedProducts.front()->name
        + "/" + axes[canonicalScale]->name;
    const std::string canonicalFeatures = canonicalDirectory + "/features.tsv";
    writeFeatureAxis(canonicalFeatures, features);
    for (std::size_t scale = 0; scale < 3; ++scale) {
        if ((impl.config.scales & scaleBits[scale]) == 0) continue;
        const std::string canonicalBarcodes = root + "/"
            + selectedProducts.front()->name + "/" + axes[scale]->name + "/barcodes.tsv";
        writeBarcodeAxis(canonicalBarcodes, *axes[scale]);
        for (const ProductDefinition *product : selectedProducts) {
            const std::string directory = root + "/" + product->name + "/" + axes[scale]->name;
            if (directory + "/features.tsv" != canonicalFeatures) {
                atomicLinkOrCopy(canonicalFeatures, directory + "/features.tsv");
            }
            if (directory + "/barcodes.tsv" != canonicalBarcodes) {
                atomicLinkOrCopy(canonicalBarcodes, directory + "/barcodes.tsv");
            }
        }
    }

    std::vector<std::uint32_t> fineIndex(
        static_cast<std::size_t>(kGridRows) * kGridColumns, 0);
    for (std::size_t index = 0; index < fine.coordinates.size(); ++index) {
        fineIndex[fine.coordinates[index]] = static_cast<std::uint32_t>(index + 1);
    }

    std::ostringstream summary;
    summary << "schema\tstar_suite.molecule_first.policy_mex.v1\n"
            << "star_suite_version\t" << impl.config.starSuiteVersion << '\n'
            << "assay\tvisium-hd\n"
            << "umi_mode\t1mm_cr\n"
            << "axes_source\tin_memory_read_cliques\n"
            << "cell_calling_order\tafter_postcollapse_integer_matrix_only\n"
            << "soft_cell_calling_allowed\tfalse\n"
            << "product\tscale\tfeatures\tbarcodes\tnnz\tmass\tmatrix_field\n";

    for (const ProductDefinition *product : selectedProducts) {
        std::vector<MatrixEntry> entries;
        for (const FinalMolecule &molecule : molecules) {
            if (molecule.policy != product->bit) continue;
            if (molecule.geneIndex >= featureIndex.size()
                || molecule.coordinateIndex >= fineIndex.size()
                || featureIndex[molecule.geneIndex] == 0
                || fineIndex[molecule.coordinateIndex] == 0) {
                throw std::runtime_error("final molecule lies outside spatial matrix axes");
            }
            MatrixEntry entry;
            entry.key = matrixKey(fineIndex[molecule.coordinateIndex],
                                  featureIndex[molecule.geneIndex]);
            entry.value = molecule.weight;
            entries.push_back(entry);
        }
        const std::vector<MatrixEntry> fineEntries = aggregateMatrixEntries(std::move(entries));
        double referenceMass = -1.0;
        for (std::size_t scale = 0; scale < 3; ++scale) {
            if ((impl.config.scales & scaleBits[scale]) == 0) continue;
            const std::vector<MatrixEntry> scaled = matrixAtScale(fineEntries, fine, *axes[scale]);
            const std::string path = root + "/" + product->name + "/"
                + axes[scale]->name + "/matrix.mtx";
            const double mass = writeMatrix(path, scaled, features.size(),
                                            axes[scale]->coordinates.size(), product->real);
            if (referenceMass < 0.0) referenceMass = mass;
            const double tolerance = 1e-9 * std::max(1.0, std::max(referenceMass, mass));
            if (std::fabs(referenceMass - mass) > tolerance) {
                throw std::runtime_error("spatial scale aggregation did not conserve mass");
            }
            summary << product->name << '\t' << axes[scale]->name << '\t'
                    << features.size() << '\t' << axes[scale]->coordinates.size() << '\t'
                    << scaled.size() << '\t' << std::setprecision(17) << mass << '\t'
                    << (product->real ? "real" : "integer") << '\n';
        }
    }
    {
        AtomicOutput output(root + "/summary.tsv");
        output.stream() << summary.str();
        output.commit();
    }
}

} // namespace

std::unique_ptr<Pipeline> Pipeline::create(const PipelineConfig &config,
                                           std::string &error)
{
    try {
        if (config.threads == 0 || config.expectedReads == 0
            || config.expectedCandidates == 0) {
            throw std::invalid_argument(
                "integrated spatial GEX requires positive thread/read/candidate capacities");
        }
        if (config.barcodeContractDirectory.empty() || config.bc1OligosPath.empty()
            || config.bc2OligosPath.empty() || config.outputDirectory.empty()) {
            throw std::invalid_argument("integrated spatial GEX input/output paths are required");
        }
        if (config.overflowPolicy == OverflowPolicy::Spill
            && config.temporaryDirectory.empty()) {
            throw std::invalid_argument("integrated spatial spill requires a temporary directory");
        }
        if (config.overflowPolicy != OverflowPolicy::Spill
            && config.spillHighWaterCandidatesPerThread != 0) {
            throw std::invalid_argument(
                "spatial spill high-water override requires overflow policy Spill");
        }
        std::unique_ptr<Impl> impl(new Impl(config));
        Capacity capacity;
        capacity.reads = config.expectedReads;
        capacity.candidates = config.expectedCandidates;
        capacity.threads = config.threads;
        if (!estimateMemory(capacity, impl->memory, error)) return std::unique_ptr<Pipeline>();
        const std::uint64_t available = availableMemoryBytes();
        std::string fitError;
        const bool fits = memoryFits(impl->memory, available, config.memoryFraction,
                                     impl->budgetBytes, fitError);
        if (!fits && config.overflowPolicy == OverflowPolicy::Fail) {
            error = fitError;
            return std::unique_ptr<Pipeline>();
        }
        if (!fits) {
            if (available == 0 || !std::isfinite(config.memoryFraction)
                || config.memoryFraction <= 0.0 || config.memoryFraction > 1.0) {
                error = fitError;
                return std::unique_ptr<Pipeline>();
            }
            impl->budgetBytes = static_cast<std::uint64_t>(
                static_cast<long double>(available) * config.memoryFraction);
        }
        if (config.overflowPolicy == OverflowPolicy::Spill) {
            if (!spillBudgetFits(impl->memory, impl->budgetBytes, error)) {
                return std::unique_ptr<Pipeline>();
            }
        }
        if (!validateBarcodeContract(config.barcodeContractDirectory, error)) {
            return std::unique_ptr<Pipeline>();
        }
        spatial_r1_decoder::Config decoderConfig;
        decoderConfig.bc1OligosPath = config.bc1OligosPath;
        decoderConfig.bc2OligosPath = config.bc2OligosPath;
        decoderConfig.gridRows = kGridRows;
        decoderConfig.gridColumns = kGridColumns;
        decoderConfig.fullStartMin = 8;
        decoderConfig.fullStartMax = 12;
        impl->decoder.reset(new spatial_r1_decoder::Decoder(decoderConfig));
        if (impl->decoder->bc1Count() != kGridColumns
            || impl->decoder->bc2Count() != kGridRows) {
            throw std::runtime_error("spatial oligo axes do not match the 3350x3350 contract");
        }
        impl->threads.resize(config.threads);
        const std::uint64_t readsPerThread =
            (config.expectedReads + config.threads - 1) / config.threads;
        const std::uint64_t candidatesPerThread =
            (config.expectedCandidates + config.threads - 1) / config.threads;
        std::uint64_t readReserve = readsPerThread;
        std::uint64_t candidateReserve = candidatesPerThread;
        if (config.overflowPolicy == OverflowPolicy::Spill) {
            const std::uint64_t bytesPerCandidate =
                sizeof(CandidateEvidence) + sizeof(ReadEvidence);
            const std::uint64_t evidenceBudget =
                impl->budgetBytes - impl->memory.accumulationFixedBytes;
            const std::uint64_t denominator =
                static_cast<std::uint64_t>(config.threads) * bytesPerCandidate;
            impl->spillCandidateHighWater = std::max<std::uint64_t>(
                1, evidenceBudget / std::max<std::uint64_t>(denominator, 1));
            impl->spillCandidateHighWater = std::min(
                impl->spillCandidateHighWater,
                std::max<std::uint64_t>(1, candidatesPerThread));
            if (config.spillHighWaterCandidatesPerThread != 0) {
                impl->spillCandidateHighWater = std::min(
                    impl->spillCandidateHighWater,
                    config.spillHighWaterCandidatesPerThread);
            }
            impl->summary.spillHighWaterCandidatesPerThread =
                impl->spillCandidateHighWater;
            impl->spillDirectory = config.temporaryDirectory + "/spatial_gex_overflow";
            readReserve = std::min(readsPerThread, impl->spillCandidateHighWater);
            candidateReserve = std::min(candidatesPerThread,
                                        impl->spillCandidateHighWater);
        }
        for (Impl::ThreadState &thread : impl->threads) {
            thread.h0.reset(kGridColumns, kGridRows);
            thread.reads.reserve(static_cast<std::size_t>(readReserve));
            thread.candidates.reserve(static_cast<std::size_t>(candidateReserve));
        }
        return std::unique_ptr<Pipeline>(new Pipeline(std::move(impl)));
    } catch (const std::exception &exception) {
        error = exception.what();
        return std::unique_ptr<Pipeline>();
    }
}

Pipeline::Pipeline(std::unique_ptr<Impl> impl) : impl_(std::move(impl)) {}
Pipeline::~Pipeline() {}

bool Pipeline::currentThreadIndex(std::uint32_t &threadIndex, std::string &error)
{
    struct Binding {
        const Impl *impl = nullptr;
        std::uint64_t generation = 0;
        std::uint32_t threadIndex = 0;
    };
    static thread_local Binding binding;
    if (binding.impl != impl_.get() || binding.generation != impl_->generation) {
        const std::uint32_t assigned = impl_->nextThreadBinding.fetch_add(
            1, std::memory_order_relaxed);
        if (assigned >= impl_->threads.size()) {
            error = "integrated spatial GEX observed more mapping threads than declared";
            return false;
        }
        binding.impl = impl_.get();
        binding.generation = impl_->generation;
        binding.threadIndex = assigned;
    }
    threadIndex = binding.threadIndex;
    return true;
}

bool Pipeline::decodeCurrentThread(const char *sequence,
                                   std::size_t sequenceLength,
                                   const char *quality,
                                   std::size_t qualityLength,
                                   std::string &error)
{
    std::uint32_t threadIndex = 0;
    if (!currentThreadIndex(threadIndex, error)) return false;
    Impl::ThreadState &thread = impl_->threads[threadIndex];
    thread.currentDecodedReady = false;
    if (!decode(threadIndex, sequence, sequenceLength, quality, qualityLength,
                thread.currentDecoded, error)) return false;
    thread.currentDecodedReady = true;
    return true;
}

bool Pipeline::appendCurrentThread(std::uint32_t geneIndex,
                                   std::uint64_t sourceOrdinal,
                                   std::string &error)
{
    std::uint32_t threadIndex = 0;
    if (!currentThreadIndex(threadIndex, error)) return false;
    Impl::ThreadState &thread = impl_->threads[threadIndex];
    if (!thread.currentDecodedReady) {
        error = "integrated spatial GEX lost current-thread decoded R1 state";
        return false;
    }
    ++thread.uniqueGeneReads;
    const bool appended = append(threadIndex, geneIndex, sourceOrdinal,
                                 thread.currentDecoded, error);
    thread.currentDecodedReady = false;
    return appended;
}

bool Pipeline::decode(std::uint32_t threadIndex, const char *sequence,
                      std::size_t sequenceLength, const char *quality,
                      std::size_t qualityLength,
                      spatial_r1_decoder::Result &result, std::string &error)
{
    if (threadIndex >= impl_->threads.size()) {
        error = "spatial decoder thread index is out of range";
        return false;
    }
    try {
        Impl::ThreadState &thread = impl_->threads[threadIndex];
        if (!impl_->decoder->decode(sequence, sequenceLength, quality, qualityLength,
                                    result, &thread.h0, error)) return false;
        ++thread.decoded;
        if (!result.candidates.empty()) ++thread.candidateReads;
        if (result.barcodeHadN) {
            ++thread.barcodeReadsWithN;
            thread.barcodeNBases += result.barcodeNCount;
        }
        if (result.barcodeDpChecked) {
            if (result.decoderAssigned) {
                ++thread.barcodeDpRecoveredReads;
            } else if (!result.candidates.empty()) {
                ++thread.barcodeDpAmbiguousReads;
            } else {
                ++thread.barcodeDpUnassignedReads;
            }
        }
        if (result.barcodeHadUnsupportedBase) ++thread.barcodeUnsupportedReads;
        if (result.rawUmiHadN) ++thread.umiReadsWithN;
        if (!result.rawUmiValid) ++thread.umiReadsWithInvalidBase;
        return true;
    } catch (const std::exception &exception) {
        error = exception.what();
        return false;
    }
}

void Pipeline::noteUniqueGene(std::uint32_t threadIndex)
{
    if (threadIndex < impl_->threads.size()) {
        ++impl_->threads[threadIndex].uniqueGeneReads;
    }
}

bool Pipeline::append(std::uint32_t threadIndex, std::uint32_t geneIndex,
                      std::uint64_t sourceOrdinal,
                      const spatial_r1_decoder::Result &decoded,
                      std::string &error)
{
    try {
        if (threadIndex >= impl_->threads.size()) {
            throw std::runtime_error("spatial accumulator thread index is out of range");
        }
        if (!decoded.rawUmiValid || decoded.candidates.empty()) return true;
        if (sourceOrdinal > std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("spatial source ordinal exceeds compact uint32 range");
        }
        if (decoded.candidates.size() > std::numeric_limits<std::uint16_t>::max()) {
            throw std::runtime_error("spatial candidate count exceeds compact uint16 range");
        }
        const std::uint64_t priorReads = impl_->totalJoinedReads.fetch_add(
            1, std::memory_order_relaxed);
        const std::uint64_t priorCandidates = impl_->totalCandidateRows.fetch_add(
            decoded.candidates.size(), std::memory_order_relaxed);
        if (priorReads >= impl_->config.expectedReads
            || decoded.candidates.size() > impl_->config.expectedCandidates
            || priorCandidates > impl_->config.expectedCandidates
                    - decoded.candidates.size()) {
            throw std::runtime_error("integrated spatial evidence exceeded declared capacity");
        }
        Impl::ThreadState &thread = impl_->threads[threadIndex];
        if (impl_->config.overflowPolicy == OverflowPolicy::Spill
            && !thread.reads.empty()
            && thread.candidates.size() + decoded.candidates.size()
                > impl_->spillCandidateHighWater) {
            spillThreadSegment(*impl_, threadIndex);
        }
        if (thread.candidates.size() > std::numeric_limits<std::uint32_t>::max()
            || decoded.candidates.size() > std::numeric_limits<std::uint32_t>::max()
                - thread.candidates.size()) {
            throw std::runtime_error("spatial per-thread candidate pool exceeds uint32 range");
        }
        if (thread.reads.size() >= std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("spatial per-thread read pool exceeds uint32 range");
        }
        ReadEvidence read;
        read.geneIndex = geneIndex;
        read.rawUmi = decoded.rawUmi;
        read.candidateBegin = static_cast<std::uint32_t>(thread.candidates.size());
        read.sourceOrdinal = static_cast<std::uint32_t>(sourceOrdinal);
        read.candidateCount = static_cast<std::uint16_t>(decoded.candidates.size());
        read.flags = decoded.decoderAssigned ? 1u : 0u;
        for (const spatial_r1_decoder::Candidate &source : decoded.candidates) {
            CandidateEvidence candidate;
            candidate.coordinateIndex = source.coordinateIndex;
            candidate.auditBits = source.auditBits;
            candidate.logSequenceLikelihood = source.logSequenceLikelihood;
            thread.candidates.push_back(candidate);
        }
        thread.reads.push_back(read);
        const std::uint64_t residentReads = impl_->residentReads.fetch_add(
            1, std::memory_order_relaxed) + 1;
        const std::uint64_t residentCandidates = impl_->residentCandidates.fetch_add(
            decoded.candidates.size(), std::memory_order_relaxed)
            + decoded.candidates.size();
        updatePeak(impl_->peakResidentReads, residentReads);
        updatePeak(impl_->peakResidentCandidates, residentCandidates);
        return true;
    } catch (const std::exception &exception) {
        error = exception.what();
        return false;
    }
}

bool Pipeline::finalize(const std::vector<std::string> &geneIds,
                        std::string &error)
{
    try {
        if (impl_->finalized) throw std::runtime_error("spatial pipeline finalized twice");
        impl_->finalized = true;
        spatial_r1_decoder::ExactH0Counts h0;
        h0.reset(kGridColumns, kGridRows);
        for (const Impl::ThreadState &thread : impl_->threads) {
            h0.add(thread.h0);
            impl_->summary.readsDecoded += thread.decoded;
            impl_->summary.readsWithCandidates += thread.candidateReads;
            impl_->summary.uniqueGeneReads += thread.uniqueGeneReads;
            impl_->summary.barcodeReadsWithN += thread.barcodeReadsWithN;
            impl_->summary.barcodeNBases += thread.barcodeNBases;
            impl_->summary.barcodeDpRecoveredReads += thread.barcodeDpRecoveredReads;
            impl_->summary.barcodeDpAmbiguousReads += thread.barcodeDpAmbiguousReads;
            impl_->summary.barcodeDpUnassignedReads += thread.barcodeDpUnassignedReads;
            impl_->summary.barcodeUnsupportedReads += thread.barcodeUnsupportedReads;
            impl_->summary.umiReadsWithN += thread.umiReadsWithN;
            impl_->summary.umiReadsWithInvalidBase += thread.umiReadsWithInvalidBase;
        }
        impl_->summary.joinedReads = impl_->totalJoinedReads.load(std::memory_order_relaxed);
        impl_->summary.candidateRows = impl_->totalCandidateRows.load(std::memory_order_relaxed);
        impl_->summary.peakResidentReads =
            impl_->peakResidentReads.load(std::memory_order_relaxed);
        impl_->summary.peakResidentCandidates =
            impl_->peakResidentCandidates.load(std::memory_order_relaxed);
        impl_->summary.exactH0Reads = h0.reads;
        if (impl_->summary.readsDecoded > impl_->config.expectedReads
            || impl_->summary.joinedReads > impl_->config.expectedReads
            || impl_->summary.candidateRows > impl_->config.expectedCandidates) {
            throw std::runtime_error("integrated spatial evidence exceeded declared capacity");
        }
        std::vector<Clique> cliques;
        std::vector<CliqueCandidate> cliqueCandidates;
        std::vector<std::uint64_t> coordinateBits(
            (static_cast<std::size_t>(kGridRows) * kGridColumns + 63) / 64, 0);
        std::vector<std::uint8_t> featureUsed(geneIds.size(), 0);
        const std::chrono::steady_clock::time_point mergeBegin =
            std::chrono::steady_clock::now();
        if (impl_->config.overflowPolicy == OverflowPolicy::Spill) {
            buildCliquesSpill(*impl_, h0, cliques, cliqueCandidates,
                              coordinateBits, featureUsed);
            impl_->summary.spillMergeSeconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - mergeBegin).count();
            impl_->summary.spillRuns = 0;
            impl_->summary.spillBytes = 0;
            for (const Impl::ThreadState &thread : impl_->threads) {
                for (const spill::Run &run : thread.spillRuns) {
                    if (impl_->summary.spillBytes
                        > std::numeric_limits<std::uint64_t>::max() - run.bytes) {
                        throw std::runtime_error("spatial spill summary arithmetic overflow");
                    }
                    ++impl_->summary.spillRuns;
                    impl_->summary.spillBytes += run.bytes;
                }
            }
        } else {
            buildCliquesMemory(*impl_, h0, cliques, cliqueCandidates,
                               coordinateBits, featureUsed);
        }
        for (Impl::ThreadState &thread : impl_->threads) {
            std::vector<ReadEvidence>().swap(thread.reads);
            std::vector<CandidateEvidence>().swap(thread.candidates);
            thread.h0.bc1.clear();
            thread.h0.bc2.clear();
        }

        std::vector<FinalMolecule> molecules;
        if ((impl_->config.products & ProductStrict) != 0) {
            std::vector<FinalMolecule> values = resolveInteger(
                cliques, cliqueCandidates, ProductStrict);
            impl_->summary.strictMolecules = values.size();
            molecules.insert(molecules.end(), values.begin(), values.end());
        }
        if ((impl_->config.products & ProductSoftExpected) != 0) {
            std::vector<FinalMolecule> values = resolveSoft(cliques, cliqueCandidates);
            for (const FinalMolecule &value : values) {
                impl_->summary.softExpectedMass += value.weight;
            }
            molecules.insert(molecules.end(), values.begin(), values.end());
        }
        if ((impl_->config.products & ProductHard) != 0) {
            std::vector<FinalMolecule> values = resolveInteger(
                cliques, cliqueCandidates, ProductHard);
            impl_->summary.hardMolecules = values.size();
            molecules.insert(molecules.end(), values.begin(), values.end());
        }
        if ((impl_->config.products & ProductGatedHard) != 0) {
            std::vector<FinalMolecule> values = resolveInteger(
                cliques, cliqueCandidates, ProductGatedHard);
            impl_->summary.gatedHardMolecules = values.size();
            molecules.insert(molecules.end(), values.begin(), values.end());
        }
        std::vector<Clique>().swap(cliques);
        std::vector<CliqueCandidate>().swap(cliqueCandidates);
        materialize(*impl_, geneIds, featureUsed, coordinateBits, molecules);

        std::ostringstream runSummary;
        runSummary << "schema\tstar_suite.spatial_gex_integrated.v1\n"
                   << "source_revision\t" << impl_->config.sourceRevision << '\n'
                   << "reads_decoded\t" << impl_->summary.readsDecoded << '\n'
                   << "reads_with_candidates\t" << impl_->summary.readsWithCandidates << '\n'
                   << "unique_gene_reads\t" << impl_->summary.uniqueGeneReads << '\n'
                   << "joined_reads\t" << impl_->summary.joinedReads << '\n'
                   << "candidate_rows\t" << impl_->summary.candidateRows << '\n'
                   << "exact_h0_reads\t" << impl_->summary.exactH0Reads << '\n'
                   << "barcode_reads_with_n\t" << impl_->summary.barcodeReadsWithN << '\n'
                   << "barcode_n_bases\t" << impl_->summary.barcodeNBases << '\n'
                   << "barcode_dp_recovered_reads\t"
                   << impl_->summary.barcodeDpRecoveredReads << '\n'
                   << "barcode_dp_ambiguous_reads\t"
                   << impl_->summary.barcodeDpAmbiguousReads << '\n'
                   << "barcode_dp_unassigned_reads\t"
                   << impl_->summary.barcodeDpUnassignedReads << '\n'
                   << "barcode_unsupported_reads\t"
                   << impl_->summary.barcodeUnsupportedReads << '\n'
                   << "umi_reads_with_n\t" << impl_->summary.umiReadsWithN << '\n'
                   << "umi_reads_with_invalid_base\t"
                   << impl_->summary.umiReadsWithInvalidBase << '\n'
                   << "read_cliques\t" << impl_->summary.readCliques << '\n'
                   << "strict_molecules\t" << impl_->summary.strictMolecules << '\n'
                   << "soft_expected_mass\t" << std::setprecision(17)
                   << impl_->summary.softExpectedMass << '\n'
                   << "hard_molecules\t" << impl_->summary.hardMolecules << '\n'
                   << "gated_hard_molecules\t" << impl_->summary.gatedHardMolecules << '\n'
                   << "estimated_peak_bytes\t" << impl_->memory.peakBytes << '\n'
                   << "memory_budget_bytes\t" << impl_->budgetBytes << '\n'
                   << "overflow_policy\t"
                   << (impl_->config.overflowPolicy == OverflowPolicy::Spill ? "Spill" : "Fail")
                   << '\n'
                   << "spill_runs\t" << impl_->summary.spillRuns << '\n'
                   << "spill_bytes\t" << impl_->summary.spillBytes << '\n'
                   << "spill_high_water_candidates_per_thread\t"
                   << impl_->summary.spillHighWaterCandidatesPerThread << '\n'
                   << "spill_merge_seconds\t" << std::setprecision(17)
                   << impl_->summary.spillMergeSeconds << '\n'
                   << "peak_resident_reads\t" << impl_->summary.peakResidentReads << '\n'
                   << "peak_resident_candidates\t"
                   << impl_->summary.peakResidentCandidates << '\n';
        {
            AtomicOutput output(impl_->config.outputDirectory + "/run_summary.tsv");
            output.stream() << runSummary.str();
            output.commit();
        }
        if (impl_->config.overflowPolicy == OverflowPolicy::Spill) {
            cleanupSpills(*impl_);
        }
        {
            AtomicOutput output(impl_->config.outputDirectory + "/RUN_COMPLETE");
            output.stream() << "star_suite.spatial_gex_integrated.v1\n";
            output.commit();
        }
        return true;
    } catch (const std::exception &exception) {
        error = exception.what();
        return false;
    }
}

const PipelineSummary &Pipeline::summary() const { return impl_->summary; }
const MemoryModel &Pipeline::memoryModel() const { return impl_->memory; }
std::uint64_t Pipeline::memoryBudgetBytes() const { return impl_->budgetBytes; }

} // namespace spatial_gex
