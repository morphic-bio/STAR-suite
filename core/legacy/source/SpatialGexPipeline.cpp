#include "SpatialGex.h"
#include "SpatialGexSpill.h"
#include "SpatialGexDownstreamSpool.h"
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
#include <sys/statvfs.h>
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
const std::uint8_t kProductBits[4] = {
    ProductStrict, ProductSoftExpected, ProductHard, ProductGatedHard
};
const char *const kProductNames[4] = {
    "strict", "soft_expected", "hard", "gated_hard"
};
const std::uint8_t kScaleBits[3] = {Scale2um, Scale8um, Scale16um};
// The downstream memory model reserves 8 GiB. This deliberately conservative
// allowance covers the resident contribution vector plus the largest
// correction/reconciliation working set derived from it.
const std::uint64_t kDownstreamWorkingBytesPerContribution = 1024;
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

class CliqueSink {
  public:
    virtual ~CliqueSink() {}
    virtual void reserve(std::uint64_t cliques, std::uint64_t candidates) = 0;
    virtual void append(std::uint32_t gene, std::uint32_t rawUmi,
                        std::uint32_t memberCount,
                        const std::vector<CliqueCandidate> &candidates) = 0;
    virtual std::uint64_t cliqueCount() const = 0;
};

class InMemoryCliqueSink : public CliqueSink {
  public:
    InMemoryCliqueSink(std::vector<Clique> &cliques,
                       std::vector<CliqueCandidate> &candidates)
        : cliques_(cliques), candidates_(candidates) {}

    void reserve(std::uint64_t cliques, std::uint64_t candidates) override
    {
        cliques_.reserve(static_cast<std::size_t>(cliques));
        candidates_.reserve(static_cast<std::size_t>(candidates));
    }

    void append(std::uint32_t gene, std::uint32_t rawUmi,
                std::uint32_t memberCount,
                const std::vector<CliqueCandidate> &values) override
    {
        if (values.empty() || values.size() > std::numeric_limits<std::uint16_t>::max()
            || candidates_.size() > std::numeric_limits<std::uint32_t>::max()
            || memberCount == 0) {
            throw std::runtime_error("invalid in-memory spatial clique");
        }
        Clique clique;
        clique.gene = gene;
        clique.rawUmi = rawUmi;
        clique.candidateBegin = static_cast<std::uint32_t>(candidates_.size());
        clique.candidateCount = static_cast<std::uint16_t>(values.size());
        clique.memberCount = memberCount;
        candidates_.insert(candidates_.end(), values.begin(), values.end());
        cliques_.push_back(clique);
    }

    std::uint64_t cliqueCount() const override { return cliques_.size(); }

  private:
    std::vector<Clique> &cliques_;
    std::vector<CliqueCandidate> &candidates_;
};

bool assignmentLess(const Assignment &left, const Assignment &right)
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

std::vector<MatrixEntry> aggregateScaledMatrixEntries(
    std::vector<MatrixEntry> values)
{
    // matrixAtScale receives fine entries in canonical (fine column, feature)
    // order. Preserve that order within equal coarse keys so floating-point
    // addition is independent of whether the complete matrix or one
    // parent-coordinate shard is being materialized.
    std::stable_sort(values.begin(), values.end(), [](const MatrixEntry &left,
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
        std::uint64_t currentSourceOrdinal = 0;
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
        std::uint64_t featureAssignedReads = 0;
        std::uint64_t featureUnassignedReads = 0;
        std::uint64_t featureAssignedReadsWithCandidates = 0;
        std::uint64_t featureAssignedReadsWithoutCandidates = 0;
        std::uint64_t featureUnassignedReadsWithCandidates = 0;
        std::uint64_t featureUnassignedReadsWithoutCandidates = 0;
        std::uint64_t flexHashH0Reads = 0;
        std::uint64_t flexHashH1Reads = 0;
        std::uint64_t flexHashDenyReads = 0;
        std::uint64_t flexAlignmentMissReads = 0;
        std::uint64_t flexAlignmentResolvedReads = 0;
        std::uint64_t flexAlignmentUnresolvedReads = 0;
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
    const spatial_r1_decoder::ExactH0Counts &h0, CliqueSink &sink,
    std::vector<std::uint64_t> &coordinateBits,
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
        std::vector<CliqueCandidate> cliqueCandidates;
        cliqueCandidates.reserve(partition.candidates.size());
        for (std::size_t index = 0; index < partition.candidates.size(); ++index) {
            CliqueCandidate candidate;
            candidate.coordinate = partition.candidates[index];
            candidate.posterior = evidence[index] / total;
            cliqueCandidates.push_back(candidate);
            coordinateBits[candidate.coordinate / 64] |=
                1ULL << (candidate.coordinate % 64);
        }
        featureUsed[reads[begin].gene] = 1;
        sink.append(reads[begin].gene, reads[begin].rawUmi,
                    static_cast<std::uint32_t>(partition.members.size()),
                    cliqueCandidates);
    }
}

void appendCompleteGroup(const std::vector<Bundle> &group,
                         const spatial_r1_decoder::ExactH0Counts &h0,
                         CliqueSink &sink,
                         std::vector<std::uint64_t> &coordinateBits,
                         std::vector<std::uint8_t> &featureUsed)
{
    if (group.empty()) return;
    if (group.front().gene >= featureUsed.size()) {
        throw std::runtime_error("spatial GeneFull index lies outside transcriptome axis");
    }
    appendCliquesForGroup(group, 0, group.size(), h0, sink,
                          coordinateBits, featureUsed);
}

void buildCliquesMemory(Pipeline::Impl &impl,
                        const spatial_r1_decoder::ExactH0Counts &h0,
                        CliqueSink &sink,
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

    sink.reserve(order.size(), impl.summary.candidateRows);
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
        appendCompleteGroup(group, h0, sink,
                            coordinateBits, featureUsed);
        begin = end;
    }
    impl.summary.readCliques = sink.cliqueCount();
}

void buildCliquesSpill(Pipeline::Impl &impl,
                       const spatial_r1_decoder::ExactH0Counts &h0,
                       CliqueSink &sink,
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
    sink.reserve(impl.summary.joinedReads, impl.summary.candidateRows);
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
        appendCompleteGroup(group, h0, sink,
                            coordinateBits, featureUsed);
    }
    impl.summary.readCliques = sink.cliqueCount();
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

class SpoolCliqueSink : public CliqueSink {
  public:
    SpoolCliqueSink(const std::string &directory,
                    const std::string &sourceRevision,
                    std::uint32_t shards, std::size_t bufferBytes,
                    std::uint8_t products)
        : products_(products)
    {
        std::string error;
        writer_ = downstream_spool::ContributionWriter::create(
            directory, sourceRevision, kGridColumns, shards, bufferBytes, error);
        if (!writer_) throw std::runtime_error(error);
    }

    void reserve(std::uint64_t, std::uint64_t) override {}

    void append(std::uint32_t gene, std::uint32_t rawUmi,
                std::uint32_t memberCount,
                const std::vector<CliqueCandidate> &values) override
    {
        if (values.empty() || values.size() > std::numeric_limits<std::uint16_t>::max()
            || cliques_ > std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("invalid spooled spatial clique");
        }
        std::size_t best = 0;
        for (std::size_t index = 1; index < values.size(); ++index) {
            if ((!nearlyEqual(values[index].posterior, values[best].posterior)
                 && values[index].posterior > values[best].posterior)
                || (nearlyEqual(values[index].posterior, values[best].posterior)
                    && coordinateLess(values[index].coordinate,
                                      values[best].coordinate))) {
                best = index;
            }
        }
        double second = 0.0;
        for (std::size_t index = 0; index < values.size(); ++index) {
            if (index != best && values[index].posterior > second) {
                second = values[index].posterior;
            }
        }
        const double margin = nearlyEqual(values[best].posterior, second)
            ? 0.0 : values[best].posterior - second;
        const bool gated = values[best].posterior >= kGateMinimumPosterior
            && margin >= kGateMinimumMargin;

        const bool retainAllCandidates =
            (products_ & ProductSoftExpected) != 0;
        for (std::size_t index = 0; index < values.size(); ++index) {
            downstream_spool::Contribution record = {};
            record.posterior = values[index].posterior;
            record.gene = gene;
            record.coordinate = values[index].coordinate;
            record.rawUmi = rawUmi;
            record.memberCount = memberCount;
            record.cliqueOrdinal = static_cast<std::uint32_t>(cliques_);
            record.candidateOrdinal = static_cast<std::uint16_t>(index);
            if (values.size() == 1) {
                record.flags |= downstream_spool::ContributionStrict;
            }
            if (index == best) {
                record.flags |= downstream_spool::ContributionHard;
                if (gated) record.flags |= downstream_spool::ContributionGatedHard;
            }
            // Soft expected counts require the complete posterior candidate
            // field. Integer policies require only records selected by their
            // policy flags. Keeping this product-aware avoids carrying losing
            // candidates through the production hard-only spool while leaving
            // the explicit diagnostic `all` path byte-for-byte unchanged.
            const bool selectedInteger =
                (((products_ & ProductStrict) != 0)
                     && ((record.flags & downstream_spool::ContributionStrict) != 0))
                || (((products_ & ProductHard) != 0)
                     && ((record.flags & downstream_spool::ContributionHard) != 0))
                || (((products_ & ProductGatedHard) != 0)
                     && ((record.flags
                          & downstream_spool::ContributionGatedHard) != 0));
            if (!retainAllCandidates && !selectedInteger) continue;
            std::string error;
            if (!writer_->append(record, error)) throw std::runtime_error(error);
        }
        ++cliques_;
    }

    std::uint64_t cliqueCount() const override { return cliques_; }

    void finish(std::vector<downstream_spool::Run> &runs,
                std::uint64_t &records, std::uint64_t &bytes)
    {
        std::string error;
        if (!writer_->finish(runs, error)) throw std::runtime_error(error);
        records = writer_->records();
        bytes = writer_->bytes();
    }

  private:
    std::unique_ptr<downstream_spool::ContributionWriter> writer_;
    std::uint8_t products_ = 0;
    std::uint64_t cliques_ = 0;
};

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

std::vector<FinalMolecule> resolveSortedIntegerAssignments(
    const std::vector<Assignment> &assignments, std::uint32_t policy)
{
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

std::vector<FinalMolecule> resolveInteger(
    const std::vector<Clique> &cliques,
    const std::vector<CliqueCandidate> &candidates, std::uint32_t policy)
{
    const std::vector<Assignment> assignments =
        integerAssignments(cliques, candidates, policy);
    return resolveSortedIntegerAssignments(assignments, policy);
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

std::vector<FinalMolecule> resolveSoftContributions(
    const std::vector<downstream_spool::Contribution> &contributions)
{
    std::vector<SoftRawSupport> raw;
    raw.reserve(contributions.size());
    for (const downstream_spool::Contribution &contribution : contributions) {
        SoftRawSupport value;
        value.gene = contribution.gene;
        value.coordinate = contribution.coordinate;
        value.rawUmi = contribution.rawUmi;
        value.correctedUmi = contribution.rawUmi;
        value.support = contribution.posterior;
        raw.push_back(value);
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
        order.reserve(end - begin);
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

    std::map<std::tuple<std::uint32_t, std::uint32_t, std::uint32_t>,
             SoftProvisional> values;
    for (const downstream_spool::Contribution &contribution : contributions) {
        const SoftRawSupport &support = findSoftSupport(
            aggregated, contribution.gene, contribution.coordinate,
            contribution.rawUmi);
        const std::tuple<std::uint32_t, std::uint32_t, std::uint32_t> key(
            contribution.coordinate, support.correctedUmi, contribution.gene);
        SoftProvisional &value = values[key];
        value.coordinate = contribution.coordinate;
        value.correctedUmi = support.correctedUmi;
        value.gene = contribution.gene;
        value.absent *= 1.0 - contribution.posterior;
        if (contribution.rawUmi == support.correctedUmi) {
            value.originalAbsent *= 1.0 - contribution.posterior;
        }
    }
    std::vector<SoftProvisional> provisional;
    provisional.reserve(values.size());
    for (const std::pair<const std::tuple<std::uint32_t, std::uint32_t, std::uint32_t>,
                         SoftProvisional> &entry : values) {
        provisional.push_back(entry.second);
    }
    std::vector<FinalMolecule> result;
    begin = 0;
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

std::vector<FinalMolecule> resolveContributionShard(
    const std::vector<downstream_spool::Contribution> &contributions,
    std::uint8_t products)
{
    std::vector<FinalMolecule> result;
    struct IntegerPolicy {
        std::uint8_t product;
        std::uint8_t flag;
    };
    const IntegerPolicy policies[] = {
        {ProductStrict, downstream_spool::ContributionStrict},
        {ProductHard, downstream_spool::ContributionHard},
        {ProductGatedHard, downstream_spool::ContributionGatedHard}
    };
    for (const IntegerPolicy &policy : policies) {
        if ((products & policy.product) == 0) continue;
        std::vector<Assignment> assignments;
        assignments.reserve(contributions.size());
        for (const downstream_spool::Contribution &contribution : contributions) {
            if ((contribution.flags & policy.flag) == 0) continue;
            Assignment assignment;
            assignment.gene = contribution.gene;
            assignment.coordinate = contribution.coordinate;
            assignment.rawUmi = contribution.rawUmi;
            assignment.count = contribution.memberCount;
            assignments.push_back(assignment);
        }
        std::sort(assignments.begin(), assignments.end(), assignmentLess);
        std::vector<FinalMolecule> values =
            resolveSortedIntegerAssignments(assignments, policy.product);
        result.insert(result.end(), values.begin(), values.end());
    }
    if ((products & ProductSoftExpected) != 0) {
        std::vector<FinalMolecule> values = resolveSoftContributions(contributions);
        result.insert(result.end(), values.begin(), values.end());
    }
    std::sort(result.begin(), result.end(), [](const FinalMolecule &left,
                                               const FinalMolecule &right) {
        return std::tie(left.policy, left.coordinateIndex, left.correctedUmi,
                        left.geneIndex)
            < std::tie(right.policy, right.coordinateIndex, right.correctedUmi,
                       right.geneIndex);
    });
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
    return aggregateScaledMatrixEntries(std::move(coarse));
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

class SpoolMaterializer {
  public:
    SpoolMaterializer(Pipeline::Impl &impl,
                      const std::vector<std::string> &geneIds,
                      const std::vector<std::uint8_t> &featureUsed,
                      const std::vector<std::uint64_t> &coordinateBits,
                      const std::string &spoolRoot)
        : impl_(impl), spoolRoot_(spoolRoot),
          outputRoot_(impl.config.outputDirectory + ".partial")
    {
        std::vector<std::pair<std::string, std::uint32_t> > featurePairs;
        for (std::size_t gene = 0; gene < featureUsed.size(); ++gene) {
            if (featureUsed[gene]) {
                featurePairs.push_back(std::make_pair(
                    geneIds[gene], static_cast<std::uint32_t>(gene)));
            }
        }
        std::sort(featurePairs.begin(), featurePairs.end());
        featureIndex_.assign(geneIds.size(), 0);
        features_.reserve(featurePairs.size());
        for (std::size_t index = 0; index < featurePairs.size(); ++index) {
            if (index > 0 && featurePairs[index - 1].first == featurePairs[index].first) {
                throw std::runtime_error("duplicate canonical gene ID in spatial feature axis");
            }
            features_.push_back(featurePairs[index].first);
            featureIndex_[featurePairs[index].second] =
                static_cast<std::uint32_t>(index + 1);
        }

        std::vector<std::uint32_t> fineCoordinates;
        for (std::size_t word = 0; word < coordinateBits.size(); ++word) {
            std::uint64_t bits = coordinateBits[word];
            while (bits != 0) {
                const unsigned offset = static_cast<unsigned>(__builtin_ctzll(bits));
                fineCoordinates.push_back(
                    static_cast<std::uint32_t>(word * 64 + offset));
                bits &= bits - 1;
            }
        }
        std::sort(fineCoordinates.begin(), fineCoordinates.end(), coordinateLess);
        if (features_.empty() || fineCoordinates.empty()) {
            throw std::runtime_error("integrated spatial read-clique universe is empty");
        }
        axes_[0] = buildAxis(fineCoordinates, 1, "square_002um");
        axes_[1] = buildAxis(fineCoordinates, 4, "square_008um");
        axes_[2] = buildAxis(fineCoordinates, 8, "square_016um");
        fineIndex_.assign(static_cast<std::size_t>(kGridRows) * kGridColumns, 0);
        for (std::size_t index = 0; index < axes_[0].coordinates.size(); ++index) {
            fineIndex_[axes_[0].coordinates[index]] =
                static_cast<std::uint32_t>(index + 1);
        }
        makeDirectories(spoolRoot_);
        for (std::size_t product = 0; product < 4; ++product) {
            if ((impl_.config.products & kProductBits[product]) == 0) continue;
            for (std::size_t scale = 0; scale < 3; ++scale) {
                if ((impl_.config.scales & kScaleBits[scale]) == 0) continue;
                makeDirectories(runDirectory(product, scale));
            }
        }
        prepareOutputAxes();
    }

    void writeShard(std::uint32_t shard,
                    const std::vector<FinalMolecule> &molecules)
    {
        for (std::size_t product = 0; product < 4; ++product) {
            if ((impl_.config.products & kProductBits[product]) == 0) continue;
            std::vector<MatrixEntry> fineEntries;
            for (const FinalMolecule &molecule : molecules) {
                if (molecule.policy != kProductBits[product]) continue;
                if (molecule.geneIndex >= featureIndex_.size()
                    || molecule.coordinateIndex >= fineIndex_.size()
                    || featureIndex_[molecule.geneIndex] == 0
                    || fineIndex_[molecule.coordinateIndex] == 0) {
                    throw std::runtime_error(
                        "spooled final molecule lies outside spatial matrix axes");
                }
                MatrixEntry entry;
                entry.key = matrixKey(fineIndex_[molecule.coordinateIndex],
                                      featureIndex_[molecule.geneIndex]);
                entry.value = molecule.weight;
                fineEntries.push_back(entry);
            }
            fineEntries = aggregateMatrixEntries(std::move(fineEntries));
            for (std::size_t scale = 0; scale < 3; ++scale) {
                if ((impl_.config.scales & kScaleBits[scale]) == 0) continue;
                const std::vector<MatrixEntry> scaled = scale == 0
                    ? fineEntries : matrixAtScale(fineEntries, axes_[0], axes_[scale]);
                if (scaled.empty()) continue;
                std::vector<downstream_spool::MatrixRecord> records;
                records.reserve(scaled.size());
                for (const MatrixEntry &entry : scaled) {
                    downstream_spool::MatrixRecord record;
                    record.key = entry.key;
                    record.value = entry.value;
                    records.push_back(record);
                }
                downstream_spool::Run run;
                std::string error;
                if (!downstream_spool::writeMatrixRun(
                        runDirectory(product, scale), impl_.config.sourceRevision,
                        shard, impl_.config.downstreamSpoolShards, 0,
                        records, run, error)) {
                    throw std::runtime_error(error);
                }
                runs_[product][scale].push_back(run);
                ++impl_.summary.downstreamMatrixRuns;
                impl_.summary.downstreamMatrixBytes += run.bytes;
            }
        }
    }

    void commit()
    {
        std::ostringstream summary;
        summary << "schema\tstar_suite.molecule_first.policy_mex.v1\n"
                << "star_suite_version\t" << impl_.config.starSuiteVersion << '\n'
                << "assay\tvisium-hd\n"
                << "umi_mode\t1mm_cr\n"
                // The canonical feature/coordinate axes still come from the
                // compact in-memory used-bitsets populated by read cliques.
                // Keep the accepted summary contract independent of whether
                // molecule resolution itself is streamed.
                << "axes_source\tin_memory_read_cliques\n"
                << "cell_calling_order\tafter_postcollapse_integer_matrix_only\n"
                << "soft_cell_calling_allowed\tfalse\n"
                << "product\tscale\tfeatures\tbarcodes\tnnz\tmass\tmatrix_field\n";
        for (std::size_t product = 0; product < 4; ++product) {
            if ((impl_.config.products & kProductBits[product]) == 0) continue;
            double referenceMass = -1.0;
            for (std::size_t scale = 0; scale < 3; ++scale) {
                if ((impl_.config.scales & kScaleBits[scale]) == 0) continue;
                const MatrixResult result = mergeMatrix(
                    product, scale, kProductBits[product] == ProductSoftExpected);
                if (referenceMass < 0.0) referenceMass = result.mass;
                const double tolerance = 1e-9
                    * std::max(1.0, std::max(referenceMass, result.mass));
                if (std::fabs(referenceMass - result.mass) > tolerance) {
                    throw std::runtime_error(
                        "spooled spatial scale aggregation did not conserve mass");
                }
                summary << kProductNames[product] << '\t' << axes_[scale].name << '\t'
                        << features_.size() << '\t' << axes_[scale].coordinates.size()
                        << '\t' << result.nnz << '\t' << std::setprecision(17)
                        << result.mass << '\t'
                        << (kProductBits[product] == ProductSoftExpected
                                ? "real" : "integer") << '\n';
            }
        }
        AtomicOutput output(outputRoot_ + "/summary.tsv");
        output.stream() << summary.str();
        output.commit();
        const std::string &finalRoot = impl_.config.outputDirectory;
        struct stat finalInfo;
        if (::stat(finalRoot.c_str(), &finalInfo) == 0) {
            if (!S_ISDIR(finalInfo.st_mode) || !directoryEmpty(finalRoot)
                || ::rmdir(finalRoot.c_str()) != 0) {
                throw std::runtime_error(
                    "spatial final output directory is not empty: " + finalRoot);
            }
        } else if (errno != ENOENT) {
            throw std::runtime_error("cannot inspect spatial final output directory: "
                                     + finalRoot);
        }
        if (std::rename(outputRoot_.c_str(), finalRoot.c_str()) != 0) {
            throw std::runtime_error("cannot commit spatial output directory "
                                     + finalRoot + ": " + std::strerror(errno));
        }
        cleanupRunDirectories();
    }

  private:
    struct MatrixResult {
        std::uint64_t nnz = 0;
        double mass = 0.0;
    };

    std::string runDirectory(std::size_t product, std::size_t scale) const
    {
        std::ostringstream path;
        path << spoolRoot_ << "/matrix.p" << product << ".s" << scale;
        return path.str();
    }

    void prepareOutputAxes()
    {
        const std::string &finalRoot = impl_.config.outputDirectory;
        struct stat info;
        if (::stat(finalRoot.c_str(), &info) == 0) {
            if (!S_ISDIR(info.st_mode) || !directoryEmpty(finalRoot)) {
                throw std::runtime_error(
                    "spatial output directory must be empty: " + finalRoot);
            }
        } else if (errno != ENOENT) {
            throw std::runtime_error("cannot inspect spatial output directory: "
                                     + finalRoot);
        }
        if (::stat(outputRoot_.c_str(), &info) == 0) {
            throw std::runtime_error("spatial partial output already exists: "
                                     + outputRoot_);
        } else if (errno != ENOENT) {
            throw std::runtime_error("cannot inspect spatial partial output: "
                                     + outputRoot_);
        }
        makeDirectories(outputRoot_);
        const std::string &root = outputRoot_;
        std::vector<std::size_t> selectedProducts;
        for (std::size_t product = 0; product < 4; ++product) {
            if ((impl_.config.products & kProductBits[product]) != 0) {
                selectedProducts.push_back(product);
                for (std::size_t scale = 0; scale < 3; ++scale) {
                    if ((impl_.config.scales & kScaleBits[scale]) != 0) {
                        makeDirectories(root + "/" + kProductNames[product]
                                        + "/" + axes_[scale].name);
                    }
                }
            }
        }
        if (selectedProducts.empty()) throw std::runtime_error("no spatial products selected");
        std::size_t canonicalScale = 0;
        while (canonicalScale < 3
               && (impl_.config.scales & kScaleBits[canonicalScale]) == 0) {
            ++canonicalScale;
        }
        const std::string canonicalDirectory = root + "/"
            + kProductNames[selectedProducts.front()] + "/" + axes_[canonicalScale].name;
        const std::string canonicalFeatures = canonicalDirectory + "/features.tsv";
        writeFeatureAxis(canonicalFeatures, features_);
        for (std::size_t scale = 0; scale < 3; ++scale) {
            if ((impl_.config.scales & kScaleBits[scale]) == 0) continue;
            const std::string canonicalBarcodes = root + "/"
                + kProductNames[selectedProducts.front()] + "/" + axes_[scale].name
                + "/barcodes.tsv";
            writeBarcodeAxis(canonicalBarcodes, axes_[scale]);
            for (std::size_t product : selectedProducts) {
                const std::string directory = root + "/" + kProductNames[product]
                    + "/" + axes_[scale].name;
                if (directory + "/features.tsv" != canonicalFeatures) {
                    atomicLinkOrCopy(canonicalFeatures, directory + "/features.tsv");
                }
                if (directory + "/barcodes.tsv" != canonicalBarcodes) {
                    atomicLinkOrCopy(canonicalBarcodes, directory + "/barcodes.tsv");
                }
            }
        }
    }

    MatrixResult mergeMatrix(std::size_t product, std::size_t scale, bool real)
    {
        struct CursorState {
            std::unique_ptr<downstream_spool::MatrixCursor> cursor;
            downstream_spool::MatrixRecord record = {};
            bool active = false;
        };
        std::vector<CursorState> cursors;
        cursors.reserve(runs_[product][scale].size());
        for (const downstream_spool::Run &run : runs_[product][scale]) {
            std::string error;
            CursorState state;
            state.cursor = downstream_spool::MatrixCursor::open(
                run, impl_.config.sourceRevision, error);
            if (!state.cursor) throw std::runtime_error(error);
            state.active = state.cursor->next(state.record, error);
            if (!state.active && !error.empty()) throw std::runtime_error(error);
            cursors.push_back(std::move(state));
        }
        struct Later {
            const std::vector<CursorState> *cursors;
            bool operator()(std::size_t left, std::size_t right) const
            {
                const downstream_spool::MatrixRecord &a = (*cursors)[left].record;
                const downstream_spool::MatrixRecord &b = (*cursors)[right].record;
                if (a.key != b.key) return a.key > b.key;
                return left > right;
            }
        };
        Later later;
        later.cursors = &cursors;
        std::priority_queue<std::size_t, std::vector<std::size_t>, Later> queue(later);
        for (std::size_t index = 0; index < cursors.size(); ++index) {
            if (cursors[index].active) queue.push(index);
        }

        std::ostringstream bodyName;
        bodyName << spoolRoot_ << "/body.p" << product << ".s" << scale << ".tmp";
        const std::string bodyPath = bodyName.str();
        std::ofstream body(bodyPath.c_str(), std::ios::binary | std::ios::trunc);
        if (!body) throw std::runtime_error("cannot create spooled matrix body");
        body << std::setprecision(17);
        MatrixResult result;
        while (!queue.empty()) {
            const std::uint64_t key = cursors[queue.top()].record.key;
            double value = 0.0;
            do {
                const std::size_t index = queue.top();
                queue.pop();
                value += cursors[index].record.value;
                std::string error;
                cursors[index].active = cursors[index].cursor->next(
                    cursors[index].record, error);
                if (!cursors[index].active && !error.empty()) {
                    throw std::runtime_error(error);
                }
                if (cursors[index].active) queue.push(index);
            } while (!queue.empty() && cursors[queue.top()].record.key == key);
            body << keyFeature(key) << ' ' << keyColumn(key) << ' ';
            if (real) {
                body << value;
            } else {
                const long long rounded = std::llround(value);
                if (rounded < 0 || !nearlyEqual(value, static_cast<double>(rounded))) {
                    throw std::runtime_error(
                        "integer spooled spatial matrix contains a non-integer value");
                }
                body << static_cast<std::uint64_t>(rounded);
            }
            body << '\n';
            ++result.nnz;
            result.mass += value;
        }
        body.flush();
        if (!body) throw std::runtime_error("cannot write spooled matrix body");
        body.close();

        const std::string matrixPath = outputRoot_ + "/"
            + kProductNames[product] + "/" + axes_[scale].name + "/matrix.mtx";
        AtomicOutput output(matrixPath);
        output.stream() << "%%MatrixMarket matrix coordinate "
                        << (real ? "real" : "integer") << " general\n"
                        << "% STAR Suite molecule-first post-collapse policy matrix\n"
                        << features_.size() << ' ' << axes_[scale].coordinates.size()
                        << ' ' << result.nnz << '\n';
        std::ifstream bodyInput(bodyPath.c_str(), std::ios::binary);
        if (!bodyInput) throw std::runtime_error("cannot reopen spooled matrix body");
        output.stream() << bodyInput.rdbuf();
        if (!bodyInput.eof() && bodyInput.fail()) {
            throw std::runtime_error("cannot read spooled matrix body");
        }
        output.commit();
        if (std::remove(bodyPath.c_str()) != 0) {
            throw std::runtime_error("cannot remove committed spooled matrix body");
        }
        return result;
    }

    void cleanupRunDirectories()
    {
        for (std::size_t product = 0; product < 4; ++product) {
            if ((impl_.config.products & kProductBits[product]) == 0) continue;
            for (std::size_t scale = 0; scale < 3; ++scale) {
                if ((impl_.config.scales & kScaleBits[scale]) == 0) continue;
                for (const downstream_spool::Run &run : runs_[product][scale]) {
                    if (std::remove(run.path.c_str()) != 0 && errno != ENOENT) {
                        throw std::runtime_error(
                            "cannot remove committed matrix run " + run.path);
                    }
                }
                runs_[product][scale].clear();
                const std::string directory = runDirectory(product, scale);
                if (::rmdir(directory.c_str()) != 0 && errno != ENOENT) {
                    throw std::runtime_error("cannot remove matrix run directory "
                                             + directory);
                }
            }
        }
    }

    Pipeline::Impl &impl_;
    std::string spoolRoot_;
    std::string outputRoot_;
    std::vector<std::string> features_;
    std::vector<std::uint32_t> featureIndex_;
    std::vector<std::uint32_t> fineIndex_;
    Axis axes_[3];
    std::vector<downstream_spool::Run> runs_[4][3];
};

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
        if (config.expectedReads > std::numeric_limits<std::uint32_t>::max()
            || config.expectedCandidates > std::numeric_limits<std::uint32_t>::max()) {
            throw std::invalid_argument(
                "integrated spatial compact capacities must fit uint32");
        }
        if (config.barcodeContractDirectory.empty() || config.bc1OligosPath.empty()
            || config.bc2OligosPath.empty() || config.outputDirectory.empty()) {
            throw std::invalid_argument("integrated spatial GEX input/output paths are required");
        }
        if (config.overflowPolicy == OverflowPolicy::Spill
            && config.temporaryDirectory.empty()) {
            throw std::invalid_argument("integrated spatial spill requires a temporary directory");
        }
        if (config.downstreamSpoolShards == 0
            || config.downstreamSpoolShards > 4096
            || config.downstreamSpoolBufferBytes == 0
            || config.downstreamSpoolBufferBytes
                > std::numeric_limits<std::size_t>::max()) {
            throw std::invalid_argument(
                "invalid integrated spatial downstream spool configuration");
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
            struct statvfs disk;
            if (::statvfs(config.temporaryDirectory.c_str(), &disk) != 0
                || disk.f_frsize == 0
                || static_cast<std::uint64_t>(disk.f_bavail)
                    > std::numeric_limits<std::uint64_t>::max()
                        / static_cast<std::uint64_t>(disk.f_frsize)) {
                error = "cannot determine integrated spatial spill disk capacity: "
                    + config.temporaryDirectory;
                return std::unique_ptr<Pipeline>();
            }
            const std::uint64_t diskAvailable =
                static_cast<std::uint64_t>(disk.f_bavail)
                * static_cast<std::uint64_t>(disk.f_frsize);
            if (impl->memory.downstreamSpoolDiskBytes > diskAvailable) {
                std::ostringstream message;
                message << "integrated spatial downstream spool requires "
                        << impl->memory.downstreamSpoolDiskBytes
                        << " temporary/output bytes but filesystem has "
                        << diskAvailable << " available at "
                        << config.temporaryDirectory;
                error = message.str();
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
                                   std::uint64_t sourceOrdinal,
                                   std::string &error)
{
    std::uint32_t threadIndex = 0;
    if (!currentThreadIndex(threadIndex, error)) return false;
    Impl::ThreadState &thread = impl_->threads[threadIndex];
    if (thread.currentDecodedReady) {
        error = "integrated spatial pipeline received a new R1 before the prior "
                "feature decision was completed";
        return false;
    }
    if (sourceOrdinal > std::numeric_limits<std::uint32_t>::max()) {
        error = "integrated spatial source ordinal exceeds compact uint32 range";
        return false;
    }
    if (!decode(threadIndex, sequence, sequenceLength, quality, qualityLength,
                thread.currentDecoded, error)) return false;
    thread.currentSourceOrdinal = sourceOrdinal;
    thread.currentDecodedReady = true;
    return true;
}

bool Pipeline::completeCurrentThread(FeatureEvidenceClass source,
                                     bool assigned,
                                     std::uint32_t geneIndex,
                                     std::uint64_t sourceOrdinal,
                                     std::string &error)
{
    std::uint32_t threadIndex = 0;
    if (!currentThreadIndex(threadIndex, error)) return false;
    Impl::ThreadState &thread = impl_->threads[threadIndex];
    if (!thread.currentDecodedReady) {
        error = "integrated spatial pipeline lost current-thread decoded R1 state";
        return false;
    }
    if (thread.currentSourceOrdinal != sourceOrdinal) {
        error = "integrated spatial feature decision does not match the current "
                "raw-R1 source ordinal";
        return false;
    }
    if ((source == FeatureEvidenceClass::FlexH0
         || source == FeatureEvidenceClass::FlexH1) && !assigned) {
        error = "integrated spatial Flex cache hit was completed without a feature";
        return false;
    }
    if (source == FeatureEvidenceClass::FlexHashDeny && assigned) {
        error = "integrated spatial Flex hash deny was completed with a feature";
        return false;
    }
    if (assigned && impl_->config.featureCount != 0
        && geneIndex >= impl_->config.featureCount) {
        error = "integrated spatial feature index is outside the declared axis";
        return false;
    }

    if (assigned) {
        ++thread.featureAssignedReads;
        ++thread.uniqueGeneReads;
        if (thread.currentDecoded.candidates.empty()) {
            ++thread.featureAssignedReadsWithoutCandidates;
        } else {
            ++thread.featureAssignedReadsWithCandidates;
        }
    } else {
        ++thread.featureUnassignedReads;
        if (thread.currentDecoded.candidates.empty()) {
            ++thread.featureUnassignedReadsWithoutCandidates;
        } else {
            ++thread.featureUnassignedReadsWithCandidates;
        }
    }
    if (source == FeatureEvidenceClass::FlexH0) {
        ++thread.flexHashH0Reads;
    } else if (source == FeatureEvidenceClass::FlexH1) {
        ++thread.flexHashH1Reads;
    } else if (source == FeatureEvidenceClass::FlexHashDeny) {
        ++thread.flexHashDenyReads;
    } else if (source == FeatureEvidenceClass::FlexAlignment) {
        ++thread.flexAlignmentMissReads;
        if (assigned) {
            ++thread.flexAlignmentResolvedReads;
        } else {
            ++thread.flexAlignmentUnresolvedReads;
        }
    }

    bool completed = true;
    if (assigned) {
        completed = append(threadIndex, geneIndex, sourceOrdinal,
                           thread.currentDecoded, error);
    }
    thread.currentDecodedReady = false;
    return completed;
}

bool Pipeline::appendCurrentThread(std::uint32_t geneIndex,
                                   std::uint64_t sourceOrdinal,
                                   std::string &error)
{
    return completeCurrentThread(FeatureEvidenceClass::Gex, true, geneIndex,
                                 sourceOrdinal, error);
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
            if (thread.currentDecodedReady) {
                throw std::runtime_error(
                    "integrated spatial finalization found an R1 without a terminal "
                    "feature decision");
            }
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
            impl_->summary.featureAssignedReads += thread.featureAssignedReads;
            impl_->summary.featureUnassignedReads += thread.featureUnassignedReads;
            impl_->summary.featureAssignedReadsWithCandidates +=
                thread.featureAssignedReadsWithCandidates;
            impl_->summary.featureAssignedReadsWithoutCandidates +=
                thread.featureAssignedReadsWithoutCandidates;
            impl_->summary.featureUnassignedReadsWithCandidates +=
                thread.featureUnassignedReadsWithCandidates;
            impl_->summary.featureUnassignedReadsWithoutCandidates +=
                thread.featureUnassignedReadsWithoutCandidates;
            impl_->summary.flexHashH0Reads += thread.flexHashH0Reads;
            impl_->summary.flexHashH1Reads += thread.flexHashH1Reads;
            impl_->summary.flexHashDenyReads += thread.flexHashDenyReads;
            impl_->summary.flexAlignmentMissReads += thread.flexAlignmentMissReads;
            impl_->summary.flexAlignmentResolvedReads +=
                thread.flexAlignmentResolvedReads;
            impl_->summary.flexAlignmentUnresolvedReads +=
                thread.flexAlignmentUnresolvedReads;
        }
        if (impl_->config.requirePairedCompletion
            && impl_->summary.featureAssignedReads
                    + impl_->summary.featureUnassignedReads
                != impl_->summary.readsDecoded) {
            throw std::runtime_error(
                "integrated spatial feature decisions do not cover every decoded R1");
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
        std::vector<std::uint64_t> coordinateBits(
            (static_cast<std::size_t>(kGridRows) * kGridColumns + 63) / 64, 0);
        std::vector<std::uint8_t> featureUsed(geneIds.size(), 0);
        std::vector<Clique> cliques;
        std::vector<CliqueCandidate> cliqueCandidates;
        std::vector<downstream_spool::Run> contributionRuns;
        std::string downstreamRoot;
        std::string contributionDirectory;
        const std::chrono::steady_clock::time_point mergeBegin =
            std::chrono::steady_clock::now();
        if (impl_->config.overflowPolicy == OverflowPolicy::Spill) {
            downstreamRoot = impl_->config.temporaryDirectory
                + "/spatial_gex_downstream";
            contributionDirectory = downstreamRoot + "/contributions";
            struct stat downstreamInfo;
            if (::stat(downstreamRoot.c_str(), &downstreamInfo) == 0) {
                if (!S_ISDIR(downstreamInfo.st_mode)
                    || !directoryEmpty(downstreamRoot)) {
                    throw std::runtime_error(
                        "spatial downstream spool directory must be empty: "
                        + downstreamRoot);
                }
            } else if (errno != ENOENT) {
                throw std::runtime_error("cannot inspect spatial downstream spool: "
                                         + downstreamRoot);
            }
            makeDirectories(contributionDirectory);
            SpoolCliqueSink sink(contributionDirectory,
                                 impl_->config.sourceRevision,
                                 impl_->config.downstreamSpoolShards,
                                 static_cast<std::size_t>(
                                     impl_->config.downstreamSpoolBufferBytes),
                                 impl_->config.products);
            buildCliquesSpill(*impl_, h0, sink,
                              coordinateBits, featureUsed);
            sink.finish(contributionRuns,
                        impl_->summary.downstreamContributionRecords,
                        impl_->summary.downstreamContributionBytes);
            impl_->summary.downstreamContributionRuns = contributionRuns.size();
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
            InMemoryCliqueSink sink(cliques, cliqueCandidates);
            buildCliquesMemory(*impl_, h0, sink,
                               coordinateBits, featureUsed);
        }
        for (Impl::ThreadState &thread : impl_->threads) {
            std::vector<ReadEvidence>().swap(thread.reads);
            std::vector<CandidateEvidence>().swap(thread.candidates);
            thread.h0.bc1.clear();
            thread.h0.bc2.clear();
        }

        if (impl_->config.overflowPolicy == OverflowPolicy::Spill) {
            if (impl_->memory.downstreamSpoolBytes
                <= impl_->memory.accumulationFixedBytes) {
                throw std::runtime_error("invalid downstream spool workspace estimate");
            }
            const std::uint64_t workspace = impl_->memory.downstreamSpoolBytes
                - impl_->memory.accumulationFixedBytes;
            const std::uint64_t maximumShardRecords = workspace
                / kDownstreamWorkingBytesPerContribution;
            if (maximumShardRecords == 0) {
                throw std::runtime_error("downstream spool workspace is empty");
            }
            const std::chrono::steady_clock::time_point materializeBegin =
                std::chrono::steady_clock::now();
            SpoolMaterializer materializer(*impl_, geneIds, featureUsed,
                                            coordinateBits, downstreamRoot);
            double resolveSeconds = 0.0;
            for (const downstream_spool::Run &run : contributionRuns) {
                impl_->summary.downstreamLargestShardRecords = std::max(
                    impl_->summary.downstreamLargestShardRecords, run.records);
                if (run.records > maximumShardRecords
                    || run.records > std::numeric_limits<std::size_t>::max()) {
                    std::ostringstream message;
                    message << "spatial downstream shard " << run.shard << " has "
                            << run.records << " contribution records; bounded workspace "
                            << "permits " << maximumShardRecords
                            << ". Increase downstream shard count.";
                    throw std::runtime_error(message.str());
                }
                const std::chrono::steady_clock::time_point resolveBegin =
                    std::chrono::steady_clock::now();
                std::string cursorError;
                std::unique_ptr<downstream_spool::ContributionCursor> cursor =
                    downstream_spool::ContributionCursor::open(
                        run, impl_->config.sourceRevision, cursorError);
                if (!cursor) throw std::runtime_error(cursorError);
                std::vector<downstream_spool::Contribution> contributions;
                contributions.reserve(static_cast<std::size_t>(run.records));
                downstream_spool::Contribution contribution = {};
                while (cursor->next(contribution, cursorError)) {
                    if (contribution.coordinate
                            >= static_cast<std::uint64_t>(kGridRows) * kGridColumns
                        || downstream_spool::shardForCoordinate(
                               contribution.coordinate, kGridColumns,
                               impl_->config.downstreamSpoolShards) != run.shard) {
                        throw std::runtime_error(
                            "spatial downstream contribution has wrong coordinate shard");
                    }
                    contributions.push_back(contribution);
                }
                if (!cursorError.empty()) throw std::runtime_error(cursorError);
                if (contributions.size() != run.records) {
                    throw std::runtime_error(
                        "spatial downstream contribution count changed while reading");
                }
                std::vector<FinalMolecule> molecules = resolveContributionShard(
                    contributions, impl_->config.products);
                resolveSeconds += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - resolveBegin).count();
                for (const FinalMolecule &molecule : molecules) {
                    if (molecule.policy == ProductStrict) {
                        ++impl_->summary.strictMolecules;
                    } else if (molecule.policy == ProductSoftExpected) {
                        impl_->summary.softExpectedMass += molecule.weight;
                    } else if (molecule.policy == ProductHard) {
                        ++impl_->summary.hardMolecules;
                    } else if (molecule.policy == ProductGatedHard) {
                        ++impl_->summary.gatedHardMolecules;
                    } else {
                        throw std::runtime_error(
                            "spatial downstream resolver returned an unknown product");
                    }
                }
                materializer.writeShard(run.shard, molecules);
            }
            impl_->summary.downstreamResolveSeconds = resolveSeconds;
            materializer.commit();
            impl_->summary.downstreamMaterializeSeconds =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - materializeBegin).count()
                - resolveSeconds;
            for (const downstream_spool::Run &run : contributionRuns) {
                if (std::remove(run.path.c_str()) != 0 && errno != ENOENT) {
                    throw std::runtime_error(
                        "cannot remove committed contribution run " + run.path
                        + ": " + std::strerror(errno));
                }
            }
            if (::rmdir(contributionDirectory.c_str()) != 0) {
                throw std::runtime_error(
                    "cannot remove downstream contribution directory "
                    + contributionDirectory + ": " + std::strerror(errno));
            }
            if (::rmdir(downstreamRoot.c_str()) != 0) {
                throw std::runtime_error("cannot remove downstream spool directory "
                                         + downstreamRoot + ": "
                                         + std::strerror(errno));
            }
        } else {
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
        }

        const char *schema = impl_->config.flexFeatureMode
            ? "star_suite.spatial_flex_integrated.v1"
            : "star_suite.spatial_gex_integrated.v1";
        std::ostringstream runSummary;
        runSummary << "schema\t" << schema << '\n'
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
                   << "estimated_downstream_spool_bytes\t"
                   << impl_->memory.downstreamSpoolBytes << '\n'
                   << "estimated_downstream_spool_disk_bytes\t"
                   << impl_->memory.downstreamSpoolDiskBytes << '\n'
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
                   << "downstream_contribution_records\t"
                   << impl_->summary.downstreamContributionRecords << '\n'
                   << "downstream_contribution_runs\t"
                   << impl_->summary.downstreamContributionRuns << '\n'
                   << "downstream_contribution_bytes\t"
                   << impl_->summary.downstreamContributionBytes << '\n'
                   << "downstream_largest_shard_records\t"
                   << impl_->summary.downstreamLargestShardRecords << '\n'
                   << "downstream_matrix_runs\t"
                   << impl_->summary.downstreamMatrixRuns << '\n'
                   << "downstream_matrix_bytes\t"
                   << impl_->summary.downstreamMatrixBytes << '\n'
                   << "downstream_resolve_seconds\t" << std::setprecision(17)
                   << impl_->summary.downstreamResolveSeconds << '\n'
                   << "downstream_materialize_seconds\t" << std::setprecision(17)
                   << impl_->summary.downstreamMaterializeSeconds << '\n'
                   << "peak_resident_reads\t" << impl_->summary.peakResidentReads << '\n'
                   << "peak_resident_candidates\t"
                   << impl_->summary.peakResidentCandidates << '\n';
        if (impl_->config.flexFeatureMode) {
            runSummary
                << "feature_axis_path\t" << impl_->config.featureAxisPath << '\n'
                << "feature_axis_sha256\t"
                << impl_->config.featureAxisSha256 << '\n'
                << "feature_cache_path\t" << impl_->config.featureCachePath << '\n'
                << "feature_cache_sha256\t"
                << impl_->config.featureCacheSha256 << '\n'
                << "feature_count\t" << impl_->config.featureCount << '\n'
                << "feature_assigned_reads\t"
                << impl_->summary.featureAssignedReads << '\n'
                << "feature_unassigned_reads\t"
                << impl_->summary.featureUnassignedReads << '\n'
                << "feature_assigned_reads_with_spatial_candidates\t"
                << impl_->summary.featureAssignedReadsWithCandidates << '\n'
                << "feature_assigned_reads_without_spatial_candidates\t"
                << impl_->summary.featureAssignedReadsWithoutCandidates << '\n'
                << "feature_unassigned_reads_with_spatial_candidates\t"
                << impl_->summary.featureUnassignedReadsWithCandidates << '\n'
                << "feature_unassigned_reads_without_spatial_candidates\t"
                << impl_->summary.featureUnassignedReadsWithoutCandidates << '\n'
                << "feature_hash_h0\t" << impl_->summary.flexHashH0Reads << '\n'
                << "feature_hash_h1\t" << impl_->summary.flexHashH1Reads << '\n'
                << "feature_hash_deny\t"
                << impl_->summary.flexHashDenyReads << '\n'
                << "feature_hash_miss\t"
                << impl_->summary.flexAlignmentMissReads << '\n'
                << "feature_alignment_resolved\t"
                << impl_->summary.flexAlignmentResolvedReads << '\n'
                << "feature_alignment_unresolved\t"
                << impl_->summary.flexAlignmentUnresolvedReads << '\n';
        }
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
            output.stream() << schema << '\n';
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
