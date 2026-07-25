#include "SpatialGexSpill.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <tuple>

namespace spatial_gex {
namespace spill {
namespace {

const std::uint64_t kMagic = 0x31584f4650475353ULL; // "SSGPFOX1"
const std::uint32_t kSchema = 1;
const std::uint32_t kReadKeyGeneUmiOrdinal = 1;
const std::uint64_t kComplete = 0x4554454c504d4f43ULL; // "COMPLETE"
const std::uint64_t kFnvOffset = 1469598103934665603ULL;
const std::uint64_t kFnvPrime = 1099511628211ULL;

#pragma pack(push, 1)
struct Header {
    std::uint64_t magic;
    std::uint32_t schema;
    std::uint32_t headerBytes;
    std::uint32_t key;
    std::uint32_t readBytes;
    std::uint32_t candidateBytes;
    std::uint32_t threadIndex;
    std::uint64_t readCount;
    std::uint64_t candidateCount;
    std::uint64_t sourceRevisionHash;
    std::uint64_t payloadBytes;
    std::uint64_t reserved;
};

struct Trailer {
    std::uint64_t checksum;
    std::uint64_t payloadBytes;
    std::uint64_t complete;
};
#pragma pack(pop)

static_assert(sizeof(Header) == 72, "spatial spill header size changed");
static_assert(sizeof(Trailer) == 24, "spatial spill trailer size changed");

std::uint64_t fnvUpdate(std::uint64_t hash, const void *data, std::size_t bytes)
{
    const unsigned char *cursor = static_cast<const unsigned char *>(data);
    for (std::size_t index = 0; index < bytes; ++index) {
        hash ^= cursor[index];
        hash *= kFnvPrime;
    }
    return hash;
}

std::uint64_t sourceHash(const std::string &sourceRevision)
{
    return fnvUpdate(kFnvOffset, sourceRevision.data(), sourceRevision.size());
}

template <typename T>
void writePayload(std::ofstream &output, const T &value, std::uint64_t &checksum)
{
    output.write(reinterpret_cast<const char *>(&value), sizeof(value));
    if (!output) throw std::runtime_error("cannot write spatial spill payload");
    checksum = fnvUpdate(checksum, &value, sizeof(value));
}

template <typename T>
void readPayload(std::ifstream &input, T &value, std::uint64_t &checksum)
{
    input.read(reinterpret_cast<char *>(&value), sizeof(value));
    if (!input) throw std::runtime_error("truncated spatial spill payload");
    checksum = fnvUpdate(checksum, &value, sizeof(value));
}

bool checkedPayloadBytes(std::uint64_t reads, std::uint64_t candidates,
                         std::uint64_t &bytes)
{
    if (reads > std::numeric_limits<std::uint64_t>::max() / sizeof(ReadEvidence)
        || candidates > std::numeric_limits<std::uint64_t>::max()
            / sizeof(CandidateEvidence)) return false;
    const std::uint64_t readBytes = reads * sizeof(ReadEvidence);
    const std::uint64_t candidateBytes = candidates * sizeof(CandidateEvidence);
    if (readBytes > std::numeric_limits<std::uint64_t>::max() - candidateBytes) {
        return false;
    }
    bytes = readBytes + candidateBytes;
    return true;
}

} // namespace

bool writeRun(const std::string &directory, const std::string &sourceRevision,
              std::uint32_t threadIndex, std::uint64_t runIndex,
              const std::vector<ReadEvidence> &reads,
              const std::vector<CandidateEvidence> &candidates,
              Run &run, std::string &error)
{
    try {
        std::uint64_t referencedCandidates = 0;
        for (const ReadEvidence &read : reads) {
            if (static_cast<std::uint64_t>(read.candidateBegin) + read.candidateCount
                    > candidates.size()
                || referencedCandidates > std::numeric_limits<std::uint64_t>::max()
                    - read.candidateCount) {
                throw std::runtime_error("spatial spill read candidate span is corrupt");
            }
            referencedCandidates += read.candidateCount;
        }
        if (referencedCandidates != candidates.size()) {
            throw std::runtime_error("spatial spill candidate pool contains orphan records");
        }

        std::uint64_t payloadBytes = 0;
        if (!checkedPayloadBytes(reads.size(), candidates.size(), payloadBytes)) {
            throw std::runtime_error("spatial spill payload size overflow");
        }
        std::vector<std::uint32_t> order(reads.size());
        for (std::size_t index = 0; index < order.size(); ++index) {
            if (index > std::numeric_limits<std::uint32_t>::max()) {
                throw std::runtime_error("spatial spill segment exceeds uint32 ordering range");
            }
            order[index] = static_cast<std::uint32_t>(index);
        }
        std::sort(order.begin(), order.end(), [&](std::uint32_t left,
                                                   std::uint32_t right) {
            const ReadEvidence &a = reads[left];
            const ReadEvidence &b = reads[right];
            return std::tie(a.geneIndex, a.rawUmi, a.sourceOrdinal)
                < std::tie(b.geneIndex, b.rawUmi, b.sourceOrdinal);
        });

        std::ostringstream name;
        name << directory << "/read.t" << std::setw(4) << std::setfill('0')
             << threadIndex << ".r" << std::setw(8) << runIndex << ".bin";
        const std::string path = name.str();
        const std::string partial = path + ".partial";
        std::ofstream output(partial.c_str(), std::ios::binary | std::ios::trunc);
        if (!output) throw std::runtime_error("cannot create spatial spill " + partial);

        Header header = {};
        header.magic = kMagic;
        header.schema = kSchema;
        header.headerBytes = sizeof(Header);
        header.key = kReadKeyGeneUmiOrdinal;
        header.readBytes = sizeof(ReadEvidence);
        header.candidateBytes = sizeof(CandidateEvidence);
        header.threadIndex = threadIndex;
        header.readCount = reads.size();
        header.candidateCount = candidates.size();
        header.sourceRevisionHash = sourceHash(sourceRevision);
        header.payloadBytes = payloadBytes;
        output.write(reinterpret_cast<const char *>(&header), sizeof(header));
        if (!output) throw std::runtime_error("cannot write spatial spill header");

        std::uint64_t checksum = kFnvOffset;
        for (std::uint32_t index : order) {
            ReadEvidence read = reads[index];
            const std::uint32_t candidateBegin = read.candidateBegin;
            read.candidateBegin = 0;
            writePayload(output, read, checksum);
            for (std::uint32_t offset = 0; offset < read.candidateCount; ++offset) {
                writePayload(output, candidates[candidateBegin + offset], checksum);
            }
        }
        Trailer trailer = {};
        trailer.checksum = checksum;
        trailer.payloadBytes = payloadBytes;
        trailer.complete = kComplete;
        output.write(reinterpret_cast<const char *>(&trailer), sizeof(trailer));
        output.close();
        if (!output || std::rename(partial.c_str(), path.c_str()) != 0) {
            const int savedErrno = errno;
            std::remove(partial.c_str());
            throw std::runtime_error("cannot commit spatial spill " + path + ": "
                                     + std::strerror(savedErrno));
        }
        run.path = path;
        run.reads = reads.size();
        run.candidates = candidates.size();
        run.bytes = sizeof(Header) + payloadBytes + sizeof(Trailer);
        return true;
    } catch (const std::exception &exception) {
        error = exception.what();
        return false;
    }
}

struct Cursor::Impl {
    Run run;
    std::ifstream input;
    Header header;
    std::uint64_t readsSeen = 0;
    std::uint64_t candidatesSeen = 0;
    std::uint64_t payloadSeen = 0;
    std::uint64_t checksum = kFnvOffset;
    bool complete = false;

    Impl(const Run &value, const std::string &sourceRevision)
        : run(value), input(value.path.c_str(), std::ios::binary), header()
    {
        if (!input) throw std::runtime_error("cannot open spatial spill " + run.path);
        input.read(reinterpret_cast<char *>(&header), sizeof(header));
        if (!input || header.magic != kMagic || header.schema != kSchema
            || header.headerBytes != sizeof(Header)
            || header.key != kReadKeyGeneUmiOrdinal
            || header.readBytes != sizeof(ReadEvidence)
            || header.candidateBytes != sizeof(CandidateEvidence)) {
            throw std::runtime_error("wrong spatial spill schema/key/layout: " + run.path);
        }
        if (header.sourceRevisionHash != sourceHash(sourceRevision)) {
            throw std::runtime_error("wrong spatial spill source revision: " + run.path);
        }
        std::uint64_t expectedPayload = 0;
        if (!checkedPayloadBytes(header.readCount, header.candidateCount,
                                 expectedPayload)
            || expectedPayload != header.payloadBytes
            || (run.reads != 0 && run.reads != header.readCount)
            || (run.candidates != 0 && run.candidates != header.candidateCount)) {
            throw std::runtime_error("inconsistent spatial spill counts: " + run.path);
        }
    }

    void finish()
    {
        if (complete) return;
        if (candidatesSeen != header.candidateCount
            || payloadSeen != header.payloadBytes) {
            throw std::runtime_error("spatial spill payload count mismatch: " + run.path);
        }
        Trailer trailer = {};
        input.read(reinterpret_cast<char *>(&trailer), sizeof(trailer));
        if (!input || trailer.checksum != checksum
            || trailer.payloadBytes != payloadSeen || trailer.complete != kComplete) {
            throw std::runtime_error("spatial spill checksum/completion mismatch: " + run.path);
        }
        char extra = 0;
        if (input.read(&extra, 1)) {
            throw std::runtime_error("spatial spill has trailing data: " + run.path);
        }
        complete = true;
    }
};

std::unique_ptr<Cursor> Cursor::open(const Run &run,
                                     const std::string &sourceRevision,
                                     std::string &error)
{
    try {
        return std::unique_ptr<Cursor>(
            new Cursor(std::unique_ptr<Impl>(new Impl(run, sourceRevision))));
    } catch (const std::exception &exception) {
        error = exception.what();
        return std::unique_ptr<Cursor>();
    }
}

Cursor::Cursor(std::unique_ptr<Impl> impl) : impl_(std::move(impl)) {}
Cursor::~Cursor() {}

bool Cursor::next(ReadEvidence &read,
                  std::vector<CandidateEvidence> &candidates,
                  std::string &error)
{
    try {
        if (impl_->readsSeen == impl_->header.readCount) {
            impl_->finish();
            return false;
        }
        readPayload(impl_->input, read, impl_->checksum);
        impl_->payloadSeen += sizeof(ReadEvidence);
        if (read.candidateBegin != 0
            || impl_->candidatesSeen + read.candidateCount
                > impl_->header.candidateCount) {
            throw std::runtime_error("corrupt spatial spill read span: " + impl_->run.path);
        }
        candidates.resize(read.candidateCount);
        for (CandidateEvidence &candidate : candidates) {
            readPayload(impl_->input, candidate, impl_->checksum);
            impl_->payloadSeen += sizeof(CandidateEvidence);
        }
        ++impl_->readsSeen;
        impl_->candidatesSeen += candidates.size();
        return true;
    } catch (const std::exception &exception) {
        error = exception.what();
        return false;
    }
}

} // namespace spill
} // namespace spatial_gex
