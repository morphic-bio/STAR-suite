#include "SpatialGexDownstreamSpool.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>

namespace spatial_gex {
namespace downstream_spool {
namespace {

const std::uint64_t kMagic = UINT64_C(0x3150574e44475353); // "SSGDNWP1"
const std::uint32_t kSchema = 1;
const std::uint32_t kContributionKey = 1;
const std::uint32_t kMatrixKey = 2;
const std::uint64_t kComplete = UINT64_C(0x4554454c504d4f43); // "COMPLETE"
const std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
const std::uint64_t kFnvPrime = UINT64_C(1099511628211);

#pragma pack(push, 1)
struct Header {
    std::uint64_t magic;
    std::uint32_t schema;
    std::uint32_t headerBytes;
    std::uint32_t kind;
    std::uint32_t recordBytes;
    std::uint32_t key;
    std::uint32_t shard;
    std::uint32_t shards;
    std::uint32_t reserved32;
    std::uint64_t recordCount;
    std::uint64_t sourceRevisionHash;
    std::uint64_t payloadBytes;
    std::uint64_t reserved64;
};

struct Trailer {
    std::uint64_t checksum;
    std::uint64_t recordCount;
    std::uint64_t payloadBytes;
    std::uint64_t complete;
};
#pragma pack(pop)

static_assert(sizeof(Header) == 72, "downstream spool header size changed");
static_assert(sizeof(Trailer) == 32, "downstream spool trailer size changed");

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

std::uint32_t recordKey(RecordKind kind)
{
    return kind == RecordKind::Contribution ? kContributionKey : kMatrixKey;
}

std::string runName(const std::string &directory, RecordKind kind,
                    std::uint32_t shard, std::uint64_t runIndex)
{
    std::ostringstream name;
    name << directory << '/'
         << (kind == RecordKind::Contribution ? "contribution" : "matrix")
         << ".s" << std::setw(4) << std::setfill('0') << shard
         << ".r" << std::setw(8) << runIndex << ".bin";
    return name.str();
}

bool directoryExists(const std::string &path)
{
    struct stat info;
    return ::stat(path.c_str(), &info) == 0 && S_ISDIR(info.st_mode);
}

template <typename T>
void validateRecord(const T &)
{
}

template <>
void validateRecord<Contribution>(const Contribution &record)
{
    const std::uint8_t known = ContributionStrict | ContributionHard
        | ContributionGatedHard;
    if (!std::isfinite(record.posterior) || record.posterior < 0.0
        || record.posterior > 1.0 || record.memberCount == 0
        || (record.flags & ~known) != 0 || record.reserved != 0) {
        throw std::runtime_error("invalid spatial downstream contribution record");
    }
}

template <>
void validateRecord<MatrixRecord>(const MatrixRecord &record)
{
    if (!std::isfinite(record.value) || record.value < 0.0) {
        throw std::runtime_error("invalid spatial downstream matrix record");
    }
}

template <typename T>
class OpenRun {
  public:
    OpenRun(const std::string &directory, const std::string &sourceRevision,
            RecordKind kind, std::uint32_t shard, std::uint32_t shards,
            std::uint64_t runIndex, std::size_t bufferBytes)
        : finalPath_(runName(directory, kind, shard, runIndex)),
          partialPath_(finalPath_ + ".partial"), kind_(kind), shard_(shard),
          shards_(shards), sourceRevisionHash_(sourceHash(sourceRevision)),
          buffer_(new char[std::max<std::size_t>(bufferBytes, 4096)])
    {
        stream_.rdbuf()->pubsetbuf(buffer_.get(),
                                   static_cast<std::streamsize>(
                                       std::max<std::size_t>(bufferBytes, 4096)));
        stream_.open(partialPath_.c_str(), std::ios::binary | std::ios::in
                                           | std::ios::out | std::ios::trunc);
        if (!stream_) {
            throw std::runtime_error("cannot create downstream spool " + partialPath_);
        }
        Header placeholder = {};
        stream_.write(reinterpret_cast<const char *>(&placeholder), sizeof(placeholder));
        if (!stream_) throw std::runtime_error("cannot write downstream spool header");
    }

    ~OpenRun()
    {
        stream_.close();
        if (!committed_) std::remove(partialPath_.c_str());
    }

    void append(const T &record)
    {
        if (finished_) throw std::runtime_error("downstream spool run already finished");
        validateRecord(record);
        stream_.write(reinterpret_cast<const char *>(&record), sizeof(record));
        if (!stream_) throw std::runtime_error("cannot write downstream spool payload");
        checksum_ = fnvUpdate(checksum_, &record, sizeof(record));
        ++records_;
    }

    Run finish()
    {
        if (finished_) throw std::runtime_error("downstream spool run finished twice");
        finished_ = true;
        if (records_ > std::numeric_limits<std::uint64_t>::max() / sizeof(T)) {
            throw std::runtime_error("downstream spool payload size overflow");
        }
        const std::uint64_t payloadBytes = records_ * sizeof(T);
        Header header = {};
        header.magic = kMagic;
        header.schema = kSchema;
        header.headerBytes = sizeof(Header);
        header.kind = static_cast<std::uint32_t>(kind_);
        header.recordBytes = sizeof(T);
        header.key = recordKey(kind_);
        header.shard = shard_;
        header.shards = shards_;
        header.recordCount = records_;
        header.sourceRevisionHash = sourceRevisionHash_;
        header.payloadBytes = payloadBytes;

        stream_.seekp(0, std::ios::beg);
        stream_.write(reinterpret_cast<const char *>(&header), sizeof(header));
        stream_.seekp(0, std::ios::end);
        Trailer trailer = {};
        trailer.checksum = checksum_;
        trailer.recordCount = records_;
        trailer.payloadBytes = payloadBytes;
        trailer.complete = kComplete;
        stream_.write(reinterpret_cast<const char *>(&trailer), sizeof(trailer));
        stream_.flush();
        if (!stream_) throw std::runtime_error("cannot finalize downstream spool payload");
        stream_.close();
        if (!stream_) {
            std::remove(partialPath_.c_str());
            throw std::runtime_error("cannot close downstream spool " + partialPath_);
        }
        if (std::rename(partialPath_.c_str(), finalPath_.c_str()) != 0) {
            const int savedErrno = errno;
            std::remove(partialPath_.c_str());
            throw std::runtime_error("cannot commit downstream spool " + finalPath_
                                     + ": " + std::strerror(savedErrno));
        }
        committed_ = true;
        Run run;
        run.path = finalPath_;
        run.kind = kind_;
        run.shard = shard_;
        run.shards = shards_;
        run.records = records_;
        run.bytes = sizeof(Header) + payloadBytes + sizeof(Trailer);
        return run;
    }

    std::uint64_t records() const { return records_; }

  private:
    std::string finalPath_;
    std::string partialPath_;
    RecordKind kind_;
    std::uint32_t shard_;
    std::uint32_t shards_;
    std::uint64_t sourceRevisionHash_;
    std::unique_ptr<char[]> buffer_;
    std::fstream stream_;
    std::uint64_t checksum_ = kFnvOffset;
    std::uint64_t records_ = 0;
    bool finished_ = false;
    bool committed_ = false;
};

template <typename T>
class CursorState {
  public:
    CursorState(const Run &run, const std::string &sourceRevision,
                RecordKind expectedKind)
        : run_(run), input_(run.path.c_str(), std::ios::binary)
    {
        if (!input_) throw std::runtime_error("cannot open downstream spool " + run.path);
        input_.read(reinterpret_cast<char *>(&header_), sizeof(header_));
        if (!input_ || header_.magic != kMagic || header_.schema != kSchema
            || header_.headerBytes != sizeof(Header)
            || header_.kind != static_cast<std::uint32_t>(expectedKind)
            || header_.recordBytes != sizeof(T)
            || header_.key != recordKey(expectedKind)
            || header_.shards == 0 || header_.shard >= header_.shards
            || header_.reserved32 != 0 || header_.reserved64 != 0) {
            throw std::runtime_error("wrong downstream spool schema/key/layout: " + run.path);
        }
        if (header_.sourceRevisionHash != sourceHash(sourceRevision)) {
            throw std::runtime_error("wrong downstream spool source revision: " + run.path);
        }
        if (header_.recordCount > std::numeric_limits<std::uint64_t>::max() / sizeof(T)
            || header_.payloadBytes != header_.recordCount * sizeof(T)
            || run.kind != expectedKind || run.shard != header_.shard
            || (run.shards != 0 && run.shards != header_.shards)
            || (run.records != 0 && run.records != header_.recordCount)) {
            throw std::runtime_error("inconsistent downstream spool counts: " + run.path);
        }
    }

    bool next(T &record)
    {
        if (recordsSeen_ == header_.recordCount) {
            finish();
            return false;
        }
        input_.read(reinterpret_cast<char *>(&record), sizeof(record));
        if (!input_) throw std::runtime_error("truncated downstream spool payload: " + run_.path);
        validateRecord(record);
        checksum_ = fnvUpdate(checksum_, &record, sizeof(record));
        ++recordsSeen_;
        return true;
    }

  private:
    void finish()
    {
        if (complete_) return;
        Trailer trailer = {};
        input_.read(reinterpret_cast<char *>(&trailer), sizeof(trailer));
        if (!input_ || trailer.checksum != checksum_
            || trailer.recordCount != recordsSeen_
            || trailer.payloadBytes != recordsSeen_ * sizeof(T)
            || trailer.complete != kComplete) {
            throw std::runtime_error("downstream spool checksum/completion mismatch: "
                                     + run_.path);
        }
        char extra = 0;
        if (input_.read(&extra, 1)) {
            throw std::runtime_error("downstream spool has trailing data: " + run_.path);
        }
        complete_ = true;
    }

    Run run_;
    std::ifstream input_;
    Header header_ = {};
    std::uint64_t recordsSeen_ = 0;
    std::uint64_t checksum_ = kFnvOffset;
    bool complete_ = false;
};

} // namespace

std::uint32_t shardForCoordinate(std::uint32_t coordinate,
                                 std::uint32_t gridColumns,
                                 std::uint32_t shards)
{
    if (gridColumns == 0 || shards == 0) {
        throw std::invalid_argument("downstream spool requires positive grid/shard counts");
    }
    const std::uint64_t row = coordinate / gridColumns;
    const std::uint64_t column = coordinate % gridColumns;
    const std::uint64_t parentColumns = (gridColumns + 7) / 8;
    std::uint64_t value = (row / 8) * parentColumns + column / 8;
    // Fixed SplitMix64 finalizer; only determinism and distribution matter.
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    value ^= value >> 31;
    return static_cast<std::uint32_t>(value % shards);
}

struct ContributionWriter::Impl {
    std::string directory;
    std::string sourceRevision;
    std::uint32_t gridColumns = 0;
    std::uint32_t shards = 0;
    std::size_t bufferBytes = 0;
    std::vector<std::unique_ptr<OpenRun<Contribution> > > writers;
    std::uint64_t totalRecords = 0;
    std::uint64_t totalBytes = 0;
    bool finished = false;
};

std::unique_ptr<ContributionWriter> ContributionWriter::create(
    const std::string &directory, const std::string &sourceRevision,
    std::uint32_t gridColumns, std::uint32_t shards,
    std::size_t bufferBytes, std::string &error)
{
    try {
        if (!directoryExists(directory) || sourceRevision.empty()
            || gridColumns == 0 || shards == 0 || shards > 4096
            || bufferBytes == 0) {
            throw std::invalid_argument("invalid downstream contribution writer configuration");
        }
        std::unique_ptr<Impl> impl(new Impl());
        impl->directory = directory;
        impl->sourceRevision = sourceRevision;
        impl->gridColumns = gridColumns;
        impl->shards = shards;
        impl->bufferBytes = bufferBytes;
        impl->writers.resize(shards);
        return std::unique_ptr<ContributionWriter>(
            new ContributionWriter(std::move(impl)));
    } catch (const std::exception &exception) {
        error = exception.what();
        return std::unique_ptr<ContributionWriter>();
    }
}

ContributionWriter::ContributionWriter(std::unique_ptr<Impl> impl)
    : impl_(std::move(impl))
{
}

ContributionWriter::~ContributionWriter() {}

bool ContributionWriter::append(const Contribution &record, std::string &error)
{
    try {
        if (impl_->finished) throw std::runtime_error("contribution writer already finished");
        const std::uint32_t shard = shardForCoordinate(
            record.coordinate, impl_->gridColumns, impl_->shards);
        if (!impl_->writers[shard]) {
            impl_->writers[shard].reset(new OpenRun<Contribution>(
                impl_->directory, impl_->sourceRevision, RecordKind::Contribution,
                shard, impl_->shards, 0, impl_->bufferBytes));
        }
        impl_->writers[shard]->append(record);
        ++impl_->totalRecords;
        return true;
    } catch (const std::exception &exception) {
        error = exception.what();
        return false;
    }
}

bool ContributionWriter::finish(std::vector<Run> &runs, std::string &error)
{
    try {
        if (impl_->finished) throw std::runtime_error("contribution writer finished twice");
        impl_->finished = true;
        runs.clear();
        runs.reserve(impl_->shards);
        std::uint64_t records = 0;
        for (std::uint32_t shard = 0; shard < impl_->shards; ++shard) {
            if (!impl_->writers[shard]) continue;
            Run run = impl_->writers[shard]->finish();
            records += run.records;
            impl_->totalBytes += run.bytes;
            runs.push_back(run);
            impl_->writers[shard].reset();
        }
        if (records != impl_->totalRecords) {
            throw std::runtime_error("downstream contribution row conservation failure");
        }
        return true;
    } catch (const std::exception &exception) {
        error = exception.what();
        return false;
    }
}

std::uint64_t ContributionWriter::records() const { return impl_->totalRecords; }
std::uint64_t ContributionWriter::bytes() const { return impl_->totalBytes; }

struct ContributionCursor::Impl {
    explicit Impl(const Run &run, const std::string &sourceRevision)
        : state(run, sourceRevision, RecordKind::Contribution) {}
    CursorState<Contribution> state;
};

std::unique_ptr<ContributionCursor> ContributionCursor::open(
    const Run &run, const std::string &sourceRevision, std::string &error)
{
    try {
        return std::unique_ptr<ContributionCursor>(
            new ContributionCursor(std::unique_ptr<Impl>(
                new Impl(run, sourceRevision))));
    } catch (const std::exception &exception) {
        error = exception.what();
        return std::unique_ptr<ContributionCursor>();
    }
}

ContributionCursor::ContributionCursor(std::unique_ptr<Impl> impl)
    : impl_(std::move(impl))
{
}

ContributionCursor::~ContributionCursor() {}

bool ContributionCursor::next(Contribution &record, std::string &error)
{
    try {
        return impl_->state.next(record);
    } catch (const std::exception &exception) {
        error = exception.what();
        return false;
    }
}

bool writeMatrixRun(const std::string &directory,
                    const std::string &sourceRevision,
                    std::uint32_t shard, std::uint32_t shards,
                    std::uint64_t runIndex,
                    const std::vector<MatrixRecord> &records,
                    Run &run, std::string &error)
{
    try {
        if (!directoryExists(directory) || sourceRevision.empty()
            || shards == 0 || shard >= shards) {
            throw std::invalid_argument("invalid downstream matrix-run configuration");
        }
        for (std::size_t index = 1; index < records.size(); ++index) {
            if (records[index].key < records[index - 1].key) {
                throw std::runtime_error("downstream matrix run is not key-sorted");
            }
        }
        OpenRun<MatrixRecord> writer(directory, sourceRevision, RecordKind::Matrix,
                                     shard, shards, runIndex, 64 * 1024);
        for (const MatrixRecord &record : records) writer.append(record);
        run = writer.finish();
        return true;
    } catch (const std::exception &exception) {
        error = exception.what();
        return false;
    }
}

struct MatrixCursor::Impl {
    explicit Impl(const Run &run, const std::string &sourceRevision)
        : state(run, sourceRevision, RecordKind::Matrix) {}
    CursorState<MatrixRecord> state;
};

std::unique_ptr<MatrixCursor> MatrixCursor::open(
    const Run &run, const std::string &sourceRevision, std::string &error)
{
    try {
        return std::unique_ptr<MatrixCursor>(
            new MatrixCursor(std::unique_ptr<Impl>(new Impl(run, sourceRevision))));
    } catch (const std::exception &exception) {
        error = exception.what();
        return std::unique_ptr<MatrixCursor>();
    }
}

MatrixCursor::MatrixCursor(std::unique_ptr<Impl> impl)
    : impl_(std::move(impl))
{
}

MatrixCursor::~MatrixCursor() {}

bool MatrixCursor::next(MatrixRecord &record, std::string &error)
{
    try {
        return impl_->state.next(record);
    } catch (const std::exception &exception) {
        error = exception.what();
        return false;
    }
}

} // namespace downstream_spool
} // namespace spatial_gex
