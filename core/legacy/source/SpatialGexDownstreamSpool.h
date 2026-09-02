#ifndef STAR_SPATIAL_GEX_DOWNSTREAM_SPOOL_H
#define STAR_SPATIAL_GEX_DOWNSTREAM_SPOOL_H

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace spatial_gex {
namespace downstream_spool {

enum ContributionFlag : std::uint8_t {
    ContributionStrict = 1u << 0,
    ContributionHard = 1u << 1,
    ContributionGatedHard = 1u << 2
};

// One record is emitted for each candidate retained by a completed read
// clique. The record is coordinate-local after the posterior has been
// normalized, so downstream correction and reconciliation may be partitioned
// without retaining the original clique.
struct Contribution {
    double posterior;
    std::uint32_t gene;
    std::uint32_t coordinate;
    std::uint32_t rawUmi;
    std::uint32_t memberCount;
    std::uint32_t cliqueOrdinal;
    std::uint16_t candidateOrdinal;
    std::uint8_t flags;
    std::uint8_t reserved;
};

struct MatrixRecord {
    std::uint64_t key;
    double value;
};

static_assert(sizeof(Contribution) == 32,
              "spatial downstream contribution size changed");
static_assert(sizeof(MatrixRecord) == 16,
              "spatial downstream matrix record size changed");

enum class RecordKind : std::uint32_t {
    Contribution = 1,
    Matrix = 2
};

struct Run {
    std::string path;
    RecordKind kind = RecordKind::Contribution;
    std::uint32_t shard = 0;
    std::uint32_t shards = 0;
    std::uint64_t records = 0;
    std::uint64_t bytes = 0;
};

// Hashes a 16-um parent coordinate so all 2-um children of an 8/16-um
// spatial bin remain in the same downstream partition.
std::uint32_t shardForCoordinate(std::uint32_t coordinate,
                                 std::uint32_t gridColumns,
                                 std::uint32_t shards);

class ContributionWriter {
  public:
    static std::unique_ptr<ContributionWriter> create(
        const std::string &directory, const std::string &sourceRevision,
        std::uint32_t gridColumns, std::uint32_t shards,
        std::size_t bufferBytes, std::string &error);
    ~ContributionWriter();
    ContributionWriter(const ContributionWriter &) = delete;
    ContributionWriter &operator=(const ContributionWriter &) = delete;

    bool append(const Contribution &record, std::string &error);
    bool finish(std::vector<Run> &runs, std::string &error);
    std::uint64_t records() const;
    std::uint64_t bytes() const;

  private:
    struct Impl;
    explicit ContributionWriter(std::unique_ptr<Impl> impl);
    std::unique_ptr<Impl> impl_;
};

class ContributionCursor {
  public:
    static std::unique_ptr<ContributionCursor> open(
        const Run &run, const std::string &sourceRevision, std::string &error);
    ~ContributionCursor();
    ContributionCursor(const ContributionCursor &) = delete;
    ContributionCursor &operator=(const ContributionCursor &) = delete;

    bool next(Contribution &record, std::string &error);

  private:
    struct Impl;
    explicit ContributionCursor(std::unique_ptr<Impl> impl);
    std::unique_ptr<Impl> impl_;
};

bool writeMatrixRun(const std::string &directory,
                    const std::string &sourceRevision,
                    std::uint32_t shard, std::uint32_t shards,
                    std::uint64_t runIndex,
                    const std::vector<MatrixRecord> &records,
                    Run &run, std::string &error);

class MatrixCursor {
  public:
    static std::unique_ptr<MatrixCursor> open(
        const Run &run, const std::string &sourceRevision, std::string &error);
    ~MatrixCursor();
    MatrixCursor(const MatrixCursor &) = delete;
    MatrixCursor &operator=(const MatrixCursor &) = delete;

    bool next(MatrixRecord &record, std::string &error);

  private:
    struct Impl;
    explicit MatrixCursor(std::unique_ptr<Impl> impl);
    std::unique_ptr<Impl> impl_;
};

} // namespace downstream_spool
} // namespace spatial_gex

#endif
