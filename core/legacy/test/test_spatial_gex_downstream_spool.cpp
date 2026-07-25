#include "SpatialGexDownstreamSpool.h"

#include <cassert>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace {

const std::uint32_t kColumns = 3350;

void makeDirectory(const std::string &path)
{
    assert(::mkdir(path.c_str(), 0755) == 0);
}

void removeRun(const spatial_gex::downstream_spool::Run &run)
{
    assert(std::remove(run.path.c_str()) == 0);
}

spatial_gex::downstream_spool::Run copyRun(
    const spatial_gex::downstream_spool::Run &source,
    const std::string &destination)
{
    std::ifstream input(source.path.c_str(), std::ios::binary);
    std::ofstream output(destination.c_str(), std::ios::binary | std::ios::trunc);
    assert(input && output);
    output << input.rdbuf();
    output.flush();
    assert(!input.bad() && output);
    spatial_gex::downstream_spool::Run result = source;
    result.path = destination;
    return result;
}

void overwriteUint32(const std::string &path, std::streamoff offset,
                     std::uint32_t value)
{
    std::fstream file(path.c_str(), std::ios::binary | std::ios::in | std::ios::out);
    assert(file);
    file.seekp(offset);
    file.write(reinterpret_cast<const char *>(&value), sizeof(value));
    file.flush();
    assert(file);
}

void overwriteUint64(const std::string &path, std::streamoff offset,
                     std::uint64_t value)
{
    std::fstream file(path.c_str(), std::ios::binary | std::ios::in | std::ios::out);
    assert(file);
    file.seekp(offset);
    file.write(reinterpret_cast<const char *>(&value), sizeof(value));
    file.flush();
    assert(file);
}

void assertRejectedAfterRead(
    const spatial_gex::downstream_spool::Run &run,
    const std::string &revision)
{
    using namespace spatial_gex::downstream_spool;
    std::string error;
    std::unique_ptr<ContributionCursor> cursor = ContributionCursor::open(
        run, revision, error);
    assert(cursor && error.empty());
    Contribution value = {};
    while (cursor->next(value, error)) {}
    assert(!error.empty());
}

spatial_gex::downstream_spool::Contribution contribution(
    std::uint32_t gene, std::uint32_t coordinate, std::uint32_t umi,
    std::uint32_t clique, double posterior, std::uint8_t flags)
{
    spatial_gex::downstream_spool::Contribution value = {};
    value.posterior = posterior;
    value.gene = gene;
    value.coordinate = coordinate;
    value.rawUmi = umi;
    value.memberCount = clique + 1;
    value.cliqueOrdinal = clique;
    value.candidateOrdinal = static_cast<std::uint16_t>(clique);
    value.flags = flags;
    return value;
}

bool same(const spatial_gex::downstream_spool::Contribution &left,
          const spatial_gex::downstream_spool::Contribution &right)
{
    return left.posterior == right.posterior && left.gene == right.gene
        && left.coordinate == right.coordinate && left.rawUmi == right.rawUmi
        && left.memberCount == right.memberCount
        && left.cliqueOrdinal == right.cliqueOrdinal
        && left.candidateOrdinal == right.candidateOrdinal
        && left.flags == right.flags && left.reserved == right.reserved;
}

} // namespace

int main()
{
    using namespace spatial_gex::downstream_spool;
    char directoryTemplate[] = "/tmp/star_spatial_downstream_spool.XXXXXX";
    const char *directoryValue = ::mkdtemp(directoryTemplate);
    assert(directoryValue != NULL);
    const std::string root(directoryValue);
    const std::string contributionsDirectory = root + "/contributions";
    const std::string corruptDirectory = root + "/corrupt";
    const std::string truncatedDirectory = root + "/truncated";
    const std::string matrixDirectory = root + "/matrix";
    makeDirectory(contributionsDirectory);
    makeDirectory(corruptDirectory);
    makeDirectory(truncatedDirectory);
    makeDirectory(matrixDirectory);

    const std::string revision = "downstream-spool-test-revision";
    const std::uint32_t shards = 7;
    const std::uint32_t parent = 16 * kColumns + 24;
    const std::uint32_t parentShard = shardForCoordinate(parent, kColumns, shards);
    for (std::uint32_t row = 0; row < 8; ++row) {
        for (std::uint32_t column = 0; column < 8; ++column) {
            assert(shardForCoordinate(parent + row * kColumns + column,
                                      kColumns, shards) == parentShard);
        }
    }

    std::vector<Contribution> input;
    input.push_back(contribution(4, parent, 17, 0, 0.75,
                                 ContributionHard | ContributionGatedHard));
    input.push_back(contribution(4, parent + 1, 17, 1, 0.25, 0));
    std::uint32_t otherCoordinate = 0;
    while (shardForCoordinate(otherCoordinate, kColumns, shards) == parentShard) {
        otherCoordinate += 8;
    }
    input.push_back(contribution(9, otherCoordinate, 3, 2, 1.0,
                                 ContributionStrict | ContributionHard
                                     | ContributionGatedHard));

    std::string error;
    std::unique_ptr<ContributionWriter> writer = ContributionWriter::create(
        contributionsDirectory, revision, kColumns, shards, 4096, error);
    assert(writer && error.empty());
    std::vector<std::vector<Contribution> > expected(shards);
    for (const Contribution &value : input) {
        assert(writer->append(value, error));
        expected[shardForCoordinate(value.coordinate, kColumns, shards)].push_back(value);
    }
    std::vector<Run> runs;
    assert(writer->finish(runs, error));
    assert(writer->records() == input.size());
    assert(writer->bytes() > input.size() * sizeof(Contribution));
    assert(runs.size() == 2);
    std::uint64_t recordsRead = 0;
    std::uint32_t priorShard = 0;
    bool firstRun = true;
    for (const Run &run : runs) {
        assert(firstRun || priorShard < run.shard);
        firstRun = false;
        priorShard = run.shard;
        assert(run.kind == RecordKind::Contribution && run.shards == shards);
        assert(run.records == expected[run.shard].size());
        std::unique_ptr<ContributionCursor> cursor = ContributionCursor::open(
            run, revision, error);
        assert(cursor && error.empty());
        Contribution value = {};
        std::size_t index = 0;
        while (cursor->next(value, error)) {
            assert(index < expected[run.shard].size());
            assert(same(value, expected[run.shard][index++]));
            ++recordsRead;
        }
        assert(error.empty() && index == expected[run.shard].size());
    }
    assert(recordsRead == input.size());
    error.clear();
    assert(!ContributionCursor::open(runs.front(), "wrong-revision", error));
    assert(!error.empty());

    std::vector<Run> layoutRuns;
    Run wrongSchema = copyRun(runs.front(), corruptDirectory + "/wrong_schema.bin");
    overwriteUint32(wrongSchema.path, 8, 99);
    layoutRuns.push_back(wrongSchema);
    error.clear();
    assert(!ContributionCursor::open(wrongSchema, revision, error));
    assert(!error.empty());

    Run wrongRecordSize = copyRun(
        runs.front(), corruptDirectory + "/wrong_record_size.bin");
    overwriteUint32(wrongRecordSize.path, 20, sizeof(Contribution) + 1);
    layoutRuns.push_back(wrongRecordSize);
    error.clear();
    assert(!ContributionCursor::open(wrongRecordSize, revision, error));
    assert(!error.empty());

    Run wrongKey = copyRun(runs.front(), corruptDirectory + "/wrong_key.bin");
    overwriteUint32(wrongKey.path, 24, 99);
    layoutRuns.push_back(wrongKey);
    error.clear();
    assert(!ContributionCursor::open(wrongKey, revision, error));
    assert(!error.empty());

    Run incomplete = copyRun(runs.front(), corruptDirectory + "/incomplete.bin");
    struct stat incompleteInfo;
    assert(::stat(incomplete.path.c_str(), &incompleteInfo) == 0);
    overwriteUint64(incomplete.path, incompleteInfo.st_size - 8, 0);
    layoutRuns.push_back(incomplete);
    assertRejectedAfterRead(incomplete, revision);

    Run trailing = copyRun(runs.front(), corruptDirectory + "/trailing.bin");
    {
        std::ofstream file(trailing.path.c_str(), std::ios::binary | std::ios::app);
        assert(file);
        file.put('x');
        file.flush();
        assert(file);
    }
    layoutRuns.push_back(trailing);
    assertRejectedAfterRead(trailing, revision);

    std::unique_ptr<ContributionWriter> corruptWriter = ContributionWriter::create(
        corruptDirectory, revision, kColumns, 1, 4096, error);
    assert(corruptWriter);
    error.clear();
    assert(corruptWriter->append(input.front(), error));
    std::vector<Run> corruptRuns;
    assert(corruptWriter->finish(corruptRuns, error) && corruptRuns.size() == 1);
    {
        std::fstream file(corruptRuns.front().path.c_str(),
                          std::ios::binary | std::ios::in | std::ios::out);
        assert(file);
        file.seekg(72 + 5);
        char value = 0;
        file.read(&value, 1);
        value ^= 1;
        file.seekp(72 + 5);
        file.write(&value, 1);
    }
    error.clear();
    std::unique_ptr<ContributionCursor> cursor = ContributionCursor::open(
        corruptRuns.front(), revision, error);
    assert(cursor);
    Contribution value = {};
    while (cursor->next(value, error)) {}
    assert(!error.empty());

    std::unique_ptr<ContributionWriter> truncatedWriter = ContributionWriter::create(
        truncatedDirectory, revision, kColumns, 1, 4096, error);
    assert(truncatedWriter);
    error.clear();
    assert(truncatedWriter->append(input.front(), error));
    std::vector<Run> truncatedRuns;
    assert(truncatedWriter->finish(truncatedRuns, error)
           && truncatedRuns.size() == 1);
    struct stat info;
    assert(::stat(truncatedRuns.front().path.c_str(), &info) == 0);
    assert(::truncate(truncatedRuns.front().path.c_str(), info.st_size - 1) == 0);
    error.clear();
    cursor = ContributionCursor::open(truncatedRuns.front(), revision, error);
    assert(cursor);
    while (cursor->next(value, error)) {}
    assert(!error.empty());

    std::vector<MatrixRecord> matrixRecords(3);
    matrixRecords[0] = {UINT64_C(1), 1.25};
    matrixRecords[1] = {UINT64_C(2), 2.5};
    matrixRecords[2] = {UINT64_C(2), 3.75};
    Run matrixRun;
    error.clear();
    assert(writeMatrixRun(matrixDirectory, revision, 2, shards, 5,
                          matrixRecords, matrixRun, error));
    std::unique_ptr<MatrixCursor> matrixCursor = MatrixCursor::open(
        matrixRun, revision, error);
    assert(matrixCursor);
    MatrixRecord matrixValue = {};
    std::size_t matrixIndex = 0;
    while (matrixCursor->next(matrixValue, error)) {
        assert(matrixIndex < matrixRecords.size());
        assert(matrixValue.key == matrixRecords[matrixIndex].key);
        assert(matrixValue.value == matrixRecords[matrixIndex].value);
        ++matrixIndex;
    }
    assert(error.empty() && matrixIndex == matrixRecords.size());
    std::swap(matrixRecords[0], matrixRecords[2]);
    Run unsortedRun;
    error.clear();
    assert(!writeMatrixRun(matrixDirectory, revision, 3, shards, 6,
                           matrixRecords, unsortedRun, error));
    assert(!error.empty());

    for (const Run &run : runs) removeRun(run);
    for (const Run &run : layoutRuns) removeRun(run);
    removeRun(corruptRuns.front());
    removeRun(truncatedRuns.front());
    removeRun(matrixRun);
    assert(::rmdir(contributionsDirectory.c_str()) == 0);
    assert(::rmdir(corruptDirectory.c_str()) == 0);
    assert(::rmdir(truncatedDirectory.c_str()) == 0);
    assert(::rmdir(matrixDirectory.c_str()) == 0);
    assert(::rmdir(root.c_str()) == 0);
    return 0;
}
