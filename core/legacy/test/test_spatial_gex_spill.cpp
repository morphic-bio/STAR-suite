#include "SpatialGexSpill.h"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace {

std::vector<char> bytes(const std::string &path)
{
    std::ifstream input(path.c_str(), std::ios::binary);
    return std::vector<char>(std::istreambuf_iterator<char>(input),
                             std::istreambuf_iterator<char>());
}

void removeRun(const spatial_gex::spill::Run &run)
{
    if (!run.path.empty()) assert(std::remove(run.path.c_str()) == 0);
}

} // namespace

int main()
{
    char directoryTemplate[] = "/tmp/star_spatial_spill_test.XXXXXX";
    const char *directoryValue = ::mkdtemp(directoryTemplate);
    assert(directoryValue != NULL);
    const std::string directory(directoryValue);
    const std::string revision = "test-revision";

    std::vector<spatial_gex::ReadEvidence> reads(3);
    reads[0] = {2, 9, 0, 20, 1, 0};
    reads[1] = {1, 8, 1, 10, 2, 0};
    reads[2] = {1, 8, 3, 11, 1, 0};
    std::vector<spatial_gex::CandidateEvidence> candidates(4);
    candidates[0] = {90, 1, -1.0};
    candidates[1] = {80, 2, -2.0};
    candidates[2] = {81, 3, -3.0};
    candidates[3] = {82, 4, -4.0};

    std::string error;
    spatial_gex::spill::Run first;
    assert(spatial_gex::spill::writeRun(directory, revision, 3, 0, reads,
                                        candidates, first, error));
    assert(first.reads == 3 && first.candidates == 4 && first.bytes > 0);
    std::unique_ptr<spatial_gex::spill::Cursor> cursor =
        spatial_gex::spill::Cursor::open(first, revision, error);
    assert(cursor);
    spatial_gex::ReadEvidence read;
    std::vector<spatial_gex::CandidateEvidence> readCandidates;
    assert(cursor->next(read, readCandidates, error));
    assert(read.geneIndex == 1 && read.rawUmi == 8 && read.sourceOrdinal == 10);
    assert(readCandidates.size() == 2 && readCandidates[0].coordinateIndex == 80);
    assert(cursor->next(read, readCandidates, error));
    assert(read.geneIndex == 1 && read.sourceOrdinal == 11);
    assert(cursor->next(read, readCandidates, error));
    assert(read.geneIndex == 2 && read.sourceOrdinal == 20);
    error.clear();
    assert(!cursor->next(read, readCandidates, error) && error.empty());

    // The run path and run ordinal are not part of the payload. Identical
    // segments therefore produce byte-identical completed runs.
    spatial_gex::spill::Run second;
    assert(spatial_gex::spill::writeRun(directory, revision, 3, 1, reads,
                                        candidates, second, error));
    assert(bytes(first.path) == bytes(second.path));

    spatial_gex::spill::Run empty;
    const std::vector<spatial_gex::ReadEvidence> noReads;
    const std::vector<spatial_gex::CandidateEvidence> noCandidates;
    assert(spatial_gex::spill::writeRun(directory, revision, 4, 0, noReads,
                                        noCandidates, empty, error));
    cursor = spatial_gex::spill::Cursor::open(empty, revision, error);
    assert(cursor);
    error.clear();
    assert(!cursor->next(read, readCandidates, error) && error.empty());

    error.clear();
    assert(!spatial_gex::spill::Cursor::open(first, "wrong-revision", error));
    assert(!error.empty());

    spatial_gex::spill::Run corrupt;
    assert(spatial_gex::spill::writeRun(directory, revision, 5, 0, reads,
                                        candidates, corrupt, error));
    {
        std::fstream file(corrupt.path.c_str(), std::ios::binary | std::ios::in
                                                   | std::ios::out);
        assert(file);
        file.seekg(73);
        char value = 0;
        file.read(&value, 1);
        value ^= 1;
        file.seekp(73);
        file.write(&value, 1);
    }
    cursor = spatial_gex::spill::Cursor::open(corrupt, revision, error);
    assert(cursor);
    error.clear();
    while (cursor->next(read, readCandidates, error)) {}
    assert(!error.empty());

    spatial_gex::spill::Run truncated;
    assert(spatial_gex::spill::writeRun(directory, revision, 6, 0, reads,
                                        candidates, truncated, error));
    struct stat info;
    assert(::stat(truncated.path.c_str(), &info) == 0 && info.st_size > 1);
    assert(::truncate(truncated.path.c_str(), info.st_size - 1) == 0);
    cursor = spatial_gex::spill::Cursor::open(truncated, revision, error);
    assert(cursor);
    error.clear();
    while (cursor->next(read, readCandidates, error)) {}
    assert(!error.empty());

    removeRun(first);
    removeRun(second);
    removeRun(empty);
    removeRun(corrupt);
    removeRun(truncated);
    assert(::rmdir(directory.c_str()) == 0);
    return 0;
}
