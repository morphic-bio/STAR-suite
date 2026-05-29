#ifndef SLAM_DUMP_H
#define SLAM_DUMP_H

#include "SlamReadBuffer.h"
#include <cstdint>
#include <string>
#include <vector>
#include <atomic>
#include <fstream>

struct SlamDumpMetadata {
    uint32_t version = 2;
    // 0=Alignments (1/nTr), 1=Uniform (1.0); stored in header flags.
    uint32_t weightMode = 0;
    uint64_t nReads = 0;
    double errorRate = 0.0;
    double convRate = 0.0;
    std::vector<std::string> geneIds;
    std::vector<std::string> geneNames;
    // Optional per-gene mask from the live STAR-SLAM quantifier. Empty means
    // no filtering; otherwise 0 entries are excluded by SlamQuant::addRead().
    std::vector<uint8_t> allowedGenes;
    std::vector<std::string> chrNames;
    std::vector<uint64_t> chrStart;
};

struct SlamWeightKey {
    uint64_t h1 = 0;
    uint64_t h2 = 0;
};

struct SlamWeightRecord {
    SlamWeightKey key;
    double weight = 0.0;
};

struct SlamWeightMetadata {
    uint32_t version = 1;
    uint32_t flags = 0;
    uint64_t nReads = 0;
    uint32_t weightMode = 0;
};

struct SlamDumpPartInfo {
    std::string path;
    uint64_t nReads = 0;
};

class SlamDumpPartWriter {
public:
    SlamDumpPartWriter(const std::string& path,
                       uint64_t maxReads,
                       std::atomic<uint64_t>* globalCounter,
                       std::string* err);

    bool writeRead(const SlamBufferedRead& read, std::string* err = nullptr);
    uint64_t written() const { return written_; }
    const std::string& path() const { return path_; }
    bool ok() const { return ok_; }
    bool limitReached() const;
    void close();

private:
    std::ofstream out_;
    std::string path_;
    uint64_t written_ = 0;
    uint64_t maxReads_ = 0;
    std::atomic<uint64_t>* globalCounter_ = nullptr;
    bool ok_ = false;
};

// Write dump file from one or more buffers.
// maxReads=0 means no limit.
bool writeSlamDump(const std::string& path,
                   const SlamDumpMetadata& meta,
                   const std::vector<const SlamReadBuffer*>& buffers,
                   uint64_t maxReads,
                   std::string* err);

// Read dump file into metadata + reads.
bool readSlamDump(const std::string& path,
                  SlamDumpMetadata* meta,
                  std::vector<SlamBufferedRead>* reads,
                  std::string* err);

// Compute a stable key for weight matching from a buffered read.
SlamWeightKey computeSlamWeightKey(const SlamBufferedRead& read);

// Write/read weight sidecar file (keyed + ordered).
bool writeSlamWeights(const std::string& path,
                      const SlamDumpMetadata& dumpMeta,
                      const std::vector<const SlamReadBuffer*>& buffers,
                      uint64_t maxReads,
                      const std::vector<double>* overrideWeights = nullptr,
                      std::string* err = nullptr);

bool readSlamWeights(const std::string& path,
                     SlamWeightMetadata* meta,
                     std::vector<SlamWeightRecord>* records,
                     std::string* err);

// Merge per-thread dump parts into a final dump file.
bool mergeSlamDumpParts(const std::string& path,
                        const SlamDumpMetadata& meta,
                        const std::vector<SlamDumpPartInfo>& parts,
                        std::string* err);

// Write weight sidecar by streaming a dump file (optionally vbGene-reweighted).
bool writeSlamWeightsFromDump(const std::string& dumpPath,
                              const std::string& weightPath,
                              const std::vector<double>* vbGenePosterior,
                              std::string* err);

#endif // SLAM_DUMP_H
