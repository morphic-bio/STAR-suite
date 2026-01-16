#ifndef SLAM_DUMP_H
#define SLAM_DUMP_H

#include "SlamReadBuffer.h"
#include <cstdint>
#include <string>
#include <vector>

struct SlamDumpMetadata {
    uint32_t version = 1;
    // 0=Alignments (1/nTr), 1=Uniform (1.0); stored in header flags.
    uint32_t weightMode = 0;
    uint64_t nReads = 0;
    double errorRate = 0.0;
    double convRate = 0.0;
    std::vector<std::string> geneIds;
    std::vector<std::string> geneNames;
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

#endif // SLAM_DUMP_H
