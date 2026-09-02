#ifndef TRIM_QC_SHARD_H
#define TRIM_QC_SHARD_H

#include <string>
#include <cstdint>
#include "TrimQc.h"

// Per-shard trim-QC exchange.
//
// TrimQcCollector is a purely additive histogram accumulator and STAR already
// merges one per chunk thread. That makes it shardable for free: each mapping
// shard dumps its own collector next to its evidence sidecar, and the gather
// merges them. No BAM is needed and no read is revisited, which is the whole
// point -- the reads are only in memory once, during mapping.
//
// The eight trim counters travel with it because writeTrimQcJson needs them and
// they are plain sums.
struct TrimQcShardCounters {
    uint64_t reads_processed = 0;
    uint64_t reads_trimmed = 0;
    uint64_t reads_too_short = 0;
    uint64_t bases_quality_trimmed = 0;
    uint64_t bases_adapter_trimmed = 0;
    uint64_t pairs_processed = 0;
    uint64_t pairs_dropped = 0;
    uint64_t pairs_kept = 0;

    void merge(const TrimQcShardCounters& o);
};

bool writeTrimQcShard(const TrimQcCollector& qc,
                      const TrimQcShardCounters& counters,
                      const std::string& path,
                      std::string* error);

bool readTrimQcShard(const std::string& path,
                     TrimQcCollector& qc,
                     TrimQcShardCounters& counters,
                     std::string* error);

#endif
