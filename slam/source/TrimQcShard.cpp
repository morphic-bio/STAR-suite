#include "TrimQcShard.h"

#include <fstream>
#include <vector>

namespace {

const uint32_t kTrimQcShardMagic = 0x51425654u;   // "TVBQ"
const uint32_t kTrimQcShardVersion = 1;

template <typename T>
void put(std::ostream& o, const T& v) {
    o.write(reinterpret_cast<const char*>(&v), sizeof(T));
}
template <typename T>
bool get(std::istream& i, T& v) {
    i.read(reinterpret_cast<char*>(&v), sizeof(T));
    return static_cast<bool>(i);
}

void putVec(std::ostream& o, const std::vector<uint64_t>& v) {
    const uint64_t n = v.size();
    put(o, n);
    if (n != 0) o.write(reinterpret_cast<const char*>(v.data()),
                        static_cast<std::streamsize>(n * sizeof(uint64_t)));
}
bool getVec(std::istream& i, std::vector<uint64_t>& v) {
    uint64_t n = 0;
    if (!get(i, n)) return false;
    v.assign(static_cast<size_t>(n), 0);
    if (n != 0) i.read(reinterpret_cast<char*>(v.data()),
                       static_cast<std::streamsize>(n * sizeof(uint64_t)));
    return static_cast<bool>(i);
}

}  // namespace

void TrimQcShardCounters::merge(const TrimQcShardCounters& o) {
    reads_processed        += o.reads_processed;
    reads_trimmed          += o.reads_trimmed;
    reads_too_short        += o.reads_too_short;
    bases_quality_trimmed  += o.bases_quality_trimmed;
    bases_adapter_trimmed  += o.bases_adapter_trimmed;
    pairs_processed        += o.pairs_processed;
    pairs_dropped          += o.pairs_dropped;
    pairs_kept             += o.pairs_kept;
}

bool writeTrimQcShard(const TrimQcCollector& qc,
                      const TrimQcShardCounters& counters,
                      const std::string& path,
                      std::string* error) {
    std::ofstream o(path.c_str(), std::ios::binary);
    if (!o) {
        if (error) *error = "cannot open " + path + " for writing";
        return false;
    }
    put(o, kTrimQcShardMagic);
    put(o, kTrimQcShardVersion);
    put(o, qc.qualityBase());
    put(o, qc.totalReads());
    put(o, counters);
    const uint64_t nMates = qc.mateCount();
    put(o, nMates);
    for (uint64_t m = 0; m < nMates; ++m) {
        const TrimQcMateStats& s = qc.mate(static_cast<size_t>(m));
        put(o, s.reads);
        putVec(o, s.lengthHist);
        putVec(o, s.posCount);
        putVec(o, s.posQualSum);
        const uint64_t nPos = s.posBaseCounts.size();
        put(o, nPos);
        for (uint64_t p = 0; p < nPos; ++p)
            o.write(reinterpret_cast<const char*>(s.posBaseCounts[p].data()),
                    static_cast<std::streamsize>(5 * sizeof(uint64_t)));
        o.write(reinterpret_cast<const char*>(s.gcHist.data()),
                static_cast<std::streamsize>(101 * sizeof(uint64_t)));
    }
    o.close();
    if (!o) {
        if (error) *error = "write failed for " + path;
        return false;
    }
    return true;
}

bool readTrimQcShard(const std::string& path,
                     TrimQcCollector& qc,
                     TrimQcShardCounters& counters,
                     std::string* error) {
    std::ifstream i(path.c_str(), std::ios::binary);
    if (!i) {
        if (error) *error = "cannot read " + path;
        return false;
    }
    uint32_t magic = 0, version = 0, qualBase = 33;
    uint64_t totalReads = 0, nMates = 0;
    if (!get(i, magic) || !get(i, version) || magic != kTrimQcShardMagic ||
        version != kTrimQcShardVersion) {
        if (error) *error = path + " is not a trim-QC shard dump";
        return false;
    }
    if (!get(i, qualBase) || !get(i, totalReads) || !get(i, counters) || !get(i, nMates)) {
        if (error) *error = "truncated header in " + path;
        return false;
    }
    std::vector<TrimQcMateStats> mates(static_cast<size_t>(nMates));
    for (uint64_t m = 0; m < nMates; ++m) {
        TrimQcMateStats& s = mates[static_cast<size_t>(m)];
        uint64_t nPos = 0;
        if (!get(i, s.reads) || !getVec(i, s.lengthHist) || !getVec(i, s.posCount) ||
            !getVec(i, s.posQualSum) || !get(i, nPos)) {
            if (error) *error = "truncated mate block in " + path;
            return false;
        }
        s.posBaseCounts.assign(static_cast<size_t>(nPos), {});
        for (uint64_t p = 0; p < nPos; ++p)
            i.read(reinterpret_cast<char*>(s.posBaseCounts[static_cast<size_t>(p)].data()),
                   static_cast<std::streamsize>(5 * sizeof(uint64_t)));
        i.read(reinterpret_cast<char*>(s.gcHist.data()),
               static_cast<std::streamsize>(101 * sizeof(uint64_t)));
        if (!i) {
            if (error) *error = "truncated mate payload in " + path;
            return false;
        }
    }
    qc.restoreFromShard(qualBase, totalReads, mates);
    return true;
}
