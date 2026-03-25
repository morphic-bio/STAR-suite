// Binary snapshot for non-Flex direct-hash bridge: exact inlineHash_ layout +
// read-accounting fields consumed by populateBridgeReadAccounting().
//
// v2: bulk I/O (block read/write), kh_resize pre-sizing, no sort on write.

#include "SoloFeature_bridgeHashSnapshot.h"
#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "Parameters.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include "SoloFeatureTypes.h"
#include "TimeFunctions.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>

#include "hash_shims_cpp_compat.h"

namespace {

static const char *nonEmptyEnv(const char *name)
{
    const char *v = std::getenv(name);
    return (v != nullptr && v[0] != '\0') ? v : nullptr;
}

constexpr char kMagic[8] = {'S', 'T', 'A', 'R', 'B', 'G', '2', '\0'};
constexpr uint32_t kVersion = 2u;

#pragma pack(push, 1)
struct SnapshotHeader {
    char magic[8];
    uint32_t version;
    uint32_t featureType;
    uint32_t cbWLsize;
    uint32_t runThreadN;
    uint64_t statsReadN;
    uint64_t totalHashEntries;
    uint32_t nImmediatePairs;
    uint32_t pinVecLen;
    uint32_t nFlagCbEntries;
    uint32_t reserved32;
};

struct HashRow {
    uint64_t key;
    uint32_t val;
    uint32_t pad;
};
#pragma pack(pop)

static_assert(sizeof(HashRow) == 16, "HashRow must be 16 bytes for bulk I/O");

static uint64_t liveHashSize(khash_t(cg_agg) *hash)
{
    if (hash == nullptr)
        return 0;
    return kh_size(hash);
}

static void writeHashSection(std::ofstream &o, khash_t(cg_agg) *hash)
{
    const uint64_t n = liveHashSize(hash);
    o.write(reinterpret_cast<const char *>(&n), sizeof(n));
    if (n == 0 || hash == nullptr)
        return;

    constexpr size_t kChunk = 1u << 18;   // 256K rows = 4 MiB per write
    std::vector<HashRow> buf;
    buf.reserve(std::min<size_t>(n, kChunk));

    for (khiter_t it = kh_begin(hash); it != kh_end(hash); ++it) {
        if (!kh_exist(hash, it))
            continue;
        buf.push_back({kh_key(hash, it), kh_val(hash, it), 0u});
        if (buf.size() == kChunk) {
            o.write(reinterpret_cast<const char *>(buf.data()),
                    static_cast<std::streamsize>(buf.size() * sizeof(HashRow)));
            buf.clear();
        }
    }
    if (!buf.empty())
        o.write(reinterpret_cast<const char *>(buf.data()),
                static_cast<std::streamsize>(buf.size() * sizeof(HashRow)));
}

static void readHashSection(std::ifstream &in, khash_t(cg_agg) *hash, Parameters &P)
{
    uint64_t n = 0;
    if (!in.read(reinterpret_cast<char *>(&n), sizeof(n)))
        exitWithError("EXITING: bridge snapshot read failed (hash section size)\n", std::cerr, P.inOut->logMain,
                      EXIT_CODE_INPUT_FILES, P);
    if (n == 0)
        return;

    kh_resize(cg_agg, hash, static_cast<khint_t>(n * 2));

    constexpr size_t kChunk = 1u << 18;
    std::vector<HashRow> buf(std::min<size_t>(n, kChunk));
    uint64_t remaining = n;
    while (remaining > 0) {
        const size_t batch = static_cast<size_t>(std::min<uint64_t>(remaining, kChunk));
        if (!in.read(reinterpret_cast<char *>(buf.data()),
                     static_cast<std::streamsize>(batch * sizeof(HashRow))))
            exitWithError("EXITING: bridge snapshot read failed (hash bulk read)\n", std::cerr, P.inOut->logMain,
                          EXIT_CODE_INPUT_FILES, P);
        for (size_t i = 0; i < batch; ++i) {
            int absent = 0;
            khiter_t it = kh_put(cg_agg, hash, buf[i].key, &absent);
            kh_val(hash, it) = buf[i].val;
        }
        remaining -= batch;
    }
}

static bool bridgeFeatureOk(int32_t ft)
{
    switch (ft) {
        case SoloFeatureTypes::Gene:
        case SoloFeatureTypes::GeneFull:
        case SoloFeatureTypes::GeneFull_Ex50pAS:
        case SoloFeatureTypes::GeneFull_ExonOverIntron:
            return true;
        default:
            return false;
    }
}

} // namespace

namespace solo_bridge_hash_snapshot {

bool replaySkipReadsEnabled(const Parameters &P)
{
    return nonEmptyEnv("STAR_SOLO_BRIDGE_HASH_SNAPSHOT_IN") != nullptr
        && nonEmptyEnv("STAR_SOLO_BRIDGE_HASH_SNAPSHOT_REPLAY_SKIP_READS") != nullptr && P.pSolo.inlineHashMode
        && !P.pSolo.flexMode && nonEmptyEnv("STAR_SOLO_NONFLEX_HASH_BRIDGE") != nullptr;
}

bool stopAfterCountEnabled(const Parameters &P)
{
    return replaySkipReadsEnabled(P) && nonEmptyEnv("STAR_SOLO_BRIDGE_SNAPSHOT_STOP_AFTER_COUNT") != nullptr;
}

} // namespace solo_bridge_hash_snapshot

void SoloFeature::bridgeHashSnapshotWrite(const char *path)
{
    if (path == nullptr || path[0] == '\0')
        return;
    if (!readFeatSum)
        return;

    uint64_t totalEntries = 0;
    for (int ii = 0; ii < P.runThreadN; ++ii) {
        if (readFeatAll[ii] != nullptr)
            totalEntries += liveHashSize(readFeatAll[ii]->inlineHash_);
    }
    totalEntries += liveHashSize(readFeatSum->inlineHash_);

    SnapshotHeader h{};
    std::memcpy(h.magic, kMagic, sizeof(kMagic));
    h.version = kVersion;
    h.featureType = static_cast<uint32_t>(featureType);
    h.cbWLsize = pSolo.cbWLsize;
    h.runThreadN = static_cast<uint32_t>(P.runThreadN);
    h.statsReadN = g_statsAll.readN;
    h.totalHashEntries = totalEntries;
    h.nImmediatePairs = static_cast<uint32_t>(readFeatSum->bridgeImmediateReadCounts_.size());
    h.pinVecLen = static_cast<uint32_t>(readFeatSum->bridgePinNreadUnique_.size());
    h.nFlagCbEntries = static_cast<uint32_t>(readFeatSum->readFlag.flagCounts.size());

    std::ofstream o(path, std::ios::binary | std::ios::trunc);
    if (!o)
        exitWithError(std::string("EXITING: cannot open STAR_SOLO_BRIDGE_HASH_SNAPSHOT_OUT file: ") + path + "\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);

    o.write(reinterpret_cast<const char *>(&h), sizeof(h));

    for (int ii = 0; ii < P.runThreadN; ++ii)
        writeHashSection(o, readFeatAll[ii] != nullptr ? readFeatAll[ii]->inlineHash_ : nullptr);
    writeHashSection(o, readFeatSum->inlineHash_);

    // bridgeImmediateReadCounts_: bulk write as flat (cb, val) pairs
    {
        struct ImmRow { uint32_t cb; uint32_t pad; uint64_t val; };
        static_assert(sizeof(ImmRow) == 16, "ImmRow must be 16 bytes");
        std::vector<ImmRow> immBuf;
        immBuf.reserve(readFeatSum->bridgeImmediateReadCounts_.size());
        for (const auto &kv : readFeatSum->bridgeImmediateReadCounts_)
            immBuf.push_back({kv.first, 0u, kv.second});
        o.write(reinterpret_cast<const char *>(immBuf.data()),
                static_cast<std::streamsize>(immBuf.size() * sizeof(ImmRow)));
    }

    // Pin vectors: bulk write
    if (h.pinVecLen > 0) {
        o.write(reinterpret_cast<const char *>(readFeatSum->bridgePinNreadUnique_.data()),
                static_cast<std::streamsize>(h.pinVecLen * sizeof(uint32_t)));
        o.write(reinterpret_cast<const char *>(readFeatSum->bridgePinNreadMulti_.data()),
                static_cast<std::streamsize>(h.pinVecLen * sizeof(uint32_t)));
    }

    // flagCountsNoCB: bulk write
    o.write(reinterpret_cast<const char *>(readFeatSum->readFlag.flagCountsNoCB.data()),
            static_cast<std::streamsize>(SoloReadFlagClass::nBits * sizeof(uint64_t)));

    // flagCounts per-CB: bulk write (cb key + nBits values)
    {
        constexpr size_t rowBytes = sizeof(uint64_t) + SoloReadFlagClass::nBits * sizeof(uint64_t);
        std::vector<char> rowBuf(rowBytes);
        for (const auto &kv : readFeatSum->readFlag.flagCounts) {
            std::memcpy(rowBuf.data(), &kv.first, sizeof(uint64_t));
            std::memcpy(rowBuf.data() + sizeof(uint64_t), kv.second.data(),
                        SoloReadFlagClass::nBits * sizeof(uint64_t));
            o.write(rowBuf.data(), static_cast<std::streamsize>(rowBytes));
        }
    }

    // stats
    o.write(reinterpret_cast<const char *>(readFeatSum->stats.V),
            static_cast<std::streamsize>(SoloReadFeatureStats::nStats * sizeof(uint64_t)));

    // cbReadCount: bulk write
    if (pSolo.cbWLyes && readFeatSum->cbReadCount.size() == pSolo.cbWLsize) {
        o.write(reinterpret_cast<const char *>(readFeatSum->cbReadCount.data()),
                static_cast<std::streamsize>(pSolo.cbWLsize * sizeof(uint32_t)));
    } else {
        std::vector<uint32_t> zeros(pSolo.cbWLsize, 0u);
        o.write(reinterpret_cast<const char *>(zeros.data()),
                static_cast<std::streamsize>(pSolo.cbWLsize * sizeof(uint32_t)));
    }

    o.close();
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... STAR_SOLO_BRIDGE_HASH_SNAPSHOT_OUT wrote " << totalEntries
                     << " bridge hash entries to " << path << endl;
}

void SoloFeature::bridgeHashSnapshotLoad(const char *path)
{
    const auto t0 = std::chrono::steady_clock::now();

    std::ifstream in(path, std::ios::binary);
    if (!in)
        exitWithError(std::string("EXITING: cannot open STAR_SOLO_BRIDGE_HASH_SNAPSHOT_IN file: ") + path + "\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);

    SnapshotHeader h{};
    if (!in.read(reinterpret_cast<char *>(&h), sizeof(h)))
        exitWithError("EXITING: bridge snapshot truncated (header)\n", std::cerr, P.inOut->logMain,
                      EXIT_CODE_INPUT_FILES, P);

    if (std::memcmp(h.magic, kMagic, sizeof(kMagic)) != 0)
        exitWithError("EXITING: bridge snapshot bad magic (not a STARBG2 file). Re-create the snapshot with the "
                      "current binary.\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    if (h.version != kVersion)
        exitWithError("EXITING: bridge snapshot version mismatch\n", std::cerr, P.inOut->logMain,
                      EXIT_CODE_INPUT_FILES, P);
    if (!bridgeFeatureOk(static_cast<int32_t>(h.featureType)))
        exitWithError("EXITING: bridge snapshot featureType not supported for replay\n", std::cerr, P.inOut->logMain,
                      EXIT_CODE_INPUT_FILES, P);
    if (static_cast<int32_t>(h.featureType) != featureType)
        exitWithError("EXITING: bridge snapshot featureType does not match current Solo feature\n", std::cerr,
                      P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    if (h.cbWLsize != pSolo.cbWLsize)
        exitWithError("EXITING: bridge snapshot cbWLsize mismatch\n", std::cerr, P.inOut->logMain,
                      EXIT_CODE_INPUT_FILES, P);
    if (h.runThreadN != static_cast<uint32_t>(P.runThreadN))
        exitWithError("EXITING: bridge snapshot runThreadN mismatch (re-run with same --runThreadN as seed)\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);

    g_statsAll.readN = h.statsReadN;
    nReadsInput = h.statsReadN + 1;

    auto resetHash = [](SoloReadFeature *rf) {
        if (rf == nullptr)
            return;
        if (rf->inlineHash_) {
            kh_destroy(cg_agg, rf->inlineHash_);
            rf->inlineHash_ = nullptr;
        }
        rf->inlineHash_ = kh_init(cg_agg);
    };

    for (int ii = 0; ii < P.runThreadN; ++ii) {
        if (readFeatAll[ii] == nullptr)
            exitWithError("EXITING: bridge snapshot replay: readFeatAll[] is null\n", std::cerr, P.inOut->logMain,
                          EXIT_CODE_INPUT_FILES, P);
        resetHash(readFeatAll[ii]);
    }
    resetHash(readFeatSum);

    for (int ii = 0; ii < P.runThreadN; ++ii)
        readHashSection(in, readFeatAll[ii]->inlineHash_, P);
    readHashSection(in, readFeatSum->inlineHash_, P);

    auto elapsed = [&t0]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };
    P.inOut->logMain << "  snapshot hash sections loaded in " << elapsed() << " s" << endl;

    // bridgeImmediateReadCounts_
    {
        struct ImmRow { uint32_t cb; uint32_t pad; uint64_t val; };
        std::vector<ImmRow> immBuf(h.nImmediatePairs);
        if (h.nImmediatePairs > 0) {
            if (!in.read(reinterpret_cast<char *>(immBuf.data()),
                         static_cast<std::streamsize>(h.nImmediatePairs * sizeof(ImmRow))))
                exitWithError("EXITING: bridge snapshot truncated (immediate counts)\n", std::cerr, P.inOut->logMain,
                              EXIT_CODE_INPUT_FILES, P);
        }
        readFeatSum->bridgeImmediateReadCounts_.clear();
        readFeatSum->bridgeImmediateReadCounts_.reserve(h.nImmediatePairs);
        for (uint32_t i = 0; i < h.nImmediatePairs; ++i)
            readFeatSum->bridgeImmediateReadCounts_[immBuf[i].cb] = immBuf[i].val;
    }

    // Pin vectors
    readFeatSum->bridgePinNreadUnique_.assign(h.pinVecLen, 0u);
    readFeatSum->bridgePinNreadMulti_.assign(h.pinVecLen, 0u);
    if (h.pinVecLen > 0) {
        if (!in.read(reinterpret_cast<char *>(readFeatSum->bridgePinNreadUnique_.data()),
                     static_cast<std::streamsize>(h.pinVecLen * sizeof(uint32_t))))
            exitWithError("EXITING: bridge snapshot truncated (pin unique)\n", std::cerr, P.inOut->logMain,
                          EXIT_CODE_INPUT_FILES, P);
        if (!in.read(reinterpret_cast<char *>(readFeatSum->bridgePinNreadMulti_.data()),
                     static_cast<std::streamsize>(h.pinVecLen * sizeof(uint32_t))))
            exitWithError("EXITING: bridge snapshot truncated (pin multi)\n", std::cerr, P.inOut->logMain,
                          EXIT_CODE_INPUT_FILES, P);
    }

    // flagCountsNoCB
    if (!in.read(reinterpret_cast<char *>(readFeatSum->readFlag.flagCountsNoCB.data()),
                 static_cast<std::streamsize>(SoloReadFlagClass::nBits * sizeof(uint64_t))))
        exitWithError("EXITING: bridge snapshot truncated (flagCountsNoCB)\n", std::cerr, P.inOut->logMain,
                      EXIT_CODE_INPUT_FILES, P);

    // flagCounts per-CB
    readFeatSum->readFlag.flagCounts.clear();
    {
        constexpr size_t rowBytes = sizeof(uint64_t) + SoloReadFlagClass::nBits * sizeof(uint64_t);
        std::vector<char> rowBuf(rowBytes);
        for (uint32_t i = 0; i < h.nFlagCbEntries; ++i) {
            if (!in.read(rowBuf.data(), static_cast<std::streamsize>(rowBytes)))
                exitWithError("EXITING: bridge snapshot truncated (flagCounts)\n", std::cerr, P.inOut->logMain,
                              EXIT_CODE_INPUT_FILES, P);
            uint64_t cb = 0;
            std::memcpy(&cb, rowBuf.data(), sizeof(uint64_t));
            auto &dst = readFeatSum->readFlag.flagCounts[cb];
            std::memcpy(dst.data(), rowBuf.data() + sizeof(uint64_t),
                        SoloReadFlagClass::nBits * sizeof(uint64_t));
        }
    }

    // stats
    if (!in.read(reinterpret_cast<char *>(readFeatSum->stats.V),
                 static_cast<std::streamsize>(SoloReadFeatureStats::nStats * sizeof(uint64_t))))
        exitWithError("EXITING: bridge snapshot truncated (stats)\n", std::cerr, P.inOut->logMain,
                      EXIT_CODE_INPUT_FILES, P);

    // cbReadCount
    if (pSolo.cbWLyes) {
        readFeatSum->cbReadCount.resize(pSolo.cbWLsize);
        if (!in.read(reinterpret_cast<char *>(readFeatSum->cbReadCount.data()),
                     static_cast<std::streamsize>(pSolo.cbWLsize * sizeof(uint32_t))))
            exitWithError("EXITING: bridge snapshot truncated (cbReadCount)\n", std::cerr, P.inOut->logMain,
                          EXIT_CODE_INPUT_FILES, P);
    }

    if (in.peek() != std::char_traits<char>::eof())
        exitWithError("EXITING: bridge snapshot has trailing garbage after expected payload\n", std::cerr,
                      P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);

    uint64_t verifyTotal = 0;
    for (int ii = 0; ii < P.runThreadN; ++ii)
        verifyTotal += liveHashSize(readFeatAll[ii]->inlineHash_);
    verifyTotal += liveHashSize(readFeatSum->inlineHash_);
    if (verifyTotal != h.totalHashEntries)
        exitWithError("EXITING: bridge snapshot hash entry count mismatch after load\n", std::cerr, P.inOut->logMain,
                      EXIT_CODE_INPUT_FILES, P);

    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... STAR_SOLO_BRIDGE_HASH_SNAPSHOT_IN replay: loaded "
                     << verifyTotal << " hash entries from " << path << " in " << elapsed() << " s" << endl;
}
