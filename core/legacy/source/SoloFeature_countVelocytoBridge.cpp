#include "SoloFeature.h"
#include "SoloMemoryProfile.h"
#include "SoloReadFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "SoloCommon.h"
#include "systemFunctions.h"
#include "ErrorWarning.h"
#include "IncludeDefine.h"
#include <algorithm>
#include <bitset>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <memory>
#include <sys/stat.h>
#include <unordered_map>
#include <vector>
#include <unistd.h>

namespace {

static int velocytoReferenceFeatureIndex(const ParametersSolo &pSolo)
{
    const int candidates[] = {
        SoloFeatureTypes::Gene,
        SoloFeatureTypes::GeneFull,
        SoloFeatureTypes::GeneFull_Ex50pAS,
        SoloFeatureTypes::GeneFull_ExonOverIntron,
    };
    for (auto featureType : candidates) {
        if (pSolo.featureInd[featureType] >= 0)
            return pSolo.featureInd[featureType];
    }
    return -1;
}

static SoloFeature *velocytoReferenceFeature(const ParametersSolo &pSolo,
                                             SoloFeature **soloFeatAll,
                                             Parameters &P)
{
    const int refIndex = velocytoReferenceFeatureIndex(pSolo);
    if (refIndex < 0) {
        exitWithError("EXITING because of fatal PARAMETERS error: Velocyto requires a gene-like Solo feature "
                          "(Gene, GeneFull, GeneFull_Ex50pAS, or GeneFull_ExonOverIntron) for CB/UMI readInfo.\n"
                          "SOLUTION: re-run STAR adding GeneFull or another gene-like feature to --soloFeatures.\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    return soloFeatAll[refIndex];
}

struct VelocytoSortedRecord {
    uint64 iread;
    uint32 iCB;
    uintUMI umi;
    std::vector<trTypeStruct> trT;
};

struct VelocytoBucketRecord {
    uint64 iread;
    uintUMI umi;
    std::vector<trTypeStruct> trT;
};

/** Merge one stream record into per-(CB,UMI) transcript vectors (intersection + mask OR). */
inline void applyVelocytoMerge(std::unordered_map<uintUMI, std::vector<trTypeStruct>> &cuMap,
                               uintUMI umi,
                               const std::vector<trTypeStruct> &trT)
{
    auto it = cuMap.find(umi);
    if (it != cuMap.end() && it->second.empty())
        return;
    if (it == cuMap.end()) {
        cuMap[umi] = trT;
        return;
    }
    const std::vector<trTypeStruct> &oldV = it->second;
    uint32 inew = 0;
    std::vector<trTypeStruct> trT1;
    trT1.reserve(oldV.size());
    for (uint32 iold = 0; iold < oldV.size(); iold++) {
        while (inew < trT.size() && oldV[iold].tr > trT[inew].tr)
            ++inew;
        if (inew == trT.size())
            break;
        if (oldV[iold].tr == trT[inew].tr) {
            trT1.push_back({trT[inew].tr, (uint8)(oldV[iold].type | trT[inew].type)});
        }
    }
    cuMap[umi] = std::move(trT1);
}

template <class T>
inline bool velocytoReadBin(std::istream &in, T &v) {
    return (bool)in.read(reinterpret_cast<char *>(&v), sizeof(T));
}

inline void velocytoWriteBin(std::ostream &out, const void *p, size_t n) {
    out.write(reinterpret_cast<const char *>(p), n);
}

inline void velocytoWriteU32(std::ostream &out, uint32_t v) {
    velocytoWriteBin(out, &v, sizeof(v));
}

inline bool velocytoSpillWriteBin(std::ostream &out, const void *p, size_t n) {
    out.write(reinterpret_cast<const char *>(p), n);
    return (bool)out;
}

inline bool velocytoSpillWriteU32(std::ostream &out, uint32_t v) {
    return velocytoSpillWriteBin(out, &v, sizeof(v));
}

static void velocytoSpillFatal(const std::string &msg, Parameters &P, int exitCode) {
    exitWithError("EXITING fatal error: Velocyto integrated-hash spill: " + msg, std::cerr, P.inOut->logMain, exitCode, P);
}

/** Read exactly sizeof(T) or exit (no silent truncate). */
template <class T>
static void velocytoSpillReadFullDie(std::istream &in, T &v, const std::string &path, const char *ctx, Parameters &P) {
    in.read(reinterpret_cast<char *>(&v), sizeof(T));
    if (in.gcount() != static_cast<std::streamsize>(sizeof(T)))
        velocytoSpillFatal("truncated reading " + std::string(ctx) + " from " + path + " (expected " + std::to_string(sizeof(T))
                               + " bytes, got " + std::to_string((unsigned long long)in.gcount()) + ")",
                           P, EXIT_CODE_INCONSISTENT_DATA);
}

static uint32_t velocytoIntegratedSpillBucketCount() {
    const char *e = std::getenv("STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS");
    if (e == nullptr || e[0] == '\0') {
        const char *lowMem = std::getenv("STAR_VELOCYTO_LOW_MEM");
        return (lowMem != nullptr && lowMem[0] != '\0' && std::strcmp(lowMem, "0") != 0) ? 4096u : 128u;
    }
    char *end = nullptr;
    unsigned long x = std::strtoul(e, &end, 10);
    if (end == e || x < 1ul)
        return 128u;
    if (x > 16384ul)
        return 16384u;
    return (uint32_t)x;
}

static bool velocytoIntegratedHashInMemory() {
    const char *e = std::getenv("STAR_VELOCYTO_INTEGRATED_HASH_INMEMORY");
    return e != nullptr && e[0] != '\0' && std::strcmp(e, "0") != 0;
}

static bool velocytoLowMemMode() {
    const char *e = std::getenv("STAR_VELOCYTO_LOW_MEM");
    return e != nullptr && e[0] != '\0' && std::strcmp(e, "0") != 0;
}

static uint32_t velocytoInitialMapReserveCap()
{
    const char *e = std::getenv("STAR_VELOCYTO_UMI_RESERVE_CAP");
    if (e == nullptr || e[0] == '\0')
        return velocytoLowMemMode() ? 64u : 2048u;
    char *end = nullptr;
    unsigned long x = std::strtoul(e, &end, 10);
    if (end == e)
        return 2048u;
    if (x > 1000000ul)
        return 1000000u;
    return static_cast<uint32_t>(x);
}

static uint32_t velocytoInitialMapReserve(uint32_t readCount)
{
    if (readCount < 32u)
        return 0u;
    uint32_t estimate = readCount / 16u;
    if (estimate < 8u)
        estimate = 8u;
    const uint32_t cap = velocytoInitialMapReserveCap();
    if (estimate > cap)
        estimate = cap;
    return estimate;
}

static void reserveVelocytoCuMaps(vector<unordered_map<uintUMI, vector<trTypeStruct>>> &cuTrTypes,
                                  const vector<uint32> &cbReadCount,
                                  uint32 nCB,
                                  Parameters &P,
                                  const char *label)
{
    uint64_t totalReserve = 0;
    uint32_t reservedCBs = 0;
    uint32_t maxReserve = 0;
    for (uint32 ii = 0; ii < nCB; ii++) {
        const uint32_t reserveN = (ii < cbReadCount.size()) ? velocytoInitialMapReserve(cbReadCount[ii]) : 0u;
        if (reserveN == 0u)
            continue;
        cuTrTypes[ii].reserve(reserveN);
        totalReserve += reserveN;
        ++reservedCBs;
        if (reserveN > maxReserve)
            maxReserve = reserveN;
    }
    P.inOut->logMain << "Velocyto UMI map reserve policy: " << label
                     << " reserved_cbs=" << reservedCBs
                     << " total_initial_buckets=" << totalReserve
                     << " max_per_cb=" << maxReserve
                     << " cap=" << velocytoInitialMapReserveCap()
                     << endl;
}

static void reserveVelocytoCuMapsRange(vector<unordered_map<uintUMI, vector<trTypeStruct>>> &cuTrTypes,
                                       const vector<uint32> &cbReadCount,
                                       uint32 cbBegin,
                                       uint32 cbEnd)
{
    for (uint32 iCB = cbBegin; iCB < cbEnd; iCB++) {
        const uint32_t reserveN = (iCB < cbReadCount.size()) ? velocytoInitialMapReserve(cbReadCount[iCB]) : 0u;
        if (reserveN == 0u)
            continue;
        cuTrTypes[iCB - cbBegin].reserve(reserveN);
    }
}

static uint32_t velocytoBucketForCb(uint32 iCB, uint32 nCB, uint32_t nBuckets)
{
    if (nBuckets <= 1u || nCB == 0u)
        return 0u;
    const uint64_t bucketWidth = (static_cast<uint64_t>(nCB) + nBuckets - 1u) / nBuckets;
    uint32_t b = static_cast<uint32_t>(static_cast<uint64_t>(iCB) / bucketWidth);
    if (b >= nBuckets)
        b = nBuckets - 1u;
    return b;
}

static void velocytoBucketRange(uint32_t bucket,
                                uint32 nCB,
                                uint32_t nBuckets,
                                uint32 &cbBegin,
                                uint32 &cbEnd)
{
    if (nBuckets <= 1u || nCB == 0u) {
        cbBegin = 0u;
        cbEnd = nCB;
        return;
    }
    const uint64_t bucketWidth = (static_cast<uint64_t>(nCB) + nBuckets - 1u) / nBuckets;
    const uint64_t begin = static_cast<uint64_t>(bucket) * bucketWidth;
    uint64_t end = begin + bucketWidth;
    if (begin > nCB)
        cbBegin = nCB;
    else
        cbBegin = static_cast<uint32>(begin);
    if (end > nCB)
        cbEnd = nCB;
    else
        cbEnd = static_cast<uint32>(end);
}

static size_t velocytoInitialCountWords(uint64 nReadsMapped, uint32 countMatStride)
{
    const uint64 defaultCap = 100000000ull;
    uint64 words = nReadsMapped * static_cast<uint64>(countMatStride) / 20ull + 16ull;
    if (words < 1048576ull)
        words = 1048576ull;
    if (words > defaultCap)
        words = defaultCap;
    return static_cast<size_t>(words);
}

} // namespace

void SoloFeature::countVelocytoStreamThreads()
{
    time_t rawTime;
    SoloFeature *readInfoSource = velocytoReferenceFeature(pSolo, soloFeatAll, P);

    nReadPerCB.resize(nCB);

    vector<unordered_map<uintUMI, vector<trTypeStruct>>> cuTrTypes(nCB);
    reserveVelocytoCuMaps(cuTrTypes, readFeatSum->cbReadCount, nCB, P, "stream_threads");

    {
        std::ostringstream extra;
        extra << "velocyto_cuTrTypes_cbs=" << nCB
              << " packedReadInfo_words=" << readInfoSource->packedReadInfo.data.size();
        soloMemoryProfileCheckpoint(P.inOut->logMain, "countVelocytoStreamThreads_maps_allocated", extra.str());
    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting: allocated arrays" << endl;

    for (int iThread = 0; iThread < P.runThreadN; iThread++) {
        fstream *streamReads = readFeatAll[iThread]->streamReads;
        streamReads->flush();
        streamReads->seekg(0, ios::beg);

        uint64 iread;
        while (*streamReads >> iread) {
            uint32_t cb = readInfoSource->getPackedCB((uint32_t)iread);
            uint32_t umi = readInfoSource->getPackedUMI((uint32_t)iread);
            uint8_t status = readInfoSource->getPackedStatus((uint32_t)iread);
            if (status != 1) {
                streamReads->ignore((uint32)-1, '\n');
                continue;
            }

            uint32 iCB = indCBwl[cb];
            nReadPerCB[iCB]++;

            if (cuTrTypes[iCB].count(umi) > 0 && cuTrTypes[iCB][umi].empty()) {
                streamReads->ignore((uint32)-1, '\n');
                continue;
            }

            uint32 nTr;
            *streamReads >> nTr;
            vector<trTypeStruct> trT(nTr);
            for (auto &tt : trT) {
                uint32 ty;
                *streamReads >> tt.tr >> ty;
                tt.type = (uint8)ty;
            }

            applyVelocytoMerge(cuTrTypes[iCB], umi, trT);
        }
    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting: finished input" << endl;

    {
        size_t umiBuckets = 0;
        size_t trEntries = 0;
        for (uint32 ii = 0; ii < nCB; ii++) {
            umiBuckets += cuTrTypes[ii].size();
            for (const auto &kv : cuTrTypes[ii]) {
                trEntries += kv.second.size();
            }
        }
        std::ostringstream extra;
        extra << "velocyto_umi_buckets=" << umiBuckets
              << " velocyto_tr_entries=" << trEntries;
        soloMemoryProfileCheckpoint(P.inOut->logMain, "countVelocytoStreamThreads_input_done", extra.str());
    }
    countVelocytoFinalizeFromCuMaps(cuTrTypes);
    soloMemoryProfileCheckpoint(P.inOut->logMain, "countVelocyto_finalize_done");
}

void SoloFeature::countVelocytoSortedReplay()
{
    time_t rawTime;
    SoloFeature *readInfoSource = velocytoReferenceFeature(pSolo, soloFeatAll, P);

    nReadPerCB.resize(nCB);

    vector<unordered_map<uintUMI, vector<trTypeStruct>>> cuTrTypes(nCB);
    reserveVelocytoCuMaps(cuTrTypes, readFeatSum->cbReadCount, nCB, P, "sorted_replay");

    // Full in-memory vector of every Velocyto-positive stream record — can dominate RSS on very large
    // inputs; 2M acceptance must budget host memory (harness enforces UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB
    // unless explicitly uncapped for dev). Stage 2 CB-bucket path avoids the global sort; see
    // countVelocytoSortedReplayCBuckets when STAR_VELOCYTO_INTEGRATED_HASH is set.
    vector<VelocytoSortedRecord> records;
    records.reserve((size_t)(nReadsMapped + 16));

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting (sorted replay): scanning streams" << endl;

    for (int iThread = 0; iThread < P.runThreadN; iThread++) {
        fstream *streamReads = readFeatAll[iThread]->streamReads;
        streamReads->flush();
        streamReads->seekg(0, ios::beg);

        uint64 iread;
        while (*streamReads >> iread) {
            uint32_t cb = readInfoSource->getPackedCB((uint32_t)iread);
            uint32_t umi = readInfoSource->getPackedUMI((uint32_t)iread);
            uint8_t status = readInfoSource->getPackedStatus((uint32_t)iread);
            if (status != 1) {
                streamReads->ignore((uint32)-1, '\n');
                continue;
            }

            uint32 iCB = indCBwl[cb];
            nReadPerCB[iCB]++;

            uint32 nTr;
            *streamReads >> nTr;
            vector<trTypeStruct> trT(nTr);
            for (auto &tt : trT) {
                uint32 ty;
                *streamReads >> tt.tr >> ty;
                tt.type = (uint8)ty;
            }

            records.push_back({iread, iCB, umi, std::move(trT)});
        }
    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting (sorted replay): collected " << records.size()
                     << " stream records\n"
                     << "RAM after Velocyto sorted-replay materialization:\n"
                     << linuxProcMemory() << flush;

    std::sort(records.begin(), records.end(), [](const VelocytoSortedRecord &a, const VelocytoSortedRecord &b) {
        if (a.iCB != b.iCB)
            return a.iCB < b.iCB;
        if (a.umi != b.umi)
            return a.umi < b.umi;
        return a.iread < b.iread;
    });

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting (sorted replay): merging " << records.size()
                     << " records\n"
                     << "RAM after Velocyto sorted-replay sort:\n"
                     << linuxProcMemory() << flush;

    for (const auto &rec : records)
        applyVelocytoMerge(cuTrTypes[rec.iCB], rec.umi, rec.trT);

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting (sorted replay): finished input" << endl;

    countVelocytoFinalizeFromCuMaps(cuTrTypes);
}

void SoloFeature::countVelocytoSortedReplayCBuckets()
{
    time_t rawTime;
    SoloFeature *readInfoSource = velocytoReferenceFeature(pSolo, soloFeatAll, P);

    nReadPerCB.resize(nCB);

    if (velocytoIntegratedHashInMemory()) {
        // Debug / A–B only: holds every Velocyto-positive record in RAM (per-CB vectors).
        vector<unordered_map<uintUMI, vector<trTypeStruct>>> cuTrTypes(nCB);
        reserveVelocytoCuMaps(cuTrTypes, readFeatSum->cbReadCount, nCB, P, "integrated_hash_inmemory");
        vector<vector<VelocytoBucketRecord>> perCB(nCB);
        for (uint32 ii = 0; ii < nCB; ii++) {
            uint32 est = readFeatSum->cbReadCount[ii];
            if (est > 0u)
                perCB[ii].reserve((size_t)est / 2u + 8u);
        }

        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime)
                         << " ... Velocyto counting (integrated hash INMEMORY=1): scanning streams" << endl;

        uint64 nRec = 0;
        for (int iThread = 0; iThread < P.runThreadN; iThread++) {
            fstream *streamReads = readFeatAll[iThread]->streamReads;
            streamReads->flush();
            streamReads->seekg(0, ios::beg);

            uint64 iread;
            while (*streamReads >> iread) {
                uint32_t cb = readInfoSource->getPackedCB((uint32_t)iread);
                uint32_t umi = readInfoSource->getPackedUMI((uint32_t)iread);
                uint8_t status = readInfoSource->getPackedStatus((uint32_t)iread);
                if (status != 1) {
                    streamReads->ignore((uint32)-1, '\n');
                    continue;
                }

                uint32 iCB = indCBwl[cb];
                nReadPerCB[iCB]++;

                uint32 nTr;
                *streamReads >> nTr;
                vector<trTypeStruct> trT(nTr);
                for (auto &tt : trT) {
                    uint32 ty;
                    *streamReads >> tt.tr >> ty;
                    tt.type = (uint8)ty;
                }

                perCB[iCB].push_back({iread, umi, std::move(trT)});
                ++nRec;
            }
        }

        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting (integrated hash INMEMORY): collected " << nRec
                         << " stream records\n"
                         << "RAM after Velocyto sorted-replay materialization:\n"
                         << linuxProcMemory() << flush;

        for (uint32 iCB = 0; iCB < nCB; iCB++) {
            vector<VelocytoBucketRecord> &bucket = perCB[iCB];
            if (bucket.empty())
                continue;

            std::sort(bucket.begin(), bucket.end(), [](const VelocytoBucketRecord &a, const VelocytoBucketRecord &b) {
                if (a.umi != b.umi)
                    return a.umi < b.umi;
                return a.iread < b.iread;
            });

            for (const auto &rec : bucket)
                applyVelocytoMerge(cuTrTypes[iCB], rec.umi, rec.trT);

            vector<VelocytoBucketRecord>().swap(bucket);
        }

        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting (integrated hash INMEMORY): finished merge" << endl;

        countVelocytoFinalizeFromCuMaps(cuTrTypes);
        return;
    }

    const uint32_t requestedBuckets = velocytoIntegratedSpillBucketCount();
    const uint32_t nBuckets = nCB == 0u ? 1u : std::min<uint32_t>(requestedBuckets, nCB);
    std::string tmpl = outputPrefix + "velocyto_integrated_spill_XXXXXX";
    std::vector<char> tpath(tmpl.begin(), tmpl.end());
    tpath.push_back('\0');
    if (mkdtemp(tpath.data()) == nullptr) {
        exitWithError("EXITING fatal error: Velocyto integrated-hash mkdtemp failed for " + tmpl, std::cerr, P.inOut->logMain,
                      EXIT_CODE_FILE_OPEN, P);
    }
    const std::string spillDir(tpath.data());

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ... Velocyto counting (integrated hash / range disk spill): scanning streams, spill_buckets=" << nBuckets
                     << " dir=" << spillDir << endl;

    std::vector<std::unique_ptr<std::ofstream>> spillOut;
    spillOut.reserve(nBuckets);
    for (uint32_t b = 0; b < nBuckets; b++) {
        std::string p = spillDir + "/bucket_" + std::to_string((unsigned long)b) + ".bin";
        spillOut.emplace_back(new std::ofstream(p, std::ios::binary | std::ios::out | std::ios::trunc));
        if (!*spillOut.back()) {
            exitWithError("EXITING fatal error: cannot open Velocyto spill file " + p, std::cerr, P.inOut->logMain, EXIT_CODE_FILE_OPEN, P);
        }
    }

    uint64 nRec = 0;
    for (int iThread = 0; iThread < P.runThreadN; iThread++) {
        fstream *streamReads = readFeatAll[iThread]->streamReads;
        streamReads->flush();
        streamReads->seekg(0, ios::beg);

        uint64 iread;
        while (*streamReads >> iread) {
            uint32_t cb = readInfoSource->getPackedCB((uint32_t)iread);
            uint32_t umi = readInfoSource->getPackedUMI((uint32_t)iread);
            uint8_t status = readInfoSource->getPackedStatus((uint32_t)iread);
            if (status != 1) {
                streamReads->ignore((uint32)-1, '\n');
                continue;
            }

            uint32 iCB = indCBwl[cb];
            nReadPerCB[iCB]++;

            uint32 nTr;
            *streamReads >> nTr;
            vector<trTypeStruct> trT(nTr);
            for (auto &tt : trT) {
                uint32 ty;
                *streamReads >> tt.tr >> ty;
                tt.type = (uint8)ty;
            }

            const uint32_t shard = velocytoBucketForCb(iCB, nCB, nBuckets);
            std::ostream &os = *spillOut[shard];
            if (!velocytoSpillWriteU32(os, iCB) || !velocytoSpillWriteBin(os, &iread, sizeof(iread)) || !velocytoSpillWriteU32(os, umi)
                || !velocytoSpillWriteU32(os, nTr)) {
                velocytoSpillFatal("write failed (disk full?) bucket=" + std::to_string((unsigned long)shard), P, EXIT_CODE_FILE_WRITE);
            }
            for (auto &tt : trT) {
                if (!velocytoSpillWriteU32(os, tt.tr) || !velocytoSpillWriteU32(os, (uint32_t)tt.type))
                    velocytoSpillFatal("write failed (disk full?) bucket=" + std::to_string((unsigned long)shard), P, EXIT_CODE_FILE_WRITE);
            }
            ++nRec;
        }
    }

    for (uint32_t b = 0; b < nBuckets; b++) {
        spillOut[b]->flush();
        if (!*spillOut[b])
            velocytoSpillFatal("flush failed bucket=" + std::to_string((unsigned long)b), P, EXIT_CODE_FILE_WRITE);
    }
    spillOut.clear();

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting (integrated hash): staged " << nRec
                     << " records to spill\n"
                     << "RAM after Velocyto integrated-hash spill (all records staged to disk):\n"
                     << linuxProcMemory() << flush;

    countVelocytoFinalizeInit();

    for (uint32_t b = 0; b < nBuckets; b++) {
        uint32 cbBegin = 0u;
        uint32 cbEnd = 0u;
        velocytoBucketRange(b, nCB, nBuckets, cbBegin, cbEnd);
        if (cbBegin >= cbEnd) {
            std::string p = spillDir + "/bucket_" + std::to_string((unsigned long)b) + ".bin";
            std::remove(p.c_str());
            continue;
        }

        std::string p = spillDir + "/bucket_" + std::to_string((unsigned long)b) + ".bin";
        struct stat st {};
        if (::stat(p.c_str(), &st) != 0) {
            unordered_map<uintUMI, vector<trTypeStruct>> empty;
            for (uint32 iCB = cbBegin; iCB < cbEnd; iCB++)
                countVelocytoFinalizeOneCb(iCB, empty);
            continue;
        }
        if (st.st_size == 0) {
            std::remove(p.c_str());
            unordered_map<uintUMI, vector<trTypeStruct>> empty;
            for (uint32 iCB = cbBegin; iCB < cbEnd; iCB++)
                countVelocytoFinalizeOneCb(iCB, empty);
            continue;
        }

        std::ifstream in(p, std::ios::binary);
        if (!in.is_open()) {
            exitWithError("EXITING fatal error: cannot read Velocyto spill file " + p, std::cerr, P.inOut->logMain, EXIT_CODE_FILE_OPEN, P);
        }

        vector<VelocytoSortedRecord> chunk;
        chunk.reserve(4096);

        for (;;) {
            uint32_t iCB = 0;
            in.read(reinterpret_cast<char *>(&iCB), sizeof(iCB));
            if (in.gcount() == 0) {
                if (in.eof())
                    break;
                velocytoSpillFatal("read failed at record boundary (non-EOF) in " + p, P, EXIT_CODE_INCONSISTENT_DATA);
            }
            if (in.gcount() != static_cast<std::streamsize>(sizeof(iCB)))
                velocytoSpillFatal("partial iCB header in " + p, P, EXIT_CODE_INCONSISTENT_DATA);

            uint64_t iread = 0;
            velocytoSpillReadFullDie(in, iread, p, "iread", P);
            uint32_t umi = 0;
            velocytoSpillReadFullDie(in, umi, p, "umi", P);
            uint32_t nTr = 0;
            velocytoSpillReadFullDie(in, nTr, p, "nTr", P);
            if (nTr > 10000000u)
                velocytoSpillFatal("absurd nTr in " + p, P, EXIT_CODE_INCONSISTENT_DATA);

            vector<trTypeStruct> trT(nTr);
            for (uint32_t i = 0; i < nTr; i++) {
                uint32_t tr = 0, ty = 0;
                velocytoSpillReadFullDie(in, tr, p, "transcript", P);
                velocytoSpillReadFullDie(in, ty, p, "type", P);
                trT[i].tr = tr;
                trT[i].type = (uint8_t)(ty & 0xFFu);
            }
            chunk.push_back({iread, iCB, (uintUMI)umi, std::move(trT)});
        }

        char trailing = 0;
        if (in.get(trailing))
            velocytoSpillFatal("trailing data after last record in " + p, P, EXIT_CODE_INCONSISTENT_DATA);
        if (in.bad())
            velocytoSpillFatal("I/O error reading " + p, P, EXIT_CODE_INCONSISTENT_DATA);

        in.close();
        std::remove(p.c_str());

        std::sort(chunk.begin(), chunk.end(), [](const VelocytoSortedRecord &a, const VelocytoSortedRecord &b) {
            if (a.iCB != b.iCB)
                return a.iCB < b.iCB;
            if (a.umi != b.umi)
                return a.umi < b.umi;
            return a.iread < b.iread;
        });

        vector<unordered_map<uintUMI, vector<trTypeStruct>>> cuTrTypes(cbEnd - cbBegin);
        reserveVelocytoCuMapsRange(cuTrTypes, readFeatSum->cbReadCount, cbBegin, cbEnd);
        for (const auto &rec : chunk)
            applyVelocytoMerge(cuTrTypes[rec.iCB - cbBegin], rec.umi, rec.trT);

        vector<VelocytoSortedRecord>().swap(chunk);
        for (uint32 iCB = cbBegin; iCB < cbEnd; iCB++)
            countVelocytoFinalizeOneCb(iCB, cuTrTypes[iCB - cbBegin]);
        vector<unordered_map<uintUMI, vector<trTypeStruct>>>().swap(cuTrTypes);
    }

    if (::rmdir(spillDir.c_str()) != 0) {
        time_t warnTime;
        time(&warnTime);
        P.inOut->logMain << timeMonthDayTime(warnTime) << " ... WARNING: could not rmdir Velocyto spill dir " << spillDir << endl;
    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting (integrated hash spill): finished merge" << endl;

    countVelocytoFinalizeFinish();
}

void SoloFeature::countVelocytoFinalizeInit()
{
    nUMIperCB.assign(nCB, 0);
    nGenePerCB.assign(nCB, 0);

    countMatStride = 4;
    countCellGeneUMI.clear();
    countCellGeneUMI.reserve(velocytoInitialCountWords(nReadsMapped, countMatStride));
    countCellGeneUMIindex.resize(nCB + 1);
    countCellGeneUMIindex[0] = 0;
}

void SoloFeature::countVelocytoFinalizeOneCb(uint32 iCB, unordered_map<uintUMI, vector<trTypeStruct>> &cuMap)
{
    map<uint32, array<uint32, 3>> geneC;
    for (auto &umi : cuMap) {
        if (umi.second.empty())
            continue;
        uint32 geneI = Trans.trGene[umi.second[0].tr];

        bool exonModel = false;
        bool intronModel = false;
        bool spanModel = true;
        bool mixedModel = false;

        for (auto &tt : umi.second) {
            if (Trans.trGene[tt.tr] != geneI) {
                geneI = (uint32)-1;
                break;
            }

            bitset<velocytoTypeGeneBits> gV(tt.type);

            mixedModel |= ((gV.test(AlignVsTranscript::Intron) && gV.test(AlignVsTranscript::Concordant))
                           || gV.test(AlignVsTranscript::ExonIntron))
                          && !gV.test(AlignVsTranscript::ExonIntronSpan);
            spanModel &= gV.test(AlignVsTranscript::ExonIntronSpan);
            exonModel |= gV.test(AlignVsTranscript::Concordant) && !gV.test(AlignVsTranscript::Intron)
                         && !gV.test(AlignVsTranscript::ExonIntron);
            intronModel |= gV.test(AlignVsTranscript::Intron) && !gV.test(AlignVsTranscript::ExonIntron)
                           && !gV.test(AlignVsTranscript::Concordant);
        }

        if (geneI + 1 == 0)
            continue;

        if (exonModel && !intronModel && !mixedModel) {
            geneC[geneI][0]++;
        } else if (spanModel || ((intronModel || mixedModel) && !exonModel)) {
            geneC[geneI][1]++;
        } else {
            geneC[geneI][2]++;
        }

        nUMIperCB[iCB]++;
    }

    countCellGeneUMIindex[iCB + 1] = countCellGeneUMIindex[iCB];

    if (nUMIperCB[iCB] == 0) {
        unordered_map<uintUMI, vector<trTypeStruct>>().swap(cuMap);
        return;
    }

    nGenePerCB[iCB] += geneC.size();
    readFeatSum->stats.V[readFeatSum->stats.yesUMIs] += nUMIperCB[iCB];
    ++readFeatSum->stats.V[readFeatSum->stats.yesCellBarcodes];

    for (auto &gg : geneC) {
        countCellGeneUMI.push_back(gg.first);
        for (uint32 ii = 0; ii < 3; ii++)
            countCellGeneUMI.push_back(gg.second[ii]);
        countCellGeneUMIindex[iCB + 1] += countMatStride;
    }

    unordered_map<uintUMI, vector<trTypeStruct>>().swap(cuMap);
}

void SoloFeature::countVelocytoFinalizeFinish()
{
    time_t rawTime;

    nReadPerCBtotal = nReadPerCB;
    nReadPerCBunique = nReadPerCB;

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Velocyto counting: finished collapsing UMIs" << endl;
    P.inOut->logMain << "RAM for solo feature " << SoloFeatureTypes::Names[featureType] << "\n"
                     << linuxProcMemory() << flush;
}

void SoloFeature::countVelocytoFinalizeFromCuMaps(vector<unordered_map<uintUMI, vector<trTypeStruct>>> &cuTrTypes)
{
    countVelocytoFinalizeInit();
    for (uint32 iCB = 0; iCB < nCB; iCB++) {
        countVelocytoFinalizeOneCb(iCB, cuTrTypes[iCB]);
    }
    countVelocytoFinalizeFinish();
}
