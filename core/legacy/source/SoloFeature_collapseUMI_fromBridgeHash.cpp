// Non-Flex Solo inline-hash bridge: collapse UMIs directly from aggregated bridge hash
// (tagless key [wlCb24][UMI24][geneFull16]) without materializeRGUFromHash() / rGeneUMI.
//
// Data flow (no global sort over the full tuple set):
// 1) Canonicalize thread-local + merged hashes, count entries per whitelist CB, then CSR-fill
//    a flat (gene, umi24, count) array with per-CB slices [cbOffsets[i], cbOffsets[i+1]).
// 2) Per CB: count/scatter into a local buffer grouped by gene (no comparison sort of the slice);
//    per gene bucket run 1MM_CR via umiArrayCorrect_CR (which sorts UMIs internally).
// 3) MultiGeneUMI_CR: collect flat (corrUmi, origUmi, geneIdx, count) rows for the CB, radix-sort
//    by corrected UMI, resolve each corrected-UMI group with vector scans (no nested umi->gene maps).

#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "SoloCommon.h"
#include "ParametersSolo.h"
#include "SoloHDMoleculeWriter.h"
#include "hash_shims_cpp_compat.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "ErrorWarning.h"
#include "IncludeDefine.h"
#include "MultiGeneUmiCr.h"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <omp.h>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

KHASH_MAP_INIT_INT(cbcount, size_t)

namespace {

const std::unordered_set<std::string>& collapseTraceBarcodeSet()
{
    static std::unordered_set<std::string> barcodes;
    static bool loaded = false;
    if (loaded) {
        return barcodes;
    }
    loaded = true;

    const char *path = std::getenv("STAR_SOLO_DEBUG_BARCODE_FILE");
    if (path == nullptr || path[0] == '\0') {
        return barcodes;
    }

    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty()) {
            barcodes.insert(line);
        }
    }
    return barcodes;
}

bool shouldTraceCollapseBarcode(const ParametersSolo &pSolo, uint32_t wlIdx)
{
    if (wlIdx >= pSolo.cbWLstr.size()) {
        return false;
    }
    const auto &debugSet = collapseTraceBarcodeSet();
    if (debugSet.empty()) {
        return false;
    }
    return debugSet.count(pSolo.cbWLstr[wlIdx]) != 0;
}

struct BridgeStageDigest {
    uint64_t records = 0;
    uint64_t totalCount = 0;
    uint64_t hashSum = 0;
    uint64_t hashXor = 0;
    uint64_t hashSum2 = 0;

    void add(uint64_t h, uint32_t count)
    {
        ++records;
        totalCount += count;
        hashSum += h;
        hashXor ^= h;
        hashSum2 += h * h + 0x9e3779b97f4a7c15ull;
    }

    void merge(const BridgeStageDigest &other)
    {
        records += other.records;
        totalCount += other.totalCount;
        hashSum += other.hashSum;
        hashXor ^= other.hashXor;
        hashSum2 += other.hashSum2;
    }
};

uint64_t bridgeSplitMix64(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ull;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
    return x ^ (x >> 31);
}

uint64_t bridgeTupleHash(uint64_t stageTag, uint32_t wlCb, uint32_t gene, uint32_t umi,
                         uint32_t corr, uint32_t count)
{
    uint64_t h = bridgeSplitMix64(stageTag);
    h ^= bridgeSplitMix64(static_cast<uint64_t>(wlCb) | (static_cast<uint64_t>(gene) << 32));
    h ^= bridgeSplitMix64(static_cast<uint64_t>(umi) | (static_cast<uint64_t>(corr) << 32));
    h ^= bridgeSplitMix64(static_cast<uint64_t>(count) | 0xd1b54a32d192ed03ull);
    return bridgeSplitMix64(h);
}

std::string formatBridgeGeneCounts(const std::vector<std::pair<uint32_t, uint32_t>> &counts,
                                   const std::vector<uint32_t> &gID)
{
    std::vector<std::pair<uint32_t, uint32_t>> ordered;
    ordered.reserve(counts.size());
    for (const auto &kv : counts) {
        const uint32_t localIdx = kv.first;
        const uint32_t geneId = localIdx < gID.size() ? gID[localIdx] : localIdx;
        ordered.push_back({geneId, kv.second});
    }
    std::sort(ordered.begin(), ordered.end());

    std::ostringstream oss;
    for (size_t i = 0; i < ordered.size(); ++i) {
        if (i != 0) {
            oss << ',';
        }
        oss << ordered[i].first << ':' << ordered[i].second;
    }
    return oss.str();
}

struct BridgeFlatSlot {
    uint32_t gene;
    uint32_t umi24;
    uint32_t count;
};

struct MgRow {
    uint32_t corr;
    uint32_t orig;
    uint32_t geneIdx;
    uint32_t count;
};

size_t bridgeHashLiveSize(khash_t(cg_agg) *hash)
{
    return hash == nullptr ? 0 : kh_size(hash);
}

void mergeBridgeHashInto(khash_t(cg_agg) *src, khash_t(cg_agg) *dst)
{
    if (src == nullptr || dst == nullptr) {
        return;
    }
    for (khiter_t it = kh_begin(src); it != kh_end(src); ++it) {
        if (!kh_exist(src, it)) {
            continue;
        }
        int absent = 0;
        khiter_t dest = kh_put(cg_agg, dst, kh_key(src, it), &absent);
        if (absent) {
            kh_val(dst, dest) = kh_val(src, it);
        } else {
            kh_val(dst, dest) += kh_val(src, it);
        }
    }
}

// Stable LSB radix on MgRow.corr (uint32). Replaces std::sort for large per-CB mgBuf.
void radixSortMgRowsByCorr(std::vector<MgRow> &rows)
{
    const size_t n = rows.size();
    if (n <= 1)
        return;
    std::vector<MgRow> tmp(n);
    bool srcIsRows = true;
    for (int pass = 0; pass < 4; ++pass) {
        const int shift = pass * 8;
        uint32_t cnt[256] = {};
        const MgRow *src = srcIsRows ? rows.data() : tmp.data();
        MgRow *dst = srcIsRows ? tmp.data() : rows.data();
        for (size_t i = 0; i < n; ++i) {
            const uint32_t b = (src[i].corr >> shift) & 0xFFu;
            cnt[b]++;
        }
        uint32_t pos[256];
        pos[0] = 0;
        for (int i = 1; i < 256; ++i)
            pos[i] = pos[i - 1] + cnt[i - 1];
        for (size_t i = 0; i < n; ++i) {
            const uint32_t b = (src[i].corr >> shift) & 0xFFu;
            dst[pos[b]++] = src[i];
        }
        srcIsRows = !srcIsRows;
    }
    if (!srcIsRows)
        rows.swap(tmp);
}

void writeBridgeDeterminismDigest(const std::string &path,
                                  const ParametersSolo &pSolo,
                                  uint32_t nCB,
                                  const std::vector<uint32_t> &indCB,
                                  size_t totalHashSize,
                                  size_t mergedAmbigHashSize,
                                  const BridgeStageDigest &preCr,
                                  const BridgeStageDigest &postCr,
                                  const BridgeStageDigest &resolved,
                                  const BridgeStageDigest &finalMatrix,
                                  const std::vector<BridgeStageDigest> *cbPreCr,
                                  const std::vector<BridgeStageDigest> *cbPostCr,
                                  const std::vector<BridgeStageDigest> *cbResolved,
                                  const std::vector<BridgeStageDigest> *cbFinalMatrix,
                                  Parameters &P)
{
    std::ofstream out(path.c_str());
    if (!out.good()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal OUTPUT FILE error: could not open STAR_SOLO_BRIDGE_DETERMINISM_OUT="
               << path << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    out << "# total_hash_entries=" << totalHashSize << "\n";
    out << "# merged_ambiguous_hash_entries=" << mergedAmbigHashSize << "\n";
    out << "# nCB=" << nCB << "\n";
    out << "stage\trecords\ttotal_count\thash_sum\thash_xor\thash_sum2\n";
    auto emit = [&](const char *stage, const BridgeStageDigest &d) {
        out << stage << '\t'
            << d.records << '\t'
            << d.totalCount << '\t'
            << d.hashSum << '\t'
            << d.hashXor << '\t'
            << d.hashSum2 << '\n';
    };
    emit("pre_cr_exact", preCr);
    emit("post_cr", postCr);
    emit("resolved_umi_gene", resolved);
    emit("final_gene_count", finalMatrix);
    out.close();

    if (cbPreCr == nullptr || cbPostCr == nullptr || cbResolved == nullptr || cbFinalMatrix == nullptr) {
        return;
    }

    const std::string cbPath = path + ".cb.tsv";
    std::ofstream cbOut(cbPath.c_str());
    if (!cbOut.good()) {
        ostringstream errOut;
        errOut << "EXITING because of fatal OUTPUT FILE error: could not open bridge per-CB digest file="
               << cbPath << "\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    cbOut << "stage\tiCB\twl_cb\tbarcode\trecords\ttotal_count\thash_sum\thash_xor\thash_sum2\n";
    auto emitCb = [&](const char *stage, const std::vector<BridgeStageDigest> &digests) {
        for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
            const BridgeStageDigest &d = digests[iCB];
            if (d.records == 0) {
                continue;
            }
            const uint32_t wlCb = indCB[iCB];
            cbOut << stage << '\t'
                  << iCB << '\t'
                  << wlCb << '\t';
            if (wlCb < pSolo.cbWLstr.size()) {
                cbOut << pSolo.cbWLstr[wlCb];
            }
            cbOut << '\t'
                  << d.records << '\t'
                  << d.totalCount << '\t'
                  << d.hashSum << '\t'
                  << d.hashXor << '\t'
                  << d.hashSum2 << '\n';
        }
    };
    emitCb("pre_cr_exact", *cbPreCr);
    emitCb("post_cr", *cbPostCr);
    emitCb("resolved_umi_gene", *cbResolved);
    emitCb("final_gene_count", *cbFinalMatrix);
}

} // namespace

void SoloFeature::collapseUMIall_fromBridgeHash()
{
    const auto t0 = std::chrono::steady_clock::now();
    auto elapsedSec = [&t0]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };

    const char *why = nullptr;
    if (pSolo.trackReadIdsForTags)
        why = "trackReadIdsForTags is enabled (needs per-read recordReadInfo)";
    else if (P.outSAMtype.empty() || P.outSAMtype[0] != "None")
        why = "--outSAMtype must be None for direct bridge collapse (or disable STAR_SOLO_NONFLEX_HASH_BRIDGE)";
    else if (pSolo.multiMap.yes.multi)
        why = "multimapper Solo output is not supported on this path (use --soloMultiMappers Unique)";
    else if (!pSolo.umiFiltering.MultiGeneUMI_CR)
        why = "--soloUMIfiltering MultiGeneUMI_CR is required";
    else if (pSolo.umiFiltering.MultiGeneUMI || pSolo.umiFiltering.MultiGeneUMI_All)
        why = "MultiGeneUMI / MultiGeneUMI_All cannot be combined with this path";
    else if (!pSolo.umiDedup.yes.CR || pSolo.umiDedup.yes.N != 1u)
        why = "exactly one --soloUMIdedup mode and it must be 1MM_CR for direct bridge collapse";

    if (why != nullptr) {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETER error: non-Flex inline-hash bridge direct collapse is not "
                  "supported for this configuration.\n"
               << "Reason: " << why << "\n"
               << "SOLUTION: use --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloMultiMappers Unique "
                  "--outSAMtype None, do not request sorted-BAM read-id tracking, or unset STAR_SOLO_NONFLEX_HASH_BRIDGE.\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    const char *hdMoleculeOutPath = std::getenv("STAR_SOLO_HD_MOLECULE_OUT");
    const bool hdMoleculeOut = hdMoleculeOutPath != nullptr && hdMoleculeOutPath[0] != '\0';
    SoloHDBarcodeMap hdBarcodeMap;
    if (hdMoleculeOut) {
        if (pSolo.umiL == 0 || pSolo.umiL > 12) {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETER error: STAR_SOLO_HD_MOLECULE_OUT requires "
                      "1 <= --soloUMIlen <= 12 on the non-Flex direct-hash bridge path.\n"
                   << "Reason: this bridge stores UMI identity in the existing 24-bit internal hash key.\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }

        const char *coordPathEnv = std::getenv("STAR_SOLO_HD_BARCODE_COORDS");
        const std::string coordPath = coordPathEnv != nullptr ? coordPathEnv : "";
        std::string hdError;
        if (!hdBarcodeMap.initialize(pSolo, coordPath, hdError)) {
            ostringstream errOut;
            errOut << "EXITING because of fatal PARAMETER error: could not initialize HD barcode coordinates for "
                      "STAR_SOLO_HD_MOLECULE_OUT.\n";
            if (!hdError.empty()) {
                errOut << "Reason: " << hdError << "\n";
            } else {
                errOut << "Reason: no whitelist barcodes resolved to s_002um row/column coordinates.\n";
            }
            errOut << "SOLUTION: use s_002um_<row>_<col> barcodes in the Solo whitelist/output whitelist, or set "
                      "STAR_SOLO_HD_BARCODE_COORDS to a TSV/CSV mapping barcode row2 col2.\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        P.inOut->logMain << "Experimental STAR-Spatial HD molecule export enabled: "
                         << hdMoleculeOutPath
                         << " coordinate_barcodes=" << hdBarcodeMap.validCount()
                         << endl;
    }

    if (!readFeatSum) {
        P.inOut->logMain << "WARNING: collapseUMIall_fromBridgeHash: readFeatSum is null" << endl;
        nCB = 0;
        return;
    }

    size_t rawThreadHashSize = 0;
    size_t threadHashCount = 0;
    for (int ii = 0; ii < P.runThreadN; ++ii) {
        if (readFeatAll[ii] && readFeatAll[ii]->inlineHash_) {
            const size_t sz = kh_size(readFeatAll[ii]->inlineHash_);
            if (sz > 0) {
                rawThreadHashSize += sz;
                ++threadHashCount;
            }
        }
    }
    const size_t mergedAmbigHashSize = bridgeHashLiveSize(readFeatSum->inlineHash_);
    const size_t rawHashSize = rawThreadHashSize + mergedAmbigHashSize;

    if (readFeatSum->inlineHash_ == nullptr) {
        readFeatSum->inlineHash_ = kh_init(cg_agg);
    }

    const size_t khMax = static_cast<size_t>(std::numeric_limits<khint_t>::max());
    const size_t resizeTarget = rawHashSize > (khMax - 1) / 2 ? khMax : rawHashSize * 2 + 1;
    if (resizeTarget > 0) {
        kh_resize(cg_agg, readFeatSum->inlineHash_, static_cast<khint_t>(resizeTarget));
    }
    for (int ii = 0; ii < P.runThreadN; ++ii) {
        if (readFeatAll[ii] && readFeatAll[ii]->inlineHash_) {
            mergeBridgeHashInto(readFeatAll[ii]->inlineHash_, readFeatSum->inlineHash_);
            kh_destroy(cg_agg, readFeatAll[ii]->inlineHash_);
            readFeatAll[ii]->inlineHash_ = nullptr;
        }
    }

    size_t totalHashSize = bridgeHashLiveSize(readFeatSum->inlineHash_);
    const size_t coalescedHashEntries = rawHashSize >= totalHashSize ? rawHashSize - totalHashSize : 0;

    if (totalHashSize == 0) {
        P.inOut->logMain << "WARNING: collapseUMIall_fromBridgeHash: no thread-local or merged bridge hash entries"
                         << endl;
        if (const char *snapPath = std::getenv("STAR_SOLO_BRIDGE_HASH_SNAPSHOT_OUT")) {
            if (snapPath[0] != '\0')
                bridgeHashSnapshotWrite(snapPath);
        }
        if (readFeatSum->inlineHash_) {
            kh_destroy(cg_agg, readFeatSum->inlineHash_);
            readFeatSum->inlineHash_ = nullptr;
        }
        nCB = 0;
        return;
    }

    const std::vector<uint32> *exactCbReadCount = nullptr;
    if (pSolo.CBmatchWL.oneExact) {
        if (readBarSum == nullptr || readBarSum->cbReadCountExact.size() < pSolo.cbWLsize) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: non-Flex direct bridge oneExact gating requires exact CB "
                      "read counts for every whitelist barcode.\n"
                   << "Reason: readBarSum->cbReadCountExact is missing or shorter than the whitelist "
                   << "(size="
                   << (readBarSum == nullptr ? 0 : readBarSum->cbReadCountExact.size())
                   << ", whitelist_size=" << pSolo.cbWLsize << ").\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
        }
        exactCbReadCount = &readBarSum->cbReadCountExact;
    }
    auto cbAllowedByOneExact = [&](uint32_t wlCb) -> bool {
        return exactCbReadCount == nullptr || (*exactCbReadCount)[wlCb] > 0;
    };

    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ... Direct bridge-hash UMI collapse (canonical hash + CSR + flat MultiGene), hash_entries="
                     << totalHashSize
                     << " raw_thread_hashes=" << threadHashCount
                     << " raw_thread_hash_entries=" << rawThreadHashSize
                     << " merged_ambiguous_hash_entries=" << mergedAmbigHashSize
                     << " coalesced_hash_entries=" << coalescedHashEntries
                     << endl;

    if (const char *snapPath = std::getenv("STAR_SOLO_BRIDGE_HASH_SNAPSHOT_OUT")) {
        if (snapPath[0] != '\0')
            bridgeHashSnapshotWrite(snapPath);
    }

    const char *determinismPathEnv = std::getenv("STAR_SOLO_BRIDGE_DETERMINISM_OUT");
    const bool bridgeDeterminismTrace = determinismPathEnv != nullptr && determinismPathEnv[0] != '\0';
    const std::string bridgeDeterminismPath = bridgeDeterminismTrace ? determinismPathEnv : "";
    const bool bridgeDeterminismPerCb = bridgeDeterminismTrace
        && std::getenv("STAR_SOLO_BRIDGE_DETERMINISM_CB") != nullptr;

    khash_t(cbcount) *cbEntryCount = kh_init(cbcount);
    kh_resize(cbcount, cbEntryCount, std::min<size_t>(totalHashSize / 8 + 64, size_t{1} << 20));
    size_t acceptedHashSize = 0;

    auto countHashEntriesPerWlCb = [&](SoloReadFeature *srcFeat, khash_t(cg_agg) *hash) {
        (void)srcFeat;
        if (hash == nullptr)
            return;
        for (khiter_t it = kh_begin(hash); it != kh_end(hash); ++it) {
            if (!kh_exist(hash, it))
                continue;
            const uint64_t key = kh_key(hash, it);
            uint32_t wlCb = 0, umi24 = 0;
            uint16_t gene16 = 0;
            unpackBridgeWlUmiGeneKey(key, &wlCb, &umi24, &gene16);
            (void)umi24;
            (void)gene16;
            if (!cbAllowedByOneExact(wlCb))
                continue;
            ++acceptedHashSize;
            int absent = 0;
            khiter_t itCb = kh_put(cbcount, cbEntryCount, wlCb, &absent);
            if (absent)
                kh_val(cbEntryCount, itCb) = 1;
            else
                ++kh_val(cbEntryCount, itCb);
        }
    };

    for (int ii = 0; ii < P.runThreadN; ++ii) {
        if (readFeatAll[ii] && readFeatAll[ii]->inlineHash_ && kh_size(readFeatAll[ii]->inlineHash_) > 0) {
            countHashEntriesPerWlCb(readFeatAll[ii], readFeatAll[ii]->inlineHash_);
        }
    }
    if (readFeatSum->inlineHash_ && kh_size(readFeatSum->inlineHash_) > 0) {
        countHashEntriesPerWlCb(readFeatSum, readFeatSum->inlineHash_);
    }

    if (acceptedHashSize == 0) {
        kh_destroy(cbcount, cbEntryCount);
        cbEntryCount = nullptr;
        P.inOut->logMain << "WARNING: collapseUMIall_fromBridgeHash: no bridge hash entries remain after "
                            "oneExact CB gating"
                         << endl;
        if (readFeatSum->inlineHash_) {
            kh_destroy(cg_agg, readFeatSum->inlineHash_);
            readFeatSum->inlineHash_ = nullptr;
        }
        nCB = 0;
        return;
    }
    if (acceptedHashSize != totalHashSize) {
        P.inOut->logMain << "Direct bridge-hash oneExact gating retained " << acceptedHashSize
                         << " of " << totalHashSize << " hash entries" << endl;
        totalHashSize = acceptedHashSize;
    }

    std::vector<uint32_t> sortedCBs;
    sortedCBs.reserve(kh_size(cbEntryCount));
    for (khiter_t it = kh_begin(cbEntryCount); it != kh_end(cbEntryCount); ++it) {
        if (!kh_exist(cbEntryCount, it))
            continue;
        sortedCBs.push_back(kh_key(cbEntryCount, it));
    }
    std::sort(sortedCBs.begin(), sortedCBs.end());

    nCB = static_cast<uint32_t>(sortedCBs.size());
    indCB.resize(nCB);
    indCBwl.assign(pSolo.cbWLsize, static_cast<uint32_t>(-1));
    for (uint32_t i = 0; i < nCB; ++i) {
        indCB[i] = sortedCBs[i];
        indCBwl[sortedCBs[i]] = i;
    }

    std::vector<uint32_t> hdRow2;
    std::vector<uint32_t> hdCol2;
    SoloHDMoleculeWriter hdWriter;
    bool hdWriteOk = true;
    std::string hdWriteError;
    if (hdMoleculeOut) {
        hdRow2.resize(nCB);
        hdCol2.resize(nCB);
        for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
            if (!hdBarcodeMap.coordForWL(indCB[iCB], hdRow2[iCB], hdCol2[iCB])) {
                ostringstream errOut;
                errOut << "EXITING because of fatal INPUT error: observed Solo whitelist barcode does not have an "
                          "HD coordinate mapping.\n"
                       << "Observed whitelist index: " << indCB[iCB] << "\n";
                if (indCB[iCB] < pSolo.cbWLstr.size()) {
                    errOut << "Barcode: " << pSolo.cbWLstr[indCB[iCB]] << "\n";
                }
                errOut << "SOLUTION: include this barcode in STAR_SOLO_HD_BARCODE_COORDS or use an s_002um output "
                          "whitelist.\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            }
        }
        if (!hdWriter.open(hdMoleculeOutPath,
                           featuresNumber,
                           SoloHDMoleculeWriter::kFlagPrededuplicated,
                           hdWriteError)) {
            ostringstream errOut;
            errOut << "EXITING because of fatal OUTPUT FILE error: could not open STAR-Spatial HD molecule stream.\n"
                   << "Reason: " << hdWriteError << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
    }

    std::vector<size_t> cbOffsets(nCB + 1, 0);
    for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
        const uint32_t wl = indCB[iCB];
        const khiter_t itCb = kh_get(cbcount, cbEntryCount, wl);
        const size_t wlCount = (itCb != kh_end(cbEntryCount)) ? kh_val(cbEntryCount, itCb) : 0;
        cbOffsets[iCB + 1] = cbOffsets[iCB] + wlCount;
    }
    kh_destroy(cbcount, cbEntryCount);
    cbEntryCount = nullptr;
    if (cbOffsets[nCB] != totalHashSize) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: bridge hash CSR sizing mismatch (expected " << totalHashSize
               << " entries, got " << cbOffsets[nCB] << ").\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
    }

    std::vector<BridgeFlatSlot> flat;
    flat.resize(totalHashSize);
    std::vector<size_t> cbWrite(nCB);
    for (uint32_t iCB = 0; iCB < nCB; ++iCB)
        cbWrite[iCB] = cbOffsets[iCB];

    auto drainHashIntoFlat = [&](SoloReadFeature *srcFeat, khash_t(cg_agg) *&hash) {
        (void)srcFeat;
        if (hash == nullptr)
            return;
        for (khiter_t it = kh_begin(hash); it != kh_end(hash); ++it) {
            if (!kh_exist(hash, it))
                continue;
            const uint64_t key = kh_key(hash, it);
            const uint32_t val = kh_val(hash, it);

            uint32_t wlCb = 0, umi24 = 0;
            uint16_t gene16 = 0;
            unpackBridgeWlUmiGeneKey(key, &wlCb, &umi24, &gene16);
            if (!cbAllowedByOneExact(wlCb))
                continue;

            const uint32_t iCB = indCBwl[wlCb];
            const size_t pos = cbWrite[iCB]++;
            flat[pos].gene = static_cast<uint32_t>(gene16);
            flat[pos].umi24 = umi24;
            flat[pos].count = val;
        }

        kh_destroy(cg_agg, hash);
        hash = nullptr;
    };

    for (int ii = 0; ii < P.runThreadN; ++ii) {
        if (readFeatAll[ii] && readFeatAll[ii]->inlineHash_ && kh_size(readFeatAll[ii]->inlineHash_) > 0) {
            drainHashIntoFlat(readFeatAll[ii], readFeatAll[ii]->inlineHash_);
        }
    }

    if (readFeatSum->inlineHash_) {
        if (kh_size(readFeatSum->inlineHash_) > 0) {
            drainHashIntoFlat(readFeatSum, readFeatSum->inlineHash_);
        } else {
            kh_destroy(cg_agg, readFeatSum->inlineHash_);
            readFeatSum->inlineHash_ = nullptr;
        }
    }

    for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
        if (cbWrite[iCB] != cbOffsets[iCB + 1]) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: bridge hash CSR fill mismatch for CB slot " << iCB << ".\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
        }
    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... CSR flat layout built (entries=" << flat.size() << ")"
                     << endl;

    countMatStride = pSolo.umiDedup.yes.N + 1;
    nUMIperCB.assign(nCB, 0);
    nGenePerCB.assign(nCB, 0);
    nReadPerCB.assign(nCB, 0);
    nReadPerCBmax = 0;
    countCellGeneUMIindex.assign(nCB + 1, 0);

    countMatMult.s = 1 + pSolo.multiMap.yes.N * pSolo.umiDedup.yes.N;
    countMatMult.m.clear();
    countMatMult.i.assign(nCB + 1, 0);

    struct ThreadScratch {
        std::vector<uint32_t> umiArray;
        std::vector<MgRow> mgBuf;
        std::vector<std::pair<uint32_t, uint32_t>> aggGene;
        std::vector<std::pair<uint32_t, uint32_t>> origAtCu;
        std::vector<multi_gene_umi_cr::GeneSupport> resolutionSupports;
        vector<uint32> geneCounts;
        std::unordered_map<uintUMI, uintUMI> emptyUmiCorr;
        std::vector<uint32_t> geneStamp;
        uint32_t stampVer;
        size_t nCbGeneSegT;
        uint32_t nReadPerCBmaxT;
        std::vector<uint32_t> geneSeenBkt;
        std::vector<uint32_t> geneDupCount;
        std::vector<uint32_t> geneWriteScratch;
        uint32_t geneBucketMark;
        std::vector<uint32_t> touchedGenes;
        std::vector<size_t> genePref;
        std::vector<BridgeFlatSlot> geneBucketBuf;
        std::vector<uint32_t> gID;
        std::vector<size_t> gBeg;
        std::vector<size_t> gEnd;
        std::vector<SoloHDMoleculeRecord> hdRecords;
        BridgeStageDigest detPreCr;
        BridgeStageDigest detPostCr;
        BridgeStageDigest detResolved;
        BridgeStageDigest detFinalMatrix;

        ThreadScratch()
            : stampVer(1),
              nCbGeneSegT(0),
              nReadPerCBmaxT(0),
              geneBucketMark(1),
              geneStamp(65536, 0),
              geneSeenBkt(65536, 0),
              geneDupCount(65536, 0),
              geneWriteScratch(65536, 0)
        {
            umiArray.reserve(umiArrayStride * 64);
            touchedGenes.reserve(8192);
            genePref.reserve(8192);
            hdRecords.reserve(1024);
        }
    };

    std::vector<uint32_t> cbRowCounts(nCB, 0);
    std::vector<std::vector<uint32_t>> cbRows(nCB);
    std::vector<BridgeStageDigest> cbDetPreCr;
    std::vector<BridgeStageDigest> cbDetPostCr;
    std::vector<BridgeStageDigest> cbDetResolved;
    std::vector<BridgeStageDigest> cbDetFinalMatrix;
    if (bridgeDeterminismPerCb) {
        cbDetPreCr.resize(nCB);
        cbDetPostCr.resize(nCB);
        cbDetResolved.resize(nCB);
        cbDetFinalMatrix.resize(nCB);
    }
    const int collapseThreads = std::max(1, P.runThreadN);
    std::vector<ThreadScratch> scratch(static_cast<size_t>(collapseThreads));

    #pragma omp parallel num_threads(collapseThreads)
    {
        ThreadScratch &ts = scratch[static_cast<size_t>(omp_get_thread_num())];

        #pragma omp for schedule(dynamic, 16)
        for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
            const uint32_t wlCb = indCB[iCB];
            const size_t cbBeg = cbOffsets[iCB];
            const size_t cbEnd = cbOffsets[iCB + 1];
            const size_t sliceLen = cbEnd - cbBeg;
            BridgeFlatSlot *const rawSlice = flat.data() + cbBeg;

            nGenePerCB[iCB] = 0;
            nUMIperCB[iCB] = 0;
            cbRowCounts[iCB] = 0;
            cbRows[iCB].clear();
            if (hdMoleculeOut) {
                ts.hdRecords.clear();
            }

            if (sliceLen == 0) {
                if (++ts.stampVer == 0u) {
                    std::fill(ts.geneStamp.begin(), ts.geneStamp.end(), 0u);
                    ts.stampVer = 1u;
                }
                continue;
            }

            uint64_t readSum = 0;
            for (size_t j = cbBeg; j < cbEnd; ++j) {
                readSum += flat[j].count;
                const uint32_t g = flat[j].gene;
                if (ts.geneStamp[g] != ts.stampVer) {
                    ts.geneStamp[g] = ts.stampVer;
                    ++ts.nCbGeneSegT;
                }
            }
            nReadPerCB[iCB] = static_cast<uint32_t>(readSum > UINT32_MAX ? UINT32_MAX : readSum);
            if (nReadPerCB[iCB] > ts.nReadPerCBmaxT)
                ts.nReadPerCBmaxT = nReadPerCB[iCB];
            if (++ts.stampVer == 0u) {
                std::fill(ts.geneStamp.begin(), ts.geneStamp.end(), 0u);
                ts.stampVer = 1u;
            }

            const uint32_t bktMark = ts.geneBucketMark++;
            if (ts.geneBucketMark == 0u) {
                std::fill(ts.geneSeenBkt.begin(), ts.geneSeenBkt.end(), 0u);
                ts.geneBucketMark = 1u;
            }
            ts.touchedGenes.clear();
            for (size_t j = 0; j < sliceLen; ++j) {
                const uint32_t g = rawSlice[j].gene;
                if (g & geneMultMark) {
                    ostringstream errOut;
                    errOut << "EXITING because of fatal ERROR: multimapper gene encountered in non-Flex bridge direct "
                              "collapse.\n"
                           << "This path requires unique-gene alignments only (--soloMultiMappers Unique).\n";
                    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
                }
                if (ts.geneSeenBkt[g] != bktMark) {
                    ts.geneSeenBkt[g] = bktMark;
                    ts.touchedGenes.push_back(g);
                    ts.geneDupCount[g] = 1;
                } else {
                    ts.geneDupCount[g]++;
                }
            }

            std::sort(ts.touchedGenes.begin(), ts.touchedGenes.end());

            ts.genePref.resize(ts.touchedGenes.size() + 1);
            ts.genePref[0] = 0;
            for (size_t i = 0; i < ts.touchedGenes.size(); ++i)
                ts.genePref[i + 1] = ts.genePref[i] + static_cast<size_t>(ts.geneDupCount[ts.touchedGenes[i]]);
            if (ts.genePref.back() != sliceLen) {
                ostringstream errOut;
                errOut << "EXITING because of fatal ERROR: bridge gene bucket size mismatch (expected " << sliceLen
                       << ", got " << ts.genePref.back() << ").\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
            }

            if (ts.geneBucketBuf.size() < sliceLen)
                ts.geneBucketBuf.resize(sliceLen);
            BridgeFlatSlot *slice = ts.geneBucketBuf.data();
            for (size_t i = 0; i < ts.touchedGenes.size(); ++i) {
                const uint32_t g = ts.touchedGenes[i];
                ts.geneWriteScratch[g] = static_cast<uint32_t>(ts.genePref[i]);
            }
            for (size_t j = 0; j < sliceLen; ++j) {
                const uint32_t g = rawSlice[j].gene;
                slice[ts.geneWriteScratch[g]++] = rawSlice[j];
            }

            ts.gID.clear();
            ts.gBeg.clear();
            ts.gEnd.clear();
            ts.gID.reserve(ts.touchedGenes.size());
            ts.gBeg.reserve(ts.touchedGenes.size());
            ts.gEnd.reserve(ts.touchedGenes.size());
            for (size_t i = 0; i < ts.touchedGenes.size(); ++i) {
                ts.gID.push_back(ts.touchedGenes[i]);
                ts.gBeg.push_back(ts.genePref[i]);
                ts.gEnd.push_back(ts.genePref[i + 1]);
            }

            const uint32_t nGenes = static_cast<uint32_t>(ts.gID.size());
            ts.mgBuf.clear();
            ts.mgBuf.reserve(sliceLen);

            for (uint32_t iG = 0; iG < nGenes; ++iG) {
                const uint32_t geneId = ts.gID[iG];
                const size_t a = ts.gBeg[iG];
                const size_t b = ts.gEnd[iG];
                if (b == a)
                    continue;

                std::sort(slice + a, slice + b,
                          [](const BridgeFlatSlot &x, const BridgeFlatSlot &y) {
                              return x.umi24 < y.umi24;
                          });

                if (ts.umiArray.size() < (b - a) * umiArrayStride)
                    ts.umiArray.resize((b - a) * umiArrayStride);

                uint32_t nU0 = 0;
                for (size_t t = a; t < b;) {
                    const uint32_t umi24 = slice[t].umi24;
                    uint64_t count = 0;
                    while (t < b && slice[t].umi24 == umi24) {
                        count += slice[t].count;
                        ++t;
                    }
                    if (count > std::numeric_limits<uint32_t>::max()) {
                        ostringstream errOut;
                        errOut << "EXITING because of fatal ERROR: bridge exact UMI count overflow for CB "
                               << iCB << ", gene " << geneId << ".\n";
                        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
                    }
                    const uint32_t count32 = static_cast<uint32_t>(count);
                    ts.umiArray[static_cast<size_t>(nU0) * umiArrayStride + 0] = umi24;
                    ts.umiArray[static_cast<size_t>(nU0) * umiArrayStride + 1] = count32;
                    ts.umiArray[static_cast<size_t>(nU0) * umiArrayStride + 2] = static_cast<uint32_t>(-1);
                    if (bridgeDeterminismTrace) {
                        const uint64_t h = bridgeTupleHash(0x50524543525f4558ull, wlCb, geneId, umi24, 0, count32);
                        ts.detPreCr.add(h, count32);
                        if (bridgeDeterminismPerCb) {
                            cbDetPreCr[iCB].add(h, count32);
                        }
                    }
                    ++nU0;
                }

                umiArrayCorrect_CR(nU0, ts.umiArray.data(), false, false, ts.emptyUmiCorr);

                for (uint64_t iu = 0; iu < static_cast<uint64_t>(nU0) * umiArrayStride; iu += umiArrayStride) {
                    MgRow row;
                    row.orig = ts.umiArray[iu + 0];
                    row.corr = ts.umiArray[iu + 2];
                    row.geneIdx = iG;
                    row.count = ts.umiArray[iu + 1];
                    if (bridgeDeterminismTrace) {
                        const uint64_t h = bridgeTupleHash(0x504f53545f43525full, wlCb, geneId,
                                                           row.orig, row.corr, row.count);
                        ts.detPostCr.add(h, row.count);
                        if (bridgeDeterminismPerCb) {
                            cbDetPostCr[iCB].add(h, row.count);
                        }
                    }
                    ts.mgBuf.push_back(row);
                }
            }

            radixSortMgRowsByCorr(ts.mgBuf);
            ts.geneCounts.assign(static_cast<size_t>(nGenes), 0);

            for (size_t p = 0; p < ts.mgBuf.size();) {
                size_t q = p + 1;
                const uint32_t cu = ts.mgBuf[p].corr;
                while (q < ts.mgBuf.size() && ts.mgBuf[q].corr == cu)
                    ++q;

                std::sort(ts.mgBuf.begin() + p, ts.mgBuf.begin() + q,
                          [](const MgRow &x, const MgRow &y) { return x.geneIdx < y.geneIdx; });

                ts.aggGene.clear();
                for (size_t i = p; i < q;) {
                    const uint32_t g0 = ts.mgBuf[i].geneIdx;
                    uint32_t s = 0;
                    while (i < q && ts.mgBuf[i].geneIdx == g0) {
                        s += ts.mgBuf[i].count;
                        ++i;
                    }
                    ts.aggGene.push_back({g0, s});
                }

                ts.origAtCu.clear();
                for (size_t i = p; i < q; ++i) {
                    if (ts.mgBuf[i].orig == cu)
                        ts.origAtCu.push_back({ts.mgBuf[i].geneIdx, ts.mgBuf[i].count});
                }
                std::sort(ts.origAtCu.begin(), ts.origAtCu.end(),
                          [](const std::pair<uint32_t, uint32_t> &x, const std::pair<uint32_t, uint32_t> &y) {
                              return x.first < y.first;
                          });
                {
                    size_t o = 0;
                    for (size_t i = 0; i < ts.origAtCu.size();) {
                        const uint32_t g0 = ts.origAtCu[i].first;
                        uint32_t s = 0;
                        while (i < ts.origAtCu.size() && ts.origAtCu[i].first == g0) {
                            s += ts.origAtCu[i].second;
                            ++i;
                        }
                        ts.origAtCu[o++] = {g0, s};
                    }
                    ts.origAtCu.resize(o);
                }

                ts.resolutionSupports.clear();
                ts.resolutionSupports.reserve(ts.aggGene.size());
                size_t originalIndex = 0;
                for (const auto &pr : ts.aggGene) {
                    while (originalIndex < ts.origAtCu.size()
                           && ts.origAtCu[originalIndex].first < pr.first) {
                        ++originalIndex;
                    }
                    multi_gene_umi_cr::GeneSupport support;
                    support.gene = pr.first;
                    support.correctedCount = pr.second;
                    support.originalAtCorrectedCount =
                        originalIndex < ts.origAtCu.size()
                            && ts.origAtCu[originalIndex].first == pr.first
                        ? ts.origAtCu[originalIndex].second : 0;
                    ts.resolutionSupports.push_back(support);
                }
                const multi_gene_umi_cr::Result resolution =
                    multi_gene_umi_cr::resolve(ts.resolutionSupports);
                const uint32_t maxg = resolution.accepted
                    ? resolution.gene : static_cast<uint32_t>(-1);

                if (shouldTraceCollapseBarcode(pSolo, indCB[iCB])) {
                    const int64_t chosenGene = (maxg + 1u == 0u) ? -1 : static_cast<int64_t>(ts.gID[maxg]);
                    P.inOut->logMain << "[GENEFULL-CR-TRACE] mode=bridge"
                                     << " cb=" << pSolo.cbWLstr[indCB[iCB]]
                                     << " corr=" << cu
                                     << " corrected={" << formatBridgeGeneCounts(ts.aggGene, ts.gID) << '}'
                                     << " orig={" << formatBridgeGeneCounts(ts.origAtCu, ts.gID) << '}'
                                     << " chosen=" << chosenGene
                                     << endl;
                }

                if (maxg + 1u != 0u) {
                    const uint32_t chosenGene = ts.gID[maxg];
                    if (bridgeDeterminismTrace) {
                        const uint64_t h = bridgeTupleHash(0x5245534f4c564544ull, wlCb, chosenGene, cu, cu, 1);
                        ts.detResolved.add(h, 1);
                        if (bridgeDeterminismPerCb) {
                            cbDetResolved[iCB].add(h, 1);
                        }
                    }
                    ts.geneCounts[maxg]++;
                    if (hdMoleculeOut) {
                        SoloHDMoleculeRecord record;
                        record.row2 = hdRow2[iCB];
                        record.col2 = hdCol2[iCB];
                        record.featureIdx = chosenGene;
                        record.umiLength = static_cast<uint16_t>(pSolo.umiL);
                        record.flags = 0;
                        if (!SoloHDMoleculeWriter::packStarUmi(
                                cu, record.umiLength, record.umiLow, record.umiHigh)) {
                            #pragma omp critical(star_hd_molecule_writer)
                            {
                                if (hdWriteOk) {
                                    hdWriteOk = false;
                                    hdWriteError = "could not pack corrected STAR UMI for HD molecule output";
                                }
                            }
                        } else {
                            ts.hdRecords.push_back(record);
                        }
                    }
                }

                p = q;
            }

            size_t cbNonzeroGenes = 0;
            for (uint32_t ig = 0; ig < nGenes; ++ig)
                if (ts.geneCounts[ig] != 0)
                    ++cbNonzeroGenes;

            cbRowCounts[iCB] = static_cast<uint32_t>(cbNonzeroGenes);
            std::vector<uint32_t> &rows = cbRows[iCB];
            rows.reserve(cbNonzeroGenes * countMatStride);
            for (uint32_t ig = 0; ig < nGenes; ++ig) {
                if (ts.geneCounts[ig] == 0)
                    continue;
                nGenePerCB[iCB]++;
                nUMIperCB[iCB] += ts.geneCounts[ig];
                if (bridgeDeterminismTrace) {
                    const uint64_t h = bridgeTupleHash(0x46494e414c4d4154ull, wlCb, ts.gID[ig],
                                                       0, 0, ts.geneCounts[ig]);
                    ts.detFinalMatrix.add(h, ts.geneCounts[ig]);
                    if (bridgeDeterminismPerCb) {
                        cbDetFinalMatrix[iCB].add(h, ts.geneCounts[ig]);
                    }
                }
                rows.push_back(ts.gID[ig]);
                rows.push_back(ts.geneCounts[ig]);
            }

            if (hdMoleculeOut && !ts.hdRecords.empty()) {
                #pragma omp critical(star_hd_molecule_writer)
                {
                    if (hdWriteOk && !hdWriter.writeRecords(ts.hdRecords, hdWriteError)) {
                        hdWriteOk = false;
                    }
                }
            }
        }
    }

    if (hdMoleculeOut) {
        if (!hdWriteOk) {
            ostringstream errOut;
            errOut << "EXITING because of fatal OUTPUT FILE error: failed while writing STAR-Spatial HD molecule "
                      "stream.\n"
                   << "Reason: " << hdWriteError << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        if (!hdWriter.close(hdWriteError)) {
            ostringstream errOut;
            errOut << "EXITING because of fatal OUTPUT FILE error: failed while finalizing STAR-Spatial HD molecule "
                      "stream.\n"
                   << "Reason: " << hdWriteError << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        }
        P.inOut->logMain << "Finished STAR-Spatial HD molecule export: records="
                         << hdWriter.recordCount()
                         << " path=" << hdMoleculeOutPath
                         << endl;
    }

    size_t nCbGeneSeg = 0;
    BridgeStageDigest detPreCr;
    BridgeStageDigest detPostCr;
    BridgeStageDigest detResolved;
    BridgeStageDigest detFinalMatrix;
    for (const ThreadScratch &ts : scratch) {
        nCbGeneSeg += ts.nCbGeneSegT;
        if (ts.nReadPerCBmaxT > nReadPerCBmax)
            nReadPerCBmax = ts.nReadPerCBmaxT;
        if (bridgeDeterminismTrace) {
            detPreCr.merge(ts.detPreCr);
            detPostCr.merge(ts.detPostCr);
            detResolved.merge(ts.detResolved);
            detFinalMatrix.merge(ts.detFinalMatrix);
        }
    }

    {
        std::vector<BridgeFlatSlot>().swap(flat);
        std::vector<size_t>().swap(cbWrite);
        std::vector<size_t>().swap(cbOffsets);
    }

    size_t finalMatSlots = 0;
    for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
        countCellGeneUMIindex[iCB] = static_cast<uint32_t>(finalMatSlots);
        finalMatSlots += static_cast<size_t>(cbRowCounts[iCB]) * countMatStride;
    }
    countCellGeneUMIindex[nCB] = static_cast<uint32_t>(finalMatSlots);
    countCellGeneUMI.assign(finalMatSlots, 0);

    for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
        const std::vector<uint32_t> &rows = cbRows[iCB];
        if (!rows.empty())
            std::copy(rows.begin(), rows.end(), countCellGeneUMI.begin() + countCellGeneUMIindex[iCB]);
        std::vector<uint32_t>().swap(cbRows[iCB]);
    }

    if (bridgeDeterminismTrace) {
        writeBridgeDeterminismDigest(bridgeDeterminismPath,
                                     pSolo,
                                     nCB,
                                     indCB,
                                     totalHashSize,
                                     mergedAmbigHashSize,
                                     detPreCr,
                                     detPostCr,
                                     detResolved,
                                     detFinalMatrix,
                                     bridgeDeterminismPerCb ? &cbDetPreCr : nullptr,
                                     bridgeDeterminismPerCb ? &cbDetPostCr : nullptr,
                                     bridgeDeterminismPerCb ? &cbDetResolved : nullptr,
                                     bridgeDeterminismPerCb ? &cbDetFinalMatrix : nullptr,
                                     P);
        P.inOut->logMain << "Wrote bridge determinism digest: " << bridgeDeterminismPath << endl;
    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ... Parallel CB collapse complete, exact matrix slots=" << finalMatSlots
                     << " collapse_threads=" << collapseThreads
                     << endl;

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Unique (CB,gene) segments (matrix row upper bound)="
                     << nCbGeneSeg << " | Finished direct bridge-hash UMI collapse, nCB=" << nCB << " wall="
                     << elapsedSec() << " s" << endl;
}
