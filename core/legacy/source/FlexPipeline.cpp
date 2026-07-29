#include "FlexPipeline.h"
#include "FlexHashScreen.h"
#include "SoloReadBarcode.h"
#include "SoloReadFeature.h"
#include "SoloRead.h"
#include "SoloFeatureTypes.h"
#include "SampleDetector.h"
#include "SequenceFuns.h"
#include "Parameters.h"
#include "ReadAlign.h"
#include "SpatialGex.h"
#include "Stats.h"
#include "ErrorWarning.h"
#include "input/CbqInputModule.h"

#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <unordered_set>
#include <unistd.h>
#include <chrono>

namespace {

static constexpr int kGzBufSize = 1 << 20;

bool gzReadLine(gzFile gz, char *buf, int maxLen) {
    char *ret = gzgets(gz, buf, maxLen);
    if (ret == Z_NULL) return false;
    int len = static_cast<int>(strlen(buf));
    if (len > 0 && buf[len - 1] == '\n') buf[len - 1] = '\0';
    return true;
}

bool copyCbqSpanToBuffer(star::input::CbqByteSpan span, char *dest, size_t capacity, size_t *lengthOut) {
    if (dest == nullptr || capacity == 0) return false;
    const size_t n = std::min(span.size, capacity - 1);
    if (n > 0 && span.data != nullptr) {
        std::memcpy(dest, span.data, n);
    }
    dest[n] = '\0';
    if (lengthOut != nullptr) *lengthOut = n;
    return span.size < capacity;
}

void copyCbqReadName(star::input::CbqByteSpan span, char *dest, size_t capacity) {
    size_t ignored = 0;
    copyCbqSpanToBuffer(span, dest, capacity, &ignored);
}

bool failSpatialFlexPipeline(Parameters &P, const std::string &detail) {
    exitWithError("EXITING because native spatial Flex pipeline processing failed: "
                      + detail + "\n",
                  std::cerr, P.inOut->logMain,
                  EXIT_CODE_INCONSISTENT_DATA, P);
    return false;
}

bool completeSpatialFlexCachedRead(Parameters &P,
                                   const char *barcodeSequence,
                                   uint32_t barcodeLength,
                                   const char *barcodeQuality,
                                   uint64_t sourceOrdinal,
                                   const FlexHashScreenDecision &decision) {
    if (!P.soloSpatialFlexIntegratedEnabled) return true;
    if (P.spatialGexPipeline == nullptr) {
        return failSpatialFlexPipeline(P, "integrated spatial accumulator is unavailable");
    }

    std::string error;
    if (!P.spatialGexPipeline->decodeCurrentThread(
            barcodeSequence, barcodeLength, barcodeQuality, barcodeLength,
            sourceOrdinal, error)) {
        return failSpatialFlexPipeline(P, "raw-R1 decoding failed: " + error);
    }

    spatial_gex::FeatureEvidenceClass source =
        spatial_gex::FeatureEvidenceClass::FlexHashDeny;
    bool assigned = false;
    uint32_t geneIndex = 0;
    if (decision.action == FlexHashScreenDecision::Keep) {
        if (decision.geneIdx15 == 0
            || (decision.cacheClass != 0 && decision.cacheClass != 1)) {
            return failSpatialFlexPipeline(P, "invalid H0/H1 cache keep decision");
        }
        source = decision.cacheClass == 0
            ? spatial_gex::FeatureEvidenceClass::FlexH0
            : spatial_gex::FeatureEvidenceClass::FlexH1;
        assigned = true;
        geneIndex = decision.geneIdx15 - 1;
    } else if (decision.action != FlexHashScreenDecision::Deny) {
        return failSpatialFlexPipeline(P, "non-terminal cache decision reached cached-read completion");
    }

    if (!P.spatialGexPipeline->completeCurrentThread(
            source, assigned, geneIndex, sourceOrdinal, error)) {
        return failSpatialFlexPipeline(P, "feature completion failed: " + error);
    }
    return true;
}

} // namespace

void *flexLaneReaderThread(void *arg) {
    FlexLaneReaderArgs *ctx = static_cast<FlexLaneReaderArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    const int laneId = ctx->laneId;
    gzFile gzR2 = ctx->gzR2;
    gzFile gzR1 = ctx->gzR1;

    char lineBuf[kFlexPipeSeqMax + 256];

    while (true) {
        ReadPacket pkt;
        pkt.laneId = static_cast<uint8_t>(laneId);
        pkt.readFilesIndex = static_cast<uint32_t>(laneId);
        pkt.eof = false;

        // R2: @name, seq, +, qual
        if (!gzReadLine(gzR2, lineBuf, sizeof(lineBuf))) break;
        {
            // Parse name line: @name [iReadAll] [readFilter] [readFilesIndex]
            const char *src = lineBuf;
            if (*src == '@') ++src;
            size_t nameLen = 0;
            while (src[nameLen] && src[nameLen] != ' ' && src[nameLen] != '\t')
                ++nameLen;
            if (nameLen >= kFlexPipeNameMax) nameLen = kFlexPipeNameMax - 1;
            std::memcpy(pkt.name, src, nameLen);
            pkt.name[nameLen] = '\0';
            // Parse readFilter from Illumina header if present
            const char *rest = src + nameLen;
            while (*rest == ' ' || *rest == '\t') ++rest;
            // iReadAll field (skip)
            while (*rest && *rest != ' ' && *rest != '\t') ++rest;
            while (*rest == ' ' || *rest == '\t') ++rest;
            // readFilter field
            if (*rest) pkt.readFilter = *rest;
        }

        if (!gzReadLine(gzR2, pkt.seq[0], kFlexPipeSeqMax)) break;
        pkt.readLen[0] = static_cast<uint32_t>(strlen(pkt.seq[0]));

        if (!gzReadLine(gzR2, lineBuf, sizeof(lineBuf))) break; // + line
        if (!gzReadLine(gzR2, pkt.qual[0], kFlexPipeSeqMax)) break;

        // R1: @name, seq, +, qual
        if (!gzReadLine(gzR1, lineBuf, sizeof(lineBuf))) break;
        if (!gzReadLine(gzR1, pkt.seq[1], kFlexPipeSeqMax)) break;
        pkt.readLen[1] = static_cast<uint32_t>(strlen(pkt.seq[1]));
        if (!gzReadLine(gzR1, lineBuf, sizeof(lineBuf))) break;
        if (!gzReadLine(gzR1, pkt.qual[1], kFlexPipeSeqMax)) break;

        pkt.iReadAll = st->iReadAllGlobal.fetch_add(1);

        st->counters.perLaneReads[laneId]++;
        st->readerQ.push(std::move(pkt));
    }

    gzclose(gzR2);
    gzclose(gzR1);

    int finishedCount = st->readersFinished.fetch_add(1) + 1;
    if (finishedCount == st->nLanes) {
        st->readerQ.close();
    }

    return nullptr;
}

void *flexLaneReaderRouterThread(void *arg) {
    FlexLaneReaderArgs *ctx = static_cast<FlexLaneReaderArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    Parameters &P = *ctx->P;
    const int laneId = ctx->laneId;
    gzFile gzR2 = ctx->gzR2;
    gzFile gzR1 = ctx->gzR1;

    auto &cache = FlexHashScreenCache::instance();
    const uint32_t sampleTagOffset = P.pSolo.sampleProbeOffset;

    char lineBuf[kFlexPipeSeqMax + 256];
    char seq0[kFlexPipeSeqMax], seq1[kFlexPipeSeqMax];
    char qual0[kFlexPipeSeqMax], qual1[kFlexPipeSeqMax];
    char name[kFlexPipeNameMax];

    while (true) {
        // R2: @name, seq, +, qual
        if (!gzReadLine(gzR2, lineBuf, sizeof(lineBuf))) break;
        {
            const char *src = lineBuf;
            if (*src == '@') ++src;
            size_t nameLen = 0;
            while (src[nameLen] && src[nameLen] != ' ' && src[nameLen] != '\t')
                ++nameLen;
            if (nameLen >= kFlexPipeNameMax) nameLen = kFlexPipeNameMax - 1;
            std::memcpy(name, src, nameLen);
            name[nameLen] = '\0';
        }
        if (!gzReadLine(gzR2, seq0, kFlexPipeSeqMax)) break;
        uint32_t readLen0 = static_cast<uint32_t>(strlen(seq0));
        if (!gzReadLine(gzR2, lineBuf, sizeof(lineBuf))) break;
        if (!gzReadLine(gzR2, qual0, kFlexPipeSeqMax)) break;

        // R1: @name, seq, +, qual
        if (!gzReadLine(gzR1, lineBuf, sizeof(lineBuf))) break;
        if (!gzReadLine(gzR1, seq1, kFlexPipeSeqMax)) break;
        uint32_t readLen1 = static_cast<uint32_t>(strlen(seq1));
        if (!gzReadLine(gzR1, lineBuf, sizeof(lineBuf))) break;
        if (!gzReadLine(gzR1, qual1, kFlexPipeSeqMax)) break;

        uint64_t iReadAll = st->iReadAllGlobal.fetch_add(1);
        st->counters.perLaneReads[laneId]++;

        FlexHashScreenDecision decision = P.soloSpatialFlexIntegratedEnabled
            ? cache.classifyReadH0H1Offset0(seq0, readLen0)
            : cache.classifyReadH0Offset0(seq0, readLen0);

        if (decision.action == FlexHashScreenDecision::Keep ||
            decision.action == FlexHashScreenDecision::Deny) {
            DecisionPacket dp;
            dp.iReadAll = iReadAll;
            dp.readFilesIndex = static_cast<uint32_t>(laneId);
            dp.eof = false;

            const uint32_t bcLen = std::min(readLen1, kFlexPipeBarcodeSeqMax - 1);
            std::memcpy(dp.barcodeSeq, seq1, bcLen);
            dp.barcodeSeq[bcLen] = '\0';
            std::memcpy(dp.barcodeQual, qual1, bcLen);
            dp.barcodeQual[bcLen] = '\0';
            dp.barcodeLen = bcLen;
            std::memcpy(dp.readName, name, kFlexPipeNameMax);

            if (sampleTagOffset + 8 <= readLen0) {
                std::memcpy(dp.sampleTagSeq, seq0 + sampleTagOffset, 8);
                dp.sampleTagLen = 8;
            } else {
                dp.sampleTagLen = 0;
            }

            if (decision.action == FlexHashScreenDecision::Keep) {
                dp.verdict = DecisionPacket::KEEP;
                dp.geneIdx15 = decision.geneIdx15;
                dp.cacheClass = decision.cacheClass;
                dp.probeRegion = decision.probeRegion;
                dp.denyReason = nullptr;
                st->counters.triageKeep.fetch_add(1);
            } else {
                dp.verdict = DecisionPacket::DENY;
                dp.geneIdx15 = 0;
                dp.cacheClass = 0;
                dp.denyReason = "NEG_PROBE_AMBIG";
                st->counters.triageDeny.fetch_add(1);
            }

            int shard = static_cast<int>(iReadAll % static_cast<uint64_t>(st->nSolo));
            st->soloQ[shard]->push(std::move(dp));
        } else {
            EnrichedPacket ep;
            std::memcpy(ep.name, name, kFlexPipeNameMax);
            std::memcpy(ep.seq[0], seq0, readLen0 + 1);
            std::memcpy(ep.seq[1], seq1, readLen1 + 1);
            std::memcpy(ep.qual[0], qual0, readLen0 + 1);
            std::memcpy(ep.qual[1], qual1, readLen1 + 1);
            ep.readLen[0] = readLen0;
            ep.readLen[1] = readLen1;
            ep.iReadAll = iReadAll;
            ep.laneId = static_cast<uint8_t>(laneId);
            ep.readFilesIndex = static_cast<uint32_t>(laneId);
            ep.readFilter = 'Y';
            ep.eof = false;
            ep.cbMatch = -1;
            ep.cbMatchIndN = 0;
            ep.umiB = 0;
            ep.detectedSampleToken = 0xFF;
            ep.hashScreenSampleIdx = 0;

            st->counters.triageMiss.fetch_add(1);
            st->alignQ.push(std::move(ep));
        }

        st->counters.readsTotal.fetch_add(1);
    }

    gzclose(gzR2);
    gzclose(gzR1);

    int finishedCount = st->readersFinished.fetch_add(1) + 1;
    if (finishedCount == st->nLanes) {
        for (int i = 0; i < st->nSolo; ++i)
            st->soloQ[i]->close();
        st->alignQ.close();
    }

    return nullptr;
}

// Shared one-time init for the sample pre-filter (thread-safe via call_once)
static std::unordered_set<uint32_t> g_samplePreFilter;
static std::once_flag g_samplePreFilterFlag;

static void buildSamplePreFilter(Parameters &P) {

    static uint8_t asciiToNib[256] = {};
    std::memset(asciiToNib, 0, sizeof(asciiToNib));
    asciiToNib['A'] = asciiToNib['a'] = 1;
    asciiToNib['C'] = asciiToNib['c'] = 2;
    asciiToNib['G'] = asciiToNib['g'] = 4;
    asciiToNib['T'] = asciiToNib['t'] = 8;

    if (P.pSolo.sampleWhitelistPath.empty() || P.pSolo.sampleWhitelistPath == "-" ||
        P.pSolo.sampleProbesPath.empty() || P.pSolo.sampleProbesPath == "-")
        return;

    std::unordered_set<std::string> userCanonicals;
    {
        std::ifstream wl(P.pSolo.sampleWhitelistPath.c_str());
        std::string line;
        while (std::getline(wl, line)) {
            std::istringstream iss(line);
            std::string id, canonical;
            if (iss >> id >> canonical) userCanonicals.insert(canonical);
        }
    }
    auto encodeTag = [&](const std::string &tag) -> uint32_t {
        if (tag.size() != 8) return 0;
        uint32_t code = 0;
        for (int i = 0; i < 8; ++i) {
            uint8_t nib = asciiToNib[static_cast<uint8_t>(tag[i])];
            if (nib == 0) return 0;
            code = (code << 4) | nib;
        }
        return code;
    };
    for (const std::string &canon : userCanonicals) {
        uint32_t c = encodeTag(canon);
        if (c) g_samplePreFilter.insert(c);
    }
    std::ifstream probeFile(P.pSolo.sampleProbesPath.c_str());
    std::string line;
    while (std::getline(probeFile, line)) {
        std::istringstream iss(line);
        std::string variant, canonical, barcodeId;
        if (!(iss >> variant >> canonical >> barcodeId)) continue;
        if (userCanonicals.find(canonical) == userCanonicals.end()) continue;
        uint32_t c = encodeTag(variant);
        if (c) g_samplePreFilter.insert(c);
    }
}

static uint8_t g_asciiToNib[256];
static std::once_flag g_asciiToNibFlag;

static void initAsciiToNib() {
    std::memset(g_asciiToNib, 0, sizeof(g_asciiToNib));
    g_asciiToNib['A'] = g_asciiToNib['a'] = 1;
    g_asciiToNib['C'] = g_asciiToNib['c'] = 2;
    g_asciiToNib['G'] = g_asciiToNib['g'] = 4;
    g_asciiToNib['T'] = g_asciiToNib['t'] = 8;
}

static void ensureAsciiToNibInit() {
    std::call_once(g_asciiToNibFlag, initAsciiToNib);
}

// Process one lane: read all FASTQ records, hash screen, inline Solo for hits, push misses to alignQ.
// All per-thread state (localBar, sampleDet, scratch buffers) is owned by the caller
// and reused across lane claims — no per-lane allocation.
static uint64_t processOneLane(
    FlexPipelineState *st, Parameters &P, int laneId,
    gzFile gzR2, gzFile gzR1,
    SoloReadFeature *readFeat, Stats *stats,
    const std::unordered_set<uint32_t> &samplePreFilter,
    SampleDetector *sampleDet, bool sampleDetReady,
    SoloReadBarcode &localBar, bool noAlign = false)
{
    auto &cache = FlexHashScreenCache::instance();
    const uint32_t sampleTagOffset = P.pSolo.sampleProbeOffset;
    const bool hasSamplePreFilter = !samplePreFilter.empty();

    uint8_t packedScratch[4];
    char dummySeq[4] = {'\0'};
    char dummyQual[4] = {'\0'};

    char lineBuf[kFlexPipeSeqMax + 256];
    char seq0[kFlexPipeSeqMax], seq1[kFlexPipeSeqMax];
    char qual0[kFlexPipeSeqMax], qual1[kFlexPipeSeqMax];
    char name[kFlexPipeNameMax];

    uint64_t nReads = 0;

    while (true) {
        if (!gzReadLine(gzR2, lineBuf, sizeof(lineBuf))) break;
        {
            const char *src = lineBuf;
            if (*src == '@') ++src;
            size_t nameLen = 0;
            while (src[nameLen] && src[nameLen] != ' ' && src[nameLen] != '\t')
                ++nameLen;
            if (nameLen >= kFlexPipeNameMax) nameLen = kFlexPipeNameMax - 1;
            std::memcpy(name, src, nameLen);
            name[nameLen] = '\0';
        }
        if (!gzReadLine(gzR2, seq0, kFlexPipeSeqMax)) break;
        uint32_t readLen0 = static_cast<uint32_t>(strlen(seq0));
        if (!gzReadLine(gzR2, lineBuf, sizeof(lineBuf))) break;
        if (!gzReadLine(gzR2, qual0, kFlexPipeSeqMax)) break;

        if (!gzReadLine(gzR1, lineBuf, sizeof(lineBuf))) break;
        if (!gzReadLine(gzR1, seq1, kFlexPipeSeqMax)) break;
        uint32_t readLen1 = static_cast<uint32_t>(strlen(seq1));
        if (!gzReadLine(gzR1, lineBuf, sizeof(lineBuf))) break;
        if (!gzReadLine(gzR1, qual1, kFlexPipeSeqMax)) break;

        uint64_t iReadAll = st->iReadAllGlobal.fetch_add(1);
        st->counters.perLaneReads[laneId]++;
        nReads++;

        bool sampleOK = true;
        if (hasSamplePreFilter && sampleTagOffset + 8 <= readLen0) {
            uint32_t tagCode = 0;
            bool encodable = true;
            for (int i = 0; i < 8; ++i) {
                uint8_t nib = g_asciiToNib[static_cast<uint8_t>(seq0[sampleTagOffset + i])];
                if (nib == 0) { encodable = false; break; }
                tagCode = (tagCode << 4) | nib;
            }
            if (!encodable || samplePreFilter.find(tagCode) == samplePreFilter.end()) {
                sampleOK = false;
            }
        }

        FlexHashScreenDecision decision;
        if (sampleOK) {
            decision = cache.classifyReadH0H1Offset0(seq0, readLen0);
        } else {
            decision.action = FlexHashScreenDecision::Pass;
        }

        if (decision.action == FlexHashScreenDecision::Keep ||
            decision.action == FlexHashScreenDecision::Deny) {

            if (P.soloSpatialFlexIntegratedEnabled) {
                if (!completeSpatialFlexCachedRead(
                        P, seq1, readLen1, qual1, iReadAll, decision)) {
                    return nReads;
                }
                if (decision.action == FlexHashScreenDecision::Keep) {
                    stats->hashScreenKeep++;
                    st->counters.triageKeep.fetch_add(1);
                } else {
                    stats->hashScreenDeny++;
                    st->counters.triageDeny.fetch_add(1);
                }
                st->counters.readsTotal.fetch_add(1);
                continue;
            }

            char *readSeqPtrs[2]  = { dummySeq, seq1 };
            char *readQualPtrs[2] = { dummyQual, qual1 };
            uint64 readLens[2]    = { 0, readLen1 };
            std::string readNameExtra;

            localBar.getCBandUMI(readSeqPtrs, readQualPtrs, readLens, readNameExtra,
                                  static_cast<uint32_t>(laneId), name);

            uint8_t detectedSampleToken = 0xFF;
            if (sampleDetReady && sampleTagOffset + 8 <= readLen0) {
                nuclPackBAM(seq0 + sampleTagOffset, reinterpret_cast<char *>(packedScratch), 8);
                uint32_t detIdx = sampleDet->detectSampleFromPackedTag(packedScratch);
                if (detIdx > 0) {
                    detectedSampleToken = static_cast<uint8_t>(detIdx & 0x1Fu);
                }
            }
            localBar.detectedSampleToken = detectedSampleToken;

            if (decision.action == FlexHashScreenDecision::Keep) {
                record_flex_hash_screen_keep(readFeat, localBar, iReadAll,
                                             decision.geneIdx15, decision.cacheClass,
                                             decision.probeRegion);
                stats->hashScreenKeep++;
                if (localBar.cbMatch < 0) stats->hashScreenKeepNoBarcode++;
                st->counters.triageKeep.fetch_add(1);
            } else {
                record_flex_hash_screen_deny(readFeat, localBar, iReadAll, "NEG_PROBE_AMBIG");
                stats->hashScreenDeny++;
                st->counters.triageDeny.fetch_add(1);
            }
        } else {
            st->counters.triageMiss.fetch_add(1);
            if (!noAlign) {
                EnrichedPacket ep;
                std::memcpy(ep.name, name, kFlexPipeNameMax);
                std::memcpy(ep.seq[0], seq0, readLen0 + 1);
                std::memcpy(ep.seq[1], seq1, readLen1 + 1);
                std::memcpy(ep.qual[0], qual0, readLen0 + 1);
                std::memcpy(ep.qual[1], qual1, readLen1 + 1);
                ep.readLen[0] = readLen0;
                ep.readLen[1] = readLen1;
                ep.iReadAll = iReadAll;
                ep.laneId = static_cast<uint8_t>(laneId);
                ep.readFilesIndex = static_cast<uint32_t>(laneId);
                ep.readFilter = 'Y';
                ep.eof = false;
                ep.cbMatch = -1;
                ep.cbMatchIndN = 0;
                ep.umiB = 0;
                ep.detectedSampleToken = 0xFF;
                ep.hashScreenSampleIdx = 0;
                st->alignQ.push(std::move(ep));
            }
        }

        st->counters.readsTotal.fetch_add(1);
    }

    gzclose(gzR2);
    gzclose(gzR1);

    return nReads;
}

static uint64_t processCbqModuleRecords(
    FlexPipelineState *st, Parameters &P, int laneId,
    const std::string &cbqPath,
    star::input::CbqInputModule &module,
    SoloReadFeature *readFeat, Stats *stats,
    const std::unordered_set<uint32_t> &samplePreFilter,
    SampleDetector *sampleDet, bool sampleDetReady,
    SoloReadBarcode &localBar, bool noAlign,
    bool deterministicReadIds,
    uint64_t globalFirst)
{
    auto &cache = FlexHashScreenCache::instance();
    const uint32_t sampleTagOffset = P.pSolo.sampleProbeOffset;
    const bool hasSamplePreFilter = !samplePreFilter.empty();

    std::string inputError;
    uint8_t packedScratch[4];
    char dummySeq[4] = {'\0'};
    char dummyQual[4] = {'\0'};

    char seq0[kFlexPipeSeqMax], seq1[kFlexPipeSeqMax];
    char qual0[kFlexPipeSeqMax], qual1[kFlexPipeSeqMax];
    char name[kFlexPipeNameMax];

    uint64_t nReads = 0;
    star::input::CbqReadBatchView batch;
    while (true) {
        inputError.clear();
        const star::input::InputStatus status = module.next_batch(&batch, &inputError);
        if (status == star::input::InputStatus::End) break;
        if (status == star::input::InputStatus::Error) {
            P.inOut->logMain << "ERROR: Flex CBQ pipeline failed while reading " << cbqPath
                             << ": " << inputError << "\n" << std::flush;
            break;
        }

        for (uint32_t i = 0; i < batch.record_count; ++i) {
            const star::input::CbqReadView &record = batch.records[i];
            if (record.segment_count < 2) {
                P.inOut->logMain << "ERROR: Flex CBQ pipeline expected paired reads in "
                                 << cbqPath << "\n" << std::flush;
                module.close();
                return nReads;
            }

            size_t readLen0Size = 0;
            size_t readLen1Size = 0;
            if (!star::input::materialize_cbq_segment_sequence_to_buffer(
                    record.segments[0], seq0, sizeof(seq0), &readLen0Size, &inputError) ||
                !star::input::materialize_cbq_segment_sequence_to_buffer(
                    record.segments[1], seq1, sizeof(seq1), &readLen1Size, &inputError)) {
                P.inOut->logMain << "ERROR: Flex CBQ sequence materialization failed in "
                                 << cbqPath << ": " << inputError << "\n" << std::flush;
                module.close();
                return nReads;
            }

            size_t qualLen0 = 0;
            size_t qualLen1 = 0;
            if (!copyCbqSpanToBuffer(record.segments[0].quality, qual0, sizeof(qual0), &qualLen0) ||
                !copyCbqSpanToBuffer(record.segments[1].quality, qual1, sizeof(qual1), &qualLen1)) {
                P.inOut->logMain << "ERROR: Flex CBQ quality length exceeds buffer in "
                                 << cbqPath << "\n" << std::flush;
                module.close();
                return nReads;
            }
            const uint32_t readLen0 = static_cast<uint32_t>(readLen0Size);
            const uint32_t readLen1 = static_cast<uint32_t>(readLen1Size);
            copyCbqReadName(record.read_name, name, sizeof(name));

            const uint64_t localOrdinal = nReads;
            uint64_t iReadAll = deterministicReadIds
                ? globalFirst + localOrdinal
                : st->iReadAllGlobal.fetch_add(1);
            st->counters.perLaneReads[laneId]++;
            nReads++;

            bool sampleOK = true;
            if (hasSamplePreFilter && sampleTagOffset + 8 <= readLen0) {
                uint32_t tagCode = 0;
                bool encodable = true;
                for (int j = 0; j < 8; ++j) {
                    uint8_t nib = g_asciiToNib[static_cast<uint8_t>(seq0[sampleTagOffset + j])];
                    if (nib == 0) { encodable = false; break; }
                    tagCode = (tagCode << 4) | nib;
                }
                if (!encodable || samplePreFilter.find(tagCode) == samplePreFilter.end()) {
                    sampleOK = false;
                }
            }

            FlexHashScreenDecision decision;
            if (sampleOK) {
                decision = cache.classifyReadH0H1Offset0(seq0, readLen0);
            } else {
                decision.action = FlexHashScreenDecision::Pass;
            }

            if (decision.action == FlexHashScreenDecision::Keep ||
                decision.action == FlexHashScreenDecision::Deny) {

                if (P.soloSpatialFlexIntegratedEnabled) {
                    if (!completeSpatialFlexCachedRead(
                            P, seq1, readLen1, qual1, iReadAll, decision)) {
                        return nReads;
                    }
                    if (decision.action == FlexHashScreenDecision::Keep) {
                        stats->hashScreenKeep++;
                        st->counters.triageKeep.fetch_add(1);
                    } else {
                        stats->hashScreenDeny++;
                        st->counters.triageDeny.fetch_add(1);
                    }
                    st->counters.readsTotal.fetch_add(1);
                    continue;
                }

                char *readSeqPtrs[2]  = { dummySeq, seq1 };
                char *readQualPtrs[2] = { dummyQual, qual1 };
                uint64 readLens[2]    = { 0, readLen1 };
                std::string readNameExtra;

                localBar.getCBandUMI(readSeqPtrs, readQualPtrs, readLens, readNameExtra,
                                      static_cast<uint32_t>(laneId), name);

                uint8_t detectedSampleToken = 0xFF;
                if (sampleDetReady && sampleTagOffset + 8 <= readLen0) {
                    nuclPackBAM(seq0 + sampleTagOffset, reinterpret_cast<char *>(packedScratch), 8);
                    uint32_t detIdx = sampleDet->detectSampleFromPackedTag(packedScratch);
                    if (detIdx > 0) {
                        detectedSampleToken = static_cast<uint8_t>(detIdx & 0x1Fu);
                    }
                }
                localBar.detectedSampleToken = detectedSampleToken;

                if (decision.action == FlexHashScreenDecision::Keep) {
                    record_flex_hash_screen_keep(readFeat, localBar, iReadAll,
                                                 decision.geneIdx15, decision.cacheClass,
                                                 decision.probeRegion);
                    stats->hashScreenKeep++;
                    if (localBar.cbMatch < 0) stats->hashScreenKeepNoBarcode++;
                    st->counters.triageKeep.fetch_add(1);
                } else {
                    record_flex_hash_screen_deny(readFeat, localBar, iReadAll, "NEG_PROBE_AMBIG");
                    stats->hashScreenDeny++;
                    st->counters.triageDeny.fetch_add(1);
                }
            } else {
                st->counters.triageMiss.fetch_add(1);
                if (!noAlign) {
                    EnrichedPacket ep;
                    std::memcpy(ep.name, name, kFlexPipeNameMax);
                    std::memcpy(ep.seq[0], seq0, readLen0 + 1);
                    std::memcpy(ep.seq[1], seq1, readLen1 + 1);
                    std::memcpy(ep.qual[0], qual0, readLen0 + 1);
                    std::memcpy(ep.qual[1], qual1, readLen1 + 1);
                    ep.readLen[0] = readLen0;
                    ep.readLen[1] = readLen1;
                    ep.iReadAll = iReadAll;
                    ep.laneId = static_cast<uint8_t>(laneId);
                    ep.readFilesIndex = static_cast<uint32_t>(laneId);
                    ep.readFilter = record.read_filter;
                    ep.eof = false;
                    ep.cbMatch = -1;
                    ep.cbMatchIndN = 0;
                    ep.umiB = 0;
                    ep.detectedSampleToken = 0xFF;
                    ep.hashScreenSampleIdx = 0;
                    st->alignQ.push(std::move(ep));
                }
            }

            st->counters.readsTotal.fetch_add(1);
        }
    }

    return nReads;
}

static star::input::InputSourcePlan makeSingleCbqLanePlan(const std::string &cbqPath) {
    std::vector<std::vector<std::string>> readFiles(1);
    readFiles[0].push_back(cbqPath);
    return star::input::make_cbq_input_source_plan(readFiles, std::vector<std::string>(), 2);
}

static uint64_t processOneCbqLane(
    FlexPipelineState *st, Parameters &P, int laneId,
    const std::string &cbqPath,
    SoloReadFeature *readFeat, Stats *stats,
    const std::unordered_set<uint32_t> &samplePreFilter,
    SampleDetector *sampleDet, bool sampleDetReady,
    SoloReadBarcode &localBar, bool noAlign = false)
{
    const star::input::InputSourcePlan plan = makeSingleCbqLanePlan(cbqPath);
    star::input::CbqInputModule module;
    std::string inputError;
    if (!module.configure(plan, &inputError) || !module.open(&inputError)) {
        P.inOut->logMain << "ERROR: Flex CBQ pipeline could not open " << cbqPath
                         << ": " << inputError << "\n" << std::flush;
        return 0;
    }
    const uint64_t nReads = processCbqModuleRecords(
        st, P, laneId, cbqPath, module, readFeat, stats, samplePreFilter,
        sampleDet, sampleDetReady, localBar, noAlign, false, 0);
    module.close();
    return nReads;
}

static uint64_t processOneCbqRange(
    FlexPipelineState *st, Parameters &P,
    const FlexCbqRangeTask &task,
    const std::string &cbqPath,
    SoloReadFeature *readFeat, Stats *stats,
    const std::unordered_set<uint32_t> &samplePreFilter,
    SampleDetector *sampleDet, bool sampleDetReady,
    SoloReadBarcode &localBar, bool noAlign)
{
    const star::input::InputSourcePlan plan = makeSingleCbqLanePlan(cbqPath);
    star::input::CbqInputModule module;
    std::string inputError;
    if (!module.configure(plan, &inputError) ||
        !module.open_range(0, task.firstRecord, task.recordCount, &inputError)) {
        P.inOut->logMain << "ERROR: Flex CBQ range pipeline could not open "
                         << cbqPath << " lane=" << task.laneId
                         << " first=" << task.firstRecord
                         << " count=" << task.recordCount
                         << ": " << inputError << "\n" << std::flush;
        return 0;
    }
    const uint64_t nReads = processCbqModuleRecords(
        st, P, task.laneId, cbqPath, module, readFeat, stats, samplePreFilter,
        sampleDet, sampleDetReady, localBar, noAlign, true, task.globalFirst);
    module.close();
    return nReads;
}

bool flexPrepareCbqRangeTasks(FlexPipelineState *state, Parameters &P,
                              int nWorkers, std::string *reason) {
    (void)P;
    if (state == nullptr || nWorkers <= 0) {
        if (reason != nullptr) {
            *reason = "invalid range task planner inputs";
        }
        return false;
    }
    state->cbqRangeTasks.clear();
    state->nextCbqRangeIdx.store(0, std::memory_order_relaxed);

    std::vector<uint64_t> laneCounts(static_cast<size_t>(state->nLanes), 0);
    uint64_t totalRecords = 0;
    for (int lane = 0; lane < state->nLanes; ++lane) {
        const std::string &cbqPath = state->laneFiles[static_cast<size_t>(lane)].r2path;
        const star::input::InputSourcePlan plan = makeSingleCbqLanePlan(cbqPath);
        star::input::CbqInputModule module;
        std::string inputError;
        if (!module.configure(plan, &inputError) ||
            !module.open_range(0, 0, std::numeric_limits<uint64_t>::max(), &inputError)) {
            if (reason != nullptr) {
                *reason = "lane " + std::to_string(lane) + ": " + inputError;
            }
            return false;
        }
        laneCounts[static_cast<size_t>(lane)] = module.current_lane_record_count();
        totalRecords += laneCounts[static_cast<size_t>(lane)];
        module.close();
    }

    if (totalRecords == 0) {
        if (reason != nullptr) {
            *reason = "no CBQ records";
        }
        return false;
    }

    const uint64_t chunkSize =
        (totalRecords + static_cast<uint64_t>(nWorkers) - 1U) /
        static_cast<uint64_t>(nWorkers);
    for (int worker = 0; worker < nWorkers; ++worker) {
        const uint64_t globalFirst = static_cast<uint64_t>(worker) * chunkSize;
        if (globalFirst >= totalRecords) {
            break;
        }
        const uint64_t globalEnd = std::min(totalRecords, globalFirst + chunkSize);
        uint64_t laneGlobalFirst = 0;
        for (int lane = 0; lane < state->nLanes; ++lane) {
            const uint64_t laneCount = laneCounts[static_cast<size_t>(lane)];
            const uint64_t laneGlobalEnd = laneGlobalFirst + laneCount;
            if (laneCount > 0 && globalFirst < laneGlobalEnd && globalEnd > laneGlobalFirst) {
                const uint64_t overlapFirst = std::max(globalFirst, laneGlobalFirst);
                const uint64_t overlapEnd = std::min(globalEnd, laneGlobalEnd);
                FlexCbqRangeTask task;
                task.laneId = lane;
                task.firstRecord = overlapFirst - laneGlobalFirst;
                task.recordCount = overlapEnd - overlapFirst;
                task.globalFirst = overlapFirst;
                state->cbqRangeTasks.push_back(task);
            }
            laneGlobalFirst = laneGlobalEnd;
            if (laneGlobalFirst >= globalEnd) {
                break;
            }
        }
    }

    if (reason != nullptr) {
        std::ostringstream out;
        out << state->cbqRangeTasks.size() << " ranges across "
            << state->nLanes << " lanes and " << totalRecords << " records";
        *reason = out.str();
    }
    return !state->cbqRangeTasks.empty();
}

void *flexLaneReaderFullThread(void *arg) {
    FlexLaneReaderArgs *ctx = static_cast<FlexLaneReaderArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    Parameters &P = *ctx->P;
    SoloReadFeature *readFeat = ctx->readFeat;
    Stats *stats = ctx->stats;
    ReadAlign *RA = ctx->RA;

    ensureAsciiToNibInit();
    std::call_once(g_samplePreFilterFlag, buildSamplePreFilter, std::ref(P));
    const auto &samplePreFilter = g_samplePreFilter;
    const bool noAlign = (P.pSolo.flexNoAlign != 0);

    // Preallocate once per thread — reused across all lane claims
    SoloReadBarcode localBar(P);
    SampleDetector *sampleDet = nullptr;
    bool sampleDetReady = false;

    auto ensureSampleDetector = [&]() {
        if (sampleDet == nullptr &&
            !P.pSolo.sampleWhitelistPath.empty() && P.pSolo.sampleWhitelistPath != "-" &&
            !P.pSolo.sampleProbesPath.empty() && P.pSolo.sampleProbesPath != "-") {
            sampleDet = new SampleDetector(P.pSolo);
            if (sampleDet->loadWhitelist(P.pSolo.sampleWhitelistPath) &&
                sampleDet->loadProbes(P.pSolo.sampleProbesPath)) {
                sampleDetReady = true;
            } else {
                delete sampleDet;
                sampleDet = nullptr;
            }
        }
    };

    const bool useCbqRangeTasks = (P.readFilesTypeN == 20 && P.cbqInputActive &&
                                   !st->cbqRangeTasks.empty());
    if (useCbqRangeTasks) {
        FlexCbqRangeTask task;
        while (st->claimNextCbqRange(&task)) {
            ensureSampleDetector();
            const std::string &cbqPath =
                st->laneFiles[static_cast<size_t>(task.laneId)].r2path;
            processOneCbqRange(st, P, task, cbqPath, readFeat, stats,
                               samplePreFilter, sampleDet, sampleDetReady,
                               localBar, noAlign);
        }
    } else {
        // Phase 1: Lane work-stealing loop — claim and process lanes until none remain
        int lane;
        while ((lane = st->claimNextLane()) >= 0) {
            // Lazy-init SampleDetector on first lane claim (threads that never
            // claim a lane skip the file I/O entirely)
            ensureSampleDetector();

            const std::string &r2path = st->laneFiles[lane].r2path;
            if (P.readFilesTypeN == 20 && P.cbqInputActive) {
                processOneCbqLane(st, P, lane, r2path, readFeat, stats,
                                  samplePreFilter, sampleDet, sampleDetReady, localBar, noAlign);
            } else {
                const std::string &r1path = st->laneFiles[lane].r1path;
                gzFile gzR2 = gzopen(r2path.c_str(), "rb");
                gzFile gzR1 = gzopen(r1path.c_str(), "rb");
                if (!gzR2 || !gzR1) {
                    if (gzR2) gzclose(gzR2);
                    if (gzR1) gzclose(gzR1);
                    continue;
                }
                gzbuffer(gzR2, kGzBufSize);
                gzbuffer(gzR1, kGzBufSize);

                processOneLane(st, P, lane, gzR2, gzR1, readFeat, stats,
                               samplePreFilter, sampleDet, sampleDetReady, localBar, noAlign);
            }
        }
    }

    // Signal that this reader is done; last one closes alignQ
    int finishedCount = st->readersFinished.fetch_add(1) + 1;
    if (finishedCount == st->nFusedThreads) {
        st->alignQ.close();
    }

    // Phase 2: Role switch — drain alignQ as an alignment worker (skipped in noAlign mode)
    if (!noAlign && RA != nullptr) {
        EnrichedPacket ep;
        while (st->alignQ.pop(ep)) {
            RA->oneReadFromPacket(ep);
        }
    }

    delete sampleDet;
    return nullptr;
}

void *flexTriageThread(void *arg) {
    FlexTriageArgs *ctx = static_cast<FlexTriageArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    Parameters &P = *ctx->P;

    auto &cache = FlexHashScreenCache::instance();
    const uint32_t sampleTagOffset = P.pSolo.sampleProbeOffset;

    ReadPacket rpkt;
    while (st->readerQ.pop(rpkt)) {

        FlexHashScreenDecision decision = P.soloSpatialFlexIntegratedEnabled
            ? cache.classifyReadH0H1Offset0(rpkt.seq[0], rpkt.readLen[0])
            : cache.classifyReadH0Offset0(rpkt.seq[0], rpkt.readLen[0]);

        if (decision.action == FlexHashScreenDecision::Keep ||
            decision.action == FlexHashScreenDecision::Deny) {
            DecisionPacket dp;
            dp.iReadAll = rpkt.iReadAll;
            dp.readFilesIndex = rpkt.readFilesIndex;
            dp.eof = false;

            // Carry raw R1 (barcode read) for deferred CB/UMI extraction
            const uint32_t bcLen = std::min(rpkt.readLen[1], kFlexPipeBarcodeSeqMax - 1);
            std::memcpy(dp.barcodeSeq, rpkt.seq[1], bcLen);
            dp.barcodeSeq[bcLen] = '\0';
            std::memcpy(dp.barcodeQual, rpkt.qual[1], bcLen);
            dp.barcodeQual[bcLen] = '\0';
            dp.barcodeLen = bcLen;
            std::memcpy(dp.readName, rpkt.name, kFlexPipeNameMax);

            // Carry 8-byte sample tag region from R2
            if (sampleTagOffset + 8 <= rpkt.readLen[0]) {
                std::memcpy(dp.sampleTagSeq, rpkt.seq[0] + sampleTagOffset, 8);
                dp.sampleTagLen = 8;
            } else {
                dp.sampleTagLen = 0;
            }

            if (decision.action == FlexHashScreenDecision::Keep) {
                dp.verdict = DecisionPacket::KEEP;
                dp.geneIdx15 = decision.geneIdx15;
                dp.cacheClass = decision.cacheClass;
                dp.probeRegion = decision.probeRegion;
                dp.denyReason = nullptr;
                st->counters.triageKeep.fetch_add(1);
            } else {
                dp.verdict = DecisionPacket::DENY;
                dp.geneIdx15 = 0;
                dp.cacheClass = 0;
                dp.denyReason = "NEG_PROBE_AMBIG";
                st->counters.triageDeny.fetch_add(1);
            }

            int shard = static_cast<int>(rpkt.iReadAll % static_cast<uint64_t>(st->nSolo));
            st->soloQ[shard]->push(std::move(dp));
        } else {
            // MISS — route to alignment queue (alignment workers do their own CB/UMI)
            EnrichedPacket ep;
            std::memcpy(ep.name, rpkt.name, kFlexPipeNameMax);
            std::memcpy(ep.seq[0], rpkt.seq[0], rpkt.readLen[0] + 1);
            std::memcpy(ep.seq[1], rpkt.seq[1], rpkt.readLen[1] + 1);
            std::memcpy(ep.qual[0], rpkt.qual[0], rpkt.readLen[0] + 1);
            std::memcpy(ep.qual[1], rpkt.qual[1], rpkt.readLen[1] + 1);
            ep.readLen[0] = rpkt.readLen[0];
            ep.readLen[1] = rpkt.readLen[1];
            ep.iReadAll = rpkt.iReadAll;
            ep.laneId = rpkt.laneId;
            ep.readFilesIndex = rpkt.readFilesIndex;
            ep.readFilter = rpkt.readFilter;
            ep.eof = false;

            ep.cbMatch = -1;
            ep.cbMatchIndN = 0;
            ep.umiB = 0;
            ep.detectedSampleToken = 0xFF;
            ep.hashScreenSampleIdx = 0;

            st->counters.triageMiss.fetch_add(1);
            st->alignQ.push(std::move(ep));
        }

        st->counters.readsTotal.fetch_add(1);
    }

    // Last triage thread to finish closes downstream queues
    int triageDone = st->triageFinished.fetch_add(1) + 1;
    if (triageDone == st->nTriage) {
        for (int i = 0; i < st->nSolo; ++i) {
            st->soloQ[i]->close();
        }
        st->alignQ.close();
    }

    return nullptr;
}

void *flexSoloConsumerThread(void *arg) {
    FlexSoloConsumerArgs *ctx = static_cast<FlexSoloConsumerArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    Parameters &P = *ctx->P;
    SoloReadFeature *readFeat = ctx->readFeat;
    Stats *stats = ctx->stats;
    const int consumerId = ctx->consumerId;

    SoloReadBarcode localBar(P);

    // Solo consumer owns sample detection (deferred from triage)
    SampleDetector *sampleDet = nullptr;
    bool sampleDetReady = false;
    if (!P.pSolo.sampleWhitelistPath.empty() && P.pSolo.sampleWhitelistPath != "-" &&
        !P.pSolo.sampleProbesPath.empty() && P.pSolo.sampleProbesPath != "-") {
        sampleDet = new SampleDetector(P.pSolo);
        if (sampleDet->loadWhitelist(P.pSolo.sampleWhitelistPath) &&
            sampleDet->loadProbes(P.pSolo.sampleProbesPath)) {
            sampleDetReady = true;
        } else {
            delete sampleDet;
            sampleDet = nullptr;
        }
    }

    uint8_t packedScratch[4];  // 8 bases in 4-bit BAM = 4 bytes

    // Scratch arrays for getCBandUMI — R2 not needed, provide dummy
    char dummySeq[4] = {'\0'};
    char dummyQual[4] = {'\0'};

    DecisionPacket dp;
    while (st->soloQ[consumerId]->pop(dp)) {

        if (P.soloSpatialFlexIntegratedEnabled) {
            FlexHashScreenDecision decision;
            decision.action = dp.verdict == DecisionPacket::KEEP
                ? FlexHashScreenDecision::Keep
                : FlexHashScreenDecision::Deny;
            decision.geneIdx15 = dp.geneIdx15;
            decision.cacheClass = dp.cacheClass;
            decision.probeRegion = dp.probeRegion;
            if (!completeSpatialFlexCachedRead(
                    P, dp.barcodeSeq, dp.barcodeLen, dp.barcodeQual,
                    dp.iReadAll, decision)) {
                break;
            }
            if (dp.verdict == DecisionPacket::KEEP) {
                stats->hashScreenKeep++;
            } else {
                stats->hashScreenDeny++;
            }
            continue;
        }

        // Deferred CB/UMI extraction from raw R1 carried in packet
        char *readSeqPtrs[2]  = { dummySeq, dp.barcodeSeq };
        char *readQualPtrs[2] = { dummyQual, dp.barcodeQual };
        uint64 readLens[2]    = { 0, dp.barcodeLen };
        std::string readNameExtra;

        localBar.getCBandUMI(readSeqPtrs, readQualPtrs, readLens, readNameExtra,
                              dp.readFilesIndex, dp.readName);

        // Deferred sample detection from the 8-byte raw tag carried in the packet
        uint8_t detectedSampleToken = 0xFF;
        if (sampleDetReady && dp.sampleTagLen == 8) {
            nuclPackBAM(dp.sampleTagSeq, reinterpret_cast<char *>(packedScratch), 8);
            uint32_t detIdx = sampleDet->detectSampleFromPackedTag(packedScratch);
            if (detIdx > 0) {
                detectedSampleToken = static_cast<uint8_t>(detIdx & 0x1Fu);
            }
        }
        localBar.detectedSampleToken = detectedSampleToken;

        if (dp.verdict == DecisionPacket::KEEP) {
            record_flex_hash_screen_keep(readFeat, localBar, dp.iReadAll,
                                         dp.geneIdx15, dp.cacheClass, dp.probeRegion);
            stats->hashScreenKeep++;
            if (localBar.cbMatch < 0) stats->hashScreenKeepNoBarcode++;
        } else {
            record_flex_hash_screen_deny(readFeat, localBar, dp.iReadAll, dp.denyReason);
            stats->hashScreenDeny++;
        }
    }

    delete sampleDet;
    return nullptr;
}

void *flexAlignWorkerThread(void *arg) {
    FlexAlignWorkerArgs *ctx = static_cast<FlexAlignWorkerArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    ReadAlign *RA = ctx->RA;

    EnrichedPacket ep;
    while (st->alignQ.pop(ep)) {
        RA->oneReadFromPacket(ep);
    }

    return nullptr;
}

void *flexStatsReporterThread(void *arg) {
    FlexStatsReporterArgs *ctx = static_cast<FlexStatsReporterArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    Parameters &P = *ctx->P;

    auto tStart = std::chrono::steady_clock::now();
    uint64_t prevTotal = 0;

    while (!st->pipelineDone.load(std::memory_order_relaxed)) {
        sleep(10);
        if (st->pipelineDone.load(std::memory_order_relaxed)) break;

        auto tNow = std::chrono::steady_clock::now();
        int elapsedSec = static_cast<int>(
            std::chrono::duration_cast<std::chrono::seconds>(tNow - tStart).count());

        uint64_t total = st->counters.readsTotal.load(std::memory_order_relaxed);
        uint64_t keep  = st->counters.triageKeep.load(std::memory_order_relaxed);
        uint64_t deny  = st->counters.triageDeny.load(std::memory_order_relaxed);
        uint64_t miss  = st->counters.triageMiss.load(std::memory_order_relaxed);

        size_t rqSize = st->readerQ.size();
        size_t aqSize = st->alignQ.size();

        std::ostringstream soloQStr;
        if (st->nSolo > 0) {
            for (int i = 0; i < st->nSolo; ++i) {
                if (i > 0) soloQStr << ",";
                soloQStr << st->soloQ[i]->size();
            }
        } else {
            soloQStr << "n/a";
        }

        uint64_t delta = total - prevTotal;
        prevTotal = total;

        P.inOut->logMain << "PIPELINE_STATS t=" << elapsedSec << "s"
                         << " total=" << total
                         << " keep=" << keep
                         << " deny=" << deny
                         << " miss=" << miss
                         << " readerQ=" << rqSize
                         << " soloQ=" << soloQStr.str()
                         << " alignQ=" << aqSize
                         << " delta10s=" << delta
                         << "\n" << std::flush;
    }

    return nullptr;
}
