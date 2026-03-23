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
#include "Stats.h"

#include <cstring>
#include <cstdlib>
#include <fstream>
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

        FlexHashScreenDecision decision = cache.classifyReadH0Offset0(seq0, readLen0);

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

void *flexLaneReaderFullThread(void *arg) {
    FlexLaneReaderArgs *ctx = static_cast<FlexLaneReaderArgs *>(arg);
    FlexPipelineState *st = ctx->state;
    Parameters &P = *ctx->P;
    const int laneId = ctx->laneId;
    gzFile gzR2 = ctx->gzR2;
    gzFile gzR1 = ctx->gzR1;
    SoloReadFeature *readFeat = ctx->readFeat;
    Stats *stats = ctx->stats;

    auto &cache = FlexHashScreenCache::instance();
    const uint32_t sampleTagOffset = P.pSolo.sampleProbeOffset;

    SoloReadBarcode localBar(P);

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

    // Build sample code set for fast pre-filter: all valid variant codes
    // across user-selected samples (nSamples * up to 8 variants each).
    std::unordered_set<uint32_t> validSampleCodes;
    if (sampleDetReady && sampleDet != nullptr) {
        // Re-load whitelist+probes to get sampleCodes — use detectSampleFromPackedTag
        // approach: pre-enumerate all valid codes by encoding each sample's variants.
        // The SampleDetector stores sampleCodes_ privately, so we build the set by
        // encoding the canonical + variants from the probe file.
        // Simpler approach: iterate all 16 possible nibble-encoded 8-mers for each
        // sample. But that's intractable. Instead, load the probe file directly.
    }

    // ASCII-to-nibble LUT for fast 8-byte sample tag encoding
    static uint8_t asciiToNib[256] = {};
    static bool nibLutInit = false;
    if (!nibLutInit) {
        std::memset(asciiToNib, 0, sizeof(asciiToNib));
        asciiToNib['A'] = asciiToNib['a'] = 1;
        asciiToNib['C'] = asciiToNib['c'] = 2;
        asciiToNib['G'] = asciiToNib['g'] = 4;
        asciiToNib['T'] = asciiToNib['t'] = 8;
        nibLutInit = true;
    }

    // Build fast sample pre-filter: hash set of all valid uint32 sample codes.
    // Encode each variant from the probe file the same way detectSampleFromPackedTag does.
    std::unordered_set<uint32_t> samplePreFilter;
    if (sampleDetReady) {
        // Re-parse the probe file to extract all variant sequences
        std::ifstream probeFile(P.pSolo.sampleProbesPath.c_str());
        if (probeFile.is_open()) {
            // Also load the whitelist canonicals into a set
            std::unordered_set<std::string> userCanonicals;
            {
                std::ifstream wl(P.pSolo.sampleWhitelistPath.c_str());
                std::string line;
                while (std::getline(wl, line)) {
                    std::istringstream iss(line);
                    std::string id, canonical;
                    if (iss >> id >> canonical) {
                        userCanonicals.insert(canonical);
                    }
                }
            }
            // Add canonical sequences
            for (const std::string &canon : userCanonicals) {
                if (canon.size() == 8) {
                    uint32_t code = 0;
                    bool ok = true;
                    for (int i = 0; i < 8; ++i) {
                        uint8_t nib = asciiToNib[static_cast<uint8_t>(canon[i])];
                        if (nib == 0) { ok = false; break; }
                        code = (code << 4) | nib;
                    }
                    if (ok) samplePreFilter.insert(code);
                }
            }
            // Add variants from probe file
            std::string line;
            while (std::getline(probeFile, line)) {
                std::istringstream iss(line);
                std::string variant, canonical, barcodeId;
                if (!(iss >> variant >> canonical >> barcodeId)) continue;
                if (variant.size() != 8 || canonical.size() != 8) continue;
                if (userCanonicals.find(canonical) == userCanonicals.end()) continue;
                uint32_t code = 0;
                bool ok = true;
                for (int i = 0; i < 8; ++i) {
                    uint8_t nib = asciiToNib[static_cast<uint8_t>(variant[i])];
                    if (nib == 0) { ok = false; break; }
                    code = (code << 4) | nib;
                }
                if (ok) samplePreFilter.insert(code);
            }
        }
    }
    const bool hasSamplePreFilter = !samplePreFilter.empty();

    uint8_t packedScratch[4];
    char dummySeq[4] = {'\0'};
    char dummyQual[4] = {'\0'};

    char lineBuf[kFlexPipeSeqMax + 256];
    char seq0[kFlexPipeSeqMax], seq1[kFlexPipeSeqMax];
    char qual0[kFlexPipeSeqMax], qual1[kFlexPipeSeqMax];
    char name[kFlexPipeNameMax];

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

        // Sample pre-filter: encode 8bp at sampleTagOffset, skip hash screen
        // if the sample tag doesn't match any user-selected sample variant.
        bool sampleOK = true;
        if (hasSamplePreFilter && sampleTagOffset + 8 <= readLen0) {
            uint32_t tagCode = 0;
            bool encodable = true;
            for (int i = 0; i < 8; ++i) {
                uint8_t nib = asciiToNib[static_cast<uint8_t>(seq0[sampleTagOffset + i])];
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
                                             decision.geneIdx15, decision.cacheClass);
                stats->hashScreenKeep++;
                st->counters.triageKeep.fetch_add(1);
            } else {
                record_flex_hash_screen_deny(readFeat, localBar, iReadAll, "NEG_PROBE_AMBIG");
                stats->hashScreenDeny++;
                st->counters.triageDeny.fetch_add(1);
            }
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
    delete sampleDet;

    int finishedCount = st->readersFinished.fetch_add(1) + 1;
    if (finishedCount == st->nLanes) {
        st->alignQ.close();
    }

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

        FlexHashScreenDecision decision = cache.classifyReadH0Offset0(
            rpkt.seq[0], rpkt.readLen[0]);

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
                                         dp.geneIdx15, dp.cacheClass);
            stats->hashScreenKeep++;
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
